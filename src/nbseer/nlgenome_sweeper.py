#!/usr/bin/env python3
"""
NLGenomeSweeper集成模块

该模块提供了NLGenomeSweeper工具的Python封装，用于在基因组中识别NBS-LRR基因。
NLGenomeSweeper基于NB-ARC结构域的存在来搜索NBS-LRR (NLR) 抗病基因。

主要功能：
1. 基因组中NBS-LRR基因的识别
2. 基于HMM模型的NB-ARC结构域搜索
3. 两轮搜索策略（通用+物种特异性）
4. 候选基因位点的BED格式输出
5. 并行处理大基因组支持

作者: NBS基因注释项目组
创建时间: 2025-06-21
"""

import os
import subprocess
import shutil
from pathlib import Path
from typing import Optional, List, Dict, Tuple, Union
from dataclasses import dataclass
import tempfile
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
import pandas as pd
from Bio import SeqIO

from ..utils.config import get_config
from ..utils.logging_setup import get_logger


@dataclass
class NLGSearchParams:
    """NLGenomeSweeper搜索参数"""
    genome_file: str
    output_dir: str
    consensus_file: Optional[str] = None
    threads: int = 4
    overwrite: bool = False
    max_intron_size: int = 1000
    evalue_threshold: float = 1e-5
    
    def __post_init__(self):
        """验证参数"""
        if not os.path.exists(self.genome_file):
            raise FileNotFoundError(f"基因组文件不存在: {self.genome_file}")
        
        if self.consensus_file and not os.path.exists(self.consensus_file):
            raise FileNotFoundError(f"共识序列文件不存在: {self.consensus_file}")
        
        if self.threads < 1:
            raise ValueError("线程数必须大于0")
        
        if self.evalue_threshold <= 0:
            raise ValueError("E-value阈值必须大于0")


@dataclass
class CandidateGene:
    """候选基因信息"""
    chromosome: str
    start: int
    end: int
    strand: str = "+"
    score: float = 0.0
    gene_type: str = "unknown"
    domain_info: Optional[Dict] = None
    
    def to_bed_line(self) -> str:
        """转换为BED格式行"""
        return f"{self.chromosome}\t{self.start}\t{self.end}\t{self.gene_type}\t{self.score}\t{self.strand}"
    
    @property
    def length(self) -> int:
        """基因长度"""
        return self.end - self.start
    
    @property
    def coordinates(self) -> str:
        """基因坐标字符串"""
        return f"{self.chromosome}:{self.start}-{self.end}"


@dataclass
class NLGSearchResult:
    """NLGenomeSweeper搜索结果"""
    candidates: List[CandidateGene]
    output_files: Dict[str, str]
    search_params: NLGSearchParams
    execution_time: float
    total_candidates: int
    
    def to_bed_file(self, output_path: str) -> None:
        """导出为BED格式文件"""
        with open(output_path, 'w') as f:
            f.write("# NLGenomeSweeper候选基因结果\n")
            f.write("# 格式: chromosome\tstart\tend\tgene_type\tscore\tstrand\n")
            for candidate in self.candidates:
                f.write(candidate.to_bed_line() + "\n")
    
    def get_summary(self) -> Dict:
        """获取结果摘要"""
        return {
            "total_candidates": self.total_candidates,
            "execution_time": self.execution_time,
            "genome_file": self.search_params.genome_file,
            "output_directory": self.search_params.output_dir,
            "threads_used": self.search_params.threads,
            "candidate_types": self._get_type_distribution()
        }
    
    def _get_type_distribution(self) -> Dict[str, int]:
        """获取候选基因类型分布"""
        type_counts = {}
        for candidate in self.candidates:
            gene_type = candidate.gene_type
            type_counts[gene_type] = type_counts.get(gene_type, 0) + 1
        return type_counts


class NLGenomeSweeperRunner:
    """
    NLGenomeSweeper工具的Python封装类
    
    该类封装了NLGenomeSweeper命令行工具的调用，提供了Python接口
    用于在基因组中识别NBS-LRR基因。
    
    主要功能：
    - 配置和执行NLGenomeSweeper搜索
    - 解析输出结果
    - 并行处理支持
    - 结果格式转换
    """
    
    def __init__(self, config_path: Optional[str] = None):
        """
        初始化NLGenomeSweeperRunner
        
        Args:
            config_path: 配置文件路径，如果为None则使用默认配置
        """
        self.config = get_config(config_path)
        self.logger = get_logger(__name__)
        
        # 获取工具路径
        self.nlg_executable = self._get_nlg_executable()
        self.tool_config = self.config.get_tool_config('nlgenome_sweeper')
        
        # 验证工具可用性
        self._validate_dependencies()
        
        self.logger.info("NLGenomeSweeperRunner初始化完成")
    
    def _get_nlg_executable(self) -> str:
        """获取NLGenomeSweeper可执行文件路径"""
        # 从配置中获取路径
        nlg_path = self.config.get("tools.nlgenome_sweeper.executable_path")
        
        if nlg_path and os.path.exists(nlg_path):
            return nlg_path
        
        # 尝试在工具目录中查找
        tools_dir = self.config.get("paths.tools_dir")
        nlg_tool_path = os.path.join(tools_dir, "NLGenomeSweeper", "NLGenomeSweeper")
        
        if os.path.exists(nlg_tool_path):
            return nlg_tool_path
        
        # 尝试在PATH中查找
        nlg_in_path = shutil.which("NLGenomeSweeper")
        if nlg_in_path:
            return nlg_in_path
        
        raise FileNotFoundError(
            "无法找到NLGenomeSweeper可执行文件。请检查配置或安装。"
        )
    
    def _validate_dependencies(self) -> None:
        """验证所需依赖工具是否可用"""
        required_tools = [
            'samtools', 'bedtools', 'blastp', 'blastx', 
            'TransDecoder.LongOrfs', 'TransDecoder.Predict',
            'muscle', 'interproscan', 'hmmbuild'
        ]
        
        missing_tools = []
        for tool in required_tools:
            if not shutil.which(tool):
                missing_tools.append(tool)
        
        if missing_tools:
            self.logger.warning(f"缺少以下工具: {', '.join(missing_tools)}")
            self.logger.warning("NLGenomeSweeper可能无法正常运行")
        else:
            self.logger.info("所有依赖工具检查通过")
    
    def run_search(self, params: NLGSearchParams) -> NLGSearchResult:
        """
        执行NLGenomeSweeper搜索
        
        Args:
            params: 搜索参数
            
        Returns:
            NLGSearchResult: 搜索结果对象
        """
        import time
        start_time = time.time()
        
        self.logger.info(f"开始NLGenomeSweeper搜索: {params.genome_file}")
        
        # 准备输出目录
        self._prepare_output_directory(params)
        
        # 构建命令
        cmd = self._build_command(params)
        
        # 执行搜索
        try:
            self._execute_command(cmd, params.output_dir)
            
            # 解析结果
            candidates = self._parse_results(params.output_dir)
            
            # 获取输出文件
            output_files = self._get_output_files(params.output_dir)
            
            execution_time = time.time() - start_time
            
            result = NLGSearchResult(
                candidates=candidates,
                output_files=output_files,
                search_params=params,
                execution_time=execution_time,
                total_candidates=len(candidates)
            )
            
            self.logger.info(f"搜索完成，找到 {len(candidates)} 个候选基因")
            return result
            
        except Exception as e:
            self.logger.error(f"NLGenomeSweeper搜索失败: {e}")
            raise
    
    def _prepare_output_directory(self, params: NLGSearchParams) -> None:
        """准备输出目录"""
        output_path = Path(params.output_dir)
        
        if output_path.exists():
            if params.overwrite:
                self.logger.info(f"清理现有输出目录: {output_path}")
                shutil.rmtree(output_path)
            else:
                raise FileExistsError(f"输出目录已存在: {output_path}")
        
        output_path.mkdir(parents=True, exist_ok=True)
        self.logger.info(f"创建输出目录: {output_path}")
    
    def _build_command(self, params: NLGSearchParams) -> List[str]:
        """构建NLGenomeSweeper命令"""
        cmd = [
            self.nlg_executable,
            "-genome", params.genome_file,
            "-outdir", str(Path(params.output_dir).parent),
            "-t", str(params.threads)
        ]
        
        if params.consensus_file:
            cmd.extend(["-consensus", params.consensus_file])
        
        if params.overwrite:
            cmd.extend(["-overwrite", "T"])
        else:
            cmd.extend(["-overwrite", "F"])
        
        self.logger.debug(f"构建命令: {' '.join(cmd)}")
        return cmd
    
    def _execute_command(self, cmd: List[str], output_dir: str) -> None:
        """执行NLGenomeSweeper命令"""
        try:
            # 设置环境变量
            env = os.environ.copy()
            
            # 执行命令
            result = subprocess.run(
                cmd,
                cwd=output_dir,
                capture_output=True,
                text=True,
                check=True,
                env=env
            )
            
            # 记录输出
            if result.stdout:
                self.logger.debug(f"NLGenomeSweeper stdout: {result.stdout}")
            if result.stderr:
                self.logger.debug(f"NLGenomeSweeper stderr: {result.stderr}")
                
        except subprocess.CalledProcessError as e:
            self.logger.error(f"NLGenomeSweeper执行失败: {e}")
            self.logger.error(f"标准错误输出: {e.stderr}")
            raise RuntimeError(f"NLGenomeSweeper执行失败: {e.stderr}")
    
    def _parse_results(self, output_dir: str) -> List[CandidateGene]:
        """解析NLGenomeSweeper结果"""
        candidates = []
        
        # 查找主要结果文件
        nlg_output_dir = os.path.join(output_dir, "NLGenomeSweeper")
        bed_files = [
            "All_candidates.bed",
            "Final_candidates.bed"
        ]
        
        for bed_file in bed_files:
            bed_path = os.path.join(nlg_output_dir, bed_file)
            if os.path.exists(bed_path):
                candidates.extend(self._parse_bed_file(bed_path))
                break
        
        return candidates
    
    def _parse_bed_file(self, bed_path: str) -> List[CandidateGene]:
        """解析BED格式文件"""
        candidates = []
        
        try:
            with open(bed_path, 'r') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    
                    # 跳过注释和空行
                    if line.startswith('#') or not line:
                        continue
                    
                    try:
                        parts = line.split('\t')
                        if len(parts) >= 3:
                            candidate = CandidateGene(
                                chromosome=parts[0],
                                start=int(parts[1]),
                                end=int(parts[2]),
                                gene_type=parts[3] if len(parts) > 3 else "unknown",
                                score=float(parts[4]) if len(parts) > 4 else 0.0,
                                strand=parts[5] if len(parts) > 5 else "+"
                            )
                            candidates.append(candidate)
                    except (ValueError, IndexError) as e:
                        self.logger.warning(f"解析BED文件第{line_num}行失败: {e}")
                        continue
                        
        except Exception as e:
            self.logger.error(f"读取BED文件失败: {e}")
            
        return candidates
    
    def _get_output_files(self, output_dir: str) -> Dict[str, str]:
        """获取输出文件路径"""
        nlg_output_dir = os.path.join(output_dir, "NLGenomeSweeper")
        
        output_files = {}
        expected_files = [
            "All_candidates.bed",
            "Final_candidates.bed", 
            "Filtered_candidates.bed",
            "All_candidates.gff3",
            "Species_specific_consensus.fa"
        ]
        
        for filename in expected_files:
            file_path = os.path.join(nlg_output_dir, filename)
            if os.path.exists(file_path):
                output_files[filename] = file_path
        
        return output_files
    
    def run_parallel_search(self, 
                          genome_files: List[str],
                          output_base_dir: str,
                          **kwargs) -> Dict[str, NLGSearchResult]:
        """
        并行处理多个基因组文件
        
        Args:
            genome_files: 基因组文件列表
            output_base_dir: 输出基础目录
            **kwargs: 其他搜索参数
            
        Returns:
            Dict[str, NLGSearchResult]: 文件名到结果的映射
        """
        results = {}
        max_workers = kwargs.get('max_workers', 2)
        
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            # 提交任务
            future_to_file = {}
            for genome_file in genome_files:
                genome_name = Path(genome_file).stem
                output_dir = os.path.join(output_base_dir, genome_name)
                
                params = NLGSearchParams(
                    genome_file=genome_file,
                    output_dir=output_dir,
                    **{k: v for k, v in kwargs.items() if k != 'max_workers'}
                )
                
                future = executor.submit(self.run_search, params)
                future_to_file[future] = genome_name
            
            # 收集结果
            for future in as_completed(future_to_file):
                genome_name = future_to_file[future]
                try:
                    result = future.result()
                    results[genome_name] = result
                    self.logger.info(f"完成基因组 {genome_name} 的处理")
                except Exception as e:
                    self.logger.error(f"处理基因组 {genome_name} 失败: {e}")
        
        return results
    
    def validate_genome_file(self, genome_file: str) -> bool:
        """
        验证基因组文件格式
        
        Args:
            genome_file: 基因组文件路径
            
        Returns:
            bool: 是否为有效的FASTA文件
        """
        try:
            # 使用Biopython验证FASTA格式
            with open(genome_file, 'r') as f:
                sequences = SeqIO.parse(f, 'fasta')
                
                # 检查是否有序列
                sequence_count = 0
                valid_nucleotides = set('ATCGN')
                
                for seq in sequences:
                    sequence_count += 1
                    
                    # 检查序列长度
                    if len(seq.seq) == 0:
                        self.logger.warning(f"发现空序列: {seq.id}")
                        continue
                    
                    # 检查序列是否包含有效的核苷酸字符
                    seq_upper = str(seq.seq).upper()
                    invalid_chars = set(seq_upper) - valid_nucleotides
                    if invalid_chars:
                        self.logger.warning(f"序列 {seq.id} 包含无效字符: {invalid_chars}")
                        return False
                
                # 必须至少有一个序列
                return sequence_count > 0
                
        except Exception as e:
            self.logger.error(f"基因组文件验证失败: {e}")
            return False


def create_nlg_runner(config_path: Optional[str] = None) -> NLGenomeSweeperRunner:
    """
    创建NLGenomeSweeperRunner实例的便捷函数
    
    Args:
        config_path: 配置文件路径
        
    Returns:
        NLGenomeSweeperRunner: 配置好的运行器实例
    """
    return NLGenomeSweeperRunner(config_path)


# 示例使用代码
if __name__ == "__main__":
    # 示例：基本使用
    runner = create_nlg_runner()
    
    # 配置搜索参数
    params = NLGSearchParams(
        genome_file="/path/to/genome.fasta",
        output_dir="/path/to/output",
        threads=4,
        overwrite=True
    )
    
    # 执行搜索
    try:
        result = runner.run_search(params)
        print(f"找到 {result.total_candidates} 个候选基因")
        print(f"执行时间: {result.execution_time:.2f} 秒")
        
        # 导出结果
        result.to_bed_file("/path/to/output/candidates.bed")
        
    except Exception as e:
        print(f"搜索失败: {e}") 