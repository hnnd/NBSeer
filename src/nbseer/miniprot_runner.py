#!/usr/bin/env python3
"""
Miniprot集成模块

该模块提供了miniprot工具的Python封装，用于基于蛋白质序列进行基因结构预测。
miniprot是一个用于将蛋白质序列映射到基因组并预测基因结构的高效工具。

主要功能：
1. 蛋白质数据库索引构建
2. 基因组区域的基因结构预测
3. GFF3格式输出解析
4. 多线程并行处理支持
5. 预测结果的筛选和排序

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
from concurrent.futures import ThreadPoolExecutor, as_completed
import pandas as pd
from Bio import SeqIO
import re

from src.utils.config import get_config
from src.utils.logging_setup import get_logger


@dataclass
class MiniprotParams:
    """Miniprot运行参数"""
    reference_genome: str
    protein_database: str
    output_dir: str
    threads: int = 4
    max_intron_size: str = "200k"
    gap_open_penalty: int = 11
    gap_extension: int = 1
    intron_open_penalty: int = 29
    frameshift_penalty: int = 23
    splice_weight: float = 1.0
    bonus_score: int = 5
    splice_model: int = 1  # 1=general, 2=mammal, 0=none
    output_format: str = "gff"  # gff, gtf, aln, trans
    min_score_ratio: float = 0.97  # --outs参数
    min_coverage: float = 0.1  # --outc参数
    max_alignments: int = 1000  # --outn参数
    overwrite: bool = False
    
    def __post_init__(self):
        """验证参数"""
        if not os.path.exists(self.reference_genome):
            raise FileNotFoundError(f"参考基因组文件不存在: {self.reference_genome}")
        
        if not os.path.exists(self.protein_database):
            raise FileNotFoundError(f"蛋白质数据库文件不存在: {self.protein_database}")
        
        if self.threads < 1:
            raise ValueError("线程数必须大于0")
        
        if not 0 < self.min_score_ratio <= 1:
            raise ValueError("最小得分比例必须在0-1之间")
        
        if not 0 < self.min_coverage <= 1:
            raise ValueError("最小覆盖度必须在0-1之间")


@dataclass
class GeneStructure:
    """基因结构信息"""
    query_id: str
    chromosome: str
    start: int
    end: int
    strand: str
    score: float
    coverage: float
    identity: float
    exons: List[Tuple[int, int]]  # [(start, end), ...]
    cds_regions: List[Tuple[int, int]]  # [(start, end), ...]
    protein_match: str
    alignment_length: int
    
    def __post_init__(self):
        """计算衍生属性"""
        if self.exons:
            self.exon_count = len(self.exons)
            self.total_exon_length = sum(end - start + 1 for start, end in self.exons)
        else:
            self.exon_count = 0
            self.total_exon_length = 0
    
    @property
    def length(self) -> int:
        """基因总长度"""
        return self.end - self.start + 1
    
    @property
    def coordinates(self) -> str:
        """基因坐标字符串"""
        return f"{self.chromosome}:{self.start}-{self.end}({self.strand})"
    
    def to_gff_line(self) -> str:
        """转换为GFF格式行"""
        attributes = f"ID={self.query_id};Score={self.score:.3f};Coverage={self.coverage:.3f};Identity={self.identity:.3f}"
        return f"{self.chromosome}\tminiprot\tgene\t{self.start}\t{self.end}\t{self.score:.1f}\t{self.strand}\t.\t{attributes}"


@dataclass
class MiniprotResult:
    """Miniprot预测结果"""
    gene_structures: List[GeneStructure]
    output_files: Dict[str, str]
    params: MiniprotParams
    execution_time: float
    total_predictions: int
    
    def get_best_predictions(self, top_n: int = 10) -> List[GeneStructure]:
        """获取得分最高的预测结果"""
        sorted_genes = sorted(self.gene_structures, key=lambda x: x.score, reverse=True)
        return sorted_genes[:top_n]
    
    def filter_by_score(self, min_score: float) -> List[GeneStructure]:
        """按得分筛选预测结果"""
        return [gene for gene in self.gene_structures if gene.score >= min_score]
    
    def filter_by_coverage(self, min_coverage: float) -> List[GeneStructure]:
        """按覆盖度筛选预测结果"""
        return [gene for gene in self.gene_structures if gene.coverage >= min_coverage]
    
    def get_summary(self) -> Dict:
        """获取结果摘要"""
        if not self.gene_structures:
            return {
                "total_predictions": 0,
                "execution_time": self.execution_time,
                "reference_genome": self.params.reference_genome,
                "protein_database": self.params.protein_database
            }
        
        scores = [gene.score for gene in self.gene_structures]
        coverages = [gene.coverage for gene in self.gene_structures]
        
        return {
            "total_predictions": self.total_predictions,
            "execution_time": self.execution_time,
            "reference_genome": self.params.reference_genome,
            "protein_database": self.params.protein_database,
            "score_stats": {
                "min": min(scores),
                "max": max(scores),
                "mean": sum(scores) / len(scores)
            },
            "coverage_stats": {
                "min": min(coverages),
                "max": max(coverages),
                "mean": sum(coverages) / len(coverages)
            },
            "chromosomes": list(set(gene.chromosome for gene in self.gene_structures))
        }
    
    def to_gff_file(self, output_path: str) -> None:
        """导出为GFF3格式文件"""
        with open(output_path, 'w') as f:
            f.write("##gff-version 3\n")
            f.write("# Miniprot基因结构预测结果\n")
            f.write(f"# 总预测数: {self.total_predictions}\n")
            f.write(f"# 参考基因组: {self.params.reference_genome}\n")
            f.write(f"# 蛋白质数据库: {self.params.protein_database}\n")
            
            for gene in self.gene_structures:
                f.write(gene.to_gff_line() + "\n")
                
                # 写入外显子信息
                for i, (start, end) in enumerate(gene.exons, 1):
                    exon_attrs = f"ID={gene.query_id}.exon{i};Parent={gene.query_id}"
                    f.write(f"{gene.chromosome}\tminiprot\texon\t{start}\t{end}\t.\t{gene.strand}\t.\t{exon_attrs}\n")


class MiniprotRunner:
    """
    Miniprot工具的Python封装类
    
    该类封装了miniprot命令行工具的调用，提供了Python接口
    用于基于蛋白质序列进行基因结构预测。
    
    主要功能：
    - 蛋白质数据库索引构建
    - 基因结构预测
    - 结果解析和筛选
    - 并行处理支持
    """
    
    def __init__(self, config_path: Optional[str] = None):
        """
        初始化MiniprotRunner
        
        Args:
            config_path: 配置文件路径，如果为None则使用默认配置
        """
        self.config = get_config(config_path)
        self.logger = get_logger(__name__)
        
        # 获取工具路径
        self.miniprot_executable = self._get_miniprot_executable()
        self.tool_config = self.config.get_tool_config('miniprot')
        
        # 验证工具可用性
        self._validate_dependencies()
        
        self.logger.info("MiniprotRunner初始化完成")
    
    def _get_miniprot_executable(self) -> str:
        """获取miniprot可执行文件路径"""
        # 从配置中获取路径
        miniprot_path = self.config.get("tools.miniprot.executable_path")
        
        if miniprot_path and os.path.exists(miniprot_path):
            return miniprot_path
        
        # 尝试在PATH中查找
        miniprot_in_path = shutil.which("miniprot")
        if miniprot_in_path:
            return miniprot_in_path
        
        raise FileNotFoundError(
            "无法找到miniprot可执行文件。请检查配置或安装。"
        )
    
    def _validate_dependencies(self) -> None:
        """验证所需依赖工具是否可用"""
        try:
            # 运行miniprot无参数来测试
            result = subprocess.run(
                [self.miniprot_executable],
                capture_output=True,
                text=True,
                timeout=10
            )
            # miniprot无参数时返回1，但这是正常的
            if result.returncode in [0, 1]:
                self.logger.info("Miniprot工具验证成功")
            else:
                self.logger.warning("Miniprot工具可能无法正常运行")
        except Exception as e:
            self.logger.warning(f"检查miniprot工具失败: {e}")
    
    def build_genome_index(self, genome_fasta: str, index_file: str) -> bool:
        """
        构建基因组索引
        
        Args:
            genome_fasta: 基因组FASTA文件路径
            index_file: 输出索引文件路径
            
        Returns:
            bool: 索引构建是否成功
        """
        self.logger.info(f"开始构建基因组索引: {genome_fasta}")
        
        try:
            # 验证基因组文件
            if not self._validate_genome_file(genome_fasta):
                self.logger.error(f"基因组文件验证失败: {genome_fasta}")
                return False
            
            cmd = [
                self.miniprot_executable,
                "-d", index_file,
                genome_fasta
            ]
            
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )
            
            if os.path.exists(index_file):
                index_size = os.path.getsize(index_file)
                self.logger.info(f"基因组索引构建成功: {index_file} (大小: {index_size / 1024 / 1024:.1f} MB)")
                return True
            else:
                self.logger.error("索引文件未生成")
                return False
                
        except subprocess.CalledProcessError as e:
            self.logger.error(f"基因组索引构建失败: {e}")
            self.logger.error(f"标准错误输出: {e.stderr}")
            return False
    
    def _validate_genome_file(self, genome_file: str) -> bool:
        """
        验证基因组文件格式
        
        Args:
            genome_file: 基因组文件路径
            
        Returns:
            bool: 是否为有效的基因组FASTA文件
        """
        try:
            # 使用Biopython验证FASTA格式
            with open(genome_file, 'r') as f:
                sequences = SeqIO.parse(f, 'fasta')
                
                sequence_count = 0
                total_length = 0
                valid_nucleotides = set('ATCGN')
                
                for seq in sequences:
                    sequence_count += 1
                    seq_len = len(seq.seq)
                    total_length += seq_len
                    
                    # 检查序列长度
                    if seq_len == 0:
                        self.logger.warning(f"发现空序列: {seq.id}")
                        continue
                    
                    # 检查序列是否包含有效的核苷酸字符
                    seq_upper = str(seq.seq).upper()
                    invalid_chars = set(seq_upper) - valid_nucleotides
                    if invalid_chars:
                        # 允许一些IUPAC核苷酸代码
                        extended_nucleotides = set('ATCGNRYSWKMBDHV')
                        extended_invalid = set(seq_upper) - extended_nucleotides
                        if extended_invalid:
                            self.logger.warning(f"序列 {seq.id} 包含非核苷酸字符: {extended_invalid}")
                            if len(extended_invalid) > 10:  # 如果无效字符太多才返回False
                                return False
                
                # 必须至少有一个序列
                self.logger.info(f"基因组文件验证通过，共 {sequence_count} 个序列，总长度 {total_length:,} bp")
                return sequence_count > 0
                
        except Exception as e:
            self.logger.error(f"基因组文件验证失败: {e}")
            return False
    
    def validate_protein_database(self, protein_file: str) -> bool:
        """
        验证蛋白质数据库文件格式
        
        Args:
            protein_file: 蛋白质文件路径
            
        Returns:
            bool: 是否为有效的蛋白质FASTA文件
        """
        try:
            # 使用Biopython验证FASTA格式
            with open(protein_file, 'r') as f:
                sequences = SeqIO.parse(f, 'fasta')
                
                sequence_count = 0
                valid_amino_acids = set('ACDEFGHIKLMNPQRSTVWYX*')
                
                for seq in sequences:
                    sequence_count += 1
                    
                    # 检查序列长度
                    if len(seq.seq) == 0:
                        self.logger.warning(f"发现空序列: {seq.id}")
                        continue
                    
                    # 检查序列是否包含有效的氨基酸字符
                    seq_upper = str(seq.seq).upper()
                    invalid_chars = set(seq_upper) - valid_amino_acids
                    if invalid_chars:
                        self.logger.warning(f"序列 {seq.id} 包含无效字符: {invalid_chars}")
                        # 对于蛋白质序列，允许一些非标准字符
                        if len(invalid_chars) > 5:  # 如果无效字符太多才返回False
                            return False
                
                # 必须至少有一个序列
                self.logger.info(f"蛋白质数据库验证通过，共 {sequence_count} 个序列")
                return sequence_count > 0
                
        except Exception as e:
            self.logger.error(f"蛋白质数据库验证失败: {e}")
            return False
    
    def direct_alignment(self, 
                        genome_fasta: str,
                        protein_fasta: str,
                        output_gff: str,
                        threads: int = 8,
                        min_score_ratio: float = 0.8,
                        min_coverage: float = 0.1,
                        max_intron_size: str = "200k") -> MiniprotResult:
        """
        直接比对模式：将蛋白质序列直接比对到基因组
        
        Args:
            genome_fasta: 基因组FASTA文件路径
            protein_fasta: 蛋白质FASTA文件路径  
            output_gff: 输出GFF文件路径
            threads: 线程数
            min_score_ratio: 最小得分比例(--outs)
            min_coverage: 最小覆盖度(--outc)
            max_intron_size: 最大内含子大小(-G)
            
        Returns:
            MiniprotResult: 预测结果对象
        """
        import time
        
        self.logger.info(f"开始直接比对: {protein_fasta} -> {genome_fasta}")
        self.logger.info(f"输出文件: {output_gff}")
        
        # 验证输入文件
        if not os.path.exists(genome_fasta):
            raise FileNotFoundError(f"基因组文件不存在: {genome_fasta}")
        if not os.path.exists(protein_fasta):
            raise FileNotFoundError(f"蛋白质文件不存在: {protein_fasta}")
        
        # 验证基因组和蛋白质文件
        self.logger.info("验证输入文件...")
        if not self._validate_genome_file(genome_fasta):
            raise ValueError(f"基因组文件验证失败: {genome_fasta}")
        if not self.validate_protein_database(protein_fasta):
            raise ValueError(f"蛋白质文件验证失败: {protein_fasta}")
        
        # 创建输出目录
        output_dir = os.path.dirname(output_gff)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)
        
        # 构建命令
        cmd = [
            self.miniprot_executable,
            "--gff",                    # 输出GFF格式
            "-t", str(threads),         # 线程数
            "-G", max_intron_size,      # 最大内含子大小
            "--outs", str(min_score_ratio),  # 最小得分比例
            "--outc", str(min_coverage),     # 最小覆盖度
            genome_fasta,               # 基因组文件
            protein_fasta               # 蛋白质文件
        ]
        
        self.logger.info(f"执行命令: {' '.join(cmd)}")
        
        # 记录开始时间
        start_time = time.time()
        
        try:
            # 执行命令并重定向输出到文件
            with open(output_gff, 'w') as f:
                result = subprocess.run(
                    cmd,
                    stdout=f,
                    stderr=subprocess.PIPE,
                    text=True,
                    check=True
                )
            
            execution_time = time.time() - start_time
            self.logger.info(f"比对完成，耗时: {execution_time:.2f}秒")
            
            # 检查输出文件
            if not os.path.exists(output_gff):
                raise RuntimeError(f"输出文件未生成: {output_gff}")
            
            gff_size = os.path.getsize(output_gff)
            self.logger.info(f"输出GFF文件: {output_gff} (大小: {gff_size / 1024:.1f} KB)")
            
            # 解析结果
            gene_structures = self._parse_gff_output(output_gff)
            
            # 创建MiniprotParams对象（用于结果记录）
            params = MiniprotParams(
                reference_genome=genome_fasta,
                protein_database=protein_fasta,
                output_dir=output_dir,
                threads=threads,
                min_score_ratio=min_score_ratio,
                min_coverage=min_coverage,
                max_intron_size=max_intron_size
            )
            
            # 创建结果对象
            result_obj = MiniprotResult(
                gene_structures=gene_structures,
                output_files={"gff": output_gff},
                params=params,
                execution_time=execution_time,
                total_predictions=len(gene_structures)
            )
            
            self.logger.info(f"解析完成，共预测 {len(gene_structures)} 个基因结构")
            
            return result_obj
            
        except subprocess.CalledProcessError as e:
            error_msg = f"miniprot执行失败: {e.stderr}"
            self.logger.error(error_msg)
            raise RuntimeError(error_msg)
        except Exception as e:
            self.logger.error(f"直接比对过程中发生错误: {str(e)}")
            raise

    def run_prediction(self, params: MiniprotParams) -> Dict:
        """执行基因结构预测"""
        import time
        start_time = time.time()
        
        self.logger.info(f"开始miniprot基因结构预测")
        self.logger.info(f"参考基因组: {params.reference_genome}")
        self.logger.info(f"蛋白质数据库: {params.protein_database}")
        
        # 准备输出目录
        self._prepare_output_directory(params)
        
        # 构建命令
        cmd = self._build_command(params)
        
        # 执行预测
        try:
            output_file = self._execute_command(cmd, params)
            execution_time = time.time() - start_time
            
            result = {
                "output_file": output_file,
                "params": params,
                "execution_time": execution_time,
                "success": True
            }
            
            self.logger.info(f"预测完成，耗时 {execution_time:.2f} 秒")
            return result
            
        except Exception as e:
            self.logger.error(f"Miniprot预测失败: {e}")
            raise
    
    def _prepare_output_directory(self, params: MiniprotParams) -> None:
        """准备输出目录"""
        output_path = Path(params.output_dir)
        
        if output_path.exists() and params.overwrite:
            self.logger.info(f"清理现有输出目录: {output_path}")
            shutil.rmtree(output_path)
        
        output_path.mkdir(parents=True, exist_ok=True)
        self.logger.info(f"创建输出目录: {output_path}")
    
    def _build_command(self, params: MiniprotParams) -> List[str]:
        """构建miniprot命令"""
        cmd = [
            self.miniprot_executable,
            "-t", str(params.threads),
            "-G", params.max_intron_size,
            "-O", str(params.gap_open_penalty),
            "-E", str(params.gap_extension),
            "-J", str(params.intron_open_penalty),
            "-F", str(params.frameshift_penalty),
            "-C", str(params.splice_weight),
            "-B", str(params.bonus_score),
            "-j", str(params.splice_model),
            f"--outs={params.min_score_ratio}",
            f"--outc={params.min_coverage}",
            f"--outn={params.max_alignments}",
            params.reference_genome,
            params.protein_database
        ]
        
        # 添加输出格式
        if params.output_format == "gff":
            cmd.append("--gff")
        elif params.output_format == "gtf":
            cmd.append("--gtf")
        elif params.output_format == "aln":
            cmd.append("--aln")
        elif params.output_format == "trans":
            cmd.append("--trans")
        
        self.logger.debug(f"构建命令: {' '.join(cmd)}")
        return cmd
    
    def _execute_command(self, cmd: List[str], params: MiniprotParams) -> str:
        """执行miniprot命令"""
        output_file = os.path.join(params.output_dir, f"miniprot_output.{params.output_format}")
        
        try:
            with open(output_file, 'w') as f:
                result = subprocess.run(
                    cmd,
                    stdout=f,
                    stderr=subprocess.PIPE,
                    text=True,
                    check=True,
                    cwd=params.output_dir
                )
            
            # 记录错误输出（如果有）
            if result.stderr:
                self.logger.debug(f"Miniprot stderr: {result.stderr}")
            
            self.logger.info(f"Miniprot执行成功，输出文件: {output_file}")
            return output_file
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Miniprot执行失败: {e}")
            self.logger.error(f"标准错误输出: {e.stderr}")
            raise RuntimeError(f"Miniprot执行失败: {e.stderr}")
    
    def _parse_results(self, output_file: str, params: MiniprotParams) -> List[GeneStructure]:
        """解析miniprot结果"""
        gene_structures = []
        
        if params.output_format == "gff":
            gene_structures = self._parse_gff_output(output_file)
        else:
            self.logger.warning(f"暂不支持解析 {params.output_format} 格式的输出")
        
        return gene_structures
    
    def _parse_gff_output(self, gff_file: str) -> List[GeneStructure]:
        """解析GFF格式输出"""
        gene_structures = []
        current_gene = None
        
        try:
            with open(gff_file, 'r') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    
                    # 跳过注释和空行
                    if line.startswith('#') or not line:
                        continue
                    
                    try:
                        parts = line.split('\t')
                        if len(parts) < 9:
                            continue
                        
                        chromosome, source, feature_type, start, end, score, strand, phase, attributes = parts
                        
                        # 解析属性
                        attr_dict = self._parse_gff_attributes(attributes)
                        
                        if feature_type == "mRNA" or feature_type == "transcript":
                            # 创建新的基因结构
                            current_gene = GeneStructure(
                                query_id=attr_dict.get("ID", f"gene_{line_num}"),
                                chromosome=chromosome,
                                start=int(start),
                                end=int(end),
                                strand=strand,
                                score=float(score) if score != '.' else 0.0,
                                coverage=float(attr_dict.get("Coverage", "0")),
                                identity=float(attr_dict.get("Identity", "0")),
                                exons=[],
                                cds_regions=[],
                                protein_match=attr_dict.get("Target", ""),
                                alignment_length=int(end) - int(start) + 1
                            )
                            gene_structures.append(current_gene)
                        
                        elif feature_type == "exon" and current_gene:
                            # 添加外显子
                            current_gene.exons.append((int(start), int(end)))
                        
                        elif feature_type == "CDS" and current_gene:
                            # 添加CDS区域
                            current_gene.cds_regions.append((int(start), int(end)))
                    
                    except (ValueError, IndexError) as e:
                        self.logger.warning(f"解析GFF文件第{line_num}行失败: {e}")
                        continue
        
        except Exception as e:
            self.logger.error(f"读取GFF文件失败: {e}")
        
        return gene_structures
    
    def _parse_gff_attributes(self, attributes: str) -> Dict[str, str]:
        """解析GFF属性字符串"""
        attr_dict = {}
        
        for attr in attributes.split(';'):
            if '=' in attr:
                key, value = attr.split('=', 1)
                attr_dict[key.strip()] = value.strip()
        
        return attr_dict
    
    def _get_output_files(self, output_dir: str) -> Dict[str, str]:
        """获取输出文件路径"""
        output_files = {}
        
        # 查找输出文件
        for filename in os.listdir(output_dir):
            file_path = os.path.join(output_dir, filename)
            if os.path.isfile(file_path):
                output_files[filename] = file_path
        
        return output_files
    
    def predict_candidate_regions(self, 
                                candidate_regions: List[Dict],  # [{"chr": str, "start": int, "end": int}, ...]
                                genome_index: str,
                                protein_database: str,
                                output_dir: str,
                                **kwargs) -> MiniprotResult:
        """
        在候选区域运行miniprot预测
        
        Args:
            candidate_regions: 候选区域列表，每个区域包含chr, start, end
            genome_index: 基因组索引文件路径
            protein_database: 蛋白质数据库文件路径
            output_dir: 输出目录
            **kwargs: 其他miniprot参数
            
        Returns:
            MiniprotResult: 预测结果
        """
        import time
        start_time = time.time()
        
        self.logger.info(f"开始在 {len(candidate_regions)} 个候选区域运行miniprot预测")
        self.logger.info(f"基因组索引: {genome_index}")
        self.logger.info(f"蛋白质数据库: {protein_database}")
        
        # 准备输出目录
        output_path = Path(output_dir)
        if output_path.exists() and kwargs.get('overwrite', False):
            shutil.rmtree(output_path)
        output_path.mkdir(parents=True, exist_ok=True)
        
        # 构建命令
        cmd = self._build_candidate_prediction_command(
            genome_index, protein_database, candidate_regions, output_dir, **kwargs
        )
        
        # 执行预测
        try:
            output_file = self._execute_candidate_prediction(cmd, output_dir, **kwargs)
            execution_time = time.time() - start_time
            
            # 解析结果
            gene_structures = self._parse_gff_output(output_file)
            
            result = MiniprotResult(
                gene_structures=gene_structures,
                output_files=self._get_output_files(output_dir),
                params=None,  # 将在后续版本中改进
                execution_time=execution_time,
                total_predictions=len(gene_structures)
            )
            
            self.logger.info(f"候选区域预测完成，找到 {len(gene_structures)} 个基因结构")
            self.logger.info(f"执行时间: {execution_time:.2f} 秒")
            
            return result
            
        except Exception as e:
            self.logger.error(f"候选区域预测失败: {e}")
            raise
    
    def _build_candidate_prediction_command(self, 
                                          genome_index: str,
                                          protein_database: str,
                                          candidate_regions: List[Dict],
                                          output_dir: str,
                                          **kwargs) -> List[str]:
        """构建候选区域预测命令"""
        
        # 转换为绝对路径，避免工作目录问题
        genome_index_abs = os.path.abspath(genome_index)
        protein_database_abs = os.path.abspath(protein_database)
        
        # 基础命令
        cmd = [
            self.miniprot_executable,
            "-t", str(kwargs.get('threads', 4)),
            "-G", kwargs.get('max_intron_size', '200k'),
            "-O", str(kwargs.get('gap_open_penalty', 11)),
            "-E", str(kwargs.get('gap_extension', 1)),
            "-J", str(kwargs.get('intron_open_penalty', 29)),
            "-F", str(kwargs.get('frameshift_penalty', 23)),
            "-C", str(kwargs.get('splice_weight', 1.0)),
            "-B", str(kwargs.get('bonus_score', 5)),
            "-j", str(kwargs.get('splice_model', 1)),
            f"--outs={kwargs.get('min_score_ratio', 0.97)}",
            f"--outc={kwargs.get('min_coverage', 0.1)}",
            f"--outn={kwargs.get('max_alignments', 1000)}",
            "--gff",  # 输出GFF格式
            genome_index_abs,  # 使用基因组索引的绝对路径
            protein_database_abs  # 蛋白质查询文件的绝对路径
        ]
        
        # 如果有候选区域，可以添加区域限制（如果miniprot支持）
        # 这里暂时不添加区域限制，因为miniprot会在整个基因组上预测
        # 后续可以通过后处理来筛选候选区域内的预测结果
        
        self.logger.debug(f"构建候选区域预测命令: {' '.join(cmd)}")
        return cmd
    
    def _execute_candidate_prediction(self, cmd: List[str], output_dir: str, **kwargs) -> str:
        """执行候选区域预测命令"""
        output_file = os.path.join(output_dir, "candidate_predictions.gff")
        
        try:
            with open(output_file, 'w') as f:
                result = subprocess.run(
                    cmd,
                    stdout=f,
                    stderr=subprocess.PIPE,
                    text=True,
                    check=True,
                    cwd=output_dir
                )
            
            # 记录错误输出（如果有）
            if result.stderr:
                self.logger.debug(f"Miniprot stderr: {result.stderr}")
            
            self.logger.info(f"候选区域预测执行成功，输出文件: {output_file}")
            return output_file
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"候选区域预测执行失败: {e}")
            self.logger.error(f"标准错误输出: {e.stderr}")
            raise RuntimeError(f"候选区域预测执行失败: {e.stderr}")
    
    def filter_predictions_by_regions(self, 
                                    gene_structures: List[GeneStructure],
                                    candidate_regions: List[Dict]) -> List[GeneStructure]:
        """
        根据候选区域筛选预测结果
        
        Args:
            gene_structures: 基因结构预测列表
            candidate_regions: 候选区域列表
            
        Returns:
            List[GeneStructure]: 筛选后的基因结构列表
        """
        filtered_structures = []
        
        for gene in gene_structures:
            for region in candidate_regions:
                # 检查基因是否与候选区域重叠
                if (gene.chromosome == region['chr'] and
                    gene.start <= region['end'] and
                    gene.end >= region['start']):
                    
                    # 计算重叠长度
                    overlap_start = max(gene.start, region['start'])
                    overlap_end = min(gene.end, region['end'])
                    overlap_length = overlap_end - overlap_start + 1
                    
                    # 如果重叠长度大于基因长度的50%，则保留
                    if overlap_length >= gene.length * 0.5:
                        filtered_structures.append(gene)
                        break  # 避免重复添加
        
        self.logger.info(f"根据候选区域筛选：{len(gene_structures)} -> {len(filtered_structures)} 个预测")
        return filtered_structures

    def run_parallel_prediction(self, 
                              genome_regions: List[Tuple[str, str, str]],  # [(chr, start, end), ...]
                              params: MiniprotParams,
                              max_workers: int = 4) -> Dict[str, MiniprotResult]:
        """
        并行处理多个基因组区域
        
        Args:
            genome_regions: 基因组区域列表 [(chromosome, start, end), ...]
            params: 基础预测参数
            max_workers: 最大并行工作线程数
            
        Returns:
            Dict[str, MiniprotResult]: 区域标识到结果的映射
        """
        results = {}
        
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            # 提交任务
            future_to_region = {}
            for i, (chromosome, start, end) in enumerate(genome_regions):
                region_id = f"{chromosome}_{start}_{end}"
                
                # 为每个区域创建独立的参数
                region_params = MiniprotParams(
                    reference_genome=params.reference_genome,
                    protein_database=params.protein_database,
                    output_dir=os.path.join(params.output_dir, f"region_{i}"),
                    threads=1,  # 每个线程使用单线程
                    **{k: v for k, v in params.__dict__.items() 
                       if k not in ['reference_genome', 'protein_database', 'output_dir', 'threads']}
                )
                
                future = executor.submit(self.run_prediction, region_params)
                future_to_region[future] = region_id
            
            # 收集结果
            for future in as_completed(future_to_region):
                region_id = future_to_region[future]
                try:
                    result = future.result()
                    results[region_id] = result
                    self.logger.info(f"完成区域 {region_id} 的预测")
                except Exception as e:
                    self.logger.error(f"处理区域 {region_id} 失败: {e}")
        
        return results


def create_miniprot_runner(config_path: Optional[str] = None) -> MiniprotRunner:
    """
    创建MiniprotRunner实例的便捷函数
    
    Args:
        config_path: 配置文件路径
        
    Returns:
        MiniprotRunner: 配置好的运行器实例
    """
    return MiniprotRunner(config_path)


# 示例使用代码
if __name__ == "__main__":
    # 示例：基本使用
    runner = create_miniprot_runner()
    
    # 配置预测参数
    params = MiniprotParams(
        reference_genome="/path/to/genome.fasta",
        protein_database="/path/to/proteins.fasta",
        output_dir="/path/to/output",
        threads=4,
        min_score_ratio=0.97,
        overwrite=True
    )
    
    # 执行预测
    try:
        result = runner.run_prediction(params)
        print(f"找到 {result['total_predictions']} 个基因结构预测")
        print(f"执行时间: {result['execution_time']:.2f} 秒")
        
        # 获取最佳预测
        best_predictions = result['gene_structures'][:5]
        for i, gene in enumerate(best_predictions, 1):
            print(f"预测 {i}: {gene.coordinates}, 得分: {gene.score:.3f}")
        
        # 导出结果
        result['gene_structures'][-1].to_gff_file("/path/to/output/predictions.gff3")
        
    except Exception as e:
        print(f"预测失败: {e}") 