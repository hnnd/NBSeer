#!/usr/bin/env python3
"""
EVM Runner for NBS Genome Annotation Pipeline
运行EVidenceModeler整合多个预测结果，生成最终基因注释
"""

import os
import sys
import json
import shutil
import logging
import subprocess
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple
from dataclasses import dataclass, asdict
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
import multiprocessing as mp
from datetime import datetime
import tempfile
import time

# 添加项目根目录到Python路径
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

from src.utils.config import ConfigManager
from src.nbseer.memory_manager import GenomeMemoryManager, PartitioningConfig


@dataclass
class EVMConfig:
    """EVM配置类"""
    genome_fasta: str
    augustus_predictions: str
    miniprot_predictions: str
    weights_file: str
    output_dir: str
    evm_binary: str
    partition_size: int = 1000000  # 1MB分区大小
    threads: int = 4
    memory_limit: str = "8G"
    debug: bool = False
    # 新增：miniprot处理结果配置
    use_processed_miniprot: bool = True
    miniprot_processed_dir: str = "results/miniprot_processed"
    include_low_quality: bool = True  # 是否包含低质量预测


@dataclass 
class EVMStats:
    """EVM统计结果"""
    total_genes_predicted: int = 0
    total_exons: int = 0
    total_introns: int = 0
    augustus_genes_used: int = 0
    miniprot_genes_used: int = 0
    miniprot_high_quality: int = 0
    miniprot_medium_quality: int = 0
    miniprot_low_quality: int = 0
    processing_time: float = 0.0
    partition_count: int = 0
    success_partitions: int = 0
    failed_partitions: int = 0


class EVMRunner:
    """EVM执行器类"""
    
    def __init__(self, config: EVMConfig, logger: Optional[logging.Logger] = None):
        """
        初始化EVM运行器
        
        Args:
            config: EVM配置
            logger: 日志记录器
        """
        self.config = config
        self.logger = logger or self._setup_logger()
        self.stats = EVMStats()
        self.checkpoints = {}
        self.temp_dir = None
        
        # 验证输入文件和工具
        self._validate_inputs()
        
    def _setup_logger(self) -> logging.Logger:
        """设置日志记录器"""
        logger = logging.getLogger("EVMRunner")
        logger.setLevel(logging.DEBUG if self.config.debug else logging.INFO)
        
        if not logger.handlers:
            handler = logging.StreamHandler()
            formatter = logging.Formatter(
                '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
            )
            handler.setFormatter(formatter)
            logger.addHandler(handler)
            
        return logger
        
    def _validate_inputs(self):
        """验证输入文件和工具"""
        self.logger.info("验证输入文件和工具...")
        
        # 检查基因组文件
        if not os.path.exists(self.config.genome_fasta):
            raise FileNotFoundError(f"基因组文件不存在: {self.config.genome_fasta}")
            
        # 检查Augustus预测结果目录
        if not os.path.exists(self.config.augustus_predictions):
            raise FileNotFoundError(f"Augustus预测目录不存在: {self.config.augustus_predictions}")
            
        # 检查miniprot预测结果
        if self.config.use_processed_miniprot:
            # 检查处理后的miniprot文件
            processed_dir = Path(self.config.miniprot_processed_dir)
            if not processed_dir.exists():
                raise FileNotFoundError(f"Miniprot处理目录不存在: {processed_dir}")
                
            required_files = [
                "miniprot_high_quality.gff3",
                "miniprot_medium_quality.gff3", 
                "miniprot_quality_report.json"
            ]
            if self.config.include_low_quality:
                required_files.append("miniprot_low_quality.gff3")
                
            for file in required_files:
                file_path = processed_dir / file
                if not file_path.exists():
                    raise FileNotFoundError(f"Miniprot处理文件不存在: {file_path}")
        else:
            # 检查原始miniprot文件
            if not os.path.exists(self.config.miniprot_predictions):
                raise FileNotFoundError(f"Miniprot预测文件不存在: {self.config.miniprot_predictions}")
            
        # 检查EVM二进制文件
        if not os.path.exists(self.config.evm_binary):
            raise FileNotFoundError(f"EVM二进制文件不存在: {self.config.evm_binary}")
            
        # 创建输出目录
        os.makedirs(self.config.output_dir, exist_ok=True)
        
        self.logger.info("输入验证完成")

    def _check_miniprot_processing_status(self) -> Dict[str, Any]:
        """检查miniprot处理状态"""
        self.logger.info("检查miniprot处理状态...")
        
        processed_dir = Path(self.config.miniprot_processed_dir)
        status = {
            "processed": False,
            "quality_report": {},
            "files": {},
            "total_predictions": 0
        }
        
        if not processed_dir.exists():
            self.logger.warning(f"Miniprot处理目录不存在: {processed_dir}")
            return status
            
        # 检查质量报告
        report_file = processed_dir / "miniprot_quality_report.json"
        if report_file.exists():
            try:
                with open(report_file, 'r') as f:
                    quality_report = json.load(f)
                    status["quality_report"] = quality_report
                    status["total_predictions"] = quality_report.get("processing_info", {}).get("total_predictions", 0)
                    
                    # 提取质量统计
                    quality_stats = quality_report.get("quality_statistics", {})
                    self.stats.miniprot_high_quality = quality_stats.get("high_quality", {}).get("count", 0)
                    self.stats.miniprot_medium_quality = quality_stats.get("medium_quality", {}).get("count", 0)
                    self.stats.miniprot_low_quality = quality_stats.get("low_quality", {}).get("count", 0)
                    
                    self.logger.info(f"Miniprot质量统计 - 高质量: {self.stats.miniprot_high_quality}, "
                                   f"中质量: {self.stats.miniprot_medium_quality}, "
                                   f"低质量: {self.stats.miniprot_low_quality}")
                    
            except Exception as e:
                self.logger.error(f"读取质量报告失败: {e}")
                return status
        
        # 检查处理后的文件
        quality_files = ["high_quality", "medium_quality", "low_quality"]
        for quality in quality_files:
            file_path = processed_dir / f"miniprot_{quality}.gff3"
            if file_path.exists():
                file_size = file_path.stat().st_size
                status["files"][quality] = {
                    "path": str(file_path),
                    "size_mb": round(file_size / (1024 * 1024), 2),
                    "exists": True
                }
            else:
                status["files"][quality] = {"exists": False}
        
        # 检查权重配置文件
        weights_file = Path("evm_config/miniprot_weights.txt")
        if weights_file.exists():
            status["weights_config"] = str(weights_file)
            
        status["processed"] = all(status["files"][q]["exists"] for q in ["high_quality", "medium_quality"])
        
        if status["processed"]:
            self.logger.info("Miniprot处理结果验证通过")
        else:
            self.logger.warning("Miniprot处理结果不完整")
            
        return status
        
    def _setup_temp_directory(self) -> str:
        """设置临时工作目录"""
        if self.temp_dir is None:
            self.temp_dir = tempfile.mkdtemp(prefix="evm_", suffix="_temp")
            self.logger.info(f"创建临时目录: {self.temp_dir}")
        return self.temp_dir
        
    def _cleanup_temp_directory(self):
        """清理临时目录"""
        if self.temp_dir and os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)
            self.logger.info(f"清理临时目录: {self.temp_dir}")
            self.temp_dir = None
            
    def _create_checkpoint(self, step: str, data: Dict[str, Any]):
        """创建检查点"""
        checkpoint_file = os.path.join(self.config.output_dir, f"checkpoint_{step}.json")
        self.checkpoints[step] = {
            'timestamp': datetime.now().isoformat(),
            'data': data
        }
        
        with open(checkpoint_file, 'w') as f:
            json.dump(self.checkpoints[step], f, indent=2)
            
        self.logger.info(f"创建检查点: {step}")
        
    def _load_checkpoint(self, step: str) -> Optional[Dict[str, Any]]:
        """加载检查点"""
        checkpoint_file = os.path.join(self.config.output_dir, f"checkpoint_{step}.json")
        
        if os.path.exists(checkpoint_file):
            try:
                with open(checkpoint_file, 'r') as f:
                    checkpoint = json.load(f)
                    self.logger.info(f"加载检查点: {step}")
                    return checkpoint.get('data')
            except Exception as e:
                self.logger.warning(f"加载检查点失败 {step}: {e}")
                
        return None
        
    def _combine_augustus_predictions(self) -> str:
        """合并所有Augustus预测结果（使用已转换的绝对坐标）"""
        self.logger.info("合并Augustus预测结果（已转换坐标）...")
        
        # 检查是否有检查点
        checkpoint_data = self._load_checkpoint("augustus_combine")
        if checkpoint_data:
            combined_file = checkpoint_data.get('combined_file')
            if combined_file and os.path.exists(combined_file):
                self.logger.info(f"使用缓存的合并文件: {combined_file}")
                return combined_file
        
        combined_file = os.path.join(self.config.output_dir, "augustus_combined.gff3")
        augustus_files = []
        
        # 优先使用已转换坐标的Augustus文件
        converted_dir = self.config.augustus_predictions.replace("augustus_nlr_candidates", "augustus_nlr_candidates_converted")
        
        if os.path.exists(converted_dir):
            self.logger.info(f"使用已转换坐标的Augustus预测结果: {converted_dir}")
            augustus_dir = Path(converted_dir)
        else:
            self.logger.warning(f"未找到转换坐标目录，使用原始目录: {self.config.augustus_predictions}")
            augustus_dir = Path(self.config.augustus_predictions)
            
        # 收集所有Augustus GFF3文件
        for gff_file in augustus_dir.glob("*.gff3"):
            augustus_files.append(str(gff_file))
            
        self.logger.info(f"找到 {len(augustus_files)} 个Augustus预测文件")
        
        if not augustus_files:
            raise ValueError("未找到Augustus预测文件")
            
        # 合并文件
        gene_count = 0
        with open(combined_file, 'w') as outf:
            outf.write("##gff-version 3\n")
            
            for gff_file in augustus_files:
                with open(gff_file, 'r') as inf:
                    for line in inf:
                        line = line.strip()
                        if line and not line.startswith('#'):
                            outf.write(line + '\n')
                            if '\tgene\t' in line:
                                gene_count += 1
                                
        self.stats.augustus_genes_used = gene_count
        self.logger.info(f"合并完成，共 {gene_count} 个Augustus基因预测")
        
        # 创建检查点
        self._create_checkpoint("augustus_combine", {"combined_file": combined_file, "gene_count": gene_count})
        
        return combined_file

    def _prepare_evm_inputs(self) -> Tuple[str, str, str]:
        """准备EVM输入文件（支持处理后的miniprot结果）"""
        self.logger.info("准备EVM输入文件...")
        
        # 检查是否有检查点
        checkpoint_data = self._load_checkpoint("evm_inputs")
        if checkpoint_data:
            augustus_file = checkpoint_data.get('augustus_file')
            miniprot_file = checkpoint_data.get('miniprot_file')
            weights_file = checkpoint_data.get('weights_file')
            
            if all(os.path.exists(f) for f in [augustus_file, miniprot_file, weights_file]):
                self.logger.info("使用缓存的EVM输入文件")
                return augustus_file, miniprot_file, weights_file
        
        # 1. 准备Augustus输入文件
        augustus_file = self._combine_augustus_predictions()
        
        # 2. 准备Miniprot输入文件
        if self.config.use_processed_miniprot:
            miniprot_file = self._combine_processed_miniprot_predictions()
        else:
            # 使用原始miniprot文件
            miniprot_file = self.config.miniprot_predictions
            
        # 3. 准备权重配置文件
        weights_file = self._generate_combined_weights()
        
        # 创建检查点
        self._create_checkpoint("evm_inputs", {
            "augustus_file": augustus_file,
            "miniprot_file": miniprot_file,
            "weights_file": weights_file
        })
        
        return augustus_file, miniprot_file, weights_file

    def _combine_processed_miniprot_predictions(self) -> str:
        """合并处理后的miniprot预测文件"""
        self.logger.info("合并处理后的miniprot预测文件...")
        
        processed_dir = Path(self.config.miniprot_processed_dir)
        combined_file = os.path.join(self.config.output_dir, "miniprot_combined.gff3")
        
        # 定义要合并的质量级别
        quality_levels = ["high_quality", "medium_quality"]
        if self.config.include_low_quality:
            quality_levels.append("low_quality")
            
        total_predictions = 0
        
        with open(combined_file, 'w') as outf:
            outf.write("##gff-version 3\n")
            outf.write("# Combined miniprot predictions from quality-filtered results\n")
            outf.write(f"# Quality levels included: {', '.join(quality_levels)}\n")
            
            for quality in quality_levels:
                input_file = processed_dir / f"miniprot_{quality}.gff3"
                
                if not input_file.exists():
                    self.logger.warning(f"质量文件不存在: {input_file}")
                    continue
                    
                self.logger.info(f"合并 {quality} 质量预测...")
                file_predictions = 0
                
                with open(input_file, 'r') as inf:
                    for line in inf:
                        line = line.strip()
                        if line and not line.startswith('#'):
                            outf.write(line + '\n')
                            if '\tmRNA\t' in line:
                                file_predictions += 1
                                
                total_predictions += file_predictions
                self.logger.info(f"{quality}: {file_predictions} 个预测")
                
        self.stats.miniprot_genes_used = total_predictions
        self.logger.info(f"Miniprot合并完成，共 {total_predictions} 个预测")
        
        return combined_file

    def _generate_combined_weights(self) -> str:
        """生成合并的权重配置文件"""
        self.logger.info("生成合并的权重配置文件...")
        
        weights_file = os.path.join(self.config.output_dir, "evm_weights.txt")
        
        with open(weights_file, 'w') as f:
            f.write("# EVM Combined Weight Configuration\n")
            f.write(f"# Generated: {datetime.now().isoformat()}\n")
            f.write("#\n")
            
            # Augustus权重
            f.write("# Augustus ab initio predictions\n")
            f.write("ABINITIO_PREDICTION\taugustus\t1.0\n")
            f.write("#\n")
            
            # Miniprot权重（基于质量）
            if self.config.use_processed_miniprot:
                f.write("# Miniprot protein alignments (quality-based weights)\n")
                f.write("PROTEIN\tminiprot_high\t1.0\n")
                f.write("PROTEIN\tminiprot_medium\t0.7\n")
                if self.config.include_low_quality:
                    f.write("PROTEIN\tminiprot_low\t0.4\n")
            else:
                f.write("# Miniprot protein alignments (default weight)\n")
                f.write("PROTEIN\tminiprot\t0.8\n")
                
        self.logger.info(f"权重配置文件生成: {weights_file}")
        return weights_file

    def run_evm_pipeline(self) -> Dict[str, Any]:
        """运行完整的EVM流程"""
        self.logger.info("开始运行EVM流程...")
        start_time = datetime.now()
        
        try:
            # 设置临时目录
            self._setup_temp_directory()
            
            # 检查miniprot处理状态
            if self.config.use_processed_miniprot:
                miniprot_status = self._check_miniprot_processing_status()
                if not miniprot_status["processed"]:
                    raise ValueError("Miniprot预处理未完成，请先运行MiniprotResultProcessor")
            
            # 准备EVM输入文件
            augustus_file, miniprot_file, weights_file = self._prepare_evm_inputs()
            
            # 检查是否使用并行处理
            if self.config.partition_size > 0:
                # 使用内存管理和并行处理
                evm_result = self._run_parallel_evm_pipeline(augustus_file, miniprot_file, weights_file)
            else:
                # 使用原始的单进程处理
                evm_result = self._run_evm_command(augustus_file, miniprot_file, weights_file)
            
            # 分析结果
            if evm_result.get('success', False):
                output_file = evm_result.get('output_file')
                if output_file and os.path.exists(output_file):
                    self._analyze_results(output_file)
                    
            # 计算处理时间
            end_time = datetime.now()
            self.stats.processing_time = (end_time - start_time).total_seconds()
            
            # 生成最终报告
            report = self._generate_report()
            
            # 集成miniprot质量报告
            if self.config.use_processed_miniprot:
                miniprot_status = self._check_miniprot_processing_status()
                report["miniprot_quality_info"] = miniprot_status["quality_report"]
            
            return report
            
        except Exception as e:
            self.logger.error(f"EVM流程执行失败: {e}")
            raise
        finally:
            # 清理临时目录
            self._cleanup_temp_directory()

    def _run_parallel_evm_pipeline(self, augustus_file: str, miniprot_file: str, weights_file: str) -> Dict[str, Any]:
        """使用GenomeMemoryManager运行并行EVM流程"""
        self.logger.info("开始并行EVM处理...")
        
        try:
            # 配置分区参数
            partition_config = PartitioningConfig(
                target_partition_size=self.config.partition_size,
                max_partition_size=self.config.partition_size * 2,
                min_partition_size=self.config.partition_size // 2,
                max_concurrent_partitions=min(self.config.threads, 4),  # 限制并发数
                memory_safety_factor=0.7  # 使用70%内存安全系数
            )
            
            # 创建输出目录用于分区处理
            partition_output_dir = Path(self.config.output_dir) / "partitions"
            partition_output_dir.mkdir(exist_ok=True)
            
            # 使用GenomeMemoryManager进行分区处理
            with GenomeMemoryManager(
                self.config.genome_fasta,
                str(partition_output_dir),
                partition_config
            ) as memory_manager:
                
                # 获取内存统计
                memory_stats = memory_manager.get_memory_stats()
                self.logger.info(f"系统内存状态: {memory_stats.memory_percent:.1f}% 已使用, "
                               f"{memory_stats.available_gb:.1f} GB 可用")
                
                # 估算特征密度
                gff_files = [augustus_file, miniprot_file]
                feature_density = memory_manager.estimate_feature_density(gff_files)
                self.logger.info(f"特征密度估算完成，涉及 {len(feature_density)} 个染色体")
                
                # 创建基因组分区
                partitions = memory_manager.create_partitions(feature_density)
                self.stats.partition_count = len(partitions)
                self.logger.info(f"创建了 {len(partitions)} 个分区")
                
                # 估算处理内存需求
                memory_estimate = memory_manager.estimate_processing_memory(
                    len(partitions), 
                    500  # 估算每个分区需要500MB
                )
                
                if not memory_estimate["memory_sufficient"]:
                    self.logger.warning("内存可能不足，将调整并发数")
                    partition_config.max_concurrent_partitions = max(1, 
                        memory_estimate["safe_concurrent_partitions"])
                
                # 并行处理分区
                results = self._process_partitions_parallel(
                    partitions, augustus_file, miniprot_file, 
                    weights_file, memory_manager
                )
                
                # 合并分区结果
                final_output = self._merge_partition_results(results)
                
                # 更新统计信息
                self.stats.success_partitions = sum(1 for r in results if r.get("success", False))
                self.stats.failed_partitions = len(results) - self.stats.success_partitions
                
                self.logger.info(f"并行处理完成: {self.stats.success_partitions}/{len(partitions)} 分区成功")
                
                return {
                    "success": self.stats.failed_partitions == 0,
                    "output_file": final_output,
                    "partition_results": results,
                    "memory_stats": memory_estimate
                }
                
        except Exception as e:
            self.logger.error(f"并行EVM处理失败: {e}")
            return {"success": False, "error": str(e)}

    def _process_partitions_parallel(self, partitions: List[Any], augustus_file: str, 
                                   miniprot_file: str, weights_file: str, 
                                   memory_manager: GenomeMemoryManager) -> List[Dict[str, Any]]:
        """并行处理基因组分区"""
        self.logger.info(f"开始并行处理 {len(partitions)} 个分区...")
        
        results = []
        max_workers = min(self.config.threads, 4)  # 限制最大并发数
        
        # 创建分区处理任务
        def process_single_partition(partition_idx: int, partition: Any) -> Dict[str, Any]:
            """处理单个分区的函数"""
            try:
                partition_id = f"partition_{partition_idx:04d}"
                self.logger.info(f"开始处理分区 {partition_id}: {partition.chromosome}:{partition.start}-{partition.end}")
                
                # 提取分区序列
                partition_fasta = memory_manager.extract_partition_sequence(partition)
                
                # 提取分区对应的注释数据
                partition_augustus = self._extract_partition_annotations(
                    augustus_file, partition, f"{partition_id}_augustus.gff3"
                )
                partition_miniprot = self._extract_partition_annotations(
                    miniprot_file, partition, f"{partition_id}_miniprot.gff3"
                )
                
                # 运行EVM处理分区
                partition_result = self._run_evm_on_partition(
                    partition_fasta, partition_augustus, partition_miniprot, 
                    weights_file, partition_id
                )
                
                # 清理临时文件
                self._cleanup_partition_temp_files([partition_fasta, partition_augustus, partition_miniprot])
                
                return {
                    "partition_id": partition_id,
                    "partition": partition,
                    "success": partition_result.get("success", False),
                    "output_file": partition_result.get("output_file"),
                    "error": partition_result.get("error")
                }
                
            except Exception as e:
                self.logger.error(f"分区 {partition_idx} 处理失败: {e}")
                return {
                    "partition_id": f"partition_{partition_idx:04d}",
                    "partition": partition,
                    "success": False,
                    "error": str(e)
                }
        
        # 使用ProcessPoolExecutor进行真正的并行处理
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            # 提交所有任务
            future_to_partition = {
                executor.submit(process_single_partition, idx, partition): idx 
                for idx, partition in enumerate(partitions)
            }
            
            # 收集结果
            for future in as_completed(future_to_partition):
                partition_idx = future_to_partition[future]
                try:
                    result = future.result()
                    results.append(result)
                    
                    if result["success"]:
                        self.logger.info(f"分区 {result['partition_id']} 处理成功")
                    else:
                        self.logger.error(f"分区 {result['partition_id']} 处理失败: {result.get('error', 'Unknown error')}")
                        
                except Exception as e:
                    self.logger.error(f"分区 {partition_idx} 执行异常: {e}")
                    results.append({
                        "partition_id": f"partition_{partition_idx:04d}",
                        "success": False,
                        "error": str(e)
                    })
        
        # 按分区顺序排序结果
        results.sort(key=lambda x: x["partition_id"])
        
        return results

    def _extract_partition_annotations(self, gff_file: str, partition: Any, output_filename: str) -> str:
        """提取分区对应的注释数据"""
        output_path = Path(self.config.output_dir) / "partitions" / output_filename
        
        try:
            with open(gff_file, 'r') as infile, open(output_path, 'w') as outfile:
                for line in infile:
                    if line.startswith('#'):
                        outfile.write(line)
                        continue
                    
                    fields = line.strip().split('\t')
                    if len(fields) >= 5:
                        chrom = fields[0]
                        start = int(fields[3])
                        end = int(fields[4])
                        
                        # 检查是否在分区范围内
                        if (chrom == partition.chromosome and 
                            start <= partition.end and end >= partition.start):
                            outfile.write(line)
            
            return str(output_path)
            
        except Exception as e:
            self.logger.error(f"提取分区注释失败: {e}")
            return str(output_path)  # 返回空文件路径

    def _run_evm_on_partition(self, partition_fasta: str, augustus_gff: str, 
                             miniprot_gff: str, weights_file: str, partition_id: str) -> Dict[str, Any]:
        """在单个分区上运行EVM"""
        try:
            output_file = Path(self.config.output_dir) / "partitions" / f"{partition_id}_evm.gff3"
            
            # 构建EVM命令
            cmd = [
                self.config.evm_binary,
                "--genome", partition_fasta,
                "--gene_predictions", augustus_gff,
                "--protein_alignments", miniprot_gff,
                "--weights", weights_file,
                "--output_file_name", str(output_file)
            ]
            
            # 执行EVM
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=3600  # 1小时超时
            )
            
            success = result.returncode == 0
            
            return {
                "success": success,
                "output_file": str(output_file) if success else None,
                "returncode": result.returncode,
                "stderr": result.stderr
            }
            
        except Exception as e:
            return {"success": False, "error": str(e)}

    def _merge_partition_results(self, results: List[Dict[str, Any]]) -> str:
        """合并分区处理结果"""
        self.logger.info("合并分区处理结果...")
        
        final_output = Path(self.config.output_dir) / "evm_output.gff3"
        
        try:
            with open(final_output, 'w') as outfile:
                # 写入GFF3头部
                outfile.write("##gff-version 3\n")
                
                # 合并所有成功的分区结果
                for result in results:
                    if result["success"] and result.get("output_file"):
                        output_file = result["output_file"]
                        if os.path.exists(output_file):
                            with open(output_file, 'r') as infile:
                                for line in infile:
                                    if not line.startswith('#'):  # 跳过头部信息
                                        outfile.write(line)
            
            self.logger.info(f"结果合并完成: {final_output}")
            return str(final_output)
            
        except Exception as e:
            self.logger.error(f"合并结果失败: {e}")
            return str(final_output)

    def _cleanup_partition_temp_files(self, file_paths: List[str]):
        """清理分区临时文件"""
        for file_path in file_paths:
            try:
                if os.path.exists(file_path):
                    os.remove(file_path)
            except Exception as e:
                self.logger.warning(f"清理临时文件失败 {file_path}: {e}")

    def _run_evm_command(self, augustus_file: str, miniprot_file: str, weights_file: str) -> Dict[str, Any]:
        """执行EVM命令"""
        self.logger.info("执行EVM命令...")
        
        output_file = os.path.join(self.config.output_dir, "evm_output.gff3")
        
        # 构建EVM命令
        cmd = [
            self.config.evm_binary,
            "--genome", self.config.genome_fasta,
            "--gene_predictions", augustus_file,
            "--protein_alignments", miniprot_file,
            "--weights", weights_file,
            "--output_file_name", output_file,
            "--segmentSize", str(self.config.partition_size),
            "--CPU", str(self.config.threads)
        ]
        
        self.logger.info(f"EVM命令: {' '.join(cmd)}")
        
        try:
            # 执行命令
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=7200  # 2小时超时
            )
            
            success = result.returncode == 0
            
            if success:
                self.logger.info("EVM执行成功")
            else:
                self.logger.error(f"EVM执行失败: {result.stderr}")
                
            return {
                "success": success,
                "output_file": output_file if success else None,
                "stdout": result.stdout,
                "stderr": result.stderr,
                "returncode": result.returncode
            }
            
        except subprocess.TimeoutExpired:
            self.logger.error("EVM执行超时")
            return {"success": False, "error": "Timeout"}
        except Exception as e:
            self.logger.error(f"EVM执行异常: {e}")
            return {"success": False, "error": str(e)}

    def _analyze_results(self, output_file: str):
        """分析EVM输出结果"""
        self.logger.info("分析EVM输出结果...")
        
        gene_count = 0
        exon_count = 0
        intron_count = 0
        
        try:
            with open(output_file, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                        
                    fields = line.strip().split('\t')
                    if len(fields) >= 3:
                        feature_type = fields[2]
                        
                        if feature_type == 'gene':
                            gene_count += 1
                        elif feature_type == 'exon':
                            exon_count += 1
                        elif feature_type == 'intron':
                            intron_count += 1
                            
            self.stats.total_genes_predicted = gene_count
            self.stats.total_exons = exon_count
            self.stats.total_introns = intron_count
            
            self.logger.info(f"EVM结果统计 - 基因: {gene_count}, 外显子: {exon_count}, 内含子: {intron_count}")
            
        except Exception as e:
            self.logger.error(f"分析结果失败: {e}")


    def _generate_report(self) -> Dict[str, Any]:
        """生成最终报告"""
        report = {
            "timestamp": datetime.now().isoformat(),
            "config": asdict(self.config),
            "statistics": asdict(self.stats),
            "checkpoints": self.checkpoints,
            "output_files": {
                "evm_output": f"{self.config.output_dir}/evm_output.gff3",
                "augustus_combined": f"{self.config.output_dir}/augustus_combined.gff3",
                "miniprot_combined": f"{self.config.output_dir}/miniprot_combined.gff3",
                "weights_file": f"{self.config.output_dir}/weights.txt"
            }
        }
        
        # 保存报告
        report_file = Path(self.config.output_dir) / "evm_integration_report.json"
        with open(report_file, 'w') as f:
            json.dump(report, f, indent=2)
            
        self.logger.info(f"报告已保存到: {report_file}")
        return report

    # ============== 新增：内存管理和分区策略方法 ==============
    
    def _get_system_memory_info(self) -> Dict[str, float]:
        """获取系统内存信息"""
        try:
            import psutil
            memory = psutil.virtual_memory()
            return {
                "total_gb": memory.total / (1024**3),
                "available_gb": memory.available / (1024**3),
                "used_gb": memory.used / (1024**3),
                "percent_used": memory.percent
            }
        except ImportError:
            self.logger.warning("psutil未安装，使用默认内存估算")
            return {
                "total_gb": 16.0,  # 默认假设16GB内存
                "available_gb": 8.0,
                "used_gb": 8.0,
                "percent_used": 50.0
            }
    
    def _calculate_optimal_partition_size(self) -> int:
        """
        基于基因组大小和可用内存计算最优分区大小
        
        Returns:
            int: 最优分区大小（bp）
        """
        self.logger.info("计算最优分区大小...")
        
        # 获取系统内存信息
        memory_info = self._get_system_memory_info()
        available_memory_gb = memory_info["available_gb"]
        
        # 获取基因组大小
        genome_size = self._get_genome_size()
        genome_size_gb = genome_size / (1024**3)
        
        self.logger.info(f"基因组大小: {genome_size_gb:.2f} GB")
        self.logger.info(f"可用内存: {available_memory_gb:.2f} GB")
        
        # 计算最优分区大小
        # 策略：确保每个分区在内存中的数据量不超过可用内存的20%
        max_partition_memory_gb = available_memory_gb * 0.2
        
        # 估算每bp数据的内存使用量（包括GFF3数据）
        # 经验值：每bp大约需要100-200字节内存（包括索引和处理开销）
        bytes_per_bp = 150
        max_partition_bp = int(max_partition_memory_gb * (1024**3) / bytes_per_bp)
        
        # 设置合理的分区大小范围
        min_partition_size = 500_000   # 最小500kb
        max_partition_size = 50_000_000  # 最大50Mb
        
        optimal_size = max(min_partition_size, min(max_partition_bp, max_partition_size))
        
        self.logger.info(f"计算得出最优分区大小: {optimal_size:,} bp ({optimal_size/1_000_000:.1f} Mb)")
        
        return optimal_size
    
    def _get_genome_size(self) -> int:
        """获取基因组总大小"""
        try:
            total_size = 0
            with open(self.config.genome_fasta, 'r') as f:
                for line in f:
                    if not line.startswith('>'):
                        total_size += len(line.strip())
            return total_size
        except Exception as e:
            self.logger.error(f"无法获取基因组大小: {e}")
            return 3_400_000_000  # 默认3.4Gb（水稻基因组大小）
    
    def _monitor_memory_usage(self) -> Dict[str, float]:
        """监控当前内存使用情况"""
        try:
            import psutil
            process = psutil.Process()
            memory_info = process.memory_info()
            system_memory = psutil.virtual_memory()
            
            return {
                "process_memory_mb": memory_info.rss / (1024**2),
                "process_memory_percent": process.memory_percent(),
                "system_available_gb": system_memory.available / (1024**3),
                "system_used_percent": system_memory.percent
            }
        except ImportError:
            return {
                "process_memory_mb": 0,
                "process_memory_percent": 0,
                "system_available_gb": 8.0,
                "system_used_percent": 50.0
            }
    
    def _create_genome_partitions(self, partition_size: int) -> List[Dict[str, Any]]:
        """
        创建基因组分区信息
        
        Args:
            partition_size: 分区大小（bp）
            
        Returns:
            List[Dict]: 分区信息列表
        """
        self.logger.info("创建基因组分区...")
        
        partitions = []
        
        # 读取基因组序列信息
        chromosomes = self._get_chromosome_info()
        
        partition_id = 0
        for chr_name, chr_length in chromosomes.items():
            # 为每个染色体创建分区
            start_pos = 1
            while start_pos <= chr_length:
                end_pos = min(start_pos + partition_size - 1, chr_length)
                
                partition = {
                    "id": partition_id,
                    "chromosome": chr_name,
                    "start": start_pos,
                    "end": end_pos,
                    "length": end_pos - start_pos + 1,
                    "partition_name": f"{chr_name}_{start_pos}-{end_pos}"
                }
                
                partitions.append(partition)
                partition_id += 1
                start_pos = end_pos + 1
        
        self.logger.info(f"创建了 {len(partitions)} 个分区")
        self.stats.partition_count = len(partitions)
        
        return partitions
    
    def _get_chromosome_info(self) -> Dict[str, int]:
        """获取染色体信息"""
        chromosomes = {}
        current_chr = None
        current_length = 0
        
        try:
            with open(self.config.genome_fasta, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        # 保存前一个染色体的长度
                        if current_chr:
                            chromosomes[current_chr] = current_length
                        
                        # 开始新染色体
                        current_chr = line[1:].split()[0]  # 取第一个空格前的部分作为染色体名
                        current_length = 0
                    else:
                        current_length += len(line)
                
                # 保存最后一个染色体
                if current_chr:
                    chromosomes[current_chr] = current_length
                    
        except Exception as e:
            self.logger.error(f"读取基因组文件失败: {e}")
            raise
        
        self.logger.info(f"检测到 {len(chromosomes)} 个染色体")
        for chr_name, length in chromosomes.items():
            self.logger.debug(f"{chr_name}: {length:,} bp")
            
        return chromosomes
    
    def _stream_process_partitions(self, partitions: List[Dict[str, Any]], 
                                 augustus_file: str, miniprot_file: str, 
                                 weights_file: str) -> Dict[str, Any]:
        """
        流式处理分区数据
        
        Args:
            partitions: 分区信息列表
            augustus_file: Augustus预测文件
            miniprot_file: Miniprot预测文件
            weights_file: 权重配置文件
            
        Returns:
            Dict: 处理结果统计
        """
        self.logger.info(f"开始流式处理 {len(partitions)} 个分区...")
        
        results = {
            "total_partitions": len(partitions),
            "successful_partitions": 0,
            "failed_partitions": 0,
            "total_genes": 0,
            "processing_times": [],
            "memory_usage": []
        }
        
        # 创建分区输出目录
        partition_dir = Path(self.config.output_dir) / "partitions"
        partition_dir.mkdir(exist_ok=True)
        
        for i, partition in enumerate(partitions):
            start_time = datetime.now()
            partition_name = partition["partition_name"]
            
            self.logger.info(f"处理分区 {i+1}/{len(partitions)}: {partition_name}")
            
            try:
                # 监控内存使用
                memory_before = self._monitor_memory_usage()
                
                # 提取分区数据
                partition_files = self._extract_partition_data(
                    partition, augustus_file, miniprot_file
                )
                
                # 运行EVM处理该分区
                partition_result = self._run_evm_partition(
                    partition, partition_files, weights_file
                )
                
                # 记录成功处理
                results["successful_partitions"] += 1
                results["total_genes"] += partition_result.get("gene_count", 0)
                
                # 监控内存使用
                memory_after = self._monitor_memory_usage()
                results["memory_usage"].append({
                    "partition": partition_name,
                    "before_mb": memory_before["process_memory_mb"],
                    "after_mb": memory_after["process_memory_mb"],
                    "system_available_gb": memory_after["system_available_gb"]
                })
                
                # 清理临时文件
                self._cleanup_partition_files(partition_files)
                
            except Exception as e:
                self.logger.error(f"处理分区 {partition_name} 失败: {e}")
                results["failed_partitions"] += 1
                
            # 记录处理时间
            processing_time = (datetime.now() - start_time).total_seconds()
            results["processing_times"].append(processing_time)
            
            self.logger.info(f"分区 {partition_name} 处理完成，耗时: {processing_time:.2f}s")
        
        self.logger.info(f"分区处理完成 - 成功: {results['successful_partitions']}, "
                        f"失败: {results['failed_partitions']}")
        
        return results
    
    def _extract_partition_data(self, partition: Dict[str, Any], 
                              augustus_file: str, miniprot_file: str) -> Dict[str, str]:
        """
        提取分区相关的数据
        
        Args:
            partition: 分区信息
            augustus_file: Augustus预测文件
            miniprot_file: Miniprot预测文件
            
        Returns:
            Dict: 分区数据文件路径
        """
        chr_name = partition["chromosome"]
        start_pos = partition["start"]
        end_pos = partition["end"]
        partition_name = partition["partition_name"]
        
        # 创建分区临时目录
        partition_temp_dir = Path(self.temp_dir) / partition_name
        partition_temp_dir.mkdir(exist_ok=True)
        
        # 提取基因组序列
        genome_file = partition_temp_dir / f"{partition_name}.fa"
        self._extract_genome_sequence(chr_name, start_pos, end_pos, genome_file)
        
        # 提取Augustus预测
        augustus_partition_file = partition_temp_dir / f"{partition_name}_augustus.gff3"
        self._extract_gff_partition(augustus_file, chr_name, start_pos, end_pos, 
                                   augustus_partition_file)
        
        # 提取Miniprot预测
        miniprot_partition_file = partition_temp_dir / f"{partition_name}_miniprot.gff3"
        self._extract_gff_partition(miniprot_file, chr_name, start_pos, end_pos,
                                   miniprot_partition_file)
        
        return {
            "genome": str(genome_file),
            "augustus": str(augustus_partition_file),
            "miniprot": str(miniprot_partition_file),
            "temp_dir": str(partition_temp_dir)
        }
    
    def _extract_genome_sequence(self, chr_name: str, start_pos: int, end_pos: int, 
                               output_file: Path):
        """提取基因组序列片段"""
        # 这里可以使用samtools faidx或者直接读取FASTA文件
        # 为简化实现，先创建一个符号链接到原始基因组文件
        # 在实际应用中应该提取具体的序列片段
        try:
            os.symlink(os.path.abspath(self.config.genome_fasta), output_file)
        except OSError:
            # 如果符号链接失败，复制文件
            shutil.copy2(self.config.genome_fasta, output_file)
    
    def _extract_gff_partition(self, input_file: str, chr_name: str, 
                             start_pos: int, end_pos: int, output_file: Path):
        """提取GFF3文件中指定区域的记录"""
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            for line in infile:
                if line.startswith('#'):
                    outfile.write(line)
                    continue
                    
                fields = line.strip().split('\t')
                if len(fields) >= 5:
                    seq_id = fields[0]
                    feature_start = int(fields[3])
                    feature_end = int(fields[4])
                    
                    # 检查是否在目标区域内
                    if (seq_id == chr_name and 
                        not (feature_end < start_pos or feature_start > end_pos)):
                        outfile.write(line)
    
    def _run_evm_partition(self, partition: Dict[str, Any], 
                          partition_files: Dict[str, str], 
                          weights_file: str) -> Dict[str, Any]:
        """运行单个分区的EVM处理"""
        partition_name = partition["partition_name"]
        
        # 构建EVM命令
        cmd = [
            self.config.evm_binary,
            "--genome", partition_files["genome"],
            "--gene_predictions", partition_files["augustus"],
            "--protein_alignments", partition_files["miniprot"],
            "--weights", weights_file,
            "--output_file_name", f"{partition_name}_evm"
        ]
        
        # 执行命令
        result = subprocess.run(cmd, capture_output=True, text=True, 
                              cwd=partition_files["temp_dir"])
        
        if result.returncode != 0:
            raise RuntimeError(f"EVM执行失败: {result.stderr}")
        
        # 分析结果
        output_file = Path(partition_files["temp_dir"]) / f"{partition_name}_evm.gff3"
        gene_count = 0
        if output_file.exists():
            with open(output_file, 'r') as f:
                for line in f:
                    if '\tgene\t' in line:
                        gene_count += 1
        
        return {
            "gene_count": gene_count,
            "output_file": str(output_file)
        }
    
    def _cleanup_partition_files(self, partition_files: Dict[str, str]):
        """清理分区临时文件"""
        temp_dir = partition_files.get("temp_dir")
        if temp_dir and os.path.exists(temp_dir):
            try:
                shutil.rmtree(temp_dir)
            except Exception as e:
                self.logger.warning(f"清理临时目录失败: {e}")
    
    def _manage_temp_files(self):
        """管理临时文件生命周期"""
        if not self.temp_dir:
            return
            
        try:
            # 获取临时目录大小
            total_size = 0
            for dirpath, dirnames, filenames in os.walk(self.temp_dir):
                for filename in filenames:
                    filepath = os.path.join(dirpath, filename)
                    total_size += os.path.getsize(filepath)
            
            size_gb = total_size / (1024**3)
            self.logger.info(f"临时文件总大小: {size_gb:.2f} GB")
            
            # 如果临时文件过大，进行清理
            if size_gb > 10.0:  # 超过10GB
                self.logger.warning("临时文件过大，开始清理...")
                self._cleanup_old_temp_files()
                
        except Exception as e:
            self.logger.error(f"管理临时文件失败: {e}")
    
    def _cleanup_old_temp_files(self):
        """清理旧的临时文件"""
        if not self.temp_dir or not os.path.exists(self.temp_dir):
            return
            
        try:
            import time
            current_time = time.time()
            
            for dirpath, dirnames, filenames in os.walk(self.temp_dir):
                for filename in filenames:
                    filepath = os.path.join(dirpath, filename)
                    file_age = current_time - os.path.getmtime(filepath)
                    
                    # 删除超过1小时的临时文件
                    if file_age > 3600:
                        os.remove(filepath)
                        self.logger.debug(f"删除旧临时文件: {filepath}")
                        
        except Exception as e:
            self.logger.error(f"清理旧临时文件失败: {e}")


def main():
    """主函数"""
    # 配置日志
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    
    # EVM配置
    config = EVMConfig(
        genome_fasta="genome/osa.fa",
        augustus_predictions="results/augustus_nlr_candidates_converted",
        miniprot_predictions="results/miniprot_full_predictions.gff",
        weights_file="evm_config/weights.txt",
        output_dir="results/evm_integration",
        evm_binary="tools/EVidenceModeler/EVidenceModeler",
        use_processed_miniprot=True,
        miniprot_processed_dir="results/miniprot_processed"
    )
    
    # 运行EVM
    runner = EVMRunner(config)
    result = runner.run_evm_pipeline()
    
    print("EVM集成完成!")
    print(f"处理结果: {json.dumps(result, indent=2)}")


if __name__ == "__main__":
    main() 