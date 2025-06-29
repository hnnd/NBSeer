"""
数据验证和基因组质量检查模块
Data Validation and Genome Quality Assessment Module

提供全面的基因组数据质量评估和验证功能。
"""

import os
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any, Union
from dataclasses import dataclass, asdict
from collections import Counter
import logging

from .file_io import FastaValidator, FastaReader, SequenceStats
from .logging_setup import get_logger


@dataclass
class QualityThresholds:
    """质量评估阈值配置"""
    min_sequence_length: int = 1000
    max_n_content: float = 0.1  # 10%
    min_gc_content: float = 0.2  # 20%
    max_gc_content: float = 0.8  # 80%
    min_genome_size: int = 1_000_000  # 1MB
    max_genome_size: int = 50_000_000_000  # 50GB
    min_n50: int = 10000
    min_sequences_for_assembly: int = 1


@dataclass
class AssemblyStats:
    """基因组组装统计信息"""
    n50: int
    n90: int
    l50: int  # 达到N50所需的最少contig数
    l90: int  # 达到N90所需的最少contig数
    largest_contig: int
    total_gaps: int
    gap_percentage: float
    contig_count: int


@dataclass
class QualityMetrics:
    """基因组质量指标"""
    sequence_stats: SequenceStats
    assembly_stats: AssemblyStats
    gc_distribution: Dict[str, float]
    n_content_distribution: Dict[str, float]
    sequence_complexity: Dict[str, float]
    quality_score: float
    quality_grade: str
    issues: List[str]
    recommendations: List[str]


class GenomeQualityAssessment:
    """基因组质量评估器"""
    
    def __init__(self, thresholds: Optional[QualityThresholds] = None, 
                 logger: Optional[logging.Logger] = None):
        self.thresholds = thresholds or QualityThresholds()
        self.logger = logger or get_logger(__name__)
        self.fasta_validator = FastaValidator(self.logger)
        self.fasta_reader = FastaReader(self.logger)
    
    def assess_genome_quality(self, file_path: Union[str, Path]) -> Optional[QualityMetrics]:
        """
        全面评估基因组质量
        
        Args:
            file_path: 基因组文件路径
            
        Returns:
            QualityMetrics对象或None
        """
        file_path = Path(file_path)
        self.logger.info(f"开始基因组质量评估: {file_path}")
        
        # 首先验证文件格式
        is_valid, error_msg, seq_count = self.fasta_validator.validate_fasta_format(file_path)
        if not is_valid:
            self.logger.error(f"文件格式验证失败: {error_msg}")
            return None
        
        # 获取基本序列统计
        sequence_stats = self.fasta_validator.get_sequence_statistics(file_path)
        if not sequence_stats:
            self.logger.error("无法获取序列统计信息")
            return None
        
        try:
            # 计算组装统计
            assembly_stats = self._calculate_assembly_stats(file_path)
            
            # 计算GC含量分布
            gc_distribution = self._calculate_gc_distribution(file_path)
            
            # 计算N含量分布
            n_content_distribution = self._calculate_n_content_distribution(file_path)
            
            # 计算序列复杂度
            sequence_complexity = self._calculate_sequence_complexity(file_path)
            
            # 评估质量分数和等级
            quality_score, quality_grade, issues, recommendations = self._evaluate_quality(
                sequence_stats, assembly_stats, gc_distribution, n_content_distribution
            )
            
            quality_metrics = QualityMetrics(
                sequence_stats=sequence_stats,
                assembly_stats=assembly_stats,
                gc_distribution=gc_distribution,
                n_content_distribution=n_content_distribution,
                sequence_complexity=sequence_complexity,
                quality_score=quality_score,
                quality_grade=quality_grade,
                issues=issues,
                recommendations=recommendations
            )
            
            self.logger.info(f"质量评估完成，质量等级: {quality_grade}，得分: {quality_score:.2f}")
            return quality_metrics
            
        except Exception as e:
            self.logger.error(f"质量评估过程中发生错误: {e}")
            return None
    
    def _calculate_assembly_stats(self, file_path: Union[str, Path]) -> AssemblyStats:
        """计算组装统计信息"""
        self.logger.info("计算组装统计信息")
        
        # 收集所有序列长度
        sequence_lengths = []
        total_gaps = 0
        
        for seq_record in self.fasta_reader.read_sequences(file_path):
            seq_len = len(seq_record.seq)
            sequence_lengths.append(seq_len)
            
            # 计算gap数量（N字符）
            seq_str = str(seq_record.seq).upper()
            total_gaps += seq_str.count('N')
        
        if not sequence_lengths:
            return AssemblyStats(0, 0, 0, 0, 0, 0, 0.0, 0)
        
        # 排序（从大到小）
        sequence_lengths.sort(reverse=True)
        
        total_length = sum(sequence_lengths)
        contig_count = len(sequence_lengths)
        largest_contig = sequence_lengths[0] if sequence_lengths else 0
        gap_percentage = total_gaps / total_length if total_length > 0 else 0.0
        
        # 计算N50和N90
        n50, l50 = self._calculate_nx(sequence_lengths, 0.5)
        n90, l90 = self._calculate_nx(sequence_lengths, 0.9)
        
        return AssemblyStats(
            n50=n50,
            n90=n90,
            l50=l50,
            l90=l90,
            largest_contig=largest_contig,
            total_gaps=total_gaps,
            gap_percentage=gap_percentage,
            contig_count=contig_count
        )
    
    def _calculate_nx(self, sorted_lengths: List[int], fraction: float) -> Tuple[int, int]:
        """
        计算Nx统计（如N50, N90）
        
        Args:
            sorted_lengths: 排序后的序列长度列表（降序）
            fraction: 分位数（0.5 for N50, 0.9 for N90）
            
        Returns:
            (Nx值, Lx值)
        """
        total_length = sum(sorted_lengths)
        target_length = total_length * fraction
        
        cumulative_length = 0
        for i, length in enumerate(sorted_lengths):
            cumulative_length += length
            if cumulative_length >= target_length:
                return length, i + 1
        
        return 0, 0
    
    def _calculate_gc_distribution(self, file_path: Union[str, Path]) -> Dict[str, float]:
        """计算GC含量分布"""
        self.logger.info("计算GC含量分布")
        
        gc_contents = []
        
        for seq_record in self.fasta_reader.read_sequences(file_path):
            seq_str = str(seq_record.seq).upper()
            seq_len = len(seq_str)
            
            if seq_len > 0:
                gc_count = seq_str.count('G') + seq_str.count('C')
                gc_content = gc_count / seq_len
                gc_contents.append(gc_content)
        
        if not gc_contents:
            return {"mean": 0.0, "std": 0.0, "min": 0.0, "max": 0.0}
        
        gc_array = np.array(gc_contents)
        
        return {
            "mean": float(np.mean(gc_array)),
            "std": float(np.std(gc_array)),
            "min": float(np.min(gc_array)),
            "max": float(np.max(gc_array)),
            "median": float(np.median(gc_array)),
            "q25": float(np.percentile(gc_array, 25)),
            "q75": float(np.percentile(gc_array, 75))
        }
    
    def _calculate_n_content_distribution(self, file_path: Union[str, Path]) -> Dict[str, float]:
        """计算N含量分布"""
        self.logger.info("计算N含量分布")
        
        n_contents = []
        
        for seq_record in self.fasta_reader.read_sequences(file_path):
            seq_str = str(seq_record.seq).upper()
            seq_len = len(seq_str)
            
            if seq_len > 0:
                n_count = seq_str.count('N')
                n_content = n_count / seq_len
                n_contents.append(n_content)
        
        if not n_contents:
            return {"mean": 0.0, "std": 0.0, "min": 0.0, "max": 0.0}
        
        n_array = np.array(n_contents)
        
        return {
            "mean": float(np.mean(n_array)),
            "std": float(np.std(n_array)),
            "min": float(np.min(n_array)),
            "max": float(np.max(n_array)),
            "median": float(np.median(n_array)),
            "sequences_with_n": sum(1 for n in n_contents if n > 0),
            "total_sequences": len(n_contents)
        }
    
    def _calculate_sequence_complexity(self, file_path: Union[str, Path]) -> Dict[str, float]:
        """计算序列复杂度"""
        self.logger.info("计算序列复杂度")
        
        complexity_scores = []
        low_complexity_count = 0
        total_sequences = 0
        
        for seq_record in self.fasta_reader.read_sequences(file_path):
            total_sequences += 1
            seq_str = str(seq_record.seq).upper()
            
            # 计算Shannon熵作为复杂度指标
            complexity = self._calculate_shannon_entropy(seq_str)
            complexity_scores.append(complexity)
            
            # 低复杂度阈值（Shannon熵 < 1.5）
            if complexity < 1.5:
                low_complexity_count += 1
        
        if not complexity_scores:
            return {"mean_complexity": 0.0, "low_complexity_fraction": 0.0}
        
        return {
            "mean_complexity": float(np.mean(complexity_scores)),
            "std_complexity": float(np.std(complexity_scores)),
            "min_complexity": float(np.min(complexity_scores)),
            "max_complexity": float(np.max(complexity_scores)),
            "low_complexity_sequences": low_complexity_count,
            "low_complexity_fraction": low_complexity_count / total_sequences,
            "total_sequences": total_sequences
        }
    
    def _calculate_shannon_entropy(self, sequence: str) -> float:
        """计算序列的Shannon熵"""
        if not sequence:
            return 0.0
        
        # 计算每个字符的频率
        char_counts = Counter(sequence)
        sequence_length = len(sequence)
        
        # 计算Shannon熵
        entropy = 0.0
        for count in char_counts.values():
            if count > 0:
                probability = count / sequence_length
                entropy -= probability * np.log2(probability)
        
        return entropy
    
    def _evaluate_quality(self, sequence_stats: SequenceStats, assembly_stats: AssemblyStats,
                         gc_distribution: Dict[str, float], n_content_distribution: Dict[str, float]) -> Tuple[float, str, List[str], List[str]]:
        """评估整体质量分数和等级"""
        issues = []
        recommendations = []
        scores = []
        
        # 1. 基因组大小评估 (权重: 15%)
        genome_size_score = self._score_genome_size(sequence_stats.total_length, issues, recommendations)
        scores.append(("genome_size", genome_size_score, 0.15))
        
        # 2. N50评估 (权重: 20%)
        n50_score = self._score_n50(assembly_stats.n50, issues, recommendations)
        scores.append(("n50", n50_score, 0.20))
        
        # 3. N含量评估 (权重: 20%)
        n_content_score = self._score_n_content(sequence_stats.n_content, n_content_distribution, issues, recommendations)
        scores.append(("n_content", n_content_score, 0.20))
        
        # 4. GC含量评估 (权重: 15%)
        gc_score = self._score_gc_content(sequence_stats.gc_content, gc_distribution, issues, recommendations)
        scores.append(("gc_content", gc_score, 0.15))
        
        # 5. 序列数量评估 (权重: 10%)
        sequence_count_score = self._score_sequence_count(sequence_stats.total_sequences, assembly_stats.contig_count, issues, recommendations)
        scores.append(("sequence_count", sequence_count_score, 0.10))
        
        # 6. 序列长度分布评估 (权重: 10%)
        length_distribution_score = self._score_length_distribution(sequence_stats, assembly_stats, issues, recommendations)
        scores.append(("length_distribution", length_distribution_score, 0.10))
        
        # 7. Gap评估 (权重: 10%)
        gap_score = self._score_gaps(assembly_stats.gap_percentage, issues, recommendations)
        scores.append(("gaps", gap_score, 0.10))
        
        # 计算加权总分
        total_score = sum(score * weight for _, score, weight in scores)
        
        # 确定质量等级
        if total_score >= 90:
            quality_grade = "优秀"
        elif total_score >= 80:
            quality_grade = "良好"
        elif total_score >= 60:
            quality_grade = "需改进"
        else:
            quality_grade = "不合格"
        
        return total_score, quality_grade, issues, recommendations
    
    def _score_genome_size(self, genome_size: int, issues: List[str], recommendations: List[str]) -> float:
        """评估基因组大小"""
        if genome_size < self.thresholds.min_genome_size:
            issues.append(f"基因组大小过小: {genome_size:,} bp < {self.thresholds.min_genome_size:,} bp")
            recommendations.append("检查是否为部分基因组或组装不完整")
            return 30.0
        elif genome_size > self.thresholds.max_genome_size:
            issues.append(f"基因组大小异常大: {genome_size:,} bp > {self.thresholds.max_genome_size:,} bp")
            recommendations.append("检查是否包含污染序列或重复组装")
            return 60.0
        else:
            return 100.0
    
    def _score_n50(self, n50: int, issues: List[str], recommendations: List[str]) -> float:
        """评估N50值"""
        if n50 < self.thresholds.min_n50:
            issues.append(f"N50值较低: {n50:,} bp < {self.thresholds.min_n50:,} bp")
            recommendations.append("考虑使用更好的组装策略或长读长测序技术")
            if n50 < 1000:
                return 20.0
            elif n50 < 5000:
                return 50.0
            else:
                return 70.0
        else:
            # N50越高分数越高，但有上限
            if n50 > 1000000:  # 1MB
                return 100.0
            elif n50 > 100000:  # 100KB
                return 95.0
            elif n50 > 50000:   # 50KB
                return 90.0
            else:
                return 85.0
    
    def _score_n_content(self, overall_n_content: float, n_distribution: Dict[str, float], 
                        issues: List[str], recommendations: List[str]) -> float:
        """评估N含量"""
        if overall_n_content > self.thresholds.max_n_content:
            issues.append(f"N含量过高: {overall_n_content:.2%} > {self.thresholds.max_n_content:.2%}")
            recommendations.append("考虑gap填充或重新组装以减少未知碱基")
            if overall_n_content > 0.2:  # 20%
                return 10.0
            elif overall_n_content > 0.15:  # 15%
                return 30.0
            else:
                return 60.0
        else:
            # N含量越低分数越高
            if overall_n_content < 0.01:  # 1%
                return 100.0
            elif overall_n_content < 0.05:  # 5%
                return 90.0
            else:
                return 80.0
    
    def _score_gc_content(self, overall_gc_content: float, gc_distribution: Dict[str, float],
                         issues: List[str], recommendations: List[str]) -> float:
        """评估GC含量"""
        if overall_gc_content < self.thresholds.min_gc_content or overall_gc_content > self.thresholds.max_gc_content:
            issues.append(f"GC含量异常: {overall_gc_content:.2%} (正常范围: {self.thresholds.min_gc_content:.2%}-{self.thresholds.max_gc_content:.2%})")
            recommendations.append("检查是否存在污染序列或特殊的基因组特征")
            return 40.0
        
        # 检查GC含量分布的标准差
        gc_std = gc_distribution.get("std", 0.0)
        if gc_std > 0.1:  # 标准差过大
            issues.append(f"GC含量分布不均匀，标准差: {gc_std:.3f}")
            recommendations.append("检查是否存在不同来源的序列混合")
            return 70.0
        
        return 100.0
    
    def _score_sequence_count(self, total_sequences: int, contig_count: int,
                            issues: List[str], recommendations: List[str]) -> float:
        """评估序列数量"""
        if contig_count < self.thresholds.min_sequences_for_assembly:
            issues.append(f"序列数量过少: {contig_count}")
            recommendations.append("检查是否为单个染色体或高质量组装")
            return 60.0
        elif contig_count > 10000:
            issues.append(f"序列片段过多: {contig_count:,}")
            recommendations.append("考虑使用scaffold或更好的组装方法减少片段化")
            return 40.0
        elif contig_count > 1000:
            return 70.0
        else:
            return 100.0
    
    def _score_length_distribution(self, sequence_stats: SequenceStats, assembly_stats: AssemblyStats,
                                 issues: List[str], recommendations: List[str]) -> float:
        """评估序列长度分布"""
        avg_length = sequence_stats.avg_length
        largest_contig = assembly_stats.largest_contig
        
        # 检查平均长度
        if avg_length < self.thresholds.min_sequence_length:
            issues.append(f"平均序列长度过短: {avg_length:.0f} bp < {self.thresholds.min_sequence_length} bp")
            recommendations.append("考虑过滤短序列或改进组装质量")
            return 50.0
        
        # 检查最大序列与平均序列的比例
        if largest_contig > 0 and avg_length > 0:
            ratio = largest_contig / avg_length
            if ratio > 1000:  # 长度分布极不均匀
                issues.append(f"序列长度分布极不均匀，最大/平均比例: {ratio:.1f}")
                recommendations.append("检查是否存在异常长的序列或组装错误")
                return 60.0
        
        return 100.0
    
    def _score_gaps(self, gap_percentage: float, issues: List[str], recommendations: List[str]) -> float:
        """评估gap比例"""
        if gap_percentage > 0.1:  # 10%
            issues.append(f"Gap比例过高: {gap_percentage:.2%}")
            recommendations.append("考虑gap填充或重新组装")
            return 30.0
        elif gap_percentage > 0.05:  # 5%
            return 70.0
        else:
            return 100.0


class DataValidator:
    """综合数据验证器"""
    
    def __init__(self, config_manager=None, logger: Optional[logging.Logger] = None):
        self.config = config_manager
        self.logger = logger or get_logger(__name__)
        self.quality_assessor = GenomeQualityAssessment(logger=self.logger)
        
    def validate_input_data(self, file_path: Union[str, Path]) -> Dict[str, Any]:
        """
        验证输入数据的完整性和质量
        
        Args:
            file_path: 输入文件路径
            
        Returns:
            验证结果字典
        """
        file_path = Path(file_path)
        self.logger.info(f"开始验证输入数据: {file_path}")
        
        result = {
            "file_path": str(file_path),
            "validation_passed": False,
            "format_valid": False,
            "quality_metrics": None,
            "errors": [],
            "warnings": [],
            "recommendations": []
        }
        
        try:
            # 1. 基本文件检查
            if not file_path.exists():
                result["errors"].append(f"文件不存在: {file_path}")
                return result
            
            # 2. 格式验证
            fasta_validator = FastaValidator(self.logger)
            is_valid, error_msg, seq_count = fasta_validator.validate_fasta_format(file_path)
            
            if not is_valid:
                result["errors"].append(f"FASTA格式验证失败: {error_msg}")
                return result
            
            result["format_valid"] = True
            result["sequence_count"] = seq_count
            
            # 3. 质量评估
            quality_metrics = self.quality_assessor.assess_genome_quality(file_path)
            if quality_metrics:
                result["quality_metrics"] = asdict(quality_metrics)
                
                # 根据质量等级确定是否通过验证
                if quality_metrics.quality_grade in ["优秀", "良好"]:
                    result["validation_passed"] = True
                elif quality_metrics.quality_grade == "需改进":
                    result["validation_passed"] = True
                    result["warnings"].extend(quality_metrics.issues)
                    result["recommendations"].extend(quality_metrics.recommendations)
                else:
                    result["validation_passed"] = False
                    result["errors"].extend(quality_metrics.issues)
                    result["recommendations"].extend(quality_metrics.recommendations)
            
            return result
            
        except Exception as e:
            error_msg = f"数据验证过程中发生错误: {e}"
            self.logger.error(error_msg)
            result["errors"].append(error_msg)
            return result


if __name__ == "__main__":
    # 测试基因组质量评估
    logger = get_logger(__name__)
    assessor = GenomeQualityAssessment(logger=logger)
    
    # 测试文件
    test_file = "data/input/test_valid.fasta"
    
    if Path(test_file).exists():
        print(f"=== 测试基因组质量评估: {test_file} ===")
        
        quality_metrics = assessor.assess_genome_quality(test_file)
        if quality_metrics:
            print(f"质量等级: {quality_metrics.quality_grade}")
            print(f"质量得分: {quality_metrics.quality_score:.2f}")
            print(f"N50: {quality_metrics.assembly_stats.n50:,} bp")
            print(f"GC含量: {quality_metrics.sequence_stats.gc_content:.2%}")
            print(f"N含量: {quality_metrics.sequence_stats.n_content:.2%}")
            
            if quality_metrics.issues:
                print(f"问题: {quality_metrics.issues}")
            if quality_metrics.recommendations:
                print(f"建议: {quality_metrics.recommendations}")
        
        print("\n=== 测试数据验证器 ===")
        validator = DataValidator(logger=logger)
        result = validator.validate_input_data(test_file)
        print(f"验证通过: {result['validation_passed']}")
        print(f"格式有效: {result['format_valid']}")
        if result['errors']:
            print(f"错误: {result['errors']}")
        if result['warnings']:
            print(f"警告: {result['warnings']}")
    else:
        print(f"测试文件不存在: {test_file}")
    
    print("\n基因组质量评估测试完成") 