"""
基因组质量检查模块

提供基因组数据质量评估功能，包括：
- 组装统计（N50/N90计算）
- GC含量分布分析
- N含量分布评估
- 序列复杂度计算
- 综合质量评分和等级评定
"""

import numpy as np
import pandas as pd
from pathlib import Path
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Union
from collections import Counter
import statistics
import math

from utils.file_io import FastaReader, SequenceStats
from utils.logging_setup import NBSLogger, get_logger


@dataclass
class AssemblyStats:
    """组装统计信息"""
    total_sequences: int
    total_length: int
    mean_length: float
    median_length: float
    n50: int
    n90: int
    longest_sequence: int
    shortest_sequence: int
    sequence_lengths: List[int] = field(default_factory=list)


@dataclass
class GCStats:
    """GC含量统计信息"""
    overall_gc_content: float
    gc_content_std: float
    gc_content_distribution: Dict[str, int] = field(default_factory=dict)
    sequence_gc_contents: List[float] = field(default_factory=list)


@dataclass
class NContentStats:
    """N含量统计信息"""
    total_n_count: int
    n_percentage: float
    sequences_with_n: int
    max_n_stretch: int
    n_stretches_distribution: Dict[int, int] = field(default_factory=dict)


@dataclass
class ComplexityStats:
    """序列复杂度统计信息"""
    mean_complexity: float
    complexity_std: float
    low_complexity_regions: int
    sequence_complexities: List[float] = field(default_factory=list)


@dataclass
class QualityAssessment:
    """基因组质量评估结果"""
    file_path: str
    assembly_stats: AssemblyStats
    gc_stats: GCStats
    n_content_stats: NContentStats
    complexity_stats: ComplexityStats
    overall_score: float
    quality_grade: str
    recommendations: List[str] = field(default_factory=list)


class GenomeQualityAssessment:
    """基因组质量评估器"""
    
    def __init__(self, logger: Optional[NBSLogger] = None):
        """
        初始化基因组质量评估器
        
        Args:
            logger: 日志记录器
        """
        if logger is None:
            nbs_logger = NBSLogger()
            self.logger = nbs_logger.setup_logging()
        else:
            self.logger = logger
        
        # 质量评估阈值
        self.quality_thresholds = {
            'excellent': {'min_score': 90, 'max_n_percent': 1.0, 'min_n50': 1000000},
            'good': {'min_score': 80, 'max_n_percent': 2.0, 'min_n50': 500000},
            'fair': {'min_score': 70, 'max_n_percent': 5.0, 'min_n50': 100000},
            'poor': {'min_score': 0, 'max_n_percent': 100.0, 'min_n50': 0}
        }
    
    def assess_genome_quality(self, fasta_file: Union[str, Path]) -> QualityAssessment:
        """
        评估基因组质量
        
        Args:
            fasta_file: FASTA文件路径
            
        Returns:
            QualityAssessment: 质量评估结果
        """
        fasta_file = Path(fasta_file)
        self.logger.info(f"开始评估基因组质量: {fasta_file}")
        
        try:
            # 读取序列数据
            sequences = self._read_sequences(fasta_file)
            
            # 计算各项统计指标
            assembly_stats = self._calculate_assembly_stats(sequences)
            gc_stats = self._calculate_gc_stats(sequences)
            n_content_stats = self._calculate_n_content_stats(sequences)
            complexity_stats = self._calculate_complexity_stats(sequences)
            
            # 计算综合评分和等级
            overall_score = self._calculate_overall_score(
                assembly_stats, gc_stats, n_content_stats, complexity_stats
            )
            quality_grade = self._determine_quality_grade(overall_score, assembly_stats, n_content_stats)
            
            # 生成建议
            recommendations = self._generate_recommendations(
                assembly_stats, gc_stats, n_content_stats, complexity_stats
            )
            
            assessment = QualityAssessment(
                file_path=str(fasta_file),
                assembly_stats=assembly_stats,
                gc_stats=gc_stats,
                n_content_stats=n_content_stats,
                complexity_stats=complexity_stats,
                overall_score=overall_score,
                quality_grade=quality_grade,
                recommendations=recommendations
            )
            
            self.logger.info(f"基因组质量评估完成，总分: {overall_score:.1f}, 等级: {quality_grade}")
            return assessment
            
        except Exception as e:
            self.logger.error(f"基因组质量评估失败: {e}")
            raise
    
    def _read_sequences(self, fasta_file: Path) -> List[Tuple[str, str]]:
        """读取序列数据"""
        sequences = []
        reader = FastaReader(self.logger)
        
        for seq_record in reader.read_sequences(fasta_file):
            sequences.append((seq_record.id, str(seq_record.seq)))
        
        self.logger.debug(f"读取了 {len(sequences)} 个序列")
        return sequences
    
    def _calculate_assembly_stats(self, sequences: List[Tuple[str, str]]) -> AssemblyStats:
        """计算组装统计信息"""
        lengths = [len(seq) for _, seq in sequences]
        lengths.sort(reverse=True)
        
        total_length = sum(lengths)
        total_sequences = len(sequences)
        
        # 计算N50和N90
        n50 = self._calculate_nx(lengths, 50)
        n90 = self._calculate_nx(lengths, 90)
        
        return AssemblyStats(
            total_sequences=total_sequences,
            total_length=total_length,
            mean_length=statistics.mean(lengths),
            median_length=statistics.median(lengths),
            n50=n50,
            n90=n90,
            longest_sequence=max(lengths),
            shortest_sequence=min(lengths),
            sequence_lengths=lengths
        )
    
    def _calculate_nx(self, sorted_lengths: List[int], x: int) -> int:
        """
        计算Nx值（如N50, N90）
        
        Args:
            sorted_lengths: 按长度降序排列的序列长度列表
            x: 百分比值（如50表示N50）
            
        Returns:
            Nx值
        """
        total_length = sum(sorted_lengths)
        target_length = total_length * x / 100
        
        cumulative_length = 0
        for length in sorted_lengths:
            cumulative_length += length
            if cumulative_length >= target_length:
                return length
        
        return sorted_lengths[-1] if sorted_lengths else 0
    
    def _calculate_gc_stats(self, sequences: List[Tuple[str, str]]) -> GCStats:
        """计算GC含量统计"""
        gc_contents = []
        total_gc = 0
        total_bases = 0
        
        for seq_id, sequence in sequences:
            sequence_upper = sequence.upper()
            gc_count = sequence_upper.count('G') + sequence_upper.count('C')
            valid_bases = len([b for b in sequence_upper if b in 'ATCG'])
            
            if valid_bases > 0:
                gc_content = (gc_count / valid_bases) * 100
                gc_contents.append(gc_content)
                total_gc += gc_count
                total_bases += valid_bases
        
        overall_gc = (total_gc / total_bases * 100) if total_bases > 0 else 0
        gc_std = statistics.stdev(gc_contents) if len(gc_contents) > 1 else 0
        
        # GC含量分布
        gc_distribution = self._create_gc_distribution(gc_contents)
        
        return GCStats(
            overall_gc_content=overall_gc,
            gc_content_std=gc_std,
            gc_content_distribution=gc_distribution,
            sequence_gc_contents=gc_contents
        )
    
    def _create_gc_distribution(self, gc_contents: List[float]) -> Dict[str, int]:
        """创建GC含量分布"""
        distribution = {
            'very_low_gc_0_30': 0,
            'low_gc_30_40': 0,
            'normal_gc_40_60': 0,
            'high_gc_60_70': 0,
            'very_high_gc_70_100': 0
        }
        
        for gc in gc_contents:
            if gc < 30:
                distribution['very_low_gc_0_30'] += 1
            elif gc < 40:
                distribution['low_gc_30_40'] += 1
            elif gc < 60:
                distribution['normal_gc_40_60'] += 1
            elif gc < 70:
                distribution['high_gc_60_70'] += 1
            else:
                distribution['very_high_gc_70_100'] += 1
        
        return distribution
    
    def _calculate_n_content_stats(self, sequences: List[Tuple[str, str]]) -> NContentStats:
        """计算N含量统计"""
        total_n_count = 0
        total_length = 0
        sequences_with_n = 0
        max_n_stretch = 0
        n_stretches = []
        
        for seq_id, sequence in sequences:
            sequence_upper = sequence.upper()
            n_count = sequence_upper.count('N')
            
            if n_count > 0:
                sequences_with_n += 1
                total_n_count += n_count
                
                # 查找N的连续片段
                current_stretch = 0
                for base in sequence_upper:
                    if base == 'N':
                        current_stretch += 1
                    else:
                        if current_stretch > 0:
                            n_stretches.append(current_stretch)
                            max_n_stretch = max(max_n_stretch, current_stretch)
                        current_stretch = 0
                
                # 处理序列末尾的N片段
                if current_stretch > 0:
                    n_stretches.append(current_stretch)
                    max_n_stretch = max(max_n_stretch, current_stretch)
            
            total_length += len(sequence)
        
        n_percentage = (total_n_count / total_length * 100) if total_length > 0 else 0
        
        # N片段分布
        n_stretches_distribution = Counter(n_stretches)
        
        return NContentStats(
            total_n_count=total_n_count,
            n_percentage=n_percentage,
            sequences_with_n=sequences_with_n,
            max_n_stretch=max_n_stretch,
            n_stretches_distribution=dict(n_stretches_distribution)
        )
    
    def _calculate_complexity_stats(self, sequences: List[Tuple[str, str]]) -> ComplexityStats:
        """计算序列复杂度统计"""
        complexities = []
        low_complexity_regions = 0
        
        for seq_id, sequence in sequences:
            complexity = self._calculate_sequence_complexity(sequence)
            complexities.append(complexity)
            
            # 检查低复杂度区域（复杂度 < 0.5）
            if complexity < 0.5:
                low_complexity_regions += 1
        
        mean_complexity = statistics.mean(complexities) if complexities else 0
        complexity_std = statistics.stdev(complexities) if len(complexities) > 1 else 0
        
        return ComplexityStats(
            mean_complexity=mean_complexity,
            complexity_std=complexity_std,
            low_complexity_regions=low_complexity_regions,
            sequence_complexities=complexities
        )
    
    def _calculate_sequence_complexity(self, sequence: str, window_size: int = 100) -> float:
        """
        计算序列复杂度（基于香农熵）
        
        Args:
            sequence: DNA序列
            window_size: 窗口大小
            
        Returns:
            复杂度值（0-1之间，1表示最高复杂度）
        """
        if len(sequence) < window_size:
            window_size = len(sequence)
        
        if window_size == 0:
            return 0.0
        
        # 计算k-mer频率（k=4）
        k = min(4, window_size)
        kmers = []
        
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i+k].upper()
            if all(base in 'ATCGN' for base in kmer):
                kmers.append(kmer)
        
        if not kmers:
            return 0.0
        
        # 计算香农熵
        kmer_counts = Counter(kmers)
        total_kmers = len(kmers)
        
        entropy = 0
        for count in kmer_counts.values():
            probability = count / total_kmers
            entropy -= probability * math.log2(probability)
        
        # 标准化熵值（最大熵为log2(4^k)）
        max_entropy = k * 2  # log2(4^k) = k * log2(4) = k * 2
        normalized_entropy = entropy / max_entropy if max_entropy > 0 else 0
        
        return min(1.0, normalized_entropy)
    
    def _calculate_overall_score(self, assembly_stats: AssemblyStats, gc_stats: GCStats,
                               n_content_stats: NContentStats, complexity_stats: ComplexityStats) -> float:
        """计算综合质量评分"""
        scores = []
        
        # 组装质量评分（40%权重）
        assembly_score = self._score_assembly_quality(assembly_stats)
        scores.append(('assembly', assembly_score, 0.4))
        
        # GC含量评分（20%权重）
        gc_score = self._score_gc_quality(gc_stats)
        scores.append(('gc', gc_score, 0.2))
        
        # N含量评分（25%权重）
        n_score = self._score_n_content(n_content_stats)
        scores.append(('n_content', n_score, 0.25))
        
        # 复杂度评分（15%权重）
        complexity_score = self._score_complexity(complexity_stats)
        scores.append(('complexity', complexity_score, 0.15))
        
        # 计算加权平均分
        weighted_score = sum(score * weight for _, score, weight in scores)
        
        self.logger.debug(f"各项评分: {[(name, f'{score:.1f}') for name, score, _ in scores]}")
        
        return weighted_score
    
    def _score_assembly_quality(self, stats: AssemblyStats) -> float:
        """评估组装质量"""
        score = 0
        
        # N50评分（50%）
        if stats.n50 >= 10000000:  # 10Mb
            score += 50
        elif stats.n50 >= 1000000:  # 1Mb
            score += 40
        elif stats.n50 >= 100000:  # 100kb
            score += 30
        elif stats.n50 >= 10000:  # 10kb
            score += 20
        else:
            score += 10
        
        # 序列数量评分（30%）
        if stats.total_sequences <= 1000:
            score += 30
        elif stats.total_sequences <= 10000:
            score += 20
        elif stats.total_sequences <= 100000:
            score += 10
        else:
            score += 5
        
        # 平均长度评分（20%）
        if stats.mean_length >= 100000:
            score += 20
        elif stats.mean_length >= 10000:
            score += 15
        elif stats.mean_length >= 1000:
            score += 10
        else:
            score += 5
        
        return min(100, score)
    
    def _score_gc_quality(self, stats: GCStats) -> float:
        """评估GC含量质量"""
        score = 100
        
        # GC含量合理性（植物基因组通常35-45%）
        if 35 <= stats.overall_gc_content <= 45:
            score -= 0
        elif 30 <= stats.overall_gc_content <= 50:
            score -= 10
        elif 25 <= stats.overall_gc_content <= 55:
            score -= 20
        else:
            score -= 30
        
        # GC含量变异性
        if stats.gc_content_std <= 5:
            score -= 0
        elif stats.gc_content_std <= 10:
            score -= 10
        elif stats.gc_content_std <= 15:
            score -= 20
        else:
            score -= 30
        
        return max(0, score)
    
    def _score_n_content(self, stats: NContentStats) -> float:
        """评估N含量质量"""
        score = 100
        
        # N含量百分比
        if stats.n_percentage <= 1:
            score -= 0
        elif stats.n_percentage <= 2:
            score -= 10
        elif stats.n_percentage <= 5:
            score -= 30
        elif stats.n_percentage <= 10:
            score -= 50
        else:
            score -= 70
        
        # 最大N片段长度
        if stats.max_n_stretch <= 100:
            score -= 0
        elif stats.max_n_stretch <= 1000:
            score -= 10
        elif stats.max_n_stretch <= 10000:
            score -= 20
        else:
            score -= 30
        
        return max(0, score)
    
    def _score_complexity(self, stats: ComplexityStats) -> float:
        """评估序列复杂度质量"""
        score = 100
        
        # 平均复杂度
        if stats.mean_complexity >= 0.8:
            score -= 0
        elif stats.mean_complexity >= 0.6:
            score -= 10
        elif stats.mean_complexity >= 0.4:
            score -= 20
        else:
            score -= 40
        
        # 低复杂度区域比例
        if hasattr(stats, 'sequence_complexities') and stats.sequence_complexities:
            low_complexity_ratio = stats.low_complexity_regions / len(stats.sequence_complexities)
            if low_complexity_ratio <= 0.1:
                score -= 0
            elif low_complexity_ratio <= 0.2:
                score -= 10
            elif low_complexity_ratio <= 0.3:
                score -= 20
            else:
                score -= 30
        
        return max(0, score)
    
    def _determine_quality_grade(self, score: float, assembly_stats: AssemblyStats,
                               n_content_stats: NContentStats) -> str:
        """确定质量等级"""
        # 基于综合评分和关键指标确定等级
        for grade, thresholds in self.quality_thresholds.items():
            if (score >= thresholds['min_score'] and
                n_content_stats.n_percentage <= thresholds['max_n_percent'] and
                assembly_stats.n50 >= thresholds['min_n50']):
                return grade
        
        return 'poor'
    
    def _generate_recommendations(self, assembly_stats: AssemblyStats, gc_stats: GCStats,
                                n_content_stats: NContentStats, complexity_stats: ComplexityStats) -> List[str]:
        """生成改进建议"""
        recommendations = []
        
        # 组装质量建议
        if assembly_stats.n50 < 100000:
            recommendations.append("建议使用长读长测序技术或scaffolding方法提高组装连续性")
        
        if assembly_stats.total_sequences > 10000:
            recommendations.append("序列片段过多，建议进行gap filling或序列合并")
        
        # GC含量建议
        if gc_stats.overall_gc_content < 30 or gc_stats.overall_gc_content > 50:
            recommendations.append("GC含量异常，建议检查数据质量或物种特异性")
        
        if gc_stats.gc_content_std > 10:
            recommendations.append("GC含量变异较大，可能存在污染序列，建议进行序列过滤")
        
        # N含量建议
        if n_content_stats.n_percentage > 5:
            recommendations.append("N含量过高，建议进行gap filling或使用更高质量的原始数据")
        
        if n_content_stats.max_n_stretch > 10000:
            recommendations.append("存在大片段N区域，建议使用长读长数据填补gaps")
        
        # 复杂度建议
        if complexity_stats.mean_complexity < 0.6:
            recommendations.append("序列复杂度较低，可能存在重复序列或简单序列，建议进行重复序列分析")
        
        if complexity_stats.low_complexity_regions > len(complexity_stats.sequence_complexities) * 0.2:
            recommendations.append("低复杂度区域较多，建议使用RepeatMasker等工具进行重复序列注释")
        
        if not recommendations:
            recommendations.append("基因组质量良好，可以进行后续分析")
        
        return recommendations
    
    def export_quality_report(self, assessment: QualityAssessment, output_file: Union[str, Path]) -> None:
        """
        导出质量评估报告
        
        Args:
            assessment: 质量评估结果
            output_file: 输出文件路径
        """
        output_file = Path(output_file)
        
        try:
            with open(output_file, 'w', encoding='utf-8') as f:
                f.write("# 基因组质量评估报告\n\n")
                f.write(f"**文件路径**: {assessment.file_path}\n")
                f.write(f"**综合评分**: {assessment.overall_score:.1f}/100\n")
                f.write(f"**质量等级**: {assessment.quality_grade.upper()}\n\n")
                
                # 组装统计
                f.write("## 组装统计\n\n")
                stats = assessment.assembly_stats
                f.write(f"- 序列总数: {stats.total_sequences:,}\n")
                f.write(f"- 总长度: {stats.total_length:,} bp\n")
                f.write(f"- 平均长度: {stats.mean_length:,.0f} bp\n")
                f.write(f"- 中位数长度: {stats.median_length:,.0f} bp\n")
                f.write(f"- N50: {stats.n50:,} bp\n")
                f.write(f"- N90: {stats.n90:,} bp\n")
                f.write(f"- 最长序列: {stats.longest_sequence:,} bp\n")
                f.write(f"- 最短序列: {stats.shortest_sequence:,} bp\n\n")
                
                # GC含量统计
                f.write("## GC含量统计\n\n")
                gc_stats = assessment.gc_stats
                f.write(f"- 总体GC含量: {gc_stats.overall_gc_content:.2f}%\n")
                f.write(f"- GC含量标准差: {gc_stats.gc_content_std:.2f}\n")
                f.write("- GC含量分布:\n")
                for category, count in gc_stats.gc_content_distribution.items():
                    f.write(f"  - {category.replace('_', ' ').title()}: {count}\n")
                f.write("\n")
                
                # N含量统计
                f.write("## N含量统计\n\n")
                n_stats = assessment.n_content_stats
                f.write(f"- N总数: {n_stats.total_n_count:,}\n")
                f.write(f"- N含量百分比: {n_stats.n_percentage:.2f}%\n")
                f.write(f"- 含N序列数: {n_stats.sequences_with_n}\n")
                f.write(f"- 最大N片段长度: {n_stats.max_n_stretch:,}\n\n")
                
                # 复杂度统计
                f.write("## 序列复杂度统计\n\n")
                comp_stats = assessment.complexity_stats
                f.write(f"- 平均复杂度: {comp_stats.mean_complexity:.3f}\n")
                f.write(f"- 复杂度标准差: {comp_stats.complexity_std:.3f}\n")
                f.write(f"- 低复杂度区域数: {comp_stats.low_complexity_regions}\n\n")
                
                # 建议
                f.write("## 改进建议\n\n")
                for i, recommendation in enumerate(assessment.recommendations, 1):
                    f.write(f"{i}. {recommendation}\n")
            
            self.logger.info(f"质量评估报告已导出至: {output_file}")
            
        except Exception as e:
            self.logger.error(f"导出质量评估报告失败: {e}")
            raise 