"""
Augustus质量评估器
Augustus Quality Assessor

本模块提供对Augustus基因预测结果的全面质量评估：
- 核苷酸级别评估 (Nucleotide-level assessment)
- 外显子级别评估 (Exon-level assessment)  
- 基因级别评估 (Gene-level assessment)
- 预测可靠性分析 (Prediction reliability analysis)
- 质量控制指标 (Quality control metrics)
"""

import os
import logging
import numpy as np
from pathlib import Path
from typing import List, Dict, Optional, Tuple, Union, Set
from dataclasses import dataclass, field
from collections import defaultdict, Counter
import json
import warnings

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis

from .augustus_output_parser import GeneStructure, AugustusParsingResult


@dataclass
class QualityMetrics:
    """质量评估指标"""
    # 核苷酸级别指标
    nucleotide_sensitivity: float = 0.0  # 敏感性 (Sn)
    nucleotide_specificity: float = 0.0  # 特异性 (Sp)  
    nucleotide_f1_score: float = 0.0     # F1分数
    nucleotide_accuracy: float = 0.0     # 准确率
    
    # 外显子级别指标
    exon_sensitivity: float = 0.0        # 外显子敏感性
    exon_specificity: float = 0.0        # 外显子特异性
    exon_f1_score: float = 0.0           # 外显子F1分数
    
    # 基因级别指标
    gene_sensitivity: float = 0.0        # 基因敏感性
    gene_specificity: float = 0.0        # 基因特异性 
    gene_f1_score: float = 0.0           # 基因F1分数
    
    # 预测质量指标
    avg_gene_score: float = 0.0          # 平均基因预测分数
    coding_potential_score: float = 0.0  # 编码潜力分数
    structural_consistency: float = 0.0   # 结构一致性分数
    
    # 统计计数
    true_positives: int = 0
    false_positives: int = 0
    false_negatives: int = 0
    true_negatives: int = 0
    
    # 详细分析
    prediction_confidence: Dict[str, float] = field(default_factory=dict)
    quality_warnings: List[str] = field(default_factory=list)


@dataclass
class ReferenceAnnotation:
    """参考注释数据"""
    genes: List[GeneStructure] = field(default_factory=list)
    coding_regions: Set[Tuple[str, int, int]] = field(default_factory=set)  # (seqname, start, end)
    exon_regions: Set[Tuple[str, int, int, str]] = field(default_factory=set)  # (seqname, start, end, strand)


class AugustusQualityAssessor:
    """Augustus质量评估器"""
    
    def __init__(self, min_overlap: float = 0.5, 
                 min_gene_length: int = 150,
                 max_intergenic_distance: int = 1000):
        """
        初始化质量评估器
        
        Args:
            min_overlap: 最小重叠比例用于匹配预测和参考基因
            min_gene_length: 最小基因长度阈值
            max_intergenic_distance: 最大基因间距离
        """
        self.min_overlap = min_overlap
        self.min_gene_length = min_gene_length
        self.max_intergenic_distance = max_intergenic_distance
        self.logger = logging.getLogger(__name__)
    
    def assess_prediction_quality(self, 
                                 predictions: AugustusParsingResult,
                                 reference: Optional[ReferenceAnnotation] = None,
                                 genome_sequence: Optional[str] = None) -> QualityMetrics:
        """
        评估预测质量
        
        Args:
            predictions: Augustus预测结果
            reference: 参考注释（如果有的话）
            genome_sequence: 基因组序列文件路径
        
        Returns:
            QualityMetrics: 质量评估指标
        """
        self.logger.info("开始Augustus预测质量评估")
        
        metrics = QualityMetrics()
        
        # 1. 内在质量评估（不需要参考注释）
        self._assess_intrinsic_quality(predictions, metrics, genome_sequence)
        
        # 2. 如果有参考注释，进行外在质量评估
        if reference:
            self._assess_extrinsic_quality(predictions, reference, metrics)
        
        # 3. 预测可靠性分析
        self._analyze_prediction_reliability(predictions, metrics)
        
        # 4. 质量控制检查
        self._quality_control_checks(predictions, metrics)
        
        self.logger.info(f"质量评估完成 - F1分数: {metrics.nucleotide_f1_score:.3f}")
        return metrics
    
    def _assess_intrinsic_quality(self, predictions: AugustusParsingResult, 
                                 metrics: QualityMetrics, 
                                 genome_sequence: Optional[str] = None):
        """评估内在质量（不依赖参考注释）"""
        genes = predictions.genes
        if not genes:
            return
        
        # 计算平均预测分数
        scores = [g.score for g in genes if g.score > 0]
        if scores:
            metrics.avg_gene_score = np.mean(scores)
        
        # 分析编码潜力
        coding_potentials = []
        structural_scores = []
        
        for gene in genes:
            # 编码潜力评估
            coding_potential = self._calculate_coding_potential(gene, genome_sequence)
            coding_potentials.append(coding_potential)
            gene.coding_potential = coding_potential
            
            # 结构一致性评估
            structural_score = self._calculate_structural_consistency(gene)
            structural_scores.append(structural_score)
        
        if coding_potentials:
            metrics.coding_potential_score = np.mean(coding_potentials)
        
        if structural_scores:
            metrics.structural_consistency = np.mean(structural_scores)
        
        self.logger.info(f"内在质量评估完成 - 平均编码潜力: {metrics.coding_potential_score:.3f}")
    
    def _assess_extrinsic_quality(self, predictions: AugustusParsingResult,
                                 reference: ReferenceAnnotation,
                                 metrics: QualityMetrics):
        """评估外在质量（基于参考注释）"""
        pred_genes = predictions.genes
        ref_genes = reference.genes
        
        if not pred_genes or not ref_genes:
            self.logger.warning("预测或参考基因列表为空，跳过外在质量评估")
            return
        
        # 1. 核苷酸级别评估
        self._assess_nucleotide_level(pred_genes, ref_genes, metrics)
        
        # 2. 外显子级别评估
        self._assess_exon_level(pred_genes, ref_genes, metrics)
        
        # 3. 基因级别评估
        self._assess_gene_level(pred_genes, ref_genes, metrics)
        
        self.logger.info(f"外在质量评估完成 - 基因F1: {metrics.gene_f1_score:.3f}")
    
    def _assess_nucleotide_level(self, pred_genes: List[GeneStructure],
                               ref_genes: List[GeneStructure],
                               metrics: QualityMetrics):
        """核苷酸级别质量评估"""
        # 构建编码区域集合
        pred_coding = set()
        ref_coding = set()
        
        for gene in pred_genes:
            for start, end in gene.cds_regions:
                for pos in range(start, end + 1):
                    pred_coding.add((gene.seqname, pos))
        
        for gene in ref_genes:
            for start, end in gene.cds_regions:
                for pos in range(start, end + 1):
                    ref_coding.add((gene.seqname, pos))
        
        # 计算TP, FP, FN
        tp = len(pred_coding.intersection(ref_coding))
        fp = len(pred_coding - ref_coding)
        fn = len(ref_coding - pred_coding)
        
        metrics.true_positives = tp
        metrics.false_positives = fp
        metrics.false_negatives = fn
        
        # 计算敏感性、特异性和F1分数
        if tp + fn > 0:
            metrics.nucleotide_sensitivity = tp / (tp + fn)
        
        if tp + fp > 0:
            metrics.nucleotide_specificity = tp / (tp + fp)
        
        if metrics.nucleotide_sensitivity + metrics.nucleotide_specificity > 0:
            metrics.nucleotide_f1_score = (2 * metrics.nucleotide_sensitivity * 
                                         metrics.nucleotide_specificity) / (
                                         metrics.nucleotide_sensitivity + 
                                         metrics.nucleotide_specificity)
        
        # 计算准确率
        total_positions = len(pred_coding.union(ref_coding))
        if total_positions > 0:
            metrics.nucleotide_accuracy = tp / total_positions
    
    def _assess_exon_level(self, pred_genes: List[GeneStructure],
                          ref_genes: List[GeneStructure],
                          metrics: QualityMetrics):
        """外显子级别质量评估"""
        pred_exons = set()
        ref_exons = set()
        
        # 收集所有外显子
        for gene in pred_genes:
            for exon in gene.exons:
                pred_exons.add((gene.seqname, exon.start, exon.end, gene.strand))
        
        for gene in ref_genes:
            for exon in gene.exons:
                ref_exons.add((gene.seqname, exon.start, exon.end, gene.strand))
        
        # 计算外显子匹配
        matched_pred = set()
        matched_ref = set()
        
        for pred_exon in pred_exons:
            for ref_exon in ref_exons:
                if self._exons_overlap(pred_exon, ref_exon, self.min_overlap):
                    matched_pred.add(pred_exon)
                    matched_ref.add(ref_exon)
        
        # 计算外显子级别指标
        tp_exons = len(matched_pred)
        fp_exons = len(pred_exons) - tp_exons
        fn_exons = len(ref_exons) - len(matched_ref)
        
        if tp_exons + fn_exons > 0:
            metrics.exon_sensitivity = tp_exons / (tp_exons + fn_exons)
        
        if tp_exons + fp_exons > 0:
            metrics.exon_specificity = tp_exons / (tp_exons + fp_exons)
        
        if metrics.exon_sensitivity + metrics.exon_specificity > 0:
            metrics.exon_f1_score = (2 * metrics.exon_sensitivity * 
                                   metrics.exon_specificity) / (
                                   metrics.exon_sensitivity + 
                                   metrics.exon_specificity)
    
    def _assess_gene_level(self, pred_genes: List[GeneStructure],
                          ref_genes: List[GeneStructure],
                          metrics: QualityMetrics):
        """基因级别质量评估"""
        matched_pred = set()
        matched_ref = set()
        
        # 基因匹配
        for i, pred_gene in enumerate(pred_genes):
            for j, ref_gene in enumerate(ref_genes):
                if self._genes_overlap(pred_gene, ref_gene, self.min_overlap):
                    matched_pred.add(i)
                    matched_ref.add(j)
        
        # 计算基因级别指标
        tp_genes = len(matched_pred)
        fp_genes = len(pred_genes) - tp_genes
        fn_genes = len(ref_genes) - len(matched_ref)
        
        if tp_genes + fn_genes > 0:
            metrics.gene_sensitivity = tp_genes / (tp_genes + fn_genes)
        
        if tp_genes + fp_genes > 0:
            metrics.gene_specificity = tp_genes / (tp_genes + fp_genes)
        
        if metrics.gene_sensitivity + metrics.gene_specificity > 0:
            metrics.gene_f1_score = (2 * metrics.gene_sensitivity * 
                                   metrics.gene_specificity) / (
                                   metrics.gene_sensitivity + 
                                   metrics.gene_specificity)
    
    def _calculate_coding_potential(self, gene: GeneStructure, 
                                  genome_sequence: Optional[str] = None) -> float:
        """计算编码潜力分数"""
        score = 0.0
        factors = 0
        
        # 1. 基因长度评估
        if gene.gene_length >= self.min_gene_length:
            score += 1.0
        else:
            score += gene.gene_length / self.min_gene_length
        factors += 1
        
        # 2. 外显子数量评估
        if gene.exon_count > 0:
            # 多外显子基因通常更可靠
            exon_score = min(1.0, gene.exon_count / 3.0)
            score += exon_score
        factors += 1
        
        # 3. CDS长度评估
        if gene.cds_length > 0:
            # CDS应该是3的倍数
            if gene.cds_length % 3 == 0:
                score += 1.0
            else:
                score += 0.5
        factors += 1
        
        # 4. GC含量评估
        if gene.gc_content > 0:
            # 合理的GC含量范围
            if 0.35 <= gene.gc_content <= 0.65:
                score += 1.0
            else:
                score += 0.5
        factors += 1
        
        # 5. 蛋白质序列评估
        if gene.protein_sequence:
            protein_score = self._evaluate_protein_sequence(gene.protein_sequence)
            score += protein_score
        factors += 1
        
        return score / factors if factors > 0 else 0.0
    
    def _evaluate_protein_sequence(self, protein_seq: str) -> float:
        """评估蛋白质序列质量"""
        if not protein_seq or len(protein_seq) < 50:
            return 0.0
        
        score = 0.0
        
        try:
            # 使用BioPython的蛋白质分析
            analysis = ProteinAnalysis(protein_seq)
            
            # 1. 检查起始密码子（M）
            if protein_seq.startswith('M'):
                score += 0.3
            
            # 2. 检查终止密码子（*）
            if protein_seq.endswith('*'):
                score += 0.2
            elif '*' not in protein_seq:
                score += 0.1  # 没有内部终止密码子
            
            # 3. 氨基酸组成分析
            aa_percent = analysis.get_amino_acids_percent()
            
            # 检查是否有合理的氨基酸分布
            common_aa = ['A', 'L', 'S', 'G', 'V', 'E', 'K', 'I', 'D', 'T']
            common_sum = sum(aa_percent.get(aa, 0) for aa in common_aa)
            
            if common_sum > 0.6:  # 常见氨基酸占60%以上
                score += 0.3
            
            # 4. 疏水性评估
            try:
                hydrophobicity = analysis.gravy()
                if -2.0 <= hydrophobicity <= 2.0:  # 合理范围
                    score += 0.2
            except:
                pass
            
        except Exception as e:
            self.logger.warning(f"蛋白质序列分析失败: {e}")
            return 0.5
        
        return min(1.0, score)
    
    def _calculate_structural_consistency(self, gene: GeneStructure) -> float:
        """计算结构一致性分数"""
        score = 0.0
        factors = 0
        
        # 1. 外显子-内含子结构合理性
        if gene.exon_count > 1:
            expected_introns = gene.exon_count - 1
            if gene.intron_count == expected_introns:
                score += 1.0
            else:
                score += 0.5
        else:
            score += 1.0  # 单外显子基因
        factors += 1
        
        # 2. 外显子长度分布
        if gene.exons:
            exon_lengths = [e.length for e in gene.exons]
            avg_length = np.mean(exon_lengths)
            
            # 合理的外显子长度范围
            if 50 <= avg_length <= 500:
                score += 1.0
            else:
                score += 0.5
        factors += 1
        
        # 3. 内含子长度检查
        if gene.introns:
            intron_lengths = [i.length for i in gene.introns]
            avg_intron_length = np.mean(intron_lengths)
            
            # 合理的内含子长度
            if 50 <= avg_intron_length <= 10000:
                score += 1.0
            else:
                score += 0.5
        factors += 1
        
        # 4. 剪接位点检查（如果有序列信息）
        splice_score = self._check_splice_sites(gene)
        score += splice_score
        factors += 1
        
        return score / factors if factors > 0 else 0.0
    
    def _check_splice_sites(self, gene: GeneStructure) -> float:
        """检查剪接位点"""
        if not gene.introns:
            return 1.0  # 单外显子基因
        
        canonical_donors = ['GT', 'GC']
        canonical_acceptors = ['AG']
        
        correct_sites = 0
        total_sites = 0
        
        for intron in gene.introns:
            total_sites += 2  # donor + acceptor
            
            if intron.donor_site in canonical_donors:
                correct_sites += 1
            
            if intron.acceptor_site in canonical_acceptors:
                correct_sites += 1
        
        return correct_sites / total_sites if total_sites > 0 else 0.0
    
    def _analyze_prediction_reliability(self, predictions: AugustusParsingResult,
                                      metrics: QualityMetrics):
        """分析预测可靠性"""
        genes = predictions.genes
        if not genes:
            return
        
        confidence_scores = {}
        
        for gene in genes:
            # 综合置信度评估
            confidence = 0.0
            factors = 0
            
            # 1. Augustus预测分数
            if gene.score > 0:
                confidence += min(1.0, gene.score / 1.0)  # 归一化到0-1
                factors += 1
            
            # 2. 编码潜力
            if gene.coding_potential > 0:
                confidence += gene.coding_potential
                factors += 1
            
            # 3. 基因长度合理性
            length_score = 1.0 if gene.gene_length >= self.min_gene_length else 0.5
            confidence += length_score
            factors += 1
            
            # 4. 结构完整性
            if gene.protein_sequence:
                structure_score = 1.0 if (gene.protein_sequence.startswith('M') and 
                                        gene.protein_sequence.endswith('*')) else 0.5
                confidence += structure_score
                factors += 1
            
            final_confidence = confidence / factors if factors > 0 else 0.0
            confidence_scores[gene.gene_id] = final_confidence
        
        metrics.prediction_confidence = confidence_scores
        
        # 计算平均置信度
        if confidence_scores:
            avg_confidence = np.mean(list(confidence_scores.values()))
            self.logger.info(f"平均预测置信度: {avg_confidence:.3f}")
    
    def _quality_control_checks(self, predictions: AugustusParsingResult,
                              metrics: QualityMetrics):
        """质量控制检查"""
        warnings = []
        genes = predictions.genes
        
        if not genes:
            warnings.append("没有预测到任何基因")
            metrics.quality_warnings = warnings
            return
        
        # 1. 检查过短基因
        short_genes = [g for g in genes if g.gene_length < self.min_gene_length]
        if short_genes:
            warnings.append(f"发现 {len(short_genes)} 个过短基因 (< {self.min_gene_length} bp)")
        
        # 2. 检查无外显子基因
        no_exon_genes = [g for g in genes if g.exon_count == 0]
        if no_exon_genes:
            warnings.append(f"发现 {len(no_exon_genes)} 个无外显子基因")
        
        # 3. 检查异常长基因
        long_genes = [g for g in genes if g.gene_length > 50000]
        if long_genes:
            warnings.append(f"发现 {len(long_genes)} 个异常长基因 (> 50kb)")
        
        # 4. 检查重叠基因
        overlaps = self._find_overlapping_genes(genes)
        if overlaps:
            warnings.append(f"发现 {len(overlaps)} 对重叠基因")
        
        # 5. 检查基因密度
        total_length = sum(g.gene_length for g in genes)
        if total_length > 0:
            gene_density = len(genes) / (total_length / 1000)  # genes per kb
            if gene_density < 0.1:
                warnings.append(f"基因密度过低: {gene_density:.3f} genes/kb")
            elif gene_density > 2.0:
                warnings.append(f"基因密度过高: {gene_density:.3f} genes/kb")
        
        metrics.quality_warnings = warnings
        
        if warnings:
            self.logger.warning(f"质量控制发现 {len(warnings)} 个问题")
            for warning in warnings:
                self.logger.warning(f"  - {warning}")
    
    def _find_overlapping_genes(self, genes: List[GeneStructure]) -> List[Tuple[str, str]]:
        """查找重叠基因"""
        overlaps = []
        
        for i, gene1 in enumerate(genes):
            for j, gene2 in enumerate(genes[i+1:], i+1):
                if (gene1.seqname == gene2.seqname and 
                    gene1.strand == gene2.strand and
                    self._ranges_overlap(gene1.start, gene1.end, 
                                       gene2.start, gene2.end)):
                    overlaps.append((gene1.gene_id, gene2.gene_id))
        
        return overlaps
    
    def _ranges_overlap(self, start1: int, end1: int, 
                       start2: int, end2: int, min_overlap: float = 0.1) -> bool:
        """检查两个区间是否重叠"""
        overlap_start = max(start1, start2)
        overlap_end = min(end1, end2)
        
        if overlap_start <= overlap_end:
            overlap_length = overlap_end - overlap_start + 1
            min_length = min(end1 - start1 + 1, end2 - start2 + 1)
            return overlap_length >= min_length * min_overlap
        
        return False
    
    def _exons_overlap(self, exon1: Tuple, exon2: Tuple, min_overlap: float) -> bool:
        """检查两个外显子是否重叠"""
        seqname1, start1, end1, strand1 = exon1
        seqname2, start2, end2, strand2 = exon2
        
        return (seqname1 == seqname2 and strand1 == strand2 and
                self._ranges_overlap(start1, end1, start2, end2, min_overlap))
    
    def _genes_overlap(self, gene1: GeneStructure, gene2: GeneStructure, 
                      min_overlap: float) -> bool:
        """检查两个基因是否重叠"""
        return (gene1.seqname == gene2.seqname and 
                gene1.strand == gene2.strand and
                self._ranges_overlap(gene1.start, gene1.end, 
                                   gene2.start, gene2.end, min_overlap))
    
    def export_quality_report(self, metrics: QualityMetrics, 
                            output_file: str, format: str = 'json'):
        """
        导出质量评估报告
        
        Args:
            metrics: 质量指标
            output_file: 输出文件路径
            format: 输出格式 ('json', 'tsv', 'txt')
        """
        if format.lower() == 'json':
            self._export_json_report(metrics, output_file)
        elif format.lower() == 'tsv':
            self._export_tsv_report(metrics, output_file)
        elif format.lower() == 'txt':
            self._export_text_report(metrics, output_file)
        else:
            raise ValueError(f"Unsupported format: {format}")
    
    def _export_json_report(self, metrics: QualityMetrics, output_file: str):
        """导出JSON格式报告"""
        report = {
            'quality_assessment': {
                'nucleotide_level': {
                    'sensitivity': metrics.nucleotide_sensitivity,
                    'specificity': metrics.nucleotide_specificity,
                    'f1_score': metrics.nucleotide_f1_score,
                    'accuracy': metrics.nucleotide_accuracy
                },
                'exon_level': {
                    'sensitivity': metrics.exon_sensitivity,
                    'specificity': metrics.exon_specificity,
                    'f1_score': metrics.exon_f1_score
                },
                'gene_level': {
                    'sensitivity': metrics.gene_sensitivity,
                    'specificity': metrics.gene_specificity,
                    'f1_score': metrics.gene_f1_score
                },
                'prediction_quality': {
                    'avg_gene_score': metrics.avg_gene_score,
                    'coding_potential_score': metrics.coding_potential_score,
                    'structural_consistency': metrics.structural_consistency
                },
                'statistics': {
                    'true_positives': metrics.true_positives,
                    'false_positives': metrics.false_positives,
                    'false_negatives': metrics.false_negatives,
                    'true_negatives': metrics.true_negatives
                },
                'prediction_confidence': metrics.prediction_confidence,
                'quality_warnings': metrics.quality_warnings
            }
        }
        
        with open(output_file, 'w', encoding='utf-8') as f:
            json.dump(report, f, indent=2, ensure_ascii=False)
        
        self.logger.info(f"质量评估报告已导出到: {output_file}")
    
    def _export_tsv_report(self, metrics: QualityMetrics, output_file: str):
        """导出TSV格式报告"""
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write("Metric\tValue\n")
            f.write(f"Nucleotide_Sensitivity\t{metrics.nucleotide_sensitivity:.4f}\n")
            f.write(f"Nucleotide_Specificity\t{metrics.nucleotide_specificity:.4f}\n")
            f.write(f"Nucleotide_F1_Score\t{metrics.nucleotide_f1_score:.4f}\n")
            f.write(f"Exon_Sensitivity\t{metrics.exon_sensitivity:.4f}\n")
            f.write(f"Exon_Specificity\t{metrics.exon_specificity:.4f}\n")
            f.write(f"Exon_F1_Score\t{metrics.exon_f1_score:.4f}\n")
            f.write(f"Gene_Sensitivity\t{metrics.gene_sensitivity:.4f}\n")
            f.write(f"Gene_Specificity\t{metrics.gene_specificity:.4f}\n")
            f.write(f"Gene_F1_Score\t{metrics.gene_f1_score:.4f}\n")
            f.write(f"Avg_Gene_Score\t{metrics.avg_gene_score:.4f}\n")
            f.write(f"Coding_Potential_Score\t{metrics.coding_potential_score:.4f}\n")
            f.write(f"Structural_Consistency\t{metrics.structural_consistency:.4f}\n")
        
        self.logger.info(f"质量评估TSV报告已导出到: {output_file}")
    
    def _export_text_report(self, metrics: QualityMetrics, output_file: str):
        """导出文本格式报告"""
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write("Augustus基因预测质量评估报告\n")
            f.write("=" * 50 + "\n\n")
            
            f.write("核苷酸级别评估:\n")
            f.write(f"  敏感性 (Sensitivity): {metrics.nucleotide_sensitivity:.4f}\n")
            f.write(f"  特异性 (Specificity): {metrics.nucleotide_specificity:.4f}\n")
            f.write(f"  F1分数: {metrics.nucleotide_f1_score:.4f}\n")
            f.write(f"  准确率: {metrics.nucleotide_accuracy:.4f}\n\n")
            
            f.write("外显子级别评估:\n")
            f.write(f"  敏感性: {metrics.exon_sensitivity:.4f}\n")
            f.write(f"  特异性: {metrics.exon_specificity:.4f}\n")
            f.write(f"  F1分数: {metrics.exon_f1_score:.4f}\n\n")
            
            f.write("基因级别评估:\n")
            f.write(f"  敏感性: {metrics.gene_sensitivity:.4f}\n")
            f.write(f"  特异性: {metrics.gene_specificity:.4f}\n")
            f.write(f"  F1分数: {metrics.gene_f1_score:.4f}\n\n")
            
            f.write("预测质量指标:\n")
            f.write(f"  平均基因分数: {metrics.avg_gene_score:.4f}\n")
            f.write(f"  编码潜力分数: {metrics.coding_potential_score:.4f}\n")
            f.write(f"  结构一致性: {metrics.structural_consistency:.4f}\n\n")
            
            if metrics.quality_warnings:
                f.write("质量警告:\n")
                for warning in metrics.quality_warnings:
                    f.write(f"  - {warning}\n")
        
        self.logger.info(f"质量评估文本报告已导出到: {output_file}")


def load_reference_annotation(gff3_file: str) -> ReferenceAnnotation:
    """
    从GFF3文件加载参考注释
    
    Args:
        gff3_file: 参考注释GFF3文件路径
    
    Returns:
        ReferenceAnnotation: 参考注释对象
    """
    from .augustus_output_parser import AugustusOutputParser
    
    parser = AugustusOutputParser()
    result = parser.parse_augustus_output(gff3_file)
    
    # 构建编码区域和外显子区域集合
    coding_regions = set()
    exon_regions = set()
    
    for gene in result.genes:
        for start, end in gene.cds_regions:
            coding_regions.add((gene.seqname, start, end))
        
        for exon in gene.exons:
            exon_regions.add((gene.seqname, exon.start, exon.end, gene.strand))
    
    return ReferenceAnnotation(
        genes=result.genes,
        coding_regions=coding_regions,
        exon_regions=exon_regions
    )


def assess_augustus_predictions(predictions_file: str,
                              reference_file: str = None,
                              genome_sequence: str = None,
                              output_dir: str = None) -> QualityMetrics:
    """
    便利函数：评估Augustus预测质量
    
    Args:
        predictions_file: Augustus预测GFF3文件
        reference_file: 参考注释GFF3文件（可选）
        genome_sequence: 基因组序列文件（可选）
        output_dir: 输出目录（可选）
    
    Returns:
        QualityMetrics: 质量评估结果
    """
    from .augustus_output_parser import AugustusOutputParser
    
    # 解析预测结果
    parser = AugustusOutputParser()
    predictions = parser.parse_augustus_output(predictions_file, genome_sequence)
    
    # 加载参考注释（如果提供）
    reference = None
    if reference_file and os.path.exists(reference_file):
        reference = load_reference_annotation(reference_file)
    
    # 质量评估
    assessor = AugustusQualityAssessor()
    metrics = assessor.assess_prediction_quality(predictions, reference, genome_sequence)
    
    # 导出报告
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        
        base_name = Path(predictions_file).stem
        
        # 导出多种格式
        assessor.export_quality_report(
            metrics, 
            os.path.join(output_dir, f"{base_name}_quality_report.json"),
            'json'
        )
        
        assessor.export_quality_report(
            metrics,
            os.path.join(output_dir, f"{base_name}_quality_report.tsv"),
            'tsv'
        )
        
        assessor.export_quality_report(
            metrics,
            os.path.join(output_dir, f"{base_name}_quality_report.txt"),
            'txt'
        )
    
    return metrics 