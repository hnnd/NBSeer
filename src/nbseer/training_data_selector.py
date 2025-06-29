#!/usr/bin/env python3
"""
训练数据选择器 - 从miniprot结果中筛选高质量基因结构用于Augustus训练
"""

import logging
from pathlib import Path
from collections import defaultdict, Counter
from typing import List, Dict, Tuple, Optional, Set
import re
from dataclasses import dataclass

@dataclass
class GeneStructure:
    """基因结构数据类"""
    gene_id: str
    chromosome: str
    start: int
    end: int
    strand: str
    identity: float
    positive: float
    score: float
    cds_regions: List[Tuple[int, int, int]]  # (start, end, phase)
    
    @property
    def length(self) -> int:
        """基因总长度"""
        return self.end - self.start + 1
    
    @property
    def exon_count(self) -> int:
        """外显子数量"""
        return len(self.cds_regions)
    
    @property
    def cds_length(self) -> int:
        """CDS总长度"""
        return sum(end - start + 1 for start, end, _ in self.cds_regions)

class TrainingDataSelector:
    """训练数据选择器类"""
    
    def __init__(self, 
                 min_identity: float = 0.6,
                 min_positive: float = 0.7,
                 min_gene_length: int = 300,
                 max_gene_length: int = 50000,
                 min_cds_length: int = 150,
                 min_exon_count: int = 1,
                 max_exon_count: int = 50,
                 min_exon_length: int = 30,
                 overlap_threshold: float = 0.8):
        """
        初始化训练数据选择器
        
        Args:
            min_identity: 最小身份相似度阈值
            min_positive: 最小正向得分阈值
            min_gene_length: 最小基因长度
            max_gene_length: 最大基因长度
            min_cds_length: 最小CDS长度
            min_exon_count: 最小外显子数量
            max_exon_count: 最大外显子数量
            min_exon_length: 最小外显子长度
            overlap_threshold: 重叠阈值（用于去重）
        """
        self.min_identity = min_identity
        self.min_positive = min_positive
        self.min_gene_length = min_gene_length
        self.max_gene_length = max_gene_length
        self.min_cds_length = min_cds_length
        self.min_exon_count = min_exon_count
        self.max_exon_count = max_exon_count
        self.min_exon_length = min_exon_length
        self.overlap_threshold = overlap_threshold
        
        self.logger = logging.getLogger(__name__)
        
        # 统计信息
        self.stats = {
            'total_parsed': 0,
            'passed_basic_filters': 0,
            'passed_structure_filters': 0,
            'passed_quality_filters': 0,
            'final_selected': 0,
            'removed_overlaps': 0
        }
    
    def parse_gff_file(self, gff_file: Path) -> List[GeneStructure]:
        """
        解析GFF文件，提取基因结构信息
        
        Args:
            gff_file: GFF文件路径
            
        Returns:
            基因结构列表
        """
        gene_structures = []
        current_gene = None
        current_cds = []
        
        self.logger.info(f"开始解析GFF文件: {gff_file}")
        
        with open(gff_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                if line_num % 1000000 == 0:
                    self.logger.info(f"已处理 {line_num:,} 行...")
                
                if line.startswith('#') or not line.strip():
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) != 9:
                    continue
                
                chrom, source, feature_type, start, end, score, strand, phase, attributes = fields
                start, end = int(start), int(end)
                
                # 解析属性
                attr_dict = self._parse_attributes(attributes)
                
                if feature_type == 'mRNA':
                    # 保存前一个基因
                    if current_gene and current_cds:
                        gene_structure = self._create_gene_structure(
                            current_gene, current_cds
                        )
                        if gene_structure:
                            gene_structures.append(gene_structure)
                            self.stats['total_parsed'] += 1
                    
                    # 开始新基因
                    current_gene = {
                        'id': attr_dict.get('ID', ''),
                        'chromosome': chrom,
                        'start': start,
                        'end': end,
                        'strand': strand,
                        'identity': float(attr_dict.get('Identity', 0)),
                        'positive': float(attr_dict.get('Positive', 0)),
                        'score': float(score) if score != '.' else 0
                    }
                    current_cds = []
                
                elif feature_type == 'CDS' and current_gene:
                    phase_val = 0 if phase == '.' else int(phase)
                    current_cds.append((start, end, phase_val))
        
        # 处理最后一个基因
        if current_gene and current_cds:
            gene_structure = self._create_gene_structure(current_gene, current_cds)
            if gene_structure:
                gene_structures.append(gene_structure)
                self.stats['total_parsed'] += 1
        
        self.logger.info(f"解析完成，共提取 {len(gene_structures)} 个基因结构")
        return gene_structures
    
    def _parse_attributes(self, attr_string: str) -> Dict[str, str]:
        """解析GFF属性字段"""
        attributes = {}
        for item in attr_string.split(';'):
            if '=' in item:
                key, value = item.split('=', 1)
                attributes[key] = value
        return attributes
    
    def _create_gene_structure(self, gene_info: Dict, cds_list: List[Tuple[int, int, int]]) -> Optional[GeneStructure]:
        """创建基因结构对象"""
        try:
            return GeneStructure(
                gene_id=gene_info['id'],
                chromosome=gene_info['chromosome'],
                start=gene_info['start'],
                end=gene_info['end'],
                strand=gene_info['strand'],
                identity=gene_info['identity'],
                positive=gene_info['positive'],
                score=gene_info['score'],
                cds_regions=sorted(cds_list)
            )
        except Exception as e:
            self.logger.warning(f"创建基因结构失败: {e}")
            return None
    
    def apply_quality_filters(self, gene_structures: List[GeneStructure]) -> List[GeneStructure]:
        """
        应用质量过滤器
        
        Args:
            gene_structures: 输入基因结构列表
            
        Returns:
            过滤后的基因结构列表
        """
        filtered = []
        
        for gene in gene_structures:
            # 基本质量过滤
            if not self._passes_basic_quality(gene):
                continue
            self.stats['passed_basic_filters'] += 1
            
            # 结构合理性过滤
            if not self._passes_structure_validation(gene):
                continue
            self.stats['passed_structure_filters'] += 1
            
            # 高级质量过滤
            if not self._passes_advanced_quality(gene):
                continue
            self.stats['passed_quality_filters'] += 1
            
            filtered.append(gene)
        
        self.logger.info(f"质量过滤后保留 {len(filtered)} 个基因结构")
        return filtered
    
    def _passes_basic_quality(self, gene: GeneStructure) -> bool:
        """基本质量检查"""
        return (
            gene.identity >= self.min_identity and
            gene.positive >= self.min_positive and
            gene.length >= self.min_gene_length and
            gene.length <= self.max_gene_length and
            gene.cds_length >= self.min_cds_length
        )
    
    def _passes_structure_validation(self, gene: GeneStructure) -> bool:
        """结构合理性检查"""
        # 外显子数量检查
        if not (self.min_exon_count <= gene.exon_count <= self.max_exon_count):
            return False
        
        # 外显子长度检查
        for start, end, _ in gene.cds_regions:
            exon_length = end - start + 1
            if exon_length < self.min_exon_length:
                return False
        
        # CDS区域不能重叠
        sorted_cds = sorted(gene.cds_regions)
        for i in range(len(sorted_cds) - 1):
            if sorted_cds[i][1] >= sorted_cds[i+1][0]:
                return False
        
        return True
    
    def _passes_advanced_quality(self, gene: GeneStructure) -> bool:
        """高级质量检查"""
        # CDS长度应该是3的倍数（对于完整基因）
        if gene.cds_length % 3 != 0:
            # 允许部分基因不是3的倍数（可能是部分序列）
            pass
        
        # 检查外显子长度分布的合理性
        exon_lengths = [end - start + 1 for start, end, _ in gene.cds_regions]
        avg_exon_length = sum(exon_lengths) / len(exon_lengths)
        
        # 过滤异常短的平均外显子长度
        if avg_exon_length < 50:
            return False
        
        return True
    
    def remove_overlapping_genes(self, gene_structures: List[GeneStructure]) -> List[GeneStructure]:
        """
        移除重叠的基因，保留质量最高的
        
        Args:
            gene_structures: 输入基因结构列表
            
        Returns:
            去重后的基因结构列表
        """
        # 按染色体分组
        chrom_genes = defaultdict(list)
        for gene in gene_structures:
            chrom_genes[gene.chromosome].append(gene)
        
        final_genes = []
        
        for chrom, genes in chrom_genes.items():
            # 按位置排序
            genes.sort(key=lambda g: (g.start, g.end))
            
            # 贪心算法选择非重叠的最优基因
            selected = []
            for gene in genes:
                # 检查是否与已选择的基因重叠
                overlaps = False
                for selected_gene in selected:
                    if self._calculate_overlap(gene, selected_gene) > self.overlap_threshold:
                        overlaps = True
                        # 如果当前基因质量更高，替换已选择的基因
                        if self._compare_gene_quality(gene, selected_gene) > 0:
                            selected.remove(selected_gene)
                            selected.append(gene)
                            self.stats['removed_overlaps'] += 1
                        break
                
                if not overlaps:
                    selected.append(gene)
            
            final_genes.extend(selected)
        
        self.logger.info(f"去重后保留 {len(final_genes)} 个基因结构")
        self.stats['final_selected'] = len(final_genes)
        return final_genes
    
    def _calculate_overlap(self, gene1: GeneStructure, gene2: GeneStructure) -> float:
        """计算两个基因的重叠度"""
        if gene1.chromosome != gene2.chromosome or gene1.strand != gene2.strand:
            return 0.0
        
        overlap_start = max(gene1.start, gene2.start)
        overlap_end = min(gene1.end, gene2.end)
        
        if overlap_start > overlap_end:
            return 0.0
        
        overlap_length = overlap_end - overlap_start + 1
        min_length = min(gene1.length, gene2.length)
        
        return overlap_length / min_length
    
    def _compare_gene_quality(self, gene1: GeneStructure, gene2: GeneStructure) -> int:
        """
        比较两个基因的质量
        
        Returns:
            1 if gene1 > gene2, -1 if gene1 < gene2, 0 if equal
        """
        # 首先比较身份相似度
        if abs(gene1.identity - gene2.identity) > 0.01:
            return 1 if gene1.identity > gene2.identity else -1
        
        # 然后比较得分
        if abs(gene1.score - gene2.score) > 1:
            return 1 if gene1.score > gene2.score else -1
        
        # 最后比较正向得分
        if abs(gene1.positive - gene2.positive) > 0.01:
            return 1 if gene1.positive > gene2.positive else -1
        
        return 0
    
    def select_training_data(self, gff_file: Path, max_genes: Optional[int] = None) -> List[GeneStructure]:
        """
        从GFF文件中选择训练数据
        
        Args:
            gff_file: 输入GFF文件路径
            max_genes: 最大基因数量限制
            
        Returns:
            选择的高质量基因结构列表
        """
        self.logger.info("开始训练数据选择流程")
        
        # 重置统计信息
        self.stats = {k: 0 for k in self.stats.keys()}
        
        # 1. 解析GFF文件
        gene_structures = self.parse_gff_file(gff_file)
        
        # 2. 应用质量过滤器
        filtered_genes = self.apply_quality_filters(gene_structures)
        
        # 3. 移除重叠基因
        final_genes = self.remove_overlapping_genes(filtered_genes)
        
        # 4. 按质量排序并限制数量
        final_genes.sort(key=lambda g: (g.identity, g.score), reverse=True)
        
        if max_genes and len(final_genes) > max_genes:
            final_genes = final_genes[:max_genes]
            self.logger.info(f"限制为前 {max_genes} 个最高质量基因")
        
        self.logger.info("训练数据选择完成")
        self._log_statistics()
        
        return final_genes
    
    def _log_statistics(self):
        """记录统计信息"""
        self.logger.info("=== 训练数据选择统计 ===")
        self.logger.info(f"总解析基因数: {self.stats['total_parsed']:,}")
        self.logger.info(f"通过基本过滤: {self.stats['passed_basic_filters']:,}")
        self.logger.info(f"通过结构过滤: {self.stats['passed_structure_filters']:,}")
        self.logger.info(f"通过质量过滤: {self.stats['passed_quality_filters']:,}")
        self.logger.info(f"移除重叠基因: {self.stats['removed_overlaps']:,}")
        self.logger.info(f"最终选择基因: {self.stats['final_selected']:,}")
        
        if self.stats['total_parsed'] > 0:
            retention_rate = self.stats['final_selected'] / self.stats['total_parsed'] * 100
            self.logger.info(f"保留率: {retention_rate:.2f}%")
    
    def generate_statistics_report(self, gene_structures: List[GeneStructure]) -> Dict:
        """
        生成训练数据统计报告
        
        Args:
            gene_structures: 基因结构列表
            
        Returns:
            统计报告字典
        """
        if not gene_structures:
            return {}
        
        # 基本统计
        total_genes = len(gene_structures)
        gene_lengths = [g.length for g in gene_structures]
        cds_lengths = [g.cds_length for g in gene_structures]
        exon_counts = [g.exon_count for g in gene_structures]
        identity_scores = [g.identity for g in gene_structures]
        
        # 染色体分布
        chrom_dist = Counter(g.chromosome for g in gene_structures)
        
        # 链方向分布
        strand_dist = Counter(g.strand for g in gene_structures)
        
        # 外显子数量分布
        exon_dist = Counter(exon_counts)
        
        report = {
            'basic_stats': {
                'total_genes': total_genes,
                'avg_gene_length': sum(gene_lengths) / total_genes,
                'avg_cds_length': sum(cds_lengths) / total_genes,
                'avg_exon_count': sum(exon_counts) / total_genes,
                'avg_identity': sum(identity_scores) / total_genes
            },
            'length_distribution': {
                'min_gene_length': min(gene_lengths),
                'max_gene_length': max(gene_lengths),
                'median_gene_length': sorted(gene_lengths)[total_genes//2],
                'min_cds_length': min(cds_lengths),
                'max_cds_length': max(cds_lengths),
                'median_cds_length': sorted(cds_lengths)[total_genes//2]
            },
            'quality_distribution': {
                'min_identity': min(identity_scores),
                'max_identity': max(identity_scores),
                'median_identity': sorted(identity_scores)[total_genes//2],
                'high_quality_count': sum(1 for i in identity_scores if i >= 0.8),
                'medium_quality_count': sum(1 for i in identity_scores if 0.6 <= i < 0.8)
            },
            'chromosome_distribution': dict(chrom_dist.most_common()),
            'strand_distribution': dict(strand_dist),
            'exon_count_distribution': dict(exon_dist.most_common(10))
        }
        
        return report
    
    def save_filtered_gff(self, gene_structures: List[GeneStructure], output_file: Path):
        """
        保存过滤后的基因结构为GFF格式
        
        Args:
            gene_structures: 基因结构列表
            output_file: 输出文件路径
        """
        self.logger.info(f"保存过滤后的GFF文件: {output_file}")
        
        with open(output_file, 'w') as f:
            f.write("##gff-version 3\n")
            
            for gene in gene_structures:
                # 写入mRNA行
                f.write(f"{gene.chromosome}\tminiprot\tmRNA\t{gene.start}\t{gene.end}\t"
                       f"{gene.score:.1f}\t{gene.strand}\t.\t"
                       f"ID={gene.gene_id};Identity={gene.identity:.4f};Positive={gene.positive:.4f}\n")
                
                # 写入CDS行
                for i, (start, end, phase) in enumerate(gene.cds_regions):
                    f.write(f"{gene.chromosome}\tminiprot\tCDS\t{start}\t{end}\t"
                           f"{gene.score:.1f}\t{gene.strand}\t{phase}\t"
                           f"Parent={gene.gene_id}\n")
        
        self.logger.info(f"已保存 {len(gene_structures)} 个基因结构到GFF文件") 