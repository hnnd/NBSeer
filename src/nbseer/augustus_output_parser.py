"""
Augustus输出解析器
Augustus Output Parser

本模块提供对Augustus预测输出的详细解析和基因结构信息提取：
- GFF3格式解析
- 基因结构分析
- 外显子/内含子边界识别
- 蛋白质序列提取
- 统计分析和质量评估
"""

import os
import re
import logging
from pathlib import Path
from typing import List, Dict, Optional, Tuple, Union
from dataclasses import dataclass, field
from collections import defaultdict
import json

try:
    import gffutils
    GFFUTILS_AVAILABLE = True
except ImportError:
    GFFUTILS_AVAILABLE = False

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


@dataclass
class ExonInfo:
    """外显子信息"""
    start: int
    end: int
    length: int
    phase: int = 0
    score: float = 0.0
    strand: str = "+"
    
    def __post_init__(self):
        if self.length <= 0:
            self.length = self.end - self.start + 1


@dataclass
class IntronInfo:
    """内含子信息"""
    start: int
    end: int
    length: int
    donor_site: str = ""  # GT, GC等
    acceptor_site: str = ""  # AG等
    strand: str = "+"
    
    def __post_init__(self):
        if self.length <= 0:
            self.length = self.end - self.start + 1


@dataclass
class GeneStructure:
    """基因结构信息"""
    gene_id: str
    seqname: str
    start: int
    end: int
    strand: str
    score: float = 0.0
    
    # 基因组件
    exons: List[ExonInfo] = field(default_factory=list)
    introns: List[IntronInfo] = field(default_factory=list)
    cds_regions: List[Tuple[int, int]] = field(default_factory=list)
    
    # 序列信息
    gene_sequence: str = ""
    cds_sequence: str = ""
    protein_sequence: str = ""
    
    # 统计信息
    gene_length: int = 0
    cds_length: int = 0
    exon_count: int = 0
    intron_count: int = 0
    
    # 质量指标
    coding_potential: float = 0.0
    gc_content: float = 0.0
    
    def __post_init__(self):
        if self.gene_length <= 0:
            self.gene_length = self.end - self.start + 1
        
        if not self.exon_count:
            self.exon_count = len(self.exons)
        
        if not self.intron_count:
            self.intron_count = len(self.introns)


@dataclass
class AugustusParsingResult:
    """Augustus解析结果"""
    genes: List[GeneStructure] = field(default_factory=list)
    total_genes: int = 0
    parsing_errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    
    # 统计信息
    statistics: Dict = field(default_factory=dict)
    
    def __post_init__(self):
        if not self.total_genes:
            self.total_genes = len(self.genes)


class AugustusOutputParser:
    """Augustus输出解析器"""
    
    def __init__(self, use_gffutils: bool = True):
        """
        初始化解析器
        
        Args:
            use_gffutils: 是否使用gffutils库进行高级解析
        """
        self.use_gffutils = use_gffutils and GFFUTILS_AVAILABLE
        self.logger = logging.getLogger(__name__)
        
        if not self.use_gffutils:
            self.logger.info("使用基本GFF3解析模式")
    
    def parse_augustus_output(self, gff3_file: str, 
                            fasta_file: str = None,
                            protein_file: str = None) -> AugustusParsingResult:
        """
        解析Augustus完整输出
        
        Args:
            gff3_file: GFF3预测文件路径
            fasta_file: 输入序列文件路径（用于提取序列）
            protein_file: 蛋白质序列文件路径
        
        Returns:
            AugustusParsingResult: 解析结果
        """
        self.logger.info(f"解析Augustus输出: {gff3_file}")
        
        # 解析GFF3文件
        if self.use_gffutils:
            genes = self._parse_with_gffutils(gff3_file)
        else:
            genes = self._parse_basic_gff3(gff3_file)
        
        # 如果提供了序列文件，提取序列信息
        if fasta_file and os.path.exists(fasta_file):
            self._extract_sequences(genes, fasta_file)
        
        # 如果提供了蛋白质文件，提取蛋白质序列
        if protein_file and os.path.exists(protein_file):
            self._extract_protein_sequences(genes, protein_file)
        
        # 计算统计信息
        statistics = self._calculate_statistics(genes)
        
        result = AugustusParsingResult(
            genes=genes,
            statistics=statistics
        )
        
        self.logger.info(f"解析完成: {len(genes)} 个基因")
        return result
    
    def _parse_basic_gff3(self, gff3_file: str) -> List[GeneStructure]:
        """基本GFF3解析"""
        genes = {}
        
        with open(gff3_file, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue
                
                seqname, source, feature, start, end, score, strand, phase, attributes = parts
                start, end = int(start), int(end)
                score = float(score) if score != '.' else 0.0
                phase = int(phase) if phase != '.' else 0
                
                # 解析属性
                attr_dict = self._parse_attributes(attributes)
                
                if feature == 'gene':
                    gene_id = attr_dict.get('ID', f'gene_{len(genes)}')
                    genes[gene_id] = {
                        'gene_id': gene_id,
                        'seqname': seqname,
                        'start': start,
                        'end': end,
                        'strand': strand,
                        'score': score,
                        'exons': [],
                        'cds_regions': []
                    }
                
                elif feature in ['exon', 'CDS']:
                    # 查找父基因
                    parent_id = attr_dict.get('Parent', '')
                    # 尝试多种父ID匹配方式
                    gene_key = None
                    for gene_id in genes.keys():
                        if parent_id == gene_id or parent_id.startswith(gene_id):
                            gene_key = gene_id
                            break
                    
                    if gene_key:
                        gene_data = genes[gene_key]
                        
                        if feature == 'exon':
                            exon = {
                                'start': start,
                                'end': end,
                                'length': end - start + 1,
                                'phase': phase,
                                'score': score,
                                'strand': strand
                            }
                            gene_data['exons'].append(exon)
                        
                        elif feature == 'CDS':
                            gene_data['cds_regions'].append((start, end))
        
        # 转换为GeneStructure对象并计算内含子信息
        gene_structures = []
        for gene_data in genes.values():
            gene = GeneStructure(
                gene_id=gene_data['gene_id'],
                seqname=gene_data['seqname'],
                start=gene_data['start'],
                end=gene_data['end'],
                strand=gene_data['strand'],
                score=gene_data['score']
            )
            
            # 添加外显子
            for exon_data in gene_data['exons']:
                exon = ExonInfo(
                    start=exon_data['start'],
                    end=exon_data['end'],
                    length=exon_data['length'],
                    phase=exon_data['phase'],
                    score=exon_data['score'],
                    strand=exon_data['strand']
                )
                gene.exons.append(exon)
            
            # 添加CDS区域
            gene.cds_regions = gene_data['cds_regions']
            
            # 更新计数
            gene.exon_count = len(gene.exons)
            gene.intron_count = 0  # 将在下面计算
            
            # 计算内含子
            if len(gene.exons) > 1:
                # 按位置排序外显子
                gene.exons.sort(key=lambda x: x.start)
                
                # 计算内含子
                for i in range(len(gene.exons) - 1):
                    intron_start = gene.exons[i].end + 1
                    intron_end = gene.exons[i + 1].start - 1
                    
                    if intron_end > intron_start:
                        intron = IntronInfo(
                            start=intron_start,
                            end=intron_end,
                            length=intron_end - intron_start + 1,
                            strand=gene.strand
                        )
                        gene.introns.append(intron)
            
            # 更新内含子计数
            gene.intron_count = len(gene.introns)
            
            gene_structures.append(gene)
        
        return gene_structures
    
    def _parse_with_gffutils(self, gff3_file: str) -> List[GeneStructure]:
        """使用gffutils进行高级解析"""
        try:
            # 创建临时数据库
            db_path = gff3_file + '.db'
            if os.path.exists(db_path):
                os.remove(db_path)
            
            db = gffutils.create_db(gff3_file, db_path, merge_strategy='create_unique')
            
            genes = []
            
            for gene_feature in db.features_of_type('gene'):
                gene = GeneStructure(
                    gene_id=gene_feature.id,
                    seqname=gene_feature.seqid,
                    start=gene_feature.start,
                    end=gene_feature.end,
                    strand=gene_feature.strand,
                    score=float(gene_feature.score) if gene_feature.score != '.' else 0.0
                )
                
                # 获取外显子
                for exon in db.children(gene_feature, featuretype='exon'):
                    exon_info = ExonInfo(
                        start=exon.start,
                        end=exon.end,
                        length=exon.end - exon.start + 1,
                        phase=int(exon.phase) if exon.phase != '.' else 0,
                        score=float(exon.score) if exon.score != '.' else 0.0,
                        strand=exon.strand
                    )
                    gene.exons.append(exon_info)
                
                # 获取CDS区域
                for cds in db.children(gene_feature, featuretype='CDS'):
                    gene.cds_regions.append((cds.start, cds.end))
                
                # 计算内含子
                if len(gene.exons) > 1:
                    gene.exons.sort(key=lambda x: x.start)
                    for i in range(len(gene.exons) - 1):
                        intron_start = gene.exons[i].end + 1
                        intron_end = gene.exons[i + 1].start - 1
                        
                        if intron_end > intron_start:
                            intron = IntronInfo(
                                start=intron_start,
                                end=intron_end,
                                length=intron_end - intron_start + 1,
                                strand=gene.strand
                            )
                            gene.introns.append(intron)
                
                genes.append(gene)
            
            # 清理临时数据库
            if os.path.exists(db_path):
                os.remove(db_path)
            
            return genes
            
        except Exception as e:
            self.logger.warning(f"gffutils解析失败，回退到基本解析: {e}")
            return self._parse_basic_gff3(gff3_file)
    
    def _parse_attributes(self, attributes_str: str) -> Dict[str, str]:
        """解析GFF3属性字符串"""
        attrs = {}
        for attr in attributes_str.split(';'):
            if '=' in attr:
                key, value = attr.split('=', 1)
                attrs[key] = value
        return attrs
    
    def _extract_sequences(self, genes: List[GeneStructure], fasta_file: str):
        """从FASTA文件中提取基因序列"""
        self.logger.info(f"从 {fasta_file} 提取基因序列")
        
        # 读取序列
        sequences = {}
        for record in SeqIO.parse(fasta_file, 'fasta'):
            sequences[record.id] = record.seq
        
        for gene in genes:
            if gene.seqname in sequences:
                seq = sequences[gene.seqname]
                
                # 提取基因序列（注意转换为0-based索引）
                gene.gene_sequence = str(seq[gene.start-1:gene.end])
                
                # 提取CDS序列
                if gene.cds_regions:
                    cds_parts = []
                    for start, end in sorted(gene.cds_regions):
                        cds_parts.append(str(seq[start-1:end]))
                    gene.cds_sequence = ''.join(cds_parts)
                    gene.cds_length = len(gene.cds_sequence)
                
                # 计算GC含量
                if gene.gene_sequence:
                    gc_count = gene.gene_sequence.count('G') + gene.gene_sequence.count('C')
                    gene.gc_content = gc_count / len(gene.gene_sequence) if gene.gene_sequence else 0
    
    def _extract_protein_sequences(self, genes: List[GeneStructure], protein_file: str):
        """从蛋白质文件中提取蛋白质序列"""
        self.logger.info(f"从 {protein_file} 提取蛋白质序列")
        
        proteins = {}
        for record in SeqIO.parse(protein_file, 'fasta'):
            # Augustus蛋白质序列ID通常包含基因ID
            gene_id = record.id.split('.')[0]  # 去掉可能的转录本后缀
            proteins[gene_id] = str(record.seq)
        
        for gene in genes:
            if gene.gene_id in proteins:
                gene.protein_sequence = proteins[gene.gene_id]
    
    def _calculate_statistics(self, genes: List[GeneStructure]) -> Dict:
        """计算基因结构统计信息"""
        if not genes:
            return {}
        
        stats = {
            'total_genes': len(genes),
            'gene_lengths': [g.gene_length for g in genes],
            'exon_counts': [g.exon_count for g in genes],
            'intron_counts': [g.intron_count for g in genes],
            'cds_lengths': [g.cds_length for g in genes if g.cds_length > 0],
            'gc_contents': [g.gc_content for g in genes if g.gc_content > 0],
        }
        
        # 计算平均值和范围
        for key in ['gene_lengths', 'exon_counts', 'intron_counts', 'cds_lengths', 'gc_contents']:
            values = stats[key]
            if values:
                stats[f'{key}_mean'] = sum(values) / len(values)
                stats[f'{key}_min'] = min(values)
                stats[f'{key}_max'] = max(values)
                stats[f'{key}_median'] = sorted(values)[len(values)//2]
        
        # 基因结构分布
        stats['single_exon_genes'] = len([g for g in genes if g.exon_count == 1])
        stats['multi_exon_genes'] = len([g for g in genes if g.exon_count > 1])
        
        # 链分布
        stats['plus_strand_genes'] = len([g for g in genes if g.strand == '+'])
        stats['minus_strand_genes'] = len([g for g in genes if g.strand == '-'])
        
        return stats
    
    def export_gene_structures(self, result: AugustusParsingResult, 
                             output_file: str, format: str = 'json'):
        """
        导出基因结构信息
        
        Args:
            result: 解析结果
            output_file: 输出文件路径
            format: 输出格式 ('json', 'tsv', 'gff3')
        """
        if format.lower() == 'json':
            self._export_json(result, output_file)
        elif format.lower() == 'tsv':
            self._export_tsv(result, output_file)
        elif format.lower() == 'gff3':
            self._export_gff3(result, output_file)
        else:
            raise ValueError(f"Unsupported format: {format}")
    
    def _export_json(self, result: AugustusParsingResult, output_file: str):
        """导出JSON格式"""
        data = {
            'summary': result.statistics,
            'genes': []
        }
        
        for gene in result.genes:
            gene_data = {
                'gene_id': gene.gene_id,
                'seqname': gene.seqname,
                'start': gene.start,
                'end': gene.end,
                'strand': gene.strand,
                'length': gene.gene_length,
                'exon_count': gene.exon_count,
                'intron_count': gene.intron_count,
                'cds_length': gene.cds_length,
                'gc_content': gene.gc_content,
                'exons': [
                    {
                        'start': e.start,
                        'end': e.end,
                        'length': e.length,
                        'phase': e.phase
                    } for e in gene.exons
                ],
                'introns': [
                    {
                        'start': i.start,
                        'end': i.end,
                        'length': i.length
                    } for i in gene.introns
                ]
            }
            
            if gene.protein_sequence:
                gene_data['protein_length'] = len(gene.protein_sequence)
            
            data['genes'].append(gene_data)
        
        with open(output_file, 'w') as f:
            json.dump(data, f, indent=2)
        
        self.logger.info(f"基因结构信息已导出到: {output_file}")
    
    def _export_tsv(self, result: AugustusParsingResult, output_file: str):
        """导出TSV格式"""
        with open(output_file, 'w') as f:
            # 写入头部
            headers = [
                'gene_id', 'seqname', 'start', 'end', 'strand', 'length',
                'exon_count', 'intron_count', 'cds_length', 'gc_content',
                'protein_length', 'score'
            ]
            f.write('\t'.join(headers) + '\n')
            
            # 写入数据
            for gene in result.genes:
                row = [
                    gene.gene_id,
                    gene.seqname,
                    str(gene.start),
                    str(gene.end),
                    gene.strand,
                    str(gene.gene_length),
                    str(gene.exon_count),
                    str(gene.intron_count),
                    str(gene.cds_length),
                    f"{gene.gc_content:.3f}",
                    str(len(gene.protein_sequence)) if gene.protein_sequence else '0',
                    f"{gene.score:.3f}"
                ]
                f.write('\t'.join(row) + '\n')
        
        self.logger.info(f"基因结构TSV文件已导出到: {output_file}")
    
    def _export_gff3(self, result: AugustusParsingResult, output_file: str):
        """导出GFF3格式（重新格式化）"""
        with open(output_file, 'w') as f:
            f.write("##gff-version 3\n")
            f.write("# Processed Augustus gene predictions\n")
            
            for gene in result.genes:
                # 基因行
                attrs = f"ID={gene.gene_id};Name={gene.gene_id}"
                f.write(f"{gene.seqname}\tAugustus\tgene\t{gene.start}\t{gene.end}\t"
                       f"{gene.score:.3f}\t{gene.strand}\t.\t{attrs}\n")
                
                # 外显子行
                for i, exon in enumerate(gene.exons):
                    exon_id = f"{gene.gene_id}.exon{i+1}"
                    attrs = f"ID={exon_id};Parent={gene.gene_id}"
                    f.write(f"{gene.seqname}\tAugustus\texon\t{exon.start}\t{exon.end}\t"
                           f"{exon.score:.3f}\t{exon.strand}\t{exon.phase}\t{attrs}\n")
        
        self.logger.info(f"基因结构GFF3文件已导出到: {output_file}")


def analyze_augustus_directory(augustus_dir: str, output_dir: str = None) -> Dict:
    """
    批量分析Augustus输出目录
    
    Args:
        augustus_dir: Augustus输出目录
        output_dir: 分析结果输出目录
    
    Returns:
        分析结果字典
    """
    logger = logging.getLogger(__name__)
    logger.info(f"分析Augustus输出目录: {augustus_dir}")
    
    if output_dir is None:
        output_dir = os.path.join(augustus_dir, "analysis")
    
    os.makedirs(output_dir, exist_ok=True)
    
    parser = AugustusOutputParser()
    all_results = []
    
    # 查找所有GFF3文件
    gff3_files = list(Path(augustus_dir).glob("*.gff3"))
    logger.info(f"找到 {len(gff3_files)} 个GFF3文件")
    
    for gff3_file in gff3_files:
        try:
            # 查找对应的序列文件
            fasta_file = gff3_file.with_suffix('.fasta')
            if not fasta_file.exists():
                fasta_file = None
            
            # 查找对应的蛋白质文件
            protein_file = gff3_file.with_suffix('.aa')
            if not protein_file.exists():
                protein_file = None
            
            # 解析
            result = parser.parse_augustus_output(
                str(gff3_file), 
                str(fasta_file) if fasta_file else None,
                str(protein_file) if protein_file else None
            )
            
            all_results.append({
                'file': gff3_file.name,
                'result': result
            })
            
            # 导出单个文件的结果
            base_name = gff3_file.stem
            parser.export_gene_structures(
                result, 
                os.path.join(output_dir, f"{base_name}_structures.json")
            )
            parser.export_gene_structures(
                result, 
                os.path.join(output_dir, f"{base_name}_structures.tsv"),
                format='tsv'
            )
            
        except Exception as e:
            logger.error(f"解析 {gff3_file} 失败: {e}")
    
    # 生成汇总报告
    summary = generate_summary_report(all_results)
    
    summary_file = os.path.join(output_dir, "augustus_analysis_summary.json")
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)
    
    logger.info(f"分析完成，结果保存在: {output_dir}")
    return summary


def generate_summary_report(results: List[Dict]) -> Dict:
    """生成汇总分析报告"""
    summary = {
        'total_files': len(results),
        'total_genes': 0,
        'files_with_genes': 0,
        'files_without_genes': 0,
        'gene_statistics': {},
        'file_details': []
    }
    
    all_gene_lengths = []
    all_exon_counts = []
    all_cds_lengths = []
    
    for item in results:
        file_name = item['file']
        result = item['result']
        
        gene_count = len(result.genes)
        summary['total_genes'] += gene_count
        
        if gene_count > 0:
            summary['files_with_genes'] += 1
        else:
            summary['files_without_genes'] += 1
        
        # 收集统计数据
        for gene in result.genes:
            all_gene_lengths.append(gene.gene_length)
            all_exon_counts.append(gene.exon_count)
            if gene.cds_length > 0:
                all_cds_lengths.append(gene.cds_length)
        
        summary['file_details'].append({
            'file': file_name,
            'genes': gene_count,
            'statistics': result.statistics
        })
    
    # 计算全局统计
    if all_gene_lengths:
        summary['gene_statistics'] = {
            'gene_length_mean': sum(all_gene_lengths) / len(all_gene_lengths),
            'gene_length_min': min(all_gene_lengths),
            'gene_length_max': max(all_gene_lengths),
            'exon_count_mean': sum(all_exon_counts) / len(all_exon_counts),
            'exon_count_min': min(all_exon_counts),
            'exon_count_max': max(all_exon_counts),
        }
        
        if all_cds_lengths:
            summary['gene_statistics'].update({
                'cds_length_mean': sum(all_cds_lengths) / len(all_cds_lengths),
                'cds_length_min': min(all_cds_lengths),
                'cds_length_max': max(all_cds_lengths),
            })
    
    return summary 