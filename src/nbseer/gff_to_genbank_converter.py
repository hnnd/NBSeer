#!/usr/bin/env python3
"""
GFF到GenBank格式转换器 - 为Augustus训练准备数据
"""

import logging
from pathlib import Path
from typing import List, Dict, Tuple, Optional
import re
from dataclasses import dataclass
from datetime import datetime
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import generic_dna

@dataclass
class GeneInfo:
    """基因信息数据类"""
    gene_id: str
    chromosome: str
    start: int
    end: int
    strand: str
    cds_regions: List[Tuple[int, int, int]]  # (start, end, phase)
    
    @property
    def length(self) -> int:
        return self.end - self.start + 1

class GFFToGenBankConverter:
    """GFF到GenBank格式转换器"""
    
    def __init__(self, 
                 flanking_length: int = 1000,
                 min_gene_length: int = 200,
                 max_gene_length: int = 50000):
        """
        初始化转换器
        
        Args:
            flanking_length: 基因上下游包含的序列长度
            min_gene_length: 最小基因长度
            max_gene_length: 最大基因长度
        """
        self.flanking_length = flanking_length
        self.min_gene_length = min_gene_length
        self.max_gene_length = max_gene_length
        
        self.logger = logging.getLogger(__name__)
        
        # 统计信息
        self.stats = {
            'total_genes': 0,
            'valid_genes': 0,
            'filtered_genes': 0,
            'sequences_written': 0,
            'total_length': 0
        }
    
    def parse_gff_file(self, gff_file: Path) -> List[GeneInfo]:
        """
        解析GFF文件，提取基因信息
        
        Args:
            gff_file: GFF文件路径
            
        Returns:
            基因信息列表
        """
        genes = []
        current_gene = None
        current_cds = []
        
        self.logger.info(f"解析GFF文件: {gff_file}")
        
        with open(gff_file, 'r') as f:
            for line in f:
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
                        gene_info = self._create_gene_info(current_gene, current_cds)
                        if gene_info:
                            genes.append(gene_info)
                            self.stats['total_genes'] += 1
                    
                    # 开始新基因
                    current_gene = {
                        'id': attr_dict.get('ID', ''),
                        'chromosome': chrom,
                        'start': start,
                        'end': end,
                        'strand': strand
                    }
                    current_cds = []
                
                elif feature_type == 'CDS' and current_gene:
                    phase_val = 0 if phase == '.' else int(phase)
                    current_cds.append((start, end, phase_val))
        
        # 处理最后一个基因
        if current_gene and current_cds:
            gene_info = self._create_gene_info(current_gene, current_cds)
            if gene_info:
                genes.append(gene_info)
                self.stats['total_genes'] += 1
        
        self.logger.info(f"解析完成，共提取 {len(genes)} 个基因")
        return genes
    
    def _parse_attributes(self, attr_string: str) -> Dict[str, str]:
        """解析GFF属性字段"""
        attributes = {}
        for item in attr_string.split(';'):
            if '=' in item:
                key, value = item.split('=', 1)
                attributes[key] = value
        return attributes
    
    def _create_gene_info(self, gene_data: Dict, cds_list: List[Tuple[int, int, int]]) -> Optional[GeneInfo]:
        """创建基因信息对象"""
        try:
            return GeneInfo(
                gene_id=gene_data['id'],
                chromosome=gene_data['chromosome'],
                start=gene_data['start'],
                end=gene_data['end'],
                strand=gene_data['strand'],
                cds_regions=sorted(cds_list)
            )
        except Exception as e:
            self.logger.warning(f"创建基因信息失败: {e}")
            return None
    
    def load_genome_sequences(self, genome_file: Path) -> Dict[str, SeqRecord]:
        """
        加载基因组序列
        
        Args:
            genome_file: 基因组FASTA文件路径
            
        Returns:
            染色体序列字典
        """
        self.logger.info(f"加载基因组序列: {genome_file}")
        
        sequences = {}
        for record in SeqIO.parse(genome_file, "fasta"):
            sequences[record.id] = record
            self.logger.debug(f"加载染色体 {record.id}, 长度: {len(record.seq):,} bp")
        
        self.logger.info(f"加载完成，共 {len(sequences)} 条染色体序列")
        return sequences
    
    def extract_gene_sequence(self, gene: GeneInfo, genome_seqs: Dict[str, SeqRecord]) -> Optional[SeqRecord]:
        """
        提取基因序列（包含上下游序列）
        
        Args:
            gene: 基因信息
            genome_seqs: 基因组序列字典
            
        Returns:
            包含基因和上下游序列的SeqRecord
        """
        if gene.chromosome not in genome_seqs:
            self.logger.warning(f"未找到染色体序列: {gene.chromosome}")
            return None
        
        chrom_seq = genome_seqs[gene.chromosome]
        chrom_length = len(chrom_seq.seq)
        
        # 计算提取区域（包含上下游序列）
        extract_start = max(1, gene.start - self.flanking_length)
        extract_end = min(chrom_length, gene.end + self.flanking_length)
        
        # 提取序列（注意：SeqRecord使用0-based索引）
        extracted_seq = chrom_seq.seq[extract_start-1:extract_end]
        
        # 创建新的SeqRecord
        seq_id = f"{gene.chromosome}_{extract_start}_{extract_end}"
        record = SeqRecord(
            extracted_seq,
            id=seq_id,
            description=f"Gene {gene.gene_id} with flanking regions"
        )
        
        # 调整CDS坐标到新序列的相对位置
        adjusted_cds = []
        for cds_start, cds_end, phase in gene.cds_regions:
            rel_start = cds_start - extract_start + 1
            rel_end = cds_end - extract_start + 1
            
            # 确保坐标在有效范围内
            if rel_start > 0 and rel_end <= len(extracted_seq):
                adjusted_cds.append((rel_start, rel_end, phase))
        
        if not adjusted_cds:
            self.logger.warning(f"基因 {gene.gene_id} 的CDS区域超出提取范围")
            return None
        
        # 添加CDS特征
        for i, (cds_start, cds_end, phase) in enumerate(adjusted_cds):
            # 创建CDS特征（注意：FeatureLocation使用0-based索引）
            if gene.strand == '+':
                location = FeatureLocation(cds_start-1, cds_end, strand=1)
            else:
                location = FeatureLocation(cds_start-1, cds_end, strand=-1)
            
            cds_feature = SeqFeature(
                location=location,
                type="CDS",
                qualifiers={
                    'gene': [gene.gene_id],
                    'note': [f"CDS {i+1} of {len(adjusted_cds)}"],
                    'codon_start': [str(phase + 1)]
                }
            )
            record.features.append(cds_feature)
        
        # 添加基因特征
        gene_start = adjusted_cds[0][0]
        gene_end = adjusted_cds[-1][1]
        
        if gene.strand == '+':
            gene_location = FeatureLocation(gene_start-1, gene_end, strand=1)
        else:
            gene_location = FeatureLocation(gene_start-1, gene_end, strand=-1)
        
        gene_feature = SeqFeature(
            location=gene_location,
            type="gene",
            qualifiers={
                'gene': [gene.gene_id],
                'note': [f"Original location: {gene.chromosome}:{gene.start}-{gene.end}"]
            }
        )
        record.features.insert(0, gene_feature)  # 基因特征放在最前面
        
        return record
    
    def validate_gene_structure(self, gene: GeneInfo, record: SeqRecord) -> bool:
        """
        验证基因结构的有效性
        
        Args:
            gene: 基因信息
            record: 序列记录
            
        Returns:
            是否通过验证
        """
        # 检查基因长度
        if not (self.min_gene_length <= gene.length <= self.max_gene_length):
            self.logger.debug(f"基因 {gene.gene_id} 长度 {gene.length} 超出范围")
            return False
        
        # 检查CDS区域
        if not gene.cds_regions:
            self.logger.debug(f"基因 {gene.gene_id} 没有CDS区域")
            return False
        
        # 检查CDS总长度是否为3的倍数（对于完整基因）
        total_cds_length = sum(end - start + 1 for start, end, _ in gene.cds_regions)
        if total_cds_length % 3 != 0:
            self.logger.debug(f"基因 {gene.gene_id} CDS长度 {total_cds_length} 不是3的倍数")
            # 注意：这里不直接返回False，因为某些基因可能是部分序列
        
        # 检查序列质量
        if len(record.seq) < self.min_gene_length:
            self.logger.debug(f"基因 {gene.gene_id} 提取序列过短")
            return False
        
        # 检查是否包含过多的N
        n_count = record.seq.upper().count('N')
        n_ratio = n_count / len(record.seq)
        if n_ratio > 0.1:  # 超过10%的N
            self.logger.debug(f"基因 {gene.gene_id} 包含过多N ({n_ratio:.2%})")
            return False
        
        return True
    
    def convert_to_genbank(self, 
                          gff_file: Path, 
                          genome_file: Path, 
                          output_file: Path) -> bool:
        """
        将GFF文件转换为GenBank格式
        
        Args:
            gff_file: 输入GFF文件
            genome_file: 基因组FASTA文件
            output_file: 输出GenBank文件
            
        Returns:
            转换是否成功
        """
        try:
            self.logger.info("开始GFF到GenBank转换")
            
            # 重置统计信息
            self.stats = {k: 0 for k in self.stats.keys()}
            
            # 1. 解析GFF文件
            genes = self.parse_gff_file(gff_file)
            
            # 2. 加载基因组序列
            genome_seqs = self.load_genome_sequences(genome_file)
            
            # 3. 转换每个基因
            valid_records = []
            
            for gene in genes:
                # 提取基因序列
                record = self.extract_gene_sequence(gene, genome_seqs)
                if not record:
                    self.stats['filtered_genes'] += 1
                    continue
                
                # 验证基因结构
                if not self.validate_gene_structure(gene, record):
                    self.stats['filtered_genes'] += 1
                    continue
                
                valid_records.append(record)
                self.stats['valid_genes'] += 1
                self.stats['total_length'] += len(record.seq)
            
            # 4. 写入GenBank文件
            if valid_records:
                with open(output_file, 'w') as f:
                    SeqIO.write(valid_records, f, "genbank")
                
                self.stats['sequences_written'] = len(valid_records)
                self.logger.info(f"成功写入 {len(valid_records)} 个基因序列到 {output_file}")
            else:
                self.logger.warning("没有有效的基因序列可写入")
                return False
            
            # 5. 记录统计信息
            self._log_conversion_stats()
            
            return True
            
        except Exception as e:
            self.logger.error(f"转换失败: {e}")
            import traceback
            traceback.print_exc()
            return False
    
    def _log_conversion_stats(self):
        """记录转换统计信息"""
        self.logger.info("=== GFF到GenBank转换统计 ===")
        self.logger.info(f"总基因数: {self.stats['total_genes']:,}")
        self.logger.info(f"有效基因数: {self.stats['valid_genes']:,}")
        self.logger.info(f"过滤基因数: {self.stats['filtered_genes']:,}")
        self.logger.info(f"写入序列数: {self.stats['sequences_written']:,}")
        self.logger.info(f"总序列长度: {self.stats['total_length']:,} bp")
        
        if self.stats['total_genes'] > 0:
            success_rate = self.stats['valid_genes'] / self.stats['total_genes'] * 100
            self.logger.info(f"转换成功率: {success_rate:.2f}%")
        
        if self.stats['valid_genes'] > 0:
            avg_length = self.stats['total_length'] / self.stats['valid_genes']
            self.logger.info(f"平均序列长度: {avg_length:.1f} bp")
    
    def split_genbank_file(self, 
                          genbank_file: Path, 
                          train_size: int = None,
                          test_size: int = 200) -> Tuple[Path, Path]:
        """
        将GenBank文件分割为训练集和测试集
        
        Args:
            genbank_file: 输入GenBank文件
            train_size: 训练集大小（None表示使用剩余所有）
            test_size: 测试集大小
            
        Returns:
            (训练集文件路径, 测试集文件路径)
        """
        self.logger.info(f"分割GenBank文件: {genbank_file}")
        
        # 读取所有记录
        records = list(SeqIO.parse(genbank_file, "genbank"))
        total_records = len(records)
        
        self.logger.info(f"总记录数: {total_records}")
        
        if total_records < test_size:
            self.logger.warning(f"记录数 ({total_records}) 少于测试集大小 ({test_size})")
            test_size = min(test_size, total_records // 2)
        
        # 随机打乱
        import random
        random.shuffle(records)
        
        # 分割
        test_records = records[:test_size]
        train_records = records[test_size:]
        
        if train_size and len(train_records) > train_size:
            train_records = train_records[:train_size]
        
        # 生成输出文件名
        base_name = genbank_file.stem
        output_dir = genbank_file.parent
        
        train_file = output_dir / f"{base_name}_train.gb"
        test_file = output_dir / f"{base_name}_test.gb"
        
        # 写入文件
        with open(train_file, 'w') as f:
            SeqIO.write(train_records, f, "genbank")
        
        with open(test_file, 'w') as f:
            SeqIO.write(test_records, f, "genbank")
        
        self.logger.info(f"训练集: {len(train_records)} 记录 -> {train_file}")
        self.logger.info(f"测试集: {len(test_records)} 记录 -> {test_file}")
        
        return train_file, test_file
    
    def generate_conversion_report(self, output_dir: Path) -> Dict:
        """
        生成转换报告
        
        Args:
            output_dir: 输出目录
            
        Returns:
            报告字典
        """
        report = {
            'conversion_stats': self.stats.copy(),
            'parameters': {
                'flanking_length': self.flanking_length,
                'min_gene_length': self.min_gene_length,
                'max_gene_length': self.max_gene_length
            },
            'generated_at': datetime.now().isoformat()
        }
        
        # 保存报告
        report_file = output_dir / "genbank_conversion_report.json"
        import json
        with open(report_file, 'w', encoding='utf-8') as f:
            json.dump(report, f, indent=2, ensure_ascii=False)
        
        self.logger.info(f"转换报告已保存: {report_file}")
        return report 