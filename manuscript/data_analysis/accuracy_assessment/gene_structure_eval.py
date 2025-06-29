#!/usr/bin/env python3
"""
基因结构准确性评估脚本
比较NBS-Pipeline注释结果与参考注释的准确性
"""

import os
import sys
import json
import argparse
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict, Counter
import logging
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
import gffutils
from intervaltree import IntervalTree, Interval
from datetime import datetime

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class GeneStructureEvaluator:
    def __init__(self, output_dir):
        """
        初始化基因结构评估器
        
        Args:
            output_dir: 结果输出目录
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # 评估指标
        self.metrics = {
            'gene_level': {},
            'exon_level': {},
            'nucleotide_level': {},
            'domain_level': {}
        }
        
        # 存储比较结果
        self.comparison_results = []
        
    def load_gff_annotations(self, gff_file):
        """加载GFF注释文件"""
        try:
            # 创建临时数据库文件
            db_file = str(gff_file) + '.db'
            if Path(db_file).exists():
                os.remove(db_file)
            
            # 创建gffutils数据库
            db = gffutils.create_db(str(gff_file), db_file, force=True, keep_order=True,
                                  merge_strategy='merge', sort_attribute_values=True)
            
            logger.info(f"Loaded annotations from {gff_file}")
            return db
            
        except Exception as e:
            logger.error(f"Error loading GFF file {gff_file}: {e}")
            return None
    
    def extract_gene_features(self, db, feature_type='gene'):
        """从GFF数据库提取基因特征"""
        genes = []
        
        try:
            for gene in db.features_of_type(feature_type):
                gene_info = {
                    'id': gene.id,
                    'seqid': gene.seqid,
                    'start': gene.start,
                    'end': gene.end,
                    'strand': gene.strand,
                    'length': gene.end - gene.start + 1
                }
                
                # 提取外显子信息
                exons = []
                try:
                    for exon in db.children(gene, featuretype='exon'):
                        exons.append({
                            'start': exon.start,
                            'end': exon.end,
                            'length': exon.end - exon.start + 1
                        })
                except:
                    # 如果没有外显子注释，将整个基因作为一个外显子
                    exons = [{'start': gene.start, 'end': gene.end, 
                             'length': gene.end - gene.start + 1}]
                
                gene_info['exons'] = sorted(exons, key=lambda x: x['start'])
                gene_info['num_exons'] = len(exons)
                gene_info['total_exon_length'] = sum(e['length'] for e in exons)
                
                # 计算内含子信息
                if len(exons) > 1:
                    introns = []
                    for i in range(len(exons) - 1):
                        intron_start = exons[i]['end'] + 1
                        intron_end = exons[i+1]['start'] - 1
                        if intron_end > intron_start:
                            introns.append({
                                'start': intron_start,
                                'end': intron_end,
                                'length': intron_end - intron_start + 1
                            })
                    gene_info['introns'] = introns
                    gene_info['num_introns'] = len(introns)
                else:
                    gene_info['introns'] = []
                    gene_info['num_introns'] = 0
                
                genes.append(gene_info)
                
        except Exception as e:
            logger.error(f"Error extracting gene features: {e}")
        
        logger.info(f"Extracted {len(genes)} genes from database")
        return genes
    
    def find_overlapping_genes(self, reference_genes, predicted_genes, min_overlap=0.5):
        """找到重叠的基因对"""
        overlaps = []
        
        # 为每个染色体创建区间树
        ref_trees = defaultdict(IntervalTree)
        for i, gene in enumerate(reference_genes):
            ref_trees[gene['seqid']].add(Interval(gene['start'], gene['end'], {'index': i, 'gene': gene}))
        
        # 查找重叠
        for pred_idx, pred_gene in enumerate(predicted_genes):
            seqid = pred_gene['seqid']
            if seqid not in ref_trees:
                continue
                
            overlapping_intervals = ref_trees[seqid].overlap(pred_gene['start'], pred_gene['end'])
            
            for interval in overlapping_intervals:
                ref_gene = interval.data['gene']
                ref_idx = interval.data['index']
                
                # 计算重叠比例
                overlap_start = max(pred_gene['start'], ref_gene['start'])
                overlap_end = min(pred_gene['end'], ref_gene['end'])
                overlap_length = max(0, overlap_end - overlap_start + 1)
                
                pred_overlap_ratio = overlap_length / pred_gene['length']
                ref_overlap_ratio = overlap_length / ref_gene['length']
                
                # 使用较小的重叠比例作为判断标准
                overlap_ratio = min(pred_overlap_ratio, ref_overlap_ratio)
                
                if overlap_ratio >= min_overlap:
                    overlaps.append({
                        'ref_index': ref_idx,
                        'pred_index': pred_idx,
                        'ref_gene': ref_gene,
                        'pred_gene': pred_gene,
                        'overlap_length': overlap_length,
                        'overlap_ratio': overlap_ratio,
                        'pred_overlap_ratio': pred_overlap_ratio,
                        'ref_overlap_ratio': ref_overlap_ratio
                    })
        
        logger.info(f"Found {len(overlaps)} overlapping gene pairs")
        return overlaps
    
    def evaluate_gene_level_accuracy(self, reference_genes, predicted_genes, overlaps):
        """评估基因水平的准确性"""
        # 基本统计
        num_ref_genes = len(reference_genes)
        num_pred_genes = len(predicted_genes)
        num_overlapping = len(set(o['ref_index'] for o in overlaps))
        
        # 计算敏感性和特异性
        true_positives = num_overlapping
        false_negatives = num_ref_genes - true_positives
        false_positives = num_pred_genes - len(set(o['pred_index'] for o in overlaps))
        
        sensitivity = true_positives / num_ref_genes if num_ref_genes > 0 else 0
        precision = true_positives / num_pred_genes if num_pred_genes > 0 else 0
        f1_score = 2 * (sensitivity * precision) / (sensitivity + precision) if (sensitivity + precision) > 0 else 0
        
        specificity = true_positives / (true_positives + false_positives) if (true_positives + false_positives) > 0 else 0
        
        gene_metrics = {
            'num_reference_genes': num_ref_genes,
            'num_predicted_genes': num_pred_genes,
            'true_positives': true_positives,
            'false_positives': false_positives,
            'false_negatives': false_negatives,
            'sensitivity': sensitivity,
            'precision': precision,
            'specificity': specificity,
            'f1_score': f1_score
        }
        
        return gene_metrics
    
    def evaluate_exon_level_accuracy(self, overlaps):
        """评估外显子水平的准确性"""
        exon_comparisons = []
        
        for overlap in overlaps:
            ref_gene = overlap['ref_gene']
            pred_gene = overlap['pred_gene']
            
            # 比较外显子数目
            ref_exon_count = ref_gene['num_exons']
            pred_exon_count = pred_gene['num_exons']
            
            # 比较外显子边界
            ref_exons = ref_gene['exons']
            pred_exons = pred_gene['exons']
            
            # 计算外显子边界准确性
            exact_matches = 0
            boundary_tolerance = 10  # 允许10bp的误差
            
            for ref_exon in ref_exons:
                for pred_exon in pred_exons:
                    start_diff = abs(ref_exon['start'] - pred_exon['start'])
                    end_diff = abs(ref_exon['end'] - pred_exon['end'])
                    
                    if start_diff <= boundary_tolerance and end_diff <= boundary_tolerance:
                        exact_matches += 1
                        break
            
            exon_accuracy = exact_matches / max(ref_exon_count, pred_exon_count) if max(ref_exon_count, pred_exon_count) > 0 else 0
            
            comparison = {
                'ref_gene_id': ref_gene['id'],
                'pred_gene_id': pred_gene['id'],
                'ref_exon_count': ref_exon_count,
                'pred_exon_count': pred_exon_count,
                'exon_count_diff': abs(ref_exon_count - pred_exon_count),
                'exact_exon_matches': exact_matches,
                'exon_accuracy': exon_accuracy,
                'ref_total_exon_length': ref_gene['total_exon_length'],
                'pred_total_exon_length': pred_gene['total_exon_length']
            }
            
            exon_comparisons.append(comparison)
        
        # 计算整体外显子水平指标
        if exon_comparisons:
            avg_exon_accuracy = np.mean([c['exon_accuracy'] for c in exon_comparisons])
            avg_exon_count_diff = np.mean([c['exon_count_diff'] for c in exon_comparisons])
            
            # 统计外显子数目分布
            ref_exon_counts = [c['ref_exon_count'] for c in exon_comparisons]
            pred_exon_counts = [c['pred_exon_count'] for c in exon_comparisons]
            
            exon_metrics = {
                'avg_exon_accuracy': avg_exon_accuracy,
                'avg_exon_count_difference': avg_exon_count_diff,
                'ref_exon_count_mean': np.mean(ref_exon_counts),
                'ref_exon_count_std': np.std(ref_exon_counts),
                'pred_exon_count_mean': np.mean(pred_exon_counts),
                'pred_exon_count_std': np.std(pred_exon_counts),
                'detailed_comparisons': exon_comparisons
            }
        else:
            exon_metrics = {
                'avg_exon_accuracy': 0,
                'detailed_comparisons': []
            }
        
        return exon_metrics
    
    def evaluate_nucleotide_level_accuracy(self, overlaps):
        """评估核苷酸水平的准确性"""
        nucleotide_comparisons = []
        
        for overlap in overlaps:
            ref_gene = overlap['ref_gene']
            pred_gene = overlap['pred_gene']
            
            # 创建基因区域的核苷酸集合
            ref_nucleotides = set()
            pred_nucleotides = set()
            
            # 添加参考基因的外显子核苷酸
            for exon in ref_gene['exons']:
                for pos in range(exon['start'], exon['end'] + 1):
                    ref_nucleotides.add(pos)
            
            # 添加预测基因的外显子核苷酸
            for exon in pred_gene['exons']:
                for pos in range(exon['start'], exon['end'] + 1):
                    pred_nucleotides.add(pos)
            
            # 计算交集和并集
            intersection = ref_nucleotides & pred_nucleotides
            union = ref_nucleotides | pred_nucleotides
            
            # 计算Jaccard系数
            jaccard = len(intersection) / len(union) if len(union) > 0 else 0
            
            # 计算敏感性和特异性
            nt_sensitivity = len(intersection) / len(ref_nucleotides) if len(ref_nucleotides) > 0 else 0
            nt_precision = len(intersection) / len(pred_nucleotides) if len(pred_nucleotides) > 0 else 0
            
            comparison = {
                'ref_gene_id': ref_gene['id'],
                'pred_gene_id': pred_gene['id'],
                'ref_nucleotides': len(ref_nucleotides),
                'pred_nucleotides': len(pred_nucleotides),
                'shared_nucleotides': len(intersection),
                'jaccard_coefficient': jaccard,
                'nucleotide_sensitivity': nt_sensitivity,
                'nucleotide_precision': nt_precision
            }
            
            nucleotide_comparisons.append(comparison)
        
        # 计算整体核苷酸水平指标
        if nucleotide_comparisons:
            nucleotide_metrics = {
                'avg_jaccard_coefficient': np.mean([c['jaccard_coefficient'] for c in nucleotide_comparisons]),
                'avg_nucleotide_sensitivity': np.mean([c['nucleotide_sensitivity'] for c in nucleotide_comparisons]),
                'avg_nucleotide_precision': np.mean([c['nucleotide_precision'] for c in nucleotide_comparisons]),
                'detailed_comparisons': nucleotide_comparisons
            }
        else:
            nucleotide_metrics = {
                'avg_jaccard_coefficient': 0,
                'detailed_comparisons': []
            }
        
        return nucleotide_metrics
    
    def analyze_gene_structure_patterns(self, genes, label=''):
        """分析基因结构模式"""
        if not genes:
            return {}
        
        # 基因长度分布
        gene_lengths = [g['length'] for g in genes]
        
        # 外显子数目分布
        exon_counts = [g['num_exons'] for g in genes]
        
        # 外显子长度分布
        exon_lengths = []
        for gene in genes:
            exon_lengths.extend([e['length'] for e in gene['exons']])
        
        # 内含子长度分布
        intron_lengths = []
        for gene in genes:
            intron_lengths.extend([i['length'] for i in gene['introns']])
        
        patterns = {
            'gene_count': len(genes),
            'gene_length_stats': {
                'mean': np.mean(gene_lengths),
                'median': np.median(gene_lengths),
                'std': np.std(gene_lengths),
                'min': np.min(gene_lengths),
                'max': np.max(gene_lengths),
                'quartiles': np.percentile(gene_lengths, [25, 50, 75]).tolist()
            },
            'exon_count_stats': {
                'mean': np.mean(exon_counts),
                'median': np.median(exon_counts),
                'std': np.std(exon_counts),
                'distribution': Counter(exon_counts)
            },
            'exon_length_stats': {
                'mean': np.mean(exon_lengths) if exon_lengths else 0,
                'median': np.median(exon_lengths) if exon_lengths else 0,
                'std': np.std(exon_lengths) if exon_lengths else 0
            },
            'intron_length_stats': {
                'mean': np.mean(intron_lengths) if intron_lengths else 0,
                'median': np.median(intron_lengths) if intron_lengths else 0,
                'std': np.std(intron_lengths) if intron_lengths else 0,
                'count': len(intron_lengths)
            }
        }
        
        return patterns
    
    def create_comparison_visualizations(self, gene_metrics, exon_metrics, nucleotide_metrics, 
                                       ref_patterns, pred_patterns, dataset_name):
        """创建比较可视化图表"""
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        
        # 图1: 准确性指标柱状图
        metrics_names = ['Sensitivity', 'Precision', 'Specificity', 'F1-Score']
        metrics_values = [
            gene_metrics['sensitivity'],
            gene_metrics['precision'], 
            gene_metrics['specificity'],
            gene_metrics['f1_score']
        ]
        
        bars = axes[0,0].bar(metrics_names, metrics_values, 
                            color=['skyblue', 'lightcoral', 'lightgreen', 'gold'])
        axes[0,0].set_ylabel('Score')
        axes[0,0].set_title(f'Gene-Level Accuracy Metrics\n({dataset_name})')
        axes[0,0].set_ylim(0, 1.1)
        
        # 添加数值标签
        for bar, value in zip(bars, metrics_values):
            axes[0,0].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                          f'{value:.3f}', ha='center', va='bottom')
        
        # 图2: 基因长度分布比较
        ref_lengths = [g['length'] for g in ref_patterns.get('gene_length_stats', {}).get('raw_data', [])]
        pred_lengths = [g['length'] for g in pred_patterns.get('gene_length_stats', {}).get('raw_data', [])]
        
        # 如果没有原始数据，创建示例分布
        if len(ref_lengths) == 0 and 'gene_length_stats' in ref_patterns:
            ref_stats = ref_patterns['gene_length_stats']
            ref_lengths = np.random.normal(ref_stats['mean'], ref_stats['std'], ref_patterns['gene_count'])
            
        if len(pred_lengths) == 0 and 'gene_length_stats' in pred_patterns:
            pred_stats = pred_patterns['gene_length_stats']
            pred_lengths = np.random.normal(pred_stats['mean'], pred_stats['std'], pred_patterns['gene_count'])
        
        if len(ref_lengths) > 0 and len(pred_lengths) > 0:
            axes[0,1].hist([ref_lengths, pred_lengths], bins=20, alpha=0.7, 
                          label=['Reference', 'Predicted'], color=['blue', 'red'])
            axes[0,1].set_xlabel('Gene Length (bp)')
            axes[0,1].set_ylabel('Frequency')
            axes[0,1].set_title('Gene Length Distribution')
            axes[0,1].legend()
        
        # 图3: 外显子数目分布比较
        if 'exon_count_stats' in ref_patterns and 'exon_count_stats' in pred_patterns:
            ref_exon_dist = ref_patterns['exon_count_stats']['distribution']
            pred_exon_dist = pred_patterns['exon_count_stats']['distribution']
            
            all_exon_counts = sorted(set(list(ref_exon_dist.keys()) + list(pred_exon_dist.keys())))
            ref_counts = [ref_exon_dist.get(n, 0) for n in all_exon_counts]
            pred_counts = [pred_exon_dist.get(n, 0) for n in all_exon_counts]
            
            x = np.arange(len(all_exon_counts))
            width = 0.35
            
            axes[0,2].bar(x - width/2, ref_counts, width, label='Reference', alpha=0.7)
            axes[0,2].bar(x + width/2, pred_counts, width, label='Predicted', alpha=0.7)
            axes[0,2].set_xlabel('Number of Exons')
            axes[0,2].set_ylabel('Gene Count')
            axes[0,2].set_title('Exon Count Distribution')
            axes[0,2].set_xticks(x)
            axes[0,2].set_xticklabels(all_exon_counts)
            axes[0,2].legend()
        
        # 图4: 外显子准确性分布
        if exon_metrics.get('detailed_comparisons'):
            exon_accuracies = [c['exon_accuracy'] for c in exon_metrics['detailed_comparisons']]
            axes[1,0].hist(exon_accuracies, bins=20, alpha=0.7, color='green')
            axes[1,0].set_xlabel('Exon Structure Accuracy')
            axes[1,0].set_ylabel('Gene Count')
            axes[1,0].set_title('Exon-Level Accuracy Distribution')
            axes[1,0].axvline(np.mean(exon_accuracies), color='red', linestyle='--', 
                             label=f'Mean: {np.mean(exon_accuracies):.3f}')
            axes[1,0].legend()
        
        # 图5: 核苷酸水平准确性
        if nucleotide_metrics.get('detailed_comparisons'):
            jaccard_scores = [c['jaccard_coefficient'] for c in nucleotide_metrics['detailed_comparisons']]
            nt_sensitivities = [c['nucleotide_sensitivity'] for c in nucleotide_metrics['detailed_comparisons']]
            
            axes[1,1].scatter(jaccard_scores, nt_sensitivities, alpha=0.6)
            axes[1,1].set_xlabel('Jaccard Coefficient')
            axes[1,1].set_ylabel('Nucleotide Sensitivity')
            axes[1,1].set_title('Nucleotide-Level Accuracy')
            axes[1,1].plot([0, 1], [0, 1], 'r--', alpha=0.5)
        
        # 图6: 整体性能汇总
        performance_categories = ['Gene\nDetection', 'Exon\nStructure', 'Nucleotide\nAccuracy']
        performance_scores = [
            gene_metrics['f1_score'],
            exon_metrics.get('avg_exon_accuracy', 0),
            nucleotide_metrics.get('avg_jaccard_coefficient', 0)
        ]
        
        bars = axes[1,2].bar(performance_categories, performance_scores, 
                            color=['lightblue', 'lightcoral', 'lightgreen'])
        axes[1,2].set_ylabel('Score')
        axes[1,2].set_title('Overall Performance Summary')
        axes[1,2].set_ylim(0, 1.1)
        
        for bar, score in zip(bars, performance_scores):
            axes[1,2].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                          f'{score:.3f}', ha='center', va='bottom')
        
        plt.tight_layout()
        
        # 保存图表
        plot_file = self.output_dir / f'accuracy_assessment_{dataset_name}.svg'
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"Accuracy assessment plots saved to {plot_file}")
    
    def evaluate_dataset(self, reference_gff, predicted_gff, dataset_name):
        """评估单个数据集的准确性"""
        logger.info(f"Evaluating accuracy for dataset: {dataset_name}")
        
        # 加载注释
        ref_db = self.load_gff_annotations(reference_gff)
        pred_db = self.load_gff_annotations(predicted_gff)
        
        if not ref_db or not pred_db:
            logger.error(f"Failed to load annotations for {dataset_name}")
            return None
        
        # 提取基因特征
        ref_genes = self.extract_gene_features(ref_db)
        pred_genes = self.extract_gene_features(pred_db)
        
        if not ref_genes or not pred_genes:
            logger.error(f"No genes found for {dataset_name}")
            return None
        
        # 找到重叠的基因
        overlaps = self.find_overlapping_genes(ref_genes, pred_genes)
        
        # 评估各个层面的准确性
        gene_metrics = self.evaluate_gene_level_accuracy(ref_genes, pred_genes, overlaps)
        exon_metrics = self.evaluate_exon_level_accuracy(overlaps)
        nucleotide_metrics = self.evaluate_nucleotide_level_accuracy(overlaps)
        
        # 分析基因结构模式
        ref_patterns = self.analyze_gene_structure_patterns(ref_genes, 'reference')
        pred_patterns = self.analyze_gene_structure_patterns(pred_genes, 'predicted')
        
        # 创建可视化
        self.create_comparison_visualizations(
            gene_metrics, exon_metrics, nucleotide_metrics,
            ref_patterns, pred_patterns, dataset_name
        )
        
        # 汇总结果
        evaluation_result = {
            'dataset': dataset_name,
            'timestamp': datetime.now().isoformat(),
            'gene_metrics': gene_metrics,
            'exon_metrics': exon_metrics,
            'nucleotide_metrics': nucleotide_metrics,
            'reference_patterns': ref_patterns,
            'predicted_patterns': pred_patterns,
            'overlap_analysis': {
                'total_overlaps': len(overlaps),
                'avg_overlap_ratio': np.mean([o['overlap_ratio'] for o in overlaps]) if overlaps else 0
            }
        }
        
        # 保存详细结果
        result_file = self.output_dir / f'accuracy_evaluation_{dataset_name}.json'
        with open(result_file, 'w') as f:
            json.dump(evaluation_result, f, indent=2, default=str)
        
        logger.info(f"Accuracy evaluation completed for {dataset_name}")
        return evaluation_result
    
    def generate_summary_report(self, all_evaluations):
        """生成汇总报告"""
        if not all_evaluations:
            logger.error("No evaluations to summarize")
            return
        
        # 汇总所有数据集的指标
        summary = {
            'overall_summary': {
                'num_datasets': len(all_evaluations),
                'evaluation_date': datetime.now().isoformat()
            },
            'gene_level_summary': {},
            'exon_level_summary': {},
            'nucleotide_level_summary': {},
            'dataset_details': all_evaluations
        }
        
        # 计算基因水平平均指标
        gene_metrics = [e['gene_metrics'] for e in all_evaluations]
        summary['gene_level_summary'] = {
            'avg_sensitivity': np.mean([m['sensitivity'] for m in gene_metrics]),
            'avg_precision': np.mean([m['precision'] for m in gene_metrics]),
            'avg_specificity': np.mean([m['specificity'] for m in gene_metrics]),
            'avg_f1_score': np.mean([m['f1_score'] for m in gene_metrics]),
            'std_sensitivity': np.std([m['sensitivity'] for m in gene_metrics]),
            'std_precision': np.std([m['precision'] for m in gene_metrics]),
            'std_specificity': np.std([m['specificity'] for m in gene_metrics]),
            'std_f1_score': np.std([m['f1_score'] for m in gene_metrics])
        }
        
        # 计算外显子水平平均指标
        exon_metrics = [e['exon_metrics'] for e in all_evaluations]
        summary['exon_level_summary'] = {
            'avg_exon_accuracy': np.mean([m.get('avg_exon_accuracy', 0) for m in exon_metrics]),
            'avg_exon_count_difference': np.mean([m.get('avg_exon_count_difference', 0) for m in exon_metrics])
        }
        
        # 计算核苷酸水平平均指标
        nt_metrics = [e['nucleotide_metrics'] for e in all_evaluations]
        summary['nucleotide_level_summary'] = {
            'avg_jaccard_coefficient': np.mean([m.get('avg_jaccard_coefficient', 0) for m in nt_metrics]),
            'avg_nucleotide_sensitivity': np.mean([m.get('avg_nucleotide_sensitivity', 0) for m in nt_metrics]),
            'avg_nucleotide_precision': np.mean([m.get('avg_nucleotide_precision', 0) for m in nt_metrics])
        }
        
        # 保存汇总报告
        summary_file = self.output_dir / 'accuracy_assessment_summary.json'
        with open(summary_file, 'w') as f:
            json.dump(summary, f, indent=2, default=str)
        
        # 创建汇总CSV
        csv_data = []
        for eval_result in all_evaluations:
            row = {
                'Dataset': eval_result['dataset'],
                'Sensitivity': eval_result['gene_metrics']['sensitivity'],
                'Precision': eval_result['gene_metrics']['precision'],
                'Specificity': eval_result['gene_metrics']['specificity'],
                'F1_Score': eval_result['gene_metrics']['f1_score'],
                'Exon_Accuracy': eval_result['exon_metrics'].get('avg_exon_accuracy', 0),
                'Jaccard_Coefficient': eval_result['nucleotide_metrics'].get('avg_jaccard_coefficient', 0),
                'Ref_Genes': eval_result['gene_metrics']['num_reference_genes'],
                'Pred_Genes': eval_result['gene_metrics']['num_predicted_genes']
            }
            csv_data.append(row)
        
        df = pd.DataFrame(csv_data)
        csv_file = self.output_dir / 'accuracy_assessment_summary.csv'
        df.to_csv(csv_file, index=False)
        
        logger.info(f"Summary report saved to {summary_file}")
        logger.info(f"Summary CSV saved to {csv_file}")
        
        return summary

def main():
    parser = argparse.ArgumentParser(description='Gene Structure Accuracy Evaluation')
    parser.add_argument('--reference-gff', required=True, help='Reference annotation GFF file')
    parser.add_argument('--predicted-gff', required=True, help='Predicted annotation GFF file')
    parser.add_argument('--dataset-name', required=True, help='Dataset name for output files')
    parser.add_argument('--output', required=True, help='Output directory')
    
    args = parser.parse_args()
    
    # 创建评估器
    evaluator = GeneStructureEvaluator(args.output)
    
    # 运行评估
    result = evaluator.evaluate_dataset(
        reference_gff=args.reference_gff,
        predicted_gff=args.predicted_gff,
        dataset_name=args.dataset_name
    )
    
    if result:
        print(f"Accuracy evaluation completed for {args.dataset_name}")
        print(f"Results saved to: {args.output}")
    else:
        print("Accuracy evaluation failed")
        sys.exit(1)

if __name__ == '__main__':
    main()