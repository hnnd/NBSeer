#!/usr/bin/env python3
"""
基因特征统计分析脚本
分析NBS基因的结构特征分布和模式
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
from collections import Counter, defaultdict
from scipy import stats
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
import logging
from datetime import datetime
import gffutils
from Bio import SeqIO

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class GeneFeatureAnalyzer:
    def __init__(self, output_dir):
        """
        初始化基因特征分析器
        
        Args:
            output_dir: 结果输出目录
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # 设置图表样式
        plt.style.use('default')
        sns.set_palette("husl")
        
        # 统计结果存储
        self.feature_stats = {}
        self.comparative_stats = {}
    
    def _json_serializer(self, obj):
        """JSON序列化自定义函数，处理numpy数据类型"""
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        return obj
    
    def load_genome_sequences(self, genome_file):
        """加载基因组序列用于GC含量分析"""
        try:
            sequences = {}
            for record in SeqIO.parse(genome_file, "fasta"):
                sequences[record.id.split()[0]] = str(record.seq)  # 只取ID的第一部分
            logger.info(f"Loaded {len(sequences)} chromosomes from {genome_file}")
            return sequences
        except Exception as e:
            logger.error(f"Error loading genome sequences: {e}")
            return {}
    
    def extract_gene_features_from_gff(self, gff_file, genome_sequences=None):
        """从GFF文件提取基因特征"""
        try:
            # 创建数据库
            db_file = str(gff_file) + '.db'
            if Path(db_file).exists():
                os.remove(db_file)
            
            db = gffutils.create_db(str(gff_file), db_file, force=True, keep_order=True,
                                  merge_strategy='create_unique', sort_attribute_values=True,
                                  id_spec={'gene': 'ID', 'mRNA': 'ID', 'exon': 'ID', 'CDS': 'ID'})
            
            genes_features = []
            
            for gene in db.features_of_type('gene'):
                # 基本信息
                gene_id = gene.id
                seqid = gene.seqid
                start = gene.start
                end = gene.end
                strand = gene.strand
                gene_length = end - start + 1
                
                # 提取外显子
                exons = []
                try:
                    for exon in db.children(gene, featuretype='exon'):
                        exons.append({
                            'start': exon.start,
                            'end': exon.end,
                            'length': exon.end - exon.start + 1
                        })
                except:
                    # 如果没有外显子，将整个基因作为单外显子
                    exons = [{'start': start, 'end': end, 'length': gene_length}]
                
                exons = sorted(exons, key=lambda x: x['start'])
                
                # 计算外显子统计
                num_exons = len(exons)
                exon_lengths = [e['length'] for e in exons]
                total_exon_length = sum(exon_lengths)
                avg_exon_length = total_exon_length / num_exons if num_exons > 0 else 0
                
                # 计算内含子
                introns = []
                if num_exons > 1:
                    for i in range(num_exons - 1):
                        intron_start = exons[i]['end'] + 1
                        intron_end = exons[i+1]['start'] - 1
                        if intron_end >= intron_start:
                            introns.append({
                                'start': intron_start,
                                'end': intron_end,
                                'length': intron_end - intron_start + 1
                            })
                
                num_introns = len(introns)
                intron_lengths = [i['length'] for i in introns]
                total_intron_length = sum(intron_lengths)
                avg_intron_length = total_intron_length / num_introns if num_introns > 0 else 0
                
                # 计算GC含量（如果有基因组序列）
                gc_content = None
                if genome_sequences and seqid in genome_sequences:
                    try:
                        sequence = genome_sequences[seqid][start-1:end]  # 转换为0-based
                        gc_count = sequence.upper().count('G') + sequence.upper().count('C')
                        gc_content = gc_count / len(sequence) if len(sequence) > 0 else 0
                    except:
                        gc_content = None
                
                # 基因密度相关（需要更多上下文信息，这里简化处理）
                gene_density = None  # 可以后续计算
                
                feature = {
                    'gene_id': gene_id,
                    'seqid': seqid,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'gene_length': gene_length,
                    'num_exons': num_exons,
                    'num_introns': num_introns,
                    'total_exon_length': total_exon_length,
                    'total_intron_length': total_intron_length,
                    'avg_exon_length': avg_exon_length,
                    'avg_intron_length': avg_intron_length,
                    'exon_lengths': exon_lengths,
                    'intron_lengths': intron_lengths,
                    'gc_content': gc_content,
                    'coding_ratio': total_exon_length / gene_length if gene_length > 0 else 0
                }
                
                genes_features.append(feature)
            
            logger.info(f"Extracted features from {len(genes_features)} genes")
            return genes_features
            
        except Exception as e:
            logger.error(f"Error extracting features from {gff_file}: {e}")
            return []
    
    def calculate_basic_statistics(self, features, dataset_name):
        """计算基本统计量"""
        if not features:
            return {}
        
        # 转换为DataFrame便于分析
        df = pd.DataFrame(features)
        
        # 数值特征
        numeric_features = [
            'gene_length', 'num_exons', 'num_introns', 
            'total_exon_length', 'total_intron_length',
            'avg_exon_length', 'avg_intron_length', 
            'gc_content', 'coding_ratio'
        ]
        
        stats_dict = {
            'dataset': dataset_name,
            'total_genes': len(features),
            'features': {}
        }
        
        for feature in numeric_features:
            if feature in df.columns:
                values = df[feature].dropna()
                if len(values) > 0:
                    stats_dict['features'][feature] = {
                        'count': len(values),
                        'mean': float(values.mean()),
                        'median': float(values.median()),
                        'std': float(values.std()),
                        'min': float(values.min()),
                        'max': float(values.max()),
                        'q25': float(values.quantile(0.25)),
                        'q75': float(values.quantile(0.75)),
                        'skewness': float(values.skew()),
                        'kurtosis': float(values.kurtosis())
                    }
        
        # 分类特征分布
        stats_dict['distributions'] = {
            'exon_count_distribution': dict(df['num_exons'].value_counts().sort_index()),
            'intron_count_distribution': dict(df['num_introns'].value_counts().sort_index()),
            'chromosome_distribution': dict(df['seqid'].value_counts())
        }
        
        return stats_dict
    
    def analyze_gene_length_distribution(self, all_features, output_prefix):
        """分析基因长度分布"""
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        colors = ['skyblue', 'lightcoral', 'lightgreen', 'gold', 'lightpink']
        
        # 准备数据
        datasets = list(all_features.keys())
        gene_lengths_by_dataset = {}
        
        for dataset in datasets:
            features = all_features[dataset]
            gene_lengths_by_dataset[dataset] = [f['gene_length'] for f in features]
        
        # 子图1: 基因长度分布直方图
        for i, (dataset, lengths) in enumerate(gene_lengths_by_dataset.items()):
            axes[0,0].hist(lengths, bins=50, alpha=0.7, label=dataset, 
                          color=colors[i % len(colors)], density=True)
        
        axes[0,0].set_xlabel('Gene Length (bp)')
        axes[0,0].set_ylabel('Density')
        axes[0,0].set_title('Gene Length Distribution')
        axes[0,0].legend()
        axes[0,0].set_xlim(0, np.percentile([l for lengths in gene_lengths_by_dataset.values() 
                                            for l in lengths], 95))
        
        # 子图2: 基因长度对数分布
        for i, (dataset, lengths) in enumerate(gene_lengths_by_dataset.items()):
            log_lengths = np.log10([l for l in lengths if l > 0])
            axes[0,1].hist(log_lengths, bins=30, alpha=0.7, label=dataset,
                          color=colors[i % len(colors)], density=True)
        
        axes[0,1].set_xlabel('Log10(Gene Length)')
        axes[0,1].set_ylabel('Density')
        axes[0,1].set_title('Gene Length Distribution (Log Scale)')
        axes[0,1].legend()
        
        # 子图3: 箱线图比较
        data_for_box = []
        labels_for_box = []
        for dataset, lengths in gene_lengths_by_dataset.items():
            data_for_box.append(lengths)
            labels_for_box.append(dataset)
        
        bp = axes[1,0].boxplot(data_for_box, labels=labels_for_box, patch_artist=True)
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
        
        axes[1,0].set_ylabel('Gene Length (bp)')
        axes[1,0].set_title('Gene Length Comparison (Boxplot)')
        axes[1,0].tick_params(axis='x', rotation=45)
        
        # 子图4: 累积分布函数
        for i, (dataset, lengths) in enumerate(gene_lengths_by_dataset.items()):
            sorted_lengths = np.sort(lengths)
            y_vals = np.arange(1, len(sorted_lengths) + 1) / len(sorted_lengths)
            axes[1,1].plot(sorted_lengths, y_vals, label=dataset, 
                          color=colors[i % len(colors)], linewidth=2)
        
        axes[1,1].set_xlabel('Gene Length (bp)')
        axes[1,1].set_ylabel('Cumulative Probability')
        axes[1,1].set_title('Gene Length Cumulative Distribution')
        axes[1,1].legend()
        axes[1,1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        # 保存图表
        plot_file = self.output_dir / f'{output_prefix}_gene_length_analysis.svg'
        plt.savefig(plot_file, format='svg', bbox_inches='tight')
        plt.close()
        
        logger.info(f"Gene length analysis plot saved to {plot_file}")
    
    def analyze_exon_intron_structure(self, all_features, output_prefix):
        """分析外显子-内含子结构"""
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        
        colors = ['skyblue', 'lightcoral', 'lightgreen', 'gold', 'lightpink']
        datasets = list(all_features.keys())
        
        # 准备数据
        exon_counts_by_dataset = {}
        exon_lengths_by_dataset = {}
        intron_lengths_by_dataset = {}
        
        for dataset in datasets:
            features = all_features[dataset]
            exon_counts_by_dataset[dataset] = [f['num_exons'] for f in features]
            
            # 收集所有外显子长度
            all_exon_lengths = []
            all_intron_lengths = []
            for f in features:
                all_exon_lengths.extend(f['exon_lengths'])
                all_intron_lengths.extend(f['intron_lengths'])
            
            exon_lengths_by_dataset[dataset] = all_exon_lengths
            intron_lengths_by_dataset[dataset] = all_intron_lengths
        
        # 子图1: 外显子数目分布
        all_exon_counts = sorted(set([c for counts in exon_counts_by_dataset.values() 
                                     for c in counts]))
        x_pos = np.arange(len(all_exon_counts))
        width = 0.8 / len(datasets)
        
        for i, (dataset, counts) in enumerate(exon_counts_by_dataset.items()):
            count_dist = Counter(counts)
            y_values = [count_dist.get(count, 0) for count in all_exon_counts]
            axes[0,0].bar(x_pos + i*width, y_values, width, 
                         label=dataset, color=colors[i % len(colors)], alpha=0.7)
        
        axes[0,0].set_xlabel('Number of Exons')
        axes[0,0].set_ylabel('Gene Count')
        axes[0,0].set_title('Exon Count Distribution')
        axes[0,0].set_xticks(x_pos + width * (len(datasets)-1) / 2)
        axes[0,0].set_xticklabels(all_exon_counts)
        axes[0,0].legend()
        
        # 子图2: 外显子长度分布
        for i, (dataset, lengths) in enumerate(exon_lengths_by_dataset.items()):
            if lengths:
                axes[0,1].hist(lengths, bins=50, alpha=0.7, label=dataset,
                              color=colors[i % len(colors)], density=True)
        
        axes[0,1].set_xlabel('Exon Length (bp)')
        axes[0,1].set_ylabel('Density')
        axes[0,1].set_title('Exon Length Distribution')
        axes[0,1].legend()
        axes[0,1].set_xlim(0, 2000)  # 限制显示范围
        
        # 子图3: 内含子长度分布
        for i, (dataset, lengths) in enumerate(intron_lengths_by_dataset.items()):
            if lengths:
                # 使用对数尺度
                log_lengths = np.log10([l for l in lengths if l > 0])
                axes[0,2].hist(log_lengths, bins=50, alpha=0.7, label=dataset,
                              color=colors[i % len(colors)], density=True)
        
        axes[0,2].set_xlabel('Log10(Intron Length)')
        axes[0,2].set_ylabel('Density')
        axes[0,2].set_title('Intron Length Distribution (Log Scale)')
        axes[0,2].legend()
        
        # 子图4: 外显子长度箱线图
        exon_data_for_box = []
        exon_labels = []
        for dataset, lengths in exon_lengths_by_dataset.items():
            if lengths:
                exon_data_for_box.append(lengths)
                exon_labels.append(dataset)
        
        if exon_data_for_box:
            bp1 = axes[1,0].boxplot(exon_data_for_box, labels=exon_labels, patch_artist=True)
            for patch, color in zip(bp1['boxes'], colors):
                patch.set_facecolor(color)
        
        axes[1,0].set_ylabel('Exon Length (bp)')
        axes[1,0].set_title('Exon Length Comparison')
        axes[1,0].tick_params(axis='x', rotation=45)
        axes[1,0].set_ylim(0, 2000)  # 限制Y轴范围
        
        # 子图5: 内含子长度箱线图
        intron_data_for_box = []
        intron_labels = []
        for dataset, lengths in intron_lengths_by_dataset.items():
            if lengths:
                intron_data_for_box.append(lengths)
                intron_labels.append(dataset)
        
        if intron_data_for_box:
            bp2 = axes[1,1].boxplot(intron_data_for_box, labels=intron_labels, patch_artist=True)
            for patch, color in zip(bp2['boxes'], colors):
                patch.set_facecolor(color)
        
        axes[1,1].set_ylabel('Intron Length (bp)')
        axes[1,1].set_title('Intron Length Comparison')
        axes[1,1].tick_params(axis='x', rotation=45)
        axes[1,1].set_yscale('log')  # 使用对数刻度
        
        # 子图6: 编码区比例分布
        for i, dataset in enumerate(datasets):
            features = all_features[dataset]
            coding_ratios = [f['coding_ratio'] for f in features if f['coding_ratio'] is not None]
            if coding_ratios:
                axes[1,2].hist(coding_ratios, bins=30, alpha=0.7, label=dataset,
                              color=colors[i % len(colors)], density=True)
        
        axes[1,2].set_xlabel('Coding Ratio (Exon Length / Gene Length)')
        axes[1,2].set_ylabel('Density')
        axes[1,2].set_title('Coding Ratio Distribution')
        axes[1,2].legend()
        
        plt.tight_layout()
        
        # 保存图表
        plot_file = self.output_dir / f'{output_prefix}_exon_intron_analysis.svg'
        plt.savefig(plot_file, format='svg', bbox_inches='tight')
        plt.close()
        
        logger.info(f"Exon-intron analysis plot saved to {plot_file}")
    
    def analyze_gc_content(self, all_features, output_prefix):
        """分析GC含量分布"""
        fig, axes = plt.subplots(1, 3, figsize=(18, 6))
        
        colors = ['skyblue', 'lightcoral', 'lightgreen', 'gold', 'lightpink']
        datasets = list(all_features.keys())
        
        # 收集GC含量数据
        gc_data_by_dataset = {}
        for dataset in datasets:
            features = all_features[dataset]
            gc_contents = [f['gc_content'] for f in features 
                          if f['gc_content'] is not None]
            gc_data_by_dataset[dataset] = gc_contents
        
        # 检查是否有GC含量数据
        has_gc_data = any(len(gc_contents) > 0 for gc_contents in gc_data_by_dataset.values())
        
        if has_gc_data:
            # 子图1: GC含量分布直方图
            for i, (dataset, gc_contents) in enumerate(gc_data_by_dataset.items()):
                if gc_contents:
                    axes[0].hist(gc_contents, bins=30, alpha=0.7, label=dataset,
                                color=colors[i % len(colors)], density=True)
            
            axes[0].set_xlabel('GC Content')
            axes[0].set_ylabel('Density')
            axes[0].set_title('GC Content Distribution')
            axes[0].legend()
            
            # 子图2: GC含量箱线图
            gc_data_for_box = []
            gc_labels = []
            for dataset, gc_contents in gc_data_by_dataset.items():
                if gc_contents:
                    gc_data_for_box.append(gc_contents)
                    gc_labels.append(dataset)
            
            if gc_data_for_box:
                bp = axes[1].boxplot(gc_data_for_box, labels=gc_labels, patch_artist=True)
                for patch, color in zip(bp['boxes'], colors):
                    patch.set_facecolor(color)
            
            axes[1].set_ylabel('GC Content')
            axes[1].set_title('GC Content Comparison')
            axes[1].tick_params(axis='x', rotation=45)
            
            # 子图3: GC含量 vs 基因长度散点图
            for i, dataset in enumerate(datasets):
                features = all_features[dataset]
                gc_contents = []
                gene_lengths = []
                
                for f in features:
                    if f['gc_content'] is not None:
                        gc_contents.append(f['gc_content'])
                        gene_lengths.append(f['gene_length'])
                
                if gc_contents and gene_lengths:
                    axes[2].scatter(gene_lengths, gc_contents, alpha=0.6, 
                                   label=dataset, color=colors[i % len(colors)])
            
            axes[2].set_xlabel('Gene Length (bp)')
            axes[2].set_ylabel('GC Content')
            axes[2].set_title('GC Content vs Gene Length')
            axes[2].legend()
            
        else:
            # 如果没有GC含量数据，显示提示信息
            for ax in axes:
                ax.text(0.5, 0.5, 'No GC content data available\n(Genome sequences not provided)',
                       ha='center', va='center', transform=ax.transAxes,
                       fontsize=12, bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgray"))
                ax.set_title('GC Content Analysis')
        
        plt.tight_layout()
        
        # 保存图表
        plot_file = self.output_dir / f'{output_prefix}_gc_content_analysis.svg'
        plt.savefig(plot_file, format='svg', bbox_inches='tight')
        plt.close()
        
        logger.info(f"GC content analysis plot saved to {plot_file}")
    
    def perform_clustering_analysis(self, all_features, output_prefix):
        """执行聚类分析"""
        # 合并所有数据集的特征
        all_genes = []
        dataset_labels = []
        
        for dataset, features in all_features.items():
            for feature in features:
                all_genes.append(feature)
                dataset_labels.append(dataset)
        
        if len(all_genes) < 10:  # 需要足够的数据点进行聚类
            logger.warning("Not enough genes for clustering analysis")
            return
        
        # 准备特征矩阵
        feature_matrix = []
        valid_indices = []
        
        for i, gene in enumerate(all_genes):
            # 选择数值特征进行聚类
            features_vec = [
                gene['gene_length'],
                gene['num_exons'],
                gene['num_introns'],
                gene['avg_exon_length'] if gene['avg_exon_length'] else 0,
                gene['avg_intron_length'] if gene['avg_intron_length'] else 0,
                gene['coding_ratio'] if gene['coding_ratio'] else 0
            ]
            
            # 检查是否有缺失值
            if all(isinstance(x, (int, float)) and not np.isnan(x) for x in features_vec):
                feature_matrix.append(features_vec)
                valid_indices.append(i)
        
        if len(feature_matrix) < 10:
            logger.warning("Not enough valid features for clustering")
            return
        
        feature_matrix = np.array(feature_matrix)
        valid_labels = [dataset_labels[i] for i in valid_indices]
        
        # 标准化特征
        scaler = StandardScaler()
        features_scaled = scaler.fit_transform(feature_matrix)
        
        # K-means聚类
        n_clusters = min(5, len(set(valid_labels)))  # 最多5个聚类
        kmeans = KMeans(n_clusters=n_clusters, random_state=42)
        cluster_labels = kmeans.fit_predict(features_scaled)
        
        # 可视化聚类结果
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        # PCA降维可视化
        from sklearn.decomposition import PCA
        pca = PCA(n_components=2)
        features_pca = pca.fit_transform(features_scaled)
        
        # 子图1: 按数据集着色
        dataset_colors = {}
        colors = ['red', 'blue', 'green', 'orange', 'purple']
        unique_datasets = list(set(valid_labels))
        for i, dataset in enumerate(unique_datasets):
            dataset_colors[dataset] = colors[i % len(colors)]
        
        for dataset in unique_datasets:
            mask = [label == dataset for label in valid_labels]
            axes[0,0].scatter(features_pca[mask, 0], features_pca[mask, 1], 
                             c=dataset_colors[dataset], label=dataset, alpha=0.6)
        
        axes[0,0].set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.2%} variance)')
        axes[0,0].set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.2%} variance)')
        axes[0,0].set_title('PCA Visualization (by Dataset)')
        axes[0,0].legend()
        
        # 子图2: 按聚类着色
        cluster_colors = ['red', 'blue', 'green', 'orange', 'purple']
        for cluster in range(n_clusters):
            mask = cluster_labels == cluster
            axes[0,1].scatter(features_pca[mask, 0], features_pca[mask, 1],
                             c=cluster_colors[cluster % len(cluster_colors)], 
                             label=f'Cluster {cluster}', alpha=0.6)
        
        axes[0,1].set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.2%} variance)')
        axes[0,1].set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.2%} variance)')
        axes[0,1].set_title('K-means Clustering Result')
        axes[0,1].legend()
        
        # 子图3: 特征重要性
        feature_names = ['Gene Length', 'Num Exons', 'Num Introns', 
                        'Avg Exon Length', 'Avg Intron Length', 'Coding Ratio']
        
        pca_components = np.abs(pca.components_)
        
        axes[1,0].bar(range(len(feature_names)), pca_components[0], alpha=0.7, label='PC1')
        axes[1,0].bar(range(len(feature_names)), pca_components[1], alpha=0.7, label='PC2')
        axes[1,0].set_xlabel('Features')
        axes[1,0].set_ylabel('Absolute Loading')
        axes[1,0].set_title('PCA Feature Loadings')
        axes[1,0].set_xticks(range(len(feature_names)))
        axes[1,0].set_xticklabels(feature_names, rotation=45)
        axes[1,0].legend()
        
        # 子图4: 聚类中心特征
        cluster_centers = scaler.inverse_transform(kmeans.cluster_centers_)
        
        x_pos = np.arange(len(feature_names))
        width = 0.8 / n_clusters
        
        for i in range(n_clusters):
            axes[1,1].bar(x_pos + i*width, cluster_centers[i], width,
                         label=f'Cluster {i}', alpha=0.7,
                         color=cluster_colors[i % len(cluster_colors)])
        
        axes[1,1].set_xlabel('Features')
        axes[1,1].set_ylabel('Feature Value')
        axes[1,1].set_title('Cluster Centers')
        axes[1,1].set_xticks(x_pos + width * (n_clusters-1) / 2)
        axes[1,1].set_xticklabels(feature_names, rotation=45)
        axes[1,1].legend()
        axes[1,1].set_yscale('log')  # 使用对数刻度
        
        plt.tight_layout()
        
        # 保存图表
        plot_file = self.output_dir / f'{output_prefix}_clustering_analysis.svg'
        plt.savefig(plot_file, format='svg', bbox_inches='tight')
        plt.close()
        
        logger.info(f"Clustering analysis plot saved to {plot_file}")
        
        # 保存聚类结果
        clustering_results = {
            'n_clusters': n_clusters,
            'cluster_assignments': cluster_labels.tolist(),
            'dataset_labels': valid_labels,
            'cluster_centers': cluster_centers.tolist(),
            'pca_explained_variance': pca.explained_variance_ratio_.tolist(),
            'feature_names': feature_names
        }
        
        cluster_file = self.output_dir / f'{output_prefix}_clustering_results.json'
        with open(cluster_file, 'w') as f:
            json.dump(clustering_results, f, indent=2, default=self._json_serializer)
        
        logger.info(f"Clustering results saved to {cluster_file}")
    
    def generate_comparative_statistics(self, all_stats):
        """生成比较统计分析"""
        if len(all_stats) < 2:
            logger.warning("Need at least 2 datasets for comparative analysis")
            return {}
        
        comparative_results = {
            'statistical_tests': {},
            'effect_sizes': {},
            'summary_comparison': {}
        }
        
        datasets = list(all_stats.keys())
        features_to_compare = ['gene_length', 'num_exons', 'avg_exon_length', 'avg_intron_length']
        
        # 两两比较
        for i in range(len(datasets)):
            for j in range(i+1, len(datasets)):
                dataset1, dataset2 = datasets[i], datasets[j]
                comparison_key = f"{dataset1}_vs_{dataset2}"
                
                comparative_results['statistical_tests'][comparison_key] = {}
                comparative_results['effect_sizes'][comparison_key] = {}
                
                for feature in features_to_compare:
                    if (feature in all_stats[dataset1]['features'] and 
                        feature in all_stats[dataset2]['features']):
                        
                        stats1 = all_stats[dataset1]['features'][feature]
                        stats2 = all_stats[dataset2]['features'][feature]
                        
                        # t检验（假设正态分布）
                        # 这里使用统计量近似，实际应该使用原始数据
                        t_stat = (stats1['mean'] - stats2['mean']) / np.sqrt(
                            (stats1['std']**2 / stats1['count']) + 
                            (stats2['std']**2 / stats2['count'])
                        )
                        
                        # Cohen's d效应量
                        pooled_std = np.sqrt(
                            ((stats1['count']-1)*stats1['std']**2 + 
                             (stats2['count']-1)*stats2['std']**2) / 
                            (stats1['count'] + stats2['count'] - 2)
                        )
                        cohens_d = (stats1['mean'] - stats2['mean']) / pooled_std
                        
                        comparative_results['statistical_tests'][comparison_key][feature] = {
                            't_statistic': float(t_stat),
                            'degrees_of_freedom': stats1['count'] + stats2['count'] - 2,
                            'mean_difference': stats1['mean'] - stats2['mean'],
                            'std_error': np.sqrt(
                                (stats1['std']**2 / stats1['count']) + 
                                (stats2['std']**2 / stats2['count'])
                            )
                        }
                        
                        comparative_results['effect_sizes'][comparison_key][feature] = {
                            'cohens_d': float(cohens_d),
                            'interpretation': self._interpret_effect_size(abs(cohens_d))
                        }
        
        # 汇总比较
        for feature in features_to_compare:
            feature_means = []
            feature_datasets = []
            
            for dataset in datasets:
                if feature in all_stats[dataset]['features']:
                    feature_means.append(all_stats[dataset]['features'][feature]['mean'])
                    feature_datasets.append(dataset)
            
            if len(feature_means) > 1:
                comparative_results['summary_comparison'][feature] = {
                    'datasets': feature_datasets,
                    'means': feature_means,
                    'overall_mean': np.mean(feature_means),
                    'coefficient_of_variation': np.std(feature_means) / np.mean(feature_means),
                    'min_dataset': feature_datasets[np.argmin(feature_means)],
                    'max_dataset': feature_datasets[np.argmax(feature_means)],
                    'range_ratio': max(feature_means) / min(feature_means) if min(feature_means) > 0 else None
                }
        
        return comparative_results
    
    def _interpret_effect_size(self, d):
        """解释效应量大小"""
        if d < 0.2:
            return "negligible"
        elif d < 0.5:
            return "small"
        elif d < 0.8:
            return "medium"
        else:
            return "large"
    
    def analyze_multiple_datasets(self, gff_files, genome_files=None, dataset_names=None):
        """分析多个数据集"""
        if dataset_names is None:
            dataset_names = [f"Dataset_{i+1}" for i in range(len(gff_files))]
        
        if genome_files is None:
            genome_files = [None] * len(gff_files)
        
        all_features = {}
        all_stats = {}
        
        # 分析每个数据集
        for i, (gff_file, genome_file, dataset_name) in enumerate(zip(gff_files, genome_files, dataset_names)):
            logger.info(f"Analyzing dataset {i+1}/{len(gff_files)}: {dataset_name}")
            
            # 加载基因组序列（如果提供）
            genome_sequences = {}
            if genome_file:
                genome_sequences = self.load_genome_sequences(genome_file)
            
            # 提取特征
            features = self.extract_gene_features_from_gff(gff_file, genome_sequences)
            
            if features:
                all_features[dataset_name] = features
                
                # 计算统计量
                stats = self.calculate_basic_statistics(features, dataset_name)
                all_stats[dataset_name] = stats
                
                # 保存单个数据集的统计结果
                stats_file = self.output_dir / f'gene_stats_{dataset_name}.json'
                with open(stats_file, 'w') as f:
                    json.dump(stats, f, indent=2, default=self._json_serializer)
        
        if not all_features:
            logger.error("No valid features extracted from any dataset")
            return
        
        # 生成比较分析
        logger.info("Generating comparative visualizations...")
        
        # 基因长度分析
        self.analyze_gene_length_distribution(all_features, 'comparative')
        
        # 外显子-内含子结构分析
        self.analyze_exon_intron_structure(all_features, 'comparative')
        
        # GC含量分析
        self.analyze_gc_content(all_features, 'comparative')
        
        # 聚类分析
        self.perform_clustering_analysis(all_features, 'comparative')
        
        # 比较统计分析
        comparative_stats = self.generate_comparative_statistics(all_stats)
        
        # 生成最终报告
        final_report = {
            'analysis_summary': {
                'datasets_analyzed': len(all_features),
                'total_genes': sum(len(features) for features in all_features.values()),
                'analysis_date': datetime.now().isoformat()
            },
            'individual_stats': all_stats,
            'comparative_analysis': comparative_stats,
            'dataset_summary': {
                dataset: {
                    'gene_count': len(features),
                    'avg_gene_length': np.mean([f['gene_length'] for f in features]),
                    'avg_exon_count': np.mean([f['num_exons'] for f in features])
                }
                for dataset, features in all_features.items()
            }
        }
        
        # 保存最终报告
        report_file = self.output_dir / 'gene_features_analysis_report.json'
        with open(report_file, 'w') as f:
            json.dump(final_report, f, indent=2, default=self._json_serializer)
        
        # 保存CSV汇总
        summary_data = []
        for dataset, features in all_features.items():
            for feature in features:
                row = {
                    'Dataset': dataset,
                    'Gene_ID': feature['gene_id'],
                    'Gene_Length': feature['gene_length'],
                    'Num_Exons': feature['num_exons'],
                    'Num_Introns': feature['num_introns'],
                    'Avg_Exon_Length': feature['avg_exon_length'],
                    'Avg_Intron_Length': feature['avg_intron_length'],
                    'Coding_Ratio': feature['coding_ratio'],
                    'GC_Content': feature['gc_content']
                }
                summary_data.append(row)
        
        df = pd.DataFrame(summary_data)
        csv_file = self.output_dir / 'gene_features_summary.csv'
        df.to_csv(csv_file, index=False)
        
        logger.info(f"Analysis completed. Report saved to {report_file}")
        logger.info(f"Summary CSV saved to {csv_file}")
        
        return final_report

def main():
    parser = argparse.ArgumentParser(description='Gene Features Statistical Analysis')
    parser.add_argument('--gff-files', nargs='+', required=True, 
                       help='GFF annotation files to analyze')
    parser.add_argument('--genome-files', nargs='*', 
                       help='Genome sequence files (optional, for GC content analysis)')
    parser.add_argument('--dataset-names', nargs='*',
                       help='Names for the datasets (optional)')
    parser.add_argument('--output', required=True, help='Output directory')
    
    args = parser.parse_args()
    
    # 创建分析器
    analyzer = GeneFeatureAnalyzer(args.output)
    
    # 运行分析
    report = analyzer.analyze_multiple_datasets(
        gff_files=args.gff_files,
        genome_files=args.genome_files,
        dataset_names=args.dataset_names
    )
    
    if report:
        print("Gene features analysis completed successfully!")
        print(f"Results saved to: {args.output}")
    else:
        print("Gene features analysis failed")
        sys.exit(1)

if __name__ == '__main__':
    main()