"""
Augustus基因预测模块
Augustus Gene Prediction Module

本模块提供对Augustus基因预测工具的封装，支持：
- 批量基因预测
- 候选区域处理
- GFF3输出解析
- 预测质量评估
- 并行处理支持
"""

import os
import subprocess
import logging
from pathlib import Path
from typing import List, Dict, Optional, Tuple, Union
from dataclasses import dataclass, field
from concurrent.futures import ProcessPoolExecutor, as_completed
import tempfile
import shutil
from datetime import datetime
import json

try:
    import gffutils
    GFFUTILS_AVAILABLE = True
except ImportError:
    GFFUTILS_AVAILABLE = False
    logging.warning("gffutils not available - using basic GFF3 parsing")

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


@dataclass
class PredictionConfig:
    """Augustus预测配置"""
    species: str = "nbs_rice_v1"
    augustus_config_path: str = "/home/wangys/opt/Augustus/config"
    output_gff3: bool = True
    alternatives_from_evidence: bool = False
    allow_hinted_splicesites: str = "atac"
    min_intron_len: int = 20
    max_intron_len: int = 50000
    strand: str = "both"  # both, forward, backward
    genemodel: str = "partial"  # complete, partial, intronless, atleastone, exactlyone
    singlestrand: bool = False
    uniqueGeneId: bool = True
    protein: bool = True
    codingseq: bool = True
    cds: bool = True
    sample: int = 100
    alternatives_from_sampling: bool = False
    maxtracks: int = 3
    temperature: float = 2.0
    additional_params: Dict[str, str] = field(default_factory=dict)


@dataclass
class PredictionResult:
    """单个预测结果"""
    sequence_id: str
    sequence_length: int
    genes_predicted: int
    genes: List[Dict] = field(default_factory=list)
    gff3_file: Optional[str] = None
    protein_file: Optional[str] = None
    cds_file: Optional[str] = None
    prediction_time: float = 0.0
    quality_metrics: Dict = field(default_factory=dict)
    warnings: List[str] = field(default_factory=list)
    errors: List[str] = field(default_factory=list)


@dataclass
class BatchPredictionResult:
    """批量预测结果"""
    total_sequences: int
    successful_predictions: int
    failed_predictions: int
    total_genes: int
    total_time: float
    results: List[PredictionResult] = field(default_factory=list)
    summary_stats: Dict = field(default_factory=dict)
    quality_report: Dict = field(default_factory=dict)


class GFF3Parser:
    """GFF3文件解析器"""
    
    def __init__(self, use_gffutils: bool = True):
        self.use_gffutils = use_gffutils and GFFUTILS_AVAILABLE
        self.logger = logging.getLogger(__name__)
    
    def parse_gff3_file(self, gff3_file: str) -> List[Dict]:
        """解析GFF3文件，提取基因结构信息"""
        if not os.path.exists(gff3_file):
            raise FileNotFoundError(f"GFF3 file not found: {gff3_file}")
        
        if self.use_gffutils:
            return self._parse_with_gffutils(gff3_file)
        else:
            return self._parse_basic(gff3_file)
    
    def _parse_basic(self, gff3_file: str) -> List[Dict]:
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
                
                # 解析属性
                attr_dict = {}
                for attr in attributes.split(';'):
                    if '=' in attr:
                        key, value = attr.split('=', 1)
                        attr_dict[key] = value
                
                if feature == 'gene':
                    gene_id = attr_dict.get('gene_id', attr_dict.get('ID', f"gene_{len(genes)}"))
                    genes[gene_id] = {
                        'gene_id': gene_id,
                        'seqname': seqname,
                        'start': int(start),
                        'end': int(end),
                        'strand': strand,
                        'length': int(end) - int(start) + 1,
                        'transcripts': []
                    }
        
        return list(genes.values())


class AugustusPredictor:
    """Augustus基因预测器"""
    
    def __init__(self, config: Optional[PredictionConfig] = None):
        """
        初始化Augustus预测器
        
        Args:
            config: 预测配置，如果为None则使用默认配置
        """
        self.config = config or PredictionConfig()
        self.logger = logging.getLogger(__name__)
        self.gff3_parser = GFF3Parser()
        
        # 验证Augustus安装
        self._validate_augustus_installation()
        
        # 验证物种模型
        self._validate_species_model()
    
    def _validate_augustus_installation(self):
        """验证Augustus安装"""
        try:
            result = subprocess.run(['augustus', '--version'], 
                                  capture_output=True, text=True, timeout=10)
            if result.returncode == 0:
                self.logger.info(f"Augustus version: {result.stdout.strip()}")
            else:
                raise RuntimeError("Augustus not properly installed")
        except (subprocess.TimeoutExpired, FileNotFoundError) as e:
            raise RuntimeError(f"Augustus not found or not working: {e}")
    
    def _validate_species_model(self):
        """验证物种模型"""
        species_path = Path(self.config.augustus_config_path) / "species" / self.config.species
        if not species_path.exists():
            available_species = list((Path(self.config.augustus_config_path) / "species").glob("*"))
            raise ValueError(f"Species model '{self.config.species}' not found. "
                           f"Available models: {[s.name for s in available_species]}")
        
        # 检查关键模型文件
        required_files = [
            f"{self.config.species}_parameters.cfg",
            f"{self.config.species}_exon_probs.pbl",
            f"{self.config.species}_intron_probs.pbl"
        ]
        
        missing_files = []
        for file in required_files:
            if not (species_path / file).exists():
                missing_files.append(file)
        
        if missing_files:
            self.logger.warning(f"Missing model files for {self.config.species}: {missing_files}")
    
    def _build_augustus_command(self, input_file: str, output_dir: str, 
                               sequence_id: str = None) -> List[str]:
        """构建Augustus命令"""
        cmd = ['augustus']
        
        # 基本参数
        cmd.extend([f'--species={self.config.species}'])
        cmd.extend([f'--AUGUSTUS_CONFIG_PATH={self.config.augustus_config_path}'])
        
        # 输出格式
        if self.config.output_gff3:
            cmd.extend(['--gff3=on'])
        else:
            cmd.extend(['--gff3=off'])
        
        # 其他参数 (使用正确的Augustus参数名)
        cmd.extend([f'--alternatives-from-evidence={str(self.config.alternatives_from_evidence).lower()}'])
        cmd.extend([f'--allow_hinted_splicesites={self.config.allow_hinted_splicesites}'])
        cmd.extend([f'--strand={self.config.strand}'])
        cmd.extend([f'--genemodel={self.config.genemodel}'])
        
        # 注意: min_intron_len 和 max_intron_len 在Augustus中可能不是命令行参数
        # 这些通常在物种配置文件中设置
        
        if self.config.singlestrand:
            cmd.extend(['--singlestrand=true'])
        
        if self.config.uniqueGeneId:
            cmd.extend(['--uniqueGeneId=true'])
        
        if self.config.protein:
            cmd.extend(['--protein=on'])
        
        if self.config.codingseq:
            cmd.extend(['--codingseq=on'])
        
        if self.config.cds:
            cmd.extend(['--cds=on'])
        
        # 添加自定义参数
        for key, value in self.config.additional_params.items():
            cmd.extend([f'--{key}={value}'])
        
        # 输入文件
        cmd.append(input_file)
        
        return cmd
    
    def predict_single_sequence(self, sequence_file: str, output_dir: str, 
                               sequence_id: str = None) -> PredictionResult:
        """
        对单个序列进行基因预测
        
        Args:
            sequence_file: 输入序列文件路径
            output_dir: 输出目录
            sequence_id: 序列ID，如果为None则从文件名推断
        
        Returns:
            PredictionResult: 预测结果
        """
        start_time = datetime.now()
        
        if sequence_id is None:
            sequence_id = Path(sequence_file).stem
        
        # 创建输出目录
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        # 准备输出文件名
        gff3_file = output_path / f"{sequence_id}.gff3"
        log_file = output_path / f"{sequence_id}.log"
        
        # 获取序列长度
        try:
            with open(sequence_file, 'r') as f:
                sequences = list(SeqIO.parse(f, 'fasta'))
                sequence_length = sum(len(seq) for seq in sequences)
        except Exception as e:
            return PredictionResult(
                sequence_id=sequence_id,
                sequence_length=0,
                genes_predicted=0,
                errors=[f"Failed to read input file: {e}"]
            )
        
        # 构建Augustus命令
        cmd = self._build_augustus_command(sequence_file, str(output_path), sequence_id)
        
        # 执行预测
        try:
            with open(gff3_file, 'w') as gff_out, \
                 open(log_file, 'w') as log_out:
                
                self.logger.info(f"Running Augustus prediction for {sequence_id}")
                self.logger.debug(f"Command: {' '.join(cmd)}")
                
                result = subprocess.run(
                    cmd,
                    stdout=gff_out,
                    stderr=log_out,
                    text=True,
                    timeout=3600  # 1小时超时
                )
                
                if result.returncode != 0:
                    with open(log_file, 'r') as f:
                        error_msg = f.read()
                    return PredictionResult(
                        sequence_id=sequence_id,
                        sequence_length=sequence_length,
                        genes_predicted=0,
                        errors=[f"Augustus failed with return code {result.returncode}: {error_msg}"]
                    )
        
        except subprocess.TimeoutExpired:
            return PredictionResult(
                sequence_id=sequence_id,
                sequence_length=sequence_length,
                genes_predicted=0,
                errors=["Augustus prediction timed out"]
            )
        except Exception as e:
            return PredictionResult(
                sequence_id=sequence_id,
                sequence_length=sequence_length,
                genes_predicted=0,
                errors=[f"Unexpected error during prediction: {e}"]
            )
        
        # 解析结果
        try:
            genes = self.gff3_parser.parse_gff3_file(str(gff3_file))
            genes_predicted = len(genes)
            
            prediction_result = PredictionResult(
                sequence_id=sequence_id,
                sequence_length=sequence_length,
                genes_predicted=genes_predicted,
                genes=genes,
                gff3_file=str(gff3_file),
                prediction_time=(datetime.now() - start_time).total_seconds()
            )
            
            self.logger.info(f"Prediction completed for {sequence_id}: {genes_predicted} genes predicted")
            
            return prediction_result
            
        except Exception as e:
            return PredictionResult(
                sequence_id=sequence_id,
                sequence_length=sequence_length,
                genes_predicted=0,
                errors=[f"Failed to parse prediction results: {e}"]
            )
    
    def predict_batch(self, input_sequences: List[str], output_dir: str, 
                     max_workers: int = None) -> BatchPredictionResult:
        """
        批量基因预测
        
        Args:
            input_sequences: 输入序列文件列表
            output_dir: 输出目录
            max_workers: 最大并行工作进程数
        
        Returns:
            BatchPredictionResult: 批量预测结果
        """
        start_time = datetime.now()
        
        if max_workers is None:
            max_workers = min(len(input_sequences), os.cpu_count())
        
        self.logger.info(f"Starting batch prediction with {max_workers} workers")
        
        results = []
        successful = 0
        failed = 0
        
        # 串行处理（避免并行问题）
        for seq_file in input_sequences:
            try:
                result = self.predict_single_sequence(seq_file, output_dir)
                results.append(result)
                
                if result.errors:
                    failed += 1
                    self.logger.error(f"Prediction failed for {seq_file}: {result.errors}")
                else:
                    successful += 1
                    self.logger.info(f"Prediction successful for {seq_file}: {result.genes_predicted} genes")
                    
            except Exception as e:
                failed += 1
                self.logger.error(f"Exception during prediction of {seq_file}: {e}")
                results.append(PredictionResult(
                    sequence_id=Path(seq_file).stem,
                    sequence_length=0,
                    genes_predicted=0,
                    errors=[str(e)]
                ))
        
        # 计算总体统计
        total_genes = sum(r.genes_predicted for r in results)
        total_time = (datetime.now() - start_time).total_seconds()
        
        batch_result = BatchPredictionResult(
            total_sequences=len(input_sequences),
            successful_predictions=successful,
            failed_predictions=failed,
            total_genes=total_genes,
            total_time=total_time,
            results=results
        )
        
        self.logger.info(f"Batch prediction completed: {successful}/{len(input_sequences)} successful, "
                        f"{total_genes} total genes predicted in {total_time:.2f}s")
        
        return batch_result
    
    def predict_candidate_regions(self, candidate_regions: List[Dict], 
                                 genome_file: str, output_dir: str,
                                 max_workers: int = None) -> BatchPredictionResult:
        """
        对候选区域进行基因预测
        
        Args:
            candidate_regions: 候选区域列表，每个区域包含seqname, start, end
            genome_file: 基因组文件路径
            output_dir: 输出目录
            max_workers: 最大并行工作进程数
        
        Returns:
            BatchPredictionResult: 批量预测结果
        """
        self.logger.info(f"Predicting genes in {len(candidate_regions)} candidate regions")
        
        # 创建临时目录存储候选区域序列
        temp_dir = Path(output_dir) / "temp_regions"
        temp_dir.mkdir(parents=True, exist_ok=True)
        
        try:
            # 提取候选区域序列
            region_files = self._extract_candidate_sequences(
                candidate_regions, genome_file, str(temp_dir)
            )
            
            # 批量预测
            batch_result = self.predict_batch(region_files, output_dir, max_workers)
            
            # 更新结果中的坐标信息
            self._update_coordinates_in_results(batch_result, candidate_regions)
            
            return batch_result
            
        finally:
            # 清理临时文件
            if temp_dir.exists():
                shutil.rmtree(temp_dir)
    
    def _extract_candidate_sequences(self, candidate_regions: List[Dict], 
                                   genome_file: str, temp_dir: str) -> List[str]:
        """提取候选区域序列"""
        # 读取基因组序列
        genome_sequences = {}
        with open(genome_file, 'r') as f:
            for record in SeqIO.parse(f, 'fasta'):
                genome_sequences[record.id] = record
        
        region_files = []
        
        for i, region in enumerate(candidate_regions):
            seqname = region['seqname']
            start = region['start']
            end = region['end']
            
            if seqname not in genome_sequences:
                self.logger.warning(f"Sequence {seqname} not found in genome file")
                continue
            
            # 提取区域序列 (转换为0-based索引)
            seq_record = genome_sequences[seqname]
            region_seq = seq_record.seq[start-1:end]
            
            # 创建新的序列记录
            region_id = f"region_{i}_{seqname}_{start}_{end}"
            region_record = SeqRecord(
                region_seq,
                id=region_id,
                description=f"Candidate region {seqname}:{start}-{end}"
            )
            
            # 保存到文件
            region_file = Path(temp_dir) / f"{region_id}.fasta"
            with open(region_file, 'w') as f:
                SeqIO.write(region_record, f, 'fasta')
            
            region_files.append(str(region_file))
        
        return region_files
    
    def _update_coordinates_in_results(self, batch_result: BatchPredictionResult, 
                                     candidate_regions: List[Dict]):
        """更新预测结果中的基因组坐标"""
        for result in batch_result.results:
            # 从sequence_id中提取区域信息
            if result.sequence_id.startswith('region_'):
                parts = result.sequence_id.split('_')
                if len(parts) >= 5:
                    region_index = int(parts[1])
                    if 0 <= region_index < len(candidate_regions):
                        region = candidate_regions[region_index]
                        region_start = region['start']
                        
                        # 更新基因坐标
                        for gene in result.genes:
                            gene['start'] += region_start - 1
                            gene['end'] += region_start - 1
                            gene['seqname'] = region['seqname']
    
    def save_batch_results(self, batch_result: BatchPredictionResult, 
                          output_file: str, format: str = 'json'):
        """
        保存批量预测结果
        
        Args:
            batch_result: 批量预测结果
            output_file: 输出文件路径
            format: 输出格式 ('json' 或 'tsv')
        """
        if format.lower() == 'json':
            self._save_results_json(batch_result, output_file)
        elif format.lower() == 'tsv':
            self._save_results_tsv(batch_result, output_file)
        else:
            raise ValueError(f"Unsupported format: {format}")
    
    def _save_results_json(self, batch_result: BatchPredictionResult, output_file: str):
        """保存JSON格式结果"""
        # 转换为可序列化的字典
        results_dict = {
            'summary': {
                'total_sequences': batch_result.total_sequences,
                'successful_predictions': batch_result.successful_predictions,
                'failed_predictions': batch_result.failed_predictions,
                'total_genes': batch_result.total_genes,
                'total_time': batch_result.total_time,
                'genes_per_sequence': batch_result.total_genes / batch_result.total_sequences if batch_result.total_sequences > 0 else 0
            },
            'predictions': []
        }
        
        for result in batch_result.results:
            pred_dict = {
                'sequence_id': result.sequence_id,
                'sequence_length': result.sequence_length,
                'genes_predicted': result.genes_predicted,
                'prediction_time': result.prediction_time,
                'genes': result.genes,
                'warnings': result.warnings,
                'errors': result.errors,
                'quality_metrics': result.quality_metrics
            }
            
            # 添加候选区域信息（如果存在）
            if hasattr(result, 'candidate_info'):
                pred_dict['candidate_info'] = result.candidate_info
            
            results_dict['predictions'].append(pred_dict)
        
        with open(output_file, 'w') as f:
            json.dump(results_dict, f, indent=2)
        
        self.logger.info(f"Results saved to JSON file: {output_file}")
    
    def _save_results_tsv(self, batch_result: BatchPredictionResult, output_file: str):
        """保存TSV格式结果"""
        with open(output_file, 'w') as f:
            # 写入头部
            headers = [
                'sequence_id', 'sequence_length', 'genes_predicted', 
                'prediction_time', 'has_errors', 'has_warnings',
                'candidate_seqname', 'candidate_start', 'candidate_end',
                'candidate_identity'
            ]
            f.write('\t'.join(headers) + '\n')
            
            # 写入数据
            for result in batch_result.results:
                row = [
                    result.sequence_id,
                    str(result.sequence_length),
                    str(result.genes_predicted),
                    f"{result.prediction_time:.2f}",
                    str(len(result.errors) > 0),
                    str(len(result.warnings) > 0)
                ]
                
                # 添加候选区域信息
                if hasattr(result, 'candidate_info') and result.candidate_info:
                    info = result.candidate_info
                    row.extend([
                        info.get('seqname', ''),
                        str(info.get('start', '')),
                        str(info.get('end', '')),
                        str(info.get('identity', ''))
                    ])
                else:
                    row.extend(['', '', '', ''])
                
                f.write('\t'.join(row) + '\n')
        
        self.logger.info(f"Results saved to TSV file: {output_file}")
 