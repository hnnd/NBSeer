#!/usr/bin/env python3
"""
NLR-Annotator集成模块

该模块提供了NLR-Annotator工具的Python封装，用于在基因组中识别NLR基因。
NLR-Annotator是一个基于Java的工具，用于识别NLR (NBS-LRR) 抗病基因。

主要功能：
1. 基因组中NLR基因的识别
2. 基于motif模式的NLR基因搜索
3. 多种输出格式支持（tabular、GFF3、BED）
4. 候选基因位点的多格式输出
5. 并行处理大基因组支持

作者: NBS基因注释项目组
创建时间: 2025-06-22
"""

import os
import subprocess
import shutil
import json
from pathlib import Path
from typing import Optional, List, Dict, Tuple, Union
from dataclasses import dataclass, asdict
import tempfile
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
import pandas as pd
from Bio import SeqIO
import time

try:
    from ..utils.config import get_config
    from ..utils.logging_setup import get_logger
    from ..utils.data_validation import FastaValidator
except ImportError:
    # Fallback for standalone usage
    def get_config(path=None):
        return {"paths": {"tools_dir": "tools"}}
    
    def get_logger(name):
        import logging
        return logging.getLogger(name)
    
    class FastaValidator:
        def __init__(self, logger):
            self.logger = logger
        
        def validate_file(self, file_path):
            return os.path.exists(file_path)


@dataclass
class NLRAnnotatorConfig:
    """NLR-Annotator配置参数"""
    input_fasta: str = ""
    output_dir: str = ""
    mot_file: str = "tools/NLR-Annotator/mot.txt"
    store_file: str = "tools/NLR-Annotator/store.txt"
    jar_file: str = "tools/NLR-Annotator/NLR-Annotator-v2.1b.jar"
    threads: int = 8
    memory_limit: str = "8G"
    sequences_per_thread: int = 1000
    overwrite: bool = False
    output_formats: List[str] = None
    
    # Legacy support
    @property
    def genome_file(self):
        return self.input_fasta
    
    @genome_file.setter
    def genome_file(self, value):
        self.input_fasta = value
    
    @property
    def memory_mb(self):
        if self.memory_limit.endswith('G'):
            return int(self.memory_limit[:-1]) * 1000
        elif self.memory_limit.endswith('M'):
            return int(self.memory_limit[:-1])
        else:
            return int(self.memory_limit)
    
    def __post_init__(self):
        """验证参数"""
        if self.input_fasta and not os.path.exists(self.input_fasta):
            raise FileNotFoundError(f"基因组文件不存在: {self.input_fasta}")
        
        if self.jar_file and not os.path.exists(self.jar_file):
            raise FileNotFoundError(f"NLR-Annotator JAR文件不存在: {self.jar_file}")
        
        if self.mot_file and not os.path.exists(self.mot_file):
            raise FileNotFoundError(f"Motif配置文件不存在: {self.mot_file}")
        
        if self.store_file and not os.path.exists(self.store_file):
            raise FileNotFoundError(f"Storage配置文件不存在: {self.store_file}")
        
        if self.threads < 1:
            raise ValueError("线程数必须大于0")
        
        if self.memory_mb < 1000:
            raise ValueError("内存限制必须至少1000MB")
        
        if self.sequences_per_thread < 1:
            raise ValueError("每线程序列数必须大于0")
        
        if self.output_formats is None:
            self.output_formats = ["tabular", "gff", "bed"]


@dataclass
class NLRGene:
    """NLR gene information"""
    gene_id: str
    sequence_id: str
    start: int
    end: int
    strand: str = "+"
    gene_type: str = "NLR"
    score: float = 0.0
    motifs: List[dict] = None
    length: int = 0
    
    def __post_init__(self):
        if self.motifs is None:
            self.motifs = []
        if self.length == 0:
            self.length = self.end - self.start + 1
    
    # Legacy compatibility properties
    @property
    def chromosome(self) -> str:
        """Legacy compatibility - chromosome is sequence_id"""
        return self.sequence_id
    
    @chromosome.setter
    def chromosome(self, value: str):
        self.sequence_id = value
    
    def to_bed_line(self) -> str:
        """Convert to BED format line"""
        return f"{self.sequence_id}\t{self.start}\t{self.end}\t{self.gene_type}\t{self.score}\t{self.strand}"
    
    def to_gff_line(self) -> str:
        """Convert to GFF format line"""
        motif_names = [m.get('name', '') for m in self.motifs] if self.motifs else []
        attributes = f"ID={self.gene_id};Type={self.gene_type};Motifs={','.join(motif_names)}"
        return f"{self.sequence_id}\tNLR-Annotator\tgene\t{self.start}\t{self.end}\t{self.score}\t{self.strand}\t.\t{attributes}"
    
    @property
    def coordinates(self) -> str:
        """Gene coordinates string"""
        return f"{self.sequence_id}:{self.start}-{self.end}"


@dataclass
class NLRAnnotatorResult:
    """NLR-Annotator analysis result"""
    genes: List[NLRGene]
    total_genes: int
    execution_time: float
    command: str
    success: bool = True
    error_message: str = ""
    
    # Legacy compatibility
    output_files: Dict[str, str] = None
    config: NLRAnnotatorConfig = None
    java_version: str = ""
    
    def __post_init__(self):
        if self.output_files is None:
            self.output_files = {}
        if self.total_genes == 0:
            self.total_genes = len(self.genes)
    
    def to_bed_file(self, output_path: str) -> None:
        """导出为BED格式文件"""
        with open(output_path, 'w') as f:
            f.write("# NLR-Annotator候选基因结果\n")
            f.write("# 格式: chromosome\tstart\tend\tgene_type\tscore\tstrand\n")
            for gene in self.genes:
                f.write(gene.to_bed_line() + "\n")
    
    def to_gff_file(self, output_path: str) -> None:
        """导出为GFF格式文件"""
        with open(output_path, 'w') as f:
            f.write("##gff-version 3\n")
            f.write("# NLR-Annotator基因预测结果\n")
            for gene in self.genes:
                f.write(gene.to_gff_line() + "\n")
    
    def to_json_file(self, output_path: str) -> None:
        """导出为JSON格式文件"""
        data = {
            "metadata": {
                "total_genes": self.total_genes,
                "execution_time": self.execution_time,
                "java_version": self.java_version,
                "genome_file": self.config.genome_file,
                "output_directory": self.config.output_dir
            },
            "genes": [asdict(gene) for gene in self.genes]
        }
        
        with open(output_path, 'w') as f:
            json.dump(data, f, indent=2, ensure_ascii=False)
    
    def get_summary(self) -> Dict:
        """获取结果摘要"""
        return {
            "total_genes": self.total_genes,
            "execution_time": self.execution_time,
            "genome_file": self.config.genome_file,
            "output_directory": self.config.output_dir,
            "threads_used": self.config.threads,
            "memory_used_mb": self.config.memory_mb,
            "gene_types": self._get_type_distribution()
        }
    
    def _get_type_distribution(self) -> Dict[str, int]:
        """获取基因类型分布"""
        type_counts = {}
        for gene in self.genes:
            gene_type = gene.gene_type
            type_counts[gene_type] = type_counts.get(gene_type, 0) + 1
        return type_counts


class NLRAnnotatorRunner:
    """
    NLR-Annotator工具的Python封装类
    
    该类封装了NLR-Annotator Java工具的调用，提供了Python接口
    用于在基因组中识别NLR基因。
    
    主要功能：
    - 配置和执行NLR-Annotator搜索
    - 解析多种输出结果格式
    - 并行处理支持
    - 结果格式转换
    """
    
    def __init__(self, config_path: Optional[str] = None):
        """
        初始化NLRAnnotatorRunner
        
        Args:
            config_path: 配置文件路径，如果为None则使用默认配置
        """
        self.config = get_config(config_path)
        self.logger = get_logger(__name__)
        self.fasta_validator = FastaValidator(self.logger)
        
        # 获取工具路径
        self.jar_file = self._get_jar_file()
        self.mot_file = self._get_mot_file()
        self.store_file = self._get_store_file()
        
        self.logger.info("NLRAnnotatorRunner初始化完成")
    
    def validate_installation(self) -> bool:
        """验证NLR-Annotator安装"""
        try:
            self._validate_installation()
            return True
        except Exception as e:
            self.logger.error(f"安装验证失败: {e}")
            return False
    
    def _get_jar_file(self) -> str:
        """获取NLR-Annotator JAR文件路径"""
        # 从配置中获取路径
        jar_path = self.config.get("tools.nlr_annotator.jar_path") if hasattr(self.config, 'get') else None
        
        if jar_path and os.path.exists(jar_path):
            return jar_path
        
        # 尝试在工具目录中查找
        tools_dir = self.config.get("paths.tools_dir", "tools") if hasattr(self.config, 'get') else "tools"
        jar_tool_path = os.path.join(tools_dir, "NLR-Annotator", "NLR-Annotator-v2.1b.jar")
        
        if os.path.exists(jar_tool_path):
            return jar_tool_path
        
        # Default fallback path
        default_path = "tools/NLR-Annotator/NLR-Annotator-v2.1b.jar"
        if os.path.exists(default_path):
            return default_path
        
        raise FileNotFoundError(
            "无法找到NLR-Annotator JAR文件。请检查配置或运行setup_nlr_annotator.py。"
        )
    
    def _get_mot_file(self) -> str:
        """获取motif配置文件路径"""
        tools_dir = self.config.get("paths.tools_dir", "tools") if hasattr(self.config, 'get') else "tools"
        mot_path = os.path.join(tools_dir, "NLR-Annotator", "mot.txt")
        
        if os.path.exists(mot_path):
            return mot_path
        
        # Default fallback path
        default_path = "tools/NLR-Annotator/mot.txt"
        if os.path.exists(default_path):
            return default_path
        
        raise FileNotFoundError(
            "无法找到mot.txt配置文件。请运行setup_nlr_annotator.py。"
        )
    
    def _get_store_file(self) -> str:
        """获取storage配置文件路径"""
        tools_dir = self.config.get("paths.tools_dir", "tools") if hasattr(self.config, 'get') else "tools"
        store_path = os.path.join(tools_dir, "NLR-Annotator", "store.txt")
        
        if os.path.exists(store_path):
            return store_path
        
        # Default fallback path
        default_path = "tools/NLR-Annotator/store.txt"
        if os.path.exists(default_path):
            return default_path
        
        raise FileNotFoundError(
            "无法找到store.txt配置文件。请运行setup_nlr_annotator.py。"
        )
    
    def _validate_installation(self) -> None:
        """验证NLR-Annotator安装"""
        self.logger.info("验证NLR-Annotator安装")
        
        # 检查Java环境
        try:
            result = subprocess.run(['java', '-version'], 
                                  capture_output=True, text=True, check=True)
            java_version = result.stderr.split('\n')[0]
            self.logger.info(f"Java环境检查通过: {java_version}")
        except (subprocess.CalledProcessError, FileNotFoundError):
            raise RuntimeError("未找到Java运行环境。请安装Java Runtime Environment (JRE) 1.8或更高版本")
        
        # 检查JAR文件
        if not os.path.exists(self.jar_file):
            raise FileNotFoundError(f"NLR-Annotator JAR文件不存在: {self.jar_file}")
        
        # 检查配置文件
        if not os.path.exists(self.mot_file):
            raise FileNotFoundError(f"Motif配置文件不存在: {self.mot_file}")
        
        if not os.path.exists(self.store_file):
            raise FileNotFoundError(f"Storage配置文件不存在: {self.store_file}")
        
        self.logger.info("NLR-Annotator安装验证完成")
    
    def validate_input(self, input_fasta: str) -> bool:
        """Validate input FASTA file"""
        self.logger.info(f"Validating input file: {input_fasta}")
        
        if hasattr(self.fasta_validator, 'validate_fasta_format'):
            is_valid, error_msg, seq_count = self.fasta_validator.validate_fasta_format(input_fasta)
            if not is_valid:
                self.logger.error(f"Input file validation failed: {error_msg}")
                return False
            self.logger.info(f"Input file validation passed, contains {seq_count} sequences")
            return True
        else:
            # Fallback validation
            if not self.fasta_validator.validate_file(input_fasta):
                self.logger.error(f"Input file does not exist or is invalid: {input_fasta}")
                return False
            self.logger.info(f"Input file validation passed: {input_fasta}")
            return True
    
    def _build_command(self, input_fasta: str, output_dir: str, output_prefix: str) -> List[str]:
        """Build NLR-Annotator command"""
        output_base = os.path.join(output_dir, output_prefix)
        
        # Use default values for command parameters
        memory_limit = "8G"
        threads = 8
        sequences_per_thread = 1000
        
        cmd = [
            'java',
            f'-Xmx{memory_limit}',
            '-jar', os.path.abspath(self.jar_file),
            '-i', os.path.abspath(input_fasta),
            '-x', os.path.abspath(self.mot_file),
            '-y', os.path.abspath(self.store_file),
            '-t', str(threads),
            '-n', str(sequences_per_thread)
        ]
        
        # Add output format parameters
        cmd.extend(['-o', f'{output_base}.txt'])
        cmd.extend(['-g', f'{output_base}.gff'])
        cmd.extend(['-b', f'{output_base}.bed'])
        
        self.logger.debug(f"Built command: {' '.join(cmd)}")
        return cmd
    
    def _parse_results(self, output_dir: str, output_prefix: str, execution_time: float) -> NLRAnnotatorResult:
        """Parse results from direct output files (fallback method)"""
        genes = []
        output_base = os.path.join(output_dir, output_prefix)
        
        # Try to parse different format files
        for ext in ['.txt', '.gff', '.bed']:
            file_path = f"{output_base}{ext}"
            if os.path.exists(file_path):
                self.logger.info(f"Found output file: {file_path}")
                # For now, just log that we found files
                # The actual parsing would need format-specific logic
        
        return NLRAnnotatorResult(
            genes=genes,
            total_genes=len(genes),
            execution_time=execution_time,
            command="fallback_parsing"
        )
    
    def run_annotation(self, input_fasta: str, output_dir: str, 
                      output_prefix: str = "nlr_results") -> NLRAnnotatorResult:
        """
        Run NLR-Annotator on input FASTA file.
        
        Args:
            input_fasta: Path to input FASTA file
            output_dir: Directory to save results
            output_prefix: Prefix for output files
            
        Returns:
            NLRAnnotatorResult object with analysis results
        """
        # Validate installation
        if not self.validate_installation():
            raise RuntimeError("NLR-Annotator installation validation failed")
            
        # Validate input
        if not self.validate_input(input_fasta):
            raise ValueError(f"Input validation failed for {input_fasta}")
            
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        
        # Build command
        cmd = self._build_command(input_fasta, output_dir, output_prefix)
        
        # Get current working directory to restore later
        original_cwd = os.getcwd()
        
        try:
            # Run command (NLR-Annotator creates temp files in current working directory)
            self.logger.info(f"Running NLR-Annotator command: {' '.join(cmd)}")
            start_time = time.time()
            
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=3600  # 1 hour timeout
            )
            
            execution_time = time.time() - start_time
            self.logger.info(f"NLR-Annotator completed in {execution_time:.2f} seconds")
            
            if result.returncode != 0:
                self.logger.error(f"NLR-Annotator failed with return code {result.returncode}")
                self.logger.error(f"STDERR: {result.stderr}")
                raise RuntimeError(f"NLR-Annotator execution failed: {result.stderr}")
            
            # Look for temporary directories created by NLR-Annotator in current working directory
            temp_dirs = [d for d in os.listdir('.') if d.startswith('temp_mast')]
            
            if not temp_dirs:
                self.logger.warning("No temp_mast directories found, checking for direct output files")
                # Try to parse direct output files if they exist
                return self._parse_results(output_dir, output_prefix, execution_time)
            
            # Process temporary files from temp_mast directories
            all_genes = []
            for temp_dir in temp_dirs:
                temp_path = os.path.join(original_cwd, temp_dir)
                if os.path.isdir(temp_path):
                    self.logger.info(f"Processing temporary directory: {temp_path}")
                    
                    # Look for TSV files in temp directory
                    tsv_files = [f for f in os.listdir(temp_path) if f.endswith('.tsv') and os.path.getsize(os.path.join(temp_path, f)) > 0]
                    
                    for tsv_file in tsv_files:
                        tsv_path = os.path.join(temp_path, tsv_file)
                        self.logger.info(f"Processing TSV file: {tsv_path}")
                        
                        # Parse motifs from TSV file
                        genes = self._parse_temp_tsv(tsv_path)
                        all_genes.extend(genes)
                        
                        # Copy TSV file to output directory with proper name
                        output_tsv = os.path.join(output_dir, f"{output_prefix}_{tsv_file}")
                        shutil.copy2(tsv_path, output_tsv)
                        self.logger.info(f"Copied {tsv_path} to {output_tsv}")
            
            # Group genes by sequence/chromosome and create final results
            grouped_genes = self._group_genes_by_sequence(all_genes)
            
            # Save results in standard formats
            self._save_results(grouped_genes, output_dir, output_prefix)
            
            # Clean up temporary directories
            for temp_dir in temp_dirs:
                temp_path = os.path.join(original_cwd, temp_dir)
                if os.path.isdir(temp_path):
                    shutil.rmtree(temp_path)
                    self.logger.info(f"Cleaned up temporary directory: {temp_path}")
            
            return NLRAnnotatorResult(
                genes=grouped_genes,
                total_genes=len(grouped_genes),
                execution_time=execution_time,
                command=' '.join(cmd)
            )
            
        except Exception as e:
            # Clean up any temp directories on error
            temp_dirs = [d for d in os.listdir('.') if d.startswith('temp_mast')]
            for temp_dir in temp_dirs:
                temp_path = os.path.join(original_cwd, temp_dir)
                if os.path.isdir(temp_path):
                    shutil.rmtree(temp_path)
                    self.logger.info(f"Cleaned up temporary directory on error: {temp_path}")
            raise
    
    def _parse_temp_tsv(self, tsv_path: str) -> List[NLRGene]:
        """Parse motifs from temporary TSV file and group into genes."""
        genes = []
        current_gene_motifs = []
        current_sequence = None
        current_start = None
        current_end = None
        
        try:
            with open(tsv_path, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    
                    parts = line.split('\t')
                    if len(parts) < 8:
                        continue
                    
                    motif_name = parts[0]
                    sequence_frame = parts[1]
                    position = parts[2]
                    motif_sequence = parts[3]
                    score = float(parts[4])
                    
                    # Extract chromosome/sequence name from frame info
                    # Format: Chr1_675000_frame+1 -> Chr1
                    sequence_id = sequence_frame.split('_')[0]
                    
                    # Extract genomic coordinates (parts 6-8)
                    if len(parts) >= 9:
                        chrom = parts[6]
                        start_pos = int(parts[7])
                        end_pos = int(parts[8])
                        strand = parts[9] if len(parts) > 9 else '+'
                    else:
                        chrom = sequence_id
                        start_pos = int(position)
                        end_pos = start_pos + len(motif_sequence) * 3  # Approximate
                        strand = '+'
                    
                    # Group motifs that are close together into genes
                    if (current_sequence != chrom or 
                        current_start is None or 
                        start_pos - current_end > 10000):  # 10kb gap threshold
                        
                        # Save previous gene if exists
                        if current_gene_motifs:
                            gene = self._create_gene_from_motifs(
                                current_gene_motifs, current_sequence, 
                                current_start, current_end
                            )
                            genes.append(gene)
                        
                        # Start new gene
                        current_gene_motifs = []
                        current_sequence = chrom
                        current_start = start_pos
                        current_end = end_pos
                    else:
                        # Extend current gene
                        current_end = max(current_end, end_pos)
                    
                    # Add motif to current gene
                    current_gene_motifs.append({
                        'name': motif_name,
                        'sequence': motif_sequence,
                        'score': score,
                        'start': start_pos,
                        'end': end_pos,
                        'strand': strand
                    })
            
            # Save last gene
            if current_gene_motifs:
                gene = self._create_gene_from_motifs(
                    current_gene_motifs, current_sequence, 
                    current_start, current_end
                )
                genes.append(gene)
                
        except Exception as e:
            self.logger.error(f"Error parsing TSV file {tsv_path}: {e}")
            
        return genes
    
    def _create_gene_from_motifs(self, motifs: List[dict], sequence: str, 
                                start: int, end: int) -> NLRGene:
        """Create NLRGene object from grouped motifs."""
        # Determine gene type based on motifs present
        motif_names = [m['name'] for m in motifs]
        
        # Common NLR motif patterns
        has_cc = any('cc' in name.lower() for name in motif_names)
        has_nbarc = any('nbarc' in name.lower() or 'nb' in name.lower() for name in motif_names)
        has_lrr = any('lrr' in name.lower() for name in motif_names)
        
        # Determine gene type
        if has_cc and has_nbarc and has_lrr:
            gene_type = "CC-NBARC-LRR"
        elif has_nbarc and has_lrr:
            gene_type = "NBARC-LRR"
        elif has_cc and has_nbarc:
            gene_type = "CC-NBARC"
        elif has_nbarc:
            gene_type = "NBARC"
        else:
            gene_type = "NLR-like"
        
        # Calculate average score
        avg_score = sum(m['score'] for m in motifs) / len(motifs) if motifs else 0.0
        
        # Create gene ID
        gene_id = f"{sequence}_{start}_{end}_{gene_type.replace('-', '_')}"
        
        return NLRGene(
            gene_id=gene_id,
            sequence_id=sequence,
            start=start,
            end=end,
            strand=motifs[0]['strand'] if motifs else '+',
            gene_type=gene_type,
            score=avg_score,
            motifs=motifs,
            length=end - start + 1
        )
    
    def _group_genes_by_sequence(self, genes: List[NLRGene]) -> List[NLRGene]:
        """Group and merge overlapping genes by sequence."""
        if not genes:
            return []
        
        # Sort genes by sequence and position
        genes.sort(key=lambda g: (g.sequence_id, g.start))
        
        grouped = []
        current_gene = genes[0]
        
        for next_gene in genes[1:]:
            # If genes are on same sequence and close/overlapping
            if (current_gene.sequence_id == next_gene.sequence_id and
                next_gene.start - current_gene.end <= 5000):  # 5kb merge threshold
                
                # Merge genes
                current_gene = NLRGene(
                    gene_id=f"{current_gene.sequence_id}_{current_gene.start}_{next_gene.end}_merged",
                    sequence_id=current_gene.sequence_id,
                    start=current_gene.start,
                    end=next_gene.end,
                    strand=current_gene.strand,
                    gene_type=f"{current_gene.gene_type}+{next_gene.gene_type}",
                    score=(current_gene.score + next_gene.score) / 2,
                    motifs=current_gene.motifs + next_gene.motifs,
                    length=next_gene.end - current_gene.start + 1
                )
            else:
                # Save current gene and start new one
                grouped.append(current_gene)
                current_gene = next_gene
        
        # Add last gene
        grouped.append(current_gene)
        
        return grouped
    
    def _save_results(self, genes: List[NLRGene], output_dir: str, output_prefix: str):
        """Save results in JSON and TSV formats."""
        # Save as JSON
        json_file = os.path.join(output_dir, f"{output_prefix}.json")
        with open(json_file, 'w') as f:
            json.dump([{
                'gene_id': g.gene_id,
                'sequence_id': g.sequence_id,
                'start': g.start,
                'end': g.end,
                'strand': g.strand,
                'gene_type': g.gene_type,
                'score': g.score,
                'length': g.length,
                'motif_count': len(g.motifs),
                'motifs': g.motifs
            } for g in genes], f, indent=2)
        
        # Save as TSV
        tsv_file = os.path.join(output_dir, f"{output_prefix}.tsv")
        with open(tsv_file, 'w') as f:
            f.write("gene_id\tsequence_id\tstart\tend\tstrand\tgene_type\tscore\tlength\tmotif_count\n")
            for g in genes:
                f.write(f"{g.gene_id}\t{g.sequence_id}\t{g.start}\t{g.end}\t{g.strand}\t{g.gene_type}\t{g.score:.6f}\t{g.length}\t{len(g.motifs)}\n")
        
        # Save as GFF3
        gff_file = os.path.join(output_dir, f"{output_prefix}.gff3")
        with open(gff_file, 'w') as f:
            f.write("##gff-version 3\n")
            for g in genes:
                attributes = f"ID={g.gene_id};gene_type={g.gene_type};score={g.score:.6f};motif_count={len(g.motifs)}"
                f.write(f"{g.sequence_id}\tNLR-Annotator\tgene\t{g.start}\t{g.end}\t{g.score:.3f}\t{g.strand}\t.\t{attributes}\n")
        
        self.logger.info(f"Results saved: {json_file}, {tsv_file}, {gff_file}")


def create_nlr_runner(config_path: Optional[str] = None) -> NLRAnnotatorRunner:
    """
    创建NLRAnnotatorRunner实例的工厂函数
    
    Args:
        config_path: 配置文件路径
        
    Returns:
        NLRAnnotatorRunner实例
    """
    return NLRAnnotatorRunner(config_path) 