"""
Data validation module for NBS annotation pipeline
NBS注释流水线的数据验证模块

提供输入文件格式验证和质量检查功能
"""

import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union
from dataclasses import dataclass
from collections import Counter

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from ..utils.logging_setup import get_logger
from ..utils.exceptions import ValidationError

logger = get_logger(__name__)


@dataclass
class ValidationResult:
    """
    Container for file validation results
    文件验证结果容器
    """
    
    is_valid: bool
    file_path: Path
    file_type: str
    file_size: int
    record_count: int
    error_message: Optional[str] = None
    warning_messages: List[str] = None
    statistics: Dict[str, Any] = None
    
    def __post_init__(self) -> None:
        """Initialize default values"""
        if self.warning_messages is None:
            self.warning_messages = []
        if self.statistics is None:
            self.statistics = {}


class DataValidator:
    """
    Comprehensive data validator for bioinformatics files
    生物信息学文件的综合数据验证器
    
    支持验证：
    - FASTA格式文件（基因组、蛋白质序列）
    - GFF/GTF格式文件
    - 文件完整性和质量检查
    """
    
    def __init__(self, strict_mode: bool = False) -> None:
        """
        Initialize data validator
        
        Args:
            strict_mode: Enable strict validation mode / 启用严格验证模式
        """
        self.strict_mode = strict_mode
        self.logger = get_logger(f"{__name__}.DataValidator")
        
        # Validation criteria
        self.min_sequence_length = 100  # Minimum sequence length
        self.max_sequence_length = 1e9   # Maximum sequence length (1Gb)
        self.max_file_size = 50 * 1024**3  # Maximum file size (50GB)
        
    def validate_file_existence(self, file_path: Union[str, Path]) -> Path:
        """
        Validate that file exists and is readable
        验证文件存在且可读
        
        Args:
            file_path: Path to file to validate
            
        Returns:
            Validated Path object
            
        Raises:
            ValidationError: If file doesn't exist or isn't readable
        """
        path = Path(file_path)
        
        if not path.exists():
            raise ValidationError(
                f"File does not exist: {path}",
                file_path=str(path),
                validation_type="existence",
            )
            
        if not path.is_file():
            raise ValidationError(
                f"Path is not a file: {path}",
                file_path=str(path),
                validation_type="file_type",
            )
            
        if not path.stat().st_size > 0:
            raise ValidationError(
                f"File is empty: {path}",
                file_path=str(path),
                validation_type="file_size",
            )
            
        # Check file size
        file_size = path.stat().st_size
        if file_size > self.max_file_size:
            raise ValidationError(
                f"File too large ({file_size / 1024**3:.1f}GB > {self.max_file_size / 1024**3:.1f}GB): {path}",
                file_path=str(path),
                validation_type="file_size",
            )
            
        return path
    
    def validate_fasta(
        self,
        file_path: Union[str, Path],
        sequence_type: str = "nucleotide",
    ) -> ValidationResult:
        """
        Validate FASTA format file
        验证FASTA格式文件
        
        Args:
            file_path: Path to FASTA file
            sequence_type: Type of sequences (nucleotide/protein)
            
        Returns:
            ValidationResult object
        """
        path = self.validate_file_existence(file_path)
        file_size = path.stat().st_size
        
        self.logger.info(f"Validating FASTA file: {path}")
        
        try:
            # Parse FASTA file
            records = list(SeqIO.parse(str(path), "fasta"))
            
            if not records:
                return ValidationResult(
                    is_valid=False,
                    file_path=path,
                    file_type="fasta",
                    file_size=file_size,
                    record_count=0,
                    error_message="No sequences found in FASTA file",
                )
            
            # Collect statistics
            sequence_lengths = [len(record.seq) for record in records]
            sequence_ids = [record.id for record in records]
            
            # Check for duplicate IDs
            id_counts = Counter(sequence_ids)
            duplicate_ids = [seq_id for seq_id, count in id_counts.items() if count > 1]
            
            # Validate sequence content
            warnings = []
            valid_sequences = 0
            
            for i, record in enumerate(records):
                seq_str = str(record.seq).upper()
                seq_len = len(seq_str)
                
                # Check sequence length
                if seq_len < self.min_sequence_length:
                    warnings.append(
                        f"Sequence {record.id} is very short ({seq_len} bp)"
                    )
                elif seq_len > self.max_sequence_length:
                    return ValidationResult(
                        is_valid=False,
                        file_path=path,
                        file_type="fasta",
                        file_size=file_size,
                        record_count=len(records),
                        error_message=f"Sequence {record.id} is too long ({seq_len} bp)",
                    )
                
                # Validate sequence composition
                if sequence_type == "nucleotide":
                    # Check for valid nucleotide characters
                    valid_chars = set("ATCGNRYSWKMBDHV-.")
                    invalid_chars = set(seq_str) - valid_chars
                    
                    if invalid_chars:
                        if self.strict_mode:
                            return ValidationResult(
                                is_valid=False,
                                file_path=path,
                                file_type="fasta",
                                file_size=file_size,
                                record_count=len(records),
                                error_message=f"Invalid nucleotide characters in {record.id}: {invalid_chars}",
                            )
                        else:
                            warnings.append(
                                f"Invalid nucleotide characters in {record.id}: {invalid_chars}"
                            )
                    
                    # Check N content
                    n_content = seq_str.count('N') / seq_len
                    if n_content > 0.5:
                        warnings.append(
                            f"High N content ({n_content:.1%}) in sequence {record.id}"
                        )
                        
                elif sequence_type == "protein":
                    # Check for valid amino acid characters
                    valid_chars = set("ACDEFGHIKLMNPQRSTVWYXZBU*-.")
                    invalid_chars = set(seq_str) - valid_chars
                    
                    if invalid_chars:
                        if self.strict_mode:
                            return ValidationResult(
                                is_valid=False,
                                file_path=path,
                                file_type="fasta",
                                file_size=file_size,
                                record_count=len(records),
                                error_message=f"Invalid amino acid characters in {record.id}: {invalid_chars}",
                            )
                        else:
                            warnings.append(
                                f"Invalid amino acid characters in {record.id}: {invalid_chars}"
                            )
                    
                    # Check X content (unknown amino acids)
                    x_content = seq_str.count('X') / seq_len
                    if x_content > 0.1:
                        warnings.append(
                            f"High X content ({x_content:.1%}) in sequence {record.id}"
                        )
                
                valid_sequences += 1
            
            # Additional validations
            if duplicate_ids:
                if self.strict_mode:
                    return ValidationResult(
                        is_valid=False,
                        file_path=path,
                        file_type="fasta",
                        file_size=file_size,
                        record_count=len(records),
                        error_message=f"Duplicate sequence IDs found: {duplicate_ids[:5]}",
                    )
                else:
                    warnings.append(f"Duplicate sequence IDs found: {len(duplicate_ids)} duplicates")
            
            # Compile statistics
            statistics = {
                "total_sequences": len(records),
                "valid_sequences": valid_sequences,
                "total_length": sum(sequence_lengths),
                "min_length": min(sequence_lengths),
                "max_length": max(sequence_lengths),
                "mean_length": sum(sequence_lengths) / len(sequence_lengths),
                "duplicate_ids": len(duplicate_ids),
                "sequence_type": sequence_type,
            }
            
            self.logger.info(f"FASTA validation completed: {len(records)} sequences, {statistics['total_length']} bp total")
            
            return ValidationResult(
                is_valid=True,
                file_path=path,
                file_type="fasta",
                file_size=file_size,
                record_count=len(records),
                warning_messages=warnings,
                statistics=statistics,
            )
            
        except Exception as e:
            error_msg = f"FASTA parsing error: {e}"
            self.logger.error(error_msg)
            return ValidationResult(
                is_valid=False,
                file_path=path,
                file_type="fasta",
                file_size=file_size,
                record_count=0,
                error_message=error_msg,
            )
    
    def validate_protein_fasta(self, file_path: Union[str, Path]) -> ValidationResult:
        """
        Validate protein FASTA file
        验证蛋白质FASTA文件
        
        Args:
            file_path: Path to protein FASTA file
            
        Returns:
            ValidationResult object
        """
        return self.validate_fasta(file_path, sequence_type="protein")
    
    def validate_genome_fasta(self, file_path: Union[str, Path]) -> ValidationResult:
        """
        Validate genome FASTA file
        验证基因组FASTA文件
        
        Args:
            file_path: Path to genome FASTA file
            
        Returns:
            ValidationResult object
        """
        return self.validate_fasta(file_path, sequence_type="nucleotide")
    
    def validate_gff(self, file_path: Union[str, Path]) -> ValidationResult:
        """
        Validate GFF/GTF format file
        验证GFF/GTF格式文件
        
        Args:
            file_path: Path to GFF file
            
        Returns:
            ValidationResult object
        """
        path = self.validate_file_existence(file_path)
        file_size = path.stat().st_size
        
        self.logger.info(f"Validating GFF file: {path}")
        
        try:
            records = []
            warnings = []
            line_count = 0
            feature_types = Counter()
            
            with open(path, 'r', encoding='utf-8') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    line_count += 1
                    
                    # Skip empty lines and comments
                    if not line or line.startswith('#'):
                        continue
                    
                    # Parse GFF line
                    fields = line.split('\t')
                    
                    if len(fields) != 9:
                        if self.strict_mode:
                            return ValidationResult(
                                is_valid=False,
                                file_path=path,
                                file_type="gff",
                                file_size=file_size,
                                record_count=0,
                                error_message=f"Invalid GFF format at line {line_num}: expected 9 fields, got {len(fields)}",
                            )
                        else:
                            warnings.append(f"Invalid GFF format at line {line_num}")
                            continue
                    
                    seqid, source, feature_type, start, end, score, strand, phase, attributes = fields
                    
                    # Validate coordinates
                    try:
                        start_pos = int(start)
                        end_pos = int(end)
                        
                        if start_pos > end_pos:
                            warnings.append(f"Invalid coordinates at line {line_num}: start > end")
                        if start_pos < 1:
                            warnings.append(f"Invalid start coordinate at line {line_num}: {start_pos}")
                            
                    except ValueError:
                        warnings.append(f"Invalid coordinates at line {line_num}: {start}-{end}")
                    
                    # Validate strand
                    if strand not in ['+', '-', '.', '?']:
                        warnings.append(f"Invalid strand at line {line_num}: {strand}")
                    
                    # Count feature types
                    feature_types[feature_type] += 1
                    records.append(fields)
            
            statistics = {
                "total_lines": line_count,
                "total_features": len(records),
                "feature_types": dict(feature_types),
                "most_common_features": feature_types.most_common(10),
            }
            
            self.logger.info(f"GFF validation completed: {len(records)} features")
            
            return ValidationResult(
                is_valid=True,
                file_path=path,
                file_type="gff",
                file_size=file_size,
                record_count=len(records),
                warning_messages=warnings,
                statistics=statistics,
            )
            
        except Exception as e:
            error_msg = f"GFF parsing error: {e}"
            self.logger.error(error_msg)
            return ValidationResult(
                is_valid=False,
                file_path=path,
                file_type="gff",
                file_size=file_size,
                record_count=0,
                error_message=error_msg,
            )
    
    def validate_multiple_files(
        self,
        file_paths: Dict[str, Path],
        file_types: Dict[str, str],
    ) -> Dict[str, ValidationResult]:
        """
        Validate multiple files
        验证多个文件
        
        Args:
            file_paths: Dictionary mapping file names to paths
            file_types: Dictionary mapping file names to types
            
        Returns:
            Dictionary of validation results
        """
        results = {}
        
        for file_name, file_path in file_paths.items():
            file_type = file_types.get(file_name, "unknown")
            
            try:
                if file_type == "genome_fasta":
                    results[file_name] = self.validate_genome_fasta(file_path)
                elif file_type == "protein_fasta":
                    results[file_name] = self.validate_protein_fasta(file_path)
                elif file_type == "gff":
                    results[file_name] = self.validate_gff(file_path)
                else:
                    # Generic file existence check
                    self.validate_file_existence(file_path)
                    results[file_name] = ValidationResult(
                        is_valid=True,
                        file_path=file_path,
                        file_type=file_type,
                        file_size=file_path.stat().st_size,
                        record_count=0,
                        warning_messages=[f"Unknown file type: {file_type}"],
                    )
                    
            except Exception as e:
                results[file_name] = ValidationResult(
                    is_valid=False,
                    file_path=file_path,
                    file_type=file_type,
                    file_size=0,
                    record_count=0,
                    error_message=str(e),
                )
        
        return results
    
    def generate_validation_report(
        self,
        validation_results: Dict[str, ValidationResult],
        output_path: Optional[Path] = None,
    ) -> str:
        """
        Generate comprehensive validation report
        生成综合验证报告
        
        Args:
            validation_results: Dictionary of validation results
            output_path: Optional path to save report
            
        Returns:
            Validation report as string
        """
        report_lines = [
            "# NBS Annotation Pipeline - Data Validation Report",
            f"Generated: {__import__('datetime').datetime.now().isoformat()}",
            f"Total files validated: {len(validation_results)}",
            "",
        ]
        
        # Summary
        valid_files = sum(1 for r in validation_results.values() if r.is_valid)
        invalid_files = len(validation_results) - valid_files
        
        report_lines.extend([
            "## Summary",
            f"- Valid files: {valid_files}",
            f"- Invalid files: {invalid_files}",
            f"- Validation success rate: {valid_files / len(validation_results) * 100:.1f}%",
            "",
        ])
        
        # File details
        report_lines.append("## File Details")
        
        for file_name, result in validation_results.items():
            status = "✅ VALID" if result.is_valid else "❌ INVALID"
            report_lines.extend([
                f"### {file_name} - {status}",
                f"- File path: `{result.file_path}`",
                f"- File type: {result.file_type}",
                f"- File size: {result.file_size / 1024**2:.1f} MB",
                f"- Record count: {result.record_count}",
            ])
            
            if result.error_message:
                report_lines.append(f"- **Error**: {result.error_message}")
            
            if result.warning_messages:
                report_lines.append("- **Warnings**:")
                for warning in result.warning_messages:
                    report_lines.append(f"  - {warning}")
            
            if result.statistics:
                report_lines.append("- **Statistics**:")
                for key, value in result.statistics.items():
                    if isinstance(value, dict):
                        report_lines.append(f"  - {key}: {len(value)} items")
                    else:
                        report_lines.append(f"  - {key}: {value}")
            
            report_lines.append("")
        
        report_text = "\n".join(report_lines)
        
        if output_path:
            with open(output_path, 'w', encoding='utf-8') as f:
                f.write(report_text)
            self.logger.info(f"Validation report saved to: {output_path}")
        
        return report_text