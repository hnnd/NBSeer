"""
文件I/O和格式验证模块
File I/O and Format Validation Module

提供文件读取、格式验证和数据处理功能，专门针对生物信息学数据格式。
"""

import os
import gzip
import hashlib
from pathlib import Path
from typing import Iterator, Dict, List, Optional, Tuple, Any, Union
from dataclasses import dataclass
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import logging

from .logging_setup import get_logger


@dataclass
class FileInfo:
    """文件信息"""
    path: Path
    size_bytes: int
    exists: bool
    readable: bool
    format_valid: bool
    error_message: Optional[str] = None


@dataclass
class SequenceStats:
    """序列统计信息"""
    total_sequences: int
    total_length: int
    min_length: int
    max_length: int
    avg_length: float
    gc_content: float
    n_content: float
    sequence_ids: List[str]


class FastaValidator:
    """FASTA文件验证器"""
    
    def __init__(self, logger: Optional[logging.Logger] = None):
        self.logger = logger or get_logger(__name__)
        
    def validate_file_format(self, file_path: Union[str, Path]) -> FileInfo:
        """
        验证文件基本信息和可读性
        
        Args:
            file_path: 文件路径
            
        Returns:
            FileInfo对象，包含文件基本信息
        """
        file_path = Path(file_path)
        
        # 检查文件是否存在
        if not file_path.exists():
            return FileInfo(
                path=file_path,
                size_bytes=0,
                exists=False,
                readable=False,
                format_valid=False,
                error_message=f"文件不存在: {file_path}"
            )
        
        # 获取文件大小
        try:
            size_bytes = file_path.stat().st_size
        except OSError as e:
            return FileInfo(
                path=file_path,
                size_bytes=0,
                exists=True,
                readable=False,
                format_valid=False,
                error_message=f"无法获取文件信息: {e}"
            )
        
        # 检查文件是否可读
        try:
            with open(file_path, 'r') as f:
                f.read(1)  # 尝试读取一个字符
            readable = True
        except (OSError, UnicodeDecodeError):
            # 可能是二进制文件或权限问题
            try:
                with open(file_path, 'rb') as f:
                    f.read(1)
                readable = True
            except OSError as e:
                return FileInfo(
                    path=file_path,
                    size_bytes=size_bytes,
                    exists=True,
                    readable=False,
                    format_valid=False,
                    error_message=f"文件不可读: {e}"
                )
        
        return FileInfo(
            path=file_path,
            size_bytes=size_bytes,
            exists=True,
            readable=readable,
            format_valid=True  # 这里只是基本检查，具体格式验证在后面
        )
    
    def validate_fasta_format(self, file_path: Union[str, Path]) -> Tuple[bool, Optional[str], int]:
        """
        验证FASTA文件格式
        
        Args:
            file_path: FASTA文件路径
            
        Returns:
            (是否有效, 错误信息, 序列数量)
        """
        file_path = Path(file_path)
        self.logger.info(f"验证FASTA文件格式: {file_path}")
        
        # 首先进行基本文件检查
        file_info = self.validate_file_format(file_path)
        if not file_info.format_valid:
            return False, file_info.error_message, 0
        
        try:
            # 检查文件是否为压缩格式
            is_gzipped = file_path.suffix.lower() == '.gz'
            
            # 尝试解析FASTA文件
            sequence_count = 0
            
            if is_gzipped:
                with gzip.open(file_path, 'rt') as handle:
                    sequences = SeqIO.parse(handle, "fasta")
                    for seq_record in sequences:
                        sequence_count += 1
                        # 检查前几个序列的基本格式
                        if sequence_count <= 5:
                            if not self._validate_sequence_record(seq_record):
                                return False, f"序列格式错误: {seq_record.id}", sequence_count
            else:
                with open(file_path, 'r') as handle:
                    sequences = SeqIO.parse(handle, "fasta")
                    for seq_record in sequences:
                        sequence_count += 1
                        # 检查前几个序列的基本格式
                        if sequence_count <= 5:
                            if not self._validate_sequence_record(seq_record):
                                return False, f"序列格式错误: {seq_record.id}", sequence_count
            
            if sequence_count == 0:
                return False, "文件中未找到有效的FASTA序列", 0
            
            self.logger.info(f"FASTA文件验证成功，包含 {sequence_count} 个序列")
            return True, None, sequence_count
            
        except Exception as e:
            error_msg = f"FASTA文件格式验证失败: {str(e)}"
            self.logger.error(error_msg)
            return False, error_msg, 0
    
    def _validate_sequence_record(self, seq_record: SeqRecord) -> bool:
        """
        验证单个序列记录的格式
        
        Args:
            seq_record: SeqRecord对象
            
        Returns:
            是否有效
        """
        # 检查序列ID是否存在
        if not seq_record.id:
            return False
        
        # 检查序列长度
        if len(seq_record.seq) == 0:
            return False
        
        # 检查序列是否包含有效字符（DNA/RNA/蛋白质）
        sequence_str = str(seq_record.seq).upper()
        
        # DNA/RNA序列检查
        dna_chars = set('ATCGRYSWKMBDHVN-')
        protein_chars = set('ACDEFGHIKLMNPQRSTVWYXBZJUO*-')
        
        seq_chars = set(sequence_str)
        
        # 如果序列字符都在DNA字符集中，或都在蛋白质字符集中，则认为有效
        if seq_chars.issubset(dna_chars) or seq_chars.issubset(protein_chars):
            return True
        
        # 如果包含非标准字符，记录警告但不认为是错误
        invalid_chars = seq_chars - dna_chars - protein_chars
        if invalid_chars:
            self.logger.warning(f"序列 {seq_record.id} 包含非标准字符: {invalid_chars}")
        
        return True
    
    def get_sequence_statistics(self, file_path: Union[str, Path]) -> Optional[SequenceStats]:
        """
        获取FASTA文件的序列统计信息
        
        Args:
            file_path: FASTA文件路径
            
        Returns:
            SequenceStats对象或None
        """
        file_path = Path(file_path)
        self.logger.info(f"计算序列统计信息: {file_path}")
        
        try:
            sequences = []
            sequence_ids = []
            total_length = 0
            total_gc = 0
            total_n = 0
            
            # 检查文件是否为压缩格式
            is_gzipped = file_path.suffix.lower() == '.gz'
            
            if is_gzipped:
                with gzip.open(file_path, 'rt') as handle:
                    for seq_record in SeqIO.parse(handle, "fasta"):
                        seq_len = len(seq_record.seq)
                        sequences.append(seq_len)
                        sequence_ids.append(seq_record.id)
                        total_length += seq_len
                        
                        # 计算GC和N含量
                        seq_str = str(seq_record.seq).upper()
                        gc_count = seq_str.count('G') + seq_str.count('C')
                        n_count = seq_str.count('N')
                        
                        total_gc += gc_count
                        total_n += n_count
            else:
                with open(file_path, 'r') as handle:
                    for seq_record in SeqIO.parse(handle, "fasta"):
                        seq_len = len(seq_record.seq)
                        sequences.append(seq_len)
                        sequence_ids.append(seq_record.id)
                        total_length += seq_len
                        
                        # 计算GC和N含量
                        seq_str = str(seq_record.seq).upper()
                        gc_count = seq_str.count('G') + seq_str.count('C')
                        n_count = seq_str.count('N')
                        
                        total_gc += gc_count
                        total_n += n_count
            
            if not sequences:
                return None
            
            stats = SequenceStats(
                total_sequences=len(sequences),
                total_length=total_length,
                min_length=min(sequences),
                max_length=max(sequences),
                avg_length=total_length / len(sequences),
                gc_content=total_gc / total_length if total_length > 0 else 0.0,
                n_content=total_n / total_length if total_length > 0 else 0.0,
                sequence_ids=sequence_ids
            )
            
            self.logger.info(f"统计完成: {stats.total_sequences} 个序列，总长度 {stats.total_length:,} bp")
            return stats
            
        except Exception as e:
            self.logger.error(f"计算序列统计信息失败: {e}")
            return None


class FastaReader:
    """FASTA文件读取器"""
    
    def __init__(self, logger: Optional[logging.Logger] = None):
        self.logger = logger or get_logger(__name__)
    
    def read_sequences(self, file_path: Union[str, Path], format: str = "fasta") -> Iterator[SeqRecord]:
        """
        流式读取序列文件
        
        Args:
            file_path: 文件路径
            format: 文件格式，默认为"fasta"
            
        Yields:
            SeqRecord对象
        """
        file_path = Path(file_path)
        
        try:
            # 检查文件是否为压缩格式
            is_gzipped = file_path.suffix.lower() == '.gz'
            
            if is_gzipped:
                with gzip.open(file_path, 'rt') as handle:
                    for seq_record in SeqIO.parse(handle, format):
                        yield seq_record
            else:
                with open(file_path, 'r') as handle:
                    for seq_record in SeqIO.parse(handle, format):
                        yield seq_record
                        
        except Exception as e:
            self.logger.error(f"读取序列文件失败: {e}")
            raise
    
    def read_sequence_by_id(self, file_path: Union[str, Path], seq_id: str, 
                           format: str = "fasta") -> Optional[SeqRecord]:
        """
        根据ID读取特定序列
        
        Args:
            file_path: 文件路径
            seq_id: 序列ID
            format: 文件格式
            
        Returns:
            SeqRecord对象或None
        """
        try:
            for seq_record in self.read_sequences(file_path, format):
                if seq_record.id == seq_id:
                    return seq_record
            return None
        except Exception as e:
            self.logger.error(f"根据ID读取序列失败: {e}")
            return None
    
    def create_sequence_index(self, file_path: Union[str, Path], 
                            format: str = "fasta") -> Optional[Dict[str, Any]]:
        """
        为FASTA文件创建索引
        
        Args:
            file_path: 文件路径
            format: 文件格式
            
        Returns:
            索引字典或None
        """
        file_path = Path(file_path)
        
        try:
            # 使用Biopython的索引功能
            if file_path.suffix.lower() == '.gz':
                # 对于压缩文件，我们需要先解压或使用其他方法
                self.logger.warning("压缩文件暂不支持索引，将使用流式读取")
                return None
            else:
                index = SeqIO.index(str(file_path), format)
                self.logger.info(f"索引创建成功，包含 {len(index)} 个序列")
                return index
                
        except Exception as e:
            self.logger.error(f"创建序列索引失败: {e}")
            return None


def calculate_file_checksum(file_path: Union[str, Path], algorithm: str = "md5") -> str:
    """
    计算文件校验和
    
    Args:
        file_path: 文件路径
        algorithm: 校验算法（md5, sha1, sha256）
        
    Returns:
        校验和字符串
    """
    file_path = Path(file_path)
    
    if algorithm.lower() == "md5":
        hash_obj = hashlib.md5()
    elif algorithm.lower() == "sha1":
        hash_obj = hashlib.sha1()
    elif algorithm.lower() == "sha256":
        hash_obj = hashlib.sha256()
    else:
        raise ValueError(f"不支持的校验算法: {algorithm}")
    
    with open(file_path, 'rb') as f:
        for chunk in iter(lambda: f.read(8192), b""):
            hash_obj.update(chunk)
    
    return hash_obj.hexdigest()


if __name__ == "__main__":
    # 测试FASTA验证器
    logger = get_logger(__name__)
    validator = FastaValidator(logger)
    
    # 测试现有的基因组文件
    test_files = [
        "genome/osa.fa",
        "db/AllResistanceGenes.fasta"
    ]
    
    for test_file in test_files:
        if Path(test_file).exists():
            print(f"\n=== 测试文件: {test_file} ===")
            
            # 验证格式
            is_valid, error_msg, seq_count = validator.validate_fasta_format(test_file)
            print(f"格式有效: {is_valid}")
            if error_msg:
                print(f"错误信息: {error_msg}")
            print(f"序列数量: {seq_count}")
            
            # 获取统计信息
            if is_valid:
                stats = validator.get_sequence_statistics(test_file)
                if stats:
                    print(f"总序列数: {stats.total_sequences}")
                    print(f"总长度: {stats.total_length:,} bp")
                    print(f"平均长度: {stats.avg_length:.1f} bp")
                    print(f"GC含量: {stats.gc_content:.2%}")
                    print(f"N含量: {stats.n_content:.2%}")
        else:
            print(f"文件不存在: {test_file}")
    
    print("\nFASTA验证器测试完成") 