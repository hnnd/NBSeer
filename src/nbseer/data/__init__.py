"""
Data processing and validation module for NBS annotation
NBS注释的数据处理和验证模块

提供数据验证、格式转换和质量控制功能
"""

from .validation import DataValidator, ValidationResult

__all__ = [
    "DataValidator", 
    "ValidationResult",
]