"""
NBS Gene Annotation Pipeline
植物NBS抗病基因重新注释流水线

A comprehensive pipeline for re-annotating plant NBS (Nucleotide-Binding Site) 
disease resistance genes using multiple computational genomics tools.

主要功能 / Main Features:
- NLR基因定位 / NLR gene localization
- 基因结构预测 / Gene structure prediction  
- 模型训练 / Model training
- 证据整合 / Evidence integration
"""

__version__ = "0.1.0"
__author__ = "NBS Annotation Team"
__email__ = "nbs@example.com"

# 核心组件导入 / Core component imports
from .pipeline.coordinator import NBSAnnotationPipeline
from .tools.base import ExternalTool, ToolResult
from .data.validation import DataValidator, ValidationResult
from .utils.config import Config, load_config
from .utils.exceptions import (
    NBSAnnotationError,
    ToolExecutionError,
    ValidationError,
    ConfigurationError,
)

__all__ = [
    "NBSAnnotationPipeline",
    "ExternalTool", 
    "ToolResult",
    "DataValidator",
    "ValidationResult", 
    "Config",
    "load_config",
    "NBSAnnotationError",
    "ToolExecutionError",
    "ValidationError", 
    "ConfigurationError",
]