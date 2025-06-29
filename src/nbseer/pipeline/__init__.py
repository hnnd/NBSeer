"""
Pipeline coordination module for NBS annotation
NBS注释流水线协调模块

提供流水线执行、管理和协调功能
"""

from .coordinator import NBSAnnotationPipeline, PipelineResults

__all__ = [
    "NBSAnnotationPipeline",
    "PipelineResults",
]