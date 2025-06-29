"""
External tool interfaces for NBS annotation pipeline
NBS注释流水线的外部工具接口

提供统一的外部工具调用接口
"""

from .base import ExternalTool, ToolResult

__all__ = [
    "ExternalTool",
    "ToolResult",
]