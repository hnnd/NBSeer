"""
Gene ID renaming configuration for post-EVM processing.
Provides user-defined prefix-based naming conventions.
"""

from typing import Dict, List, Optional, Tuple
import re


class GeneNamingConfig:
    """Configuration for gene ID generation with user-defined prefixes."""
    
    DEFAULT_PREFIX = "NBS"  # 默认前缀
    
    def __init__(self, prefix: str = None):
        """
        Initialize gene naming configuration.
        
        Args:
            prefix: User-defined prefix for gene IDs (default: "NBS")
        """
        self.prefix = prefix or self.DEFAULT_PREFIX
    
    def get_prefix(self) -> str:
        """Get the configured prefix."""
        return self.prefix
    
    def generate_gene_id(self, gene_number: int, feature_type: str = "gene") -> str:
        """
        Generate gene ID following the pattern: {prefix}_{number}.{suffix}
        
        Args:
            gene_number: Sequential gene number
            feature_type: Feature type (gene, mRNA, etc.)
            
        Returns:
            Formatted gene ID (e.g., "OsNP_0001.g", "OsNP_0001.t")
        """
        base_id = f"{self.prefix}_{gene_number:04d}"
        
        # 根据特征类型添加后缀
        if feature_type.lower() in ['gene']:
            return f"{base_id}.g"
        elif feature_type.lower() in ['mrna', 'transcript']:
            return f"{base_id}.t"
        else:
            # 对于其他特征类型（如exon, CDS），返回基础ID
            return base_id
    
    def validate_prefix(self, prefix: str) -> bool:
        """
        Validate prefix format.
        
        Args:
            prefix: User-provided prefix
            
        Returns:
            True if prefix is valid, False otherwise
        """
        # 基本验证：不能为空，不能包含特殊字符
        if not prefix:
            return False
        
        # 允许字母、数字和下划线
        return re.match(r'^[A-Za-z0-9_]+$', prefix) is not None
    
    @classmethod
    def create_from_prefix(cls, prefix: str) -> 'GeneNamingConfig':
        """
        Create configuration from user-provided prefix.
        
        Args:
            prefix: User-defined prefix
            
        Returns:
            GeneNamingConfig instance
        """
        config = cls(prefix)
        if not config.validate_prefix(prefix):
            raise ValueError(f"Invalid prefix: {prefix}. Only alphanumeric characters and underscores are allowed.")
        return config
    
    def __str__(self) -> str:
        return f"GeneNamingConfig(prefix='{self.prefix}')"
    
    def __repr__(self) -> str:
        return self.__str__()