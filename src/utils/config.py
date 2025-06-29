"""
配置管理模块
Configuration Management Module

提供统一的配置文件读取和管理功能。
"""

import os
import yaml
from pathlib import Path
from typing import Any, Dict, Optional, Union
from dataclasses import dataclass


@dataclass
class ProjectPaths:
    """项目路径配置"""
    project_root: Path
    data_dir: Path
    input_dir: Path
    reference_dir: Path
    output_dir: Path
    results_dir: Path
    intermediate_dir: Path
    final_dir: Path
    temp_dir: Path
    log_dir: Path
    tools_dir: Path


@dataclass
class ToolConfig:
    """工具配置"""
    executable: str
    version_required: str
    default_params: Dict[str, Any]


@dataclass
class DatabaseConfig:
    """数据库配置"""
    path: str
    format: str
    description: str


class ConfigManager:
    """配置管理器"""
    
    def __init__(self, config_file: Optional[Union[str, Path]] = None):
        """
        初始化配置管理器
        
        Args:
            config_file: 配置文件路径，默认为项目根目录下的config/config.yaml
        """
        self._config_data = None
        self._project_root = self._find_project_root()
        
        if config_file is None:
            config_file = self._project_root / "config" / "config.yaml"
        
        self.config_file = Path(config_file)
        self._load_config()
        
    def _find_project_root(self) -> Path:
        """
        自动查找项目根目录
        
        Returns:
            项目根目录路径
        """
        current = Path(__file__).parent
        
        # 向上查找，直到找到包含pyproject.toml或setup.py的目录
        while current != current.parent:
            if (current / "pyproject.toml").exists() or (current / "setup.py").exists():
                return current
            current = current.parent
        
        # 如果没找到，返回当前文件的父目录的父目录
        return Path(__file__).parent.parent
    
    def _load_config(self) -> None:
        """加载配置文件"""
        if not self.config_file.exists():
            raise FileNotFoundError(f"配置文件不存在: {self.config_file}")
        
        try:
            with open(self.config_file, 'r', encoding='utf-8') as f:
                self._config_data = yaml.safe_load(f)
        except yaml.YAMLError as e:
            raise ValueError(f"配置文件格式错误: {e}")
        except Exception as e:
            raise RuntimeError(f"读取配置文件失败: {e}")
    
    def get(self, key: str, default: Any = None) -> Any:
        """
        获取配置项
        
        Args:
            key: 配置键，支持点号分隔的嵌套键，如 'tools.miniprot.executable'
            default: 默认值
            
        Returns:
            配置值
        """
        if self._config_data is None:
            return default
        
        keys = key.split('.')
        value = self._config_data
        
        try:
            for k in keys:
                value = value[k]
            return value
        except (KeyError, TypeError):
            return default
    
    def set(self, key: str, value: Any) -> None:
        """
        设置配置项
        
        Args:
            key: 配置键
            value: 配置值
        """
        if self._config_data is None:
            self._config_data = {}
        
        keys = key.split('.')
        current = self._config_data
        
        # 创建嵌套字典结构
        for k in keys[:-1]:
            if k not in current:
                current[k] = {}
            current = current[k]
        
        current[keys[-1]] = value
    
    def save(self, file_path: Optional[Union[str, Path]] = None) -> None:
        """
        保存配置到文件
        
        Args:
            file_path: 保存路径，默认为原配置文件
        """
        if file_path is None:
            file_path = self.config_file
        
        file_path = Path(file_path)
        file_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(file_path, 'w', encoding='utf-8') as f:
            yaml.dump(self._config_data, f, default_flow_style=False, 
                     allow_unicode=True, sort_keys=False)
    
    def get_paths(self) -> ProjectPaths:
        """
        获取项目路径配置
        
        Returns:
            ProjectPaths对象
        """
        # 设置项目根目录
        project_root = self._project_root
        if self.get('paths.project_root'):
            project_root = Path(self.get('paths.project_root'))
        
        return ProjectPaths(
            project_root=project_root,
            data_dir=project_root / self.get('paths.data_dir', 'data'),
            input_dir=project_root / self.get('paths.input_dir', 'data/input'),
            reference_dir=project_root / self.get('paths.reference_dir', 'data/reference'),
            output_dir=project_root / self.get('paths.output_dir', 'output'),
            results_dir=project_root / self.get('paths.results_dir', 'results'),
            intermediate_dir=project_root / self.get('paths.intermediate_dir', 'results/intermediate'),
            final_dir=project_root / self.get('paths.final_dir', 'results/final'),
            temp_dir=project_root / self.get('paths.temp_dir', 'temp'),
            log_dir=project_root / self.get('paths.log_dir', 'logs'),
            tools_dir=project_root / self.get('paths.tools_dir', 'tools')
        )
    
    def get_tool_config(self, tool_name: str) -> Optional[ToolConfig]:
        """
        获取工具配置
        
        Args:
            tool_name: 工具名称
            
        Returns:
            ToolConfig对象或None
        """
        tool_data = self.get(f'tools.{tool_name}')
        if not tool_data:
            return None
        
        # 处理可执行文件路径
        executable = tool_data.get('executable', '')
        if not os.path.isabs(executable):
            executable = str(self._project_root / executable)
        
        return ToolConfig(
            executable=executable,
            version_required=tool_data.get('version_required', ''),
            default_params=tool_data.get('default_params', {})
        )
    
    def get_database_config(self, db_name: str) -> Optional[DatabaseConfig]:
        """
        获取数据库配置
        
        Args:
            db_name: 数据库名称
            
        Returns:
            DatabaseConfig对象或None
        """
        db_data = self.get(f'databases.{db_name}')
        if not db_data:
            return None
        
        # 处理数据库文件路径
        db_path = db_data.get('path', '')
        if not os.path.isabs(db_path):
            db_path = str(self._project_root / db_path)
        
        return DatabaseConfig(
            path=db_path,
            format=db_data.get('format', ''),
            description=db_data.get('description', '')
        )
    
    def get_analysis_params(self, analysis_type: str) -> Dict[str, Any]:
        """
        获取分析参数
        
        Args:
            analysis_type: 分析类型 (localization, prediction, training, integration)
            
        Returns:
            分析参数字典
        """
        return self.get(f'analysis.{analysis_type}', {})
    
    def get_logging_config(self) -> Dict[str, Any]:
        """
        获取日志配置
        
        Returns:
            日志配置字典
        """
        return self.get('logging', {
            'level': 'INFO',
            'console_output': True,
            'file_output': True,
            'max_file_size_mb': 10,
            'backup_count': 5
        })
    
    def get_performance_config(self) -> Dict[str, Any]:
        """
        获取性能配置
        
        Returns:
            性能配置字典
        """
        return self.get('performance', {
            'parallel': {'max_processes': 4, 'chunk_size': 1000},
            'memory': {'max_memory_gb': 16, 'temp_cleanup': True}
        })
    
    def validate_config(self) -> Dict[str, Any]:
        """
        验证配置文件的完整性
        
        Returns:
            验证结果字典，包含错误和警告信息
        """
        errors = []
        warnings = []
        
        # 检查必需的配置项
        required_sections = ['project', 'paths', 'tools', 'databases']
        for section in required_sections:
            if not self.get(section):
                errors.append(f"缺少必需的配置节: {section}")
        
        # 检查工具可执行文件
        tools = ['nlgenomesweeper', 'miniprot', 'augustus', 'evm']
        for tool in tools:
            tool_config = self.get_tool_config(tool)
            if tool_config:
                if not Path(tool_config.executable).exists():
                    warnings.append(f"工具可执行文件不存在: {tool_config.executable}")
        
        # 检查数据库文件
        databases = ['resistance_genes', 'reference_genome']
        for db in databases:
            db_config = self.get_database_config(db)
            if db_config:
                if not Path(db_config.path).exists():
                    warnings.append(f"数据库文件不存在: {db_config.path}")
        
        return {
            'valid': len(errors) == 0,
            'errors': errors,
            'warnings': warnings
        }
    
    def create_directories(self) -> None:
        """创建项目所需的目录"""
        paths = self.get_paths()
        
        directories = [
            paths.data_dir, paths.input_dir, paths.reference_dir,
            paths.output_dir, paths.results_dir, paths.intermediate_dir,
            paths.final_dir, paths.temp_dir, paths.log_dir
        ]
        
        for directory in directories:
            directory.mkdir(parents=True, exist_ok=True)
    
    @property
    def project_root(self) -> Path:
        """获取项目根目录"""
        return self._project_root
    
    @property
    def config_data(self) -> Dict[str, Any]:
        """获取完整配置数据"""
        return self._config_data or {}


# 全局配置管理器实例
_global_config = None


def get_config(config_file: Optional[Union[str, Path]] = None) -> ConfigManager:
    """
    获取全局配置管理器实例
    
    Args:
        config_file: 配置文件路径
        
    Returns:
        ConfigManager实例
    """
    global _global_config
    
    if _global_config is None or config_file is not None:
        _global_config = ConfigManager(config_file)
    
    return _global_config


def reload_config(config_file: Optional[Union[str, Path]] = None) -> ConfigManager:
    """
    重新加载配置
    
    Args:
        config_file: 配置文件路径
        
    Returns:
        新的ConfigManager实例
    """
    global _global_config
    _global_config = ConfigManager(config_file)
    return _global_config


if __name__ == "__main__":
    # 测试配置管理器
    try:
        config = ConfigManager()
        
        print("=== 项目信息 ===")
        print(f"项目名称: {config.get('project.name')}")
        print(f"项目版本: {config.get('project.version')}")
        print(f"项目根目录: {config.project_root}")
        
        print("\n=== 路径配置 ===")
        paths = config.get_paths()
        print(f"数据目录: {paths.data_dir}")
        print(f"输出目录: {paths.output_dir}")
        print(f"日志目录: {paths.log_dir}")
        
        print("\n=== 工具配置 ===")
        miniprot_config = config.get_tool_config('miniprot')
        if miniprot_config:
            print(f"miniprot: {miniprot_config.executable}")
            print(f"参数: {miniprot_config.default_params}")
        
        print("\n=== 配置验证 ===")
        validation = config.validate_config()
        print(f"配置有效: {validation['valid']}")
        if validation['errors']:
            print(f"错误: {validation['errors']}")
        if validation['warnings']:
            print(f"警告: {validation['warnings']}")
        
        print("\n配置管理器测试完成")
        
    except Exception as e:
        print(f"配置管理器测试失败: {e}") 