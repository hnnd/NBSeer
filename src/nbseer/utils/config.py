"""
Configuration management for NBS annotation pipeline
NBS注释流水线的配置管理模块

提供配置文件加载、验证和管理功能
"""

import os
from pathlib import Path
from typing import Any, Dict, List, Optional, Union
from dataclasses import dataclass
import yaml

from .exceptions import ConfigurationError
from .logging_setup import get_logger

logger = get_logger(__name__)


@dataclass
class ToolConfig:
    """
    Configuration for external tools
    外部工具配置
    """
    
    executable: str
    parameters: Dict[str, Any]
    timeout: int
    additional_config: Dict[str, Any]
    
    def __post_init__(self) -> None:
        """Validate tool configuration after initialization"""
        if not self.executable:
            raise ConfigurationError("Tool executable cannot be empty")
        if self.timeout <= 0:
            raise ConfigurationError("Tool timeout must be positive")
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert ToolConfig to dictionary"""
        result = {
            "executable": self.executable,
            "parameters": self.parameters,
            "timeout": self.timeout,
        }
        result.update(self.additional_config)
        return result


class Config:
    """
    Main configuration class for NBS annotation pipeline
    NBS注释流水线的主要配置类
    
    提供配置加载、验证、访问和管理功能
    """
    
    def __init__(self, config_data: Dict[str, Any], config_path: Optional[Path] = None) -> None:
        """
        Initialize configuration
        
        Args:
            config_data: Configuration dictionary
            config_path: Path to configuration file
        """
        self.config_path = config_path
        self._config_data = config_data
        
        # Resolve relative paths in configuration
        self._resolve_relative_paths()
        
        self._validate_config()
        
        # Set up environment variables if specified
        self._setup_environment()
        
        logger.info(f"Configuration loaded from: {config_path or 'dictionary'}")
    
    def _validate_config(self) -> None:
        """
        Validate configuration structure and values
        验证配置结构和值
        """
        required_sections = ["tools", "io", "pipeline", "logging"]
        
        for section in required_sections:
            if section not in self._config_data:
                raise ConfigurationError(
                    f"Required configuration section missing: {section}",
                    config_path=str(self.config_path) if self.config_path else None,
                    config_key=section,
                )
        
        # Validate tools section
        tools_config = self._config_data.get("tools", {})
        required_tools = ["nlr_annotator", "miniprot", "augustus", "evm"]
        
        for tool in required_tools:
            if tool not in tools_config:
                raise ConfigurationError(
                    f"Required tool configuration missing: {tool}",
                    config_path=str(self.config_path) if self.config_path else None,
                    config_key=f"tools.{tool}",
                )
        
        # Validate I/O paths
        io_config = self._config_data.get("io", {})
        for key, path_str in io_config.items():
            if key.endswith("_dir") and path_str:
                path = Path(path_str)
                # Create directories if they don't exist
                if not path.exists():
                    try:
                        path.mkdir(parents=True, exist_ok=True)
                        logger.info(f"Created directory: {path}")
                    except OSError as e:
                        logger.warning(f"Could not create directory {path}: {e}")
    
    def _resolve_relative_paths(self) -> None:
        """
        Resolve relative paths in configuration to absolute paths
        将配置中的相对路径解析为绝对路径
        """
        if not self.config_path:
            # If no config path, assume we're in project root
            project_root = Path.cwd()
        else:
            # Config file is in config/ subdirectory, so project root is parent
            project_root = self.config_path.parent.parent
        
        # Define path keys that need resolution
        path_keys = [
            ("tools", "nlr_annotator", "jar_path"),
            ("tools", "nlr_annotator", "mot_file"),
            ("tools", "nlr_annotator", "store_file"),
            ("tools", "miniprot", "executable"),
            ("tools", "augustus", "executable"),
            ("tools", "augustus", "config_dir"),
            ("tools", "augustus", "scripts_path"),
            ("tools", "evm", "executable_dir"),
            ("post_processing", "gene_renaming", "tool_path"),
        ]
        
        # Resolve tool paths
        for path_key in path_keys:
            self._resolve_path_in_config(path_key, project_root)
        
        # Resolve environment variable paths
        self._resolve_environment_paths(project_root)
    
    def _resolve_path_in_config(self, path_key: tuple, project_root: Path) -> None:
        """
        Resolve a specific path in configuration
        解析配置中的特定路径
        
        Args:
            path_key: Tuple of keys to navigate to the path value
            project_root: Project root directory
        """
        current = self._config_data
        
        # Navigate to the parent of the target key
        for key in path_key[:-1]:
            if key not in current:
                return  # Path doesn't exist in config
            current = current[key]
        
        # Get the final key and value
        final_key = path_key[-1]
        if final_key not in current:
            return  # Path doesn't exist in config
        
        path_value = current[final_key]
        if not isinstance(path_value, str):
            return  # Not a string path
        
        # Convert relative path to absolute path
        path_obj = Path(path_value)
        if not path_obj.is_absolute():
            absolute_path = project_root / path_obj
            current[final_key] = str(absolute_path)
            logger.debug(f"Resolved path {'.'.join(path_key)}: {path_value} -> {absolute_path}")

    def _resolve_environment_paths(self, project_root: Path) -> None:
        """
        Resolve relative paths in environment variables
        解析环境变量中的相对路径
        
        Args:
            project_root: Project root directory
        """
        env_config = self._config_data.get("environment", {})
        
        # Resolve environment variable paths
        env_vars = env_config.get("variables", {})
        for var_name, var_value in env_vars.items():
            if isinstance(var_value, str):
                path_obj = Path(var_value)
                if not path_obj.is_absolute():
                    absolute_path = project_root / path_obj
                    env_vars[var_name] = str(absolute_path)
                    logger.debug(f"Resolved environment variable {var_name}: {var_value} -> {absolute_path}")
        
        # Resolve PATH extensions
        path_extensions = env_config.get("path_extensions", [])
        for i, path_value in enumerate(path_extensions):
            if isinstance(path_value, str):
                path_obj = Path(path_value)
                if not path_obj.is_absolute():
                    absolute_path = project_root / path_obj
                    path_extensions[i] = str(absolute_path)
                    logger.debug(f"Resolved path extension: {path_value} -> {absolute_path}")

    def _setup_environment(self) -> None:
        """
        Set up environment variables from configuration
        从配置设置环境变量
        """
        env_config = self._config_data.get("environment", {})
        
        # Set environment variables
        env_vars = env_config.get("variables", {})
        for var_name, var_value in env_vars.items():
            if var_value:
                os.environ[var_name] = str(var_value)
                logger.debug(f"Set environment variable: {var_name}={var_value}")
        
        # Extend PATH
        path_extensions = env_config.get("path_extensions", [])
        if path_extensions:
            current_path = os.environ.get("PATH", "")
            new_paths = [p for p in path_extensions if Path(p).exists()]
            if new_paths:
                extended_path = ":".join(new_paths + [current_path])
                os.environ["PATH"] = extended_path
                logger.debug(f"Extended PATH with: {new_paths}")
    
    def get(self, key: str, default: Any = None) -> Any:
        """
        Get configuration value using dot notation
        使用点表示法获取配置值
        
        Args:
            key: Configuration key (e.g., 'tools.augustus.parameters')
            default: Default value if key not found
            
        Returns:
            Configuration value
        """
        keys = key.split('.')
        value = self._config_data
        
        try:
            for k in keys:
                value = value[k]
            return value
        except (KeyError, TypeError):
            return default
    
    def get_tool_config(self, tool_name: str) -> ToolConfig:
        """
        Get configuration for specific tool
        获取特定工具的配置
        
        Args:
            tool_name: Name of the tool
            
        Returns:
            ToolConfig object
        """
        tool_data = self.get(f"tools.{tool_name}")
        
        if not tool_data:
            raise ConfigurationError(
                f"Tool configuration not found: {tool_name}",
                config_path=str(self.config_path) if self.config_path else None,
                config_key=f"tools.{tool_name}",
            )
        
        return ToolConfig(
            executable=tool_data.get("executable", tool_name),
            parameters=tool_data.get("parameters", {}),
            timeout=tool_data.get("timeout", 3600),
            additional_config={k: v for k, v in tool_data.items() 
                             if k not in ["executable", "parameters", "timeout"]},
        )
    
    def get_io_config(self) -> Dict[str, Any]:
        """
        Get I/O configuration
        获取I/O配置
        
        Returns:
            I/O configuration dictionary
        """
        return self.get("io", {})
    
    def get_pipeline_config(self) -> Dict[str, Any]:
        """
        Get pipeline configuration
        获取流水线配置
        
        Returns:
            Pipeline configuration dictionary
        """
        return self.get("pipeline", {})
    
    def get_validation_config(self) -> Dict[str, Any]:
        """
        Get validation configuration
        获取验证配置
        
        Returns:
            Validation configuration dictionary
        """
        return self.get("validation", {})
    
    def get_logging_config(self) -> Dict[str, Any]:
        """
        Get logging configuration
        获取日志配置
        
        Returns:
            Logging configuration dictionary
        """
        return self.get("logging", {})
    
    def get_output_config(self) -> Dict[str, Any]:
        """
        Get output configuration
        获取输出配置
        
        Returns:
            Output configuration dictionary
        """
        return self.get("output", {})
    
    def get_quality_control_config(self) -> Dict[str, Any]:
        """
        Get quality control configuration
        获取质量控制配置
        
        Returns:
            Quality control configuration dictionary
        """
        return self.get("quality_control", {})
    
    def set(self, key: str, value: Any) -> None:
        """
        Set configuration value using dot notation
        使用点表示法设置配置值
        
        Args:
            key: Configuration key
            value: Value to set
        """
        keys = key.split('.')
        target = self._config_data
        
        # Navigate to the parent of the target key
        for k in keys[:-1]:
            if k not in target:
                target[k] = {}
            target = target[k]
        
        # Set the value
        target[keys[-1]] = value
        logger.debug(f"Set configuration: {key} = {value}")
    
    def update(self, updates: Dict[str, Any]) -> None:
        """
        Update configuration with new values
        用新值更新配置
        
        Args:
            updates: Dictionary of updates using dot notation keys
        """
        for key, value in updates.items():
            self.set(key, value)
    
    def to_dict(self) -> Dict[str, Any]:
        """
        Get configuration as dictionary
        获取配置字典
        
        Returns:
            Configuration dictionary
        """
        return self._config_data.copy()
    
    def save(self, output_path: Path) -> None:
        """
        Save configuration to YAML file
        保存配置到YAML文件
        
        Args:
            output_path: Path to save configuration
        """
        with open(output_path, 'w', encoding='utf-8') as f:
            yaml.dump(self._config_data, f, default_flow_style=False, 
                     allow_unicode=True, indent=2)
        
        logger.info(f"Configuration saved to: {output_path}")
    
    def merge(self, other_config: 'Config') -> 'Config':
        """
        Merge with another configuration
        与另一个配置合并
        
        Args:
            other_config: Configuration to merge with
            
        Returns:
            New merged configuration
        """
        merged_data = self._deep_merge(self._config_data, other_config._config_data)
        return Config(merged_data)
    
    @staticmethod
    def _deep_merge(dict1: Dict[str, Any], dict2: Dict[str, Any]) -> Dict[str, Any]:
        """
        Deep merge two dictionaries
        深度合并两个字典
        """
        result = dict1.copy()
        
        for key, value in dict2.items():
            if key in result and isinstance(result[key], dict) and isinstance(value, dict):
                result[key] = Config._deep_merge(result[key], value)
            else:
                result[key] = value
        
        return result
    
    def create_environment_profile(self) -> Dict[str, str]:
        """
        Create environment profile for external tools
        为外部工具创建环境配置文件
        
        Returns:
            Environment variables dictionary
        """
        env_profile = os.environ.copy()
        
        # Add configuration-specific environment variables
        env_config = self.get("environment", {})
        env_vars = env_config.get("variables", {})
        
        for var_name, var_value in env_vars.items():
            if var_value:
                env_profile[var_name] = str(var_value)
        
        # Extend PATH
        path_extensions = env_config.get("path_extensions", [])
        if path_extensions:
            current_path = env_profile.get("PATH", "")
            new_paths = [p for p in path_extensions if Path(p).exists()]
            if new_paths:
                env_profile["PATH"] = ":".join(new_paths + [current_path])
        
        return env_profile
    
    def validate_tool_availability(self) -> Dict[str, bool]:
        """
        Check if all configured tools are available
        检查所有配置的工具是否可用
        
        Returns:
            Dictionary mapping tool names to availability status
        """
        import shutil
        
        availability = {}
        tools_config = self.get("tools", {})
        
        for tool_name, tool_config in tools_config.items():
            # Handle different tool configuration patterns
            if tool_name == "evm":
                # EVM uses executable_dir instead of executable
                executable_dir = tool_config.get("executable_dir")
                if executable_dir:
                    evm_dir = Path(executable_dir)
                    # Check if EVM directory exists and contains required scripts
                    if evm_dir.exists():
                        # Check for essential EVM scripts
                        scripts = tool_config.get("scripts", {})
                        partition_script = evm_dir / scripts.get("partition", "EvmUtils/partition_EVM_inputs.pl")
                        execute_script = evm_dir / scripts.get("execute", "EvmUtils/execute_EVM_in_parallel.pl")
                        
                        if partition_script.exists() and execute_script.exists():
                            availability[tool_name] = True
                        else:
                            availability[tool_name] = False
                            logger.warning(f"EVM scripts not found in {executable_dir}")
                    else:
                        availability[tool_name] = False
                        logger.warning(f"EVM directory not found: {executable_dir}")
                else:
                    availability[tool_name] = False
                    logger.warning(f"EVM executable_dir not configured")
            elif tool_name == "nlr_annotator":
                # NLR-Annotator uses jar_path
                jar_path = tool_config.get("jar_path")
                executable = tool_config.get("executable", "java")
                
                # Check if Java is available
                java_available = shutil.which(executable) is not None
                # Check if JAR file exists
                jar_available = jar_path and Path(jar_path).exists()
                
                availability[tool_name] = java_available and jar_available
                if not java_available:
                    logger.warning(f"Java executable not found: {executable}")
                if not jar_available:
                    logger.warning(f"NLR-Annotator JAR not found: {jar_path}")
            else:
                # Standard tools with executable field
                executable = tool_config.get("executable", tool_name)
                
                # Check if executable exists in PATH
                if shutil.which(executable):
                    availability[tool_name] = True
                else:
                    # Check if it's an absolute path
                    if Path(executable).exists():
                        availability[tool_name] = True
                    else:
                        availability[tool_name] = False
                        logger.warning(f"Tool not found: {tool_name} ({executable})")
        
        return availability
    
    def __str__(self) -> str:
        """String representation of configuration"""
        return f"Config(path={self.config_path}, sections={list(self._config_data.keys())})"
    
    def __repr__(self) -> str:
        """Detailed string representation"""
        return f"Config(config_path='{self.config_path}', config_data={self._config_data})"


def load_config(config_path: Optional[Union[str, Path]] = None) -> Config:
    """
    Load configuration from YAML file or use default
    从YAML文件加载配置或使用默认配置
    
    Args:
        config_path: Path to configuration file. If None, uses default config.
        
    Returns:
        Config object
        
    Raises:
        ConfigurationError: If configuration file cannot be loaded or is invalid
    """
    if config_path is None:
        # Use default configuration file
        config_path = Path(__file__).parent.parent.parent.parent / "config" / "default.yaml"
    
    config_path = Path(config_path)
    
    if not config_path.exists():
        raise ConfigurationError(
            f"Configuration file does not exist: {config_path}",
            config_path=str(config_path),
        )
    
    try:
        with open(config_path, 'r', encoding='utf-8') as f:
            config_data = yaml.safe_load(f)
        
        if not config_data:
            raise ConfigurationError(
                f"Configuration file is empty: {config_path}",
                config_path=str(config_path),
            )
        
        return Config(config_data, config_path)
        
    except yaml.YAMLError as e:
        raise ConfigurationError(
            f"YAML parsing error in {config_path}: {e}",
            config_path=str(config_path),
        ) from e
    except Exception as e:
        raise ConfigurationError(
            f"Error loading configuration from {config_path}: {e}",
            config_path=str(config_path),
        ) from e


def load_config_with_overrides(
    base_config_path: Optional[Union[str, Path]] = None,
    override_config_path: Optional[Union[str, Path]] = None,
    env_overrides: Optional[Dict[str, Any]] = None,
) -> Config:
    """
    Load configuration with overrides
    加载带覆盖选项的配置
    
    Args:
        base_config_path: Path to base configuration file
        override_config_path: Path to override configuration file
        env_overrides: Environment-based overrides
        
    Returns:
        Config object with overrides applied
    """
    # Load base configuration
    config = load_config(base_config_path)
    
    # Apply override file if specified
    if override_config_path:
        override_path = Path(override_config_path)
        if override_path.exists():
            override_config = load_config(override_path)
            config = config.merge(override_config)
            logger.info(f"Applied configuration overrides from: {override_path}")
    
    # Apply environment overrides
    if env_overrides:
        config.update(env_overrides)
        logger.info(f"Applied {len(env_overrides)} environment overrides")
    
    return config


def create_config_template(output_path: Path) -> None:
    """
    Create a configuration template file
    创建配置模板文件
    
    Args:
        output_path: Path to save template
    """
    template_config = {
        "project": {
            "name": "nbs-gene-annotation",
            "version": "0.1.0",
            "description": "NBS gene annotation pipeline configuration",
        },
        "io": {
            "input_dir": "data/input",
            "output_dir": "data/output",
            "interim_dir": "data/interim",
            "log_dir": "logs",
        },
        "tools": {
            "nlr_annotator": {
                "executable": "java",
                "jar_path": "/path/to/NLR-Annotator.jar",
                "parameters": {},
                "timeout": 3600,
            },
            "miniprot": {
                "executable": "miniprot",
                "parameters": {"threads": 8},
                "timeout": 7200,
            },
            "augustus": {
                "executable": "augustus",
                "config_dir": "/opt/augustus/config",
                "parameters": {"genemodel": "complete"},
                "timeout": 1800,
            },
            "evm": {
                "executable_dir": "/opt/EVidenceModeler",
                "parameters": {"cpu": 8},
                "timeout": 14400,
            },
        },
        "pipeline": {
            "parallel": {"max_workers": 4},
            "memory": {"max_memory_gb": 32},
        },
        "logging": {
            "level": "INFO",
            "format": "%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        },
    }
    
    with open(output_path, 'w', encoding='utf-8') as f:
        yaml.dump(template_config, f, default_flow_style=False, 
                 allow_unicode=True, indent=2)
    
    logger.info(f"Configuration template saved to: {output_path}")