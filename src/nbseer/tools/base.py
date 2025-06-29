"""
Base classes for external tool interfaces
外部工具接口基类

提供统一的外部工具调用接口和结果处理框架
"""

import shutil
import subprocess
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any, Dict, List, Optional, Union
from dataclasses import dataclass
from datetime import datetime

from ..utils.exceptions import ToolExecutionError, ValidationError
from ..utils.logging_setup import get_logger

logger = get_logger(__name__)


@dataclass
class ToolResult:
    """
    Standardized result object for external tool execution
    外部工具执行的标准化结果对象
    """
    
    tool_name: str
    command: str
    return_code: int
    stdout: str
    stderr: str
    execution_time: float
    input_files: List[Path]
    output_files: List[Path]
    success: bool
    metadata: Dict[str, Any]
    
    def __post_init__(self) -> None:
        """Validate and process result after initialization"""
        self.success = self.return_code == 0
        
    def _convert_paths_to_strings(self, obj: Any) -> Any:
        """Recursively convert Path objects to strings for JSON serialization"""
        if isinstance(obj, Path):
            return str(obj)
        elif isinstance(obj, dict):
            return {key: self._convert_paths_to_strings(value) for key, value in obj.items()}
        elif isinstance(obj, list):
            return [self._convert_paths_to_strings(item) for item in obj]
        elif isinstance(obj, tuple):
            return tuple(self._convert_paths_to_strings(item) for item in obj)
        else:
            return obj
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert result to dictionary"""
        return {
            "tool_name": self.tool_name,
            "command": self.command,
            "return_code": self.return_code,
            "success": self.success,
            "execution_time": self.execution_time,
            "input_files": [str(f) for f in self.input_files],
            "output_files": [str(f) for f in self.output_files],
            "metadata": self._convert_paths_to_strings(self.metadata),
            "timestamp": datetime.now().isoformat(),
        }


class ExternalTool(ABC):
    """
    Abstract base class for external bioinformatics tools
    外部生物信息学工具的抽象基类
    
    所有外部工具接口都应继承此类并实现必要的抽象方法
    """
    
    def __init__(
        self,
        tool_name: str,
        executable_path: Optional[Path] = None,
        config: Optional[Dict[str, Any]] = None,
    ) -> None:
        """
        Initialize external tool interface
        
        Args:
            tool_name: Name of the tool / 工具名称
            executable_path: Path to tool executable / 工具可执行文件路径
            config: Tool-specific configuration / 工具特定配置
        """
        self.tool_name = tool_name
        self.config = config or {}
        
        # Determine executable path from various sources
        if executable_path:
            self.executable_path = Path(executable_path)
        elif self.config.get("executable"):
            self.executable_path = Path(self.config["executable"])
        else:
            self.executable_path = self._find_executable()
        self.logger = get_logger(f"{__name__}.{tool_name}")
        
        # Validate tool installation
        self._validate_installation()
        
    def _find_executable(self) -> Path:
        """
        Find tool executable in system PATH
        在系统PATH中查找工具可执行文件
        """
        executable = shutil.which(self.tool_name)
        if not executable:
            raise ToolExecutionError(
                tool_name=self.tool_name,
                command="which",
                return_code=1,
                stderr=f"Tool '{self.tool_name}' not found in PATH",
            )
        return Path(executable)
    
    def _validate_installation(self) -> None:
        """
        Validate that the tool is properly installed
        验证工具是否正确安装
        """
        if not self.executable_path.exists():
            raise ToolExecutionError(
                tool_name=self.tool_name,
                command="existence_check",
                return_code=1,
                stderr=f"Executable not found: {self.executable_path}",
            )
            
        if not self.executable_path.is_file():
            raise ToolExecutionError(
                tool_name=self.tool_name,
                command="file_check",
                return_code=1,
                stderr=f"Path is not a file: {self.executable_path}",
            )
            
        # Check if executable is actually executable
        if not self.executable_path.stat().st_mode & 0o111:
            raise ToolExecutionError(
                tool_name=self.tool_name,
                command="permission_check",
                return_code=1,
                stderr=f"File is not executable: {self.executable_path}",
            )
    
    @abstractmethod
    def prepare_input(self, **kwargs: Any) -> Dict[str, Path]:
        """
        Prepare input files and parameters for tool execution
        为工具执行准备输入文件和参数
        
        Returns:
            Dictionary mapping parameter names to file paths
        """
        pass
    
    @abstractmethod
    def build_command(
        self,
        input_files: Dict[str, Path],
        output_dir: Path,
        **kwargs: Any,
    ) -> List[str]:
        """
        Build command line for tool execution
        构建工具执行的命令行
        
        Args:
            input_files: Prepared input files
            output_dir: Output directory
            **kwargs: Additional parameters
            
        Returns:
            Command as list of strings
        """
        pass
    
    @abstractmethod
    def parse_output(self, output_dir: Path) -> Dict[str, Any]:
        """
        Parse tool output files
        解析工具输出文件
        
        Args:
            output_dir: Directory containing output files
            
        Returns:
            Parsed results as dictionary
        """
        pass
    
    def validate_input(self, **kwargs: Any) -> None:
        """
        Validate input parameters and files
        验证输入参数和文件
        
        Can be overridden by subclasses for specific validation
        """
        # Basic validation - check required files exist
        for key, value in kwargs.items():
            if isinstance(value, (str, Path)):
                path = Path(value)
                if key.endswith('_file') or key.endswith('_path'):
                    if not path.exists():
                        raise ValidationError(
                            f"Input file does not exist: {path}",
                            file_path=str(path),
                            validation_type="file_existence",
                        )
    
    def execute(
        self,
        output_dir: Path,
        timeout: Optional[int] = None,
        **kwargs: Any,
    ) -> ToolResult:
        """
        Execute the external tool
        执行外部工具
        
        Args:
            output_dir: Directory for output files
            timeout: Execution timeout in seconds
            **kwargs: Tool-specific parameters
            
        Returns:
            ToolResult object containing execution results
        """
        start_time = datetime.now()
        
        try:
            # Validate inputs
            self.validate_input(**kwargs)
            
            # Prepare inputs
            input_files = self.prepare_input(**kwargs)
            
            # Build command
            command = self.build_command(input_files, output_dir, **kwargs)
            
            # Ensure output directory exists
            output_dir.mkdir(parents=True, exist_ok=True)
            
            # Log execution start
            self.logger.info(
                f"Starting {self.tool_name} execution: {' '.join(command)}"
            )
            
            # Execute command
            result = subprocess.run(
                command,
                cwd=output_dir,
                capture_output=True,
                text=True,
                timeout=timeout,
            )
            
            execution_time = (datetime.now() - start_time).total_seconds()
            
            # Get output files
            output_files = list(output_dir.glob("*"))
            
            # Create result object
            tool_result = ToolResult(
                tool_name=self.tool_name,
                command=" ".join(command),
                return_code=result.returncode,
                stdout=result.stdout,
                stderr=result.stderr,
                execution_time=execution_time,
                input_files=list(input_files.values()),
                output_files=output_files,
                success=result.returncode == 0,
                metadata={"config": self.config, "kwargs": kwargs},
            )
            
            # Log result
            if tool_result.success:
                self.logger.info(
                    f"{self.tool_name} completed successfully in "
                    f"{execution_time:.2f}s"
                )
            else:
                self.logger.error(
                    f"{self.tool_name} failed with return code "
                    f"{result.returncode}"
                )
                self.logger.error(f"STDERR: {result.stderr}")
            
            return tool_result
            
        except subprocess.TimeoutExpired as e:
            execution_time = (datetime.now() - start_time).total_seconds()
            raise ToolExecutionError(
                tool_name=self.tool_name,
                command=" ".join(command) if 'command' in locals() else "unknown",
                return_code=-1,
                stderr=f"Command timed out after {timeout} seconds",
            ) from e
            
        except Exception as e:
            execution_time = (datetime.now() - start_time).total_seconds()
            self.logger.error(f"{self.tool_name} execution failed: {e}")
            
            if isinstance(e, (ToolExecutionError, ValidationError)):
                raise
            else:
                raise ToolExecutionError(
                    tool_name=self.tool_name,
                    command=" ".join(command) if 'command' in locals() else "unknown",
                    return_code=-1,
                    stderr=str(e),
                ) from e
    
    def get_version(self) -> str:
        """
        Get tool version information
        获取工具版本信息
        """
        try:
            # Most tools support --version or -v
            for version_flag in ["--version", "-v", "-version"]:
                try:
                    result = subprocess.run(
                        [str(self.executable_path), version_flag],
                        capture_output=True,
                        text=True,
                        timeout=10,
                    )
                    if result.returncode == 0:
                        return result.stdout.strip()
                except (subprocess.TimeoutExpired, subprocess.CalledProcessError):
                    continue
                    
            return "Unknown version"
            
        except Exception as e:
            self.logger.warning(f"Could not determine version for {self.tool_name}: {e}")
            return "Unknown version"
    
    def check_dependencies(self) -> Dict[str, bool]:
        """
        Check if tool dependencies are available
        检查工具依赖是否可用
        
        Returns:
            Dictionary mapping dependency names to availability status
        """
        # Base implementation - can be overridden by subclasses
        return {"executable": self.executable_path.exists()}
    
    def __str__(self) -> str:
        """String representation of the tool"""
        return f"{self.tool_name} ({self.executable_path})"
    
    def __repr__(self) -> str:
        """Detailed string representation"""
        return (
            f"{self.__class__.__name__}("
            f"tool_name='{self.tool_name}', "
            f"executable_path='{self.executable_path}'"
            f")"
        )