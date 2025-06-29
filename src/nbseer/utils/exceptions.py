"""
Exception classes for NBS annotation pipeline
NBS注释流水线异常类定义

定义了流水线中使用的所有自定义异常类型
"""

from typing import Any, Dict, Optional


class NBSAnnotationError(Exception):
    """
    Base exception for NBS annotation pipeline
    NBS注释流水线的基础异常类
    """
    
    def __init__(
        self,
        message: str,
        details: Optional[Dict[str, Any]] = None,
        component: Optional[str] = None,
    ) -> None:
        """
        Initialize NBS annotation error
        
        Args:
            message: Error message / 错误消息
            details: Additional error details / 错误详细信息
            component: Component that caused the error / 引起错误的组件
        """
        super().__init__(message)
        self.details = details or {}
        self.component = component
        
    def __str__(self) -> str:
        """String representation of the error"""
        base_msg = super().__str__()
        
        if self.component:
            base_msg = f"[{self.component}] {base_msg}"
            
        if self.details:
            detail_str = ", ".join(f"{k}={v}" for k, v in self.details.items())
            base_msg = f"{base_msg} ({detail_str})"
            
        return base_msg


class ToolExecutionError(NBSAnnotationError):
    """
    Exception raised when external tool execution fails
    外部工具执行失败时抛出的异常
    """
    
    def __init__(
        self,
        tool_name: str,
        command: str,
        return_code: int,
        stdout: str = "",
        stderr: str = "",
        **kwargs: Any,
    ) -> None:
        """
        Initialize tool execution error
        
        Args:
            tool_name: Name of the tool that failed / 失败的工具名称
            command: Command that was executed / 执行的命令
            return_code: Exit code of the command / 命令退出码
            stdout: Standard output / 标准输出
            stderr: Standard error / 标准错误
        """
        message = f"Tool '{tool_name}' failed with exit code {return_code}"
        details = {
            "command": command,
            "return_code": return_code,
            "stdout": stdout[:500] + "..." if len(stdout) > 500 else stdout,
            "stderr": stderr[:500] + "..." if len(stderr) > 500 else stderr,
        }
        super().__init__(message, details, tool_name, **kwargs)
        
        self.tool_name = tool_name
        self.command = command
        self.return_code = return_code
        self.stdout = stdout
        self.stderr = stderr


class ValidationError(NBSAnnotationError):
    """
    Exception raised when data validation fails
    数据验证失败时抛出的异常
    """
    
    def __init__(
        self,
        message: str,
        file_path: Optional[str] = None,
        line_number: Optional[int] = None,
        validation_type: Optional[str] = None,
        **kwargs: Any,
    ) -> None:
        """
        Initialize validation error
        
        Args:
            message: Validation error message / 验证错误消息
            file_path: Path to the file that failed validation / 验证失败的文件路径
            line_number: Line number where validation failed / 验证失败的行号
            validation_type: Type of validation that failed / 失败的验证类型
        """
        details = {}
        if file_path:
            details["file_path"] = file_path
        if line_number:
            details["line_number"] = line_number
        if validation_type:
            details["validation_type"] = validation_type
            
        super().__init__(message, details, "DataValidator", **kwargs)
        
        self.file_path = file_path
        self.line_number = line_number
        self.validation_type = validation_type


class ConfigurationError(NBSAnnotationError):
    """
    Exception raised when configuration is invalid
    配置无效时抛出的异常
    """
    
    def __init__(
        self,
        message: str,
        config_path: Optional[str] = None,
        config_key: Optional[str] = None,
        **kwargs: Any,
    ) -> None:
        """
        Initialize configuration error
        
        Args:
            message: Configuration error message / 配置错误消息
            config_path: Path to the configuration file / 配置文件路径
            config_key: Configuration key that caused the error / 引起错误的配置键
        """
        details = {}
        if config_path:
            details["config_path"] = config_path
        if config_key:
            details["config_key"] = config_key
            
        super().__init__(message, details, "ConfigManager", **kwargs)
        
        self.config_path = config_path
        self.config_key = config_key


class ModelTrainingError(NBSAnnotationError):
    """
    Exception raised when model training fails
    模型训练失败时抛出的异常
    """
    
    def __init__(
        self,
        message: str,
        model_type: Optional[str] = None,
        training_data_path: Optional[str] = None,
        **kwargs: Any,
    ) -> None:
        """
        Initialize model training error
        
        Args:
            message: Training error message / 训练错误消息
            model_type: Type of model being trained / 正在训练的模型类型
            training_data_path: Path to training data / 训练数据路径
        """
        details = {}
        if model_type:
            details["model_type"] = model_type
        if training_data_path:
            details["training_data_path"] = training_data_path
            
        super().__init__(message, details, "ModelTrainer", **kwargs)
        
        self.model_type = model_type
        self.training_data_path = training_data_path


class PipelineError(NBSAnnotationError):
    """
    Exception raised when pipeline execution fails
    流水线执行失败时抛出的异常
    """
    
    def __init__(
        self,
        message: str,
        pipeline_stage: Optional[str] = None,
        input_files: Optional[list] = None,
        **kwargs: Any,
    ) -> None:
        """
        Initialize pipeline error
        
        Args:
            message: Pipeline error message / 流水线错误消息
            pipeline_stage: Stage where the error occurred / 发生错误的阶段
            input_files: Input files being processed / 正在处理的输入文件
        """
        details = {}
        if pipeline_stage:
            details["pipeline_stage"] = pipeline_stage
        if input_files:
            details["input_files"] = input_files
            
        super().__init__(message, details, "PipelineCoordinator", **kwargs)
        
        self.pipeline_stage = pipeline_stage
        self.input_files = input_files


class DataProcessingError(NBSAnnotationError):
    """
    Exception raised when data processing fails
    数据处理失败时抛出的异常
    """
    
    def __init__(
        self,
        message: str,
        data_type: Optional[str] = None,
        processing_step: Optional[str] = None,
        **kwargs: Any,
    ) -> None:
        """
        Initialize data processing error
        
        Args:
            message: Data processing error message / 数据处理错误消息
            data_type: Type of data being processed / 正在处理的数据类型
            processing_step: Processing step that failed / 失败的处理步骤
        """
        details = {}
        if data_type:
            details["data_type"] = data_type
        if processing_step:
            details["processing_step"] = processing_step
            
        super().__init__(message, details, "DataProcessor", **kwargs)
        
        self.data_type = data_type
        self.processing_step = processing_step