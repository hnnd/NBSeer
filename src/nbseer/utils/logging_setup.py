"""
Logging setup and configuration for NBS annotation pipeline
NBS注释流水线的日志设置和配置

提供统一的日志配置和管理功能
"""

import logging
import logging.handlers
import sys
from pathlib import Path
from typing import Optional, Union
import json
from datetime import datetime


class JsonFormatter(logging.Formatter):
    """
    JSON formatter for structured logging
    结构化日志的JSON格式化器
    """
    
    def format(self, record: logging.LogRecord) -> str:
        """Format log record as JSON"""
        log_entry = {
            "timestamp": datetime.fromtimestamp(record.created).isoformat(),
            "level": record.levelname,
            "logger": record.name,
            "message": record.getMessage(),
            "module": record.module,
            "function": record.funcName,
            "line": record.lineno,
        }
        
        # Add exception info if present
        if record.exc_info:
            log_entry["exception"] = self.formatException(record.exc_info)
        
        # Add extra fields
        for key, value in record.__dict__.items():
            if key not in ["name", "msg", "args", "levelname", "levelno", "pathname", 
                          "filename", "module", "lineno", "funcName", "created", 
                          "msecs", "relativeCreated", "thread", "threadName", 
                          "processName", "process", "getMessage", "exc_info", 
                          "exc_text", "stack_info"]:
                log_entry[key] = value
        
        return json.dumps(log_entry, ensure_ascii=False)


def setup_logging(
    level: str = "INFO",
    log_file: Optional[Union[str, Path]] = None,
    console_output: bool = True,
    json_format: bool = False,
    max_file_size: int = 100 * 1024 * 1024,  # 100MB
    backup_count: int = 5,
) -> None:
    """
    Set up logging configuration for the application
    为应用程序设置日志配置
    
    Args:
        level: Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
        log_file: Path to log file (optional)
        console_output: Whether to output to console
        json_format: Whether to use JSON formatting
        max_file_size: Maximum log file size in bytes
        backup_count: Number of backup log files to keep
    """
    # Convert level string to logging level
    numeric_level = getattr(logging, level.upper(), logging.INFO)
    
    # Clear any existing handlers
    root_logger = logging.getLogger()
    root_logger.handlers.clear()
    root_logger.setLevel(numeric_level)
    
    # Set up formatters
    if json_format:
        formatter = JsonFormatter()
    else:
        formatter = logging.Formatter(
            "%(asctime)s - %(name)s - %(levelname)s - %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S"
        )
    
    # Console handler
    if console_output:
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(numeric_level)
        console_handler.setFormatter(formatter)
        root_logger.addHandler(console_handler)
    
    # File handler
    if log_file:
        log_path = Path(log_file)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Use rotating file handler to prevent huge log files
        file_handler = logging.handlers.RotatingFileHandler(
            log_path,
            maxBytes=max_file_size,
            backupCount=backup_count,
            encoding='utf-8'
        )
        file_handler.setLevel(numeric_level)
        file_handler.setFormatter(formatter)
        root_logger.addHandler(file_handler)
    
    # Set specific logger levels to reduce noise
    logging.getLogger("urllib3").setLevel(logging.WARNING)
    logging.getLogger("requests").setLevel(logging.WARNING)
    logging.getLogger("Bio").setLevel(logging.WARNING)


def get_logger(name: str) -> logging.Logger:
    """
    Get a logger instance with the specified name
    获取指定名称的日志记录器实例
    
    Args:
        name: Logger name (usually __name__)
        
    Returns:
        Logger instance
    """
    return logging.getLogger(name)


def log_function_call(func):
    """
    Decorator to log function calls
    装饰器，用于记录函数调用
    
    Args:
        func: Function to decorate
        
    Returns:
        Decorated function
    """
    def wrapper(*args, **kwargs):
        logger = get_logger(func.__module__)
        logger.debug(f"Calling {func.__name__} with args={args}, kwargs={kwargs}")
        
        try:
            result = func(*args, **kwargs)
            logger.debug(f"{func.__name__} completed successfully")
            return result
        except Exception as e:
            logger.error(f"{func.__name__} failed with error: {e}")
            raise
    
    return wrapper


def log_execution_time(func):
    """
    Decorator to log function execution time
    装饰器，用于记录函数执行时间
    
    Args:
        func: Function to decorate
        
    Returns:
        Decorated function
    """
    def wrapper(*args, **kwargs):
        import time
        
        logger = get_logger(func.__module__)
        start_time = time.time()
        
        try:
            result = func(*args, **kwargs)
            execution_time = time.time() - start_time
            logger.info(f"{func.__name__} completed in {execution_time:.2f} seconds")
            return result
        except Exception as e:
            execution_time = time.time() - start_time
            logger.error(f"{func.__name__} failed after {execution_time:.2f} seconds: {e}")
            raise
    
    return wrapper


class LogCapture:
    """
    Context manager to capture log messages for testing
    用于测试的日志消息捕获上下文管理器
    """
    
    def __init__(self, logger_name: str, level: str = "INFO"):
        """
        Initialize log capture
        
        Args:
            logger_name: Name of logger to capture
            level: Minimum level to capture
        """
        self.logger_name = logger_name
        self.level = getattr(logging, level.upper())
        self.handler = None
        self.logs = []
    
    def __enter__(self):
        """Start capturing logs"""
        self.handler = logging.Handler()
        self.handler.setLevel(self.level)
        self.handler.emit = lambda record: self.logs.append(record)
        
        logger = logging.getLogger(self.logger_name)
        logger.addHandler(self.handler)
        
        return self.logs
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Stop capturing logs"""
        if self.handler:
            logger = logging.getLogger(self.logger_name)
            logger.removeHandler(self.handler)


def configure_pipeline_logging(
    pipeline_id: str,
    output_dir: Path,
    level: str = "INFO",
    console_output: bool = True,
) -> Path:
    """
    Configure logging for a specific pipeline run
    为特定的流水线运行配置日志
    
    Args:
        pipeline_id: Unique pipeline identifier
        output_dir: Pipeline output directory
        level: Logging level
        console_output: Whether to output to console
        
    Returns:
        Path to log file
    """
    log_dir = output_dir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    
    log_file = log_dir / f"{pipeline_id}.log"
    
    setup_logging(
        level=level,
        log_file=log_file,
        console_output=console_output,
        json_format=False,
    )
    
    # Log pipeline start
    logger = get_logger(__name__)
    logger.info(f"Pipeline logging configured: {pipeline_id}")
    logger.info(f"Log file: {log_file}")
    
    return log_file


# Set up basic logging on module import
if not logging.getLogger().handlers:
    setup_logging(level="INFO", console_output=True)