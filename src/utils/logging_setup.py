"""
多级日志系统配置模块
Multi-level Logging System Configuration

提供统一的日志配置和管理功能，支持多种输出格式和级别。
"""

import logging
import logging.handlers
import os
import sys
from datetime import datetime
from pathlib import Path
from typing import Optional, Dict, Any


class NBSLogger:
    """NBS基因注释项目专用日志器"""
    
    def __init__(self, name: str = "nbs_annotation"):
        self.name = name
        self.logger = logging.getLogger(name)
        self._configured = False
        
    def setup_logging(
        self,
        log_level: str = "INFO",
        log_dir: Optional[str] = None,
        console_output: bool = True,
        file_output: bool = True,
        max_file_size: int = 10 * 1024 * 1024,  # 10MB
        backup_count: int = 5
    ) -> logging.Logger:
        """
        配置多级日志系统
        
        Args:
            log_level: 日志级别 (DEBUG, INFO, WARNING, ERROR, CRITICAL)
            log_dir: 日志文件目录，默认为项目根目录下的logs/
            console_output: 是否输出到控制台
            file_output: 是否输出到文件
            max_file_size: 单个日志文件最大大小（字节）
            backup_count: 保留的备份文件数量
            
        Returns:
            配置好的Logger对象
        """
        if self._configured:
            return self.logger
            
        # 设置日志级别
        level = getattr(logging, log_level.upper(), logging.INFO)
        self.logger.setLevel(level)
        
        # 清除现有的handlers
        self.logger.handlers.clear()
        
        # 创建格式化器
        formatters = self._create_formatters()
        
        # 控制台输出
        if console_output:
            console_handler = self._create_console_handler(formatters['console'])
            self.logger.addHandler(console_handler)
        
        # 文件输出
        if file_output:
            if log_dir is None:
                # 自动确定项目根目录
                current_dir = Path(__file__).parent
                project_root = current_dir.parent.parent
                log_dir = project_root / "logs"
            
            log_dir = Path(log_dir)
            log_dir.mkdir(exist_ok=True)
            
            # 创建不同级别的文件处理器
            file_handlers = self._create_file_handlers(
                log_dir, formatters['file'], max_file_size, backup_count
            )
            
            for handler in file_handlers:
                self.logger.addHandler(handler)
        
        self._configured = True
        self.logger.info(f"日志系统初始化完成 - 级别: {log_level}")
        return self.logger
    
    def _create_formatters(self) -> Dict[str, logging.Formatter]:
        """创建不同的格式化器"""
        return {
            'console': logging.Formatter(
                '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                datefmt='%Y-%m-%d %H:%M:%S'
            ),
            'file': logging.Formatter(
                '%(asctime)s - %(name)s - %(levelname)s - %(filename)s:%(lineno)d - %(funcName)s() - %(message)s',
                datefmt='%Y-%m-%d %H:%M:%S'
            ),
            'debug': logging.Formatter(
                '%(asctime)s - %(name)s - %(levelname)s - %(pathname)s:%(lineno)d - %(funcName)s() - %(message)s',
                datefmt='%Y-%m-%d %H:%M:%S.%f'
            )
        }
    
    def _create_console_handler(self, formatter: logging.Formatter) -> logging.StreamHandler:
        """创建控制台处理器"""
        handler = logging.StreamHandler(sys.stdout)
        handler.setFormatter(formatter)
        return handler
    
    def _create_file_handlers(
        self,
        log_dir: Path,
        formatter: logging.Formatter,
        max_file_size: int,
        backup_count: int
    ) -> list:
        """创建文件处理器"""
        handlers = []
        
        # 通用日志文件（所有级别）
        general_log = log_dir / "nbs_annotation.log"
        general_handler = logging.handlers.RotatingFileHandler(
            general_log, maxBytes=max_file_size, backupCount=backup_count
        )
        general_handler.setFormatter(formatter)
        general_handler.setLevel(logging.DEBUG)
        handlers.append(general_handler)
        
        # 错误日志文件（只记录WARNING及以上）
        error_log = log_dir / "nbs_annotation_error.log"
        error_handler = logging.handlers.RotatingFileHandler(
            error_log, maxBytes=max_file_size, backupCount=backup_count
        )
        error_handler.setFormatter(formatter)
        error_handler.setLevel(logging.WARNING)
        handlers.append(error_handler)
        
        # 调试日志文件（只在DEBUG级别时创建）
        if self.logger.level <= logging.DEBUG:
            debug_log = log_dir / "nbs_annotation_debug.log"
            debug_handler = logging.handlers.RotatingFileHandler(
                debug_log, maxBytes=max_file_size, backupCount=backup_count
            )
            debug_handler.setFormatter(formatter)
            debug_handler.setLevel(logging.DEBUG)
            handlers.append(debug_handler)
        
        return handlers
    
    def get_logger(self) -> logging.Logger:
        """获取配置好的logger实例"""
        if not self._configured:
            self.setup_logging()
        return self.logger


# 全局日志器实例
_global_logger = None


def setup_logging(
    log_level: str = "INFO",
    log_dir: Optional[str] = None,
    console_output: bool = True,
    file_output: bool = True,
    **kwargs
) -> logging.Logger:
    """
    设置全局日志系统的便捷函数
    
    Args:
        log_level: 日志级别
        log_dir: 日志目录
        console_output: 控制台输出
        file_output: 文件输出
        **kwargs: 其他配置参数
        
    Returns:
        配置好的Logger对象
    """
    global _global_logger
    
    if _global_logger is None:
        _global_logger = NBSLogger("nbs_annotation")
    
    return _global_logger.setup_logging(
        log_level=log_level,
        log_dir=log_dir,
        console_output=console_output,
        file_output=file_output,
        **kwargs
    )


def get_logger(name: Optional[str] = None) -> logging.Logger:
    """
    获取日志器实例
    
    Args:
        name: 日志器名称，默认使用全局日志器
        
    Returns:
        Logger对象
    """
    if name is None:
        global _global_logger
        if _global_logger is None:
            _global_logger = NBSLogger("nbs_annotation")
            _global_logger.setup_logging()
        return _global_logger.get_logger()
    else:
        return logging.getLogger(name)


def log_function_call(func):
    """
    装饰器：自动记录函数调用
    
    Usage:
        @log_function_call
        def my_function(arg1, arg2):
            pass
    """
    def wrapper(*args, **kwargs):
        logger = get_logger()
        logger.debug(f"调用函数: {func.__name__} - args: {args}, kwargs: {kwargs}")
        try:
            result = func(*args, **kwargs)
            logger.debug(f"函数 {func.__name__} 执行成功")
            return result
        except Exception as e:
            logger.error(f"函数 {func.__name__} 执行失败: {str(e)}")
            raise
    return wrapper


def log_execution_time(func):
    """
    装饰器：记录函数执行时间
    
    Usage:
        @log_execution_time
        def slow_function():
            pass
    """
    def wrapper(*args, **kwargs):
        logger = get_logger()
        start_time = datetime.now()
        logger.info(f"开始执行函数: {func.__name__}")
        
        try:
            result = func(*args, **kwargs)
            end_time = datetime.now()
            duration = (end_time - start_time).total_seconds()
            logger.info(f"函数 {func.__name__} 执行完成，耗时: {duration:.2f}秒")
            return result
        except Exception as e:
            end_time = datetime.now()
            duration = (end_time - start_time).total_seconds()
            logger.error(f"函数 {func.__name__} 执行失败，耗时: {duration:.2f}秒，错误: {str(e)}")
            raise
    return wrapper


if __name__ == "__main__":
    # 测试日志系统
    logger = setup_logging(log_level="DEBUG")
    
    logger.debug("这是一条调试信息")
    logger.info("这是一条信息")
    logger.warning("这是一条警告")
    logger.error("这是一条错误")
    logger.critical("这是一条严重错误")
    
    print("日志系统测试完成") 