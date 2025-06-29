#!/usr/bin/env python3
"""
Augustus坐标转换模块

将Augustus重新预测的NBS基因结果从相对坐标转换为基因组绝对坐标，
为EVM证据整合做准备。

Author: NBS Annotation Pipeline
Date: 2025-01-21
"""

import os
import re
import glob
import logging
from pathlib import Path
from typing import List, Dict, Tuple, Optional, NamedTuple
from dataclasses import dataclass
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import csv
from tqdm import tqdm

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class RegionInfo(NamedTuple):
    """区域信息数据类"""
    chromosome: str
    start: int
    end: int
    
    @classmethod
    def from_string(cls, region_str: str) -> 'RegionInfo':
        """从字符串解析区域信息 (例如: Chr1:687993-694473)"""
        match = re.match(r'(\w+):(\d+)-(\d+)', region_str)
        if not match:
            raise ValueError(f"Invalid region format: {region_str}")
        
        chromosome = match.group(1)
        start = int(match.group(2))
        end = int(match.group(3))
        
        return cls(chromosome, start, end)

@dataclass
class ConversionResult:
    """坐标转换结果"""
    input_file: str
    output_file: str
    region_info: RegionInfo
    features_converted: int
    conversion_success: bool
    error_message: Optional[str] = None

@dataclass
class ConversionStats:
    """批量转换统计"""
    total_files: int
    successful_conversions: int
    failed_conversions: int
    total_features_converted: int
    conversion_results: List[ConversionResult]

class AugustusCoordinateConverter:
    """
    Augustus坐标转换器
    
    将Augustus重新预测结果从相对坐标转换为基因组绝对坐标
    """
    
    def __init__(self, augustus_dir: str = "results/augustus_nlr_candidates"):
        """
        初始化转换器
        
        Args:
            augustus_dir: Augustus预测结果目录
        """
        self.augustus_dir = Path(augustus_dir)
        self.logger = logging.getLogger(__name__)
        
        # 验证输入目录
        if not self.augustus_dir.exists():
            raise FileNotFoundError(f"Augustus directory not found: {augustus_dir}")
    
    def extract_region_from_filename(self, filename: str) -> RegionInfo:
        """
        从文件名提取区域信息
        
        Args:
            filename: GFF3文件名 (例如: Chr1_nlr1.gff3)
            
        Returns:
            RegionInfo: 解析的区域信息
        """
        # 读取文件获取区域信息
        filepath = self.augustus_dir / filename
        
        try:
            with open(filepath, 'r') as f:
                for line in f:
                    line = line.strip()
                    
                    # 查找以##sequence-region开头的行
                    if line.startswith('##sequence-region'):
                        # 格式: ##sequence-region Chr1:687993-694473 1 6481
                        parts = line.split()
                        if len(parts) >= 2:
                            region_str = parts[1]
                            return RegionInfo.from_string(region_str)
                    
                    # 查找Augustus预测区域信息注释行
                    # 格式: # ----- prediction on sequence number 1 (length = 7876, name = Chr10:1640748-1648623) -----
                    elif 'prediction on sequence number' in line and 'name =' in line:
                        # 提取name部分
                        match = re.search(r'name\s*=\s*([^)]+)', line)
                        if match:
                            region_str = match.group(1).strip()
                            return RegionInfo.from_string(region_str)
                    
                    # 如果到达数据行（非注释行），停止查找
                    elif not line.startswith('#') and line:
                        break
            
            # 如果没有找到region信息，尝试从特征行推断
            self.logger.warning(f"No region info found in comments for {filename}, attempting to parse from feature lines")
            
            with open(filepath, 'r') as f:
                for line in f:
                    line = line.strip()
                    # 跳过注释行
                    if line.startswith('#') or not line:
                        continue
                    
                    # 解析第一个特征行
                    fields = line.split('\t')
                    if len(fields) >= 9:
                        seqid = fields[0]  # 应该包含区域信息，如 Chr10:1640748-1648623
                        try:
                            return RegionInfo.from_string(seqid)
                        except ValueError:
                            # 如果seqid不是区域格式，继续查找
                            continue
            
            # 最后尝试从文件名推断
            self.logger.error(f"Cannot extract region information from {filename}")
            raise ValueError(f"Cannot determine region coordinates for {filename}")
            
        except Exception as e:
            self.logger.error(f"Error extracting region from {filename}: {e}")
            raise
    
    def convert_coordinates(self, relative_start: int, relative_end: int, 
                          region_info: RegionInfo) -> Tuple[int, int]:
        """
        将相对坐标转换为绝对坐标
        
        Args:
            relative_start: 相对起始坐标
            relative_end: 相对结束坐标
            region_info: 区域信息
            
        Returns:
            Tuple[int, int]: (绝对起始坐标, 绝对结束坐标)
        """
        # 转换公式: absolute = region_start + relative - 1
        absolute_start = region_info.start + relative_start - 1
        absolute_end = region_info.start + relative_end - 1
        
        return absolute_start, absolute_end
    
    def validate_coordinates(self, start: int, end: int, 
                           region_info: RegionInfo) -> bool:
        """
        验证转换后坐标的有效性
        
        Args:
            start: 绝对起始坐标
            end: 绝对结束坐标
            region_info: 区域信息
            
        Returns:
            bool: 坐标是否有效
        """
        # 检查坐标是否在区域范围内
        if start < region_info.start or end > region_info.end:
            return False
            
        # 检查起始坐标是否小于结束坐标
        if start > end:
            return False
            
        return True
    
    def convert_gff3_file(self, input_file: str, output_file: str) -> ConversionResult:
        """
        转换单个GFF3文件的坐标
        
        Args:
            input_file: 输入GFF3文件路径
            output_file: 输出GFF3文件路径
            
        Returns:
            ConversionResult: 转换结果
        """
        input_path = Path(input_file)
        output_path = Path(output_file)
        
        try:
            # 提取区域信息
            region_info = self.extract_region_from_filename(input_path.name)
            self.logger.info(f"Converting {input_file}: region {region_info}")
            
            # 确保输出目录存在
            output_path.parent.mkdir(parents=True, exist_ok=True)
            
            features_converted = 0
            
            with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
                for line in infile:
                    line = line.strip()
                    
                    # 保持注释行不变（除了sequence-region行）
                    if line.startswith('#'):
                        if line.startswith('##sequence-region'):
                            # 更新sequence-region行使用绝对坐标
                            parts = line.split()
                            if len(parts) >= 4:
                                # 格式: ##sequence-region Chr1:687993-694473 1 6481
                                # 转换为: ##sequence-region Chr1 687993 694473
                                chromosome = region_info.chromosome
                                start = region_info.start
                                end = region_info.end
                                outfile.write(f"##sequence-region {chromosome} {start} {end}\n")
                            else:
                                outfile.write(line + '\n')
                        else:
                            outfile.write(line + '\n')
                        continue
                    
                    # 处理特征行
                    if line:
                        fields = line.split('\t')
                        if len(fields) >= 9:  # 标准GFF3格式有9个字段
                            # 字段索引: 0=seqid, 1=source, 2=type, 3=start, 4=end, 5=score, 6=strand, 7=phase, 8=attributes
                            
                            # 更新seqid为染色体名
                            fields[0] = region_info.chromosome
                            
                            # 转换坐标
                            try:
                                relative_start = int(fields[3])
                                relative_end = int(fields[4])
                                
                                absolute_start, absolute_end = self.convert_coordinates(
                                    relative_start, relative_end, region_info
                                )
                                
                                # 验证坐标
                                if not self.validate_coordinates(absolute_start, absolute_end, region_info):
                                    self.logger.warning(
                                        f"Invalid coordinates after conversion: {absolute_start}-{absolute_end} "
                                        f"for region {region_info}"
                                    )
                                
                                fields[3] = str(absolute_start)
                                fields[4] = str(absolute_end)
                                
                                # 确保source字段为Augustus
                                fields[1] = "Augustus"
                                
                                features_converted += 1
                                
                            except ValueError as e:
                                self.logger.error(f"Error converting coordinates in line: {line}, error: {e}")
                                continue
                            
                            # 写入转换后的行
                            outfile.write('\t'.join(fields) + '\n')
            
            return ConversionResult(
                input_file=str(input_path),
                output_file=str(output_path),
                region_info=region_info,
                features_converted=features_converted,
                conversion_success=True
            )
            
        except Exception as e:
            error_msg = f"Error converting {input_file}: {e}"
            self.logger.error(error_msg)
            
            return ConversionResult(
                input_file=str(input_path),
                output_file=str(output_path),
                region_info=None,
                features_converted=0,
                conversion_success=False,
                error_message=error_msg
            )
    
    def get_augustus_files(self) -> List[str]:
        """
        获取所有Augustus预测文件
        
        Returns:
            List[str]: GFF3文件路径列表
        """
        pattern = str(self.augustus_dir / "Chr*_nlr*.gff3")
        files = glob.glob(pattern)
        files.sort()  # 按文件名排序
        
        self.logger.info(f"Found {len(files)} Augustus prediction files")
        return files
    
    def batch_convert(self, output_dir: str = "results/augustus_nlr_candidates_converted",
                     max_workers: int = 4) -> ConversionStats:
        """
        批量转换所有Augustus预测文件
        
        Args:
            output_dir: 输出目录
            max_workers: 最大并行工作线程数
            
        Returns:
            ConversionStats: 批量转换统计结果
        """
        # 获取所有文件
        input_files = self.get_augustus_files()
        
        if not input_files:
            self.logger.warning("No Augustus prediction files found")
            return ConversionStats(0, 0, 0, 0, [])
        
        # 创建输出目录
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        # 准备转换任务
        conversion_tasks = []
        for input_file in input_files:
            input_path = Path(input_file)
            output_file = output_path / input_path.name
            conversion_tasks.append((input_file, str(output_file)))
        
        # 批量处理
        results = []
        successful_conversions = 0
        failed_conversions = 0
        total_features = 0
        
        self.logger.info(f"Starting batch conversion of {len(conversion_tasks)} files using {max_workers} workers")
        
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            # 提交所有任务
            future_to_task = {
                executor.submit(self.convert_gff3_file, input_file, output_file): (input_file, output_file)
                for input_file, output_file in conversion_tasks
            }
            
            # 处理完成的任务
            with tqdm(total=len(conversion_tasks), desc="Converting files") as pbar:
                for future in as_completed(future_to_task):
                    input_file, output_file = future_to_task[future]
                    
                    try:
                        result = future.result()
                        results.append(result)
                        
                        if result.conversion_success:
                            successful_conversions += 1
                            total_features += result.features_converted
                        else:
                            failed_conversions += 1
                            
                    except Exception as e:
                        self.logger.error(f"Exception in converting {input_file}: {e}")
                        results.append(ConversionResult(
                            input_file=input_file,
                            output_file=output_file,
                            region_info=None,
                            features_converted=0,
                            conversion_success=False,
                            error_message=str(e)
                        ))
                        failed_conversions += 1
                    
                    pbar.update(1)
        
        # 创建统计结果
        stats = ConversionStats(
            total_files=len(conversion_tasks),
            successful_conversions=successful_conversions,
            failed_conversions=failed_conversions,
            total_features_converted=total_features,
            conversion_results=results
        )
        
        self.logger.info(f"Batch conversion completed: {successful_conversions} successful, {failed_conversions} failed")
        
        return stats
    
    def generate_conversion_report(self, stats: ConversionStats, 
                                 report_file: str = "results/augustus_coordinate_conversion_report.json") -> None:
        """
        生成转换统计报告
        
        Args:
            stats: 转换统计结果
            report_file: 报告文件路径
        """
        report_data = {
            "conversion_summary": {
                "total_files": stats.total_files,
                "successful_conversions": stats.successful_conversions,
                "failed_conversions": stats.failed_conversions,
                "success_rate": stats.successful_conversions / stats.total_files if stats.total_files > 0 else 0,
                "total_features_converted": stats.total_features_converted
            },
            "conversion_details": []
        }
        
        for result in stats.conversion_results:
            detail = {
                "input_file": result.input_file,
                "output_file": result.output_file,
                "conversion_success": result.conversion_success,
                "features_converted": result.features_converted
            }
            
            if result.region_info:
                detail["region_info"] = {
                    "chromosome": result.region_info.chromosome,
                    "start": result.region_info.start,
                    "end": result.region_info.end
                }
            
            if result.error_message:
                detail["error_message"] = result.error_message
                
            report_data["conversion_details"].append(detail)
        
        # 保存报告
        report_path = Path(report_file)
        report_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(report_path, 'w') as f:
            json.dump(report_data, f, indent=2)
        
        self.logger.info(f"Conversion report saved to {report_file}")
    
    def export_conversion_summary(self, stats: ConversionStats,
                                tsv_file: str = "results/augustus_coordinate_conversion_summary.tsv") -> None:
        """
        导出转换摘要为TSV格式
        
        Args:
            stats: 转换统计结果
            tsv_file: TSV文件路径
        """
        tsv_path = Path(tsv_file)
        tsv_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(tsv_path, 'w', newline='') as f:
            writer = csv.writer(f, delimiter='\t')
            
            # 写入标题行
            writer.writerow([
                "input_file", "output_file", "chromosome", "region_start", "region_end",
                "features_converted", "conversion_success", "error_message"
            ])
            
            # 写入数据行
            for result in stats.conversion_results:
                row = [
                    result.input_file,
                    result.output_file,
                    result.region_info.chromosome if result.region_info else "",
                    result.region_info.start if result.region_info else "",
                    result.region_info.end if result.region_info else "",
                    result.features_converted,
                    result.conversion_success,
                    result.error_message or ""
                ]
                writer.writerow(row)
        
        self.logger.info(f"Conversion summary exported to {tsv_file}")

def main():
    """主函数 - 命令行工具入口"""
    import argparse
    
    parser = argparse.ArgumentParser(description="Convert Augustus coordinates from relative to absolute")
    parser.add_argument("--augustus-dir", default="results/augustus_nlr_candidates",
                       help="Augustus prediction results directory")
    parser.add_argument("--output-dir", default="results/augustus_nlr_candidates_converted",
                       help="Output directory for converted files")
    parser.add_argument("--max-workers", type=int, default=4,
                       help="Maximum number of parallel workers")
    parser.add_argument("--report-file", default="results/augustus_coordinate_conversion_report.json",
                       help="Conversion report file")
    parser.add_argument("--summary-file", default="results/augustus_coordinate_conversion_summary.tsv",
                       help="Conversion summary TSV file")
    
    args = parser.parse_args()
    
    # 创建转换器
    converter = AugustusCoordinateConverter(args.augustus_dir)
    
    # 执行批量转换
    stats = converter.batch_convert(args.output_dir, args.max_workers)
    
    # 生成报告
    converter.generate_conversion_report(stats, args.report_file)
    converter.export_conversion_summary(stats, args.summary_file)
    
    # 打印摘要
    print(f"\nConversion Summary:")
    print(f"Total files: {stats.total_files}")
    print(f"Successful conversions: {stats.successful_conversions}")
    print(f"Failed conversions: {stats.failed_conversions}")
    print(f"Success rate: {stats.successful_conversions/stats.total_files*100:.1f}%")
    print(f"Total features converted: {stats.total_features_converted}")
    
    if stats.failed_conversions > 0:
        print(f"\nFailed files:")
        for result in stats.conversion_results:
            if not result.conversion_success:
                print(f"  {result.input_file}: {result.error_message}")

if __name__ == "__main__":
    main() 