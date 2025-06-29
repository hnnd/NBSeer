#!/usr/bin/env python3
"""
Runtime Benchmark Testing Script
For testing NBSeer performance across different genome sizes
"""

import os
import sys
import time
import json
import psutil
import subprocess
import argparse
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
import logging
import yaml
import shutil

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class RuntimeBenchmark:
    def __init__(self, manuscript_data_dir, output_dir, pipeline_script):
        """
        初始化运行时间基准测试
        
        Args:
            manuscript_data_dir: manuscript/data目录路径
            output_dir: 结果输出目录
            pipeline_script: NBS流水线主脚本路径
        """
        self.manuscript_data_dir = Path(manuscript_data_dir)
        self.output_dir = Path(output_dir)
        self.pipeline_script = Path(pipeline_script)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # 测试数据集配置
        self.datasets = {
            'tair10': {
                'name': 'Arabidopsis thaliana',
                'genome': self.manuscript_data_dir / 'genome' / 'tair10.fa',
                'expected_size_mb': 120,
                'description': 'Small genome for quick testing'
            },
            'osa': {
                'name': 'Oryza sativa',
                'genome': self.manuscript_data_dir / 'genome' / 'osa.fa', 
                'expected_size_mb': 380,
                'description': 'Medium-sized genome'
            },
            'CM334': {
                'name': 'Capsicum annuum',
                'genome': self.manuscript_data_dir / 'genome' / 'CM334.fa',
                'expected_size_mb': 3000,
                'description': 'Large genome for scalability testing'
            }
        }
        
        self.protein_db = self.manuscript_data_dir / 'prgdb' / 'prg_nbs.fasta'
        self.results = []
        
        # 配置文件路径
        self.config_template = Path('config/default.yaml')
        self.temp_config_dir = self.output_dir / 'temp_configs'
        self.temp_config_dir.mkdir(parents=True, exist_ok=True)
        
    def get_genome_size(self, genome_path):
        """获取基因组文件大小（MB）"""
        try:
            size_bytes = genome_path.stat().st_size
            size_mb = size_bytes / (1024 * 1024)
            return size_mb
        except Exception as e:
            logger.error(f"Error getting genome size for {genome_path}: {e}")
            return 0
    
    def monitor_resources(self, process):
        """监控进程资源使用情况"""
        max_memory_mb = 0
        cpu_times = []
        memory_usage = []
        timestamps = []
        
        try:
            ps_process = psutil.Process(process.pid)
            start_time = time.time()
            
            while process.poll() is None:
                try:
                    # 获取内存使用情况 (MB)
                    memory_info = ps_process.memory_info()
                    current_memory_mb = memory_info.rss / (1024 * 1024)
                    max_memory_mb = max(max_memory_mb, current_memory_mb)
                    
                    # 获取CPU使用情况
                    cpu_percent = ps_process.cpu_percent()
                    
                    # 记录时间序列数据
                    current_time = time.time() - start_time
                    timestamps.append(current_time)
                    memory_usage.append(current_memory_mb)
                    cpu_times.append(cpu_percent)
                    
                    time.sleep(1)  # 每秒采样一次
                    
                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    break
                    
        except Exception as e:
            logger.warning(f"Resource monitoring error: {e}")
        
        return {
            'max_memory_mb': max_memory_mb,
            'avg_cpu_percent': sum(cpu_times) / len(cpu_times) if cpu_times else 0,
            'memory_timeline': memory_usage,
            'cpu_timeline': cpu_times,
            'timestamps': timestamps
        }
    
    def create_threaded_config(self, num_threads, config_name):
        """
        创建一个临时配置文件,将所有线程/CPU参数设置为指定值
        
        Args:
            num_threads: 要设置的线程数
            config_name: 配置文件名称
        
        Returns:
            临时配置文件路径
        """
        if not self.config_template.exists():
            logger.error(f"Configuration template not found: {self.config_template}")
            return None
        
        # 读取配置模板
        with open(self.config_template, 'r', encoding='utf-8') as f:
            config = yaml.safe_load(f)
        
        # 更新所有线程/CPU相关参数
        if 'tools' in config:
            # NLR-Annotator线程数
            if 'nlr_annotator' in config['tools'] and 'parameters' in config['tools']['nlr_annotator']:
                config['tools']['nlr_annotator']['parameters']['threads'] = num_threads
                logger.info(f"Set NLR-Annotator threads to {num_threads}")
            
            # Miniprot线程数
            if 'miniprot' in config['tools'] and 'parameters' in config['tools']['miniprot']:
                config['tools']['miniprot']['parameters']['threads'] = num_threads
                logger.info(f"Set Miniprot threads to {num_threads}")
            
            # Augustus训练CPU数
            if 'augustus' in config['tools'] and 'training' in config['tools']['augustus']:
                if 'miniprot_training' in config['tools']['augustus']['training']:
                    config['tools']['augustus']['training']['miniprot_training']['cpus'] = num_threads
                    logger.info(f"Set Augustus training CPUs to {num_threads}")
            
            # EVidenceModeler CPU数
            if 'evm' in config['tools'] and 'parameters' in config['tools']['evm']:
                config['tools']['evm']['parameters']['cpu'] = num_threads
                logger.info(f"Set EVidenceModeler CPUs to {num_threads}")
        
        # 更新流水线并行处理参数
        if 'pipeline' in config and 'parallel' in config['pipeline']:
            config['pipeline']['parallel']['max_workers'] = min(num_threads, 8)  # 限制最大工作线程数
            logger.info(f"Set pipeline max_workers to {min(num_threads, 8)}")
        
        # 保存临时配置文件
        temp_config_file = self.temp_config_dir / f'{config_name}.yaml'
        with open(temp_config_file, 'w', encoding='utf-8') as f:
            yaml.dump(config, f, default_flow_style=False, allow_unicode=True)
        
        logger.info(f"Created threaded config: {temp_config_file}")
        return temp_config_file
    
    def run_pipeline_benchmark(self, dataset_key, num_threads=8, repeat=3):
        """运行单个数据集的基准测试"""
        dataset = self.datasets[dataset_key]
        genome_path = dataset['genome']
        
        if not genome_path.exists():
            logger.error(f"Genome file not found: {genome_path}")
            return None
        
        if not self.protein_db.exists():
            logger.error(f"Protein database not found: {self.protein_db}")
            return None
        
        logger.info(f"Starting benchmark for {dataset['name']} ({dataset_key})")
        
        # 获取实际基因组大小
        actual_size_mb = self.get_genome_size(genome_path)
        
        # 多次运行取平均值
        run_results = []
        
        for run_idx in range(repeat):
            logger.info(f"  Run {run_idx + 1}/{repeat}")
            
            # 创建临时输出目录
            temp_output = self.output_dir / f"temp_{dataset_key}_run{run_idx}"
            temp_output.mkdir(exist_ok=True)
            
            # 创建线程配置文件
            config_name = f"{dataset_key}_threads{num_threads}_run{run_idx}"
            temp_config = self.create_threaded_config(num_threads, config_name)
            
            if temp_config is None:
                logger.error(f"Failed to create config for {dataset_key} run {run_idx}")
                continue
            
            # 构建命令
            cmd = [
                'python', '-m', 'src.nbseer.main',
                '--config', str(temp_config),
                '--genome', str(genome_path),
                '--proteins', str(self.protein_db),
                '--output', str(temp_output),
                '--threads', str(num_threads),
                '--enable-training',
                '--training-species-name', 'nbs_training_species'
            ]
            
            # 运行流水线并监控资源
            start_time = time.time()
            
            try:
                process = subprocess.Popen(
                    cmd,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    cwd=self.pipeline_script.parent.parent.parent  # 项目根目录
                )
                
                # 监控资源使用
                resource_stats = self.monitor_resources(process)
                
                # 等待完成
                stdout, stderr = process.communicate()
                end_time = time.time()
                
                runtime_seconds = end_time - start_time
                runtime_minutes = runtime_seconds / 60
                runtime_hours = runtime_seconds / 3600
                
                # 检查运行结果
                success = process.returncode == 0
                
                if not success:
                    logger.error(f"Pipeline failed for {dataset_key} run {run_idx}")
                    logger.error(f"Error: {stderr.decode()}")
                
                run_result = {
                    'dataset': dataset_key,
                    'run_index': run_idx,
                    'runtime_seconds': runtime_seconds,
                    'runtime_minutes': runtime_minutes,
                    'runtime_hours': runtime_hours,
                    'success': success,
                    'max_memory_mb': resource_stats['max_memory_mb'],
                    'avg_cpu_percent': resource_stats['avg_cpu_percent'],
                    'num_threads': num_threads,
                    'genome_size_mb': actual_size_mb,
                    'timestamp': datetime.now().isoformat()
                }
                
                run_results.append(run_result)
                
                # 清理临时文件
                subprocess.run(['rm', '-rf', str(temp_output)], check=False)
                
                # 清理临时配置文件
                if temp_config and temp_config.exists():
                    temp_config.unlink()
                
            except Exception as e:
                logger.error(f"Error running pipeline: {e}")
                run_results.append({
                    'dataset': dataset_key,
                    'run_index': run_idx,
                    'success': False,
                    'error': str(e),
                    'timestamp': datetime.now().isoformat()
                })
        
        # 计算平均值（仅成功的运行）
        successful_runs = [r for r in run_results if r.get('success', False)]
        
        if successful_runs:
            avg_result = {
                'dataset': dataset_key,
                'species_name': dataset['name'],
                'genome_size_mb': actual_size_mb,
                'expected_size_mb': dataset['expected_size_mb'],
                'successful_runs': len(successful_runs),
                'total_runs': repeat,
                'avg_runtime_seconds': sum(r['runtime_seconds'] for r in successful_runs) / len(successful_runs),
                'avg_runtime_minutes': sum(r['runtime_minutes'] for r in successful_runs) / len(successful_runs),
                'avg_runtime_hours': sum(r['runtime_hours'] for r in successful_runs) / len(successful_runs),
                'std_runtime_seconds': pd.Series([r['runtime_seconds'] for r in successful_runs]).std(),
                'avg_max_memory_mb': sum(r['max_memory_mb'] for r in successful_runs) / len(successful_runs),
                'std_max_memory_mb': pd.Series([r['max_memory_mb'] for r in successful_runs]).std(),
                'avg_cpu_percent': sum(r['avg_cpu_percent'] for r in successful_runs) / len(successful_runs),
                'num_threads': num_threads,
                'benchmark_timestamp': datetime.now().isoformat()
            }
            
            return avg_result, run_results
        else:
            logger.error(f"No successful runs for {dataset_key}")
            return None, run_results
    
    def run_parallelization_test(self, dataset_key='tair10', thread_counts=[1, 2, 4, 8, 16]):
        """测试并行化效果"""
        logger.info(f"Running parallelization test on {dataset_key}")
        
        parallel_results = []
        baseline_time = None
        
        for threads in thread_counts:
            logger.info(f"Testing with {threads} threads")
            
            result, _ = self.run_pipeline_benchmark(dataset_key, num_threads=threads, repeat=2)
            
            if result:
                if baseline_time is None and threads == 1:
                    baseline_time = result['avg_runtime_seconds']
                
                speedup = baseline_time / result['avg_runtime_seconds'] if baseline_time else 1.0
                efficiency = speedup / threads if threads > 0 else 0
                
                parallel_result = {
                    'threads': threads,
                    'runtime_seconds': result['avg_runtime_seconds'],
                    'speedup': speedup,
                    'efficiency': efficiency,
                    'memory_mb': result['avg_max_memory_mb']
                }
                
                parallel_results.append(parallel_result)
        
        return parallel_results
    
    def generate_benchmark_report(self, results, parallel_results=None):
        """生成基准测试报告"""
        # 创建DataFrame
        df = pd.DataFrame(results)
        
        # 生成汇总报告
        report = {
            'benchmark_summary': {
                'total_datasets': len(results),
                'total_successful': len([r for r in results if r.get('successful_runs', 0) > 0]),
                'benchmark_date': datetime.now().isoformat()
            },
            'performance_metrics': results,
            'parallelization_analysis': parallel_results if parallel_results else []
        }
        
        # 保存JSON报告
        report_file = self.output_dir / 'runtime_benchmark_report.json'
        with open(report_file, 'w') as f:
            json.dump(report, f, indent=2)
        
        # 保存CSV数据
        csv_file = self.output_dir / 'runtime_benchmark_data.csv'
        df.to_csv(csv_file, index=False)
        
        logger.info(f"Benchmark report saved to {report_file}")
        logger.info(f"Benchmark data saved to {csv_file}")
        
        return report
    
    def create_visualizations(self, results, parallel_results=None):
        """创建可视化图表"""
        # 设置图表样式
        plt.style.use('default')
        sns.set_palette("husl")
        
        # 图1: 运行时间 vs 基因组大小
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        # 提取数据
        successful_results = [r for r in results if r.get('successful_runs', 0) > 0]
        genome_sizes = [r['genome_size_mb'] for r in successful_results]
        runtimes_hours = [r['avg_runtime_hours'] for r in successful_results]
        memory_usage = [r['avg_max_memory_mb'] for r in successful_results]
        species_names = [r['species_name'] for r in successful_results]
        
        # 子图1: 运行时间 vs 基因组大小
        axes[0,0].scatter(genome_sizes, runtimes_hours, s=100, alpha=0.7)
        for i, name in enumerate(species_names):
            axes[0,0].annotate(name.split()[0], (genome_sizes[i], runtimes_hours[i]), 
                              xytext=(5, 5), textcoords='offset points')
        axes[0,0].set_xlabel('Genome Size (MB)')
        axes[0,0].set_ylabel('Runtime (Hours)')
        axes[0,0].set_title('Runtime vs Genome Size')
        axes[0,0].grid(True, alpha=0.3)
        
        # 子图2: 内存使用 vs 基因组大小
        axes[0,1].scatter(genome_sizes, memory_usage, s=100, alpha=0.7, color='orange')
        for i, name in enumerate(species_names):
            axes[0,1].annotate(name.split()[0], (genome_sizes[i], memory_usage[i]),
                              xytext=(5, 5), textcoords='offset points')
        axes[0,1].set_xlabel('Genome Size (MB)')
        axes[0,1].set_ylabel('Peak Memory Usage (MB)')
        axes[0,1].set_title('Memory Usage vs Genome Size')
        axes[0,1].grid(True, alpha=0.3)
        
        # 子图3: 并行化效果（如果有数据）
        if parallel_results:
            threads = [r['threads'] for r in parallel_results]
            speedups = [r['speedup'] for r in parallel_results]
            efficiencies = [r['efficiency'] for r in parallel_results]
            
            ax3 = axes[1,0]
            ax3.plot(threads, speedups, 'o-', label='Actual Speedup', linewidth=2)
            ax3.plot(threads, threads, '--', label='Ideal Speedup', alpha=0.5)
            ax3.set_xlabel('Number of Threads')
            ax3.set_ylabel('Speedup Factor')
            ax3.set_title('Parallelization Speedup')
            ax3.legend()
            ax3.grid(True, alpha=0.3)
            
            # 效率图
            ax4 = axes[1,1]
            ax4.plot(threads, efficiencies, 'o-', color='red', linewidth=2)
            ax4.set_xlabel('Number of Threads')
            ax4.set_ylabel('Parallel Efficiency')
            ax4.set_title('Parallelization Efficiency')
            ax4.grid(True, alpha=0.3)
            ax4.set_ylim(0, 1.1)
        else:
            axes[1,0].text(0.5, 0.5, 'No parallelization data', 
                          ha='center', va='center', transform=axes[1,0].transAxes)
            axes[1,1].text(0.5, 0.5, 'No efficiency data',
                          ha='center', va='center', transform=axes[1,1].transAxes)
        
        plt.tight_layout()
        
        # 保存图表
        plot_file = self.output_dir / 'runtime_benchmark_plots.png'
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"Benchmark plots saved to {plot_file}")
    
    def cleanup_temp_configs(self):
        """清理临时配置文件"""
        try:
            if self.temp_config_dir.exists():
                shutil.rmtree(self.temp_config_dir)
                logger.info("Cleaned up temporary configuration files")
        except Exception as e:
            logger.warning(f"Error cleaning up temp configs: {e}")
    
    def run_full_benchmark(self, include_parallelization=True):
        """运行完整的基准测试"""
        logger.info("Starting full runtime benchmark analysis")
        
        try:
            all_results = []
            
            # 运行所有数据集的基准测试
            for dataset_key in self.datasets.keys():
                logger.info(f"Processing dataset: {dataset_key}")
                
                result, detailed_results = self.run_pipeline_benchmark(dataset_key, num_threads=8, repeat=3)
                
                if result:
                    all_results.append(result)
                    
                    # 保存详细结果
                    detailed_file = self.output_dir / f'detailed_results_{dataset_key}.json'
                    with open(detailed_file, 'w') as f:
                        json.dump(detailed_results, f, indent=2)
            
            # 并行化测试（在最小的数据集上）
            parallel_results = None
            if include_parallelization and all_results:
                logger.info("Starting parallelization test with properly configured thread counts")
                parallel_results = self.run_parallelization_test('tair10', [1, 2, 4, 8, 16])
            
            # 生成报告和可视化
            report = self.generate_benchmark_report(all_results, parallel_results)
            self.create_visualizations(all_results, parallel_results)
            
            logger.info("Runtime benchmark analysis completed")
            return report
            
        finally:
            # 确保清理临时配置文件
            self.cleanup_temp_configs()

def main():
    parser = argparse.ArgumentParser(description='NBS-Pipeline Runtime Benchmark Analysis')
    parser.add_argument('--manuscript-data', required=True, help='Path to manuscript/data directory')
    parser.add_argument('--output', required=True, help='Output directory for results')
    parser.add_argument('--pipeline-script', required=True, help='Path to NBS pipeline main script')
    parser.add_argument('--skip-parallelization', action='store_true', 
                       help='Skip parallelization testing')
    
    args = parser.parse_args()
    
    # 创建基准测试对象
    benchmark = RuntimeBenchmark(
        manuscript_data_dir=args.manuscript_data,
        output_dir=args.output,
        pipeline_script=args.pipeline_script
    )
    
    # 运行基准测试
    report = benchmark.run_full_benchmark(include_parallelization=not args.skip_parallelization)
    
    print("Benchmark completed successfully!")
    print(f"Results saved to: {args.output}")

if __name__ == '__main__':
    main()