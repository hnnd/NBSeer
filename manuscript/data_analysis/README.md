# NBSeer Manuscript Data Analysis Script Suite

## Overview

This is a complete data analysis script suite prepared for NBSeer project publication. This suite contains all necessary tools for performance benchmarking, accuracy assessment, statistical analysis, and manuscript figure generation.

## 脚本组织结构

```
data_analysis/
├── README.md                           # 本文件
├── run_analysis_pipeline.sh            # 主控脚本
├── analysis_tasks.md                   # 详细任务清单
├── benchmark_scripts/                  # 性能基准测试
│   └── runtime_analysis.py
├── accuracy_assessment/                # 准确性评估
│   └── gene_structure_eval.py
├── statistical_analysis/               # 统计分析
│   └── gene_features_stats.py
├── visualization/                      # 可视化工具
├── reports/                           # 报告生成
│   └── manuscript_data_generator.py
└── results/                           # 分析结果输出
```

## 快速开始

### 1. 环境准备

确保您的系统已安装必要的依赖：

```bash
# Python 依赖
pip install pandas numpy matplotlib seaborn scipy scikit-learn
pip install biopython gffutils psutil intervaltree

# 系统工具
# 确保有足够的计算资源（建议64核心，128GB内存）
```

### 2. 数据准备

确保您的数据目录结构如下：

```
manuscript/data/
├── genome/
│   ├── tair10.fa          # 拟南芥基因组
│   ├── osa.fa             # 水稻基因组
│   └── CM334.fa           # 辣椒基因组
└── prgdb/
    └── prg_nbs.fasta      # NBS蛋白参考数据库
```

### 3. 运行完整分析

```bash
# 进入分析目录
cd manuscript/data_analysis

# 运行完整分析流程
./run_analysis_pipeline.sh

# 或者查看帮助信息
./run_analysis_pipeline.sh --help
```

### 4. 快速模式运行

如果您想要快速测试，可以使用：

```bash
# 快速模式（减少重复次数，跳过并行化测试）
./run_analysis_pipeline.sh --quick

# 跳过耗时的性能测试
./run_analysis_pipeline.sh --skip-runtime

# 只生成论文图表
./run_analysis_pipeline.sh --skip-runtime --skip-accuracy --skip-features
```

## 详细脚本说明

### 1. 性能基准测试 (`runtime_analysis.py`)

**功能**: 测试NBS-Pipeline在不同基因组大小上的运行性能

**输入**:
- 基因组FASTA文件
- 蛋白质数据库
- NBS-Pipeline主脚本

**输出**:
- 运行时间统计
- 内存使用分析
- 并行化效果评估
- 性能可视化图表

**使用示例**:
```bash
python runtime_analysis.py \
    --manuscript-data ../data \
    --output ./results/runtime_benchmarks \
    --pipeline-script ../../src/nbs_annotation/main.py
```

### 2. 准确性评估 (`gene_structure_eval.py`)

**功能**: 评估NBS-Pipeline注释结果与参考注释的准确性

**输入**:
- 参考注释GFF文件
- 预测注释GFF文件
- 数据集名称

**输出**:
- 基因水平准确性指标
- 外显子结构准确性
- 核苷酸水平准确性
- 详细比较可视化

**使用示例**:
```bash
python gene_structure_eval.py \
    --reference-gff reference_annotations/tair10_nbs_genes.gff3 \
    --predicted-gff ../results_tair10/evidence_integration/final_annotations.gff \
    --dataset-name tair10 \
    --output ./results/accuracy_assessment
```

### 3. 统计分析 (`gene_features_stats.py`)

**功能**: 分析NBS基因的结构特征和分布模式

**输入**:
- 多个GFF注释文件
- 对应的基因组序列文件（可选）
- 数据集名称

**输出**:
- 基因长度分布分析
- 外显子-内含子结构统计
- GC含量分析
- 聚类分析结果

**使用示例**:
```bash
python gene_features_stats.py \
    --gff-files ../results_tair10/evidence_integration/final_annotations.gff \
                ../results_osa/evidence_integration/final_annotations.gff \
    --genome-files ../data/genome/tair10.fa ../data/genome/osa.fa \
    --dataset-names Arabidopsis Rice \
    --output ./results/gene_features
```

### 4. 论文数据生成 (`manuscript_data_generator.py`)

**功能**: 整合所有分析结果，生成符合期刊要求的图表和数据表

**输入**:
- 所有分析结果目录
- 图表输出目录
- 表格输出目录

**输出**:
- Figure 1: Pipeline Workflow (TIFF, 300 DPI)
- Figure 2: Accuracy Assessment (TIFF, 300 DPI)
- Figure 3: Performance Benchmarks (TIFF, 300 DPI)
- Table 1: Dataset Summary (CSV + LaTeX)
- 补充材料图表

**使用示例**:
```bash
python manuscript_data_generator.py \
    --analysis-results ./results \
    --figures-output ../figures/processed \
    --tables-output ../tables/formatted
```

## 输出文件说明

### 图表文件
- **格式**: TIFF, 300 DPI (符合Bioinformatics期刊要求)
- **位置**: `manuscript/figures/processed/`
- **命名**: `Figure[N]_[Description].tiff`

### 数据表文件
- **格式**: CSV (数据) + LaTeX (排版)
- **位置**: `manuscript/tables/formatted/`
- **命名**: `Table[N]_[Description].[csv|tex]`

### 分析结果
- **位置**: `manuscript/data_analysis/results/`
- **内容**: 详细的分析数据、统计结果、中间文件

## 故障排除

### 常见问题

1. **内存不足错误**
   ```bash
   # 减少并行线程数
   ./run_analysis_pipeline.sh --threads 4
   ```

2. **缺少参考注释文件**
   ```bash
   # 跳过准确性评估
   ./run_analysis_pipeline.sh --skip-accuracy
   ```

3. **NBS-Pipeline未安装**
   ```bash
   # 指定Pipeline脚本路径
   ./run_analysis_pipeline.sh -p /path/to/nbs/pipeline/main.py
   ```

4. **Python依赖缺失**
   ```bash
   # 安装必要的包
   pip install -r requirements.txt  # 如果有requirements文件
   # 或者手动安装
   pip install pandas numpy matplotlib seaborn scipy scikit-learn biopython gffutils psutil intervaltree
   ```

### 调试模式

如果遇到问题，可以单独运行各个脚本进行调试：

```bash
# 单独测试性能分析
python benchmark_scripts/runtime_analysis.py --help

# 单独测试准确性评估
python accuracy_assessment/gene_structure_eval.py --help

# 单独测试特征分析
python statistical_analysis/gene_features_stats.py --help
```

## 自定义分析

如果您需要修改分析参数或添加新的分析功能：

1. **修改脚本参数**: 直接编辑相应的Python脚本
2. **添加新数据集**: 更新`run_analysis_pipeline.sh`中的数据集配置
3. **自定义可视化**: 修改`manuscript_data_generator.py`中的图表生成函数

## 性能要求

### 推荐硬件配置
- **CPU**: 64核心或以上
- **内存**: 128GB或以上
- **存储**: 2TB高速SSD
- **预计运行时间**: 8-24小时（取决于数据集大小）

### 最小硬件配置
- **CPU**: 16核心
- **内存**: 32GB
- **存储**: 500GB
- **预计运行时间**: 24-48小时

## 技术支持

如果您在使用过程中遇到问题：

1. 检查日志文件（位于结果目录中）
2. 确认输入数据格式正确
3. 验证所有依赖已正确安装
4. 查看具体错误信息并根据上述故障排除指南解决

## 版本信息

- **脚本版本**: 1.0.0
- **兼容的NBS-Pipeline版本**: 0.1.0+
- **Python版本要求**: 3.11+
- **最后更新**: 2024-12-27

---

**注意**: 这些脚本是为学术论文发表准备的，确保您有足够的计算资源和时间来完成所有分析。建议先在小数据集上测试，确认流程正常后再运行完整分析。