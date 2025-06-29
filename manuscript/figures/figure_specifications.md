# 图表制作指南和规格要求

## Bioinformatics期刊图表要求

### 技术规格
- **分辨率**: 350 DPI (位图) / 1200 DPI (线条图)
- **格式**: TIFF, EPS, 或高质量PDF
- **色彩**: 支持彩色，但需考虑色盲友好
- **字体**: 清晰可读，最小10pt
- **尺寸**: 单栏宽度85mm，双栏宽度175mm

---

## Figure 1: Pipeline Workflow架构图

### 设计要求
- **类型**: 流程图/架构图
- **尺寸**: 双栏宽度 (175mm)
- **布局**: 垂直流程，从上到下

### 内容元素
```
输入数据 (顶部)
├── 基因组序列 (FASTA)
└── 参考蛋白序列 (FASTA)
        ↓
阶段1: 输入验证 (Input Validation)
├── 序列格式检查
├── 完整性验证  
└── 质量评估
        ↓
阶段2: NLR基因定位 (NLR Localization)
├── NLR-Annotator扫描
├── 保守域识别
└── 候选区域提取
        ↓
阶段3: 蛋白质比对 (Protein Alignment)  
├── miniprot比对
├── 剪接位点预测
└── 基因模型生成
        ↓
阶段4: 基因结构预测 (Gene Prediction)
├── Augustus物种训练
├── 从头基因预测
└── 结构优化
        ↓
阶段5: 证据整合 (Evidence Integration)
├── EVidenceModeler整合
├── 权重优化
└── 共识模型生成
        ↓
输出结果 (底部)
├── 最终注释 (GFF3)
├── 基因序列 (FASTA)
└── 质量报告 (JSON/HTML)
```

### 颜色方案
- 输入/输出: 浅蓝色 (#E3F2FD)
- 处理阶段: 渐变绿色 (#C8E6C9 to #4CAF50)
- 工具名称: 深灰色 (#424242)
- 连接线: 深蓝色 (#1976D2)

### 制作脚本
```python
# figures/scripts/create_figure1_workflow.py
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyBboxPatch
```

---

## Figure 2: 注释质量评估

### 设计要求
- **类型**: 组合图表 (3个子图)
- **尺寸**: 双栏宽度 (175mm)
- **布局**: 2×2网格，空余位置放图例

### 子图A: 准确性对比柱状图
```python
# 数据结构示例
tools = ['NBS-Pipeline', 'MAKER2', 'Augustus', 'GeneMark-ES', 'BRAKER2']
sensitivity = [0.85, 0.72, 0.68, 0.65, 0.78]  # 示例数据
specificity = [0.91, 0.81, 0.76, 0.73, 0.84]  # 示例数据
f1_score = [0.88, 0.76, 0.72, 0.69, 0.81]     # 示例数据
```

### 子图B: 基因结构准确性箱线图
```python
# 外显子边界准确性分布
# Y轴: 准确性百分比 (0-100%)
# X轴: 不同工具
# 数据: 每个工具在多个基因上的表现分布
```

### 子图C: 多物种性能热图
```python
# 数据矩阵: 物种 × 评估指标
species = ['Arabidopsis', 'Rice', 'Pepper', 'Tomato', 'Soybean']
metrics = ['Sensitivity', 'Specificity', 'F1-Score', 'Domain Accuracy']
# 颜色映射: 0-1区间，白色到深绿色
```

### 制作脚本位置
- `figures/scripts/create_figure2_quality.py`

---

## Figure 3: 性能基准测试

### 设计要求
- **类型**: 散点图和折线图组合
- **尺寸**: 双栏宽度 (175mm)  
- **布局**: 2×2网格

### 子图A: 运行时间vs基因组大小
```python
genome_sizes = [120, 380, 900, 2100]  # MB
runtime_hours = [0.5, 2.1, 8.3, 24.7]  # 示例数据
# 对数坐标，拟合趋势线
```

### 子图B: 内存使用情况
```python
# 峰值内存使用 (GB) vs 基因组大小 (MB)
# 包含error bars显示变异范围
```

### 子图C: 并行化效果
```python
cpu_cores = [1, 2, 4, 8, 16, 32]
speedup = [1.0, 1.8, 3.2, 5.1, 7.3, 8.9]  # 实际数据待测试
# 理论线性加速 vs 实际测量结果
```

### 子图D: 各阶段耗时分布
```python
# 堆积柱状图或饼图
stages = ['NLR Localization', 'Protein Alignment', 'Gene Prediction', 'Evidence Integration']
time_percentages = [15, 35, 40, 10]  # 百分比
```

### 制作脚本位置
- `figures/scripts/create_figure3_performance.py`

---

## Table 1: 数据集汇总和结果统计

### 表格结构
```
| Species | Genome Size (Mb) | Known NBS | Predicted | Validated | Sensitivity | Specificity | F1-Score |
|---------|------------------|-----------|-----------|-----------|-------------|-------------|----------|
| A.thaliana | 120 | 165 | 172 | 158 | 0.86 | 0.92 | 0.89 |
| O.sativa | 380 | 341 | 356 | 325 | 0.84 | 0.91 | 0.87 |
| C.annuum | 900 | 278 | 289 | 267 | 0.88 | 0.92 | 0.90 |
| S.lycopersicum | 780 | 298 | 312 | 287 | 0.85 | 0.92 | 0.88 |
| G.max | 1100 | 331 | 347 | 318 | 0.83 | 0.92 | 0.87 |
```

### 格式要求
- 使用LaTeX表格格式
- 数值对齐，保留2-3位小数
- 添加表注解释缩写和方法

---

## Supplementary Figures

### Supplementary Figure S1: 详细案例研究
- 展示具体基因的注释过程
- 多证据对比和整合结果
- 与已知注释的详细比较

### Supplementary Figure S2: 参数敏感性分析
- 关键参数对结果影响的热图
- 参数优化建议的可视化

### Supplementary Figure S3: 质量控制统计
- 假阳性/假阴性案例分析
- 质量分数分布直方图

---

## 图表制作工作流

### 1. 数据准备阶段
```bash
# 创建数据准备脚本
mkdir -p figures/data_prep
# 从分析结果提取图表数据
python data_analysis/extract_figure_data.py
```

### 2. 图表生成阶段
```bash
# 批量生成所有图表
python figures/scripts/generate_all_figures.py
```

### 3. 质量检查阶段
```bash
# 检查图表规格和质量
python figures/scripts/validate_figures.py
```

### 4. 格式转换阶段
```bash
# 转换为期刊要求格式
python figures/scripts/convert_to_journal_format.py
```

---

## 制作时间表

| 图表 | 数据依赖 | 预计制作时间 | 负责人 |
|------|----------|--------------|--------|
| Figure 1 | 无 | 1天 | [待指定] |
| Figure 2 | 准确性评估完成 | 2天 | [待指定] |
| Figure 3 | 性能测试完成 | 2天 | [待指定] |
| Table 1 | 所有验证完成 | 1天 | [待指定] |
| Supp Figs | 深度分析完成 | 2天 | [待指定] |

## 质量检查清单

- [ ] 分辨率符合期刊要求
- [ ] 字体大小清晰可读
- [ ] 颜色对色盲友好
- [ ] 图例完整清晰
- [ ] 坐标轴标签正确
- [ ] 统计显著性标注
- [ ] 文件格式正确