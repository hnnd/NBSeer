# NBS-Pipeline 论文数据分析汇总报告

**生成时间**: 2025-06-29 11:08:24
**分析脚本版本**: 1.0.0

## 分析配置

- **数据目录**: /home/wangys/data/work/nbs/manuscript/data
- **输出目录**: /home/wangys/data/work/nbs/manuscript/data_analysis/results
- **并行线程数**: 8
- **快速模式**: false

## 执行的分析阶段

- **性能基准测试**: ✓ 完成
- **准确性评估**: ✓ 完成
- **特征统计分析**: ✓ 完成
- **论文数据生成**: ✓ 完成

## 输出文件结构

```
/home/wangys/data/work/nbs/manuscript/data_analysis/results/
├── runtime_benchmarks/          # 性能基准测试结果
├── accuracy_assessment/         # 准确性评估结果
├── gene_features/              # 基因特征分析结果
└── analysis_summary.md         # 本报告

/home/wangys/data/work/nbs/manuscript/figures/processed/        # 论文图表
├── Figure1_Pipeline_Workflow.svg
├── Figure2_Accuracy_Assessment.svg
├── Figure3_Performance_Benchmarks.svg
└── FigureS*_*.tiff             # 补充图表

/home/wangys/data/work/nbs/manuscript/tables/formatted/         # 论文数据表
├── Table1_Dataset_Summary.csv
├── Table1_Dataset_Summary.tex
└── *.csv                       # 其他数据表
```

## 下一步操作

1. 检查各分析阶段的输出结果
2. 将图表和数据表整合到论文中
3. 根据需要调整图表格式和内容
4. 准备投稿材料

## 注意事项

- 所有图表已按Bioinformatics期刊要求格式化 (300 DPI TIFF)
- 表格提供CSV和LaTeX两种格式
- 如有错误或需要调整，请重新运行相应分析阶段

---
**分析完成**: 2025-06-29 11:08:24
