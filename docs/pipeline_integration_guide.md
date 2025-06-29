# NBS注释流水线中的Augustus训练集成

本文档描述了如何在NBS基因注释流水线中使用集成的Augustus训练功能。

## 概述

Augustus训练功能已完全集成到NBS注释流水线中，提供以下特性：

- **自动集成**: 训练在蛋白质比对和基因预测之间自动执行
- **命令行控制**: 通过命令行参数启用和配置训练
- **配置文件支持**: 支持通过配置文件预设训练参数
- **检查点支持**: 训练结果被缓存，避免重复训练
- **质量控制**: 自动验证训练数据质量和模型有效性

## 快速开始

### 基本用法

```bash
# 启用训练的完整流水线
python -m nbs_annotation.main \
    --genome genome.fa \
    --proteins proteins.fa \
    --enable-training \
    --training-species-name nbs_rice_v2

# 仅运行Augustus训练阶段
python -m nbs_annotation.main \
    --genome genome.fa \
    --proteins proteins.fa \
    --stage augustus_training \
    --enable-training \
    --training-species-name nbs_rice_v2
```

### 高级配置

```bash
# 自定义训练参数
python -m nbs_annotation.main \
    --genome genome.fa \
    --proteins proteins.fa \
    --enable-training \
    --training-species-name custom_nbs_model \
    --training-quality medium \
    --training-min-genes 30 \
    --training-optimization-rounds 2 \
    --training-timeout 360 \
    --output results/
```

## 流水线集成

### 流水线阶段

Augustus训练作为新阶段插入到现有流水线中：

1. **NLR基因定位** (`nlr_localization`)
2. **蛋白质比对** (`protein_alignment`) 
3. **Augustus训练** (`augustus_training`) - 新增
4. **基因结构预测** (`gene_prediction`)
5. **证据整合** (`evidence_integration`)

### 自动工作流

当启用训练时，流水线会：

1. 运行Miniprot蛋白质比对
2. 提取高质量比对结果用于训练
3. 执行Augustus模型训练
4. 验证训练的模型
5. 在基因预测中使用新训练的模型

### 检查点机制

- **训练缓存**: 训练完成后创建 `training_completed.json` 标记文件
- **模型重用**: 后续运行会自动检测并使用已训练的模型  
- **增量运行**: 支持从任意阶段恢复流水线执行

## 命令行选项

### 训练控制选项

| 选项 | 类型 | 默认值 | 描述 |
|------|------|--------|------|
| `--enable-training` | flag | false | 启用Augustus模型训练 |
| `--training-species-name` | string | auto | 新训练模型的物种名称 |
| `--training-quality` | choice | high | 训练数据质量级别 (high/medium/low/all) |
| `--training-min-genes` | int | 20 | 所需的最少训练基因数 |
| `--training-optimization-rounds` | int | 1 | 优化轮数 |
| `--training-timeout` | int | 240 | 训练超时时间（分钟） |
| `--skip-training-optimization` | flag | false | 跳过参数优化 |

### 阶段选择

```bash
# 运行所有阶段（包括训练）
--stage all

# 仅运行Augustus训练
--stage augustus_training

# 运行特定阶段序列
--stage nlr_localization
--stage protein_alignment  
--stage augustus_training
--stage gene_prediction
--stage evidence_integration
```

## 配置文件集成

### YAML配置

在 `config/default.yaml` 中预设训练参数：

```yaml
tools:
  augustus:
    training:
      miniprot_training:
        enabled: false  # 通过命令行启用
        species_name: "nbs_rice_miniprot"
        quality_filter: "high"
        min_training_genes: 20
        min_identity_threshold: 0.95
        max_frameshifts: 0
        max_stop_codons: 0
        flanking_dna_length: 4000
        optimization_rounds: 1
        cpus: 4
        timeout_minutes: 240
        backup_existing_model: true
        create_training_report: true
```

### 配置覆盖

命令行参数会自动覆盖配置文件设置：

```bash
# 命令行参数覆盖配置文件
python -m nbs_annotation.main \
    --enable-training \
    --training-quality medium \
    --training-min-genes 30
```

等价于在配置中设置：
```yaml
tools.augustus.training.miniprot_training.enabled: true
tools.augustus.training.miniprot_training.quality_filter: "medium"  
tools.augustus.training.miniprot_training.min_training_genes: 30
```

## 使用场景

### 场景1: 标准流水线带训练

```bash
# 一次性运行完整流水线并训练自定义模型
python -m nbs_annotation.main \
    --genome rice_genome.fa \
    --proteins resistance_proteins.fa \
    --enable-training \
    --training-species-name rice_nb_resistance \
    --output results/rice_annotation/
```

**流程**:
1. NLR基因定位
2. 蛋白质比对 → 生成高质量训练数据
3. Augustus训练 → 创建 `rice_nb_resistance` 模型
4. 基因预测 → 使用新训练的模型
5. 证据整合 → 最终注释结果

### 场景2: 分阶段训练和预测

```bash
# 第一步：运行到训练完成
python -m nbs_annotation.main \
    --genome genome.fa \
    --proteins proteins.fa \
    --stage augustus_training \
    --enable-training \
    --training-species-name custom_model

# 第二步：使用训练的模型进行预测
python -m nbs_annotation.main \
    --genome genome.fa \
    --proteins proteins.fa \
    --augustus-model custom_model \
    --resume-from gene_prediction
```

### 场景3: 批量模型训练

```bash
# 训练多个质量级别的模型
for quality in high medium low; do
    python -m nbs_annotation.main \
        --genome genome.fa \
        --proteins proteins.fa \
        --stage augustus_training \
        --enable-training \
        --training-species-name "nbs_${quality}_model" \
        --training-quality "${quality}"
done
```

### 场景4: 测试和验证

```bash
# 快速训练测试（跳过优化）
python -m nbs_annotation.main \
    --genome test_genome.fa \
    --proteins test_proteins.fa \
    --enable-training \
    --training-species-name test_model \
    --training-min-genes 5 \
    --skip-training-optimization \
    --training-timeout 30 \
    --dry-run
```

## 输出文件

### 训练相关输出

```
results/
├── augustus_training/
│   ├── training_completed.json          # 训练完成标记
│   ├── nbs_species_training_data.gff3   # 处理后的训练数据
│   ├── nbs_species_training.log         # 训练日志
│   └── nbs_species_miniprot_training_report.md  # 详细报告
├── gene_prediction/
│   ├── predictions.gff                  # 使用训练模型的预测结果
│   └── predictions_genome_coords.gff    # 基因组坐标版本
└── pipeline_results.json               # 完整流水线结果
```

### Augustus模型文件

```
$AUGUSTUS_CONFIG_PATH/species/nbs_species_name/
├── nbs_species_name_parameters.cfg      # 模型参数
├── nbs_species_name_exon_probs.pbl      # 外显子概率
├── nbs_species_name_intron_probs.pbl    # 内含子概率
├── nbs_species_name_igenic_probs.pbl    # 基因间区概率
└── ... (其他模型文件)
```

## 质量控制

### 训练数据验证

流水线自动验证：
- **最少基因数**: 检查训练基因数量是否满足要求
- **序列质量**: 验证基因组和GFF文件格式
- **比对质量**: 确保Miniprot结果质量符合标准

### 模型验证

训练完成后自动检查：
- **模型文件完整性**: 验证所有必需文件存在
- **预测功能测试**: 运行简单序列预测测试
- **参数合理性**: 检查模型参数范围

### 错误处理

- **训练失败**: 自动回退到默认Augustus模型
- **数据不足**: 提供详细错误信息和建议
- **超时处理**: 支持自定义超时时间

## 性能优化

### 资源配置

```bash
# 高性能配置
python -m nbs_annotation.main \
    --enable-training \
    --threads 16 \
    --memory 64 \
    --training-optimization-rounds 3 \
    --training-timeout 480
```

### 快速模式

```bash
# 快速训练模式
python -m nbs_annotation.main \
    --enable-training \
    --training-quality medium \
    --training-min-genes 10 \
    --skip-training-optimization \
    --training-timeout 60
```

## 故障排除

### 常见问题

#### 1. 训练数据不足
```
错误: Insufficient training genes: 5 < 20
解决: 降低 --training-min-genes 或使用更低质量数据 --training-quality medium
```

#### 2. Augustus配置问题
```
错误: autoAugTrain.pl script not found
解决: 检查Augustus安装和AUGUSTUS_CONFIG_PATH环境变量
```

#### 3. 训练超时
```
错误: Training timed out after 240 minutes
解决: 增加 --training-timeout 或减少 --training-optimization-rounds
```

#### 4. 内存不足
```
错误: Training failed with memory error
解决: 增加 --memory 或减少 --training-min-genes
```

### 调试模式

```bash
# 启用详细日志
python -m nbs_annotation.main \
    --enable-training \
    -vv \
    --log-file training_debug.log
```

### 手动训练

如果流水线训练失败，可以手动运行：

```bash
# 使用独立训练脚本
python src/nbs_annotation/train_augustus_cli.py \
    --species manual_model \
    --quality high \
    --genome genome.fa \
    --verbose
```

## 高级主题

### 自定义训练数据

```python
# 使用Python API自定义训练流程
from nbs_annotation.augustus_miniprot_trainer import AugustusMiniprotTrainer, MiniprotTrainingConfig

config = MiniprotTrainingConfig(
    species_name="custom_nbs",
    miniprot_gff_file="custom_miniprot_results.gff3",
    quality_filter="high",
    optimization_rounds=2
)

trainer = AugustusMiniprotTrainer(config)
result = trainer.train_augustus_model()
```

### 与其他工具集成

```bash
# 训练后立即用于其他预测工具
python -m nbs_annotation.main --enable-training --training-species-name new_model
augustus --species=new_model --gff3=on other_genome.fa > predictions.gff3
```

### 模型版本管理

```bash
# 为不同数据集训练不同版本
python -m nbs_annotation.main \
    --enable-training \
    --training-species-name "nbs_v$(date +%Y%m%d)" \
    --training-quality high
```

## 最佳实践

### 1. 训练数据质量
- 优先使用 `--training-quality high`
- 确保至少50-100个高质量训练基因
- 验证基因组质量和完整性

### 2. 资源管理
- 根据数据大小调整 `--training-timeout`
- 使用适当的CPU数量 (`--threads`)
- 监控内存使用情况

### 3. 迭代改进
- 从快速训练开始验证流程
- 逐步增加优化轮数
- 比较不同质量级别的结果

### 4. 版本控制
- 使用有意义的模型名称
- 保存训练配置和日志
- 记录模型性能指标

## 参考链接

- [Augustus训练详细指南](augustus_training_guide.md)
- [Miniprot结果处理](miniprot_processing_guide.md)
- [EVM证据整合](evm_integration_guide.md)
- [流水线配置参考](pipeline_configuration.md)