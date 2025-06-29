# 基因ID重命名工具使用指南
# Gene ID Renaming Tool Usage Guide

## 概述 / Overview

这是一个用于EVM后处理的基因ID重命名工具，支持用户自定义前缀的基因ID命名规则。

This is a post-EVM processing tool for gene ID renaming that supports user-defined prefix naming conventions.

## 特性 / Features

- ✅ **用户指定前缀**: 支持任意用户定义的前缀，如 `NBS`, `CaCM334`, `AtCol0` 等
- ✅ **分层的ID格式**: 基因使用 `前缀_编号.g`，转录本使用 `前缀_编号.t` 格式
- ✅ **自动基因创建**: 为每个mRNA自动创建对应的基因特征
- ✅ **独立工具**: 不集成在EVM流程中，作为独立的后处理步骤
- ✅ **批量处理**: 支持单文件和多文件批量处理
- ✅ **完整性保持**: 保持GFF3文件的层次关系和完整性
- ✅ **可追溯性**: 生成ID映射文件和详细统计信息

## 安装和配置 / Installation and Configuration

工具已集成到NBS基因注释流水线中，无需额外安装。

The tool is integrated into the NBS gene annotation pipeline, no additional installation required.

## 使用方法 / Usage

### 1. 命令行使用 / Command Line Usage

#### 基本用法 / Basic Usage

```bash
# 使用默认前缀 "NBS"
python src/nbs_annotation/post_evm_renamer.py -i evm_output.gff3

# 使用自定义前缀 "CaCM334"
python src/nbs_annotation/post_evm_renamer.py -i evm_output.gff3 -p CaCM334

# 指定输出目录
python src/nbs_annotation/post_evm_renamer.py -i evm_output.gff3 -p CaCM334 -o /path/to/output

# 启用详细日志
python src/nbs_annotation/post_evm_renamer.py -i evm_output.gff3 -p CaCM334 -v
```

#### 批量处理 / Batch Processing

```bash
# 批量重命名多个文件
python src/nbs_annotation/post_evm_renamer.py -i file1.gff3 file2.gff3 file3.gff3 -p MyPrefix

# 使用通配符（需要shell展开）
python src/nbs_annotation/post_evm_renamer.py -i results/evm_*.gff3 -p CaCM334
```

### 2. Python API使用 / Python API Usage

```python
from src.nbs_annotation.post_evm_renamer import PostEVMRenamer

# 创建重命名器
renamer = PostEVMRenamer(prefix='CaCM334', output_dir='./renamed_output')

# 重命名单个文件
result = renamer.rename_evm_output('evm_output.gff3')

if result['success']:
    print(f"重命名成功: {result['statistics']['total_ids_renamed']} 个基因ID")
    print(f"输出文件: {result['output_file']}")
    print(f"映射文件: {result['mapping_file']}")
else:
    print(f"重命名失败: {result['error']}")

# 批量重命名
files = ['file1.gff3', 'file2.gff3', 'file3.gff3']
batch_result = renamer.batch_rename(files, prefix='MyPrefix')
print(f"批量处理完成: {batch_result['successful']}/{batch_result['total_files']} 成功")
```

### 3. 便捷方法 / Convenience Method

```python
from src.nbs_annotation.gene_id_renamer import GeneIDRenamer

# 一行代码完成重命名
id_mapping = GeneIDRenamer.rename_evm_output(
    'evm_output.gff3', 
    'renamed_output.gff3', 
    prefix='CaCM334'
)
print(f"重命名了 {len(id_mapping)} 个基因ID")
```

## 输出文件 / Output Files

每次重命名操作会生成以下文件：

Each renaming operation generates the following files:

1. **重命名后的GFF3文件** / **Renamed GFF3 file**: `{input_name}_{prefix}_renamed.gff3`
   - 包含重命名后的基因注释
   - Contains renamed gene annotations

2. **ID映射文件** / **ID mapping file**: `{input_name}_{prefix}_mapping.txt`
   - 原始ID到新ID的映射关系
   - Mapping from original IDs to new IDs

3. **统计信息文件** / **Statistics file**: `{input_name}_{prefix}_stats.json`
   - 详细的重命名统计信息
   - Detailed renaming statistics

4. **日志文件** / **Log file**: `gene_renaming.log`
   - 操作日志记录
   - Operation log records

## 示例 / Examples

### 示例1：辣椒CM334基因重命名 / Example 1: Pepper CM334 Gene Renaming

```bash
# 输入文件：evm_output.gff3
# 希望前缀：CaCM334
python src/nbs_annotation/post_evm_renamer.py -i evm_output.gff3 -p CaCM334

# 输出：
# - evm_output_CaCM334_renamed.gff3 (重命名后的注释文件)
# - evm_output_CaCM334_mapping.txt  (ID映射文件)
# - evm_output_CaCM334_stats.json   (统计信息)
```

**输入示例** (evm_output.gff3):
```gff3
Chr1    EVM    mRNA    1000    3000    .    +    .    ID=evm.model.Chr1.1;Name=evm.model.Chr1.1
Chr1    EVM    exon    1000    1500    .    +    .    ID=evm.model.Chr1.1.exon1;Parent=evm.model.Chr1.1
Chr5    EVM    mRNA    5000    7000    .    -    .    ID=evm.model.Chr5.1;Name=evm.model.Chr5.1
Chr5    EVM    exon    5000    5500    .    -    .    ID=evm.model.Chr5.1.exon1;Parent=evm.model.Chr5.1
```

**输出示例** (evm_output_CaCM334_renamed.gff3):
```gff3
##gff-version 3
Chr1    EVM    gene    1000    3000    .    +    .    ID=CaCM334_0001.g
Chr5    EVM    gene    5000    7000    .    -    .    ID=CaCM334_0002.g
Chr1    EVM    mRNA    1000    3000    .    +    .    ID=CaCM334_0001.t;Name=CaCM334_0001.t
Chr1    EVM    exon    1000    1500    .    +    .    ID=evm.model.Chr1.1.exon1;Parent=CaCM334_0001.t
Chr5    EVM    mRNA    5000    7000    .    -    .    ID=CaCM334_0002.t;Name=CaCM334_0002.t
Chr5    EVM    exon    5000    5500    .    -    .    ID=evm.model.Chr5.1.exon1;Parent=CaCM334_0002.t
```

**映射文件示例** (evm_output_CaCM334_mapping.txt):
```
# Gene ID Mapping File
# Original_ID	New_ID
evm.model.Chr1.1	CaCM334_0001.t
evm.model.Chr5.1	CaCM334_0002.t
```

### 示例2：水稻日本晴基因重命名 / Example 2: Rice Nipponbare Gene Renaming

```bash
python src/nbs_annotation/post_evm_renamer.py -i rice_evm_output.gff3 -p OsNB -o rice_renamed
```

输出ID格式：
- 基因: `OsNB_0001.g`, `OsNB_0002.g`, `OsNB_0003.g` ...
- 转录本: `OsNB_0001.t`, `OsNB_0002.t`, `OsNB_0003.t` ...

### 示例3：批量处理 / Example 3: Batch Processing

```bash
# 处理多个染色体的EVM输出文件
python src/nbs_annotation/post_evm_renamer.py \
    -i chr1_evm.gff3 chr2_evm.gff3 chr3_evm.gff3 \
    -p CaCM334 \
    -o batch_output
```

## 配置选项 / Configuration Options

### 预设前缀 / Preset Prefixes

配置文件 `config/default.yaml` 中包含常用的预设前缀：

The configuration file `config/default.yaml` contains common preset prefixes:

```yaml
post_processing:
  gene_renaming:
    preset_prefixes:
      rice_nipponbare: "OsNB"      # 水稻日本晴
      rice_93_11: "Os93"           # 水稻93-11
      pepper_cm334: "CaCM334"      # 辣椒CM334
      pepper_zunla: "CaZunla"      # 辣椒Zunla-1
      arabidopsis_col0: "AtCol0"   # 拟南芥Col-0
      tomato_heinz: "SlHeinz"      # 番茄Heinz
```

### 前缀验证规则 / Prefix Validation Rules

- 只允许字母、数字和下划线：`[A-Za-z0-9_]+`
- 不能为空
- 示例有效前缀：`NBS`, `CaCM334`, `Rice_NB`, `Test_123`
- 示例无效前缀：`Ca-CM334`, `Test.Prefix`, ``, `Prefix with spaces`

## 与EVM流程集成 / Integration with EVM Pipeline

该工具设计为EVM流程的**独立后处理步骤**，建议的工作流程：

This tool is designed as an **independent post-processing step** for the EVM pipeline. Recommended workflow:

1. **运行EVM集成** / **Run EVM Integration**:
   ```bash
   python src/nbs_annotation/evm_runner.py
   ```

2. **验证EVM输出** / **Validate EVM Output**:
   检查 `results/evm_integration/evm_output.gff3`

3. **运行基因ID重命名** / **Run Gene ID Renaming**:
   ```bash
   python src/nbs_annotation/post_evm_renamer.py \
       -i results/evm_integration/evm_output.gff3 \
       -p CaCM334 \
       -o results/renamed_genes
   ```

4. **验证重命名结果** / **Validate Renaming Results**:
   检查输出文件和统计信息

## 故障排除 / Troubleshooting

### 常见问题 / Common Issues

1. **输入文件不存在**
   ```
   Error: Input file does not exist: /path/to/file.gff3
   ```
   **解决方案**: 检查文件路径是否正确

2. **无效的前缀**
   ```
   Error: Invalid prefix: Ca-CM334. Only alphanumeric characters and underscores are allowed.
   ```
   **解决方案**: 使用只包含字母、数字和下划线的前缀

3. **没有找到EVM ID**
   ```
   Warning: No EVM IDs found in input file.
   ```
   **解决方案**: 确认输入文件是EVM的输出文件

4. **权限错误**
   ```
   Error: Permission denied: /path/to/output
   ```
   **解决方案**: 检查输出目录的写权限

### 调试模式 / Debug Mode

启用详细日志以获得更多调试信息：

Enable verbose logging for more debugging information:

```bash
python src/nbs_annotation/post_evm_renamer.py -i evm_output.gff3 -p CaCM334 -v
```

## 最佳实践 / Best Practices

1. **备份原始文件**: 重命名前始终备份原始EVM输出文件
2. **验证输出**: 重命名后检查输出文件的完整性
3. **保存映射文件**: ID映射文件对于后续分析很重要
4. **使用有意义的前缀**: 选择能反映物种、品种或项目的前缀
5. **批量处理**: 对于多个文件使用批量处理功能
6. **检查日志**: 查看日志文件了解处理详情

## 性能指标 / Performance Metrics

- **处理速度**: 约1000个基因ID/秒
- **内存使用**: 典型GFF3文件 < 100MB
- **支持文件大小**: 最大50GB（配置可调）
- **并发处理**: 支持多文件并行处理

## 更新日志 / Change Log

### v1.0.0 (2025-06-26)
- ✅ 初始版本发布
- ✅ 支持用户定义前缀
- ✅ 独立的后处理工具
- ✅ 命令行和Python API
- ✅ 批量处理功能
- ✅ 完整的文档和测试

## 联系支持 / Support

如有问题或建议，请联系开发团队或提交Issue。

For questions or suggestions, please contact the development team or submit an issue.