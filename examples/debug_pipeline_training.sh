#!/bin/bash
"""
Debug version of pipeline training example

调试版本的流水线训练示例脚本
"""

set -e  # Exit on any error

echo "🧪 调试Augustus训练集成"
echo "========================"

# 设置正确的工作目录和Python路径
cd /data/usdata/wangys/work/nbs
export PYTHONPATH="/data/usdata/wangys/work/nbs/src"

# 检查文件
echo "🔍 检查文件..."
echo "当前目录: $(pwd)"
echo "Python路径: $PYTHONPATH"

# 检查输入文件是否存在
if [[ ! -f "genome/osa.fa" ]]; then
    echo "❌ 基因组文件不存在: genome/osa.fa"
    echo "可用文件:"
    ls -la genome/ || echo "genome目录不存在"
    exit 1
fi

if [[ ! -f "db/AllResistanceGenes.fasta" ]]; then
    echo "❌ 蛋白质文件不存在: db/AllResistanceGenes.fasta" 
    echo "可用文件:"
    ls -la db/ || echo "db目录不存在"
    exit 1
fi

echo "✅ 输入文件检查通过"
echo ""

# 测试1: 检查命令行帮助
echo "🧪 测试1: 检查命令行帮助..."
python -m nbs_annotation.main --help | grep -A 5 "Augustus Training" || echo "未找到Augustus Training选项"
echo ""

# 测试2: Dry run测试
echo "🧪 测试2: Dry run测试..."
python -m nbs_annotation.main \
    --genome genome/osa.fa \
    --proteins db/AllResistanceGenes.fasta \
    --enable-training \
    --training-species-name debug_test_model \
    --training-quality medium \
    --training-min-genes 5 \
    --skip-training-optimization \
    --training-timeout 30 \
    --dry-run

echo ""

# 测试3: 验证输入文件
echo "🧪 测试3: 验证输入文件..."
python -m nbs_annotation.main \
    --genome genome/osa.fa \
    --proteins db/AllResistanceGenes.fasta \
    --validate-only

echo ""

# 测试4: 检查工具可用性
echo "🧪 测试4: 检查工具可用性..."
python -m nbs_annotation.main \
    --check-tools

echo ""

# 测试5: 运行单个阶段 - NLR定位
echo "🧪 测试5: 运行NLR定位阶段..."
OUTPUT_DIR="results/debug_test"
mkdir -p "$OUTPUT_DIR"

python -m nbs_annotation.main \
    --genome genome/osa.fa \
    --stage nlr_localization \
    --output "$OUTPUT_DIR" \
    --verbose

echo "NLR定位结果:"
ls -la "$OUTPUT_DIR/nlr_localization/" || echo "NLR定位输出目录不存在"
echo ""

# 测试6: 运行蛋白质比对阶段
echo "🧪 测试6: 运行蛋白质比对阶段..."
python -m nbs_annotation.main \
    --genome genome/osa.fa \
    --proteins db/AllResistanceGenes.fasta \
    --stage protein_alignment \
    --output "$OUTPUT_DIR" \
    --verbose

echo "蛋白质比对结果:"
ls -la "$OUTPUT_DIR/protein_alignment/" || echo "蛋白质比对输出目录不存在"
echo ""

# 测试7: 检查是否有高质量Miniprot结果
echo "🧪 测试7: 检查Miniprot质量过滤结果..."
if [[ -d "$OUTPUT_DIR/protein_alignment/filtered" ]]; then
    echo "发现过滤结果目录:"
    ls -la "$OUTPUT_DIR/protein_alignment/filtered/"
    
    # 检查高质量结果
    if [[ -f "$OUTPUT_DIR/protein_alignment/filtered/miniprot_high_quality.gff3" ]]; then
        echo "✅ 发现高质量Miniprot结果"
        echo "高质量结果统计:"
        grep -c "^[^#]" "$OUTPUT_DIR/protein_alignment/filtered/miniprot_high_quality.gff3" || echo "0 条记录"
    else
        echo "❌ 未发现高质量Miniprot结果"
    fi
else
    echo "❌ 未发现过滤结果目录"
fi
echo ""

# 测试8: 尝试运行Augustus训练（如果有训练数据）
if [[ -f "$OUTPUT_DIR/protein_alignment/filtered/miniprot_high_quality.gff3" ]]; then
    echo "🧪 测试8: 运行Augustus训练阶段..."
    python -m nbs_annotation.main \
        --genome genome/osa.fa \
        --proteins db/AllResistanceGenes.fasta \
        --stage augustus_training \
        --output "$OUTPUT_DIR" \
        --enable-training \
        --training-species-name debug_trained_model \
        --training-quality high \
        --training-min-genes 2 \
        --skip-training-optimization \
        --training-timeout 60 \
        --verbose
    
    echo "Augustus训练结果:"
    ls -la "$OUTPUT_DIR/augustus_training/" || echo "Augustus训练输出目录不存在"
else
    echo "⏭️  跳过Augustus训练测试 - 缺少训练数据"
fi

echo ""
echo "🎉 调试测试完成!"
echo "查看结果目录: $OUTPUT_DIR"