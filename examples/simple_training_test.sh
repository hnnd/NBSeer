#!/bin/bash
# Simple Augustus Training Test
# 简单的Augustus训练测试

set -e

echo "🧪 简单Augustus训练测试"
echo "====================="

# 设置环境
cd /data/usdata/wangys/work/nbs
export PYTHONPATH="/data/usdata/wangys/work/nbs/src"

echo "当前目录: $(pwd)"
echo "Python路径: $PYTHONPATH"
echo ""

# 测试1: 检查basic功能
echo "🔍 测试1: 检查基本功能..."
python -m nbs_annotation.main --help | head -5
echo ""

# 测试2: Dry run with training
echo "🔍 测试2: Dry run测试带训练选项..."
python -m nbs_annotation.main \
    --genome genome/osa.fa \
    --proteins db/AllResistanceGenes.fasta \
    --enable-training \
    --training-species-name simple_test_model \
    --training-quality high \
    --training-min-genes 3 \
    --skip-training-optimization \
    --training-timeout 30 \
    --dry-run
echo ""

# 测试3: 检查工具可用性
echo "🔍 测试3: 检查工具可用性..."
python -m nbs_annotation.main --check-tools
echo ""

# 测试4: 验证配置文件中的训练设置
echo "🔍 测试4: 检查配置文件..."
echo "Augustus训练配置:"
grep -A 10 "miniprot_training:" config/default.yaml
echo ""

# 测试5: 尝试直接使用训练CLI
echo "🔍 测试5: 测试独立训练CLI..."
python src/nbs_annotation/train_augustus_cli.py --help | head -10
echo ""

# 测试6: 检查是否存在现有的Miniprot结果
echo "🔍 测试6: 检查现有的测试数据..."
if [[ -f "data/test/test_output/protein_alignment/filtered/miniprot_high_quality.gff3" ]]; then
    echo "✅ 发现现有的高质量Miniprot结果"
    echo "记录数量:"
    grep -c "^[^#]" data/test/test_output/protein_alignment/filtered/miniprot_high_quality.gff3 || echo "0条记录"
    
    echo ""
    echo "🧪 测试Augustus训练（使用现有数据）..."
    python src/nbs_annotation/train_augustus_cli.py \
        --species simple_test_model_cli \
        --genome genome/osa.fa \
        --miniprot-file data/test/test_output/protein_alignment/filtered/miniprot_high_quality.gff3 \
        --quality high \
        --min-genes 2 \
        --no-optimization \
        --timeout 60 \
        --working-dir results/simple_training_test \
        --verbose
else
    echo "❌ 未发现现有的Miniprot结果"
    echo "可用的测试文件:"
    find data/test/ -name "*.gff*" -o -name "*.fa" | head -10
fi

echo ""
echo "🎉 简单测试完成!"