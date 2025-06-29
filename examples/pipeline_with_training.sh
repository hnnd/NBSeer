#!/bin/bash
"""
Example: NBS Annotation Pipeline with Augustus Training
使用Augustus训练的NBS注释流水线示例

This script demonstrates how to run the complete NBS annotation pipeline
with integrated Augustus model training using high-quality Miniprot results.

Author: NBS Annotation Pipeline
Date: 2025-06-25
"""

set -e  # Exit on any error

echo "NBS注释流水线 - Augustus训练集成示例"
echo "=================================================="

# 配置参数
GENOME_FILE="genome/osa.fa"
PROTEIN_FILE="db/AllResistanceGenes.fasta"
OUTPUT_DIR="results/example_with_training"
SPECIES_NAME="nbs_rice_example_$(date +%Y%m%d)"

echo "📋 配置信息:"
echo "  基因组文件: $GENOME_FILE"
echo "  蛋白质文件: $PROTEIN_FILE"
echo "  输出目录: $OUTPUT_DIR"
echo "  模型名称: $SPECIES_NAME"
echo ""

# 检查输入文件
echo "🔍 检查输入文件..."
if [[ ! -f "$GENOME_FILE" ]]; then
    echo "❌ 基因组文件不存在: $GENOME_FILE"
    exit 1
fi

if [[ ! -f "$PROTEIN_FILE" ]]; then
    echo "❌ 蛋白质文件不存在: $PROTEIN_FILE" 
    exit 1
fi

echo "✅ 输入文件检查通过"
echo ""

# 示例1: 完整流水线带训练（推荐）
echo "🚀 示例1: 运行完整流水线带Augustus训练"
echo "----------------------------------------"

python -m nbs_annotation.main \
    --genome "$GENOME_FILE" \
    --proteins "$PROTEIN_FILE" \
    --output "$OUTPUT_DIR/full_pipeline" \
    --enable-training \
    --training-species-name "$SPECIES_NAME" \
    --training-quality high \
    --training-min-genes 10 \
    --training-optimization-rounds 1 \
    --training-timeout 120 \
    --threads 4 \
    --verbose

echo "✅ 完整流水线完成!"
echo ""

# 示例2: 分阶段运行
echo "🔄 示例2: 分阶段运行流水线"
echo "------------------------"

STAGE_OUTPUT="$OUTPUT_DIR/staged_pipeline"

echo "  阶段1: NLR基因定位"
python -m nbs_annotation.main \
    --genome "$GENOME_FILE" \
    --stage nlr_localization \
    --output "$STAGE_OUTPUT" \
    --verbose

echo "  阶段2: 蛋白质比对"
python -m nbs_annotation.main \
    --genome "$GENOME_FILE" \
    --proteins "$PROTEIN_FILE" \
    --stage protein_alignment \
    --output "$STAGE_OUTPUT" \
    --verbose

echo "  阶段3: Augustus训练"
python -m nbs_annotation.main \
    --genome "$GENOME_FILE" \
    --proteins "$PROTEIN_FILE" \
    --stage augustus_training \
    --output "$STAGE_OUTPUT" \
    --enable-training \
    --training-species-name "${SPECIES_NAME}_staged" \
    --training-quality medium \
    --skip-training-optimization \
    --verbose

echo "  阶段4: 基因预测（使用训练的模型）"
python -m nbs_annotation.main \
    --genome "$GENOME_FILE" \
    --stage gene_prediction \
    --output "$STAGE_OUTPUT" \
    --augustus-model "${SPECIES_NAME}_staged" \
    --nlr-candidates "$STAGE_OUTPUT/nlr_localization/NLR_candidates.gff" \
    --verbose

echo "  阶段5: 证据整合"
python -m nbs_annotation.main \
    --genome "$GENOME_FILE" \
    --stage evidence_integration \
    --output "$STAGE_OUTPUT" \
    --verbose

echo "✅ 分阶段流水线完成!"
echo ""

# 示例3: 快速测试模式
echo "⚡ 示例3: 快速测试模式"
echo "--------------------"

python -m nbs_annotation.main \
    --genome "$GENOME_FILE" \
    --proteins "$PROTEIN_FILE" \
    --output "$OUTPUT_DIR/quick_test" \
    --enable-training \
    --training-species-name "${SPECIES_NAME}_quick" \
    --training-quality medium \
    --training-min-genes 5 \
    --skip-training-optimization \
    --training-timeout 30 \
    --threads 2 \
    --verbose

echo "✅ 快速测试完成!"
echo ""

# 示例4: 批量训练不同质量级别
echo "📊 示例4: 批量训练不同质量模型"
echo "-----------------------------"

for quality in high medium low; do
    echo "  训练质量级别: $quality"
    python -m nbs_annotation.main \
        --genome "$GENOME_FILE" \
        --proteins "$PROTEIN_FILE" \
        --stage augustus_training \
        --output "$OUTPUT_DIR/batch_training" \
        --enable-training \
        --training-species-name "nbs_${quality}_$(date +%H%M)" \
        --training-quality "$quality" \
        --training-min-genes 3 \
        --skip-training-optimization \
        --training-timeout 20 \
        --verbose
done

echo "✅ 批量训练完成!"
echo ""

# 检查结果
echo "📋 检查结果文件"
echo "---------------"

echo "输出目录结构:"
find "$OUTPUT_DIR" -name "*.gff*" -o -name "*.json" -o -name "*.md" | head -20

echo ""
echo "训练的Augustus模型:"
if [[ -n "$AUGUSTUS_CONFIG_PATH" ]]; then
    find "$AUGUSTUS_CONFIG_PATH/species" -name "*${SPECIES_NAME}*" -type d 2>/dev/null || echo "  (需要检查AUGUSTUS_CONFIG_PATH)"
else
    echo "  (需要设置AUGUSTUS_CONFIG_PATH环境变量)"
fi

echo ""
echo "🎉 所有示例运行完成!"
echo "📖 查看详细结果:"
echo "  - 完整流水线: $OUTPUT_DIR/full_pipeline/pipeline_results.json"
echo "  - 分阶段结果: $STAGE_OUTPUT/pipeline_results.json"
echo "  - 训练报告: $OUTPUT_DIR/*/augustus_training/*_report.md"
echo ""
echo "🔧 使用训练的模型:"
echo "  augustus --species=$SPECIES_NAME --gff3=on genome.fa"
echo ""