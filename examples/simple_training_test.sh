#!/bin/bash
# Simple Augustus Training Test
# Basic Augustus training test script

set -e

echo "🧪 Simple Augustus Training Test"
echo "================================"

# Set up environment
cd /data/usdata/wangys/work/nbs
export PYTHONPATH="/data/usdata/wangys/work/nbs/src"

echo "Current directory: $(pwd)"
echo "Python path: $PYTHONPATH"
echo ""

# Test 1: Check basic functionality
echo "🔍 Test 1: Checking basic functionality..."
python -m nbs_annotation.main --help | head -5
echo ""

# Test 2: Dry run with training
echo "🔍 Test 2: Dry run test with training options..."
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

# Test 3: Check tool availability
echo "🔍 Test 3: Checking tool availability..."
python -m nbs_annotation.main --check-tools
echo ""

# Test 4: Validate training settings in config file
echo "🔍 Test 4: Checking configuration file..."
echo "Augustus training configuration:"
grep -A 10 "miniprot_training:" config/default.yaml
echo ""

# Test 5: Try using training CLI directly
echo "🔍 Test 5: Testing standalone training CLI..."
python src/nbs_annotation/train_augustus_cli.py --help | head -10
echo ""

# Test 6: Check if existing Miniprot results exist
echo "🔍 Test 6: Checking existing test data..."
if [[ -f "data/test/test_output/protein_alignment/filtered/miniprot_high_quality.gff3" ]]; then
    echo "✅ Found existing high-quality Miniprot results"
    echo "Number of records:"
    grep -c "^[^#]" data/test/test_output/protein_alignment/filtered/miniprot_high_quality.gff3 || echo "0 records"
    
    echo ""
    echo "🧪 Testing Augustus training (using existing data)..."
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
    echo "❌ No existing Miniprot results found"
    echo "Available test files:"
    find data/test/ -name "*.gff*" -o -name "*.fa" | head -10
fi

echo ""
echo "🎉 Simple test completed!"