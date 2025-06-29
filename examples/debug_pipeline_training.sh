#!/bin/bash
"""
Debug version of pipeline training example

Debugging version of the pipeline training example script
"""

set -e  # Exit on any error

echo "🧪 Debug Augustus Training Integration"
echo "========================================"

# Set up correct working directory and Python path
cd /data/usdata/wangys/work/nbs
export PYTHONPATH="/data/usdata/wangys/work/nbs/src"

# Check files
echo "🔍 Checking files..."
echo "Current directory: $(pwd)"
echo "Python path: $PYTHONPATH"

# Check if input files exist
if [[ ! -f "genome/osa.fa" ]]; then
    echo "❌ Genome file does not exist: genome/osa.fa"
    echo "Available files:"
    ls -la genome/ || echo "genome directory does not exist"
    exit 1
fi

if [[ ! -f "db/AllResistanceGenes.fasta" ]]; then
    echo "❌ Protein file does not exist: db/AllResistanceGenes.fasta" 
    echo "Available files:"
    ls -la db/ || echo "db directory does not exist"
    exit 1
fi

echo "✅ Input file validation passed"
echo ""

# Test 1: Check command line help
echo "🧪 Test 1: Checking command line help..."
python -m nbseer.main --help | grep -A 5 "Augustus Training" || echo "Augustus Training options not found"
echo ""

# Test 2: Dry run test
echo "🧪 Test 2: Dry run test..."
python -m nbseer.main \
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

# Test 3: Validate input files
echo "🧪 Test 3: Validating input files..."
python -m nbseer.main \
    --genome genome/osa.fa \
    --proteins db/AllResistanceGenes.fasta \
    --validate-only

echo ""

# Test 4: Check tool availability
echo "🧪 Test 4: Checking tool availability..."
python -m nbseer.main \
    --check-tools

echo ""

# Test 5: Run single stage - NLR localization
echo "🧪 Test 5: Running NLR localization stage..."
OUTPUT_DIR="results/debug_test"
mkdir -p "$OUTPUT_DIR"

python -m nbseer.main \
    --genome genome/osa.fa \
    --stage nlr_localization \
    --output "$OUTPUT_DIR" \
    --verbose

echo "NLR localization results:"
ls -la "$OUTPUT_DIR/nlr_localization/" || echo "NLR localization output directory does not exist"
echo ""

# Test 6: Run protein alignment stage
echo "🧪 Test 6: Running protein alignment stage..."
python -m nbseer.main \
    --genome genome/osa.fa \
    --proteins db/AllResistanceGenes.fasta \
    --stage protein_alignment \
    --output "$OUTPUT_DIR" \
    --verbose

echo "Protein alignment results:"
ls -la "$OUTPUT_DIR/protein_alignment/" || echo "Protein alignment output directory does not exist"
echo ""

# Test 7: Check if there are high-quality Miniprot results
echo "🧪 Test 7: Checking Miniprot quality filtering results..."
if [[ -d "$OUTPUT_DIR/protein_alignment/filtered" ]]; then
    echo "Found filtering results directory:"
    ls -la "$OUTPUT_DIR/protein_alignment/filtered/"
    
    # Check high-quality results
    if [[ -f "$OUTPUT_DIR/protein_alignment/filtered/miniprot_high_quality.gff3" ]]; then
        echo "✅ Found high-quality Miniprot results"
        echo "High-quality results statistics:"
        grep -c "^[^#]" "$OUTPUT_DIR/protein_alignment/filtered/miniprot_high_quality.gff3" || echo "0 records"
    else
        echo "❌ No high-quality Miniprot results found"
    fi
else
    echo "❌ No filtering results directory found"
fi
echo ""

# Test 8: Try running Augustus training (if training data is available)
if [[ -f "$OUTPUT_DIR/protein_alignment/filtered/miniprot_high_quality.gff3" ]]; then
    echo "🧪 Test 8: Running Augustus training stage..."
    python -m nbseer.main \
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
    
    echo "Augustus training results:"
    ls -la "$OUTPUT_DIR/augustus_training/" || echo "Augustus training output directory does not exist"
else
    echo "⏭️  Skipping Augustus training test - missing training data"
fi

echo ""
echo "🎉 Debug test completed!"
echo "View results directory: $OUTPUT_DIR"