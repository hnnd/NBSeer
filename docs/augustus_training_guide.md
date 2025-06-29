# Augustus Training with High-Quality Miniprot Results

This guide explains how to use the enhanced Augustus training functionality that leverages high-quality Miniprot protein alignment results to train custom Augustus species models for NBS gene prediction.

## Overview

The Augustus training system integrates high-quality Miniprot protein alignment results with the `autoAugTrain.pl` script to create optimized gene prediction models specifically for NBS resistance genes. This approach provides several advantages:

- **High-quality training data**: Uses only high-confidence protein alignments
- **Automated workflow**: Streamlined process from Miniprot results to trained models
- **Quality filtering**: Removes frameshifts, stop codons, and low-identity alignments
- **Comprehensive validation**: Tests trained models and generates detailed reports
- **Flexible configuration**: Supports various quality levels and training parameters

## Quick Start

### Basic Usage

Train a new Augustus model using high-quality Miniprot results:

```bash
# Basic training with default parameters
python src/nbs_annotation/train_augustus_cli.py --species nbs_rice_v2 --quality high

# Training with custom parameters
python src/nbs_annotation/train_augustus_cli.py \
    --species nbs_rice_v3 \
    --quality medium \
    --cpus 8 \
    --timeout 360 \
    --verbose
```

### Using the Python API

```python
from nbs_annotation.augustus_miniprot_trainer import (
    AugustusMiniprotTrainer, 
    MiniprotTrainingConfig
)

# Configure training
config = MiniprotTrainingConfig(
    species_name="nbs_rice_custom",
    quality_filter="high",
    optimization_rounds=1,
    cpus=4
)

# Run training
trainer = AugustusMiniprotTrainer(config)
result = trainer.train_augustus_model()

if result.training_success:
    print(f"Training successful! Model: {result.species_name}")
    print(f"Model directory: {result.model_directory}")
else:
    print(f"Training failed: {result.error_message}")
```

## Prerequisites

### Software Requirements

1. **Augustus**: Properly installed with scripts
   ```bash
   # Check Augustus installation
   augustus --help
   which autoAugTrain.pl
   ```

2. **Environment Variables**:
   ```bash
   export AUGUSTUS_CONFIG_PATH="/path/to/augustus/config"
   export PATH="$PATH:/path/to/augustus/scripts"
   ```

3. **Required Files**:
   - Genome FASTA file (e.g., `genome/osa.fa`)
   - High-quality Miniprot results (auto-detected in standard locations)

### Data Requirements

The training system expects Miniprot results in GFF3 format with quality filtering already applied. The pipeline automatically looks for these files:

- `data/test/test_output/protein_alignment/filtered/miniprot_high_quality.gff3`
- `data/output/protein_alignment/filtered/miniprot_high_quality.gff3`
- `results/protein_alignment/filtered/miniprot_high_quality.gff3`

## Configuration Options

### Quality Levels

| Quality Level | Identity Threshold | Max Frameshifts | Max Stop Codons | Description |
|---------------|-------------------|-----------------|------------------|-------------|
| `high`        | ≥95%              | 0               | 0                | Perfect alignments only |
| `medium`      | ≥85%              | 0               | 0                | Good quality alignments |
| `low`         | ≥70%              | 1               | 1                | Lower quality alignments |
| `all`         | Any               | Any             | Any              | All available alignments |

### Training Parameters

#### Basic Parameters
- `--species`: Species name for the Augustus model (required)
- `--genome`: Path to genome FASTA file
- `--quality`: Quality level of Miniprot results to use
- `--working-dir`: Directory for training files
- `--min-genes`: Minimum number of training genes required

#### Advanced Parameters
- `--flanking-dna`: Length of flanking DNA sequence (default: 4000 bp)
- `--optimization-rounds`: Number of parameter optimization rounds
- `--cpus`: Number of CPU cores to use
- `--timeout`: Training timeout in minutes

#### Quality Thresholds
- `--min-identity`: Minimum identity threshold (0.0-1.0)
- `--max-frameshifts`: Maximum frameshifts allowed
- `--max-stop-codons`: Maximum stop codons allowed

## Training Workflow

The training process consists of several automated steps:

### 1. Data Preparation
- Finds and validates Miniprot results
- Filters training data based on quality criteria
- Converts to Augustus-compatible GFF3 format
- Adds gene features and validates structure

### 2. Model Training
- Creates new Augustus species configuration
- Runs `autoAugTrain.pl` with specified parameters
- Performs parameter optimization (if requested)
- Validates trained model functionality

### 3. Quality Assessment
- Tests model prediction capabilities
- Generates comprehensive training report
- Validates all model files are present
- Provides usage recommendations

## Output Files

After successful training, you'll find these files in your working directory:

### Model Files
```
/path/to/augustus/config/species/your_species_name/
├── your_species_name_parameters.cfg      # Model parameters
├── your_species_name_exon_probs.pbl      # Exon probabilities
├── your_species_name_intron_probs.pbl    # Intron probabilities
├── your_species_name_igenic_probs.pbl    # Intergenic probabilities
└── ... (other model files)
```

### Training Files
```
results/augustus_training/
├── your_species_name_training_data.gff3  # Processed training data
├── your_species_name_training.log        # Training log
├── your_species_name_miniprot_training_report.md  # Detailed report
└── ... (intermediate files)
```

## Using Trained Models

Once training is complete, use your custom model with Augustus:

### Basic Prediction
```bash
augustus --species=your_species_name --gff3=on genome.fa > predictions.gff3
```

### Advanced Prediction
```bash
augustus \
    --species=your_species_name \
    --gff3=on \
    --uniqueGeneId=true \
    --alternatives-from-evidence=true \
    genome.fa > predictions.gff3
```

### Integration with EVM
Add your trained model to EVM configuration:
```bash
# In EVM weights file
ABINITIO_PREDICTION    augustus_your_species_name    2
```

## Best Practices

### Training Data Quality
1. **Use high-quality Miniprot results**: Start with `--quality high`
2. **Sufficient training genes**: Aim for at least 50-100 training genes
3. **Representative data**: Ensure training data covers diverse gene structures
4. **Quality filtering**: Remove problematic alignments (frameshifts, stop codons)

### Training Parameters
1. **Start simple**: Begin with default parameters
2. **Optimize gradually**: Increase optimization rounds for better accuracy
3. **Resource allocation**: Use appropriate CPU and timeout settings
4. **Backup models**: Always backup existing models before training

### Validation and Testing
1. **Check training reports**: Review detailed training reports
2. **Test predictions**: Validate model on known sequences
3. **Compare performance**: Compare with existing models
4. **Iterative improvement**: Refine training data and parameters

## Troubleshooting

### Common Issues

#### Training Data Problems
```
Error: Insufficient training genes: X < Y
```
**Solution**: Lower `--min-genes` threshold or use lower quality data

#### Augustus Configuration Issues
```
Error: autoAugTrain.pl script not found
```
**Solution**: Check Augustus installation and set correct `--augustus-scripts` path

#### Resource Limitations
```
Error: Training timed out after X minutes
```
**Solution**: Increase `--timeout` or reduce `--optimization-rounds`

### Debug Mode
For detailed troubleshooting, run with debug output:
```bash
python src/nbs_annotation/train_augustus_cli.py \
    --species debug_model \
    --debug \
    --log-file debug.log
```

### Log Analysis
Check training logs for detailed information:
```bash
# View training log
less results/augustus_training/your_species_name_training.log

# Check for errors
grep -i error results/augustus_training/your_species_name_training.log

# Monitor progress
tail -f results/augustus_training/your_species_name_training.log
```

## Configuration File Integration

Add Augustus training configuration to `config/default.yaml`:

```yaml
tools:
  augustus:
    training:
      miniprot_training:
        enabled: true
        species_name: "nbs_rice_miniprot"
        quality_filter: "high"
        min_training_genes: 20
        optimization_rounds: 1
        timeout_minutes: 240
```

## Performance Considerations

### Training Time
- **Without optimization**: 10-30 minutes
- **With optimization**: 2-6 hours
- **Large datasets**: 6+ hours

### Resource Usage
- **Memory**: 2-8 GB RAM
- **CPU**: Scales with `--cpus` parameter
- **Storage**: 100-500 MB for model files

### Optimization Tips
1. Use fewer optimization rounds for testing
2. Start with smaller, high-quality datasets
3. Use multiple CPU cores for faster training
4. Monitor resource usage during training

## Advanced Usage

### Custom Training Data
```python
# Prepare custom training data
trainer = AugustusMiniprotTrainer(config)
custom_gff = trainer.prepare_training_data("custom_miniprot.gff3")

# Use custom data for training
config.miniprot_gff_file = custom_gff
result = trainer.train_augustus_model()
```

### Batch Training
```bash
# Train multiple models with different parameters
for quality in high medium low; do
    python src/nbs_annotation/train_augustus_cli.py \
        --species "nbs_rice_${quality}" \
        --quality "${quality}" \
        --no-optimization
done
```

### Integration with Pipeline
```python
# Integrate with existing pipeline
from nbs_annotation.augustus_miniprot_trainer import AugustusMiniprotTrainer

def train_custom_augustus_model(species_name, miniprot_results):
    """Train Augustus model as part of annotation pipeline."""
    config = MiniprotTrainingConfig(
        species_name=species_name,
        miniprot_gff_file=miniprot_results,
        quality_filter="high"
    )
    
    trainer = AugustusMiniprotTrainer(config)
    result = trainer.train_augustus_model()
    
    return result.species_name if result.training_success else None
```

## API Reference

### Core Classes

#### `MiniprotTrainingConfig`
Configuration dataclass for Augustus training parameters.

#### `AugustusMiniprotTrainer`
Main class for executing Augustus training workflow.

#### `MiniprotTrainingResult`
Result object containing training outcomes and file paths.

### Key Methods

#### `train_augustus_model()`
Execute complete training workflow.

#### `prepare_training_data(miniprot_gff_file)`
Process Miniprot results into Augustus training format.

#### `validate_trained_model()`
Test trained model functionality.

## See Also

- [Augustus Documentation](http://bioinf.uni-greifswald.de/augustus/)
- [autoAugTrain.pl Manual](http://bioinf.uni-greifswald.de/augustus/binaries/tutorial-cgp/)
- [Miniprot Results Processing](miniprot_processing_guide.md)
- [EVM Integration Guide](evm_integration_guide.md)