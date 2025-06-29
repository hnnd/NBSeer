# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is **NBSeer** - a sophisticated bioinformatics tool for accurately annotating plant NBS (Nucleotide-Binding Site) disease resistance genes. NBSeer integrates multiple specialized tools to provide comprehensive gene structure prediction and annotation.

### Core Capabilities
- Automated re-annotation of plant NBS genes from genomic sequences
- Integration of multiple evidence sources (ab initio, protein alignment, RNA-seq)
- Custom Augustus training for improved gene structure prediction
- Parallel processing support for large-scale analyses
- Comprehensive quality control and validation

## Quick Start Commands

### Initial Setup (Required)
```bash
# Install all required bioinformatics tools
./setup_tools.sh

# Set up environment variables
source setup_env.sh

# Verify tool installation
./verify_tools.sh
```

### Running the Pipeline
```bash
# Basic annotation run
python -m src.nbseer.main --config config/example_config.yaml

# Development mode with verbose output
python -m src.nbseer.main --config config/example_config.yaml --verbose --debug
```

### Development Commands
```bash
# Run all tests
pytest

# Run specific test categories
pytest tests/unit/
pytest tests/integration/
pytest tests/functional/

# Code formatting and linting
ruff format .
ruff check .

# Type checking
mypy src/
```

## Architecture Overview

### Pipeline Stages
1. **Input Processing** - Genome sequence and annotation parsing
2. **NLR Detection** - NLR-Annotator identifies candidate regions
3. **Gene Structure Prediction** - Augustus with custom training models
4. **Protein Evidence** - miniprot alignment of known NBS proteins
5. **Evidence Integration** - EVidenceModeler combines all evidence
6. **Post-processing** - Quality filtering and output formatting

### Key Tool Integration
- **NLR-Annotator v2.1b**: Initial NBS gene localization
- **Augustus v3.5.0**: Gene structure prediction with custom species training
- **miniprot v0.17**: High-accuracy protein-to-genome alignment
- **EVidenceModeler**: Weighted evidence integration
- **ParaFly**: Parallel command execution

### Data Flow
```
Genome FASTA → NLR-Annotator → Candidate Regions →
Augustus Training → Gene Prediction → miniprot Alignment →
EVidenceModeler → Final Annotations → GFF3/FASTA Output
```

## Key Files and Directories

### Core Source Code
- `src/nbseer/main.py` - Main pipeline entry point
- `src/nbseer/pipeline/` - Core pipeline stages
- `src/nbseer/tools/` - Tool-specific wrappers
- `src/nbseer/config/` - Configuration management
- `src/nbseer/utils/` - Shared utilities

### Configuration
- `config/example_config.yaml` - Template configuration
- `evm_config/` - EVidenceModeler weight configurations
- `pyproject.toml` - Python project dependencies and metadata

### Tool Integration
- `tools/` - Integrated bioinformatics software
- `setup_tools.sh` - Automated tool installation script
- `verify_tools.sh` - Tool validation script

### Testing and Examples
- `tests/` - Comprehensive test suite (unit, integration, functional)
- `examples/` - Usage examples and test datasets

## Configuration Management

### Main Config Structure (YAML)
```yaml
input:
  genome: path/to/genome.fasta
  proteins: path/to/proteins.fasta
  
output:
  directory: results/
  prefix: species_name

tools:
  augustus_species: custom_trained_model
  evm_weights: config/weights.txt
  
parameters:
  min_gene_length: 300
  max_intron_length: 10000
```

### Environment Variables
- Set via `source setup_env.sh`
- `AUGUSTUS_CONFIG_PATH` - Augustus configuration directory
- `EVM_HOME` - EVidenceModeler installation path
- Tool-specific PATH additions

## Development Guidelines

### Code Organization
- Pipeline stages are modular and can be run independently
- Tool wrappers provide consistent interfaces to external software
- Configuration is centralized and validation is enforced
- Comprehensive logging at multiple verbosity levels

### Testing Strategy
- Unit tests for individual functions and classes
- Integration tests for tool interactions
- Functional tests with real biological data
- Use pytest fixtures for common test data

### Error Handling
- Custom exception classes for different error types
- Graceful handling of tool failures with informative messages
- Automatic cleanup of temporary files
- Progress tracking and resumption capabilities

## Common Issues and Solutions

### Tool Installation Problems
- Run `./verify_tools.sh` to diagnose installation issues
- Check `TOOLS_INSTALLATION_REPORT.md` for known issues
- Augustus training requires species-specific parameter files

### Pipeline Execution Issues
- Ensure sufficient disk space for intermediate files
- Check memory requirements for large genomes (>1GB)
- Validate input file formats before running
- Use `--debug` flag for detailed error messages

### Performance Optimization
- Use `--threads` parameter for parallel processing
- Adjust EVidenceModeler chunk sizes for memory constraints
- Consider splitting large genomes into chromosomes
- Monitor temp directory usage during execution

## Dependencies and Requirements

### System Requirements
- Linux/Unix environment (tested on Ubuntu 20.04+)
- Python 3.11+
- 8GB+ RAM recommended
- 50GB+ disk space for large genomes

### Python Dependencies
- Managed via pyproject.toml with uv lock file
- Key packages: BioPython, pandas, PyYAML, pytest
- Development dependencies include ruff, mypy, black

### External Tools
- All tools are automatically installed via setup_tools.sh
- Source installations with specific versions for reproducibility
- Custom Augustus training models for plant species