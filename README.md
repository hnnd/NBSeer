# NBSeer: Intelligent Plant NBS Gene Annotation Tool

[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Code style: ruff](https://img.shields.io/badge/code%20style-ruff-000000.svg)](https://github.com/astral-sh/ruff)

NBSeer is an intelligent computational tool for accurately annotating plant NBS (Nucleotide-Binding Site) disease resistance genes using multiple bioinformatics tools and evidence integration.

## 🧬 Overview

Plant NBS disease resistance genes are crucial for plant immunity but are challenging to annotate accurately due to their complex structure, high copy numbers, and rapid evolution. NBSeer combines multiple state-of-the-art tools to achieve accurate gene structure prediction and annotation.

### ✨ Key Features

- 🔍 **NLR Gene Localization**: Uses NLR-Annotator to identify candidate regions
- 🧬 **Gene Structure Prediction**: Employs Augustus with custom training
- 🔗 **Protein Alignment**: Utilizes miniprot for protein-based evidence
- 📊 **Evidence Integration**: Combines multiple predictions using EVidenceModeler
- ✅ **Quality Control**: Comprehensive validation and quality assessment
- 🚀 **Parallel Processing**: Optimized for large-scale genomic data
- 🏷️ **Gene ID Renaming**: Species-specific systematic naming
- ⚙️ **Automated Installation**: One-click setup for all dependencies

## 🚀 Quick Start

### 1. Installation

```bash
# Create conda env
conda create -n nbseer python=3.11
conda activate nbseer

# Clone the repository
git clone https://github.com/hnnd/nbseer.git
cd nbseer

# Install dependencies
./install_system_deps.sh

# Automated installation
./setup_tools.sh

# Load environment
source setup_env.sh

# Verify installation
python -m src.nbseer.main --check-tools
```

### 2. Basic Usage

```bash
# Run full pipeline
python -m src.nbseer.main \
    --genome genome.fa \
    --proteins data/prgdb/prg_nbs.fasta \
    --output results/ \
    --threads 16

# Use species-specific gene naming
python -m src.nbseer.post_evm_renamer \
    --input results/final_annotations.gff3 \
    --output results/renamed_annotations.gff3 \
    --species rice_nipponbare
```

### 3. Output

```
results/
├── final_annotations.gff3      # Final annotations
├── gene_sequences.fasta        # Gene sequences  
├── protein_sequences.fasta     # Protein sequences
├── pipeline_results.json       # Statistics
├── quality_report.html         # Quality report
└── renamed_annotations.gff3    # Renamed annotations
```

## 🛠️ System Requirements

### Hardware
- **CPU**: 16+ cores recommended
- **Memory**: 32GB+ RAM
- **Storage**: 100GB+ available space

### Software
- **OS**: Linux (Ubuntu 20.04+) or macOS
- **Python**: 3.11+
- **Java**: 8+ (for NLR-Annotator)

## 🔧 Integrated Tools

| Tool | Version | Purpose | Status |
|------|---------|---------|--------|
| **Augustus** | 3.5.0 | Gene prediction | ✅ Auto-install |
| **miniprot** | 0.17 | Protein alignment | ✅ Auto-install |
| **EVidenceModeler** | master | Evidence integration | ✅ Auto-install |
| **ParaFly** | master | Parallel execution | ✅ Auto-install |
| **NLR-Annotator** | 2.1b | NLR identification | ✅ Auto-install |

## 📊 Performance

| Genome Size | Processing Time | Memory Usage | Predicted Genes |
|-------------|-----------------|--------------|-----------------|
| 120 Mb | ~0.5 hours | ~8 GB | ~200-500 |
| 380 Mb | ~1 hours | ~24 GB | ~500-1000 |
| 2.8 Gb | ~4 hours | ~64 GB | ~500-1000 |


## Citation
If you use NBSeer in your research, please cite:

```bibtex
@software{nbseer_2025,
  title={NBSeer: Intelligent Plant NBS Gene Annotation Tool},
  author={LIU Qian, LIU Pingbo, DAI Liangying, CHEN Wu, LIU Tianbo, LU Yaoxiong, WANG Yunsheng},
  year={2025},
  version={0.1.0},
  url={https://github.com/hnnd/nbseer}
}
```

## 📄 License

This project is licensed under the MIT License.

---

**Version**: 0.1.0 | **Updated**: 2025-06-29 