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
# Clone the repository
git clone https://github.com/hnnd/nbseer.git
cd nbseer

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
    --proteins proteins.fa \
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

## 📚 Documentation

| Document | Description | Audience |
|----------|-------------|----------|
| **[📋 Complete Documentation](./NBS_ANNOTATION_PIPELINE_DOCUMENTATION.md)** | Complete technical documentation | Developers & Advanced Users |
| **[📖 Project Overview](./PROJECT_OVERVIEW.md)** | Documentation navigation & overview | All Users |
| **[⚙️ Installation Guide](./INSTALLATION_GUIDE.md)** | Detailed installation guide | System Administrators |
| **[🔧 Configuration Guide](./CONFIG_UPDATE_GUIDE.md)** | Configuration & customization | Advanced Users |

### Quick Navigation

- 🆕 **New Users**: Start with this README → [Installation Guide](./INSTALLATION_GUIDE.md)
- 🔧 **System Admins**: [Installation Guide](./INSTALLATION_GUIDE.md) → [Tools Setup](./TOOLS_SETUP_GUIDE.md)
- 👨‍💻 **Developers**: [Complete Documentation](./NBS_ANNOTATION_PIPELINE_DOCUMENTATION.md) → Source Code (`src/`)
- 🧬 **Bioinformaticians**: [Augustus Training Guide](./docs/augustus_training_guide.md) → [Gene ID Renaming](./docs/gene_id_renaming_guide.md)

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
| 100 Mb | ~2 hours | ~8 GB | ~200-500 |
| 500 Mb | ~8 hours | ~24 GB | ~1000-2000 |
| 3 Gb | ~24 hours | ~64 GB | ~5000-10000 |

## 🧪 Example Workflows

### Workflow 1: Complete Pipeline
```bash
# Step 1: Run annotation pipeline
python -m src.nbseer.main \
    --genome rice_genome.fa \
    --proteins rice_proteins.fa \
    --output rice_results/ \
    --threads 24

# Step 2: Rename genes with species prefix
python -m src.nbseer.post_evm_renamer \
    --input rice_results/final_annotations.gff3 \
    --output rice_results/rice_nbs_genes.gff3 \
    --species rice_nipponbare

# Step 3: Generate summary report
python -c "
import json
with open('rice_results/pipeline_results.json') as f:
    stats = json.load(f)
    print(f'Total genes predicted: {stats[\"total_genes\"]}')
    print(f'Average gene length: {stats[\"average_gene_length\"]:.1f} bp')
"
```

### Workflow 2: Custom Training
```bash
# Train Augustus model for your species
python -m src.nbseer.train_augustus_cli \
    --genome your_genome.fa \
    --proteins your_proteins.fa \
    --species-name your_species \
    --training-type miniprot \
    --output training_results/

# Use custom model for annotation
python -m src.nbseer.main \
    --genome your_genome.fa \
    --proteins your_proteins.fa \
    --augustus-species your_species \
    --output results/
```

### Workflow 3: Quality Control Focus
```bash
# High-stringency annotation
python -m src.nbseer.main \
    --genome genome.fa \
    --proteins proteins.fa \
    --quality-filter strict \
    --min-identity 0.98 \
    --min-gene-length 500 \
    --output high_quality_results/
```

## 🤝 Support & Contribution

### Getting Help
- 📖 **Documentation**: [Complete Documentation](./NBS_ANNOTATION_PIPELINE_DOCUMENTATION.md)
- 🐛 **Issues**: [GitHub Issues](https://github.com/hnnd/nbseer/issues)
- 💬 **Discussions**: [GitHub Discussions](https://github.com/hnnd/nbseer/discussions)

### Contributing
We welcome contributions! Please see our [Contributing Guidelines](./CONTRIBUTING.md) for details.

### Citation
If you use NBSeer in your research, please cite:

```bibtex
@software{nbseer_2025,
  title={NBSeer: Intelligent Plant NBS Gene Annotation Tool},
  author={NBSeer Development Team},
  year={2025},
  version={0.1.0},
  url={https://github.com/hnnd/nbseer}
}
```

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

**🔗 Quick Links**:
[📋 Complete Docs](./NBS_ANNOTATION_PIPELINE_DOCUMENTATION.md) | 
[📖 Project Overview](./PROJECT_OVERVIEW.md) | 
[⚙️ Installation](./INSTALLATION_GUIDE.md) | 
[🔧 Configuration](./CONFIG_UPDATE_GUIDE.md) | 
[🧬 Augustus Training](./docs/augustus_training_guide.md)

---

**Version**: 0.1.0 | **Updated**: 2025-06-29 | **Team**: NBSeer Development Team