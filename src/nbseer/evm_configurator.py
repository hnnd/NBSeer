#!/usr/bin/env python3
"""
EVidenceModeler Configuration Manager

This module provides the EVMConfigurator class for managing EVidenceModeler (EVM)
configuration, evidence weights, and execution parameters for gene annotation pipelines.
"""

import os
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass, field
import json
import subprocess
import tempfile
from datetime import datetime

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@dataclass
class EVMWeights:
    """Configuration for evidence weights in EVM"""
    protein_weights: Dict[str, int] = field(default_factory=dict)
    abinitio_weights: Dict[str, int] = field(default_factory=dict)
    transcript_weights: Dict[str, int] = field(default_factory=dict)
    other_weights: Dict[str, int] = field(default_factory=dict)
    
    def to_weights_file_content(self) -> str:
        """Generate weights.txt file content"""
        lines = []
        
        # Add protein evidence weights
        for source, weight in self.protein_weights.items():
            lines.append(f"PROTEIN\t{source}\t{weight}")
        
        # Add ab initio prediction weights
        for source, weight in self.abinitio_weights.items():
            lines.append(f"ABINITIO_PREDICTION\t{source}\t{weight}")
        
        # Add transcript evidence weights
        for source, weight in self.transcript_weights.items():
            lines.append(f"TRANSCRIPT\t{source}\t{weight}")
        
        # Add other prediction weights
        for source, weight in self.other_weights.items():
            lines.append(f"OTHER_PREDICTION\t{source}\t{weight}")
        
        return "\n".join(lines) + "\n"


@dataclass
class EVMPartitionConfig:
    """Configuration for EVM partitioning strategy"""
    segment_size: int = 100000  # Size of genome segments
    overlap_size: int = 10000   # Overlap between segments
    min_intergenic_length: int = 1000  # Minimum intergenic length
    
    def validate(self) -> bool:
        """Validate partition configuration"""
        if self.segment_size <= 0:
            raise ValueError("Segment size must be positive")
        if self.overlap_size < 0:
            raise ValueError("Overlap size cannot be negative")
        if self.overlap_size >= self.segment_size:
            raise ValueError("Overlap size must be less than segment size")
        if self.min_intergenic_length < 0:
            raise ValueError("Minimum intergenic length cannot be negative")
        return True


@dataclass
class EVMExecutionConfig:
    """Configuration for EVM execution parameters"""
    cpu_cores: int = 4
    memory_limit: str = "8G"
    temp_dir: Optional[str] = None
    output_prefix: str = "evm"
    debug_mode: bool = False
    
    def __post_init__(self):
        if self.temp_dir is None:
            self.temp_dir = tempfile.gettempdir()


@dataclass
class EVMInputFiles:
    """Configuration for EVM input files"""
    genome_file: str
    gene_predictions: List[str] = field(default_factory=list)
    protein_alignments: List[str] = field(default_factory=list)
    transcript_alignments: List[str] = field(default_factory=list)
    
    def validate_files(self) -> bool:
        """Validate that all input files exist and are readable"""
        all_files = [self.genome_file] + self.gene_predictions + \
                   self.protein_alignments + self.transcript_alignments
        
        for file_path in all_files:
            if not os.path.exists(file_path):
                raise FileNotFoundError(f"Input file not found: {file_path}")
            if not os.access(file_path, os.R_OK):
                raise PermissionError(f"Cannot read input file: {file_path}")
        
        return True


class EVMConfigurator:
    """
    EVidenceModeler Configuration Manager
    
    This class manages the configuration of EVidenceModeler (EVM) for gene annotation,
    including evidence weights, partitioning strategies, and execution parameters.
    """
    
    def __init__(self, 
                 output_dir: str = "evm_output",
                 config_name: str = "evm_config"):
        """
        Initialize EVMConfigurator
        
        Args:
            output_dir: Directory for EVM output files
            config_name: Base name for configuration files
        """
        self.output_dir = Path(output_dir)
        self.config_name = config_name
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize default configurations
        self.weights = EVMWeights()
        self.partition_config = EVMPartitionConfig()
        self.execution_config = EVMExecutionConfig()
        self.input_files = None
        
        # Set default evidence weights for NBS annotation pipeline
        self.set_default_weights()
        
        logger.info(f"EVMConfigurator initialized with output directory: {self.output_dir}")
    
    def set_default_weights(self):
        """Set default evidence weights for NBS gene annotation"""
        # Protein evidence (higher weight as it's more reliable)
        self.weights.protein_weights["miniprot"] = 10
        
        # Ab initio predictions (lower weight as they're less reliable)
        self.weights.abinitio_weights["augustus"] = 8
        
        logger.info("Set default evidence weights: miniprot=10, augustus=8")
    
    def add_protein_evidence(self, source: str, weight: int):
        """Add protein evidence source with weight"""
        self.weights.protein_weights[source] = weight
        logger.info(f"Added protein evidence: {source} with weight {weight}")
    
    def add_abinitio_evidence(self, source: str, weight: int):
        """Add ab initio prediction source with weight"""
        self.weights.abinitio_weights[source] = weight
        logger.info(f"Added ab initio evidence: {source} with weight {weight}")
    
    def add_transcript_evidence(self, source: str, weight: int):
        """Add transcript evidence source with weight"""
        self.weights.transcript_weights[source] = weight
        logger.info(f"Added transcript evidence: {source} with weight {weight}")
    
    def configure_partitioning(self, 
                             segment_size: int = 100000,
                             overlap_size: int = 10000,
                             min_intergenic_length: int = 1000):
        """Configure genome partitioning parameters"""
        self.partition_config = EVMPartitionConfig(
            segment_size=segment_size,
            overlap_size=overlap_size,
            min_intergenic_length=min_intergenic_length
        )
        self.partition_config.validate()
        logger.info(f"Configured partitioning: segment={segment_size}, overlap={overlap_size}")
    
    def configure_execution(self,
                          cpu_cores: int = 4,
                          memory_limit: str = "8G",
                          temp_dir: Optional[str] = None,
                          debug_mode: bool = False):
        """Configure EVM execution parameters"""
        self.execution_config = EVMExecutionConfig(
            cpu_cores=cpu_cores,
            memory_limit=memory_limit,
            temp_dir=temp_dir,
            debug_mode=debug_mode
        )
        logger.info(f"Configured execution: cores={cpu_cores}, memory={memory_limit}")
    
    def set_input_files(self,
                       genome_file: str,
                       gene_predictions: List[str] = None,
                       protein_alignments: List[str] = None,
                       transcript_alignments: List[str] = None):
        """Set input files for EVM"""
        self.input_files = EVMInputFiles(
            genome_file=genome_file,
            gene_predictions=gene_predictions or [],
            protein_alignments=protein_alignments or [],
            transcript_alignments=transcript_alignments or []
        )
        self.input_files.validate_files()
        logger.info(f"Set input files: genome={genome_file}, "
                   f"predictions={len(self.input_files.gene_predictions)}, "
                   f"proteins={len(self.input_files.protein_alignments)}")
    
    def create_weights_file(self, output_path: Optional[str] = None) -> str:
        """Create weights.txt file for EVM"""
        if output_path is None:
            output_path = self.output_dir / "weights.txt"
        
        weights_content = self.weights.to_weights_file_content()
        
        with open(output_path, 'w') as f:
            f.write(weights_content)
        
        logger.info(f"Created weights file: {output_path}")
        return str(output_path)
    
    def validate_gff3_sources(self, gff3_file: str) -> Dict[str, int]:
        """
        Validate that GFF3 file sources match configured weights
        
        Returns:
            Dictionary of source names and their counts
        """
        sources = {}
        
        try:
            with open(gff3_file, 'r') as f:
                for line in f:
                    if line.startswith('#') or not line.strip():
                        continue
                    
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:
                        source = parts[1]  # Second column is source
                        sources[source] = sources.get(source, 0) + 1
        
        except Exception as e:
            logger.error(f"Error validating GFF3 file {gff3_file}: {e}")
            raise
        
        return sources
    
    def validate_input_consistency(self) -> bool:
        """
        Validate that input files are consistent with configured weights
        """
        if not self.input_files:
            raise ValueError("Input files not configured")
        
        # Check protein alignments
        for protein_file in self.input_files.protein_alignments:
            sources = self.validate_gff3_sources(protein_file)
            for source in sources:
                if source not in self.weights.protein_weights:
                    logger.warning(f"Protein source '{source}' not in weights configuration")
        
        # Check gene predictions
        for prediction_file in self.input_files.gene_predictions:
            sources = self.validate_gff3_sources(prediction_file)
            for source in sources:
                if source not in self.weights.abinitio_weights:
                    logger.warning(f"Ab initio source '{source}' not in weights configuration")
        
        return True
    
    def generate_partition_command(self) -> str:
        """Generate command for partitioning EVM inputs"""
        if not self.input_files:
            raise ValueError("Input files not configured")
        
        cmd_parts = [
            "EvmUtils/partition_EVM_inputs.pl",
            f"--genome {self.input_files.genome_file}",
            f"--segmentSize {self.partition_config.segment_size}",
            f"--overlapSize {self.partition_config.overlap_size}",
            f"--partition_listing {self.output_dir}/partitions_list.out"
        ]
        
        # Add gene predictions
        for pred_file in self.input_files.gene_predictions:
            cmd_parts.append(f"--gene_predictions {pred_file}")
        
        # Add protein alignments
        for prot_file in self.input_files.protein_alignments:
            cmd_parts.append(f"--protein_alignments {prot_file}")
        
        # Add transcript alignments
        for trans_file in self.input_files.transcript_alignments:
            cmd_parts.append(f"--transcript_alignments {trans_file}")
        
        return " \\\n    ".join(cmd_parts)
    
    def generate_evm_commands(self, weights_file: str) -> str:
        """Generate command for creating EVM execution commands"""
        if not self.input_files:
            raise ValueError("Input files not configured")
        
        cmd_parts = [
            "EvmUtils/write_EVM_commands.pl",
            f"--genome {self.input_files.genome_file}",
            f"--weights {weights_file}",
            f"--output_file_name {self.execution_config.output_prefix}.out",
            f"--partitions {self.output_dir}/partitions_list.out"
        ]
        
        # Add gene predictions
        for pred_file in self.input_files.gene_predictions:
            cmd_parts.append(f"--gene_predictions {pred_file}")
        
        # Add protein alignments
        for prot_file in self.input_files.protein_alignments:
            cmd_parts.append(f"--protein_alignments {prot_file}")
        
        # Add transcript alignments
        for trans_file in self.input_files.transcript_alignments:
            cmd_parts.append(f"--transcript_alignments {trans_file}")
        
        cmd_parts.append(f"> {self.output_dir}/commands.list")
        
        return " \\\n    ".join(cmd_parts)
    
    def generate_parallel_execution_command(self) -> str:
        """Generate command for parallel EVM execution"""
        return f"ParaFly -c {self.output_dir}/commands.list -CPU {self.execution_config.cpu_cores}"
    
    def generate_recombine_command(self) -> str:
        """Generate command for recombining EVM outputs"""
        cmd_parts = [
            "EvmUtils/recombine_EVM_partial_outputs.pl",
            f"--partitions {self.output_dir}/partitions_list.out",
            f"--output_file_name {self.execution_config.output_prefix}.out"
        ]
        
        return " \\\n    ".join(cmd_parts)
    
    def generate_gff3_conversion_command(self) -> str:
        """Generate command for converting EVM output to GFF3"""
        if not self.input_files:
            raise ValueError("Input files not configured")
        
        cmd_parts = [
            "EvmUtils/convert_EVM_outputs_to_GFF3.pl",
            f"--partitions {self.output_dir}/partitions_list.out",
            f"--output {self.execution_config.output_prefix}.out",
            f"--genome {self.input_files.genome_file}"
        ]
        
        return " \\\n    ".join(cmd_parts)
    
    def save_configuration(self, config_file: Optional[str] = None) -> str:
        """Save complete configuration to JSON file"""
        if config_file is None:
            config_file = self.output_dir / f"{self.config_name}.json"
        
        config_data = {
            "created_at": datetime.now().isoformat(),
            "weights": {
                "protein": self.weights.protein_weights,
                "abinitio": self.weights.abinitio_weights,
                "transcript": self.weights.transcript_weights,
                "other": self.weights.other_weights
            },
            "partition_config": {
                "segment_size": self.partition_config.segment_size,
                "overlap_size": self.partition_config.overlap_size,
                "min_intergenic_length": self.partition_config.min_intergenic_length
            },
            "execution_config": {
                "cpu_cores": self.execution_config.cpu_cores,
                "memory_limit": self.execution_config.memory_limit,
                "temp_dir": self.execution_config.temp_dir,
                "output_prefix": self.execution_config.output_prefix,
                "debug_mode": self.execution_config.debug_mode
            },
            "input_files": {
                "genome_file": self.input_files.genome_file if self.input_files else None,
                "gene_predictions": self.input_files.gene_predictions if self.input_files else [],
                "protein_alignments": self.input_files.protein_alignments if self.input_files else [],
                "transcript_alignments": self.input_files.transcript_alignments if self.input_files else []
            }
        }
        
        with open(config_file, 'w') as f:
            json.dump(config_data, f, indent=2)
        
        logger.info(f"Saved configuration to: {config_file}")
        return str(config_file)
    
    def load_configuration(self, config_file: str):
        """Load configuration from JSON file"""
        with open(config_file, 'r') as f:
            config_data = json.load(f)
        
        # Load weights
        weights_data = config_data.get("weights", {})
        self.weights = EVMWeights(
            protein_weights=weights_data.get("protein", {}),
            abinitio_weights=weights_data.get("abinitio", {}),
            transcript_weights=weights_data.get("transcript", {}),
            other_weights=weights_data.get("other", {})
        )
        
        # Load partition config
        partition_data = config_data.get("partition_config", {})
        self.partition_config = EVMPartitionConfig(**partition_data)
        
        # Load execution config
        exec_data = config_data.get("execution_config", {})
        self.execution_config = EVMExecutionConfig(**exec_data)
        
        # Load input files
        input_data = config_data.get("input_files", {})
        if input_data.get("genome_file"):
            self.input_files = EVMInputFiles(**input_data)
        
        logger.info(f"Loaded configuration from: {config_file}")
    
    def generate_full_pipeline_script(self, script_path: Optional[str] = None) -> str:
        """Generate a complete shell script for EVM pipeline execution"""
        if script_path is None:
            script_path = self.output_dir / "run_evm_pipeline.sh"
        
        # Create weights file first
        weights_file = self.create_weights_file()
        
        # Generate script content
        script_content = f"""#!/bin/bash
# EVidenceModeler Pipeline Script
# Generated on {datetime.now().isoformat()}
# Configuration: {self.config_name}

set -e  # Exit on any error

echo "Starting EVM pipeline at $(date)"

# Configuration
OUTPUT_DIR="{self.output_dir}"
WEIGHTS_FILE="{weights_file}"
CPU_CORES={self.execution_config.cpu_cores}

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "Step 1: Partitioning genome and evidence files..."
{self.generate_partition_command()}

echo "Step 2: Generating EVM commands..."
{self.generate_evm_commands(weights_file)}

echo "Step 3: Running EVM in parallel..."
{self.generate_parallel_execution_command()}

echo "Step 4: Recombining partial outputs..."
{self.generate_recombine_command()}

echo "Step 5: Converting to GFF3 format..."
{self.generate_gff3_conversion_command()}

echo "Step 6: Finalizing output..."
find . -regex ".*{self.execution_config.output_prefix}.out.gff3" -exec cat {{}} \\; | \\
    bedtools sort -i - > "$OUTPUT_DIR/EVM.all.gff"

echo "EVM pipeline completed successfully at $(date)"
echo "Final output: $OUTPUT_DIR/EVM.all.gff"
"""
        
        with open(script_path, 'w') as f:
            f.write(script_content)
        
        # Make script executable
        os.chmod(script_path, 0o755)
        
        logger.info(f"Generated pipeline script: {script_path}")
        return str(script_path)
    
    def get_summary(self) -> Dict[str, Any]:
        """Get a summary of the current configuration"""
        return {
            "output_directory": str(self.output_dir),
            "weights_summary": {
                "protein_sources": len(self.weights.protein_weights),
                "abinitio_sources": len(self.weights.abinitio_weights),
                "transcript_sources": len(self.weights.transcript_weights),
                "total_evidence_types": (len(self.weights.protein_weights) + 
                                       len(self.weights.abinitio_weights) + 
                                       len(self.weights.transcript_weights))
            },
            "partition_config": {
                "segment_size": self.partition_config.segment_size,
                "overlap_size": self.partition_config.overlap_size
            },
            "execution_config": {
                "cpu_cores": self.execution_config.cpu_cores,
                "memory_limit": self.execution_config.memory_limit
            },
            "input_files_configured": self.input_files is not None
        }


def create_default_nbs_config(output_dir: str = "evm_output") -> EVMConfigurator:
    """
    Create a default EVMConfigurator for NBS gene annotation pipeline
    
    Args:
        output_dir: Directory for EVM output files
        
    Returns:
        Configured EVMConfigurator instance
    """
    configurator = EVMConfigurator(output_dir=output_dir, config_name="nbs_evm_config")
    
    # Configure for typical NBS annotation pipeline
    configurator.configure_partitioning(
        segment_size=100000,  # 100kb segments
        overlap_size=10000,   # 10kb overlap
        min_intergenic_length=1000
    )
    
    configurator.configure_execution(
        cpu_cores=8,
        memory_limit="16G",
        debug_mode=False
    )
    
    logger.info("Created default NBS EVM configuration")
    return configurator


if __name__ == "__main__":
    # Example usage
    configurator = create_default_nbs_config("test_evm_output")
    
    # Print configuration summary
    summary = configurator.get_summary()
    print("EVM Configuration Summary:")
    print(json.dumps(summary, indent=2))
    
    # Save configuration
    config_file = configurator.save_configuration()
    print(f"Configuration saved to: {config_file}") 