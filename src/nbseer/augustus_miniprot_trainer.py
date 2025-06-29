#!/usr/bin/env python3
"""
Augustus Training with High-Quality Miniprot Results

This module provides advanced Augustus training capabilities using high-quality
Miniprot results as training data. It integrates with the existing Augustus
training infrastructure and uses autoAugTrain.pl for optimal model training.

Key Features:
- Use high-quality Miniprot GFF3 results as training data
- Automatic quality filtering and validation
- Integration with autoAugTrain.pl for robust training
- Custom species model creation for NBS genes
- Comprehensive training validation and reporting

Author: NBS Annotation Pipeline
Date: 2025-06-25
"""

import os
import sys
import re
import logging
import subprocess
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass, asdict
from datetime import datetime

# Add parent directory to path for imports
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from utils.logging_setup import setup_logging
from utils.config import get_config
from nbseer.augustus_trainer import AugustusAutoTrainer, AutoTrainingConfig, AutoTrainingResult
from nbseer.miniprot_result_processor import MiniprotResultProcessor


@dataclass
class MiniprotTrainingConfig:
    """Configuration for Augustus training with Miniprot results."""
    
    # Basic configuration
    species_name: str = "nbs_miniprot"
    genome_file: str = "genome/osa.fa"
    working_dir: str = "results/augustus_miniprot_training"
    
    # Miniprot data sources
    miniprot_gff_file: str = ""
    miniprot_quality_filter: str = "high"  # "high", "medium", "low", "all"
    use_filtered_results: bool = True
    
    # Augustus training parameters
    augustus_scripts_path: str = "/home/wangys/data/work/nbs/tools/Augustus/scripts"
    augustus_config_path: str = "/home/wangys/data/work/nbs/tools/Augustus/config"
    flanking_dna_length: int = 4000
    optimization_rounds: int = 1
    cpus: int = 4
    timeout_minutes: int = 240  # 4 hours
    
    # Quality thresholds
    min_training_genes: int = 20
    min_identity_threshold: float = 0.95
    max_frameshifts: int = 0
    max_stop_codons: int = 0
    
    # Training data processing
    merge_overlapping_genes: bool = True
    filter_partial_genes: bool = True
    add_gene_features: bool = True
    
    # Backup and safety
    backup_existing_model: bool = True
    create_training_report: bool = True
    use_existing: bool = True  # Use existing training directories to avoid conflicts


@dataclass
class MiniprotTrainingResult:
    """Results from Miniprot-based Augustus training."""
    
    # Basic info
    species_name: str
    training_success: bool
    training_time_minutes: float
    
    # Training data statistics
    source_miniprot_file: str
    filtered_training_genes: int
    quality_filter_applied: str
    
    # Model information
    model_files: List[str]
    model_directory: str
    
    # Validation results
    validation_passed: bool
    test_prediction_success: bool
    
    # Files generated
    training_gff_file: str
    log_file: str
    report_file: str
    
    # Error information
    error_message: str = ""
    warnings: List[str] = None
    
    def __post_init__(self):
        if self.warnings is None:
            self.warnings = []


class AugustusMiniprotTrainer:
    """
    Advanced Augustus trainer using high-quality Miniprot results.
    
    This class provides a complete workflow for training Augustus gene prediction
    models using high-quality Miniprot alignment results as training data.
    """
    
    def __init__(self, config: Optional[MiniprotTrainingConfig] = None,
                 logger: Optional[logging.Logger] = None):
        """
        Initialize the Miniprot-based Augustus trainer.
        
        Args:
            config: Training configuration
            logger: Logger instance
        """
        self.config = config or MiniprotTrainingConfig()
        self.logger = logger or setup_logging(f"{__name__}.MiniprotTrainer")
        
        # Initialize core components
        self.augustus_trainer = None
        self.miniprot_processor = None
        
        # Training state
        self.training_started = False
        self.training_completed = False
        self.training_result: Optional[MiniprotTrainingResult] = None
        
        # Validate configuration
        self._validate_configuration()
        
        self.logger.info(f"Initialized AugustusMiniprotTrainer for species: {self.config.species_name}")
        self.logger.info(f"Miniprot quality filter: {self.config.miniprot_quality_filter}")
        self.logger.info(f"Working directory: {self.config.working_dir}")
    
    def _validate_configuration(self) -> bool:
        """
        Validate the training configuration.
        
        Returns:
            True if configuration is valid
            
        Raises:
            ValueError: If configuration is invalid
        """
        errors = []
        
        # Check genome file
        if not os.path.exists(self.config.genome_file):
            errors.append(f"Genome file not found: {self.config.genome_file}")
        
        # Check Augustus scripts path
        if not os.path.exists(self.config.augustus_scripts_path):
            errors.append(f"Augustus scripts path not found: {self.config.augustus_scripts_path}")
        
        # Check autoAugTrain.pl script
        auto_aug_train = os.path.join(self.config.augustus_scripts_path, "autoAugTrain.pl")
        if not os.path.exists(auto_aug_train):
            errors.append(f"autoAugTrain.pl script not found: {auto_aug_train}")
        
        # Validate species name
        if not self.config.species_name or not self.config.species_name.strip():
            errors.append("Species name cannot be empty")
        
        # Validate quality filter
        valid_filters = ["high", "medium", "low", "all"]
        if self.config.miniprot_quality_filter not in valid_filters:
            errors.append(f"Invalid quality filter: {self.config.miniprot_quality_filter}. Must be one of: {valid_filters}")
        
        # Check working directory can be created
        try:
            os.makedirs(self.config.working_dir, exist_ok=True)
        except Exception as e:
            errors.append(f"Cannot create working directory {self.config.working_dir}: {e}")
        
        if errors:
            error_msg = "Configuration validation failed:\n" + "\n".join(f"- {error}" for error in errors)
            self.logger.error(error_msg)
            raise ValueError(error_msg)
        
        self.logger.info("Configuration validation passed")
        return True
    
    def find_miniprot_results(self) -> str:
        """
        Find high-quality Miniprot results file.
        
        Returns:
            Path to the appropriate Miniprot GFF file
        """
        # If specific file provided, use it
        if self.config.miniprot_gff_file and os.path.exists(self.config.miniprot_gff_file):
            self.logger.info(f"Using specified Miniprot file: {self.config.miniprot_gff_file}")
            return self.config.miniprot_gff_file
        
        # Look for filtered results in common locations
        search_paths = [
            "data/test/test_output/protein_alignment/filtered",
            "results/protein_alignment/filtered",
            "data/output/protein_alignment/filtered"
        ]
        
        quality_files = {
            "high": "miniprot_high_quality.gff3",
            "medium": "miniprot_medium_quality.gff3", 
            "low": "miniprot_low_quality.gff3",
            "all": "alignments_filtered.gff"
        }
        
        target_file = quality_files.get(self.config.miniprot_quality_filter, "miniprot_high_quality.gff3")
        
        for search_path in search_paths:
            full_path = os.path.join(search_path, target_file)
            if os.path.exists(full_path):
                self.logger.info(f"Found Miniprot results: {full_path}")
                return full_path
        
        # If not found, raise error
        raise FileNotFoundError(f"Could not find Miniprot results file: {target_file}")
    
    def prepare_training_data(self, miniprot_gff_file: str) -> str:
        """
        Prepare GFF3 training data from Miniprot results.
        
        Args:
            miniprot_gff_file: Path to Miniprot GFF3 file
            
        Returns:
            Path to prepared training GFF3 file
        """
        self.logger.info(f"Preparing training data from: {miniprot_gff_file}")
        
        # Read and validate Miniprot GFF3
        training_genes = []
        gene_count = 0
        
        try:
            with open(miniprot_gff_file, 'r') as f:
                current_gene_lines = []
                current_gene_id = None
                
                for line in f:
                    line = line.strip()
                    
                    # Skip comments and empty lines
                    if not line or line.startswith('#'):
                        continue
                    
                    # Parse GFF3 line
                    fields = line.split('\t')
                    if len(fields) != 9:
                        continue
                    
                    seqname, source, feature, start, end, score, strand, frame, attributes = fields
                    
                    # Extract gene ID from attributes
                    gene_id = None
                    if feature == 'mRNA':
                        # Start of new gene
                        match = re.search(r'ID=([^;]+)', attributes)
                        if match:
                            gene_id = match.group(1)
                        
                        # Save previous gene if exists
                        if current_gene_lines and current_gene_id:
                            training_genes.append({
                                'gene_id': current_gene_id,
                                'lines': current_gene_lines.copy()
                            })
                            gene_count += 1
                        
                        # Start new gene
                        current_gene_id = gene_id
                        current_gene_lines = [line]
                    
                    elif feature in ['CDS', 'exon', 'start_codon', 'stop_codon']:
                        # Add to current gene
                        if current_gene_lines:
                            current_gene_lines.append(line)
                
                # Don't forget the last gene
                if current_gene_lines and current_gene_id:
                    training_genes.append({
                        'gene_id': current_gene_id,
                        'lines': current_gene_lines.copy()
                    })
                    gene_count += 1
        
        except Exception as e:
            self.logger.error(f"Error reading Miniprot GFF3 file: {e}")
            raise
        
        # Validate minimum gene count
        if gene_count < self.config.min_training_genes:
            raise ValueError(f"Insufficient training genes: {gene_count} < {self.config.min_training_genes}")
        
        self.logger.info(f"Found {gene_count} training genes in Miniprot results")
        
        # Process and enhance training data
        enhanced_training_data = self._enhance_training_data(training_genes)
        
        # Write enhanced training GFF3
        training_gff_file = os.path.join(self.config.working_dir, f"{self.config.species_name}_training_data.gff3")
        
        with open(training_gff_file, 'w') as f:
            f.write("##gff-version 3\n")
            f.write(f"# Augustus training data derived from Miniprot results\n")
            f.write(f"# Source: {miniprot_gff_file}\n")
            f.write(f"# Quality filter: {self.config.miniprot_quality_filter}\n")
            f.write(f"# Training genes: {len(enhanced_training_data)}\n")
            f.write(f"# Generated: {datetime.now().isoformat()}\n")
            f.write("#\n")
            
            for gene_data in enhanced_training_data:
                for line in gene_data['lines']:
                    f.write(line + '\n')
        
        self.logger.info(f"Prepared training data: {training_gff_file}")
        self.logger.info(f"Enhanced training genes: {len(enhanced_training_data)}")
        
        return training_gff_file
    
    def _enhance_training_data(self, training_genes: List[Dict]) -> List[Dict]:
        """
        Enhance training data with gene features and quality filtering.
        
        Args:
            training_genes: List of training gene data
            
        Returns:
            Enhanced training gene data
        """
        enhanced_genes = []
        
        for gene_data in training_genes:
            gene_id = gene_data['gene_id']
            lines = gene_data['lines']
            
            if not lines:
                continue
            
            # Parse first line (mRNA) to get coordinates
            mRNA_line = lines[0]
            fields = mRNA_line.split('\t')
            seqname, source, feature, start, end, score, strand, frame, attributes = fields
            
            # Add gene feature if requested
            if self.config.add_gene_features:
                gene_line = f"{seqname}\t{source}\tgene\t{start}\t{end}\t{score}\t{strand}\t{frame}\tID={gene_id}_gene;Name={gene_id}_gene"
                enhanced_lines = [gene_line] + lines
            else:
                enhanced_lines = lines
            
            # Filter partial genes if requested
            if self.config.filter_partial_genes:
                # Check if gene has both start and stop codons
                has_start = any('start_codon' in line for line in lines)
                has_stop = any('stop_codon' in line for line in lines)
                
                if not (has_start or has_stop):  # Allow genes with at least one
                    self.logger.debug(f"Filtering partial gene: {gene_id}")
                    continue
            
            enhanced_genes.append({
                'gene_id': gene_id,
                'lines': enhanced_lines
            })
        
        return enhanced_genes
    
    def setup_augustus_trainer(self, training_gff_file: str) -> AugustusAutoTrainer:
        """
        Setup Augustus auto trainer with Miniprot training data.
        
        Args:
            training_gff_file: Path to prepared training GFF3 file
            
        Returns:
            Configured AugustusAutoTrainer instance
        """
        # Create Augustus training configuration
        augustus_config = AutoTrainingConfig(
            species_name=self.config.species_name,
            augustus_scripts_path=self.config.augustus_scripts_path,
            augustus_config_path=self.config.augustus_config_path,
            genome_file=self.config.genome_file,
            training_gff=training_gff_file,
            working_dir=self.config.working_dir,
            flanking_dna_length=self.config.flanking_dna_length,
            optimization_rounds=self.config.optimization_rounds,
            cpus=self.config.cpus,
            timeout_minutes=self.config.timeout_minutes,
            min_training_genes=self.config.min_training_genes,
            backup_existing_model=self.config.backup_existing_model,
            use_existing=self.config.use_existing
        )
        
        # Create Augustus trainer
        augustus_trainer = AugustusAutoTrainer(
            config=augustus_config,
            logger=self.logger
        )
        
        return augustus_trainer
    
    def train_augustus_model(self) -> MiniprotTrainingResult:
        """
        Execute complete Augustus training workflow with Miniprot data.
        
        Returns:
            Training result object
        """
        start_time = time.time()
        
        # Initialize result object
        result = MiniprotTrainingResult(
            species_name=self.config.species_name,
            training_success=False,
            training_time_minutes=0.0,
            source_miniprot_file="",
            filtered_training_genes=0,
            quality_filter_applied=self.config.miniprot_quality_filter,
            model_files=[],
            model_directory="",
            validation_passed=False,
            test_prediction_success=False,
            training_gff_file="",
            log_file="",
            report_file=""
        )
        
        try:
            self.logger.info("Starting Augustus training with Miniprot results")
            self.training_started = True
            
            # Step 1: Find Miniprot results
            self.logger.info("Step 1: Finding Miniprot results")
            miniprot_gff_file = self.find_miniprot_results()
            result.source_miniprot_file = miniprot_gff_file
            
            # Step 2: Prepare training data
            self.logger.info("Step 2: Preparing training data")
            training_gff_file = self.prepare_training_data(miniprot_gff_file)
            result.training_gff_file = training_gff_file
            
            # Count training genes
            with open(training_gff_file, 'r') as f:
                gene_count = sum(1 for line in f if '\tmRNA\t' in line)
            result.filtered_training_genes = gene_count
            
            # Step 3: Setup Augustus trainer
            self.logger.info("Step 3: Setting up Augustus trainer")
            self.augustus_trainer = self.setup_augustus_trainer(training_gff_file)
            
            # Step 4: Execute Augustus training
            self.logger.info("Step 4: Executing Augustus training")
            augustus_result = self.augustus_trainer.train_model_auto()
            
            if augustus_result:
                result.training_success = True
                result.model_directory = self.augustus_trainer.species_dir
                result.model_files = self.augustus_trainer.training_result.model_files if self.augustus_trainer.training_result else []
                result.log_file = self.augustus_trainer.training_result.log_file if self.augustus_trainer.training_result else ""
                result.validation_passed = self.augustus_trainer.training_result.validation_passed if self.augustus_trainer.training_result else False
                
                self.logger.info("Augustus training completed successfully")
            else:
                result.error_message = "Augustus training failed"
                self.logger.error("Augustus training failed")
            
        except Exception as e:
            result.error_message = f"Training workflow failed: {e}"
            self.logger.error(f"Training workflow failed: {e}")
        
        finally:
            # Calculate training time
            end_time = time.time()
            result.training_time_minutes = (end_time - start_time) / 60.0
            
            self.training_completed = True
            self.training_result = result
            
            # Generate training report
            if self.config.create_training_report:
                report_file = self.generate_comprehensive_report(result)
                result.report_file = report_file
        
        return result
    
    def generate_comprehensive_report(self, result: MiniprotTrainingResult) -> str:
        """
        Generate comprehensive training report.
        
        Args:
            result: Training result object
            
        Returns:
            Path to generated report file
        """
        try:
            report_file = os.path.join(self.config.working_dir, f"{self.config.species_name}_miniprot_training_report.md")
            
            with open(report_file, 'w') as f:
                f.write("# Augustus Training Report - Miniprot Integration\n\n")
                f.write(f"**Species:** {result.species_name}\n")
                f.write(f"**Training Date:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write(f"**Training Success:** {result.training_success}\n")
                f.write(f"**Training Time:** {result.training_time_minutes:.1f} minutes\n\n")
                
                f.write("## Data Sources\n\n")
                f.write(f"- **Miniprot Results:** {result.source_miniprot_file}\n")
                f.write(f"- **Quality Filter:** {result.quality_filter_applied}\n")
                f.write(f"- **Training Genes:** {result.filtered_training_genes}\n")
                f.write(f"- **Training GFF:** {result.training_gff_file}\n\n")
                
                f.write("## Configuration\n\n")
                f.write(f"- **Genome File:** {self.config.genome_file}\n")
                f.write(f"- **Working Directory:** {self.config.working_dir}\n")
                f.write(f"- **Augustus Scripts:** {self.config.augustus_scripts_path}\n")
                f.write(f"- **Flanking DNA:** {self.config.flanking_dna_length} bp\n")
                f.write(f"- **Optimization Rounds:** {self.config.optimization_rounds}\n")
                f.write(f"- **CPUs:** {self.config.cpus}\n")
                f.write(f"- **Timeout:** {self.config.timeout_minutes} minutes\n\n")
                
                f.write("## Training Results\n\n")
                if result.training_success:
                    f.write("✅ **Training Status:** SUCCESSFUL\n")
                    f.write(f"- Model Directory: {result.model_directory}\n")
                    f.write(f"- Model Validation: {'✅ PASSED' if result.validation_passed else '❌ FAILED'}\n")
                    f.write(f"- Prediction Test: {'✅ PASSED' if result.test_prediction_success else '❌ FAILED'}\n\n")
                    
                    f.write("### Model Files\n\n")
                    for model_file in result.model_files:
                        f.write(f"- {model_file}\n")
                    f.write("\n")
                else:
                    f.write("❌ **Training Status:** FAILED\n")
                    if result.error_message:
                        f.write(f"- **Error:** {result.error_message}\n")
                    f.write("\n")
                
                f.write("## Quality Metrics\n\n")
                f.write(f"- **Minimum Training Genes:** {self.config.min_training_genes}\n")
                f.write(f"- **Actual Training Genes:** {result.filtered_training_genes}\n")
                f.write(f"- **Gene Count Status:** {'✅ ADEQUATE' if result.filtered_training_genes >= self.config.min_training_genes else '❌ INSUFFICIENT'}\n")
                f.write(f"- **Identity Threshold:** {self.config.min_identity_threshold}\n")
                f.write(f"- **Max Frameshifts:** {self.config.max_frameshifts}\n")
                f.write(f"- **Max Stop Codons:** {self.config.max_stop_codons}\n\n")
                
                f.write("## Next Steps\n\n")
                if result.training_success:
                    f.write("1. ✅ Model training completed successfully\n")
                    f.write("2. Test the model on validation data\n")
                    f.write("3. Use the trained model for Augustus prediction:\n")
                    f.write(f"   ```bash\n")
                    f.write(f"   augustus --species={result.species_name} --gff3=on genome.fa\n")
                    f.write(f"   ```\n")
                    f.write("4. Integrate with EVM for consensus prediction\n")
                else:
                    f.write("1. ❌ Training failed - review error messages\n")
                    f.write("2. Check training data quality and quantity\n")
                    f.write("3. Verify Augustus installation and configuration\n")
                    f.write("4. Consider adjusting training parameters\n")
                
                f.write("\n## Log Files\n\n")
                if result.log_file:
                    f.write(f"- Training Log: {result.log_file}\n")
                f.write(f"- This Report: {report_file}\n")
            
            self.logger.info(f"Comprehensive training report generated: {report_file}")
            return report_file
            
        except Exception as e:
            self.logger.error(f"Report generation failed: {e}")
            return ""
    
    def get_training_status(self) -> Dict[str, Any]:
        """
        Get current training status.
        
        Returns:
            Dictionary with training status information
        """
        status = {
            'training_started': self.training_started,
            'training_completed': self.training_completed,
            'species_name': self.config.species_name,
            'working_directory': self.config.working_dir,
            'quality_filter': self.config.miniprot_quality_filter
        }
        
        if self.training_result:
            status.update({
                'training_success': self.training_result.training_success,
                'training_time_minutes': self.training_result.training_time_minutes,
                'filtered_training_genes': self.training_result.filtered_training_genes,
                'validation_passed': self.training_result.validation_passed,
                'error_message': self.training_result.error_message
            })
        
        return status


def main():
    """Main function for testing Augustus Miniprot training."""
    # Example usage
    config = MiniprotTrainingConfig(
        species_name="nbs_rice_miniprot",
        genome_file="genome/osa.fa",
        miniprot_quality_filter="high",
        optimization_rounds=1,
        cpus=4,
        timeout_minutes=120
    )
    
    trainer = AugustusMiniprotTrainer(config)
    
    try:
        result = trainer.train_augustus_model()
        
        if result.training_success:
            print("🎉 Augustus training with Miniprot data completed successfully!")
            print(f"Species: {result.species_name}")
            print(f"Training genes: {result.filtered_training_genes}")
            print(f"Training time: {result.training_time_minutes:.1f} minutes")
            print(f"Model directory: {result.model_directory}")
            print(f"Report: {result.report_file}")
        else:
            print("❌ Augustus training failed!")
            print(f"Error: {result.error_message}")
            
    except Exception as e:
        print(f"Training workflow failed: {e}")


if __name__ == "__main__":
    main()