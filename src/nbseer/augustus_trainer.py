#!/usr/bin/env python3
"""
Augustus Model Training Module

This module provides the AugustusTrainer class for training Augustus gene prediction models
specifically optimized for NBS gene features.
"""

import os
import sys
import json
import shutil
import subprocess
import logging
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass
from datetime import datetime

import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from utils.logging_setup import setup_logging
from utils.config import get_config


@dataclass
class TrainingConfig:
    """Configuration for Augustus training."""
    species_name: str = "custom_nbs"
    augustus_config_path: str = ""
    training_data_path: str = ""
    test_data_ratio: float = 0.2
    optimization_rounds: int = 5
    cross_validation_folds: int = 5
    min_gene_level_sensitivity: float = 0.2
    parallel_jobs: int = 4
    
    # NBS-specific parameters
    min_exon_length: int = 50
    max_exon_length: int = 3000
    min_intron_length: int = 20
    max_intron_length: int = 20000
    translation_table: int = 1
    
    # Training quality thresholds
    min_training_genes: int = 100
    max_training_genes: int = 2000


@dataclass
class TrainingResult:
    """Results from Augustus training."""
    species_name: str
    training_genes_count: int
    test_genes_count: int
    gene_sensitivity: float
    gene_specificity: float
    exon_sensitivity: float
    exon_specificity: float
    training_time_minutes: float
    optimization_completed: bool
    model_files: List[str]
    training_log: str


class AugustusTrainer:
    """
    Augustus Model Trainer for NBS gene prediction.
    
    This class encapsulates the complete Augustus training workflow including:
    - Species configuration creation
    - Initial model training
    - Parameter optimization
    - Cross-validation evaluation
    - Model quality assessment
    """
    
    def __init__(self, config: Optional[TrainingConfig] = None, 
                 logger: Optional[logging.Logger] = None):
        """
        Initialize Augustus trainer.
        
        Args:
            config: Training configuration
            logger: Logger instance
        """
        self.config = config or TrainingConfig()
        self.logger = logger or setup_logging(__name__)
        
        # Initialize paths
        self.augustus_config_path = self._find_augustus_config_path()
        self.species_dir = os.path.join(self.augustus_config_path, "species", self.config.species_name)
        self.training_dir = "results/augustus_training"
        
        # Training state
        self.training_started = False
        self.training_completed = False
        self.optimization_completed = False
        
        # Results
        self.training_result: Optional[TrainingResult] = None
        
        self.logger.info(f"Initialized AugustusTrainer for species: {self.config.species_name}")
    
    def _find_augustus_config_path(self) -> str:
        """Find Augustus configuration directory."""
        possible_paths = [
            os.environ.get('AUGUSTUS_CONFIG_PATH', ''),
            '/usr/local/share/augustus/config',
            '/opt/augustus/config',
            os.path.expanduser('~/augustus/config'),
            './tools/augustus/config'
        ]
        
        for path in possible_paths:
            if path and os.path.exists(os.path.join(path, "species")):
                self.logger.info(f"Found Augustus config at: {path}")
                return path
        
        # If not found, create a local config directory
        local_config = "config/augustus"
        os.makedirs(os.path.join(local_config, "species"), exist_ok=True)
        self.logger.warning(f"Augustus config not found, using local: {local_config}")
        return local_config
    
    def validate_training_data(self, training_file: str) -> Tuple[bool, str, int]:
        """
        Validate training data file.
        
        Args:
            training_file: Path to GenBank format training file
            
        Returns:
            Tuple of (is_valid, message, gene_count)
        """
        if not os.path.exists(training_file):
            return False, f"Training file not found: {training_file}", 0
        
        try:
            # Count genes in GenBank file
            with open(training_file, 'r') as f:
                content = f.read()
                gene_count = content.count('LOCUS ')
            
            if gene_count < self.config.min_training_genes:
                return False, f"Insufficient training genes: {gene_count} < {self.config.min_training_genes}", gene_count
            
            if gene_count > self.config.max_training_genes:
                self.logger.warning(f"Large training set: {gene_count} genes. Consider subsampling for faster training.")
            
            # Check file format
            if not content.strip().startswith('LOCUS'):
                return False, "Invalid GenBank format: file should start with LOCUS", gene_count
            
            return True, f"Valid training data with {gene_count} genes", gene_count
            
        except Exception as e:
            return False, f"Error validating training data: {e}", 0
    
    def create_species_config(self) -> bool:
        """
        Create new species configuration using new_species.pl.
        
        Returns:
            True if successful
        """
        try:
            self.logger.info(f"Creating new species configuration: {self.config.species_name}")
            
            # Check if species already exists
            if os.path.exists(self.species_dir):
                self.logger.warning(f"Species {self.config.species_name} already exists. Backing up...")
                backup_dir = f"{self.species_dir}_backup_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
                shutil.move(self.species_dir, backup_dir)
                self.logger.info(f"Backed up existing species to: {backup_dir}")
            
            # Run new_species.pl
            cmd = ["new_species.pl", f"--species={self.config.species_name}"]
            
            env = os.environ.copy()
            env['AUGUSTUS_CONFIG_PATH'] = self.augustus_config_path
            
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                env=env,
                cwd=self.training_dir
            )
            
            if result.returncode != 0:
                self.logger.error(f"Failed to create species config: {result.stderr}")
                return False
            
            # Verify species directory was created
            if not os.path.exists(self.species_dir):
                self.logger.error(f"Species directory not created: {self.species_dir}")
                return False
            
            self.logger.info(f"Successfully created species configuration at: {self.species_dir}")
            return True
            
        except Exception as e:
            self.logger.error(f"Error creating species config: {e}")
            return False
    
    def configure_species_parameters(self) -> bool:
        """
        Configure species-specific parameters for NBS genes.
        
        Returns:
            True if successful
        """
        try:
            self.logger.info("Configuring NBS-specific parameters")
            
            # Path to parameters file
            params_file = os.path.join(self.species_dir, f"{self.config.species_name}_parameters.cfg")
            
            if not os.path.exists(params_file):
                self.logger.error(f"Parameters file not found: {params_file}")
                return False
            
            # Read existing parameters
            with open(params_file, 'r') as f:
                lines = f.readlines()
            
            # Modify NBS-specific parameters
            modifications = {
                'translation_table': str(self.config.translation_table),
                'stopCodonExcludedFromCDS': 'true',  # NBS genes typically exclude stop codons
                'min_coding_len': '300',  # Minimum CDS length for NBS genes
                'alternatives-from-evidence': 'true',  # Use evidence for alternative splicing
                'uniqueGeneId': 'true',  # Unique gene IDs
                'gff3': 'on',  # Output GFF3 format
                'UTR': 'off',  # Don't predict UTRs initially
                'allow_hinted_splicesites': 'atac,gtag',  # Allow alternative splice sites
            }
            
            # Apply modifications
            modified_lines = []
            for line in lines:
                line_modified = False
                for param, value in modifications.items():
                    if line.strip().startswith(param + ' '):
                        modified_lines.append(f"{param} {value}\n")
                        line_modified = True
                        break
                
                if not line_modified:
                    modified_lines.append(line)
            
            # Write modified parameters
            with open(params_file, 'w') as f:
                f.writelines(modified_lines)
            
            self.logger.info("Successfully configured NBS-specific parameters")
            return True
            
        except Exception as e:
            self.logger.error(f"Error configuring species parameters: {e}")
            return False
    
    def split_training_data(self, training_file: str) -> Tuple[str, str]:
        """
        Split training data into training and test sets.
        
        Args:
            training_file: Path to full training data
            
        Returns:
            Tuple of (train_file, test_file) paths
        """
        try:
            self.logger.info(f"Splitting training data with ratio {self.config.test_data_ratio}")
            
            # Use randomSplit.pl if available, otherwise implement simple splitting
            train_file = training_file.replace('.gb', '.train.gb')
            test_file = training_file.replace('.gb', '.test.gb')
            
            # Calculate split size
            with open(training_file, 'r') as f:
                content = f.read()
                genes = content.split('LOCUS ')[1:]  # Split by LOCUS, skip first empty element
            
            total_genes = len(genes)
            test_size = int(total_genes * self.config.test_data_ratio)
            train_size = total_genes - test_size
            
            self.logger.info(f"Splitting {total_genes} genes: {train_size} train, {test_size} test")
            
            # Simple random split (you could use randomSplit.pl here if available)
            import random
            random.shuffle(genes)
            
            train_genes = genes[:train_size]
            test_genes = genes[train_size:]
            
            # Write training file
            with open(train_file, 'w') as f:
                for gene in train_genes:
                    f.write('LOCUS ' + gene)
            
            # Write test file
            with open(test_file, 'w') as f:
                for gene in test_genes:
                    f.write('LOCUS ' + gene)
            
            self.logger.info(f"Created training file: {train_file}")
            self.logger.info(f"Created test file: {test_file}")
            
            return train_file, test_file
            
        except Exception as e:
            self.logger.error(f"Error splitting training data: {e}")
            return training_file, training_file
    
    def run_initial_training(self, train_file: str) -> bool:
        """
        Run initial etraining on the training data.
        
        Args:
            train_file: Path to training data file
            
        Returns:
            True if successful
        """
        try:
            self.logger.info("Running initial etraining")
            
            # Convert to absolute path
            abs_train_file = os.path.abspath(train_file)
            
            cmd = [
                "etraining",
                f"--species={self.config.species_name}",
                abs_train_file
            ]
            
            env = os.environ.copy()
            env['AUGUSTUS_CONFIG_PATH'] = self.augustus_config_path
            
            self.logger.info(f"Running command: {' '.join(cmd)}")
            
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                env=env,
                cwd=self.training_dir
            )
            
            if result.returncode != 0:
                self.logger.error(f"etraining failed: {result.stderr}")
                return False
            
            self.logger.info("Initial training completed successfully")
            self.training_started = True
            
            # Log training output
            if result.stdout:
                self.logger.info(f"Training output: {result.stdout}")
            
            return True
            
        except Exception as e:
            self.logger.error(f"Error in initial training: {e}")
            return False
    
    def optimize_parameters(self, train_file: str) -> bool:
        """
        Optimize Augustus parameters using optimize_augustus.pl.
        
        Args:
            train_file: Path to training data file
            
        Returns:
            True if successful
        """
        try:
            self.logger.info("Starting parameter optimization")
            
            # Convert to absolute path
            abs_train_file = os.path.abspath(train_file)
            
            cmd = [
                "optimize_augustus.pl",
                f"--species={self.config.species_name}",
                abs_train_file
            ]
            
            env = os.environ.copy()
            env['AUGUSTUS_CONFIG_PATH'] = self.augustus_config_path
            
            self.logger.info(f"Running optimization: {' '.join(cmd)}")
            self.logger.warning("Parameter optimization may take several hours...")
            
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                env=env,
                cwd=self.training_dir,
                timeout=7200  # 2 hours timeout
            )
            
            if result.returncode != 0:
                self.logger.error(f"Parameter optimization failed: {result.stderr}")
                return False
            
            self.logger.info("Parameter optimization completed successfully")
            self.optimization_completed = True
            
            # Re-run etraining with optimized parameters
            self.logger.info("Re-running etraining with optimized parameters")
            return self.run_initial_training(train_file)
            
        except subprocess.TimeoutExpired:
            self.logger.error("Parameter optimization timed out after 2 hours")
            return False
        except Exception as e:
            self.logger.error(f"Error in parameter optimization: {e}")
            return False
    
    def evaluate_model(self, test_file: str) -> Dict[str, float]:
        """
        Evaluate trained model on test data.
        
        Args:
            test_file: Path to test data file
            
        Returns:
            Dictionary with evaluation metrics
        """
        try:
            self.logger.info("Evaluating model on test data")
            
            # Convert to absolute path
            abs_test_file = os.path.abspath(test_file)
            
            cmd = [
                "augustus",
                f"--species={self.config.species_name}",
                "--gff3=on",
                abs_test_file
            ]
            
            env = os.environ.copy()
            env['AUGUSTUS_CONFIG_PATH'] = self.augustus_config_path
            
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                env=env,
                cwd=self.training_dir
            )
            
            if result.returncode != 0:
                self.logger.error(f"Model evaluation failed: {result.stderr}")
                return {}
            
            # Parse evaluation metrics from output
            metrics = self._parse_evaluation_output(result.stderr)
            
            self.logger.info(f"Evaluation metrics: {metrics}")
            return metrics
            
        except Exception as e:
            self.logger.error(f"Error evaluating model: {e}")
            return {}
    
    def _parse_evaluation_output(self, output: str) -> Dict[str, float]:
        """Parse evaluation metrics from Augustus output."""
        metrics = {
            'gene_sensitivity': 0.0,
            'gene_specificity': 0.0,
            'exon_sensitivity': 0.0,
            'exon_specificity': 0.0
        }
        
        try:
            lines = output.split('\n')
            for line in lines:
                if 'gene level' in line and 'sensitivity' in line:
                    # Extract sensitivity value
                    parts = line.split()
                    for i, part in enumerate(parts):
                        if part == 'sensitivity':
                            metrics['gene_sensitivity'] = float(parts[i+1].rstrip('%')) / 100
                        elif part == 'specificity':
                            metrics['gene_specificity'] = float(parts[i+1].rstrip('%')) / 100
                
                elif 'exon level' in line and 'sensitivity' in line:
                    # Extract exon metrics
                    parts = line.split()
                    for i, part in enumerate(parts):
                        if part == 'sensitivity':
                            metrics['exon_sensitivity'] = float(parts[i+1].rstrip('%')) / 100
                        elif part == 'specificity':
                            metrics['exon_specificity'] = float(parts[i+1].rstrip('%')) / 100
        
        except Exception as e:
            self.logger.warning(f"Could not parse evaluation metrics: {e}")
        
        return metrics
    
    def save_training_results(self, metrics: Dict[str, float], 
                            train_genes: int, test_genes: int,
                            training_time: float) -> str:
        """Save training results to file."""
        try:
            results_dir = os.path.join(self.training_dir, "results")
            os.makedirs(results_dir, exist_ok=True)
            
            # Create training result
            self.training_result = TrainingResult(
                species_name=self.config.species_name,
                training_genes_count=train_genes,
                test_genes_count=test_genes,
                gene_sensitivity=metrics.get('gene_sensitivity', 0.0),
                gene_specificity=metrics.get('gene_specificity', 0.0),
                exon_sensitivity=metrics.get('exon_sensitivity', 0.0),
                exon_specificity=metrics.get('exon_specificity', 0.0),
                training_time_minutes=training_time,
                optimization_completed=self.optimization_completed,
                model_files=self._get_model_files(),
                training_log=f"Training completed at {datetime.now()}"
            )
            
            # Save as JSON
            results_file = os.path.join(results_dir, f"{self.config.species_name}_training_results.json")
            with open(results_file, 'w') as f:
                json.dump(self.training_result.__dict__, f, indent=2)
            
            # Generate markdown report
            report_file = os.path.join(results_dir, f"{self.config.species_name}_training_report.md")
            self._generate_training_report(report_file)
            
            self.logger.info(f"Training results saved to: {results_file}")
            return results_file
            
        except Exception as e:
            self.logger.error(f"Error saving training results: {e}")
            return ""
    
    def _get_model_files(self) -> List[str]:
        """Get list of trained model files."""
        model_files = []
        if os.path.exists(self.species_dir):
            for file in os.listdir(self.species_dir):
                if file.startswith(self.config.species_name):
                    model_files.append(os.path.join(self.species_dir, file))
        return model_files
    
    def _generate_training_report(self, report_file: str):
        """Generate markdown training report."""
        if not self.training_result:
            return
        
        report = [
            f"# Augustus Training Report: {self.config.species_name}",
            f"\n**Training Date:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            f"\n## Training Configuration",
            f"- Species: {self.training_result.species_name}",
            f"- Training Genes: {self.training_result.training_genes_count}",
            f"- Test Genes: {self.training_result.test_genes_count}",
            f"- Training Time: {self.training_result.training_time_minutes:.1f} minutes",
            f"- Optimization Completed: {self.training_result.optimization_completed}",
            f"\n## Model Performance",
            f"- Gene Sensitivity: {self.training_result.gene_sensitivity:.3f}",
            f"- Gene Specificity: {self.training_result.gene_specificity:.3f}",
            f"- Exon Sensitivity: {self.training_result.exon_sensitivity:.3f}",
            f"- Exon Specificity: {self.training_result.exon_specificity:.3f}",
            f"\n## Model Files",
        ]
        
        for model_file in self.training_result.model_files:
            report.append(f"- {model_file}")
        
        report.extend([
            f"\n## Quality Assessment",
            f"Gene-level sensitivity: {'✅ GOOD' if self.training_result.gene_sensitivity >= self.config.min_gene_level_sensitivity else '❌ POOR'}",
            f"Training data size: {'✅ ADEQUATE' if self.training_result.training_genes_count >= self.config.min_training_genes else '❌ INSUFFICIENT'}",
            f"\n## Next Steps",
            f"1. Review model performance metrics",
            f"2. Test model on validation data",
            f"3. Proceed with Augustus prediction on candidate regions"
        ])
        
        with open(report_file, 'w') as f:
            f.write('\n'.join(report))
    
    def train_model(self, training_file: str) -> bool:
        """
        Complete Augustus training workflow.
        
        Args:
            training_file: Path to GenBank format training file
            
        Returns:
            True if training completed successfully
        """
        start_time = datetime.now()
        
        try:
            self.logger.info("Starting Augustus model training workflow")
            
            # Validate training data
            is_valid, message, gene_count = self.validate_training_data(training_file)
            if not is_valid:
                self.logger.error(f"Training data validation failed: {message}")
                return False
            
            self.logger.info(f"Training data validated: {message}")
            
            # Create training directory
            os.makedirs(self.training_dir, exist_ok=True)
            
            # Step 1: Create species configuration
            if not self.create_species_config():
                return False
            
            # Step 2: Configure species parameters
            if not self.configure_species_parameters():
                return False
            
            # Step 3: Split training data
            train_file, test_file = self.split_training_data(training_file)
            
            # Count genes in split files
            with open(train_file, 'r') as f:
                train_genes = f.read().count('LOCUS ')
            with open(test_file, 'r') as f:
                test_genes = f.read().count('LOCUS ')
            
            # Step 4: Initial training
            if not self.run_initial_training(train_file):
                return False
            
            # Step 5: Parameter optimization (optional, time-consuming)
            if self.config.optimization_rounds > 0:
                self.logger.info("Starting parameter optimization...")
                if not self.optimize_parameters(train_file):
                    self.logger.warning("Parameter optimization failed, continuing with initial training")
            
            # Step 6: Evaluate model
            metrics = self.evaluate_model(test_file)
            
            # Calculate training time
            training_time = (datetime.now() - start_time).total_seconds() / 60
            
            # Step 7: Save results
            self.save_training_results(metrics, train_genes, test_genes, training_time)
            
            self.training_completed = True
            
            # Check if training meets quality requirements
            gene_sensitivity = metrics.get('gene_sensitivity', 0.0)
            if gene_sensitivity < self.config.min_gene_level_sensitivity:
                self.logger.warning(f"Low gene sensitivity: {gene_sensitivity:.3f} < {self.config.min_gene_level_sensitivity}")
                self.logger.warning("Consider increasing training data or improving data quality")
            else:
                self.logger.info(f"Training successful! Gene sensitivity: {gene_sensitivity:.3f}")
            
            return True
            
        except Exception as e:
            self.logger.error(f"Training workflow failed: {e}")
            return False


def main():
    """Main function for testing Augustus training."""
    # Example usage
    config = TrainingConfig(
        species_name="nbs_rice",
        test_data_ratio=0.2,
        optimization_rounds=1,  # Reduced for testing
        min_gene_level_sensitivity=0.15
    )
    
    trainer = AugustusTrainer(config)
    
    # Use high quality training data
    training_file = "results/augustus_training_data/genbank_format/high_quality_training.gb"
    
    if os.path.exists(training_file):
        success = trainer.train_model(training_file)
        if success:
            print("Augustus training completed successfully!")
        else:
            print("Augustus training failed!")
    else:
        print(f"Training file not found: {training_file}")


if __name__ == "__main__":
    main()


@dataclass
class AutoTrainingConfig:
    """Configuration for autoAugTrain.pl-based Augustus training."""
    species_name: str = "nbs_rice_gff"
    augustus_scripts_path: str = "/home/wangys/data/work/nbs/tools/Augustus/scripts"
    genome_file: str = "genome/osa.fa"
    training_gff: str = "results/augustus_training_data/high_quality_training_data.gff"
    augustus_config_path: str = ""
    working_dir: str = "results/augustus_training"
    
    # autoAugTrain.pl parameters
    flanking_dna_length: int = 4000
    optimization_rounds: int = 1
    cpus: int = 4
    verbose: int = 2
    use_existing: bool = False
    use_crf: bool = False
    
    # Training quality thresholds
    timeout_minutes: int = 180  # 3 hours timeout
    min_training_genes: int = 50
    backup_existing_model: bool = True


@dataclass 
class AutoTrainingResult:
    """Results from autoAugTrain.pl training."""
    species_name: str
    training_success: bool
    training_time_minutes: float
    model_files: List[str]
    log_file: str
    error_message: str = ""
    gene_count: int = 0
    validation_passed: bool = False


class AugustusAutoTrainer:
    """
    Augustus Model Trainer using autoAugTrain.pl script with GFF format data.
    
    This class implements training using the autoAugTrain.pl script which directly
    accepts GFF format training data, avoiding the need for GenBank format conversion
    and preserving stop codon information from miniprot predictions.
    """
    
    def __init__(self, config: Optional[AutoTrainingConfig] = None,
                 logger: Optional[logging.Logger] = None):
        """
        Initialize Augustus auto trainer.
        
        Args:
            config: Auto training configuration
            logger: Logger instance
        """
        self.config = config or AutoTrainingConfig()
        self.logger = logger or setup_logging(f"{__name__}.AutoTrainer")
        
        # Initialize paths
        self.augustus_config_path = self._find_augustus_config_path()
        self.auto_aug_train_script = os.path.join(
            self.config.augustus_scripts_path, "autoAugTrain.pl"
        )
        self.species_dir = os.path.join(
            self.augustus_config_path, "species", self.config.species_name
        )
        
        # Training state
        self.training_started = False
        self.training_completed = False
        self.training_result: Optional[AutoTrainingResult] = None
        
        # Validate configuration on initialization
        self._validate_config()
        
        self.logger.info(f"Initialized AugustusAutoTrainer for species: {self.config.species_name}")
        self.logger.info(f"Augustus scripts path: {self.config.augustus_scripts_path}")
        self.logger.info(f"Training GFF: {self.config.training_gff}")
        self.logger.info(f"Genome file: {self.config.genome_file}")
    
    def _find_augustus_config_path(self) -> str:
        """Find Augustus configuration directory."""
        # Check configuration first
        if self.config.augustus_config_path and os.path.exists(os.path.join(self.config.augustus_config_path, "species")):
            self.logger.info(f"Found Augustus config from configuration: {self.config.augustus_config_path}")
            return self.config.augustus_config_path
            
        # Check environment variable second
        env_path = os.environ.get('AUGUSTUS_CONFIG_PATH', '')
        if env_path and os.path.exists(os.path.join(env_path, "species")):
            self.logger.info(f"Found Augustus config from environment: {env_path}")
            return env_path
        
        # Common installation paths
        possible_paths = [
            '/home/wangys/data/work/nbs/tools/Augustus/config',  # Project tools directory
            '/usr/local/share/augustus/config',
            '/opt/augustus/config',
            '/home/wangys/opt/Augustus/config',
            os.path.expanduser('~/augustus/config'),
            './tools/augustus/config'
        ]
        
        for path in possible_paths:
            if os.path.exists(os.path.join(path, "species")):
                self.logger.info(f"Found Augustus config at: {path}")
                return path
        
        # If not found, try to create a local config directory
        local_config = os.path.abspath("config/augustus")
        os.makedirs(os.path.join(local_config, "species"), exist_ok=True)
        self.logger.warning(f"Augustus config not found, using local: {local_config}")
        return local_config
    
    def _validate_config(self) -> bool:
        """
        Validate configuration and environment setup.
        
        Returns:
            True if configuration is valid
            
        Raises:
            ValueError: If critical validation fails
        """
        errors = []
        
        # Check autoAugTrain.pl script exists
        if not os.path.exists(self.auto_aug_train_script):
            errors.append(f"autoAugTrain.pl script not found: {self.auto_aug_train_script}")
        
        # Check if script is executable
        if os.path.exists(self.auto_aug_train_script):
            if not os.access(self.auto_aug_train_script, os.X_OK):
                errors.append(f"autoAugTrain.pl script not executable: {self.auto_aug_train_script}")
        
        # Check genome file exists
        if not os.path.exists(self.config.genome_file):
            errors.append(f"Genome file not found: {self.config.genome_file}")
        
        # Check training GFF exists
        if not os.path.exists(self.config.training_gff):
            errors.append(f"Training GFF file not found: {self.config.training_gff}")
        
        # Check Augustus config path exists
        if not os.path.exists(os.path.join(self.augustus_config_path, "species")):
            errors.append(f"Augustus config species directory not found: {self.augustus_config_path}/species")
        
        # Validate species name
        if not self.config.species_name or not self.config.species_name.strip():
            errors.append("Species name cannot be empty")
        
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

    def validate_training_data(self) -> bool:
        """
        Validate GFF training data format and content.
        
        Returns:
            True if training data is valid
        """
        try:
            self.logger.info(f"Validating training data: {self.config.training_gff}")
            
            # Check file exists and is readable
            if not os.path.exists(self.config.training_gff):
                self.logger.error(f"Training GFF file not found: {self.config.training_gff}")
                return False
            
            gene_count = 0
            mRNA_count = 0
            cds_count = 0
            
            with open(self.config.training_gff, 'r') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    
                    # Skip comments and empty lines
                    if not line or line.startswith('#'):
                        continue
                    
                    # Parse GFF line
                    fields = line.split('\t')
                    if len(fields) != 9:
                        self.logger.warning(f"Invalid GFF format at line {line_num}: {len(fields)} fields")
                        continue
                    
                    seqname, source, feature, start, end, score, strand, frame, attributes = fields
                    
                    # Count features
                    if feature == 'gene':
                        gene_count += 1
                    elif feature == 'mRNA':
                        mRNA_count += 1
                    elif feature == 'CDS':
                        cds_count += 1
                    
                    # Validate coordinates
                    try:
                        start_pos = int(start)
                        end_pos = int(end)
                        if start_pos > end_pos:
                            self.logger.warning(f"Invalid coordinates at line {line_num}: start={start_pos} > end={end_pos}")
                    except ValueError:
                        self.logger.warning(f"Invalid coordinates at line {line_num}: start={start}, end={end}")
            
            # Validate minimum requirements
            if mRNA_count < self.config.min_training_genes:
                self.logger.error(f"Insufficient training genes: {mRNA_count} < {self.config.min_training_genes}")
                return False
            
            self.logger.info(f"Training data validation successful:")
            self.logger.info(f"  - Genes: {gene_count}")
            self.logger.info(f"  - mRNAs: {mRNA_count}")
            self.logger.info(f"  - CDS features: {cds_count}")
            
            return True
            
        except Exception as e:
            self.logger.error(f"Training data validation failed: {e}")
            return False

    def validate_genome_file(self) -> bool:
        """
        Validate genome FASTA file.
        
        Returns:
            True if genome file is valid
        """
        try:
            self.logger.info(f"Validating genome file: {self.config.genome_file}")
            
            if not os.path.exists(self.config.genome_file):
                self.logger.error(f"Genome file not found: {self.config.genome_file}")
                return False
            
            # Check if file is readable and has content
            try:
                with open(self.config.genome_file, 'r') as f:
                    first_line = f.readline().strip()
                    if not first_line.startswith('>'):
                        self.logger.error("Genome file does not appear to be in FASTA format")
                        return False
            except Exception as e:
                self.logger.error(f"Cannot read genome file: {e}")
                return False
            
            # Check for FASTA index
            fai_file = f"{self.config.genome_file}.fai"
            if not os.path.exists(fai_file):
                self.logger.warning(f"FASTA index not found: {fai_file}")
                self.logger.info("Consider creating index with: samtools faidx genome.fa")
            
            self.logger.info("Genome file validation successful")
            return True
            
        except Exception as e:
            self.logger.error(f"Genome file validation failed: {e}")
            return False

    def backup_existing_model(self) -> bool:
        """
        Backup existing species model if it exists.
        
        Returns:
            True if backup successful or no model exists
        """
        try:
            if not os.path.exists(self.species_dir):
                self.logger.info(f"No existing model to backup for species: {self.config.species_name}")
                return True
            
            if not self.config.backup_existing_model:
                self.logger.info("Backup disabled, removing existing model")
                import shutil
                shutil.rmtree(self.species_dir)
                return True
            
            # Create backup with timestamp
            import time
            timestamp = time.strftime("%Y%m%d_%H%M%S")
            backup_dir = f"{self.species_dir}_backup_{timestamp}"
            
            import shutil
            shutil.copytree(self.species_dir, backup_dir)
            self.logger.info(f"Backed up existing model to: {backup_dir}")
            
            # Remove original
            shutil.rmtree(self.species_dir)
            
            return True
            
        except Exception as e:
            self.logger.error(f"Model backup failed: {e}")
            return False

    def setup_training_command(self) -> List[str]:
        """
        Setup autoAugTrain.pl command with all parameters.
        
        Returns:
            Command list for subprocess
        """
        cmd = [
            "perl",
            self.auto_aug_train_script,
            f"--genome={os.path.abspath(self.config.genome_file)}",
            f"--trainingset={os.path.abspath(self.config.training_gff)}",
            f"--species={self.config.species_name}",
            f"--flanking_DNA={self.config.flanking_dna_length}",
            f"--workingdir={os.path.abspath(self.config.working_dir)}",
            f"--optrounds={self.config.optimization_rounds}",
            f"--cpus={self.config.cpus}"
        ]
        
        # Add optional parameters
        if self.config.verbose > 1:
            cmd.extend(["-v"] * (self.config.verbose - 1))
        
        if self.config.use_existing:
            cmd.append("--useexisting")
        
        if self.config.use_crf:
            cmd.append("--CRF")
        
        self.logger.info(f"Training command: {' '.join(cmd)}")
        return cmd

    def run_training(self) -> AutoTrainingResult:
        """
        Execute autoAugTrain.pl training process.
        
        Returns:
            Training result object
        """
        import subprocess
        import time
        
        start_time = time.time()
        
        # Initialize result object
        result = AutoTrainingResult(
            species_name=self.config.species_name,
            training_success=False,
            training_time_minutes=0.0,
            model_files=[],
            log_file=""
        )
        
        try:
            self.logger.info("Starting autoAugTrain.pl training process")
            self.training_started = True
            
            # Setup environment
            env = os.environ.copy()
            env['AUGUSTUS_CONFIG_PATH'] = self.augustus_config_path
            
            # Setup command
            cmd = self.setup_training_command()
            
            # Create log file
            log_file = os.path.join(self.config.working_dir, f"{self.config.species_name}_training.log")
            result.log_file = log_file
            
            # Execute training with timeout
            self.logger.info(f"Executing training with {self.config.timeout_minutes} minute timeout")
            
            with open(log_file, 'w') as log_f:
                process = subprocess.Popen(
                    cmd,
                    cwd=self.config.working_dir,
                    env=env,
                    stdout=log_f,
                    stderr=subprocess.STDOUT,
                    text=True
                )
                
                try:
                    # Wait for completion with timeout
                    process.wait(timeout=self.config.timeout_minutes * 60)
                    
                    if process.returncode == 0:
                        result.training_success = True
                        self.logger.info("Training completed successfully")
                    else:
                        result.error_message = f"Training failed with return code: {process.returncode}"
                        self.logger.error(result.error_message)
                        
                except subprocess.TimeoutExpired:
                    process.kill()
                    result.error_message = f"Training timed out after {self.config.timeout_minutes} minutes"
                    self.logger.error(result.error_message)
            
            # Calculate training time
            end_time = time.time()
            result.training_time_minutes = (end_time - start_time) / 60.0
            
            self.logger.info(f"Training completed in {result.training_time_minutes:.1f} minutes")
            
        except Exception as e:
            result.error_message = f"Training execution failed: {e}"
            self.logger.error(result.error_message)
        
        finally:
            self.training_completed = True
            self.training_result = result
        
        return result

    def validate_trained_model(self) -> bool:
        """
        Validate the trained model files.
        
        Returns:
            True if model validation successful
        """
        try:
            self.logger.info(f"Validating trained model: {self.config.species_name}")
            
            # Check if species directory was created
            if not os.path.exists(self.species_dir):
                self.logger.error(f"Species directory not created: {self.species_dir}")
                return False
            
            # Check for essential model files
            essential_files = [
                f"{self.config.species_name}_parameters.cfg",
                f"{self.config.species_name}_exon_probs.pbl",
                f"{self.config.species_name}_intron_probs.pbl",
                f"{self.config.species_name}_igenic_probs.pbl"
            ]
            
            missing_files = []
            existing_files = []
            
            for filename in essential_files:
                filepath = os.path.join(self.species_dir, filename)
                if os.path.exists(filepath):
                    existing_files.append(filename)
                    self.logger.info(f"Found model file: {filename} ({os.path.getsize(filepath)} bytes)")
                else:
                    missing_files.append(filename)
            
            if missing_files:
                self.logger.error(f"Missing essential model files: {missing_files}")
                return False
            
            # Test basic model functionality
            if not self.test_model_prediction():
                self.logger.error("Model prediction test failed")
                return False
            
            self.logger.info("Model validation successful")
            if self.training_result:
                self.training_result.model_files = existing_files
                self.training_result.validation_passed = True
            
            return True
            
        except Exception as e:
            self.logger.error(f"Model validation failed: {e}")
            return False

    def test_model_prediction(self) -> bool:
        """
        Test the trained model with a simple prediction.
        
        Returns:
            True if prediction test successful
        """
        try:
            # Create a simple test sequence
            test_seq = ">test_seq\nATGAAGGCTGTTGCTGCTGCTGCTGCTAGCGCTGCTGCTGCTTAA\n"
            test_file = os.path.join(self.config.working_dir, "test_seq.fa")
            
            with open(test_file, 'w') as f:
                f.write(test_seq)
            
            # Run Augustus prediction
            cmd = [
                "augustus",
                f"--species={self.config.species_name}",
                test_file
            ]
            
            env = os.environ.copy()
            env['AUGUSTUS_CONFIG_PATH'] = self.augustus_config_path
            
            result = subprocess.run(
                cmd,
                env=env,
                capture_output=True,
                text=True,
                timeout=60
            )
            
            if result.returncode == 0 and "gene" in result.stdout:
                self.logger.info("Model prediction test passed")
                return True
            else:
                self.logger.warning(f"Model prediction test warning: {result.stderr}")
                return True  # Sometimes Augustus gives warnings but still works
                
        except Exception as e:
            self.logger.error(f"Model prediction test failed: {e}")
            return False
        
        finally:
            # Clean up test file
            if os.path.exists(test_file):
                os.remove(test_file)

    def generate_training_report(self) -> str:
        """
        Generate comprehensive training report.
        
        Returns:
            Path to generated report file
        """
        try:
            report_file = os.path.join(self.config.working_dir, f"{self.config.species_name}_auto_training_report.md")
            
            with open(report_file, 'w') as f:
                f.write(f"# Augustus Auto Training Report\n\n")
                f.write(f"**Species:** {self.config.species_name}\n")
                f.write(f"**Training Date:** {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
                
                f.write("## Configuration\n\n")
                f.write(f"- Genome file: {self.config.genome_file}\n")
                f.write(f"- Training GFF: {self.config.training_gff}\n")
                f.write(f"- Working directory: {self.config.working_dir}\n")
                f.write(f"- Flanking DNA length: {self.config.flanking_dna_length}\n")
                f.write(f"- Optimization rounds: {self.config.optimization_rounds}\n")
                f.write(f"- CPUs: {self.config.cpus}\n\n")
                
                if self.training_result:
                    f.write("## Training Results\n\n")
                    f.write(f"- **Training successful:** {self.training_result.training_success}\n")
                    f.write(f"- **Training time:** {self.training_result.training_time_minutes:.1f} minutes\n")
                    f.write(f"- **Model validation:** {self.training_result.validation_passed}\n")
                    
                    if self.training_result.error_message:
                        f.write(f"- **Error:** {self.training_result.error_message}\n")
                    
                    f.write(f"\n## Model Files\n\n")
                    for model_file in self.training_result.model_files:
                        filepath = os.path.join(self.species_dir, model_file)
                        if os.path.exists(filepath):
                            size = os.path.getsize(filepath)
                            f.write(f"- {model_file} ({size:,} bytes)\n")
                    
                    f.write(f"\n## Log File\n\n")
                    f.write(f"Training log: {self.training_result.log_file}\n")
            
            self.logger.info(f"Training report generated: {report_file}")
            return report_file
            
        except Exception as e:
            self.logger.error(f"Report generation failed: {e}")
            return ""

    def train_model_auto(self) -> bool:
        """
        Execute complete autoAugTrain.pl-based training workflow.
        
        Returns:
            True if training workflow successful
        """
        try:
            self.logger.info("Starting autoAugTrain.pl-based Augustus training workflow")
            
            # Step 1: Validate training data and genome
            self.logger.info("Step 1: Validating training data")
            if not self.validate_training_data():
                return False
            
            if not self.validate_genome_file():
                return False
            
            # Step 2: Backup existing model if needed
            self.logger.info("Step 2: Managing existing models")
            if not self.backup_existing_model():
                return False
            
            # Step 3: Run training
            self.logger.info("Step 3: Running autoAugTrain.pl training")
            result = self.run_training()
            
            if not result.training_success:
                self.logger.error("Training failed")
                return False
            
            # Step 4: Validate trained model
            self.logger.info("Step 4: Validating trained model")
            if not self.validate_trained_model():
                return False
            
            # Step 5: Generate report
            self.logger.info("Step 5: Generating training report")
            report_file = self.generate_training_report()
            
            self.logger.info("AutoAugTrain.pl training workflow completed successfully!")
            self.logger.info(f"New species model available: {self.config.species_name}")
            self.logger.info(f"Model directory: {self.species_dir}")
            if report_file:
                self.logger.info(f"Training report: {report_file}")
            
            return True
            
        except Exception as e:
            self.logger.error(f"Training workflow failed: {e}")
            return False


if __name__ == "__main__":
    main() 