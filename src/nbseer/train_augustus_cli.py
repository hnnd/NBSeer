#!/usr/bin/env python3
"""
Augustus Training CLI Tool

Command-line interface for training Augustus models using high-quality
Miniprot results. This tool provides an easy way to train custom Augustus
species models for NBS gene prediction.

Usage:
    python train_augustus_cli.py [options]

Examples:
    # Train with high-quality Miniprot results
    python train_augustus_cli.py --species nbs_rice_v2 --quality high
    
    # Train with medium quality, more CPUs
    python train_augustus_cli.py --species nbs_rice_v3 --quality medium --cpus 8
    
    # Quick training for testing (no optimization)
    python train_augustus_cli.py --species test_model --quality high --no-optimization

Author: NBS Annotation Pipeline
Date: 2025-06-25
"""

import os
import sys
import argparse
import logging
from pathlib import Path

# Add parent directory to path for imports
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from utils.logging_setup import setup_logging
from utils.config import get_config
from nbs_annotation.augustus_miniprot_trainer import (
    AugustusMiniprotTrainer, 
    MiniprotTrainingConfig
)


def setup_argument_parser():
    """Setup command line argument parser."""
    parser = argparse.ArgumentParser(
        description="Train Augustus models using high-quality Miniprot results",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s --species nbs_rice_v2 --quality high
  %(prog)s --species nbs_rice_v3 --quality medium --cpus 8 --timeout 360
  %(prog)s --species test_model --quality high --no-optimization --verbose
        """
    )
    
    # Basic configuration
    parser.add_argument(
        "--species", "-s",
        required=True,
        help="Species name for the Augustus model (e.g., nbs_rice_v2)"
    )
    
    parser.add_argument(
        "--genome", "-g",
        default="genome/osa.fa",
        help="Path to genome FASTA file (default: genome/osa.fa)"
    )
    
    parser.add_argument(
        "--miniprot-file", "-m",
        help="Path to specific Miniprot GFF3 file (auto-detected if not provided)"
    )
    
    parser.add_argument(
        "--quality", "-q",
        choices=["high", "medium", "low", "all"],
        default="high",
        help="Quality level of Miniprot results to use (default: high)"
    )
    
    parser.add_argument(
        "--working-dir", "-w",
        default="results/augustus_training",
        help="Working directory for training files (default: results/augustus_training)"
    )
    
    # Training parameters
    parser.add_argument(
        "--min-genes",
        type=int,
        default=20,
        help="Minimum number of training genes required (default: 20)"
    )
    
    parser.add_argument(
        "--flanking-dna",
        type=int,
        default=4000,
        help="Length of flanking DNA sequence (default: 4000)"
    )
    
    parser.add_argument(
        "--optimization-rounds",
        type=int,
        default=1,
        help="Number of optimization rounds (default: 1, use 0 to skip)"
    )
    
    parser.add_argument(
        "--no-optimization",
        action="store_true",
        help="Skip parameter optimization (faster but less accurate)"
    )
    
    parser.add_argument(
        "--cpus",
        type=int,
        default=4,
        help="Number of CPU cores to use (default: 4)"
    )
    
    parser.add_argument(
        "--timeout",
        type=int,
        default=240,
        help="Training timeout in minutes (default: 240)"
    )
    
    # Quality thresholds
    parser.add_argument(
        "--min-identity",
        type=float,
        default=0.95,
        help="Minimum identity threshold for training data (default: 0.95)"
    )
    
    parser.add_argument(
        "--max-frameshifts",
        type=int,
        default=0,
        help="Maximum frameshifts allowed in training data (default: 0)"
    )
    
    parser.add_argument(
        "--max-stop-codons",
        type=int,
        default=0,
        help="Maximum stop codons allowed in training data (default: 0)"
    )
    
    # Processing options
    parser.add_argument(
        "--no-backup",
        action="store_true",
        help="Don't backup existing model (overwrites without backup)"
    )
    
    parser.add_argument(
        "--no-report",
        action="store_true",
        help="Don't generate training report"
    )
    
    parser.add_argument(
        "--no-gene-features",
        action="store_true",
        help="Don't add gene features to training data"
    )
    
    parser.add_argument(
        "--allow-partial-genes",
        action="store_true",
        help="Allow partial genes (missing start/stop codons) in training data"
    )
    
    # Output and logging
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Verbose output"
    )
    
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Debug mode (very verbose)"
    )
    
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Quiet mode (minimal output)"
    )
    
    parser.add_argument(
        "--log-file",
        help="Log file path (default: auto-generated)"
    )
    
    # Augustus configuration
    parser.add_argument(
        "--augustus-scripts",
        default="/home/wangys/data/work/nbs/tools/Augustus/scripts",
        help="Path to Augustus scripts directory"
    )
    
    return parser


def validate_arguments(args):
    """Validate command line arguments."""
    errors = []
    
    # Check genome file
    if not os.path.exists(args.genome):
        errors.append(f"Genome file not found: {args.genome}")
    
    # Check Augustus scripts
    if not os.path.exists(args.augustus_scripts):
        errors.append(f"Augustus scripts directory not found: {args.augustus_scripts}")
    
    # Check autoAugTrain.pl
    auto_aug_train = os.path.join(args.augustus_scripts, "autoAugTrain.pl")
    if not os.path.exists(auto_aug_train):
        errors.append(f"autoAugTrain.pl not found: {auto_aug_train}")
    
    # Check miniprot file if specified
    if args.miniprot_file and not os.path.exists(args.miniprot_file):
        errors.append(f"Miniprot file not found: {args.miniprot_file}")
    
    # Validate species name
    if not args.species.replace('_', '').replace('-', '').isalnum():
        errors.append("Species name should contain only alphanumeric characters, underscores, and hyphens")
    
    # Validate numerical parameters
    if args.min_genes < 1:
        errors.append("Minimum genes must be at least 1")
    
    if args.flanking_dna < 0:
        errors.append("Flanking DNA length cannot be negative")
    
    if args.optimization_rounds < 0:
        errors.append("Optimization rounds cannot be negative")
    
    if args.cpus < 1:
        errors.append("Number of CPUs must be at least 1")
    
    if args.timeout < 1:
        errors.append("Timeout must be at least 1 minute")
    
    if not 0.0 <= args.min_identity <= 1.0:
        errors.append("Minimum identity must be between 0.0 and 1.0")
    
    if args.max_frameshifts < 0:
        errors.append("Maximum frameshifts cannot be negative")
    
    if args.max_stop_codons < 0:
        errors.append("Maximum stop codons cannot be negative")
    
    return errors


def create_training_config(args):
    """Create training configuration from command line arguments."""
    # Handle optimization rounds
    optimization_rounds = 0 if args.no_optimization else args.optimization_rounds
    
    config = MiniprotTrainingConfig(
        # Basic configuration
        species_name=args.species,
        genome_file=args.genome,
        working_dir=args.working_dir,
        miniprot_gff_file=args.miniprot_file or "",
        miniprot_quality_filter=args.quality,
        
        # Augustus training parameters
        augustus_scripts_path=args.augustus_scripts,
        flanking_dna_length=args.flanking_dna,
        optimization_rounds=optimization_rounds,
        cpus=args.cpus,
        timeout_minutes=args.timeout,
        
        # Quality thresholds
        min_training_genes=args.min_genes,
        min_identity_threshold=args.min_identity,
        max_frameshifts=args.max_frameshifts,
        max_stop_codons=args.max_stop_codons,
        
        # Training data processing
        filter_partial_genes=not args.allow_partial_genes,
        add_gene_features=not args.no_gene_features,
        
        # Backup and reporting
        backup_existing_model=not args.no_backup,
        create_training_report=not args.no_report
    )
    
    return config


def setup_logging_from_args(args):
    """Setup logging based on command line arguments."""
    if args.debug:
        level = "DEBUG"
    elif args.verbose:
        level = "INFO"
    elif args.quiet:
        level = "WARNING"
    else:
        level = "INFO"
    
    logger = setup_logging("augustus_training_cli", level=level)
    
    if args.log_file:
        # Add file handler if specified
        file_handler = logging.FileHandler(args.log_file)
        file_handler.setLevel(getattr(logging, level))
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    
    return logger


def print_configuration_summary(config, logger):
    """Print configuration summary."""
    logger.info("Augustus Training Configuration:")
    logger.info("=" * 40)
    logger.info(f"Species name: {config.species_name}")
    logger.info(f"Genome file: {config.genome_file}")
    logger.info(f"Quality filter: {config.miniprot_quality_filter}")
    logger.info(f"Working directory: {config.working_dir}")
    logger.info(f"Min training genes: {config.min_training_genes}")
    logger.info(f"Flanking DNA: {config.flanking_dna_length} bp")
    logger.info(f"Optimization rounds: {config.optimization_rounds}")
    logger.info(f"CPUs: {config.cpus}")
    logger.info(f"Timeout: {config.timeout_minutes} minutes")
    logger.info(f"Min identity: {config.min_identity_threshold}")
    logger.info("=" * 40)


def print_results_summary(result, logger):
    """Print training results summary."""
    logger.info("Training Results:")
    logger.info("=" * 30)
    
    if result.training_success:
        logger.info("🎉 Training Status: SUCCESS")
        logger.info(f"📊 Training time: {result.training_time_minutes:.1f} minutes")
        logger.info(f"📈 Training genes: {result.filtered_training_genes}")
        logger.info(f"🔬 Quality filter: {result.quality_filter_applied}")
        logger.info(f"✅ Validation: {'PASSED' if result.validation_passed else 'FAILED'}")
        logger.info(f"📁 Model directory: {result.model_directory}")
        
        if result.training_gff_file:
            logger.info(f"📄 Training GFF: {result.training_gff_file}")
        
        if result.log_file:
            logger.info(f"📝 Log file: {result.log_file}")
        
        if result.report_file:
            logger.info(f"📋 Report: {result.report_file}")
        
        logger.info("\n🚀 Next Steps:")
        logger.info(f"Use your trained model with:")
        logger.info(f"augustus --species={result.species_name} --gff3=on genome.fa")
        
    else:
        logger.error("❌ Training Status: FAILED")
        logger.error(f"💥 Error: {result.error_message}")
        
        if result.warnings:
            logger.warning("⚠️  Warnings:")
            for warning in result.warnings:
                logger.warning(f"  - {warning}")
    
    logger.info("=" * 30)


def main():
    """Main function."""
    parser = setup_argument_parser()
    args = parser.parse_args()
    
    # Setup logging
    logger = setup_logging_from_args(args)
    
    try:
        logger.info("🧬 Augustus Training with Miniprot Results")
        logger.info("=" * 50)
        
        # Validate arguments
        errors = validate_arguments(args)
        if errors:
            logger.error("Configuration validation failed:")
            for error in errors:
                logger.error(f"  - {error}")
            return 1
        
        # Create training configuration
        config = create_training_config(args)
        print_configuration_summary(config, logger)
        
        # Initialize trainer
        logger.info("🚀 Initializing Augustus trainer...")
        trainer = AugustusMiniprotTrainer(config, logger)
        
        # Run training
        logger.info("🎯 Starting Augustus training workflow...")
        logger.info("This may take several hours depending on configuration...")
        
        result = trainer.train_augustus_model()
        
        # Print results
        print_results_summary(result, logger)
        
        # Return appropriate exit code
        return 0 if result.training_success else 1
        
    except KeyboardInterrupt:
        logger.warning("⚠️  Training interrupted by user")
        return 130
    except Exception as e:
        logger.error(f"💥 Training failed with exception: {e}")
        if args.debug:
            import traceback
            logger.error(traceback.format_exc())
        return 1


if __name__ == "__main__":
    sys.exit(main())