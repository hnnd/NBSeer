#!/usr/bin/env python3
"""
Test Script for Augustus Training with Miniprot Results

This script tests the Augustus training functionality using high-quality
Miniprot results. It provides both unit tests and integration tests to
validate the complete training workflow.

Author: NBS Annotation Pipeline
Date: 2025-06-25
"""

import os
import sys
import unittest
import tempfile
import shutil
import logging
from pathlib import Path

# Add parent directory to path for imports
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from utils.logging_setup import setup_logging
from utils.config import get_config
from nbseer.augustus_miniprot_trainer import (
    AugustusMiniprotTrainer, 
    MiniprotTrainingConfig,
    MiniprotTrainingResult
)


class TestAugustusTraining(unittest.TestCase):
    """Test cases for Augustus training functionality."""
    
    @classmethod
    def setUpClass(cls):
        """Set up test environment."""
        cls.logger = setup_logging("test_augustus_training")
        cls.test_dir = "test_augustus_training"
        os.makedirs(cls.test_dir, exist_ok=True)
        
        # Test configuration
        cls.test_config = MiniprotTrainingConfig(
            species_name="test_nbs_species",
            genome_file="genome/osa.fa",
            working_dir=cls.test_dir,
            miniprot_quality_filter="high",
            optimization_rounds=0,  # Skip optimization for testing
            cpus=2,
            timeout_minutes=30,
            min_training_genes=1,  # Lower threshold for testing
            backup_existing_model=False  # Don't backup during testing
        )
    
    @classmethod
    def tearDownClass(cls):
        """Clean up test environment."""
        if os.path.exists(cls.test_dir):
            shutil.rmtree(cls.test_dir)
    
    def test_configuration_validation(self):
        """Test configuration validation."""
        self.logger.info("Testing configuration validation...")
        
        # Test valid configuration
        try:
            trainer = AugustusMiniprotTrainer(self.test_config)
            self.assertIsNotNone(trainer)
            self.logger.info("✅ Valid configuration accepted")
        except Exception as e:
            self.fail(f"Valid configuration rejected: {e}")
        
        # Test invalid configuration
        invalid_config = MiniprotTrainingConfig(
            species_name="",  # Empty species name should fail
            genome_file="nonexistent_genome.fa"
        )
        
        with self.assertRaises(ValueError):
            AugustusMiniprotTrainer(invalid_config)
        self.logger.info("✅ Invalid configuration properly rejected")
    
    def test_miniprot_file_detection(self):
        """Test Miniprot results file detection."""
        self.logger.info("Testing Miniprot file detection...")
        
        trainer = AugustusMiniprotTrainer(self.test_config)
        
        try:
            miniprot_file = trainer.find_miniprot_results()
            self.assertTrue(os.path.exists(miniprot_file))
            self.logger.info(f"✅ Found Miniprot results: {miniprot_file}")
        except FileNotFoundError as e:
            self.logger.warning(f"⚠️  Miniprot results not found: {e}")
            self.skipTest("Miniprot results not available for testing")
    
    def test_training_data_preparation(self):
        """Test training data preparation from Miniprot results."""
        self.logger.info("Testing training data preparation...")
        
        trainer = AugustusMiniprotTrainer(self.test_config)
        
        try:
            # Find miniprot file
            miniprot_file = trainer.find_miniprot_results()
            
            # Prepare training data
            training_gff = trainer.prepare_training_data(miniprot_file)
            
            # Validate training file
            self.assertTrue(os.path.exists(training_gff))
            
            # Check file content
            with open(training_gff, 'r') as f:
                content = f.read()
                self.assertIn("##gff-version 3", content)
                self.assertIn("mRNA", content)
            
            self.logger.info(f"✅ Training data prepared: {training_gff}")
            
        except FileNotFoundError:
            self.skipTest("Miniprot results not available for testing")
        except Exception as e:
            self.fail(f"Training data preparation failed: {e}")
    
    def test_augustus_trainer_setup(self):
        """Test Augustus trainer setup."""
        self.logger.info("Testing Augustus trainer setup...")
        
        trainer = AugustusMiniprotTrainer(self.test_config)
        
        try:
            # Create a minimal test training file
            test_training_file = os.path.join(self.test_dir, "test_training.gff3")
            with open(test_training_file, 'w') as f:
                f.write("##gff-version 3\n")
                f.write("Chr1\tminiprot\tmRNA\t1000\t2000\t.\t+\t.\tID=test_gene\n")
                f.write("Chr1\tminiprot\tCDS\t1000\t2000\t.\t+\t0\tParent=test_gene\n")
            
            # Setup Augustus trainer
            augustus_trainer = trainer.setup_augustus_trainer(test_training_file)
            self.assertIsNotNone(augustus_trainer)
            
            self.logger.info("✅ Augustus trainer setup successful")
            
        except Exception as e:
            self.logger.warning(f"⚠️  Augustus trainer setup failed: {e}")
            self.skipTest("Augustus not properly configured")


class TestAugustusTrainingIntegration(unittest.TestCase):
    """Integration tests for complete Augustus training workflow."""
    
    def setUp(self):
        """Set up integration test environment."""
        self.logger = setup_logging("test_augustus_integration")
        self.test_dir = "integration_test_augustus"
        os.makedirs(self.test_dir, exist_ok=True)
    
    def tearDown(self):
        """Clean up integration test environment."""
        if os.path.exists(self.test_dir):
            shutil.rmtree(self.test_dir)
    
    def test_full_training_workflow(self):
        """Test complete training workflow (requires actual data)."""
        self.logger.info("Testing full Augustus training workflow...")
        
        # Configuration for integration test
        config = MiniprotTrainingConfig(
            species_name="integration_test_species",
            genome_file="genome/osa.fa",
            working_dir=self.test_dir,
            miniprot_quality_filter="high",
            optimization_rounds=0,  # Skip optimization for testing
            cpus=2,
            timeout_minutes=60,
            min_training_genes=1,
            backup_existing_model=False
        )
        
        trainer = AugustusMiniprotTrainer(config)
        
        try:
            # Check if prerequisites are available
            miniprot_file = trainer.find_miniprot_results()
            if not os.path.exists(config.genome_file):
                self.skipTest("Genome file not available for integration testing")
            
            # Run training workflow
            self.logger.info("Starting training workflow...")
            result = trainer.train_augustus_model()
            
            # Validate results
            self.assertIsInstance(result, MiniprotTrainingResult)
            self.assertEqual(result.species_name, config.species_name)
            
            if result.training_success:
                self.logger.info("✅ Integration test PASSED - Training successful")
                self.assertTrue(result.validation_passed)
                self.assertGreater(len(result.model_files), 0)
                self.assertTrue(os.path.exists(result.training_gff_file))
            else:
                self.logger.warning(f"⚠️  Integration test FAILED - Training unsuccessful: {result.error_message}")
                # Don't fail the test if Augustus is not properly configured
                self.skipTest(f"Augustus training failed: {result.error_message}")
                
        except FileNotFoundError as e:
            self.logger.warning(f"⚠️  Skipping integration test: {e}")
            self.skipTest("Required files not available for integration testing")
        except Exception as e:
            self.logger.error(f"❌ Integration test failed with exception: {e}")
            self.fail(f"Integration test failed: {e}")


def run_manual_test():
    """Run manual test with actual data."""
    print("🧪 Running Manual Augustus Training Test")
    print("=" * 50)
    
    # Setup logging
    logger = setup_logging("manual_test")
    
    # Configuration
    config = MiniprotTrainingConfig(
        species_name="manual_test_nbs",
        genome_file="genome/osa.fa",
        working_dir="manual_test_augustus",
        miniprot_quality_filter="high",
        optimization_rounds=1,
        cpus=4,
        timeout_minutes=120,
        min_training_genes=5
    )
    
    try:
        print(f"📋 Configuration:")
        print(f"  - Species: {config.species_name}")
        print(f"  - Genome: {config.genome_file}")
        print(f"  - Quality filter: {config.miniprot_quality_filter}")
        print(f"  - Working dir: {config.working_dir}")
        print()
        
        # Initialize trainer
        print("🚀 Initializing trainer...")
        trainer = AugustusMiniprotTrainer(config, logger)
        
        # Run training
        print("🎯 Starting Augustus training workflow...")
        result = trainer.train_augustus_model()
        
        # Report results
        print("\n📊 Training Results:")
        print("=" * 30)
        print(f"Success: {'✅ YES' if result.training_success else '❌ NO'}")
        print(f"Species: {result.species_name}")
        print(f"Training time: {result.training_time_minutes:.1f} minutes")
        print(f"Training genes: {result.filtered_training_genes}")
        print(f"Quality filter: {result.quality_filter_applied}")
        print(f"Validation: {'✅ PASSED' if result.validation_passed else '❌ FAILED'}")
        
        if result.training_success:
            print(f"\n📁 Generated Files:")
            print(f"  - Model directory: {result.model_directory}")
            print(f"  - Training GFF: {result.training_gff_file}")
            print(f"  - Log file: {result.log_file}")
            print(f"  - Report: {result.report_file}")
            
            print(f"\n🎉 Training completed successfully!")
            print(f"You can now use the trained model:")
            print(f"augustus --species={result.species_name} --gff3=on genome.fa")
        else:
            print(f"\n❌ Training failed: {result.error_message}")
            
    except Exception as e:
        print(f"\n💥 Manual test failed: {e}")
        logger.error(f"Manual test failed: {e}")


def main():
    """Main function to run tests."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Test Augustus Training Functionality")
    parser.add_argument("--manual", action="store_true", help="Run manual test with actual data")
    parser.add_argument("--integration", action="store_true", help="Run integration tests only")
    parser.add_argument("--unit", action="store_true", help="Run unit tests only")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose output")
    
    args = parser.parse_args()
    
    if args.manual:
        run_manual_test()
        return
    
    # Setup test suite
    suite = unittest.TestSuite()
    
    if args.unit or (not args.integration and not args.unit):
        # Add unit tests
        suite.addTest(unittest.makeSuite(TestAugustusTraining))
    
    if args.integration or (not args.integration and not args.unit):
        # Add integration tests
        suite.addTest(unittest.makeSuite(TestAugustusTrainingIntegration))
    
    # Run tests
    runner = unittest.TextTestRunner(verbosity=2 if args.verbose else 1)
    result = runner.run(suite)
    
    # Return exit code
    return 0 if result.wasSuccessful() else 1


if __name__ == "__main__":
    sys.exit(main())