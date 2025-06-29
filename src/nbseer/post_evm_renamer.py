#!/usr/bin/env python3
"""
Post-EVM Gene ID Renaming Tool
EVM后基因ID重命名工具

用户指定前缀的基因ID重命名功能，独立于EVM流程运行
User-defined prefix gene ID renaming functionality, runs independently from EVM pipeline
"""

import os
import sys
import argparse
import logging
from pathlib import Path
from typing import Dict, Any
import json

# Add project root to Python path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

from src.nbseer.gene_id_renamer import GeneIDRenamer
from src.nbseer.species_config import GeneNamingConfig


class PostEVMRenamer:
    """Post-EVM gene ID renaming tool."""
    
    def __init__(self, prefix: str = None, output_dir: str = None):
        """
        Initialize the post-EVM renamer.
        
        Args:
            prefix: User-defined prefix for gene IDs (default: "NBS")
            output_dir: Output directory for renamed files
        """
        self.prefix = prefix or "NBS"
        self.output_dir = Path(output_dir) if output_dir else Path("./renamed_output")
        self.output_dir.mkdir(exist_ok=True)
        
        # Setup logger
        self.logger = self._setup_logger()
        
    def _setup_logger(self) -> logging.Logger:
        """Setup logger for the renamer."""
        logger = logging.getLogger("PostEVMRenamer")
        logger.setLevel(logging.INFO)
        
        if not logger.handlers:
            # Console handler
            console_handler = logging.StreamHandler()
            console_formatter = logging.Formatter(
                '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
            )
            console_handler.setFormatter(console_formatter)
            logger.addHandler(console_handler)
            
            # File handler
            log_file = self.output_dir / "gene_renaming.log"
            file_handler = logging.FileHandler(log_file)
            file_handler.setFormatter(console_formatter)
            logger.addHandler(file_handler)
            
        return logger
    
    def validate_input_file(self, input_file: str) -> bool:
        """
        Validate input GFF3 file.
        
        Args:
            input_file: Path to input GFF3 file
            
        Returns:
            True if valid, False otherwise
        """
        if not os.path.exists(input_file):
            self.logger.error(f"Input file does not exist: {input_file}")
            return False
        
        if not input_file.lower().endswith(('.gff', '.gff3')):
            self.logger.warning(f"Input file may not be GFF3 format: {input_file}")
        
        # Check if file contains EVM IDs
        evm_ids_found = False
        try:
            with open(input_file, 'r') as f:
                for i, line in enumerate(f):
                    if i > 100:  # Check first 100 lines
                        break
                    if 'evm.model' in line or 'evm.TU' in line:
                        evm_ids_found = True
                        break
        except Exception as e:
            self.logger.error(f"Error reading input file: {e}")
            return False
        
        if not evm_ids_found:
            self.logger.warning("No EVM IDs found in input file. Make sure this is EVM output.")
        
        return True
    
    def rename_evm_output(self, input_file: str, output_prefix: str = None) -> Dict[str, Any]:
        """
        Rename EVM output gene IDs.
        
        Args:
            input_file: Path to EVM output GFF3 file
            output_prefix: Optional prefix override for this specific run
            
        Returns:
            Dictionary with renaming results and statistics
        """
        self.logger.info(f"Starting gene ID renaming for: {input_file}")
        
        # Validate input
        if not self.validate_input_file(input_file):
            return {"success": False, "error": "Input validation failed"}
        
        # Use provided prefix or default
        prefix = output_prefix or self.prefix
        
        # Validate prefix
        config = GeneNamingConfig(prefix)
        if not config.validate_prefix(prefix):
            self.logger.error(f"Invalid prefix: {prefix}")
            return {"success": False, "error": f"Invalid prefix: {prefix}"}
        
        # Generate output file paths
        input_path = Path(input_file)
        output_file = self.output_dir / f"{input_path.stem}_{prefix}_renamed.gff3"
        mapping_file = self.output_dir / f"{input_path.stem}_{prefix}_mapping.txt"
        stats_file = self.output_dir / f"{input_path.stem}_{prefix}_stats.json"
        
        try:
            # Create renamer
            renamer = GeneIDRenamer(prefix)
            
            # Get chromosome summary before renaming
            chr_summary = renamer.get_chromosome_summary(input_file)
            self.logger.info(f"Found chromosomes: {list(chr_summary.keys())}")
            
            # Perform renaming
            id_mapping = renamer.rename_gff3_file(str(input_file), str(output_file))
            
            # Save mapping file
            renamer.save_mapping_file(str(mapping_file))
            
            # Get statistics
            mapping_summary = renamer.get_mapping_summary()
            
            # Create detailed statistics
            stats = {
                "input_file": str(input_file),
                "output_file": str(output_file),
                "mapping_file": str(mapping_file),
                "prefix_used": prefix,
                "total_ids_renamed": len(id_mapping),
                "chromosome_summary": chr_summary,
                "mapping_summary": mapping_summary,
                "first_few_mappings": dict(list(id_mapping.items())[:10]),  # First 10 mappings as examples
            }
            
            # Save statistics
            with open(stats_file, 'w') as f:
                json.dump(stats, f, indent=2)
            
            self.logger.info(f"Gene ID renaming completed successfully:")
            self.logger.info(f"  - Prefix used: {prefix}")
            self.logger.info(f"  - Total IDs renamed: {len(id_mapping)}")
            self.logger.info(f"  - Output file: {output_file}")
            self.logger.info(f"  - Mapping file: {mapping_file}")
            self.logger.info(f"  - Statistics file: {stats_file}")
            
            return {
                "success": True,
                "output_file": str(output_file),
                "mapping_file": str(mapping_file),
                "stats_file": str(stats_file),
                "statistics": stats
            }
            
        except Exception as e:
            self.logger.error(f"Gene ID renaming failed: {e}")
            return {"success": False, "error": str(e)}
    
    def batch_rename(self, input_files: list, prefix: str = None) -> Dict[str, Any]:
        """
        Rename multiple EVM output files.
        
        Args:
            input_files: List of input GFF3 files
            prefix: Optional prefix override
            
        Returns:
            Dictionary with batch renaming results
        """
        self.logger.info(f"Starting batch renaming for {len(input_files)} files")
        
        results = {}
        total_success = 0
        total_failed = 0
        
        for input_file in input_files:
            self.logger.info(f"Processing: {input_file}")
            result = self.rename_evm_output(input_file, prefix)
            
            results[input_file] = result
            if result["success"]:
                total_success += 1
            else:
                total_failed += 1
        
        # Create batch summary
        batch_summary = {
            "total_files": len(input_files),
            "successful": total_success,
            "failed": total_failed,
            "prefix_used": prefix or self.prefix,
            "output_directory": str(self.output_dir),
            "results": results
        }
        
        # Save batch summary
        batch_summary_file = self.output_dir / f"batch_rename_summary_{prefix or self.prefix}.json"
        with open(batch_summary_file, 'w') as f:
            json.dump(batch_summary, f, indent=2)
        
        self.logger.info(f"Batch renaming completed:")
        self.logger.info(f"  - Total files: {len(input_files)}")
        self.logger.info(f"  - Successful: {total_success}")
        self.logger.info(f"  - Failed: {total_failed}")
        self.logger.info(f"  - Summary file: {batch_summary_file}")
        
        return batch_summary


def main():
    """Main function for command-line usage."""
    parser = argparse.ArgumentParser(
        description="Post-EVM Gene ID Renaming Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Rename with default prefix "NBS"
  python post_evm_renamer.py -i evm_output.gff3
  
  # Rename with custom prefix "CaCM334"
  python post_evm_renamer.py -i evm_output.gff3 -p CaCM334
  
  # Batch rename multiple files
  python post_evm_renamer.py -i file1.gff3 file2.gff3 file3.gff3 -p MyPrefix
  
  # Specify output directory
  python post_evm_renamer.py -i evm_output.gff3 -p CaCM334 -o /path/to/output
        """
    )
    
    parser.add_argument(
        '-i', '--input',
        nargs='+',
        required=True,
        help='Input EVM output GFF3 file(s)'
    )
    
    parser.add_argument(
        '-p', '--prefix',
        default='NBS',
        help='Gene ID prefix (default: NBS). Examples: NBS, CaCM334, AtCol0'
    )
    
    parser.add_argument(
        '-o', '--output-dir',
        default='./renamed_output',
        help='Output directory for renamed files (default: ./renamed_output)'
    )
    
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Enable verbose logging'
    )
    
    args = parser.parse_args()
    
    # Setup logging level
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Create renamer
    renamer = PostEVMRenamer(prefix=args.prefix, output_dir=args.output_dir)
    
    # Process files
    if len(args.input) == 1:
        # Single file
        result = renamer.rename_evm_output(args.input[0])
        if not result["success"]:
            print(f"Error: {result['error']}")
            sys.exit(1)
    else:
        # Batch processing
        result = renamer.batch_rename(args.input)
        if result["failed"] > 0:
            print(f"Warning: {result['failed']} files failed to process")
    
    print("Gene ID renaming completed successfully!")


if __name__ == "__main__":
    main()