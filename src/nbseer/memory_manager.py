#!/usr/bin/env python3
"""
Memory Management Module for Large Genome Processing

This module provides comprehensive memory management and partitioning strategies
for processing large genomic datasets in the EVM pipeline.

Author: Assistant
Date: 2025-06-24
"""

import os
import sys
import gc
import psutil
import logging
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple, Generator, Any
from pathlib import Path
import subprocess
import json
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
import tempfile

@dataclass
class MemoryStats:
    """Memory usage statistics and thresholds."""
    total_memory: int  # Total system memory in bytes
    available_memory: int  # Available memory in bytes
    used_memory: int  # Used memory in bytes
    memory_percent: float  # Memory usage percentage
    warning_threshold: float = 80.0  # Warning threshold percentage
    danger_threshold: float = 90.0  # Danger threshold percentage
    
    @property
    def is_warning(self) -> bool:
        """Check if memory usage is at warning level."""
        return self.memory_percent >= self.warning_threshold
    
    @property
    def is_danger(self) -> bool:
        """Check if memory usage is at danger level."""
        return self.memory_percent >= self.danger_threshold
    
    @property
    def available_gb(self) -> float:
        """Available memory in GB."""
        return self.available_memory / (1024**3)
    
    @property
    def total_gb(self) -> float:
        """Total memory in GB."""
        return self.total_memory / (1024**3)

@dataclass
class GenomePartition:
    """Information about a genome partition."""
    chromosome: str
    partition_id: str
    start_pos: int
    end_pos: int
    size: int
    feature_count: int = 0
    file_path: Optional[str] = None
    processed: bool = False
    
    @property
    def region_str(self) -> str:
        """Get region string for samtools."""
        return f"{self.chromosome}:{self.start_pos}-{self.end_pos}"
    
    @property
    def size_mb(self) -> float:
        """Partition size in MB."""
        return self.size / (1024**2)

@dataclass
class PartitioningConfig:
    """Configuration for genome partitioning."""
    target_partition_size: int = 10 * 1024 * 1024  # 10MB default
    min_partition_size: int = 1 * 1024 * 1024  # 1MB minimum
    max_partition_size: int = 100 * 1024 * 1024  # 100MB maximum
    overlap_size: int = 10000  # 10kb overlap between partitions
    memory_safety_factor: float = 0.7  # Use 70% of available memory
    max_concurrent_partitions: int = 4  # Maximum parallel partitions
    feature_density_threshold: int = 100  # Features per MB for dense regions

class GenomeMemoryManager:
    """
    Comprehensive memory management for large genome processing.
    
    Provides intelligent partitioning, memory monitoring, and resource management
    for processing large genomic datasets efficiently.
    """
    
    def __init__(self, 
                 genome_file: str,
                 output_dir: str,
                 config: Optional[PartitioningConfig] = None,
                 logger: Optional[logging.Logger] = None):
        """
        Initialize the GenomeMemoryManager.
        
        Args:
            genome_file: Path to the genome FASTA file
            output_dir: Output directory for partitions and temporary files
            config: Partitioning configuration
            logger: Logger instance
        """
        self.genome_file = Path(genome_file)
        self.output_dir = Path(output_dir)
        self.config = config or PartitioningConfig()
        self.logger = logger or self._setup_logger()
        
        # Initialize attributes
        self.partitions: List[GenomePartition] = []
        self.chromosome_info: Dict[str, Dict] = {}
        self.temp_dir: Optional[Path] = None
        
        # Validate inputs
        self._validate_inputs()
        
        # Setup output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize genome index
        self._initialize_genome_index()
    
    def _setup_logger(self) -> logging.Logger:
        """Setup logger for memory manager."""
        logger = logging.getLogger(f"{__name__}.GenomeMemoryManager")
        if not logger.handlers:
            handler = logging.StreamHandler()
            formatter = logging.Formatter(
                '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
            )
            handler.setFormatter(formatter)
            logger.addHandler(handler)
            logger.setLevel(logging.INFO)
        return logger
    
    def _validate_inputs(self) -> None:
        """Validate input files and directories."""
        if not self.genome_file.exists():
            raise FileNotFoundError(f"Genome file not found: {self.genome_file}")
        
        if not self.genome_file.is_file():
            raise ValueError(f"Genome file is not a regular file: {self.genome_file}")
        
        # Check if genome file is readable
        try:
            with open(self.genome_file, 'r') as f:
                f.read(1)
        except Exception as e:
            raise ValueError(f"Cannot read genome file: {e}")
    
    def _initialize_genome_index(self) -> None:
        """Initialize samtools faidx index for the genome."""
        index_file = f"{self.genome_file}.fai"
        
        if not os.path.exists(index_file):
            self.logger.info(f"Creating samtools index for {self.genome_file}")
            try:
                subprocess.run(
                    ["samtools", "faidx", str(self.genome_file)],
                    check=True,
                    capture_output=True,
                    text=True
                )
            except subprocess.CalledProcessError as e:
                raise RuntimeError(f"Failed to create genome index: {e.stderr}")
        
        # Load chromosome information
        self._load_chromosome_info(index_file)
    
    def _load_chromosome_info(self, index_file: str) -> None:
        """Load chromosome information from faidx index."""
        try:
            with open(index_file, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:
                        chrom = parts[0]
                        length = int(parts[1])
                        self.chromosome_info[chrom] = {
                            'length': length,
                            'size_mb': length / (1024**2)
                        }
        except Exception as e:
            raise RuntimeError(f"Failed to load chromosome info: {e}")
        
        self.logger.info(f"Loaded {len(self.chromosome_info)} chromosomes")
        total_size = sum(info['length'] for info in self.chromosome_info.values())
        self.logger.info(f"Total genome size: {total_size / (1024**3):.2f} GB")
    
    def get_memory_stats(self) -> MemoryStats:
        """Get current memory statistics."""
        memory = psutil.virtual_memory()
        return MemoryStats(
            total_memory=memory.total,
            available_memory=memory.available,
            used_memory=memory.used,
            memory_percent=memory.percent
        )
    
    def monitor_memory(self, interval: float = 1.0) -> Generator[MemoryStats, None, None]:
        """
        Monitor memory usage continuously.
        
        Args:
            interval: Monitoring interval in seconds
            
        Yields:
            MemoryStats objects with current memory information
        """
        while True:
            stats = self.get_memory_stats()
            yield stats
            time.sleep(interval)
    
    def estimate_feature_density(self, gff_files: List[str]) -> Dict[str, float]:
        """
        Estimate feature density per chromosome from GFF files.
        
        Args:
            gff_files: List of GFF3 files to analyze
            
        Returns:
            Dictionary mapping chromosome to features per MB
        """
        feature_counts = {}
        
        for gff_file in gff_files:
            if not os.path.exists(gff_file):
                continue
                
            try:
                with open(gff_file, 'r') as f:
                    for line in f:
                        if line.startswith('#'):
                            continue
                        
                        parts = line.strip().split('\t')
                        if len(parts) >= 9:
                            chrom = parts[0]
                            feature_counts[chrom] = feature_counts.get(chrom, 0) + 1
            except Exception as e:
                self.logger.warning(f"Error reading {gff_file}: {e}")
        
        # Calculate density (features per MB)
        density = {}
        for chrom, count in feature_counts.items():
            if chrom in self.chromosome_info:
                size_mb = self.chromosome_info[chrom]['size_mb']
                density[chrom] = count / size_mb if size_mb > 0 else 0
        
        return density
    
    def calculate_optimal_partition_size(self, 
                                       feature_density: Optional[Dict[str, float]] = None) -> int:
        """
        Calculate optimal partition size based on available memory and feature density.
        
        Args:
            feature_density: Feature density per chromosome (features per MB)
            
        Returns:
            Optimal partition size in bytes
        """
        memory_stats = self.get_memory_stats()
        available_memory = memory_stats.available_memory * self.config.memory_safety_factor
        
        # Base partition size on available memory
        # Aim for partitions that use ~5% of available memory each
        base_partition_size = int(available_memory * 0.05)
        
        # Adjust based on feature density if provided
        if feature_density:
            avg_density = sum(feature_density.values()) / len(feature_density)
            if avg_density > self.config.feature_density_threshold:
                # Reduce partition size for dense regions
                base_partition_size = int(base_partition_size * 0.5)
        
        # Apply constraints
        partition_size = max(self.config.min_partition_size, 
                           min(base_partition_size, self.config.max_partition_size))
        
        self.logger.info(f"Calculated optimal partition size: {partition_size / (1024**2):.1f} MB")
        return partition_size
    
    def create_partitions(self, 
                         feature_density: Optional[Dict[str, float]] = None) -> List[GenomePartition]:
        """
        Create genome partitions based on chromosome sizes and feature density.
        
        Args:
            feature_density: Feature density per chromosome
            
        Returns:
            List of GenomePartition objects
        """
        self.partitions = []
        partition_size = self.calculate_optimal_partition_size(feature_density)
        
        partition_id = 0
        for chrom, info in self.chromosome_info.items():
            chrom_length = info['length']
            
            # Adjust partition size for this chromosome based on density
            current_partition_size = partition_size
            if feature_density and chrom in feature_density:
                density = feature_density[chrom]
                if density > self.config.feature_density_threshold:
                    current_partition_size = int(partition_size * 0.5)
            
            # Create partitions for this chromosome
            start = 1
            while start <= chrom_length:
                end = min(start + current_partition_size - 1, chrom_length)
                
                partition = GenomePartition(
                    chromosome=chrom,
                    partition_id=f"{chrom}_{start}-{end}",
                    start_pos=start,
                    end_pos=end,
                    size=end - start + 1
                )
                
                self.partitions.append(partition)
                partition_id += 1
                
                # Move to next partition with overlap
                start = end - self.config.overlap_size + 1
                if start > chrom_length:
                    break
        
        self.logger.info(f"Created {len(self.partitions)} partitions")
        self._log_partition_summary()
        return self.partitions
    
    def _log_partition_summary(self) -> None:
        """Log summary of created partitions."""
        if not self.partitions:
            return
        
        total_size = sum(p.size for p in self.partitions)
        avg_size = total_size / len(self.partitions)
        
        self.logger.info(f"Partition summary:")
        self.logger.info(f"  Total partitions: {len(self.partitions)}")
        self.logger.info(f"  Average size: {avg_size / (1024**2):.1f} MB")
        self.logger.info(f"  Total coverage: {total_size / (1024**3):.2f} GB")
        
        # Log per-chromosome statistics
        chrom_stats = {}
        for partition in self.partitions:
            chrom = partition.chromosome
            if chrom not in chrom_stats:
                chrom_stats[chrom] = {'count': 0, 'total_size': 0}
            chrom_stats[chrom]['count'] += 1
            chrom_stats[chrom]['total_size'] += partition.size
        
        for chrom, stats in chrom_stats.items():
            self.logger.info(f"  {chrom}: {stats['count']} partitions, "
                           f"{stats['total_size'] / (1024**2):.1f} MB")
    
    def extract_partition_sequence(self, partition: GenomePartition) -> str:
        """
        Extract sequence for a specific partition using samtools.
        
        Args:
            partition: GenomePartition object
            
        Returns:
            Path to the extracted sequence file
        """
        if not self.temp_dir:
            self.temp_dir = Path(tempfile.mkdtemp(prefix="genome_partitions_"))
        
        output_file = self.temp_dir / f"{partition.partition_id}.fasta"
        
        try:
            cmd = [
                "samtools", "faidx", 
                str(self.genome_file),
                partition.region_str
            ]
            
            with open(output_file, 'w') as f:
                subprocess.run(cmd, check=True, stdout=f, stderr=subprocess.PIPE)
            
            partition.file_path = str(output_file)
            return str(output_file)
            
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Failed to extract partition {partition.partition_id}: {e}")
    
    def estimate_processing_memory(self, 
                                 num_partitions: int,
                                 avg_features_per_partition: int = 1000) -> Dict[str, Any]:
        """
        Estimate memory requirements for processing partitions.
        
        Args:
            num_partitions: Number of partitions to process
            avg_features_per_partition: Average features per partition
            
        Returns:
            Memory estimation dictionary
        """
        # Base memory estimates (rough approximations)
        base_memory_per_partition = 50 * 1024 * 1024  # 50MB base
        memory_per_feature = 1024  # 1KB per feature
        
        estimated_memory_per_partition = (
            base_memory_per_partition + 
            (avg_features_per_partition * memory_per_feature)
        )
        
        total_estimated_memory = estimated_memory_per_partition * num_partitions
        
        memory_stats = self.get_memory_stats()
        available_memory = memory_stats.available_memory * self.config.memory_safety_factor
        
        # Calculate safe number of concurrent partitions
        safe_concurrent = max(1, int(available_memory / estimated_memory_per_partition))
        safe_concurrent = min(safe_concurrent, self.config.max_concurrent_partitions)
        
        return {
            'estimated_memory_per_partition': estimated_memory_per_partition,
            'total_estimated_memory': total_estimated_memory,
            'available_memory': available_memory,
            'safe_concurrent_partitions': safe_concurrent,
            'memory_sufficient': total_estimated_memory <= available_memory,
            'estimated_memory_per_partition_mb': estimated_memory_per_partition / (1024**2),
            'total_estimated_memory_gb': total_estimated_memory / (1024**3),
            'available_memory_gb': available_memory / (1024**3)
        }
    
    def process_partitions_with_memory_management(self,
                                                processing_func,
                                                *args,
                                                **kwargs) -> List[Any]:
        """
        Process partitions with intelligent memory management.
        
        Args:
            processing_func: Function to process each partition
            *args, **kwargs: Arguments for processing function
            
        Returns:
            List of processing results
        """
        if not self.partitions:
            raise ValueError("No partitions created. Call create_partitions() first.")
        
        # Estimate memory requirements
        avg_features = kwargs.get('avg_features_per_partition', 1000)
        memory_estimate = self.estimate_processing_memory(
            len(self.partitions), avg_features
        )
        
        self.logger.info(f"Memory estimation: {memory_estimate}")
        
        if not memory_estimate['memory_sufficient']:
            self.logger.warning("Estimated memory requirements exceed available memory")
            self.logger.warning("Consider reducing partition size or feature density")
        
        # Process partitions with concurrency control
        max_workers = memory_estimate['safe_concurrent_partitions']
        results = []
        
        self.logger.info(f"Processing {len(self.partitions)} partitions with "
                        f"{max_workers} concurrent workers")
        
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            # Submit all tasks
            future_to_partition = {
                executor.submit(processing_func, partition, *args, **kwargs): partition
                for partition in self.partitions
            }
            
            # Process completed tasks
            for future in as_completed(future_to_partition):
                partition = future_to_partition[future]
                try:
                    result = future.result()
                    results.append(result)
                    partition.processed = True
                    
                    # Monitor memory during processing
                    memory_stats = self.get_memory_stats()
                    if memory_stats.is_warning:
                        self.logger.warning(f"Memory usage at {memory_stats.memory_percent:.1f}%")
                        if memory_stats.is_danger:
                            self.logger.error(f"Memory usage critical: {memory_stats.memory_percent:.1f}%")
                            # Force garbage collection
                            gc.collect()
                    
                except Exception as e:
                    self.logger.error(f"Error processing partition {partition.partition_id}: {e}")
                    raise
        
        self.logger.info(f"Completed processing {len(results)} partitions")
        return results
    
    def save_partition_info(self, output_file: str) -> None:
        """
        Save partition information to JSON file.
        
        Args:
            output_file: Path to output JSON file
        """
        partition_data = {
            'config': {
                'target_partition_size': self.config.target_partition_size,
                'min_partition_size': self.config.min_partition_size,
                'max_partition_size': self.config.max_partition_size,
                'overlap_size': self.config.overlap_size,
                'memory_safety_factor': self.config.memory_safety_factor,
                'max_concurrent_partitions': self.config.max_concurrent_partitions
            },
            'genome_info': {
                'genome_file': str(self.genome_file),
                'chromosomes': self.chromosome_info
            },
            'partitions': [
                {
                    'chromosome': p.chromosome,
                    'partition_id': p.partition_id,
                    'start_pos': p.start_pos,
                    'end_pos': p.end_pos,
                    'size': p.size,
                    'size_mb': p.size_mb,
                    'feature_count': p.feature_count,
                    'processed': p.processed
                }
                for p in self.partitions
            ],
            'summary': {
                'total_partitions': len(self.partitions),
                'total_size': sum(p.size for p in self.partitions),
                'average_size': sum(p.size for p in self.partitions) / len(self.partitions) if self.partitions else 0,
                'processed_count': sum(1 for p in self.partitions if p.processed)
            }
        }
        
        with open(output_file, 'w') as f:
            json.dump(partition_data, f, indent=2)
        
        self.logger.info(f"Saved partition information to {output_file}")
    
    def cleanup_temp_files(self) -> None:
        """Clean up temporary files and directories."""
        if self.temp_dir and self.temp_dir.exists():
            import shutil
            try:
                shutil.rmtree(self.temp_dir)
                self.logger.info(f"Cleaned up temporary directory: {self.temp_dir}")
            except Exception as e:
                self.logger.warning(f"Failed to clean up temporary directory: {e}")
    
    def __enter__(self):
        """Context manager entry."""
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit with cleanup."""
        self.cleanup_temp_files()


def main():
    """Example usage of GenomeMemoryManager."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Genome Memory Manager Demo")
    parser.add_argument("--genome", required=True, help="Path to genome FASTA file")
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument("--gff-files", nargs="*", help="GFF files for feature density estimation")
    
    args = parser.parse_args()
    
    # Setup logging
    logging.basicConfig(level=logging.INFO)
    
    try:
        with GenomeMemoryManager(args.genome, args.output) as manager:
            # Get memory stats
            memory_stats = manager.get_memory_stats()
            print(f"Memory Stats: {memory_stats.memory_percent:.1f}% used, "
                  f"{memory_stats.available_gb:.1f} GB available")
            
            # Estimate feature density if GFF files provided
            feature_density = None
            if args.gff_files:
                feature_density = manager.estimate_feature_density(args.gff_files)
                print(f"Feature density: {feature_density}")
            
            # Create partitions
            partitions = manager.create_partitions(feature_density)
            
            # Estimate processing memory
            memory_estimate = manager.estimate_processing_memory(len(partitions))
            print(f"Memory estimate: {memory_estimate}")
            
            # Save partition info
            manager.save_partition_info(os.path.join(args.output, "partition_info.json"))
            
    except Exception as e:
        print(f"Error: {e}")
        return 1
    
    return 0


if __name__ == "__main__":
    sys.exit(main()) 