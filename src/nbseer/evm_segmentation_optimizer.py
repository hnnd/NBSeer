#!/usr/bin/env python3
"""
EVM Segmentation Optimizer

This module provides tools for optimizing EVidenceModeler segmentation parameters
and overlap handling strategies for efficient processing of large genomes.
"""

import os
import logging
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass, field
from datetime import datetime
import subprocess
import re

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@dataclass
class GenomeStatistics:
    """Statistics about the genome for segmentation optimization"""
    total_length: int = 0
    num_chromosomes: int = 0
    chromosome_lengths: Dict[str, int] = field(default_factory=dict)
    longest_gene_length: int = 0
    gene_density: float = 0.0  # genes per kb
    gc_content: float = 0.0
    repeat_content: float = 0.0
    
    def get_average_chromosome_length(self) -> float:
        """Get average chromosome length"""
        if not self.chromosome_lengths:
            return 0.0
        return sum(self.chromosome_lengths.values()) / len(self.chromosome_lengths)
    
    def get_chromosome_size_distribution(self) -> Dict[str, int]:
        """Get chromosome size distribution"""
        if not self.chromosome_lengths:
            return {}
        
        sorted_lengths = sorted(self.chromosome_lengths.values(), reverse=True)
        return {
            "largest": sorted_lengths[0] if sorted_lengths else 0,
            "smallest": sorted_lengths[-1] if sorted_lengths else 0,
            "median": sorted_lengths[len(sorted_lengths)//2] if sorted_lengths else 0
        }


@dataclass
class SegmentationStrategy:
    """Configuration for genome segmentation strategy"""
    name: str
    segment_size: int
    overlap_size: int
    min_intergenic_length: int = 1000
    max_segments_per_chromosome: int = 100
    adaptive_sizing: bool = False
    
    def get_overlap_percentage(self) -> float:
        """Calculate overlap as percentage of segment size"""
        return (self.overlap_size / self.segment_size) * 100
    
    def validate(self) -> bool:
        """Validate segmentation parameters"""
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
class SegmentationResult:
    """Results from segmentation analysis"""
    strategy: SegmentationStrategy
    estimated_segments: int
    estimated_runtime_hours: float
    memory_usage_gb: float
    parallel_efficiency: float
    quality_score: float
    warnings: List[str] = field(default_factory=list)
    
    def get_summary(self) -> Dict[str, Any]:
        """Get summary of segmentation results"""
        return {
            "strategy_name": self.strategy.name,
            "segment_size": self.strategy.segment_size,
            "overlap_size": self.strategy.overlap_size,
            "overlap_percentage": self.strategy.get_overlap_percentage(),
            "estimated_segments": self.estimated_segments,
            "estimated_runtime_hours": self.estimated_runtime_hours,
            "memory_usage_gb": self.memory_usage_gb,
            "parallel_efficiency": self.parallel_efficiency,
            "quality_score": self.quality_score,
            "warnings": self.warnings
        }


class EVMSegmentationOptimizer:
    """
    EVM Segmentation Parameter Optimizer
    
    This class analyzes genome characteristics and optimizes segmentation
    parameters for efficient EVM processing.
    """
    
    def __init__(self, output_dir: str = "segmentation_analysis"):
        """
        Initialize the segmentation optimizer
        
        Args:
            output_dir: Directory for analysis outputs
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        self.genome_stats = GenomeStatistics()
        self.strategies = []
        self.analysis_results = []
        
        logger.info(f"EVMSegmentationOptimizer initialized with output directory: {self.output_dir}")
    
    def analyze_genome_statistics(self, genome_file: str, gff_files: List[str] = None) -> GenomeStatistics:
        """
        Analyze genome statistics for segmentation optimization
        
        Args:
            genome_file: Path to genome FASTA file
            gff_files: Optional list of GFF files for gene analysis
            
        Returns:
            GenomeStatistics object with analyzed data
        """
        logger.info(f"Analyzing genome statistics from {genome_file}")
        
        # Analyze genome FASTA
        self._analyze_fasta_file(genome_file)
        
        # Analyze gene features if GFF files provided
        if gff_files:
            self._analyze_gff_files(gff_files)
        
        # Save statistics
        self._save_genome_statistics()
        
        logger.info(f"Genome analysis completed: {self.genome_stats.total_length:,} bp, "
                   f"{self.genome_stats.num_chromosomes} chromosomes")
        
        return self.genome_stats
    
    def _analyze_fasta_file(self, genome_file: str):
        """Analyze FASTA file for basic statistics"""
        try:
            with open(genome_file, 'r') as f:
                current_chr = None
                current_length = 0
                total_length = 0
                gc_count = 0
                
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        # Save previous chromosome
                        if current_chr is not None:
                            self.genome_stats.chromosome_lengths[current_chr] = current_length
                            total_length += current_length
                        
                        # Start new chromosome
                        current_chr = line[1:].split()[0]  # Get chromosome name
                        current_length = 0
                    else:
                        # Count sequence
                        current_length += len(line)
                        gc_count += line.upper().count('G') + line.upper().count('C')
                
                # Save last chromosome
                if current_chr is not None:
                    self.genome_stats.chromosome_lengths[current_chr] = current_length
                    total_length += current_length
                
                self.genome_stats.total_length = total_length
                self.genome_stats.num_chromosomes = len(self.genome_stats.chromosome_lengths)
                self.genome_stats.gc_content = (gc_count / total_length) * 100 if total_length > 0 else 0
                
        except Exception as e:
            logger.error(f"Error analyzing FASTA file: {e}")
            raise
    
    def _analyze_gff_files(self, gff_files: List[str]):
        """Analyze GFF files for gene statistics"""
        gene_lengths = []
        gene_count = 0
        
        for gff_file in gff_files:
            try:
                with open(gff_file, 'r') as f:
                    for line in f:
                        if line.startswith('#') or not line.strip():
                            continue
                        
                        parts = line.strip().split('\t')
                        if len(parts) >= 5 and parts[2] == 'gene':
                            try:
                                start = int(parts[3])
                                end = int(parts[4])
                                gene_length = end - start + 1
                                gene_lengths.append(gene_length)
                                gene_count += 1
                            except ValueError:
                                continue
                                
            except Exception as e:
                logger.warning(f"Error analyzing GFF file {gff_file}: {e}")
                continue
        
        if gene_lengths:
            self.genome_stats.longest_gene_length = max(gene_lengths)
            self.genome_stats.gene_density = (gene_count / self.genome_stats.total_length) * 1000  # genes per kb
        
        logger.info(f"Analyzed {gene_count} genes, longest: {self.genome_stats.longest_gene_length:,} bp")
    
    def _save_genome_statistics(self):
        """Save genome statistics to file"""
        stats_file = self.output_dir / "genome_statistics.json"
        
        stats_data = {
            "total_length": self.genome_stats.total_length,
            "num_chromosomes": self.genome_stats.num_chromosomes,
            "chromosome_lengths": self.genome_stats.chromosome_lengths,
            "longest_gene_length": self.genome_stats.longest_gene_length,
            "gene_density": self.genome_stats.gene_density,
            "gc_content": self.genome_stats.gc_content,
            "repeat_content": self.genome_stats.repeat_content,
            "analysis_timestamp": datetime.now().isoformat()
        }
        
        with open(stats_file, 'w') as f:
            json.dump(stats_data, f, indent=2)
        
        logger.info(f"Genome statistics saved to {stats_file}")
    
    def generate_segmentation_strategies(self) -> List[SegmentationStrategy]:
        """
        Generate multiple segmentation strategies based on genome characteristics
        
        Returns:
            List of segmentation strategies to evaluate
        """
        strategies = []
        
        # Get basic genome metrics
        total_length = self.genome_stats.total_length
        longest_gene = self.genome_stats.longest_gene_length
        avg_chr_length = self.genome_stats.get_average_chromosome_length()
        
        # Strategy 1: Conservative (small segments, large overlap)
        conservative_segment = max(100000, longest_gene * 2)
        conservative_overlap = max(25000, longest_gene)
        strategies.append(SegmentationStrategy(
            name="Conservative",
            segment_size=conservative_segment,
            overlap_size=conservative_overlap,
            min_intergenic_length=1000
        ))
        
        # Strategy 2: Balanced (medium segments, medium overlap)
        balanced_segment = max(200000, longest_gene * 3)
        balanced_overlap = max(50000, longest_gene)
        strategies.append(SegmentationStrategy(
            name="Balanced",
            segment_size=balanced_segment,
            overlap_size=balanced_overlap,
            min_intergenic_length=1000
        ))
        
        # Strategy 3: Aggressive (large segments, small overlap)
        aggressive_segment = max(500000, longest_gene * 4)
        aggressive_overlap = max(100000, longest_gene)
        strategies.append(SegmentationStrategy(
            name="Aggressive",
            segment_size=aggressive_segment,
            overlap_size=aggressive_overlap,
            min_intergenic_length=1000
        ))
        
        # Strategy 4: Memory-optimized (very small segments)
        memory_segment = max(50000, longest_gene * 1.5)
        memory_overlap = max(15000, longest_gene // 2)
        strategies.append(SegmentationStrategy(
            name="Memory-Optimized",
            segment_size=int(memory_segment),
            overlap_size=int(memory_overlap),
            min_intergenic_length=500
        ))
        
        # Strategy 5: Speed-optimized (larger segments)
        speed_segment = max(1000000, longest_gene * 5)
        speed_overlap = max(200000, longest_gene)
        strategies.append(SegmentationStrategy(
            name="Speed-Optimized",
            segment_size=speed_segment,
            overlap_size=speed_overlap,
            min_intergenic_length=2000
        ))
        
        # Strategy 6: Adaptive (varies by chromosome size)
        if avg_chr_length > 0:
            adaptive_segment = min(int(avg_chr_length / 10), 500000)
            adaptive_overlap = min(adaptive_segment // 4, 100000)
            strategies.append(SegmentationStrategy(
                name="Adaptive",
                segment_size=max(adaptive_segment, longest_gene * 2),
                overlap_size=max(adaptive_overlap, longest_gene),
                adaptive_sizing=True
            ))
        
        # Validate all strategies
        for strategy in strategies:
            try:
                strategy.validate()
            except ValueError as e:
                logger.warning(f"Invalid strategy {strategy.name}: {e}")
                continue
        
        self.strategies = strategies
        logger.info(f"Generated {len(strategies)} segmentation strategies")
        
        return strategies
    
    def evaluate_segmentation_strategy(self, strategy: SegmentationStrategy) -> SegmentationResult:
        """
        Evaluate a segmentation strategy for performance and efficiency
        
        Args:
            strategy: Segmentation strategy to evaluate
            
        Returns:
            SegmentationResult with performance metrics
        """
        # Calculate estimated number of segments
        total_segments = 0
        for chr_length in self.genome_stats.chromosome_lengths.values():
            if strategy.adaptive_sizing:
                # Adaptive sizing based on chromosome length
                adaptive_segment_size = min(strategy.segment_size, chr_length // 5)
                adaptive_segment_size = max(adaptive_segment_size, self.genome_stats.longest_gene_length * 2)
                segments_per_chr = max(1, (chr_length - strategy.overlap_size) // 
                                     (adaptive_segment_size - strategy.overlap_size))
            else:
                segments_per_chr = max(1, (chr_length - strategy.overlap_size) // 
                                     (strategy.segment_size - strategy.overlap_size))
            
            total_segments += segments_per_chr
        
        # Estimate runtime (based on segment count and size)
        base_time_per_segment = 0.1  # hours per segment (rough estimate)
        estimated_runtime = total_segments * base_time_per_segment
        
        # Estimate memory usage
        memory_per_segment = (strategy.segment_size / 1000000) * 0.5  # GB per Mb
        max_parallel_segments = min(total_segments, 8)  # Assume max 8 parallel processes
        estimated_memory = memory_per_segment * max_parallel_segments
        
        # Calculate parallel efficiency
        parallel_efficiency = min(1.0, 8 / max(1, total_segments))
        
        # Calculate quality score (balance of speed, memory, and completeness)
        overlap_ratio = strategy.overlap_size / strategy.segment_size
        segment_size_score = 1.0 - min(1.0, strategy.segment_size / 2000000)  # Prefer smaller segments
        overlap_score = min(1.0, overlap_ratio * 4)  # Prefer reasonable overlap
        efficiency_score = parallel_efficiency
        
        quality_score = (segment_size_score + overlap_score + efficiency_score) / 3
        
        # Generate warnings
        warnings = []
        if strategy.overlap_size < self.genome_stats.longest_gene_length:
            warnings.append(f"Overlap size ({strategy.overlap_size:,}) smaller than longest gene ({self.genome_stats.longest_gene_length:,})")
        
        if total_segments > 1000:
            warnings.append(f"Very high segment count ({total_segments}) may impact performance")
        
        if estimated_memory > 64:
            warnings.append(f"High memory usage ({estimated_memory:.1f} GB) may require large memory systems")
        
        if strategy.get_overlap_percentage() > 50:
            warnings.append(f"High overlap percentage ({strategy.get_overlap_percentage():.1f}%) increases computational overhead")
        
        return SegmentationResult(
            strategy=strategy,
            estimated_segments=total_segments,
            estimated_runtime_hours=estimated_runtime,
            memory_usage_gb=estimated_memory,
            parallel_efficiency=parallel_efficiency,
            quality_score=quality_score,
            warnings=warnings
        )
    
    def optimize_segmentation_parameters(self) -> SegmentationResult:
        """
        Find the optimal segmentation parameters for the genome
        
        Returns:
            Best segmentation result
        """
        logger.info("Optimizing segmentation parameters...")
        
        if not self.strategies:
            self.generate_segmentation_strategies()
        
        # Evaluate all strategies
        results = []
        for strategy in self.strategies:
            result = self.evaluate_segmentation_strategy(strategy)
            results.append(result)
            
            logger.info(f"Strategy '{strategy.name}': "
                       f"{result.estimated_segments} segments, "
                       f"{result.estimated_runtime_hours:.1f}h runtime, "
                       f"{result.memory_usage_gb:.1f}GB memory, "
                       f"quality={result.quality_score:.2f}")
        
        self.analysis_results = results
        
        # Find best strategy (highest quality score with reasonable constraints)
        valid_results = [r for r in results if len(r.warnings) <= 2 and r.memory_usage_gb <= 32]
        if not valid_results:
            valid_results = results  # Fall back to all results if none meet constraints
        
        best_result = max(valid_results, key=lambda r: r.quality_score)
        
        logger.info(f"Optimal strategy: '{best_result.strategy.name}' "
                   f"(quality score: {best_result.quality_score:.2f})")
        
        return best_result
    
    def save_optimization_report(self, best_result: SegmentationResult) -> str:
        """
        Save detailed optimization report
        
        Args:
            best_result: Best segmentation result
            
        Returns:
            Path to saved report
        """
        report_file = self.output_dir / "segmentation_optimization_report.json"
        
        report_data = {
            "analysis_timestamp": datetime.now().isoformat(),
            "genome_statistics": {
                "total_length": self.genome_stats.total_length,
                "num_chromosomes": self.genome_stats.num_chromosomes,
                "longest_gene_length": self.genome_stats.longest_gene_length,
                "gene_density": self.genome_stats.gene_density,
                "gc_content": self.genome_stats.gc_content
            },
            "best_strategy": best_result.get_summary(),
            "all_strategies": [result.get_summary() for result in self.analysis_results],
            "recommendations": self._generate_recommendations(best_result)
        }
        
        with open(report_file, 'w') as f:
            json.dump(report_data, f, indent=2)
        
        # Also create human-readable report
        text_report = self._create_text_report(best_result)
        text_file = self.output_dir / "segmentation_optimization_report.txt"
        with open(text_file, 'w') as f:
            f.write(text_report)
        
        logger.info(f"Optimization report saved to {report_file}")
        return str(report_file)
    
    def _generate_recommendations(self, best_result: SegmentationResult) -> List[str]:
        """Generate recommendations based on analysis"""
        recommendations = []
        
        strategy = best_result.strategy
        
        # Segment size recommendations
        if strategy.segment_size < 100000:
            recommendations.append("Consider increasing segment size if memory allows for better performance")
        elif strategy.segment_size > 1000000:
            recommendations.append("Large segment size may cause memory issues on smaller systems")
        
        # Overlap recommendations
        overlap_ratio = strategy.get_overlap_percentage()
        if overlap_ratio < 10:
            recommendations.append("Consider increasing overlap to ensure gene completeness")
        elif overlap_ratio > 30:
            recommendations.append("High overlap may increase computational overhead unnecessarily")
        
        # Performance recommendations
        if best_result.estimated_segments > 500:
            recommendations.append("High segment count - consider using distributed computing for better performance")
        
        if best_result.memory_usage_gb > 16:
            recommendations.append("High memory usage - ensure adequate RAM or consider memory-optimized strategy")
        
        # Genome-specific recommendations
        if self.genome_stats.gene_density > 0.1:  # High gene density
            recommendations.append("High gene density detected - smaller segments may improve accuracy")
        
        if len(best_result.warnings) > 0:
            recommendations.append("Review warnings and consider adjusting parameters if needed")
        
        return recommendations
    
    def _create_text_report(self, best_result: SegmentationResult) -> str:
        """Create human-readable text report"""
        lines = [
            "EVM Segmentation Optimization Report",
            "=" * 50,
            f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            "",
            "Genome Statistics:",
            f"  Total Length: {self.genome_stats.total_length:,} bp",
            f"  Chromosomes: {self.genome_stats.num_chromosomes}",
            f"  Longest Gene: {self.genome_stats.longest_gene_length:,} bp",
            f"  Gene Density: {self.genome_stats.gene_density:.3f} genes/kb",
            f"  GC Content: {self.genome_stats.gc_content:.1f}%",
            "",
            "Recommended Strategy:",
            f"  Name: {best_result.strategy.name}",
            f"  Segment Size: {best_result.strategy.segment_size:,} bp",
            f"  Overlap Size: {best_result.strategy.overlap_size:,} bp",
            f"  Overlap Percentage: {best_result.strategy.get_overlap_percentage():.1f}%",
            f"  Min Intergenic Length: {best_result.strategy.min_intergenic_length:,} bp",
            "",
            "Performance Estimates:",
            f"  Estimated Segments: {best_result.estimated_segments:,}",
            f"  Estimated Runtime: {best_result.estimated_runtime_hours:.1f} hours",
            f"  Memory Usage: {best_result.memory_usage_gb:.1f} GB",
            f"  Parallel Efficiency: {best_result.parallel_efficiency:.1%}",
            f"  Quality Score: {best_result.quality_score:.2f}",
            ""
        ]
        
        if best_result.warnings:
            lines.extend([
                "Warnings:",
                *[f"  - {warning}" for warning in best_result.warnings],
                ""
            ])
        
        recommendations = self._generate_recommendations(best_result)
        if recommendations:
            lines.extend([
                "Recommendations:",
                *[f"  - {rec}" for rec in recommendations],
                ""
            ])
        
        lines.extend([
            "All Evaluated Strategies:",
            "-" * 30
        ])
        
        for result in self.analysis_results:
            lines.extend([
                f"Strategy: {result.strategy.name}",
                f"  Segment Size: {result.strategy.segment_size:,} bp",
                f"  Overlap Size: {result.strategy.overlap_size:,} bp",
                f"  Segments: {result.estimated_segments:,}",
                f"  Runtime: {result.estimated_runtime_hours:.1f}h",
                f"  Memory: {result.memory_usage_gb:.1f}GB",
                f"  Quality: {result.quality_score:.2f}",
                ""
            ])
        
        return "\n".join(lines)


def optimize_for_nbs_genome(genome_file: str = None, 
                           gff_files: List[str] = None,
                           output_dir: str = "nbs_segmentation_analysis") -> SegmentationResult:
    """
    Convenience function to optimize segmentation for NBS genome annotation
    
    Args:
        genome_file: Path to genome FASTA file
        gff_files: List of GFF files with gene annotations
        output_dir: Output directory for analysis
        
    Returns:
        Best segmentation result for NBS annotation
    """
    optimizer = EVMSegmentationOptimizer(output_dir)
    
    # Use default values if files not provided
    if genome_file and os.path.exists(genome_file):
        optimizer.analyze_genome_statistics(genome_file, gff_files)
    else:
        # Use typical rice genome statistics as defaults
        optimizer.genome_stats = GenomeStatistics(
            total_length=373245519,  # Rice genome size
            num_chromosomes=12,
            longest_gene_length=150000,  # Typical long gene
            gene_density=0.08,  # genes per kb
            gc_content=43.5
        )
        logger.info("Using default rice genome statistics for optimization")
    
    # Generate and evaluate strategies
    best_result = optimizer.optimize_segmentation_parameters()
    
    # Save report
    optimizer.save_optimization_report(best_result)
    
    return best_result


if __name__ == "__main__":
    # Example usage
    result = optimize_for_nbs_genome()
    print(f"Optimal segmentation strategy: {result.strategy.name}")
    print(f"Segment size: {result.strategy.segment_size:,} bp")
    print(f"Overlap size: {result.strategy.overlap_size:,} bp")
    print(f"Quality score: {result.quality_score:.2f}") 