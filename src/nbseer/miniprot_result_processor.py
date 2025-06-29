#!/usr/bin/env python3
"""
Miniprot Result Processor for EVM Integration

This module processes miniprot prediction results, performs quality-based
classification and filtering, and generates EVM-compatible output files
with appropriate weights based on prediction quality.

Author: NBS Annotation Pipeline
Date: 2025-06-24
"""

import re
import json
import logging
from pathlib import Path
from dataclasses import dataclass, asdict
from typing import Dict, List, Optional, Tuple, Set
from collections import defaultdict
import gzip


@dataclass
class MiniprotPrediction:
    """Data structure for a single miniprot prediction"""
    # Basic information
    prediction_id: str
    chromosome: str
    start: int
    end: int
    strand: str
    
    # Quality metrics from mRNA feature
    identity: float
    positive: float
    rank: int
    target_info: str
    
    # Quality metrics from PAF header
    frameshift_count: int = 0
    stop_codon_count: int = 0
    
    # Classification
    quality_class: str = ""
    weight: float = 0.0
    
    # Raw GFF3 lines
    paf_header: str = ""
    gff3_lines: List[str] = None
    
    def __post_init__(self):
        if self.gff3_lines is None:
            self.gff3_lines = []
        if not (0.0 <= self.identity <= 1.0):
            raise ValueError(f"Identity score {self.identity} must be between 0.0 and 1.0")
        if not (0.0 <= self.positive <= 1.0):
            raise ValueError(f"Positive score {self.positive} must be between 0.0 and 1.0")
        if self.frameshift_count < 0:
            raise ValueError(f"Frameshift count {self.frameshift_count} cannot be negative")
        if self.stop_codon_count < 0:
            raise ValueError(f"Stop codon count {self.stop_codon_count} cannot be negative")


@dataclass
class QualityStats:
    """Statistics for quality assessment"""
    total_predictions: int = 0
    high_quality: int = 0
    medium_quality: int = 0
    low_quality: int = 0
    filtered_out: int = 0
    
    # Detailed filtering reasons
    low_identity: int = 0
    frameshift_issues: int = 0
    stop_codon_issues: int = 0
    combined_issues: int = 0
    
    def to_dict(self) -> Dict:
        return asdict(self)


@dataclass
class QualityThresholds:
    """Configurable quality thresholds"""
    high_identity_min: float = 0.9
    medium_identity_min: float = 0.7
    low_identity_min: float = 0.5
    
    high_weight: float = 1.0
    medium_weight: float = 0.7
    low_weight: float = 0.4
    
    max_frameshift: int = 0
    max_stop_codons: int = 0


class MiniprotResultProcessor:
    """
    Process miniprot prediction results for EVM integration
    
    This class handles:
    1. Parsing GFF3 files with PAF headers
    2. Extracting quality metrics
    3. Quality-based classification and filtering
    4. Generating quality-separated output files
    5. Creating EVM weight configurations
    """
    
    def __init__(self, 
                 input_gff3: str,
                 output_dir: str = "results",
                 config_dir: str = "evm_config",
                 thresholds: Optional[QualityThresholds] = None,
                 logger: Optional[logging.Logger] = None):
        """
        Initialize the processor
        
        Args:
            input_gff3: Path to input miniprot GFF3 file
            output_dir: Directory for output files
            config_dir: Directory for EVM configuration files
            thresholds: Quality thresholds configuration
            logger: Logger instance
        """
        self.input_gff3 = Path(input_gff3)
        self.output_dir = Path(output_dir)
        self.config_dir = Path(config_dir)
        self.thresholds = thresholds or QualityThresholds()
        
        # Setup logging
        self.logger = logger or self._setup_logger()
        
        # Create output directories
        self.output_dir.mkdir(exist_ok=True)
        self.config_dir.mkdir(exist_ok=True)
        
        # Initialize data structures
        self.predictions: Dict[str, MiniprotPrediction] = {}
        self.paf_data: Dict[str, Dict] = {}  # target_id -> PAF info
        self.stats = QualityStats()
        
        self.logger.info(f"Initialized MiniprotResultProcessor")
        self.logger.info(f"Input: {self.input_gff3}")
        self.logger.info(f"Output directory: {self.output_dir}")
        self.logger.info(f"Config directory: {self.config_dir}")
    
    def _setup_logger(self) -> logging.Logger:
        """Setup default logger"""
        logger = logging.getLogger("MiniprotProcessor")
        if not logger.handlers:
            handler = logging.StreamHandler()
            formatter = logging.Formatter(
                '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
            )
            handler.setFormatter(formatter)
            logger.addHandler(handler)
            logger.setLevel(logging.INFO)
        return logger
    
    def process_miniprot_results(self) -> QualityStats:
        """
        Main processing pipeline
        
        Returns:
            QualityStats: Processing statistics
        """
        self.logger.info("Starting miniprot results processing...")
        
        try:
            # Step 1: Parse GFF3 file
            self.logger.info("Step 1: Parsing GFF3 file...")
            self._parse_gff3_file()
            
            # Step 2: Extract quality metrics
            self.logger.info("Step 2: Extracting quality metrics...")
            self._extract_quality_metrics()
            
            # Step 3: Classify predictions by quality
            self.logger.info("Step 3: Classifying predictions by quality...")
            self._classify_predictions()
            
            # Step 4: Generate quality-separated output files
            self.logger.info("Step 4: Generating output files...")
            self._generate_output_files()
            
            # Step 5: Generate EVM weight configuration
            self.logger.info("Step 5: Generating EVM weight configuration...")
            self._generate_evm_weights()
            
            # Step 6: Generate quality report
            self.logger.info("Step 6: Generating quality report...")
            self._generate_quality_report()
            
            self.logger.info("Processing completed successfully!")
            self._log_statistics()
            
            return self.stats
            
        except Exception as e:
            self.logger.error(f"Error during processing: {e}")
            raise
    
    def _parse_gff3_file(self):
        """Parse the GFF3 file and extract PAF headers and features"""
        self.logger.info(f"Parsing GFF3 file: {self.input_gff3}")
        
        current_paf_data = {}
        current_prediction_lines = []
        current_prediction_id = None
        
        with self._open_file(self.input_gff3) as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line:
                    continue
                
                if line.startswith('##PAF'):
                    # Save previous prediction if exists
                    if current_prediction_id and current_prediction_lines:
                        self._save_current_prediction(
                            current_prediction_id, 
                            current_paf_data, 
                            current_prediction_lines
                        )
                    
                    # Parse new PAF header
                    current_paf_data = self._parse_paf_header(line)
                    current_prediction_lines = [line]
                    current_prediction_id = None
                    
                elif line.startswith('##'):
                    # Other headers, add to current prediction
                    if current_prediction_lines:
                        current_prediction_lines.append(line)
                
                elif not line.startswith('#'):
                    # GFF3 feature line
                    fields = line.split('\t')
                    if len(fields) >= 9:
                        feature_type = fields[2]
                        if feature_type == 'mRNA':
                            # Extract prediction ID from mRNA line
                            attributes = self._parse_attributes(fields[8])
                            current_prediction_id = attributes.get('ID', '')
                        
                        current_prediction_lines.append(line)
                
                # Progress reporting
                if line_num % 1000000 == 0:
                    self.logger.info(f"Processed {line_num:,} lines...")
        
        # Save last prediction
        if current_prediction_id and current_prediction_lines:
            self._save_current_prediction(
                current_prediction_id, 
                current_paf_data, 
                current_prediction_lines
            )
        
        self.logger.info(f"Parsed {len(self.predictions):,} predictions")
    
    def _parse_paf_header(self, paf_line: str) -> Dict:
        """Parse PAF header line to extract quality metrics"""
        # Example PAF line:
        # ##PAF   2767407 283     0       277     +       Chr8    28443022        25419984        25421274        6998
        # 34      0       AS:i:1098       ms:i:1228       np:i:249        fs:i:0  st:i:0  da:i:51 do:i:126
        
        paf_data = {}
        
        # Extract fs (frameshift) and st (stop codon) counts
        fs_match = re.search(r'fs:i:(\d+)', paf_line)
        st_match = re.search(r'st:i:(\d+)', paf_line)
        
        paf_data['frameshift_count'] = int(fs_match.group(1)) if fs_match else 0
        paf_data['stop_codon_count'] = int(st_match.group(1)) if st_match else 0
        
        # Extract target ID (first field after ##PAF)
        fields = paf_line.split()
        if len(fields) > 1:
            paf_data['target_id'] = fields[1]
        
        return paf_data
    
    def _parse_attributes(self, attr_string: str) -> Dict[str, str]:
        """Parse GFF3 attributes string"""
        attributes = {}
        for item in attr_string.split(';'):
            if '=' in item:
                key, value = item.split('=', 1)
                attributes[key] = value
        return attributes
    
    def _save_current_prediction(self, prediction_id: str, paf_data: Dict, gff3_lines: List[str]):
        """Save current prediction data"""
        if not prediction_id:
            return
        
        # Find mRNA line to extract quality metrics
        mrna_line = None
        for line in gff3_lines:
            if '\tmRNA\t' in line:
                mrna_line = line
                break
        
        if not mrna_line:
            return
        
        # Parse mRNA line
        fields = mrna_line.split('\t')
        if len(fields) < 9:
            return
        
        attributes = self._parse_attributes(fields[8])
        
        # Create prediction object
        prediction = MiniprotPrediction(
            prediction_id=prediction_id,
            chromosome=fields[0],
            start=int(fields[3]),
            end=int(fields[4]),
            strand=fields[6],
            identity=float(attributes.get('Identity', '0.0')),
            positive=float(attributes.get('Positive', '0.0')),
            rank=int(attributes.get('Rank', '1')),
            target_info=attributes.get('Target', ''),
            frameshift_count=paf_data.get('frameshift_count', 0),
            stop_codon_count=paf_data.get('stop_codon_count', 0),
            gff3_lines=gff3_lines.copy()
        )
        
        self.predictions[prediction_id] = prediction
    
    def _extract_quality_metrics(self):
        """Extract and validate quality metrics for all predictions"""
        self.logger.info("Extracting quality metrics...")
        
        for pred_id, prediction in self.predictions.items():
            # Quality metrics are already extracted during parsing
            # Just validate ranges
            if not (0.0 <= prediction.identity <= 1.0):
                self.logger.warning(f"Invalid identity score for {pred_id}: {prediction.identity}")
            
            if not (0.0 <= prediction.positive <= 1.0):
                self.logger.warning(f"Invalid positive score for {pred_id}: {prediction.positive}")
    
    def _classify_predictions(self):
        """Classify predictions based on quality thresholds"""
        self.logger.info("Classifying predictions by quality...")
        
        for pred_id, prediction in self.predictions.items():
            self.stats.total_predictions += 1
            
            # Check if should be filtered out
            should_filter = (
                prediction.identity < self.thresholds.low_identity_min or
                prediction.frameshift_count > self.thresholds.max_frameshift or
                prediction.stop_codon_count > self.thresholds.max_stop_codons
            )
            
            if should_filter:
                prediction.quality_class = "filtered"
                prediction.weight = 0.0
                self.stats.filtered_out += 1
                
                # Count filtering reasons
                if prediction.identity < self.thresholds.low_identity_min:
                    self.stats.low_identity += 1
                if prediction.frameshift_count > self.thresholds.max_frameshift:
                    self.stats.frameshift_issues += 1
                if prediction.stop_codon_count > self.thresholds.max_stop_codons:
                    self.stats.stop_codon_issues += 1
                if (prediction.frameshift_count > 0 and prediction.stop_codon_count > 0):
                    self.stats.combined_issues += 1
                
            elif prediction.identity >= self.thresholds.high_identity_min:
                prediction.quality_class = "high"
                prediction.weight = self.thresholds.high_weight
                self.stats.high_quality += 1
                
            elif prediction.identity >= self.thresholds.medium_identity_min:
                prediction.quality_class = "medium"
                prediction.weight = self.thresholds.medium_weight
                self.stats.medium_quality += 1
                
            else:
                prediction.quality_class = "low"
                prediction.weight = self.thresholds.low_weight
                self.stats.low_quality += 1
    
    def _generate_output_files(self):
        """Generate quality-separated GFF3 output files"""
        self.logger.info("Generating quality-separated GFF3 files...")
        
        # Group predictions by quality class
        quality_groups = defaultdict(list)
        for prediction in self.predictions.values():
            if prediction.quality_class != "filtered":
                quality_groups[prediction.quality_class].append(prediction)
        
        # Write separate files for each quality class
        for quality_class, predictions in quality_groups.items():
            output_file = self.output_dir / f"miniprot_{quality_class}_quality.gff3"
            
            self.logger.info(f"Writing {len(predictions):,} {quality_class} quality predictions to {output_file}")
            
            with open(output_file, 'w') as f:
                # Write GFF3 header
                f.write("##gff-version 3\n")
                f.write(f"# Miniprot {quality_class} quality predictions\n")
                f.write(f"# Generated by MiniprotResultProcessor\n")
                f.write(f"# Total predictions: {len(predictions):,}\n")
                f.write(f"# Quality class: {quality_class}\n")
                f.write(f"# Weight: {predictions[0].weight if predictions else 0.0}\n")
                f.write("\n")
                
                # Write predictions
                for prediction in predictions:
                    # Write only the GFF3 feature lines (skip PAF headers)
                    for line in prediction.gff3_lines:
                        if not line.startswith('##PAF'):
                            f.write(line + '\n')
                    f.write('\n')  # Separate predictions
    
    def _generate_evm_weights(self):
        """Generate EVM weight configuration file"""
        weights_file = self.config_dir / "miniprot_weights.txt"
        
        self.logger.info(f"Generating EVM weights file: {weights_file}")
        
        with open(weights_file, 'w') as f:
            f.write("# EVM Weight Configuration for Miniprot Predictions\n")
            f.write("# Generated by MiniprotResultProcessor\n")
            f.write(f"# Processing date: {self._get_timestamp()}\n")
            f.write("#\n")
            f.write("# Format: ABINITIO_PREDICTION <program> <weight>\n")
            f.write("#         PROTEIN <program> <weight>\n")
            f.write("#\n")
            
            # Count predictions by quality
            quality_counts = defaultdict(int)
            for prediction in self.predictions.values():
                if prediction.quality_class != "filtered":
                    quality_counts[prediction.quality_class] += 1
            
            # Write weights for each quality class
            if quality_counts['high'] > 0:
                f.write(f"PROTEIN\tminiprot_high\t{self.thresholds.high_weight}\n")
                f.write(f"# High quality: {quality_counts['high']:,} predictions (Identity >= {self.thresholds.high_identity_min})\n")
            
            if quality_counts['medium'] > 0:
                f.write(f"PROTEIN\tminiprot_medium\t{self.thresholds.medium_weight}\n")
                f.write(f"# Medium quality: {quality_counts['medium']:,} predictions ({self.thresholds.medium_identity_min} <= Identity < {self.thresholds.high_identity_min})\n")
            
            if quality_counts['low'] > 0:
                f.write(f"PROTEIN\tminiprot_low\t{self.thresholds.low_weight}\n")
                f.write(f"# Low quality: {quality_counts['low']:,} predictions ({self.thresholds.low_identity_min} <= Identity < {self.thresholds.medium_identity_min})\n")
            
            f.write(f"#\n")
            f.write(f"# Filtered out: {self.stats.filtered_out:,} predictions\n")
            f.write(f"# Total processed: {self.stats.total_predictions:,} predictions\n")
    
    def _generate_quality_report(self):
        """Generate comprehensive quality assessment report"""
        report_file = self.output_dir / "miniprot_quality_report.json"
        
        self.logger.info(f"Generating quality report: {report_file}")
        
        # Calculate percentages
        total = self.stats.total_predictions
        report_data = {
            "processing_info": {
                "input_file": str(self.input_gff3),
                "processing_date": self._get_timestamp(),
                "total_predictions": total
            },
            "quality_thresholds": asdict(self.thresholds),
            "quality_statistics": {
                "total_predictions": total,
                "high_quality": {
                    "count": self.stats.high_quality,
                    "percentage": round(self.stats.high_quality / total * 100, 2) if total > 0 else 0
                },
                "medium_quality": {
                    "count": self.stats.medium_quality,
                    "percentage": round(self.stats.medium_quality / total * 100, 2) if total > 0 else 0
                },
                "low_quality": {
                    "count": self.stats.low_quality,
                    "percentage": round(self.stats.low_quality / total * 100, 2) if total > 0 else 0
                },
                "filtered_out": {
                    "count": self.stats.filtered_out,
                    "percentage": round(self.stats.filtered_out / total * 100, 2) if total > 0 else 0
                }
            },
            "filtering_details": {
                "low_identity": self.stats.low_identity,
                "frameshift_issues": self.stats.frameshift_issues,
                "stop_codon_issues": self.stats.stop_codon_issues,
                "combined_issues": self.stats.combined_issues
            },
            "output_files": {
                "high_quality_gff3": str(self.output_dir / "miniprot_high_quality.gff3"),
                "medium_quality_gff3": str(self.output_dir / "miniprot_medium_quality.gff3"),
                "low_quality_gff3": str(self.output_dir / "miniprot_low_quality.gff3"),
                "evm_weights": str(self.config_dir / "miniprot_weights.txt"),
                "quality_report": str(report_file)
            }
        }
        
        with open(report_file, 'w') as f:
            json.dump(report_data, f, indent=2)
    
    def _get_timestamp(self) -> str:
        """Get current timestamp string"""
        from datetime import datetime
        return datetime.now().isoformat()
    
    def _log_statistics(self):
        """Log processing statistics"""
        total = self.stats.total_predictions
        
        self.logger.info("="*60)
        self.logger.info("MINIPROT PROCESSING STATISTICS")
        self.logger.info("="*60)
        self.logger.info(f"Total predictions processed: {total:,}")
        self.logger.info(f"High quality (Identity ≥ {self.thresholds.high_identity_min}): {self.stats.high_quality:,} ({self.stats.high_quality/total*100:.1f}%)")
        self.logger.info(f"Medium quality ({self.thresholds.medium_identity_min} ≤ Identity < {self.thresholds.high_identity_min}): {self.stats.medium_quality:,} ({self.stats.medium_quality/total*100:.1f}%)")
        self.logger.info(f"Low quality ({self.thresholds.low_identity_min} ≤ Identity < {self.thresholds.medium_identity_min}): {self.stats.low_quality:,} ({self.stats.low_quality/total*100:.1f}%)")
        self.logger.info(f"Filtered out: {self.stats.filtered_out:,} ({self.stats.filtered_out/total*100:.1f}%)")
        self.logger.info("")
        self.logger.info("Filtering breakdown:")
        self.logger.info(f"  Low identity (< {self.thresholds.low_identity_min}): {self.stats.low_identity:,}")
        self.logger.info(f"  Frameshift issues (fs > {self.thresholds.max_frameshift}): {self.stats.frameshift_issues:,}")
        self.logger.info(f"  Stop codon issues (st > {self.thresholds.max_stop_codons}): {self.stats.stop_codon_issues:,}")
        self.logger.info(f"  Combined issues: {self.stats.combined_issues:,}")
        self.logger.info("="*60)

    def _open_file(self, filepath: Path):
        """Open file, handling both regular and gzipped files"""
        if filepath.suffix == '.gz':
            return gzip.open(filepath, 'rt')
        else:
            return open(filepath, 'r')


def main():
    """Main function for command-line usage"""
    import argparse
    
    parser = argparse.ArgumentParser(description="Process miniprot results for EVM integration")
    parser.add_argument("input_gff3", help="Input miniprot GFF3 file")
    parser.add_argument("--output-dir", default="results", help="Output directory")
    parser.add_argument("--config-dir", default="evm_config", help="EVM config directory")
    parser.add_argument("--high-threshold", type=float, default=0.9, help="High quality threshold")
    parser.add_argument("--medium-threshold", type=float, default=0.7, help="Medium quality threshold")
    parser.add_argument("--low-threshold", type=float, default=0.5, help="Low quality threshold")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose logging")
    
    args = parser.parse_args()
    
    # Setup logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(level=log_level)
    
    # Create thresholds
    thresholds = QualityThresholds(
        high_identity_min=args.high_threshold,
        medium_identity_min=args.medium_threshold,
        low_identity_min=args.low_threshold
    )
    
    # Process results
    processor = MiniprotResultProcessor(
        input_gff3=args.input_gff3,
        output_dir=args.output_dir,
        config_dir=args.config_dir,
        thresholds=thresholds
    )
    
    stats = processor.process_miniprot_results()
    
    print(f"\nProcessing completed!")
    print(f"Total predictions: {stats.total_predictions:,}")
    print(f"High quality: {stats.high_quality:,}")
    print(f"Medium quality: {stats.medium_quality:,}")
    print(f"Low quality: {stats.low_quality:,}")
    print(f"Filtered out: {stats.filtered_out:,}")


if __name__ == "__main__":
    main()