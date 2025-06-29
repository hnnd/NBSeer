"""
Miniprot tool interface
Miniprot工具接口

提供蛋白质比对功能的工具接口
"""

from pathlib import Path
from typing import Any, Dict, List

from .base import ExternalTool, ToolResult
from ..utils.logging_setup import get_logger

logger = get_logger(__name__)


class MiniprotTool(ExternalTool):
    """
    Interface for miniprot tool
    Miniprot工具接口
    
    Miniprot用于将蛋白质序列比对到基因组序列，生成基于蛋白质的基因模型
    """
    
    def __init__(self, config: Dict[str, Any] = None) -> None:
        """
        Initialize miniprot tool interface
        
        Args:
            config: Tool configuration dictionary
        """
        super().__init__(
            tool_name="miniprot",
            config=config,
        )
    
    def prepare_input(self, **kwargs: Any) -> Dict[str, Path]:
        """
        Prepare input files for miniprot
        
        Args:
            genome_file: Path to genome FASTA file
            protein_file: Path to protein FASTA file
            
        Returns:
            Dictionary of prepared input files
        """
        genome_file = kwargs.get("genome_file")
        protein_file = kwargs.get("protein_file")
        
        if not genome_file:
            raise ValueError("genome_file is required for miniprot")
        if not protein_file:
            raise ValueError("protein_file is required for miniprot")
        
        genome_path = Path(genome_file)
        protein_path = Path(protein_file)
        
        self.validate_input(genome_file=genome_path, protein_file=protein_path)
        
        return {
            "genome": genome_path,
            "protein": protein_path,
        }
    
    def build_command(
        self,
        input_files: Dict[str, Path],
        output_dir: Path,
        **kwargs: Any,
    ) -> List[str]:
        """
        Build miniprot command
        
        Args:
            input_files: Prepared input files
            output_dir: Output directory
            **kwargs: Additional parameters
            
        Returns:
            Command as list of strings
        """
        genome_file = input_files["genome"]
        protein_file = input_files["protein"]
        output_file = output_dir / "alignments.gff"
        
        # Basic command structure
        command = [str(self.executable_path)]
        
        # Add parameters from config
        parameters = self.config.get("parameters", {})
        
        # Thread count
        threads = parameters.get("threads", 8)
        command.extend(["-t", str(threads)])
        
        # Output format - use --gff for GFF3 format and ensure PAF headers are included
        output_format = parameters.get("output_format", "gff")
        if output_format == "gff":
            command.append("--gff")
        elif output_format == "gtf":
            command.append("--gtf")
        
        # Score threshold for output
        outs = parameters.get("outs", 0.99)
        command.extend(["--outs", str(outs)])
        
        # Genome preset
        genome_preset = parameters.get("genome_preset")
        if genome_preset:
            command.extend(["-G", genome_preset])
        
        # Note: Removed -O (splicing_penalty), -s (min_dp_score), -c (min_chain_score) 
        # parameters as requested to use miniprot default settings
        
        # Maximum gap length
        max_gap = parameters.get("max_gap")
        if max_gap:
            command.extend(["-e", str(max_gap)])
        
        # Maximum intron length
        max_intron = parameters.get("max_intron")
        if max_intron:
            command.extend(["-G", str(max_intron)])
        
        # Add any additional parameters from kwargs
        for key, value in kwargs.items():
            if key not in ["genome_file", "protein_file"] and value is not None:
                if isinstance(value, bool):
                    if value:
                        command.append(f"--{key}")
                else:
                    command.extend([f"--{key}", str(value)])
        
        # Add input files (use absolute paths to avoid cwd issues)
        command.append(str(genome_file.resolve()))
        command.append(str(protein_file.resolve()))
        
        # Note: Output will be captured from stdout and written to file in execute method
        return command
    
    def _apply_quality_filter(self, output_file: Path) -> Path:
        """
        Apply quality filter to remove predictions with premature stops and frameshifts
        应用质量过滤，移除有提前终止和移码突变的预测
        """
        try:
            # Import the result processor
            import sys
            sys.path.append(str(Path(__file__).parent.parent))
            from miniprot_result_processor import MiniprotResultProcessor, QualityThresholds
            
            # Create filtered output file
            filtered_file = output_file.parent / "alignments_filtered.gff"
            
            # Set strict filtering thresholds
            thresholds = QualityThresholds(
                high_identity_min=0.95,     # 高质量：95%以上同一性
                medium_identity_min=0.85,   # 中等质量：85-95%同一性  
                low_identity_min=0.7,       # 低质量：70-85%同一性
                max_frameshift=0,           # 不允许移码突变
                max_stop_codons=0,          # 不允许提前终止
                high_weight=1.0,
                medium_weight=0.8,
                low_weight=0.5
            )
            
            # Process and filter results
            processor = MiniprotResultProcessor(
                input_gff3=str(output_file),
                output_dir=str(output_file.parent / "filtered"),
                thresholds=thresholds,
                logger=self.logger
            )
            
            stats = processor.process_miniprot_results()
            
            # Combine high and medium quality results into final filtered file
            with open(filtered_file, 'w') as outf:
                outf.write("##gff-version 3\n")
                outf.write("# Filtered miniprot alignments (no frameshifts or premature stops)\n")
                outf.write(f"# Original predictions: {stats.total_predictions}\n")
                outf.write(f"# High quality: {stats.high_quality}\n") 
                outf.write(f"# Medium quality: {stats.medium_quality}\n")
                outf.write(f"# Filtered out: {stats.filtered_out} (frameshifts: {stats.frameshift_issues}, stops: {stats.stop_codon_issues})\n")
                outf.write("\n")
                
                # Write high quality predictions
                high_file = output_file.parent / "filtered" / "miniprot_high_quality.gff3"
                if high_file.exists():
                    with open(high_file, 'r') as f:
                        for line in f:
                            if not line.startswith('#'):
                                outf.write(line)
                
                # Write medium quality predictions  
                medium_file = output_file.parent / "filtered" / "miniprot_medium_quality.gff3"
                if medium_file.exists():
                    with open(medium_file, 'r') as f:
                        for line in f:
                            if not line.startswith('#'):
                                outf.write(line)
            
            self.logger.info(f"Quality filtering completed:")
            self.logger.info(f"  Original: {stats.total_predictions} predictions")
            self.logger.info(f"  High quality: {stats.high_quality}")
            self.logger.info(f"  Medium quality: {stats.medium_quality}")
            self.logger.info(f"  Filtered out: {stats.filtered_out} (frameshifts: {stats.frameshift_issues}, stops: {stats.stop_codon_issues})")
            self.logger.info(f"  Final filtered file: {filtered_file}")
            
            return filtered_file
            
        except Exception as e:
            self.logger.warning(f"Quality filtering failed: {e}")
            self.logger.warning("Using original unfiltered results")
            return output_file

    def execute(self, output_dir: Path, timeout: int = None, **kwargs: Any) -> ToolResult:
        """
        Execute miniprot with shell redirection
        
        Args:
            output_dir: Directory for output files
            timeout: Execution timeout in seconds
            **kwargs: Tool-specific parameters
            
        Returns:
            ToolResult object
        """
        import subprocess
        from datetime import datetime
        
        start_time = datetime.now()
        
        try:
            # Validate inputs
            self.validate_input(**kwargs)
            
            # Prepare inputs
            input_files = self.prepare_input(**kwargs)
            
            # Build command
            command_list = self.build_command(input_files, output_dir, **kwargs)
            
            # Convert to shell command with redirection
            output_file = output_dir / "alignments.gff"
            command_without_redirect = []
            for item in command_list:
                if item not in [">", str(output_file)]:
                    command_without_redirect.append(item)
            
            # Ensure output directory exists
            output_dir.mkdir(parents=True, exist_ok=True)
            
            # Log execution start
            self.logger.info(f"Starting {self.tool_name} execution")
            
            # Display command parameters on screen
            cmd_display = " ".join(command_without_redirect)
            print(f"\n{'='*80}")
            print(f"🧬 MINIPROT COMMAND")
            print(f"{'='*80}")
            print(f"Command: {cmd_display}")
            print(f"Output:  {output_file}")
            print(f"{'='*80}\n")
            
            # Execute command with shell redirection
            with open(output_file, 'w') as f:
                result = subprocess.run(
                    command_without_redirect,
                    stdout=f,
                    stderr=subprocess.PIPE,
                    text=True,
                    timeout=timeout or self.config.get("timeout", 7200),
                    cwd=output_dir,
                )
            
            execution_time = (datetime.now() - start_time).total_seconds()
            
            # Apply quality filtering to remove frameshifts and premature stops
            final_output_file = output_file
            if result.returncode == 0 and output_file.exists():
                # Check if filtering is enabled in config
                enable_filtering = self.config.get("parameters", {}).get("enable_quality_filter", True)
                
                if enable_filtering:
                    self.logger.info("Applying quality filter to remove frameshifts and premature stops...")
                    try:
                        final_output_file = self._apply_quality_filter(output_file)
                    except Exception as e:
                        self.logger.warning(f"Quality filtering failed: {e}, using original results")
                        final_output_file = output_file
            
            # Get output files
            output_files = [final_output_file] if final_output_file.exists() else []
            
            # Create result object
            final_command = " ".join(command_without_redirect) + f" > {output_file}"
            if final_output_file != output_file:
                final_command += f" | quality_filter -> {final_output_file}"
                
            tool_result = ToolResult(
                tool_name=self.tool_name,
                command=final_command,
                return_code=result.returncode,
                stdout="",  # Redirected to file
                stderr=result.stderr,
                execution_time=execution_time,
                input_files=list(input_files.values()),
                output_files=output_files,
                success=result.returncode == 0,
                metadata={"config": self.config, "kwargs": kwargs},
            )
            
            # Log result
            if tool_result.success:
                self.logger.info(f"{self.tool_name} completed successfully in {execution_time:.2f}s")
            else:
                self.logger.error(f"{self.tool_name} failed with return code {result.returncode}")
                self.logger.error(f"STDERR: {result.stderr}")
            
            return tool_result
            
        except Exception as e:
            execution_time = (datetime.now() - start_time).total_seconds()
            self.logger.error(f"{self.tool_name} execution failed: {e}")
            raise
    
    def parse_output(self, output_dir: Path) -> Dict[str, Any]:
        """
        Parse miniprot output files
        
        Args:
            output_dir: Directory containing output files
            
        Returns:
            Parsed results dictionary
        """
        results = {
            "alignment_count": 0,
            "output_files": [],
            "alignments": [],
        }
        
        # Look for GFF output file
        gff_file = output_dir / "alignments.gff"
        if gff_file.exists():
            results["output_files"].append(str(gff_file))
            
            # Count alignments
            alignment_count = 0
            try:
                with open(gff_file, 'r') as f:
                    for line in f:
                        if not line.startswith('#') and line.strip():
                            fields = line.strip().split('\t')
                            if len(fields) >= 3:
                                feature_type = fields[2]
                                if feature_type in ['mRNA', 'transcript']:
                                    alignment_count += 1
                                    
                                    # Extract alignment info
                                    if len(fields) >= 9:
                                        alignment_info = {
                                            "seqid": fields[0],
                                            "source": fields[1],
                                            "type": fields[2],
                                            "start": int(fields[3]),
                                            "end": int(fields[4]),
                                            "score": fields[5],
                                            "strand": fields[6],
                                            "phase": fields[7],
                                            "attributes": fields[8],
                                        }
                                        results["alignments"].append(alignment_info)
                
                results["alignment_count"] = alignment_count
                
            except Exception as e:
                logger.warning(f"Error parsing miniprot output: {e}")
        
        return results
    
    def check_dependencies(self) -> Dict[str, bool]:
        """Check miniprot dependencies"""
        dependencies = super().check_dependencies()
        
        # Additional checks specific to miniprot
        try:
            import subprocess
            result = subprocess.run(
                [str(self.executable_path), "--version"],
                capture_output=True,
                text=True,
                timeout=10,
            )
            dependencies["version_check"] = result.returncode == 0
        except Exception:
            dependencies["version_check"] = False
        
        return dependencies