"""
EVidenceModeler (EVM) tool interface
EVidenceModeler工具接口

提供证据整合功能的工具接口
"""

from pathlib import Path
from typing import Any, Dict, List

from .base import ExternalTool, ToolResult
from ..utils.logging_setup import get_logger

logger = get_logger(__name__)


class EVMTool(ExternalTool):
    """
    Interface for EVidenceModeler tool
    EVidenceModeler工具接口
    
    EVM用于整合多种基因预测证据，生成最终的一致性注释
    """
    
    def __init__(self, config: Dict[str, Any] = None) -> None:
        """
        Initialize EVM tool interface
        
        Args:
            config: Tool configuration dictionary
        """
        self.config = config or {}
        
        # EVM is typically a collection of Perl scripts
        self.executable_dir = Path(self.config.get("executable_dir", "/home/wangys/data/work/nbs/tools/EVidenceModeler"))
        
        # Set up script paths
        self.scripts = {
            "main": self.executable_dir / "EVidenceModeler",
            "partition": self.executable_dir / "EvmUtils/partition_EVM_inputs.pl",
            "write_commands": self.executable_dir / "EvmUtils/write_EVM_commands.pl",
            "execute_commands": self.executable_dir / "EvmUtils/execute_EVM_commands.pl",
            "recombine": self.executable_dir / "EvmUtils/recombine_EVM_partial_outputs.pl",
            "convert": self.executable_dir / "EvmUtils/convert_EVM_outputs_to_GFF3.pl",
        }
        
        # Use perl as the executable
        super().__init__(
            tool_name="evm",
            executable_path=self._find_perl_executable(),
            config=config,
        )
    
    def _find_perl_executable(self) -> Path:
        """Find Perl executable"""
        import shutil
        perl_path = shutil.which("perl")
        if not perl_path:
            raise RuntimeError("Perl executable not found")
        return Path(perl_path)
    
    def prepare_input(self, **kwargs: Any) -> Dict[str, Path]:
        """
        Prepare input files for EVM - only integrating in NLR candidate regions
        为EVM准备输入文件 - 仅在NLR候选区域进行整合
        
        Args:
            prediction_files: List of prediction files to integrate
            genome_file: Path to genome FASTA file
            nlr_candidates_file: Path to NLR candidates file (optional, for region filtering)
            weights_file: Path to weights configuration file (optional)
            
        Returns:
            Dictionary of prepared input files
        """
        prediction_files = kwargs.get("prediction_files", [])
        genome_file = kwargs.get("genome_file")
        nlr_candidates_file = kwargs.get("nlr_candidates_file")
        
        if not prediction_files:
            raise ValueError("prediction_files list is required for EVM")
        if not genome_file:
            raise ValueError("genome_file is required for EVM")
        
        # Validate input files
        genome_path = Path(genome_file)
        self.validate_input(genome_file=genome_path)
        
        prediction_paths = []
        for pred_file in prediction_files:
            pred_path = Path(pred_file)
            self.validate_input(prediction_file=pred_path)
            prediction_paths.append(pred_path)
        
        input_files = {
            "genome": genome_path,
            "predictions": prediction_paths,
        }
        
        # Add NLR candidates file for region filtering
        if nlr_candidates_file:
            nlr_path = Path(nlr_candidates_file)
            if nlr_path.exists():
                input_files["nlr_candidates"] = nlr_path
                self.logger.info(f"Using NLR candidates for region filtering: {nlr_path}")
        
        # Handle weights file
        weights_file = kwargs.get("weights_file")
        if weights_file:
            weights_path = Path(weights_file)
            self.validate_input(weights_file=weights_path)
            input_files["weights"] = weights_path
        
        return input_files
    
    def create_weights_file(self, output_dir: Path) -> Path:
        """
        Create weights configuration file for EVM
        
        Args:
            output_dir: Output directory
            
        Returns:
            Path to created weights file
        """
        weights_file = output_dir / "weights.txt"
        
        # Set equal weights for Augustus and Miniprot predictions (both as gene predictions)
        weights = self.config.get("weights", {
            "ABINITIO_PREDICTION": 1,  # Augustus ab initio gene predictions
            "OTHER_PREDICTION": 1,     # Miniprot-based gene predictions
        })
        
        with open(weights_file, 'w') as f:
            f.write("# EVM weights configuration\n")
            f.write("# Equal weights for Augustus (ab initio) and Miniprot (evidence-based) gene predictions\n")
            f.write("ABINITIO_PREDICTION\tAUGUSTUS\t1\n")
            f.write("OTHER_PREDICTION\tminiprot\t1\n")
        
        return weights_file
    
    def _parse_nlr_regions(self, nlr_file: Path, extend_bp: int = 5000) -> List[Dict]:
        """
        Parse NLR candidate regions from GFF file and extend them for complete gene structures
        从GFF文件解析NLR候选区域并延长以包含完整基因结构
        
        Args:
            nlr_file: NLR candidates GFF file
            extend_bp: Number of base pairs to extend upstream and downstream
        """
        regions = []
        
        with open(nlr_file, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                    
                fields = line.strip().split('\t')
                if len(fields) >= 5:
                    original_start = int(fields[3])
                    original_end = int(fields[4])
                    
                    # Extend regions by 5kb upstream and downstream
                    extended_start = max(1, original_start - extend_bp)
                    extended_end = original_end + extend_bp
                    
                    region = {
                        'seqid': fields[0],
                        'start': extended_start,
                        'end': extended_end,
                        'strand': fields[6] if len(fields) > 6 else '.',
                        'original_start': original_start,
                        'original_end': original_end,
                    }
                    regions.append(region)
        
        self.logger.info(f"Parsed {len(regions)} NLR candidate regions, extended by ±{extend_bp}bp for complete gene structures")
        return regions

    def _filter_predictions_by_nlr_regions(self, pred_file: Path, nlr_regions: List[Dict], 
                                          output_dir: Path, file_suffix: str) -> Path:
        """
        Filter prediction file to only include features in NLR regions
        过滤预测文件，仅包含NLR区域内的特征
        """
        filtered_file = output_dir / f"nlr_filtered_{file_suffix}.gff"
        
        features_in_regions = 0
        total_features = 0
        
        with open(filtered_file, 'w') as out_f:
            out_f.write("##gff-version 3\n")
            out_f.write(f"# Filtered predictions in NLR candidate regions from {pred_file.name}\n")
            
            with open(pred_file, 'r') as in_f:
                for line in in_f:
                    if line.startswith('#'):
                        continue
                        
                    line = line.strip()
                    if not line:
                        continue
                        
                    fields = line.split('\t')
                    if len(fields) >= 5:
                        total_features += 1
                        seqid = fields[0]
                        start = int(fields[3])
                        end = int(fields[4])
                        
                        # Check if feature overlaps with any NLR region
                        for region in nlr_regions:
                            if (seqid == region['seqid'] and 
                                not (end < region['start'] or start > region['end'])):
                                # Feature overlaps with NLR region
                                out_f.write(line + '\n')
                                features_in_regions += 1
                                break
        
        self.logger.info(f"Filtered {file_suffix}: {features_in_regions}/{total_features} features in NLR regions")
        return filtered_file

    def create_evm_inputs(self, input_files: Dict[str, Path], output_dir: Path) -> Dict[str, Path]:
        """
        Create EVM-specific input files, filtered to NLR candidate regions
        创建EVM专用输入文件，过滤到NLR候选区域
        """
        evm_inputs = {}
        
        # Genome file
        evm_inputs["genome"] = input_files["genome"]
        
        # Create weights file
        weights_file = self.create_weights_file(output_dir)
        evm_inputs["weights"] = weights_file
        
        # Parse NLR regions if available - extend regions significantly for miniprot integration
        nlr_regions = []
        if "nlr_candidates" in input_files:
            nlr_regions = self._parse_nlr_regions(input_files["nlr_candidates"], extend_bp=50000)  # Extend to 50kb for better coverage
        
        # Separate and filter Augustus and Miniprot files (both as gene predictions)
        augustus_file = None
        miniprot_file = None
        
        for pred_file in input_files.get("predictions", []):
            file_path_str = str(pred_file).lower()
            if ("augustus" in file_path_str or "predictions" in file_path_str or 
                "gene_prediction" in file_path_str):
                augustus_file = pred_file
            elif ("miniprot" in file_path_str or "alignments" in file_path_str or 
                  "protein_alignment" in file_path_str):
                miniprot_file = pred_file
        
        # Filter files to NLR regions if NLR regions are available
        if nlr_regions:
            gene_prediction_files = []
            
            if augustus_file:
                filtered_augustus = self._filter_predictions_by_nlr_regions(
                    augustus_file, nlr_regions, output_dir, "augustus"
                )
                gene_prediction_files.append(filtered_augustus)
            
            if miniprot_file:
                # Convert Miniprot to gene prediction format and filter
                converted_miniprot = self._convert_miniprot_to_gene_predictions(
                    miniprot_file, output_dir
                )
                filtered_miniprot = self._filter_predictions_by_nlr_regions(
                    converted_miniprot, nlr_regions, output_dir, "miniprot"
                )
                gene_prediction_files.append(filtered_miniprot)
            
            # Combine multiple gene prediction files if needed
            if len(gene_prediction_files) == 1:
                evm_inputs["gene_predictions"] = gene_prediction_files[0]
            elif len(gene_prediction_files) > 1:
                combined_file = self._combine_gene_prediction_files(gene_prediction_files, output_dir)
                evm_inputs["gene_predictions"] = combined_file
        else:
            # No NLR filtering - use original files
            self.logger.warning("No NLR candidates provided - using full genome predictions")
            gene_prediction_files = []
            
            if augustus_file:
                gene_prediction_files.append(augustus_file)
            
            if miniprot_file:
                converted_miniprot = self._convert_miniprot_to_gene_predictions(
                    miniprot_file, output_dir
                )
                gene_prediction_files.append(converted_miniprot)
            
            # Combine files
            if len(gene_prediction_files) == 1:
                evm_inputs["gene_predictions"] = gene_prediction_files[0]
            elif len(gene_prediction_files) > 1:
                combined_file = self._combine_gene_prediction_files(gene_prediction_files, output_dir)
                evm_inputs["gene_predictions"] = combined_file
            
        return evm_inputs

    def _convert_miniprot_to_gene_predictions(self, miniprot_file: Path, output_dir: Path) -> Path:
        """
        Convert Miniprot alignment output to gene prediction format for EVM
        将Miniprot比对输出转换为EVM的基因预测格式
        
        This now prioritizes filtered/processed miniprot results over raw output
        """
        converted_file = output_dir / "miniprot_gene_predictions.gff"
        
        # Check if we have processed miniprot results to use instead
        processed_miniprot_dir = Path("results/miniprot_processed")
        if processed_miniprot_dir.exists():
            self.logger.info("Found processed miniprot directory, using filtered results")
            return self._convert_processed_miniprot_to_gene_predictions(processed_miniprot_dir, converted_file)
        
        # Fallback to converting raw miniprot file
        self.logger.warning("No processed miniprot found, using raw miniprot output")
        
        with open(converted_file, 'w') as out_f:
            out_f.write("##gff-version 3\n")
            out_f.write("# Miniprot alignments converted to gene predictions for EVM\n")
            
            with open(miniprot_file, 'r') as in_f:
                for line in in_f:
                    if line.startswith('#') or not line.strip():
                        continue
                        
                    fields = line.strip().split('\t')
                    if len(fields) >= 9:
                        # Filter out low quality predictions here if no processed results
                        attrs = fields[8]
                        
                        # Skip predictions with frameshifts or stop codons for raw data
                        if "Frameshift=" in attrs or "StopCodon=" in attrs:
                            continue
                            
                        # Skip very low identity matches
                        if "Identity=" in attrs:
                            try:
                                identity_str = [a for a in attrs.split(';') if a.startswith('Identity=')][0]
                                identity = float(identity_str.split('=')[1])
                                if identity < 0.3:  # Skip < 30% identity
                                    continue
                            except:
                                pass
                        
                        # Change source to indicate gene prediction origin
                        fields[1] = "miniprot"
                        
                        # Keep all feature types (gene, mRNA, CDS, exon)
                        # EVM will understand these as gene predictions
                        out_f.write('\t'.join(fields) + '\n')
        
        self.logger.info(f"Converted Miniprot alignments to gene predictions: {converted_file}")
        return converted_file
    
    def _convert_processed_miniprot_to_gene_predictions(self, processed_dir: Path, output_file: Path) -> Path:
        """
        Convert processed/filtered miniprot results to gene predictions for EVM
        """
        self.logger.info("Using processed miniprot results for EVM integration")
        
        # Use high and medium quality predictions only
        quality_files = [
            processed_dir / "miniprot_high_quality.gff3",
            processed_dir / "miniprot_medium_quality.gff3"
        ]
        
        with open(output_file, 'w') as out_f:
            out_f.write("##gff-version 3\n")
            out_f.write("# Processed Miniprot predictions converted to gene predictions for EVM\n")
            out_f.write("# Using high and medium quality predictions only\n")
            
            for quality_file in quality_files:
                if quality_file.exists():
                    self.logger.info(f"Adding predictions from: {quality_file}")
                    with open(quality_file, 'r') as in_f:
                        for line in in_f:
                            if line.startswith('#') or not line.strip():
                                continue
                                
                            fields = line.strip().split('\t')
                            if len(fields) >= 9:
                                # Change source to indicate gene prediction origin
                                fields[1] = "miniprot"
                                out_f.write('\t'.join(fields) + '\n')
                else:
                    self.logger.warning(f"Quality file not found: {quality_file}")
        
        self.logger.info(f"Converted processed Miniprot to gene predictions: {output_file}")
        return output_file

    def _combine_gene_prediction_files(self, gene_files: List[Path], output_dir: Path) -> Path:
        """
        Combine multiple gene prediction files into a single file for EVM
        将多个基因预测文件合并为EVM的单个文件
        """
        combined_file = output_dir / "combined_gene_predictions.gff"
        
        with open(combined_file, 'w') as out_f:
            out_f.write("##gff-version 3\n")
            out_f.write("# Combined gene predictions from multiple sources for EVM\n")
            
            for gene_file in gene_files:
                if gene_file.exists():
                    source_name = gene_file.stem
                    out_f.write(f"# === Predictions from {source_name} ===\n")
                    
                    with open(gene_file, 'r') as in_f:
                        for line in in_f:
                            if not line.startswith('#'):
                                out_f.write(line)
        
        self.logger.info(f"Combined {len(gene_files)} gene prediction files: {combined_file}")
        return combined_file

    def build_evm_command_pipeline(
        self,
        evm_inputs: Dict[str, Path],
        output_dir: Path,
        **kwargs: Any,
    ) -> List[str]:
        """
        Build complete EVM command pipeline
        构建完整的EVM命令流水线
        
        Returns:
            List of shell commands for EVM execution
        """
        genome_file = evm_inputs["genome"]
        weights_file = evm_inputs["weights"]
        
        # Parameters
        parameters = self.config.get("parameters", {})
        segment_size = parameters.get("segmentSize", 100000)
        overlap_size = parameters.get("overlapSize", 10000)
        
        commands = []
        
        # Step 1: Partition inputs
        partition_cmd = [
            "perl", str(self.scripts["partition"]),
            "--genome", str(genome_file),
            "--segmentSize", str(segment_size),
            "--overlapSize", str(overlap_size),
            "--partition_listing", str(output_dir / "partitions_list.out")
        ]
        
        # Add gene predictions (now includes both Augustus and Miniprot)
        if "gene_predictions" in evm_inputs:
            partition_cmd.extend(["--gene_predictions", str(evm_inputs["gene_predictions"])])
        
        commands.append(" ".join(partition_cmd))
        
        # Step 2: Write EVM commands
        write_cmd = [
            "perl", str(self.scripts["write_commands"]),
            "--genome", str(genome_file),
            "--weights", str(weights_file),
            "--partition_listing", str(output_dir / "partitions_list.out"),
            "--output_file_name", "evm.out"
        ]
        
        # Add gene predictions
        if "gene_predictions" in evm_inputs:
            write_cmd.extend(["--gene_predictions", str(evm_inputs["gene_predictions"])])
            
        commands.append(" ".join(write_cmd))
        
        # Step 3: Execute EVM commands (simplified single-threaded)
        execute_cmd = [
            "perl", str(self.scripts["execute_commands"]),
            str(output_dir / "partitions_list.out")
        ]
        commands.append(" ".join(execute_cmd))
        
        # Step 4: Recombine outputs
        recombine_cmd = [
            "perl", str(self.scripts["recombine"]),
            "--partition_listing", str(output_dir / "partitions_list.out"),
            "--output_file_name", "evm.out"
        ]
        commands.append(" ".join(recombine_cmd))
        
        # Step 5: Convert to GFF3
        convert_cmd = [
            "perl", str(self.scripts["convert"]),
            "--partition_listing", str(output_dir / "partitions_list.out"),
            "--output_file_name", "evm.out",
            "--genome", str(genome_file)
        ]
        commands.append(" ".join(convert_cmd))
        
        return commands

    def build_simple_evm_command(
        self,
        evm_inputs: Dict[str, Path],
        output_dir: Path,
        **kwargs: Any,
    ) -> str:
        """
        Build simple single-command EVM execution
        构建简单的单命令EVM执行
        
        Returns:
            Shell command string for EVM execution
        """
        genome_file = evm_inputs.get("genome")
        weights_file = evm_inputs.get("weights")
        
        if not genome_file:
            raise ValueError(f"Missing genome file in evm_inputs: {evm_inputs}")
        if not weights_file:
            raise ValueError(f"Missing weights file in evm_inputs: {evm_inputs}")
        
        # Build EVM command
        cmd_parts = [str(self.scripts["main"])]
        cmd_parts.extend(["--genome", str(Path(genome_file).resolve())])
        cmd_parts.extend(["--weights", str(Path(weights_file).resolve())])
        
        # Add required sample_id parameter
        sample_id = kwargs.get("sample_id", "sample1")
        cmd_parts.extend(["--sample_id", str(sample_id)])
        
        # Add gene predictions (includes both Augustus and Miniprot as gene predictions)
        if "gene_predictions" in evm_inputs:
            cmd_parts.extend(["--gene_predictions", str(Path(evm_inputs["gene_predictions"]).resolve())])
        
        # Parameters
        parameters = self.config.get("parameters", {})
        segment_size = parameters.get("segmentSize", 100000)
        overlap_size = parameters.get("overlapSize", 10000)
        
        cmd_parts.extend(["--segmentSize", str(segment_size)])
        cmd_parts.extend(["--overlapSize", str(overlap_size)])
        
        # Output
        cmd_parts.extend([">", str(output_dir / "evm.out")])
        
        return " ".join(cmd_parts)

    def build_command(
        self,
        input_files: Dict[str, Path],
        output_dir: Path,
        **kwargs: Any,
    ) -> List[str]:
        """
        Build EVM command for execution
        
        Args:
            input_files: Prepared input files
            output_dir: Output directory
            **kwargs: Additional parameters
            
        Returns:
            Command as list of strings
        """
        # Create EVM-specific inputs
        evm_inputs = self.create_evm_inputs(input_files, output_dir)
        
        # Try simple command first
        cmd_string = self.build_simple_evm_command(evm_inputs, output_dir, **kwargs)
        
        # Return as shell command
        return ["bash", "-c", cmd_string]
    
    def execute(self, output_dir: Path, timeout: int = None, **kwargs: Any) -> ToolResult:
        """
        Execute EVM pipeline using real EVM tool
        
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
            output_dir.mkdir(parents=True, exist_ok=True)
            
            # Check EVM availability
            if not self._check_evm_availability():
                raise RuntimeError(f"EVM not found at: {self.executable_dir}")
            
            # Execute real EVM
            return self._execute_real_evm(input_files, output_dir, timeout, **kwargs)
            
        except Exception as e:
            execution_time = (datetime.now() - start_time).total_seconds()
            self.logger.error(f"EVM execution failed: {e}")
            raise

    def _check_evm_availability(self) -> bool:
        """Check if real EVM is available"""
        evm_main = self.scripts["main"]
        partition_script = self.scripts["partition"]
        
        available = (evm_main.exists() and 
                    partition_script.exists() and 
                    self.executable_path.exists())
        
        if available:
            self.logger.info(f"EVM found at: {self.executable_dir}")
        else:
            self.logger.info(f"EVM not found at: {self.executable_dir}")
            self.logger.info(f"Main script exists: {evm_main.exists()}")
            self.logger.info(f"Partition script exists: {partition_script.exists()}")
            
        return available

    def _execute_real_evm(self, input_files: Dict[str, Path], output_dir: Path, 
                         timeout: int = None, **kwargs: Any) -> ToolResult:
        """Execute real EVM command"""
        import subprocess
        from datetime import datetime
        
        start_time = datetime.now()
        
        # Build command
        command_list = self.build_command(input_files, output_dir, **kwargs)
        
        self.logger.info(f"Executing real EVM command: {' '.join(command_list)}")
        
        # Execute command
        result = subprocess.run(
            command_list,
            cwd=output_dir,
            capture_output=True,
            text=True,
            timeout=timeout or self.config.get("timeout", 3600),
        )
        
        execution_time = (datetime.now() - start_time).total_seconds()
        
        # Check for output files
        output_file = output_dir / "evm.out"
        gff_output = output_dir / "final_annotations.gff"
        
        # Convert EVM output to standard GFF if needed
        if output_file.exists() and not gff_output.exists():
            self._convert_evm_output_to_gff(output_file, gff_output)
        
        output_files = [f for f in [output_file, gff_output] if f.exists()]
        
        tool_result = ToolResult(
            tool_name=self.tool_name,
            command=" ".join(command_list),
            return_code=result.returncode,
            stdout=result.stdout,
            stderr=result.stderr,
            execution_time=execution_time,
            input_files=input_files.get("predictions", []),
            output_files=output_files,
            success=result.returncode == 0,
            metadata={"config": self.config, "kwargs": kwargs, "mode": "real_evm"},
        )
        
        if tool_result.success:
            self.logger.info(f"Real EVM completed successfully in {execution_time:.2f}s")
        else:
            self.logger.error(f"Real EVM failed with return code {result.returncode}")
            self.logger.error(f"STDERR: {result.stderr}")
            raise RuntimeError(f"EVM execution failed: {result.stderr}")
        
        return tool_result

    def _convert_evm_output_to_gff(self, evm_out: Path, gff_out: Path) -> None:
        """Convert EVM output to standard GFF3 format"""
        self.logger.info(f"Converting EVM output to GFF3: {evm_out} -> {gff_out}")
        
        # Check for the proper EVM GFF3 output first
        output_dir = evm_out.parent
        evm_gff3_file = None
        
        # Look for pattern: sample*.EVM.gff3
        for gff3_file in output_dir.glob("*.EVM.gff3"):
            evm_gff3_file = gff3_file
            break
            
        if evm_gff3_file and evm_gff3_file.exists():
            self.logger.info(f"Found EVM GFF3 output: {evm_gff3_file}")
            # Copy the proper GFF3 file
            import shutil
            shutil.copy2(evm_gff3_file, gff_out)
            self.logger.info(f"Copied EVM GFF3 output to: {gff_out}")
            return
        
        # Fallback to converting evm.out if no GFF3 found
        self.logger.warning(f"No EVM GFF3 file found, trying to convert: {evm_out}")
        with open(gff_out, 'w') as out_f:
            out_f.write("##gff-version 3\n")
            out_f.write("# EVM output converted to GFF3\n")
            
            if evm_out.exists():
                with open(evm_out, 'r') as in_f:
                    for line in in_f:
                        line = line.strip()
                        if line and not line.startswith('#') and '\t' in line:
                            # Only write lines that look like GFF format
                            fields = line.split('\t')
                            if len(fields) >= 8:
                                out_f.write(line + '\n')
    
    def filter_miniprot_to_genes(self, miniprot_file: Path, output_dir: Path) -> Path:
        """
        Filter Miniprot results to gene-level predictions for EVM
        过滤Miniprot结果为基因级别预测用于EVM整合
        
        Args:
            miniprot_file: Original Miniprot GFF file
            output_dir: Output directory
            
        Returns:
            Path to filtered GFF file
        """
        filtered_file = output_dir / "miniprot_genes_only.gff"
        
        with open(filtered_file, 'w') as out_f:
            out_f.write("##gff-version 3\n")
            out_f.write("# Miniprot protein alignments filtered to gene predictions\n")
            
            if miniprot_file.exists():
                with open(miniprot_file, 'r') as in_f:
                    current_gene = []
                    for line in in_f:
                        if line.startswith('#'):
                            continue
                            
                        line = line.strip()
                        if not line:
                            continue
                            
                        fields = line.split('\t')
                        if len(fields) >= 9:
                            feature_type = fields[2]
                            
                            # Convert CDS/exon to gene predictions
                            if feature_type in ['CDS', 'exon']:
                                # Change source to PROTEIN for EVM
                                fields[1] = "PROTEIN"
                                fields[2] = "gene"  # Convert to gene feature
                                
                                # Update attributes
                                attrs = fields[8]
                                if "ID=" not in attrs:
                                    attrs = f"ID=miniprot_gene_{fields[0]}_{fields[3]}_{fields[4]};{attrs}"
                                    fields[8] = attrs
                                
                                out_f.write('\t'.join(fields) + '\n')
        
        return filtered_file

    def _parse_gff_features(self, gff_file: Path) -> Dict[str, List[Dict]]:
        """
        Parse GFF file and group features by gene with proper hierarchy handling
        解析GFF文件并按基因分组特征，正确处理层级关系
        """
        genes = {}
        id_to_gene_map = {}  # Map feature IDs to their gene group
        
        # First pass: identify all features and their relationships
        all_features = []
        with open(gff_file, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                    
                fields = line.strip().split('\t')
                if len(fields) >= 9:
                    feature_type = fields[2]
                    attrs = fields[8]
                    
                    # Extract ID and Parent from attributes
                    feature_id = None
                    parent_id = None
                    
                    for attr in attrs.split(';'):
                        attr = attr.strip()
                        if attr.startswith('ID='):
                            feature_id = attr.split('=')[1]
                        elif attr.startswith('Parent='):
                            parent_id = attr.split('=')[1]
                    
                    feature = {
                        'seqid': fields[0],
                        'source': fields[1],
                        'type': fields[2],
                        'start': int(fields[3]),
                        'end': int(fields[4]),
                        'score': fields[5],
                        'strand': fields[6],
                        'phase': fields[7],
                        'attributes': fields[8],
                        'feature_id': feature_id,
                        'parent_id': parent_id,
                        'original_line': line.strip()
                    }
                    all_features.append(feature)
        
        # Second pass: determine gene groupings
        for feature in all_features:
            gene_id = None
            
            if feature['type'] == 'gene':
                # This is a gene feature
                gene_id = feature['feature_id'] or f"gene_{feature['seqid']}_{feature['start']}_{feature['end']}"
                id_to_gene_map[feature['feature_id']] = gene_id
            else:
                # Non-gene feature - find its gene
                if feature['parent_id']:
                    # Trace parent hierarchy to find gene
                    current_parent = feature['parent_id']
                    while current_parent and current_parent not in id_to_gene_map:
                        # Look for parent feature
                        parent_feature = next((f for f in all_features if f['feature_id'] == current_parent), None)
                        if parent_feature:
                            if parent_feature['type'] == 'gene':
                                gene_id = parent_feature['feature_id'] or f"gene_{parent_feature['seqid']}_{parent_feature['start']}_{parent_feature['end']}"
                                id_to_gene_map[current_parent] = gene_id
                                break
                            else:
                                current_parent = parent_feature['parent_id']
                        else:
                            break
                    
                    gene_id = id_to_gene_map.get(current_parent)
                
                # For features like mRNA that have no parent or no gene parent, treat each as a separate gene
                if not gene_id:
                    if feature['type'] in ['mRNA', 'transcript']:
                        # Use the mRNA/transcript ID as the gene ID
                        gene_id = feature['feature_id'] or f"gene_{feature['seqid']}_{feature['start']}_{feature['end']}"
                    else:
                        # For other features without clear gene grouping, create gene ID based on coordinates
                        gene_id = f"gene_{feature['seqid']}_{feature['start']}_{feature['end']}"
            
            # Add feature to gene group
            if gene_id:
                feature['gene_id'] = gene_id
                if gene_id not in genes:
                    genes[gene_id] = []
                genes[gene_id].append(feature)
        
        # Log feature type statistics
        total_features_by_type = {}
        for gene_id, features in genes.items():
            feature_types = [f['type'] for f in features]
            for ftype in feature_types:
                total_features_by_type[ftype] = total_features_by_type.get(ftype, 0) + 1
            self.logger.debug(f"Gene {gene_id}: {len(features)} features - {set(feature_types)}")
        
        self.logger.info(f"Parsed GFF features by type: {total_features_by_type}")
        return genes

    def _merge_overlapping_genes(self, augustus_genes: Dict, miniprot_genes: Dict) -> List[Dict]:
        """
        Merge overlapping genes from Augustus and Miniprot
        合并Augustus和Miniprot重叠的基因
        """
        merged_genes = []
        used_augustus = set()
        used_miniprot = set()
        
        # Find overlapping genes
        for aug_id, aug_features in augustus_genes.items():
            if aug_id in used_augustus:
                continue
                
            aug_gene = next((f for f in aug_features if f['type'] == 'gene'), None)
            if not aug_gene:
                continue
                
            best_overlap = None
            best_overlap_ratio = 0
            
            for mini_id, mini_features in miniprot_genes.items():
                if mini_id in used_miniprot:
                    continue
                    
                mini_gene = next((f for f in mini_features if f['type'] == 'gene'), None)
                if not mini_gene:
                    continue
                
                # Check if genes overlap
                if (aug_gene['seqid'] == mini_gene['seqid'] and 
                    aug_gene['strand'] == mini_gene['strand']):
                    
                    overlap_start = max(aug_gene['start'], mini_gene['start'])
                    overlap_end = min(aug_gene['end'], mini_gene['end'])
                    
                    if overlap_start <= overlap_end:
                        overlap_len = overlap_end - overlap_start + 1
                        aug_len = aug_gene['end'] - aug_gene['start'] + 1
                        mini_len = mini_gene['end'] - mini_gene['start'] + 1
                        
                        overlap_ratio = overlap_len / min(aug_len, mini_len)
                        
                        if overlap_ratio > 0.3 and overlap_ratio > best_overlap_ratio:
                            best_overlap = mini_id
                            best_overlap_ratio = overlap_ratio
            
            # Merge genes if overlap found
            if best_overlap:
                merged_gene = self._merge_gene_structures(
                    aug_features, miniprot_genes[best_overlap]
                )
                merged_genes.append(merged_gene)
                used_augustus.add(aug_id)
                used_miniprot.add(best_overlap)
            else:
                # Keep Augustus gene alone
                merged_genes.append({'augustus': aug_features, 'miniprot': []})
                used_augustus.add(aug_id)
        
        # Add remaining Miniprot genes
        for mini_id, mini_features in miniprot_genes.items():
            if mini_id not in used_miniprot:
                merged_genes.append({'augustus': [], 'miniprot': mini_features})
        
        return merged_genes

    def _merge_gene_structures(self, augustus_features: List[Dict], miniprot_features: List[Dict]) -> Dict:
        """
        Merge gene structures from Augustus and Miniprot
        合并Augustus和Miniprot的基因结构
        """
        aug_gene = next((f for f in augustus_features if f['type'] == 'gene'), None)
        mini_gene = next((f for f in miniprot_features if f['type'] == 'gene'), None)
        
        # Use Augustus structure as primary, supplement with Miniprot
        if aug_gene:
            # Extend gene boundaries if Miniprot provides wider coverage
            if mini_gene:
                start = min(aug_gene['start'], mini_gene['start'])
                end = max(aug_gene['end'], mini_gene['end'])
            else:
                start = aug_gene['start']
                end = aug_gene['end']
            
            return {
                'augustus': augustus_features,
                'miniprot': miniprot_features,
                'merged_start': start,
                'merged_end': end,
                'primary_source': 'augustus'
            }
        elif mini_gene:
            return {
                'augustus': [],
                'miniprot': miniprot_features,
                'merged_start': mini_gene['start'],
                'merged_end': mini_gene['end'],
                'primary_source': 'miniprot'
            }
        else:
            return {'augustus': augustus_features, 'miniprot': miniprot_features}

    def _filter_genes_by_nlr_regions(self, genes_dict: Dict[str, List[Dict]], 
                                    nlr_regions: List[Dict]) -> Dict[str, List[Dict]]:
        """
        Filter genes dictionary to only include genes in NLR regions
        过滤基因字典，仅包含NLR区域内的基因
        """
        filtered_genes = {}
        
        for gene_id, features in genes_dict.items():
            # Find gene feature to get coordinates
            gene_feature = next((f for f in features if f['type'] == 'gene'), None)
            
            # If no explicit gene feature, use the representative feature (e.g., mRNA or first feature)
            if not gene_feature:
                # Try to find mRNA/transcript feature first
                gene_feature = next((f for f in features if f['type'] in ['mRNA', 'transcript']), None)
                # If still no representative feature, use the first feature
                if not gene_feature and features:
                    gene_feature = features[0]
            
            if not gene_feature:
                continue
                
            # Check if gene overlaps with any NLR region
            for region in nlr_regions:
                if (gene_feature['seqid'] == region['seqid'] and 
                    not (gene_feature['end'] < region['start'] or 
                         gene_feature['start'] > region['end'])):
                    # Gene overlaps with NLR region
                    filtered_genes[gene_id] = features
                    break
        
        return filtered_genes

    def _create_integrated_evm_output(self, input_files: Dict[str, Path], output_file: Path) -> None:
        """
        Create integrated EVM output with proper gene structure merging
        创建带有正确基因结构合并的EVM整合输出
        
        Args:
            input_files: Input prediction files
            output_file: Output GFF file
        """
        prediction_files = input_files.get("predictions", [])
        
        augustus_genes = {}
        miniprot_genes = {}
        
        # Parse NLR regions for filtering - use same extension as main pipeline
        nlr_regions = []
        if "nlr_candidates" in input_files:
            nlr_regions = self._parse_nlr_regions(input_files["nlr_candidates"], extend_bp=50000)  # Use same extension as main pipeline
        
        # Parse input files and filter by NLR regions
        for pred_file in prediction_files:
            if not pred_file.exists():
                self.logger.warning(f"Prediction file not found: {pred_file}")
                continue
            
            file_path_str = str(pred_file).lower()
            self.logger.info(f"Processing prediction file: {pred_file}")
            
            # Check file content to determine type (more reliable than path-based detection)
            is_augustus_file = False
            is_miniprot_file = False
            
            try:
                with open(pred_file, 'r') as f:
                    # Read first few lines to check source
                    for _ in range(10):  # Check first 10 lines
                        line = f.readline()
                        if not line:
                            break
                        if line.startswith('#'):
                            continue
                        fields = line.strip().split('\t')
                        if len(fields) >= 2:
                            source = fields[1].lower()
                            if source == "augustus":
                                is_augustus_file = True
                                break
                            elif source == "miniprot":
                                is_miniprot_file = True
                                break
            except Exception as e:
                self.logger.warning(f"Error reading file {pred_file}: {e}")
                
            # Fallback to path-based detection if content detection fails
            if not is_augustus_file and not is_miniprot_file:
                if ("augustus" in file_path_str or "predictions" in file_path_str or 
                    "gene_prediction" in file_path_str):
                    is_augustus_file = True
                elif ("miniprot" in file_path_str or "alignments" in file_path_str or 
                      "protein_alignment" in file_path_str):
                    is_miniprot_file = True
            
            if is_augustus_file:
                self.logger.info(f"Identified as Augustus file: {pred_file}")
                all_augustus_genes = self._parse_gff_features(pred_file)
                # Filter to NLR regions only
                if nlr_regions:
                    augustus_genes = self._filter_genes_by_nlr_regions(all_augustus_genes, nlr_regions)
                    self.logger.info(f"Filtered Augustus genes: {len(augustus_genes)}/{len(all_augustus_genes)} in NLR regions")
                else:
                    augustus_genes = all_augustus_genes
                    self.logger.info(f"Parsed {len(augustus_genes)} Augustus genes (no NLR filtering)")
                    
            elif is_miniprot_file:
                self.logger.info(f"Identified as Miniprot file: {pred_file}")
                all_miniprot_genes = self._parse_gff_features(pred_file)
                # Filter to NLR regions only
                if nlr_regions:
                    miniprot_genes = self._filter_genes_by_nlr_regions(all_miniprot_genes, nlr_regions)
                    self.logger.info(f"Filtered Miniprot genes: {len(miniprot_genes)}/{len(all_miniprot_genes)} in NLR regions")
                else:
                    miniprot_genes = all_miniprot_genes
                    self.logger.info(f"Parsed {len(miniprot_genes)} Miniprot genes (no NLR filtering)")
            else:
                self.logger.warning(f"Unknown file type: {pred_file}")
        
        self.logger.info(f"Augustus genes: {len(augustus_genes)}, Miniprot genes: {len(miniprot_genes)}")
        
        # Merge overlapping genes
        merged_genes = self._merge_overlapping_genes(augustus_genes, miniprot_genes)
        self.logger.info(f"After merging: {len(merged_genes)} gene structures")
        
        # Write output
        with open(output_file, 'w') as out_f:
            out_f.write("##gff-version 3\n")
            out_f.write("# EVM integrated annotations\n")
            out_f.write("# Merged Augustus and Miniprot predictions\n")
            out_f.write(f"# Input: {len(augustus_genes)} Augustus + {len(miniprot_genes)} Miniprot genes\n")
            out_f.write(f"# Output: {len(merged_genes)} integrated gene structures\n")
            
            gene_counter = 1
            written_genes = 0
            written_features_by_type = {}
            
            for merged_gene in merged_genes:
                augustus_features = merged_gene.get('augustus', [])
                miniprot_features = merged_gene.get('miniprot', [])
                primary_source = merged_gene.get('primary_source', 'augustus')
                
                # Determine which features to write and integration strategy
                if augustus_features and miniprot_features:
                    # Both sources available - true integration
                    if primary_source == 'augustus':
                        features_to_write = augustus_features
                        evidence_type = "INTEGRATED_ABINITIO_PRIMARY"
                    else:
                        features_to_write = miniprot_features  
                        evidence_type = "INTEGRATED_OTHER_PRIMARY"
                elif augustus_features:
                    # Augustus only
                    features_to_write = augustus_features
                    evidence_type = "ABINITIO_PREDICTION"
                elif miniprot_features:
                    # Miniprot only (as gene prediction)
                    features_to_write = miniprot_features
                    evidence_type = "OTHER_PREDICTION"
                else:
                    continue
                
                # Sort features by hierarchy: gene → mRNA/transcript → CDS/exon
                feature_order = {'gene': 0, 'mRNA': 1, 'transcript': 1, 'CDS': 2, 'exon': 2, 'five_prime_UTR': 2, 'three_prime_UTR': 2}
                sorted_features = sorted(features_to_write, key=lambda f: (feature_order.get(f['type'], 99), f['start']))
                
                # Write all features (gene, mRNA, CDS, exon) in correct order
                gene_written = False
                current_gene_id = f"evm_gene_{gene_counter}"
                current_mrna_id = f"evm_mRNA_{gene_counter}"
                
                for feature in sorted_features:
                    fields = [
                        feature['seqid'],
                        'EVM',  # Change source to EVM
                        feature['type'],
                        str(feature['start']),
                        str(feature['end']),
                        feature['score'],
                        feature['strand'],
                        feature['phase'],
                        feature['attributes']
                    ]
                    
                    # Update attributes for EVM with proper hierarchy
                    if feature['type'] == 'gene':
                        attrs = f"ID={current_gene_id};evidence_type={evidence_type}"
                        if augustus_features and miniprot_features:
                            attrs += ";merged=true;augustus_support=true;miniprot_support=true"
                        elif augustus_features:
                            attrs += ";augustus_support=true"
                        elif miniprot_features:
                            attrs += ";miniprot_support=true"
                        fields[8] = attrs
                        gene_written = True
                        
                    elif feature['type'] in ['mRNA', 'transcript']:
                        # Update mRNA parent reference
                        fields[8] = f"ID={current_mrna_id};Parent={current_gene_id}"
                        
                    elif feature['type'] in ['CDS', 'exon', 'five_prime_UTR', 'three_prime_UTR']:
                        # Update CDS/exon parent reference
                        original_attrs = feature['attributes']
                        # Preserve original ID if exists, but update Parent
                        original_id = None
                        for attr in original_attrs.split(';'):
                            if attr.startswith('ID='):
                                original_id = attr.split('=')[1]
                                break
                        
                        if original_id:
                            fields[8] = f"ID={original_id};Parent={current_mrna_id}"
                        else:
                            fields[8] = f"Parent={current_mrna_id}"
                    
                    out_f.write('\t'.join(fields) + '\n')
                    
                    # Track written feature types
                    ftype = feature['type']
                    written_features_by_type[ftype] = written_features_by_type.get(ftype, 0) + 1
                
                if gene_written:
                    written_genes += 1
                    gene_counter += 1
            
            self.logger.info(f"Wrote {written_genes} integrated genes to output")
            self.logger.info(f"Written features by type: {written_features_by_type}")
        
        self.logger.info(f"EVM integration completed: {len(merged_genes)} integrated gene structures")
    
    def parse_output(self, output_dir: Path) -> Dict[str, Any]:
        """
        Parse EVM output files
        
        Args:
            output_dir: Directory containing output files
            
        Returns:
            Parsed results dictionary
        """
        results = {
            "final_gene_count": 0,
            "output_files": [],
            "final_genes": [],
        }
        
        # Look for final annotations
        gff_file = output_dir / "final_annotations.gff"
        if gff_file.exists():
            results["output_files"].append(str(gff_file))
            
            # Count final genes
            gene_count = 0
            try:
                with open(gff_file, 'r') as f:
                    for line in f:
                        if not line.startswith('#') and line.strip():
                            fields = line.strip().split('\t')
                            if len(fields) >= 3:
                                feature_type = fields[2]
                                if feature_type == 'gene':
                                    gene_count += 1
                                    
                                    # Extract gene info
                                    if len(fields) >= 9:
                                        gene_info = {
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
                                        results["final_genes"].append(gene_info)
                
                results["final_gene_count"] = gene_count
                
            except Exception as e:
                logger.warning(f"Error parsing EVM output: {e}")
        
        return results
    
    def check_dependencies(self) -> Dict[str, bool]:
        """Check EVM dependencies"""
        dependencies = super().check_dependencies()
        
        # Check Perl
        dependencies["perl"] = self.executable_path.exists()
        
        # Check EVM directory
        dependencies["evm_directory"] = self.executable_dir.exists()
        
        # Check EVM scripts
        for script_name, script_path in self.scripts.items():
            dependencies[f"script_{script_name}"] = script_path.exists()
        
        # Overall EVM availability
        dependencies["evm_available"] = self._check_evm_availability()
        
        return dependencies