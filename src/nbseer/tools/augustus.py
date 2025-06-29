"""
Augustus tool interface
Augustus工具接口

提供基因结构预测功能的工具接口
"""

from pathlib import Path
from typing import Any, Dict, List

from .base import ExternalTool, ToolResult
from ..utils.logging_setup import get_logger

logger = get_logger(__name__)


class AugustusTool(ExternalTool):
    """
    Interface for Augustus gene prediction tool
    Augustus基因预测工具接口
    
    Augustus用于从头基因预测，特别是针对NLR基因候选区域
    """
    
    def __init__(self, config: Dict[str, Any] = None) -> None:
        """
        Initialize Augustus tool interface
        
        Args:
            config: Tool configuration dictionary
        """
        super().__init__(
            tool_name="augustus",
            config=config,
        )
        
        # Set Augustus config path if provided
        config_dir = self.config.get("config_dir")
        if config_dir:
            import os
            os.environ["AUGUSTUS_CONFIG_PATH"] = str(config_dir)
    
    def prepare_input(self, **kwargs: Any) -> Dict[str, Path]:
        """
        Prepare input files for Augustus
        
        Args:
            candidates_file: Path to NLR candidates file (FASTA or GFF)
            genome_file: Path to full genome file (required for GFF input)
            output_dir: Output directory for temporary files
            
        Returns:
            Dictionary of prepared input files
        """
        candidates_file = kwargs.get("candidates_file")
        genome_file = kwargs.get("genome_file")
        output_dir = kwargs.get("output_dir")
        
        if not candidates_file:
            raise ValueError("candidates_file is required for Augustus")
        
        candidates_path = Path(candidates_file)
        self.validate_input(candidates_file=candidates_path)
        
        # Check if candidates file is empty (no actual candidates)
        has_candidates = self._has_candidates(candidates_path)
        
        if not has_candidates:
            input_files = {
                "candidates": candidates_path,
                "has_candidates": False,
            }
            return input_files
        
        # Determine if input is GFF or FASTA
        is_gff = self._is_gff_file(candidates_path)
        
        if is_gff and genome_file:
            # Extract FASTA sequences from GFF coordinates
            genome_path = Path(genome_file)
            self.validate_input(genome_file=genome_path)
            
            if not output_dir:
                raise ValueError("output_dir is required when converting GFF to FASTA")
                
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
            
            fasta_file = self._extract_sequences_from_gff(
                candidates_path, genome_path, output_dir
            )
            
            input_files = {
                "candidates": fasta_file,
                "has_candidates": True,
                "genome": genome_path,
                "original_gff": candidates_path,
            }
        else:
            # Input is already FASTA or no genome provided
            input_files = {
                "candidates": candidates_path,
                "has_candidates": has_candidates,
            }
            
            if genome_file:
                genome_path = Path(genome_file)
                self.validate_input(genome_file=genome_path)
                input_files["genome"] = genome_path
        
        return input_files
    
    def _has_candidates(self, candidates_file: Path) -> bool:
        """
        Check if candidates file contains actual gene candidates
        检查候选文件是否包含实际的基因候选
        """
        try:
            with open(candidates_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith('#'):
                        fields = line.split('\t')
                        # Accept any GFF feature type (gene, NBSLRR, etc.)
                        if len(fields) >= 3:
                            return True
            return False
        except Exception:
            return False
    
    def _is_gff_file(self, file_path: Path) -> bool:
        """
        Check if file is in GFF format
        检查文件是否为GFF格式
        """
        try:
            with open(file_path, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith('#'):
                        fields = line.split('\t')
                        # GFF has 9 tab-separated fields
                        if len(fields) == 9:
                            return True
                        else:
                            return False
            return False
        except Exception:
            return False
    
    def _extract_sequences_from_gff(
        self, 
        gff_file: Path, 
        genome_file: Path, 
        output_dir: Path
    ) -> Path:
        """
        Extract FASTA sequences from GFF coordinates
        从GFF坐标提取FASTA序列
        """
        output_fasta = output_dir / "candidate_sequences.fasta"
        
        # Read genome sequences
        genome_sequences = {}
        with open(genome_file, 'r') as f:
            current_seq = None
            current_data = []
            
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_seq:
                        genome_sequences[current_seq] = ''.join(current_data)
                    current_seq = line[1:].split()[0]  # Get sequence name
                    current_data = []
                else:
                    current_data.append(line)
            
            if current_seq:
                genome_sequences[current_seq] = ''.join(current_data)
        
        # Extract candidate regions
        with open(output_fasta, 'w') as fasta_out:
            with open(gff_file, 'r') as gff_in:
                for line_num, line in enumerate(gff_in, 1):
                    line = line.strip()
                    if line and not line.startswith('#'):
                        fields = line.split('\t')
                        if len(fields) >= 5:
                            seqid = fields[0]
                            start = int(fields[3]) - 1  # Convert to 0-based
                            end = int(fields[4])
                            
                            if seqid in genome_sequences:
                                sequence = genome_sequences[seqid][start:end]
                                candidate_id = f"{seqid}_{start+1}_{end}"
                                
                                fasta_out.write(f">{candidate_id}\n")
                                # Write sequence in 80-character lines
                                for i in range(0, len(sequence), 80):
                                    fasta_out.write(sequence[i:i+80] + '\n')
                            else:
                                self.logger.warning(f"Sequence {seqid} not found in genome file")
        
        return output_fasta
    
    def build_command(
        self,
        input_files: Dict[str, Path],
        output_dir: Path,
        **kwargs: Any,
    ) -> List[str]:
        """
        Build Augustus command
        
        Args:
            input_files: Prepared input files
            output_dir: Output directory
            **kwargs: Additional parameters
            
        Returns:
            Command as list of strings
        """
        candidates_file = input_files["candidates"]
        output_file = output_dir / "predictions.gff"
        
        # Basic command structure
        command = [str(self.executable_path)]
        
        # Species model
        model = kwargs.get("model") or self.config.get("default_species", "rice")
        command.append(f"--species={model}")
        
        # Parameters from config
        parameters = self.config.get("parameters", {})
        
        # Gene model
        genemodel = parameters.get("genemodel", "complete")
        command.append(f"--genemodel={genemodel}")
        
        # Strand
        if parameters.get("singlestrand", False):
            command.append("--singlestrand=true")
        else:
            strand = parameters.get("strand", "both")
            command.append(f"--strand={strand}")
        
        # Output format
        if parameters.get("gff3", True):
            command.append("--gff3=on")
        
        # Unique gene IDs
        if parameters.get("uniqueGeneId", True):
            command.append("--uniqueGeneId=true")
        
        # Add any additional parameters from kwargs
        for key, value in kwargs.items():
            if key not in ["candidates_file", "genome_file", "model", "output_dir", "augustus_model", "training_model"] and value is not None:
                if isinstance(value, bool):
                    if value:
                        command.append(f"--{key}")
                else:
                    command.extend([f"--{key}", str(value)])
        
        # Input file
        command.append(str(candidates_file))
        
        # Output redirection
        command.extend([">", str(output_file)])
        
        return command
    
    def execute(self, output_dir: Path, timeout: int = None, **kwargs: Any) -> ToolResult:
        """
        Execute Augustus with shell redirection
        
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
            input_files = self.prepare_input(output_dir=output_dir, **kwargs)
            
            # Check if we have any candidates to process
            if not input_files.get("has_candidates", True):
                self.logger.info("No candidates found, creating empty prediction file")
                output_dir.mkdir(parents=True, exist_ok=True)
                output_file = output_dir / "predictions.gff"
                with open(output_file, 'w') as f:
                    f.write("##gff-version 3\n")
                    f.write("# No gene predictions - no NLR candidates found\n")
                
                execution_time = (datetime.now() - start_time).total_seconds()
                path_values = [v for v in input_files.values() if isinstance(v, Path)]
                return ToolResult(
                    tool_name=self.tool_name,
                    command="# Skipped - no candidates",
                    return_code=0,
                    stdout="",
                    stderr="",
                    execution_time=execution_time,
                    input_files=path_values,
                    output_files=[output_file],
                    success=True,
                    metadata={"config": self.config, "kwargs": kwargs, "skipped": True},
                )
            
            # Build command
            command_list = self.build_command(input_files, output_dir, **kwargs)
            
            # Convert to shell command with redirection
            output_file = output_dir / "predictions.gff"
            command_without_redirect = []
            for item in command_list:
                if item not in [">", str(output_file)]:
                    command_without_redirect.append(item)
            
            # Ensure output directory exists
            output_dir.mkdir(parents=True, exist_ok=True)
            
            # Log execution start
            self.logger.info(f"Starting {self.tool_name} execution")
            self.logger.info(f"Command: {' '.join(command_without_redirect)}")
            
            # Execute command with shell redirection
            with open(output_file, 'w') as f:
                result = subprocess.run(
                    command_without_redirect,
                    stdout=f,
                    stderr=subprocess.PIPE,
                    text=True,
                    timeout=timeout or self.config.get("timeout", 1800),
                    cwd=output_dir,
                )
            
            execution_time = (datetime.now() - start_time).total_seconds()
            
            # Get output files
            output_files = [output_file] if output_file.exists() else []
            
            # Convert coordinates if needed
            if result.returncode == 0 and output_file.exists():
                converted_gff = self.convert_coordinates_to_genome(output_file, input_files)
                if converted_gff != output_file:
                    output_files.append(converted_gff)
            
            # Create result object
            tool_result = ToolResult(
                tool_name=self.tool_name,
                command=" ".join(command_without_redirect) + f" > {output_file}",
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
    
    def convert_coordinates_to_genome(self, gff_file: Path, input_files: Dict[str, Path]) -> Path:
        """
        Convert Augustus predictions from candidate sequence coordinates to genome coordinates
        将Augustus预测结果从候选序列坐标转换为基因组坐标
        
        Args:
            gff_file: Augustus output GFF file with relative coordinates
            input_files: Input files containing coordinate mapping info
            
        Returns:
            Path to converted GFF file
        """
        if "original_gff" not in input_files:
            # No coordinate conversion needed - input was already genome coordinates
            return gff_file
            
        original_gff = input_files["original_gff"]
        converted_gff = gff_file.parent / "predictions_genome_coords.gff"
        
        # Parse original candidate coordinates
        candidate_coords = {}
        try:
            with open(original_gff, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith('#'):
                        fields = line.split('\t')
                        if len(fields) >= 5:
                            seqid = fields[0]
                            start = int(fields[3])  # 1-based
                            end = int(fields[4])
                            candidate_id = f"{seqid}_{start}_{end}"
                            candidate_coords[candidate_id] = {
                                "seqid": seqid,
                                "offset": start - 1  # Convert to 0-based offset
                            }
        except Exception as e:
            self.logger.warning(f"Error reading original GFF coordinates: {e}")
            return gff_file
        
        # Convert Augustus predictions
        try:
            with open(gff_file, 'r') as in_f, open(converted_gff, 'w') as out_f:
                for line in in_f:
                    if line.startswith('#'):
                        out_f.write(line)
                        continue
                        
                    line = line.strip()
                    if not line:
                        out_f.write('\n')
                        continue
                        
                    fields = line.split('\t')
                    if len(fields) >= 9:
                        candidate_seqid = fields[0]
                        
                        # Find matching candidate coordinates
                        coord_info = None
                        for cand_id, coord_data in candidate_coords.items():
                            if candidate_seqid in cand_id or cand_id.startswith(candidate_seqid):
                                coord_info = coord_data
                                break
                        
                        if coord_info:
                            # Convert coordinates
                            fields[0] = coord_info["seqid"]  # Use original chromosome/scaffold name
                            fields[3] = str(int(fields[3]) + coord_info["offset"])  # Adjust start
                            fields[4] = str(int(fields[4]) + coord_info["offset"])  # Adjust end
                            
                            # Update attributes to include original coordinates
                            attrs = fields[8]
                            if "orig_coords=" not in attrs:
                                attrs += f";orig_coords={candidate_seqid}"
                                fields[8] = attrs
                        
                        out_f.write('\t'.join(fields) + '\n')
                    else:
                        out_f.write(line + '\n')
                        
            self.logger.info(f"Converted coordinates for {len(candidate_coords)} candidates")
            return converted_gff
            
        except Exception as e:
            self.logger.error(f"Error converting coordinates: {e}")
            return gff_file

    def parse_output(self, output_dir: Path) -> Dict[str, Any]:
        """
        Parse Augustus output files
        
        Args:
            output_dir: Directory containing output files
            
        Returns:
            Parsed results dictionary
        """
        results = {
            "gene_count": 0,
            "output_files": [],
            "genes": [],
        }
        
        # Look for GFF output file
        gff_file = output_dir / "predictions.gff"
        if gff_file.exists():
            results["output_files"].append(str(gff_file))
            
            # Count genes
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
                                        results["genes"].append(gene_info)
                
                results["gene_count"] = gene_count
                
            except Exception as e:
                logger.warning(f"Error parsing Augustus output: {e}")
        
        return results
    
    def check_dependencies(self) -> Dict[str, bool]:
        """Check Augustus dependencies"""
        dependencies = super().check_dependencies()
        
        # Check Augustus config path
        import os
        config_path = os.environ.get("AUGUSTUS_CONFIG_PATH")
        if config_path:
            dependencies["config_path"] = Path(config_path).exists()
        else:
            dependencies["config_path"] = False
        
        # Check species models
        if config_path:
            species_dir = Path(config_path) / "species"
            dependencies["species_models"] = species_dir.exists()
        else:
            dependencies["species_models"] = False
        
        # Check version
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