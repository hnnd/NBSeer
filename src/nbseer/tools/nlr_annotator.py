"""
NLR-Annotator tool interface
NLR-Annotator工具接口

提供NLR基因定位功能的工具接口
"""

import shutil
import subprocess
from pathlib import Path
from typing import Any, Dict, List

from .base import ExternalTool, ToolResult
from ..utils.logging_setup import get_logger

logger = get_logger(__name__)


class NLRAnnotatorTool(ExternalTool):
    """
    Interface for NLR-Annotator tool
    NLR-Annotator工具接口
    
    NLR-Annotator用于识别植物基因组中的NLR基因候选区域
    """
    
    def __init__(self, config: Dict[str, Any] = None) -> None:
        """
        Initialize NLR-Annotator tool interface
        
        Args:
            config: Tool configuration dictionary
        """
        self.config = config or {}
        
        # NLR-Annotator is typically run via Java
        java_executable = self.config.get("executable", "java")
        super().__init__(
            tool_name="nlr_annotator",
            executable_path=self._find_java_executable(java_executable),
            config=config,
        )
        
        # Get JAR file path
        self.jar_path = Path(self.config.get("jar_path", "NLR-Annotator.jar"))
        if not self.jar_path.exists():
            logger.warning(f"NLR-Annotator JAR not found at: {self.jar_path}")
        
        # Get configuration file paths
        jar_dir = self.jar_path.parent
        self.mot_file = Path(self.config.get("mot_file", jar_dir / "mot.txt"))
        self.store_file = Path(self.config.get("store_file", jar_dir / "store.txt"))
        
        if not self.mot_file.exists():
            logger.warning(f"mot.txt configuration file not found at: {self.mot_file}")
        if not self.store_file.exists():
            logger.warning(f"store.txt configuration file not found at: {self.store_file}")
    
    def _find_java_executable(self, java_cmd: str) -> Path:
        """Find Java executable"""
        java_path = shutil.which(java_cmd)
        if not java_path:
            raise RuntimeError(f"Java executable not found: {java_cmd}")
        return Path(java_path)
    
    def prepare_input(self, **kwargs: Any) -> Dict[str, Path]:
        """
        Prepare input files for NLR-Annotator
        
        Args:
            genome_file: Path to genome FASTA file
            
        Returns:
            Dictionary of prepared input files
        """
        genome_file = kwargs.get("genome_file")
        if not genome_file:
            raise ValueError("genome_file is required for NLR-Annotator")
        
        genome_path = Path(genome_file)
        self.validate_input(genome_file=genome_path)
        
        return {"genome": genome_path}
    
    def build_command(
        self,
        input_files: Dict[str, Path],
        output_dir: Path,
        **kwargs: Any,
    ) -> List[str]:
        """
        Build NLR-Annotator command
        
        Args:
            input_files: Prepared input files
            output_dir: Output directory
            **kwargs: Additional parameters
            
        Returns:
            Command as list of strings
        """
        genome_file = input_files["genome"]
        
        # Basic command structure
        command = [str(self.executable_path)]
        
        # Add memory settings
        default_args = self.config.get("default_args", ["-Xmx8g"])
        if default_args:
            command.extend(default_args)
        
        # Add JAR file
        command.extend(["-jar", str(self.jar_path)])
        
        # Add input genome - use absolute path
        command.extend(["-i", str(genome_file.resolve())])
        
        # Add configuration files
        command.extend(["-x", str(self.mot_file)])
        command.extend(["-y", str(self.store_file)])
        
        # Add output file (GFF format) - output_dir should already be absolute
        output_gff = output_dir / "NLR_candidates.gff"
        command.extend(["-g", str(output_gff)])
        
        # Add parameters from config
        parameters = self.config.get("parameters", {})
        
        # Add thread count
        threads = parameters.get("threads", 1)
        command.extend(["-t", str(threads)])
        
        # Add number of sequences per thread
        seqs_per_thread = parameters.get("seqs_per_thread", 1000)
        command.extend(["-n", str(seqs_per_thread)])
        
        # Add distance parameters
        if "distanceWithinMotifCombination" in parameters:
            command.extend(["-distanceWithinMotifCombination", str(parameters["distanceWithinMotifCombination"])])
        
        if "distanceForElongating" in parameters:
            command.extend(["-distanceForElongating", str(parameters["distanceForElongating"])])
            
        if "distanceBetweenMotifCombinations" in parameters:
            command.extend(["-distanceBetweenMotifCombinations", str(parameters["distanceBetweenMotifCombinations"])])
        
        # Add any additional parameters (filter out parameters not for NLR-Annotator)
        nlr_specific_params = {
            "c", "o", "m", "a", "f", "b"  # Known NLR-Annotator parameters
        }
        for key, value in kwargs.items():
            if key not in ["genome_file", "augustus_model", "training_model"] and value is not None:
                if key in nlr_specific_params:
                    if isinstance(value, bool):
                        if value:
                            command.append(f"-{key}")
                    else:
                        command.extend([f"-{key}", str(value)])
        
        return command
    
    def parse_output(self, output_dir: Path) -> Dict[str, Any]:
        """
        Parse NLR-Annotator output files
        
        Args:
            output_dir: Directory containing output files
            
        Returns:
            Parsed results dictionary
        """
        results = {
            "candidate_count": 0,
            "output_files": [],
            "candidates": [],
        }
        
        # Look for the specific output GFF file
        gff_file = output_dir / "NLR_candidates.gff"
        if gff_file.exists():
            gff_files = [gff_file]
        else:
            # Fallback to search for any GFF files
            gff_files = list(output_dir.glob("*.gff*"))
        
        results["output_files"] = [str(f) for f in gff_files]
        
        # Count candidates from GFF files
        total_candidates = 0
        for gff_file in gff_files:
            try:
                with open(gff_file, 'r') as f:
                    for line in f:
                        if not line.startswith('#') and line.strip():
                            fields = line.strip().split('\t')
                            if len(fields) >= 3 and fields[2] in ['gene', 'NBSLRR', 'NBARC', 'CC-NBARC-LRR', 'NBARC-LRR', 'TIR-NBARC-LRR']:
                                total_candidates += 1
                                
                                # Extract basic info
                                if len(fields) >= 9:
                                    candidate_info = {
                                        "seqid": fields[0],
                                        "start": int(fields[3]),
                                        "end": int(fields[4]),
                                        "strand": fields[6],
                                        "attributes": fields[8],
                                    }
                                    results["candidates"].append(candidate_info)
            except Exception as e:
                logger.warning(f"Error parsing GFF file {gff_file}: {e}")
        
        results["candidate_count"] = total_candidates
        
        # Look for additional output files
        txt_files = list(output_dir.glob("*.txt"))
        if txt_files:
            results["output_files"].extend([str(f) for f in txt_files])
        
        return results
    
    def check_dependencies(self) -> Dict[str, bool]:
        """Check NLR-Annotator dependencies"""
        dependencies = super().check_dependencies()
        
        # Check Java
        dependencies["java"] = self.executable_path.exists()
        
        # Check JAR file
        dependencies["jar_file"] = self.jar_path.exists()
        
        # Check configuration files
        dependencies["mot_file"] = self.mot_file.exists()
        dependencies["store_file"] = self.store_file.exists()
        
        # Check Java version
        try:
            result = subprocess.run(
                [str(self.executable_path), "-version"],
                capture_output=True,
                text=True,
                timeout=10,
            )
            dependencies["java_version"] = result.returncode == 0
        except Exception:
            dependencies["java_version"] = False
        
        return dependencies