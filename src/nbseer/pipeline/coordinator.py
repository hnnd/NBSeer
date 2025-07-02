"""
Pipeline coordinator for NBS gene annotation
NBS基因注释流水线协调器

统一管理和协调整个NBS基因注释流水线的执行过程
"""

import json
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Union
from dataclasses import dataclass, asdict

from ..utils.config import Config, load_config
from ..utils.logging_setup import get_logger
from ..utils.exceptions import PipelineError, ConfigurationError
from ..data.validation import DataValidator, ValidationResult
from ..tools.nlr_annotator import NLRAnnotatorTool
from ..tools.miniprot import MiniprotTool
from ..tools.augustus import AugustusTool
from ..tools.evm import EVMTool
from ..augustus_miniprot_trainer import AugustusMiniprotTrainer, MiniprotTrainingConfig

logger = get_logger(__name__)


@dataclass
class PipelineResults:
    """
    Container for pipeline execution results
    流水线执行结果容器
    """
    
    pipeline_id: str
    start_time: datetime
    end_time: Optional[datetime]
    status: str  # running, completed, failed
    stages_completed: List[str]
    input_files: Dict[str, Path]
    output_files: Dict[str, Path]
    stage_results: Dict[str, Dict[str, Any]]
    errors: List[str]
    metadata: Dict[str, Any]
    
    def _convert_paths_to_strings(self, obj: Any) -> Any:
        """Recursively convert Path objects to strings for JSON serialization"""
        if isinstance(obj, Path):
            return str(obj)
        elif isinstance(obj, dict):
            return {key: self._convert_paths_to_strings(value) for key, value in obj.items()}
        elif isinstance(obj, list):
            return [self._convert_paths_to_strings(item) for item in obj]
        elif isinstance(obj, tuple):
            return tuple(self._convert_paths_to_strings(item) for item in obj)
        elif hasattr(obj, '__dict__'):
            # Handle objects like ValidationResult that have Path attributes
            return self._convert_paths_to_strings(obj.__dict__)
        else:
            return obj
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert results to dictionary for serialization"""
        result = asdict(self)
        # Convert Path objects to strings
        result["input_files"] = {k: str(v) for k, v in self.input_files.items()}
        result["output_files"] = {k: str(v) for k, v in self.output_files.items()}
        # Convert datetime objects to ISO format
        result["start_time"] = self.start_time.isoformat()
        if self.end_time:
            result["end_time"] = self.end_time.isoformat()
        # Recursively convert any Path objects in stage_results and metadata
        result["stage_results"] = self._convert_paths_to_strings(result["stage_results"])
        result["metadata"] = self._convert_paths_to_strings(result["metadata"])
        return result
    
    def save_to_file(self, output_path: Path) -> None:
        """Save results to JSON file"""
        with open(output_path, 'w', encoding='utf-8') as f:
            json.dump(self.to_dict(), f, indent=2, ensure_ascii=False)


class NBSAnnotationPipeline:
    """
    Main coordinator class for NBS gene annotation pipeline
    NBS基因注释流水线的主要协调器类
    
    管理整个流水线的执行，包括：
    1. NLR基因定位 / NLR gene localization
    2. 蛋白质比对 / Protein alignment  
    3. 基因结构预测 / Gene structure prediction
    4. 模型训练 / Model training
    5. 证据整合 / Evidence integration
    """
    
    def __init__(
        self,
        config: Optional[Union[str, Path, Config]] = None,
        output_base_dir: Optional[Path] = None,
        pipeline_id: Optional[str] = None,
    ) -> None:
        """
        Initialize pipeline coordinator
        
        Args:
            config: Configuration file path or Config object / 配置文件路径或配置对象
            output_base_dir: Base directory for all outputs / 所有输出的基础目录
            pipeline_id: Unique identifier for this pipeline run / 此次流水线运行的唯一标识符
        """
        self.pipeline_id = pipeline_id or self._generate_pipeline_id()
        self.logger = get_logger(f"{__name__}.{self.pipeline_id}")
        
        # Load configuration
        if isinstance(config, Config):
            self.config = config
        else:
            self.config = load_config(config)
            
        # Set up directories 
        self.output_base_dir = output_base_dir or Path("output")
        # Use output directory directly for checkpoint support
        if pipeline_id:
            # If pipeline_id is explicitly provided, use subdirectory
            self.pipeline_output_dir = self.output_base_dir / self.pipeline_id
        else:
            # Use output directory directly for stable checkpoint paths
            self.pipeline_output_dir = self.output_base_dir
        self.pipeline_output_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize data validator
        self.validator = DataValidator()
        
        # Initialize tools (lazy loading)
        self._tools: Dict[str, Any] = {}
        
        # Initialize results
        self.results = PipelineResults(
            pipeline_id=self.pipeline_id,
            start_time=datetime.now(),
            end_time=None,
            status="initialized",
            stages_completed=[],
            input_files={},
            output_files={},
            stage_results={},
            errors=[],
            metadata={"config": self.config.to_dict()},
        )
        
        self.logger.info(f"Initialized pipeline {self.pipeline_id}")
    
    def _generate_pipeline_id(self) -> str:
        """Generate unique pipeline identifier"""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        return f"nbs_pipeline_{timestamp}"
    
    def _get_tool(self, tool_name: str) -> Any:
        """
        Get tool instance (lazy loading)
        获取工具实例（延迟加载）
        """
        if tool_name not in self._tools:
            tool_config = self.config.get_tool_config(tool_name)
            # Convert ToolConfig dataclass to dictionary for tool constructors
            tool_config_dict = tool_config.to_dict()
            
            if tool_name == "nlr_annotator":
                self._tools[tool_name] = NLRAnnotatorTool(config=tool_config_dict)
            elif tool_name == "miniprot":
                self._tools[tool_name] = MiniprotTool(config=tool_config_dict)
            elif tool_name == "augustus":
                self._tools[tool_name] = AugustusTool(config=tool_config_dict)
            elif tool_name == "evm":
                self._tools[tool_name] = EVMTool(config=tool_config_dict)
            else:
                raise ConfigurationError(
                    f"Unknown tool: {tool_name}",
                    config_key="tools",
                )
                
        return self._tools[tool_name]
    
    def validate_inputs(
        self,
        genome_path: Path,
        protein_path: Path,
        **kwargs: Any,
    ) -> Dict[str, ValidationResult]:
        """
        Validate all input files
        验证所有输入文件
        
        Args:
            genome_path: Path to genome FASTA file / 基因组FASTA文件路径
            protein_path: Path to protein FASTA file / 蛋白质FASTA文件路径
            
        Returns:
            Dictionary of validation results
        """
        self.logger.info("Validating input files...")
        
        validation_results = {}
        
        # Validate genome file
        try:
            validation_results["genome"] = self.validator.validate_fasta(genome_path)
        except Exception as e:
            self.logger.error(f"Genome validation failed: {e}")
            validation_results["genome"] = ValidationResult(
                is_valid=False,
                error_message=str(e),
                file_path=genome_path,
            )
        
        # Validate protein file
        try:
            validation_results["protein"] = self.validator.validate_protein_fasta(protein_path)
        except Exception as e:
            self.logger.error(f"Protein validation failed: {e}")
            validation_results["protein"] = ValidationResult(
                is_valid=False,
                error_message=str(e),
                file_path=protein_path,
            )
        
        # Check overall validation status
        all_valid = all(result.is_valid for result in validation_results.values())
        
        if not all_valid:
            failed_files = [
                str(result.file_path) 
                for result in validation_results.values() 
                if not result.is_valid
            ]
            raise PipelineError(
                f"Input validation failed for files: {', '.join(failed_files)}",
                pipeline_stage="validation",
                input_files=failed_files,
            )
        
        self.logger.info("Input validation completed successfully")
        return validation_results
    
    def run_nlr_localization(
        self,
        genome_path: Path,
        **kwargs: Any,
    ) -> Dict[str, Any]:
        """
        Run NLR gene localization stage
        运行NLR基因定位阶段
        
        Args:
            genome_path: Path to genome FASTA file
            
        Returns:
            Results from NLR localization
        """
        self.logger.info("Starting NLR gene localization...")
        stage_name = "nlr_localization"
        
        try:
            # Get NLR annotator tool
            nlr_tool = self._get_tool("nlr_annotator")
            
            # Create stage output directory
            stage_output_dir = self.pipeline_output_dir / stage_name
            
            # Execute NLR annotation
            tool_result = nlr_tool.execute(
                output_dir=stage_output_dir,
                genome_file=genome_path,
                **kwargs,
            )
            
            if not tool_result.success:
                raise PipelineError(
                    f"NLR localization failed: {tool_result.stderr}",
                    pipeline_stage=stage_name,
                )
            
            # Parse results
            parsed_results = nlr_tool.parse_output(stage_output_dir)
            
            # Create empty GFF file if no candidates were found
            candidates_file = stage_output_dir / "NLR_candidates.gff"
            if not candidates_file.exists() or parsed_results.get('candidate_count', 0) == 0:
                self.logger.info("No NLR candidates found, creating empty GFF file")
                candidates_file.parent.mkdir(parents=True, exist_ok=True)
                with open(candidates_file, 'w') as f:
                    f.write("##gff-version 3\n")
                    f.write("# No NLR candidates found in genome\n")
            
            # Update pipeline results
            self.results.stages_completed.append(stage_name)
            self.results.stage_results[stage_name] = {
                "tool_result": tool_result.to_dict(),
                "parsed_results": parsed_results,
            }
            self.results.output_files[f"{stage_name}_candidates"] = candidates_file
            
            self.logger.info(f"NLR localization completed. Found {parsed_results.get('candidate_count', 0)} candidates")
            return parsed_results
            
        except Exception as e:
            error_msg = f"NLR localization failed: {e}"
            self.logger.error(error_msg)
            self.results.errors.append(error_msg)
            raise PipelineError(error_msg, pipeline_stage=stage_name) from e
    
    def run_protein_alignment(
        self,
        genome_path: Path,
        protein_path: Path,
        **kwargs: Any,
    ) -> Dict[str, Any]:
        """
        Run protein alignment stage using miniprot
        使用miniprot运行蛋白质比对阶段
        
        Args:
            genome_path: Path to genome FASTA file
            protein_path: Path to protein FASTA file
            
        Returns:
            Results from protein alignment
        """
        self.logger.info("Starting protein alignment...")
        stage_name = "protein_alignment"
        
        try:
            # Get miniprot tool
            miniprot_tool = self._get_tool("miniprot")
            
            # Create stage output directory
            stage_output_dir = self.pipeline_output_dir / stage_name
            
            # Execute protein alignment
            # Filter out Augustus/training-specific parameters that miniprot doesn't understand
            miniprot_kwargs = {k: v for k, v in kwargs.items() 
                             if k not in ["augustus_model", "training_model", "training_species_name", 
                                        "training_quality", "training_min_genes", "training_optimization_rounds",
                                        "training_timeout", "enable_training"]}
            
            tool_result = miniprot_tool.execute(
                output_dir=stage_output_dir,
                genome_file=genome_path,
                protein_file=protein_path,
                **miniprot_kwargs,
            )
            
            if not tool_result.success:
                raise PipelineError(
                    f"Protein alignment failed: {tool_result.stderr}",
                    pipeline_stage=stage_name,
                )
            
            # Parse results
            parsed_results = miniprot_tool.parse_output(stage_output_dir)
            
            # Update pipeline results
            self.results.stages_completed.append(stage_name)
            self.results.stage_results[stage_name] = {
                "tool_result": tool_result.to_dict(),
                "parsed_results": parsed_results,
            }
            
            # Use the actual output file from tool_result (filtered if available)
            actual_output_file = tool_result.output_files[0] if tool_result.output_files else stage_output_dir / "alignments.gff"
            self.results.output_files[f"{stage_name}_alignments"] = actual_output_file
            
            self.logger.info(f"Protein alignment completed. Generated {parsed_results.get('alignment_count', 0)} alignments")
            return parsed_results
            
        except Exception as e:
            error_msg = f"Protein alignment failed: {e}"
            self.logger.error(error_msg)
            self.results.errors.append(error_msg)
            raise PipelineError(error_msg, pipeline_stage=stage_name) from e
    
    def run_gene_prediction(
        self,
        nlr_candidates_path: Path,
        genome_path: Path,
        training_model: Optional[str] = None,
        **kwargs: Any,
    ) -> Dict[str, Any]:
        """
        Run gene structure prediction using Augustus
        使用Augustus运行基因结构预测
        
        Args:
            nlr_candidates_path: Path to NLR candidate regions
            genome_path: Path to genome FASTA file
            training_model: Name of Augustus training model to use
            
        Returns:
            Results from gene prediction
        """
        self.logger.info("Starting gene structure prediction...")
        stage_name = "gene_prediction"
        
        try:
            # Get Augustus tool
            augustus_tool = self._get_tool("augustus")
            
            # Create stage output directory
            stage_output_dir = self.pipeline_output_dir / stage_name
            
            # Execute gene prediction
            # Get Augustus model from kwargs or config
            augustus_config = self.config.get("tools", {}).get("augustus", {})
            default_model = augustus_config.get("default_species", "rice")
            model = kwargs.get("augustus_model") or training_model or default_model
            
            # Filter out training-specific parameters that don't belong in gene prediction
            prediction_kwargs = {k: v for k, v in kwargs.items() 
                                if k not in ["training_species_name", "training_quality", 
                                           "training_min_genes", "training_optimization_rounds",
                                           "training_timeout", "enable_training"]}
            
            tool_result = augustus_tool.execute(
                output_dir=stage_output_dir,
                candidates_file=nlr_candidates_path,
                genome_file=genome_path,
                model=model,
                **prediction_kwargs,
            )
            
            if not tool_result.success:
                raise PipelineError(
                    f"Gene prediction failed: {tool_result.stderr}",
                    pipeline_stage=stage_name,
                )
            
            # Parse results
            parsed_results = augustus_tool.parse_output(stage_output_dir)
            
            # Update pipeline results
            self.results.stages_completed.append(stage_name)
            self.results.stage_results[stage_name] = {
                "tool_result": tool_result.to_dict(),
                "parsed_results": parsed_results,
            }
            self.results.output_files[f"{stage_name}_predictions"] = (
                stage_output_dir / "predictions.gff"
            )
            
            self.logger.info(f"Gene prediction completed. Predicted {parsed_results.get('gene_count', 0)} genes")
            return parsed_results
            
        except Exception as e:
            error_msg = f"Gene prediction failed: {e}"
            self.logger.error(error_msg)
            self.results.errors.append(error_msg)
            raise PipelineError(error_msg, pipeline_stage=stage_name) from e
    
    def run_evidence_integration(
        self,
        prediction_files: List[Path],
        **kwargs: Any,
    ) -> Dict[str, Any]:
        """
        Run evidence integration using EVM
        使用EVM运行证据整合
        
        Args:
            prediction_files: List of prediction files to integrate
            
        Returns:
            Results from evidence integration
        """
        self.logger.info("Starting evidence integration...")
        stage_name = "evidence_integration"
        
        try:
            # Get EVM tool
            evm_tool = self._get_tool("evm")
            
            # Create stage output directory
            stage_output_dir = self.pipeline_output_dir / stage_name
            
            # Execute evidence integration
            tool_result = evm_tool.execute(
                output_dir=stage_output_dir,
                prediction_files=prediction_files,
                **kwargs,
            )
            
            if not tool_result.success:
                raise PipelineError(
                    f"Evidence integration failed: {tool_result.stderr}",
                    pipeline_stage=stage_name,
                )
            
            # Parse results
            parsed_results = evm_tool.parse_output(stage_output_dir)
            
            # Update pipeline results
            self.results.stages_completed.append(stage_name)
            self.results.stage_results[stage_name] = {
                "tool_result": tool_result.to_dict(),
                "parsed_results": parsed_results,
            }
            self.results.output_files[f"{stage_name}_final"] = (
                stage_output_dir / "final_annotations.gff"
            )
            
            self.logger.info(f"Evidence integration completed. Final annotation count: {parsed_results.get('final_gene_count', 0)}")
            return parsed_results
            
        except Exception as e:
            error_msg = f"Evidence integration failed: {e}"
            self.logger.error(error_msg)
            self.results.errors.append(error_msg)
            raise PipelineError(error_msg, pipeline_stage=stage_name) from e
    
    def run_augustus_training(
        self,
        miniprot_results_path: Optional[Path] = None,
        training_species_name: Optional[str] = None,
        **kwargs: Any,
    ) -> Dict[str, Any]:
        """
        Run Augustus model training using high-quality Miniprot results
        使用高质量Miniprot结果运行Augustus模型训练
        
        Args:
            miniprot_results_path: Path to Miniprot GFF3 results (auto-detected if None)
            training_species_name: Name for the new Augustus species model
            
        Returns:
            Results from Augustus training
        """
        self.logger.info("Starting Augustus model training...")
        stage_name = "augustus_training"
        
        try:
            # Get training configuration from pipeline config
            augustus_config = self.config.get("tools", {}).get("augustus", {})
            training_config = augustus_config.get("training", {}).get("miniprot_training", {})
            
            # Skip training if disabled
            if not training_config.get("enabled", False):
                self.logger.info("Augustus training is disabled in configuration")
                return {"training_skipped": True, "reason": "disabled_in_config"}
            
            # Get training parameters
            species_name = training_species_name or training_config.get("species_name", "nbs_trained_model")
            quality_filter = training_config.get("quality_filter", "high")
            min_training_genes = training_config.get("min_training_genes", 20)
            
            # Create stage output directory
            stage_output_dir = self.pipeline_output_dir / stage_name
            stage_output_dir.mkdir(parents=True, exist_ok=True)
            
            # Create training configuration
            miniprot_training_config = MiniprotTrainingConfig(
                species_name=species_name,
                genome_file=str(self.results.input_files.get("genome", "")),
                working_dir=str(stage_output_dir),
                miniprot_gff_file=str(miniprot_results_path) if miniprot_results_path else "",
                miniprot_quality_filter=quality_filter,
                min_training_genes=min_training_genes,
                augustus_scripts_path=augustus_config.get("scripts_path", "/home/wangys/opt/Augustus/scripts"),
                augustus_config_path=augustus_config.get("config_dir", "/home/wangys/opt/Augustus/config"),
                flanking_dna_length=training_config.get("flanking_dna_length", 4000),
                optimization_rounds=training_config.get("optimization_rounds", 1),
                cpus=training_config.get("cpus", 4),
                timeout_minutes=training_config.get("timeout_minutes", 240),
                min_identity_threshold=training_config.get("min_identity_threshold", 0.95),
                max_frameshifts=training_config.get("max_frameshifts", 0),
                max_stop_codons=training_config.get("max_stop_codons", 0),
                backup_existing_model=training_config.get("backup_existing_model", True),
                create_training_report=training_config.get("create_training_report", True)
            )
            
            # Initialize Augustus trainer
            trainer = AugustusMiniprotTrainer(miniprot_training_config, self.logger)
            
            # Execute training
            self.logger.info(f"Training Augustus model: {species_name}")
            self.logger.info(f"Quality filter: {quality_filter}")
            self.logger.info(f"Working directory: {stage_output_dir}")
            
            training_result = trainer.train_augustus_model()
            
            # Check training success
            if not training_result.training_success:
                raise PipelineError(
                    f"Augustus training failed: {training_result.error_message}",
                    pipeline_stage=stage_name,
                )
            
            # Prepare results
            parsed_results = {
                "training_success": training_result.training_success,
                "species_name": training_result.species_name,
                "training_time_minutes": training_result.training_time_minutes,
                "filtered_training_genes": training_result.filtered_training_genes,
                "quality_filter_applied": training_result.quality_filter_applied,
                "model_directory": training_result.model_directory,
                "validation_passed": training_result.validation_passed,
                "training_gff_file": training_result.training_gff_file,
                "log_file": training_result.log_file,
                "report_file": training_result.report_file,
                "model_files": training_result.model_files
            }
            
            # Update pipeline results
            self.results.stages_completed.append(stage_name)
            self.results.stage_results[stage_name] = {
                "training_result": parsed_results,
                "config_used": miniprot_training_config.__dict__
            }
            
            # Store trained model info for later use
            self.results.metadata["trained_augustus_model"] = training_result.species_name
            
            # Output files
            if training_result.training_gff_file:
                self.results.output_files[f"{stage_name}_training_data"] = Path(training_result.training_gff_file)
            if training_result.report_file:
                self.results.output_files[f"{stage_name}_report"] = Path(training_result.report_file)
            if training_result.log_file:
                self.results.output_files[f"{stage_name}_log"] = Path(training_result.log_file)
            
            self.logger.info(f"Augustus training completed successfully!")
            self.logger.info(f"Trained model: {training_result.species_name}")
            self.logger.info(f"Training genes: {training_result.filtered_training_genes}")
            self.logger.info(f"Training time: {training_result.training_time_minutes:.1f} minutes")
            self.logger.info(f"Model validation: {'PASSED' if training_result.validation_passed else 'FAILED'}")
            
            return parsed_results
            
        except Exception as e:
            error_msg = f"Augustus training failed: {e}"
            self.logger.error(error_msg)
            self.results.errors.append(error_msg)
            raise PipelineError(error_msg, pipeline_stage=stage_name) from e
    
    def _check_stage_completion(self, stage_name: str) -> bool:
        """
        Check if a pipeline stage has already been completed
        检查流水线阶段是否已完成
        
        Args:
            stage_name: Name of the stage to check
            
        Returns:
            True if stage is completed, False otherwise
        """
        expected_outputs = {
            "nlr_localization": ["NLR_candidates.gff"],
            "protein_alignment": ["alignments.gff"],
            "augustus_training": [],  # Training doesn't produce standard output files
            "gene_prediction": ["predictions.gff", "predictions_genome_coords.gff"],
            "evidence_integration": ["final_annotations.gff"],
        }
        
        if stage_name not in expected_outputs:
            return False
            
        stage_dir = self.pipeline_output_dir / stage_name
        if not stage_dir.exists():
            return False
            
        # Check if all expected output files exist and are non-empty
        for output_file in expected_outputs[stage_name]:
            file_path = stage_dir / output_file
            if not file_path.exists() or file_path.stat().st_size == 0:
                return False
                
        return True

    def run_full_pipeline(
        self,
        genome_path: Path,
        protein_path: Path,
        **kwargs: Any,
    ) -> PipelineResults:
        """
        Run the complete NBS annotation pipeline with checkpoint support
        运行完整的NBS注释流水线并支持检查点
        
        Args:
            genome_path: Path to genome FASTA file
            protein_path: Path to protein FASTA file
            
        Returns:
            Complete pipeline results
        """
        self.logger.info(f"Starting full NBS annotation pipeline: {self.pipeline_id}")
        self.results.status = "running"
        self.results.input_files = {
            "genome": genome_path,
            "protein": protein_path,
        }
        
        try:
            # Stage 1: Validate inputs
            self.logger.info("Stage 1: Input validation")
            validation_results = self.validate_inputs(genome_path, protein_path, **kwargs)
            self.results.stage_results["validation"] = validation_results
            
            # Stage 2: NLR gene localization
            stage_name = "nlr_localization"
            if self._check_stage_completion(stage_name):
                self.logger.info(f"Stage 2: {stage_name} - SKIPPED (already completed)")
                nlr_candidates_file = self.pipeline_output_dir / stage_name / "NLR_candidates.gff"
                self.results.output_files[f"{stage_name}_candidates"] = nlr_candidates_file
                self.results.stages_completed.append(stage_name)
            else:
                self.logger.info("Stage 2: NLR gene localization")
                nlr_results = self.run_nlr_localization(genome_path, **kwargs)
                nlr_candidates_file = self.results.output_files["nlr_localization_candidates"]
            
            # Stage 3: Protein alignment
            stage_name = "protein_alignment"
            if self._check_stage_completion(stage_name):
                self.logger.info(f"Stage 3: {stage_name} - SKIPPED (already completed)")
                # Check for filtered output first, then fall back to original
                filtered_file = self.pipeline_output_dir / stage_name / "alignments_filtered.gff"
                original_file = self.pipeline_output_dir / stage_name / "alignments.gff"
                
                if filtered_file.exists():
                    alignment_file = filtered_file
                    self.logger.info(f"Using filtered miniprot results: {filtered_file}")
                else:
                    alignment_file = original_file
                    self.logger.info(f"Using original miniprot results: {original_file}")
                
                self.results.output_files[f"{stage_name}_alignments"] = alignment_file
                self.results.stages_completed.append(stage_name)
            else:
                self.logger.info("Stage 3: Protein alignment")
                alignment_results = self.run_protein_alignment(genome_path, protein_path, **kwargs)
                alignment_file = self.results.output_files["protein_alignment_alignments"]
            
            # Stage 4: Augustus model training (optional)
            augustus_model_to_use = None
            stage_name = "augustus_training"
            
            # Check if training is enabled and should be performed
            augustus_config = self.config.get("tools", {}).get("augustus", {})
            training_config = augustus_config.get("training", {}).get("miniprot_training", {})
            
            if training_config.get("enabled", False):
                # Check if we should run training or skip
                training_completed = False
                try:
                    # Check if training was already completed
                    training_results_file = self.pipeline_output_dir / stage_name / "training_completed.json"
                    if training_results_file.exists():
                        self.logger.info(f"Stage 4a: {stage_name} - SKIPPED (already completed)")
                        # Load previous training results to get model name
                        import json
                        with open(training_results_file, 'r') as f:
                            training_data = json.load(f)
                            augustus_model_to_use = training_data.get("species_name")
                        training_completed = True
                        self.results.stages_completed.append(stage_name)
                    else:
                        self.logger.info("Stage 4a: Augustus model training")
                        
                        # Get high-quality Miniprot results for training
                        miniprot_training_data = None
                        
                        # Look for high-quality filtered results
                        if "protein_alignment_alignments" in self.results.output_files:
                            alignment_file = self.results.output_files["protein_alignment_alignments"]
                            # Check if this is a filtered file or if high-quality version exists
                            high_quality_file = alignment_file.parent / "filtered" / "miniprot_high_quality.gff3"
                            if high_quality_file.exists():
                                miniprot_training_data = high_quality_file
                                self.logger.info(f"Using high-quality Miniprot data for training: {high_quality_file}")
                            else:
                                self.logger.warning("High-quality Miniprot data not found, attempting to use regular alignments")
                                miniprot_training_data = alignment_file
                        
                        if miniprot_training_data and miniprot_training_data.exists():
                            # Run Augustus training
                            # Filter out parameters that are already explicitly specified
                            training_kwargs = {k: v for k, v in kwargs.items() 
                                             if k not in ["training_species_name"]}
                            training_results = self.run_augustus_training(
                                miniprot_results_path=miniprot_training_data,
                                training_species_name=kwargs.get("training_species_name"),
                                **training_kwargs
                            )
                            
                            if training_results.get("training_success", False):
                                augustus_model_to_use = training_results["species_name"]
                                training_completed = True
                                
                                # Save training completion marker
                                training_results_file.parent.mkdir(parents=True, exist_ok=True)
                                with open(training_results_file, 'w') as f:
                                    json.dump(training_results, f, indent=2)
                                    
                                self.logger.info(f"Augustus training completed! New model: {augustus_model_to_use}")
                            else:
                                self.logger.warning("Augustus training failed, will use default model for prediction")
                        else:
                            self.logger.warning("No suitable Miniprot data found for training, skipping Augustus training")
                            
                except Exception as e:
                    self.logger.error(f"Augustus training failed: {e}")
                    self.logger.warning("Continuing with default Augustus model")
            else:
                self.logger.info("Stage 4a: Augustus training - DISABLED in configuration")
            
            # Stage 4b: Gene structure prediction  
            stage_name = "gene_prediction"
            if self._check_stage_completion(stage_name):
                self.logger.info(f"Stage 4b: {stage_name} - SKIPPED (already completed)")
                prediction_file = self.pipeline_output_dir / stage_name / "predictions.gff"
                # Check for genome coordinate version (preferred for EVM)
                genome_coords_file = self.pipeline_output_dir / stage_name / "predictions_genome_coords.gff"
                if genome_coords_file.exists():
                    prediction_file = genome_coords_file
                    self.logger.info(f"Using genome coordinate file for EVM: {genome_coords_file}")
                else:
                    self.logger.info(f"Using original predictions file: {prediction_file}")
                self.results.output_files[f"{stage_name}_predictions"] = prediction_file
                self.results.stages_completed.append(stage_name)
            else:
                self.logger.info("Stage 4b: Gene structure prediction")
                # Use trained model if available, otherwise use specified or default model
                training_model = augustus_model_to_use or kwargs.get("augustus_model") or kwargs.get("training_model")
                if augustus_model_to_use:
                    self.logger.info(f"Using newly trained Augustus model: {augustus_model_to_use}")
                elif training_model:
                    self.logger.info(f"Using specified Augustus model: {training_model}")
                else:
                    self.logger.info("Using default Augustus model from configuration")
                    
                prediction_results = self.run_gene_prediction(
                    nlr_candidates_file, genome_path, training_model=training_model, **kwargs
                )
                # Use genome coordinate version if available
                prediction_file = self.results.output_files["gene_prediction_predictions"]
                genome_coords_file = prediction_file.parent / "predictions_genome_coords.gff"
                if genome_coords_file.exists():
                    prediction_file = genome_coords_file
                    self.results.output_files["gene_prediction_predictions"] = prediction_file
                    self.logger.info(f"Updated to use genome coordinate file: {genome_coords_file}")
            
            # Stage 5: Evidence integration
            stage_name = "evidence_integration"
            if self._check_stage_completion(stage_name):
                self.logger.info(f"Stage 5: {stage_name} - SKIPPED (already completed)")
                final_file = self.pipeline_output_dir / stage_name / "final_annotations.gff"
                self.results.output_files[f"{stage_name}_final"] = final_file
                self.results.stages_completed.append(stage_name)
            else:
                self.logger.info("Stage 5: Evidence integration")
                self.logger.info(f"EVM input files:")
                self.logger.info(f"  - Augustus predictions: {prediction_file}")
                self.logger.info(f"  - Miniprot alignments: {alignment_file}")
                self.logger.info(f"  - NLR candidates: {nlr_candidates_file}")
                integration_results = self.run_evidence_integration(
                    [alignment_file, prediction_file], 
                    genome_file=genome_path,
                    nlr_candidates_file=nlr_candidates_file,
                    **kwargs
                )
            
            # Pipeline completed successfully
            self.results.status = "completed"
            self.results.end_time = datetime.now()
            
            # Save results
            results_file = self.pipeline_output_dir / "pipeline_results.json"
            self.results.save_to_file(results_file)
            
            self.logger.info(f"Pipeline {self.pipeline_id} completed successfully")
            return self.results
            
        except Exception as e:
            self.results.status = "failed"
            self.results.end_time = datetime.now()
            error_msg = f"Pipeline {self.pipeline_id} failed: {e}"
            self.logger.error(error_msg)
            self.results.errors.append(error_msg)
            
            # Save error results
            results_file = self.pipeline_output_dir / "pipeline_results.json"
            self.results.save_to_file(results_file)
            
            raise PipelineError(error_msg, pipeline_stage="full_pipeline") from e
    
    def run_pipeline_from_stage(
        self,
        start_stage: str,
        genome_path: Path,
        protein_path: Path,
        **kwargs: Any,
    ) -> PipelineResults:
        """
        Resume pipeline execution from a specific stage
        从特定阶段恢复流水线执行
        
        Args:
            start_stage: Stage to start from
            genome_path: Path to genome FASTA file
            protein_path: Path to protein FASTA file
            
        Returns:
            Pipeline results from resume point
        """
        # Define stage execution order
        stage_order = [
            "nlr_localization", 
            "protein_alignment", 
            "augustus_training", 
            "gene_prediction", 
            "evidence_integration"
        ]
        
        if start_stage not in stage_order:
            raise PipelineError(f"Invalid start stage: {start_stage}")
        
        self.logger.info(f"Resuming pipeline from stage: {start_stage}")
        
        # Find the starting index
        start_idx = stage_order.index(start_stage)
        
        # Initialize results
        self.results.status = "running"
        self.results.input_files = {
            "genome": genome_path,
            "protein": protein_path,
        }
        
        try:
            # Stage 1: Always validate inputs
            self.logger.info("Stage 0: Input validation")
            validation_results = self.validate_inputs(genome_path, protein_path, **kwargs)
            self.results.stage_results["validation"] = validation_results
            
            # Check for existing outputs from previous stages
            nlr_candidates_file = None
            alignment_file = None
            augustus_model_to_use = None
            
            # Initialize variables for required files
            if start_idx > 0:  # If not starting from beginning
                # Look for NLR candidates from previous run
                nlr_candidates_file = self.pipeline_output_dir / "nlr_localization" / "NLR_candidates.gff"
                if not nlr_candidates_file.exists():
                    raise PipelineError(f"Missing prerequisite: NLR candidates file not found at {nlr_candidates_file}")
                self.results.output_files["nlr_localization_candidates"] = nlr_candidates_file
                if "nlr_localization" not in self.results.stages_completed:
                    self.results.stages_completed.append("nlr_localization")
            
            if start_idx > 1:  # If starting after protein alignment
                # Look for protein alignment results
                filtered_file = self.pipeline_output_dir / "protein_alignment" / "alignments_filtered.gff"
                original_file = self.pipeline_output_dir / "protein_alignment" / "alignments.gff"
                
                if filtered_file.exists():
                    alignment_file = filtered_file
                elif original_file.exists():
                    alignment_file = original_file
                else:
                    raise PipelineError(f"Missing prerequisite: Protein alignment file not found")
                
                self.results.output_files["protein_alignment_alignments"] = alignment_file
                if "protein_alignment" not in self.results.stages_completed:
                    self.results.stages_completed.append("protein_alignment")
            
            # Execute stages from start_stage onwards
            for stage in stage_order[start_idx:]:
                if stage == "nlr_localization":
                    self.logger.info("Stage 1: NLR gene localization")
                    nlr_results = self.run_nlr_localization(genome_path, **kwargs)
                    nlr_candidates_file = self.results.output_files["nlr_localization_candidates"]
                
                elif stage == "protein_alignment":
                    self.logger.info("Stage 2: Protein alignment")
                    alignment_results = self.run_protein_alignment(genome_path, protein_path, **kwargs)
                    alignment_file = self.results.output_files["protein_alignment_alignments"]
                
                elif stage == "augustus_training":
                    # Check if training is enabled
                    augustus_config = self.config.get("tools", {}).get("augustus", {})
                    training_config = augustus_config.get("training", {}).get("miniprot_training", {})
                    
                    if training_config.get("enabled", False):
                        self.logger.info("Stage 3: Augustus model training")
                        training_results = self.run_augustus_training(**kwargs)
                        augustus_model_to_use = training_results.get("species_name")
                    else:
                        self.logger.info("Stage 3: Augustus training - SKIPPED (disabled)")
                
                elif stage == "gene_prediction":
                    self.logger.info("Stage 4: Gene prediction")
                    if not nlr_candidates_file:
                        raise PipelineError("NLR candidates file required for gene prediction")
                    
                    prediction_results = self.run_gene_prediction(
                        nlr_candidates_path=nlr_candidates_file,
                        training_model=augustus_model_to_use,
                        **kwargs
                    )
                    prediction_file = self.results.output_files["gene_prediction_predictions"]
                
                elif stage == "evidence_integration":
                    self.logger.info("Stage 5: Evidence integration")
                    if not alignment_file or not nlr_candidates_file:
                        raise PipelineError("Missing prerequisite files for evidence integration")
                    
                    # Get prediction file
                    prediction_file = self.results.output_files.get("gene_prediction_predictions")
                    if not prediction_file or not prediction_file.exists():
                        # Look for existing prediction file
                        prediction_file = self.pipeline_output_dir / "gene_prediction" / "predictions_genome_coords.gff"
                        if not prediction_file.exists():
                            prediction_file = self.pipeline_output_dir / "gene_prediction" / "predictions.gff"
                        if not prediction_file.exists():
                            raise PipelineError("Gene prediction file required for evidence integration")
                    
                    integration_results = self.run_evidence_integration(
                        [alignment_file, prediction_file], 
                        genome_file=genome_path,
                        nlr_candidates_file=nlr_candidates_file,
                        **kwargs
                    )
            
            # Pipeline completed successfully
            self.results.status = "completed"
            self.results.end_time = datetime.now()
            
            # Save results
            results_file = self.pipeline_output_dir / "pipeline_results.json"
            self.results.save_to_file(results_file)
            
            self.logger.info(f"Pipeline resumed from {start_stage} and completed successfully")
            return self.results
            
        except Exception as e:
            self.results.status = "failed"
            self.results.end_time = datetime.now()
            error_msg = f"Pipeline resume from {start_stage} failed: {e}"
            self.logger.error(error_msg)
            self.results.errors.append(error_msg)
            
            # Save error results
            results_file = self.pipeline_output_dir / "pipeline_results.json"
            self.results.save_to_file(results_file)
            
            raise PipelineError(error_msg, pipeline_stage=f"resume_from_{start_stage}") from e

    def get_results(self) -> PipelineResults:
        """Get current pipeline results"""
        return self.results