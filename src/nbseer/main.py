"""
Main entry point for NBSeer - Plant NBS gene annotation tool

Provides command-line interface and main execution logic
"""

import argparse
import sys
import time
from pathlib import Path

from .data.validation import DataValidator
from .pipeline.coordinator import NBSAnnotationPipeline
from .utils.config import load_config, load_config_with_overrides
from .utils.exceptions import NBSAnnotationError
from .utils.logging_setup import get_logger, setup_logging


def setup_argument_parser() -> argparse.ArgumentParser:
    """
    Set up command line argument parser
    设置命令行参数解析器
    
    Returns:
        Configured ArgumentParser
    """
    parser = argparse.ArgumentParser(
        description="NBS Gene Annotation Pipeline - 植物NBS抗病基因重新注释流水线",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples / 示例:
  # Run full pipeline / 运行完整流水线
  %(prog)s --genome genome.fa --proteins proteins.fa --output results/

  # Run with custom configuration / 使用自定义配置运行
  %(prog)s --genome genome.fa --proteins proteins.fa --config custom.yaml

  # Run specific stage only / 仅运行特定阶段
  %(prog)s --genome genome.fa --stage nlr_localization

  # Validate input files only / 仅验证输入文件
  %(prog)s --genome genome.fa --proteins proteins.fa --validate-only

  # Check tool availability / 检查工具可用性
  %(prog)s --check-tools
        """,
    )

    # Input files / 输入文件
    input_group = parser.add_argument_group("Input Files", "输入文件")
    input_group.add_argument(
        "--genome", "-g",
        type=Path,
        help="Path to genome FASTA file / 基因组FASTA文件路径",
    )
    input_group.add_argument(
        "--proteins", "-p",
        type=Path,
        help="Path to protein FASTA file / 蛋白质FASTA文件路径",
    )
    input_group.add_argument(
        "--nlr-candidates",
        type=Path,
        help="Path to pre-computed NLR candidates file / 预计算的NLR候选基因文件路径",
    )

    # Output options / 输出选项
    output_group = parser.add_argument_group("Output Options", "输出选项")
    output_group.add_argument(
        "--output", "-o",
        type=Path,
        default=Path("output"),
        help="Output directory / 输出目录 (default: output)",
    )
    output_group.add_argument(
        "--pipeline-id",
        type=str,
        help="Unique pipeline identifier / 唯一流水线标识符",
    )

    # Configuration / 配置
    config_group = parser.add_argument_group("Configuration", "配置")
    config_group.add_argument(
        "--config", "-c",
        type=Path,
        help="Path to configuration file / 配置文件路径",
    )
    config_group.add_argument(
        "--override-config",
        type=Path,
        help="Path to override configuration file / 覆盖配置文件路径",
    )

    # Pipeline control / 流水线控制
    pipeline_group = parser.add_argument_group("Pipeline Control", "流水线控制")
    pipeline_group.add_argument(
        "--stage",
        choices=["nlr_localization", "protein_alignment", "augustus_training", 
                "gene_prediction", "evidence_integration", "all"],
        default="all",
        help="Pipeline stage to run / 要运行的流水线阶段 (default: all)",
    )
    pipeline_group.add_argument(
        "--resume-from",
        choices=["nlr_localization", "protein_alignment", "augustus_training", 
                "gene_prediction", "evidence_integration"],
        help="Resume pipeline from specific stage / 从特定阶段恢复流水线",
    )
    pipeline_group.add_argument(
        "--skip-validation",
        action="store_true",
        help="Skip input file validation / 跳过输入文件验证",
    )

    # Tool options / 工具选项
    tool_group = parser.add_argument_group("Tool Options", "工具选项")
    tool_group.add_argument(
        "--augustus-model",
        type=str,
        default=None,
        help="Augustus species model / Augustus物种模型 (default: from config)",
    )
    tool_group.add_argument(
        "--threads", "-t",
        type=int,
        default=8,
        help="Number of threads to use / 使用的线程数 (default: 8)",
    )
    tool_group.add_argument(
        "--memory",
        type=int,
        default=32,
        help="Maximum memory in GB / 最大内存GB (default: 32)",
    )
    
    # Augustus training options / Augustus训练选项
    training_group = parser.add_argument_group("Augustus Training", "Augustus训练")
    training_group.add_argument(
        "--enable-training",
        action="store_true",
        help="Enable Augustus model training with Miniprot results / 启用基于Miniprot结果的Augustus模型训练",
    )
    training_group.add_argument(
        "--training-species-name",
        type=str,
        help="Name for the new trained Augustus species model / 新训练的Augustus物种模型名称",
    )
    training_group.add_argument(
        "--training-quality",
        choices=["high", "medium", "low", "all"],
        default="high",
        help="Quality level of Miniprot results for training / 用于训练的Miniprot结果质量级别 (default: high)",
    )
    training_group.add_argument(
        "--training-min-genes",
        type=int,
        default=20,
        help="Minimum number of training genes required / 所需的最少训练基因数 (default: 20)",
    )
    training_group.add_argument(
        "--training-optimization-rounds",
        type=int,
        default=1,
        help="Number of optimization rounds for training / 训练的优化轮数 (default: 1)",
    )
    training_group.add_argument(
        "--training-timeout",
        type=int,
        default=240,
        help="Training timeout in minutes / 训练超时时间（分钟） (default: 240)",
    )
    training_group.add_argument(
        "--skip-training-optimization",
        action="store_true",
        help="Skip parameter optimization during training (faster) / 跳过训练期间的参数优化（更快）",
    )

    # Validation and debugging / 验证和调试
    debug_group = parser.add_argument_group("Validation & Debugging", "验证和调试")
    debug_group.add_argument(
        "--validate-only",
        action="store_true",
        help="Only validate input files, don't run pipeline / 仅验证输入文件，不运行流水线",
    )
    debug_group.add_argument(
        "--check-tools",
        action="store_true",
        help="Check tool availability and exit / 检查工具可用性并退出",
    )
    debug_group.add_argument(
        "--verbose", "-v",
        action="count",
        default=0,
        help="Increase verbosity / 增加详细程度 (use -vv for debug)",
    )
    debug_group.add_argument(
        "--log-file",
        type=Path,
        help="Path to log file / 日志文件路径",
    )
    debug_group.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be done without executing / 显示将要执行的操作但不实际执行",
    )

    return parser


def validate_arguments(args: argparse.Namespace) -> None:
    """
    Validate command line arguments
    验证命令行参数
    
    Args:
        args: Parsed command line arguments
        
    Raises:
        SystemExit: If validation fails
    """
    errors = []

    # Check required inputs for pipeline execution
    if not args.check_tools and not args.validate_only:
        if args.stage in ["all", "nlr_localization", "protein_alignment"] and not args.genome:
            errors.append("Genome file (--genome) is required")

        if args.stage in ["all", "protein_alignment"] and not args.proteins:
            errors.append("Protein file (--proteins) is required")
        
        # Special validation for --resume-from
        if args.resume_from:
            if args.stage != "all":
                errors.append("--resume-from can only be used with --stage all")
            if not args.genome:
                errors.append("Genome file (--genome) is required when using --resume-from")
            if not args.proteins:
                errors.append("Protein file (--proteins) is required when using --resume-from")

    # Check file existence
    for file_arg, file_path in [
        ("--genome", args.genome),
        ("--proteins", args.proteins),
        ("--nlr-candidates", args.nlr_candidates),
        ("--config", args.config),
        ("--override-config", args.override_config),
    ]:
        if file_path and not file_path.exists():
            errors.append(f"File not found for {file_arg}: {file_path}")

    # Check output directory
    if args.output:
        try:
            args.output.mkdir(parents=True, exist_ok=True)
        except OSError as e:
            errors.append(f"Cannot create output directory {args.output}: {e}")

    # Check numeric arguments
    if args.threads <= 0:
        errors.append("Number of threads must be positive")

    if args.memory <= 0:
        errors.append("Memory limit must be positive")

    if errors:
        print("❌ Argument validation failed:", file=sys.stderr)
        for error in errors:
            print(f"  - {error}", file=sys.stderr)
        sys.exit(1)


def setup_logging_from_args(args: argparse.Namespace) -> None:
    """
    Set up logging based on command line arguments
    根据命令行参数设置日志
    
    Args:
        args: Parsed command line arguments
    """
    # Determine log level
    if args.verbose == 0:
        log_level = "INFO"
    elif args.verbose == 1:
        log_level = "DEBUG"
    else:
        log_level = "DEBUG"  # Very verbose

    # Set up logging
    setup_logging(
        level=log_level,
        log_file=args.log_file,
        console_output=True,
    )


def check_tool_availability(config_path: Path | None = None) -> bool:
    """
    Check availability of all required tools
    检查所有必需工具的可用性
    
    Args:
        config_path: Path to configuration file
        
    Returns:
        True if all tools are available, False otherwise
    """
    logger = get_logger(__name__)

    try:
        config = load_config(config_path)
        availability = config.validate_tool_availability()

        print("🔧 Tool Availability Check / 工具可用性检查:")
        print("=" * 50)

        all_available = True
        for tool_name, is_available in availability.items():
            status = "✅ Available" if is_available else "❌ Not Found"
            print(f"{tool_name:20} {status}")
            if not is_available:
                all_available = False

        print("=" * 50)
        if all_available:
            print("✅ All tools are available / 所有工具都可用")
        else:
            print("❌ Some tools are missing / 某些工具缺失")
            print("Please install missing tools or update configuration")

        return all_available

    except Exception as e:
        logger.error(f"Tool availability check failed: {e}")
        return False


def validate_input_files(args: argparse.Namespace) -> bool:
    """
    Validate input files
    验证输入文件
    
    Args:
        args: Parsed command line arguments
        
    Returns:
        True if validation passes, False otherwise
    """
    logger = get_logger(__name__)

    try:
        validator = DataValidator(strict_mode=False)

        print("📁 Input File Validation / 输入文件验证:")
        print("=" * 50)

        validation_results = {}
        all_valid = True

        # Validate genome file
        if args.genome:
            print(f"Validating genome file: {args.genome}")
            result = validator.validate_genome_fasta(args.genome)
            validation_results["genome"] = result

            if result.is_valid:
                print(f"  ✅ Valid genome file ({result.record_count} sequences)")
            else:
                print(f"  ❌ Invalid genome file: {result.error_message}")
                all_valid = False

        # Validate protein file
        if args.proteins:
            print(f"Validating protein file: {args.proteins}")
            result = validator.validate_protein_fasta(args.proteins)
            validation_results["proteins"] = result

            if result.is_valid:
                print(f"  ✅ Valid protein file ({result.record_count} sequences)")
            else:
                print(f"  ❌ Invalid protein file: {result.error_message}")
                all_valid = False

        # Validate NLR candidates file
        if args.nlr_candidates:
            print(f"Validating NLR candidates file: {args.nlr_candidates}")
            result = validator.validate_gff(args.nlr_candidates)
            validation_results["nlr_candidates"] = result

            if result.is_valid:
                print(f"  ✅ Valid GFF file ({result.record_count} features)")
            else:
                print(f"  ❌ Invalid GFF file: {result.error_message}")
                all_valid = False

        print("=" * 50)
        if all_valid:
            print("✅ All input files are valid / 所有输入文件都有效")
        else:
            print("❌ Some input files are invalid / 某些输入文件无效")

        # Save validation report
        if validation_results:
            report_path = args.output / "input_validation_report.md"
            report = validator.generate_validation_report(validation_results, report_path)
            print(f"📄 Validation report saved to: {report_path}")

        return all_valid

    except Exception as e:
        logger.error(f"Input validation failed: {e}")
        return False


def run_pipeline(args: argparse.Namespace) -> bool:
    """
    Run the NBS annotation pipeline
    运行NBS注释流水线
    
    Args:
        args: Parsed command line arguments
        
    Returns:
        True if pipeline succeeds, False otherwise
    """
    logger = get_logger(__name__)

    try:
        # Load configuration
        env_overrides = {
            "pipeline.parallel.max_workers": args.threads,
            "pipeline.memory.max_memory_gb": args.memory,
        }
        
        # Add Augustus training overrides if enabled
        if args.enable_training:
            env_overrides.update({
                "tools.augustus.training.miniprot_training.enabled": True,
                "tools.augustus.training.miniprot_training.quality_filter": args.training_quality,
                "tools.augustus.training.miniprot_training.min_training_genes": args.training_min_genes,
                "tools.augustus.training.miniprot_training.optimization_rounds": 0 if args.skip_training_optimization else args.training_optimization_rounds,
                "tools.augustus.training.miniprot_training.timeout_minutes": args.training_timeout,
            })
            
            if args.training_species_name:
                env_overrides["tools.augustus.training.miniprot_training.species_name"] = args.training_species_name
        
        config = load_config_with_overrides(
            base_config_path=args.config,
            override_config_path=args.override_config,
            env_overrides=env_overrides,
        )

        # Initialize pipeline with absolute output path
        pipeline = NBSAnnotationPipeline(
            config=config,
            output_base_dir=args.output.resolve(),
            pipeline_id=args.pipeline_id,
        )

        print(f"🚀 Starting NBS annotation pipeline: {pipeline.pipeline_id}")
        print(f"📂 Output directory: {pipeline.pipeline_output_dir}")

        if args.dry_run:
            print("🔍 Dry run mode - showing planned execution:")
            print(f"  - Genome file: {args.genome}")
            print(f"  - Protein file: {args.proteins}")
            print(f"  - Pipeline stage: {args.stage}")
            if args.resume_from:
                print(f"  - Resume from stage: {args.resume_from}")
            print(f"  - Augustus model: {args.augustus_model}")
            print(f"  - Augustus training enabled: {args.enable_training}")
            if args.enable_training:
                print(f"  - Training species name: {args.training_species_name}")
                print(f"  - Training quality: {args.training_quality}")
                print(f"  - Training min genes: {args.training_min_genes}")
                print(f"  - Training optimization rounds: {args.training_optimization_rounds}")
                print(f"  - Training timeout: {args.training_timeout} minutes")
            print(f"  - Threads: {args.threads}")
            print(f"  - Memory: {args.memory}GB")
            return True

        # Run pipeline based on stage
        start_time = time.time()

        if args.stage == "all":
            # Run full pipeline or resume from specific stage
            if args.resume_from:
                # Resume from specific stage - run remaining stages
                results = pipeline.run_pipeline_from_stage(
                    start_stage=args.resume_from,
                    genome_path=args.genome,
                    protein_path=args.proteins,
                    augustus_model=args.augustus_model,
                    training_species_name=args.training_species_name,
                )
            else:
                # Run full pipeline from beginning
                results = pipeline.run_full_pipeline(
                    genome_path=args.genome,
                    protein_path=args.proteins,
                    augustus_model=args.augustus_model,
                    training_species_name=args.training_species_name,
                )
        # Run specific stage
        elif args.stage == "nlr_localization":
            results = pipeline.run_nlr_localization(genome_path=args.genome)
        elif args.stage == "protein_alignment":
            results = pipeline.run_protein_alignment(
                genome_path=args.genome,
                protein_path=args.proteins,
            )
        elif args.stage == "augustus_training":
            # Run Augustus training as standalone stage
            results = pipeline.run_augustus_training(
                training_species_name=args.training_species_name,
            )
        elif args.stage == "gene_prediction":
            if not args.nlr_candidates:
                raise NBSAnnotationError("NLR candidates file required for gene prediction stage")
            results = pipeline.run_gene_prediction(
                nlr_candidates_path=args.nlr_candidates,
                training_model=args.augustus_model,
            )
        elif args.stage == "evidence_integration":
            # For evidence integration, we need existing prediction files
            prediction_files = []
            
            # Look for Augustus predictions
            augustus_file = args.output / "gene_prediction" / "predictions_genome_coords.gff"
            if not augustus_file.exists():
                augustus_file = args.output / "gene_prediction" / "predictions.gff"
            if augustus_file.exists():
                prediction_files.append(augustus_file)
            
            # Look for Miniprot alignments
            miniprot_file = args.output / "protein_alignment" / "alignments.gff"
            if miniprot_file.exists():
                prediction_files.append(miniprot_file)
            
            # Look for NLR candidates
            nlr_file = args.nlr_candidates or (args.output / "nlr_localization" / "NLR_candidates.gff")
            
            if not prediction_files:
                raise NBSAnnotationError("No prediction files found for evidence integration. Run previous stages first.")
            
            results = pipeline.run_evidence_integration(
                prediction_files=prediction_files,
                genome_file=args.genome,
                nlr_candidates_file=nlr_file if nlr_file and nlr_file.exists() else None,
            )
        else:
            raise NBSAnnotationError(f"Stage not implemented: {args.stage}")

        execution_time = time.time() - start_time

        # Report results
        print("=" * 50)
        print("✅ Pipeline completed successfully!")
        print(f"⏱️  Execution time: {execution_time:.1f} seconds")

        if hasattr(results, "stages_completed"):
            print(f"📊 Stages completed: {', '.join(results.stages_completed)}")

        if hasattr(results, "output_files"):
            print("📁 Output files:")
            for file_type, file_path in results.output_files.items():
                print(f"  - {file_type}: {file_path}")

        return True

    except NBSAnnotationError as e:
        logger.error(f"Pipeline error: {e}")
        print(f"❌ Pipeline failed: {e}", file=sys.stderr)
        return False
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        print(f"❌ Unexpected error: {e}", file=sys.stderr)
        return False


def main() -> int:
    """
    Main entry point
    主要入口点
    
    Returns:
        Exit code (0 for success, 1 for failure)
    """
    # Parse command line arguments
    parser = setup_argument_parser()
    args = parser.parse_args()

    # Set up logging early
    setup_logging_from_args(args)
    logger = get_logger(__name__)

    try:
        # Validate arguments
        validate_arguments(args)

        # Handle special modes
        if args.check_tools:
            success = check_tool_availability(args.config)
            return 0 if success else 1

        if args.validate_only:
            success = validate_input_files(args)
            return 0 if success else 1

        # Validate inputs unless skipped
        if not args.skip_validation:
            validation_success = validate_input_files(args)
            if not validation_success:
                print("❌ Input validation failed. Use --skip-validation to bypass.")
                return 1

        # Run pipeline
        success = run_pipeline(args)
        return 0 if success else 1

    except KeyboardInterrupt:
        print("\n⚠️  Pipeline interrupted by user", file=sys.stderr)
        logger.info("Pipeline interrupted by user")
        return 1
    except Exception as e:
        print(f"❌ Fatal error: {e}", file=sys.stderr)
        logger.error(f"Fatal error: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
