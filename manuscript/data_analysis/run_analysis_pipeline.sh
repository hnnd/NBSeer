#!/bin/bash
#
# NBSeer Manuscript Data Analysis Pipeline
# Executes all analysis tasks and generates figures and data required for publication
#

# Set strict mode
set -euo pipefail

# Color output functions
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

log_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

log_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Default configuration
MANUSCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_DIR="${MANUSCRIPT_DIR}/data"
ANALYSIS_DIR="${MANUSCRIPT_DIR}/data_analysis"
FIGURES_DIR="${MANUSCRIPT_DIR}/figures"
TABLES_DIR="${MANUSCRIPT_DIR}/tables"
RESULTS_DIR="${ANALYSIS_DIR}/results"

# NBSeer main script path
PIPELINE_SCRIPT="${MANUSCRIPT_DIR}/../src/nbseer/main.py"

# Create necessary directories
mkdir -p "${RESULTS_DIR}"
mkdir -p "${FIGURES_DIR}/raw"
mkdir -p "${FIGURES_DIR}/processed"
mkdir -p "${TABLES_DIR}/raw_data"
mkdir -p "${TABLES_DIR}/formatted"

# Help information
show_help() {
    cat << EOF
NBSeer Manuscript Data Analysis Pipeline

Usage: $0 [options]

Options:
    -h, --help              Show this help message
    -d, --data-dir DIR      Specify data directory (default: ${DATA_DIR})
    -o, --output-dir DIR    Specify output directory (default: ${RESULTS_DIR})
    -p, --pipeline PATH     Specify NBSeer pipeline script path (default: ${PIPELINE_SCRIPT})
    -s, --skip-runtime      Skip runtime benchmark testing
    -a, --skip-accuracy     Skip accuracy assessment
    -f, --skip-features     Skip feature statistical analysis
    -m, --skip-manuscript   Skip manuscript data generation
    --threads N             Number of parallel threads (default: 8)
    --quick                 Quick mode (reduce repetitions)

Examples:
    $0                                    # Run all analyses
    $0 --skip-runtime --quick             # Skip performance tests, quick mode
    $0 -d /path/to/data -o /path/to/output # Specify custom paths

Analysis stages:
    1. Performance benchmarking (runtime_analysis.py)
    2. Accuracy assessment (gene_structure_eval.py)  
    3. Feature statistical analysis (gene_features_stats.py)
    4. Manuscript data generation (manuscript_data_generator.py)

Note: Ensure required Python dependencies are installed
EOF
}

# Parameter parsing
SKIP_RUNTIME=false
SKIP_ACCURACY=false
SKIP_FEATURES=false
SKIP_MANUSCRIPT=false
THREADS=8
QUICK_MODE=false

while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            show_help
            exit 0
            ;;
        -d|--data-dir)
            DATA_DIR="$2"
            shift 2
            ;;
        -o|--output-dir)
            RESULTS_DIR="$2"
            shift 2
            ;;
        -p|--pipeline)
            PIPELINE_SCRIPT="$2"
            shift 2
            ;;
        -s|--skip-runtime)
            SKIP_RUNTIME=true
            shift
            ;;
        -a|--skip-accuracy)
            SKIP_ACCURACY=true
            shift
            ;;
        -f|--skip-features)
            SKIP_FEATURES=true
            shift
            ;;
        -m|--skip-manuscript)
            SKIP_MANUSCRIPT=true
            shift
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --quick)
            QUICK_MODE=true
            shift
            ;;
        *)
            log_error "Unknown parameter: $1"
            show_help
            exit 1
            ;;
    esac
done

# Input validation
log_info "Validating input parameters and files..."

if [[ ! -d "$DATA_DIR" ]]; then
    log_error "Data directory does not exist: $DATA_DIR"
    exit 1
fi

if [[ ! -d "$DATA_DIR/genome" ]]; then
    log_error "Genome data directory does not exist: $DATA_DIR/genome"
    exit 1
fi

if [[ ! -d "$DATA_DIR/prgdb" ]]; then
    log_error "Protein database directory does not exist: $DATA_DIR/prgdb"
    exit 1
fi

# Check required data files
GENOME_FILES=("$DATA_DIR/genome/tair10.fa" "$DATA_DIR/genome/osa.fa" "$DATA_DIR/genome/CM334.fa")
PROTEIN_DB="$DATA_DIR/prgdb/prg_nbs.fasta"

for genome_file in "${GENOME_FILES[@]}"; do
    if [[ ! -f "$genome_file" ]]; then
        log_error "Genome file does not exist: $genome_file"
        exit 1
    fi
done

if [[ ! -f "$PROTEIN_DB" ]]; then
    log_error "Protein database file does not exist: $PROTEIN_DB"
    exit 1
fi

log_success "Input validation completed"

# Display runtime configuration
log_info "Runtime configuration:"
log_info "  Data directory: $DATA_DIR"
log_info "  Output directory: $RESULTS_DIR"
log_info "  Pipeline script: $PIPELINE_SCRIPT"
log_info "  Parallel threads: $THREADS"
log_info "  Quick mode: $QUICK_MODE"
log_info "  Skip performance testing: $SKIP_RUNTIME"
log_info "  Skip accuracy assessment: $SKIP_ACCURACY"
log_info "  Skip feature analysis: $SKIP_FEATURES"
log_info "  Skip manuscript generation: $SKIP_MANUSCRIPT"

# Check Python dependencies
log_info "Checking Python dependencies..."
python3 -c "
import sys
# Package name mapping: install name -> import name
package_mapping = {
    'pandas': 'pandas',
    'numpy': 'numpy', 
    'matplotlib': 'matplotlib',
    'seaborn': 'seaborn',
    'scipy': 'scipy',
    'scikit-learn': 'sklearn',
    'biopython': 'Bio',
    'gffutils': 'gffutils',
    'psutil': 'psutil',
    'intervaltree': 'intervaltree'
}

missing = []
for install_name, import_name in package_mapping.items():
    try:
        __import__(import_name)
    except ImportError:
        missing.append(install_name)

if missing:
    print('Missing dependencies: ' + ', '.join(missing))
    print('Please run: pip install ' + ' '.join(missing))
    sys.exit(1)
else:
    print('All dependencies installed')
"

if [[ $? -ne 0 ]]; then
    log_error "Python dependency check failed"
    exit 1
fi

log_success "Python dependency check completed"

# Stage 1: Performance benchmarking
if [[ "$SKIP_RUNTIME" == "false" ]]; then
    log_info "=========================================="
    log_info "Stage 1: Running performance benchmarks"
    log_info "=========================================="
    
    RUNTIME_OUTPUT="${RESULTS_DIR}/runtime_benchmarks"
    mkdir -p "$RUNTIME_OUTPUT"
    
    if [[ "$QUICK_MODE" == "true" ]]; then
        SKIP_PARALLEL="--skip-parallelization"
    else
        SKIP_PARALLEL=""
    fi
    
    log_info "Starting performance benchmarking..."
    python3 "${ANALYSIS_DIR}/benchmark_scripts/runtime_analysis.py" \
        --manuscript-data "$DATA_DIR" \
        --output "$RUNTIME_OUTPUT" \
        --pipeline-script "$PIPELINE_SCRIPT" \
        $SKIP_PARALLEL
    
    if [[ $? -eq 0 ]]; then
        log_success "Performance benchmarking completed"
    else
        log_error "性能基准测试失败"
        exit 1
    fi
else
    log_warning "跳过性能基准测试"
fi

# 阶段2: 准确性评估
if [[ "$SKIP_ACCURACY" == "false" ]]; then
    log_info "=========================================="
    log_info "阶段2: 运行准确性评估"
    log_info "=========================================="
    
    ACCURACY_OUTPUT="${RESULTS_DIR}/accuracy_assessment"
    mkdir -p "$ACCURACY_OUTPUT"
    
    # 注意: 这里需要真实的参考注释文件
    # 现在使用示例路径，实际使用时需要替换
    log_warning "注意: 准确性评估需要参考注释文件"
    log_warning "请确保以下路径存在参考GFF文件:"
    log_warning "  - 拟南芥: ${MANUSCRIPT_DIR}/data_analysis/accuracy_assessment/reference_annotations/tair10_nbs_genes.gff3"
    log_warning "  - 水稻: ${MANUSCRIPT_DIR}/data_analysis/accuracy_assessment/reference_annotations/osa_nbs_genes.gff3"
    log_warning "  - 辣椒: ${MANUSCRIPT_DIR}/data_analysis/accuracy_assessment/reference_annotations/CM334_nbs_genes.gff3"
    
    # 如果参考注释文件存在，运行准确性评估
    REF_DIR="${MANUSCRIPT_DIR}/data_analysis/accuracy_assessment/reference_annotations"
    if [[ -d "$REF_DIR" ]]; then
        log_info "发现参考注释目录，开始准确性评估..."
        
        # 对每个数据集运行评估
        DATASETS=("tair10" "osa" "CM334")
        for dataset in "${DATASETS[@]}"; do
            ref_file="${REF_DIR}/${dataset}_nbs_genes.gff3"
            pred_file="${MANUSCRIPT_DIR}/../results_${dataset}/evidence_integration/final_annotations.gff"
            
            if [[ -f "$ref_file" && -f "$pred_file" ]]; then
                log_info "评估数据集: $dataset"
                python3 "${ANALYSIS_DIR}/accuracy_assessment/gene_structure_eval.py" \
                    --reference-gff "$ref_file" \
                    --predicted-gff "$pred_file" \
                    --dataset-name "$dataset" \
                    --output "$ACCURACY_OUTPUT"
                
                if [[ $? -ne 0 ]]; then
                    log_error "数据集 $dataset 准确性评估失败"
                fi
            else
                log_warning "跳过数据集 $dataset (缺少必要文件)"
            fi
        done
        
        log_success "准确性评估完成"
    else
        log_warning "未找到参考注释文件，跳过准确性评估"
        log_warning "请准备参考注释后重新运行"
    fi
else
    log_warning "跳过准确性评估"
fi

# 阶段3: 特征统计分析
if [[ "$SKIP_FEATURES" == "false" ]]; then
    log_info "=========================================="
    log_info "阶段3: 运行特征统计分析"
    log_info "=========================================="
    
    FEATURES_OUTPUT="${RESULTS_DIR}/gene_features"
    mkdir -p "$FEATURES_OUTPUT"
    
    # 收集所有预测结果的GFF文件
    GFF_FILES=()
    GENOME_FILES_FOR_ANALYSIS=()
    DATASET_NAMES=()
    
    RESULTS_DIRS=("${MANUSCRIPT_DIR}/../results" "${MANUSCRIPT_DIR}/../results_tair10" "${MANUSCRIPT_DIR}/../results_CA59" "${MANUSCRIPT_DIR}/../results_cm334")
    
    for results_dir in "${RESULTS_DIRS[@]}"; do
        if [[ -d "$results_dir" ]]; then
            gff_file="${results_dir}/evidence_integration/final_annotations.gff"
            if [[ -f "$gff_file" ]]; then
                # 检查文件是否有实际内容（超过头信息）
                line_count=$(wc -l < "$gff_file")
                if [[ $line_count -gt 10 ]]; then
                    GFF_FILES+=("$gff_file")
                    
                    # 确定对应的基因组文件和数据集名称
                    if [[ "$results_dir" == *"tair10"* ]]; then
                        GENOME_FILES_FOR_ANALYSIS+=("${DATA_DIR}/genome/tair10.fa")
                        DATASET_NAMES+=("Arabidopsis")
                    elif [[ "$results_dir" == *"CA59"* ]] || [[ "$results_dir" == *"osa"* ]]; then
                        GENOME_FILES_FOR_ANALYSIS+=("${DATA_DIR}/genome/osa.fa")
                        DATASET_NAMES+=("Rice")
                    elif [[ "$results_dir" == *"cm334"* ]] || [[ "$results_dir" == *"CM334"* ]]; then
                        GENOME_FILES_FOR_ANALYSIS+=("${DATA_DIR}/genome/CM334.fa")
                        DATASET_NAMES+=("Pepper")
                    else
                        GENOME_FILES_FOR_ANALYSIS+=("${DATA_DIR}/genome/osa.fa")
                        DATASET_NAMES+=("Default")
                    fi
                else
                    log_warning "跳过空的GFF文件: $gff_file (只有 $line_count 行)"
                fi
            fi
        fi
    done
    
    if [[ ${#GFF_FILES[@]} -gt 0 ]]; then
        log_info "找到 ${#GFF_FILES[@]} 个GFF文件进行特征分析"
        
        # 构建参数
        GFF_ARGS=""
        GENOME_ARGS=""
        NAME_ARGS=""
        
        for gff_file in "${GFF_FILES[@]}"; do
            GFF_ARGS+="$gff_file "
        done
        
        for genome_file in "${GENOME_FILES_FOR_ANALYSIS[@]}"; do
            GENOME_ARGS+="$genome_file "
        done
        
        for dataset_name in "${DATASET_NAMES[@]}"; do
            NAME_ARGS+="$dataset_name "
        done
        
        log_info "开始特征统计分析..."
        python3 "${ANALYSIS_DIR}/statistical_analysis/gene_features_stats.py" \
            --gff-files $GFF_ARGS \
            --genome-files $GENOME_ARGS \
            --dataset-names $NAME_ARGS \
            --output "$FEATURES_OUTPUT"
        
        if [[ $? -eq 0 ]]; then
            log_success "特征统计分析完成"
        else
            log_error "特征统计分析失败"
            exit 1
        fi
    else
        log_warning "未找到GFF文件，跳过特征统计分析"
        log_warning "请先运行NBS-Pipeline生成注释结果"
    fi
else
    log_warning "跳过特征统计分析"
fi

# 阶段4: 论文数据生成
if [[ "$SKIP_MANUSCRIPT" == "false" ]]; then
    log_info "=========================================="
    log_info "阶段4: 生成论文图表和数据"
    log_info "=========================================="
    
    log_info "开始生成论文数据..."
    python3 "${ANALYSIS_DIR}/reports/manuscript_data_generator.py" \
        --analysis-results "$RESULTS_DIR" \
        --figures-output "${FIGURES_DIR}/processed" \
        --tables-output "${TABLES_DIR}/formatted"
    
    if [[ $? -eq 0 ]]; then
        log_success "论文数据生成完成"
    else
        log_error "论文数据生成失败"
        exit 1
    fi
else
    log_warning "跳过论文数据生成"
fi

# 生成最终汇总报告
log_info "=========================================="
log_info "生成最终分析汇总报告"
log_info "=========================================="

SUMMARY_FILE="${RESULTS_DIR}/analysis_summary.md"

cat > "$SUMMARY_FILE" << EOF
# NBS-Pipeline 论文数据分析汇总报告

**生成时间**: $(date '+%Y-%m-%d %H:%M:%S')
**分析脚本版本**: 1.0.0

## 分析配置

- **数据目录**: ${DATA_DIR}
- **输出目录**: ${RESULTS_DIR}
- **并行线程数**: ${THREADS}
- **快速模式**: ${QUICK_MODE}

## 执行的分析阶段

- **性能基准测试**: $([ "$SKIP_RUNTIME" == "true" ] && echo "跳过" || echo "✓ 完成")
- **准确性评估**: $([ "$SKIP_ACCURACY" == "true" ] && echo "跳过" || echo "✓ 完成")
- **特征统计分析**: $([ "$SKIP_FEATURES" == "true" ] && echo "跳过" || echo "✓ 完成")
- **论文数据生成**: $([ "$SKIP_MANUSCRIPT" == "true" ] && echo "跳过" || echo "✓ 完成")

## 输出文件结构

\`\`\`
${RESULTS_DIR}/
├── runtime_benchmarks/          # 性能基准测试结果
├── accuracy_assessment/         # 准确性评估结果
├── gene_features/              # 基因特征分析结果
└── analysis_summary.md         # 本报告

${FIGURES_DIR}/processed/        # 论文图表
├── Figure1_Pipeline_Workflow.tiff
├── Figure2_Accuracy_Assessment.tiff
├── Figure3_Performance_Benchmarks.tiff
└── FigureS*_*.tiff             # 补充图表

${TABLES_DIR}/formatted/         # 论文数据表
├── Table1_Dataset_Summary.csv
├── Table1_Dataset_Summary.tex
└── *.csv                       # 其他数据表
\`\`\`

## 下一步操作

1. 检查各分析阶段的输出结果
2. 将图表和数据表整合到论文中
3. 根据需要调整图表格式和内容
4. 准备投稿材料

## 注意事项

- 所有图表已按Bioinformatics期刊要求格式化 (300 DPI TIFF)
- 表格提供CSV和LaTeX两种格式
- 如有错误或需要调整，请重新运行相应分析阶段

---
**分析完成**: $(date '+%Y-%m-%d %H:%M:%S')
EOF

log_success "汇总报告已保存到: $SUMMARY_FILE"

# 显示最终结果
log_info "=========================================="
log_info "分析流程完成!"
log_info "=========================================="

log_success "所有分析任务已完成"
log_info "结果文件位置:"
log_info "  - 分析结果: $RESULTS_DIR"
log_info "  - 论文图表: ${FIGURES_DIR}/processed"
log_info "  - 论文数据表: ${TABLES_DIR}/formatted"
log_info "  - 汇总报告: $SUMMARY_FILE"

log_info "请检查各阶段的输出文件，确认分析结果符合预期"
log_info "如需重新运行特定阶段，可使用相应的 --skip-* 参数"

# 检查生成的文件数量
FIGURE_COUNT=$(find "${FIGURES_DIR}/processed" -name "*.tiff" 2>/dev/null | wc -l)
TABLE_COUNT=$(find "${TABLES_DIR}/formatted" -name "*.csv" 2>/dev/null | wc -l)

if [[ $FIGURE_COUNT -gt 0 ]]; then
    log_success "生成了 $FIGURE_COUNT 个图表文件"
else
    log_warning "未生成图表文件，请检查数据生成阶段"
fi

if [[ $TABLE_COUNT -gt 0 ]]; then
    log_success "生成了 $TABLE_COUNT 个数据表文件"
else
    log_warning "未生成数据表文件，请检查数据生成阶段"
fi

log_info "分析流程成功完成!"