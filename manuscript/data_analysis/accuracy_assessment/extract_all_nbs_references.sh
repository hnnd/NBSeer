#!/bin/bash

# 批量提取所有基因组的NBS参考基因

set -e  # 遇到错误时退出

# 设置颜色输出
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${GREEN}=== NBS Reference Genes Extraction Pipeline ===${NC}"
echo "Starting extraction for all genomes..."

# 工作目录
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="/data/usdata/wangys/work/nbs"
DATA_DIR="$PROJECT_ROOT/manuscript/data/genome"
OUTPUT_DIR="$SCRIPT_DIR/reference_genes"

# 创建输出目录
mkdir -p "$OUTPUT_DIR"

# 检查输入文件
echo -e "${YELLOW}Checking input files...${NC}"

# 定义数据集
declare -A DATASETS=(
    ["arabidopsis"]="tair10"
    ["rice"]="osa" 
    ["pepper"]="CM334"
)

# 检查所有需要的文件是否存在
all_files_exist=true
for dataset in "${!DATASETS[@]}"; do
    file_prefix="${DATASETS[$dataset]}"
    genome_file="$DATA_DIR/${file_prefix}.fa"
    gff_file="$DATA_DIR/${file_prefix}.gff"
    
    if [[ ! -f "$genome_file" ]]; then
        echo -e "${RED}Error: Genome file not found: $genome_file${NC}"
        all_files_exist=false
    fi
    
    if [[ ! -f "$gff_file" ]]; then
        echo -e "${RED}Error: GFF file not found: $gff_file${NC}"
        all_files_exist=false
    fi
done

if [[ "$all_files_exist" != true ]]; then
    echo -e "${RED}Missing required files. Please check the data directory.${NC}"
    exit 1
fi

echo -e "${GREEN}All input files found.${NC}"

# 运行提取流程
extraction_script="$SCRIPT_DIR/extract_nbs_reference_genes.py"

if [[ ! -f "$extraction_script" ]]; then
    echo -e "${RED}Error: Extraction script not found: $extraction_script${NC}"
    exit 1
fi

echo -e "${YELLOW}Starting NBS gene extraction...${NC}"

# 处理每个数据集
for dataset in "${!DATASETS[@]}"; do
    file_prefix="${DATASETS[$dataset]}"
    genome_file="$DATA_DIR/${file_prefix}.fa"
    gff_file="$DATA_DIR/${file_prefix}.gff"
    dataset_output_dir="$OUTPUT_DIR/$dataset"
    
    echo -e "${GREEN}Processing $dataset ($file_prefix)...${NC}"
    
    # 创建数据集特定的输出目录
    mkdir -p "$dataset_output_dir"
    
    # 运行提取脚本
    if python "$extraction_script" \
        --genome "$genome_file" \
        --gff "$gff_file" \
        --output "$dataset_output_dir" \
        --dataset-name "$dataset"; then
        echo -e "${GREEN}✓ $dataset extraction completed successfully${NC}"
    else
        echo -e "${RED}✗ $dataset extraction failed${NC}"
        continue
    fi
    
    # 显示结果摘要
    if [[ -f "$dataset_output_dir/${dataset}_nbs_extraction_report.json" ]]; then
        nbs_count=$(grep -o '"nbs_genes_identified": [0-9]*' "$dataset_output_dir/${dataset}_nbs_extraction_report.json" | grep -o '[0-9]*')
        echo -e "  NBS genes identified: ${GREEN}$nbs_count${NC}"
    fi
done

echo -e "${GREEN}=== Extraction Summary ===${NC}"

# 生成汇总报告
summary_file="$OUTPUT_DIR/extraction_summary.txt"
echo "NBS Reference Genes Extraction Summary" > "$summary_file"
echo "Generated on: $(date)" >> "$summary_file"
echo "=======================================" >> "$summary_file"
echo "" >> "$summary_file"

total_nbs_genes=0

for dataset in "${!DATASETS[@]}"; do
    dataset_output_dir="$OUTPUT_DIR/$dataset"
    report_file="$dataset_output_dir/${dataset}_nbs_extraction_report.json"
    
    if [[ -f "$report_file" ]]; then
        nbs_count=$(grep -o '"nbs_genes_identified": [0-9]*' "$report_file" | grep -o '[0-9]*')
        total_genes=$(grep -o '"total_genes_parsed": [0-9]*' "$report_file" | grep -o '[0-9]*')
        proteins=$(grep -o '"proteins_translated": [0-9]*' "$report_file" | grep -o '[0-9]*')
        
        echo "$dataset:" >> "$summary_file"
        echo "  Total genes parsed: $total_genes" >> "$summary_file"
        echo "  Proteins translated: $proteins" >> "$summary_file"
        echo "  NBS genes identified: $nbs_count" >> "$summary_file"
        echo "  NBS percentage: $(echo "scale=2; $nbs_count * 100 / $total_genes" | bc -l)%" >> "$summary_file"
        echo "" >> "$summary_file"
        
        total_nbs_genes=$((total_nbs_genes + nbs_count))
        
        echo -e "${dataset}: ${GREEN}$nbs_count${NC} NBS genes (${YELLOW}$(echo "scale=1; $nbs_count * 100 / $total_genes" | bc -l)%${NC})"
    else
        echo -e "${RED}$dataset: Report file not found${NC}"
        echo "$dataset: FAILED" >> "$summary_file"
        echo "" >> "$summary_file"
    fi
done

echo "" >> "$summary_file"
echo "Total NBS genes identified: $total_nbs_genes" >> "$summary_file"

echo -e "\nTotal NBS genes identified: ${GREEN}$total_nbs_genes${NC}"
echo -e "Detailed summary saved to: ${YELLOW}$summary_file${NC}"

# 创建用于准确性评估的文件列表
reference_list_file="$OUTPUT_DIR/nbs_reference_files.txt"
echo "# NBS Reference Gene Files for Accuracy Assessment" > "$reference_list_file"
echo "# Format: dataset_name:gff_file:gene_list_file" >> "$reference_list_file"
echo "# Generated on: $(date)" >> "$reference_list_file"

for dataset in "${!DATASETS[@]}"; do
    dataset_output_dir="$OUTPUT_DIR/$dataset"
    gff_file="$dataset_output_dir/${dataset}_nbs_genes.gff"
    list_file="$dataset_output_dir/${dataset}_nbs_gene_list.txt"
    
    if [[ -f "$gff_file" ]] && [[ -f "$list_file" ]]; then
        echo "$dataset:$gff_file:$list_file" >> "$reference_list_file"
    fi
done

echo -e "Reference files list: ${YELLOW}$reference_list_file${NC}"

echo -e "\n${GREEN}=== NBS Reference Extraction Completed Successfully! ===${NC}"
echo -e "Output directory: ${YELLOW}$OUTPUT_DIR${NC}"
echo ""
echo "Next steps for accuracy assessment:"
echo "1. Review the extracted NBS gene lists"
echo "2. Use the reference GFF files for comparison with pipeline predictions"
echo "3. Run the accuracy assessment analysis"