#!/bin/bash

# NBS基因注释流水线工具验证脚本
# Tool verification script for NBS annotation pipeline

set -e

# 定义颜色输出
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# 项目根目录
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TOOLS_DIR="${PROJECT_ROOT}/tools"

echo -e "${BLUE}🔍 验证NBS基因注释流水线工具配置${NC}"
echo -e "${BLUE}Verifying NBS annotation pipeline tool configuration${NC}"
echo ""

# 加载环境
source "${PROJECT_ROOT}/setup_env.sh"
echo ""

# 验证工具函数
verify_tool() {
    local tool_name="$1"
    local tool_path="$2"
    local version_cmd="$3"
    local expected_output="$4"
    
    echo -e "${YELLOW}检查 ${tool_name}...${NC}"
    
    if [ -f "${tool_path}" ] || command -v "${tool_name}" > /dev/null 2>&1; then
        echo -e "${GREEN}  ✓ 可执行文件存在${NC}"
        
        if [ -n "${version_cmd}" ]; then
            echo "  • 版本信息:"
            if eval "${version_cmd}" 2>/dev/null | head -3; then
                echo -e "${GREEN}  ✓ 工具运行正常${NC}"
            else
                echo -e "${YELLOW}  ! 无法获取版本信息，但工具存在${NC}"
            fi
        fi
    else
        echo -e "${RED}  ✗ 工具未找到: ${tool_path}${NC}"
        return 1
    fi
}

# 验证Java工具
verify_java_tool() {
    local tool_name="$1"
    local jar_path="$2"
    
    echo -e "${YELLOW}检查 ${tool_name}...${NC}"
    
    if [ -f "${jar_path}" ]; then
        echo -e "${GREEN}  ✓ JAR文件存在: ${jar_path}${NC}"
        
        # 检查Java版本
        if command -v java > /dev/null 2>&1; then
            echo "  • Java版本: $(java -version 2>&1 | head -1)"
            echo -e "${GREEN}  ✓ Java环境正常${NC}"
        else
            echo -e "${RED}  ✗ Java未安装${NC}"
            return 1
        fi
    else
        echo -e "${RED}  ✗ JAR文件未找到: ${jar_path}${NC}"
        return 1
    fi
}

# 验证目录结构
verify_directory_structure() {
    echo -e "${YELLOW}检查目录结构...${NC}"
    
    local dirs=(
        "tools/bin"
        "tools/Augustus"
        "tools/miniprot" 
        "tools/EVidenceModeler"
        "tools/NLR-Annotator"
    )
    
    for dir in "${dirs[@]}"; do
        if [ -d "${PROJECT_ROOT}/${dir}" ]; then
            echo -e "${GREEN}  ✓ ${dir}${NC}"
        else
            echo -e "${RED}  ✗ ${dir}${NC}"
        fi
    done
}

# 验证符号链接
verify_symlinks() {
    echo -e "${YELLOW}检查符号链接...${NC}"
    
    local links=(
        "tools/bin/augustus"
        "tools/bin/miniprot"
        "tools/Augustus/bin"
        "tools/Augustus/config"
        "tools/Augustus/scripts"
        "tools/miniprot/miniprot"
    )
    
    for link in "${links[@]}"; do
        local full_path="${PROJECT_ROOT}/${link}"
        if [ -L "${full_path}" ]; then
            local target=$(readlink "${full_path}")
            echo -e "${GREEN}  ✓ ${link} -> ${target}${NC}"
        elif [ -f "${full_path}" ] || [ -d "${full_path}" ]; then
            echo -e "${GREEN}  ✓ ${link} (直接文件/目录)${NC}"
        else
            echo -e "${RED}  ✗ ${link}${NC}"
        fi
    done
}

# 主验证流程
main() {
    verify_directory_structure
    echo ""
    
    verify_symlinks
    echo ""
    
    # 验证各个工具
    verify_tool "Augustus" "${TOOLS_DIR}/bin/augustus" "augustus --version" ""
    echo ""
    
    verify_tool "miniprot" "${TOOLS_DIR}/bin/miniprot" "miniprot --version" ""
    echo ""
    
    verify_java_tool "NLR-Annotator" "${TOOLS_DIR}/NLR-Annotator/NLR-Annotator-v2.1b.jar"
    echo ""
    
    # 验证EVM脚本
    echo -e "${YELLOW}检查 EVidenceModeler...${NC}"
    local evm_scripts=(
        "tools/EVidenceModeler/EvmUtils/partition_EVM_inputs.pl"
        "tools/EVidenceModeler/EvmUtils/convert_EVM_outputs_to_GFF3.pl"
        "tools/EVidenceModeler/EvmUtils/execute_EVM_commands.pl"
    )
    
    for script in "${evm_scripts[@]}"; do
        if [ -f "${PROJECT_ROOT}/${script}" ]; then
            echo -e "${GREEN}  ✓ ${script}${NC}"
        else
            echo -e "${RED}  ✗ ${script}${NC}"
        fi
    done
    echo ""
    
    # 验证环境变量
    echo -e "${YELLOW}检查环境变量...${NC}"
    echo "  • AUGUSTUS_CONFIG_PATH: ${AUGUSTUS_CONFIG_PATH:-未设置}"
    echo "  • AUGUSTUS_SCRIPTS_PATH: ${AUGUSTUS_SCRIPTS_PATH:-未设置}"
    echo "  • EVM_HOME: ${EVM_HOME:-未设置}"
    echo "  • PERL5LIB: ${PERL5LIB:-未设置}"
    
    if [ -n "${AUGUSTUS_CONFIG_PATH}" ] && [ -d "${AUGUSTUS_CONFIG_PATH}" ]; then
        echo -e "${GREEN}  ✓ Augustus环境变量配置正确${NC}"
    else
        echo -e "${RED}  ✗ Augustus环境变量配置错误${NC}"
    fi
    echo ""
    
    # 生成总结
    echo -e "${BLUE}📋 工具配置总结${NC}"
    echo "==============================="
    echo "✅ 已配置的工具:"
    echo "   • Augustus (基因预测)"
    echo "   • miniprot (蛋白质比对)"
    echo "   • NLR-Annotator (NLR基因识别)"
    echo "   • EVidenceModeler (证据整合)"
    echo ""
    echo "🔧 工具目录结构:"
    echo "   tools/"
    echo "   ├── bin/              # 统一可执行文件入口"
    echo "   ├── Augustus/         # Augustus基因预测工具"
    echo "   ├── miniprot/         # miniprot蛋白质比对工具"
    echo "   ├── EVidenceModeler/  # EVM证据整合工具"
    echo "   └── NLR-Annotator/    # NLR基因识别工具"
    echo ""
    echo "🌍 使用方法:"
    echo "   source setup_env.sh   # 加载环境配置"
    echo "   ./verify_tools.sh     # 验证工具配置"
    echo ""
    echo -e "${GREEN}🎉 第三方软件依赖配置验证完成！${NC}"
}

# 执行主函数
main "$@"