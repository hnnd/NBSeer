#!/bin/bash

# NBS基因注释流水线第三方软件自动安装脚本
# Automated setup script for third-party software dependencies

set -e  # Exit on any error

# 定义颜色输出
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# 项目根目录
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TOOLS_DIR="${PROJECT_ROOT}/tools"
TEMP_DIR="${TOOLS_DIR}/temp"  # 临时下载目录

# 软件版本定义
AUGUSTUS_VERSION="3.5.0"
MINIPROT_VERSION="0.17"
EVM_VERSION="master"
PARAFLY_VERSION="master"
NLR_ANNOTATOR_VERSION="2.1b"

# 下载URLs
AUGUSTUS_URL="https://github.com/Gaius-Augustus/Augustus/archive/refs/tags/v${AUGUSTUS_VERSION}.tar.gz"
MINIPROT_URL="https://github.com/lh3/miniprot/releases/download/v${MINIPROT_VERSION}/miniprot-${MINIPROT_VERSION}_x64-linux.tar.bz2"
EVM_URL="https://github.com/EVidenceModeler/EVidenceModeler/archive/refs/heads/master.tar.gz"
PARAFLY_URL="https://github.com/ParaFly/ParaFly/archive/refs/tags/v${PARAFLY_VERSION}.tar.gz"
NLR_ANNOTATOR_URL="https://github.com/steuernb/NLR-Annotator/raw/refs/heads/NLR-Annotator-2/NLR-Annotator-v2.1b.jar"

echo -e "${BLUE}🔧 NBS基因注释流水线第三方软件自动安装${NC}"
echo -e "${BLUE}Automated setup for NBS annotation pipeline dependencies${NC}"
echo -e "${CYAN}软件版本 / Software versions:${NC}"
echo -e "  • Augustus: ${AUGUSTUS_VERSION}"
echo -e "  • miniprot: ${MINIPROT_VERSION}"
echo -e "  • EVidenceModeler: ${EVM_VERSION}"
echo -e "  • ParaFly: ${PARAFLY_VERSION}"
echo -e "  • NLR-Annotator: ${NLR_ANNOTATOR_VERSION}"
echo ""

# 检测系统和依赖
check_system_dependencies() {
    echo -e "${YELLOW}🔍 检查系统依赖...${NC}"
    
    local missing_deps=()
    
    # 检查基本工具
    for cmd in wget curl tar make gcc g++ perl python3 java; do
        if ! command -v "$cmd" &> /dev/null; then
            missing_deps+=("$cmd")
        fi
    done
    
    # 检查开发库
    if ! ldconfig -p | grep -q "libz.so"; then
        missing_deps+=("zlib-dev")
    fi
    
    # 检查GSL库
    if ! ldconfig -p | grep -q "libgsl.so"; then
        missing_deps+=("gsl-dev")
    fi
    
    # 检查SQLite库
    if ! ldconfig -p | grep -q "libsqlite3.so"; then
        missing_deps+=("sqlite3-dev")
    fi
    
    if [ ${#missing_deps[@]} -ne 0 ]; then
        echo -e "${RED}❌ 缺少必要依赖: ${missing_deps[*]}${NC}"
        echo -e "${YELLOW}请先安装缺少的依赖：${NC}"
        echo -e "  Ubuntu/Debian: sudo apt-get install wget curl tar build-essential perl python3 default-jre zlib1g-dev libgsl-dev libsqlite3-dev libmysql++-dev liblpsolve55-dev libbamtools-dev"
        echo -e "  CentOS/RHEL: sudo yum install wget curl tar gcc gcc-c++ make perl python3 java zlib-devel gsl-devel sqlite-devel mysql++-devel lpsolve-devel bamtools-devel"
        exit 1
    fi
    
    echo -e "${GREEN}✓ 系统依赖检查完成${NC}"
}

# 创建tools目录结构
create_tools_structure() {
    echo -e "${YELLOW}📁 创建tools目录结构...${NC}"
    
    mkdir -p "${TOOLS_DIR}/bin"
    mkdir -p "${TOOLS_DIR}/lib"
    mkdir -p "${TEMP_DIR}"
    mkdir -p "${TOOLS_DIR}/augustus"
    mkdir -p "${TOOLS_DIR}/miniprot"
    mkdir -p "${TOOLS_DIR}/evidencemodeler"
    mkdir -p "${TOOLS_DIR}/nlr-annotator"
    
    echo -e "${GREEN}✓ Tools目录结构创建完成${NC}"
}

# 下载文件函数
download_file() {
    local url="$1"
    local output_file="$2"
    local description="$3"
    
    echo -e "  📥 下载 ${description}..."
    
    if [ -f "${output_file}" ]; then
        echo -e "    ℹ️  文件已存在，跳过下载"
        return 0
    fi
    
    # 尝试使用wget，如果失败则使用curl
    if command -v wget &> /dev/null; then
        wget -q --show-progress -O "${output_file}" "${url}"
    elif command -v curl &> /dev/null; then
        curl -L -o "${output_file}" "${url}" --progress-bar
    else
        echo -e "${RED}❌ 无法找到wget或curl，无法下载文件${NC}"
        exit 1
    fi
    
    if [ $? -eq 0 ]; then
        echo -e "    ✓ 下载完成"
    else
        echo -e "${RED}❌ 下载失败: ${description}${NC}"
        exit 1
    fi
}

# 安装Augustus
install_augustus() {
    echo -e "${YELLOW}🧬 安装Augustus...${NC}"
    
    local augustus_archive="${TEMP_DIR}/augustus-${AUGUSTUS_VERSION}.tar.gz"
    local augustus_dir="${TOOLS_DIR}/augustus"
    
    # 检查是否已安装
    if [ -f "${augustus_dir}/bin/augustus" ]; then
        echo -e "${GREEN}  ✓ Augustus已安装${NC}"
        return 0
    fi
    
    # 下载
    download_file "${AUGUSTUS_URL}" "${augustus_archive}" "Augustus ${AUGUSTUS_VERSION}"
    
    # 在临时目录解压和编译
    echo -e "  📦 解压Augustus..."
    cd "${TEMP_DIR}"
    tar -xzf "augustus-${AUGUSTUS_VERSION}.tar.gz"
    
    # 编译
    echo -e "  🔨 编译Augustus..."
    cd "Augustus-${AUGUSTUS_VERSION}"
    
    # 设置编译选项
    export AUGUSTUS_CONFIG_PATH="${augustus_dir}/config"
    
    # 禁用可选依赖的编译选项
    echo "# Augustus simplified configuration" > common.mk.local
    echo "COMPGENEPRED = false" >> common.mk.local
    echo "SQLITE = false" >> common.mk.local
    echo "MYSQL = false" >> common.mk.local
    echo "BAMTOOLS = false" >> common.mk.local
    echo "GSL = false" >> common.mk.local
    echo "LPSOLVE = false" >> common.mk.local
    echo "ZIPINPUT = false" >> common.mk.local
    
    # 编译核心Augustus
    make clean &> /dev/null || true
    make augustus COMPGENEPRED=false SQLITE=false MYSQL=false BAMTOOLS=false GSL=false LPSOLVE=false ZIPINPUT=false -j$(nproc) &> augustus_build.log
    
    if [ $? -ne 0 ]; then
        echo -e "${RED}❌ Augustus编译失败${NC}"
        echo -e "查看编译日志: ${TEMP_DIR}/Augustus-${AUGUSTUS_VERSION}/augustus_build.log"
        exit 1
    fi
    
    # 直接安装到augustus目录
    echo -e "  📋 安装Augustus到 ${augustus_dir}..."
    cp -r bin "${augustus_dir}/"
    cp -r config "${augustus_dir}/"
    cp -r scripts "${augustus_dir}/"
    [ -d docs ] && cp -r docs "${augustus_dir}/"
    
    # 创建可执行文件链接
    ln -sf "../augustus/bin/augustus" "${TOOLS_DIR}/bin/augustus"
    ln -sf "../augustus/bin/etraining" "${TOOLS_DIR}/bin/etraining""
    
    # 设置权限
    chmod +x "${augustus_dir}/bin/"*
    chmod +x "${augustus_dir}/scripts/"*.pl 2>/dev/null || true
    
    echo -e "${GREEN}  ✓ Augustus安装完成${NC}"
    cd "${PROJECT_ROOT}"
}

# 安装miniprot
install_miniprot() {
    echo -e "${YELLOW}🧬 安装miniprot...${NC}"
    
    local miniprot_archive="${TEMP_DIR}/miniprot-${MINIPROT_VERSION}.tar.bz2"
    local miniprot_dir="${TOOLS_DIR}/miniprot"
    
    # 检查是否已安装
    if [ -f "${miniprot_dir}/bin/miniprot" ]; then
        echo -e "${GREEN}  ✓ miniprot已安装${NC}"
        return 0
    fi
    
    # 下载
    download_file "${MINIPROT_URL}" "${miniprot_archive}" "miniprot ${MINIPROT_VERSION}"
    
    # 解压到临时目录
    echo -e "  📦 解压miniprot..."
    cd "${TEMP_DIR}"
    tar -xjf "miniprot-${MINIPROT_VERSION}.tar.bz2"
    
    # 创建miniprot目录结构并安装
    echo -e "  📋 安装miniprot到 ${miniprot_dir}..."
    mkdir -p "${miniprot_dir}/bin"
    mkdir -p "${miniprot_dir}/doc"
    
    # 复制可执行文件
    cp "miniprot-${MINIPROT_VERSION}_x64-linux/miniprot" "${miniprot_dir}/bin/"
    
    # 复制文档（如果存在）
    if [ -f "miniprot-${MINIPROT_VERSION}_x64-linux/README.md" ]; then
        cp "miniprot-${MINIPROT_VERSION}_x64-linux/README.md" "${miniprot_dir}/doc/"
    fi
    
    # 创建可执行文件链接
    ln -sf "../miniprot/bin/miniprot" "${TOOLS_DIR}/bin/miniprot"
    
    # 设置权限
    chmod +x "${miniprot_dir}/bin/miniprot"
    
    echo -e "${GREEN}  ✓ miniprot安装完成${NC}"
    cd "${PROJECT_ROOT}"
}

# 安装ParaFly (EVM依赖)
install_parafly() {
    local evm_dir="$1"
    local parafly_dir="${evm_dir}/plugins/ParaFly"
    
    echo -e "  📦 安装ParaFly..."
    
    # 检查是否已安装
    if [ -f "${parafly_dir}/bin/ParaFly" ]; then
        echo -e "    ✓ ParaFly已安装"
        return 0
    fi
    
    # 下载ParaFly
    local parafly_archive="${TEMP_DIR}/parafly-${PARAFLY_VERSION}.tar.gz"
    download_file "${PARAFLY_URL}" "${parafly_archive}" "ParaFly ${PARAFLY_VERSION}"
    
    # 解压并编译
    cd "${TEMP_DIR}"
    tar -xzf "parafly-${PARAFLY_VERSION}.tar.gz"
    if [ "${PARAFLY_VERSION}" = "master" ]; then
        cd "ParaFly-master"
    else
        cd "ParaFly-${PARAFLY_VERSION}"
    fi
    
    # 检查是否已预编译
    if [ -f "bin/ParaFly" ]; then
        echo -e "    ✓ 使用预编译的ParaFly"
        # 安装到EVM的plugins目录
        mkdir -p "${parafly_dir}/bin"
        cp bin/ParaFly "${parafly_dir}/bin/"
        chmod +x "${parafly_dir}/bin/ParaFly"
    else
        # 编译ParaFly
        echo -e "    🔨 编译ParaFly..."
        if make -j$(nproc) > parafly_build.log 2>&1; then
            echo -e "    ✓ ParaFly编译成功"
            # 安装到EVM的plugins目录
            mkdir -p "${parafly_dir}/bin"
            cp ParaFly "${parafly_dir}/bin/"
            chmod +x "${parafly_dir}/bin/ParaFly"
        else
            echo -e "    ⚠️  ParaFly编译失败，尝试简单编译..."
            # 尝试简单的编译方法
            if g++ -o ParaFly src/ParaFly.cpp -std=c++11 -pthread 2>/dev/null; then
                mkdir -p "${parafly_dir}/bin"
                cp ParaFly "${parafly_dir}/bin/"
                chmod +x "${parafly_dir}/bin/ParaFly"
            else
                echo -e "    ❌ ParaFly编译失败，跳过安装"
                cd "${PROJECT_ROOT}"
                return 1
            fi
        fi
    fi
    
    # 创建符号链接到EVM目录
    ln -sf "plugins/ParaFly/bin/ParaFly" "${evm_dir}/ParaFly" 2>/dev/null || true
    
    echo -e "    ✓ ParaFly安装完成"
    cd "${PROJECT_ROOT}"
}

# 安装EVidenceModeler
install_evm() {
    echo -e "${YELLOW}🧬 安装EVidenceModeler...${NC}"
    
    local evm_archive="${TEMP_DIR}/evm-${EVM_VERSION}.tar.gz"
    local evm_dir="${TOOLS_DIR}/evidencemodeler"
    
    # 检查是否已安装
    if [ -f "${evm_dir}/EvmUtils/partition_EVM_inputs.pl" ]; then
        echo -e "${GREEN}  ✓ EVidenceModeler已安装${NC}"
        return 0
    fi
    
    # 下载
    download_file "${EVM_URL}" "${evm_archive}" "EVidenceModeler ${EVM_VERSION}"
    
    # 解压到临时目录
    echo -e "  📦 解压EVidenceModeler..."
    cd "${TEMP_DIR}"
    tar -xzf "evm-${EVM_VERSION}.tar.gz"
    
    # 直接安装到evidencemodeler目录
    echo -e "  📋 安装EVidenceModeler到 ${evm_dir}..."
    if [ "${EVM_VERSION}" = "master" ]; then
        cp -r "EVidenceModeler-master/"* "${evm_dir}/"
    else
        cp -r "EVidenceModeler-${EVM_VERSION}/"* "${evm_dir}/"
    fi
    
    # 设置权限
    find "${evm_dir}" -name "*.pl" -exec chmod +x {} \;
    chmod +x "${evm_dir}/EVidenceModeler" 2>/dev/null || true
    
    # 安装ParaFly
    install_parafly "${evm_dir}"
    
    # 创建版本信息文件
    echo "EVidenceModeler ${EVM_VERSION}" > "${evm_dir}/VERSION"
    echo "ParaFly ${PARAFLY_VERSION}" >> "${evm_dir}/VERSION"
    echo "Installed on: $(date)" >> "${evm_dir}/VERSION"
    
    echo -e "${GREEN}  ✓ EVidenceModeler安装完成${NC}"
    cd "${PROJECT_ROOT}"
}

# 安装NLR-Annotator
install_nlr_annotator() {
    echo -e "${YELLOW}🧬 安装NLR-Annotator...${NC}"
    
    local nlr_dir="${TOOLS_DIR}/nlr-annotator"
    local nlr_jar="${nlr_dir}/NLR-Annotator-v${NLR_ANNOTATOR_VERSION}.jar"
    
    # 检查是否已安装
    if [ -f "${nlr_jar}" ]; then
        echo -e "${GREEN}  ✓ NLR-Annotator已安装${NC}"
        return 0
    fi
    
    # 创建目录结构
    mkdir -p "${nlr_dir}/bin"
    mkdir -p "${nlr_dir}/config"
    mkdir -p "${nlr_dir}/doc"
    
    # 下载主程序
    echo -e "  📋 安装NLR-Annotator到 ${nlr_dir}..."
    download_file "${NLR_ANNOTATOR_URL}" "${nlr_jar}" "NLR-Annotator ${NLR_ANNOTATOR_VERSION}"
    
    # 下载配置文件
    local mot_url="https://raw.githubusercontent.com/steuernb/NLR-Annotator/NLR-Annotator-2/src/mot.txt"
    local store_url="https://raw.githubusercontent.com/steuernb/NLR-Annotator/NLR-Annotator-2/src/store.txt"
    
    download_file "${mot_url}" "${nlr_dir}/config/mot.txt" "NLR-Annotator mot.txt"
    download_file "${store_url}" "${nlr_dir}/config/store.txt" "NLR-Annotator store.txt"
    
    # 创建启动脚本
    cat > "${nlr_dir}/bin/nlr-annotator" << EOF
#!/bin/bash
# NLR-Annotator wrapper script
NLR_DIR="\$(cd "\$(dirname "\${BASH_SOURCE[0]}")" && pwd)/.."
java -jar "\${NLR_DIR}/NLR-Annotator-v${NLR_ANNOTATOR_VERSION}.jar" "\$@"
EOF
    
    chmod +x "${nlr_dir}/bin/nlr-annotator"
    
    # 创建版本信息
    echo "NLR-Annotator v${NLR_ANNOTATOR_VERSION}" > "${nlr_dir}/VERSION"
    echo "Installed on: $(date)" >> "${nlr_dir}/VERSION"
    
    # 创建README
    cat > "${nlr_dir}/README.md" << EOF
# NLR-Annotator v${NLR_ANNOTATOR_VERSION}

## Usage
\`\`\`bash
# Using wrapper script
\${TOOLS_DIR}/nlr-annotator/bin/nlr-annotator [options]

# Using java directly
java -jar \${TOOLS_DIR}/nlr-annotator/NLR-Annotator-v${NLR_ANNOTATOR_VERSION}.jar [options]
\`\`\`

## Files
- \`NLR-Annotator-v${NLR_ANNOTATOR_VERSION}.jar\`: Main program
- \`config/mot.txt\`: Motif configuration file
- \`config/store.txt\`: Store configuration file
- \`bin/nlr-annotator\`: Wrapper script
EOF
    
    echo -e "${GREEN}  ✓ NLR-Annotator安装完成${NC}"
}

# 验证安装
verify_installations() {
    echo -e "${YELLOW}✅ 验证安装...${NC}"
    
    local all_good=true
    
    # 验证Augustus
    if [ -f "${TOOLS_DIR}/bin/augustus" ] && [ -x "${TOOLS_DIR}/bin/augustus" ]; then
        local augustus_version=$("${TOOLS_DIR}/bin/augustus" --version 2>&1 | head -n1 || echo "unknown")
        echo -e "${GREEN}  ✓ Augustus: ${augustus_version}${NC}"
    else
        echo -e "${RED}  ❌ Augustus安装失败${NC}"
        all_good=false
    fi
    
    # 验证miniprot
    if [ -f "${TOOLS_DIR}/bin/miniprot" ] && [ -x "${TOOLS_DIR}/bin/miniprot" ]; then
        local miniprot_version=$("${TOOLS_DIR}/bin/miniprot" --version 2>&1 || echo "unknown")
        echo -e "${GREEN}  ✓ miniprot: ${miniprot_version}${NC}"
    else
        echo -e "${RED}  ❌ miniprot安装失败${NC}"
        all_good=false
    fi
    
    # 验证EVidenceModeler
    if [ -f "${TOOLS_DIR}/evidencemodeler/EvmUtils/partition_EVM_inputs.pl" ]; then
        echo -e "${GREEN}  ✓ EVidenceModeler: ${EVM_VERSION}${NC}"
        
        # 验证ParaFly
        if [ -f "${TOOLS_DIR}/evidencemodeler/plugins/ParaFly/bin/ParaFly" ]; then
            echo -e "${GREEN}  ✓ ParaFly: ${PARAFLY_VERSION}${NC}"
        else
            echo -e "${YELLOW}  ⚠️  ParaFly未安装（可选组件）${NC}"
        fi
    else
        echo -e "${RED}  ❌ EVidenceModeler安装失败${NC}"
        all_good=false
    fi
    
    # 验证NLR-Annotator
    if [ -f "${TOOLS_DIR}/nlr-annotator/NLR-Annotator-v${NLR_ANNOTATOR_VERSION}.jar" ]; then
        echo -e "${GREEN}  ✓ NLR-Annotator: v${NLR_ANNOTATOR_VERSION}${NC}"
    else
        echo -e "${RED}  ❌ NLR-Annotator安装失败${NC}"
        all_good=false
    fi
    
    if [ "$all_good" = true ]; then
        echo -e "${GREEN}🎉 所有软件安装验证成功！${NC}"
    else
        echo -e "${RED}❌ 部分软件安装失败${NC}"
        exit 1
    fi
}

# 创建环境配置脚本
create_env_script() {
    echo -e "${YELLOW}🌍 创建环境配置脚本...${NC}"
    
    cat > "${PROJECT_ROOT}/setup_env.sh" << 'EOF'
#!/bin/bash
# NBS基因注释流水线环境配置脚本

# 获取项目根目录
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TOOLS_DIR="${PROJECT_ROOT}/tools"

# 添加tools/bin到PATH
export PATH="${TOOLS_DIR}/bin:$PATH"

# Augustus相关环境变量
export AUGUSTUS_CONFIG_PATH="${TOOLS_DIR}/augustus/config"
export AUGUSTUS_SCRIPTS_PATH="${TOOLS_DIR}/augustus/scripts"

# EVM相关环境变量
export PERL5LIB="${TOOLS_DIR}/evidencemodeler:$PERL5LIB"
export EVM_HOME="${TOOLS_DIR}/evidencemodeler"

echo "✅ NBS基因注释流水线环境已配置"
echo "🔧 工具目录: ${TOOLS_DIR}"
echo "📋 可用工具:"
echo "   • Augustus: $(which augustus 2>/dev/null || echo '未找到')"
echo "   • miniprot: $(which miniprot 2>/dev/null || echo '未找到')"
echo "   • NLR-Annotator: ${TOOLS_DIR}/nlr-annotator/NLR-Annotator-v2.1b.jar"
echo "   • EVidenceModeler: ${TOOLS_DIR}/evidencemodeler/"
EOF

    chmod +x "${PROJECT_ROOT}/setup_env.sh"
    echo -e "${GREEN}  ✓ 环境配置脚本已创建: setup_env.sh${NC}"
}

# 生成工具配置报告
generate_setup_report() {
    echo -e "${YELLOW}📋 生成安装报告...${NC}"
    
    cat > "${PROJECT_ROOT}/TOOLS_INSTALLATION_REPORT.md" << EOF
# NBS注释流水线工具安装报告

**安装时间**: $(date)
**安装目录**: ${TOOLS_DIR}

## 已安装软件

### Augustus ${AUGUSTUS_VERSION}
- 安装路径: ${TOOLS_DIR}/augustus/
- 可执行文件: ${TOOLS_DIR}/bin/augustus
- 配置目录: ${TOOLS_DIR}/augustus/config/
- 脚本目录: ${TOOLS_DIR}/augustus/scripts/

### miniprot ${MINIPROT_VERSION}
- 安装路径: ${TOOLS_DIR}/miniprot/
- 可执行文件: ${TOOLS_DIR}/bin/miniprot

### EVidenceModeler ${EVM_VERSION}
- 安装路径: ${TOOLS_DIR}/evidencemodeler/
- 主要脚本: EvmUtils/partition_EVM_inputs.pl

### ParaFly ${PARAFLY_VERSION}
- 安装路径: ${TOOLS_DIR}/evidencemodeler/plugins/ParaFly/
- 可执行文件: ParaFly
- 描述: EVM并行执行依赖

### NLR-Annotator ${NLR_ANNOTATOR_VERSION}
- 安装路径: ${TOOLS_DIR}/nlr-annotator/
- JAR文件: NLR-Annotator-v${NLR_ANNOTATOR_VERSION}.jar
- 配置文件: mot.txt, store.txt

## 使用方法

1. 加载环境配置:
   \`\`\`bash
   source setup_env.sh
   \`\`\`

2. 验证工具安装:
   \`\`\`bash
   ./verify_tools.sh
   \`\`\`

3. 运行流水线:
   \`\`\`bash
   python -m src.nbs_annotation.main --help
   \`\`\`

## 目录结构

\`\`\`
tools/
├── bin/                    # 统一可执行文件入口
│   ├── augustus           # Augustus可执行文件链接
│   └── miniprot           # miniprot可执行文件链接
├── augustus/              # Augustus基因预测工具
│   ├── bin/
│   ├── config/
│   └── scripts/
├── miniprot/              # miniprot蛋白质比对工具
│   └── bin/
├── evidencemodeler/       # EVM证据整合工具
│   ├── EvmUtils/
│   ├── plugins/
│   │   └── ParaFly/       # ParaFly并行执行工具
│   └── PerlLib/
├── nlr-annotator/         # NLR基因识别工具
│   ├── bin/
│   └── config/
└── temp/                  # 临时下载文件目录
\`\`\`
EOF

    echo -e "${GREEN}  ✓ 安装报告已生成: TOOLS_INSTALLATION_REPORT.md${NC}"
}

# 清理函数
cleanup_downloads() {
    echo -e "${YELLOW}🧹 清理下载文件...${NC}"
    
    read -p "是否删除下载的源文件以节省空间？[y/N] " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        rm -rf "${TEMP_DIR}"
        echo -e "${GREEN}  ✓ 临时文件已清理${NC}"
    else
        echo -e "  ℹ️  保留临时文件在: ${TEMP_DIR}"
    fi
}

# 主执行函数
main() {
    echo "开始自动安装..."
    
    check_system_dependencies
    create_tools_structure
    
    install_augustus
    install_miniprot
    install_evm
    install_nlr_annotator
    
    verify_installations
    create_env_script
    generate_setup_report
    
    echo ""
    echo -e "${GREEN}🎉 NBS基因注释流水线第三方软件安装完成！${NC}"
    echo ""
    echo -e "${BLUE}下一步操作：${NC}"
    echo "1. 加载环境配置: source setup_env.sh"
    echo "2. 验证工具安装: ./verify_tools.sh"
    echo "3. 查看安装报告: cat TOOLS_INSTALLATION_REPORT.md"
    echo ""
    
    cleanup_downloads
}

# 检查是否以root权限运行
if [ "$EUID" -eq 0 ]; then
    echo -e "${YELLOW}⚠️  警告：正在以root权限运行，建议使用普通用户权限${NC}"
    read -p "是否继续？[y/N] " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        exit 1
    fi
fi

# 执行主函数
main "$@"
