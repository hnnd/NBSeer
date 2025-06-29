#!/bin/bash

# NBS基因注释流水线系统依赖安装脚本
# System dependencies installer for NBS annotation pipeline

set -e

# 定义颜色输出
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

echo -e "${BLUE}🔧 NBS基因注释流水线系统依赖安装${NC}"
echo -e "${BLUE}Installing system dependencies for NBS annotation pipeline${NC}"
echo ""

# 检测操作系统
detect_os() {
    if [ -f /etc/os-release ]; then
        . /etc/os-release
        OS=$NAME
        VERSION=$VERSION_ID
    elif type lsb_release >/dev/null 2>&1; then
        OS=$(lsb_release -si)
        VERSION=$(lsb_release -sr)
    elif [ -f /etc/redhat-release ]; then
        OS="Red Hat Enterprise Linux"
        VERSION=$(cat /etc/redhat-release | sed 's/.*release \([0-9]\+\).*/\1/')
    else
        OS=$(uname -s)
        VERSION=$(uname -r)
    fi
    
    echo -e "${YELLOW}检测到操作系统: ${OS} ${VERSION}${NC}"
}

# Ubuntu/Debian系统安装
install_ubuntu_deps() {
    echo -e "${YELLOW}📦 安装Ubuntu/Debian依赖...${NC}"
    
    # 更新包列表
    echo "更新包列表..."
    sudo apt-get update -q
    
    # 安装基本依赖
    echo "安装基本开发工具..."
    sudo apt-get install -y \
        wget \
        curl \
        tar \
        bzip2 \
        build-essential \
        gcc \
        g++ \
        make \
        cmake \
        perl \
        python3 \
        python3-pip \
        default-jre \
        default-jdk \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libncurses5-dev \
        libssl-dev \
        libffi-dev \
        libgsl-dev \
        libsqlite3-dev \
        libmysql++-dev \
        liblpsolve55-dev \
        libbamtools-dev \
        git
    
    echo -e "${GREEN}✓ Ubuntu/Debian依赖安装完成${NC}"
}

# CentOS/RHEL/Fedora系统安装
install_redhat_deps() {
    echo -e "${YELLOW}📦 安装CentOS/RHEL/Fedora依赖...${NC}"
    
    # 检测包管理器
    if command -v dnf &> /dev/null; then
        PKG_MGR="dnf"
    elif command -v yum &> /dev/null; then
        PKG_MGR="yum"
    else
        echo -e "${RED}❌ 未找到包管理器 (yum/dnf)${NC}"
        exit 1
    fi
    
    # 安装基本依赖
    echo "安装基本开发工具..."
    sudo $PKG_MGR groupinstall -y "Development Tools"
    sudo $PKG_MGR install -y \
        wget \
        curl \
        tar \
        bzip2 \
        gcc \
        gcc-c++ \
        make \
        cmake \
        perl \
        python3 \
        python3-pip \
        java-11-openjdk \
        java-11-openjdk-devel \
        zlib-devel \
        bzip2-devel \
        xz-devel \
        ncurses-devel \
        openssl-devel \
        libffi-devel \
        gsl-devel \
        sqlite-devel \
        mysql++-devel \
        lpsolve-devel \
        bamtools-devel \
        git
    
    echo -e "${GREEN}✓ CentOS/RHEL/Fedora依赖安装完成${NC}"
}

# macOS系统安装
install_macos_deps() {
    echo -e "${YELLOW}📦 安装macOS依赖...${NC}"
    
    # 检查是否安装了Homebrew
    if ! command -v brew &> /dev/null; then
        echo "安装Homebrew..."
        /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
    fi
    
    # 安装基本依赖
    echo "安装基本开发工具..."
    brew install \
        wget \
        curl \
        tar \
        bzip2 \
        gcc \
        make \
        cmake \
        perl \
        python3 \
        openjdk \
        zlib \
        bzip2 \
        xz \
        ncurses \
        openssl \
        libffi \
        git
    
    echo -e "${GREEN}✓ macOS依赖安装完成${NC}"
}

# 验证安装
verify_installation() {
    echo -e "${YELLOW}✅ 验证系统依赖安装...${NC}"
    
    local missing=()
    
    # 检查基本命令
    for cmd in wget curl tar make gcc g++ perl python3 java; do
        if ! command -v "$cmd" &> /dev/null; then
            missing+=("$cmd")
        else
            echo -e "${GREEN}  ✓ $cmd${NC}"
        fi
    done
    
    # 检查Python包管理器
    if command -v pip3 &> /dev/null; then
        echo -e "${GREEN}  ✓ pip3${NC}"
    else
        missing+=("pip3")
    fi
    
    # 检查Java版本
    if command -v java &> /dev/null; then
        java_version=$(java -version 2>&1 | head -n1)
        echo -e "${GREEN}  ✓ Java: ${java_version}${NC}"
    fi
    
    # 检查开发库
    if ldconfig -p 2>/dev/null | grep -q "libz.so" || [ -f /usr/lib/libz.dylib ] || [ -f /opt/homebrew/lib/libz.dylib ]; then
        echo -e "${GREEN}  ✓ zlib development library${NC}"
    else
        missing+=("zlib-dev")
    fi
    
    if [ ${#missing[@]} -eq 0 ]; then
        echo -e "${GREEN}🎉 所有系统依赖验证成功！${NC}"
        return 0
    else
        echo -e "${RED}❌ 缺少依赖: ${missing[*]}${NC}"
        return 1
    fi
}

# 安装Python依赖
install_python_deps() {
    echo -e "${YELLOW}🐍 安装Python依赖...${NC}"
    
    # 升级pip
    echo "升级pip..."
    python3 -m pip install --upgrade pip
    
    # 安装常用的生物信息学Python包
    echo "安装Python生物信息学包..."
    pip3 install --user \
        biopython \
        pandas \
        numpy \
        matplotlib \
        seaborn \
        pyyaml \
        requests
    
    echo -e "${GREEN}✓ Python依赖安装完成${NC}"
}

# 主函数
main() {
    detect_os
    
    case "$OS" in
        "Ubuntu"*|"Debian"*)
            install_ubuntu_deps
            ;;
        "CentOS"*|"Red Hat"*|"Fedora"*)
            install_redhat_deps
            ;;
        "Darwin")
            install_macos_deps
            ;;
        *)
            echo -e "${YELLOW}⚠️  未识别的操作系统: $OS${NC}"
            echo -e "${YELLOW}请手动安装以下依赖:${NC}"
            echo "  - wget/curl"
            echo "  - tar, bzip2"
            echo "  - gcc, g++, make, cmake"
            echo "  - perl, python3, java"
            echo "  - zlib, bzip2, xz, ncurses development libraries"
            exit 1
            ;;
    esac
    
    # 安装Python依赖
    install_python_deps
    
    # 验证安装
    if verify_installation; then
        echo ""
        echo -e "${GREEN}🎉 系统依赖安装完成！${NC}"
        echo ""
        echo -e "${BLUE}下一步操作：${NC}"
        echo "1. 运行第三方软件安装: ./setup_tools.sh"
        echo "2. 加载环境配置: source setup_env.sh"
        echo "3. 验证工具安装: ./verify_tools.sh"
    else
        echo ""
        echo -e "${RED}❌ 系统依赖安装不完整，请检查错误信息${NC}"
        exit 1
    fi
}

# 检查是否有root权限
if [ "$EUID" -eq 0 ]; then
    echo -e "${YELLOW}⚠️  以root权限运行${NC}"
else
    echo -e "${YELLOW}ℹ️  需要sudo权限来安装系统包${NC}"
fi

# 执行主函数
main "$@"