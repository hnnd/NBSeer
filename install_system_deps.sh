#!/bin/bash

# NBSeer System Dependencies Installer
# Install system dependencies and Python environment using uv

set -e

# Define color output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

echo -e "${BLUE}🔧 NBSeer System Dependencies Installation${NC}"
echo -e "${BLUE}Installing system dependencies and Python environment with uv${NC}"
echo ""

# Detect operating system
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
    
    echo -e "${YELLOW}Detected OS: ${OS} ${VERSION}${NC}"
}

# Ubuntu/Debian system installation
install_ubuntu_deps() {
    echo -e "${YELLOW}📦 Installing Ubuntu/Debian dependencies...${NC}"
    
    # Update package list
    echo "Updating package list..."
    sudo apt-get update -q
    
    # Install basic dependencies
    echo "Installing basic development tools..."
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
    
    echo -e "${GREEN}✓ Ubuntu/Debian dependencies installed${NC}"
}

# CentOS/RHEL/Fedora system installation
install_redhat_deps() {
    echo -e "${YELLOW}📦 Installing CentOS/RHEL/Fedora dependencies...${NC}"
    
    # Detect package manager
    if command -v dnf &> /dev/null; then
        PKG_MGR="dnf"
    elif command -v yum &> /dev/null; then
        PKG_MGR="yum"
    else
        echo -e "${RED}❌ Package manager not found (yum/dnf)${NC}"
        exit 1
    fi
    
    # Install basic dependencies
    echo "Installing basic development tools..."
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
    
    echo -e "${GREEN}✓ CentOS/RHEL/Fedora dependencies installed${NC}"
}

# macOS system installation
install_macos_deps() {
    echo -e "${YELLOW}📦 Installing macOS dependencies...${NC}"
    
    # Check if Homebrew is installed
    if ! command -v brew &> /dev/null; then
        echo "Installing Homebrew..."
        /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
    fi
    
    # Install basic dependencies
    echo "Installing basic development tools..."
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
    
    echo -e "${GREEN}✓ macOS dependencies installed${NC}"
}

# Verify installation
verify_installation() {
    echo -e "${YELLOW}✅ Verifying system dependencies installation...${NC}"
    
    local missing=()
    
    # Check basic commands
    for cmd in wget curl tar make gcc g++ perl python3 java uv; do
        if ! command -v "$cmd" &> /dev/null; then
            missing+=("$cmd")
        else
            echo -e "${GREEN}  ✓ $cmd${NC}"
        fi
    done
    
    # Check if virtual environment exists
    if [ -d ".venv" ]; then
        echo -e "${GREEN}  ✓ Python virtual environment (.venv)${NC}"
    else
        echo -e "${YELLOW}  ⚠ Virtual environment not found${NC}"
    fi
    
    # Check Java version
    if command -v java &> /dev/null; then
        java_version=$(java -version 2>&1 | head -n1)
        echo -e "${GREEN}  ✓ Java: ${java_version}${NC}"
    fi
    
    # Check development libraries
    if ldconfig -p 2>/dev/null | grep -q "libz.so" || [ -f /usr/lib/libz.dylib ] || [ -f /opt/homebrew/lib/libz.dylib ]; then
        echo -e "${GREEN}  ✓ zlib development library${NC}"
    else
        missing+=("zlib-dev")
    fi
    
    if [ ${#missing[@]} -eq 0 ]; then
        echo -e "${GREEN}🎉 All system dependencies verified successfully!${NC}"
        return 0
    else
        echo -e "${RED}❌ Missing dependencies: ${missing[*]}${NC}"
        return 1
    fi
}

# Install uv and setup Python environment
install_uv_and_python_env() {
    echo -e "${YELLOW}🐍 Installing uv and setting up Python environment...${NC}"
    
    # Check if uv is already installed
    if command -v uv &> /dev/null; then
        echo "uv is already installed"
        uv --version
    else
        echo "Installing uv..."
        curl -LsSf https://astral.sh/uv/install.sh | sh
        export PATH="$HOME/.cargo/bin:$PATH"
    fi
    
    # Check if we're in a project directory with pyproject.toml
    if [ -f "pyproject.toml" ]; then
        echo "Creating virtual environment and installing dependencies with uv..."
        uv venv .venv
        echo "Virtual environment created at .venv"
        
        echo "Installing project dependencies..."
        uv pip install -e .
        
        echo "Installing additional bioinformatics tools..."
        uv pip install \
            gffutils \
            intervaltree \
            psutil
    else
        echo "No pyproject.toml found. Creating basic Python environment..."
        uv venv .venv
        echo "Virtual environment created at .venv"
        
        echo "Installing basic dependencies..."
        uv pip install \
            biopython \
            pandas \
            numpy \
            matplotlib \
            seaborn \
            pyyaml \
            requests \
            gffutils \
            intervaltree \
            psutil
    fi
    
    echo -e "${GREEN}✓ Python environment setup with uv completed${NC}"
}

# Main function
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
            echo -e "${YELLOW}⚠️  Unrecognized operating system: $OS${NC}"
            echo -e "${YELLOW}Please manually install the following dependencies:${NC}"
            echo "  - wget/curl"
            echo "  - tar, bzip2"
            echo "  - gcc, g++, make, cmake"
            echo "  - perl, python3, java"
            echo "  - zlib, bzip2, xz, ncurses development libraries"
            exit 1
            ;;
    esac
    
    # Install uv and Python environment
    install_uv_and_python_env
    
    # Verify installation
    if verify_installation; then
        echo ""
        echo -e "${GREEN}🎉 System dependencies installation completed!${NC}"
        echo ""
        echo -e "${BLUE}Next steps:${NC}"
        echo "1. Activate Python environment: source .venv/bin/activate"
        echo "2. Install bioinformatics tools: ./setup_tools.sh"
        echo "3. Load environment configuration: source setup_env.sh"
        echo "4. Verify tools installation: ./verify_tools.sh"
    else
        echo ""
        echo -e "${RED}❌ System dependencies installation incomplete, please check error messages${NC}"
        exit 1
    fi
}

# Check if running with root privileges
if [ "$EUID" -eq 0 ]; then
    echo -e "${YELLOW}⚠️  Running as root${NC}"
else
    echo -e "${YELLOW}ℹ️  Sudo privileges required for system package installation${NC}"
fi

# Execute main function
main "$@"