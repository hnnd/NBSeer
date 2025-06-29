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
