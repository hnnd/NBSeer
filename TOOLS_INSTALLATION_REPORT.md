# NBS注释流水线工具安装报告

**安装时间**: Sun 29 Jun 2025 05:11:40 PM CST
**安装目录**: /home/wangys/data/work/nbs/tools

## 已安装软件

### Augustus 3.5.0
- 安装路径: /home/wangys/data/work/nbs/tools/augustus/
- 可执行文件: /home/wangys/data/work/nbs/tools/bin/augustus
- 配置目录: /home/wangys/data/work/nbs/tools/augustus/config/
- 脚本目录: /home/wangys/data/work/nbs/tools/augustus/scripts/

### miniprot 0.17
- 安装路径: /home/wangys/data/work/nbs/tools/miniprot/
- 可执行文件: /home/wangys/data/work/nbs/tools/bin/miniprot

### EVidenceModeler master
- 安装路径: /home/wangys/data/work/nbs/tools/evidencemodeler/
- 主要脚本: EvmUtils/partition_EVM_inputs.pl

### ParaFly master
- 安装路径: /home/wangys/data/work/nbs/tools/evidencemodeler/plugins/ParaFly/
- 可执行文件: ParaFly
- 描述: EVM并行执行依赖

### NLR-Annotator 2.1b
- 安装路径: /home/wangys/data/work/nbs/tools/nlr-annotator/
- JAR文件: NLR-Annotator-v2.1b.jar
- 配置文件: mot.txt, store.txt

## 使用方法

1. 加载环境配置:
   ```bash
   source setup_env.sh
   ```

2. 验证工具安装:
   ```bash
   ./verify_tools.sh
   ```

3. 运行流水线:
   ```bash
   python -m src.nbs_annotation.main --help
   ```

## 目录结构

```
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
```
