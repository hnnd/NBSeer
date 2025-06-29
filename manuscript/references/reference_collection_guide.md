# 参考文献收集和管理指南

## 文献检索策略

### 核心关键词组合
1. **植物免疫基础**
   - "plant immunity" + "NBS genes"
   - "plant disease resistance" + "NLR genes"  
   - "innate immunity" + "plant pathogens"

2. **基因注释方法**
   - "gene annotation" + "bioinformatics"
   - "genome annotation" + "pipeline"
   - "ab initio gene prediction"
   - "evidence-based gene annotation"

3. **NBS基因特异性研究**
   - "NBS gene annotation" 
   - "NLR gene prediction"
   - "resistance gene analogs"
   - "plant R gene identification"

4. **计算工具和算法**
   - "Augustus gene prediction"
   - "EVidenceModeler"
   - "protein alignment genome"
   - "miniprot spliced alignment"

### 主要数据库和期刊

#### 高影响因子期刊 (优先检索)
- **Nature** (IF: 64.8)
- **Science** (IF: 56.9) 
- **Nature Plants** (IF: 15.8)
- **The Plant Cell** (IF: 12.0)
- **Plant Physiology** (IF: 7.4)
- **Bioinformatics** (IF: 5.8)
- **Genome Research** (IF: 6.2)

#### 专业期刊
- **BMC Bioinformatics** (IF: 3.3)
- **BMC Genomics** (IF: 4.0)
- **Plant Methods** (IF: 5.1)
- **Frontiers in Plant Science** (IF: 5.6)
- **New Phytologist** (IF: 10.3)

#### 数据库检索
- **PubMed/MEDLINE**: 生物医学文献主库
- **Web of Science**: 综合学术数据库
- **Google Scholar**: 补充检索
- **bioRxiv**: 预印本文献
- **arXiv**: 计算机科学预印本

---

## 文献分类体系

### 1. 背景和综述类 (Introduction用)
- [ ] 植物免疫系统综述
- [ ] NBS/NLR基因家族进化
- [ ] 植物-病原菌互作机制
- [ ] 作物抗病基因挖掘现状

**目标数量**: 15-20篇

### 2. 方法学和工具类 (Methods用)
- [ ] 基因注释算法和流水线
- [ ] Augustus基因预测方法
- [ ] EVidenceModeler证据整合
- [ ] 蛋白质比对算法
- [ ] 序列分析工具

**目标数量**: 20-25篇

### 3. NBS基因专门研究 (全文用)
- [ ] NBS基因结构特征
- [ ] NBS基因聚类和分布
- [ ] NBS基因功能验证
- [ ] NBS基因注释工具比较

**目标数量**: 15-20篇

### 4. 基因组学参考数据 (Methods和Results用)
- [ ] 模式植物基因组注释
- [ ] 作物基因组项目
- [ ] 基因组质量评估标准
- [ ] 参考数据库和资源

**目标数量**: 10-15篇

### 5. 评估和验证方法 (Results用)
- [ ] 基因注释质量评估
- [ ] 生物信息学工具基准测试
- [ ] 统计分析方法
- [ ] 性能评估指标

**目标数量**: 10-12篇

---

## 文献质量筛选标准

### 必须包含的文献
- **工具开发原始论文**: Augustus, EVidenceModeler, miniprot等
- **NBS基因经典综述**: 近5年内高引用综述
- **基准数据集论文**: TAIR, RAP-DB等参考基因组
- **方法学重要论文**: 基因注释评估标准

### 优先选择标准
1. **时效性**: 优先选择2018年后发表的文献
2. **影响因子**: 优先选择IF>3.0的期刊文献
3. **引用次数**: 优先选择高引用文献
4. **相关性**: 与本研究方法和结果直接相关

### 排除标准
- 发表时间超过10年且非经典论文
- 影响因子过低 (<1.5) 且非专业领域期刊
- 与研究主题相关性较低
- 质量不高的会议论文和技术报告

---

## 文献管理工具推荐

### 1. Zotero (推荐)
- **优点**: 免费、浏览器集成、协作功能
- **插件**: Better BibTeX for Zotero
- **同步**: 云端同步，多设备访问

### 2. Mendeley
- **优点**: PDF管理强大、社交功能
- **缺点**: 存储空间限制

### 3. EndNote
- **优点**: 功能全面、机构支持
- **缺点**: 需要付费许可

### 推荐工作流
```
文献检索 → Zotero收集 → PDF获取 → 笔记整理 → BibTeX导出 → LaTeX集成
```

---

## 引用格式和规范

### Bioinformatics期刊引用格式
```
期刊文章:
Author,A.B. and Author,C.D. (year) Title of the article. Journal Name, vol, pages.

书籍:
Author,A.B. (year) Book Title. Publisher, Place.

网络资源:
Author,A.B. (year) Title. URL (访问日期).
```

### 引用管理检查清单
- [ ] 所有引用在文本中有对应标记
- [ ] 参考文献格式统一
- [ ] DOI和URL链接有效
- [ ] 作者姓名拼写正确
- [ ] 期刊名称缩写标准

---

## 文献阅读和笔记策略

### 快速筛选方法
1. **标题+摘要**: 判断相关性
2. **图表浏览**: 了解主要结果
3. **结论部分**: 把握核心贡献
4. **方法概述**: 评估技术细节

### 深度阅读重点
- **Introduction**: 背景知识和研究现状
- **Methods**: 技术细节和参数设置
- **Results**: 关键数据和性能指标
- **Discussion**: 局限性和改进建议

### 笔记记录模板
```markdown
## 文献信息
- 标题: [标题]
- 作者: [作者]
- 期刊: [期刊名] (年份)
- DOI: [DOI链接]

## 主要贡献
- 

## 方法特点
- 

## 关键结果
- 

## 相关性评估
- 与本研究关系: [高/中/低]
- 引用价值: [引言/方法/结果/讨论]

## 引用要点
- 
```

---

## 时间安排和里程碑

### 第1周: 基础文献收集
- [ ] 检索和下载核心工具论文 (20篇)
- [ ] 收集NBS基因综述文献 (15篇)
- [ ] 整理基因组参考数据论文 (10篇)

### 第2周: 方法学文献补充
- [ ] 深入检索基因注释方法论文 (25篇)
- [ ] 收集评估和验证相关文献 (12篇)
- [ ] 补充最新技术发展文献 (10篇)

### 第3周: 文献整理和分类
- [ ] 建立Zotero文献库
- [ ] 按主题分类整理
- [ ] 制作文献阅读优先级列表

### 第4周: 精读和笔记
- [ ] 精读核心文献 (30篇)
- [ ] 整理引用要点
- [ ] 完善参考文献列表

### 质量控制检查点
- 总文献数量: 80-100篇
- 高质量文献比例: >60%
- 近5年文献比例: >70%
- 直接相关文献比例: >80%