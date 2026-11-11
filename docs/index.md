# DNAKit

**DNAKit: A Comprehensive Toolkit for Efficient DNA Research**

<p align="center">
  <a href="https://pypi.org/project/dnakit/">
    <img src="https://img.shields.io/pypi/v/dnakit?include_prereleases=true&label=PyPI&logo=pypi&cacheSeconds=300" alt="PyPI version">
  </a>
  <a href="https://github.com/mapengsen/DNAKit/tree/main/packaging/bioconda">
    <img src="https://img.shields.io/badge/Bioconda-recipe-43B02A?logo=anaconda&logoColor=white" alt="Bioconda recipe">
  </a>
  <a href="https://github.com/mapengsen/DNAKit/tree/main/galaxy/dnakit">
    <img src="https://img.shields.io/badge/Galaxy-wrapper-2C3143?logo=galaxy&logoColor=white" alt="Galaxy wrapper">
  </a>
  <a href="https://github.com/mapengsen/DNAKit/tree/main/packaging/guix">
    <img src="https://img.shields.io/badge/GNU%20Guix-package-A42E2B?logo=gnu&logoColor=white" alt="GNU Guix package">
  </a>
  <a href="llm_tools/mcp/">
    <img src="https://img.shields.io/badge/MCP-server-1A1A1A?logo=modelcontextprotocol&logoColor=white" alt="MCP server">
  </a>
  <a href="llm_tools/skills/">
    <img src="https://img.shields.io/badge/Skills-agent%20workflows-D97757?logo=data:image/svg%2bxml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCI+PHBhdGggZmlsbD0id2hpdGUiIGQ9Ik0xMyAxIDMgMTRoN2wtMS41IDlMMjAgOS41aC03LjJ6Ii8+PC9zdmc+" alt="Agent Skills">
  </a>
  <a href="https://dnakit.mapengsen.asia">
    <img src="https://img.shields.io/badge/Web%20App-online-2EA44F?logo=googlechrome&logoColor=white" alt="Web app">
  </a>
  <a href="gui/">
    <img src="https://img.shields.io/badge/GUI-Windows-0078D6?logo=data:image/svg%2bxml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCI+PHBhdGggZmlsbD0id2hpdGUiIGQ9Ik0wIDBoMTEuNHYxMS40SDB6TTEyLjYgMEgyNHYxMS40SDEyLjZ6TTAgMTIuNmgxMS40VjI0SDB6TTEyLjYgMTIuNkgyNFYyNEgxMi42eiIvPjwvc3ZnPg==" alt="GUI for Windows">
  </a>
  <a href="gui/">
    <img src="https://img.shields.io/badge/GUI-macOS-000000?logo=apple&logoColor=white" alt="GUI for macOS">
  </a>
  <a href="gui/">
    <img src="https://img.shields.io/badge/GUI-Linux-FCC624?logo=linux&logoColor=black" alt="GUI for Linux">
  </a>
</p>

DNAKit 是面向 DNA 分析的统一工具包。它按同一套规范整合了 284 个现有软件、Web 服务和 Agent 工具集的功能，提供 6,698 个操作；研究者和 AI Agent 用同一套函数，就能完成从数据下载到模型训练的整个流程。

<div class="dnakit-stats">
  <div><strong>6,698</strong><span>操作</span></div>
  <div><strong>12</strong><span>DNA 分析领域</span></div>
  <div><strong>148</strong><span>文件格式</span></div>
  <div><strong>107</strong><span>公共数据库</span></div>
  <div><strong>48</strong><span>DNA 基础模型</span></div>
  <div><strong>7</strong><span>使用方式</span></div>
</div>

# 为什么用 DNAKit

- **一套规范，不用逐个适配工具：** 所有操作使用相同的命名、数据格式、返回结果和错误提示，不同软件之间不再需要自己写衔接代码。
- **一次安装，七种用法：** [Python API](llm_tools/python_api.md)、命令行、[网页版](https://dnakit.mapengsen.asia)、[桌面版](gui.md)、[MCP](llm_tools/mcp.md)、[Skills](llm_tools/skills.md) 和内置 [Agent](llm_tools/agent.md) 共用同一套实现，结果一致。
- **为 AI Agent 准备好的 DNA 工具：** 6,120 个操作可通过 MCP 调用，另有 575 个工作流 Skills；同一个请求在 Codex、Claude Code 和 Gemini CLI 中得到相同结果。
- **不写代码也能用：** 网页版打开浏览器即可使用，无需安装；桌面版（`dnakit gui`）支持 Windows、macOS 和 Linux。
- **每一步都可追溯：** 每次操作对应一次明确的函数调用，输入、参数、后端和参考数据库都会记录下来，便于检查和复现。
- **可以自己扩展：** 写一个带类型注解和文档的 Python 函数即可注册为工具；也可以用自然语言描述需求，由 Agent 生成工具，通过检查后注册。

# 能做什么

| 领域 | 操作数 | 例子 |
|:---|---:|:---|
| 数据与数据库 | 692 | FASTA、BAM、VCF、BED 等 148 种格式的读写与转换；NCBI、Ensembl、ENCODE、ClinVar 等 107 个数据库的查询与下载 |
| 序列分析 | 251 | 反向互补、转录与翻译、ORF 查找、k-mer 与序列组成、DNA 描述符 |
| 比对与同源 | 328 | 序列比对、相似度与距离、BLAST 同源搜索 |
| 测序读段与组装 | 566 | 读段修剪与质控、UMI 处理、读段比对、基因组组装 |
| 变异与临床 | 907 | 变异检测、VEP 注释、ClinVar 临床意义、结构变异 |
| 群体遗传与 GWAS | 680 | 群体遗传统计、单倍型与连锁不平衡、全基因组关联分析 |
| 基因与表达 | 1,009 | 基因注释与同源基因、表达与富集分析、单细胞与空间转录组 |
| 调控与表观基因组 | 666 | motif 扫描、峰分析、甲基化、Hi-C 接触矩阵 |
| 进化与微生物组 | 452 | 系统发育树、dN/dS、泛基因组、宏基因组与物种分类 |
| 设计与工程 | 389 | 引物设计、熔解温度与热力学、限制性酶切、Golden Gate 组装、CRISPR 设计 |
| 序列机器学习 | 195 | 数据集去重、聚类与划分，指纹与特征，性质预测模型训练 |
| DNA 基础模型 | 233 | AlphaGenome、Evo 2 等 48 个模型的表征提取、微调与预测 |

DNAKit 还提供生成序列的评价（有效性、唯一性、多样性、新颖性）和可视化（序列文字图、位置高亮、线性与环状 DNA 图、序列比对图）。另有 330 个通用操作（如基础设施和文献检索）不归入上述领域。完整清单见[功能树总览](api/features/function_tree.md)。

# 快速开始

```bash
pip install dnakit
```

```python
from dnakit.ops import reverse_complement

print(reverse_complement("AACG"))  # CGTT
```

```bash
pip install "dnakit[agent]"
# 注册为 MCP 服务器
claude mcp add dnakit -- dnakit-mcp     # Claude Code
codex mcp add dnakit -- dnakit-mcp      # Codex
# 安装 Skills（或 --agent claude-code、--agent cursor）
dnakit skills install dnakit --agent codex
```

更多安装方式见[安装](installation.md)，Agent 的用法见 [MCP](llm_tools/mcp.md) 和 [Skills](llm_tools/skills.md)。

# 支持与社区

**Github仓库：**[github.com/mapengsen/DNAKit](https://github.com/mapengsen/DNAKit)

# 引用 DNAKit

DNAKit 尚未发布论文或 DOI；当前版本的中文引用信息保存在仓库 `README_CN.md` 中，可引用为：

```text
DNAKit contributors. DNAKit 0.2.0, 2026.
```

引用 DNAKit 不能替代对实际算法、参数集、后端和参考数据库的引用。项目使用的方法、论文、数据库和网站统一列在[致谢与主要来源](acknowledgements.md#acknowledgements)中；结果中的 `Provenance`、`BackendInfo`、`ReferenceLibrary` 和 `RunManifest` 用于保留这些信息。
