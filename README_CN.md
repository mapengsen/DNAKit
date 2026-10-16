<p align="right">
  <strong>简体中文</strong> | <a href="https://github.com/mapengsen/DNAKit/blob/main/README.md">English</a>
</p>

<p align="center">
  <img alt="DNAKit logo" src="https://raw.githubusercontent.com/mapengsen/DNAKit/main/docs/assets/images/DNAKit-icon.png" width="40%">
</p>

<p align="center">
  <a href="https://github.com/mapengsen/DNAKit/actions/workflows/ci.yml">
    <img src="https://github.com/mapengsen/DNAKit/actions/workflows/ci.yml/badge.svg?branch=main&event=push" alt="CI status">
  </a>
  <a href="https://mapengsen.github.io/DNAKit/">
    <img src="https://img.shields.io/badge/docs-GitHub%20Pages-4051b5" alt="Documentation">
  </a>
  <a href="https://www.python.org/downloads/">
    <img src="https://img.shields.io/badge/Python-3.10%2B-3776AB" alt="Python 3.10+">
  </a>
  <a href="https://github.com/mapengsen/DNAKit/blob/main/LICENSE">
    <img src="https://img.shields.io/github/license/mapengsen/DNAKit" alt="MIT license">
  </a>
</p>

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
</p>

<p align="center">
  <a href="https://github.com/mapengsen/DNAKit/blob/main/docs/llm_tools/mcp.md">
    <img src="https://img.shields.io/badge/MCP-server-1A1A1A?logo=modelcontextprotocol&logoColor=white" alt="MCP server">
  </a>
  <a href="https://github.com/mapengsen/DNAKit/blob/main/docs/llm_tools/skills.md">
    <img src="https://img.shields.io/badge/Skills-agent%20workflows-D97757?logo=data:image/svg%2bxml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCI+PHBhdGggZmlsbD0id2hpdGUiIGQ9Ik0xMyAxIDMgMTRoN2wtMS41IDlMMjAgOS41aC03LjJ6Ii8+PC9zdmc+" alt="Agent Skills">
  </a>
</p>

<p align="center">
  <a href="https://dnakit.mapengsen.asia">
    <img src="https://img.shields.io/badge/Web%20App-online-2EA44F?logo=googlechrome&logoColor=white" alt="Web app">
  </a>
  <a href="https://github.com/mapengsen/DNAKit/blob/main/docs/gui.md">
    <img src="https://img.shields.io/badge/GUI-Windows-0078D6?logo=data:image/svg%2bxml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCI+PHBhdGggZmlsbD0id2hpdGUiIGQ9Ik0wIDBoMTEuNHYxMS40SDB6TTEyLjYgMEgyNHYxMS40SDEyLjZ6TTAgMTIuNmgxMS40VjI0SDB6TTEyLjYgMTIuNkgyNFYyNEgxMi42eiIvPjwvc3ZnPg==" alt="GUI for Windows">
  </a>
  <a href="https://github.com/mapengsen/DNAKit/blob/main/docs/gui.md">
    <img src="https://img.shields.io/badge/GUI-macOS-000000?logo=apple&logoColor=white" alt="GUI for macOS">
  </a>
  <a href="https://github.com/mapengsen/DNAKit/blob/main/docs/gui.md">
    <img src="https://img.shields.io/badge/GUI-Linux-FCC624?logo=linux&logoColor=black" alt="GUI for Linux">
  </a>
</p>

<p align="center">
  <strong>DNAkit: A Comprehensive Toolkit for Efficient DNA Research</strong>
</p>

# DNAKit 是什么？

DNAKit 是面向 DNA 序列的可复现 Python 工具包，定位类似“DNA 领域的 RDKit”。它覆盖标准化、文件与注释格式、基础操作、描述符、模式扫描、热力学、指纹、相似度、聚类与数据划分、综合评价、分子生物学模拟和可视化，并可选用 DNA 基础模型直接预测性质或提取 rep 后进行 k-means 聚类等各种各样的操作。

它提供：

1. **数据准备：** DNA 序列与注释文件读写、数据库下载与查询、合法性检查。
2. **数据处理：** 序列标准化与编辑、转录与翻译、去重与聚类、数据集划分。
3. **数据分析：** DNA 描述符、表征与指纹、序列与模式搜索、理化及热力学性质、深度学习性质预测。
4. **数据评价：** 有效性、唯一性、多样性与新颖性评价，以及相似度、距离和序列比对。
5. **可视化：** 序列文字图、位置高亮、线性与环状 DNA 图、序列比对图。
6. **大语言模型工具：** 提供 [Python API](docs/llm_tools/python_api.md)、[MCP 接入](docs/llm_tools/mcp.md)、[Skills 工作流](docs/llm_tools/skills.md)与内置 [Agent](docs/llm_tools/agent.md)，支持通用大模型调用与自然语言任务执行。
7. **图形界面与网页版：** [网页版](https://dnakit.mapengsen.asia)无需安装，打开浏览器即可使用；[桌面版](docs/gui.md)（`dnakit gui`）支持 Windows、macOS 和 Linux，提供序列与数据集工具、全部函数的检索与运行以及 Skills 目录。

# 网页版与桌面版

不写代码也能使用 DNAKit。网页版和桌面版使用同一套界面：单序列工具（描述符、ORF、motif、序列比较）、
数据集工具（去重、训练/验证/测试集划分、格式转换、描述符特征表）、全部 Agent 可调用函数（附等价的
Python 和 MCP 调用），以及 Skills 目录。

**网页版：** 浏览器打开 [https://dnakit.mapengsen.asia](https://dnakit.mapengsen.asia)，无需安装。网页运行在共享 CPU 上，单次调用限时，
不执行写文件、下载模型或调用外部程序的函数。

**桌面版：** 在自己电脑的窗口中运行同一界面，没有上述限制。

| **系统** | <img src="https://img.shields.io/badge/-Windows%2010%20%7C%2011-0078D6?logo=data:image/svg%2bxml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCI+PHBhdGggZmlsbD0id2hpdGUiIGQ9Ik0wIDBoMTEuNHYxMS40SDB6TTEyLjYgMEgyNHYxMS40SDEyLjZ6TTAgMTIuNmgxMS40VjI0SDB6TTEyLjYgMTIuNkgyNFYyNEgxMi42eiIvPjwvc3ZnPg==" alt="Windows"> | <img src="https://img.shields.io/badge/-macOS-000000?logo=apple&logoColor=white" alt="macOS"> | <img src="https://img.shields.io/badge/-Linux-FCC624?logo=linux&logoColor=black" alt="Linux"> |
|:---|:---:|:---:|:---:|
| **1. 安装** | `pip install "dnakit[gui]"` | `pip install "dnakit[gui]"` | `pip install "dnakit[gui-qt]"` |
| **2. 打开桌面版** | `dnakit gui` | `dnakit gui` | `dnakit gui` |

在 PyPI 发布下一个版本之前，可以从 GitHub 安装：`pip install "dnakit[gui] @ git+https://github.com/mapengsen/DNAKit"`（Linux 用 `gui-qt`）。桌面版在 Windows 上
使用 WebView2，在 macOS 上使用 WebKit；Linux 上使用 GTK WebKit 或 Qt 后端，没有显示器的服务器会改用
浏览器打开。详见[图形界面说明](https://github.com/mapengsen/DNAKit/blob/main/docs/gui.md)。

<p align="center">
  <img alt="DNAKit 桌面版" src="https://raw.githubusercontent.com/mapengsen/DNAKit/main/docs/assets/images/gui-window.png" width="80%">
</p>

# 安装与快速入门

Pypi安装当前版本：

```bash
pip install dnakit==0.2.0
```

如需 VEP、ClinVar、dbSNP、gnomAD、dN/dS 和 Golden Gate
功能，可执行 `pip install "dnakit[external-tools]"` 安装可选后端。

如需让支持 MCP 的 Agent 调用 DNAKit，可执行 `pip install "dnakit[agent]"`，然后在 Agent
中配置 `dnakit-mcp` 命令。

当前源码中的内置 `run_agent()` 还支持按需生成工具、优化说明或兼容修复实现、本地文件结果缓存与断点续跑，以及检索用户指定的知识库目录（RAG）。配置和示例见 [Agent 教程](docs/llm_tools/agent.md#agent-actions)；工具生成与优化要求可写源码检出及显式副作用授权。

当前源码还包含统一的 DNA／基因组学 Skill 目录。先用 `python -m pip install -e "."`
安装本仓库，再用 `dnakit skills search "变异注释"` 搜索，或用
`dnakit skills install dnakit --agent codex` 安装统一入口。
Claude Code、Cursor、MCP 和引用方式见 [Skills 教程](docs/llm_tools/skills.md)。
DNAKit 代码保持 MIT；收录技能保留[各自的许可证和作者署名](THIRD_PARTY_NOTICES.md)。

当前源码还新增了[育种、基因识别、甲基化、DNA 存储与纳米结构接口](docs/api/features/38_genomics_extensions.md)，以及 [外部分析与 MCP 连接](docs/api/features/39_external_genomics.md)。
从本地安装所需依赖，例如 `python -m pip install -e ".[bio-analysis,genomics-genes,agent]"`。
各 DNA 方向的原始项目、作者与尚缺能力见[来源 FAQ](docs/faq/genomics_sources.md)。[335 个补充项目审核](docs/faq/user_catalog.md)逐项记录去重结果、新增接口及未接入原因。

如果需要参与开发，克隆仓库后通过 Conda 创建完整环境：

```bash
conda env create -f environment-dev.yml
conda activate dnakit-dev
```

# 文档

完整文档位于 [DNAKit 文档首页](https://mapengsen.github.io/DNAKit/)，所有功能模块可在 [功能树总览](https://mapengsen.github.io/DNAKit/api/features/function_tree/) 中查看，常见问题见 [FAQ](https://mapengsen.github.io/DNAKit/faq/)。

# 支持与社区

**Github仓库：**

[github.com/mapengsen/DNAKit](https://github.com/mapengsen/DNAKit)

如果有问题、意见或建议，可以先查看：

- [常见问题](https://mapengsen.github.io/DNAKit/faq/)
- [GitHub Issues](https://github.com/mapengsen/DNAKit/issues)

# 更多发布平台

仓库已提供 [Bioconda 与 GNU Guix 配方](https://github.com/mapengsen/DNAKit/tree/main/packaging)
以及 [Galaxy Tool Shed 包装器](https://github.com/mapengsen/DNAKit/tree/main/galaxy/dnakit)。这些文件已固定对应
`0.1.1`；在各平台审核通过前，仍请使用上面的 PyPI 安装方式。

# 更新日志

详细版本更新记录请查看 [CHANGELOG.md](https://github.com/mapengsen/DNAKit/blob/main/CHANGELOG.md)。

# 引用

DNAKit 的项目论文尚未发表。当前版本可引用为：

```text
DNAKit contributors. DNAKit 0.2.0, 2026. https://github.com/mapengsen/DNAKit
```

正式论文发表后，此处将更新为论文引用信息。使用具体算法、数据库或可选后端时，还应引用其对应来源。

当前工作区还提供原生 DNA 编码、生成评价、统计和约束优化，以及可选格式编解码与比对/NGS 专业后端；用法按任务归入[功能教程](docs/api/features/function_tree.md)，依赖与实现约定见[生物信息学说明](docs/faq/bioinformatics.md)。使用 `python -m pip install -e ".[bio]"` 安装。
