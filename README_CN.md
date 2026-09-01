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
  <strong>DNAkit: A Comprehensive Toolkit for Efficient DNA Research</strong>
</p>

# DNAKit 是什么？

DNAKit 是面向 DNA 序列的可复现 Python 工具包，定位类似“DNA 领域的 RDKit”。它覆盖标准化、文件与注释格式、基础操作、描述符、模式扫描、热力学、指纹、相似度、聚类与数据划分、综合评价、分子生物学模拟和可视化，并可选用 DNA 基础模型直接预测性质或提取 rep 后进行 k-means 聚类等各种各样的操作。

# 安装与快速入门

Pypi安装当前版本：

```bash
pip install dnakit==0.1.3
```

如需 VEP、ClinVar、dbSNP、gnomAD、dN/dS 和 Golden Gate
功能，可执行 `pip install "dnakit[external-tools]"` 安装可选后端。

如需让支持 MCP 的 Agent 调用 DNAKit，可执行 `pip install "dnakit[agent]"`，然后在 Agent
中配置 `dnakit-mcp` 命令。

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
DNAKit contributors. DNAKit 0.1.3, 2026. https://github.com/mapengsen/DNAKit
```

正式论文发表后，此处将更新为论文引用信息。使用具体算法、数据库或可选后端时，还应引用其对应来源。
