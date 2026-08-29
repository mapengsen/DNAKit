<p align="center">
  <img alt="DNAKit logo" src="https://raw.githubusercontent.com/mapengsen/DNAKit/main/docs/assets/images/DNAKit-icon.png" width="40%">
</p>

<p align="center">
  <a href="https://pypi.org/project/dnakit/">
    <img src="https://img.shields.io/pypi/v/dnakit?include_prereleases=true&amp;label=PyPI&amp;cacheSeconds=300" alt="PyPI version">
  </a>
  <a href="https://github.com/mapengsen/DNAKit/actions/workflows/ci.yml">
    <img src="https://github.com/mapengsen/DNAKit/actions/workflows/ci.yml/badge.svg?branch=main&amp;event=push" alt="CI status">
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
  <strong>DNAkit: A Comprehensive Toolkit for Efficient DNA Research</strong>
</p>

# DNAKit 是什么？

DNAKit 是面向 DNA 序列的可复现 Python 工具包，定位类似“DNA 领域的 RDKit”。它覆盖标准化、文件与注释格式、基础操作、描述符、模式扫描、热力学、指纹、相似度、聚类与数据划分、综合评价、分子生物学模拟和可视化，并可选用 DNA 基础模型提取 rep 后进行 k-means 聚类等各种各样的操作。

# 安装与快速入门

安装当前版本：

```bash
pip install dnakit==0.1.1
```

如果需要参与开发，克隆仓库后通过 Conda 创建完整环境：

```bash
conda env create -f environment-dev.yml
conda activate dnakit-dev
```

随后可以阅读 [Python 快速入门指南](https://mapengsen.github.io/DNAKit/quickstart/)。

# 文档

完整文档位于 [DNAKit 文档首页](https://mapengsen.github.io/DNAKit/)，所有功能模块可在 [功能树总览](https://mapengsen.github.io/DNAKit/api/features/function_tree/) 中查看，常见问题见 [FAQ](https://mapengsen.github.io/DNAKit/faq/)。

# 支持与社区

**Github仓库：**

[github.com/mapengsen/DNAKit](https://github.com/mapengsen/DNAKit)

如果有问题、意见或建议，可以先查看：

- [常见问题](https://mapengsen.github.io/DNAKit/faq/)
- [GitHub Issues](https://github.com/mapengsen/DNAKit/issues)

# 更新日志

详细版本更新记录请查看 [CHANGELOG.md](https://github.com/mapengsen/DNAKit/blob/main/CHANGELOG.md)。

# 引用

DNAKit 的项目论文尚未发表。当前版本可引用为：

```text
DNAKit contributors. DNAKit 0.1.1, 2026. https://github.com/mapengsen/DNAKit
```

正式论文发表后，此处将更新为论文引用信息。使用具体算法、数据库或可选后端时，还应引用其对应来源。
