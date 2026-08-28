<p align="center">
  <img alt="DNAKit logo" src="https://raw.githubusercontent.com/mapengsen/DNAKit/main/docs/assets/images/DNAKit-icon.png" width="40%">
</p>

<p align="center">
  <a href="https://pypi.org/project/dnakit/">
    <img src="https://img.shields.io/pypi/v/dnakit?include_prereleases&amp;label=PyPI" alt="PyPI version">
  </a>
  <a href="https://github.com/mapengsen/DNAKit/actions/workflows/ci.yml">
    <img src="https://github.com/mapengsen/DNAKit/actions/workflows/ci.yml/badge.svg" alt="CI status">
  </a>
  <a href="https://pypi.org/project/dnakit/">
    <img src="https://img.shields.io/pypi/pyversions/dnakit" alt="Supported Python versions">
  </a>
  <a href="https://github.com/mapengsen/DNAKit/blob/main/LICENSE">
    <img src="https://img.shields.io/github/license/mapengsen/DNAKit" alt="MIT license">
  </a>
</p>

# DNAKit 是什么？

DNAKit 是面向 DNA 序列的可复现 Python 工具包，定位类似“DNA 领域的 RDKit”。它覆盖标准化、文件与注释格式、基础操作、描述符、模式扫描、热力学、指纹、相似度、聚类与数据划分、综合评价、分子生物学模拟和可视化，并可选用 DNA 基础模型提取 rep 后进行 k-means 聚类；不集成启动子活性、表达量、TF 结合强度或 CRISPR 编辑效率等任务型深度学习预测模型。

# 安装与快速入门

安装当前开发预览版：

```bash
pip install dnakit==0.1.0.dev0
```

如果需要参与开发，克隆仓库后通过 Conda 创建完整环境：

```bash
conda env create -f environment-dev.yml
conda activate dnakit-dev
```

随后可以阅读 [Python 快速入门指南](https://github.com/mapengsen/DNAKit/blob/main/docs/quickstart.md)。

更详细的依赖、可选后端和安装说明见 [安装文档](https://github.com/mapengsen/DNAKit/blob/main/docs/installation.md)。

# 文档

完整文档位于 [DNAKit 文档首页](https://github.com/mapengsen/DNAKit/blob/main/docs/index.md) 和仓库的 [`docs`](https://github.com/mapengsen/DNAKit/tree/main/docs) 目录。

所有功能模块可在 [功能树总览](https://github.com/mapengsen/DNAKit/blob/main/docs/api/features/function_tree.md) 中查看，常见问题见 [FAQ](https://github.com/mapengsen/DNAKit/blob/main/docs/faq.md)。

# 支持与社区

如果有问题、意见或建议，可以先查看：

- [常见问题](https://github.com/mapengsen/DNAKit/blob/main/docs/faq.md)
- [贡献指南](https://github.com/mapengsen/DNAKit/blob/main/CONTRIBUTING.md)
- [GitHub Issues](https://github.com/mapengsen/DNAKit/issues)

如果发现缺陷或希望增加功能，请按照 [贡献指南](https://github.com/mapengsen/DNAKit/blob/main/CONTRIBUTING.md) 中的要求记录问题和复现信息。

DNAKit 当前处于开发预览阶段，公开接口在正式版本发布前可能调整。
