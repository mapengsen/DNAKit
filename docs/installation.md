# 安装

DNAKit 要求 Python 3.10 或更高版本。

## 从 PyPI 安装

安装 `0.1.1`：

```bash
pip install dnakit==0.1.1
```

验证安装：

```bash
python -c "import dnakit; print(dnakit.__version__)"
```

## 其他平台

项目仓库已提供以下平台适配文件：

- [Bioconda 配方](https://github.com/mapengsen/DNAKit/tree/main/packaging/bioconda)；
- [GNU Guix 包定义](https://github.com/mapengsen/DNAKit/tree/main/packaging/guix)；
- [Galaxy Tool Shed 包装器](https://github.com/mapengsen/DNAKit/tree/main/galaxy/dnakit)。

这些适配文件固定对应 `0.1.1`，正式安装命令需要等相应平台审核并发布后才能使用。

## 开发环境

在项目根目录创建并激活完整 Conda 开发环境：

```bash
conda env create -f environment-dev.yml
conda activate dnakit-dev
```

开发环境包含测试、类型检查、文档和构建工具。可选功能的额外依赖见
[`pyproject.toml`](https://github.com/mapengsen/DNAKit/blob/main/pyproject.toml)。

# 仓库地址

[github.com/mapengsen/DNAKit](https://github.com/mapengsen/DNAKit)
