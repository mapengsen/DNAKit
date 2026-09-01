# 安装

DNAKit 要求 Python 3.10 或更高版本。

## 从 PyPI 安装

安装 `0.1.3`：

```bash
pip install dnakit==0.1.3
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

## 可选生物信息后端

安装 VEP、ClinVar、dbSNP、gnomAD、dN/dS 和
Golden Gate 功能所需的可选后端：

```bash
python -m pip install "dnakit[external-tools]"
```

默认的 `pip install dnakit` 不会安装或加载该后端。

## Agent 与 MCP

如需让 Codex、Claude、Cursor 等支持 MCP 的 Agent 调用 DNAKit：

```bash
python -m pip install "dnakit[agent]"
```

安装后使用 `dnakit-mcp` 启动本地 MCP 服务，具体配置见 [Agent 与 MCP 工具](agent_tools.md)。

# 仓库地址

[github.com/mapengsen/DNAKit](https://github.com/mapengsen/DNAKit)
