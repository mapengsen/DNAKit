# 安装

## 从 PyPI 安装

安装 `0.1.4`：

```bash
pip install dnakit==0.1.4
```

地址为：[pypi.org/project/dnakit](https://pypi.org/project/dnakit/)

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

## 深度学习性质预测

安装 Enformer 轨道以及 NT Revised/Genomic Benchmarks 任务 checkpoint 推理依赖：

```bash
python -m pip install "dnakit[neural,neural-enformer]"
```

27 个任务 checkpoint 不随安装包分发，预测时会自动下载当前任务缺失的权重，已有文件直接复用。默认缓存目录为 `./ckpt/enformer-benchmarks/`，也可在调用时使用 `checkpoint_dir` 或
`checkpoint_path` 指定。缓存规则见[权重缓存与自动下载](api/features/23_deep_learning_property_prediction.md#checkpoint-cache)，其他模型的依赖和边界见该页对应说明。

## Agent 与 MCP

如需让 Codex、Claude、Cursor 等支持 MCP 的 Agent 调用 DNAKit：

```bash
python -m pip install "dnakit[agent]"
```

安装后使用 `dnakit-mcp` 启动本地 MCP 服务，具体配置见 [大模型工具 → MCP](llm_tools/mcp.md)。

# 仓库地址

[github.com/mapengsen/DNAKit](https://github.com/mapengsen/DNAKit)

## 生物信息学 Python API 可选依赖

当前工作区可使用 `python -m pip install -e ".[bio]"` 安装全部 Python 分析依赖，或仅安装 `.[bio-analysis]`。MAFFT、BLAST、HMMER、BWA、samtools、bcftools 等可执行程序需单独安装。详见[依赖与实现方式](api/features/27_bioinformatics.md#dependencies)。
