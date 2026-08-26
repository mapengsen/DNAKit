# 教程

仓库中的 Notebook 只调用真实公开 API，使用固定小型输入或临时目录，不下载数据库、不访问网络。当前 Notebook 不调用条件 NUPACK adapter；公开二级/三维结构示例在对应 API 页面和 `examples/analyze_dna_structures.py`。

## Notebook 目录

| Notebook | 内容 | 外部条件 |
| --- | --- | --- |
| `notebooks/00_skeleton_check.ipynb` | 核心对象、FASTA/JSON I/O、描述符、基础指纹/相似度、去重、固定 seed 划分和 SVG | 仅核心+docs |
| `notebooks/01_advanced_workflow.ipynb` | 模式/热力学、固定 16 维指纹、sketch、层次/阈值聚类、leakage、参考库 novelty/memorization、synthesis-risk、引物性质/设计、HTML/图片导出和运行摘要 | 图片单元需要 `viz` extra；Primer3 仅在已安装时由示例显式执行；不包含 NUPACK 能力 |

执行全部 Notebook：

```bash
PYTHONNOUSERSITE=1 python -m pytest --nbmake notebooks
```

## 配置工作流

`examples/advanced_workflow.yml` 是固定输入、白名单步骤的本地 pipeline。它读取同目录的 `fixed_demo.fasta`，依次执行标准化、验证、描述符、指纹、去重、划分、写出和报告，并生成带 SHA-256 的 `run-manifest.json`。为避免在源码目录生成工件，先把两个示例文件复制到临时目录：

```bash
DNAKIT_WORK_DIR="$(mktemp -d)"
cp examples/fixed_demo.fasta examples/advanced_workflow.yml "$DNAKIT_WORK_DIR/"
PYTHONNOUSERSITE=1 dnakit workflow \
  "$DNAKIT_WORK_DIR/advanced_workflow.yml" \
  --no-progress
```

配置内的 `output_dir` 相对配置文件解析；CLI 只接受 `--dry-run`、`--resume` 和 `--no-progress`，没有 `--output-dir` 或 `--cache-dir`。兼容/开发入口是 `python -m dnakit.cli.workflow run CONFIG`。未知步骤/字段、输出路径逃逸和超过资源上限的请求都会被拒绝。完整 schema 和输出见[示例页面](../examples/index.md)。

## 推荐阅读顺序

1. [快速入门](../quickstart.md)：理解不可变对象、显式策略和结果审计。
2. `00_skeleton_check.ipynb`：跑通最小端到端流程。
3. `01_advanced_workflow.ipynb`：连接高级模块和评价。
4. [API 参考](../api/index.md)：按模块查看完整参数。
5. [验证与 benchmark](../validation.md)：核对数值模型、容差和环境边界。

## 教程不会声称的内容

- Primer3 adapter 对照不等于独立热力学算法验证；
- 没有真实 NUPACK 二级结构数值验证；
- novelty 不脱离参考库定义；
- synthesis-risk、PCR 和 assembly 不是实验成功概率；
- benchmark 只比较本机已配置的 DNAKit/Biopython 对等任务，不是跨机器或任意外部工具排名；
- 本地 wheel/sdist 构建不等于 TestPyPI/PyPI 发布。
