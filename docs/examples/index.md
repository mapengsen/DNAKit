# 示例

## 固定数据夹具

仓库提供两份固定夹具：

- `examples/fixed_demo.fasta`：三条短 DNA 序列，其中两条完全相同；
- `examples/fixed_demo_expected.json`：人工核验的预期结果。

文档站的[固定输入演示](../demo/index.md)使用同一份数据。演示是只读页面，不包含 DNA 序列输入框、上传控件或在线计算后端。

## 可执行 CLI 工作流

以下命令依次完成格式转换、精确去重和可复现随机划分；需要输入顺序无关时可将 CLI 的 `--method` 改为 `hash`：

```bash
DNAKIT_DEMO_DIR="$(mktemp -d)"

dnakit convert \
  examples/fixed_demo.fasta \
  "$DNAKIT_DEMO_DIR/input.json" \
  --output-format json \
  --no-progress

dnakit deduplicate \
  examples/fixed_demo.fasta \
  "$DNAKIT_DEMO_DIR/nonredundant.fasta" \
  --equivalence exact

dnakit split \
  "$DNAKIT_DEMO_DIR/nonredundant.fasta" \
  "$DNAKIT_DEMO_DIR/split" \
  --ratios train=0.5,test=0.5 \
  --seed 7

find "$DNAKIT_DEMO_DIR" -maxdepth 2 -type f -print
```

预期结果：转换得到 3 条记录，精确去重保留 2 条记录，划分得到 1 条训练记录和 1 条测试记录，并生成划分分配信息。命令写入临时目录，不会修改输入夹具。

## Python 工作流

[快速入门](../quickstart.md)提供可直接运行的内存工作流，覆盖读取、标准化、序列操作、描述符、指纹、相似度、去重、划分和写出。主要模块的公开入口见 [API 参考](../api/index.md)。

## 严格 YAML 工作流

仓库文件 `examples/advanced_workflow.yml` 使用 schema `dnakit-workflow-v1`，只允许以下 8 种步骤：`normalize`、`validate`、`descriptors`、`fingerprint`、`deduplicate`、`split`、`write`、`report`。根字段为 `schema_version`、`run_id`、`input`、`output_dir`、`seed`、`error_policy`、`overwrite`、`limits` 和 `steps`；每步只有 `id`、`operation`、`input`、`params`。

以下命令先复制固定输入和配置，因此所有输出都在临时目录：

```bash
DNAKIT_WORK_DIR="$(mktemp -d)"
cp examples/fixed_demo.fasta examples/advanced_workflow.yml "$DNAKIT_WORK_DIR/"
PYTHONNOUSERSITE=1 dnakit workflow \
  "$DNAKIT_WORK_DIR/advanced_workflow.yml" \
  --no-progress
find "$DNAKIT_WORK_DIR/workflow-output" -maxdepth 3 -type f -print
```

成功后生成 `splits/train.fasta`、`splits/test.fasta`、`report.html`、`run-manifest.json` 和专用目录 marker。manifest 保存可复现的兼容命令、配置/输入 SHA-256、resolved config、seed、软件版本、每步参数/状态/耗时/摘要/工件。

只验证计划且不写文件：

```bash
DNAKIT_DRY_RUN_DIR="$(mktemp -d)"
cp examples/fixed_demo.fasta examples/advanced_workflow.yml "$DNAKIT_DRY_RUN_DIR/"
PYTHONNOUSERSITE=1 dnakit workflow \
  "$DNAKIT_DRY_RUN_DIR/advanced_workflow.yml" \
  --dry-run \
  --no-progress
```

`--resume` 会重算内存步骤，仅在配置 SHA、原始输入 SHA、步骤输入状态 SHA 和输出路径/大小/SHA 全部匹配时跳过 `write`/`report`。它不是通用缓存。

## SVG 可视化示例

```python
from pathlib import Path
from tempfile import TemporaryDirectory

from dnakit import DNASequence, Gap
from dnakit.visualization import Highlight, plot_sequence, save_svg

sequence = DNASequence(["ACGT", Gap(100), "TGCA"])
artifact = plot_sequence(
    sequence,
    highlights=[Highlight(1, 4, label="left region")],
)
with TemporaryDirectory() as directory:
    saved = save_svg(artifact, Path(directory) / "sequence.svg")

print(saved.target_artifact.sha256)
```

已知 Gap 显示为 `[100 bp]`，未知 Gap 显示为 `[… bp]`，两者都不会为绘图而展开成大量字符。`save_svg()` 默认拒绝覆盖，且返回字节数、SHA-256 和目标 artifact 信息。

高级 API 端到端示例见 `notebooks/01_advanced_workflow.ipynb`：它实际运行热力学、固定 16 维热力学指纹、sketch/index、阈值/层次聚类、代表选择、leakage、reference-scoped novelty/memorization、synthesis-risk 和引物性质/设计。Primer3 只在环境中已安装时由示例显式执行；Notebook 不执行条件 NUPACK adapter。下载结构的 PDB/DBN/DSSR 本地分析由 `examples/analyze_dna_structures.py` 完成。
