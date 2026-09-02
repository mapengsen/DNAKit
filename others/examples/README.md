# DNAKit examples

本目录保存可复现示例的输入和预期工件。

当前包含人工可核验的固定演示夹具和一个严格工作流：

- `fixed_demo.fasta`：三条短 DNA 序列；
- `fixed_demo_expected.json`：固定预期结果。
- `advanced_workflow.yml`：`dnakit-workflow-v1` 固定输入端到端配置。

文档中的 CLI 和 Python 工作流直接使用这些固定输入；预期工件用于人工核验，不会被示例命令覆盖。可执行命令见 `docs/examples/index.md`，端到端 Notebook 见 `others/notebooks/00_skeleton_check.ipynb` 和 `others/notebooks/01_advanced_workflow.ipynb`。
