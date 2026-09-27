# pip 基因组计算集成：本轮变更

版本说明：新增 `genotype_pca()` 和 `average_nucleotide_identity()`，前者独立实现明确的二倍体剂量 PCA，后者通过 PyFastANI 进程内调用 FastANI。新增 `genomics` extra，一并安装 NumPy、Pyrodigal、pysam、PyFastANI；现有基因预测、关系矩阵、育种预测和 BAM 统计可复用。新增代码提供阶段进度，并完成独立代码审查。尚未发布到 PyPI。

43 项外部分析已逐项分类；本轮落地两项新 API，其余迁移建议并非已经实现。原有 `run_analysis()` 的输入与输出约定继续适用。

## 验证

- 新接口 80 项测试、相关套件 251 项测试通过；PCA 另在 Python 3.14.6 通过 51 项测试。
- 真实 NumPy、PyFastANI、Pyrodigal、pysam 的四项人工数据调用通过。
- wheel 与新运行时代码一致，在全新 Python 3.11 环境安装 `genomics` extra 后，从源码目录之外、清空 PATH/PYTHONPATH 执行验证成功。
- 两项新 API 经 DNAKit 工具目录实际 JSON 调用通过，可由统一 Agent/MCP 入口发现。
- Ruff、严格类型检查、中英文 MkDocs 严格构建通过；CI 已加入 `genomics` 依赖与真实后端验证。
- 审查修复：ANI 自比属于近似估计，验证容差设为 0.01 个百分点；进程守卫证据准确限定为拦截 `subprocess.Popen`。

完整机器可读记录：[pip-genomics-validation.json](../../registry/pip-genomics-validation.json)。上述验证不代表真实生物样本准确性或所有平台已验证。

## 执行示例

在 Python ≥3.10 的环境、项目根目录安装当前源码及组合依赖：

```bash
python -m pip install -e ".[genomics]"
```

运行 PCA、ANI、CDS 预测和 BAM 统计的真实计算验证，显示进度并写入新 JSON 文件：

```bash
python others/validation/pip_genomics.py --output results/pip_genomics.json
```

运行无需输入文件的基因型 PCA 示例：

```bash
python -c "from dnakit.genomics import genotype_pca; print(genotype_pca([[0, 0], [1, 2], [2, 1]]))"
```

## 本轮全部新增／修改文件

共 25 个源码、配置、文档与验证记录文件。只列本轮改动；此前工作区中的其他修改以及临时安装／构建文件不计入。

| 状态 | 文件 |
| --- | --- |
| 修改 | [pyproject.toml](../../pyproject.toml) |
| 修改 | [src/dnakit/genomics/__init__.py](../../src/dnakit/genomics/__init__.py) |
| 修改 | [docs/api/features/38_genomics_extensions.md](../../docs/api/features/38_genomics_extensions.md) |
| 修改 | [docs/api/features/38_genomics_extensions.en.md](../../docs/api/features/38_genomics_extensions.en.md) |
| 修改 | [docs/faq/api/38_genomics_extensions.md](../../docs/faq/api/38_genomics_extensions.md) |
| 修改 | [docs/faq/api/38_genomics_extensions.en.md](../../docs/faq/api/38_genomics_extensions.en.md) |
| 修改 | [docs/api/features/39_external_genomics.md](../../docs/api/features/39_external_genomics.md) |
| 修改 | [docs/api/features/39_external_genomics.en.md](../../docs/api/features/39_external_genomics.en.md) |
| 修改 | [docs/installation.md](../../docs/installation.md) |
| 修改 | [docs/installation.en.md](../../docs/installation.en.md) |
| 修改 | [docs/api/features/function_tree.md](../../docs/api/features/function_tree.md) |
| 修改 | [docs/api/features/function_tree.en.md](../../docs/api/features/function_tree.en.md) |
| 修改 | [docs/faq.md](../../docs/faq.md) |
| 修改 | [docs/faq.en.md](../../docs/faq.en.md) |
| 修改 | [.github/workflows/ci.yml](../../.github/workflows/ci.yml) |
| 修改 | [mkdocs.yml](../../mkdocs.yml) |
| 新增 | [src/dnakit/genomics/ani.py](../../src/dnakit/genomics/ani.py) |
| 新增 | [src/dnakit/genomics/population.py](../../src/dnakit/genomics/population.py) |
| 新增 | [others/tests/unit/test_genomics_ani.py](../../others/tests/unit/test_genomics_ani.py) |
| 新增 | [others/tests/unit/test_genomics_population.py](../../others/tests/unit/test_genomics_population.py) |
| 新增 | [others/validation/pip_genomics.py](../../others/validation/pip_genomics.py) |
| 新增 | [docs/faq/external_analysis_migration.md](../../docs/faq/external_analysis_migration.md) |
| 新增 | [docs/faq/external_analysis_migration.en.md](../../docs/faq/external_analysis_migration.en.md) |
| 新增 | [registry/pip-genomics-validation.json](../../registry/pip-genomics-validation.json) |
| 新增 | [others/docs/pip_genomics_changes.md](../../others/docs/pip_genomics_changes.md) |
