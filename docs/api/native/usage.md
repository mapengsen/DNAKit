# 新增函数的调用约定与示例

[484 个函数完整索引](index.md)按用途列出全部新增接口。每个函数有独立条目、当前签名、参数和返回值说明。这里说明通用调用方式，并给出可复制的完整示例。

## 安装与使用入口

这些修改尚未发布到 PyPI。在 DNAKit 源码根目录安装：

```bash
python -m pip install -e "."
```

Python 直接使用 `import dnakit`，不需要初始化 ToolUniverse，也不要求先安装 Skill。科学软件、模型权重、API 凭据和参考数据仍按具体任务准备。来源和许可见[第三方声明](../../acknowledgements.md#third-party-notices)。

## 离线查找函数与参数

下面的代码只读取包内接口目录，不调用数据库或模型：

```python
import dnakit

matches = dnakit.databases.search("gnomad", limit=5, verbose=False)
print([item["name"] for item in matches])
info = dnakit.databases.describe("gnomad_get_sv_by_region", verbose=False)
print(info["input_schema"])
print(len(dnakit.databases.list_operations(verbose=False)))  # 本版本为 450
```

`list_operations()` 中的 450 个操作包括兼容分发入口 `execute_tool`；另外 5 个目录和调用函数构成 `dnakit.databases` 的 455 个公开函数。不能把 450 理解为 450 种独立科学算法。

CLI 提供相同的发现入口：

```bash
dnakit database search "gnomad"
dnakit database info gnomad_get_sv_by_region
```

网页也可以直接搜索 `gnomad_get_sv_by_region`，跳转到[对应函数](variants.md#dnakit.databases.gnomad_get_sv_by_region)。

## 直接查询数据库

以下示例会访问远程数据库。输出取决于服务和查询时间，文档不把示例当作实时结果记录。

```python
import dnakit

gene = dnakit.databases.ensembl_lookup_gene(
    id="BRCA1", species="homo_sapiens",
)
print(gene)
```

需要动态选择函数时，使用参数字典；下面三种调用形式等价，选择一种运行即可：

```python
import dnakit

params = {"id": "BRCA1", "species": "homo_sapiens"}
result = dnakit.databases.call("ensembl_lookup_gene", params)

# 等价形式：
# result = dnakit.databases.execute(
#     {"name": "ensembl_lookup_gene", "arguments": params}
# )
# result = dnakit.databases.execute_tool(
#     tool_name="ensembl_lookup_gene", arguments=params
# )
print(result)
```

`describe()` 的 `input_schema` 是 DNAKit 当前接受的参数约束；`config.parameter` 保留提供方原参数。生成包装的可选参数默认 `None` 时通常表示省略该参数，不代表上游服务的默认值也是 `None`。别名冲突、条件组合、取值范围仍要遵循 Schema。

返回值保留不同提供方的结构，不统一包装为 `result["data"]`。提供方结果 Schema 是结构参考；新增查询分支、分页截断和真实返回字段还需核对实际结果。`call(..., max_result_bytes=...)` 控制序列化结果上限，默认 20 MB。错误结果或超限会抛出异常。

不同服务使用不同坐标和基因组版本，不能把一个接口的坐标约定直接套到另一个接口。ChIP-Atlas 的 `requires_data`、`requires_submission` 状态说明仍需要数据或提交任务，详见[表达与调控接口](regulation.md)。

## 从检索结果下载序列

本例依次使用原有且已扩展的 `entrez`、新增的 `accessions` 和新增的 NCBI 序列下载函数。查询会联网，下载会写入当前目录；实际命中可能变化。

```python
import dnakit

found = dnakit.search.ncbi.entrez(
    "nuccore", organism="Escherichia coli", strain="K-12",
    seq_type="complete_genome", retmax=3,
)
accessions = dnakit.search.ncbi.accessions(found.metadata["uids"])
if accessions:
    downloaded = dnakit.download.ncbi.sequence(
        accessions[0], "ecoli_k12.fasta", format="fasta",
    )
    print(downloaded.files[0].path)
    print(downloaded.metadata["sequence_length"])
    print(downloaded.manifest_path)
else:
    print("没有匹配序列，请核对查询条件。")
```

ENA 使用 INSDC 编号，例如 `dnakit.download.ena.sequence("U00096.3", "ecoli_ena.fasta")`。这与数据库模块的 `ncbi_get_sequence` 不同：下载命名空间把文件和来源清单保存到磁盘，数据库包装按对应服务返回记录。格式转换依赖与坐标说明见[下载教程](../features/15_download.md#accession-sequences)。

## ENCODE 查询与批量下载预览

```python
import json
import dnakit

experiments = json.loads(dnakit.encode.search_experiments(
    assay_title="ATAC-seq", organism="Homo sapiens", limit=3,
))
print(experiments)

preview = json.loads(dnakit.encode.batch_download(
    download_dir="downloads/encode", assay_title="ATAC-seq",
    assembly="GRCh38", file_format="bed", limit=3, dry_run=True,
))
print(preview)
```

`dry_run=True` 查询待下载文件但不执行文件下载；设为 `False` 才下载。ENCODE 的 20 个函数见[专页](encode.md)。多数返回 JSON 文本，引用和表格导出也可返回 BibTeX、RIS、CSV、TSV，不能对这些非 JSON 文本使用 `json.loads()`。追踪、凭据和文件下载会有相应的本地写入。

## 执行随包坐标转换脚本

此例在当前 Python 环境运行随包脚本，把 BED 0-based 半开坐标转换成 GFF 1-based 闭区间；它不是基因组版本转换。

```python
import dnakit

print(dnakit.workflows.scripts("genomic-coordinates", verbose=False))
path = dnakit.workflows.resource_path(
    "genomic-coordinates", "scripts/convert_coords.py", verbose=False,
)
print(path)
result = dnakit.workflows.run(
    "genomic-coordinates", "scripts/convert_coords.py",
    arguments=["--from", "bed", "--to", "gff", "chr1", "999", "1000", "--format", "json"],
    timeout_seconds=60, verbose=False,
)
if result["returncode"] != 0:
    raise RuntimeError(result["output"])
print(result["output"])
```

`run_support()` 运行已打包的支持流程，允许的名称为 `clinical-variant-reporter`、`clinical-trial-finder`、`omics-target-evidence-mapper`、`gwas-lookup`、`fine-mapping`。先通过 `arguments=["--help"]` 查看各自参数；脚本导入的科学依赖仍需准备。参数、资源校验和返回值详见[工作流接口](execution.md)。

## 独立位点遗传模拟

下面使用一个显式的合成位点，不访问远程服务。父方为 A/A、母方为 C/C，关闭突变后每个后代在该位点得到 A/C。

```python
import dnakit

locus = {
    "chromosome": "1", "position": 1, "ref": "A", "alt": "C",
    "dominance": "additive", "effect_size": 1,
}
father = {"id": "father", "sex": "Male", "loci": {"L1": {**locus, "alleles": ["A", "A"]}}}
mother = {"id": "mother", "sex": "Female", "loci": {"L1": {**locus, "alleles": ["C", "C"]}}}
traits = {"traits": {"dosage": {"loci": [dict(locus, id="L1", effect=1)]}}}
result = dnakit.workflows.simulate_inheritance(
    father, mother, traits, {"diseases": {}},
    offspring_count=2, mutation_rate=0, seed=13, verbose=False,
)
print([child["loci"]["L1"]["alleles"] for child in result["offspring"]])
# [['A', 'C'], ['A', 'C']]
```

该函数模拟独立位点的分离与突变，不模拟染色体交叉，不用于推断真实临床结局。[完整参数](analysis.md#dnakit.workflows.simulate_inheritance)。

## 从已有报告生成 PDF

安装报告渲染依赖：

```bash
python -m pip install -e ".[reports]"
```

下面的示例使用临时目录，展示真实的输入、输出和来源记录：

```python
from pathlib import Path
from tempfile import TemporaryDirectory
import dnakit

with TemporaryDirectory() as directory:
    root = Path(directory)
    report = root / "report.md"
    report.write_text("# Results\nUser-supplied report content.\n", encoding="utf-8")
    result = dnakit.workflows.write_wes_report(report, root / "pdf", verbose=False)
    print(result["output"])
    print(result["manifest"])
    print(Path(result["output"]).exists())  # True
```

实际保存时使用自己的输出目录。该函数渲染已有内容，不增加医学解释；`language` 只切换英文或西班牙文标签，不翻译正文。标题、Logo、字体和覆盖选项见[完整参数](execution.md#dnakit.workflows.write_wes_report)。

## MCP 与验证范围

同一函数也可由 MCP 的工具搜索、参数描述和调用入口使用，详见[MCP 教程](../../llm_tools/mcp.md)。Python、CLI、MCP 是同一能力的不同入口，统计时不重复计数。

文档校验检查全部 484 个函数是否有唯一条目、签名及参数是否与源码一致，并构建中英文网页。文档覆盖与离线示例验证不等于所有远程服务、模型和完整科学流程都已端到端实测。
