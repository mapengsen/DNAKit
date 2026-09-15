# Skill 统一接口

如果需要直接使用函数，请查看[新增 484 个函数的完整索引](../api/native/index.md)和[调用示例](../api/native/usage.md)。全部函数已按下载、查询、分析、模型和工程化用途归类，每项都有独立参数与返回值说明。

DNAKit 为整个技能目录提供统一的搜索、安装和执行入口。集合封装中的数据库访问、共享脚本、报告与调用代码迁入 DNAKit；Biopython、samtools、AlphaGenome 等实际科学软件或模型继续使用原实现。

当前修改尚未发布到 PyPI。在源码根目录安装本版本：

```bash
python -m pip install -e "."
dnakit skills search "变异注释" --limit 5
dnakit skills info variant-annotation
dnakit skills doctor variant-annotation
dnakit skills install variant-annotation --agent codex
```

发布包含这些功能的版本后，可以用 `python -m pip install dnakit` 替代源码安装。Skill 名称去掉集合前缀，旧名称保留为别名；真实软件名称继续保留。例如 `tooluniverse-sequence-retrieval` 对应 `sequence-retrieval`，`virulencefinder` 和 `alphagenome-atlas-website-links` 保持原名。

## 数据库函数

先搜索所需接口，再查看参数：

```bash
dnakit database search "gene lookup"
dnakit database info ensembl_lookup_gene
dnakit database call ensembl_lookup_gene '{"id":"BRCA1","species":"homo_sapiens"}'
```

Python 直接使用能力函数，无需初始化 ToolUniverse 或加载全部工具：

```python
import dnakit

gene = dnakit.databases.ensembl_lookup_gene(
    id="BRCA1",                # 也接受原参数 gene_id
    species="homo_sapiens",
)
variants = dnakit.databases.gnomad_get_sv_by_region(
    chrom="17", start=43044295, end=43070295,  # end 也可写成 stop
)
print(gene)
print(variants)
```

`dnakit.databases.describe()` 的 `input_schema` 描述 DNAKit 接受的参数，`config.parameter` 保留原服务参数。别名出现冲突、参数不合法、提供方报错或结果超出上限时抛出异常。结果保留提供方的字典、列表和字段结构，不能假定每个接口都有相同的 `data` 层。

可选参数会实际参与处理。基因名称可按接口要求解析为 ID，坐标别名会转换后调用；Open Targets 的 `size` 限制当前返回页，并标记剩余结果，GWAS 的 `p_value_threshold` 对当前页实施过滤。当前页结果不能被视为完整数据库检索结果。

## 随包脚本

用户不需要克隆原技能仓库来定位已打包的脚本：

```bash
dnakit skills scripts genomic-coordinates
dnakit skills run genomic-coordinates --script scripts/convert_coords.py -- \
  --from bed --to gff chr1 999 1000 --format json
```

Python 调用同一份脚本，保留原来的可选参数：

```python
import dnakit

result = dnakit.workflows.run(
    "genomic-coordinates",
    "scripts/convert_coords.py",
    arguments=["--from", "bed", "--to", "gff", "chr1", "999", "1000", "--format", "json"],
    timeout_seconds=60,
)
print(result["returncode"])
print(result["output"])
```

该示例把 BED 的 0-based 半开区间转换成 GFF 的 1-based 闭区间；不是参考基因组版本转换。CLI 转发标准输入、输出和退出码；Python 接口限制执行时间和捕获输出大小。脚本及其辅助资源在执行前核对已登记的 SHA-256。流程依赖的科学软件仍须位于运行环境中。

## ENCODE、组合流程与报告

ENCODE 的查询、下载和本地追踪使用同一个 DNAKit 包：

```python
import json
import dnakit

experiments = json.loads(dnakit.encode.search_experiments(
    assay_title="ATAC-seq", organism="Homo sapiens", limit=3,
))
print(experiments)
```

ENCODE 函数保留原接口选项和 JSON 文本返回形式。本地追踪及凭据使用 `~/.dnakit-encode/`；持久化凭据存储需要 `dnakit[encode-credentials]`。源代码许可为 AGPL-3.0-only，详见第三方声明。

nf-core 等组合流程需要的已许可支持脚本可通过 `dnakit skills support-run NAME -- ...` 或 `dnakit.workflows.run_support()` 执行。支持步骤包括 `clinical-variant-reporter`、`clinical-trial-finder` 、`omics-target-evidence-mapper`、`gwas-lookup` 和 `fine-mapping`。这些函数运行原流程，不保证其结论已经通过临床验证。

DNAKit 自己的报告函数接收已有 Markdown 内容，生成 PDF 和来源哈希，不自动生成新的医学解释：

```bash
python -m pip install -e ".[reports]"
dnakit skills report --input wes_report.md --output report-pdf --language en
```

对应函数为 `dnakit.workflows.write_wes_report()`，提供语言标签、标题、左右 Logo、字体和覆盖选项。没有复制上游标为专有许可的英文或西班牙文报告实现。

## MCP 与验证范围

按 [MCP 教程](mcp.md)安装 `dnakit[agent]` 后，Agent 使用 `search_dnakit_tools` → `describe_dnakit_tool` → `call_dnakit_tool` 发现并调用上述函数。脚本执行、模型计算、文件写入等操作带有相应的操作类型；目录增加函数不会把全部定义放入模型的初始上下文。

每个 Skill 的 `DNAKIT.md`、运行依赖和验证状态可以通过 `dnakit skills info ID --json` 查看。`registry/migration.json` 记录逐项接入方式，`registry/adaptations.json` 绑定维护模板和上游版本。函数签名、文件完整性、单元测试、真实服务调用及完整科研流程属于不同验证层级；未经执行的流程继续标为 `not_run`。

移入的代码保留原作者、许可及原始文件哈希。DNAKit 新增的参数转换、执行入口和报告函数另行记录，不把上游模型或代码改称原创。RCSB 坐标文件接口依据 [RCSB 官方下载说明](https://www.rcsb.org/docs/programmatic-access/file-download-services)，CELLxGENE 有界查询依据 [官方 Python API](https://chanzuckerberg.github.io/cellxgene-census/_autosummary/cellxgene_census.get_anndata.html)。

ChIP-Atlas 的基因区域查询先返回 `requires_data` 状态及峰文件下载链接；提供核对过的 `peak_file` 后才计算区域重叠。原富集接口只提供提交说明，标为 `requires_submission`，不能当成已算出的统计结果。区域重叠也不能直接证明 motif 被破坏，详见 [ChIP-Atlas 官方 API](https://chip-atlas.org/agents)。
