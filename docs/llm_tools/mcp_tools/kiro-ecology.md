# AWS Kiro Life Sciences: ecology (`kiro-ecology`)

DNAKit 的 MCP 工具索引收录了这个服务器的 **3** 个 DNA/基因相关工具（服务器共列出 7 个，2026-10-16 读取）。工具在用户自己安装服务器后运行，DNAKit 本身不包含它们。调用方式见 [调用索引中的工具](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found)。下面的工具说明是上游英文原文。

| 项目 | 内容 |
| --- | --- |
| 源码 | [aws-samples/sample-kiro-power-life-sciences @ `95b67e7`](https://github.com/aws-samples/sample-kiro-power-life-sciences/tree/95b67e763df43e94e5bcfd8a230a58553a888cfc) |
| 许可证 | MIT-0 in root LICENSE.txt; MIT declared by component pyproject.toml |
| 启动命令 | `life-sciences-ecology` |
| 传输方式 | stdio |
| 安装 | Prepare the pinned aws-samples/sample-kiro-power-life-sciences source. Separately install life-sciences-common and the life-sciences-ecology component in its Python &gt;=3.10 environment, then set executable_path to the installed life-sciences-ecology console script. No uvx or automatic install is used. |
| 副作用 | Selected tools query external biological resources; molecular-analysis and workflow tools can submit sequences or access local result files. Named credentials are forwarded only for providers that need them. |

| 工具 | 类别 | 简介 |
| --- | --- | --- |
| [`bold_search`](#bold_search) | `sequence` | Searches BOLD Systems DNA barcode records (e.g. COI sequences) by taxon name. |
| [`genbank_env_search`](#genbank_env_search) | `sequence` | Searches GenBank environmental sample sequences by keyword via NCBI E-utilities. |
| [`mgnify_search`](#mgnify_search) | `sequencing_data` | Searches EBI MGnify metagenomics studies by keyword. |

## `bold_search` {#bold_search}

Searches BOLD Systems DNA barcode records (e.g. COI sequences) by taxon name.

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`taxon`
- **关键词：**dna barcode, coi, barcoding, taxon, species identification

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `taxon` | string | 是 | Taxon name. |

??? note "上游说明（英文原文）"

    Search BOLD barcode records by taxon name.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-ecology',
    'bold_search',
    {'taxon': '...'},
)
```

## `genbank_env_search` {#genbank_env_search}

Searches GenBank environmental sample sequences by keyword via NCBI E-utilities.

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`query`
- **关键词：**environmental sequences, genbank, nucleotide, env_sample, metabarcoding

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `query` | string | 是 | Search keyword. |

??? note "上游说明（英文原文）"

    Search GenBank environmental sequences by keyword.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-ecology',
    'genbank_env_search',
    {'query': '...'},
)
```

## `mgnify_search` {#mgnify_search}

Searches EBI MGnify metagenomics studies by keyword.

- **类别：**`sequencing_data`（Sequencing datasets, runs, experiments and read files）
- **必填参数：**`keyword`
- **关键词：**metagenomics, microbiome, amplicon, 16s, ebi

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `keyword` | string | 是 | Search keyword. |

??? note "上游说明（英文原文）"

    Search MGnify metagenomics studies by keyword.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-ecology',
    'mgnify_search',
    {'keyword': '...'},
)
```
