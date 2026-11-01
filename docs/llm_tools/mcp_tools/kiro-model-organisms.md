# AWS Kiro Life Sciences: model-organisms (`kiro-model-organisms`)

DNAKit 的 MCP 工具索引收录了这个服务器的 **5** 个 DNA/基因相关工具（服务器共列出 5 个，2026-10-16 读取）。工具在用户自己安装服务器后运行，DNAKit 本身不包含它们。调用方式见 [调用索引中的工具](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found)。下面的工具说明是上游英文原文。

| 项目 | 内容 |
| --- | --- |
| 源码 | [aws-samples/sample-kiro-power-life-sciences @ `95b67e7`](https://github.com/aws-samples/sample-kiro-power-life-sciences/tree/95b67e763df43e94e5bcfd8a230a58553a888cfc) |
| 许可证 | MIT-0 in root LICENSE.txt; MIT declared by component pyproject.toml |
| 启动命令 | `life-sciences-model-organisms` |
| 传输方式 | stdio |
| 安装 | Prepare the pinned aws-samples/sample-kiro-power-life-sciences source. Separately install life-sciences-common and the life-sciences-model-organisms component in its Python &gt;=3.10 environment, then set executable_path to the installed life-sciences-model-organisms console script. No uvx or automatic install is used. |
| 副作用 | Selected tools query external biological resources; molecular-analysis and workflow tools can submit sequences or access local result files. Named credentials are forwarded only for providers that need them. |

| 工具 | 类别 | 简介 |
| --- | --- | --- |
| [`flybase_gene`](#flybase_gene) | `genome_annotation` | Looks up a Drosophila melanogaster gene by symbol in FlyBase. |
| [`mgi_gene`](#mgi_gene) | `genome_annotation` | Looks up a mouse gene by symbol in Mouse Genome Informatics (MGI). |
| [`sgd_gene`](#sgd_gene) | `genome_annotation` | Looks up a Saccharomyces cerevisiae gene by name in the Saccharomyces Genome Database (SGD). |
| [`wormbase_gene`](#wormbase_gene) | `genome_annotation` | Looks up a Caenorhabditis elegans gene by name in WormBase. |
| [`zfin_gene`](#zfin_gene) | `genome_annotation` | Looks up a zebrafish (Danio rerio) gene by name in ZFIN. |

## `flybase_gene` {#flybase_gene}

Looks up a Drosophila melanogaster gene by symbol in FlyBase.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`symbol`
- **关键词：**drosophila, fruit fly, model organism, fbgn, gene symbol

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `symbol` | string | 是 | Gene symbol (e.g. 'dpp', 'Notch'). |

??? note "上游说明（英文原文）"

    Lookup a gene by symbol in FlyBase (Drosophila).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-model-organisms',
    'flybase_gene',
    {'symbol': '...'},
)
```

## `mgi_gene` {#mgi_gene}

Looks up a mouse gene by symbol in Mouse Genome Informatics (MGI).

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`symbol`
- **关键词：**mouse, mus musculus, mgi id, model organism, gene symbol

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `symbol` | string | 是 | Gene symbol (e.g. 'Trp53', 'Brca1'). |

??? note "上游说明（英文原文）"

    Lookup a gene by symbol in MGI (mouse).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-model-organisms',
    'mgi_gene',
    {'symbol': '...'},
)
```

## `sgd_gene` {#sgd_gene}

Looks up a Saccharomyces cerevisiae gene by name in the Saccharomyces Genome Database (SGD).

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`name`
- **关键词：**yeast, saccharomyces cerevisiae, model organism, orf, gene name

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `name` | string | 是 | Gene name (e.g. 'CDC28', 'RAD51'). |

??? note "上游说明（英文原文）"

    Lookup a gene by name in SGD (yeast).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-model-organisms',
    'sgd_gene',
    {'name': '...'},
)
```

## `wormbase_gene` {#wormbase_gene}

Looks up a Caenorhabditis elegans gene by name in WormBase.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`name`
- **关键词：**c. elegans, nematode, worm, model organism, gene name

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `name` | string | 是 | Gene name (e.g. 'unc-86', 'lin-12'). |

??? note "上游说明（英文原文）"

    Lookup a gene by name in WormBase (C. elegans).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-model-organisms',
    'wormbase_gene',
    {'name': '...'},
)
```

## `zfin_gene` {#zfin_gene}

Looks up a zebrafish (Danio rerio) gene by name in ZFIN.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`name`
- **关键词：**zebrafish, danio rerio, model organism, zdb id, gene name

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `name` | string | 是 | Gene name (e.g. 'shha', 'tp53'). |

??? note "上游说明（英文原文）"

    Lookup a gene by name in ZFIN (zebrafish).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-model-organisms',
    'zfin_gene',
    {'name': '...'},
)
```
