# AWS Kiro Life Sciences: epigenomics (`kiro-epigenomics`)

DNAKit 的 MCP 工具索引收录了这个服务器的 **3** 个 DNA/基因相关工具（服务器共列出 3 个，2026-10-16 读取）。工具在用户自己安装服务器后运行，DNAKit 本身不包含它们。调用方式见 [调用索引中的工具](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found)。下面的工具说明是上游英文原文。

| 项目 | 内容 |
| --- | --- |
| 源码 | [aws-samples/sample-kiro-power-life-sciences @ `95b67e7`](https://github.com/aws-samples/sample-kiro-power-life-sciences/tree/95b67e763df43e94e5bcfd8a230a58553a888cfc) |
| 许可证 | MIT-0 in root LICENSE.txt; MIT declared by component pyproject.toml |
| 启动命令 | `life-sciences-epigenomics` |
| 传输方式 | stdio |
| 安装 | Prepare the pinned aws-samples/sample-kiro-power-life-sciences source. Separately install life-sciences-common and the life-sciences-epigenomics component in its Python &gt;=3.10 environment, then set executable_path to the installed life-sciences-epigenomics console script. No uvx or automatic install is used. |
| 副作用 | Selected tools query external biological resources; molecular-analysis and workflow tools can submit sequences or access local result files. Named credentials are forwarded only for providers that need them. |

| 工具 | 类别 | 简介 |
| --- | --- | --- |
| [`ihec_search`](#ihec_search) | `regulation_epigenome` | Searches International Human Epigenome Consortium (IHEC) epigenome datasets by tissue type. |
| [`methbase_search`](#methbase_search) | `regulation_epigenome` | Searches MethBase DNA methylomes by species and optional tissue. |
| [`roadmap_search`](#roadmap_search) | `regulation_epigenome` | Searches Roadmap Epigenomics data files by tissue and optional histone mark (e.g. H3K4me3). |

## `ihec_search` {#ihec_search}

Searches International Human Epigenome Consortium (IHEC) epigenome datasets by tissue type.

- **类别：**`regulation_epigenome`（Regulatory elements, TF binding, chromatin, methylation and 3D genome）
- **必填参数：**`tissue`
- **关键词：**reference epigenome, histone, methylation, chip-seq, tissue

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `tissue` | string | 是 | Tissue type to search for. |

??? note "上游说明（英文原文）"

    Search IHEC datasets by tissue type.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-epigenomics',
    'ihec_search',
    {'tissue': '...'},
)
```

## `methbase_search` {#methbase_search}

Searches MethBase DNA methylomes by species and optional tissue.

- **类别：**`regulation_epigenome`（Regulatory elements, TF binding, chromatin, methylation and 3D genome）
- **必填参数：**`species`
- **关键词：**dna methylation, methylome, bisulfite, cpg, wgbs

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `species` | string | 是 | Species name. |
| `tissue` | string | 否 | Tissue type. 默认 `` |

??? note "上游说明（英文原文）"

    Search MethBase methylomes by species and optional tissue.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-epigenomics',
    'methbase_search',
    {'species': '...'},
)
```

## `roadmap_search` {#roadmap_search}

Searches Roadmap Epigenomics data files by tissue and optional histone mark (e.g. H3K4me3).

- **类别：**`regulation_epigenome`（Regulatory elements, TF binding, chromatin, methylation and 3D genome）
- **必填参数：**`tissue`
- **关键词：**histone mark, h3k4me3, chromatin, epigenome, tissue

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `tissue` | string | 是 | Tissue type to search for. |
| `mark` | string | 否 | Histone mark (e.g. 'H3K4me3'). 默认 `` |

??? note "上游说明（英文原文）"

    Search Roadmap Epigenomics data by tissue and optional histone mark.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-epigenomics',
    'roadmap_search',
    {'tissue': '...'},
)
```
