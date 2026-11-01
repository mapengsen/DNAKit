# 搜索第三方 MCP 服务器里的 DNA 与基因工具

适用读者：想找到只存在于第三方 MCP 服务器中的功能的 DNAKit 用户和 Agent 开发者。

`list_mcp_servers()` 列出 49 个第三方 MCP 服务器，但以前 Agent 只能看到每个服务器的名称和能力标签，服务器里有哪些工具，要先安装并启动它才知道。现在 DNAKit 随包提供这些工具的离线索引，只保留处理 DNA 或基因的工具，`search_mcp_tools()` 和 Agent 服务器的 `search_dnakit_tools` 不用启动任何服务器就能搜到它们。

每个收录工具的参数和通过 DNAKit 调用的写法见 [第三方 MCP 工具目录](../../llm_tools/mcp_tools/index.md)。

## 索引的建立方式

- **列出工具。** 2026-10-16，49 个服务器都按固定提交安装，并请求 `tools/list`。46 个实际启动后列出（stdio 或官方托管的 HTTP 端点）。`ncbi-database`、`seqeron-sequence` 和 `grebi` 从固定提交的源码读取：第一个在该提交下 MCP 握手失败，本机没有 .NET，grebi 的托管端点返回 404。
- **隐藏的工具。** 有两个服务器在 `tools/list` 里只暴露一个分发工具。ToolUniverse 默认以精简模式启动（5 个工具），2,922 个底层工具通过 `execute_tool` 调用；genomi 的 75 个能力操作通过 `genomi.invoke` 调用。这两份注册表都从固定版本的包中读取，并在索引里记下所需的分发调用方式。
- **筛选。** 列出的 4,005 个工具逐个阅读了说明和输入 schema，只保留主要输入或输出是 DNA 或基因的工具，共 1,402 个。其余 2,603 个连同原因一起记录：服务器自身管理、文献检索、药物与化学、只涉及蛋白质、非遗传的临床数据、通用本体查询、单细胞的细胞层面步骤、通用平台管理，或与主题无关。
- 保留的工具分为十类：`sequence`、`genome_annotation`、`variant`、`association`、`expression`、`regulation_epigenome`、`sequencing_data`、`prediction_model`、`molecular_tools` 和 `breeding_population`。

索引是一份快照。调用工具仍需安装对应服务器；`call_mcp_tool` 会按已安装服务器的实时 schema 校验参数。

## 1) `dnakit.integrations.mcp_tool_index_summary()`

- **作用：** 统计列出、保留和排除的 MCP 工具数，按类别和按服务器分别给出。
- **处理方法：** 读取随包的索引文件，不启动服务器，不联网。
- **API：** `dnakit.integrations.mcp_tool_index_summary(verbose[可选])`。
- **输入：** 除 `verbose` 外无。
- **输出：** 字典，包含 `catalog_servers`、`listed_servers`、`servers_with_indexed_tools`、`listed_tools`、`indexed_tools`、`excluded_tools`、`indexed_by_category`、`excluded_by_reason`、`categories` 和 `servers`。`servers` 的每一项有 `server_id`、`listing`、`listed_date`、`source_commit`、`listed_tools`、`indexed_tools` 和 `note`。
- **示例代码：**

```python
from dnakit.integrations import mcp_tool_index_summary

summary = mcp_tool_index_summary(verbose=False)
print(summary["catalog_servers"], summary["listed_tools"], summary["indexed_tools"], summary["excluded_tools"])
for category, count in summary["indexed_by_category"].items():
    print(category, count)
for server in summary["servers"]:
    if server["server_id"] in {"addgene", "genomi", "galaxy-mcp", "tooluniverse"}:
        print(server["server_id"], server["listing"], server["listed_tools"], server["indexed_tools"])
```

- **示例输出：**

```text
49 4005 1402 2603
sequence 116
genome_annotation 360
variant 305
association 209
expression 77
regulation_epigenome 72
sequencing_data 122
prediction_model 66
molecular_tools 29
breeding_population 46
tooluniverse live-stdio+source-catalog 2927 845
genomi live-stdio+source-catalog 112 57
galaxy-mcp live-stdio 48 0
addgene live-stdio 4 4
```

## 2) `dnakit.integrations.search_mcp_tools()`

- **作用：** 按功能搜索第三方 MCP 服务器中的 DNA 与基因工具。
- **处理方法：** 在本地用 Okapi BM25 对索引中的工具排序，检索内容包括工具名、摘要、关键词、类别、上游说明和服务器名称。查询为空且给出 `server_id` 或 `category` 时，按索引顺序列出这些工具。
- **API：** `dnakit.integrations.search_mcp_tools(query, server_id[可选], category[可选], limit[可选], verbose[可选])`。
- **输入：** `query` 为英文文本（最多 1,000 个字符）。`server_id` 为目录中的 id，`category` 为十个类别之一，`limit` 为 1 到 100 的整数（默认 10）。
- **输出：** 字典，包含 `query`、`count`、`total_matches`、`indexed_tools` 和 `tools`。每个工具有 `server_id`、`tool_name`、`category`、`summary`、`required_arguments`、`listing`、`server_requires`（环境变量和安装说明）以及 `call`，`call` 说明如何用 `call_mcp_tool` 调用它。需要 API key 的 ToolUniverse 工具还带有 `tool_requires`。
- **示例代码：**

```python
from dnakit.integrations import search_mcp_tools

result = search_mcp_tools("polygenic risk score", limit=3, verbose=False)
print(result["count"], result["total_matches"], result["indexed_tools"])
for tool in result["tools"]:
    print(tool["server_id"], tool["tool_name"], tool["category"], tool["required_arguments"])
```

- **示例输出：**

```text
3 170 1402
locus polygenic_risk variant []
genechat calculate_prs variant []
gwas-mcp calculate_prs variant ['sumstats_path', 'target_plink_prefix']
```

列出一个服务器保留的工具：

- **示例代码：**

```python
from dnakit.integrations import search_mcp_tools

listed = search_mcp_tools("", server_id="gnomad", limit=3, verbose=False)
print(listed["total_matches"])
for tool in listed["tools"]:
    print(tool["tool_name"], "|", tool["summary"])
```

- **示例输出：**

```text
11
get_clinvar_variant_info | Retrieves ClinVar information for a variant ID from a gnomAD dataset (v2/v3/v4) on GRCh37 or GRCh38.
get_copy_number_variant_info | Retrieves copy number variant (CNV) information by CNV ID from gnomAD v4 (GRCh38).
get_gene_info | Retrieves gnomAD v4 gene information by Ensembl gene ID or gene symbol.
```

## 3) `dnakit.integrations.describe_mcp_tool()`

- **作用：** 调用前查看某个已索引工具的完整说明、输入 schema 和服务器信息。
- **处理方法：** 在索引中查找该工具。列出过但被排除的工具会报错，错误中给出排除类别和原因。
- **API：** `dnakit.integrations.describe_mcp_tool(server_id, tool_name, verbose[可选])`。
- **输入：** `search_mcp_tools` 返回的 `server_id` 和 `tool_name`。
- **输出：** 搜索结果中的字段，加上 `description`、`keywords`、`input_schema`、`annotations` 和 `server`（名称、仓库、许可证、固定提交、命令、参数、传输方式、安装说明、环境变量和对外影响）。
- **示例代码：**

```python
from dnakit.integrations import describe_mcp_tool

tool = describe_mcp_tool("tooluniverse", "DNA_reverse_complement", verbose=False)
print(tool["category"], tool["listing"])
print(tool["call"])
print(sorted(tool["input_schema"]["properties"]))
print(tool["server"]["license"], tool["server"]["command"])
```

- **示例输出：**

```text
sequence source-catalog
{'function': 'dnakit.integrations.call_mcp_tool', 'server_id': 'tooluniverse', 'tool_name': 'execute_tool', 'arguments': {'tool_name': 'DNA_reverse_complement', 'arguments': "the tool's own arguments (see input_schema)"}}
['operation', 'sequence']
Apache-2.0 tooluniverse
```

被排除的工具会说明原因：

- **示例代码：**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.integrations import describe_mcp_tool

try:
    describe_mcp_tool("galaxy-mcp", "run_tool", verbose=False)
except ConfigurationError as error:
    print(error.code, error.context["excluded_as"])
    print(error.context["reason"])
```

- **示例输出：**

```text
MCP_INDEX_TOOL_NOT_INDEXED generic_platform
Galaxy tool discovery/execution is generic platform management, not a DNA function itself.
```

## 4) 调用搜到的工具

每条结果的 `call` 给出 `call_mcp_tool` 的参数。多数服务器直接用工具自己的名字调用。ToolUniverse 和 genomi 则调用分发工具（`execute_tool` 或 `genomi.invoke`），把工具名和它的参数嵌套在里面。安装服务器后，以 `bioinformatics-mcp` 为例：

```python
from dnakit.integrations import call_mcp_tool, search_mcp_tools

hit = search_mcp_tools("reverse complement", server_id="bioinformatics-mcp", limit=1, verbose=False)["tools"][0]
result = call_mcp_tool(
    hit["server_id"],
    hit["call"]["tool_name"],
    {"sequence": "ATGGCC"},
    executable_path="/path/to/venv/bin/bioinformatics-mcp",
    verbose=False,
)
print(result["status"], result["content"][0]["text"])
```

2026-10-16 实际运行了三种调用方式：直接调用（`bioinformatics-mcp` 的 `reverse_complement`）、genomi 操作（通过 `genomi.invoke` 调用 `sequence.translate`）和 ToolUniverse 工具（通过 `execute_tool` 调用 `DNA_reverse_complement`），三者都对 `ATGGCC` 返回了反向互补序列或翻译结果。

## 5) 在 Agent 服务器（`dnakit-mcp`）中

`search_dnakit_tools` 新增参数 `mcp_limit`（默认 5，最多 50）。未设置 `category` 且查询非空时，返回结果多一个 `mcp_tools` 字段：同一查询下 `search_mcp_tools` 的前 `mcp_limit` 条结果。DNAKit 函数的 `tools` 列表不变。Agent 通过 `call_dnakit_tool` 调用 `dnakit.integrations.call_mcp_tool` 来使用搜到的工具，可以先用 `dnakit.integrations.describe_mcp_tool` 查看 schema。设置 `mcp_limit=0` 可不返回 MCP 工具。无论 `mode` 是什么，MCP 结果都用 BM25 排序。

`integration_inventory()` 也新增了 `mcp_listed_servers`、`mcp_listed_tools` 和 `mcp_indexed_dna_gene_tools`。

## 6) 各服务器数量

"列出"是从服务器读到的全部工具，"保留"是索引中的 DNA/基因子集。

| 服务器 | 列出方式 | 列出 | 保留 |
| --- | --- | ---: | ---: |
| `tooluniverse` | 实测 + 注册表 | 2,927 | 845 |
| `genomi` | 实测 + 注册表 | 112 | 57 |
| `plant-genomics` | 实测（stdio） | 50 | 40 |
| `bioinformatics-mcp` | 实测（stdio） | 70 | 33 |
| `gwas-mcp` | 实测（stdio） | 44 | 33 |
| `opentargets` | 实测（stdio） | 68 | 30 |
| `seqeron-sequence` | 源码 | 35 | 30 |
| `gigwa-mcp` | 实测（stdio） | 35 | 27 |
| `alphagenome-jolab` | 实测（stdio） | 24 | 23 |
| `locus` | 实测（stdio） | 22 | 21 |
| `rewire-genomics` | 实测（stdio） | 24 | 20 |
| `chorus` | 实测（stdio） | 24 | 17 |
| `kiro-genomics` | 实测（stdio） | 18 | 17 |
| `just-prs` | 实测（stdio） | 20 | 16 |
| `gget` | 实测（stdio） | 20 | 15 |
| `grebi` | 源码 | 33 | 15 |
| `mavedb-link` | 实测（HTTP） | 15 | 13 |
| `dnallm` | 实测（stdio） | 16 | 11 |
| `gnomad` | 实测（stdio） | 12 | 11 |
| `ucsc-genome` | 实测（stdio） | 12 | 11 |
| `biothings` | 实测（stdio） | 19 | 10 |
| `genechat` | 实测（stdio） | 10 | 10 |
| `encode-toolkit` | 实测（stdio） | 20 | 9 |
| `biocontext-kb` | 实测（HTTP） | 49 | 7 |
| `ensembl` | 实测（HTTP） | 7 | 7 |
| `uniprot-provenance` | 实测（stdio） | 41 | 7 |
| `brapi` | 实测（HTTP） | 22 | 6 |
| `juicebox` | 实测（stdio） | 17 | 6 |
| `biomcp` | 实测（stdio） | 7 | 5 |
| `kiro-model-organisms` | 实测（stdio） | 5 | 5 |
| `sequence-simulation` | 实测（stdio） | 6 | 5 |
| `addgene` | 实测（stdio） | 4 | 4 |
| `bioinformatics-ncbi-uniprot` | 实测（stdio） | 10 | 4 |
| `kiro-agriculture` | 实测（stdio） | 4 | 4 |
| `ncbi-database` | 源码 | 4 | 4 |
| `ncbi-sra` | 实测（stdio） | 6 | 4 |
| `10x-genomics` | 实测（stdio） | 27 | 3 |
| `kiro-clinical` | 实测（stdio） | 10 | 3 |
| `kiro-ecology` | 实测（stdio） | 7 | 3 |
| `kiro-epigenomics` | 实测（stdio） | 3 | 3 |
| `smartapi-biothings` | 实测（stdio） | 5 | 3 |
| `just-dna-agents` | 实测（stdio） | 6 | 2 |
| `opentargets-platform` | 实测（HTTP） | 5 | 2 |
| `ebi-tools` | 实测（HTTP） | 4 | 1 |
| `galaxy-mcp` | 实测（stdio） | 48 | 0 |
| `genepattern` | 实测（stdio） | 52 | 0 |
| `kiro-pipelines` | 实测（stdio） | 8 | 0 |
| `pubmed` | 实测（HTTP） | 11 | 0 |
| `pubmed-anthropic` | 实测（HTTP） | 7 | 0 |

## 范围

- 筛选按书面规则逐个工具进行。处于边界的工具，其原因以 "borderline:" 开头；对被排除的工具调用 `describe_mcp_tool` 可看到它的类别和原因。
- schema 是固定提交下的快照。上游服务器可能增加、改名或删除工具；`list_mcp_tools(server_id)` 总是显示已安装服务器当前提供的工具。
- ToolUniverse 和 genomi 的条目来自它们的注册表，而不是 `tools/list`；其中 208 个 ToolUniverse 工具只有在设置了对应 API key（或某个 LLM key）后才会出现（见 `tool_requires`）。grebi 的查询模板工具取决于具体部署。
- 有些服务器需要 `mcp<2` 才能启动（`genechat`、`gwas-mcp`、`encode-toolkit`、`ucsc-genome`、`bioinformatics-ncbi-uniprot`），目录中的 `installation` 说明已写明。
- `others/tests/unit/test_mcp_tool_index.py` 检查：每个目录服务器都有列出记录，保留和排除的工具数加起来等于列出数，搜索过滤条件，分发调用的格式，以及 Agent 服务器返回的 `mcp_tools` 字段。
