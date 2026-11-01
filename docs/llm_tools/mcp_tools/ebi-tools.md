# EBI Tools MCP (QuentinCody) (`ebi-tools`)

DNAKit 的 MCP 工具索引收录了这个服务器的 **1** 个 DNA/基因相关工具（服务器共列出 4 个，2026-10-16 读取）。工具在用户自己安装服务器后运行，DNAKit 本身不包含它们。调用方式见 [调用索引中的工具](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found)。下面的工具说明是上游英文原文。

| 项目 | 内容 |
| --- | --- |
| 源码 | [QuentinCody/ebi-tools-mcp-server @ `32437e6`](https://github.com/QuentinCody/ebi-tools-mcp-server/tree/32437e6a1ea9dee56bba4a98d502606b503c2cee) |
| 许可证 | Not specified in repository snapshot; remote service only, no upstream code bundled |
| 启动命令 | — |
| 传输方式 | http |
| 安装 | Select transport='http' for the repository author's published Cloudflare Worker MCP service. No standalone Node stdio entry point is present; no mcp-remote/npx helper is required by DNAKit. |
| 副作用 | Selected tools can submit sequences to EMBL-EBI Job Dispatcher and retrieve job results through a community-hosted service. |

| 工具 | 类别 | 简介 |
| --- | --- | --- |
| [`ebi_tools_execute`](#ebi_tools_execute) | `sequence` | Runs JavaScript against the EBI Job Dispatcher API to submit, poll and fetch jobs for BLAST, FASTA, Clustal Omega, MAFFT, MUSCLE, EMBOSS, InterProScan and simi… |

## `ebi_tools_execute` {#ebi_tools_execute}

Runs JavaScript against the EBI Job Dispatcher API to submit, poll and fetch jobs for BLAST, FASTA, Clustal Omega, MAFFT, MUSCLE, EMBOSS, InterProScan and similar tools.

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`code`
- **关键词：**blast, multiple sequence alignment, clustal omega, mafft, backtranseq, ebi job dispatcher

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `code` | string | 是 | JavaScript code to execute. Use api.get/api.post for API calls. The last expression or explicit return value becomes the result. Example: return await api.get('/dataset/tissueSiteDetail') |
| `workspace` | string | 否 | Shared workspace id — stage into a cross-server workspace DO so other servers' datasets can be JOINed in one query. Omit for per-server staging. |

??? note "上游说明（英文原文）"

    Execute JavaScript code against the EBI Job Dispatcher API (18 endpoints). Code runs in a sandboxed V8 isolate with:
    - api.get(path, params) — make GET requests (path params auto-interpolated from params object)
    - api.post(path, body, params) — make POST requests
    - api.getAll(path, params, opts) — fetch EVERY page of a paged endpoint (avoids silent under-counting). opts: \{strategy:'offset'\|'page'\|'cursor', pageSize, offsetParam, limitParam, max, maxPages, itemsField\}. Returns \{items, count, pages, total_available, completeness\} (or auto-stages if large).
    - searchSpec(query) — search the API catalog
    - listCategories() — list endpoint categories
    - getEndpoint(path) — get full endpoint docs
    - console logging (log, warn, error, info) — captured output

    The last expression or return value is the result.

    API REFERENCE (EBI Job Dispatcher, 18 of 18 endpoints):
    reverse_translation:
    GET /emboss_backtranseq/parameters — List parameters for EMBOSS Backtranseq (reverse translation)
    GET /emboss_backtranseq/parameterdetails/\{parameterId\} — Get details for a specific Backtranseq parameter (parameterId)
    sequence_similarity:
    GET /ncbiblast/parameters — List parameters for NCBI BLAST
    GET /ncbiblast/parameterdetails/\{parameterId\} — Get details for a BLAST parameter (e.g. available databases) (parameterId)
    GET /fasta/parameters — List parameters for FASTA similarity search
    GET /psisearch/parameters — List parameters for PSI-Search (iterative BLAST)
    alignment:
    GET /clustalo/parameters — List parameters for Clustal Omega multiple alignment
    GET /clustalo/parameterdetails/\{parameterId\} — Get details for a Clustal Omega parameter (parameterId)
    GET /mafft/parameters — List parameters for MAFFT multiple alignment
    GET /muscle/parameters — List parameters for MUSCLE multiple alignment
    GET /tcoffee/parameters — List parameters for T-Coffee multiple alignment
    GET /emboss_water/parameters — List parameters for EMBOSS Water (local pairwise alignment)
    protein_domains:
    GET /iprscan5/parameters — List parameters for InterProScan 5 (protein domain annotation)
    GET /iprscan5/parameterdetails/\{parameterId\} — Get details for an InterProScan parameter (parameterId)
    protein_features:
    GET /phobius/parameters — List parameters for Phobius (signal peptide + TM prediction)
    GET /radar/parameters — List parameters for RADAR internal repeat detection
    results:
    GET /ncbiblast/resulttypes/\{jobId\} — List available result types for a completed BLAST job (jobId)
    GET /iprscan5/resulttypes/\{jobId\} — List available result types for a completed InterProScan job (jobId)

    STAGING: Large responses (&gt;30KB) are auto-staged into SQLite. When this happens, api.get/api.post returns \{__staged: true, data_access_id, schema, tables_created, total_rows, message\}. Scalar properties from the original response (.count, .total, .meta) are preserved on the staged object.

    When staging occurs:
    1. Check result.__staged === true
    2. Read any preserved scalars (result.count, result.total, etc.)
    3. Return the staging metadata — the caller will use ebi_tools_query_data with the data_access_id to explore the data with SQL

    DO NOT try to access .results, .data, .entries, .items on a staged response — those arrays were replaced by SQLite tables.

    For advanced use: api.query(data_access_id, sql) and db.queryStaged(data_access_id, sql) are available to query staged data within the same execution (returns \{results, row_count\}, max 1000 rows, SELECT only). This is useful when you need to aggregate or filter staged data before returning.

    SCRATCHPAD: db.stage(data, tableName?) stages any array/object into SQLite and returns \{data_access_id, tables_created, total_rows\}. Use this to persist computed or filtered results for SQL queries without re-entering the context window. Example: const filtered = await api.query(id, "SELECT * WHERE score &gt; 0.8"); const saved = await db.stage(filtered.results, "top_hits");

    TYPED STAGING: db.stage(da \[truncated\]

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'ebi-tools',
    'ebi_tools_execute',
    {'code': '...'},
)
```
