# 通用数据库接口

本类收录 `databases` 模块的 **8** 个函数。其中 5 个（`cache`、`download`、`export`、`fetch`、`list_versions`）始终返回 `status: "not_implemented"`：它们各自需要 DNAKit 未定义的指定数据库、发布版本目录或缓存布局，且每条消息都指向可用的替代方法。`query` 按名称调用一个内置数据库操作，`map_ids` 通过 MyGene.info 映射基因标识符；二者都要访问远程服务，其中 `map_ids` 在安装可选包 `mygene` 之前返回 `not_implemented`。`search` 列出匹配的内置操作。示例输出均为真实的离线运行结果；`input_summary` 简写为 `{...}`。

## 1) cache

- **作用：** 未实现：数据库适配器没有定义任何数据库结果缓存。

- **处理方法：** 不做任何计算。数据库适配器不保存缓存，因此这个函数没有可填充或检查的内容。

- **API：** `dnakit.databases.cache(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：任意值（仅在 `input_summary` 中汇总）；`options` 不起作用。

- **输出：** 始终返回 `status: "not_implemented"` 的封装，含 `required_dependency`（“具有明确布局和过期策略的数据库结果缓存”）和 `message`（数据库适配器不保存缓存，这个函数没有可填充或检查的内容），没有 `result`；正常调用不会抛出异常。

- **示例代码：**

```python
import json
from dnakit.databases import cache

result = cache("uniprot")
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "databases",
  "operation": "cache",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a database-result cache with a defined layout and expiry policy",
  "message": "The database adapters keep no cache; there is nothing for this function to populate or inspect."
}
```

## 2) download

- **作用：** 未实现：通用下载需要指定的数据库版本；请使用具体的下载或获取函数。

- **处理方法：** 不做任何计算，也不下载任何内容。各内置数据库之间没有通用下载；可用 `dnakit.databases.search` 找到具体操作，再由 `query` 调用。

- **API：** `dnakit.databases.download(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：任意值（仅在 `input_summary` 中汇总）；`options` 不起作用。

- **输出：** 始终返回 `status: "not_implemented"` 的封装，含 `required_dependency`（“指定的数据库版本和网络访问”）和 `message`（不存在通用的数据库下载；请用 dnakit.databases.search 查找具体操作），没有 `result`；正常调用不会抛出异常。

- **示例代码：**

```python
import json
from dnakit.databases import download

result = download({"database": "clinvar", "release": "2024-01"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "databases",
  "operation": "download",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a named database release and network access",
  "message": "No generic database download exists; use dnakit.databases.search to find the specific operation."
}
```

## 3) export

- **作用：** 未实现：不存在通用的数据库内容导出。

- **处理方法：** 不做任何计算，也不写任何文件。如需保存数据库调用的结果，请把 `query` 的结果交给 `io_formats` 模块的函数。

- **API：** `dnakit.databases.export(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：任意值（仅在 `input_summary` 中汇总）；`options` 不起作用。

- **输出：** 始终返回 `status: "not_implemented"` 的封装，含 `required_dependency`（“数据库快照和目标格式”）和 `message`（不存在通用的数据库导出；请用 io 函数导出查询结果），没有 `result`；正常调用不会抛出异常。

- **示例代码：**

```python
import json
from dnakit.databases import export

result = export({"database": "clinvar"}, options={"format": "tsv"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "databases",
  "operation": "export",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a database snapshot and a target format",
  "message": "No generic database export exists; export the result of a query with the io functions."
}
```

## 4) fetch

- **作用：** 未实现：`fetch` 没有指明数据库；请使用 `query` 或具体的 `fetch_*` 函数。

- **处理方法：** 不做任何计算，也不发出请求。请改用带内置操作名的 `dnakit.databases.query`，或 `fetch_ensembl_gene_info` 等具体函数。

- **API：** `dnakit.databases.fetch(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：任意值（仅在 `input_summary` 中汇总）；`options` 不起作用。

- **输出：** 始终返回 `status: "not_implemented"` 的封装，含 `required_dependency`（“指定的数据库和网络访问”）和 `message`（fetch 没有指明数据库；请改为用操作名调用 dnakit.databases.query），没有 `result`；正常调用不会抛出异常。

- **示例代码：**

```python
import json
from dnakit.databases import fetch

result = fetch("ENSG00000012048")
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "databases",
  "operation": "fetch",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a named database and network access",
  "message": "fetch does not name a database; call dnakit.databases.query with an operation name instead."
}
```

## 5) list_versions

- **作用：** 未实现：DNAKit 不保存数据库发布版本的目录。

- **处理方法：** 不做任何计算。没有任何 DNAKit 函数列出数据库发布版本；各提供方操作返回的是远程服务当前提供的版本。

- **API：** `dnakit.databases.list_versions(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：任意值（仅在 `input_summary` 中汇总）；`options` 不起作用。

- **输出：** 始终返回 `status: "not_implemented"` 的封装，含 `required_dependency`（“该数据库的发布版本目录”）和 `message`（没有任何 DNAKit 函数列出数据库发布版本），没有 `result`；正常调用不会抛出异常。

- **示例代码：**

```python
import json
from dnakit.databases import list_versions

result = list_versions("clinvar")
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "databases",
  "operation": "list_versions",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a release catalogue of the database",
  "message": "No DNAKit function lists database releases."
}
```

## 6) map_ids

- **作用：** 通过委托给 `dnakit.databases.map_gene_identifiers`，把基因标识符映射到 MyGene.info 记录。

- **处理方法：** 标识符会发送到 MyGene.info 服务，因此需要网络访问和可选包 `mygene`。缺少 `mygene` 时函数返回 `status: "not_implemented"`，`required_dependency` 为 `"pip install mygene"`，`message` 为导入错误信息（示例即此情形）。标识符不会在本地映射。

- **API：** `dnakit.databases.map_ids(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：标识符列表或以逗号分隔的字符串（空列表或其他类型会抛 `InputFormatError`）。`options`：`scopes`（标识符类型，默认 `"entrezgene"`）、`species`（默认 `"human"`）、`fields`（默认 `"symbol,name,taxid"`）。

- **输出：** 封装包含 `status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`（`provenance.implementation` 为 `dnakit.databases.map_gene_identifiers`），前提是已安装 `mygene` 且服务有响应；`result` 为 `map_gene_identifiers` 的文档，每个标识符在 `results` 下占一行。没有 `mygene` 时，封装改为 `status: "not_implemented"` 并含 `required_dependency` 和 `message`，如示例所示（示例在没有 `mygene` 的环境中运行，未发出任何请求）。

- **示例代码：**

```python
import json
from dnakit.databases import map_ids

result = map_ids("1017, 7157", options={"scopes": "entrezgene", "species": "human"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "databases",
  "operation": "map_ids",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "pip install mygene",
  "message": "[BACKEND_UNAVAILABLE] mygene could not be imported: No module named 'mygene' Hint: Install the optional dependency with: python -m pip install \"dnakit[bio-mygene]\""
}
```

## 7) query

- **作用：** 按名称通过 `dnakit.databases.call` 调用一个内置的 DNAKit 数据库操作。

- **处理方法：** 操作名必须是内置操作之一（见 `dnakit.databases.list_operations` 和 `search`）；参数按该操作的输入模式校验，然后分派给提供方实现，由其访问远程数据库，因此成功调用需要网络。未知操作会抛 `ConfigurationError`（`UNKNOWN_DATABASE_OPERATION`）；提供方失败会抛 `BackendExecutionError`（`DATABASE_PROVIDER_ERROR`）。`query` 本身不解释也不缓存返回的文档。示例必须能离线运行，因此下面没有展示成功的远程调用，而是展示未知操作的错误。

- **API：** `dnakit.databases.query(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：操作名（参数放在 `options["arguments"]`），或映射 `{"name": operation, "arguments": {...}}`（其他键会抛 `InputFormatError`）；`arguments` 必须是映射（否则抛 `ConfigurationError`）。`options`：`arguments`、`max_result_bytes`（默认 20,000,000，限制结果大小）。

- **输出：** 封装包含 `status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`（`provenance.implementation` 为 `dnakit.databases.call`）；`result` 为提供方返回的 JSON 文档，字段随操作而异。示例打印的是未知操作的错误（`ConfigurationError` 的 `code` 和 `message`），而不是封装。

- **示例代码：**

```python
import json
from dnakit.databases import query
from dnakit.exceptions import ConfigurationError

try:
    query("definitely_unknown_operation")
except ConfigurationError as exc:
    print(json.dumps({"error": type(exc).__name__, "code": exc.code, "message": exc.message}, indent=2))
    # result = query({"name": "pdb_get_structure", "arguments": {"pdb_id": "1ABC"}})  # needs network
```

- **示例输出：**

```json
{
  "error": "ConfigurationError",
  "code": "UNKNOWN_DATABASE_OPERATION",
  "message": "Unknown DNAKit database operation: definitely_unknown_operation"
}
```

## 8) search

- **作用：** 列出名称或描述与查询匹配的内置数据库和网络 API 操作。它搜索的是操作的本地目录（目前 450 项），不是数据库本身，也不发起网络请求。

- **处理方法：** `dnakit.databases.search` 就是 `dnakit.databases._runtime.search`，一个普通的关键词匹配函数，不走共享的候选运行时（`src/dnakit/databases/_candidates.py` 里通用的 `search(input, options)` 占位函数被它覆盖，不是公开的那个函数）。查询先转为小写并拆分成互不相同的词元（`\w+`）。对每个操作，把 `name + upstream_name + description` 转为小写，得分为在其中作为子串出现的查询词元个数（因此像 `a` 这样的短词元几乎匹配所有操作）。得分大于 0 的操作按得分降序、再按名称排序，返回前 `limit` 项。没有词干处理、同义词处理或相关性模型。要查看单个操作的参数用 `dnakit.databases.describe(name)`，要运行它用 `query` 或 `call`（这些可能访问远程服务）。

- **API：** `dnakit.databases.search(query[required], limit[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `query`：非空字符串。`limit`：1 到 100 的整数，默认 10（布尔值会被拒绝）。`verbose`：是否打印进度信息，默认 `True`。`query` 或 `limit` 不合法时抛出 `ConfigurationError`（代码 `CONFIGURATION_ERROR`，消息 `Provide a nonempty query and limit between 1 and 100.`）。

- **输出：** 与候选运行时的函数不同，它不返回封装（没有 `status`、`provenance`、`input_summary`），而是返回映射列表，最佳匹配在前，每项含 `name`（DNAKit 操作名）、`description` 和 `upstream_name`（原始 API 操作名）。无匹配时列表为空。

- **示例代码：**

```python
import json
from dnakit.databases import search

result = search("clinvar variant", limit=2)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
[
  {
    "name": "clin_var_get_clinical_significance",
    "description": "Get clinical significance information for a variant from ClinVar. Returns pathogenicity classification and clinical interpretations.",
    "upstream_name": "ClinVar_get_clinical_significance"
  },
  {
    "name": "clin_var_get_variant_details",
    "description": "Get variant summary information from ClinVar by variant ID. Returns accession, title, genes, clinical significance, review status, chromosome location, and variation name. For the same variant, this provides a different field view than ClinVar_get_clinical_significance.",
    "upstream_name": "ClinVar_get_variant_details"
  }
]
```
