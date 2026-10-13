# 服务后端与本地 Toil 作业

这一组函数位于 `dnakit.integrations`，是基于可选依赖的服务客户端和本地作业适配器，用于基因、变异、序列和 API 查询。[Toil](https://github.com/DataBiosphere/toil) 在本地运行两个串联的作业，计算反向互补序列和 GC 含量；[smartapi-mcp](https://github.com/biothings/smartapi-mcp) 搜索 SmartAPI 注册表并读取 OpenAPI 规范；[biomcp](https://github.com/genomoncology/biomcp) 通过命令行查询基因和变异；[galaxy-mcp](https://github.com/galaxyproject/galaxy-mcp) 与 [bioblend](https://github.com/galaxyproject/bioblend) 读取 Galaxy 工具列表；[clawbio](https://github.com/ClawBio/ClawBio) 的 `analyze-fasta` 技能统计一个 FASTA 文件。

这一组不做序列比对、变异检测、参考数据下载，也不运行 Galaxy 工具或工作流。网络函数会访问公开服务，结果随上游数据库变化；本地函数只处理调用时显式给出的序列或文件。所有依赖都是可选的，只有调用对应函数时才会导入。

## 1) `dnakit.integrations.compute_toil_sequence_stats()`

- **作用：** 对一组 DNA 序列计算反向互补序列和 GC 含量。
- **处理方法：** 在本地 Toil 作业图中运行两个串联作业。第一个作业用 Biopython 的 `reverse_complement` 生成反向互补序列；第二个作业接收第一个作业的输出，用 `gc_fraction(ambiguous="ignore")` 计算 GC 比例，即 G 与 C 的个数除以全部字符数（包括 N 等模糊碱基）。输入先转为大写。Toil 使用 `single_machine` 批处理系统，临时目录在函数返回时删除。
- **API：** `dnakit.integrations.compute_toil_sequence_stats(sequences, *, max_sequences=100, max_total_bases=1000000)`。
- **输入：** `sequences` 是 `{名称: 序列}` 映射或 `DNA` 对象，只接受 IUPAC 核苷酸字符。`max_sequences` 最大为 1000，`max_total_bases` 最大为 100000000，超过时报错。
- **输出：** 含 `sequence_count`、`sequences`（每项含 `id`、`length`、`reverse_complement`、`gc_fraction`）和 `provenance` 的字典，顺序与输入一致。
- **示例代码：**

```python
from dnakit.integrations import compute_toil_sequence_stats

result = compute_toil_sequence_stats({"first": "ATGCGCAATT", "second": "GGGCCCAAAT"})
for row in result["sequences"]:
    print(row["id"], row["reverse_complement"], row["gc_fraction"])
```

- **示例输出：**

```text
first AATTGCGCAT 0.4
second ATTTGGGCCC 0.6
```

## 2) `dnakit.integrations.search_smartapi_registry()`

- **作用：** 按关键词搜索 SmartAPI 注册表中的 API，并按相关性排序。
- **处理方法：** 先把关键词转为小写再请求注册表，因为注册表的相关性区分大小写，"MyGene" 与 "mygene" 返回的结果不同。随后用 smartapi-mcp 的 `rank_apis` 按关键词与标题、描述、标签的匹配程度排序，并去掉得分为 0 的 API。
- **API：** `dnakit.integrations.search_smartapi_registry(query, *, limit=10)`。
- **输入：** 非空查询字符串，最长 200 个字符；`limit` 取 1 到 100。
- **输出：** 含 `query`、`registry_query`、`registry_hits`（注册表返回的条目数）、`matches`（每项含 `smartapi_id`、`name`、`title`、`description`、`tags`、`score`）和 `provenance` 的字典。
- **示例代码：**

```python
from dnakit.integrations import search_smartapi_registry

result = search_smartapi_registry("mygene", limit=3)
for match in result["matches"]:
    print(match["smartapi_id"], match["title"])
```

- **示例输出：**

```text
59dce17363dce279d389100834e43648 MyGene.info API
a9baab984f6af0021aaece6a163e0cf4 MyGene.info API (for test)
```

## 3) `dnakit.integrations.describe_smartapi_api()`

- **作用：** 读取一个 SmartAPI 登记的 OpenAPI 规范，并列出其中的 HTTP 操作。
- **处理方法：** 用 smartapi-mcp 的 `load_api_spec` 下载并校验规范（拒绝外部 `$ref`，限制规范大小），用 `get_base_server_url` 取得服务地址，再按路径排序、按 GET、PUT、POST、DELETE 等顺序列出操作。
- **API：** `dnakit.integrations.describe_smartapi_api(smartapi_id, *, max_operations=200)`。
- **输入：** 32 位小写十六进制的 SmartAPI ID；`max_operations` 取 1 到 1000。
- **输出：** 含 `title`、`version`、`openapi`、`server_url`、`operation_count`、`truncated`、`operations`（每项含 `method`、`path`、`operation_id`、`summary`）和 `provenance` 的字典。
- **示例代码：**

```python
from dnakit.integrations import describe_smartapi_api

spec = describe_smartapi_api("59dce17363dce279d389100834e43648")
print(spec["title"], spec["server_url"], spec["operation_count"])
for operation in spec["operations"][:2]:
    print(operation["method"], operation["path"])
```

- **示例输出：**

```text
MyGene.info API https://mygene.info/v3 6
POST /gene
GET /gene/{id}
```

## 4) `dnakit.integrations.get_biomcp_gene()`

- **作用：** 按人类基因符号查询基因的标识符、位置和摘要。
- **处理方法：** 调用 `biomcp -j get gene <symbol>`。命令不经过 shell，参数固定，输出限制为 1 MB，超时 60 秒。biomcp 没有 Python API，因此这里使用它的 JSON 命令行输出。可执行文件先从 `PATH` 查找，找不到时查找与当前 Python 同目录的 `biomcp`。
- **API：** `dnakit.integrations.get_biomcp_gene(symbol)`。
- **输入：** 基因符号，例如 `BRCA1`、`HLA-A`、`C9orf72`，最长 32 个字符。
- **输出：** 含 `symbol`、`name`、`entrez_id`、`ensembl_id`、`uniprot_id`、`omim_id`、`gene_type`、`location`、`coordinate`、`genome_build`、`aliases`、`clinical_diseases`、`summary` 和 `provenance` 的字典。坐标和摘要由 biomcp 原样返回，不在本地重新计算。
- **示例代码：**

```python
from dnakit.integrations import get_biomcp_gene

gene = get_biomcp_gene("BRCA1")
print(gene["symbol"], gene["entrez_id"], gene["location"], gene["genome_build"])
```

- **示例输出：**

```text
BRCA1 672 17 GRCh38
```

## 5) `dnakit.integrations.search_biomcp_genes()`

- **作用：** 按文本搜索基因记录。
- **处理方法：** 调用 `biomcp -j search gene -q <query> --limit <n>`，并解析返回的 `results`。
- **API：** `dnakit.integrations.search_biomcp_genes(query, *, limit=10)`。
- **输入：** 非空查询字符串，最长 200 个字符，且不能以 `-` 开头；`limit` 取 1 到 50。
- **输出：** 含 `query`、`total`、`genes`（每项含 `symbol`、`name`、`entrez_id`、`uniprot_id`、`omim_id`、`coordinate`）和 `provenance` 的字典。
- **示例代码：**

```python
from dnakit.integrations import search_biomcp_genes

result = search_biomcp_genes("BRCA1", limit=2)
for gene in result["genes"]:
    print(gene["symbol"], gene["entrez_id"], gene["uniprot_id"])
```

- **示例输出：**

```text
BRCA1 672 P38398
```

## 6) `dnakit.integrations.get_biomcp_variant()`

- **作用：** 按 rsID 或基因组 HGVS 写法查询变异注释。
- **处理方法：** 调用 `biomcp -j get variant <id>`。`cadd_score` 不是有限数字时返回 `None`。
- **API：** `dnakit.integrations.get_biomcp_variant(variant_id)`。
- **输入：** rsID（如 `rs334`）或基因组 HGVS（如 `chr7:g.140753336A>T`），最长 128 个字符，不能包含空白，也不能以 `-` 开头。
- **输出：** 含 `id`、`rsid`、`gene`、`hgvs_c`、`hgvs_p`、`consequence`、`significance`、`transcript`、`genome_build`、`cadd_score`、`sift_pred`、`polyphen_pred` 和 `provenance` 的字典。
- **示例代码：**

```python
from dnakit.integrations import get_biomcp_variant

variant = get_biomcp_variant("rs334")
print(variant["gene"], variant["hgvs_p"], variant["significance"])
```

- **示例输出：**

```text
HBB p.Glu7Val Pathogenic
```

## 7) `dnakit.integrations.list_galaxy_tools()`

- **作用：** 列出一个 Galaxy 服务器上的工具。
- **处理方法：** 使用 galaxy-mcp 的 `GalaxyInstance`，它禁止跟随重定向，并拒绝带用户名、密码、查询串或片段的地址；请求通过 bioblend 的 `get_tools()` 完成。API key 放在 `x-api-key` 请求头中，不会出现在输出里。工具按 `id` 排序。bioblend 会先下载完整列表，`limit` 在之后截取。
- **API：** `dnakit.integrations.list_galaxy_tools(url, *, api_key=None, limit=50)`。
- **输入：** Galaxy 服务器的 HTTP(S) 地址。`api_key` 可以省略，公共服务器的部分接口允许匿名访问。`limit` 取 1 到 500。
- **输出：** 含 `url`（规范化后的地址）、`total_tools`、`returned`、`tools`（每项含 `id`、`name`、`version`、`description`）和 `provenance` 的字典。
- **示例代码：**

```python
from dnakit.integrations import list_galaxy_tools

result = list_galaxy_tools("https://usegalaxy.org", limit=2)
print(result["url"], result["returned"], result["total_tools"] > 1000)
```

- **示例输出：**

```text
https://usegalaxy.org/ 2 True
```

## 8) `dnakit.integrations.show_galaxy_tool()`

- **作用：** 查看一个 Galaxy 工具的基本信息。
- **处理方法：** 通过 bioblend 的 `show_tool(tool_id)` 请求工具详情，使用与上一个函数相同的重定向和凭据规则。
- **API：** `dnakit.integrations.show_galaxy_tool(url, *, tool_id, api_key=None)`。
- **输入：** 服务器地址和工具 ID，例如 `Cut1`。工具 ID 不能包含 `..`。
- **输出：** 含 `url`、`tool`（含 `id`、`name`、`version`、`description`、`model_class`）和 `provenance` 的字典。
- **示例代码：**

```python
from dnakit.integrations import show_galaxy_tool

result = show_galaxy_tool("https://usegalaxy.org", tool_id="Cut1")
print(result["tool"]["id"], result["tool"]["name"], result["tool"]["version"])
```

- **示例输出：**

```text
Cut1 Cut 1.0.2
```

## 9) `dnakit.integrations.analyze_fasta_with_clawbio()`

- **作用：** 用 clawbio 的 `analyze-fasta` 技能统计一个核酸 FASTA 文件。
- **处理方法：** 从已安装的 clawbio 包中加载 `skills/analyze-fasta/analyze_fasta.py`，调用其中的 `analyze_fasta()`，并整理每条序列的长度、GC 百分比和碱基组成。输出中去掉了分析时间戳和文件的绝对路径。如果 clawbio 报告错误，或第一条记录被判断为蛋白质，函数会报错。文件大小超过 `max_bytes` 时直接拒绝。
- **API：** `dnakit.integrations.analyze_fasta_with_clawbio(path, *, max_bytes=10000000)`。
- **输入：** 已存在的 FASTA 文件路径；`max_bytes` 取 1 到 100000000。
- **输出：** 含 `file`、`sequence_type`、`sequences`（每项含 `index`、`id`、`length_bp`、`gc_percent`、`composition`）、`summary`（含 `total_sequences`、`total_residues`、`min_length`、`max_length`、`n50`、`avg_gc_percent`）和 `provenance` 的字典。
- **示例文件：** `sequences.fa`

```text
>seq1 test
ATGCGCAA
TT
>seq2
GGGCCCAAAT
```

- **示例代码：**

```python
from dnakit.integrations import analyze_fasta_with_clawbio

result = analyze_fasta_with_clawbio("sequences.fa")
print(result["summary"]["total_sequences"], result["summary"]["n50"], result["summary"]["avg_gc_percent"])
for row in result["sequences"]:
    print(row["id"], row["length_bp"], row["gc_percent"])
```

- **示例输出：**

```text
2 10 50.0
seq1 10 40.0
seq2 10 60.0
```

## 安装可选依赖

```bash
python -m pip install "dnakit[bio-toil]"          # toil 与 biopython，compute_toil_sequence_stats
python -m pip install "dnakit[bio-smartapi]"      # smartapi-mcp，SmartAPI 搜索与规范读取
python -m pip install "dnakit[bio-biomcp]"        # biomcp-cli，biomcp 命令
python -m pip install "dnakit[bio-galaxy]"        # galaxy-mcp 与 bioblend，Galaxy 工具列表
python -m pip install "dnakit[bio-clawbio]"       # clawbio 与 biopython，analyze_fasta_with_clawbio
```

clawbio 的依赖较多，安装时间明显长于其他包。各包的项目地址：[toil](https://github.com/DataBiosphere/toil)、[smartapi-mcp](https://github.com/biothings/smartapi-mcp)、[biomcp](https://github.com/genomoncology/biomcp)、[galaxy-mcp](https://github.com/galaxyproject/galaxy-mcp)、[bioblend](https://github.com/galaxyproject/bioblend)、[clawbio](https://github.com/ClawBio/ClawBio)。

## 验证方式

测试文件是 `others/tests/unit/test_service_backends.py`，共 30 项：27 项离线测试，3 项实时测试（需要 `DNAKIT_LIVE_TESTS=1`）。主要检查：

- Toil 结果与 Biopython 的直接调用一致，并与手工计算一致。例如 ATGCGCAATT 的反向互补是 AATTGCGCAT，GC 为 0.4。
- SmartAPI 搜索使用本地模拟的 HTTP 响应。排序结果与直接调用 `rank_apis` 一致，包括得分。搜索词不匹配的 API 被排除。
- biomcp 用一个只记录参数并输出固定 JSON 的 POSIX 脚本代替，检查固定的命令参数（HGVS 字符串作为单个参数传递）、JSON 解析、错误 JSON 和非 JSON 输出。实时测试将结果与直接运行 biomcp 的结果比较。
- Galaxy 工具列表使用本机回环地址上的桩服务器测试。检查排序、截断、API key 只出现在请求头中，以及带用户名或查询串的地址被拒绝。实时测试访问 usegalaxy.org。
- clawbio 的统计结果与手工计算的长度、GC 百分比和碱基组成一致，并与直接调用脚本的结果一致。
- 缺少可选包时，每个模块都会给出包含 `bio-*` 额外名称的 `BackendUnavailableError`。biomcp 没有 Python 模块，它的缺失测试改为替换可执行文件查找函数。

## 当前范围

这一批覆盖 5 个上游包的 9 个函数：

- Toil：只封装一个两步作业图。集群、云后端和 CWL、WDL 工作流接口未封装。
- smartapi-mcp：封装注册表搜索和规范读取。服务器构建和 MCP 工具注册未封装。
- biomcp：封装基因、变异查询和基因搜索。药物、疾病、文献、临床试验和 `enrich` 等命令未封装。
- galaxy-mcp：只读取工具列表和工具详情。启动工具、上传数据和导入工作流会修改服务器状态，未封装。
- clawbio：只封装 `analyze-fasta` 技能。其他技能需要外部程序或特定数据，未封装。
- `list_galaxy_tools()` 会下载服务器上的完整工具列表，再截取 `limit` 条。这一步的内存用量没有上限，对大型服务器应注意。
