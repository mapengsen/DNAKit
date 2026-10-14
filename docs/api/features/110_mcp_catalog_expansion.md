# MCP 目录：新增 17 个第三方服务器

适用读者：想接入更多 DNA 与基因 MCP 服务器的 DNAKit 用户。

`list_mcp_servers()` 现在返回 49 个条目，比之前多 17 个。DNAKit 不内置这些服务器。每个条目记录如何启动一个单独安装的服务器：命令、环境变量、许可证、固定的上游提交，以及工具对外部世界会做什么。收录只说明上游给出了启动方式，不代表服务正常运行，也不代表结果经过科学验证。

## 核实方式

- 每个上游仓库都克隆后，在一个固定提交上阅读。许可证以 LICENSE 文件为准；README 或 `pyproject.toml` 与之不一致的候选没有收录。
- 工具名、环境变量和启动命令来自源码，而不只是 README。
- 12 个服务器在干净环境里通过 stdio 实际启动并查询了工具：`genomi`、`galaxy-mcp`、`gnomad`、`addgene`、`locus`、`ncbi-sra`、`bioinformatics-mcp`、`rewire-genomics`、`bioinformatics-ncbi-uniprot`、`uniprot-provenance`、`smartapi-biothings` 和 `just-prs`。`tool_names` 里的每个名字都在列出的工具中。另外 5 个需要 .NET、Node.js 18、Python 3.14 或 torch，本机没有，所以它们的 `startup_check` 字段写明没有启动过。
- 不保存任何凭据：机密变量列在 `secret_environment`，只有你设置后才会传入。

## 1) `dnakit.integrations.list_mcp_servers()`

- **作用：** 列出 DNAKit 可以连接的 MCP 服务器，以及安装和启动每一个所需的信息。
- **处理方法：** 读取本地目录文件，不启动任何服务器，也不读取凭据。新条目增加了字段 `tool_names`（工具名已对照源码核实）、`startup_check`（2026-10-14 实际运行了什么），必要时还有 `limitations`。
- **API：** `dnakit.integrations.list_mcp_servers(verbose[可选])`。
- **输入：** 除 `verbose` 外无。
- **输出：** 字典列表。每项包含 `id`、`name`、`command`、`arguments`、`repository`、`sources`、`license`、`attribution`、`installation`、环境变量相关字段、`required_files`、`capabilities`、`effects`、`source_commit`、`verified_date` 和 `transports`。
- **示例代码：**

```python
from dnakit.integrations import list_mcp_servers

servers = list_mcp_servers(verbose=False)
print(len(servers))
new = {"genomi", "galaxy-mcp", "gnomad", "addgene", "locus"}
for server in servers:
    if server["id"] in new:
        print(server["id"], server["license"], server["command"], len(server["tool_names"]))
```

- **示例输出：**

```text
49
genomi Apache-2.0 genomi 8
galaxy-mcp MIT galaxy-mcp 12
gnomad Apache-2.0 python 12
addgene MIT addgene-mcp 4
locus MIT locus-mcp 22
```

再看另一个条目，了解凭据和使用限制是怎样记录的。

- **示例代码：**

```python
from dnakit.integrations import list_mcp_servers

server = next(s for s in list_mcp_servers(verbose=False) if s["id"] == "alphagenome-jolab")
print(server["secret_environment"])
print(server["startup_check"])
print(server["effects"])
```

- **示例输出：**

```text
['ALPHAGENOME_API_KEY']
2026-10-14: not started here: needs Node.js >=18 (16 installed) and an AlphaGenome API key; tool names were checked against the source.
Sends variant coordinates and alleles through a Python subprocess to the Google DeepMind AlphaGenome API (Atlas lookups or live inference) using the caller's API key; no local files are written by the tools. The upstream README states the API is free for non-commercial use only, so commercial use is restricted and requires accepting DeepMind's AlphaGenome Terms of Use and rate limits.
```

## 2) 新增的 17 个条目

| ID | 上游 | 许可证 | 启动时列出的工具数 | 说明 |
| --- | --- | --- | ---: | --- |
| `seqeron-sequence` | AlexanderV/Seqeron | MIT | — | 源码运行，需要 .NET 10 SDK；该仓库另有十个服务器，没有登记。 |
| `smartapi-biothings` | biothings/smartapi-mcp | Apache-2.0 | 5 | 以 `--facade on` 启动；设置了 SOCKS 代理时需安装 `socksio`。 |
| `just-prs` | dna-seq/just-prs-mcp | MIT | 20 | 需要 Python 3.13 或更高；扩展模式的上传工具无法运行，因为不会传入令牌。 |
| `dnallm` | zhangtaolab/DNALLM | MIT | — | 需要 torch；安装包不含配置文件，必须设置 `workdir`。 |
| `just-dna-agents` | dna-seq/just-dna-agents | MIT | — | 需要 Python 3.14 和 git 安装；README 里说已发布到 PyPI 并不属实。 |
| `genomi` | exon-research/genomi | Apache-2.0 | 37 | 未发布到 PyPI，需从源码安装；各项能力通过 `genomi.invoke` 调用。 |
| `galaxy-mcp` | galaxyproject/galaxy-mcp | MIT | 48 | 需要 `GALAXY_URL` 和 `GALAXY_API_KEY`；工具偏通用生信，PyPI 版本落后于固定提交。 |
| `gnomad` | koido/gnomad-mcp | Apache-2.0 | 12 | 源码运行，需要 Python 3.13 或更高以及 `workdir`。 |
| `addgene` | longevity-genie/addgene-mcp | MIT | 4 | 通过 Scrapy 子进程抓取 addgene.org；只接受 GenBank 和 SnapGene 格式。 |
| `alphagenome-jolab` | taehojo/alphagenome-mcp | MIT | — | 需要 Node.js 18 或更高以及 AlphaGenome API 密钥；该 API 仅限非商业免费使用。 |
| `locus` | mattsmcnulty/locus | MIT | 22 | 需要自建本地 DuckDB 库并设置 `LOCUS_DATA_DIR`；它下载的 AlphaMissense 数据是 CC-BY-NC 4.0（非商业）。 |
| `sequence-simulation` | coreymhudson/mcp-sequence-simulation | BSD-2-Clause | — | 需要 Node.js 并先构建（`npm install && npm run build`）。 |
| `ncbi-sra` | CSI-Genomics-and-Data-Analytics-Core/ncbi-sra-mcp-server | MIT | 6 | 用 `python -c` 启动，避开脚本向标准输出打印的横幅。 |
| `bioinformatics-mcp` | raonyguimaraes/bioinformatics-mcp | MIT | 70 | 未发布到 PyPI，需从 git 安装。 |
| `rewire-genomics` | rewire-bio/genomics-mcp | MIT | 24 | 仅支持 Python 3.12；`pyBigWig` 需 `--no-binary` 和 libcurl。 |
| `bioinformatics-ncbi-uniprot` | AiAgentKarl/bioinformatics-mcp-server | MIT | 10 | 安装时需固定 `mcp<2`；在 mcp 2.x 下服务器无法导入。 |
| `uniprot-provenance` | smaniches/uniprot-mcp | Apache-2.0 | 41 | 以蛋白为主；41 个工具中有 15 个写入 `tool_names`。 |

安装服务器后，`list_mcp_tools(server_id)` 会连接并返回它的工具模式，`call_mcp_tool(server_id, tool_name, arguments)` 调用其中一个工具。源码运行的条目需要 `workdir`；以 `dotnet`、`node` 或 `python` 启动的条目也可以传入 `executable_path`。

## 3) 没有收录的候选

| 候选 | 原因 |
| --- | --- |
| BACH-AI-Tools/biotools-mcp-server | BLAST、结构、比对和系统发育工具返回写死或随机的结果，README 却把它们写成真实分析。 |
| dna-seq/ensembl-mcp | LICENSE 文件是 Apache-2.0，README 徽章写 MIT，pyproject 未声明。其中 10 个 Skills 也因此暂不收录。 |
| longevity-genie/alphagenome-mcp | LICENSE 文件是 MIT，README 写 Apache-2.0，pyproject 未声明。 |
| longevity-genie/benchling-mcp | 每个工具都需要 Benchling 商业租户，没有公开的免费使用方式。 |
| sian-davies/mcp-ncbi-blast | Gradio 应用，只在本地 HTTP 端口提供 MCP；没有 stdio 入口，也没有核实过的 https 地址。 |
| bio-mcp/bio-mcp-blast | 没有 LICENSE 文件，只有 pyproject 写了 MIT。 |
| jvogan/motif | 两个工具只打开 HTML 工作台，不返回序列分析结果；npm 包是私有的。 |
| cafferychen777/ChatSpatial | 空间转录组分析，除一个拷贝数步骤外没有 DNA 或基因层面的工具。 |
| fastmcp-me/uniprot-mcp | josefdc/Uniprot-MCP 的镜像，只有 5 个蛋白工具。 |
| AWS Power Life Sciences MCP | 同一个仓库，已作为 `kiro-*` 条目收录。 |

## 当前范围

- 目录只给出启动方式。DNAKit 不安装这些服务器，也不为它们的结果担保。
- 有些服务器会调用外部服务（Ensembl、NCBI、GWAS Catalog、Addgene、AlphaGenome）。调用工具前请先读该条目的 `effects`。
- `others/tests/unit/test_mcp_catalog_expansion.py` 检查固定提交、启动命令、机密变量和未收录名单。
