# MCP catalog: 17 more third-party servers

Written for: DNAKit users who want to connect to more DNA and gene MCP servers.

`list_mcp_servers()` now returns 49 entries, 17 more than before. DNAKit does not bundle these servers. Each entry records how to start a separately installed server: the command, the environment variables, the licence, the pinned upstream commit and what the tools do to the outside world. Catalog inclusion confirms an upstream start-up recipe, not that the service is healthy or scientifically validated.

## How entries were checked

- Each upstream repository was cloned and read at one pinned commit. The licence is the one in the LICENSE file; when the README or `pyproject.toml` disagreed, the candidate was left out.
- Tool names, environment variables and the launch command come from the source code, not only the README.
- 12 servers were started over stdio in a clean environment and queried for their tools: `genomi`, `galaxy-mcp`, `gnomad`, `addgene`, `locus`, `ncbi-sra`, `bioinformatics-mcp`, `rewire-genomics`, `bioinformatics-ncbi-uniprot`, `uniprot-provenance`, `smartapi-biothings` and `just-prs`. Every name in `tool_names` was among the listed tools. The other five need .NET, Node.js 18, Python 3.14 or a torch build that this machine lacks, so their `startup_check` field says they were not started.
- Credentials are never stored: secret variables are listed in `secret_environment` and are passed only when you set them.

## 1) `dnakit.integrations.list_mcp_servers()`

- **Purpose：** List the MCP servers DNAKit can connect to, with the facts needed to install and start each one.
- **Method：** Reads the local catalog file; it starts no server and reads no credentials. The new entries add the fields `tool_names` (names checked against the source), `startup_check` (what was run on 2026-10-14) and, where relevant, `limitations`.
- **API：** `dnakit.integrations.list_mcp_servers(verbose[optional])`.
- **Input：** None besides `verbose`.
- **Output：** A list of dicts. Each has `id`, `name`, `command`, `arguments`, `repository`, `sources`, `license`, `attribution`, `installation`, the environment fields, `required_files`, `capabilities`, `effects`, `source_commit`, `verified_date` and `transports`.
- **Example code：**

```python
from dnakit.integrations import list_mcp_servers

servers = list_mcp_servers(verbose=False)
print(len(servers))
new = {"genomi", "galaxy-mcp", "gnomad", "addgene", "locus"}
for server in servers:
    if server["id"] in new:
        print(server["id"], server["license"], server["command"], len(server["tool_names"]))
```

- **Example output：**

```text
49
genomi Apache-2.0 genomi 8
galaxy-mcp MIT galaxy-mcp 12
gnomad Apache-2.0 python 12
addgene MIT addgene-mcp 4
locus MIT locus-mcp 22
```

Another entry shows how credentials and restrictions are recorded.

- **Example code：**

```python
from dnakit.integrations import list_mcp_servers

server = next(s for s in list_mcp_servers(verbose=False) if s["id"] == "alphagenome-jolab")
print(server["secret_environment"])
print(server["startup_check"])
print(server["effects"])
```

- **Example output：**

```text
['ALPHAGENOME_API_KEY']
2026-10-14: not started here: needs Node.js >=18 (16 installed) and an AlphaGenome API key; tool names were checked against the source.
Sends variant coordinates and alleles through a Python subprocess to the Google DeepMind AlphaGenome API (Atlas lookups or live inference) using the caller's API key; no local files are written by the tools. The upstream README states the API is free for non-commercial use only, so commercial use is restricted and requires accepting DeepMind's AlphaGenome Terms of Use and rate limits.
```

## 2) The 17 new entries

| ID | Upstream | Licence | Tools listed at start | Notes |
| --- | --- | --- | ---: | --- |
| `seqeron-sequence` | AlexanderV/Seqeron | MIT | — | Source-based; needs the .NET 10 SDK. The repository has ten more servers that are not catalogued. |
| `smartapi-biothings` | biothings/smartapi-mcp | Apache-2.0 | 5 | Started with `--facade on`; with a SOCKS proxy set, install `socksio`. |
| `just-prs` | dna-seq/just-prs-mcp | MIT | 20 | Needs Python 3.13 or newer. The upload tool of the extended mode cannot run, because no token is passed. |
| `dnallm` | zhangtaolab/DNALLM | MIT | — | Needs a torch build; the wheel lacks its config file, so set `workdir`. |
| `just-dna-agents` | dna-seq/just-dna-agents | MIT | — | Needs Python 3.14 and a git install; the README's PyPI claim is wrong. |
| `genomi` | exon-research/genomi | Apache-2.0 | 37 | Not on PyPI; install from source. Capability operations go through `genomi.invoke`. |
| `galaxy-mcp` | galaxyproject/galaxy-mcp | MIT | 48 | Needs `GALAXY_URL` and `GALAXY_API_KEY`. Tools are general bioinformatics; PyPI lags the pinned commit. |
| `gnomad` | koido/gnomad-mcp | Apache-2.0 | 12 | Source-based; needs Python 3.13 or newer and `workdir`. |
| `addgene` | longevity-genie/addgene-mcp | MIT | 4 | Scrapes addgene.org through a Scrapy subprocess; accepts GenBank and SnapGene formats only. |
| `alphagenome-jolab` | taehojo/alphagenome-mcp | MIT | — | Needs Node.js 18 or newer and an AlphaGenome API key; the API is free for non-commercial use only. |
| `locus` | mattsmcnulty/locus | MIT | 22 | Needs a local DuckDB store you build and `LOCUS_DATA_DIR`; the AlphaMissense data it downloads is CC-BY-NC 4.0. |
| `sequence-simulation` | coreymhudson/mcp-sequence-simulation | BSD-2-Clause | — | Needs Node.js and a build (`npm install && npm run build`). |
| `ncbi-sra` | CSI-Genomics-and-Data-Analytics-Core/ncbi-sra-mcp-server | MIT | 6 | Started with `python -c` to avoid the banner the script prints to stdout. |
| `bioinformatics-mcp` | raonyguimaraes/bioinformatics-mcp | MIT | 70 | Not on PyPI; install from git. |
| `rewire-genomics` | rewire-bio/genomics-mcp | MIT | 24 | Python 3.12 only; `pyBigWig` needs `--no-binary` and libcurl. |
| `bioinformatics-ncbi-uniprot` | AiAgentKarl/bioinformatics-mcp-server | MIT | 10 | Install with `mcp<2`; under mcp 2.x the server fails to import. |
| `uniprot-provenance` | smaniches/uniprot-mcp | Apache-2.0 | 41 | Protein-centred; 15 of its 41 tools are listed in `tool_names`. |

After installing a server, `list_mcp_tools(server_id)` connects and returns its schemas, and `call_mcp_tool(server_id, tool_name, arguments)` calls one tool. Source-based entries need `workdir`; entries that start with `dotnet`, `node` or `python` can also take `executable_path`.

## 3) Candidates left out

| Candidate | Reason |
| --- | --- |
| BACH-AI-Tools/biotools-mcp-server | BLAST, structure, alignment and phylogeny tools return hard-coded or random results while the README presents them as real analysis. |
| dna-seq/ensembl-mcp | The LICENSE file is Apache-2.0, the README badge says MIT and pyproject declares none. Its 10 Skills are held back for the same reason. |
| longevity-genie/alphagenome-mcp | The LICENSE file is MIT, the README says Apache-2.0 and pyproject declares none. |
| longevity-genie/benchling-mcp | Every tool needs a commercial Benchling tenant; no free access is documented. |
| sian-davies/mcp-ncbi-blast | A Gradio app that serves MCP only on a local HTTP port; no stdio entry point and no verified https endpoint. |
| bio-mcp/bio-mcp-blast | No LICENSE file; only pyproject says MIT. |
| jvogan/motif | Its two tools only open an HTML workbench and return no sequence results; the npm package is private. |
| cafferychen777/ChatSpatial | Spatial transcriptomics analysis; no DNA or gene-level tools beyond a copy-number step. |
| fastmcp-me/uniprot-mcp | A mirror of josefdc/Uniprot-MCP with five protein-only tools. |
| AWS Power Life Sciences MCP | Already in the catalog as the `kiro-*` entries (same repository). |

## Scope

- The catalog gives start-up recipes only. DNAKit does not install the servers and does not vouch for their results.
- Several servers call outside services (Ensembl, NCBI, GWAS Catalog, Addgene, AlphaGenome). Read each entry's `effects` before calling a tool.
- `others/tests/unit/test_mcp_catalog_expansion.py` checks the pinned commits, launch commands, secrets and the left-out list.
