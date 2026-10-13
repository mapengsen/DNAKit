# Service backends and a local Toil job graph

The functions in this group live in `dnakit.integrations`. They are adapters over optional packages for gene, variant, sequence, and API lookups. [Toil](https://github.com/DataBiosphere/toil) runs two chained jobs locally to compute reverse complements and GC content. [smartapi-mcp](https://github.com/biothings/smartapi-mcp) searches the SmartAPI registry and loads OpenAPI specs. [biomcp](https://github.com/genomoncology/biomcp) queries genes and variants through its command line. [galaxy-mcp](https://github.com/galaxyproject/galaxy-mcp) and [bioblend](https://github.com/galaxyproject/bioblend) read Galaxy tool lists. The `analyze-fasta` skill from [clawbio](https://github.com/ClawBio/ClawBio) summarizes one FASTA file.

This group does not perform sequence alignment, variant calling, reference downloads, or Galaxy tool and workflow runs. Network functions call public services, and their results change with the upstream data. Local functions only process the sequences or files passed to them. All dependencies are optional and are imported only when the matching function is called.

## 1) `dnakit.integrations.compute_toil_sequence_stats()`

- **Purpose:** Compute reverse complements and GC content for a set of DNA sequences.
- **Method:** Run two chained jobs in a local Toil job graph. The first job uses Biopython's `reverse_complement` to build each reverse complement. The second job, a follow-on that receives the first job's output, uses `gc_fraction(ambiguous="ignore")` to compute the GC fraction: the count of G and C divided by all characters, including ambiguous bases such as N. Input is uppercased first. Toil uses the `single_machine` batch system, and its temporary directory is removed when the function returns.
- **API:** `dnakit.integrations.compute_toil_sequence_stats(sequences, *, max_sequences=100, max_total_bases=1000000)`.
- **Input:** `sequences` is a `{name: sequence}` mapping or a `DNA` object, and accepts only IUPAC nucleotide characters. `max_sequences` is at most 1000, and `max_total_bases` is at most 100000000. Larger inputs raise an error.
- **Output:** A dictionary with `sequence_count`, `sequences` (each item has `id`, `length`, `reverse_complement`, and `gc_fraction`), and `provenance`. The order matches the input.
- **Example code:**

```python
from dnakit.integrations import compute_toil_sequence_stats

result = compute_toil_sequence_stats({"first": "ATGCGCAATT", "second": "GGGCCCAAAT"})
for row in result["sequences"]:
    print(row["id"], row["reverse_complement"], row["gc_fraction"])
```

- **Example output:**

```text
first AATTGCGCAT 0.4
second ATTTGGGCCC 0.6
```

## 2) `dnakit.integrations.search_smartapi_registry()`

- **Purpose:** Search the SmartAPI registry by keyword and rank the matching APIs.
- **Method:** The keyword is lowercased before the registry request, because the registry's relevance depends on case: "MyGene" and "mygene" return different hits. The results are then ranked with smartapi-mcp's `rank_apis`, which weights the keyword matches in each API's title, description, and tags. APIs with a score of 0 are removed.
- **API:** `dnakit.integrations.search_smartapi_registry(query, *, limit=10)`.
- **Input:** A non-empty query string of at most 200 characters. `limit` is from 1 to 100.
- **Output:** A dictionary with `query`, `registry_query`, `registry_hits` (the number of registry entries returned), `matches` (each item has `smartapi_id`, `name`, `title`, `description`, `tags`, and `score`), and `provenance`.
- **Example code:**

```python
from dnakit.integrations import search_smartapi_registry

result = search_smartapi_registry("mygene", limit=3)
for match in result["matches"]:
    print(match["smartapi_id"], match["title"])
```

- **Example output:**

```text
59dce17363dce279d389100834e43648 MyGene.info API
a9baab984f6af0021aaece6a163e0cf4 MyGene.info API (for test)
```

## 3) `dnakit.integrations.describe_smartapi_api()`

- **Purpose:** Load one SmartAPI-registered OpenAPI spec and list its HTTP operations.
- **Method:** smartapi-mcp's `load_api_spec` downloads and validates the spec. It rejects external `$ref` values and limits the spec size. `get_base_server_url` gets the service address. Operations are sorted by path, then by method in the order GET, PUT, POST, DELETE.
- **API:** `dnakit.integrations.describe_smartapi_api(smartapi_id, *, max_operations=200)`.
- **Input:** A 32-character lowercase hexadecimal SmartAPI ID. `max_operations` is from 1 to 1000.
- **Output:** A dictionary with `title`, `version`, `openapi`, `server_url`, `operation_count`, `truncated`, `operations` (each item has `method`, `path`, `operation_id`, and `summary`), and `provenance`.
- **Example code:**

```python
from dnakit.integrations import describe_smartapi_api

spec = describe_smartapi_api("59dce17363dce279d389100834e43648")
print(spec["title"], spec["server_url"], spec["operation_count"])
for operation in spec["operations"][:2]:
    print(operation["method"], operation["path"])
```

- **Example output:**

```text
MyGene.info API https://mygene.info/v3 6
POST /gene
GET /gene/{id}
```

## 4) `dnakit.integrations.get_biomcp_gene()`

- **Purpose:** Look up a human gene symbol for its identifiers, location, and summary.
- **Method:** Run `biomcp -j get gene <symbol>`. The command does not use a shell, its arguments are fixed, its output is limited to 1 MB, and it times out after 60 seconds. biomcp has no Python API, so this function uses its JSON command-line output. The executable is looked up on `PATH` first, then next to the current Python executable.
- **API:** `dnakit.integrations.get_biomcp_gene(symbol)`.
- **Input:** A gene symbol such as `BRCA1`, `HLA-A`, or `C9orf72`, up to 32 characters.
- **Output:** A dictionary with `symbol`, `name`, `entrez_id`, `ensembl_id`, `uniprot_id`, `omim_id`, `gene_type`, `location`, `coordinate`, `genome_build`, `aliases`, `clinical_diseases`, `summary`, and `provenance`. Coordinates and the summary are returned as biomcp provides them and are not recomputed locally.
- **Example code:**

```python
from dnakit.integrations import get_biomcp_gene

gene = get_biomcp_gene("BRCA1")
print(gene["symbol"], gene["entrez_id"], gene["location"], gene["genome_build"])
```

- **Example output:**

```text
BRCA1 672 17 GRCh38
```

## 5) `dnakit.integrations.search_biomcp_genes()`

- **Purpose:** Search gene records by free text.
- **Method:** Run `biomcp -j search gene -q <query> --limit <n>` and parse the returned `results`.
- **API:** `dnakit.integrations.search_biomcp_genes(query, *, limit=10)`.
- **Input:** A non-empty query string of at most 200 characters that does not start with `-`. `limit` is from 1 to 50.
- **Output:** A dictionary with `query`, `total`, `genes` (each item has `symbol`, `name`, `entrez_id`, `uniprot_id`, `omim_id`, and `coordinate`), and `provenance`.
- **Example code:**

```python
from dnakit.integrations import search_biomcp_genes

result = search_biomcp_genes("BRCA1", limit=2)
for gene in result["genes"]:
    print(gene["symbol"], gene["entrez_id"], gene["uniprot_id"])
```

- **Example output:**

```text
BRCA1 672 P38398
```

## 6) `dnakit.integrations.get_biomcp_variant()`

- **Purpose:** Look up variant annotations by rsID or genomic HGVS notation.
- **Method:** Run `biomcp -j get variant <id>`. `cadd_score` is returned as `None` when it is not a finite number.
- **API:** `dnakit.integrations.get_biomcp_variant(variant_id)`.
- **Input:** An rsID such as `rs334`, or a genomic HGVS string such as `chr7:g.140753336A>T`. It can be up to 128 characters, must not contain whitespace, and must not start with `-`.
- **Output:** A dictionary with `id`, `rsid`, `gene`, `hgvs_c`, `hgvs_p`, `consequence`, `significance`, `transcript`, `genome_build`, `cadd_score`, `sift_pred`, `polyphen_pred`, and `provenance`.
- **Example code:**

```python
from dnakit.integrations import get_biomcp_variant

variant = get_biomcp_variant("rs334")
print(variant["gene"], variant["hgvs_p"], variant["significance"])
```

- **Example output:**

```text
HBB p.Glu7Val Pathogenic
```

## 7) `dnakit.integrations.list_galaxy_tools()`

- **Purpose:** List the tools on a Galaxy server.
- **Method:** Use galaxy-mcp's `GalaxyInstance`, which does not follow redirects and rejects addresses that contain a username, password, query string, or fragment. The request goes through bioblend's `get_tools()`. The API key is sent in the `x-api-key` request header and never appears in the output. Tools are sorted by `id`. bioblend downloads the complete list first, and `limit` is applied afterward.
- **API:** `dnakit.integrations.list_galaxy_tools(url, *, api_key=None, limit=50)`.
- **Input:** The HTTP(S) address of a Galaxy server. `api_key` can be omitted, because public servers allow anonymous access to some endpoints. `limit` is from 1 to 500.
- **Output:** A dictionary with `url` (the normalized address), `total_tools`, `returned`, `tools` (each item has `id`, `name`, `version`, and `description`), and `provenance`.
- **Example code:**

```python
from dnakit.integrations import list_galaxy_tools

result = list_galaxy_tools("https://usegalaxy.org", limit=2)
print(result["url"], result["returned"], result["total_tools"] > 1000)
```

- **Example output:**

```text
https://usegalaxy.org/ 2 True
```

## 8) `dnakit.integrations.show_galaxy_tool()`

- **Purpose:** Show the basic metadata of one Galaxy tool.
- **Method:** Request the tool details through bioblend's `show_tool(tool_id)`, with the same redirect and credential rules as the previous function.
- **API:** `dnakit.integrations.show_galaxy_tool(url, *, tool_id, api_key=None)`.
- **Input:** A server address and a tool ID such as `Cut1`. The tool ID must not contain `..`.
- **Output:** A dictionary with `url`, `tool` (with `id`, `name`, `version`, `description`, and `model_class`), and `provenance`.
- **Example code:**

```python
from dnakit.integrations import show_galaxy_tool

result = show_galaxy_tool("https://usegalaxy.org", tool_id="Cut1")
print(result["tool"]["id"], result["tool"]["name"], result["tool"]["version"])
```

- **Example output:**

```text
Cut1 Cut 1.0.2
```

## 9) `dnakit.integrations.analyze_fasta_with_clawbio()`

- **Purpose:** Summarize one nucleotide FASTA file with clawbio's `analyze-fasta` skill.
- **Method:** Load `skills/analyze-fasta/analyze_fasta.py` from the installed clawbio package, call its `analyze_fasta()`, and reshape the results into per-sequence length, GC percentage, and base composition. The analysis timestamp and the absolute file path are removed from the output. The function raises an error when clawbio reports a problem or when the first record is classified as protein. Files larger than `max_bytes` are rejected before parsing.
- **API:** `dnakit.integrations.analyze_fasta_with_clawbio(path, *, max_bytes=10000000)`.
- **Input:** The path to an existing FASTA file. `max_bytes` is from 1 to 100000000.
- **Output:** A dictionary with `file`, `sequence_type`, `sequences` (each item has `index`, `id`, `length_bp`, `gc_percent`, and `composition`), `summary` (with `total_sequences`, `total_residues`, `min_length`, `max_length`, `n50`, and `avg_gc_percent`), and `provenance`.
- **Example file:** `sequences.fa`

```text
>seq1 test
ATGCGCAA
TT
>seq2
GGGCCCAAAT
```

- **Example code:**

```python
from dnakit.integrations import analyze_fasta_with_clawbio

result = analyze_fasta_with_clawbio("sequences.fa")
print(result["summary"]["total_sequences"], result["summary"]["n50"], result["summary"]["avg_gc_percent"])
for row in result["sequences"]:
    print(row["id"], row["length_bp"], row["gc_percent"])
```

- **Example output:**

```text
2 10 50.0
seq1 10 40.0
seq2 10 60.0
```

## Installing the optional dependencies

```bash
python -m pip install "dnakit[bio-toil]"          # toil and biopython, compute_toil_sequence_stats
python -m pip install "dnakit[bio-smartapi]"      # smartapi-mcp, SmartAPI search and spec loading
python -m pip install "dnakit[bio-biomcp]"        # biomcp-cli, the biomcp command
python -m pip install "dnakit[bio-galaxy]"        # galaxy-mcp and bioblend, Galaxy tool listing
python -m pip install "dnakit[bio-clawbio]"       # clawbio and biopython, analyze_fasta_with_clawbio
```

clawbio has many dependencies, so its installation takes noticeably longer than the other packages. Project links: [toil](https://github.com/DataBiosphere/toil), [smartapi-mcp](https://github.com/biothings/smartapi-mcp), [biomcp](https://github.com/genomoncology/biomcp), [galaxy-mcp](https://github.com/galaxyproject/galaxy-mcp), [bioblend](https://github.com/galaxyproject/bioblend), [clawbio](https://github.com/ClawBio/ClawBio).

## How these functions are verified

The test file is `others/tests/unit/test_service_backends.py`. It has 30 tests: 27 offline tests and 3 live tests, which run only when `DNAKIT_LIVE_TESTS=1` is set. The tests check the following:

- The Toil results match direct Biopython calls and hand calculations. For example, the reverse complement of ATGCGCAATT is AATTGCGCAT, and its GC fraction is 0.4.
- SmartAPI search runs against a locally mocked HTTP response. The order and scores match a direct call to `rank_apis`. APIs that share no query term with the search are excluded.
- biomcp is replaced by a POSIX script that records its arguments and prints fixed JSON. The tests check the fixed argument vector (an HGVS string is passed as one argument), JSON parsing, error JSON, and non-JSON output. The live test compares the results with a direct run of biomcp.
- Galaxy tool listing runs against a stub server on the loopback address. The tests check sorting, truncation, that the API key appears only in the request header, and that addresses with a username or query string are rejected. The live test queries usegalaxy.org.
- clawbio statistics match hand-computed lengths, GC percentages, and base composition, and they also match a direct call to the packaged script.
- When an optional package is missing, each module raises a `BackendUnavailableError` that names the `bio-*` extra. biomcp has no Python module, so its missing-tool test replaces the executable lookup instead.

## Current scope

This batch covers 9 functions from 5 upstream packages:

- Toil: only one two-job graph is wrapped. Cluster and cloud backends, and the CWL and WDL workflow interfaces, are not wrapped.
- smartapi-mcp: registry search and spec loading are wrapped. Server building and MCP tool registration are not.
- biomcp: gene and variant lookup and gene search are wrapped. Drug, disease, literature, trial, and `enrich` commands are not.
- galaxy-mcp: only tool lists and tool details are read. Running tools, uploading data, and importing workflows change server state, so they are not wrapped.
- clawbio: only the `analyze-fasta` skill is wrapped. Other skills need external programs or specific data, so they are not wrapped.
- `list_galaxy_tools()` downloads the server's complete tool list and then keeps `limit` items. Memory use for this download is not capped, so large servers need care.
