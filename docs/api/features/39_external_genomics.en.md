# External genomics analyses and MCP

Use one Python API to discover registered capabilities, run external analyses and connect to third-party MCP services. See the [parameter FAQ](../../faq/api/39_external_genomics.md) for installation, full signatures and execution details, and the [source FAQ](../../faq/genomics_sources.md) for authors, papers and licenses.

For scientific dependencies installed together through pip, start with the [`genomics` extra and tutorial](38_genomics_extensions.md#pip-genomics). It provides in-process genotype PCA, ANI and prokaryotic CDS prediction. The `run_analysis()` adapters on this page still use their registered external executables; input and algorithm scopes are documented separately. See the [43-analysis migration FAQ](../../faq/external_analysis_migration.md).

## 1) Discover projects and executable analyses

- **Purpose:** Find relevant projects and registered DNAKit analyses, including their required inputs and backends.
- **Method:** Match keywords against local catalog records and filter by domain; `describe_analysis` retrieves the definition for one analysis ID. Discovery does not run an analysis.
- **API:** `dnakit.integrations.search_projects(query[optional], domain[optional], verbose[optional])`; `dnakit.integrations.list_analyses(query[optional], domain[optional], verbose[optional])`; `dnakit.integrations.describe_analysis(analysis[required], verbose[optional])`.
- **Input:** A project keyword, domain or analysis ID, such as `plink2`, `population` or `plink2-pca`.
- **Output:** Returns a project list, available analyses, and analysis-definition mappings.
- **Example code:**

```python
from dnakit.integrations import search_projects, list_analyses, describe_analysis

projects = search_projects("plink2", verbose=False)
analyses = list_analyses("plink2", verbose=False)
definition = describe_analysis("plink2-pca", verbose=False)
print([project["id"] for project in projects["projects"]])
print([analysis["id"] for analysis in analyses])
print(definition["backend"], sorted(definition["inputs"]))
```

- **Example result:** Offline results from the current catalog.

```text
['plink2']
['plink2-pca', 'plink2-king', 'plink2-ld-prune', 'plink2-gwas']
plink2 ['genotypes']
```

## 2) External software

- **Purpose:** Plan a registered analysis and execute it once the backend and inputs are ready.
- **Method:** Validate input files, parameters and the output directory, then build argument vectors from registered templates. `plan_analysis` returns the plan; `run_analysis` invokes the backend, checks output files and records the run.
- **API:** `dnakit.integrations.plan_analysis(analysis[required], inputs[required], output_directory[required], parameters[optional], threads[optional], executable_paths[optional], verbose[optional])`; `dnakit.integrations.run_analysis(analysis[required], inputs[required], output_directory[required], parameters[optional], threads[optional], executable_paths[optional], timeout_seconds[optional], max_artifact_bytes[optional], verbose[optional])`.
- **Input:** An analysis ID, input file mapping and a new output directory. The example requires existing `data/cohort.bed`, `data/cohort.bim` and `data/cohort.fam` files.
- **Output:** `plan_analysis` returns an execution plan; execution returns status, steps, and artifact paths.
- **Example code:**

```python
from dnakit.integrations import plan_analysis, run_analysis

inputs = {"genotypes": "data/cohort"}
plan = plan_analysis(
    "plink2-pca", inputs, "results/cohort-pca",
    parameters={"components": 3}, threads=4, verbose=False,
)
print(plan["analysis"], plan["status"])
print(plan["parameters"])
print(plan["steps"][0]["expected_outputs"])

# With PLINK 2 installed and the inputs and plan checked:
# result = run_analysis(
#     "plink2-pca", inputs, "results/cohort-pca",
#     parameters={"components": 3}, threads=4,
# )
```

- **Example result:** Plan output when the input files exist; this step neither runs PLINK nor writes result files.

```text
plink2-pca planned
{'components': 3}
['pca.eigenvec', 'pca.eigenval']
```

<details markdown="1">
<summary>Registered analyses and backends</summary>

<!-- generated:analyses:start -->
| Analysis ID / 分析 ID | Capability / 功能 | Backend / 后端 | Validation / 验证 |
| --- | --- | --- | --- |
| `plink2-pca` | Population PCA | plink2 | passed_synthetic_fixture |
| `plink2-king` | KING robust kinship | plink2 | passed_synthetic_fixture |
| `plink2-ld-prune` | LD pruning | plink2 | passed_synthetic_fixture |
| `plink2-gwas` | GWAS with covariates | plink2 | passed_synthetic_fixture |
| `gcta-grm` | GCTA GRM | gcta64 | not_run |
| `gcta-reml` | GCTA GREML heritability | gcta64 | not_run |
| `minigraph-pangenome` | Reference plus assembly graph | minigraph | not_run |
| `hifiasm-assembly` | HiFi assembly | hifiasm | not_run |
| `flye-assembly` | Nanopore repeat-graph assembly | flye | not_run |
| `sniffles-sv` | Long-read structural variants | sniffles | not_run |
| `modkit-pileup` | Modified-base pileup | modkit | not_run |
| `mosdepth-coverage` | BAM sequencing depth | mosdepth | not_run |
| `fastani` | Average nucleotide identity | fastANI | not_run |
| `kraken2-classify` | Taxonomic read classification | kraken2 | not_run |
| `amrfinder-nucleotide` | AMRFinderPlus nucleotide annotation | amrfinder | not_run |
| `prokka-annotation` | Prokaryotic annotation | prokka | not_run |
| `barrnap-rrna` | Genomic rRNA annotation | barrnap | not_run |
| `fimo-motifs` | Motif occurrence scanning | fimo | not_run |
| `cas-offinder` | CRISPR off-target enumeration | cas-offinder | not_run |
| `bakta-annotation` | Bacterial genome annotation with Bakta | bakta | not_run |
| `busco-proteins` | Protein completeness assessment with BUSCO | busco | not_run |
| `busco-genome` | Genome completeness assessment with BUSCO | busco | not_run |
| `augustus-genes` | Species-trained gene prediction with AUGUSTUS | augustus | not_run |
| `diamond-blastp` | Protein similarity search with DIAMOND | diamond | passed_synthetic_smoke |
| `diamond-blastx` | Translated DNA similarity search with DIAMOND | diamond | passed_synthetic_smoke |
| `eggnog-proteins` | Protein functional annotation with eggNOG-mapper | emapper.py | not_run |
| `miniprot-align` | Protein-to-genome alignment with miniprot | miniprot | passed_synthetic_smoke |
| `liftoff-annotation` | Reference gene annotation transfer with Liftoff | liftoff | not_run |
| `checkm2-predict` | Genome completeness and contamination with CheckM2 | checkm2 | not_run |
| `quast-assembly` | Reference-based assembly evaluation with QUAST | quast.py | not_run |
| `vsearch-dereplicate` | Nucleotide sequence dereplication with VSEARCH | vsearch | passed_synthetic_smoke |
| `vsearch-cluster` | Nucleotide sequence clustering with VSEARCH | vsearch | passed_synthetic_smoke |
| `mmseqs-protein-search` | Protein similarity search with MMseqs2 | mmseqs | passed_synthetic_smoke |
| `mmseqs-protein-cluster` | Protein sequence clustering with MMseqs2 | mmseqs | passed_synthetic_smoke |
| `abricate-screen` | Assembly gene screening with ABRicate | abricate | not_run |
| `mlst-typing` | Multilocus sequence typing with MLST | mlst | not_run |
| `genomad-mobile-elements` | Virus and plasmid detection with geNomad | genomad | not_run |
| `repeatmasker-library` | Custom-library repeat masking with RepeatMasker | RepeatMasker | not_run |
| `mageck-test` | CRISPR screen enrichment testing with MAGeCK | mageck | not_run |
| `mageck-count` | CRISPR guide read counting with MAGeCK | mageck | not_run |
| `crispresso2-amplicon` | Amplicon editing outcome analysis with CRISPResso2 | CRISPResso | not_run |
| `drugz-unpaired` | Unpaired chemical-genetic screen analysis with DrugZ | drugz.py | passed_synthetic_smoke |
| `gemma-lmm` | Linear mixed-model association with GEMMA | gemma | not_run |
<!-- generated:analyses:end -->

</details>

## 3) Third-party MCP

- **Purpose:** Discover registered services, then select and call third-party capabilities through their live tool schemas.
- **Method:** `list_mcp_servers` reads the local catalog; `list_mcp_tools` connects to a service and retrieves schemas; `call_mcp_tool` validates arguments and invokes one tool.
- **API:** `dnakit.integrations.list_mcp_servers(verbose[optional])`; `dnakit.integrations.list_mcp_tools(server_id[required], transport[optional], executable_path[optional], workdir[optional])`; `dnakit.integrations.call_mcp_tool(server_id[required], tool_name[required], arguments[required], transport[optional], executable_path[optional], workdir[optional])`. Connection timeouts, environment variables and response limits are covered in the [parameter FAQ](../../faq/api/39_external_genomics.md#mcp-connection).
- **Input:** A catalog service ID. Live calls also require a configured service and a tool name and arguments from its current schema.
- **Output:** Returns lists of MCP servers and tools, with service IDs, commands, or tool definitions.
- **Example code:**

```python
from dnakit.integrations import list_mcp_servers, list_mcp_tools, call_mcp_tool

server = next(
    item for item in list_mcp_servers(verbose=False) if item["id"] == "biomcp"
)
print(server["id"], server["command"])
print(server["arguments"])

# With BioMCP configured, retrieve live tools; this starts the service:
# tools = list_mcp_tools("biomcp")
# Prepare a name and arguments from the returned inputSchema, then call:
# result = call_mcp_tool("biomcp", tool_name, arguments)
```

- **Example result:** Offline service metadata; no service or remote tool was started.

```text
biomcp biomcp
['serve']
```

<details markdown="1">
<summary>Registered MCP services and upstream configuration</summary>

<!-- generated:mcp:start -->
| Server ID / 服务 ID | Upstream / 上游 | Configuration / 配置 |
| --- | --- | --- |
| `gigwa-mcp` | [Gigwa MCP](https://github.com/gkanogiannis/Gigwa-MCP) | Install gigwa-mcp separately into a Python >=3.10 environment; configure GIGWA_URL and optional paired GIGWA_USER/GIGWA_PASS. |
| `chorus` | [Chorus](https://github.com/pinellolab/chorus) | Install the pinellolab/chorus source project and provision its oracle environments, model weights and references separately. The PyPI package named chorus is unrelated. Set executable_path to the installed chorus-mcp. |
| `mavedb-link` | [MaveDB-Link](https://github.com/berntpopp/mavedb-link) | For stdio, install the source project's dependencies in Python >=3.12 and set workdir to that checkout and executable_path to its Python. The author-documented HTTP service is available with transport='http'. |
| `biomcp` | [BioMCP](https://github.com/genomoncology/biomcp) | Install genomoncology/biomcp separately; its PyPI distribution is biomcp-cli, not the unrelated biomcp distribution. The installed executable is biomcp. |
| `juicebox` | [Juicebox MCP](https://github.com/aidenlab/juicebox-mcp) | Prepare the official source checkout and npm dependencies separately, then set workdir to the checkout and executable_path to Node if necessary. |
| `genechat` | [GeneChat MCP](https://github.com/natecostello/genechat-mcp) | Install and initialize GeneChat and its local genome/annotation databases separately. GENECHAT_CONFIG optionally selects an existing TOML configuration. |
| `genepattern` | [GenePattern MCP](https://github.com/genepattern/genepattern-mcp) | Install the official source checkout requirements separately with Python >=3.11. Set workdir to the checkout, executable_path to that Python, and GENEPATTERN_KEY to an existing API token. |
| `brapi` | [BrAPI MCP](https://github.com/cyanheads/brapi-mcp-server) | Separately prepare the pinned @cyanheads/brapi-mcp-server release/source with its dist build and dependencies; set workdir to the package root and executable_path to Node >=24. The bridge never runs npm, npx, bunx or an installer. |
| `ensembl` | [Ensembl MCP (cyanheads)](https://github.com/cyanheads/ensembl-mcp-server) | Separately prepare the pinned @cyanheads/ensembl-mcp-server build and dependencies; set workdir to the package root and executable_path to Node >=24. |
| `biothings` | [BioThings MCP](https://github.com/longevity-genie/biothings-mcp) | Install the pinned longevity-genie/biothings-mcp project separately in Python >=3.10; use its installed biothings-mcp entry point. |
| `gget` | [gget MCP](https://github.com/longevity-genie/gget-mcp) | Install gget-mcp and its separately specified gget dependencies in Python >=3.10; the verified command is gget-mcp stdio. |
| `ucsc-genome` | [UCSC Genome MCP](https://github.com/hlydecker/ucsc-genome-mcp) | Prepare the pinned source checkout and its requirements in Python >=3.10, set workdir to the checkout and executable_path to that Python. Use the existing hyphenated ucsc-genome-mcp.py script; upstream README/packaging also mention nonexistent underscore filenames. |
| `opentargets` | [Open Targets MCP](https://github.com/nickzren/opentargets-mcp) | Install nickzren/opentargets-mcp separately in Python >=3.10; set executable_path to its interpreter. Uses the verified Python module entry point with explicit stdio. |
| `plant-genomics` | [Plant Genomics MCP](https://github.com/musharna/plant-genomics-mcp) | Install the pinned plant-genomics-mcp Python distribution and dependencies separately. The plant-genomics-mcp entry point selects stdio; no HTTP token is needed for this local transport. |
| `biocontext-kb` | [BioContextAI Knowledgebase MCP](https://github.com/biocontext-ai/knowledgebase-mcp) | Install biocontext_kb separately in Python >=3.12. DNAKit fixes MCP_ENVIRONMENT=DEVELOPMENT for the verified stdio entry point. |
| `grebi` | [GrEBI](https://github.com/EBISPOT/GrEBI) | Select transport='http' for the author-documented development service. Local deployment is a separate Java/API plus graph-database and HPC data-pipeline task; no local startup recipe is claimed here. |
| `ebi-tools` | [EBI Tools MCP (QuentinCody)](https://github.com/QuentinCody/ebi-tools-mcp-server) | Select transport='http' for the repository author's published Cloudflare Worker MCP service. No standalone Node stdio entry point is present; no mcp-remote/npx helper is required by DNAKit. |
| `gwas-mcp` | [GWAS MCP (muslus)](https://github.com/muslus/gwas-mcp) | Prepare the pinned muslus/gwas-mcp source and its dependencies in an isolated Python >=3.10 environment, then point executable_path at that interpreter. Avoid assuming that the separately published gwas-mcp distribution has identical code. |
| `ncbi-database` | [NCBI Database MCP (hpend2373)](https://github.com/hpend2373/NCBI-Database-MCP) | Prepare the pinned NCBI-Database-MCP checkout and its requirements, set workdir to the checkout and executable_path to its Python. Uses the README-supported standard MCP script rather than executing its shell installer/launcher. |
| `pubmed` | [PubMed MCP (cyanheads)](https://github.com/cyanheads/pubmed-mcp-server) | Prepare the pinned @cyanheads/pubmed-mcp-server build and dependencies separately, then set workdir to the package root and executable_path to Node >=24. |
| `encode-toolkit` | [ENCODE Toolkit](https://github.com/ammawla/encode-toolkit) | Install the Python encode-toolkit distribution separately and set executable_path to its interpreter. DNAKit invokes the verified Python module directly; it does not use the npm shim which can bootstrap uv/Python. |
| `tooluniverse` | [ToolUniverse](https://github.com/mims-harvard/ToolUniverse) | Install the pinned ToolUniverse distribution and required provider dependencies separately. Use the installed tooluniverse console script; unlike the README uvx --refresh example, this bridge does not refresh or install anything. |
| `kiro-genomics` | [AWS Kiro Life Sciences: genomics](https://github.com/aws-samples/sample-kiro-power-life-sciences) | Prepare the pinned aws-samples/sample-kiro-power-life-sciences source. Separately install life-sciences-common and the life-sciences-genomics component in its Python >=3.10 environment, then set executable_path to the installed life-sciences-genomics console script. No uvx or automatic install is used. |
| `kiro-model-organisms` | [AWS Kiro Life Sciences: model-organisms](https://github.com/aws-samples/sample-kiro-power-life-sciences) | Prepare the pinned aws-samples/sample-kiro-power-life-sciences source. Separately install life-sciences-common and the life-sciences-model-organisms component in its Python >=3.10 environment, then set executable_path to the installed life-sciences-model-organisms console script. No uvx or automatic install is used. |
| `kiro-epigenomics` | [AWS Kiro Life Sciences: epigenomics](https://github.com/aws-samples/sample-kiro-power-life-sciences) | Prepare the pinned aws-samples/sample-kiro-power-life-sciences source. Separately install life-sciences-common and the life-sciences-epigenomics component in its Python >=3.10 environment, then set executable_path to the installed life-sciences-epigenomics console script. No uvx or automatic install is used. |
| `kiro-agriculture` | [AWS Kiro Life Sciences: agriculture](https://github.com/aws-samples/sample-kiro-power-life-sciences) | Prepare the pinned aws-samples/sample-kiro-power-life-sciences source. Separately install life-sciences-common and the life-sciences-agriculture component in its Python >=3.10 environment, then set executable_path to the installed life-sciences-agriculture console script. No uvx or automatic install is used. |
| `kiro-ecology` | [AWS Kiro Life Sciences: ecology](https://github.com/aws-samples/sample-kiro-power-life-sciences) | Prepare the pinned aws-samples/sample-kiro-power-life-sciences source. Separately install life-sciences-common and the life-sciences-ecology component in its Python >=3.10 environment, then set executable_path to the installed life-sciences-ecology console script. No uvx or automatic install is used. |
| `kiro-clinical` | [AWS Kiro Life Sciences: clinical](https://github.com/aws-samples/sample-kiro-power-life-sciences) | Prepare the pinned aws-samples/sample-kiro-power-life-sciences source. Separately install life-sciences-common and the life-sciences-clinical component in its Python >=3.10 environment, then set executable_path to the installed life-sciences-clinical console script. No uvx or automatic install is used. |
| `kiro-pipelines` | [AWS Kiro Life Sciences: pipelines](https://github.com/aws-samples/sample-kiro-power-life-sciences) | Prepare the pinned aws-samples/sample-kiro-power-life-sciences source. Separately install life-sciences-common and the life-sciences-pipelines component in its Python >=3.10 environment, then set executable_path to the installed life-sciences-pipelines console script. No uvx or automatic install is used. |
| `pubmed-anthropic` | [PubMed (Anthropic Life Sciences hosted connector)](https://github.com/anthropics/life-sciences) | Select transport=http. The provider-published endpoint supports anonymous MCP initialization and tools/list as verified on 2026-09-20. No local installer is run. |
| `opentargets-platform` | [Open Targets Platform official MCP](https://github.com/anthropics/life-sciences) | Select transport=http. The provider-published endpoint supports anonymous MCP initialization and tools/list as verified on 2026-09-20. No local installer is run. |
| `10x-genomics` | [10x Genomics Cloud MCP](https://github.com/10XGenomics/txg-mcp) | Separately prepare and unpack the upstream 10x Genomics MCP bundle, including its server/, node_modules/ and platform bin/ files. Set workdir to that unpacked bundle, executable_path to installed Node.js >=18, and the named TXG_CLI_ACCESS_TOKEN environment variable. DNAKit does not download a bundle or TXG binary. |
<!-- generated:mcp:end -->

</details>

## 4) Stateful workflows

- **Purpose:** Execute successive MCP calls in one session while preserving service state, such as a loaded model.
- **Method:** Refresh schemas, validate explicit arguments and invoke each step in order. A remote error stops the sequence and returns the executed steps and error status.
- **API:** `dnakit.integrations.call_mcp_workflow(server_id[required], calls[required], transport[optional], executable_path[optional], workdir[optional], timeout_seconds[optional])`; see the [parameter FAQ](../../faq/api/39_external_genomics.md#mcp-workflow) for other connection options.
- **Input:** A service ID and a list of calls, each with a `tool_name` and `arguments` from live schemas. Supply actual load and prediction tool names and arguments to the function below.
- **Output:** Returns a mapping of tool-call results from the same session.
- **Example code:**

```python
from dnakit.integrations import call_mcp_workflow


def load_then_predict(load_tool, load_arguments, predict_tool, predict_arguments):
    calls = [
        {"tool_name": load_tool, "arguments": load_arguments},
        {"tool_name": predict_tool, "arguments": predict_arguments},
    ]
    result = call_mcp_workflow("chorus", calls, timeout_seconds=300)
    print(result["status"], result["requested_steps"], result["executed_steps"])
    return result
```

- **Example result:** Output form when both calls succeed. Chorus was not connected here; actual results depend on the service.

```text
completed 2 2
```

## 5) DNAKit MCP exposure

- **Purpose:** Expose genomics and external-tool APIs to MCP clients, with input schemas and execution metadata available before invocation.
- **Method:** DNAKit builds a tool catalog, signatures and JSON schemas from public functions. Clients discover a tool and then invoke it. The example reads execution metadata from the same catalog offline.
- **API:** Server entry point `dnakit-mcp`; MCP tools `search_dnakit_tools(query[required], category[optional], limit[optional])` and `describe_dnakit_tool(tool_name[required])`; local catalog entry point `dnakit.tools.default_tool_registry(verbose[optional])`.
- **Input:** The keywords `genomics` or `integrations`, and a full tool name such as `dnakit.integrations.run_analysis`.
- **Output:** Returns DNAKit MCP tool descriptions or execution results.
- **Example code:**

```python
from dnakit.tools import default_tool_registry

registry = default_tool_registry(verbose=False)
info = registry.describe("dnakit.integrations.run_analysis")
print(info["name"])
print(info["effect"])
```

- **Example result:** Offline tool metadata; no analysis was executed.

```text
dnakit.integrations.run_analysis
external_process
```
