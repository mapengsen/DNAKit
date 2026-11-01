# just-dna-agents MCP (`just-dna-agents`)

DNAKit's MCP tool index keeps **2** DNA/gene tools of this server (6 listed, read on 2026-10-16). They run once you install the server; DNAKit does not bundle them. See [calling a tool that was found](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found).

| Item | Value |
| --- | --- |
| Source | [dna-seq/just-dna-agents @ `2896f4f`](https://github.com/dna-seq/just-dna-agents/tree/2896f4fe7b2f2e1b7a48dda2f721eab87123e770) |
| License | MIT |
| Command | `just-dna-agents-mcp serve --transport stdio` |
| Transports | stdio |
| Installation | Install the just-dna-agents-mcp and just-dna-agents packages from the pinned GitHub repository subdirectories (neither name resolves on PyPI although the README says the package is published) into a Python &gt;=3.14 environment; optional JUST_DNA_AGENTS_MCP_OUTPUT_DIR and JUST_DNA_AGENTS_MCP_RESOLVE_WITH_ENSEMBL settings control compilation. |
| Effects | validate_spec reads a caller-given spec directory; compile_module writes parquet files to the output directory (default current directory) and, with Ensembl resolution enabled, may use an Ensembl variation reference whose parquet cache the repository's resolver downloads from the Hugging Face dataset just-dna-seq/ensembl_variations; the compiler internals live in the separate just-dna-compiler package and were not inspected. |

| Tool | Category | Summary |
| --- | --- | --- |
| [`compile_module`](#compile_module) | `variant` | Validates a module spec directory and compiles it into parquet files (SNP weights, annotations, studies, PGS and PGx tables), optionally resolving rsIDs via En… |
| [`validate_spec`](#validate_spec) | `variant` | Validates a module spec directory (module_spec.yaml, variants.csv, studies.csv) and returns errors, warnings and stats. |

## `compile_module` {#compile_module}

Validates a module spec directory and compiles it into parquet files (SNP weights, annotations, studies, PGS and PGx tables), optionally resolving rsIDs via Ensembl.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `spec_dir`
- **Keywords:** snp, rsid, polygenic score, pharmacogenomics, parquet, annotation module, genotype interpretation

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `spec_dir` | string | yes | Path to the module spec directory. |
| `output_dir` | string \| null | no | Output directory for parquet files. If not provided, uses the configured default output directory. Default `None` |
| `compression` | string | no | Parquet compression codec (zstd, snappy, lz4, gzip). Default `zstd` |

??? note "Upstream description"

    Compile a module spec directory into deployable parquet files.

    Validates the spec, optionally resolves missing rsid/position via
    Ensembl DuckDB, and produces the composed parquet artifact in the
    output directory: the SNP core (weights/annotations/studies) when the
    module carries variants, plus one parquet per 0.4 table kind present
    (pgs, diplotypes, pharm_variants, the binning kinds, ...).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'just-dna-agents',
    'compile_module',
    {'spec_dir': '...'},
)
```

## `validate_spec` {#validate_spec}

Validates a module spec directory (module_spec.yaml, variants.csv, studies.csv) and returns errors, warnings and stats.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `spec_dir`
- **Keywords:** variants csv, snp weights, validation, annotation module, studies

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `spec_dir` | string | yes | Path to the module spec directory containing module_spec.yaml, variants.csv, and studies.csv. |

??? note "Upstream description"

    Validate a module spec directory.

    Checks YAML structure, CSV row validity, cross-row consistency,
    and weight/state directionality. Returns a validation result dict
    with 'valid', 'errors', 'warnings', and 'stats' keys.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'just-dna-agents',
    'validate_spec',
    {'spec_dir': '...'},
)
```
