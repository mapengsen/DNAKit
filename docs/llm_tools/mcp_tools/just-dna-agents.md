# just-dna-agents MCP (`just-dna-agents`)

DNAKit 的 MCP 工具索引收录了这个服务器的 **2** 个 DNA/基因相关工具（服务器共列出 6 个，2026-10-16 读取）。工具在用户自己安装服务器后运行，DNAKit 本身不包含它们。调用方式见 [调用索引中的工具](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found)。下面的工具说明是上游英文原文。

| 项目 | 内容 |
| --- | --- |
| 源码 | [dna-seq/just-dna-agents @ `2896f4f`](https://github.com/dna-seq/just-dna-agents/tree/2896f4fe7b2f2e1b7a48dda2f721eab87123e770) |
| 许可证 | MIT |
| 启动命令 | `just-dna-agents-mcp serve --transport stdio` |
| 传输方式 | stdio |
| 安装 | Install the just-dna-agents-mcp and just-dna-agents packages from the pinned GitHub repository subdirectories (neither name resolves on PyPI although the README says the package is published) into a Python &gt;=3.14 environment; optional JUST_DNA_AGENTS_MCP_OUTPUT_DIR and JUST_DNA_AGENTS_MCP_RESOLVE_WITH_ENSEMBL settings control compilation. |
| 副作用 | validate_spec reads a caller-given spec directory; compile_module writes parquet files to the output directory (default current directory) and, with Ensembl resolution enabled, may use an Ensembl variation reference whose parquet cache the repository's resolver downloads from the Hugging Face dataset just-dna-seq/ensembl_variations; the compiler internals live in the separate just-dna-compiler package and were not inspected. |

| 工具 | 类别 | 简介 |
| --- | --- | --- |
| [`compile_module`](#compile_module) | `variant` | Validates a module spec directory and compiles it into parquet files (SNP weights, annotations, studies, PGS and PGx tables), optionally resolving rsIDs via En… |
| [`validate_spec`](#validate_spec) | `variant` | Validates a module spec directory (module_spec.yaml, variants.csv, studies.csv) and returns errors, warnings and stats. |

## `compile_module` {#compile_module}

Validates a module spec directory and compiles it into parquet files (SNP weights, annotations, studies, PGS and PGx tables), optionally resolving rsIDs via Ensembl.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`spec_dir`
- **关键词：**snp, rsid, polygenic score, pharmacogenomics, parquet, annotation module, genotype interpretation

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `spec_dir` | string | 是 | Path to the module spec directory. |
| `output_dir` | string \| null | 否 | Output directory for parquet files. If not provided, uses the configured default output directory. 默认 `None` |
| `compression` | string | 否 | Parquet compression codec (zstd, snappy, lz4, gzip). 默认 `zstd` |

??? note "上游说明（英文原文）"

    Compile a module spec directory into deployable parquet files.

    Validates the spec, optionally resolves missing rsid/position via
    Ensembl DuckDB, and produces the composed parquet artifact in the
    output directory: the SNP core (weights/annotations/studies) when the
    module carries variants, plus one parquet per 0.4 table kind present
    (pgs, diplotypes, pharm_variants, the binning kinds, ...).

安装服务器后通过 DNAKit 调用：

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

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`spec_dir`
- **关键词：**variants csv, snp weights, validation, annotation module, studies

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `spec_dir` | string | 是 | Path to the module spec directory containing module_spec.yaml, variants.csv, and studies.csv. |

??? note "上游说明（英文原文）"

    Validate a module spec directory.

    Checks YAML structure, CSV row validity, cross-row consistency,
    and weight/state directionality. Returns a validation result dict
    with 'valid', 'errors', 'warnings', and 'stats' keys.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'just-dna-agents',
    'validate_spec',
    {'spec_dir': '...'},
)
```
