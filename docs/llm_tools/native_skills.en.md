# Unified Skill APIs

For direct function usage, see the [complete 484-function index](../api/native/index.md) and [calling examples](../api/native/usage.md). All functions are grouped by retrieval, querying, analysis, models and engineering, with individual parameters and returns.

DNAKit provides search, installation and execution across the skill catalog. Collection-specific database adapters and shared workflow code are bundled in DNAKit with their original attribution. Scientific programs, libraries and trained models, including samtools, Biopython and AlphaGenome, keep their original implementations and requirements.

These changes are not yet published to PyPI. Install the current checkout:

```bash
python -m pip install -e "."
dnakit skills search "variant annotation" --limit 5
dnakit skills info variant-annotation
dnakit skills doctor variant-annotation
dnakit skills install variant-annotation --agent codex
```

After a release containing these changes, `python -m pip install dnakit` can replace the checkout installation. Collection prefixes are removed selectively; old names remain aliases. For example, `tooluniverse-sequence-retrieval` becomes `sequence-retrieval`, while `virulencefinder` and `alphagenome-atlas-website-links` keep their names.

## Database functions

Discover the operation and inspect its parameters before calling it:

```bash
dnakit database search "gene lookup"
dnakit database info ensembl_lookup_gene
dnakit database call ensembl_lookup_gene '{"id":"BRCA1","species":"homo_sapiens"}'
```

```python
import dnakit

gene = dnakit.databases.ensembl_lookup_gene(id="BRCA1", species="homo_sapiens")
variants = dnakit.databases.gnomad_get_sv_by_region(
    chrom="17", start=43044295, end=43070295,
)
```

There is no ToolUniverse initialization or tool-loading step. `describe()` returns DNAKit's accepted `input_schema` together with the original service schema in `config.parameter`. Conflicting aliases, invalid arguments, provider failures and oversized results raise exceptions. Results preserve the provider's dictionary/list structure; do not assume every function has the same `data` envelope.

Optional parameters have explicit behavior: identifiers and coordinate aliases are normalized; Open Targets `size` limits the current returned page and indicates truncation; GWAS `p_value_threshold` filters the current page. Neither option establishes exhaustive database coverage.

## Bundled scripts

Resolve scripts from the installed package without cloning an upstream collection:

```bash
dnakit skills scripts genomic-coordinates
dnakit skills run genomic-coordinates --script scripts/convert_coords.py -- \
  --from bed --to gff chr1 999 1000 --format json
```

```python
import dnakit

result = dnakit.workflows.run(
    "genomic-coordinates", "scripts/convert_coords.py",
    arguments=["--from", "bed", "--to", "gff", "chr1", "999", "1000", "--format", "json"],
    timeout_seconds=60,
)
print(result["returncode"], result["output"])
```

This converts a BED interval to GFF coordinates; it does not change genome assembly. The CLI forwards standard input/output and exit status. The Python API bounds execution time and captured output. Bundled scripts and helper resources are checked against recorded SHA-256 hashes before execution. Actual scientific dependencies must be available in the execution environment.

## ENCODE, composite workflows and reports

```python
import json
import dnakit

experiments = json.loads(dnakit.encode.search_experiments(
    assay_title="ATAC-seq", organism="Homo sapiens", limit=3,
))
```

ENCODE functions preserve their options and JSON-text return type. Each call owns and closes its client resources. Tracking and fallback credentials use `~/.dnakit-encode/`; persistent credential storage requires `dnakit[encode-credentials]`. The adapted ENCODE implementation retains AGPL-3.0-only licensing.

Composite workflows can call `dnakit skills support-run NAME -- ...` or `dnakit.workflows.run_support()`. Bundled supporting stages are `clinical-variant-reporter`, `clinical-trial-finder`, `omics-target-evidence-mapper`, `gwas-lookup` and `fine-mapping`. Their scientific outputs still need appropriate validation.

DNAKit's independently written PDF renderer accepts existing Markdown content:

```bash
python -m pip install -e ".[reports]"
dnakit skills report --input wes_report.md --output report-pdf --language en
```

`dnakit.workflows.write_wes_report()` supports language labels, title, logos, fonts and overwrite options, and records input/output hashes. It does not generate new medical interpretations. No proprietary upstream English/Spanish report implementation was copied.

## MCP and validation

With `dnakit[agent]` and the [MCP configuration](mcp.md), agents discover these functions through `search_dnakit_tools`, inspect them with `describe_dnakit_tool`, and execute them with `call_dnakit_tool`. File writes, external processes and model execution retain their corresponding operation classifications. The initial context does not contain every function definition.

Read each skill's `DNAKIT.md` and inspect `dnakit skills info ID --json`. `registry/migration.json` records integration details; `registry/adaptations.json` binds maintained templates to source revisions. File checks, signature checks, unit tests, live service calls and complete scientific workflow validation are distinct. Unexecuted workflows remain `not_run`.

For ChIP-Atlas, `ch_ipatlas_get_peak_data(gene=...)` resolves a GRCh38 locus and returns a track retrieval plan with `status="requires_data"`. Supply the downloaded, verified BED as `peak_file` to compute overlap. The inherited enrichment adapter returns submission instructions, marked `requires_submission`, rather than computed statistics. A track URL or peak overlap does not establish motif disruption. These distinctions follow the [official ChIP-Atlas API](https://chip-atlas.org/agents).

Original authors, licenses and file hashes are retained. Coordinate-file retrieval follows [RCSB documentation](https://www.rcsb.org/docs/programmatic-access/file-download-services), and bounded Census queries follow the [CELLxGENE API](https://chanzuckerberg.github.io/cellxgene-census/_autosummary/cellxgene_census.get_anndata.html).
