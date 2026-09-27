# DNAKit Skills: search, install, and use

DNAKit collects DNA/genomics workflows with offline discovery, selective installation, provenance, and MCP access. A Skill contains instructions and supporting resources. Its analysis software, reference data, and APIs remain separate dependencies.

Start with the **`dnakit` discovery skill**. Describe your task; the Agent searches the catalog, reads a suitable workflow, and uses the available tools to perform it.

## Install DNAKit and the discovery skill

This feature requires a DNAKit release containing the skill catalog. From the current source checkout:

```bash
python -m pip install -e "."
dnakit skills stats
```

After a release containing this feature is published, users can obtain it with `python -m pip install -U dnakit`. Discovery and installation do not require a model API or MCP.

In your analysis project, choose one client:

```bash
dnakit skills install dnakit --agent codex
# Or:
dnakit skills install dnakit --agent claude-code
# Or:
dnakit skills install dnakit --agent cursor
```

| Client | Project directory | User directory with `--global` |
| --- | --- | --- |
| Codex | `.agents/skills/dnakit/` | `~/.agents/skills/dnakit/` |
| Claude Code | `.claude/skills/dnakit/` | `~/.claude/skills/dnakit/` |
| Cursor | `.cursor/skills/dnakit/` | `~/.cursor/skills/dnakit/` |

See the official [Codex](https://learn.chatgpt.com/docs/build-skills), [Claude Code](https://code.claude.com/docs/en/skills), and [Cursor](https://cursor.com/docs/skills) documentation. Start a new client session after installation, and make the DNAKit Python environment accessible to the client.

For example, ask Codex:

```text
$dnakit Convert BED interval chr1 999 1000 to GFF coordinates.
Find a suitable skill, check dependencies, run the conversion, and explain the conventions.
```

The expected route is discovery of `genomic-coordinates`, reading its instructions, preparing its runtime, and running its conversion script. The result is `chr1 1000 1000`, a 1 bp interval.

## Plan a complex task with `dnakit-plan`

The catalog includes `dnakit-plan`, a Codex skill for checking DNAKit MCP capabilities, resolving consequential questions with the user, and saving an approved Markdown plan for a later session. Its canonical source is `skills/dnakit-plan/`; this repository links `.agents/skills/dnakit-plan/` to that directory so local edits are immediately visible to Codex. To use it in another project, install the cataloged version:

```bash
dnakit skills install dnakit-plan --agent codex
# Add --global to use it in every project.
```

Open a new Codex session after installation. Installing the skill does not configure DNAKit MCP; follow the [MCP guide](mcp.md) to connect it. Codex sessions started at this repository root can use the project-local link without installing it again.

Invoke `$dnakit-plan` explicitly when planning must happen; installing the skill alone does not force every request through planning.

First follow the [MCP tutorial](mcp.md) and use `/mcp` in Codex to check that the DNAKit tools are available. Then provide only the task goal, for example:

```text
$dnakit-plan Download the hg38 reference genome and prepare training data.
```

The skill requires Codex to search DNAKit MCP for the needed tools, inspect their schemas, and record the exact tool names, key arguments, inputs, outputs, and checks for each step in the Markdown plan. It also defines consequential questions, saving after confirmation, and stopping before execution. Answer its questions and confirm the final plan; in a mode that permits file writes, Codex saves it under `plans/` without another workflow prompt. In a later session, simply ask Codex to execute `plans/<task-name>.md`. `/plan` is an optional extra planning mode; if it prevents file writes, switch back to regular mode after confirming the plan. A skill guides this process; a program-level confirmation gate is needed when execution before approval must be impossible.

## Find and install an individual workflow

```bash
dnakit skills search "variant annotation" --limit 5
dnakit skills search "变异注释" --limit 5
dnakit skills search "GWAS" --category gwas --limit 5
dnakit skills categories
dnakit skills info variant-annotation
dnakit skills doctor variant-annotation
dnakit skills install variant-annotation --agent codex
```

Search uses local keyword matching and curated Chinese query aliases. It makes no model API calls; summaries and workflow text returned to an Agent still consume conversation tokens. Skill discovery finds workflow instructions; tool discovery finds callable functions.

`doctor` probes declared Python modules, commands, and environment-variable presence without printing secret values. Read the workflow for complete version, R package, reference-data, API, and model/data-license requirements. Passing local checks is not execution validation.

Installation includes the selected skill's supporting files and licenses. Existing unmanaged directories are protected. Identical installations are unchanged; replacing a modified DNAKit installation requires explicit `--force`. Analysis dependencies are not installed by this command.

## Run an offline example

The upstream `genomic-coordinates` skill declares Python 3.11+:

```bash
dnakit skills install genomic-coordinates --destination ./analysis-skills
python3 analysis-skills/genomic-coordinates/scripts/convert_coords.py \
  --from bed --to gff chr1 999 1000 --format json
```

The output has `start=1000`, `end=1000`, and `length=1`: BED is 0-based half-open, while GFF is 1-based inclusive. Converting conventions is distinct from mapping between reference genome builds.

## Use MCP

After following the [MCP tutorial](mcp.md), an Agent can call:

```text
search_dnakit_skills(query="variant annotation", limit=5)
read_dnakit_skill(skill_id="variant-annotation")
```

Use a relative path from the response's `files` list as `resource` to read a supporting document. Continue with `next_offset` when provided. Reading text through MCP does not install scripts on the client; install the selected folder in the execution environment when scripts are needed.

The full catalog now has DNAKit integration guides and maintained adaptations.
Collection adapters and bundled scripts use the [unified DNAKit APIs](#unified-skill-api),
while actual scientific software and models retain their dependencies.
The `sequence-retrieval` workflow uses native NCBI/ENA adapters. FASTA, NCBI GenBank and ENA EMBL use the base package;
cross-format conversion requires `dnakit[bio-formats]`. See
[sequence search](../api/features/22_database_query_sequences.md) and
[sequence downloads](../api/features/15_download.md#accession-sequences).

```bash
dnakit skills search "sequence retrieval" --limit 5
dnakit skills install sequence-retrieval --agent codex
dnakit skills read sequence-retrieval --file examples.md
```

The old `tooluniverse-sequence-retrieval` and `public-db-lookup` identifiers remain
aliases. Existing client installation folders are not automatically moved or
deleted; remove obsolete copies when appropriate to avoid duplicate discovery.
Use each skill's documentation and `skills info` runtime metadata to check its status.

## Unified Skill APIs {#unified-skill-api}

Every bundled skill now includes `DNAKIT.md` with its DNAKit APIs, packaged script entrypoints and retained scientific dependencies. Use `dnakit skills scripts ID` and `dnakit skills run ID --script FILE -- <arguments>` without an upstream checkout.

For direct function usage, open the relevant tutorial: [queries](../api/features/22_database_query_discovery.md), [downloads](../api/features/15_download.md), [variant annotation](../api/features/35_variant_annotation.md), [variant effect scores](../api/features/36_variant_effect_scores.md), [custom model training](../api/features/23_prediction_custom_training.md), [regulatory tracks](../api/features/23_prediction_regulatory_tracks.md), [protein structure](../api/features/23_prediction_protein_features.md), or [workflows and reports](../api/features/14_engineering.md). Each guide includes parameters, returns and examples.

DNAKit provides search, installation and execution across the skill catalog. Collection-specific database adapters and shared workflow code are bundled in DNAKit with their original attribution. Scientific programs, libraries and trained models, including samtools, Biopython and AlphaGenome, keep their original implementations and requirements.

Collection prefixes are removed selectively; old names remain aliases. For example, `tooluniverse-sequence-retrieval` becomes `sequence-retrieval`, while `virulencefinder` and `alphagenome-atlas-website-links` keep their names.

### Database functions

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

### Bundled scripts

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

### ENCODE, composite workflows and reports

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

### MCP and validation

With `dnakit[agent]` and the [MCP configuration](mcp.md), agents discover these functions through `search_dnakit_tools`, inspect them with `describe_dnakit_tool`, and execute them with `call_dnakit_tool`. File writes, external processes and model execution retain their corresponding operation classifications. The initial context does not contain every function definition.

Read each skill's `DNAKIT.md` and inspect `dnakit skills info ID --json`. `registry/migration.json` records integration details; `registry/adaptations.json` binds maintained templates to source revisions. File checks, signature checks, unit tests, live service calls and complete scientific workflow validation are distinct. Unexecuted workflows remain `not_run`.

For ChIP-Atlas, `ch_ipatlas_get_peak_data(gene=...)` resolves a GRCh38 locus and returns a track retrieval plan with `status="requires_data"`. Supply the downloaded, verified BED as `peak_file` to compute overlap. The inherited enrichment adapter returns submission instructions, marked `requires_submission`, rather than computed statistics. A track URL or peak overlap does not establish motif disruption. These distinctions follow the [official ChIP-Atlas API](https://chip-atlas.org/agents).

Original authors, licenses and file hashes are retained. Coordinate-file retrieval follows [RCSB documentation](https://www.rcsb.org/docs/programmatic-access/file-download-services), and bounded Census queries follow the [CELLxGENE API](https://chanzuckerberg.github.io/cellxgene-census/_autosummary/cellxgene_census.get_anndata.html).
