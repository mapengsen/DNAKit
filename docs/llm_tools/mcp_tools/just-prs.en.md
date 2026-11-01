# just-prs MCP (`just-prs`)

DNAKit's MCP tool index keeps **16** DNA/gene tools of this server (20 listed, read on 2026-10-16). They run once you install the server; DNAKit does not bundle them. See [calling a tool that was found](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found).

| Item | Value |
| --- | --- |
| Source | [dna-seq/just-prs-mcp @ `a231ce8`](https://github.com/dna-seq/just-prs-mcp/tree/a231ce82970b8fabdc2e163b3e77661b140cb22d) |
| License | MIT |
| Command | `just-prs-mcp stdio` |
| Transports | stdio |
| Installation | pip install just-prs-mcp (PyPI, version 0.5.0 at verification) into a Python &gt;=3.13 environment (the optional 'reference' extra adds pgenlib on non-Windows); default mode is essentials and no API key is needed. Optionally set PRS_MCP_MODE (essentials or extended), PRS_MCP_CACHE_DIR, PRS_MCP_DEFAULT_GENOME_BUILD, PRS_MCP_DEFAULT_PANEL, PRS_MCP_DUCKDB_MEMORY_LIMIT. |
| Effects | Queries the PGS Catalog and downloads scoring files, reference data and sample genomes (Zenodo) into a local cache, reads local VCF or array files, writes normalized Parquet files, job state, logs and trait HTML reports to the cache directory; extended mode adds bulk downloads and a HuggingFace upload tool that needs an HF token, which this entry deliberately does not expose. |

| Tool | Category | Summary |
| --- | --- | --- |
| [`absolute_risk`](#absolute_risk) | `variant` | Estimates absolute disease risk for a PGS score from a PRS z-score, sex and trait prevalence/effect-size data. |
| [`assess_quality`](#assess_quality) | `variant` | Classifies the quality of a PRS result from variant match rate, AUROC and percentile, returning a label and interpretation. |
| [`best_performance`](#best_performance) | `variant` | Looks up the best published evaluation of a PGS score, returning effect sizes (OR/HR/beta) and AUROC/C-index. |
| [`compare_genomes`](#compare_genomes) | `variant` | Compares saved by-trait PRS reports across two or more genomes, ranking by median in-scope percentile per trait. |
| [`compute_prs`](#compute_prs) | `variant` | Computes a polygenic risk score for one genome VCF against one PGS Catalog model, with optional reference restoration and ancestry. |
| [`compute_prs_batch`](#compute_prs_batch) | `variant` | Starts a job computing PRS for one genome VCF against many PGS IDs using a DuckDB engine and cached normalization. |
| [`compute_prs_by_trait`](#compute_prs_by_trait) | `variant` | Computes all PGS scores linked to an EFO/MONDO trait for one or more genome VCFs and saves an HTML trait report. |
| [`download_sample_genome`](#download_sample_genome) | `variant` | Downloads a public sample WGS VCF from Zenodo (or any Zenodo record) and optionally normalizes it to Parquet. |
| [`normalize_vcf`](#normalize_vcf) | `variant` | Normalizes a VCF into a quality-filtered genotype Parquet (chr prefix stripped, rsid, GT-derived genotype) as a background job. |
| [`percentile`](#percentile) | `variant` | Estimates the population percentile and z-score of a computed PRS using reference panels, theoretical distribution or AUROC fallback. |
| [`plot_trait_panel`](#plot_trait_panel) | `variant` | Writes an interactive HTML trait report (Vega-Lite panel, model table) from saved PRS result files and returns its link. |
| [`score_info`](#score_info) | `variant` | Returns cleaned metadata for a single PGS Catalog score by its PGS ID. |
| [`search_scores`](#search_scores) | `variant` | Searches PGS Catalog polygenic scores by free text over ID, name and trait, optionally filtered by genome build. |
| [`search_traits`](#search_traits) | `variant` | Searches PGS Catalog traits by term and returns counts or lists of directly associated and child-trait PGS IDs. |
| [`trait_info`](#trait_info) | `variant` | Fetches a PGS Catalog trait by EFO or MONDO ID along with its associated PGS IDs. |
| [`vcf_metainfo`](#vcf_metainfo) | `variant` | Detects a genome VCF's build, WGS vs array vs gVCF input type and inferred super-population ancestry. |

## `absolute_risk` {#absolute_risk}

Estimates absolute disease risk for a PGS score from a PRS z-score, sex and trait prevalence/effect-size data.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `pgs_id`, `z_score`
- **Keywords:** polygenic risk score, pgs catalog, z-score, prevalence, lifetime risk, prs
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `pgs_id` | string | yes |  |
| `z_score` | number | yes |  |
| `sex` | string \| null | no | Default `None` |

??? note "Upstream description"

    Estimate absolute disease risk from a PRS z-score and population prevalence.

    Joins the score's trait to prevalence + effect-size data. ``z_score`` is
    the PRS in standard deviations from the population mean. Raises if the
    required prevalence / effect-size data is unavailable for this score.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'just-prs',
    'absolute_risk',
    {'pgs_id': '...', 'z_score': 0.0},
)
```

## `assess_quality` {#assess_quality}

Classifies the quality of a PRS result from variant match rate, AUROC and percentile, returning a label and interpretation.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `match_rate`
- **Keywords:** polygenic risk score, match rate, auroc, prs quality, percentile, prs
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `match_rate` | number | yes |  |
| `auroc` | number \| null | no | Default `None` |
| `percentile` | number \| null | no | Default `None` |
| `percentile_method` | string \| null | no | Default `None` |
| `reliable` | boolean | no | Default `True` |
| `caveat` | string | no | Default `` |

??? note "Upstream description"

    Classify and interpret a PRS result's quality (pure logic — no I/O).

    ``match_rate`` is the fraction of scoring variants matched (0-1). Returns
    a quality label/color and a human-readable interpretation combining match
    rate, AUROC, and (optionally) the result percentile. Pass
    ``percentile_method`` / ``reliable`` / ``caveat`` from the ``percentile``
    tool so the summary describes how the percentile was actually derived and
    echoes any low-coverage caveat.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'just-prs',
    'assess_quality',
    {'match_rate': 0.0},
)
```

## `best_performance` {#best_performance}

Looks up the best published evaluation of a PGS score, returning effect sizes (OR/HR/beta) and AUROC/C-index.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `pgs_id`
- **Keywords:** pgs catalog, polygenic score, auroc, odds ratio, effect size, evaluation
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `pgs_id` | string | yes |  |

??? note "Upstream description"

    Look up the best evaluation performance for a score (largest sample, EUR-preferred).

    Returns parsed effect sizes (OR/HR/Beta) and classification metrics
    (AUROC/C-index) plus pre-formatted display strings.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'just-prs',
    'best_performance',
    {'pgs_id': '...'},
)
```

## `compare_genomes` {#compare_genomes}

Compares saved by-trait PRS reports across two or more genomes, ranking by median in-scope percentile per trait.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `result_paths`
- **Keywords:** polygenic risk score, genome comparison, family, percentile, trait panel, prs
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `result_paths` | array of string | yes |  |
| `genome_labels` | array \| null | no | Default `None` |
| `model_scope` | string | no | Default `usable` |
| `include_prompt` | boolean | no | Default `True` |
| `assistant` | string | no | Default `other` |
| `ancestry` | string | no | Default `EUR` |

??? note "Upstream description"

    Compare saved by-trait PRS reports across two or more genomes.

    ``result_paths`` — the ``result_paths`` list from ``compute_prs_by_trait``
    (or any saved report JSONs). At least two paths are required. Files are
    grouped by trait, so a family comparison across several traits is one call.

    Rankings use just-prs ``summarize_trait_rows``: the **median** in-scope
    percentile is the headline (not a single "best" model). High percentile is
    not automatically good or bad — interpret direction from the trait.

    When every path is the same trait, ``include_prompt=True`` (default) attaches
    the same multi-sample LLM prompt as ``prs prompt`` / the UI Ask-AI buttons.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'just-prs',
    'compare_genomes',
    {'result_paths': []},
)
```

## `compute_prs` {#compute_prs}

Computes a polygenic risk score for one genome VCF against one PGS Catalog model, with optional reference restoration and ancestry.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `vcf_path`, `pgs_id`
- **Keywords:** polygenic risk score, vcf, genotype, pgs catalog, reference allele, wgs, ancestry
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `vcf_path` | string | yes |  |
| `pgs_id` | string | yes |  |
| `genome_build` | string \| null | no | Default `None` |
| `genotypes_path` | string \| null | no | Default `None` |
| `attach_performance` | boolean | no | Default `False` |
| `reference_restoration` | string | no | Default `auto` |
| `genotype_input_mode` | string | no | Default `auto` |
| `infer_ancestry` | boolean | no | Default `False` |

??? note "Upstream description"

    Compute a polygenic risk score for one genome against one PGS model.

    Pass the raw VCF directly; a separate ``normalize_vcf`` call is not needed.
    If this VCF has a fresh normalized cache, it is reused automatically;
    ``genotypes_path`` is the explicit Parquet fast path. For a trait panel —
    one genome or several — use
    ``compute_prs_by_trait(samples=\[...\])``, then
    ``compare_genomes`` / ``build_prs_prompt`` / ``plot_trait_panel``.

    **Reference restoration (F15/F22).** ``reference_restoration`` fills a scoring
    variant's missing reference allele from the precomputed reference-allele
    universe (pulled from HuggingFace on first use) so a locus the sample is
    confidently hom-ref at scores as dose-0 instead of being dropped as
    unscorable-absent — recovering the ~50% coverage gap on genome-wide WGS
    scores. Default ``"auto"`` detects the input: whole-universe restoration for
    WGS, the matching chip scope for a genotyping array, off otherwise. Pass
    ``"off"``, ``"wgs"``, or a chip id (e.g. ``"gsa_v3"``) to force a scope. The
    result's ``variants_ref_resolved_panel`` / ``variants_ref_resolved_fasta``
    counters show how many loci were recovered.

    Set ``infer_ancestry=True`` to also infer and attach the sample's genetic
    ancestry (``sample_ancestry``). ``genotype_input_mode`` (``"auto"`` /
    ``"variant_only"`` / ``"all_sites"``) controls how absent loci are read.

    Set ``attach_performance=True`` to embed the score's best published
    performance (effect sizes, AUROC/C-index, evaluation ancestry) on the
    result in the same call. The result also carries ``detected_genome_build`` /
    ``build_mismatch`` from the VCF (F4).

    **Recommended follow-up:** ``percentile`` then ``absolute_risk`` (disease
    traits), or ``build_prs_prompt(kind='score')`` for a reusable LLM write-up.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'just-prs',
    'compute_prs',
    {'vcf_path': '...', 'pgs_id': '...'},
)
```

## `compute_prs_batch` {#compute_prs_batch}

Starts a job computing PRS for one genome VCF against many PGS IDs using a DuckDB engine and cached normalization.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `vcf_path`, `pgs_ids`
- **Keywords:** polygenic risk score, vcf, batch scoring, pgs catalog, genotype, duckdb
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `vcf_path` | string | yes |  |
| `pgs_ids` | array of string | yes |  |
| `genome_build` | string \| null | no | Default `None` |
| `attach_performance` | boolean | no | Default `False` |
| `reference_restoration` | string | no | Default `auto` |
| `genotype_input_mode` | string | no | Default `auto` |

??? note "Upstream description"

    Compute PRS for one genome against many PGS models.

    Pass the raw VCF directly; a separate normalization call is not required.
    The tool transparently creates or reuses a collision-safe normalized cache,
    then uses the memory-safe DuckDB engine with spill-to-disk. It reuses the
    scoring caches across scores. This is still **one genome × many PGS IDs**.
    For a trait panel across people, use ``compute_prs_by_trait(samples=\[...\])``
    then ``compare_genomes`` / ``build_prs_prompt``.
    Set ``attach_performance=True`` to embed each score's best published
    performance on its result in the same pass.

    ``reference_restoration`` (default ``"auto"``) recovers absent-hom-ref
    coverage on WGS/array input — see ``compute_prs``; the reference-allele
    universe is resolved once and reused across every score. ``genotype_input_mode``
    (``"auto"`` / ``"variant_only"`` / ``"all_sites"``) controls how absent loci
    are read.

    Returns a ``PRSJob`` immediately. When done, ``result`` is the list of
    successful ``PRSResult`` objects. If status is queued or running, poll
    ``prs_job_status(job_id, wait_seconds=25)``.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'just-prs',
    'compute_prs_batch',
    {'vcf_path': '...', 'pgs_ids': []},
)
```

## `compute_prs_by_trait` {#compute_prs_by_trait}

Computes all PGS scores linked to an EFO/MONDO trait for one or more genome VCFs and saves an HTML trait report.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `trait_id`, `samples`
- **Keywords:** polygenic risk score, efo, mondo, trait, vcf, percentile, ancestry

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `trait_id` | string | yes |  |
| `samples` | array of string | yes |  |
| `genome_build` | string \| null | no | Default `None` |
| `include_children` | boolean | no | Default `False` |
| `limit` | integer \| null | no | Default `None` |
| `interpret` | boolean | no | Default `True` |
| `include_prompt` | boolean | no | Default `False` |
| `superpopulation` | string | no | Default `EUR` |
| `panel` | string \| null | no | Default `None` |
| `top_n` | integer \| null | no | Default `None` |
| `profile` | string | no | Default `all` |
| `min_match_rate` | number \| null | no | Default `None` |
| `min_auroc` | number \| null | no | Default `None` |
| `build` | string \| null | no | Default `None` |
| `ancestry` | string \| null | no | Default `None` |
| `reference_restoration` | string | no | Default `auto` |
| `genotype_input_mode` | string | no | Default `auto` |

??? note "Upstream description"

    Compute the PRS scores associated with a trait for one or more genomes.

    REQUIRED REPORT DELIVERY: when html_link is set, your reply MUST include it as a markdown link (for example \[Open the interactive trait report\](&lt;html_link&gt;)) and explicitly tell the user to open it. This is the primary user-facing result, not an implementation detail. Do not skip the link, merely say the report was saved, fetch it, or paste the HTML. html_link is an HTTP URL when available and a file:// URI on local stdio.

    Pass a raw .vcf/.vcf.gz path directly. Do NOT call normalize_vcf first as a prerequisite: multi-score tools transparently create or reuse a safe normalized cache, and single-score tools read the VCF directly. normalize_vcf is only an optional explicit step for custom filters. If list_genomes already shows a normalized Parquet, it can be passed directly.

    ``trait_id`` may be an EFO or MONDO identifier. ``samples`` is always a
    list (one genome is fine): each entry is ``Label=path`` or a bare path;
    ``.parquet`` is a normalized genotype source. Every genome is auto-saved;
    the return carries ``reports``, ``result_paths``, and (when models have
    percentiles) the primary ``html_link`` plus transport-specific
    ``html_url`` / ``html_path`` — the just-prs trait HTML report (same page
    as ``prs plot trait -o *.html``).
    Pass ``result_paths`` to ``compare_genomes`` (2+), ``plot_trait_panel``,
    or ``build_prs_prompt``.

    **Profile (curation):** ``profile="all"`` (default) returns the raw
    associated panel. ``profile="curated"`` is the interpreted *shortlist* —
    it drops toy scores (&lt;10 variants), scores with no performance evidence,
    scores below the C_wt coverage floor, and de-dups score families.

    ``interpret=True`` (default) attaches percentiles/quality per model.
    Set ``interpret=False`` for raw scores only — the report then says so.
    ``reference_restoration`` defaults to ``"auto"``. ``superpopulation``
    may be a 1000G code or ``"auto"``. Set ``include_prompt=True`` to also
    attach the just-prs LLM prompt covering every genome in this call.

    This tool returns a PRSJob immediately. If status is queued or running, call prs_job_status(job_id, wait_seconds=25) until done or failed. Do not start a second identical compute — the server reuses the in-flight job. When done, html_link and result (TraitPRSReports) are on the job. Native multi-score work is isolated in a child process so a DuckDB/Polars crash fails the batch without killing the MCP transport.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'just-prs',
    'compute_prs_by_trait',
    {'trait_id': '...', 'samples': []},
)
```

## `download_sample_genome` {#download_sample_genome}

Downloads a public sample WGS VCF from Zenodo (or any Zenodo record) and optionally normalizes it to Parquet.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** none
- **Keywords:** wgs, vcf, zenodo, sample genome, personal genome, download

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sample` | string | no | Default `anton` |
| `output_dir` | string \| null | no | Default `None` |
| `record_url` | string \| null | no | Default `None` |
| `filename` | string \| null | no | Default `None` |
| `auto_normalize` | boolean | no | Default `True` |
| `force` | boolean | no | Default `False` |

??? note "Upstream description"

    Download a public sample WGS VCF from Zenodo to try PRS without your own data.

    Two whole-genome sequencing (WGS) datasets open-sourced by the
    just-dna-lite project are pre-configured — download both when you want a
    multi-person comparison:

    - ``sample="anton"`` — Anton Kulaga's genome (~482 MB, CC0 public domain,
    Zenodo record 18370498, file: antonkulaga.vcf).
    - ``sample="livia"`` — Livia Zaharia's genome (~349 MB, CC-BY-4.0,
    Zenodo record 19487816, file: SIMHIFQTILQ.hard-filtered.vcf.gz).

    Pass ``record_url`` (e.g. 'https://zenodo.org/records/18370498') to fetch
    any other Zenodo record, and ``filename`` to pick a specific file when a
    record has several.

    The downloaded VCF lands under ``&lt;cache_dir&gt;/samples/`` (or ``output_dir``)
    and is a drop-in path for ``normalize_vcf`` / ``compute_prs`` /
    ``compute_prs_by_trait`` (including ``samples=\["Anton=&lt;path&gt;", ...\]``).

    ``auto_normalize`` defaults to ``True``: the download is normalized to a
    reusable Parquet in the same call, so ``data`` carries both ``path`` (the
    raw VCF) and ``normalized_path`` — a one-call, compute-ready genotype
    source with no separate ``normalize_vcf`` round-trip. (Normalization is
    idempotent, so a re-download of an already-staged sample is cheap.) Pass
    ``auto_normalize=False`` to fetch the raw VCF only.

    Use ``list_genomes`` to see which genomes have already been downloaded
    and/or normalized.

    Idempotent: if the target VCF already exists with the size Zenodo
    reports, the ~hundreds-of-MB download is skipped and the cached file is
    reused; likewise a present Parquet skips re-normalization. ``data`` echoes
    ``reused_cache`` (download skipped) and ``downloaded_bytes`` (bytes
    actually transferred, 0 on a cache hit) so the caller can tell a cache
    hit from a fresh fetch. Pass ``force=True`` to re-download/re-normalize
    regardless.

    Returns a ``PRSJob`` immediately. When done, ``result`` is an ``OpResult``
    whose ``data`` carries the local ``path``. Poll ``prs_job_status`` if
    status is still running — the download is hundreds of MB.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'just-prs',
    'download_sample_genome',
    {},
)
```

## `normalize_vcf` {#normalize_vcf}

Normalizes a VCF into a quality-filtered genotype Parquet (chr prefix stripped, rsid, GT-derived genotype) as a background job.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `vcf_path`
- **Keywords:** vcf, genotype, parquet, quality filter, rsid, normalization, depth

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `vcf_path` | string | yes |  |
| `output_path` | string \| null | no | Default `None` |
| `pass_filters` | array \| null | no | Default `None` |
| `min_depth` | integer \| null | no | Default `None` |
| `min_qual` | number \| null | no | Default `None` |
| `sex` | string \| null | no | Default `None` |
| `genome_build` | string \| null | no | Default `None` |
| `force` | boolean | no | Default `False` |

??? note "Upstream description"

    Normalize a VCF to a quality-filtered genotype Parquet.

    This tool is **optional preprocessing**, not a prerequisite for scoring:
    ``compute_prs``, ``compute_prs_batch``, and ``compute_prs_by_trait`` accept
    raw VCF paths directly. Use normalization when you intentionally want a
    reusable Parquet for many later analyses or custom quality filters. Before
    normalizing, ``list_genomes`` can show whether a cached Parquet already exists.

    Strips the chr prefix, renames id→rsid, computes genotype from GT, applies
    optional quality filters (FILTER allow-list, min DP, min QUAL), and writes
    zstd-compressed Parquet. The output is a drop-in genotype source for
    ``compute_prs`` / ``compute_prs_batch`` (``genotypes_path``) and for
    ``compute_prs_by_trait`` (a ``samples`` entry ending in ``.parquet``), so
    each VCF is normalized once and reused. Normalize every genome you plan
    to compare; scoring still happens per genome.

    Idempotent: if the target Parquet already exists it is reused and the
    (slow) normalization is skipped — ``reused_cache=True`` flags the hit.
    Custom filters (``pass_filters`` / ``min_depth`` / ``min_qual`` / ``sex``)
    always re-run, since the cached Parquet may not reflect them. Pass
    ``force=True`` to re-normalize unconditionally.

    Normalization is the slow step (seconds to minutes depending on VCF size).
    Returns a ``PRSJob`` immediately; poll ``prs_job_status`` when status is
    running. When done, ``result`` is a ``NormalizeResult``.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'just-prs',
    'normalize_vcf',
    {'vcf_path': '...'},
)
```

## `percentile` {#percentile}

Estimates the population percentile and z-score of a computed PRS using reference panels, theoretical distribution or AUROC fallback.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `prs_score`, `pgs_id`
- **Keywords:** polygenic risk score, percentile, 1000 genomes, superpopulation, z-score, reference panel
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `prs_score` | number | yes |  |
| `pgs_id` | string | yes |  |
| `superpopulation` | string | no | Default `EUR` |
| `panel` | string \| null | no | Default `None` |
| `weight_mass_coverage` | number \| null | no | Default `None` |

??? note "Upstream description"

    Estimate the population percentile (0-100) for a computed PRS value.

    Uses the 3-tier fallback: precomputed reference-panel distributions
    (best), then a theoretical distribution, then an AUROC approximation.
    ``superpopulation`` is a 1000G code (AFR/AMR/EAS/EUR/SAS) and should
    match **that genome's** ancestry — when comparing people, call this
    once per score×genome (or rely on ``compute_prs_by_trait(interpret=True)``).
    Pass ``weight_mass_coverage`` (C_wt) from ``compute_prs`` so a deflated
    low-coverage percentile is flagged ``reliable=False`` with a caveat
    instead of presented as authoritative. Also returns the true z-score and
    reference mean/std used, so absolute risk can be computed without
    inverting the percentile.

    **Important next step:** for disease traits, feed the returned
    ``z_score`` directly into ``absolute_risk`` to get the concrete
    lifetime probability and risk ratio vs the population average. This is
    more informative than the percentile alone.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'just-prs',
    'percentile',
    {'prs_score': 0.0, 'pgs_id': '...'},
)
```

## `plot_trait_panel` {#plot_trait_panel}

Writes an interactive HTML trait report (Vega-Lite panel, model table) from saved PRS result files and returns its link.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `result_paths`
- **Keywords:** polygenic risk score, html report, visualization, trait panel, vega-lite

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `result_paths` | array of string | yes |  |
| `output_path` | string \| null | no | Default `None` |
| `genome_labels` | array \| null | no | Default `None` |
| `ancestry` | string | no | Default `EUR` |
| `model_scope` | string | no | Default `usable` |
| `panel` | string \| null | no | Default `None` |

??? note "Upstream description"

    Write the just-prs trait HTML report and return a URL or file path.

    REQUIRED REPORT DELIVERY: when html_link is set, your reply MUST include it as a markdown link (for example \[Open the interactive trait report\](&lt;html_link&gt;)) and explicitly tell the user to open it. This is the primary user-facing result, not an implementation detail. Do not skip the link, merely say the report was saved, fetch it, or paste the HTML. html_link is an HTTP URL when available and a file:// URI on local stdio.

    This is the MCP equivalent of ``prs plot trait ... -o intel_o_family.html``:
    ``just_prs.viz.plot_trait_scores`` + ``save_trait_report``. The page has
    the Vega-Lite panel, median cards, model table, and Ask-AI buttons.
    ``result_paths`` is always a list (one genome is fine; two or more
    overlay people). On Streamable HTTP ``html_url`` is served from this
    same origin at ``/reports/&lt;file&gt;`` with no extra config.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'just-prs',
    'plot_trait_panel',
    {'result_paths': []},
)
```

## `score_info` {#score_info}

Returns cleaned metadata for a single PGS Catalog score by its PGS ID.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `pgs_id`
- **Keywords:** pgs catalog, polygenic score, metadata, pgs id
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `pgs_id` | string | yes |  |

??? note "Upstream description"

    Get cleaned metadata for a single PGS score by its ID (e.g. 'PGS000001').

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'just-prs',
    'score_info',
    {'pgs_id': '...'},
)
```

## `search_scores` {#search_scores}

Searches PGS Catalog polygenic scores by free text over ID, name and trait, optionally filtered by genome build.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `query`
- **Keywords:** pgs catalog, polygenic score, trait search, grch38, grch37
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string | yes |  |
| `genome_build` | string \| null | no | Default `None` |
| `limit` | integer | no | Default `25` |

??? note "Upstream description"

    Search the PGS Catalog for polygenic scores by free text.

    Case-insensitive substring match across PGS ID, score name, reported
    trait, and EFO trait. Optionally filter to a genome build (GRCh37 /
    GRCh38, harmonized cross-build scores included). Returns up to ``limit``
    matches with their key metadata.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'just-prs',
    'search_scores',
    {'query': '...'},
)
```

## `search_traits` {#search_traits}

Searches PGS Catalog traits by term and returns counts or lists of directly associated and child-trait PGS IDs.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `term`
- **Keywords:** pgs catalog, trait, efo, polygenic score, synonyms
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `term` | string | yes |  |
| `limit` | integer | no | Default `25` |
| `include_pgs_ids` | boolean | no | Default `False` |

??? note "Upstream description"

    Search the PGS Catalog REST API for traits by term.

    Upstream matching is exact-substring over labels and synonyms, so this
    wrapper retries a few punctuation/order variants when the first query is
    empty. By default, results include counts of directly associated PGS IDs
    and child-trait PGS IDs; set ``include_pgs_ids`` for the full arrays.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'just-prs',
    'search_traits',
    {'term': '...'},
)
```

## `trait_info` {#trait_info}

Fetches a PGS Catalog trait by EFO or MONDO ID along with its associated PGS IDs.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** none
- **Keywords:** pgs catalog, efo, mondo, trait, polygenic score
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `trait_id` | string \| null | no | Default `None` |
| `efo_id` | string \| null | no | Default `None` |

??? note "Upstream description"

    Fetch a trait by ontology ID (EFO or MONDO) with its associated PGS IDs.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'just-prs',
    'trait_info',
    {},
)
```

## `vcf_metainfo` {#vcf_metainfo}

Detects a genome VCF's build, WGS vs array vs gVCF input type and inferred super-population ancestry.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** none
- **Keywords:** vcf, genome build, genotyping array, gvcf, ancestry, superpopulation
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `vcf_path` | string \| null | no | Default `None` |
| `genotypes_path` | string \| null | no | Default `None` |
| `genome_build` | string \| null | no | Default `None` |
| `panel` | string \| null | no | Default `None` |
| `resolution` | string | no | Default `superpop` |
| `infer_population` | boolean | no | Default `True` |

??? note "Upstream description"

    Detect a genome's build, input type (WGS/array/gVCF), and genetic ancestry.

    One call that pairs **WGS-vs-array/gVCF detection** with **population
    (super-population) inference**, so before scoring you already know: the
    genome build (read from the VCF header), whether the input is whole-genome
    sequencing, a consumer genotyping array (and which chip), or a gVCF /
    all-sites callset, and the sample's inferred super-population.

    Pass ``genotypes_path`` (a normalized Parquet from ``normalize_vcf`` /
    ``download_sample_genome``) for the fast path, or a raw ``vcf_path``.
    Genome-build detection reads the VCF header, so it is only available on the
    ``vcf_path`` branch.

    The result carries ready-to-use recommendations:
    ``recommended_reference_restoration`` (feed to ``compute_prs`` /
    ``compute_prs_batch`` / ``compute_prs_by_trait`` as ``reference_restoration``
    — though those tools default to ``"auto"`` and derive the same thing) and
    ``recommended_superpopulation`` (feed as ``superpopulation``). Set
    ``infer_population=False`` to skip ancestry inference (it pulls a small model
    from HuggingFace on first use).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'just-prs',
    'vcf_metainfo',
    {},
)
```
