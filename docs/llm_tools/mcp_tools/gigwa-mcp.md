# Gigwa MCP (`gigwa-mcp`)

DNAKit 的 MCP 工具索引收录了这个服务器的 **27** 个 DNA/基因相关工具（服务器共列出 35 个，2026-10-16 读取）。工具在用户自己安装服务器后运行，DNAKit 本身不包含它们。调用方式见 [调用索引中的工具](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found)。下面的工具说明是上游英文原文。

| 项目 | 内容 |
| --- | --- |
| 源码 | [gkanogiannis/Gigwa-MCP @ `d6ccdd1`](https://github.com/gkanogiannis/Gigwa-MCP/tree/d6ccdd1a5eb19b8a0cc224b96ad459868d9c4aaf) |
| 许可证 | Apache-2.0 |
| 启动命令 | `gigwa-mcp` |
| 传输方式 | stdio |
| 安装 | Install gigwa-mcp separately into a Python &gt;=3.10 environment; configure GIGWA_URL and optional paired GIGWA_USER/GIGWA_PASS. |
| 必需环境变量 | `GIGWA_URL` |
| 副作用 | Tools may upload data, modify Gigwa, or write local result files; inspect each tool before invocation. |

| 工具 | 类别 | 简介 |
| --- | --- | --- |
| [`audit_import_quality`](#audit_import_quality) | `breeding_population` | Scans Gigwa variant sets for genotype-encoding artifacts by sampling genotypes and flagging abnormal heterozygosity, completeness or monomorphism. |
| [`count_variants`](#count_variants) | `variant` | Counts variants in a Gigwa variant set server-side, filtered by genomic region, MAF and missing-data fraction. |
| [`diversity_by_group`](#diversity_by_group) | `breeding_population` | Computes per-group He, Ho, Fis, MAF, percent polymorphic and rarefied allelic richness from genotypes. |
| [`diversity_core_collection`](#diversity_core_collection) | `breeding_population` | Selects a core collection of accessions that maximises captured marker-allele diversity (Core Hunter style). |
| [`diversity_fst`](#diversity_fst) | `breeding_population` | Computes pairwise Weir and Cockerham Fst between user-defined groups of samples from genotypes. |
| [`diversity_kinship`](#diversity_kinship) | `breeding_population` | Computes a VanRaden genomic relationship matrix from alt-allele dosage and reports most-related pairs. |
| [`diversity_pca`](#diversity_pca) | `breeding_population` | Runs PCA on the alt-allele dosage matrix to describe population structure and flag outlier samples. |
| [`diversity_structure`](#diversity_structure) | `breeding_population` | Clusters samples into genetic groups via PCA and K-means, choosing K by pseudo-F. |
| [`diversity_summary`](#diversity_summary) | `breeding_population` | Computes per-marker MAF, expected/observed heterozygosity and PIC plus dataset means. |
| [`diversity_tree`](#diversity_tree) | `breeding_population` | Builds a UPGMA Newick tree of accessions from pairwise IBS allele-sharing distance. |
| [`export_genotypes`](#export_genotypes) | `breeding_population` | Exports a Gigwa variant set or filtered subset to VCF, PLINK, Flapjack or other genotype formats. |
| [`filter_individuals_by_metadata`](#filter_individuals_by_metadata) | `breeding_population` | Selects individuals whose stored germplasm metadata match field/value filters, returning IDs for genotype export. |
| [`get_germplasm_metadata`](#get_germplasm_metadata) | `breeding_population` | Fetches server-stored per-individual germplasm attributes for a Gigwa database and writes them to CSV. |
| [`import_dartseq`](#import_dartseq) | `breeding_population` | Converts DArTseq SNP/Silico-DArT xlsx reports to VCF with correct genotype calling and imports them into a Gigwa database. |
| [`import_metadata`](#import_metadata) | `breeding_population` | Imports a per-individual metadata TSV into an existing Gigwa database as searchable attributes. |
| [`import_vcf`](#import_vcf) | `breeding_population` | Imports a VCF file into a Gigwa database/project/run, optionally replacing existing project data. |
| [`list_content`](#list_content) | `breeding_population` | Lists Gigwa databases, projects and runs with their exact BrAPI variantSetDbIds. |
| [`list_metadata_values`](#list_metadata_values) | `breeding_population` | Lists individual-metadata fields available for a Gigwa database and their distinct values. |
| [`list_sequences`](#list_sequences) | `breeding_population` | Lists reference sequences (chromosomes/contigs) present in a Gigwa variant set. |
| [`list_variant_sets`](#list_variant_sets) | `breeding_population` | Lists every Gigwa variant set (run) with its BrAPI variantSetDbId and variant/callset counts. |
| [`map_dartseq_to_reference`](#map_dartseq_to_reference) | `breeding_population` | Aligns DArTseq SNP marker tag sequences to a reference genome with minimap2 to infer chromosome, position and strand. |
| [`qc_call_rate`](#qc_call_rate) | `breeding_population` | Computes per-sample and per-marker call rate (missingness) for a variant set and flags low-quality ones. |
| [`qc_duplicate_accessions`](#qc_duplicate_accessions) | `breeding_population` | Detects duplicate or clonal accessions by pairwise IBS similarity and groups them into duplicate sets. |
| [`qc_heterozygosity`](#qc_heterozygosity) | `breeding_population` | Computes per-sample observed heterozygosity and flags outliers suggesting contamination or inbreeding. |
| [`qc_maf_filter`](#qc_maf_filter) | `breeding_population` | Reports per-marker MAF and missingness and counts markers that would be filtered, without modifying data. |
| [`search_callsets`](#search_callsets) | `breeding_population` | Retrieves BrAPI callsets for a run and writes per-sample metadata (names, IDs, attributes) to CSV. |
| [`search_variants`](#search_variants) | `variant` | Searches variants matching region, MAF and missing-data filters and writes their id, chrom, pos, ref, alt to CSV. |

## `audit_import_quality` {#audit_import_quality}

Scans Gigwa variant sets for genotype-encoding artifacts by sampling genotypes and flagging abnormal heterozygosity, completeness or monomorphism.

- **类别：**`breeding_population`（Breeding genotyping, markers, population genetics and strain resources）
- **必填参数：**无
- **关键词：**genotype qc, heterozygosity, allele matrix, brapi, import artifacts

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant_set_db_id` | string \| null | 否 | BrAPI variantSetDbId identifying the run (MODULE§project§run) -- copy the exact string from list_variant_sets / list_content, never assemble one by hand: the middle segment is a numeric project index, not the project's name, and a wrong guess fails with an opaque HTTP 500 rather than a clear error. 默认 `None` |
| `max_markers` | integer | 否 | Cap analysis to the first N markers in canonical Gigwa search order; omit to use all. 默认 `1000` |
| `max_samples` | integer | 否 | Cap the number of samples/callsets sampled (allelematrix path). 默认 `300` |
| `het_threshold` | number | 否 | Mean observed-heterozygosity above which a run is flagged BROKEN (mis-called heterozygotes). 默认 `0.6` |
| `complete_call_rate` | number | 否 | Call-rate above which a run is flagged as suspiciously complete (no missing data). 默认 `0.999` |
| `monomorphic_threshold` | number | 否 | Monomorphic-marker fraction above which a run is flagged for low informativeness. 默认 `0.9` |
| `output_dir` | string \| null | 否 | Directory for the output CSV(s) (default ./gigwa_results/&lt;module&gt;/). 默认 `None` |

??? note "上游说明（英文原文）"

    Scan a Gigwa instance for databases imported with genotype-encoding artifacts.

    With no ``variant_set_db_id`` this audits **every run on the instance**; pass one to
    audit a single variant set. For each run it pulls a bounded genotype sample (up to
    ``max_markers`` markers × ``max_samples`` callsets) via paged BrAPI
    ``search/allelematrix`` — cheap and constant-cost regardless of how large the variant
    set is, so it is safe to run across a whole production instance without exporting
    multi-GB VCFs. The aggregate genotype-class fractions it needs are estimated tightly
    from the sample (a *true* zero hom-alt class stays zero; a rare-but-real one shows up).
    It flags two import failure modes plus two weaker signals:

    - **BROKEN** — cohort mean Ho above ``het_threshold`` (DArT 2-row mis-call), or
    homozygous-alt genotypes far below their HWE expectation given the alt-allele
    frequency (lost hom-alt class; the HWE test avoids false positives on low-MAF /
    mostly-monomorphic panels where near-zero hom-alt is genuine).
    - **SUSPECT** — call rate above ``complete_call_rate`` (no missing data, often missing
    forced to 0/0), monomorphic fraction above ``monomorphic_threshold``, or AD/DP depth
    fields present but uniformly zero (a VCF synthesised from genotype calls with
    fabricated depth/likelihoods — the same converter often miscalls GT too).

    Writes ``import_quality_scan.csv`` (one row per run) under ``output_dir`` (default
    ``./gigwa_results/``) and returns a summary ranked worst-first. Read-only — it never
    modifies Gigwa.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gigwa-mcp',
    'audit_import_quality',
    {},
)
```

## `count_variants` {#count_variants}

Counts variants in a Gigwa variant set server-side, filtered by genomic region, MAF and missing-data fraction.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`variant_set_db_id`
- **关键词：**variant count, maf, region filter, brapi, missing data

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant_set_db_id` | string | 是 | BrAPI variantSetDbId identifying the run (MODULE§project§run) -- copy the exact string from list_variant_sets / list_content, never assemble one by hand: the middle segment is a numeric project index, not the project's name, and a wrong guess fails with an opaque HTTP 500 rather than a clear error. |
| `reference_name` | string \| null | 否 | Chromosome/contig name to restrict the search to (see list_sequences). 默认 `None` |
| `start` | integer \| null | 否 | Region start position, 1-based inclusive. 默认 `None` |
| `end` | integer \| null | 否 | Region end position, 1-based inclusive. 默认 `None` |
| `min_maf` | number \| null | 否 | Minimum minor-allele frequency (0-1). 默认 `None` |
| `max_maf` | number \| null | 否 | Maximum minor-allele frequency (0-1). 默认 `None` |
| `max_missing_data` | number \| null | 否 | Maximum per-variant missing-data fraction (0-1). 默认 `None` |

??? note "上游说明（英文原文）"

    Count variants matching filters, computed server-side (nothing is downloaded).

    Fast way to size a query before pulling data. Filter by genomic region
    (``reference_name`` + optional ``start``/``end``, from ``list_sequences``), minor-
    allele frequency (``min_maf``/``max_maf``) and/or ``max_missing_data`` (0–1 fraction).
    With no filters this returns the total variant count of the set. ``variant_set_db_id``
    is a BrAPI variantSetDbId (from ``list_variant_sets`` / ``list_content``).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gigwa-mcp',
    'count_variants',
    {'variant_set_db_id': '...'},
)
```

## `diversity_by_group` {#diversity_by_group}

Computes per-group He, Ho, Fis, MAF, percent polymorphic and rarefied allelic richness from genotypes.

- **类别：**`breeding_population`（Breeding genotyping, markers, population genetics and strain resources）
- **必填参数：**`variant_set_db_id`
- **关键词：**heterozygosity, fis, allelic richness, population genetics, nei

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant_set_db_id` | string | 是 | BrAPI variantSetDbId identifying the run (MODULE§project§run) -- copy the exact string from list_variant_sets / list_content, never assemble one by hand: the middle segment is a numeric project index, not the project's name, and a wrong guess fails with an opaque HTTP 500 rather than a clear error. |
| `groups_json` | string \| null | 否 | JSON object mapping each group name to a list of accession names/ids. 默认 `None` |
| `metadata_tsv` | string \| null | 否 | Path to a metadata TSV (import_metadata format) used to define groups. 默认 `None` |
| `group_column` | string \| null | 否 | Column in the metadata TSV holding the group/population label. 默认 `None` |
| `id_column` | string | 否 | Column in the metadata TSV holding the individual/accession id (default 'individual'). 默认 `individual` |
| `max_markers` | integer \| null | 否 | Cap analysis to the first N markers in canonical Gigwa search order; omit to use all. 默认 `None` |
| `method` | string | 否 | Genotype source: 'vcf' (full export, cached) or 'allelematrix' (paged, server-side subset). 默认 `vcf` |
| `region` | string \| null | 否 | Restrict analysis to a genomic window: 'chrom' or 'chrom:start-end' (1-based). 默认 `None` |
| `output_dir` | string \| null | 否 | Directory for the output CSV(s) (default ./gigwa_results/&lt;module&gt;/). 默认 `None` |

??? note "上游说明（英文原文）"

    Per-population diversity: He, Ho, Fis, MAF, % polymorphic, allelic richness.

    Define groups the same way as ``diversity_fst`` — either ``groups_json``
    ``\{group: \[names\]\}`` or ``metadata_tsv`` + ``group_column``. For each group computes
    n, % polymorphic markers, mean MAF, Nei's He, observed Ho, Fis (1−Ho/He), mean
    observed allelic richness, and **rarefied** allelic richness (rarefied to the
    smallest group's gene-copy count so unequal group sizes are comparable). Writes
    ``diversity_by_group.csv``.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gigwa-mcp',
    'diversity_by_group',
    {'variant_set_db_id': '...'},
)
```

## `diversity_core_collection` {#diversity_core_collection}

Selects a core collection of accessions that maximises captured marker-allele diversity (Core Hunter style).

- **类别：**`breeding_population`（Breeding genotyping, markers, population genetics and strain resources）
- **必填参数：**`variant_set_db_id`
- **关键词：**core collection, germplasm, allelic diversity, genebank

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant_set_db_id` | string | 是 | BrAPI variantSetDbId identifying the run (MODULE§project§run) -- copy the exact string from list_variant_sets / list_content, never assemble one by hand: the middle segment is a numeric project index, not the project's name, and a wrong guess fails with an opaque HTTP 500 rather than a clear error. |
| `size` | integer \| null | 否 | Explicit core-collection size (number of accessions); overrides fraction. 默认 `None` |
| `fraction` | number | 否 | Core-collection size as a fraction of all accessions (default 0.1). 默认 `0.1` |
| `max_markers` | integer \| null | 否 | Cap analysis to the first N markers in canonical Gigwa search order; omit to use all. 默认 `None` |
| `method` | string | 否 | Genotype source: 'vcf' (full export, cached) or 'allelematrix' (paged, server-side subset). 默认 `vcf` |
| `region` | string \| null | 否 | Restrict analysis to a genomic window: 'chrom' or 'chrom:start-end' (1-based). 默认 `None` |
| `output_dir` | string \| null | 否 | Directory for the output CSV(s) (default ./gigwa_results/&lt;module&gt;/). 默认 `None` |

??? note "上游说明（英文原文）"

    Select a core collection that maximises captured allelic diversity.

    Greedy allele-coverage selection (Core-Hunter style): repeatedly add the accession
    that contributes the most not-yet-captured marker-alleles. Pick the core ``size``
    directly, or as ``fraction`` of all accessions (default 10%). Writes
    ``core_collection.csv`` (rank, accession, cumulative allele coverage) and reports the
    fraction of total allelic diversity the core captures.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gigwa-mcp',
    'diversity_core_collection',
    {'variant_set_db_id': '...'},
)
```

## `diversity_fst` {#diversity_fst}

Computes pairwise Weir and Cockerham Fst between user-defined groups of samples from genotypes.

- **类别：**`breeding_population`（Breeding genotyping, markers, population genetics and strain resources）
- **必填参数：**`variant_set_db_id`
- **关键词：**fst, weir cockerham, population differentiation, groups

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant_set_db_id` | string | 是 | BrAPI variantSetDbId identifying the run (MODULE§project§run) -- copy the exact string from list_variant_sets / list_content, never assemble one by hand: the middle segment is a numeric project index, not the project's name, and a wrong guess fails with an opaque HTTP 500 rather than a clear error. |
| `groups_json` | string \| null | 否 | JSON object mapping each group name to a list of accession names/ids. 默认 `None` |
| `metadata_tsv` | string \| null | 否 | Path to a metadata TSV (import_metadata format) used to define groups. 默认 `None` |
| `group_column` | string \| null | 否 | Column in the metadata TSV holding the group/population label. 默认 `None` |
| `id_column` | string | 否 | Column in the metadata TSV holding the individual/accession id (default 'individual'). 默认 `individual` |
| `max_markers` | integer \| null | 否 | Cap analysis to the first N markers in canonical Gigwa search order; omit to use all. 默认 `None` |
| `method` | string | 否 | Genotype source: 'vcf' (full export, cached) or 'allelematrix' (paged, server-side subset). 默认 `vcf` |
| `region` | string \| null | 否 | Restrict analysis to a genomic window: 'chrom' or 'chrom:start-end' (1-based). 默认 `None` |
| `output_dir` | string \| null | 否 | Directory for the output CSV(s) (default ./gigwa_results/&lt;module&gt;/). 默认 `None` |

??? note "上游说明（英文原文）"

    Pairwise Weir &amp; Cockerham Fst between groups of samples.

    Define the groups one of two ways:
    - ``groups_json`` — a JSON object mapping each group name to a list of accession
    names (or callset ids), e.g. ``\{"north": \["112","156"\], "south": \["11","42"\]\}``.
    - ``metadata_tsv`` + ``group_column`` — read groups from a metadata TSV (the same
    file format used by ``import_metadata``), keyed on ``id_column`` (default
    ``individual``) and grouped by ``group_column``.

    Writes ``fst_pairwise.csv`` with the Fst for every group pair. (Server-side BrAPI
    attributes are not used for grouping — that endpoint is unavailable on the target
    Gigwa 2.12 build.)

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gigwa-mcp',
    'diversity_fst',
    {'variant_set_db_id': '...'},
)
```

## `diversity_kinship` {#diversity_kinship}

Computes a VanRaden genomic relationship matrix from alt-allele dosage and reports most-related pairs.

- **类别：**`breeding_population`（Breeding genotyping, markers, population genetics and strain resources）
- **必填参数：**`variant_set_db_id`
- **关键词：**kinship, grm, vanraden, relatedness, inbreeding

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant_set_db_id` | string | 是 | BrAPI variantSetDbId identifying the run (MODULE§project§run) -- copy the exact string from list_variant_sets / list_content, never assemble one by hand: the middle segment is a numeric project index, not the project's name, and a wrong guess fails with an opaque HTTP 500 rather than a clear error. |
| `max_markers` | integer \| null | 否 | Cap analysis to the first N markers in canonical Gigwa search order; omit to use all. 默认 `None` |
| `top_pairs` | integer | 否 | How many most-related sample pairs to report. 默认 `15` |
| `method` | string | 否 | Genotype source: 'vcf' (full export, cached) or 'allelematrix' (paged, server-side subset). 默认 `vcf` |
| `region` | string \| null | 否 | Restrict analysis to a genomic window: 'chrom' or 'chrom:start-end' (1-based). 默认 `None` |
| `output_dir` | string \| null | 否 | Directory for the output CSV(s) (default ./gigwa_results/&lt;module&gt;/). 默认 `None` |

??? note "上游说明（英文原文）"

    VanRaden genomic relationship (kinship) matrix.

    Computes G = ZZ'/(2 Σp(1-p)) from alt dosage. Writes the full matrix as
    ``kinship_matrix.csv`` (samples × samples) and reports the most-related pairs
    and the diagonal (self-relationship / inbreeding) range. For large sets pass
    ``method="allelematrix"`` + ``max_markers`` to avoid a full VCF export.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gigwa-mcp',
    'diversity_kinship',
    {'variant_set_db_id': '...'},
)
```

## `diversity_pca` {#diversity_pca}

Runs PCA on the alt-allele dosage matrix to describe population structure and flag outlier samples.

- **类别：**`breeding_population`（Breeding genotyping, markers, population genetics and strain resources）
- **必填参数：**`variant_set_db_id`
- **关键词：**pca, population structure, dosage, outliers

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant_set_db_id` | string | 是 | BrAPI variantSetDbId identifying the run (MODULE§project§run) -- copy the exact string from list_variant_sets / list_content, never assemble one by hand: the middle segment is a numeric project index, not the project's name, and a wrong guess fails with an opaque HTTP 500 rather than a clear error. |
| `n_components` | integer | 否 | Number of principal components to compute. 默认 `10` |
| `max_markers` | integer \| null | 否 | Cap analysis to the first N markers in canonical Gigwa search order; omit to use all. 默认 `None` |
| `outlier_sd` | number | 否 | Flag points more than this many standard deviations from the mean. 默认 `6.0` |
| `metadata_tsv` | string \| null | 否 | Path to a metadata TSV (import_metadata format) used to define groups. 默认 `None` |
| `group_column` | string \| null | 否 | Column in the metadata TSV holding the group/population label. 默认 `None` |
| `id_column` | string | 否 | Column in the metadata TSV holding the individual/accession id (default 'individual'). 默认 `individual` |
| `method` | string | 否 | Genotype source: 'vcf' (full export, cached) or 'allelematrix' (paged, server-side subset). 默认 `vcf` |
| `region` | string \| null | 否 | Restrict analysis to a genomic window: 'chrom' or 'chrom:start-end' (1-based). 默认 `None` |
| `output_dir` | string \| null | 否 | Directory for the output CSV(s) (default ./gigwa_results/&lt;module&gt;/). 默认 `None` |

??? note "上游说明（英文原文）"

    Principal component analysis of population structure.

    Runs PCA on the alt-allele dosage matrix (monomorphic markers dropped, missing
    mean-imputed, Patterson scaling). Writes ``pca_coords.csv`` (per-sample PC
    coordinates) and reports variance explained plus any PC1/PC2 outlier samples
    (beyond ``outlier_sd`` SD). Pass ``metadata_tsv`` + ``group_column`` to add a
    ``group`` column (population label per sample) for colouring the PC plot. For large
    sets pass ``method="allelematrix"`` + ``max_markers`` to avoid a full VCF export.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gigwa-mcp',
    'diversity_pca',
    {'variant_set_db_id': '...'},
)
```

## `diversity_structure` {#diversity_structure}

Clusters samples into genetic groups via PCA and K-means, choosing K by pseudo-F.

- **类别：**`breeding_population`（Breeding genotyping, markers, population genetics and strain resources）
- **必填参数：**`variant_set_db_id`
- **关键词：**population structure, k-means, admixture, clusters

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant_set_db_id` | string | 是 | BrAPI variantSetDbId identifying the run (MODULE§project§run) -- copy the exact string from list_variant_sets / list_content, never assemble one by hand: the middle segment is a numeric project index, not the project's name, and a wrong guess fails with an opaque HTTP 500 rather than a clear error. |
| `k_min` | integer | 否 | Smallest number of clusters (K) to evaluate. 默认 `2` |
| `k_max` | integer | 否 | Largest number of clusters (K) to evaluate. 默认 `10` |
| `max_markers` | integer \| null | 否 | Cap analysis to the first N markers in canonical Gigwa search order; omit to use all. 默认 `None` |
| `method` | string | 否 | Genotype source: 'vcf' (full export, cached) or 'allelematrix' (paged, server-side subset). 默认 `vcf` |
| `region` | string \| null | 否 | Restrict analysis to a genomic window: 'chrom' or 'chrom:start-end' (1-based). 默认 `None` |
| `output_dir` | string \| null | 否 | Directory for the output CSV(s) (default ./gigwa_results/&lt;module&gt;/). 默认 `None` |

??? note "上游说明（英文原文）"

    Lightweight population-structure clustering (PCA + K-means, in-Python).

    Reduces the alt-dosage matrix with PCA (Patterson scaling), then runs K-means for
    K in ``k_min..k_max`` and picks the K with the highest pseudo-F (Calinski-Harabasz)
    between/within variance ratio — a clear maximum when groups are well separated.
    Writes ``structure_clusters.csv`` (sample, assigned cluster at the best K, PC coords)
    and reports the chosen K with cluster sizes. (No external ADMIXTURE binary — computed
    entirely in Python, consistent with the rest of the analysis layer.)

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gigwa-mcp',
    'diversity_structure',
    {'variant_set_db_id': '...'},
)
```

## `diversity_summary` {#diversity_summary}

Computes per-marker MAF, expected/observed heterozygosity and PIC plus dataset means.

- **类别：**`breeding_population`（Breeding genotyping, markers, population genetics and strain resources）
- **必填参数：**`variant_set_db_id`
- **关键词：**maf, pic, heterozygosity, marker diversity

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant_set_db_id` | string | 是 | BrAPI variantSetDbId identifying the run (MODULE§project§run) -- copy the exact string from list_variant_sets / list_content, never assemble one by hand: the middle segment is a numeric project index, not the project's name, and a wrong guess fails with an opaque HTTP 500 rather than a clear error. |
| `max_markers` | integer \| null | 否 | Cap analysis to the first N markers in canonical Gigwa search order; omit to use all. 默认 `None` |
| `method` | string | 否 | Genotype source: 'vcf' (full export, cached) or 'allelematrix' (paged, server-side subset). 默认 `vcf` |
| `region` | string \| null | 否 | Restrict analysis to a genomic window: 'chrom' or 'chrom:start-end' (1-based). 默认 `None` |
| `output_dir` | string \| null | 否 | Directory for the output CSV(s) (default ./gigwa_results/&lt;module&gt;/). 默认 `None` |

??? note "上游说明（英文原文）"

    Per-marker diversity statistics (MAF, He, Ho, PIC) and dataset means.

    He is Nei's gene diversity (1 - Σpᵢ²), Ho is observed heterozygosity, PIC is
    polymorphism information content. Writes ``diversity_markers.csv``. For large sets
    pass ``method="allelematrix"`` + ``max_markers`` to sample server-side.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gigwa-mcp',
    'diversity_summary',
    {'variant_set_db_id': '...'},
)
```

## `diversity_tree` {#diversity_tree}

Builds a UPGMA Newick tree of accessions from pairwise IBS allele-sharing distance.

- **类别：**`breeding_population`（Breeding genotyping, markers, population genetics and strain resources）
- **必填参数：**`variant_set_db_id`
- **关键词：**upgma, ibs, newick, phylogeny, dendrogram

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant_set_db_id` | string | 是 | BrAPI variantSetDbId identifying the run (MODULE§project§run) -- copy the exact string from list_variant_sets / list_content, never assemble one by hand: the middle segment is a numeric project index, not the project's name, and a wrong guess fails with an opaque HTTP 500 rather than a clear error. |
| `max_markers` | integer \| null | 否 | Cap analysis to the first N markers in canonical Gigwa search order; omit to use all. 默认 `5000` |
| `method` | string | 否 | Genotype source: 'vcf' (full export, cached) or 'allelematrix' (paged, server-side subset). 默认 `vcf` |
| `region` | string \| null | 否 | Restrict analysis to a genomic window: 'chrom' or 'chrom:start-end' (1-based). 默认 `None` |
| `output_dir` | string \| null | 否 | Directory for the output CSV(s) (default ./gigwa_results/&lt;module&gt;/). 默认 `None` |

??? note "上游说明（英文原文）"

    UPGMA dendrogram of accessions from IBS allele-sharing distance (Newick).

    Builds a pairwise IBS similarity matrix, converts to distance (1 − IBS), and writes a
    UPGMA tree as ``tree.nwk`` (standard Newick, loadable in FigTree / iTOL / ape). Marker
    subsampling (``max_markers``) keeps it tractable on large sets.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gigwa-mcp',
    'diversity_tree',
    {'variant_set_db_id': '...'},
)
```

## `export_genotypes` {#export_genotypes}

Exports a Gigwa variant set or filtered subset to VCF, PLINK, Flapjack or other genotype formats.

- **类别：**`breeding_population`（Breeding genotyping, markers, population genetics and strain resources）
- **必填参数：**`variant_set_db_id`, `output_path`
- **关键词：**vcf, plink, flapjack, genotype export, maf filter

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant_set_db_id` | string | 是 | BrAPI variantSetDbId identifying the run (MODULE§project§run) -- copy the exact string from list_variant_sets / list_content, never assemble one by hand: the middle segment is a numeric project index, not the project's name, and a wrong guess fails with an opaque HTTP 500 rather than a clear error. |
| `output_path` | string | 是 | Destination file path for the export. |
| `format` | string | 否 | Export format name, e.g. VCF (default), PLINK, FLAPJACK or VCF.gz; see list_export_formats for what this instance offers (and each format's type/ploidy restrictions). 默认 `VCF` |
| `region` | string \| null | 否 | Restrict analysis to a genomic window: 'chrom' or 'chrom:start-end' (1-based). 默认 `None` |
| `selected_variant_types` | string \| null | 否 | Restrict the export to these variant types, ';'-joined (e.g. 'SNP' or 'SNP;INDEL'); omit for all types. 默认 `None` |
| `min_maf` | number \| null | 否 | Minimum minor-allele frequency (0-1). 默认 `None` |
| `max_maf` | number \| null | 否 | Maximum minor-allele frequency (0-1). 默认 `None` |
| `max_missing_data` | number \| null | 否 | Maximum per-variant missing-data fraction (0-1). 默认 `None` |
| `individuals` | array \| null | 否 | Individual-level identifiers to include in the export (e.g. from filter_individuals_by_metadata), not sample/callset ids; omit for all. Gigwa resolves each individual to all of its samples/callsets across runs server-side, including when an individual has more than one, so no manual sample mapping or dedup is needed. 默认 `None` |
| `metadata_fields` | array \| null | 否 | Individual metadata columns to embed in the export (from get_germplasm_metadata); omit for none. 默认 `None` |
| `keep_on_server` | boolean | 否 | Also leave a copy of the export in the user's Gigwa temp-output area after downloading it here. 默认 `False` |
| `wait` | boolean | 否 | Block until the job finishes (True, default) or return immediately once it's kicked off (False) -- an import returns a progress token to poll with get_import_progress, an export returns a download URL to poll with get_export_progress and retrieve with fetch_export_file. 默认 `True` |
| `timeout` | number | 否 | Maximum seconds to wait for the export to complete. 默认 `1800.0` |

??? note "上游说明（英文原文）"

    Export a variant set — or a filtered/selected subset of it — to a file.

    With none of the filter/selection parameters set, exports the whole set via Gigwa's
    plain per-format export (``format`` one of ``VCF`` (default), ``PLINK`` or
    ``FLAPJACK``; availability varies by build — check ``list_export_formats`` for what
    this instance actually offers, and each format's variant-type/ploidy restrictions).

    Passing any of ``region``, ``selected_variant_types``, ``min_maf``/``max_maf``,
    ``max_missing_data``, ``individuals`` or ``metadata_fields`` instead drives Gigwa's
    selection-aware export — the same endpoint the Gigwa web UI uses for a filtered
    download — which additionally accepts any format the server advertises, including the
    bgzipped ``"VCF.gz"``. ``keep_on_server`` leaves a copy in the user's Gigwa temp-output
    area after this downloads it (default False).

    For large sets this can take a while; raise ``timeout`` (seconds), or set
    ``wait=False`` to return immediately once the export is kicked off instead of blocking
    for the whole thing. A server that queues the work returns a download URL; check progress
    with ``get_export_progress`` and, once it reports complete, retrieve the file with
    ``fetch_export_file`` (that URL, plus ``output_path``). ``wait=False`` always goes
    through the selection-aware endpoint, even with no filters set. Some Gigwa builds
    return the completed bytes immediately; in that case they are written to
    ``output_path`` and no progress/fetch step is needed.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gigwa-mcp',
    'export_genotypes',
    {'variant_set_db_id': '...', 'output_path': '...'},
)
```

## `filter_individuals_by_metadata` {#filter_individuals_by_metadata}

Selects individuals whose stored germplasm metadata match field/value filters, returning IDs for genotype export.

- **类别：**`breeding_population`（Breeding genotyping, markers, population genetics and strain resources）
- **必填参数：**`variant_set_db_id`, `filters_json`
- **关键词：**germplasm, individuals, metadata filter, selection

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant_set_db_id` | string | 是 | BrAPI variantSetDbId identifying the run (MODULE§project§run) -- copy the exact string from list_variant_sets / list_content, never assemble one by hand: the middle segment is a numeric project index, not the project's name, and a wrong guess fails with an opaque HTTP 500 rather than a clear error. |
| `filters_json` | string | 是 | JSON object mapping each metadata field name to a list of acceptable values, e.g. \{"GroupK4": \["cA"\]\} (see list_metadata_values for field/value names). Multiple fields AND together; multiple values for one field OR together. |

??? note "上游说明（英文原文）"

    Select individuals whose stored metadata matches the given field/value filters.

    ``filters_json`` is a JSON object mapping each metadata field name (see
    ``list_metadata_values``) to a list of acceptable values, e.g. ``\{"GroupK4": \["cA"\]\}``.
    Multiple fields combine with AND; multiple values for one field combine with OR.
    Returns the matching **individual**-level identifiers, ready to pass directly to
    ``export_genotypes``'s ``individuals`` parameter. Diversity ``groups_json`` instead
    matches analysis sample names; obtain those from ``get_germplasm_metadata``'s
    ``sample_name`` column. Gigwa's export resolves each individual to all of its samples
    across runs server-side (verified against the Gigwa server source), including when an
    individual has more than one sample.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gigwa-mcp',
    'filter_individuals_by_metadata',
    {'variant_set_db_id': '...', 'filters_json': '...'},
)
```

## `get_germplasm_metadata` {#get_germplasm_metadata}

Fetches server-stored per-individual germplasm attributes for a Gigwa database and writes them to CSV.

- **类别：**`breeding_population`（Breeding genotyping, markers, population genetics and strain resources）
- **必填参数：**`variant_set_db_id`
- **关键词：**germplasm, accession attributes, passport data, metadata

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant_set_db_id` | string | 是 | BrAPI variantSetDbId identifying the run (MODULE§project§run) -- copy the exact string from list_variant_sets / list_content, never assemble one by hand: the middle segment is a numeric project index, not the project's name, and a wrong guess fails with an opaque HTTP 500 rather than a clear error. |
| `output_dir` | string \| null | 否 | Directory for the output CSV(s) (default ./gigwa_results/&lt;module&gt;/). 默认 `None` |

??? note "上游说明（英文原文）"

    Fetch server-stored per-individual metadata (germplasm attributes) for a database.

    Reads the attributes already stored in Gigwa (imported earlier via ``import_metadata``
    or a BrAPI source) for the module of ``variant_set_db_id`` and writes
    ``germplasm_metadata.csv`` — one row per accession, with attribute columns.

    **Which column to join on.** ``germplasm_name`` is Gigwa's *individual* id, which on
    real data does **not** match the sample names the analysis tools use. Join on the
    ``sample_name`` column instead: it is the name ``diversity_fst`` /
    ``diversity_by_group`` match against (via ``metadata_tsv``, which expects **tab**-
    separated input — convert the CSV first). Where one individual carries several samples
    the column holds them ``;``-joined and needs splitting first. It is blank when the
    callsets could not be read.

    For the sample/callset level itself — raw ``callSetName``, ``sampleDbId`` and per-sample
    ``additionalInfo`` — use ``search_callsets`` instead.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gigwa-mcp',
    'get_germplasm_metadata',
    {'variant_set_db_id': '...'},
)
```

## `import_dartseq` {#import_dartseq}

Converts DArTseq SNP/Silico-DArT xlsx reports to VCF with correct genotype calling and imports them into a Gigwa database.

- **类别：**`breeding_population`（Breeding genotyping, markers, population genetics and strain resources）
- **必填参数：**`module`, `project`, `run`
- **关键词：**dartseq, silicodart, genotyping, vcf, markers

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `module` | string | 是 | Target Gigwa database (module) name. |
| `project` | string | 是 | Target project name within the database. |
| `run` | string | 是 | Target run name within the project. |
| `snp_xlsx` | string \| null | 否 | Path to a DArTseq SNP xlsx report. 默认 `None` |
| `silico_xlsx` | string \| null | 否 | Path to a Silico-DArT xlsx report. 默认 `None` |
| `technology` | string | 否 | Free-text genotyping technology label (e.g. 'DArTseq', 'WGS', 'GBS'). 默认 `DArTseq` |
| `ploidy` | integer | 否 | Sample ploidy (default 2). 默认 `2` |
| `skip_monomorphic` | boolean | 否 | Drop non-variant (monomorphic) markers during import. 默认 `False` |
| `clear_project_data` | boolean | 否 | Replace any existing data in the project before importing. 默认 `False` |
| `reference_fasta` | string \| null | 否 | Path to a reference genome FASTA or a prebuilt minimap2 .mmi index, for genome-anchoring. 默认 `None` |
| `positions_csv` | string \| null | 否 | Path to a dartseq_positions.csv (from map_dartseq_to_reference) to reuse instead of re-aligning. 默认 `None` |
| `min_mapq` | integer | 否 | Minimum mapping quality for a tag to count as uniquely mapped. 默认 `20` |
| `wait` | boolean | 否 | Block until the job finishes (True, default) or return immediately once it's kicked off (False) -- an import returns a progress token to poll with get_import_progress, an export returns a download URL to poll with get_export_progress and retrieve with fetch_export_file. 默认 `True` |

??? note "上游说明（英文原文）"

    Import DArTseq data from xlsx report(s) into Gigwa.

    Converts the DArTseq SNP and/or Silico-DArT xlsx report(s) to a standard VCF —
    doing the 2-row genotype calling in Python (so reference homozygotes are not
    mis-imported as heterozygous, as Gigwa's built-in DArT parser does) — and
    uploads it to create/append a database (``module``), ``project`` and ``run``.

    Provide at least one of ``snp_xlsx`` / ``silico_xlsx`` (absolute paths). SNP
    and Silico use different allele models; importing both into the *same* run is
    unusual — prefer separate runs unless you specifically intend to combine them.

    If ``reference_fasta`` is given (a reference genome FASTA *or* a prebuilt
    minimap2 ``.mmi`` index — an ``.mmi`` is loaded directly with no re-indexing,
    preferred for large genomes), the SNP markers' tag sequences are aligned to it
    and uniquely-mapped markers (mapq ≥ ``min_mapq``) are imported genome-anchored
    (real chromosome/position); the rest stay on an ``Unmapped`` contig. Without it,
    all markers go on ``Unmapped``.

    ``positions_csv`` reuses a mapping already produced by
    ``map_dartseq_to_reference`` (its ``dartseq_positions.csv``) instead of
    re-aligning — much faster when you've already inspected the mapping. Provide
    either ``reference_fasta`` or ``positions_csv``, not both.

    Set ``clear_project_data=True`` to replace any existing data in the project,
    ``skip_monomorphic=True`` to drop non-variant markers, and ``wait=False`` to
    return immediately with a progress token instead of blocking until done.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gigwa-mcp',
    'import_dartseq',
    {'module': '...', 'project': '...', 'run': '...'},
)
```

## `import_metadata` {#import_metadata}

Imports a per-individual metadata TSV into an existing Gigwa database as searchable attributes.

- **类别：**`breeding_population`（Breeding genotyping, markers, population genetics and strain resources）
- **必填参数：**`tsv_path`, `module`
- **关键词：**germplasm, metadata import, tsv, individuals

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `tsv_path` | string | 是 | Path to the metadata TSV file. |
| `module` | string | 是 | Target Gigwa database (module) name. |
| `metadata_type` | string | 否 | Metadata entity type / id-column name (default 'individual'). 默认 `individual` |
| `validate_first` | boolean | 否 | Validate the metadata file before importing. 默认 `True` |

??? note "上游说明（英文原文）"

    Import individual metadata (per-individual attributes) into an existing Gigwa database.

    The file is a TSV whose first column header equals ``metadata_type``
    (``individual`` for individual metadata) and whose values match the
    individual/sample names already present in the database. Remaining columns
    become searchable attributes. By default the file is validated first; set
    ``validate_first=False`` to skip that check.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gigwa-mcp',
    'import_metadata',
    {'tsv_path': '...', 'module': '...'},
)
```

## `import_vcf` {#import_vcf}

Imports a VCF file into a Gigwa database/project/run, optionally replacing existing project data.

- **类别：**`breeding_population`（Breeding genotyping, markers, population genetics and strain resources）
- **必填参数：**`vcf_path`, `module`, `project`, `run`
- **关键词：**vcf, genotype import, gbs, wgs, ploidy

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `vcf_path` | string | 是 | Path to the VCF file (.vcf or .vcf.gz) to import. |
| `module` | string | 是 | Target Gigwa database (module) name. |
| `project` | string | 是 | Target project name within the database. |
| `run` | string | 是 | Target run name within the project. |
| `technology` | string \| null | 否 | Free-text genotyping technology label (e.g. 'DArTseq', 'WGS', 'GBS'). 默认 `None` |
| `ploidy` | integer | 否 | Sample ploidy (default 2). 默认 `2` |
| `skip_monomorphic` | boolean | 否 | Drop non-variant (monomorphic) markers during import. 默认 `False` |
| `clear_project_data` | boolean | 否 | Replace any existing data in the project before importing. 默认 `False` |
| `wait` | boolean | 否 | Block until the job finishes (True, default) or return immediately once it's kicked off (False) -- an import returns a progress token to poll with get_import_progress, an export returns a download URL to poll with get_export_progress and retrieve with fetch_export_file. 默认 `True` |

??? note "上游说明（英文原文）"

    Import a VCF file (``.vcf`` or ``.vcf.gz``) into Gigwa.

    Uploads the VCF to create/append a database (``module``), ``project`` and
    ``run``. ``technology`` is optional free-text (e.g. 'WGS', 'GBS'). Use
    ``clear_project_data=True`` to replace existing project data and ``wait=False``
    to return a progress token instead of blocking.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gigwa-mcp',
    'import_vcf',
    {'vcf_path': '...', 'module': '...', 'project': '...', 'run': '...'},
)
```

## `list_content` {#list_content}

Lists Gigwa databases, projects and runs with their exact BrAPI variantSetDbIds.

- **类别：**`breeding_population`（Breeding genotyping, markers, population genetics and strain resources）
- **必填参数：**无
- **关键词：**variant sets, databases, projects, runs, brapi

??? note "上游说明（英文原文）"

    List the databases, projects and runs currently hosted on the Gigwa server.

    Each run's exact ``variant_set_db_id`` is printed alongside it — use that string
    verbatim in other tools, do not assemble one by hand. In particular, the middle
    segment is a numeric project index (``1``, ``2``, ...), not the project's name (e.g.
    project 'refNB' still has id segment '1') — using the name there is a common mistake
    and fails every downstream call with an opaque HTTP 500 rather than a clear error.
    ``list_variant_sets`` gives the same ids in a flatter list, if preferred.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gigwa-mcp',
    'list_content',
    {},
)
```

## `list_metadata_values` {#list_metadata_values}

Lists individual-metadata fields available for a Gigwa database and their distinct values.

- **类别：**`breeding_population`（Breeding genotyping, markers, population genetics and strain resources）
- **必填参数：**`variant_set_db_id`
- **关键词：**germplasm, metadata fields, filter values

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant_set_db_id` | string | 是 | BrAPI variantSetDbId identifying the run (MODULE§project§run) -- copy the exact string from list_variant_sets / list_content, never assemble one by hand: the middle segment is a numeric project index, not the project's name, and a wrong guess fails with an opaque HTTP 500 rather than a clear error. |

??? note "上游说明（英文原文）"

    List the individual-metadata fields available for a database and their distinct values.

    Backed by the endpoint the Gigwa web UI itself uses to populate its metadata-based
    selection filters, so it works even on builds where ``get_germplasm_metadata``'s BrAPI
    fallback returns no attributes. Use this to discover field names/values before calling
    ``filter_individuals_by_metadata``.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gigwa-mcp',
    'list_metadata_values',
    {'variant_set_db_id': '...'},
)
```

## `list_sequences` {#list_sequences}

Lists reference sequences (chromosomes/contigs) present in a Gigwa variant set.

- **类别：**`breeding_population`（Breeding genotyping, markers, population genetics and strain resources）
- **必填参数：**`variant_set_db_id`
- **关键词：**chromosomes, contigs, reference sequences, region filter

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant_set_db_id` | string | 是 | BrAPI variantSetDbId identifying the run (MODULE§project§run) -- copy the exact string from list_variant_sets / list_content, never assemble one by hand: the middle segment is a numeric project index, not the project's name, and a wrong guess fails with an opaque HTTP 500 rather than a clear error. |

??? note "上游说明（英文原文）"

    List the reference sequences (chromosomes/contigs) available in a variant set.

    Use this to discover valid ``reference_name`` values for the region filters on
    ``count_variants`` / ``search_variants`` / the QC &amp; diversity tools.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gigwa-mcp',
    'list_sequences',
    {'variant_set_db_id': '...'},
)
```

## `list_variant_sets` {#list_variant_sets}

Lists every Gigwa variant set (run) with its BrAPI variantSetDbId and variant/callset counts.

- **类别：**`breeding_population`（Breeding genotyping, markers, population genetics and strain resources）
- **必填参数：**无
- **关键词：**variant sets, brapi, callsets, runs

??? note "上游说明（英文原文）"

    List every variant set (run) with its exact BrAPI variantSetDbId.

    The other tools take a ``variant_set_db_id``; this returns those ids directly (plus
    name and variant/callset counts when the server provides them), complementing the
    human-readable database/project/run view from ``list_content``.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gigwa-mcp',
    'list_variant_sets',
    {},
)
```

## `map_dartseq_to_reference` {#map_dartseq_to_reference}

Aligns DArTseq SNP marker tag sequences to a reference genome with minimap2 to infer chromosome, position and strand.

- **类别：**`breeding_population`（Breeding genotyping, markers, population genetics and strain resources）
- **必填参数：**`snp_xlsx`, `reference_fasta`
- **关键词：**dartseq, minimap2, marker positions, reference genome, alignment

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `snp_xlsx` | string | 是 | Path to a DArTseq SNP xlsx report. |
| `reference_fasta` | string | 是 | Path to a reference genome FASTA or a prebuilt minimap2 .mmi index, for genome-anchoring. |
| `output_dir` | string \| null | 否 | Directory for the output CSV(s) (default ./gigwa_results/&lt;module&gt;/). 默认 `None` |
| `min_mapq` | integer | 否 | Minimum mapping quality for a tag to count as uniquely mapped. 默认 `20` |
| `preset` | string | 否 | minimap2 preset (default 'sr' for short reads). 默认 `sr` |
| `backend` | string | 否 | Aligner backend: 'auto' (minimap2 CLI if available, else mappy), 'cli', or 'mappy'. 默认 `auto` |

??? note "上游说明（英文原文）"

    Guess genomic positions for DArTseq SNP markers by aligning their tag sequences.

    Aligns each marker's ~69 bp ``AlleleSequence`` tag to ``reference_fasta`` (a
    reference genome FASTA, or a prebuilt minimap2 ``.mmi`` index) and reports the
    inferred chromosome, position and strand of each SNP. Writes
    ``dartseq_positions.csv`` (allele_id, chrom, pos, strand, mapq, ref, alt, status).
    The result can be passed to ``import_dartseq`` (``reference_fasta=``) to import
    the data genome-anchored instead of on an ``Unmapped`` contig.

    ``backend``: ``"auto"`` uses the minimap2 CLI when available (streams over
    multi-part indexes → bounded RAM, best for large multi-gigabase genomes),
    falling back to the in-process ``mappy`` binding. Markers are classified
    ``unique`` (mapq ≥ ``min_mapq``), ``multi`` (ambiguous), or ``unmapped``.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gigwa-mcp',
    'map_dartseq_to_reference',
    {'snp_xlsx': '...', 'reference_fasta': '...'},
)
```

## `qc_call_rate` {#qc_call_rate}

Computes per-sample and per-marker call rate (missingness) for a variant set and flags low-quality ones.

- **类别：**`breeding_population`（Breeding genotyping, markers, population genetics and strain resources）
- **必填参数：**`variant_set_db_id`
- **关键词：**call rate, missingness, genotype qc, markers

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant_set_db_id` | string | 是 | BrAPI variantSetDbId identifying the run (MODULE§project§run) -- copy the exact string from list_variant_sets / list_content, never assemble one by hand: the middle segment is a numeric project index, not the project's name, and a wrong guess fails with an opaque HTTP 500 rather than a clear error. |
| `min_sample_call_rate` | number | 否 | Flag samples with call rate below this (0-1). 默认 `0.5` |
| `min_marker_call_rate` | number | 否 | Flag markers with call rate below this (0-1). 默认 `0.5` |
| `max_markers` | integer \| null | 否 | Cap analysis to the first N markers in canonical Gigwa search order; omit to use all. 默认 `None` |
| `method` | string | 否 | Genotype source: 'vcf' (full export, cached) or 'allelematrix' (paged, server-side subset). 默认 `vcf` |
| `region` | string \| null | 否 | Restrict analysis to a genomic window: 'chrom' or 'chrom:start-end' (1-based). 默认 `None` |
| `output_dir` | string \| null | 否 | Directory for the output CSV(s) (default ./gigwa_results/&lt;module&gt;/). 默认 `None` |

??? note "上游说明（英文原文）"

    Per-sample and per-marker call rate (missingness) QC for a variant set.

    Flags samples/markers below the given thresholds. Writes
    ``call_rate_samples.csv`` and ``call_rate_markers.csv`` and returns a summary
    with the overall call rate and the worst offenders. ``variant_set_db_id`` is a
    BrAPI variantSetDbId (from ``list_content`` / BrAPI variantsets). For large
    production sets pass ``method="allelematrix"`` with ``max_markers`` (e.g. 20000)
    to estimate from a server-side marker subset instead of a full VCF export.
    ``region`` (``"chrom"`` or ``"chrom:start-end"``, 1-based; from ``list_sequences``)
    restricts the analysis to one genomic window — available on every QC/diversity tool.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gigwa-mcp',
    'qc_call_rate',
    {'variant_set_db_id': '...'},
)
```

## `qc_duplicate_accessions` {#qc_duplicate_accessions}

Detects duplicate or clonal accessions by pairwise IBS similarity and groups them into duplicate sets.

- **类别：**`breeding_population`（Breeding genotyping, markers, population genetics and strain resources）
- **必填参数：**`variant_set_db_id`
- **关键词：**ibs, duplicates, genebank, clones, accessions

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant_set_db_id` | string | 是 | BrAPI variantSetDbId identifying the run (MODULE§project§run) -- copy the exact string from list_variant_sets / list_content, never assemble one by hand: the middle segment is a numeric project index, not the project's name, and a wrong guess fails with an opaque HTTP 500 rather than a clear error. |
| `similarity_threshold` | number | 否 | IBS similarity (0-1) at/above which accessions are grouped as duplicates. 默认 `0.95` |
| `max_markers` | integer \| null | 否 | Cap analysis to the first N markers in canonical Gigwa search order; omit to use all. 默认 `5000` |
| `method` | string | 否 | Genotype source: 'vcf' (full export, cached) or 'allelematrix' (paged, server-side subset). 默认 `vcf` |
| `region` | string \| null | 否 | Restrict analysis to a genomic window: 'chrom' or 'chrom:start-end' (1-based). 默认 `None` |
| `output_dir` | string \| null | 否 | Directory for the output CSV(s) (default ./gigwa_results/&lt;module&gt;/). 默认 `None` |

??? note "上游说明（英文原文）"

    Detect duplicate / clonal accessions via pairwise identity-by-state (IBS).

    Computes IBS allele-sharing similarity between every pair of samples and groups
    pairs at or above ``similarity_threshold`` into duplicate sets — the core
    genebank "cleaning" check for mislabelled duplicates and clones. By default
    subsamples to the first ``max_markers`` in canonical Gigwa search order for speed (set to 0/None to
    use all). Writes ``duplicate_pairs.csv`` and ``duplicate_groups.csv``. For large
    sets pass ``method="allelematrix"`` to fetch the marker subset without a full export.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gigwa-mcp',
    'qc_duplicate_accessions',
    {'variant_set_db_id': '...'},
)
```

## `qc_heterozygosity` {#qc_heterozygosity}

Computes per-sample observed heterozygosity and flags outliers suggesting contamination or inbreeding.

- **类别：**`breeding_population`（Breeding genotyping, markers, population genetics and strain resources）
- **必填参数：**`variant_set_db_id`
- **关键词：**heterozygosity, genotype qc, contamination, outliers

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant_set_db_id` | string | 是 | BrAPI variantSetDbId identifying the run (MODULE§project§run) -- copy the exact string from list_variant_sets / list_content, never assemble one by hand: the middle segment is a numeric project index, not the project's name, and a wrong guess fails with an opaque HTTP 500 rather than a clear error. |
| `outlier_sd` | number | 否 | Flag points more than this many standard deviations from the mean. 默认 `3.0` |
| `max_markers` | integer \| null | 否 | Cap analysis to the first N markers in canonical Gigwa search order; omit to use all. 默认 `None` |
| `method` | string | 否 | Genotype source: 'vcf' (full export, cached) or 'allelematrix' (paged, server-side subset). 默认 `vcf` |
| `region` | string \| null | 否 | Restrict analysis to a genomic window: 'chrom' or 'chrom:start-end' (1-based). 默认 `None` |
| `output_dir` | string \| null | 否 | Directory for the output CSV(s) (default ./gigwa_results/&lt;module&gt;/). 默认 `None` |

??? note "上游说明（英文原文）"

    Per-sample observed heterozygosity QC, flagging outliers.

    High Ho relative to the cohort suggests contamination or off-types; very low
    Ho suggests selfed/inbred or duplicated material. Flags samples more than
    ``outlier_sd`` standard deviations from the mean. Writes
    ``heterozygosity_samples.csv``. For large sets pass ``method="allelematrix"`` +
    ``max_markers`` to avoid a full VCF export.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gigwa-mcp',
    'qc_heterozygosity',
    {'variant_set_db_id': '...'},
)
```

## `qc_maf_filter` {#qc_maf_filter}

Reports per-marker MAF and missingness and counts markers that would be filtered, without modifying data.

- **类别：**`breeding_population`（Breeding genotyping, markers, population genetics and strain resources）
- **必填参数：**`variant_set_db_id`
- **关键词：**maf, missingness, marker filtering, monomorphic

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant_set_db_id` | string | 是 | BrAPI variantSetDbId identifying the run (MODULE§project§run) -- copy the exact string from list_variant_sets / list_content, never assemble one by hand: the middle segment is a numeric project index, not the project's name, and a wrong guess fails with an opaque HTTP 500 rather than a clear error. |
| `maf_threshold` | number | 否 | Minor-allele-frequency threshold below which markers are flagged. 默认 `0.05` |
| `max_missing` | number | 否 | Maximum per-marker missing-data fraction (0-1) before a marker is flagged. 默认 `0.5` |
| `max_markers` | integer \| null | 否 | Cap analysis to the first N markers in canonical Gigwa search order; omit to use all. 默认 `None` |
| `method` | string | 否 | Genotype source: 'vcf' (full export, cached) or 'allelematrix' (paged, server-side subset). 默认 `vcf` |
| `region` | string \| null | 否 | Restrict analysis to a genomic window: 'chrom' or 'chrom:start-end' (1-based). 默认 `None` |
| `output_dir` | string \| null | 否 | Directory for the output CSV(s) (default ./gigwa_results/&lt;module&gt;/). 默认 `None` |

??? note "上游说明（英文原文）"

    Report markers that would be filtered by MAF / missingness (no changes applied).

    Computes per-marker minor-allele frequency and missing rate, and counts how many
    markers are monomorphic, below ``maf_threshold``, or above ``max_missing`` missing.
    Writes ``marker_filter_stats.csv``. This is a report only — it does not modify Gigwa.
    For large sets pass ``method="allelematrix"`` + ``max_markers`` to sample server-side.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gigwa-mcp',
    'qc_maf_filter',
    {'variant_set_db_id': '...'},
)
```

## `search_callsets` {#search_callsets}

Retrieves BrAPI callsets for a run and writes per-sample metadata (names, IDs, attributes) to CSV.

- **类别：**`breeding_population`（Breeding genotyping, markers, population genetics and strain resources）
- **必填参数：**`variant_set_db_id`
- **关键词：**callsets, samples, brapi, sample metadata

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant_set_db_id` | string | 是 | BrAPI variantSetDbId identifying the run (MODULE§project§run) -- copy the exact string from list_variant_sets / list_content, never assemble one by hand: the middle segment is a numeric project index, not the project's name, and a wrong guess fails with an opaque HTTP 500 rather than a clear error. |
| `output_dir` | string \| null | 否 | Directory for the output CSV(s) (default ./gigwa_results/&lt;module&gt;/). 默认 `None` |

??? note "上游说明（英文原文）"

    Dump per-sample (callset) metadata for a run: names + ``additionalInfo`` attributes.

    Retrieves the run's callsets via BrAPI ``search/callsets`` and writes
    ``sample_metadata.csv`` (one row per sample) with ``sample_name``, ``callSetName``,
    ``sampleDbId``, ``callSetDbId`` and every attribute present in the callsets'
    ``additionalInfo`` (e.g. ICARDA_IG, SeedID, Country, Latitude, Longitude, SiteCode,
    PopulationType).

    ``sample_name`` is the same name the analysis tools use for the sample, so the file
    joins to their outputs; ``callSetName`` keeps the server's raw label, which on some
    builds is the more human-readable of the two.

    This is the sample/callset-level counterpart to ``get_germplasm_metadata``: use it
    when the germplasm (accession) level exposes no attributes but the samples do.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gigwa-mcp',
    'search_callsets',
    {'variant_set_db_id': '...'},
)
```

## `search_variants` {#search_variants}

Searches variants matching region, MAF and missing-data filters and writes their id, chrom, pos, ref, alt to CSV.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`variant_set_db_id`
- **关键词：**variant search, maf, region filter, brapi, snp list

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant_set_db_id` | string | 是 | BrAPI variantSetDbId identifying the run (MODULE§project§run) -- copy the exact string from list_variant_sets / list_content, never assemble one by hand: the middle segment is a numeric project index, not the project's name, and a wrong guess fails with an opaque HTTP 500 rather than a clear error. |
| `reference_name` | string \| null | 否 | Chromosome/contig name to restrict the search to (see list_sequences). 默认 `None` |
| `start` | integer \| null | 否 | Region start position, 1-based inclusive. 默认 `None` |
| `end` | integer \| null | 否 | Region end position, 1-based inclusive. 默认 `None` |
| `min_maf` | number \| null | 否 | Minimum minor-allele frequency (0-1). 默认 `None` |
| `max_maf` | number \| null | 否 | Maximum minor-allele frequency (0-1). 默认 `None` |
| `max_missing_data` | number \| null | 否 | Maximum per-variant missing-data fraction (0-1). 默认 `None` |
| `max_variants` | integer | 否 | Maximum number of matching variants to retrieve. 默认 `100000` |
| `output_dir` | string \| null | 否 | Directory for the output CSV(s) (default ./gigwa_results/&lt;module&gt;/). 默认 `None` |

??? note "上游说明（英文原文）"

    Search variants matching filters server-side and write the matching list to CSV.

    Same filters as ``count_variants`` (region / MAF / missing-data). Returns variant
    metadata only (id, chrom, pos, ref, alt) — no genotypes are fetched — and writes
    ``variant_search.csv``. Use ``count_variants`` first to size the result; ``max_variants``
    caps how many are retrieved. For downstream genotype analysis on a filtered subset, use
    the ``region``/``min_maf`` options on the QC/diversity tools instead.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gigwa-mcp',
    'search_variants',
    {'variant_set_db_id': '...'},
)
```
