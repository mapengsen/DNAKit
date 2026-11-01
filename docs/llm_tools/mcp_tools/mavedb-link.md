# MaveDB-Link (`mavedb-link`)

DNAKit 的 MCP 工具索引收录了这个服务器的 **13** 个 DNA/基因相关工具（服务器共列出 15 个，2026-10-16 读取）。工具在用户自己安装服务器后运行，DNAKit 本身不包含它们。调用方式见 [调用索引中的工具](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found)。下面的工具说明是上游英文原文。

| 项目 | 内容 |
| --- | --- |
| 源码 | [berntpopp/mavedb-link @ `1ae150f`](https://github.com/berntpopp/mavedb-link/tree/1ae150ff3f041045708001403be13f908d58d077) |
| 许可证 | MIT |
| 启动命令 | `python mcp_server.py` |
| 传输方式 | stdio, http |
| 安装 | For stdio, install the source project's dependencies in Python &gt;=3.12 and set workdir to that checkout and executable_path to its Python. The author-documented HTTP service is available with transport='http'. |
| 副作用 | Reads local MaveDB mirrors or contacts the MaveDB API; HTTP mode sends tool arguments to the author's hosted service. |

| 工具 | 类别 | 简介 |
| --- | --- | --- |
| [`find_variant`](#find_variant) | `variant` | Finds one variant (VRS id, MaveDB URN or HGVS) across every MaveDB score set with scores and calibrated classes. |
| [`get_classified_variants`](#get_classified_variants) | `variant` | Returns a score set's variants in a calibrated functional class (abnormal PS3 / normal BS3) with scores and HGVS. |
| [`get_collection`](#get_collection) | `variant` | Returns a curated MaveDB collection by URN with its description and paged member experiment and score-set URNs. |
| [`get_experiment`](#get_experiment) | `variant` | Returns a MaveDB experiment record by URN: title, description, child score sets, keywords and publications. |
| [`get_gene_score_sets`](#get_gene_score_sets) | `variant` | Resolves an HGNC gene symbol to its identity and every published MaveDB score set targeting that gene. |
| [`get_hgvs_validation`](#get_hgvs_validation) | `variant` | Validates an HGVS variant string with MaveDB's validator, returning validity and the upstream reason if invalid. |
| [`get_mapped_variants`](#get_mapped_variants) | `variant` | Returns genome-mapped GA4GH VRS alleles and ClinGen Allele IDs for a score set's variants, paged. |
| [`get_score_distribution`](#get_score_distribution) | `variant` | Summarizes a score set's variant score distribution (quartiles, histogram) and locates a given score's percentile and class. |
| [`get_score_set`](#get_score_set) | `variant` | Returns a MaveDB score-set record by URN: title, target genes with external IDs, publications, license and variant count. |
| [`get_variant_score`](#get_variant_score) | `variant` | Returns the functional score and calibrated class for one variant via a variant URN or score-set URN plus HGVS. |
| [`get_variant_scores`](#get_variant_scores) | `variant` | Returns a score set's paged per-variant table with HGVS (nt/splice/protein), scores and derived classifications. |
| [`search_experiments`](#search_experiments) | `variant` | Searches MaveDB experiments by text and target gene/organism/type/author facets, reranking target-gene matches. |
| [`search_score_sets`](#search_score_sets) | `variant` | Searches MaveDB score sets by free text and target gene, organism, target type and author facets. |

## `find_variant` {#find_variant}

Finds one variant (VRS id, MaveDB URN or HGVS) across every MaveDB score set with scores and calibrated classes.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`variant`
- **关键词：**mave, deep mutational scanning, hgvs, vrs, clingen allele, functional score
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant` | string | 是 | The variant to look up, in ANY of three forms (auto-detected): a GA4GH VRS allele id ('ga4gh:VA…'), a full variant URN ('urn:mavedb:00000001-a-1\#2'), or a bare HGVS string ('p.Asp2723His', 'NM_000059.4:c.8167G&gt;A' — pass gene_symbol= alongside so it can resolve). vrs_id/variant_urn/hgvs are accepted as aliases. |
| `gene_symbol` | string \| null | 否 | HGNC gene symbol that scopes an hgvs= lookup (required when the HGVS is not in the mirror and must be resolved live). Ignored unless hgvs= is set. `gene` is accepted as a compatibility alias. 默认 `None` |
| `only_current` | boolean | 否 | Keep only current genome mappings (default true). 默认 `True` |
| `enrich` | boolean | 否 | Attach each hit's score + classifications (default true). 默认 `True` |
| `limit` | integer | 否 | Max cross-dataset hits (default 25). 默认 `25` |
| `offset` | integer | 否 | Rows to skip for forward paging (default 0). 默认 `0` |
| `response_mode` | string | 否 | Verbosity: minimal\|compact\|standard\|full (default compact). 可选值：`minimal`, `compact`, `standard`, `full` 默认 `compact` |

??? note "上游说明（英文原文）"

    Find ONE variant across EVERY MaveDB score set — the same variant's functional measurements wherever it was assayed (the cross-dataset rollup for 'every assay that measured this variant'). Pass variant= in ANY form (auto-detected): a VRS id (ga4gh:VA…), a variant URN ('urn:mavedb:…-a-1\#2', resolved to its VRS internally so you do NOT map it first — chain straight from get_variant_score), or a bare HGVS string (+ optional gene_symbol=), resolved to its VRS internally — no map-first round-trip. Each hit carries its score_set_urn, variant_urn, ClinGen Allele ID, and (when enrich=true, default) the score + calibrated classifications. ClinGen Allele IDs are not accepted upstream; pass the variant URN instead. Paged via offset/limit. Signature: find_variant(variant, gene_symbol=, only_current=, enrich=, limit=, offset=, response_mode=).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'mavedb-link',
    'find_variant',
    {'variant': '...'},
)
```

## `get_classified_variants` {#get_classified_variants}

Returns a score set's variants in a calibrated functional class (abnormal PS3 / normal BS3) with scores and HGVS.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`urn`
- **关键词：**mave, acmg, ps3, bs3, functional classification, calibration
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `urn` | string | 是 | A MaveDB score-set URN that carries a functional-classification calibration ('urn:mavedb:00000013-a-1'). Uncalibrated sets yield not_found -- confirm via get_score_set (score_calibrations present). |
| `classification` | string \| null | 否 | Filter to one functional class (omit for all). 默认 `None` |
| `calibration_urn` | string \| null | 否 | A specific calibration URN; omit to use the primary calibration. 默认 `None` |
| `limit` | integer | 否 | Max variants (default 100). 默认 `100` |
| `offset` | integer | 否 | Rows to skip for forward paging (default 0). 默认 `0` |
| `response_mode` | string | 否 | Verbosity: minimal\|compact\|standard\|full (default compact). 可选值：`minimal`, `compact`, `standard`, `full` 默认 `compact` |

??? note "上游说明（英文原文）"

    Return a score set's variants grouped into a calibrated functional class — e.g. every 'abnormal' (PS3) or 'normal' (BS3) variant — without paging the whole table. Resolves the primary calibration (or a given calibration_urn); filter with classification=abnormal\|normal\|not_specified (omit for all). Each variant carries its score, HGVS, class label, and ACMG criterion. Paged via offset/limit. Signature: get_classified_variants(urn, classification=, calibration_urn=, limit=, offset=, response_mode=).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'mavedb-link',
    'get_classified_variants',
    {'urn': '...'},
)
```

## `get_collection` {#get_collection}

Returns a curated MaveDB collection by URN with its description and paged member experiment and score-set URNs.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`urn`
- **关键词：**mave, collection, score sets, dataset grouping, consortium
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `urn` | string | 是 | A MaveDB collection URN ('urn:mavedb:collection-&lt;uuid&gt;'). Obtain one from a member score set's official_collections (get_score_set at standard/full). |
| `limit` | integer | 否 | Max member score sets (default 100). 默认 `100` |
| `offset` | integer | 否 | Members to skip for forward paging (default 0). 默认 `0` |
| `response_mode` | string | 否 | Verbosity: minimal\|compact\|standard\|full (default compact). 可选值：`minimal`, `compact`, `standard`, `full` 默认 `compact` |

??? note "上游说明（英文原文）"

    Return a curated MaveDB collection by URN: its name, description, badge, num_score_sets/num_experiments, and the member experiment and score-set URNs. Collections group related datasets (e.g. by gene, consortium, or theme). The member lists are PAGED (limit/offset, truncated/next_offset) so large collections stay light. Open a member with get_score_set. Signature: get_collection(urn, limit=, offset=, response_mode=).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'mavedb-link',
    'get_collection',
    {'urn': '...'},
)
```

## `get_experiment` {#get_experiment}

Returns a MaveDB experiment record by URN: title, description, child score sets, keywords and publications.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`urn`
- **关键词：**mave, experiment, deep mutational scanning, score sets, assay
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `urn` | string | 是 | A MaveDB experiment URN ('urn:mavedb:00000001-a'). Groups one or more score sets; find one via a score set's experiment_urn. |
| `response_mode` | string | 否 | Verbosity: minimal\|compact\|standard\|full (default compact). 可选值：`minimal`, `compact`, `standard`, `full` 默认 `compact` |

??? note "上游说明（英文原文）"

    Return a MaveDB experiment record by URN (urn:mavedb:00000001-a): title, description, parent experiment set, the child score-set URNs, keywords, and publications. An experiment groups one or more score sets from one assay context. Open a child dataset with get_score_set. Signature: get_experiment(urn, response_mode=).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'mavedb-link',
    'get_experiment',
    {'urn': '...'},
)
```

## `get_gene_score_sets` {#get_gene_score_sets}

Resolves an HGNC gene symbol to its identity and every published MaveDB score set targeting that gene.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`gene_symbol`
- **关键词：**mave, hgnc, gene symbol, score sets, variant effect, deep mutational scanning
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene_symbol` | string | 是 | An HGNC gene symbol (e.g. BRCA1, TP53, PTEN). |
| `limit` | integer | 否 | Max score sets (default 20). 默认 `20` |
| `offset` | integer | 否 | Score sets to skip for paging (default 0). 默认 `0` |
| `response_mode` | string | 否 | Verbosity: minimal\|compact\|standard\|full (default compact). 可选值：`minimal`, `compact`, `standard`, `full` 默认 `compact` |

??? note "上游说明（英文原文）"

    Resolve an HGNC gene symbol (e.g. BRCA1) to its gene identity (name, HGNC id, location) AND every published MaveDB score set that targets it — the COMPLETE set, unioned from HGNC resolution and the target-name facet and deduped by URN (see the coverage block). The fastest complete way to find all MAVE data for a gene. Paged via offset/limit. Signature: get_gene_score_sets(gene_symbol, limit=, offset=, response_mode=).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'mavedb-link',
    'get_gene_score_sets',
    {'gene_symbol': '...'},
)
```

## `get_hgvs_validation` {#get_hgvs_validation}

Validates an HGVS variant string with MaveDB's validator, returning validity and the upstream reason if invalid.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`variant`
- **关键词：**hgvs, validation, variant nomenclature, transcript accession, reference base
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant` | string | 是 | An HGVS string (accession-prefixed recommended). |
| `response_mode` | string | 否 | Verbosity: minimal\|compact\|standard\|full (default compact). 可选值：`minimal`, `compact`, `standard`, `full` 默认 `compact` |

??? note "上游说明（英文原文）"

    Validate an HGVS variant string against MaveDB's validator. Returns \{variant, valid, message\}: a valid string -&gt; valid=true; a parseable-but-wrong one -&gt; valid=false WITH the upstream reason (e.g. reference-base disagreement, missing transcript accession) so you can fix it before a lookup fails. Not a normalizer. Signature: get_hgvs_validation(variant, response_mode=).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'mavedb-link',
    'get_hgvs_validation',
    {'variant': '...'},
)
```

## `get_mapped_variants` {#get_mapped_variants}

Returns genome-mapped GA4GH VRS alleles and ClinGen Allele IDs for a score set's variants, paged.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`urn`
- **关键词：**vrs, ga4gh, clingen allele, genome mapping, mave
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `urn` | string | 是 | A MaveDB score-set URN ('urn:mavedb:00000001-a-1'). Find one via search_score_sets or get_gene_score_sets. |
| `current_only` | boolean | 否 | Keep only the current mapping per variant (default true). 默认 `True` |
| `limit` | integer | 否 | Max mapped variants (default 50). 默认 `50` |
| `offset` | integer | 否 | Rows to skip for forward paging (default 0). 默认 `0` |
| `response_mode` | string | 否 | Verbosity: minimal\|compact\|standard\|full (default compact). 可选值：`minimal`, `compact`, `standard`, `full` 默认 `compact` |

??? note "上游说明（英文原文）"

    Return the genome-mapped GA4GH VRS alleles for a score set's variants (urn:mavedb:...-a-1), each with its source variant URN, variant_index, VRS allele id, ClinGen Allele ID, and current flag — the bridge from assay coordinates to reference-genome/clinical coordinates. Rows are ordered numerically by variant_index to match get_variant_scores, but JOIN on variant_urn/variant_index (do not zip by row — some variants are unmapped). current_only (default true) collapses the current/superseded pair to one row per variant. Paged via offset/limit. Signature: get_mapped_variants(urn, current_only=, limit=, offset=, response_mode=).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'mavedb-link',
    'get_mapped_variants',
    {'urn': '...'},
)
```

## `get_score_distribution` {#get_score_distribution}

Summarizes a score set's variant score distribution (quartiles, histogram) and locates a given score's percentile and class.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`urn`
- **关键词：**mave, score distribution, histogram, percentile, calibration thresholds
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `urn` | string | 是 | A MaveDB score-set URN ('urn:mavedb:00000001-a-1'). Find one via search_score_sets or get_gene_score_sets. |
| `score` | number \| null | 否 | A score to locate within the distribution (percentile + class). 默认 `None` |
| `response_mode` | string | 否 | Verbosity: minimal\|compact\|standard\|full (default compact). 可选值：`minimal`, `compact`, `standard`, `full` 默认 `compact` |

??? note "上游说明（英文原文）"

    Summarise a score set's score distribution server-side (MaveDB has no stats endpoint): n, min/max, mean, median, quartiles, stdev, and a 10-bin histogram — a compact summary INSTEAD of paging the whole table. Pass score= to locate that value (its percentile + calibrated classification). Carries the calibration thresholds when present. Signature: get_score_distribution(urn, score=, response_mode=).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'mavedb-link',
    'get_score_distribution',
    {'urn': '...'},
)
```

## `get_score_set` {#get_score_set}

Returns a MaveDB score-set record by URN: title, target genes with external IDs, publications, license and variant count.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`urn`
- **关键词：**mave, score set, target gene, metadata, deep mutational scanning
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `urn` | string | 是 | A MaveDB score-set URN ('urn:mavedb:00000001-a-1'). Find one via search_score_sets or get_gene_score_sets. |
| `response_mode` | string | 否 | Verbosity: minimal\|compact\|standard\|full (default compact). 可选值：`minimal`, `compact`, `standard`, `full` 默认 `compact` |

??? note "上游说明（英文原文）"

    Return a MaveDB score-set record by URN: title, description, target gene(s) with external IDs, parent experiment, publications, the per-record license, variant count, and (standard/full) dataset columns and method/abstract text. Get the actual scores with get_variant_scores. Signature: get_score_set(urn, response_mode=).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'mavedb-link',
    'get_score_set',
    {'urn': '...'},
)
```

## `get_variant_score` {#get_variant_score}

Returns the functional score and calibrated class for one variant via a variant URN or score-set URN plus HGVS.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`urn`
- **关键词：**mave, functional score, hgvs, variant urn, calibrated class
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `urn` | string | 是 | EITHER a full variant URN ('urn:mavedb:00000001-a-1\#2', resolved directly) OR a score-set URN ('urn:mavedb:00000001-a-1', used WITH hgvs=). The first example is a full variant URN so a bare urn= call resolves; pass hgvs= alongside when urn= is a score-set URN. |
| `hgvs` | string \| null | 否 | HGVS string (hgvs_nt or hgvs_pro) to match when urn is a score-set URN; omit when urn is a full variant URN. Accession-prefix insensitive: a bare 'c.8168A&gt;G' resolves a stored 'ENST00000380152.8:c.8168A&gt;G' and vice-versa. Note: many SGE sets leave hgvs_pro null, so a p. form may not match there — use the c. form. 默认 `None` |
| `response_mode` | string | 否 | Verbosity: minimal\|compact\|standard\|full (default compact). 可选值：`minimal`, `compact`, `standard`, `full` 默认 `compact` |

??? note "上游说明（英文原文）"

    Look up the functional score for ONE variant without paging the whole table. Pass a full variant URN (urn:mavedb:...-a-1\#2) to resolve it directly, OR a score-set URN plus hgvs= (e.g. 'c.8168A&gt;G' or 'p.Arg1699Trp') to scan that score set's table for the matching row(s). BOTH paths return the same shape: \{urn, query, resolved_by, match_count, variants:\[\{variant_urn, variant_index, hgvs_nt, hgvs_pro, score, classifications?\}\], calibrations?\}. mapped_variants are opt-in (standard/full) and current-only unless full. The fast path for 'what is the score (and calibrated class) for this variant?'. Signature: get_variant_score(urn, hgvs=, response_mode=).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'mavedb-link',
    'get_variant_score',
    {'urn': '...'},
)
```

## `get_variant_scores` {#get_variant_scores}

Returns a score set's paged per-variant table with HGVS (nt/splice/protein), scores and derived classifications.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`urn`
- **关键词：**mave, score table, hgvs, deep mutational scanning, functional score
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `urn` | string | 是 | A MaveDB score-set URN ('urn:mavedb:00000001-a-1'). Find one via search_score_sets or get_gene_score_sets. |
| `offset` | integer | 否 | Rows to skip in the score table (default 0); `start` is accepted as a start alias for compatibility with MaveDB's upstream scores endpoint. 默认 `0` |
| `limit` | integer | 否 | Max score rows (default 100). 默认 `100` |
| `drop_na_columns` | boolean | 否 | Drop columns that are entirely NA (default false). 默认 `False` |
| `response_mode` | string | 否 | Verbosity: minimal\|compact\|standard\|full (default compact). 可选值：`minimal`, `compact`, `standard`, `full` 默认 `compact` |

??? note "上游说明（英文原文）"

    Return the quantitative variant-by-variant score table for a score set (urn:mavedb:...-a-1) as parsed rows: each row carries accession, variant_index (numeric join key for get_mapped_variants), HGVS (hgvs_nt/hgvs_splice/hgvs_pro), the numeric score (+ derived classification when calibrated), and score-set specific columns. Paged via offset/limit; NA values become null. Page forward with offset=next_offset. For large pulls use response_mode='minimal' to drop the HGVS columns to \{accession, variant_index, score, classification\} (token-safe); a full ~1000-row page at standard can exceed the MCP token cap and be auto-saved to a tool-results file — page via offset or use get_score_distribution for summary statistics instead. The top-level calibrations threshold ladder is record-level data shipped once (the first page, start=0) or at full — not re-sent on every forward page; the per-row classification rides on every page. Signature: get_variant_scores(urn, offset=, limit=, drop_na_columns=, response_mode=).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'mavedb-link',
    'get_variant_scores',
    {'urn': '...'},
)
```

## `search_experiments` {#search_experiments}

Searches MaveDB experiments by text and target gene/organism/type/author facets, reranking target-gene matches.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**无
- **关键词：**mave, experiment search, target gene, organism, author
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `text` | string \| null | 否 | Free-text query over gene/target, title, and abstract. 默认 `None` |
| `targets` | array \| null | 否 | Filter to score sets whose target gene is one of these HGNC symbols (server-side facet). 默认 `None` |
| `target_organism_names` | array \| null | 否 | Filter to these target organisms (client-side, null-inclusive). Use full scientific names. 默认 `None` |
| `target_types` | array \| null | 否 | Filter to these MaveDB target categories (client-side, null-inclusive). Closed set: protein_coding \| regulatory \| other_noncoding (the exact values the runtime accepts; an unlisted value is invalid_input, never a silent-empty result). 默认 `None` |
| `authors` | array \| null | 否 | Filter to score sets with these author name substrings (case-insensitive). 默认 `None` |
| `published` | boolean | 否 | Restrict to published records (default true). 默认 `True` |
| `limit` | integer | 否 | Max hits (default 25). 默认 `25` |
| `offset` | integer | 否 | Rows to skip for forward paging (default 0). 默认 `0` |
| `response_mode` | string | 否 | Verbosity: minimal\|compact\|standard\|full (default compact). 可选值：`minimal`, `compact`, `standard`, `full` 默认 `compact` |

??? note "上游说明（英文原文）"

    Search MaveDB experiments by free text (author facet, plus target facets — targets/target_organism_names/target_types — derived from the score-set search and grouped by parent experiment). A gene-token query re-ranks experiments whose score sets target the gene above abstract namesakes (reranked_by:'target_gene'). Returns experiment hits \{urn, score_set_urns, ...\} plus a pagination block (paged client-side). Use search_score_sets when you want datasets/scores; use this for the assay-context grouping. Signature: search_experiments(text=, targets=, target_organism_names=, target_types=, authors=, published=, limit=, offset=, response_mode=).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'mavedb-link',
    'search_experiments',
    {},
)
```

## `search_score_sets` {#search_score_sets}

Searches MaveDB score sets by free text and target gene, organism, target type and author facets.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**无
- **关键词：**mave, score set search, target gene, deep mutational scanning, facets
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `text` | string \| null | 否 | Free-text query over gene/target, title, and abstract. 默认 `None` |
| `targets` | array \| null | 否 | Filter to score sets whose target gene is one of these HGNC symbols (server-side facet). 默认 `None` |
| `target_organism_names` | array \| null | 否 | Filter to these target organisms (client-side, null-inclusive). Use full scientific names. 默认 `None` |
| `target_types` | array \| null | 否 | Filter to these MaveDB target categories (client-side, null-inclusive). Closed set: protein_coding \| regulatory \| other_noncoding (the exact values the runtime accepts; an unlisted value is invalid_input, never a silent-empty result). 默认 `None` |
| `authors` | array \| null | 否 | Filter to score sets with these author name substrings (case-insensitive). 默认 `None` |
| `facet_mode` | string | 否 | 'inclusive' (default; keep unknown-metadata records) or 'strict' (drop them). 可选值：`inclusive`, `strict` 默认 `inclusive` |
| `published` | boolean | 否 | Restrict to published records (default true). 默认 `True` |
| `limit` | integer | 否 | Max hits (default 25). 默认 `25` |
| `offset` | integer | 否 | Rows to skip for forward paging (default 0). 默认 `0` |
| `response_mode` | string | 否 | Verbosity: minimal\|compact\|standard\|full (default compact). 可选值：`minimal`, `compact`, `standard`, `full` 默认 `compact` |

??? note "上游说明（英文原文）"

    Search MaveDB score sets (the datasets carrying scored variants) by free text and facets — target gene(s), organism, target type, and author. Returns score-set hits \{urn, title, num_variants, targets, license, ...\} plus a pagination block. Organism/target-type facets are null-inclusive by default (records with unknown metadata are KEPT and _meta.facet_excluded reports the drops); pass facet_mode='strict' to also drop unknown-metadata records. This is the MaveDB front door. Signature: search_score_sets(text=, targets=, target_organism_names=, target_types=, authors=, facet_mode=, published=, limit=, offset=, response_mode=).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'mavedb-link',
    'search_score_sets',
    {},
)
```
