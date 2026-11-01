# MaveDB-Link (`mavedb-link`)

DNAKit's MCP tool index keeps **13** DNA/gene tools of this server (15 listed, read on 2026-10-16). They run once you install the server; DNAKit does not bundle them. See [calling a tool that was found](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found).

| Item | Value |
| --- | --- |
| Source | [berntpopp/mavedb-link @ `1ae150f`](https://github.com/berntpopp/mavedb-link/tree/1ae150ff3f041045708001403be13f908d58d077) |
| License | MIT |
| Command | `python mcp_server.py` |
| Transports | stdio, http |
| Installation | For stdio, install the source project's dependencies in Python &gt;=3.12 and set workdir to that checkout and executable_path to its Python. The author-documented HTTP service is available with transport='http'. |
| Effects | Reads local MaveDB mirrors or contacts the MaveDB API; HTTP mode sends tool arguments to the author's hosted service. |

| Tool | Category | Summary |
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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `variant`
- **Keywords:** mave, deep mutational scanning, hgvs, vrs, clingen allele, functional score
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `variant` | string | yes | The variant to look up, in ANY of three forms (auto-detected): a GA4GH VRS allele id ('ga4gh:VA…'), a full variant URN ('urn:mavedb:00000001-a-1\#2'), or a bare HGVS string ('p.Asp2723His', 'NM_000059.4:c.8167G&gt;A' — pass gene_symbol= alongside so it can resolve). vrs_id/variant_urn/hgvs are accepted as aliases. |
| `gene_symbol` | string \| null | no | HGNC gene symbol that scopes an hgvs= lookup (required when the HGVS is not in the mirror and must be resolved live). Ignored unless hgvs= is set. `gene` is accepted as a compatibility alias. Default `None` |
| `only_current` | boolean | no | Keep only current genome mappings (default true). Default `True` |
| `enrich` | boolean | no | Attach each hit's score + classifications (default true). Default `True` |
| `limit` | integer | no | Max cross-dataset hits (default 25). Default `25` |
| `offset` | integer | no | Rows to skip for forward paging (default 0). Default `0` |
| `response_mode` | string | no | Verbosity: minimal\|compact\|standard\|full (default compact). One of `minimal`, `compact`, `standard`, `full` Default `compact` |

??? note "Upstream description"

    Find ONE variant across EVERY MaveDB score set — the same variant's functional measurements wherever it was assayed (the cross-dataset rollup for 'every assay that measured this variant'). Pass variant= in ANY form (auto-detected): a VRS id (ga4gh:VA…), a variant URN ('urn:mavedb:…-a-1\#2', resolved to its VRS internally so you do NOT map it first — chain straight from get_variant_score), or a bare HGVS string (+ optional gene_symbol=), resolved to its VRS internally — no map-first round-trip. Each hit carries its score_set_urn, variant_urn, ClinGen Allele ID, and (when enrich=true, default) the score + calibrated classifications. ClinGen Allele IDs are not accepted upstream; pass the variant URN instead. Paged via offset/limit. Signature: find_variant(variant, gene_symbol=, only_current=, enrich=, limit=, offset=, response_mode=).

Call it through DNAKit once the server is installed:

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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `urn`
- **Keywords:** mave, acmg, ps3, bs3, functional classification, calibration
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `urn` | string | yes | A MaveDB score-set URN that carries a functional-classification calibration ('urn:mavedb:00000013-a-1'). Uncalibrated sets yield not_found -- confirm via get_score_set (score_calibrations present). |
| `classification` | string \| null | no | Filter to one functional class (omit for all). Default `None` |
| `calibration_urn` | string \| null | no | A specific calibration URN; omit to use the primary calibration. Default `None` |
| `limit` | integer | no | Max variants (default 100). Default `100` |
| `offset` | integer | no | Rows to skip for forward paging (default 0). Default `0` |
| `response_mode` | string | no | Verbosity: minimal\|compact\|standard\|full (default compact). One of `minimal`, `compact`, `standard`, `full` Default `compact` |

??? note "Upstream description"

    Return a score set's variants grouped into a calibrated functional class — e.g. every 'abnormal' (PS3) or 'normal' (BS3) variant — without paging the whole table. Resolves the primary calibration (or a given calibration_urn); filter with classification=abnormal\|normal\|not_specified (omit for all). Each variant carries its score, HGVS, class label, and ACMG criterion. Paged via offset/limit. Signature: get_classified_variants(urn, classification=, calibration_urn=, limit=, offset=, response_mode=).

Call it through DNAKit once the server is installed:

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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `urn`
- **Keywords:** mave, collection, score sets, dataset grouping, consortium
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `urn` | string | yes | A MaveDB collection URN ('urn:mavedb:collection-&lt;uuid&gt;'). Obtain one from a member score set's official_collections (get_score_set at standard/full). |
| `limit` | integer | no | Max member score sets (default 100). Default `100` |
| `offset` | integer | no | Members to skip for forward paging (default 0). Default `0` |
| `response_mode` | string | no | Verbosity: minimal\|compact\|standard\|full (default compact). One of `minimal`, `compact`, `standard`, `full` Default `compact` |

??? note "Upstream description"

    Return a curated MaveDB collection by URN: its name, description, badge, num_score_sets/num_experiments, and the member experiment and score-set URNs. Collections group related datasets (e.g. by gene, consortium, or theme). The member lists are PAGED (limit/offset, truncated/next_offset) so large collections stay light. Open a member with get_score_set. Signature: get_collection(urn, limit=, offset=, response_mode=).

Call it through DNAKit once the server is installed:

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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `urn`
- **Keywords:** mave, experiment, deep mutational scanning, score sets, assay
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `urn` | string | yes | A MaveDB experiment URN ('urn:mavedb:00000001-a'). Groups one or more score sets; find one via a score set's experiment_urn. |
| `response_mode` | string | no | Verbosity: minimal\|compact\|standard\|full (default compact). One of `minimal`, `compact`, `standard`, `full` Default `compact` |

??? note "Upstream description"

    Return a MaveDB experiment record by URN (urn:mavedb:00000001-a): title, description, parent experiment set, the child score-set URNs, keywords, and publications. An experiment groups one or more score sets from one assay context. Open a child dataset with get_score_set. Signature: get_experiment(urn, response_mode=).

Call it through DNAKit once the server is installed:

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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `gene_symbol`
- **Keywords:** mave, hgnc, gene symbol, score sets, variant effect, deep mutational scanning
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_symbol` | string | yes | An HGNC gene symbol (e.g. BRCA1, TP53, PTEN). |
| `limit` | integer | no | Max score sets (default 20). Default `20` |
| `offset` | integer | no | Score sets to skip for paging (default 0). Default `0` |
| `response_mode` | string | no | Verbosity: minimal\|compact\|standard\|full (default compact). One of `minimal`, `compact`, `standard`, `full` Default `compact` |

??? note "Upstream description"

    Resolve an HGNC gene symbol (e.g. BRCA1) to its gene identity (name, HGNC id, location) AND every published MaveDB score set that targets it — the COMPLETE set, unioned from HGNC resolution and the target-name facet and deduped by URN (see the coverage block). The fastest complete way to find all MAVE data for a gene. Paged via offset/limit. Signature: get_gene_score_sets(gene_symbol, limit=, offset=, response_mode=).

Call it through DNAKit once the server is installed:

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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `variant`
- **Keywords:** hgvs, validation, variant nomenclature, transcript accession, reference base
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `variant` | string | yes | An HGVS string (accession-prefixed recommended). |
| `response_mode` | string | no | Verbosity: minimal\|compact\|standard\|full (default compact). One of `minimal`, `compact`, `standard`, `full` Default `compact` |

??? note "Upstream description"

    Validate an HGVS variant string against MaveDB's validator. Returns \{variant, valid, message\}: a valid string -&gt; valid=true; a parseable-but-wrong one -&gt; valid=false WITH the upstream reason (e.g. reference-base disagreement, missing transcript accession) so you can fix it before a lookup fails. Not a normalizer. Signature: get_hgvs_validation(variant, response_mode=).

Call it through DNAKit once the server is installed:

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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `urn`
- **Keywords:** vrs, ga4gh, clingen allele, genome mapping, mave
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `urn` | string | yes | A MaveDB score-set URN ('urn:mavedb:00000001-a-1'). Find one via search_score_sets or get_gene_score_sets. |
| `current_only` | boolean | no | Keep only the current mapping per variant (default true). Default `True` |
| `limit` | integer | no | Max mapped variants (default 50). Default `50` |
| `offset` | integer | no | Rows to skip for forward paging (default 0). Default `0` |
| `response_mode` | string | no | Verbosity: minimal\|compact\|standard\|full (default compact). One of `minimal`, `compact`, `standard`, `full` Default `compact` |

??? note "Upstream description"

    Return the genome-mapped GA4GH VRS alleles for a score set's variants (urn:mavedb:...-a-1), each with its source variant URN, variant_index, VRS allele id, ClinGen Allele ID, and current flag — the bridge from assay coordinates to reference-genome/clinical coordinates. Rows are ordered numerically by variant_index to match get_variant_scores, but JOIN on variant_urn/variant_index (do not zip by row — some variants are unmapped). current_only (default true) collapses the current/superseded pair to one row per variant. Paged via offset/limit. Signature: get_mapped_variants(urn, current_only=, limit=, offset=, response_mode=).

Call it through DNAKit once the server is installed:

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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `urn`
- **Keywords:** mave, score distribution, histogram, percentile, calibration thresholds
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `urn` | string | yes | A MaveDB score-set URN ('urn:mavedb:00000001-a-1'). Find one via search_score_sets or get_gene_score_sets. |
| `score` | number \| null | no | A score to locate within the distribution (percentile + class). Default `None` |
| `response_mode` | string | no | Verbosity: minimal\|compact\|standard\|full (default compact). One of `minimal`, `compact`, `standard`, `full` Default `compact` |

??? note "Upstream description"

    Summarise a score set's score distribution server-side (MaveDB has no stats endpoint): n, min/max, mean, median, quartiles, stdev, and a 10-bin histogram — a compact summary INSTEAD of paging the whole table. Pass score= to locate that value (its percentile + calibrated classification). Carries the calibration thresholds when present. Signature: get_score_distribution(urn, score=, response_mode=).

Call it through DNAKit once the server is installed:

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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `urn`
- **Keywords:** mave, score set, target gene, metadata, deep mutational scanning
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `urn` | string | yes | A MaveDB score-set URN ('urn:mavedb:00000001-a-1'). Find one via search_score_sets or get_gene_score_sets. |
| `response_mode` | string | no | Verbosity: minimal\|compact\|standard\|full (default compact). One of `minimal`, `compact`, `standard`, `full` Default `compact` |

??? note "Upstream description"

    Return a MaveDB score-set record by URN: title, description, target gene(s) with external IDs, parent experiment, publications, the per-record license, variant count, and (standard/full) dataset columns and method/abstract text. Get the actual scores with get_variant_scores. Signature: get_score_set(urn, response_mode=).

Call it through DNAKit once the server is installed:

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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `urn`
- **Keywords:** mave, functional score, hgvs, variant urn, calibrated class
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `urn` | string | yes | EITHER a full variant URN ('urn:mavedb:00000001-a-1\#2', resolved directly) OR a score-set URN ('urn:mavedb:00000001-a-1', used WITH hgvs=). The first example is a full variant URN so a bare urn= call resolves; pass hgvs= alongside when urn= is a score-set URN. |
| `hgvs` | string \| null | no | HGVS string (hgvs_nt or hgvs_pro) to match when urn is a score-set URN; omit when urn is a full variant URN. Accession-prefix insensitive: a bare 'c.8168A&gt;G' resolves a stored 'ENST00000380152.8:c.8168A&gt;G' and vice-versa. Note: many SGE sets leave hgvs_pro null, so a p. form may not match there — use the c. form. Default `None` |
| `response_mode` | string | no | Verbosity: minimal\|compact\|standard\|full (default compact). One of `minimal`, `compact`, `standard`, `full` Default `compact` |

??? note "Upstream description"

    Look up the functional score for ONE variant without paging the whole table. Pass a full variant URN (urn:mavedb:...-a-1\#2) to resolve it directly, OR a score-set URN plus hgvs= (e.g. 'c.8168A&gt;G' or 'p.Arg1699Trp') to scan that score set's table for the matching row(s). BOTH paths return the same shape: \{urn, query, resolved_by, match_count, variants:\[\{variant_urn, variant_index, hgvs_nt, hgvs_pro, score, classifications?\}\], calibrations?\}. mapped_variants are opt-in (standard/full) and current-only unless full. The fast path for 'what is the score (and calibrated class) for this variant?'. Signature: get_variant_score(urn, hgvs=, response_mode=).

Call it through DNAKit once the server is installed:

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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `urn`
- **Keywords:** mave, score table, hgvs, deep mutational scanning, functional score
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `urn` | string | yes | A MaveDB score-set URN ('urn:mavedb:00000001-a-1'). Find one via search_score_sets or get_gene_score_sets. |
| `offset` | integer | no | Rows to skip in the score table (default 0); `start` is accepted as a start alias for compatibility with MaveDB's upstream scores endpoint. Default `0` |
| `limit` | integer | no | Max score rows (default 100). Default `100` |
| `drop_na_columns` | boolean | no | Drop columns that are entirely NA (default false). Default `False` |
| `response_mode` | string | no | Verbosity: minimal\|compact\|standard\|full (default compact). One of `minimal`, `compact`, `standard`, `full` Default `compact` |

??? note "Upstream description"

    Return the quantitative variant-by-variant score table for a score set (urn:mavedb:...-a-1) as parsed rows: each row carries accession, variant_index (numeric join key for get_mapped_variants), HGVS (hgvs_nt/hgvs_splice/hgvs_pro), the numeric score (+ derived classification when calibrated), and score-set specific columns. Paged via offset/limit; NA values become null. Page forward with offset=next_offset. For large pulls use response_mode='minimal' to drop the HGVS columns to \{accession, variant_index, score, classification\} (token-safe); a full ~1000-row page at standard can exceed the MCP token cap and be auto-saved to a tool-results file — page via offset or use get_score_distribution for summary statistics instead. The top-level calibrations threshold ladder is record-level data shipped once (the first page, start=0) or at full — not re-sent on every forward page; the per-row classification rides on every page. Signature: get_variant_scores(urn, offset=, limit=, drop_na_columns=, response_mode=).

Call it through DNAKit once the server is installed:

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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** none
- **Keywords:** mave, experiment search, target gene, organism, author
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `text` | string \| null | no | Free-text query over gene/target, title, and abstract. Default `None` |
| `targets` | array \| null | no | Filter to score sets whose target gene is one of these HGNC symbols (server-side facet). Default `None` |
| `target_organism_names` | array \| null | no | Filter to these target organisms (client-side, null-inclusive). Use full scientific names. Default `None` |
| `target_types` | array \| null | no | Filter to these MaveDB target categories (client-side, null-inclusive). Closed set: protein_coding \| regulatory \| other_noncoding (the exact values the runtime accepts; an unlisted value is invalid_input, never a silent-empty result). Default `None` |
| `authors` | array \| null | no | Filter to score sets with these author name substrings (case-insensitive). Default `None` |
| `published` | boolean | no | Restrict to published records (default true). Default `True` |
| `limit` | integer | no | Max hits (default 25). Default `25` |
| `offset` | integer | no | Rows to skip for forward paging (default 0). Default `0` |
| `response_mode` | string | no | Verbosity: minimal\|compact\|standard\|full (default compact). One of `minimal`, `compact`, `standard`, `full` Default `compact` |

??? note "Upstream description"

    Search MaveDB experiments by free text (author facet, plus target facets — targets/target_organism_names/target_types — derived from the score-set search and grouped by parent experiment). A gene-token query re-ranks experiments whose score sets target the gene above abstract namesakes (reranked_by:'target_gene'). Returns experiment hits \{urn, score_set_urns, ...\} plus a pagination block (paged client-side). Use search_score_sets when you want datasets/scores; use this for the assay-context grouping. Signature: search_experiments(text=, targets=, target_organism_names=, target_types=, authors=, published=, limit=, offset=, response_mode=).

Call it through DNAKit once the server is installed:

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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** none
- **Keywords:** mave, score set search, target gene, deep mutational scanning, facets
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `text` | string \| null | no | Free-text query over gene/target, title, and abstract. Default `None` |
| `targets` | array \| null | no | Filter to score sets whose target gene is one of these HGNC symbols (server-side facet). Default `None` |
| `target_organism_names` | array \| null | no | Filter to these target organisms (client-side, null-inclusive). Use full scientific names. Default `None` |
| `target_types` | array \| null | no | Filter to these MaveDB target categories (client-side, null-inclusive). Closed set: protein_coding \| regulatory \| other_noncoding (the exact values the runtime accepts; an unlisted value is invalid_input, never a silent-empty result). Default `None` |
| `authors` | array \| null | no | Filter to score sets with these author name substrings (case-insensitive). Default `None` |
| `facet_mode` | string | no | 'inclusive' (default; keep unknown-metadata records) or 'strict' (drop them). One of `inclusive`, `strict` Default `inclusive` |
| `published` | boolean | no | Restrict to published records (default true). Default `True` |
| `limit` | integer | no | Max hits (default 25). Default `25` |
| `offset` | integer | no | Rows to skip for forward paging (default 0). Default `0` |
| `response_mode` | string | no | Verbosity: minimal\|compact\|standard\|full (default compact). One of `minimal`, `compact`, `standard`, `full` Default `compact` |

??? note "Upstream description"

    Search MaveDB score sets (the datasets carrying scored variants) by free text and facets — target gene(s), organism, target type, and author. Returns score-set hits \{urn, title, num_variants, targets, license, ...\} plus a pagination block. Organism/target-type facets are null-inclusive by default (records with unknown metadata are KEPT and _meta.facet_excluded reports the drops); pass facet_mode='strict' to also drop unknown-metadata records. This is the MaveDB front door. Signature: search_score_sets(text=, targets=, target_organism_names=, target_types=, authors=, facet_mode=, published=, limit=, offset=, response_mode=).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'mavedb-link',
    'search_score_sets',
    {},
)
```
