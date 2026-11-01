# Locus (`locus`)

DNAKit's MCP tool index keeps **21** DNA/gene tools of this server (22 listed, read on 2026-10-16). They run once you install the server; DNAKit does not bundle them. See [calling a tool that was found](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found).

| Item | Value |
| --- | --- |
| Source | [mattsmcnulty/locus @ `24b0465`](https://github.com/mattsmcnulty/locus/tree/24b046548e9e976fc7e28395b86ff262e6d3c343) |
| License | MIT |
| Command | `locus-mcp` |
| Transports | stdio |
| Installation | Install separately into a Python &gt;=3.11 environment from the pinned source (pip install git+https://github.com/mattsmcnulty/locus.git@24b046548e9e976fc7e28395b86ff262e6d3c343; not the unrelated PyPI package named locus), then build the local DuckDB genome store with the upstream pipeline (locus download all; locus pipeline); set LOCUS_DATA_DIR (or LOCUS_DB_PATH) to the data directory because the default location is derived from the source tree. |
| Required environment | `LOCUS_DATA_DIR` |
| Effects | Most tools run read-only queries against the local DuckDB store; ask_about and variants_in_study send rsIDs to the Ensembl REST API and write temporary marker files under the Locus work directory, variants_in_study queries the GWAS Catalog, and literature_for sends the query to NCBI PubMed E-utilities. |

| Tool | Category | Summary |
| --- | --- | --- |
| [`allele_frequency`](#allele_frequency) | `variant` | Returns gnomAD allele frequencies for ClinVar- or AlphaMissense-flagged variants at a genomic position or region. |
| [`ancestry`](#ancestry) | `variant` | Returns estimated continental ancestry proportions and PCA placement of the loaded genome versus 1000 Genomes and HGDP. |
| [`ask_about`](#ask_about) | `variant` | Reports this genome's genotypes at listed rsIDs or its carried GWAS associations for a named condition or trait. |
| [`carrier_status`](#carrier_status) | `variant` | Reports which curated recessive conditions the genome carries one pathogenic copy of, plus a not-assessed list. |
| [`clinical_findings`](#clinical_findings) | `variant` | Lists ClinVar pathogenic/likely-pathogenic or other-significance variants in the loaded genome, filterable by gene. |
| [`gene_disease_validity`](#gene_disease_validity) | `association` | Returns ClinGen expert-curated gene-disease validity classifications (Definitive to Refuted) for a gene symbol. |
| [`genome_overview`](#genome_overview) | `variant` | Summarizes the loaded genome: variant counts, genome build and ClinVar/gnomAD/pharmacogenomics annotation coverage. |
| [`gwas_associations`](#gwas_associations) | `association` | Lists genome-wide-significant GWAS Catalog risk alleles carried by this genome, optionally filtered by trait. |
| [`lookup_variant_by_rsid`](#lookup_variant_by_rsid) | `variant` | Returns the genome's genotype and ClinVar/gnomAD annotations for a dbSNP rsID. |
| [`lookup_variants_in_gene`](#lookup_variants_in_gene) | `variant` | Lists the genome's variants within a gene by HGNC symbol, paginated, using SnpEff annotation. |
| [`lookup_variants_in_region`](#lookup_variants_in_region) | `variant` | Lists the genome's variants in a chr:start-end region or at a single position, paginated. |
| [`pharmacogenomics`](#pharmacogenomics) | `variant` | Returns PharmCAT star-allele diplotypes, metabolizer phenotypes and CPIC/DPWG drug guidance, filterable by gene or drug. |
| [`polygenic_risk`](#polygenic_risk) | `variant` | Returns ancestry-matched polygenic risk percentiles for a curated trait set (CAD, LDL, T2D, Lp(a)). |
| [`predicted_damaging`](#predicted_damaging) | `variant` | Lists rare missense variants unclassified by ClinVar that AlphaMissense predicts pathogenic, optionally by gene. |
| [`run_sql`](#run_sql) | `variant` | Runs a read-only SELECT over the genome database tables (variants, cnv, sv, pgx, associations, pgs_scores), capped at 200 rows. |
| [`secondary_findings`](#secondary_findings) | `variant` | Reports pathogenic/likely-pathogenic ClinVar variants in the 84 ACMG SF v3.3 medically actionable genes. |
| [`structural_variants`](#structural_variants) | `variant` | Returns copy-number and structural variant events overlapping a region, with estimated copy number. |
| [`traits`](#traits) | `variant` | Reports single-SNP wellness traits, the HLA-B*57:01 proxy and mtDNA maternal haplogroup, filterable by category. |
| [`variant_dossier`](#variant_dossier) | `variant` | Returns all evidence for one rsID in this genome: zygosity, consequence, ClinVar, AlphaMissense, gnomAD, ACMG and GWAS hits. |
| [`variants_in_study`](#variants_in_study) | `association` | Pulls the variants a GWAS study (by PMID) reported via the GWAS Catalog and genotypes them in this genome. |
| [`whats_new`](#whats_new) | `variant` | Returns the changelog since the last refresh, e.g. carried ClinVar variants newly reclassified, filterable by date or tier. |

## `allele_frequency` {#allele_frequency}

Returns gnomAD allele frequencies for ClinVar- or AlphaMissense-flagged variants at a genomic position or region.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `region`
- **Keywords:** gnomad, allele frequency, rare variant, genomic region, clinvar

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `region` | string | yes |  |

??? note "Upstream description"

    How common is a variant? Returns gnomAD allele frequencies for variants at a position/region
    ('chr1:55052000' or 'chr1:55050000-55060000'). Low gnomad_af = rare. Frequencies are carried for
    the variants where rarity is informative — those ClinVar has classified or AlphaMissense calls
    pathogenic; other variants return a null gnomad_af rather than a frequency.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'locus',
    'allele_frequency',
    {'region': '...'},
)
```

## `ancestry` {#ancestry}

Returns estimated continental ancestry proportions and PCA placement of the loaded genome versus 1000 Genomes and HGDP.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** none
- **Keywords:** ancestry, 1000 genomes, hgdp, pca, admixture, personal genome

??? note "Upstream description"

    Estimated biogeographic ancestry: continental proportions (k-NN over 1000 Genomes + HGDP) and the
    PCA placement. Continental ancestry is reliable; finer/admixed breakdowns are approximate.
    Requires `locus ancestry` to have run.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'locus',
    'ancestry',
    {},
)
```

## `ask_about` {#ask_about}

Reports this genome's genotypes at listed rsIDs or its carried GWAS associations for a named condition or trait.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `query`
- **Keywords:** rsid, genotype, gwas, trait, personal genome, clinvar

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string | yes |  |
| `limit` | integer | no | Default `50` |

??? note "Upstream description"

    On-demand: given a condition/trait name OR a paper's rsIDs (space/comma-separated), report what
    THIS genome carries. rsIDs → live genotype lookup at those positions (hom-ref-aware) plus any
    ClinVar/AlphaMissense we already have; a trait name → your carried GWAS associations for it. This is
    the 'keep up with new studies' tool: paste rsIDs from a new paper to see your genotypes. Weak
    single-variant evidence — informational, not diagnostic.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'locus',
    'ask_about',
    {'query': '...'},
)
```

## `carrier_status` {#carrier_status}

Reports which curated recessive conditions the genome carries one pathogenic copy of, plus a not-assessed list.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** none
- **Keywords:** carrier screening, recessive, pathogenic, family planning, heterozygous

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `limit` | integer | no | Default `100` |

??? note "Upstream description"

    Carrier status for common recessive conditions — which ones this person carries ONE
    pathogenic copy of. This is about their children, not their own health: a carrier is
    typically unaffected. It matters when both partners carry the same condition (1-in-4 risk per
    pregnancy), so frame results around family planning and a genetic counselor, never as a
    personal diagnosis. 'likely_affected' = two pathogenic copies (phase unknown from a VCF, so a
    compound-het is presumed, not proven).

    IMPORTANT: always report the `not_assessed` list alongside any result. This is a curated panel
    of VCF-assessable conditions, not a clinical carrier screen — and several of the most
    important screens (SMN1/spinal muscular atrophy, FMR1/Fragile X) CANNOT be answered from this
    data at all. An empty `hits` list is NOT a negative carrier screen; do not present it as one.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'locus',
    'carrier_status',
    {},
)
```

## `clinical_findings` {#clinical_findings}

Lists ClinVar pathogenic/likely-pathogenic or other-significance variants in the loaded genome, filterable by gene.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** none
- **Keywords:** clinvar, pathogenic, clinical significance, personal genome, hgnc

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene` | string | no | Default `` |
| `significance` | string | no | Default `` |
| `limit` | integer | no | Default `100` |
| `offset` | integer | no | Default `0` |

??? note "Upstream description"

    ClinVar clinical findings in this genome. With no args, returns pathogenic / likely-pathogenic
    variants. Optionally filter by `gene` (HGNC symbol) and/or `significance` substring
    (e.g. 'pathogenic', 'risk_factor', 'drug_response'). Always confirm health-relevant hits clinically.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'locus',
    'clinical_findings',
    {},
)
```

## `gene_disease_validity` {#gene_disease_validity}

Returns ClinGen expert-curated gene-disease validity classifications (Definitive to Refuted) for a gene symbol.

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `gene`
- **Keywords:** clingen, gene-disease validity, definitive, curation, mendelian

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene` | string | yes |  |

??? note "Upstream description"

    Is a GENE itself established as causing disease? Returns ClinGen's expert-curated
    gene-disease validity classifications for a gene (Definitive → Strong → Moderate → Limited →
    Disputed → Refuted), strongest first. Use this to judge whether a gene *matters* — distinct
    from whether a specific variant is pathogenic (that's `clinical_findings`/`variant_dossier`).

    A "Definitive" link means the gene-disease relationship is real; it does NOT mean the person's
    particular variant is pathogenic. An empty result means ClinGen hasn't curated this gene, not
    that the gene is benign. Reads a local file — nothing about the genome leaves the machine.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'locus',
    'gene_disease_validity',
    {'gene': '...'},
)
```

## `genome_overview` {#genome_overview}

Summarizes the loaded genome: variant counts, genome build and ClinVar/gnomAD/pharmacogenomics annotation coverage.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** none
- **Keywords:** personal genome, variant counts, genome build, annotation summary

??? note "Upstream description"

    Summarize the loaded genome: variant counts, genome build, and how much is
    annotated (ClinVar, gnomAD, pharmacogenomics). Call this first to see what's available.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'locus',
    'genome_overview',
    {},
)
```

## `gwas_associations` {#gwas_associations}

Lists genome-wide-significant GWAS Catalog risk alleles carried by this genome, optionally filtered by trait.

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** none
- **Keywords:** gwas catalog, risk allele, trait, p-value, personal genome

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `trait` | string | no | Default `` |
| `limit` | integer | no | Default `100` |

??? note "Upstream description"

    GWAS Catalog risk alleles THIS genome carries (genome-wide significant, p&lt;5e-8), optionally
    filtered by `trait` substring (e.g. 'type 2 diabetes', 'height'). WEAK / EXPLORATORY: single hits
    with tiny effects — do NOT sum them or read them like a calibrated score (use polygenic_risk for
    that). Ordered by significance. Requires `locus gwas` to have run.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'locus',
    'gwas_associations',
    {},
)
```

## `lookup_variant_by_rsid` {#lookup_variant_by_rsid}

Returns the genome's genotype and ClinVar/gnomAD annotations for a dbSNP rsID.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `rsid`
- **Keywords:** dbsnp, rsid, genotype, clinvar, gnomad

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `rsid` | string | yes |  |

??? note "Upstream description"

    Look up variant(s) by dbSNP rsID (e.g. 'rs1799853'). Returns the genotype and
    any ClinVar/gnomAD annotations carried for that variant.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'locus',
    'lookup_variant_by_rsid',
    {'rsid': '...'},
)
```

## `lookup_variants_in_gene` {#lookup_variants_in_gene}

Lists the genome's variants within a gene by HGNC symbol, paginated, using SnpEff annotation.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `gene`
- **Keywords:** hgnc, gene variants, snpeff, personal genome, pagination

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene` | string | yes |  |
| `limit` | integer | no | Default `100` |
| `offset` | integer | no | Default `0` |

??? note "Upstream description"

    List variants in a gene by HGNC symbol (e.g. 'BRCA1'). Paginated — pass offset to page.
    `total` is the full match count. Requires the SnpEff annotation step to have run.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'locus',
    'lookup_variants_in_gene',
    {'gene': '...'},
)
```

## `lookup_variants_in_region` {#lookup_variants_in_region}

Lists the genome's variants in a chr:start-end region or at a single position, paginated.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `region`
- **Keywords:** genomic region, coordinates, personal genome, vcf, interval

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `region` | string | yes |  |
| `limit` | integer | no | Default `200` |
| `offset` | integer | no | Default `0` |

??? note "Upstream description"

    List variants in a genomic region. `region` is 'chr7:117480000-117670000' (1-based, chr-prefixed)
    or a single position 'chr7:117559590'. Paginated.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'locus',
    'lookup_variants_in_region',
    {'region': '...'},
)
```

## `pharmacogenomics` {#pharmacogenomics}

Returns PharmCAT star-allele diplotypes, metabolizer phenotypes and CPIC/DPWG drug guidance, filterable by gene or drug.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** none
- **Keywords:** pharmcat, star allele, diplotype, metabolizer, cpic, cyp2c19

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene` | string | no | Default `` |
| `drug` | string | no | Default `` |

??? note "Upstream description"

    Pharmacogenomic results from PharmCAT: star-allele diplotypes, metabolizer phenotypes, and
    CPIC/DPWG drug guidance. Filter by `gene` (e.g. 'CYP2C19') or `drug` (e.g. 'clopidogrel').
    Some genes report 'Unknown'/'No Result' — notably CYP2D6, whose calling needs copy-number data
    from raw reads that Locus does not process. Treat those as "not determined", not "normal".
    Requires the PharmCAT annotation step to have run.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'locus',
    'pharmacogenomics',
    {},
)
```

## `polygenic_risk` {#polygenic_risk}

Returns ancestry-matched polygenic risk percentiles for a curated trait set (CAD, LDL, T2D, Lp(a)).

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** none
- **Keywords:** polygenic risk score, percentile, coronary artery disease, type 2 diabetes, ancestry

??? note "Upstream description"

    Polygenic (aggregate) risk scores for the curated set (coronary artery disease, LDL
    cholesterol, type 2 diabetes, Lp(a)) — call it to see what is actually scored, reported as
    an ancestry-matched percentile where available. Percentiles are only meaningful within the matched
    ancestry; these are research-grade estimates, not diagnoses. Requires `locus ancestry` to have run.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'locus',
    'polygenic_risk',
    {},
)
```

## `predicted_damaging` {#predicted_damaging}

Lists rare missense variants unclassified by ClinVar that AlphaMissense predicts pathogenic, optionally by gene.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** none
- **Keywords:** alphamissense, missense, pathogenicity prediction, rare variant, vus

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene` | string | no | Default `` |
| `limit` | integer | no | Default `100` |

??? note "Upstream description"

    Rare, predicted-damaging missense variants from AlphaMissense — the 'ClinVar is silent' set:
    variants ClinVar has NOT classified but AlphaMissense scores as likely pathogenic. Optionally
    filter by gene. Use this when ClinVar returns nothing for a gene of interest. These are
    computational predictions, not clinical assertions — far weaker evidence than a ClinVar
    pathogenic call, and most are benign in reality. Rarity filtering is best-effort: gnomAD
    frequencies are only carried for some variants, so a hit is not guaranteed to be rare.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'locus',
    'predicted_damaging',
    {},
)
```

## `run_sql` {#run_sql}

Runs a read-only SELECT over the genome database tables (variants, cnv, sv, pgx, associations, pgs_scores), capped at 200 rows.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `query`
- **Keywords:** sql, duckdb, variants table, cnv, power query, personal genome

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string | yes |  |

??? note "Upstream description"

    Run a read-only SELECT against the genome database for power queries. Tables: variants
    (chrom,pos,ref,alt,rsid,gt,filter,gene,consequence,clnsig,clndn,clnrevstat,gnomad_af,gnomad_af_grpmax,
    am_pathogenicity,am_class), cnv, sv, pgx_genes, pgx_drugs, meta, traits, associations (carried GWAS
    hits), pgs_scores, ancestry_global, ancestry_pca, watch_findings (the refresh changelog),
    watch_seen_ids, sources. Mutating statements are rejected; results capped at 200 rows.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'locus',
    'run_sql',
    {'query': '...'},
)
```

## `secondary_findings` {#secondary_findings}

Reports pathogenic/likely-pathogenic ClinVar variants in the 84 ACMG SF v3.3 medically actionable genes.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** none
- **Keywords:** acmg, secondary findings, actionable genes, pathogenic, clinvar

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `limit` | integer | no | Default `100` |

??? note "Upstream description"

    ACMG SF v3.3 secondary (incidental) findings: pathogenic / likely-pathogenic ClinVar variants
    in the 84 medically-actionable genes ACMG recommends reporting (cancer, cardiac, metabolic).
    ACMG's recessive genes are reported only when two P/LP variants are present, so single carriers
    are excluded by design. Empty is the common, reassuring result — but it means "no P/LP ClinVar
    variant in these 84 genes", NOT a clean bill of health: it cannot see variants ClinVar hasn't
    classified, non-coding or structural variants, or any gene off this list. Confirm hits clinically.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'locus',
    'secondary_findings',
    {},
)
```

## `structural_variants` {#structural_variants}

Returns copy-number and structural variant events overlapping a region, with estimated copy number.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `region`
- **Keywords:** cnv, structural variant, copy number, deletion, duplication

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `region` | string | yes |  |
| `limit` | integer | no | Default `100` |

??? note "Upstream description"

    Copy-number (CNV) and structural (SV) events overlapping a region ('chr1:1000000-2000000').
    For CNVs, `cn` is the estimated copy number (2 = normal diploid).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'locus',
    'structural_variants',
    {'region': '...'},
)
```

## `traits` {#traits}

Reports single-SNP wellness traits, the HLA-B*57:01 proxy and mtDNA maternal haplogroup, filterable by category.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** none
- **Keywords:** snp traits, haplogroup, mtdna, hla-b*57:01, lactose, wellness

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `category` | string | no | Default `` |

??? note "Upstream description"

    Single-SNP traits/wellness (lactose, caffeine, alcohol flush, earwax, eye color, muscle type),
    the HLA-B*57:01 abacavir-hypersensitivity screening proxy, and the mtDNA maternal-lineage
    haplogroup. Filter by `category`: 'wellness' \| 'pharmacogenomic' \| 'maternal lineage' (the
    haplogroup lives under the last one). Call with no filter to see everything. Informational,
    not diagnostic. Requires `locus traits` to have run.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'locus',
    'traits',
    {},
)
```

## `variant_dossier` {#variant_dossier}

Returns all evidence for one rsID in this genome: zygosity, consequence, ClinVar, AlphaMissense, gnomAD, ACMG and GWAS hits.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `rsid`
- **Keywords:** rsid, zygosity, clinvar, alphamissense, gnomad, variant interpretation

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `rsid` | string | yes |  |

??? note "Upstream description"

    Everything this genome knows about ONE variant, in a single call: genotype/zygosity, gene,
    consequence, the ClinVar call and its review status, the AlphaMissense prediction, gnomAD
    frequency, whether the gene is ACMG-actionable or on the carrier panel, any carried GWAS
    associations, and a literature link.

    Prefer this over assembling a picture from separate lookups. The evidence types routinely
    DISAGREE and only mean something together — an AlphaMissense "pathogenic" call looks alarming
    until you see ClinVar calls it benign and 40% of people carry it. Precedence: ClinVar (with its
    review status) outranks AlphaMissense, which is a prediction and never overrides a benign call;
    a high gnomAD frequency argues strongly against a severe effect regardless of prediction. If
    the fields conflict, say so rather than reporting the scariest one. A variant not found is
    almost certainly homozygous-reference, not 'ruled out'.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'locus',
    'variant_dossier',
    {'rsid': '...'},
)
```

## `variants_in_study` {#variants_in_study}

Pulls the variants a GWAS study (by PMID) reported via the GWAS Catalog and genotypes them in this genome.

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `pmid`
- **Keywords:** gwas catalog, pmid, study variants, genotype, risk allele

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `pmid` | string | yes |  |

??? note "Upstream description"

    Given the PubMed ID of a GWAS study, pull the variants it reported (via the GWAS Catalog) and
    show which ones THIS genome carries — live, hom-ref-aware genotyping. Answers 'this new paper
    found N variants, which do I have?'. Weak single-variant evidence — informational, not diagnostic.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'locus',
    'variants_in_study',
    {'pmid': '...'},
)
```

## `whats_new` {#whats_new}

Returns the changelog since the last refresh, e.g. carried ClinVar variants newly reclassified, filterable by date or tier.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** none
- **Keywords:** clinvar reclassification, changelog, refresh, personal genome, pathogenic

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `since` | string | no | Default `` |
| `tier` | string | no | Default `` |

??? note "Upstream description"

    What changed about THIS genome since the last `locus refresh` — the deterministic
    changelog (e.g. ClinVar variants you carry that were newly classified pathogenic or
    reclassified). Ranked strongest-first. Optionally filter by `since` (ISO date like
    '2026-06-01') or `tier` ('strong'\|'moderate'\|'weak'\|'info'). Empty if refresh hasn't run.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'locus',
    'whats_new',
    {},
)
```
