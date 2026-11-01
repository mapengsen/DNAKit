# Open Targets MCP (`opentargets`)

DNAKit's MCP tool index keeps **30** DNA/gene tools of this server (68 listed, read on 2026-10-16). They run once you install the server; DNAKit does not bundle them. See [calling a tool that was found](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found).

| Item | Value |
| --- | --- |
| Source | [nickzren/opentargets-mcp @ `ecb92f8`](https://github.com/nickzren/opentargets-mcp/tree/ecb92f8effdd6f8ea8f6a8c8de2fedb23334c126) |
| License | MIT |
| Command | `python -m opentargets_mcp.server --transport stdio` |
| Transports | stdio |
| Installation | Install nickzren/opentargets-mcp separately in Python &gt;=3.10; set executable_path to its interpreter. Uses the verified Python module entry point with explicit stdio. |
| Effects | Queries Open Targets evidence; results are research evidence and do not establish treatment suitability. |

| Tool | Category | Summary |
| --- | --- | --- |
| [`get_credible_set_by_id`](#get_credible_set_by_id) | `association` | Fetches details of an Open Targets study-locus credible set by ID, including lead variant and fine-mapped variants. |
| [`get_credible_set_colocalisation`](#get_credible_set_colocalisation) | `association` | Retrieves colocalisation results (e.g. with eQTL or other GWAS) for a study-locus credible set. |
| [`get_credible_sets`](#get_credible_sets) | `association` | Queries Open Targets credible sets filtered by study-locus IDs, study IDs, variant IDs, study types or genomic regions. |
| [`get_disease_associated_targets`](#get_disease_associated_targets) | `association` | Lists targets (genes) associated with a disease (EFO ID) with overall and per-datatype association scores. |
| [`get_drug_pharmacogenomics`](#get_drug_pharmacogenomics) | `variant` | Retrieves pharmacogenomics evidence linking genetic variants and genotypes to response for a ChEMBL drug. |
| [`get_similar_targets`](#get_similar_targets) | `association` | Identifies targets (genes) whose disease-association profiles are similar to a seed target. |
| [`get_studies_by_disease`](#get_studies_by_disease) | `association` | Lists GWAS and molecular QTL studies linked to one or more diseases, optionally including indirect ontology children. |
| [`get_study_credible_sets`](#get_study_credible_sets) | `association` | Retrieves fine-mapped credible sets (loci) for a GWAS or QTL study. |
| [`get_study_info`](#get_study_info) | `association` | Retrieves metadata, trait, cohort and sample details for a GWAS study. |
| [`get_target_alternative_genes`](#get_target_alternative_genes) | `genome_annotation` | Lists alternate gene symbols, synonyms and database cross-references for an Ensembl gene target. |
| [`get_target_associated_diseases`](#get_target_associated_diseases) | `association` | Lists diseases linked to an Ensembl gene target with association scores. |
| [`get_target_depmap_essentiality`](#get_target_depmap_essentiality) | `genome_annotation` | Fetches DepMap CRISPR knockout gene essentiality scores across cancer cell lines for a gene. |
| [`get_target_disease_biomarkers`](#get_target_disease_biomarkers) | `association` | Inspects evidence for biomarker annotations (e.g. variants, expression) linking a target gene to a disease. |
| [`get_target_disease_evidence`](#get_target_disease_evidence) | `association` | Retrieves evidence strings linking an Ensembl target to an EFO disease, filterable by data source. |
| [`get_target_expression`](#get_target_expression) | `expression` | Returns baseline RNA and protein expression profiles for a target gene across tissues. |
| [`get_target_genetic_constraint`](#get_target_genetic_constraint) | `variant` | Fetches gnomAD genetic constraint metrics (e.g. pLI, LOEUF, missense Z) for a target gene. |
| [`get_target_hallmarks`](#get_target_hallmarks) | `genome_annotation` | Returns COSMIC cancer hallmark annotations and attributes associated with a target gene. |
| [`get_target_homologues`](#get_target_homologues) | `genome_annotation` | Retrieves cross-species orthologues and paralogues for a target gene. |
| [`get_target_info`](#get_target_info) | `genome_annotation` | Retrieves core identity details for a target gene: symbol, name, biotype and genomic location. |
| [`get_target_mouse_phenotypes`](#get_target_mouse_phenotypes) | `association` | Retrieves mouse model knockout phenotypes (MGI/IMPC) associated with a target gene. |
| [`get_target_pathways_and_go_terms`](#get_target_pathways_and_go_terms) | `genome_annotation` | Lists Reactome pathway memberships and Gene Ontology annotations for a target gene. |
| [`get_targets_batch`](#get_targets_batch) | `genome_annotation` | Retrieves core records for multiple targets by Ensembl gene IDs in one request. |
| [`get_variant_credible_sets`](#get_variant_credible_sets) | `association` | Lists fine-mapped credible sets that include a specific variant, filterable by study type. |
| [`get_variant_evidences`](#get_variant_evidences) | `association` | Retrieves evidence strings linking a variant to target genes or diseases, filterable by data source. |
| [`get_variant_info`](#get_variant_info) | `variant` | Retrieves core metadata, allele frequencies and functional consequence annotations for a variant. |
| [`get_variant_intervals`](#get_variant_intervals) | `regulation_epigenome` | Retrieves enhancer-to-gene (E2G) regulatory interval predictions overlapping a variant. |
| [`get_variant_pharmacogenomics`](#get_variant_pharmacogenomics) | `variant` | Retrieves pharmacogenomic annotations (drug response) for a variant. |
| [`get_variant_protein_coordinates`](#get_variant_protein_coordinates) | `variant` | Retrieves protein-level consequences and coordinates for a variant (amino acid change, position). |
| [`graphql_batch_query`](#graphql_batch_query) | `association` | Executes one Open Targets GraphQL query against many variable sets concurrently, keyed by a field. |
| [`graphql_query`](#graphql_query) | `association` | Executes a raw GraphQL query with variables against the Open Targets Platform API. |

## `get_credible_set_by_id` {#get_credible_set_by_id}

Fetches details of an Open Targets study-locus credible set by ID, including lead variant and fine-mapped variants.

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `study_locus_id`
- **Keywords:** fine-mapping, credible set, gwas, study locus, l2g
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `study_locus_id` | string | yes |  |

??? note "Upstream description"

    Fetch detailed information for a specific study locus credible set.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets',
    'get_credible_set_by_id',
    {'study_locus_id': '...'},
)
```

## `get_credible_set_colocalisation` {#get_credible_set_colocalisation}

Retrieves colocalisation results (e.g. with eQTL or other GWAS) for a study-locus credible set.

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `study_locus_id`
- **Keywords:** colocalisation, coloc, eqtl, credible set, gwas
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `study_locus_id` | string | yes |  |
| `page_index` | integer | no | Default `0` |
| `page_size` | integer | no | Default `10` |

??? note "Upstream description"

    Retrieve colocalisation results for a credible set.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets',
    'get_credible_set_colocalisation',
    {'study_locus_id': '...'},
)
```

## `get_credible_sets` {#get_credible_sets}

Queries Open Targets credible sets filtered by study-locus IDs, study IDs, variant IDs, study types or genomic regions.

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** none
- **Keywords:** fine-mapping, credible set, gwas, qtl, genomic region
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `study_locus_ids` | array \| null | no | Default `None` |
| `study_ids` | array \| null | no | Default `None` |
| `variant_ids` | array \| null | no | Default `None` |
| `study_types` | array \| null | no | Default `None` |
| `regions` | array \| null | no | Default `None` |
| `page_index` | integer | no | Default `0` |
| `page_size` | integer | no | Default `10` |
| `fields` | array \| null | no | Default `None` |

??? note "Upstream description"

    Query credible sets with flexible filtering options.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets',
    'get_credible_sets',
    {},
)
```

## `get_disease_associated_targets` {#get_disease_associated_targets}

Lists targets (genes) associated with a disease (EFO ID) with overall and per-datatype association scores.

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `efo_id`
- **Keywords:** target-disease, efo, association score, genes, genetic evidence
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `efo_id` | string | yes |  |
| `fields` | array \| null | no | Default `None` |
| `page_index` | integer | no | Default `0` |
| `page_size` | integer | no | Default `10` |

??? note "Upstream description"

    List targets associated with a disease, including evidence scores.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets',
    'get_disease_associated_targets',
    {'efo_id': '...'},
)
```

## `get_drug_pharmacogenomics` {#get_drug_pharmacogenomics}

Retrieves pharmacogenomics evidence linking genetic variants and genotypes to response for a ChEMBL drug.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `chembl_id`
- **Keywords:** pharmacogenomics, pharmgkb, genotype, drug response, chembl
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `chembl_id` | string | yes |  |
| `page_index` | integer | no | Default `0` |
| `page_size` | integer | no | Default `10` |

??? note "Upstream description"

    Retrieve pharmacogenomics data linking genetic variants to drug response.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets',
    'get_drug_pharmacogenomics',
    {'chembl_id': '...'},
)
```

## `get_similar_targets` {#get_similar_targets}

Identifies targets (genes) whose disease-association profiles are similar to a seed target.

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `entity_id`
- **Keywords:** similar genes, association profile, target, embeddings
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `entity_id` | string | yes |  |
| `threshold` | number \| null | no | Default `0.5` |
| `size` | integer | no | Default `10` |

??? note "Upstream description"

    Identify targets with similar association profiles to the seed target.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets',
    'get_similar_targets',
    {'entity_id': '...'},
)
```

## `get_studies_by_disease` {#get_studies_by_disease}

Lists GWAS and molecular QTL studies linked to one or more diseases, optionally including indirect ontology children.

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `disease_ids`
- **Keywords:** gwas, study, efo, qtl, gwas catalog
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `disease_ids` | array of string | yes |  |
| `enable_indirect` | boolean | no | Default `False` |
| `study_id` | string \| null | no | Default `None` |
| `page_index` | integer | no | Default `0` |
| `page_size` | integer | no | Default `10` |
| `fields` | array \| null | no | Default `None` |

??? note "Upstream description"

    List studies linked to one or more diseases.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets',
    'get_studies_by_disease',
    {'disease_ids': []},
)
```

## `get_study_credible_sets` {#get_study_credible_sets}

Retrieves fine-mapped credible sets (loci) for a GWAS or QTL study.

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `study_id`
- **Keywords:** fine-mapping, credible set, gwas, loci
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `study_id` | string | yes |  |
| `page_index` | integer | no | Default `0` |
| `page_size` | integer | no | Default `10` |
| `fields` | array \| null | no | Default `None` |

??? note "Upstream description"

    Retrieve fine-mapped credible sets for a study.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets',
    'get_study_credible_sets',
    {'study_id': '...'},
)
```

## `get_study_info` {#get_study_info}

Retrieves metadata, trait, cohort and sample details for a GWAS study.

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `study_id`
- **Keywords:** gwas, cohort, sample size, trait, ancestry
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `study_id` | string | yes |  |
| `fields` | array \| null | no | Default `None` |

??? note "Upstream description"

    Retrieve metadata and cohort details for a GWAS study.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets',
    'get_study_info',
    {'study_id': '...'},
)
```

## `get_target_alternative_genes` {#get_target_alternative_genes}

Lists alternate gene symbols, synonyms and database cross-references for an Ensembl gene target.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `ensembl_id`
- **Keywords:** gene symbol, synonyms, xrefs, ensembl id
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `ensembl_id` | string | yes |  |

??? note "Upstream description"

    List alternate gene symbols and database cross-references for a target.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets',
    'get_target_alternative_genes',
    {'ensembl_id': '...'},
)
```

## `get_target_associated_diseases` {#get_target_associated_diseases}

Lists diseases linked to an Ensembl gene target with association scores.

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `ensembl_id`
- **Keywords:** target-disease, association score, efo, genetic evidence
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `ensembl_id` | string | yes |  |
| `fields` | array \| null | no | Default `None` |
| `page_index` | integer | no | Default `0` |
| `page_size` | integer | no | Default `10` |

??? note "Upstream description"

    List diseases linked to a target with association scores.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets',
    'get_target_associated_diseases',
    {'ensembl_id': '...'},
)
```

## `get_target_depmap_essentiality` {#get_target_depmap_essentiality}

Fetches DepMap CRISPR knockout gene essentiality scores across cancer cell lines for a gene.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `ensembl_id`
- **Keywords:** depmap, crispr screen, gene essentiality, cell lines, dependency
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `ensembl_id` | string | yes |  |

??? note "Upstream description"

    Fetch DepMap CRISPR essentiality scores across cell lines.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets',
    'get_target_depmap_essentiality',
    {'ensembl_id': '...'},
)
```

## `get_target_disease_biomarkers` {#get_target_disease_biomarkers}

Inspects evidence for biomarker annotations (e.g. variants, expression) linking a target gene to a disease.

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `ensembl_id`, `efo_id`
- **Keywords:** biomarker, target-disease, evidence, efo
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `ensembl_id` | string | yes |  |
| `efo_id` | string | yes |  |
| `size` | integer | no | Default `10` |
| `cursor` | string \| null | no | Default `None` |
| `fields` | array \| null | no | Default `None` |

??? note "Upstream description"

    Inspect evidence for biomarker annotations linking a target and disease.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets',
    'get_target_disease_biomarkers',
    {'ensembl_id': '...', 'efo_id': '...'},
)
```

## `get_target_disease_evidence` {#get_target_disease_evidence}

Retrieves evidence strings linking an Ensembl target to an EFO disease, filterable by data source.

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `ensembl_id`, `efo_id`
- **Keywords:** evidence, target-disease, genetic association, datasource
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `ensembl_id` | string | yes |  |
| `efo_id` | string | yes |  |
| `datasource_ids` | array \| null | no | Default `None` |
| `size` | integer | no | Default `10` |
| `cursor` | string \| null | no | Default `None` |
| `fields` | array \| null | no | Default `None` |

??? note "Upstream description"

    Retrieve evidence strings linking a target to a disease.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets',
    'get_target_disease_evidence',
    {'ensembl_id': '...', 'efo_id': '...'},
)
```

## `get_target_expression` {#get_target_expression}

Returns baseline RNA and protein expression profiles for a target gene across tissues.

- **Category:** `expression` (Gene-level expression and transcript abundance)
- **Required arguments:** `ensembl_id`
- **Keywords:** tissue expression, rna, gtex, hpa, baseline
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `ensembl_id` | string | yes |  |
| `fields` | array \| null | no | Default `None` |

??? note "Upstream description"

    Return RNA and protein expression profiles for a target across tissues.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets',
    'get_target_expression',
    {'ensembl_id': '...'},
)
```

## `get_target_genetic_constraint` {#get_target_genetic_constraint}

Fetches gnomAD genetic constraint metrics (e.g. pLI, LOEUF, missense Z) for a target gene.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `ensembl_id`
- **Keywords:** gnomad, pli, loeuf, loss-of-function, constraint
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `ensembl_id` | string | yes |  |

??? note "Upstream description"

    Fetch genetic constraint metrics (gnomAD) for a target.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets',
    'get_target_genetic_constraint',
    {'ensembl_id': '...'},
)
```

## `get_target_hallmarks` {#get_target_hallmarks}

Returns COSMIC cancer hallmark annotations and attributes associated with a target gene.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `ensembl_id`
- **Keywords:** cancer hallmarks, cosmic, gene function, oncogene
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `ensembl_id` | string | yes |  |

??? note "Upstream description"

    Return cancer hallmark annotations associated with a target.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets',
    'get_target_hallmarks',
    {'ensembl_id': '...'},
)
```

## `get_target_homologues` {#get_target_homologues}

Retrieves cross-species orthologues and paralogues for a target gene.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `ensembl_id`
- **Keywords:** orthologs, paralogs, homology, species
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `ensembl_id` | string | yes |  |

??? note "Upstream description"

    Retrieve cross-species homologues for a target gene.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets',
    'get_target_homologues',
    {'ensembl_id': '...'},
)
```

## `get_target_info` {#get_target_info}

Retrieves core identity details for a target gene: symbol, name, biotype and genomic location.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `ensembl_id`
- **Keywords:** gene symbol, biotype, genomic location, ensembl id
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `ensembl_id` | string | yes |  |
| `fields` | array \| null | no | Default `None` |

??? note "Upstream description"

    Retrieve core identity details for a target gene.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets',
    'get_target_info',
    {'ensembl_id': '...'},
)
```

## `get_target_mouse_phenotypes` {#get_target_mouse_phenotypes}

Retrieves mouse model knockout phenotypes (MGI/IMPC) associated with a target gene.

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `ensembl_id`
- **Keywords:** mouse knockout, mgi, impc, phenotype, model organism
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `ensembl_id` | string | yes |  |
| `page_index` | integer | no | Default `0` |
| `page_size` | integer | no | Default `10` |

??? note "Upstream description"

    Retrieve mouse knockout phenotypes associated with a target.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets',
    'get_target_mouse_phenotypes',
    {'ensembl_id': '...'},
)
```

## `get_target_pathways_and_go_terms` {#get_target_pathways_and_go_terms}

Lists Reactome pathway memberships and Gene Ontology annotations for a target gene.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `ensembl_id`
- **Keywords:** reactome, go annotations, pathways, gene function
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `ensembl_id` | string | yes |  |
| `page_index` | integer | no | Default `0` |
| `page_size` | integer | no | Default `10` |

??? note "Upstream description"

    List pathway memberships and Gene Ontology annotations for a target.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets',
    'get_target_pathways_and_go_terms',
    {'ensembl_id': '...'},
)
```

## `get_targets_batch` {#get_targets_batch}

Retrieves core records for multiple targets by Ensembl gene IDs in one request.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `ensembl_ids`
- **Keywords:** ensembl id, genes, batch, target info
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `ensembl_ids` | array of string | yes |  |

??? note "Upstream description"

    Retrieve multiple targets in a single request.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets',
    'get_targets_batch',
    {'ensembl_ids': []},
)
```

## `get_variant_credible_sets` {#get_variant_credible_sets}

Lists fine-mapped credible sets that include a specific variant, filterable by study type.

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `variant_id`
- **Keywords:** fine-mapping, credible set, gwas, qtl, variant
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `variant_id` | string | yes |  |
| `study_types` | array \| null | no | Default `None` |
| `page_index` | integer | no | Default `0` |
| `page_size` | integer | no | Default `10` |

??? note "Upstream description"

    List credible sets that include a specific variant.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets',
    'get_variant_credible_sets',
    {'variant_id': '...'},
)
```

## `get_variant_evidences` {#get_variant_evidences}

Retrieves evidence strings linking a variant to target genes or diseases, filterable by data source.

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `variant_id`
- **Keywords:** variant-to-gene, evidence, variant-disease, l2g
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `variant_id` | string | yes |  |
| `datasource_ids` | array \| null | no | Default `None` |
| `size` | integer | no | Default `10` |
| `cursor` | string \| null | no | Default `None` |
| `fields` | array \| null | no | Default `None` |

??? note "Upstream description"

    Retrieve evidence strings linking a variant to targets or diseases.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets',
    'get_variant_evidences',
    {'variant_id': '...'},
)
```

## `get_variant_info` {#get_variant_info}

Retrieves core metadata, allele frequencies and functional consequence annotations for a variant.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `variant_id`
- **Keywords:** variant annotation, rsid, allele frequency, vep, consequence
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `variant_id` | string | yes |  |
| `fields` | array \| null | no | Default `None` |

??? note "Upstream description"

    Retrieve core metadata and functional annotations for a variant.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets',
    'get_variant_info',
    {'variant_id': '...'},
)
```

## `get_variant_intervals` {#get_variant_intervals}

Retrieves enhancer-to-gene (E2G) regulatory interval predictions overlapping a variant.

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** `variant_id`
- **Keywords:** enhancer, e2g, regulatory element, variant-to-gene
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `variant_id` | string | yes |  |
| `page_index` | integer | no | Default `0` |
| `page_size` | integer | no | Default `10` |

??? note "Upstream description"

    Retrieve enhancer-to-gene (E2G) predictions overlapping a variant.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets',
    'get_variant_intervals',
    {'variant_id': '...'},
)
```

## `get_variant_pharmacogenomics` {#get_variant_pharmacogenomics}

Retrieves pharmacogenomic annotations (drug response) for a variant.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `variant_id`
- **Keywords:** pharmacogenomics, pharmgkb, drug response, genotype
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `variant_id` | string | yes |  |
| `page_index` | integer | no | Default `0` |
| `page_size` | integer | no | Default `10` |

??? note "Upstream description"

    Retrieve pharmacogenomic annotations for a variant.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets',
    'get_variant_pharmacogenomics',
    {'variant_id': '...'},
)
```

## `get_variant_protein_coordinates` {#get_variant_protein_coordinates}

Retrieves protein-level consequences and coordinates for a variant (amino acid change, position).

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `variant_id`
- **Keywords:** amino acid change, protein consequence, missense, variant effect
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `variant_id` | string | yes |  |
| `page_index` | integer | no | Default `0` |
| `page_size` | integer | no | Default `10` |

??? note "Upstream description"

    Retrieve protein-level consequences for a variant.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets',
    'get_variant_protein_coordinates',
    {'variant_id': '...'},
)
```

## `graphql_batch_query` {#graphql_batch_query}

Executes one Open Targets GraphQL query against many variable sets concurrently, keyed by a field.

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `query_string`, `variables_list`
- **Keywords:** graphql, batch query, target-disease, genetic evidence
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query_string` | string | yes |  |
| `variables_list` | array of object | yes |  |
| `key_field` | string \| null | no | Default `None` |
| `operation_name` | string \| null | no | Default `None` |
| `max_concurrency` | integer | no | Default `4` |

??? note "Upstream description"

    Execute one GraphQL query against many variable sets.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets',
    'graphql_batch_query',
    {'query_string': '...', 'variables_list': []},
)
```

## `graphql_query` {#graphql_query}

Executes a raw GraphQL query with variables against the Open Targets Platform API.

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `query_string`
- **Keywords:** graphql, raw query, target-disease, genetic evidence
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query_string` | string | yes |  |
| `variables` | object \| null | no | Default `None` |
| `operation_name` | string \| null | no | Default `None` |

??? note "Upstream description"

    ADVANCED: Execute a raw GraphQL query against Open Targets.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets',
    'graphql_query',
    {'query_string': '...'},
)
```
