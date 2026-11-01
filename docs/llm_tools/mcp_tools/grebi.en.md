# GrEBI (`grebi`)

DNAKit's MCP tool index keeps **15** DNA/gene tools of this server (33 listed, read on 2026-10-16). They run once you install the server; DNAKit does not bundle them. See [calling a tool that was found](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found).

| Item | Value |
| --- | --- |
| Source | [EBISPOT/GrEBI @ `74aef79`](https://github.com/EBISPOT/GrEBI/tree/74aef794aab9b3d0e7a52d01c15c48a86d043e09) |
| License | Apache-2.0 |
| Command | — |
| Transports | http |
| Installation | Select transport='http' for the author-documented development service. Local deployment is a separate Java/API plus graph-database and HPC data-pipeline task; no local startup recipe is claimed here. |
| Effects | Sends queries to the EMBL-EBI development knowledge-graph service. |

| Tool | Category | Summary |
| --- | --- | --- |
| [`disease_to_genes`](#disease_to_genes) | `association` | Returns genes associated with a disease ID (e.g. MONDO) from Monarch and GWAS Catalog data in the GrEBI graph, with source datasets. |
| [`gene_to_diseases`](#gene_to_diseases) | `association` | Returns diseases associated with a gene ID from Monarch and GWAS Catalog data in the GrEBI graph, with source datasets. |
| [`gwas_by_cell_type`](#gwas_by_cell_type) | `association` | Returns GWAS Catalog associations (SNP, gene, allele, p-value, study) for traits linked to a Cell Ontology cell type. |
| [`gwas_by_disease_location`](#gwas_by_disease_location) | `association` | Returns GWAS Catalog associations (SNP, gene, allele, p-value, study) for diseases located in a given anatomical location. |
| [`gwas_by_drugs_indicated_for_disease`](#gwas_by_drugs_indicated_for_disease) | `association` | Returns GWAS Catalog associations (SNP, gene, allele, p-value) for traits involving drugs indicated for a disease. |
| [`gwas_by_gene_and_cell_type`](#gwas_by_gene_and_cell_type) | `association` | Returns GWAS Catalog associations mapped to a gene for traits linked to a given cell type. |
| [`gwas_by_gene_and_disease`](#gwas_by_gene_and_disease) | `association` | Returns GWAS Catalog associations mapped to a gene for traits within a given disease context. |
| [`gwas_by_gene_and_location`](#gwas_by_gene_and_location) | `association` | Returns GWAS Catalog associations mapped to a gene for traits linked to an anatomical location. |
| [`gwas_by_location`](#gwas_by_location) | `association` | Returns GWAS Catalog associations (SNP, gene, allele, p-value, study) for traits linked to an anatomical location or its parts. |
| [`gwas_by_pathway`](#gwas_by_pathway) | `association` | Returns GWAS Catalog SNP associations for genes in a Reactome pathway, with alleles, traits and p-values. |
| [`gwas_studies_by_cell_type`](#gwas_studies_by_cell_type) | `association` | Returns GWAS Catalog studies annotated with traits linked to a cell type, with sample sizes, association counts and summary statistics availability. |
| [`gwas_studies_by_disease_location`](#gwas_studies_by_disease_location) | `association` | Returns GWAS Catalog studies for diseases located in an anatomical location, with study metadata. |
| [`gwas_studies_by_drugs_indicated_for_disease`](#gwas_studies_by_drugs_indicated_for_disease) | `association` | Returns GWAS Catalog studies whose traits involve drugs indicated for a disease, with study metadata. |
| [`gwas_trait_to_mouse_models_via_embeddings`](#gwas_trait_to_mouse_models_via_embeddings) | `association` | Finds mouse gene models with phenotypes semantically similar to a GWAS trait whose human orthologs have GWAS associations with that trait. |
| [`mouse_gene_to_opentargets`](#mouse_gene_to_opentargets) | `association` | Returns Open Targets disease associations, with scores, for the human ortholog of a mouse gene. |

## `disease_to_genes` {#disease_to_genes}

Returns genes associated with a disease ID (e.g. MONDO) from Monarch and GWAS Catalog data in the GrEBI graph, with source datasets.

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `graph`, `disease_id`, `sortBy`, `sortDir`, `pageNum`, `pageSize`
- **Keywords:** gene-disease association, mondo, monarch, gwas catalog, disease genes

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `graph` | any | yes | One of `ebi_monarch_xspecies` |
| `disease_id` | string | yes | Disease |
| `sortBy` | any | yes | One of `disease`, `gene`, `from_datasources` |
| `sortDir` | any | yes | One of `asc`, `desc` |
| `pageNum` | integer | yes | Page number (0-based) |
| `pageSize` | integer | yes | Number of results per page |

??? note "Upstream description"

    Search for gene disease associations

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'grebi',
    'disease_to_genes',
    {'graph': 'ebi_monarch_xspecies',
     'disease_id': '...',
     'sortBy': 'disease',
     'sortDir': 'asc',
     'pageNum': 0,
     'pageSize': 0},
)
```

## `gene_to_diseases` {#gene_to_diseases}

Returns diseases associated with a gene ID from Monarch and GWAS Catalog data in the GrEBI graph, with source datasets.

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `graph`, `gene_id`, `sortBy`, `sortDir`, `pageNum`, `pageSize`
- **Keywords:** gene-disease association, hgnc, monarch, gwas catalog, disease

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `graph` | any | yes | One of `ebi_monarch_xspecies` |
| `gene_id` | string | yes | Gene |
| `sortBy` | any | yes | One of `disease`, `gene`, `from_datasources` |
| `sortDir` | any | yes | One of `asc`, `desc` |
| `pageNum` | integer | yes | Page number (0-based) |
| `pageSize` | integer | yes | Number of results per page |

??? note "Upstream description"

    Search for diseases associated with a gene

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'grebi',
    'gene_to_diseases',
    {'graph': 'ebi_monarch_xspecies',
     'gene_id': '...',
     'sortBy': 'disease',
     'sortDir': 'asc',
     'pageNum': 0,
     'pageSize': 0},
)
```

## `gwas_by_cell_type` {#gwas_by_cell_type}

Returns GWAS Catalog associations (SNP, gene, allele, p-value, study) for traits linked to a Cell Ontology cell type.

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `graph`, `cell_type_id`, `sortBy`, `sortDir`, `pageNum`, `pageSize`
- **Keywords:** gwas catalog, snp, cell ontology, trait, p-value, mapped gene

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `graph` | any | yes | One of `ebi_monarch_xspecies` |
| `cell_type_id` | string | yes | Cell Type |
| `sortBy` | any | yes | One of `cell_type`, `trait`, `snp`, `effect_allele`, `gene`, `study_accession`, `p_value`, `pubmed_id`, `reported_trait`, `chromosome`, `position`, `risk_allele_frequency` |
| `sortDir` | any | yes | One of `asc`, `desc` |
| `pageNum` | integer | yes | Page number (0-based) |
| `pageSize` | integer | yes | Number of results per page |

??? note "Upstream description"

    Search GWAS Catalog by cell type

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'grebi',
    'gwas_by_cell_type',
    {'graph': 'ebi_monarch_xspecies',
     'cell_type_id': '...',
     'sortBy': 'cell_type',
     'sortDir': 'asc',
     'pageNum': 0,
     'pageSize': 0},
)
```

## `gwas_by_disease_location` {#gwas_by_disease_location}

Returns GWAS Catalog associations (SNP, gene, allele, p-value, study) for diseases located in a given anatomical location.

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `graph`, `location_id`, `sortBy`, `sortDir`, `pageNum`, `pageSize`
- **Keywords:** gwas catalog, snp, anatomy, disease location, p-value

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `graph` | any | yes | One of `ebi_monarch_xspecies`, `impc_x_gwas` |
| `location_id` | string | yes | Anatomical Location |
| `sortBy` | any | yes | One of `location`, `disease`, `snp`, `effect_allele`, `gene`, `study_accession`, `p_value`, `pubmed_id`, `reported_trait`, `chromosome`, `position`, `risk_allele_frequency` |
| `sortDir` | any | yes | One of `asc`, `desc` |
| `pageNum` | integer | yes | Page number (0-based) |
| `pageSize` | integer | yes | Number of results per page |

??? note "Upstream description"

    Search GWAS Catalog by disease location

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'grebi',
    'gwas_by_disease_location',
    {'graph': 'ebi_monarch_xspecies',
     'location_id': '...',
     'sortBy': 'location',
     'sortDir': 'asc',
     'pageNum': 0,
     'pageSize': 0},
)
```

## `gwas_by_drugs_indicated_for_disease` {#gwas_by_drugs_indicated_for_disease}

Returns GWAS Catalog associations (SNP, gene, allele, p-value) for traits involving drugs indicated for a disease.

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `graph`, `disease_id`, `sortBy`, `sortDir`, `pageNum`, `pageSize`
- **Keywords:** gwas catalog, snp, drug response, chebi, pharmacogenomics

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `graph` | any | yes | One of `ebi_monarch_xspecies` |
| `disease_id` | string | yes | Disease |
| `sortBy` | any | yes | One of `disease`, `trait`, `drug`, `snp`, `effect_allele`, `gene`, `study_accession`, `p_value`, `pubmed_id`, `reported_trait`, `chromosome`, `position` |
| `sortDir` | any | yes | One of `asc`, `desc` |
| `pageNum` | integer | yes | Page number (0-based) |
| `pageSize` | integer | yes | Number of results per page |

??? note "Upstream description"

    Search GWAS Catalog by drugs or drug classes indicated for a disease

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'grebi',
    'gwas_by_drugs_indicated_for_disease',
    {'graph': 'ebi_monarch_xspecies',
     'disease_id': '...',
     'sortBy': 'disease',
     'sortDir': 'asc',
     'pageNum': 0,
     'pageSize': 0},
)
```

## `gwas_by_gene_and_cell_type` {#gwas_by_gene_and_cell_type}

Returns GWAS Catalog associations mapped to a gene for traits linked to a given cell type.

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `graph`, `gene_id`, `cell_type_id`, `sortBy`, `sortDir`, `pageNum`, `pageSize`
- **Keywords:** gwas catalog, gene, snp, cell type, trait

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `graph` | any | yes | One of `ebi_monarch_xspecies` |
| `gene_id` | string | yes | Gene |
| `cell_type_id` | string | yes | Cell Type |
| `sortBy` | any | yes | One of `cell_type`, `trait`, `snp`, `effect_allele`, `gene`, `study_accession`, `p_value`, `pubmed_id`, `reported_trait`, `chromosome`, `position`, `risk_allele_frequency` |
| `sortDir` | any | yes | One of `asc`, `desc` |
| `pageNum` | integer | yes | Page number (0-based) |
| `pageSize` | integer | yes | Number of results per page |

??? note "Upstream description"

    Search GWAS Catalog by gene and cell type

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'grebi',
    'gwas_by_gene_and_cell_type',
    {'graph': 'ebi_monarch_xspecies',
     'gene_id': '...',
     'cell_type_id': '...',
     'sortBy': 'cell_type',
     'sortDir': 'asc',
     'pageNum': 0,
     'pageSize': 0},
)
```

## `gwas_by_gene_and_disease` {#gwas_by_gene_and_disease}

Returns GWAS Catalog associations mapped to a gene for traits within a given disease context.

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `graph`, `gene_id`, `disease_id`, `sortBy`, `sortDir`, `pageNum`, `pageSize`
- **Keywords:** gwas catalog, gene, snp, disease, trait

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `graph` | any | yes | One of `ebi_monarch_xspecies` |
| `gene_id` | string | yes | Gene |
| `disease_id` | string | yes | Disease |
| `sortBy` | any | yes | One of `disease`, `trait`, `snp`, `effect_allele`, `gene`, `study_accession`, `p_value`, `pubmed_id`, `reported_trait`, `chromosome`, `position`, `risk_allele_frequency` |
| `sortDir` | any | yes | One of `asc`, `desc` |
| `pageNum` | integer | yes | Page number (0-based) |
| `pageSize` | integer | yes | Number of results per page |

??? note "Upstream description"

    Search GWAS Catalog by gene and disease context

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'grebi',
    'gwas_by_gene_and_disease',
    {'graph': 'ebi_monarch_xspecies',
     'gene_id': '...',
     'disease_id': '...',
     'sortBy': 'disease',
     'sortDir': 'asc',
     'pageNum': 0,
     'pageSize': 0},
)
```

## `gwas_by_gene_and_location` {#gwas_by_gene_and_location}

Returns GWAS Catalog associations mapped to a gene for traits linked to an anatomical location.

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `graph`, `gene_id`, `location_id`, `sortBy`, `sortDir`, `pageNum`, `pageSize`
- **Keywords:** gwas catalog, gene, snp, anatomy, trait

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `graph` | any | yes | One of `ebi_monarch_xspecies` |
| `gene_id` | string | yes | Gene |
| `location_id` | string | yes | Anatomical Location |
| `sortBy` | any | yes | One of `location`, `trait`, `snp`, `effect_allele`, `gene`, `study_accession`, `p_value`, `pubmed_id`, `reported_trait`, `chromosome`, `position`, `risk_allele_frequency` |
| `sortDir` | any | yes | One of `asc`, `desc` |
| `pageNum` | integer | yes | Page number (0-based) |
| `pageSize` | integer | yes | Number of results per page |

??? note "Upstream description"

    Search GWAS Catalog by gene and anatomical location

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'grebi',
    'gwas_by_gene_and_location',
    {'graph': 'ebi_monarch_xspecies',
     'gene_id': '...',
     'location_id': '...',
     'sortBy': 'location',
     'sortDir': 'asc',
     'pageNum': 0,
     'pageSize': 0},
)
```

## `gwas_by_location` {#gwas_by_location}

Returns GWAS Catalog associations (SNP, gene, allele, p-value, study) for traits linked to an anatomical location or its parts.

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `graph`, `location_id`, `sortBy`, `sortDir`, `pageNum`, `pageSize`
- **Keywords:** gwas catalog, snp, uberon, anatomy, p-value

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `graph` | any | yes | One of `ebi_monarch_xspecies` |
| `location_id` | string | yes | Anatomical Location |
| `sortBy` | any | yes | One of `location`, `trait`, `snp`, `effect_allele`, `gene`, `study_accession`, `p_value`, `pubmed_id`, `reported_trait`, `chromosome`, `position`, `risk_allele_frequency` |
| `sortDir` | any | yes | One of `asc`, `desc` |
| `pageNum` | integer | yes | Page number (0-based) |
| `pageSize` | integer | yes | Number of results per page |

??? note "Upstream description"

    Search GWAS Catalog by anatomical location

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'grebi',
    'gwas_by_location',
    {'graph': 'ebi_monarch_xspecies',
     'location_id': '...',
     'sortBy': 'location',
     'sortDir': 'asc',
     'pageNum': 0,
     'pageSize': 0},
)
```

## `gwas_by_pathway` {#gwas_by_pathway}

Returns GWAS Catalog SNP associations for genes in a Reactome pathway, with alleles, traits and p-values.

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `graph`, `pathway_id`, `sortBy`, `sortDir`, `pageNum`, `pageSize`
- **Keywords:** gwas catalog, reactome, pathway genes, snp, p-value

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `graph` | any | yes | One of `ebi_monarch_xspecies` |
| `pathway_id` | string | yes | Reactome Pathway |
| `sortBy` | any | yes | One of `pathway`, `gene`, `snp`, `effect_allele`, `trait`, `study_accession`, `p_value`, `pubmed_id`, `reported_trait`, `chromosome`, `position`, `risk_allele_frequency` |
| `sortDir` | any | yes | One of `asc`, `desc` |
| `pageNum` | integer | yes | Page number (0-based) |
| `pageSize` | integer | yes | Number of results per page |

??? note "Upstream description"

    Search GWAS Catalog by Reactome pathway

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'grebi',
    'gwas_by_pathway',
    {'graph': 'ebi_monarch_xspecies',
     'pathway_id': '...',
     'sortBy': 'pathway',
     'sortDir': 'asc',
     'pageNum': 0,
     'pageSize': 0},
)
```

## `gwas_studies_by_cell_type` {#gwas_studies_by_cell_type}

Returns GWAS Catalog studies annotated with traits linked to a cell type, with sample sizes, association counts and summary statistics availability.

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `graph`, `cell_type_id`, `sortBy`, `sortDir`, `pageNum`, `pageSize`
- **Keywords:** gwas study, cell type, summary statistics, gcst, sample size

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `graph` | any | yes | One of `ebi_monarch_xspecies` |
| `cell_type_id` | string | yes | Cell Type |
| `sortBy` | any | yes | One of `cell_type`, `trait`, `study`, `study_accession`, `reported_trait`, `first_author`, `journal`, `publication_date`, `pubmed_id`, `initial_sample_size`, `replication_sample_size`, `association_count` |
| `sortDir` | any | yes | One of `asc`, `desc` |
| `pageNum` | integer | yes | Page number (0-based) |
| `pageSize` | integer | yes | Number of results per page |

??? note "Upstream description"

    Search GWAS Catalog studies by cell type

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'grebi',
    'gwas_studies_by_cell_type',
    {'graph': 'ebi_monarch_xspecies',
     'cell_type_id': '...',
     'sortBy': 'cell_type',
     'sortDir': 'asc',
     'pageNum': 0,
     'pageSize': 0},
)
```

## `gwas_studies_by_disease_location` {#gwas_studies_by_disease_location}

Returns GWAS Catalog studies for diseases located in an anatomical location, with study metadata.

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `graph`, `location_id`, `sortBy`, `sortDir`, `pageNum`, `pageSize`
- **Keywords:** gwas study, anatomy, disease location, summary statistics, gcst

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `graph` | any | yes | One of `ebi_monarch_xspecies` |
| `location_id` | string | yes | Anatomical Location |
| `sortBy` | any | yes | One of `location`, `disease`, `study`, `study_accession`, `reported_trait`, `first_author`, `journal`, `publication_date`, `pubmed_id`, `initial_sample_size`, `replication_sample_size`, `association_count` |
| `sortDir` | any | yes | One of `asc`, `desc` |
| `pageNum` | integer | yes | Page number (0-based) |
| `pageSize` | integer | yes | Number of results per page |

??? note "Upstream description"

    Search GWAS Catalog studies by disease location

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'grebi',
    'gwas_studies_by_disease_location',
    {'graph': 'ebi_monarch_xspecies',
     'location_id': '...',
     'sortBy': 'location',
     'sortDir': 'asc',
     'pageNum': 0,
     'pageSize': 0},
)
```

## `gwas_studies_by_drugs_indicated_for_disease` {#gwas_studies_by_drugs_indicated_for_disease}

Returns GWAS Catalog studies whose traits involve drugs indicated for a disease, with study metadata.

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `graph`, `disease_id`, `sortBy`, `sortDir`, `pageNum`, `pageSize`
- **Keywords:** gwas study, drug response, disease, gcst, summary statistics

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `graph` | any | yes | One of `ebi_monarch_xspecies` |
| `disease_id` | string | yes | Disease |
| `sortBy` | any | yes | One of `disease`, `trait`, `drug`, `study`, `study_accession`, `reported_trait`, `first_author`, `journal`, `publication_date`, `pubmed_id`, `initial_sample_size`, `replication_sample_size` |
| `sortDir` | any | yes | One of `asc`, `desc` |
| `pageNum` | integer | yes | Page number (0-based) |
| `pageSize` | integer | yes | Number of results per page |

??? note "Upstream description"

    Search GWAS Catalog studies by drugs or drug classes indicated for a disease

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'grebi',
    'gwas_studies_by_drugs_indicated_for_disease',
    {'graph': 'ebi_monarch_xspecies',
     'disease_id': '...',
     'sortBy': 'disease',
     'sortDir': 'asc',
     'pageNum': 0,
     'pageSize': 0},
)
```

## `gwas_trait_to_mouse_models_via_embeddings` {#gwas_trait_to_mouse_models_via_embeddings}

Finds mouse gene models with phenotypes semantically similar to a GWAS trait whose human orthologs have GWAS associations with that trait.

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `graph`, `trait_id`, `n_embedding_candidates`, `n_top_mp_hits`, `sortBy`, `sortDir`, `pageNum`, `pageSize`
- **Keywords:** gwas, mouse model, impc, orthologs, embedding similarity, mp phenotype

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `graph` | any | yes | One of `ebi_monarch_xspecies` |
| `trait_id` | string | yes | GWAS Trait (EFO term) |
| `n_embedding_candidates` | string | yes | Number of closest embedding matches to consider (across all ontologies) |
| `n_top_mp_hits` | string | yes | Number of top MP phenotype matches to expand into mouse models |
| `sortBy` | any | yes | One of `trait`, `mp_phenotype`, `embedding_distance`, `mouse_gene`, `human_gene`, `snp`, `effect_allele`, `study_accession`, `p_value`, `or_or_beta` |
| `sortDir` | any | yes | One of `asc`, `desc` |
| `pageNum` | integer | yes | Page number (0-based) |
| `pageSize` | integer | yes | Number of results per page |

??? note "Upstream description"

    GWAS trait to mouse models via embedding similarity

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'grebi',
    'gwas_trait_to_mouse_models_via_embeddings',
    {'graph': 'ebi_monarch_xspecies',
     'trait_id': '...',
     'n_embedding_candidates': '...',
     'n_top_mp_hits': '...',
     'sortBy': 'trait',
     'sortDir': 'asc',
     'pageNum': 0,
     'pageSize': 0},
)
```

## `mouse_gene_to_opentargets` {#mouse_gene_to_opentargets}

Returns Open Targets disease associations, with scores, for the human ortholog of a mouse gene.

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `graph`, `mouse_gene_id`, `min_score`, `sortBy`, `sortDir`, `pageNum`, `pageSize`
- **Keywords:** open targets, mouse gene, human ortholog, target-disease, impc

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `graph` | any | yes | One of `ebi_monarch_xspecies`, `impc_x_gwas` |
| `mouse_gene_id` | string | yes | Mouse Gene |
| `min_score` | string | yes | Minimum Score |
| `sortBy` | any | yes | One of `mouse_gene`, `human_gene`, `disease`, `score` |
| `sortDir` | any | yes | One of `asc`, `desc` |
| `pageNum` | integer | yes | Page number (0-based) |
| `pageSize` | integer | yes | Number of results per page |

??? note "Upstream description"

    Search OpenTargets disease associations by mouse gene

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'grebi',
    'mouse_gene_to_opentargets',
    {'graph': 'ebi_monarch_xspecies',
     'mouse_gene_id': '...',
     'min_score': '...',
     'sortBy': 'mouse_gene',
     'sortDir': 'asc',
     'pageNum': 0,
     'pageSize': 0},
)
```
