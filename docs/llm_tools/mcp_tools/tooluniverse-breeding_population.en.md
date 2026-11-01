# ToolUniverse: breeding_population (`tooluniverse`)

DNAKit's MCP tool index keeps **845** DNA/gene tools of this server (2927 listed, read on 2026-10-16). They run once you install the server; DNAKit does not bundle them. See [calling a tool that was found](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found).

This page lists the 14 tools in category `breeding_population`; other categories: [`sequence`](tooluniverse-sequence.md) (44), [`genome_annotation`](tooluniverse-genome_annotation.md) (253), [`variant`](tooluniverse-variant.md) (155), [`association`](tooluniverse-association.md) (138), [`expression`](tooluniverse-expression.md) (56), [`regulation_epigenome`](tooluniverse-regulation_epigenome.md) (54), [`sequencing_data`](tooluniverse-sequencing_data.md) (91), [`prediction_model`](tooluniverse-prediction_model.md) (17), [`molecular_tools`](tooluniverse-molecular_tools.md) (23).

| Item | Value |
| --- | --- |
| Source | [mims-harvard/ToolUniverse @ `9959ebe`](https://github.com/mims-harvard/ToolUniverse/tree/9959ebee7638773314362f1746980a15f5b17e5b) |
| License | Apache-2.0 |
| Command | `tooluniverse` |
| Transports | stdio |
| Installation | Install the pinned ToolUniverse distribution and required provider dependencies separately. Use the installed tooluniverse console script; unlike the README uvx --refresh example, this bridge does not refresh or install anything. |
| Effects | Starts the compact ToolUniverse MCP catalog. Selected tools may contact services, load models, execute programs or write files according to their provider definitions. |

| Tool | Category | Summary |
| --- | --- | --- |
| [`GeneNetwork_get_dataset_info`](#GeneNetwork_get_dataset_info) | `breeding_population` | Gets GeneNetwork dataset metadata: full name, data type, scale, tissue and access status. |
| [`GeneNetwork_list_datasets`](#GeneNetwork_list_datasets) | `breeding_population` | Lists datasets (tissues, platforms) available for a GeneNetwork genetic cross population such as BXD. |
| [`GeneNetwork_list_groups`](#GeneNetwork_list_groups) | `breeding_population` | Lists GeneNetwork genetic cross populations (e.g. BXD, HXB recombinant inbred strains) for a species. |
| [`GeneNetwork_list_species`](#GeneNetwork_list_species) | `breeding_population` | Lists all species available in GeneNetwork (mouse, rat, Arabidopsis, barley, etc.). |
| [`IGSR_search_populations`](#IGSR_search_populations) | `breeding_population` | Searches 1000 Genomes populations and superpopulations (AFR, AMR, EAS, EUR, SAS). |
| [`MouseMine_search_alleles`](#MouseMine_search_alleles) | `breeding_population` | Searches MouseMine for mouse alleles and mutants with allele types and attributes. |
| [`OpenTargets_get_biological_mouse_models_by_ensemblID`](#OpenTargets_get_biological_mouse_models_by_ensemblID) | `breeding_population` | Retrieves mouse models with allelic composition and genetic background for a target Ensembl gene ID. |
| [`PopGen_fst`](#PopGen_fst) | `breeding_population` | Calculates Weir-Cockerham Fst between two populations from allele frequencies and sample sizes. |
| [`PopGen_haplotype_count`](#PopGen_haplotype_count) | `breeding_population` | Estimates distinct haplotypes after recombination over generations for a multi-SNP region. |
| [`PopGen_hwe_test`](#PopGen_hwe_test) | `breeding_population` | Tests Hardy-Weinberg equilibrium from observed genotype counts with chi-square. |
| [`PopGen_inbreeding`](#PopGen_inbreeding) | `breeding_population` | Computes inbreeding coefficient F for a pedigree type over generations. |
| [`RGD_get_strain`](#RGD_get_strain) | `breeding_population` | Returns an RGD rat strain record with type, genetics, origin and background. |
| [`RGD_get_strain_annotations`](#RGD_get_strain_annotations) | `breeding_population` | Returns curated disease, phenotype and variant-trait annotations for an RGD rat strain. |
| [`RGD_search_strains`](#RGD_search_strains) | `breeding_population` | Searches the RGD rat strain catalog by keyword and strain type. |

## `GeneNetwork_get_dataset_info` {#GeneNetwork_get_dataset_info}

Gets GeneNetwork dataset metadata: full name, data type, scale, tissue and access status.

- **Category:** `breeding_population` (Breeding genotyping, markers, population genetics and strain resources)
- **Required arguments:** `dataset_name`
- **Keywords:** genenetwork, dataset metadata, bxd, expression profiling
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `dataset_name` | string | yes | Dataset name/abbreviation. Examples: 'HC_M2_0606_P' (Hippocampus Consortium M430v2), 'BR_M2_1106_R' (Brain mRNA). |

??? note "Upstream description"

    Get metadata about a specific GeneNetwork dataset including its full name, data type, data scale, tissue, and public access status. Returns information about expression profiling datasets from genetic crosses (BXD, LXS, etc.) covering various tissues and platforms. Use GeneNetwork_list_datasets to find available dataset names first.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GeneNetwork_get_dataset_info',
     'arguments': {'dataset_name': '...'}},
)
```

## `GeneNetwork_list_datasets` {#GeneNetwork_list_datasets}

Lists datasets (tissues, platforms) available for a GeneNetwork genetic cross population such as BXD.

- **Category:** `breeding_population` (Breeding genotyping, markers, population genetics and strain resources)
- **Required arguments:** `group`
- **Keywords:** genenetwork, bxd, recombinant inbred, systems genetics, datasets
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `group` | string | yes | Genetic cross group name in lowercase (e.g., 'bxd', 'bxh', 'hxb', 'bnlx-srjr-f2'). Obtain from GeneNetwork_list_groups. |

??? note "Upstream description"

    List all datasets available for a specific genetic cross population in GeneNetwork. Datasets represent different tissues, time points, or measurement technologies for the population (e.g., 'BXD Brain mRNA', 'BXD Liver Metabolites', 'BXD Phenotype'). Returns dataset names, abbreviations, data types, and timestamps. Use GeneNetwork_list_groups to find group names (e.g., 'bxd' for BXD cross).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GeneNetwork_list_datasets',
     'arguments': {'group': '...'}},
)
```

## `GeneNetwork_list_groups` {#GeneNetwork_list_groups}

Lists GeneNetwork genetic cross populations (e.g. BXD, HXB recombinant inbred strains) for a species.

- **Category:** `breeding_population` (Breeding genotyping, markers, population genetics and strain resources)
- **Required arguments:** `species`
- **Keywords:** genenetwork, recombinant inbred, genetic cross, qtl mapping
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `species` | string | yes | Species short name (e.g., 'mouse', 'rat', 'human', 'arabidopsis'). Get from GeneNetwork_list_species. |

??? note "Upstream description"

    List genetic cross populations (groups) available in GeneNetwork for a given species. These groups represent specific genetic crosses (e.g., BXD = C57BL/6J x DBA/2J recombinant inbred strains, BXH = BXH RI strains, HXB = Hybrid strains) used for QTL mapping and systems genetics. Each group can have hundreds of datasets measuring gene expression, metabolites, and phenotypes across the population. Returns group names, genetic types (riset, intercross, etc.), and IDs.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GeneNetwork_list_groups',
     'arguments': {'species': '...'}},
)
```

## `GeneNetwork_list_species` {#GeneNetwork_list_species}

Lists all species available in GeneNetwork (mouse, rat, Arabidopsis, barley, etc.).

- **Category:** `breeding_population` (Breeding genotyping, markers, population genetics and strain resources)
- **Required arguments:** none
- **Keywords:** genenetwork, species, systems genetics, model organisms
- **Read-only:** yes (declared by the server)

??? note "Upstream description"

    List all species available in GeneNetwork, the open systems genetics database for behavioral, neural, and physiological phenotypes. GeneNetwork contains genetic cross data for 13+ model organisms including mouse (BXD, BXH, AXB, etc.), rat, Arabidopsis, human, zebrafish, and barley. Returns species names and taxonomy IDs for use with other GeneNetwork tools.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GeneNetwork_list_species', 'arguments': {}},
)
```

## `IGSR_search_populations` {#IGSR_search_populations}

Searches 1000 Genomes populations and superpopulations (AFR, AMR, EAS, EUR, SAS).

- **Category:** `breeding_population` (Breeding genotyping, markers, population genetics and strain resources)
- **Required arguments:** none
- **Keywords:** 1000 genomes, populations, superpopulation, ancestry
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | Operation type (fixed: search_populations). One of `search_populations` Default `search_populations` |
| `superpopulation` | string \| null | no | Filter by superpopulation code: AFR (African), AMR (Admixed American), EAS (East Asian), EUR (European), SAS (South Asian). Leave empty for all. |
| `query` | string \| null | no | Free-text search query to find populations by name or description (e.g., 'Chinese', 'Yoruba', 'British'). |
| `limit` | integer | no | Maximum number of results to return. Default `25` |

??? note "Upstream description"

    Search 1000 Genomes Project populations from the International Genome Sample Resource (IGSR). The 1000 Genomes Project characterized genetic variation across 26 populations grouped into 5 superpopulations: AFR (African), AMR (Admixed American), EAS (East Asian), EUR (European), SAS (South Asian). Returns population codes, names, geographic coordinates, sample counts, and superpopulation assignments. Use this to find populations for allele frequency comparisons, population stratification analysis, or ancestry-specific variant studies.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'IGSR_search_populations', 'arguments': {}},
)
```

## `MouseMine_search_alleles` {#MouseMine_search_alleles}

Searches MouseMine for mouse alleles and mutants with allele types and attributes.

- **Category:** `breeding_population` (Breeding genotyping, markers, population genetics and strain resources)
- **Required arguments:** `q`
- **Keywords:** mouse alleles, knockout, targeted allele, mgi
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `q` | string | yes | Allele search query: gene symbol (e.g., 'Brca1'), allele symbol (e.g., 'Brca1&lt;tm1Whl&gt;'), or keyword. |
| `size` | integer | no | Number of results to return (default: 10). Default `10` |

??? note "Upstream description"

    Search MouseMine specifically for mouse alleles and mutants. Returns allele symbols, names, MGI identifiers, allele types (Targeted, Endonuclease-mediated, Spontaneous, etc.), and attribute strings (Null/knockout, etc.) from Mouse Genome Informatics (MGI). Useful for finding knockout/knockin alleles, spontaneous mutations, and transgenic alleles for specific genes. For broader searches, use MouseMine_search instead.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'MouseMine_search_alleles', 'arguments': {'q': '...'}},
)
```

## `OpenTargets_get_biological_mouse_models_by_ensemblID` {#OpenTargets_get_biological_mouse_models_by_ensemblID}

Retrieves mouse models with allelic composition and genetic background for a target Ensembl gene ID.

- **Category:** `breeding_population` (Breeding genotyping, markers, population genetics and strain resources)
- **Required arguments:** `ensemblId`
- **Keywords:** mouse model, allele, knockout, mgi, genetic background
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `ensemblId` | string | yes | The Ensembl ID of the target. |

??? note "Upstream description"

    Retrieve biological mouse models, including allelic compositions and genetic backgrounds, for a specific target.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'OpenTargets_get_biological_mouse_models_by_ensemblID',
     'arguments': {'ensemblId': '...'}},
)
```

## `PopGen_fst` {#PopGen_fst}

Calculates Weir-Cockerham Fst between two populations from allele frequencies and sample sizes.

- **Category:** `breeding_population` (Breeding genotyping, markers, population genetics and strain resources)
- **Required arguments:** `p1`, `p2`, `n1`, `n2`
- **Keywords:** fst, weir-cockerham, allele frequency, differentiation
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | Operation type One of `fst` Default `fst` |
| `p1` | number | yes | Allele frequency in population 1 \[0, 1\] |
| `p2` | number | yes | Allele frequency in population 2 \[0, 1\] |
| `n1` | integer | yes | Sample size (number of individuals) for population 1 |
| `n2` | integer | yes | Sample size (number of individuals) for population 2 |

??? note "Upstream description"

    Calculate the Weir-Cockerham (1984) theta Fst estimator between two populations given allele frequencies and sample sizes, including the n_c sample-size correction (theta = (MSP - MSG) / (MSP + (n_c - 1) * MSG)) so the result reflects allele-frequency divergence rather than sample size. Fst quantifies genetic differentiation: &lt;0.05 little, 0.05-0.15 moderate, 0.15-0.25 great, &gt;0.25 very great differentiation. Useful for measuring population structure and divergence at individual loci.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'PopGen_fst',
     'arguments': {'p1': 0.0, 'p2': 0.0, 'n1': 0, 'n2': 0}},
)
```

## `PopGen_haplotype_count` {#PopGen_haplotype_count}

Estimates distinct haplotypes after recombination over generations for a multi-SNP region.

- **Category:** `breeding_population` (Breeding genotyping, markers, population genetics and strain resources)
- **Required arguments:** `n_snps`
- **Keywords:** haplotype, recombination, snp
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | Operation type One of `haplotype_count` Default `haplotype_count` |
| `n_snps` | integer | yes | Number of biallelic SNP loci in the region |
| `generations` | integer | no | Number of generations of random mating (default 1) |
| `recomb_rate` | number | no | Per-adjacent-interval recombination rate \[0, 1\] per generation (default 0.5) |

??? note "Upstream description"

    Estimate distinct haplotypes after recombination across generations for a region with multiple biallelic SNPs. Reports maximum theoretical haplotypes (2^n), expected crossovers per gamete, probability of at least one recombination event, and estimated haplotype diversity based on LD decay. Useful for understanding recombination impact on genetic variation.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'PopGen_haplotype_count', 'arguments': {'n_snps': 0}},
)
```

## `PopGen_hwe_test` {#PopGen_hwe_test}

Tests Hardy-Weinberg equilibrium from observed genotype counts with chi-square.

- **Category:** `breeding_population` (Breeding genotyping, markers, population genetics and strain resources)
- **Required arguments:** `obs_AA`, `obs_Aa`, `obs_aa`
- **Keywords:** hwe, genotype counts, chi-square
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | Operation type One of `hwe_test` Default `hwe_test` |
| `obs_AA` | integer | yes | Observed count of homozygous dominant (AA) genotypes |
| `obs_Aa` | integer | yes | Observed count of heterozygous (Aa) genotypes |
| `obs_aa` | integer | yes | Observed count of homozygous recessive (aa) genotypes |

??? note "Upstream description"

    Test Hardy-Weinberg equilibrium for a biallelic locus given observed genotype counts (AA, Aa, aa). Computes allele frequencies, expected genotype counts under HWE, chi-square statistic with 1 df, and p-value. Reports whether the locus is in HWE (p &gt; 0.05) and, if not, whether heterozygote excess or deficit is observed. Useful for detecting selection, inbreeding, assortative mating, or genotyping error at individual loci.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'PopGen_hwe_test',
     'arguments': {'obs_AA': 0, 'obs_Aa': 0, 'obs_aa': 0}},
)
```

## `PopGen_inbreeding` {#PopGen_inbreeding}

Computes inbreeding coefficient F for a pedigree type over generations.

- **Category:** `breeding_population` (Breeding genotyping, markers, population genetics and strain resources)
- **Required arguments:** `pedigree`
- **Keywords:** inbreeding, pedigree, coefficient f
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | Operation type One of `inbreeding` Default `inbreeding` |
| `pedigree` | string | yes | Pedigree type (e.g., 'half-sib', 'first-cousin', 'full-sib') |
| `generations` | integer | no | Number of successive inbreeding generations (default 1) |

??? note "Upstream description"

    Compute inbreeding coefficient F for a known pedigree type over one or more generations. Supports: self, full-sib, half-sib, first-cousin, double-first-cousin, uncle-niece, aunt-nephew, half-first-cousin, second-cousin. Uses the recurrence F_t = base_F + (1-base_F)*F_\{t-1\} for iterated inbreeding. Reports F per generation, cumulative F, and retained heterozygosity.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'PopGen_inbreeding', 'arguments': {'pedigree': '...'}},
)
```

## `RGD_get_strain` {#RGD_get_strain}

Returns an RGD rat strain record with type, genetics, origin and background.

- **Category:** `breeding_population` (Breeding genotyping, markers, population genetics and strain resources)
- **Required arguments:** `symbol`
- **Keywords:** rat strain, inbred, congenic, model organism
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `symbol` | string | yes | Exact rat strain symbol (case-insensitive), e.g. 'SHR' (spontaneously hypertensive rat), 'BN' (Brown Norway), 'SS' (Dahl salt-sensitive). |

??? note "Upstream description"

    Get the full Rat Genome Database (RGD) record for a single rat strain by its exact strain symbol (e.g. 'SHR', 'BN', 'SS', 'WKY', 'GK'). Returns the strain RGD ID, name, full strain name, substrain, strain type, genetics formula, origin, background strain, and description. Use RGD_search_strains first if you only know a keyword.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'RGD_get_strain', 'arguments': {'symbol': '...'}},
)
```

## `RGD_get_strain_annotations` {#RGD_get_strain_annotations}

Returns curated disease, phenotype and variant-trait annotations for an RGD rat strain.

- **Category:** `breeding_population` (Breeding genotyping, markers, population genetics and strain resources)
- **Required arguments:** `symbol`
- **Keywords:** rat strain, disease model, phenotype
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `symbol` | string | yes | Exact rat strain symbol (case-insensitive), e.g. 'SHR', 'BN', 'SS'. |
| `category` | string \| null | no | Optional annotation category filter: 'disease', 'phenotype', 'strain', 'pathway', etc. Omit to return all categories. |

??? note "Upstream description"

    Get curated disease, phenotype, and variant-trait ontology annotations for a rat strain from the Rat Genome Database (RGD), by strain symbol (e.g. 'SHR', 'BN', 'SS'). This exposes RGD's rat disease-MODEL annotations - for example SHR is annotated to 'Left Ventricular Hypertrophy' (disease, MODEL: spontaneous) and 'arterial blood pressure trait' (phenotype). Each annotation includes the ontology term, accession (DOID/VT/CMO/MP/RS), category, qualifier, and evidence code. Optional 'category' filter: disease, phenotype, strain, pathway, biological_process, molecular_function, cellular_component.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'RGD_get_strain_annotations',
     'arguments': {'symbol': '...'}},
)
```

## `RGD_search_strains` {#RGD_search_strains}

Searches the RGD rat strain catalog by keyword and strain type.

- **Category:** `breeding_population` (Breeding genotyping, markers, population genetics and strain resources)
- **Required arguments:** none
- **Keywords:** rat strain, transgenic, inbred
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string \| null | no | Keyword to match against strain symbol, name, genetics, and description (case-insensitive substring). Example: 'hypertensive' or 'Dahl'. |
| `strain_type` | string \| null | no | Filter by strain type, e.g. 'inbred', 'congenic', 'transgenic', 'recombinant_inbred', 'consomic'. |
| `limit` | integer \| null | no | Maximum number of strains to return (default 25, max 200). |

??? note "Upstream description"

    Search the Rat Genome Database (RGD) rat strain catalog by keyword and/or strain type. Rat inbred/congenic/transgenic strains are RGD's flagship rat-specific resource (not in Alliance of Genome Resources or MGI). Provide 'query' (a keyword matched against symbol, name, full strain name, genetics, and description) and/or 'strain_type' (one of: inbred, congenic, mutant, transgenic, recombinant_inbred, outbred, consomic, coisogenic, segregating_inbred, hybrid, wild, conplastic, advanced_intercross_line, transchromosomal). At least one of the two must be supplied. Returns matching strains with their RGD ID, symbol, name, type, genetics, and origin.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'RGD_search_strains', 'arguments': {}},
)
```
