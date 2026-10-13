# G7 database functions: gget, MyGene, MyVariant and biothings_client

This page documents the G7 database wrappers. Each function is exported from `dnakit.databases`, for example `from dnakit.databases import gget_cosmic_local`.

- Local function: `gget_cosmic_local` reads a COSMIC file on disk and needs no network or account.
- Network functions: the rest query public databases. Default unit tests replace the client, so no real request is made. The example outputs below are real runs and need access to the corresponding service.
- Input files are declared with the `Example file` label; the other examples run as written.

## Existing and unwrapped functions

These are already wrapped in `gene_backends.py` and are marked as duplicates here: gget `search`, `info`, `seq`, `mutate`; MyGene `query`, `querymany`, `getgene`; MyVariant `getvariants`, `query`, `format_hgvs`, `get_hgvs_from_vcf`; the biothings_client HGVS helpers. `findgenes` is a deprecated alias.

The other unwrapped functions need an external program (`diamond`, `muscle`, and the DIAMOND path of `elm`), a large download (`ref`, `virus`), a GPU (`alphafold`), an extra package (`cellxgene` needs `cellxgene_census`), or are plotting, CLI, cache or configuration helpers. Protein-level `g2p` and `pdb` are not DNA-level. The full status for every inventory key is in the G7 ledger.

## 1) `gget_cosmic_local`

- **Purpose:** Search a local COSMIC mutation table by gene name, mutation ID, CDS or AA change. Nothing is downloaded and no login is needed.
- **Method:** Calls the local branch of `gget.cosmic`. When `cosmic_project` is omitted, the project is inferred as gget does: a file name containing `CancerMutationCensus_AllData` gives `cancer`, otherwise `other`. No match returns an empty list; spaces in column names become underscores.
- **API:** `gget_cosmic_local(searchterm[required], cosmic_tsv_path[required], cosmic_project[optional], limit[optional], verbose[optional])`
- **Input:** Search term; path to a local COSMIC TSV; `cosmic_project` one of cancer, cancer_example, census, resistance, cell_line, genome_screen, targeted_screen, other; `limit` from 1 to 10000.
- **Output:** Dictionary with `count`, `records` (in file order), `cosmic_project` and `provenance`. Commercial use of COSMIC is subject to its licence terms.

**Example file:** `CancerMutationCensus_AllData_mini.tsv`

```text
GENE_NAME	ACCESSION_NUMBER	LEGACY_MUTATION_ID	Mutation CDS	Mutation AA	GENOMIC_MUTATION_ID	Primary site
EGFR	NM_005228	COSM6223	c.2573T>G	p.L858R	COSV55365	lung
BRAF	NM_004333	COSM476	c.1799T>A	p.V600E	COSV55410	skin
EGFR	NM_005228	COSM6224	c.2236_2250del	p.E746_A750del	COSV55366	lung
```

**Example code:**

```python
from dnakit.databases import gget_cosmic_local

result = gget_cosmic_local("EGFR", "CancerMutationCensus_AllData_mini.tsv", verbose=False)
print(result["count"], [row["Mutation_AA"] for row in result["records"]])
```

**Example output:**

```text
2 ['p.L858R', 'p.E746_A750del']
```

## 2) `gget_blast`

- **Purpose:** Submit one nucleotide or protein sequence to NCBI BLAST through gget and return the hits.
- **Method:** Before submission the sequence is upper-cased and checked against the alphabet, and `program`, `database`, `limit` and `expect` are validated. A remote job can take minutes, so the example stops at validation; the network path is tested with a replacement client.
- **API:** `gget_blast(sequence[required], program[optional], database[optional], limit[optional], expect[optional], verbose[optional])`
- **Input:** Sequence string up to 20000 characters; `program` one of default, blastn, blastp, blastx, tblastn, tblastx; `expect` between 1e-30 and 1000.
- **Output:** Dictionary with `query_length`, `program`, `database`, `expect`, `count` and `hits`.

**Example code:**

```python
from dnakit.databases import gget_blast
from dnakit.exceptions import ConfigurationError

try:
    gget_blast("ACGTACGT", program="blastz", verbose=False)
except ConfigurationError as error:
    print(error)
```

**Example output:**

```text
[CONFIGURATION_ERROR] program must be one of default, blastn, blastp, blastx, tblastn, tblastx; got 'blastz'.
```

## 3) `gget_blat`

- **Purpose:** Locate a DNA or protein sequence in a chosen UCSC assembly with UCSC BLAT.
- **Method:** When `seqtype` is omitted, gget decides between nucleotide and protein from the letters. The wrapper only validates arguments and converts results to plain dictionaries. The example sequence comes from hg38 chr17:7668402-7668461, so the first hit should match exactly.
- **API:** `gget_blat(sequence[required], seqtype[optional], assembly[optional], verbose[optional])`
- **Input:** Sequence string; `seqtype` one of default, DNA, protein, translated RNA, translated DNA; `assembly` a UCSC assembly name such as human.
- **Output:** Dictionary with `query_length`, `seqtype`, `assembly`, `count` and `hits`. Each hit has `chromosome`, `start`, `end`, `strand` and `%_matched`.

**Example code:**

```python
from dnakit.databases import gget_blat

sequence = "CACCCCTCAGACACACAGGTGGCAGCAAAGTTTTATTGTAAAATAAGAGATCGATATAAA"
result = gget_blat(sequence, seqtype="DNA", verbose=False)
top = result["hits"][0]
print(result["query_length"], result["count"], top["chromosome"], top["start"], top["end"], top["%_matched"])
```

**Example output:**

```text
60 4 chr17 7668402 7668461 100.0
```

## 4) `gget_opentargets`

- **Purpose:** Query OpenTargets for diseases, drugs, expression, interactions or pharmacogenetic records of one Ensembl gene.
- **Method:** Accepts only Ensembl gene IDs (ENS...G) and a `resource` from the supported list; calls the OpenTargets GraphQL endpoint through gget.
- **API:** `gget_opentargets(ensembl_id[required], resource[optional], limit[optional], verbose[optional])`
- **Input:** Ensembl gene ID such as ENSG00000157764 (BRAF); `resource` one of depmap, diseases, drugs, expression, interactions, pharmacogenetics, tractability.
- **Output:** Dictionary with `ensembl_id`, `resource`, `count` and `records`. Record fields come from OpenTargets with dotted names.

**Example code:**

```python
from dnakit.databases import gget_opentargets

result = gget_opentargets("ENSG00000157764", resource="drugs", limit=2, verbose=False)
print(result["count"], [row["drug.name"] for row in result["records"]])
```

**Example output:**

```text
2 ['BELVARAFENIB', 'PLIXORAFENIB']
```

## 5) `gget_archs4`

- **Purpose:** Find genes most co-expressed with a gene of interest in ARCHS4, or read its tissue expression atlas.
- **Method:** `which` is correlation (co-expression) or tissue (tissue atlas); `species` is human or mouse. With `ensembl=True` the input must be an Ensembl gene ID.
- **API:** `gget_archs4(gene[required], ensembl[optional], which[optional], gene_count[optional], species[optional], verbose[optional])`
- **Input:** Gene symbol or Ensembl ID; `gene_count` from 1 to 1000.
- **Output:** Dictionary with `gene`, `which`, `species`, `count` and `records`.

**Example code:**

```python
from dnakit.databases import gget_archs4

result = gget_archs4("SOX2", which="correlation", gene_count=3, verbose=False)
print(result["count"], [row["gene_symbol"] for row in result["records"]])
```

**Example output:**

```text
3 ['SOX21', 'PTPRZ1', 'SALL3']
```

## 6) `gget_bgee`

- **Purpose:** Retrieve orthologs of one gene, or expression calls for several genes, from Bgee.
- **Method:** `record_type='orthologs'` accepts exactly one Ensembl gene ID; `record_type='expression'` accepts up to 100 IDs.
- **API:** `gget_bgee(gene_id[required], record_type[optional], verbose[optional])`
- **Input:** One Ensembl gene ID or a list of IDs; `record_type` is orthologs or expression.
- **Output:** Dictionary with `gene_ids`, `record_type`, `count` and `records`.

**Example code:**

```python
from dnakit.databases import gget_bgee

result = gget_bgee("ENSG00000139618", record_type="orthologs", verbose=False)
print(result["count"], result["records"][0]["gene_name"])
```

**Example output:**

```text
131 BRCA2
```

## 7) `gget_enrichr`

- **Purpose:** Run Enrichr gene-set enrichment on a gene list and return the statistics for each pathway or term.
- **Method:** Up to 500 genes, order kept; an optional background list of up to 50000 genes. The first record is the server's top-ranked term.
- **API:** `gget_enrichr(genes[required], database[required], species[optional], background_list[optional], verbose[optional])`
- **Input:** Gene symbols; `database` an Enrichr library name such as `KEGG_2021_Human`; `species` human or mouse.
- **Output:** Dictionary with `gene_count`, `database`, `species`, `count` and `records`.

**Example code:**

```python
from dnakit.databases import gget_enrichr

result = gget_enrichr(["TP53", "MDM2", "CDKN1A", "BAX"], database="KEGG_2021_Human", verbose=False)
print(result["gene_count"], result["count"], result["records"][0]["path_name"])
```

**Example output:**

```text
4 71 Melanoma
```

## 8) `gget_specificity`

- **Purpose:** Retrieve gene-level tissue or cell specificity statistics from 8cubeDB.
- **Method:** The gene list is validated and passed to gget's 8cubeDB endpoint; records keep the server order.
- **API:** `gget_specificity(gene_list[required], verbose[optional])`
- **Input:** List of gene symbols or IDs, up to 500.
- **Output:** Dictionary with `gene_count`, `count` and `records`.

**Example code:**

```python
from dnakit.databases import gget_specificity

result = gget_specificity(["SOX2", "PAX6"], verbose=False)
print(result["gene_count"], result["count"], result["records"][0]["gene_name"])
```

**Example output:**

```text
2 124 Pax6
```

## 9) `gget_gene_expression`

- **Purpose:** Retrieve normalized expression values from 8cubeDB for one analysis level and partition design.
- **Method:** `analysis_level` and `analysis_type` are passed to 8cubeDB unchanged (for example Across_tissues and Sex); only their characters are checked.
- **API:** `gget_gene_expression(gene_list[required], analysis_level[required], analysis_type[required], verbose[optional])`
- **Input:** Gene list up to 500; analysis level and partition design strings.
- **Output:** Dictionary with `gene_count`, `analysis_level`, `analysis_type`, `count` and `records`.

**Example code:**

```python
from dnakit.databases import gget_gene_expression

result = gget_gene_expression(
    ["SOX2", "PAX6"], analysis_level="Across_tissues", analysis_type="Sex", verbose=False
)
print(result["count"], [row["gene_name"] for row in result["records"]])
```

**Example output:**

```text
2 ['Pax6', 'Sox2']
```

## 10) `gget_psi_block`

- **Purpose:** Retrieve psi-block specificity scores from 8cubeDB to compare how specific each gene is across partitions.
- **Method:** Same validation as `gget_gene_expression`; calls the 8cubeDB psi-block endpoint.
- **API:** `gget_psi_block(gene_list[required], analysis_level[required], analysis_type[required], verbose[optional])`
- **Input:** Gene list up to 500; analysis level and partition design strings.
- **Output:** Dictionary with `gene_count`, `analysis_level`, `analysis_type`, `count` and `records`.

**Example code:**

```python
from dnakit.databases import gget_psi_block

result = gget_psi_block(
    ["SOX2", "PAX6"], analysis_level="Across_tissues", analysis_type="Sex", verbose=False
)
print(result["count"], [row["gene_name"] for row in result["records"]])
```

**Example output:**

```text
2 ['Pax6', 'Sox2']
```

## 11) `mygene_getgenes`

- **Purpose:** Batch-fetch MyGene.info gene records for Entrez or Ensembl IDs.
- **Method:** Calls `mygene.MyGeneInfo().getgenes`. Records follow the input order; IDs MyGene does not recognize are listed in `not_found` rather than dropped.
- **API:** `mygene_getgenes(gene_ids[required], fields[optional], species[optional], verbose[optional])`
- **Input:** List of Entrez or Ensembl IDs, up to 1000; `fields` a comma-separated field list; `species` human or an NCBI taxon ID.
- **Output:** Dictionary with `requested`, `found`, `not_found`, `records` and `provenance`.

**Example code:**

```python
from dnakit.databases import mygene_getgenes

result = mygene_getgenes(["1017", "999999999", "7157"], fields="symbol,name", verbose=False)
print(result["found"], [row["symbol"] for row in result["records"]], result["not_found"])
```

**Example output:**

```text
2 ['CDK2', 'TP53'] ['999999999']
```

## 12) `myvariant_querymany`

- **Purpose:** Batch-query variant annotations from MyVariant.info by HGVS ID or rsID.
- **Method:** Calls `myvariant.MyVariantInfo().querymany`. One input can return several hits; output follows the input order. `matched_terms` counts inputs with at least one hit and `not_found` lists inputs without one.
- **API:** `myvariant_querymany(query_terms[required], scopes[optional], fields[optional], assembly[optional], verbose[optional])`
- **Input:** List of HGVS IDs or rsIDs, up to 1000; `assembly` hg19 or hg38.
- **Output:** Dictionary with `requested`, `matched_terms`, `hit_count`, `not_found`, `records` and `provenance`.

**Example code:**

```python
from dnakit.databases import myvariant_querymany

result = myvariant_querymany(["rs58991260"], scopes="dbsnp.rsid", fields="dbsnp.rsid", assembly="hg38", verbose=False)
print(result["matched_terms"], result["records"][0]["_id"], result["not_found"])
```

**Example output:**

```text
1 chr1:g.218458480G>A []
```

## 13) `biothings_client_get_homologs`

- **Purpose:** Find orthologous genes in other species with `utils.gene.get_homologs` from biothings_client.
- **Method:** Reads homologene annotations first, then fetches the homolog records in batch. The joining logic is biothings_client's; the species filter applies to its result.
- **API:** `biothings_client_get_homologs(gene_list[required], fields[optional], species[optional], verbose[optional])`
- **Input:** Entrez or Ensembl gene IDs, up to 500; `species` is `all` or comma-separated NCBI taxon IDs.
- **Output:** Dictionary with `requested`, `results` and `provenance`. Each gene lists its homologs grouped by taxon ID.

**Example code:**

```python
from dnakit.databases import biothings_client_get_homologs

result = biothings_client_get_homologs(["7157"], fields="symbol", species="10090", verbose=False)
homologs = result["results"][0]["homologs"]
print(result["requested"], [(h["taxid"], g["symbol"]) for h in homologs for g in h["genes"]])
```

**Example output:**

```text
1 [(10090, 'Trp53')]
```
