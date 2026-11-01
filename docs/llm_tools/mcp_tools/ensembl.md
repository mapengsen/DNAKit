# Ensembl MCP (cyanheads) (`ensembl`)

DNAKit 的 MCP 工具索引收录了这个服务器的 **7** 个 DNA/基因相关工具（服务器共列出 7 个，2026-10-16 读取）。工具在用户自己安装服务器后运行，DNAKit 本身不包含它们。调用方式见 [调用索引中的工具](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found)。下面的工具说明是上游英文原文。

| 项目 | 内容 |
| --- | --- |
| 源码 | [cyanheads/ensembl-mcp-server @ `2ff6cd7`](https://github.com/cyanheads/ensembl-mcp-server/tree/2ff6cd793801fea4158e6c0814bd364f774a2743) |
| 许可证 | Apache-2.0 |
| 启动命令 | `node dist/index.js` |
| 传输方式 | stdio, http |
| 安装 | Separately prepare the pinned @cyanheads/ensembl-mcp-server build and dependencies; set workdir to the package root and executable_path to Node &gt;=24. |
| 副作用 | Sends selected gene, sequence, variant and homology queries to Ensembl REST or the configured mirror. |

| 工具 | 类别 | 简介 |
| --- | --- | --- |
| [`ensembl_get_homology`](#ensembl_get_homology) | `genome_annotation` | Finds orthologs and/or paralogs of a gene (by symbol+species or Ensembl ID), with homology type, percent identity and taxonomy level. |
| [`ensembl_get_sequence`](#ensembl_get_sequence) | `sequence` | Fetches a paged window of genomic, cDNA, CDS or protein sequence for an Ensembl ID or a chr:start-end region, with optional flanking and reverse strand. |
| [`ensembl_get_xrefs`](#ensembl_get_xrefs) | `genome_annotation` | Retrieves cross-database references for an Ensembl gene or transcript ID, optionally filtered to one database such as HGNC, EntrezGene, RefSeq or UniProt. |
| [`ensembl_list_species`](#ensembl_list_species) | `genome_annotation` | Lists Ensembl species with display name, common name, genome assembly, taxon ID and division, with name filtering and paging. |
| [`ensembl_lookup_gene`](#ensembl_lookup_gene) | `genome_annotation` | Resolves a gene symbol or Ensembl ID (single or batch of 20) to its stable ID, genomic location, biotype, description and transcripts, handling aliases. |
| [`ensembl_predict_variant`](#ensembl_predict_variant) | `variant` | Runs Ensembl VEP on an HGVS, region+allele or rsID variant and returns consequence terms, impact, affected transcripts and colocated known variants. |
| [`ensembl_query_region`](#ensembl_query_region) | `genome_annotation` | Finds genes, transcripts, exons, variants or regulatory elements overlapping a chromosomal region (up to 5 Mb), with biotype filtering. |

## `ensembl_get_homology` {#ensembl_get_homology}

Finds orthologs and/or paralogs of a gene (by symbol+species or Ensembl ID), with homology type, percent identity and taxonomy level.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**无
- **关键词：**ortholog, paralog, homolog, cross-species, conservation, compara
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `symbol` | string | 否 | Gene symbol in the source species (e.g. BRCA2, TP53). Species defaults to homo_sapiens; set species for other organisms. Cannot be combined with id. |
| `id` | string | 否 | Ensembl stable gene ID (e.g. ENSG00000139618). A version suffix (ENSG00000139618.17) is accepted and dropped before the query. Use ensembl_lookup_gene to get the stable ID from a symbol. Cannot be combined with symbol. |
| `species` | string | 否 | Source species (the species the query gene belongs to) in Ensembl internal format. Default is homo_sapiens. Use ensembl_list_species with nameContains (a genus or common name, e.g. mouse) to find valid values. 默认 `homo_sapiens` |
| `target_species` | string | 否 | Filter to homologs in a single target species (e.g. mus_musculus for mouse). Omit to return homologs across all available species. Use ensembl_list_species with nameContains (a genus or common name) to find valid values. |
| `type` | string | 否 | Type of homologs to return. orthologues: genes related by speciation (cross-species equivalents). paralogues: genes related by duplication (within or across species). all: both orthologs and paralogs. 可选值：`orthologues`, `paralogues`, `all` 默认 `orthologues` |
| `max_results` | integer | 否 | Maximum number of homologs to return. Broad orthology queries (e.g. BRCA2 across all species) can return 150+ homologs; the default keeps responses focused. Set to 0 to return every homolog uncapped. totalCount always reports the true number available before this cap. 默认 `25` |

??? note "上游说明（英文原文）"

    Find orthologs and/or paralogs of a gene across species. Returns each homolog's stable ID, species, homology type (ortholog_one2one, ortholog_one2many, paralog_many2many, etc.), perc_id (percent identity), perc_pos (percent positives), and taxonomy level. Essential for cross-species research — for example, "what is the mouse equivalent of human TP53?" or "how conserved is BRCA2 across mammals?". Provide either symbol + species or a stable gene ID. Target species can be filtered to a single species or left open to return all available homologs.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'ensembl',
    'ensembl_get_homology',
    {},
)
```

## `ensembl_get_sequence` {#ensembl_get_sequence}

Fetches a paged window of genomic, cDNA, CDS or protein sequence for an Ensembl ID or a chr:start-end region, with optional flanking and reverse strand.

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`id`
- **关键词：**dna sequence, cdna, cds, genomic region, flanking, reverse strand
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `id` | string | 是 | Ensembl stable ID (ENSG…, ENST…, ENSP…) or a genomic region for region mode — an id containing a colon is a region. A version suffix on a stable ID (ENST00000380152.8) is dropped before the request. Region accepts species:chr:start-end (e.g. homo_sapiens:13:32315086-32400268) or a bare chr:start-end (e.g. 13:32315086-32400268) when the species field is set. Coordinates are whole numbers; start..e… |
| `type` | string | 否 | Sequence type to retrieve. genomic: full genomic DNA including introns (default). cdna: spliced transcript sequence (requires ENST… ID). cds: coding sequence only, no UTRs (requires ENST… ID with coding transcript). protein: amino acid sequence (requires ENST… or ENSP… ID). Region ids are genomic-only — request cdna, cds, or protein from a transcript or protein stable ID. A protein ID (ENSP…) ser… 可选值：`genomic`, `cdna`, `cds`, `protein` 默认 `genomic` |
| `species` | string | 否 | Species in Ensembl internal format (e.g. homo_sapiens). Required for a bare chr:start-end region. For the species:chr:start-end form the species in id takes precedence and this field is ignored. Unused for stable ID lookups — Ensembl infers species from the ID prefix. |
| `expand_5prime` | integer | 否 | Number of base pairs to extend upstream (5' direction) of the requested feature. Default 0. Only applies to genomic sequences and region queries; on a :-1 region it extends the higher coordinate. 默认 `0` |
| `expand_3prime` | integer | 否 | Number of base pairs to extend downstream (3' direction) of the requested feature. Default 0. Only applies to genomic sequences and region queries; on a :-1 region it extends the lower coordinate. 默认 `0` |
| `offset` | integer | 否 | 0-based character offset where the returned window starts, counted in the resolved sequence (including any expand_5prime/expand_3prime flank). Default 0. Pass nextOffset from a truncated response to fetch the following window; an offset at or past the end returns an empty window. 默认 `0` |
| `max_length` | integer | 否 | Maximum number of characters in the returned window. Default 10000. Set to 0 to return everything from offset to the end, uncapped. 默认 `10000` |

??? note "上游说明（英文原文）"

    Fetch the DNA, cDNA, CDS, or protein sequence for a gene, transcript, protein, or genomic region. Returns a window of the sequence — the first 10,000 characters by default — with its stable ID, molecule type, and full length. When more follows the window, truncated is true and nextOffset is the offset to request next; walking nextOffset reconstructs the whole sequence, and max_length 0 returns everything from offset to the end. The type parameter selects which sequence is fetched: genomic (default, includes introns), cdna (spliced transcript), cds (coding sequence only), protein. For region mode, set id to a region — either species:chr:start-end (e.g. homo_sapiens:13:32315086-32400268) or a bare chr:start-end with species set (e.g. id 13:32315086-32400268, species homo_sapiens), spanning at most 10,000,000 bases; a :-1 suffix (13:32315086-32315100:-1) returns the reverse-strand sequence. Regions return genomic DNA only. Protein sequences require a transcript or protein stable ID (ENST…/ENSP…), not a gene ID — use ensembl_lookup_gene with expand_transcripts=true to get the canonical transcript ID first.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'ensembl',
    'ensembl_get_sequence',
    {'id': '...'},
)
```

## `ensembl_get_xrefs` {#ensembl_get_xrefs}

Retrieves cross-database references for an Ensembl gene or transcript ID, optionally filtered to one database such as HGNC, EntrezGene, RefSeq or UniProt.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`id`
- **关键词：**cross reference, id mapping, hgnc, entrez, refseq, omim
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `id` | string | 是 | Ensembl stable gene ID (ENSG…) or transcript ID (ENST…). A version suffix on an ENS… ID (ENSG00000139618.17) is dropped before the request. Use ensembl_lookup_gene to get the stable ID from a gene symbol. xrefs/id returns the full cross-reference set (56+ entries for well-annotated genes like BRCA2). |
| `dbname` | string | 否 | Filter to a specific external database by its Ensembl internal name. Examples: HGNC (HGNC gene ID), Uniprot_gn (UniProt gene name), EntrezGene (NCBI Gene ID), MIM_GENE (OMIM disease gene), RefSeq_mRNA (NCBI RefSeq transcript), Reactome (pathway IDs), GO (Gene Ontology terms). Omit to return all available xrefs. |

??? note "上游说明（英文原文）"

    Retrieve cross-database references for a gene or feature — HGNC, UniProt, EntrezGene, OMIM, RefSeq, Reactome, and others. Returns each xref with its database name, primary ID, display ID, and description. The dbname filter narrows to specific databases; omit to return all xrefs. IDs returned here chain to protein (pubchem via UniProt), literature (pubmed via PubMed IDs), disease (OMIM via MIM_GENE), and pathway (Reactome) resources. Requires an Ensembl stable ID — use ensembl_lookup_gene to get the ENSG… ID first. Common dbname values: HGNC, Uniprot_gn, EntrezGene, MIM_GENE, RefSeq_mRNA, RefSeq_peptide, Reactome, GO (Gene Ontology), ChEMBL.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'ensembl',
    'ensembl_get_xrefs',
    {'id': '...'},
)
```

## `ensembl_list_species` {#ensembl_list_species}

Lists Ensembl species with display name, common name, genome assembly, taxon ID and division, with name filtering and paging.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**无
- **关键词：**species, genome assembly, taxon id, division, discovery
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `division` | string | 否 | Filter to a specific Ensembl division. EnsemblVertebrates includes human, mouse, zebrafish, and other vertebrates. EnsemblPlants covers crop and model plant genomes. EnsemblFungi, EnsemblMetazoa, EnsemblProtists cover non-vertebrate model organisms. Omit to return the endpoint default division (vertebrates). 可选值：`EnsemblVertebrates`, `EnsemblPlants`, `EnsemblFungi`, `EnsemblMetazoa`, `EnsemblProtists` |
| `nameContains` | string | 否 | Case-insensitive substring filter applied locally after fetching. Matches against species name, display name, and common name. Example: "sapiens" matches homo_sapiens; "mouse" matches mus_musculus. |
| `max_results` | integer | 否 | Maximum number of species to return, counted from offset. Default 50. Set to 0 to return every matching species from offset to the end, uncapped. totalCount always reports the full number of matches before this cap. 默认 `50` |
| `offset` | integer | 否 | 0-based position in the name-sorted matches where the returned page starts. Default 0. Pass nextOffset from a previous response to fetch the following page; an offset at or past totalCount returns no species. 默认 `0` |

??? note "上游说明（英文原文）"

    List species supported by Ensembl with display name, common name, assembly, taxon ID, and division. Required discovery step — species names like homo_sapiens are opaque to non-biologists and are the input format every other Ensembl tool expects. Use nameContains to find a species by partial name match; filter by division to select one. With no division, lists the endpoint default division — the vertebrates (~356 species on the default GRCh38 endpoint). Species are sorted by internal name and returned max_results at a time (default 50); totalCount reports every match, and when more follow, nextOffset is the offset for the next page.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'ensembl',
    'ensembl_list_species',
    {},
)
```

## `ensembl_lookup_gene` {#ensembl_lookup_gene}

Resolves a gene symbol or Ensembl ID (single or batch of 20) to its stable ID, genomic location, biotype, description and transcripts, handling aliases.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**无
- **关键词：**gene symbol, ensembl id, genomic location, biotype, transcripts, alias
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `symbol` | string | 否 | Gene symbol to look up (e.g. BRCA2, TP53, EGFR). Species defaults to homo_sapiens; set species for other organisms. Case-insensitive in most species. A synonym or former symbol (HER2, TB1) resolves too; one shared by several genes returns ambiguous with the candidates instead of gene. |
| `id` | string | 否 | Ensembl stable gene ID (e.g. ENSG00000139618). A version suffix on an ENS… ID (ENSG00000139618.19) is dropped before the lookup, which returns the current version. Species is not required for ID lookup. |
| `species` | string | 否 | Species in Ensembl internal format: lowercase scientific name with underscores (e.g. homo_sapiens, mus_musculus, danio_rerio). Optional for symbol lookups — defaults to homo_sapiens; set it for other organisms. Use ensembl_list_species with nameContains (a genus or common name) to find valid values. |
| `ids` | array of string | 否 | Batch lookup: up to 20 Ensembl stable IDs (ENSG…, ENST…). A version suffix on an ENS… ID is dropped before the lookup. Returns a succeeded/failed split; each item carries its query — the ID as submitted, version included — and each failed item a reason and recovery. Provide exactly one of symbol, id, ids, or symbols. |
| `symbols` | array of string | 否 | Batch lookup: up to 20 gene symbols. Species defaults to homo_sapiens; set species for other organisms. Returns a succeeded/failed split; each item carries its submitted query. Synonyms and former symbols resolve with matchedAlias; one shared by several genes fails with reason ambiguous and its candidates. Provide exactly one of symbol, id, ids, or symbols. |
| `expand_transcripts` | boolean | 否 | When true, include the full transcript list in the response. Each transcript has its ID, biotype, canonical flag, and coordinates. Default is false to keep responses compact. 默认 `False` |

??? note "上游说明（英文原文）"

    Resolve a gene by symbol + species (or by stable ID) to its Ensembl ID, genomic location (chr:start-end:strand), biotype, description, and transcript list. Entry point for most workflows — the stable ID and coordinates returned here are inputs to other tools. Accepts both symbol lookup (BRCA2 + homo_sapiens) and direct ID lookup (ENSG00000139618). A symbol that is not a current gene name is matched against synonyms and former symbols (HER2 → ERBB2), reported as matchedAlias; when it maps to several genes the result is ambiguous, listing the candidates to retry by id. Supports batch lookup of up to 20 IDs or symbols in one call via the ids or symbols field. Provide exactly one of symbol, id, ids, or symbols. For symbol lookups species defaults to homo_sapiens (override for other organisms); for ID lookups species is not needed. Use ensembl_list_species with nameContains (a genus or common name) to find valid species names.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'ensembl',
    'ensembl_lookup_gene',
    {},
)
```

## `ensembl_predict_variant` {#ensembl_predict_variant}

Runs Ensembl VEP on an HGVS, region+allele or rsID variant and returns consequence terms, impact, affected transcripts and colocated known variants.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`variant`
- **关键词：**vep, variant consequence, hgvs, rsid, impact, clinical significance
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant` | string | 是 | Variant in one of three formats: (1) HGVS notation — transcript-relative: ENST00000380152.8:c.2T&gt;A; genomic: 13:g.32316462T&gt;A; (2) Region+allele: chr:start:end:strand/allele — e.g. 1:65568:65568:1/T (strand is 1 for forward or -1 for reverse); (3) dbSNP rsID — e.g. rs334. Ensembl normalizes chromosome names; canonical vertebrate output omits the "chr" prefix, though a chr-prefixed name is also ac… |
| `species` | string | 否 | Species in Ensembl internal format. Default is homo_sapiens. For non-human variants, set the appropriate species (e.g. mus_musculus for mouse). Use ensembl_list_species with nameContains (a genus or common name) to find valid values. 默认 `homo_sapiens` |
| `max_transcript_consequences` | integer | 否 | Maximum transcript consequences to return per VEP record. High-impact variants can affect 60+ transcripts; the default keeps the response focused on the top consequences. Set to 0 to return every transcript consequence uncapped. transcriptConsequencesTotal on each record always reports the true pre-cap count. 默认 `10` |
| `max_pubmed_ids_per_variant` | integer | 否 | Maximum PubMed IDs to return per colocated known variant. Well-studied variants (e.g. rs334) cite 100+ papers; the default trims each list. Set to 0 to return every PubMed ID uncapped. pubmedTotal on each colocated variant reports the true pre-cap count. Ignored when include_all_colocated_pubmed is true. 默认 `10` |
| `include_all_colocated_pubmed` | boolean | 否 | When true, return every PubMed ID for each colocated variant, overriding max_pubmed_ids_per_variant. Default false to keep responses compact. 默认 `False` |

??? note "上游说明（英文原文）"

    Predict the functional consequences of a sequence variant using the Ensembl Variant Effect Predictor (VEP). Accepts three input formats: HGVS notation (transcript-relative, e.g. ENST00000380152.8:c.2T&gt;A, or genomic, e.g. 13:g.32316462T&gt;A); region+allele (chr:start:end:strand/allele, e.g. 1:65568:65568:1/T); and a dbSNP rsID (e.g. rs334). Returns the most severe consequence term, affected transcripts and genes, impact level (HIGH/MODERATE/LOW/MODIFIER), and any colocated known variants with clinical significance. HGVS input: provide the full notation including transcript version for best results. Region+allele input: Ensembl normalizes chromosome names and canonical vertebrate output omits the chr prefix (a chr-prefixed name is also accepted). By default the response caps transcript consequences (max_transcript_consequences) and per-variant PubMed IDs (max_pubmed_ids_per_variant) to keep large VEP results compact — well-studied variants like rs334 otherwise carry 60+ consequences and 100+ citations. Truthful totals are always reported; set a cap to 0 (or include_all_colocated_pubmed=true) to retrieve the full set.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'ensembl',
    'ensembl_predict_variant',
    {'variant': '...'},
)
```

## `ensembl_query_region` {#ensembl_query_region}

Finds genes, transcripts, exons, variants or regulatory elements overlapping a chromosomal region (up to 5 Mb), with biotype filtering.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`region`
- **关键词：**genomic region, overlap, locus, regulatory feature, exon, coordinates
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `species` | string | 否 | Species in Ensembl internal format (e.g. homo_sapiens, mus_musculus). Required for a bare chr:start-end region; for the species:chr:start-end form the species in region takes precedence and this field is ignored. Use ensembl_list_species with nameContains (a genus or common name, e.g. mouse) to find valid values. |
| `region` | string | 是 | Genomic region as chr:start-end with species set (e.g. 13:32315086-32400268), or species:chr:start-end (e.g. homo_sapiens:13:32315086-32400268), with start at or below end. Ensembl serves at most 5,000,000 bases per region; split a larger area into smaller windows. Ensembl normalizes chromosome names and canonical vertebrate output omits the chr prefix (13, not chr13); a chr-prefixed name like ch… |
| `feature` | array of string | 否 | Feature types to retrieve — at least one. Default is gene only. Requesting variation in a large region can match tens of thousands of features. Include variation only for targeted small regions (single gene loci or smaller). 默认 `\['gene'\]` |
| `biotype` | string | 否 | Optional biotype filter (e.g. protein_coding, lncRNA, SNV). Applied server-side by Ensembl. Not all feature types support biotype filtering. |
| `max_results` | integer | 否 | Maximum number of features to return. A gene-length region can hold tens of thousands of variation features; the default keeps the response compact. Set to 0 to return every feature uncapped. totalCount always reports the true number found before this cap. 默认 `100` |

??? note "上游说明（英文原文）"

    Find genomic features overlapping a chromosomal region: genes, transcripts, variants, regulatory elements, or exons. Returns each feature with its stable ID, type, location, biotype, and name, plus the genome assembly the coordinates are on. Useful for "what's in this locus?" and for seeding follow-up lookups. Region format is chr:start-end with species set (e.g. 13:32315086-32400268 for the BRCA2 locus), or species:chr:start-end (e.g. homo_sapiens:13:32315086-32400268), spanning at most 5,000,000 bases. Ensembl normalizes chromosome names and canonical vertebrate output omits the chr prefix (13, not chr13); a chr-prefixed name like chr13 is also accepted. The feature parameter defaults to gene only — requesting variation in an 85 kb region matches 44,000+ entries. Explicitly include variation, regulatory, transcript, or exon only when needed. The response returns up to max_results features (default 100) while totalCount always reports the full count; set max_results to 0 for every feature, or query a smaller region to see a different slice. Exon rows carry the parent transcript ID, so the same exon appears once per transcript it belongs to.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'ensembl',
    'ensembl_query_region',
    {'region': '...'},
)
```
