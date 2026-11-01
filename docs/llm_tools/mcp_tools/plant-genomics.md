# Plant Genomics MCP (`plant-genomics`)

DNAKit 的 MCP 工具索引收录了这个服务器的 **40** 个 DNA/基因相关工具（服务器共列出 50 个，2026-10-16 读取）。工具在用户自己安装服务器后运行，DNAKit 本身不包含它们。调用方式见 [调用索引中的工具](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found)。下面的工具说明是上游英文原文。

| 项目 | 内容 |
| --- | --- |
| 源码 | [musharna/plant-genomics-mcp @ `c689ffd`](https://github.com/musharna/plant-genomics-mcp/tree/c689ffd64f4110aff85f32cec2ed5681f6a4853f) |
| 许可证 | MIT |
| 启动命令 | `plant-genomics-mcp` |
| 传输方式 | stdio |
| 安装 | Install the pinned plant-genomics-mcp Python distribution and dependencies separately. The plant-genomics-mcp entry point selects stdio; no HTTP token is needed for this local transport. |
| 副作用 | Queries plant genomics services and may submit selected sequences to NCBI BLAST. Configure the named NCBI contact-email variable before BLAST use. |

| 工具 | 类别 | 简介 |
| --- | --- | --- |
| [`analyze_locus_synth`](#analyze_locus_synth) | `genome_annotation` | Resolves a plant locus via Ensembl Plants and fans out to xrefs, UniProt, Europe PMC and QuickGO, reconciling cross-source names. |
| [`arabidopsis_natural_variation`](#arabidopsis_natural_variation) | `variant` | Fetches 1001 Genomes natural-variation SNPs and their predicted effects across 1135 Arabidopsis accessions for a locus. |
| [`aragwas_associations`](#aragwas_associations) | `association` | Fetches AraGWAS significant SNP-phenotype associations overlapping an Arabidopsis locus, with score, MAF and study. |
| [`atted_coexpression`](#atted_coexpression) | `expression` | Fetches top co-expressed gene neighbours for a plant locus from ATTED-II with z-scores. |
| [`bar_efp_expression`](#bar_efp_expression) | `expression` | Fetches BAR world-eFP expression for an Arabidopsis locus across about 36 ecotypes with replicate values. |
| [`bar_gene_summary`](#bar_gene_summary) | `genome_annotation` | Fetches the TAIR/Araport11 curator summary, NCBI Gene ID and cross-database aliases for an Arabidopsis locus from BAR. |
| [`batch_atted_coexpression`](#batch_atted_coexpression) | `expression` | Batch version of ATTED-II co-expression lookup for up to 50 plant loci. |
| [`batch_bar_gene_summary`](#batch_bar_gene_summary) | `genome_annotation` | Batch fetch of BAR TAIR/Araport11 gene summaries and aliases for up to 50 Arabidopsis loci. |
| [`batch_ensembl_plants_lookup_locus`](#batch_ensembl_plants_lookup_locus) | `genome_annotation` | Batch lookup of Ensembl Plants gene metadata for up to 50 loci via POST /lookup/id. |
| [`batch_get_gene_xrefs`](#batch_get_gene_xrefs) | `genome_annotation` | Batch fetch of Ensembl Plants cross-database references for up to 50 plant loci. |
| [`batch_gramene_homologs`](#batch_gramene_homologs) | `genome_annotation` | Batch fetch of Gramene compara orthologs/paralogs for up to 50 plant loci. |
| [`batch_kegg_pathways`](#batch_kegg_pathways) | `genome_annotation` | Batch fetch of KEGG pathway memberships for up to 50 Arabidopsis loci. |
| [`batch_locus_go_annotations`](#batch_locus_go_annotations) | `genome_annotation` | Batch fetch of QuickGO Gene Ontology annotations for up to 50 plant loci via UniProt resolution. |
| [`batch_phytozome_lookup_locus`](#batch_phytozome_lookup_locus) | `genome_annotation` | Batch fetch of Phytozome BioMart gene records (chromosome, start/end, strand, description) for up to 50 loci. |
| [`batch_resolve_locus_to_uniprot`](#batch_resolve_locus_to_uniprot) | `genome_annotation` | Batch resolution of up to 50 plant loci to canonical UniProtKB records. |
| [`biological_context_synth`](#biological_context_synth) | `genome_annotation` | Resolves a locus and gathers Gramene homologs, KEGG pathways, STRING partners and ATTED-II coexpression with a consensus partner ranking. |
| [`blast_sequence`](#blast_sequence) | `sequence` | Runs an NCBI BLAST URLAPI search (blastn/blastp/blastx/tblastn/tblastx) and returns parsed top hits. |
| [`consensus_homologs`](#consensus_homologs) | `genome_annotation` | Builds a cross-source homology consensus for a locus by combining Gramene homologs and NCBI BLAST hits. |
| [`ensembl_plants_lookup_locus`](#ensembl_plants_lookup_locus) | `genome_annotation` | Fetches Ensembl Plants metadata for a plant locus identifier (coordinates, biotype, description). |
| [`ensembl_region_query`](#ensembl_region_query) | `genome_annotation` | Lists genes, transcripts, CDS or exons overlapping a genomic interval in an Ensembl Plants genome. |
| [`find_homologs_synth`](#find_homologs_synth) | `sequence` | Runs BLAST for a protein or nucleotide sequence and annotates ranked hits with UniProt records. |
| [`gene_report`](#gene_report) | `genome_annotation` | Generates a Markdown gene dossier for a plant locus from Ensembl Plants, UniProt, KEGG, STRING, Europe PMC and QuickGO. |
| [`get_gene_xrefs`](#get_gene_xrefs) | `genome_annotation` | Fetches Ensembl Plants cross-database references (UniProt, NCBI Gene, TAIR, etc.) for a plant locus. |
| [`get_sequence`](#get_sequence) | `sequence` | Fetches a locus's genomic, CDS, cDNA or protein sequence from Ensembl Plants. |
| [`go_enrichment`](#go_enrichment) | `genome_annotation` | Runs g:Profiler GO and KEGG over-representation analysis on a plant gene list with optional background. |
| [`gramene_homologs`](#gramene_homologs) | `genome_annotation` | Fetches orthologs and paralogs for a plant locus from Gramene compara with gene-tree IDs. |
| [`jaspar_motif`](#jaspar_motif) | `regulation_epigenome` | Fetches one JASPAR TF binding profile by matrix ID, including the position-frequency matrix, TF class and consensus. |
| [`kegg_pathways`](#kegg_pathways) | `genome_annotation` | Fetches KEGG pathway memberships for an Arabidopsis locus. |
| [`locus_gene_rifs`](#locus_gene_rifs) | `genome_annotation` | Fetches curated GeneRIF functional statements with PubMed IDs for an Arabidopsis locus from ThaleMine. |
| [`locus_go_annotations`](#locus_go_annotations) | `genome_annotation` | Fetches QuickGO Gene Ontology annotations for a plant locus, grouped by aspect. |
| [`locus_plant_ontology`](#locus_plant_ontology) | `genome_annotation` | Fetches Plant Ontology, Trait Ontology and PECO annotations for a plant locus from Planteome. |
| [`locus_variants`](#locus_variants) | `variant` | Lists EVA/dbSNP germline variants (SNPs, indels) overlapping a locus's genomic span via Ensembl. |
| [`orthodb_orthologs`](#orthodb_orthologs) | `genome_annotation` | Resolves a locus to its OrthoDB ortholog group at Viridiplantae level and lists member genes by organism. |
| [`panther_family`](#panther_family) | `genome_annotation` | Fetches PANTHER family/subfamily classification, curated GO terms, protein class and pathways for a locus. |
| [`phytozome_lookup_locus`](#phytozome_lookup_locus) | `genome_annotation` | Fetches a Phytozome BioMart gene record (chromosome, start/end, strand, description) for a locus. |
| [`plantcyc_locus_info`](#plantcyc_locus_info) | `genome_annotation` | Fetches PlantCyc metabolic annotation for a locus: enzymes, catalyzed reactions and pathways. |
| [`resolve_locus_to_uniprot`](#resolve_locus_to_uniprot) | `genome_annotation` | Resolves a plant locus to its canonical UniProtKB record (accession, gene names, sequence length). |
| [`tair_locus_info`](#tair_locus_info) | `genome_annotation` | Fetches the TAIR curator summary, Araport11 description, NCBI Gene ID and aliases for an Arabidopsis locus via BAR. |
| [`tf_binding_motifs`](#tf_binding_motifs) | `regulation_epigenome` | Fetches JASPAR transcription-factor DNA binding motifs for a locus, confirmed by UniProt accession, with consensus and assay type. |
| [`vep_annotate`](#vep_annotate) | `variant` | Predicts a variant's molecular consequences per overlapping transcript with Ensembl VEP, including SIFT/PolyPhen. |

## `analyze_locus_synth` {#analyze_locus_synth}

Resolves a plant locus via Ensembl Plants and fans out to xrefs, UniProt, Europe PMC and QuickGO, reconciling cross-source names.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`locus`
- **关键词：**plant gene, ensembl plants, xrefs, go, gene dossier
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `locus` | string | 是 | Locus name, e.g. AT1G01010 |
| `organism` | string \| integer | 否 | Plant organism — accepts canonical slug (arabidopsis_thaliana), scientific or common name, or NCBI taxid 默认 `arabidopsis_thaliana` |

??? note "上游说明（英文原文）"

    Synthesis: one-call equivalent of the analyze_locus prompt. Resolves a locus through Ensembl Plants, then fans out to xrefs, UniProt, Europe PMC, and QuickGO in parallel. Returns a SynthesisEnvelope with per-step status and a reconciled summary flagging cross-source name/accession disagreements.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'plant-genomics',
    'analyze_locus_synth',
    {'locus': '...'},
)
```

## `arabidopsis_natural_variation` {#arabidopsis_natural_variation}

Fetches 1001 Genomes natural-variation SNPs and their predicted effects across 1135 Arabidopsis accessions for a locus.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`locus`
- **关键词：**1001 genomes, snp, natural variation, accessions, arabidopsis
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `locus` | string | 是 | Arabidopsis AGI locus, e.g. AT1G01060 (a bare AGI is transcript-scoped to .1) |
| `organism` | string \| integer | 否 | Arabidopsis only (the 1001 Genomes panel is A. thaliana) 默认 `arabidopsis_thaliana` |

??? note "上游说明（英文原文）"

    Fetch 1001 Genomes natural-variation SNP effects for an Arabidopsis locus (tools.1001genomes.org; free, no key) — the variation observed across 1135 resequenced natural accessions. Returns per-SNP effect rows (chromosome, position, accession id, effect, impact, amino-acid change, transcript) plus the gene's genomic span. variant_count is the true row total even when capped. ARABIDOPSIS-ONLY — any other organism raises OrganismNotSupported. Defaults to arabidopsis_thaliana.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'plant-genomics',
    'arabidopsis_natural_variation',
    {'locus': '...'},
)
```

## `aragwas_associations` {#aragwas_associations}

Fetches AraGWAS significant SNP-phenotype associations overlapping an Arabidopsis locus, with score, MAF and study.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`locus`
- **关键词：**aragwas, gwas, snp, phenotype, maf, arabidopsis
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `locus` | string | 是 | Arabidopsis AGI locus, e.g. AT1G01060 |
| `organism` | string \| integer | 否 | Arabidopsis only (the 1001 Genomes panel is A. thaliana) 默认 `arabidopsis_thaliana` |

??? note "上游说明（英文原文）"

    Fetch AraGWAS genome-wide association study hits for an Arabidopsis locus (aragwas.1001genomes.org; free, no key). Returns each significant SNP association overlapping the gene with effect size (score), minor-allele frequency, the SNP's predicted molecular effect (impact, amino-acid change), and the phenotype/study it came from. association_count is the true total even when page-capped. ARABIDOPSIS-ONLY — any other organism raises OrganismNotSupported. Defaults to arabidopsis_thaliana.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'plant-genomics',
    'aragwas_associations',
    {'locus': '...'},
)
```

## `atted_coexpression` {#atted_coexpression}

Fetches top co-expressed gene neighbours for a plant locus from ATTED-II with z-scores.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**`locus`
- **关键词：**atted-ii, coexpression, plant, z-score, gene network
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `locus` | string | 是 | Plant locus, e.g. AT1G01010 (Arabidopsis) or Os01g0100100 (rice) |
| `organism` | string \| integer | 否 | Plant organism — accepts canonical slug (arabidopsis_thaliana), scientific or common name, or NCBI taxid 默认 `arabidopsis_thaliana` |
| `top_n` | integer | 否 | 默认 `25` |

??? note "上游说明（英文原文）"

    Fetch co-expressed gene neighbors from ATTED-II (atted.jp, API v5) for a plant locus. Returns top_n neighbors with target locus + NCBI Entrez gene ID + z-score (higher = stronger coexpression). The ATTED-II release (e.g. Ath-u.c4-0 for Arabidopsis, Osa-u.c1-0 for rice) is resolved per-organism; wheat, sorghum, barley, poplar, and brachypodium have no published release and raise OrganismNotSupported. Pairs with string_interactions to surface high-confidence functional partners (interactors that are also coexpressed).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'plant-genomics',
    'atted_coexpression',
    {'locus': '...'},
)
```

## `bar_efp_expression` {#bar_efp_expression}

Fetches BAR world-eFP expression for an Arabidopsis locus across about 36 ecotypes with replicate values.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**`locus`
- **关键词：**efp, ecotypes, microarray, arabidopsis, natural variation
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `locus` | string | 是 | Arabidopsis AGI locus, e.g. AT1G01010 |

??? note "上游说明（英文原文）"

    Fetch BAR/eFP world-map natural-variation expression for an Arabidopsis locus. Wraps the world-eFP view at /microarray_gene_expression/world_efp/arabidopsis/\{locus\} — returns expression across ~36 ecotypes (Bay-0, Col-0, Cvi-1, Ler-2, ...) with per-replicate values, control samples, collection lat/lng, and a per-ecotype mean computed client-side. Arabidopsis only. BAR is keyless and a Global Core Biodata Resource (2023).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'plant-genomics',
    'bar_efp_expression',
    {'locus': '...'},
)
```

## `bar_gene_summary` {#bar_gene_summary}

Fetches the TAIR/Araport11 curator summary, NCBI Gene ID and cross-database aliases for an Arabidopsis locus from BAR.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`locus`
- **关键词：**tair, araport11, thalemine, gene aliases, arabidopsis
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `locus` | string | 是 | Arabidopsis AGI locus, e.g. AT1G01010 |

??? note "上游说明（英文原文）"

    Fetch the BAR (Bio-Analytic Resource, U Toronto) merged ThaleMine + GAIA-aliases summary for an Arabidopsis locus. Returns the TAIR curator summary + Araport11 computational description from /thalemine/gene_information/ together with the NCBI Gene ID and cross-DB aliases (RefSeq, UniProt, TIGR locus-model IDs) from /gaia/aliases/. Arabidopsis only — ThaleMine carries taxon 3702 plus yeast/human for ortholog cross-reference. BAR is keyless and a Global Core Biodata Resource (2023); replaces the v0.9 subscription-gated tair_locus_info stub for the curator-summary use case.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'plant-genomics',
    'bar_gene_summary',
    {'locus': '...'},
)
```

## `batch_atted_coexpression` {#batch_atted_coexpression}

Batch version of ATTED-II co-expression lookup for up to 50 plant loci.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**`loci`
- **关键词：**atted-ii, coexpression, batch, plant
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `loci` | array of string | 是 |  |
| `organism` | string \| integer | 否 | Plant organism — accepts canonical slug (arabidopsis_thaliana), scientific or common name, or NCBI taxid 默认 `arabidopsis_thaliana` |
| `top_n` | integer | 否 | 默认 `25` |

??? note "上游说明（英文原文）"

    Batch version of atted_coexpression. Up to 50 loci per call.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'plant-genomics',
    'batch_atted_coexpression',
    {'loci': []},
)
```

## `batch_bar_gene_summary` {#batch_bar_gene_summary}

Batch fetch of BAR TAIR/Araport11 gene summaries and aliases for up to 50 Arabidopsis loci.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`loci`
- **关键词：**tair, araport11, gene aliases, batch, arabidopsis
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `loci` | array of string | 是 | List of locus identifiers (1–50). Successes land in results\[locus\]; PlantGenomicsError failures in errors\[locus\]. |

??? note "上游说明（英文原文）"

    Batch variant of bar_gene_summary. Fans out per-locus BAR ThaleMine + GAIA-aliases calls in parallel (up to 50 loci). Each results\[locus\] is the full single-locus payload (curator summary, computational description, NCBI Gene ID, cross-DB aliases). Arabidopsis only.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'plant-genomics',
    'batch_bar_gene_summary',
    {'loci': []},
)
```

## `batch_ensembl_plants_lookup_locus` {#batch_ensembl_plants_lookup_locus}

Batch lookup of Ensembl Plants gene metadata for up to 50 loci via POST /lookup/id.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`loci`
- **关键词：**ensembl plants, gene coordinates, biotype, batch
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `loci` | array of string | 是 | List of locus identifiers (1–50). Successes land in results\[locus\]; PlantGenomicsError failures in errors\[locus\]. |
| `organism` | string \| integer | 否 | Plant organism — accepts canonical slug (arabidopsis_thaliana), scientific or common name, or NCBI taxid 默认 `arabidopsis_thaliana` |

??? note "上游说明（英文原文）"

    Batch variant of ensembl_plants_lookup_locus. Uses Ensembl's native POST /lookup/id endpoint — one HTTP round-trip for up to 50 loci, materially cheaper than N parallel GETs. Successes in results\[\] with the same shape as the single-locus tool. Retries 429/5xx via the shared `_http` helper (Retry-After capped at 60 s). Misses (loci with no record) still land in `errors\[\]` with the `\[NotFoundError\]` prefix; the whole batch only fails when the retry budget is exhausted.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'plant-genomics',
    'batch_ensembl_plants_lookup_locus',
    {'loci': []},
)
```

## `batch_get_gene_xrefs` {#batch_get_gene_xrefs}

Batch fetch of Ensembl Plants cross-database references for up to 50 plant loci.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`loci`
- **关键词：**xrefs, uniprot, entrez, batch, ensembl plants
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `loci` | array of string | 是 | List of locus identifiers (1–50). Successes land in results\[locus\]; PlantGenomicsError failures in errors\[locus\]. |
| `organism` | string \| integer | 否 | Plant organism — accepts canonical slug (arabidopsis_thaliana), scientific or common name, or NCBI taxid 默认 `arabidopsis_thaliana` |

??? note "上游说明（英文原文）"

    Batch variant of get_gene_xrefs. Fans out per-locus xref lookups over Ensembl Plants in parallel (up to 50 loci). Each results\[locus\] is the full single-locus shape (count + xrefs\[\] + by_db rollup).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'plant-genomics',
    'batch_get_gene_xrefs',
    {'loci': []},
)
```

## `batch_gramene_homologs` {#batch_gramene_homologs}

Batch fetch of Gramene compara orthologs/paralogs for up to 50 plant loci.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`loci`
- **关键词：**gramene, orthologs, paralogs, gene tree, batch
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `loci` | array of string | 是 | List of locus identifiers (max 50) |
| `homology_type` | string | 否 | 可选值：`ortholog`, `paralog`, `all` 默认 `ortholog` |

??? note "上游说明（英文原文）"

    Batch version of gramene_homologs. Up to 50 loci per call; shares the homology_type filter across all loci. Returns the standard batch envelope (count + results dict + errors dict).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'plant-genomics',
    'batch_gramene_homologs',
    {'loci': []},
)
```

## `batch_kegg_pathways` {#batch_kegg_pathways}

Batch fetch of KEGG pathway memberships for up to 50 Arabidopsis loci.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`loci`
- **关键词：**kegg, pathways, batch, arabidopsis
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `loci` | array of string | 是 |  |
| `organism` | string \| integer | 否 | Plant organism — only arabidopsis_thaliana is supported in v1.1.0; other plants raise OrganismNotSupported until an Entrez bridge lands 默认 `arabidopsis_thaliana` |

??? note "上游说明（英文原文）"

    Batch version of kegg_pathways. Up to 50 loci per call. v1.1.0: only arabidopsis_thaliana resolves — KEGG uses NCBI Entrez Gene IDs for other plants and our cross-backend locus contract can't produce those yet, so a non-ath organism= raises OrganismNotSupported before any HTTP fan-out.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'plant-genomics',
    'batch_kegg_pathways',
    {'loci': []},
)
```

## `batch_locus_go_annotations` {#batch_locus_go_annotations}

Batch fetch of QuickGO Gene Ontology annotations for up to 50 plant loci via UniProt resolution.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`loci`
- **关键词：**go annotations, quickgo, batch, gene function
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `loci` | array of string | 是 | List of locus identifiers (1–50). Successes land in results\[locus\]; PlantGenomicsError failures in errors\[locus\]. |
| `organism` | string \| integer | 否 | Plant organism — accepts canonical slug (arabidopsis_thaliana), scientific or common name, or NCBI taxid 默认 `arabidopsis_thaliana` |
| `limit` | integer | 否 | Max annotations per locus from QuickGO (1–100, default 50) 默认 `50` |

??? note "上游说明（英文原文）"

    Batch variant of locus_go_annotations. Two-stage fanout — each locus is resolved to UniProt and then queried in QuickGO. Per-locus NotFoundError from either stage lands in errors\[\] with the typed prefix preserved. Capped at 50 loci.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'plant-genomics',
    'batch_locus_go_annotations',
    {'loci': []},
)
```

## `batch_phytozome_lookup_locus` {#batch_phytozome_lookup_locus}

Batch fetch of Phytozome BioMart gene records (chromosome, start/end, strand, description) for up to 50 loci.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`loci`
- **关键词：**phytozome, biomart, gene coordinates, batch
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `loci` | array of string | 是 | List of locus identifiers (1–50). Successes land in results\[locus\]; PlantGenomicsError failures in errors\[locus\]. |
| `organism` | string \| integer | 否 | Plant organism — accepts canonical slug (arabidopsis_thaliana), scientific or common name, or NCBI taxid 默认 `arabidopsis_thaliana` |

??? note "上游说明（英文原文）"

    Batch variant of phytozome_lookup_locus. Fans out per-locus BioMart queries in parallel (up to 50 loci). Each results\[locus\] is the full single-locus row (organism_name, gene_name, chromosome, start/end/strand, description).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'plant-genomics',
    'batch_phytozome_lookup_locus',
    {'loci': []},
)
```

## `batch_resolve_locus_to_uniprot` {#batch_resolve_locus_to_uniprot}

Batch resolution of up to 50 plant loci to canonical UniProtKB records.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`loci`
- **关键词：**uniprot, id mapping, swiss-prot, batch
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `loci` | array of string | 是 | List of locus identifiers (1–50). Successes land in results\[locus\]; PlantGenomicsError failures in errors\[locus\]. |
| `organism` | string \| integer | 否 | Plant organism — accepts canonical slug (arabidopsis_thaliana), scientific or common name, or NCBI taxid 默认 `arabidopsis_thaliana` |

??? note "上游说明（英文原文）"

    Batch variant of resolve_locus_to_uniprot. Fans out per-locus UniProtKB searches in parallel (up to 50 loci). Each results\[locus\] is the full single-locus record (primaryAccession + uniProtkbId + entryType + geneNames + organism + sequenceLength + web_url + …).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'plant-genomics',
    'batch_resolve_locus_to_uniprot',
    {'loci': []},
)
```

## `biological_context_synth` {#biological_context_synth}

Resolves a locus and gathers Gramene homologs, KEGG pathways, STRING partners and ATTED-II coexpression with a consensus partner ranking.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`locus`
- **关键词：**homologs, kegg, coexpression, functional partners, plant gene
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `locus` | string | 是 |  |
| `organism` | string \| integer | 否 | Plant organism — accepts canonical slug (arabidopsis_thaliana), scientific or common name, or NCBI taxid 默认 `arabidopsis_thaliana` |
| `top_n` | integer | 否 | 默认 `10` |

??? note "上游说明（英文原文）"

    Synthesis: one-call equivalent of the biological_context prompt. Resolves UniProt accession, then fans out to Gramene homologs, KEGG pathways, STRING-DB partners, and ATTED-II coexpression in parallel. Adds a consensus_partners ranking that merges STRING + ATTED scores.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'plant-genomics',
    'biological_context_synth',
    {'locus': '...'},
)
```

## `blast_sequence` {#blast_sequence}

Runs an NCBI BLAST URLAPI search (blastn/blastp/blastx/tblastn/tblastx) and returns parsed top hits.

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`sequence`
- **关键词：**ncbi blast, similarity search, core_nt, e-value, rid
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `sequence` | string | 是 | Raw or FASTA-formatted query sequence. |
| `program` | string | 否 | BLAST program — default blastp. 可选值：`blastn`, `blastp`, `blastx`, `tblastn`, `tblastx` 默认 `blastp` |
| `database` | string | 否 | NCBI BLAST database slug (e.g. swissprot, core_nt, refseq_protein). Defaults to swissprot for protein programs and core_nt for nucleotide programs. |
| `hitlist_size` | integer | 否 | Max hits to return (default 10). 默认 `10` |
| `expect` | number | 否 | E-value threshold (default 10). 默认 `10.0` |
| `megablast` | boolean | 否 | Enable megablast (blastn only). Default false. 默认 `False` |
| `poll_interval` | number | 否 | Seconds between polls. Clamped up to NCBI's per-RID 60s floor. 默认 `60.0` |
| `max_wait` | number | 否 | Max seconds to wait for the search to finish before raising NotFoundError with the RID preserved (default 600). 默认 `600.0` |

??? note "上游说明（英文原文）"

    Run a BLAST sequence-similarity search against NCBI BLAST URLAPI. Async Put/Get under the hood — submits the query, polls the RID (honoring NCBI's per-RID 60s floor), and returns the parsed top hits + raw text report excerpt. Programs: blastn / blastp / blastx / tblastn / tblastx. Database defaults to swissprot for protein programs, core_nt for nucleotide. Emits notifications/progress on each poll. Long searches (&gt;10 min) raise \[NotFoundError\] with the RID preserved so the client can re-poll. Set PLANT_GENOMICS_MCP_NCBI_EMAIL to identify the request per NCBI etiquette.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'plant-genomics',
    'blast_sequence',
    {'sequence': '...'},
)
```

## `consensus_homologs` {#consensus_homologs}

Builds a cross-source homology consensus for a locus by combining Gramene homologs and NCBI BLAST hits.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`locus`
- **关键词：**homologs, gramene, blast, orthologs, consensus
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `locus` | string | 是 |  |
| `organism` | string \| integer | 否 | Plant organism — accepts canonical slug (arabidopsis_thaliana), scientific or common name, or NCBI taxid 默认 `arabidopsis_thaliana` |
| `top_n` | integer | 否 | 默认 `10` |

??? note "上游说明（英文原文）"

    Synthesis: cross-source homology consensus. Resolves UniProt + FASTA sequence, then runs Gramene homology calls and NCBI BLAST in parallel. Dedupes hits by normalized locus token and scores by n_sources * mean_identity — Gramene contributes identity=1.0, BLAST contributes pident/100.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'plant-genomics',
    'consensus_homologs',
    {'locus': '...'},
)
```

## `ensembl_plants_lookup_locus` {#ensembl_plants_lookup_locus}

Fetches Ensembl Plants metadata for a plant locus identifier (coordinates, biotype, description).

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`locus`
- **关键词：**ensembl plants, gene coordinates, biotype, tair
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `locus` | string | 是 | e.g. AT1G01010 (Arabidopsis), Os01g0100100 (rice) |
| `organism` | string \| integer | 否 | Plant organism — accepts canonical slug (arabidopsis_thaliana), scientific or common name, or NCBI taxid 默认 `arabidopsis_thaliana` |

??? note "上游说明（英文原文）"

    Fetch metadata for a plant locus identifier from Ensembl Plants. Defaults to arabidopsis_thaliana; pass organism= for other plant species (oryza_sativa, zea_mays, ...). Locus is the TAIR-style identifier (e.g. AT1G01010 for Arabidopsis NAC001).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'plant-genomics',
    'ensembl_plants_lookup_locus',
    {'locus': '...'},
)
```

## `ensembl_region_query` {#ensembl_region_query}

Lists genes, transcripts, CDS or exons overlapping a genomic interval in an Ensembl Plants genome.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`region`, `start`, `end`
- **关键词：**genomic interval, overlap, qtl interval, exons, chromosome
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `region` | string | 是 | seq-region name (chromosome / contig), e.g. '1' or 'Chr1' |
| `start` | integer | 是 | 1-based start |
| `end` | integer | 是 | 1-based inclusive end |
| `organism` | string \| integer | 否 | Plant organism — accepts canonical slug (arabidopsis_thaliana), scientific or common name, or NCBI taxid 默认 `arabidopsis_thaliana` |
| `feature` | string | 否 | Feature type to return 可选值：`gene`, `transcript`, `cds`, `exon` 默认 `gene` |

??? note "上游说明（英文原文）"

    List features overlapping a genomic interval via Ensembl Plants /overlap/region. region is the seq-region name (chromosome / contig, e.g. '1'); start and end are 1-based inclusive. feature is one of gene / transcript / cds / exon (default gene). Answers 'what genes are in this QTL interval / assembly window' without a per-locus lookup. Ensembl caps the span — oversized regions error. Defaults to arabidopsis_thaliana; pass organism= for other species.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'plant-genomics',
    'ensembl_region_query',
    {'region': '...', 'start': 0, 'end': 0},
)
```

## `find_homologs_synth` {#find_homologs_synth}

Runs BLAST for a protein or nucleotide sequence and annotates ranked hits with UniProt records.

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`sequence`
- **关键词：**blast, homologs, uniprot, similarity search
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `sequence` | string | 是 | Query sequence (protein or nucleotide) |
| `program` | string | 否 | 可选值：`blastn`, `blastp`, `blastx`, `tblastn`, `tblastx` 默认 `blastp` |
| `top_n` | integer | 否 | 默认 `10` |

??? note "上游说明（英文原文）"

    Synthesis: one-call equivalent of the find_homologs prompt. Runs BLAST then resolves UniProt-shaped subject accessions via the batch UniProt helper. Returns ranked hits each annotated with their UniProt record (or null if subject_id is not a UniProt accession).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'plant-genomics',
    'find_homologs_synth',
    {'sequence': '...'},
)
```

## `gene_report` {#gene_report}

Generates a Markdown gene dossier for a plant locus from Ensembl Plants, UniProt, KEGG, STRING, Europe PMC and QuickGO.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`locus`
- **关键词：**gene dossier, plant gene, kegg, go, markdown report
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `locus` | string | 是 | Locus name, e.g. AT1G01010 |
| `organism` | string \| integer | 否 | Plant organism — accepts canonical slug (arabidopsis_thaliana), scientific or common name, or NCBI taxid 默认 `arabidopsis_thaliana` |
| `top_n` | integer | 否 | Caps GO terms, pathways, interactors, xrefs, and papers per section 默认 `10` |

??? note "上游说明（英文原文）"

    Synthesis: one-shot 'tell me about this gene' dossier. Resolves a locus through Ensembl Plants + UniProt, then fans out to cross-references, KEGG pathways, STRING interactors, Europe PMC literature, and QuickGO GO terms. Returns a SynthesisEnvelope whose result.markdown is a rendered Markdown gene dossier (the headline output) alongside a structured result.sections mirror. Any single backend failure degrades that section to an 'Unavailable' note; the rest of the dossier still renders.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'plant-genomics',
    'gene_report',
    {'locus': '...'},
)
```

## `get_gene_xrefs` {#get_gene_xrefs}

Fetches Ensembl Plants cross-database references (UniProt, NCBI Gene, TAIR, etc.) for a plant locus.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`locus`
- **关键词：**xrefs, uniprot, entrez, tair, ensembl plants
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `locus` | string | 是 | e.g. AT1G01010 (Arabidopsis), Os01g0100100 (rice) |
| `organism` | string \| integer | 否 | Plant organism — accepts canonical slug (arabidopsis_thaliana), scientific or common name, or NCBI taxid 默认 `arabidopsis_thaliana` |

??? note "上游说明（英文原文）"

    Fetch cross-database references (UniProt, NCBI Gene, TAIR, ArrayExpress, …) for a plant locus from Ensembl Plants. Defaults to arabidopsis_thaliana; pass organism= for other Ensembl Plants species. Returns count + raw xref list + a by_db rollup keyed on Ensembl's dbname (e.g. 'Uniprot_gn', 'EntrezGene') for fast lookup of a single foreign identifier.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'plant-genomics',
    'get_gene_xrefs',
    {'locus': '...'},
)
```

## `get_sequence` {#get_sequence}

Fetches a locus's genomic, CDS, cDNA or protein sequence from Ensembl Plants.

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`locus`
- **关键词：**genomic sequence, cds, cdna, fasta, ensembl plants
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `locus` | string | 是 | e.g. AT1G01010 (Arabidopsis), Os01g0100100 (rice) |
| `organism` | string \| integer | 否 | Plant organism — accepts canonical slug (arabidopsis_thaliana), scientific or common name, or NCBI taxid 默认 `arabidopsis_thaliana` |
| `seq_type` | string | 否 | Sequence type to fetch 可选值：`genomic`, `cds`, `cdna`, `protein` 默认 `protein` |

??? note "上游说明（英文原文）"

    Fetch a locus's sequence from Ensembl Plants. seq_type is one of genomic / cds / cdna / protein (default protein — the canonical-transcript product). Closes the lookup → fetch → BLAST loop: feed the returned `sequence` straight to blast_sequence (protein for blastp, cds/cdna for blastn). Defaults to arabidopsis_thaliana; pass organism= for other plant species.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'plant-genomics',
    'get_sequence',
    {'locus': '...'},
)
```

## `go_enrichment` {#go_enrichment}

Runs g:Profiler GO and KEGG over-representation analysis on a plant gene list with optional background.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`loci`
- **关键词：**g:profiler, enrichment, gene set, kegg, go
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `loci` | array of string | 是 | Query gene set, e.g. \['AT2G46830', 'AT1G01060', ...\] |
| `organism` | string \| integer | 否 | Plant organism — accepts canonical slug (arabidopsis_thaliana), scientific or common name, or NCBI taxid 默认 `arabidopsis_thaliana` |
| `sources` | array of string | 否 | Annotation sources to test (default: all four) |
| `background` | array of string | 否 | Optional custom statistical background gene set |
| `user_threshold` | number | 否 | Significance cutoff, g:SCS-corrected (default 0.05) 默认 `0.05` |
| `top_n` | integer | 否 | Max terms returned, sorted by p-value (1–200, default 50) 默认 `50` |

??? note "上游说明（英文原文）"

    GO + KEGG over-representation analysis for a gene LIST via g:Profiler g:GOSt (biit.cs.ut.ee/gprofiler; free, no API key). Unlike locus_go_annotations (one locus → its terms), this answers 'what is my gene SET enriched for?' — the dominant question for a differential-expression or co-expression cluster. loci is the query gene list (e.g. AT-codes for Arabidopsis, RAP-DB IDs for rice). sources defaults to GO:BP/GO:MF/GO:CC + KEGG; user_threshold is the g:SCS-corrected significance cutoff (default 0.05). Optional background sets a custom statistical domain (default: all annotated genes). Returns enriched\[\] (term_id/name/p_value/intersection_size/…, capped at top_n by p-value) plus unmapped\[\] — query loci g:Profiler could not recognize, surfaced so a locus-namespace mismatch is visible. Defaults to arabidopsis_thaliana; pass organism= for any of the 12 species.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'plant-genomics',
    'go_enrichment',
    {'loci': []},
)
```

## `gramene_homologs` {#gramene_homologs}

Fetches orthologs and paralogs for a plant locus from Gramene compara with gene-tree IDs.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`locus`
- **关键词：**gramene, orthologs, paralogs, gene tree
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `locus` | string | 是 | e.g. AT1G01010 (Arabidopsis), Os01g0100100 (rice) |
| `homology_type` | string | 否 | Filter on homology kind 可选值：`ortholog`, `paralog`, `all` 默认 `ortholog` |
| `limit` | integer | 否 | Max homolog rows to return. 'total' always reports the true pre-cap count and 'truncated' says whether the cap bit. 默认 `100` |

??? note "上游说明（英文原文）"

    Fetch orthologs and paralogs for a plant locus from Gramene compara (data.gramene.org v69). Default homology_type='ortholog'; pass 'paralog' for in-species duplicates or 'all' for everything. Returns target_locus + homology category (type) + shared gene_tree_id per hit. The fl=homology projection does not carry per-row taxon, identity, or protein ID; pair with resolve_locus_to_uniprot for protein-level enrichment and with blast_sequence for sequence similarity discovery.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'plant-genomics',
    'gramene_homologs',
    {'locus': '...'},
)
```

## `jaspar_motif` {#jaspar_motif}

Fetches one JASPAR TF binding profile by matrix ID, including the position-frequency matrix, TF class and consensus.

- **类别：**`regulation_epigenome`（Regulatory elements, TF binding, chromatin, methylation and 3D genome）
- **必填参数：**`matrix_id`
- **关键词：**jaspar, pfm, tf binding, motif, sequence logo
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `matrix_id` | string | 是 | JASPAR profile id, e.g. MA0570.1 or MA0570 (latest version) |

??? note "上游说明（英文原文）"

    Fetch one JASPAR binding profile by matrix id, including its raw position-frequency matrix (PFM: per-base count vectors keyed A/C/G/T) plus TF class/family, assay type, source species, UniProt accessions, PubMed refs, IUPAC consensus, and the sequence-logo URL. The drill-down companion to tf_binding_motifs, which returns the derived consensus but not the matrix. Accepts a versioned id (MA0570.1) or a bare base id (MA0570, which resolves to the newest version). Unknown ids raise a typed NotFoundError.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'plant-genomics',
    'jaspar_motif',
    {'matrix_id': '...'},
)
```

## `kegg_pathways` {#kegg_pathways}

Fetches KEGG pathway memberships for an Arabidopsis locus.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`locus`
- **关键词：**kegg, pathways, arabidopsis, metabolism
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `locus` | string | 是 | Arabidopsis AGI locus, e.g. AT1G01010 (case preserved verbatim — KEGG v118+ is case-sensitive) |
| `organism` | string \| integer | 否 | Plant organism — only arabidopsis_thaliana is supported in v1.1.0; other plants raise OrganismNotSupported until an Entrez bridge lands 默认 `arabidopsis_thaliana` |

??? note "上游说明（英文原文）"

    Fetch KEGG pathway memberships for an Arabidopsis locus from rest.kegg.jp. Returns a list of pathway IDs + names + KEGG category classes the locus participates in. Pairs with locus_go_annotations for the GO-level functional view. Multi-organism caveat (v1.1.0): the organism= field accepts any plant in the matrix for symmetry with the other backends, but only arabidopsis_thaliana resolves — KEGG uses NCBI Entrez Gene IDs for rice/maize/etc. and our cross-backend locus contract can't produce those yet, so any other organism raises OrganismNotSupported before any HTTP call. KEGG v118+ is case-sensitive on the locus: pass AGI loci as uppercase.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'plant-genomics',
    'kegg_pathways',
    {'locus': '...'},
)
```

## `locus_gene_rifs` {#locus_gene_rifs}

Fetches curated GeneRIF functional statements with PubMed IDs for an Arabidopsis locus from ThaleMine.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`locus`
- **关键词：**generif, gene function, thalemine, pubmed, arabidopsis
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `locus` | string | 是 | AGI locus, e.g. AT5G11260 (HY5) |
| `organism` | string \| integer | 否 | Plant organism — accepts canonical slug (arabidopsis_thaliana), scientific or common name, or NCBI taxid. ThaleMine supports Arabidopsis only. 默认 `arabidopsis_thaliana` |

??? note "上游说明（英文原文）"

    Fetch curated GeneRIF functional statements for an Arabidopsis locus from ThaleMine (free, no key). A GeneRIF is a one-sentence, manually curated statement of what the gene does, each anchored to the PubMed ID of the publication that demonstrated it — dense, directly citable functional context that GO terms (locus_go_annotations) and raw abstracts (locus_literature) do not provide. Well-studied genes have many: HY5 (AT5G11260) has 114. Upstream order is preserved because ThaleMine supplies no meaningful ranking, so `truncated` means later statements were cut, not that they were less relevant. found=false means the gene exists but has no GeneRIF; an unknown locus raises a typed NotFoundError. Arabidopsis only.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'plant-genomics',
    'locus_gene_rifs',
    {'locus': '...'},
)
```

## `locus_go_annotations` {#locus_go_annotations}

Fetches QuickGO Gene Ontology annotations for a plant locus, grouped by aspect.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`locus`
- **关键词：**go annotations, quickgo, gene function, molecular function
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `locus` | string | 是 | e.g. AT1G01010 (Arabidopsis), Os01g0100100 (rice) |
| `organism` | string \| integer | 否 | Plant organism — accepts canonical slug (arabidopsis_thaliana), scientific or common name, or NCBI taxid 默认 `arabidopsis_thaliana` |
| `limit` | integer | 否 | Max annotations from QuickGO (1–100, default 50) 默认 `50` |

??? note "上游说明（英文原文）"

    Fetch Gene Ontology annotations for a plant locus from QuickGO (EBI). Free, no API key. The locus is first resolved to a UniProt accession via the same logic as resolve_locus_to_uniprot, then QuickGO is queried by geneProductId. Returns annotations\[\] with goId/goName/goAspect/qualifier/evidence + a by_aspect rollup (\{molecular_function: \[\{goId, goName\}, ...\], biological_process: \[...\], cellular_component: \[...\]\}) deduped on goId so the high-level term set is one read away.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'plant-genomics',
    'locus_go_annotations',
    {'locus': '...'},
)
```

## `locus_plant_ontology` {#locus_plant_ontology}

Fetches Plant Ontology, Trait Ontology and PECO annotations for a plant locus from Planteome.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`locus`
- **关键词：**planteome, trait ontology, plant ontology, peco, gene annotation
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `locus` | string | 是 | e.g. AT1G01010 (Arabidopsis), Os01g0100100 (rice) |
| `organism` | string \| integer | 否 | Plant organism — accepts canonical slug (arabidopsis_thaliana), scientific or common name, or NCBI taxid 默认 `arabidopsis_thaliana` |
| `limit` | integer | 否 | Max annotations from Planteome (1–200, default 100) 默认 `100` |

??? note "上游说明（英文原文）"

    Fetch Plant Ontology (PO) + Trait Ontology (TO) + experimental-condition (PECO) annotations for a plant locus from Planteome (browser.planteome.org, AmiGO2/GOlr; free, no API key). Complements locus_go_annotations: QuickGO serves GO (species-agnostic), Planteome serves the plant-specific ontologies — PO (anatomy + developmental stage), TO (traits). The locus is matched across Planteome's searchable bioentity fields and filtered by the organism's NCBI taxon. Returns annotations\[\] (term_id / term_name / ontology / aspect / evidence / reference) + a by_ontology rollup (\{PO: \[\{term_id, term_name\}, ...\], TO: \[...\], PECO: \[...\]\}) deduped on term_id. Coverage is strong for arabidopsis, rice, maize, grape, soybean, tomato; other organisms return an empty list, not an error. Defaults to arabidopsis_thaliana; pass organism= for other species.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'plant-genomics',
    'locus_plant_ontology',
    {'locus': '...'},
)
```

## `locus_variants` {#locus_variants}

Lists EVA/dbSNP germline variants (SNPs, indels) overlapping a locus's genomic span via Ensembl.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`locus`
- **关键词：**snp, indel, eva, dbsnp, consequence
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `locus` | string | 是 | e.g. AT1G01010 (Arabidopsis), Os01g0100100 (rice RAP-DB) |
| `organism` | string \| integer | 否 | Plant organism — accepts canonical slug (arabidopsis_thaliana), scientific or common name, or NCBI taxid 默认 `arabidopsis_thaliana` |
| `limit` | integer | 否 | Max variant rows to return. 'variant_count' always reports the true pre-cap total and 'truncated' says whether the cap bit. 默认 `500` |

??? note "上游说明（英文原文）"

    List natural (germline) variants overlapping a locus's genomic span via Ensembl (rest.ensembl.org; free, no key). Resolves the locus → gene coordinates, then returns EVA/dbSNP-sourced SNPs and indels with id, source, consequence class, alleles, and clinical significance. variant_count is the true overlap total; the variant list is capped for payload size with truncated flagged. Opens the variation axis (distinct from get_sequence / ensembl_region_query). Works for all 12 organisms. Defaults to arabidopsis_thaliana; pass organism= for other species.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'plant-genomics',
    'locus_variants',
    {'locus': '...'},
)
```

## `orthodb_orthologs` {#orthodb_orthologs}

Resolves a locus to its OrthoDB ortholog group at Viridiplantae level and lists member genes by organism.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`locus`
- **关键词：**orthodb, orthologs, ortholog group, viridiplantae
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `locus` | string | 是 | e.g. AT1G01060 (Arabidopsis), Os01g0100100 (rice RAP-DB) |
| `organism` | string \| integer | 否 | Plant organism — accepts canonical slug (arabidopsis_thaliana), scientific or common name, or NCBI taxid. Validated and echoed only: it does NOT scope the OrthoDB search, which keys on the locus id at the Viridiplantae level 默认 `arabidopsis_thaliana` |
| `limit` | integer | 否 | Max ortholog member rows to return. 'member_count' always reports the true pre-cap total and 'truncated' says whether the cap bit. 默认 `100` |

??? note "上游说明（英文原文）"

    Resolve a locus to its OrthoDB ortholog group and cross-species member genes (data.orthodb.org; free, no key). Searches at the Viridiplantae level, then returns the group metadata (name, evolutionary rate) and member genes grouped by organism (organism, gene id, description). organism_count is the true cluster total; the member list is capped with truncated flagged. found=false when the locus maps to no ortholog group. Works for all 12 organisms. NOTE: unlike the other locus tools, organism= does NOT scope the search — the group is resolved from the locus id alone at the Viridiplantae level, and organism is only validated and echoed back. Passing a mismatched organism therefore still returns the locus's real group.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'plant-genomics',
    'orthodb_orthologs',
    {'locus': '...'},
)
```

## `panther_family` {#panther_family}

Fetches PANTHER family/subfamily classification, curated GO terms, protein class and pathways for a locus.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`locus`
- **关键词：**panther, gene family, go, protein class, pathways
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `locus` | string | 是 | e.g. AT1G01060 (Arabidopsis), Os01g0100100 (rice RAP-DB) |
| `organism` | string \| integer | 否 | Plant organism — accepts canonical slug (arabidopsis_thaliana), scientific or common name, or NCBI taxid 默认 `arabidopsis_thaliana` |

??? note "上游说明（英文原文）"

    Fetch the PANTHER protein-family classification for a locus (pantherdb.org; free, no key). Returns the PANTHER family and subfamily (id + name) plus curated GO terms grouped by aspect (molecular_function / biological_process / cellular_component), the PANTHER protein class, and pathways. found=false when PANTHER cannot classify the locus. Complements the sequence-homology tools (gramene_homologs / consensus_homologs) with an evolutionary-family view. Works for all 12 organisms. Defaults to arabidopsis_thaliana; pass organism= for other species.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'plant-genomics',
    'panther_family',
    {'locus': '...'},
)
```

## `phytozome_lookup_locus` {#phytozome_lookup_locus}

Fetches a Phytozome BioMart gene record (chromosome, start/end, strand, description) for a locus.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`locus`
- **关键词：**phytozome, biomart, gene coordinates, soybean
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `locus` | string | 是 | e.g. AT1G01010 (Arabidopsis), Glyma.01G000100 (soybean) |
| `organism` | string \| integer | 否 | Plant organism — accepts canonical slug (arabidopsis_thaliana), scientific or common name, or NCBI taxid 默认 `arabidopsis_thaliana` |

??? note "上游说明（英文原文）"

    Fetch a gene record from Phytozome BioMart (phytozome-next.jgi.doe.gov). Defaults to arabidopsis_thaliana; pass organism= for other Phytozome proteomes (slug, scientific/common name, or NCBI taxid — e.g. glycine_max, sorghum_bicolor). Locus is the source-genome gene name (e.g. AT1G01010, Glyma.01G000100). Returns organism_name, gene_name, chromosome, gene_start, gene_end, strand, description.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'plant-genomics',
    'phytozome_lookup_locus',
    {'locus': '...'},
)
```

## `plantcyc_locus_info` {#plantcyc_locus_info}

Fetches PlantCyc metabolic annotation for a locus: enzymes, catalyzed reactions and pathways.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`locus`
- **关键词：**plantcyc, metabolic pathway, enzymes, reactions, pmn
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `locus` | string | 是 | e.g. AT3G51240 (Arabidopsis), Os11g0530600 (rice RAP-DB) |
| `organism` | string \| integer | 否 | Plant organism — accepts canonical slug (arabidopsis_thaliana), scientific or common name, or NCBI taxid 默认 `arabidopsis_thaliana` |

??? note "上游说明（英文原文）"

    Fetch metabolic annotation for a locus from PlantCyc / the Plant Metabolic Network (pmn.plantcyc.org; free BioCyc web-services API, no key). Walks gene → enzyme → catalyzed reactions → PlantCyc pathways in the organism's PGDB, returning enzymes\[\] + reactions\[\] (id/name) + pathways\[\] (id/name) — the metabolic-pathway view KEGG and GO don't provide. A non-enzymatic gene (e.g. a transcription factor) returns found=false with empty lists, not an error. reaction_count / pathway_count report true totals even when the lists are capped. 11 organisms have a PGDB (arabidopsis, rice, maize, soybean, grape, poplar, tomato, barley, sorghum, medicago, brachypodium); wheat is not yet mapped. Defaults to arabidopsis_thaliana (AraCyc, the best-curated); pass organism= for other species.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'plant-genomics',
    'plantcyc_locus_info',
    {'locus': '...'},
)
```

## `resolve_locus_to_uniprot` {#resolve_locus_to_uniprot}

Resolves a plant locus to its canonical UniProtKB record (accession, gene names, sequence length).

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`locus`
- **关键词：**uniprot, id mapping, swiss-prot, trembl
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `locus` | string | 是 | e.g. AT1G01010 (Arabidopsis), Os01g0100100 (rice) |
| `organism` | string \| integer | 否 | Plant organism — accepts canonical slug (arabidopsis_thaliana), scientific or common name, or NCBI taxid 默认 `arabidopsis_thaliana` |

??? note "上游说明（英文原文）"

    Resolve a plant locus to its canonical UniProtKB record. Prefers reviewed (Swiss-Prot) entries; falls back to unreviewed (TrEMBL) when no curated record exists (common for non-Arabidopsis plants). organism accepts a canonical slug, scientific/common name, or NCBI taxid (default arabidopsis_thaliana; e.g. oryza_sativa, zea_mays). Returns primaryAccession, uniProtkbId, entryType, recommendedName, geneNames, organism, taxonId, sequenceLength, web_url. This is the protein-side entry point — pair with InterPro / AlphaFold / Reactome / structural-bio tools.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'plant-genomics',
    'resolve_locus_to_uniprot',
    {'locus': '...'},
)
```

## `tair_locus_info` {#tair_locus_info}

Fetches the TAIR curator summary, Araport11 description, NCBI Gene ID and aliases for an Arabidopsis locus via BAR.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`locus`
- **关键词：**tair, araport11, gene summary, arabidopsis
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `locus` | string | 是 | Arabidopsis AGI locus, e.g. AT1G01010 |

??? note "上游说明（英文原文）"

    Fetch the TAIR curator-vetted Arabidopsis locus summary. Served via BAR/ThaleMine (U Toronto, Global Core Biodata Resource 2023) since TAIR's free per-locus REST API is gated behind a paid Phoenix Bioinformatics subscription. Returns TAIR curator summary + Araport11 computational description + NCBI Gene ID + cross-DB aliases (RefSeq, UniProt, TIGR locus-model IDs). Arabidopsis only. Alias of bar_gene_summary.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'plant-genomics',
    'tair_locus_info',
    {'locus': '...'},
)
```

## `tf_binding_motifs` {#tf_binding_motifs}

Fetches JASPAR transcription-factor DNA binding motifs for a locus, confirmed by UniProt accession, with consensus and assay type.

- **类别：**`regulation_epigenome`（Regulatory elements, TF binding, chromatin, methylation and 3D genome）
- **必填参数：**`locus`
- **关键词：**jaspar, tf binding, motif, cis-regulatory, consensus
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `locus` | string | 是 | e.g. AT2G46830 (Arabidopsis CCA1), Os01g0100100 (rice RAP-DB) |
| `organism` | string \| integer | 否 | Plant organism — accepts canonical slug (arabidopsis_thaliana), scientific or common name, or NCBI taxid 默认 `arabidopsis_thaliana` |

??? note "上游说明（英文原文）"

    Fetch curated transcription-factor DNA binding motifs for a locus from JASPAR (jaspar.elixir.no; free, no key) — the cis-regulatory view. Resolves the locus → UniProt accession + gene symbol, searches JASPAR by symbol scoped to the organism's taxid, then CONFIRMS each candidate by matching the accession against the profile's uniprot_ids. Returns per motif the JASPAR matrix id, TF class/family, assay type (SELEX / ChIP-seq / PBM / DAP-seq), an IUPAC consensus derived from the position-frequency matrix (e.g. CACGTG, the G-box/ABRE core), motif length, PubMed refs, and an SVG sequence-logo URL. IMPORTANT: JASPAR's name search is fuzzy, so name-similarity hits belonging to a DIFFERENT gene are returned separately in name_only_matches and must NOT be attributed to this locus; only `motifs` is UniProt-confirmed. found=false means the gene has no curated profile (not a TF, or its family is unprofiled for that species) — a normal outcome, not an error. Use jaspar_motif to retrieve the raw matrix for any matrix_id. Coverage is Arabidopsis-heavy (1236 profiles) and thin elsewhere (maize 131, soybean 91, wheat 58, tomato 51, rice 10; Brachypodium and sorghum have none). Defaults to arabidopsis_thaliana; pass organism= for other species.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'plant-genomics',
    'tf_binding_motifs',
    {'locus': '...'},
)
```

## `vep_annotate` {#vep_annotate}

Predicts a variant's molecular consequences per overlapping transcript with Ensembl VEP, including SIFT/PolyPhen.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`region`, `allele`
- **关键词：**vep, consequence, sift, polyphen, impact
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `region` | string | 是 | Ensembl region chr:start-end:strand, e.g. '1:10000-10000:1' |
| `allele` | string | 是 | Alternate allele, e.g. 'C' (or 'A/C', an insertion, etc.) |
| `organism` | string \| integer | 否 | Plant organism — accepts canonical slug (arabidopsis_thaliana), scientific or common name, or NCBI taxid 默认 `arabidopsis_thaliana` |

??? note "上游说明（英文原文）"

    Predict a variant's molecular consequences with Ensembl VEP (rest.ensembl.org; free, no key). Variant-first (not locus-first): supply an Ensembl region (chr:start-end:strand, e.g. '1:10000-10000:1') and an alternate allele (e.g. 'C'); returns the most-severe consequence plus one row per overlapping transcript (consequence terms, IMPACT, and SIFT/PolyPhen when the variant is coding-missense). found=false when Ensembl reports no overlapping feature. Works for all 12 organisms. Defaults to arabidopsis_thaliana; pass organism= for other species.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'plant-genomics',
    'vep_annotate',
    {'region': '...', 'allele': '...'},
)
```
