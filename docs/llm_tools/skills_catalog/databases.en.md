# Genomic databases

**43** Skills. A Skill is a set of instructions that guides an Agent through an analysis; the software, reference data and APIs it uses are separate dependencies. See [Skills](../skills.md) for searching, installing and using them; the Skills tab of the web demo shows each Skill in full.

| Skill | Source | Summary |
| --- | --- | --- |
| [`article-data-fetcher`](#article-data-fetcher) | clawbio | Given an article DOI or PubMed ID, discover and download the genomics data files deposited by the authors (VCF, FASTA, H5AD, CSV, JSON, BAM, etc.) from public… |
| [`batch-downloads`](#batch-downloads) | bioskills | Download large datasets from NCBI efficiently using EPost, history server, batching, rate limiting, and retry logic. Use when bulk-fetching tens of thousands o… |
| [`biobakery-databases`](#biobakery-databases) | pku-metagenomics | Install and pin bioBakery reference databases for MetaPhlAn and HUMAnN (ChocoPhlAn markers, UniRef gene families, MetaCyc pathways, utility mapping). Use this… |
| [`biomart-queries`](#biomart-queries) | bioskills | Bulk-query Ensembl BioMart (and other BioMart instances) for cross-database ID mapping, gene/transcript/exon coordinates, and ortholog tables. Use when batch-c… |
| [`bioservices-multi-database`](#bioservices-multi-database) | sciagent | Unified Python interface to 40+ bioinformatics web services: UniProt proteins, KEGG pathways, ChEMBL/ChEBI/PubChem, BLAST, cross-database ID mapping, GO annota… |
| [`blast-searches`](#blast-searches) | bioskills | Run remote BLAST searches against NCBI servers using Biopython Bio.Blast.NCBIWWW. Use when identifying unknown sequences, finding homologs, picking the correct… |
| [`cbioportal-database`](#cbioportal-database) | sciagent | Cancer genomics (TCGA et al.) via cBioPortal REST API. Retrieve somatic mutations, CNAs, expression, clinical data (survival/stage/treatment) across thousands… |
| [`cosmic-database`](#cosmic-database) | sciagent | Query COSMIC for cancer somatic mutations, gene census, mutational signatures, drug resistance variants. REST API v3.1 supports gene/sample/variant queries; fr… |
| [`depmap-crispr-essentiality`](#depmap-crispr-essentiality) | sciagent | DepMap CRISPR gene effect (Chronos) analysis: sign convention for essentiality, per-gene NaN-safe Spearman correlation, data loading/alignment. For general NaN… |
| [`ena-database`](#ena-database) | sciagent | ENA REST API for sequences, reads, assemblies, and annotations. Portal API search, Browser API retrieval (XML/FASTA/EMBL), file reports for FASTQ/BAM URLs, tax… |
| [`encode-ccres-database`](#encode-ccres-database) | deepmind | Query the ENCODE Registry of cis-Regulatory Elements (cCREs) via the SCREEN GraphQL API, or make custom queries to the ENCODE Portal REST API for experiments a… |
| [`encode-database`](#encode-database) | sciagent | ENCODE Portal REST API for regulatory genomics: TF ChIP-seq, ATAC-seq/DNase-seq peaks, histone marks, and RNA-seq across 1000+ cell types. Search experiments b… |
| [`ensembl-rest`](#ensembl-rest) | bioskills | Query the Ensembl REST API for gene/transcript/protein lookup, sequence retrieval, comparative genomics (Compara), variant effect prediction (VEP), regulatory… |
| [`entrez-fetch`](#entrez-fetch) | bioskills | Retrieve records from NCBI databases using Biopython Bio.Entrez (EFetch, ESummary). Use when downloading sequences, fetching GenBank/GenPept records, getting d… |
| [`entrez-link`](#entrez-link) | bioskills | Find cross-database references between NCBI databases using Biopython Bio.Entrez (ELink). Use when navigating gene to protein/structure, sequence to publicatio… |
| [`entrez-search`](#entrez-search) | bioskills | Search NCBI databases using Biopython Bio.Entrez (ESearch, EInfo, EGQuery, ESpell). Use when finding records by keyword, building reproducible field-qualified… |
| [`eqtl-catalogue-region-fetch`](#eqtl-catalogue-region-fetch) | clawbio | Fetch a region of cis-eQTL summary statistics from EBI eQTL Catalogue v7+ via tabix-on-FTP. Use when an agent needs eQTL beta / SE / p-value for every variant… |
| [`gene-database`](#gene-database) | sciagent | NCBI Gene via E-utilities: curated records across 1M+ taxa. Official symbols, aliases, RefSeq IDs, summaries, coordinates, GO, interactions. Use for gene ID re… |
| [`geo-data`](#geo-data) | bioskills | Query and download from NCBI Gene Expression Omnibus (GEO) and EMBL-EBI's BioStudies/ArrayExpress mirror. Use when finding expression datasets, navigating Supe… |
| [`gget-genomic-databases`](#gget-genomic-databases) | sciagent | Unified CLI/Python interface to 20+ genomic databases. Gene lookups (Ensembl search/info/seq), BLAST/BLAT, AlphaFold, Enrichr enrichment, OpenTargets disease/d… |
| [`gigwa-explore-instance`](#gigwa-explore-instance) | gigwa-mcp | Get an overview of an entire Gigwa instance and flag anything that needs attention. Use this when a user asks what data is available on a Gigwa server, wants t… |
| [`gtdb`](#gtdb) | pku-metagenomics | Use the Genome Taxonomy Database (GTDB) public releases — taxonomy tables, trees, metadata, and representative genomes. Use this skill when: (1) Choosing or ci… |
| [`gwas-catalog-region-fetch`](#gwas-catalog-region-fetch) | clawbio | Fetch a region of GWAS summary statistics from the NHGRI-EBI GWAS Catalog harmonised collection via tabix-on-FTP. Use when an agent needs GWAS beta / SE / p-va… |
| [`gwas-database`](#gwas-database) | sciagent | NHGRI-EBI GWAS Catalog REST API for SNP-trait associations from published GWAS. Query studies, associations, variants, traits, genes, summary stats. Build PRS… |
| [`gwas-fetch`](#gwas-fetch) | statgen-skills | Fetch harmonised GWAS summary statistics from the GWAS Catalog by accession (GCST*). Resolves the study's metadata (genome build, FTP URL, sample size, ancestr… |
| [`jaspar-database`](#jaspar-database) | sciagent | JASPAR 2024 TF binding profiles via REST API and pyJASPAR. Retrieve PFMs/PWMs by TF name, JASPAR ID, species, or structural class. Scan DNA for TFBS; browse by… |
| [`local-blast`](#local-blast) | bioskills | Build local BLAST databases and run searches using NCBI BLAST+ command-line tools. Use when running &gt;50 queries, building custom databases with -parse_seqids a… |
| [`mgnify`](#mgnify) | pku-metagenomics | Use EMBL-EBI MGnify for analysed microbiome studies and biome genome catalogues (including UHGG / human-gut and other biomes). Use this skill when: (1) Browsin… |
| [`model-organism-genetics`](#model-organism-genetics) | tooluniverse | Cross-species genetic analysis using model organism databases (MGI mouse, ZFIN zebrafish, FlyBase fruit fly, WormBase worm, SGD yeast, RGD rat, GBIF taxonomy).… |
| [`monarch-database`](#monarch-database) | sciagent | Monarch Initiative knowledge graph REST API for disease-gene-phenotype associations and cross-species orthology. MONDO disease-to-gene/phenotype, HP phenotype… |
| [`mouse-phenome-database`](#mouse-phenome-database) | sciagent | Retrieve mouse phenotype data from the Jackson Laboratory Mouse Phenome Database (MPD) via its REST API. Browse 520+ projects, look up per-project measure meta… |
| [`ncbi-datasets-cli`](#ncbi-datasets-cli) | bioskills | Download genome assemblies, gene records, and ortholog data from NCBI using the modern Datasets v2 CLI (replaces assembly_summary.txt scraping and many EFetch… |
| [`okg-search`](#okg-search) | statgen-skills | Query the statgen-analysis OKG for datasets matching a trait, EFO term, seed-dataset, or curated category, with controlled-vocab ancestry filtering (EUR/AFR/EA… |
| [`ortholog-inference`](#ortholog-inference) | bioskills | Pull pre-computed ortholog calls from public databases (OrthoDB, Ensembl Compara, OMA browser, eggNOG, PANTHER, KEGG Orthology, HomoloGene) via their REST APIs… |
| [`pdb-database`](#pdb-database) | user335-deepmind | Use when you want to search for or download experimentally-determined 3D structures for biomolecules (proteins, nucleic acids, bound ligands). Supports searchi… |
| [`population-genetics-1000genomes`](#population-genetics-1000genomes) | tooluniverse | Population genetics using the 1000 Genomes Project (IGSR) — superpopulation/population search, sample metadata, variant frequencies across AFR/AMR/EAS/EUR/SAS,… |
| [`progenomes`](#progenomes) | pku-metagenomics | proGenomes4 — nearly 2 million consistently QC'd and annotated high-quality prokaryotic genomes with habitat metadata, pan-genomes, MGE/BGC annotations. Use wh… |
| [`regulomedb-database`](#regulomedb-database) | sciagent | Query RegulomeDB v2 GET REST API to score variants for regulatory function and retrieve overlapping evidence (TF binding, histone marks, DNase peaks, footprint… |
| [`remap-database`](#remap-database) | sciagent | Query ReMap 2022 TF ChIP-seq peak database via REST API and BED downloads. Retrieve TF peaks overlapping a region (chr:start-end), peaks near a gene, TFs by sp… |
| [`sequence-retrieval`](#sequence-retrieval) | tooluniverse | Retrieve DNA, RNA and protein sequences with DNAKit from NCBI and ENA. Use for accession/version lookup, organism and gene disambiguation, transcript isoforms,… |
| [`sra-data`](#sra-data) | bioskills | Download raw sequencing reads from NCBI SRA using sra-tools (prefetch, fasterq-dump, vdb-validate) or the ENA mirror. Use when pulling FASTQ for SRR/ERR/DRR ac… |
| [`ucsc-genome-browser`](#ucsc-genome-browser) | sciagent | Query UCSC Genome Browser REST API for DNA sequences, tracks, gene models, and conservation across 100+ assemblies. Retrieve sequence by region, list/fetch BED… |
| [`unibind-database`](#unibind-database) | deepmind | Queries the UniBind database for experimentally validated transcription factor (TF) binding sites. Use when retrieving direct TF-DNA interaction datasets, down… |

## `article-data-fetcher` {#article-data-fetcher}

Given an article DOI or PubMed ID, discover and download the genomics data files deposited by the authors (VCF, FASTA, H5AD, CSV, JSON, BAM, etc.) from public repositories such as GEO, ENA, Zenodo, Figshare, Dryad, and OSF.

- **Source:** [ClawBio; ClawBio contributors @ `1fcecb7`](https://github.com/ClawBio/ClawBio/tree/1fcecb7e9531210baabe6d5cdbdd08de1ff3ab8d/skills/article-data-fetcher)
- **License:** `MIT`
- **Runs with:** commands `python3`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info article-data-fetcher
dnakit skills install article-data-fetcher --agent claude-code   # or --agent codex / cursor
```

## `batch-downloads` {#batch-downloads}

Download large datasets from NCBI efficiently using EPost, history server, batching, rate limiting, and retry logic. Use when bulk-fetching tens of thousands of sequences, pulling all results of a large ESearch, designing reproducible pipelines, comparing E-utilities to NCBI Datasets v2 CLI, or implementing checksum-validated downloads. Encodes WebEnv TTL (~8h), EPost 200-ID limit, retmax caps, parallelization design, and integrity verification.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/database-access/batch-downloads)
- **License:** `MIT`
- **Runs with:** Python modules `Bio.Entrez`, `dnakit`
- **Also known as:** `bio-batch-downloads`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info batch-downloads
dnakit skills install batch-downloads --agent claude-code   # or --agent codex / cursor
```

## `biobakery-databases` {#biobakery-databases}

Install and pin bioBakery reference databases for MetaPhlAn and HUMAnN (ChocoPhlAn markers, UniRef gene families, MetaCyc pathways, utility mapping). Use this skill when: (1) Downloading MetaPhlAn/HUMAnN databases, (2) Choosing UniRef90 vs UniRef50, (3) Aligning MetaPhlAn and HUMAnN DB versions, (4) Troubleshooting missing index paths. Upstream: https://github.com/biobakery/MetaPhlAn and https://github.com/biobakery/humann. For Kraken indexes use kraken2. For eggNOG/CARD/GTDB use those skills. For routing use public-databases.

- **Source:** [Copyright (c) 2026 PKU-EMBL @ `e0e19a6`](https://github.com/PKU-EMBL/Metagenomics-Skills/tree/e0e19a600d9916736c0713e7659dec3782c51258/skills/biobakery-databases)
- **License:** `MIT`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info biobakery-databases
dnakit skills install biobakery-databases --agent claude-code   # or --agent codex / cursor
```

## `biomart-queries` {#biomart-queries}

Bulk-query Ensembl BioMart (and other BioMart instances) for cross-database ID mapping, gene/transcript/exon coordinates, and ortholog tables. Use when batch-converting Ensembl IDs to other namespaces (HGNC, RefSeq, UniProt, Entrez), pulling gene coordinate tables for thousands of genes, building ortholog wide-tables across species, or replacing slow Ensembl REST loops with one-shot bulk export. Encodes BioMart's XML query format, R biomaRt vs Python pybiomart trade-off, mart-vs-dataset hierarchy, and the URL endpoint that's BioMart-specific (separate from rest.ensembl.org).

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/database-access/biomart-queries)
- **License:** `MIT`
- **Also known as:** `bio-biomart-queries`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info biomart-queries
dnakit skills install biomart-queries --agent claude-code   # or --agent codex / cursor
```

## `bioservices-multi-database` {#bioservices-multi-database}

Unified Python interface to 40+ bioinformatics web services: UniProt proteins, KEGG pathways, ChEMBL/ChEBI/PubChem, BLAST, cross-database ID mapping, GO annotations, PPI. For deep single-DB queries use dedicated tools (gget for Ensembl, pubchempy for PubChem); bioservices excels at cross-database workflows.

- **Source:** [jaechang-hits; credited K-Dense-AI and Biomni contributors @ `fe505ca`](https://github.com/jaechang-hits/SciAgent-Skills/tree/fe505cae14d20b6c33be2e49666425be98f005bb/skills/genomics-bioinformatics/databases/bioservices-multi-database)
- **License:** `CC-BY-4.0`
- **Also known as:** `bioservices`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info bioservices-multi-database
dnakit skills install bioservices-multi-database --agent claude-code   # or --agent codex / cursor
```

## `blast-searches` {#blast-searches}

Run remote BLAST searches against NCBI servers using Biopython Bio.Blast.NCBIWWW. Use when identifying unknown sequences, finding homologs, picking the correct BLAST program (blastn/blastp/blastx/tblastn/tblastx/psiblast/megablast/dc-megablast), interpreting Karlin-Altschul E-values, avoiding the max_target_seqs trap (Shah 2019), choosing composition-based statistics, or limiting searches by organism. Covers RID lifecycle, database choice (nt/nr/refseq_select/swissprot), word-size and CBS taxonomy.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/database-access/blast-searches)
- **License:** `MIT`
- **Runs with:** Python modules `Bio.Blast.NCBIWWW`, `dnakit`
- **Also known as:** `bio-blast-searches`, `blast-search`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info blast-searches
dnakit skills install blast-searches --agent claude-code   # or --agent codex / cursor
```

## `cbioportal-database` {#cbioportal-database}

Cancer genomics (TCGA et al.) via cBioPortal REST API. Retrieve somatic mutations, CNAs, expression, clinical data (survival/stage/treatment) across thousands of studies. Use for TMB, oncoprints, survival analysis. For population frequencies use gnomad-database; for drug-gene interactions use opentargets-database.

- **Source:** [jaechang-hits; credited K-Dense-AI and Biomni contributors @ `fe505ca`](https://github.com/jaechang-hits/SciAgent-Skills/tree/fe505cae14d20b6c33be2e49666425be98f005bb/skills/genomics-bioinformatics/databases/cbioportal-database)
- **License:** `CC-BY-4.0`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info cbioportal-database
dnakit skills install cbioportal-database --agent claude-code   # or --agent codex / cursor
```

## `cosmic-database` {#cosmic-database}

Query COSMIC for cancer somatic mutations, gene census, mutational signatures, drug resistance variants. REST API v3.1 supports gene/sample/variant queries; free registration. For germline use clinvar-database; for drug-target data use opentargets-database or chembl-database-bioactivity.

- **Source:** [jaechang-hits; credited K-Dense-AI and Biomni contributors @ `fe505ca`](https://github.com/jaechang-hits/SciAgent-Skills/tree/fe505cae14d20b6c33be2e49666425be98f005bb/skills/genomics-bioinformatics/databases/cosmic-database)
- **License:** `CC-BY-4.0`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info cosmic-database
dnakit skills install cosmic-database --agent claude-code   # or --agent codex / cursor
```

## `depmap-crispr-essentiality` {#depmap-crispr-essentiality}

DepMap CRISPR gene effect (Chronos) analysis: sign convention for essentiality, per-gene NaN-safe Spearman correlation, data loading/alignment. For general NaN-safe correlation see nan-safe-correlation; for quality filtering see degenerate-input-filtering.

- **Source:** [jaechang-hits; credited K-Dense-AI and Biomni contributors @ `fe505ca`](https://github.com/jaechang-hits/SciAgent-Skills/tree/fe505cae14d20b6c33be2e49666425be98f005bb/skills/genomics-bioinformatics/databases/depmap-crispr-essentiality)
- **License:** `CC-BY-4.0`
- **Also known as:** `depmap`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info depmap-crispr-essentiality
dnakit skills install depmap-crispr-essentiality --agent claude-code   # or --agent codex / cursor
```

## `ena-database` {#ena-database}

ENA REST API for sequences, reads, assemblies, and annotations. Portal API search, Browser API retrieval (XML/FASTA/EMBL), file reports for FASTQ/BAM URLs, taxonomy, cross-refs. For multi-DB Python use bioservices; for NCBI-only use pubmed-database or Biopython Entrez.

- **Source:** [jaechang-hits; credited K-Dense-AI and Biomni contributors @ `fe505ca`](https://github.com/jaechang-hits/SciAgent-Skills/tree/fe505cae14d20b6c33be2e49666425be98f005bb/skills/genomics-bioinformatics/databases/ena-database)
- **License:** `CC-BY-4.0`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info ena-database
dnakit skills install ena-database --agent claude-code   # or --agent codex / cursor
```

## `encode-ccres-database` {#encode-ccres-database}

Query the ENCODE Registry of cis-Regulatory Elements (cCREs) via the SCREEN GraphQL API, or make custom queries to the ENCODE Portal REST API for experiments and files (ChIP-seq peaks, etc.). Use when you want to query regulatory annotations or raw experimental data across human cell types.

- **Source:** [Copyright 2026 Google LLC @ `28b8482`](https://github.com/google-deepmind/science-skills/tree/28b8482603a420708c8896f6fe5e06c276d9933d/skills/encode_ccres_database)
- **License:** `Apache-2.0 AND CC-BY-4.0`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info encode-ccres-database
dnakit skills install encode-ccres-database --agent claude-code   # or --agent codex / cursor
```

## `encode-database` {#encode-database}

ENCODE Portal REST API for regulatory genomics: TF ChIP-seq, ATAC-seq/DNase-seq peaks, histone marks, and RNA-seq across 1000+ cell types. Search experiments by assay/biosample/target; download BED/bigWig; retrieve SCREEN cCREs by region or gene. Use to annotate variants with regulatory tracks, find open chromatin in a cell type, or fetch peak files for ChIP/ATAC analysis. For regulatory variant scoring use regulomedb-database; for GWAS associations use gwas-database.

- **Source:** [jaechang-hits; credited K-Dense-AI and Biomni contributors @ `fe505ca`](https://github.com/jaechang-hits/SciAgent-Skills/tree/fe505cae14d20b6c33be2e49666425be98f005bb/skills/genomics-bioinformatics/databases/encode-database)
- **License:** `CC-BY-4.0`
- **Also known as:** `download-encode`, `search-encode`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info encode-database
dnakit skills install encode-database --agent claude-code   # or --agent codex / cursor
```

## `ensembl-rest` {#ensembl-rest}

Query the Ensembl REST API for gene/transcript/protein lookup, sequence retrieval, comparative genomics (Compara), variant effect prediction (VEP), regulatory features, and cross-species ortholog/paralog calls. Use when pulling Ensembl-native data (Ensembl Gene IDs, version-pinned releases, archive endpoints for reproducibility), gene/transcript/exon structure with stable IDs, or VEP for variant annotation. Encodes the 15 req/sec rate limit, archive (e110.rest.ensembl.org) for reproducibility, Ensembl divisions (vertebrates / plants / fungi / metazoa / bacteria), and the symbol-vs-ID stability problem.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/database-access/ensembl-rest)
- **License:** `MIT`
- **Runs with:** Python modules `dnakit`, `requests`
- **Also known as:** `bio-ensembl-rest`, `ensembl-annotation`, `ensembl-database`, `query-ensembl`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info ensembl-rest
dnakit skills install ensembl-rest --agent claude-code   # or --agent codex / cursor
```

## `entrez-fetch` {#entrez-fetch}

Retrieve records from NCBI databases using Biopython Bio.Entrez (EFetch, ESummary). Use when downloading sequences, fetching GenBank/GenPept records, getting document summaries, parsing nested XML, navigating GI deprecation, choosing between rettype+retmode combinations, and parsing into Biopython SeqRecord/SwissProt objects. Covers nucleotide, protein, gene, pubmed, sra, gds, taxonomy, snp, clinvar.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/database-access/entrez-fetch)
- **License:** `MIT`
- **Runs with:** Python modules `Bio.Entrez`, `dnakit`
- **Also known as:** `bio-entrez-fetch`, `ncbi-sequence-fetch`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info entrez-fetch
dnakit skills install entrez-fetch --agent claude-code   # or --agent codex / cursor
```

## `entrez-link` {#entrez-link}

Find cross-database references between NCBI databases using Biopython Bio.Entrez (ELink). Use when navigating gene to protein/structure, sequence to publication, PubMed to GEO, BioProject to SRA runs, or discovering all link relationships for a record. Covers linkname semantics, cmd= variants, asymmetric link warnings, neighbor_history for &gt;200 input IDs, and per-database link tables.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/database-access/entrez-link)
- **License:** `MIT`
- **Runs with:** Python modules `Bio.Entrez`, `dnakit`
- **Also known as:** `bio-entrez-link`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info entrez-link
dnakit skills install entrez-link --agent claude-code   # or --agent codex / cursor
```

## `entrez-search` {#entrez-search}

Search NCBI databases using Biopython Bio.Entrez (ESearch, EInfo, EGQuery, ESpell). Use when finding records by keyword, building reproducible field-qualified queries, navigating the Entrez Query Translator, exploiting the history server for large result sets, handling retmax caps, or interpreting weekly index lag. Covers PubMed, Nucleotide, Protein, Gene, SRA, GEO, Assembly, Taxonomy, ClinVar, dbSNP.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/database-access/entrez-search)
- **License:** `MIT`
- **Runs with:** Python modules `Bio.Entrez`, `dnakit`
- **Also known as:** `bio-entrez-search`, `ncbi-edirect`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info entrez-search
dnakit skills install entrez-search --agent claude-code   # or --agent codex / cursor
```

## `eqtl-catalogue-region-fetch` {#eqtl-catalogue-region-fetch}

Fetch a region of cis-eQTL summary statistics from EBI eQTL Catalogue v7+
via tabix-on-FTP. Use when an agent needs eQTL beta / SE / p-value for
every variant in a window around a gene's TSS for one specific dataset
(study × tissue × quantification method). Input: dataset_id, chromosome,
start, end, optional molecular_trait_id. Output: harmonised TSV slice.

- **Source:** [Aviv Madar; ClawBio contributors @ `1fcecb7`](https://github.com/ClawBio/ClawBio/tree/1fcecb7e9531210baabe6d5cdbdd08de1ff3ab8d/skills/eqtl-catalogue-region-fetch)
- **License:** `MIT`
- **Runs with:** commands `python3`, `tabix`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info eqtl-catalogue-region-fetch
dnakit skills install eqtl-catalogue-region-fetch --agent claude-code   # or --agent codex / cursor
```

## `gene-database` {#gene-database}

NCBI Gene via E-utilities: curated records across 1M+ taxa. Official symbols, aliases, RefSeq IDs, summaries, coordinates, GO, interactions. Use for gene ID resolution and cross-species function queries. For sequences use Ensembl; for expression use geo-database.

- **Source:** [jaechang-hits; credited K-Dense-AI and Biomni contributors @ `fe505ca`](https://github.com/jaechang-hits/SciAgent-Skills/tree/fe505cae14d20b6c33be2e49666425be98f005bb/skills/genomics-bioinformatics/databases/gene-database)
- **License:** `CC-BY-4.0`
- **Also known as:** `alterlab-gene-db`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info gene-database
dnakit skills install gene-database --agent claude-code   # or --agent codex / cursor
```

## `geo-data` {#geo-data}

Query and download from NCBI Gene Expression Omnibus (GEO) and EMBL-EBI's BioStudies/ArrayExpress mirror. Use when finding expression datasets, navigating SuperSeries vs SubSeries, choosing between series-matrix (submitter-normalized) and raw supplementary files, downloading via GEOparse (Python) or GEOquery (R/Bioconductor), linking GEO to SRA for raw reads, or distinguishing GSE/GSM/GPL/GDS record types. Encodes the SuperSeries trap, the series-matrix normalization-trust caveat, GEOmetadb deprecation, ArrayExpress migration to BioStudies, and processed-vs-raw decision matrix.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/database-access/geo-data)
- **License:** `MIT`
- **Also known as:** `bio-geo-data`, `geo-database`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info geo-data
dnakit skills install geo-data --agent claude-code   # or --agent codex / cursor
```

## `gget-genomic-databases` {#gget-genomic-databases}

Unified CLI/Python interface to 20+ genomic databases. Gene lookups (Ensembl search/info/seq), BLAST/BLAT, AlphaFold, Enrichr enrichment, OpenTargets disease/drug, CELLxGENE single-cell, cBioPortal/COSMIC cancer, ARCHS4 expression. Spans genomics, proteomics, disease. For batch/advanced BLAST use biopython; for multi-DB Python SDK use bioservices.

- **Source:** [jaechang-hits; credited K-Dense-AI and Biomni contributors @ `fe505ca`](https://github.com/jaechang-hits/SciAgent-Skills/tree/fe505cae14d20b6c33be2e49666425be98f005bb/skills/genomics-bioinformatics/databases/gget-genomic-databases)
- **License:** `CC-BY-4.0`
- **Also known as:** `gget`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info gget-genomic-databases
dnakit skills install gget-genomic-databases --agent claude-code   # or --agent codex / cursor
```

## `gigwa-explore-instance` {#gigwa-explore-instance}

Get an overview of an entire Gigwa instance and flag anything that needs attention. Use this when a user asks what data is available on a Gigwa server, wants to enumerate databases/projects/runs, or has inherited an instance and needs a health check before analysis. Confirms connectivity, server version and authenticated user; lists all databases, projects and runs with their exact BrAPI variantSetDbIds and variant/call-set counts; and runs an instance-wide import-quality audit that ranks runs BROKEN/SUSPECT/OK so badly imported data is caught early. Trigger on 'what is on this Gigwa server', 'give me an overview of the instance', 'list the databases and runs', or 'which runs were imported badly'. Keywords: Gigwa, instance overview, server info, list content, variant sets, variantSetDbId, import audit, data catalog.

- **Source:** [Copyright 2026 Anestis Gkanogiannis; Gigwa-MCP contributors @ `d6ccdd1`](https://github.com/gkanogiannis/Gigwa-MCP/tree/d6ccdd1a5eb19b8a0cc224b96ad459868d9c4aaf/skills/gigwa-explore-instance)
- **License:** `Apache-2.0`
- **Runs with:** Python modules `gigwa_mcp`; commands `gigwa-mcp`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info gigwa-explore-instance
dnakit skills install gigwa-explore-instance --agent claude-code   # or --agent codex / cursor
```

## `gtdb` {#gtdb}

Use the Genome Taxonomy Database (GTDB) public releases — taxonomy tables, trees, metadata, and representative genomes. Use this skill when: (1) Choosing or citing a GTDB release, (2) Downloading bac120/ar53 taxonomy or reps, (3) Distinguishing GTDB data files from the GTDB-Tk runtime package, (4) Preparing reference sets for MAG comparison. Upstream: https://gtdb.ecogenomic.org/ and https://github.com/Ecogenomics/GTDBTk. For classifying your MAGs, use gtdbtk. For biome catalogues, use mgnify. For routing, use public-databases.

- **Source:** [Copyright (c) 2026 PKU-EMBL @ `e0e19a6`](https://github.com/PKU-EMBL/Metagenomics-Skills/tree/e0e19a600d9916736c0713e7659dec3782c51258/skills/gtdb)
- **License:** `MIT`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info gtdb
dnakit skills install gtdb --agent claude-code   # or --agent codex / cursor
```

## `gwas-catalog-region-fetch` {#gwas-catalog-region-fetch}

Fetch a region of GWAS summary statistics from the NHGRI-EBI GWAS Catalog
harmonised collection via tabix-on-FTP. Use when an agent needs GWAS beta /
SE / p-value for every variant in a window for one specific study (GCST
accession). Input: accession, chromosome, start, end. Output: harmonised
TSV slice in canonical format.

- **Source:** [Aviv Madar; ClawBio contributors @ `1fcecb7`](https://github.com/ClawBio/ClawBio/tree/1fcecb7e9531210baabe6d5cdbdd08de1ff3ab8d/skills/gwas-catalog-region-fetch)
- **License:** `MIT`
- **Runs with:** commands `python3`, `tabix`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info gwas-catalog-region-fetch
dnakit skills install gwas-catalog-region-fetch --agent claude-code   # or --agent codex / cursor
```

## `gwas-database` {#gwas-database}

NHGRI-EBI GWAS Catalog REST API for SNP-trait associations from published GWAS. Query studies, associations, variants, traits, genes, summary stats. Build PRS candidates, analyze pleiotropy, fetch stats for Manhattan plots. No auth.

- **Source:** [jaechang-hits; credited K-Dense-AI and Biomni contributors @ `fe505ca`](https://github.com/jaechang-hits/SciAgent-Skills/tree/fe505cae14d20b6c33be2e49666425be98f005bb/skills/genomics-bioinformatics/databases/gwas-database)
- **License:** `CC-BY-4.0`
- **Also known as:** `alterlab-gwas`, `gwas-catalog`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info gwas-database
dnakit skills install gwas-database --agent claude-code   # or --agent codex / cursor
```

## `gwas-fetch` {#gwas-fetch}

Fetch harmonised GWAS summary statistics from the GWAS Catalog by accession (GCST*). Resolves the study's metadata (genome build, FTP URL, sample size, ancestry) via the OKG when available, falling back to the GWAS Catalog REST API. Caches downloads to a configurable directory and writes a sidecar provenance manifest with SHA-256 + KG node IDs. Use when a user asks to download a GWAS Catalog study, fetch sumstats for a GCST accession, or pull sumstats for an analysis pipeline. Does NOT lift coordinates — chain into the `liftover` skill for that.

- **Source:** [Copyright 2026 Stephen Dorn; statgen-skills contributors @ `dc5f0c6`](https://github.com/svdorn/statgen-skills/tree/dc5f0c60d2995aebfe58b3f8a07b46f27562f4e0/gwas-fetch)
- **License:** `MIT`
- **Runs with:** commands `python3`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info gwas-fetch
dnakit skills install gwas-fetch --agent claude-code   # or --agent codex / cursor
```

## `jaspar-database` {#jaspar-database}

JASPAR 2024 TF binding profiles via REST API and pyJASPAR. Retrieve PFMs/PWMs by TF name, JASPAR ID, species, or structural class. Scan DNA for TFBS; browse by taxon (human, mouse) or TF family (bHLH, zinc finger). Use for motif enrichment input, TFBS scanning, and regulatory sequence analysis. For ChIP-seq peak motif discovery use homer-motif-analysis; for regulatory variant scoring use regulomedb-database.

- **Source:** [jaechang-hits; credited K-Dense-AI and Biomni contributors @ `fe505ca`](https://github.com/jaechang-hits/SciAgent-Skills/tree/fe505cae14d20b6c33be2e49666425be98f005bb/skills/genomics-bioinformatics/databases/jaspar-database)
- **License:** `CC-BY-4.0`
- **Also known as:** `jaspar-motifs`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info jaspar-database
dnakit skills install jaspar-database --agent claude-code   # or --agent codex / cursor
```

## `local-blast` {#local-blast}

Build local BLAST databases and run searches using NCBI BLAST+ command-line tools. Use when running &gt;50 queries, building custom databases with -parse_seqids and -taxid, downloading prebuilt NCBI databases via update_blastdb.pl, choosing -task variants (megablast/dc-megablast/blastn/blastn-short), tuning soft/hard masking, scaling threads, or extracting hits with blastdbcmd. Encodes BLAST v5 vs v4 database format, taxonomy filtering, makeblastdb pitfalls.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/database-access/local-blast)
- **License:** `MIT`
- **Runs with:** commands `BLAST+`
- **Also known as:** `bio-local-blast`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info local-blast
dnakit skills install local-blast --agent claude-code   # or --agent codex / cursor
```

## `mgnify` {#mgnify}

Use EMBL-EBI MGnify for analysed microbiome studies and biome genome catalogues (including UHGG / human-gut and other biomes). Use this skill when: (1) Browsing or citing MGnify analyses, (2) Downloading UHGG or other MGnify Genomes catalogues, (3) Searching genes/MAGs against catalogues, (4) Building comparison sets for MAG novelty. Upstream: https://www.ebi.ac.uk/metagenomics/ and https://github.com/EBI-Metagenomics. For raw reads use sra-ena. For taxonomy releases use gtdb. For routing use public-databases.

- **Source:** [Copyright (c) 2026 PKU-EMBL @ `e0e19a6`](https://github.com/PKU-EMBL/Metagenomics-Skills/tree/e0e19a600d9916736c0713e7659dec3782c51258/skills/mgnify)
- **License:** `MIT`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info mgnify
dnakit skills install mgnify --agent claude-code   # or --agent codex / cursor
```

## `model-organism-genetics` {#model-organism-genetics}

Cross-species genetic analysis using model organism databases (MGI mouse, ZFIN zebrafish, FlyBase fruit fly, WormBase worm, SGD yeast, RGD rat, GBIF taxonomy). Maps human genes to orthologs, retrieves phenotype/expression/functional data, assesses gene function conservation, and identifies the best animal model for studying a human gene or disease.

- **Source:** [ToolUniverse contributors, mims-harvard @ `7a0ceb2`](https://github.com/mims-harvard/ToolUniverse/tree/7a0ceb2f3a2ae28594e55ec27344e839df8bc9f7/skills/tooluniverse-model-organism-genetics)
- **License:** `Apache-2.0`
- **Also known as:** `tooluniverse-model-organism-genetics`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info model-organism-genetics
dnakit skills install model-organism-genetics --agent claude-code   # or --agent codex / cursor
```

## `monarch-database` {#monarch-database}

Monarch Initiative knowledge graph REST API for disease-gene-phenotype associations and cross-species orthology. MONDO disease-to-gene/phenotype, HP phenotype profiles, cross-species comparisons. Use for rare disease gene prioritization and phenotype-based candidate ranking. For GWAS use gwas-database; for clinical pathogenicity use clinvar-database.

- **Source:** [jaechang-hits; credited K-Dense-AI and Biomni contributors @ `fe505ca`](https://github.com/jaechang-hits/SciAgent-Skills/tree/fe505cae14d20b6c33be2e49666425be98f005bb/skills/genomics-bioinformatics/databases/monarch-database)
- **License:** `CC-BY-4.0`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info monarch-database
dnakit skills install monarch-database --agent claude-code   # or --agent codex / cursor
```

## `mouse-phenome-database` {#mouse-phenome-database}

Retrieve mouse phenotype data from the Jackson Laboratory Mouse Phenome Database (MPD) via its REST API. Browse 520+ projects, look up per-project measure metadata, pull strain-level means (raw or LS-mean adjusted) and per-animal values, find measures by MP/VT ontology terms, and resolve strain nomenclature or gene coordinates. Use for QTL support, cross-strain comparison, mouse model selection, and ontology-driven phenotype discovery. Use monarch-database for disease-gene-phenotype knowledge graphs; ensembl-database for mouse genome annotations.

- **Source:** [jaechang-hits; credited K-Dense-AI and Biomni contributors @ `fe505ca`](https://github.com/jaechang-hits/SciAgent-Skills/tree/fe505cae14d20b6c33be2e49666425be98f005bb/skills/genomics-bioinformatics/databases/mouse-phenome-database)
- **License:** `CC-BY-4.0`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info mouse-phenome-database
dnakit skills install mouse-phenome-database --agent claude-code   # or --agent codex / cursor
```

## `ncbi-datasets-cli` {#ncbi-datasets-cli}

Download genome assemblies, gene records, and ortholog data from NCBI using the modern Datasets v2 CLI (replaces assembly_summary.txt scraping and many EFetch workflows). Use when bulk-pulling genome assemblies, gene metadata across species, ortholog sets, or BLAST databases; when E-utilities are too slow for genome-scale work; or when automatic checksum verification, parallel download, and clean accession-driven retrieval are required. Encodes the JSON-lines output format, dataformat conversion, --dehydrated for cloud workflows, and when Datasets is/isn't the right tool.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/database-access/ncbi-datasets-cli)
- **License:** `MIT`
- **Also known as:** `bio-ncbi-datasets-cli`, `ncbi-datasets`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info ncbi-datasets-cli
dnakit skills install ncbi-datasets-cli --agent claude-code   # or --agent codex / cursor
```

## `okg-search` {#okg-search}

Query the statgen-analysis OKG for datasets matching a trait, EFO term, seed-dataset, or curated category, with controlled-vocab ancestry filtering (EUR/AFR/EAS/AMR/SAS/admixed/unknown), minimum-N power threshold, and sumstats-availability filter. Outputs a TSV ready to feed into `gwas-fetch GCST(X)` for batch chaining (e.g. genetic correlation with related immune traits). Auto-installs psycopg on first use and connects via the OKG_DSN env var; no MCP-server roundtrip overhead. Use when a user asks "what other GWAS exist for trait X", "what immune-trait studies are in the graph", "find sibling datasets like this one", or wants to bulk-pull a panel of GCSTs for rg / cross-trait analyses.

- **Source:** [Copyright 2026 Stephen Dorn; statgen-skills contributors @ `dc5f0c6`](https://github.com/svdorn/statgen-skills/tree/dc5f0c60d2995aebfe58b3f8a07b46f27562f4e0/okg-search)
- **License:** `MIT`
- **Runs with:** Python modules `psycopg`; commands `python3`; environment `OKG_DSN`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info okg-search
dnakit skills install okg-search --agent claude-code   # or --agent codex / cursor
```

## `ortholog-inference` {#ortholog-inference}

Pull pre-computed ortholog calls from public databases (OrthoDB, Ensembl Compara, OMA browser, eggNOG, PANTHER, KEGG Orthology, HomoloGene) via their REST APIs. Use when orthologs are already curated upstream, when the question is "what is the X ortholog of Y" rather than "how to infer orthology de novo", when batch-mapping gene IDs across species, or when comparing the resources for consensus calls. Encodes confidence-level semantics, 1:1 vs 1:many vs many:many, HomoloGene deprecation, and when to defect to de novo computation.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/database-access/ortholog-inference)
- **License:** `MIT`
- **Runs with:** Python modules `dnakit`, `requests`
- **Also known as:** `bio-ortholog-inference`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info ortholog-inference
dnakit skills install ortholog-inference --agent claude-code   # or --agent codex / cursor
```

## `pdb-database` {#pdb-database}

Use when you want to search for or download experimentally-determined 3D structures for biomolecules (proteins, nucleic acids, bound ligands). Supports searching by sequence similarity, structure similarity, chemical and other attributes. Also use to get metadata about biomolecular structure experiments.

- **Source:** [Copyright 2026 Google LLC @ `6883275`](https://github.com/google-deepmind/science-skills/tree/68832757cbbf941c620b71df5756cf6e5cc287b0/skills/pdb_database)
- **License:** `Apache-2.0 AND CC-BY-4.0`
- **Runs with:** Python modules `polite_http`; commands `uv`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info pdb-database
dnakit skills install pdb-database --agent claude-code   # or --agent codex / cursor
```

## `population-genetics-1000genomes` {#population-genetics-1000genomes}

Population genetics using the 1000 Genomes Project (IGSR) — superpopulation/population search, sample metadata, variant frequencies across AFR/AMR/EAS/EUR/SAS, ancestry-specific analyses. Use for ancestry comparison, population-aware allele frequency lookups, and 1000-Genomes-cohort-specific analyses (distinct from gnomAD which has different sample composition).

- **Source:** [ToolUniverse contributors, mims-harvard @ `7a0ceb2`](https://github.com/mims-harvard/ToolUniverse/tree/7a0ceb2f3a2ae28594e55ec27344e839df8bc9f7/skills/tooluniverse-population-genetics-1000genomes)
- **License:** `Apache-2.0`
- **Also known as:** `onekgpd`, `tooluniverse-population-genetics-1000genomes`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info population-genetics-1000genomes
dnakit skills install population-genetics-1000genomes --agent claude-code   # or --agent codex / cursor
```

## `progenomes` {#progenomes}

proGenomes4 — nearly 2 million consistently QC'd and annotated high-quality prokaryotic genomes with habitat metadata, pan-genomes, MGE/BGC annotations. Use when: (1) Building large microbial genome reference sets, (2) Comparative genomics / pan-genome analyses, (3) Selecting habitat-filtered genome cohorts, (4) Distinguishing proGenomes from GTDB reps or MGnify biome catalogues. Portal: https://progenomes.embl.de/. For taxonomy releases use gtdb. For biome MAG catalogues use mgnify. For routing use public-databases.

- **Source:** [Copyright (c) 2026 PKU-EMBL @ `e0e19a6`](https://github.com/PKU-EMBL/Metagenomics-Skills/tree/e0e19a600d9916736c0713e7659dec3782c51258/skills/progenomes)
- **License:** `MIT`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info progenomes
dnakit skills install progenomes --agent claude-code   # or --agent codex / cursor
```

## `regulomedb-database` {#regulomedb-database}

Query RegulomeDB v2 GET REST API to score variants for regulatory function and retrieve overlapping evidence (TF binding, histone marks, DNase peaks, footprints, motifs, eQTLs, chromatin state). Scores range 1a (strongest) to 7 (none). Use for GWAS hit prioritization, regulatory variant annotation, cis-regulatory discovery. Use clinvar-database for pathogenicity; gwas-database for trait associations.

- **Source:** [jaechang-hits; credited K-Dense-AI and Biomni contributors @ `fe505ca`](https://github.com/jaechang-hits/SciAgent-Skills/tree/fe505cae14d20b6c33be2e49666425be98f005bb/skills/genomics-bioinformatics/databases/regulomedb-database)
- **License:** `CC-BY-4.0`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info regulomedb-database
dnakit skills install regulomedb-database --agent claude-code   # or --agent codex / cursor
```

## `remap-database` {#remap-database}

Query ReMap 2022 TF ChIP-seq peak database via REST API and BED downloads. Retrieve TF peaks overlapping a region (chr:start-end), peaks near a gene, TFs by species, peaks filtered by biotype (promoter, enhancer), and BED files for a TF-cell type pair. Use for TF co-occupancy, regulatory annotation, and TF binding atlases. Use jaspar-database for PWM motifs; encode-database for ENCODE tracks.

- **Source:** [jaechang-hits; credited K-Dense-AI and Biomni contributors @ `fe505ca`](https://github.com/jaechang-hits/SciAgent-Skills/tree/fe505cae14d20b6c33be2e49666425be98f005bb/skills/genomics-bioinformatics/databases/remap-database)
- **License:** `CC-BY-4.0`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info remap-database
dnakit skills install remap-database --agent claude-code   # or --agent codex / cursor
```

## `sequence-retrieval` {#sequence-retrieval}

Retrieve DNA, RNA and protein sequences with DNAKit from NCBI and ENA. Use for accession/version lookup, organism and gene disambiguation, transcript isoforms, reference-versus-submission comparisons, and FASTA/GenBank/EMBL downloads with source records.

- **Source:** [ToolUniverse contributors, mims-harvard @ `7a0ceb2`](https://github.com/mims-harvard/ToolUniverse/tree/7a0ceb2f3a2ae28594e55ec27344e839df8bc9f7/skills/tooluniverse-sequence-retrieval)
- **License:** `Apache-2.0`
- **Also known as:** `public-db-lookup`, `tooluniverse-sequence-retrieval`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info sequence-retrieval
dnakit skills install sequence-retrieval --agent claude-code   # or --agent codex / cursor
```

## `sra-data` {#sra-data}

Download raw sequencing reads from NCBI SRA using sra-tools (prefetch, fasterq-dump, vdb-validate) or the ENA mirror. Use when pulling FASTQ for SRR/ERR/DRR accessions, deciding between SRA-direct, ENA mirror, or AWS/GCP cloud mirror (STRIDES), handling --include-technical for 10x and other single-cell records, validating with MD5/vdb-validate, navigating SRR/SRX/SRS/SRP/PRJNA hierarchy, or finding accessions via pysradb. Encodes SRA cloud-egress economics, the fasterq-dump uncompressed-scratch trap, and the --max-size default that silently truncates large prefetches.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/database-access/sra-data)
- **License:** `MIT`
- **Runs with:** commands `sra-tools`
- **Also known as:** `bio-sra-data`, `sra-ena`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info sra-data
dnakit skills install sra-data --agent claude-code   # or --agent codex / cursor
```

## `ucsc-genome-browser` {#ucsc-genome-browser}

Query UCSC Genome Browser REST API for DNA sequences, tracks, gene models, and conservation across 100+ assemblies. Retrieve sequence by region, list/fetch BED/bigWig tracks, chromosome sizes, RefSeq/GENCODE gene structures, PhyloP/PhastCons scores. Use for UCSC annotations; Ensembl REST API for Ensembl gene IDs and VEP variant annotation.

- **Source:** [jaechang-hits; credited K-Dense-AI and Biomni contributors @ `fe505ca`](https://github.com/jaechang-hits/SciAgent-Skills/tree/fe505cae14d20b6c33be2e49666425be98f005bb/skills/genomics-bioinformatics/databases/ucsc-genome-browser)
- **License:** `CC-BY-4.0`
- **Also known as:** `ucsc-browser`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info ucsc-genome-browser
dnakit skills install ucsc-genome-browser --agent claude-code   # or --agent codex / cursor
```

## `unibind-database` {#unibind-database}

Queries the UniBind database for experimentally validated transcription factor (TF) binding sites. Use when retrieving direct TF-DNA interaction datasets, downloading binding site coordinates (BED/FASTA) for local analysis, or listing available datasets by species, cell line, or TF name. Don't use to query specific intervals, locations, genes, motif models or expression data.

- **Source:** [Copyright 2026 Google LLC @ `28b8482`](https://github.com/google-deepmind/science-skills/tree/28b8482603a420708c8896f6fe5e06c276d9933d/skills/unibind_database)
- **License:** `Apache-2.0 AND CC-BY-4.0`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info unibind-database
dnakit skills install unibind-database --agent claude-code   # or --agent codex / cursor
```
