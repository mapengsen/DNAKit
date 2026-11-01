# DNA sequence processing

**24** Skills. A Skill is a set of instructions that guides an Agent through an analysis; the software, reference data and APIs it uses are separate dependencies. See [Skills](../skills.md) for searching, installing and using them; the Skills tab of the web demo shows each Skill in full.

| Skill | Source | Summary |
| --- | --- | --- |
| [`analyze-fasta`](#analyze-fasta) | clawbio | Analyze a single FASTA file (nucleotide or protein), compute sequence-level metrics (GC, ORFs, MW, pI, GRAVY, secondary-structure fractions) with Biopython, an… |
| [`batch-processing`](#batch-processing) | bioskills | Process many sequence files in batch (count, merge, split, convert, summarize) with memory-safe streaming and on-disk indexing using Biopython, pysam, or pyfas… |
| [`codon-usage`](#codon-usage) | bioskills | Analyze codon usage and calculate CAI (Codon Adaptation Index), RSCU, and Nc with Biopython, and produce naive max-CAI codon-optimized sequences. Use when scor… |
| [`compressed-files`](#compressed-files) | bioskills | Read, write, and index compressed sequence files (gzip, bzip2, xz, BGZF) with Biopython and bgzip/samtools. Use when working with .gz, .bz2, or .bgz sequence f… |
| [`data-visualization-sequence-logos`](#data-visualization-sequence-logos) | bioskills | Build sequence logos from aligned DNA, RNA, or protein motifs using ggseqlogo (R), Logomaker (Python), or WebLogo with explicit bits vs probability encoding, b… |
| [`erast`](#erast) | pku-metagenomics | ERAST — efficient retrieval-augmented homology search combining sequence language-model embeddings with a billion-scale vector database (protein and nucleotide… |
| [`fasta-database-curator`](#fasta-database-curator) | omics | Curate and validate FASTA or FAA databases. Use when standardizing headers, merging references, deduplicating sequences, converting GenBank files, or preparing… |
| [`fastq-quality`](#fastq-quality) | bioskills | Work with FASTQ quality scores using Biopython - access Phred scores, filter and trim by quality, compute per-position profiles, and convert between Sanger/Phr… |
| [`filter-sequences`](#filter-sequences) | bioskills | Filter and select sequences by criteria (length, ID, GC content, N content, motifs, patterns, description) using Biopython, streaming so large files never load… |
| [`format-conversion`](#format-conversion) | bioskills | Convert between sequence file formats (FASTA, FASTQ, GenBank, EMBL, Stockholm) and re-encode FASTQ quality offsets using Biopython Bio.SeqIO. Use when changing… |
| [`lexicmap`](#lexicmap) | pku-metagenomics | LexicMap — efficient nucleotide alignment of moderate-length queries (genes, plasmids, long reads; typically &gt;250 bp) against up to millions of prokaryotic gen… |
| [`motif-search`](#motif-search) | bioskills | Find sequence motifs, degenerate IUPAC patterns, and transcription-factor binding sites in DNA/RNA using Biopython and regex, including position weight matrix… |
| [`paired-end-fastq`](#paired-end-fastq) | bioskills | Handle paired-end FASTQ files (R1/R2) using Biopython while keeping mates synchronized. Use when working with Illumina paired reads, synchronizing pairs, filte… |
| [`phylign`](#phylign) | pku-metagenomics | Phylign — BLAST-like alignment of genes, plasmids, or reads against phylogeny-compressed microbial genome collections (e.g. pre-2019 661k bacteria) on ordinary… |
| [`read-sequences`](#read-sequences) | bioskills | Read biological sequence files (FASTA, FASTQ, GenBank, EMBL, ABI, SFF) with Biopython Bio.SeqIO, choosing between streaming, in-memory, and on-disk-indexed acc… |
| [`reverse-complement`](#reverse-complement) | bioskills | Generate reverse complements and complements of DNA/RNA sequences using Biopython, including IUPAC ambiguity codes, gapped alignments, and minus-strand feature… |
| [`scikit-bio`](#scikit-bio) | sciagent | Python library for biology: sequence manipulation (DNA/RNA/protein), pairwise/multiple alignment, phylogenetic trees (NJ, UPGMA), diversity (Shannon, Faith PD,… |
| [`seq-objects`](#seq-objects) | bioskills | Create and manipulate Seq, MutableSeq, and SeqRecord objects using Biopython. Use when creating sequences from strings, modifying sequence data in-place, build… |
| [`sequence-analysis`](#sequence-analysis) | tooluniverse | Biological sequence analysis — gene/protein sequence retrieval (NCBI, Ensembl, UniProt), nucleotide/protein search, ortholog discovery, and FASTQ QC + alignmen… |
| [`sequence-properties`](#sequence-properties) | bioskills | Calculate nucleotide and protein sequence properties (GC content, GC skew, molecular weight, melting temperature, isoelectric point, instability, hydropathy) w… |
| [`sequence-slicing`](#sequence-slicing) | bioskills | Slice, extract, and concatenate biological sequences and annotated records using Biopython. Use when extracting subsequences by position, splicing exons into a… |
| [`sequence-statistics`](#sequence-statistics) | bioskills | Calculate assembly and sequence statistics (N50/L50, auN, NG50/NGA50, length distribution, GC content with ambiguity handling, summary reports) using Biopython… |
| [`transcription-translation`](#transcription-translation) | bioskills | Transcribe DNA to RNA and translate to protein using Biopython, with NCBI codon-table selection, CDS validation, and six-frame ORF finding. Use when converting… |
| [`write-sequences`](#write-sequences) | bioskills | Write biological sequences to files (FASTA, FASTQ, GenBank, EMBL) using Biopython Bio.SeqIO. Use when saving sequences, creating new sequence files, or outputt… |

## `analyze-fasta` {#analyze-fasta}

Analyze a single FASTA file (nucleotide or protein), compute sequence-level metrics (GC, ORFs, MW, pI, GRAVY, secondary-structure fractions) with Biopython, and write a Markdown report plus structured JSON for downstream chaining.

- **Source:** [Santiago Rodriguez Salinas; ClawBio contributors @ `1fcecb7`](https://github.com/ClawBio/ClawBio/tree/1fcecb7e9531210baabe6d5cdbdd08de1ff3ab8d/skills/analyze-fasta)
- **License:** `MIT`
- **Runs with:** commands `python3`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info analyze-fasta
dnakit skills install analyze-fasta --agent claude-code   # or --agent codex / cursor
```

## `batch-processing` {#batch-processing}

Process many sequence files in batch (count, merge, split, convert, summarize) with memory-safe streaming and on-disk indexing using Biopython, pysam, or pyfastx. Use when iterating over a directory of FASTA/FASTQ files, merging or splitting datasets, building random access across many or huge files, or automating per-file operations without exhausting RAM.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/sequence-io/batch-processing)
- **License:** `MIT`
- **Runs with:** Python modules `Bio.SeqIO`, `dnakit`
- **Also known as:** `bio-batch-processing`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info batch-processing
dnakit skills install batch-processing --agent claude-code   # or --agent codex / cursor
```

## `codon-usage` {#codon-usage}

Analyze codon usage and calculate CAI (Codon Adaptation Index), RSCU, and Nc with Biopython, and produce naive max-CAI codon-optimized sequences. Use when scoring a gene's codon bias against a host, optimizing a CDS for heterologous expression, or studying synonymous codon selection.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/sequence-manipulation/codon-usage)
- **License:** `MIT`
- **Runs with:** Python modules `Bio.SeqUtils.CodonAdaptationIndex`, `dnakit`
- **Also known as:** `bio-codon-usage`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info codon-usage
dnakit skills install codon-usage --agent claude-code   # or --agent codex / cursor
```

## `compressed-files` {#compressed-files}

Read, write, and index compressed sequence files (gzip, bzip2, xz, BGZF) with Biopython and bgzip/samtools. Use when working with .gz, .bz2, or .bgz sequence files, when random access into a compressed FASTA/FASTQ is needed, or when SeqIO.index/faidx/tabix rejects a plain .gz. Covers the BGZF-vs-gzip seekability asymmetry, the 'rt'-not-'rb' handle trap, virtual offsets, and gzip-to-BGZF conversion.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/sequence-io/compressed-files)
- **License:** `MIT`
- **Also known as:** `bio-compressed-files`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info compressed-files
dnakit skills install compressed-files --agent claude-code   # or --agent codex / cursor
```

## `data-visualization-sequence-logos` {#data-visualization-sequence-logos}

Build sequence logos from aligned DNA, RNA, or protein motifs using ggseqlogo (R), Logomaker (Python), or WebLogo with explicit bits vs probability encoding, background-frequency correction, custom alphabets, and multi-logo stacking. Use when visualizing motif PWMs (TF binding, splice sites, CRISPR spacers), aligned-position composition, or comparing two motif sets.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/data-visualization/sequence-logos)
- **License:** `MIT`
- **Also known as:** `bio-data-visualization-sequence-logos`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info data-visualization-sequence-logos
dnakit skills install data-visualization-sequence-logos --agent claude-code   # or --agent codex / cursor
```

## `erast` {#erast}

ERAST — efficient retrieval-augmented homology search combining sequence language-model embeddings with a billion-scale vector database (protein and nucleotide). Use when: (1) Classical BLAST/MMseqs/Foldseek scale poorly, (2) You need vector-DB nearest-neighbor search plus reranking, (3) Accessing the public ERAST web DB or running the local encode→predict pipeline. Upstream: https://github.com/TencentAILabHealthcare/ERAST. Web DB: https://ai4s.tencent.com/erast. Route via homology-search.

- **Source:** [Copyright (c) 2026 PKU-EMBL @ `e0e19a6`](https://github.com/PKU-EMBL/Metagenomics-Skills/tree/e0e19a600d9916736c0713e7659dec3782c51258/skills/erast)
- **License:** `MIT`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info erast
dnakit skills install erast --agent claude-code   # or --agent codex / cursor
```

## `fasta-database-curator` {#fasta-database-curator}

Curate and validate FASTA or FAA databases. Use when standardizing headers, merging references, deduplicating sequences, converting GenBank files, or preparing BLAST, MMseqs2, and HMM inputs.

- **Source:** [Copyright (c) 2026 Frederik Schulz @ `7a4eb5b`](https://github.com/fmschulz/omics-skills/tree/7a4eb5bace3f874c3be2ec697e9d2e695c374637/skills/bio-fasta-database-curator)
- **License:** `MIT`
- **Also known as:** `bio-fasta-database-curator`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info fasta-database-curator
dnakit skills install fasta-database-curator --agent claude-code   # or --agent codex / cursor
```

## `fastq-quality` {#fastq-quality}

Work with FASTQ quality scores using Biopython - access Phred scores, filter and trim by quality, compute per-position profiles, and convert between Sanger/Phred+33, Solexa, and Illumina/Phred+64 encodings. Use when analyzing read quality, filtering or trimming low-quality bases, generating quality reports, or deciding which FASTQ quality encoding a file uses before parsing.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/sequence-io/fastq-quality)
- **License:** `MIT`
- **Runs with:** Python modules `Bio.SeqIO`, `dnakit`
- **Also known as:** `bio-fastq-quality`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info fastq-quality
dnakit skills install fastq-quality --agent claude-code   # or --agent codex / cursor
```

## `filter-sequences` {#filter-sequences}

Filter and select sequences by criteria (length, ID, GC content, N content, motifs, patterns, description) using Biopython, streaming so large files never load into RAM. Use when subsetting a FASTA/FASTQ file, removing unwanted or low-quality records, or selecting records by specific criteria. Use the paired-end-fastq skill instead whenever the input is paired R1/R2 reads.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/sequence-io/filter-sequences)
- **License:** `MIT`
- **Runs with:** Python modules `Bio.SeqIO`, `dnakit`
- **Also known as:** `bio-filter-sequences`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info filter-sequences
dnakit skills install filter-sequences --agent claude-code   # or --agent codex / cursor
```

## `format-conversion` {#format-conversion}

Convert between sequence file formats (FASTA, FASTQ, GenBank, EMBL, Stockholm) and re-encode FASTQ quality offsets using Biopython Bio.SeqIO. Use when changing a file format for a downstream tool, fixing FASTQ quality encoding (Phred+33 vs Phred+64 vs Solexa), or when a conversion risks silently dropping annotations or quality scores.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/sequence-io/format-conversion)
- **License:** `MIT`
- **Runs with:** Python modules `Bio.SeqIO`, `dnakit`
- **Also known as:** `bio-format-conversion`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info format-conversion
dnakit skills install format-conversion --agent claude-code   # or --agent codex / cursor
```

## `lexicmap` {#lexicmap}

LexicMap — efficient nucleotide alignment of moderate-length queries (genes, plasmids, long reads; typically &gt;250 bp) against up to millions of prokaryotic genomes. Use when: (1) Searching genes/plasmids vs huge genome databases, (2) Scaling beyond BLAST/minimap2 for AllTheBacteria / GTDB-scale refs, (3) Epidemiology/ecology queries needing fast genome hits. Upstream: https://github.com/shenwei356/LexicMap. Not a substitute for short-read→contig coverage mappers (strobealign/bowtie2) used in binning.

- **Source:** [Copyright (c) 2026 PKU-EMBL @ `e0e19a6`](https://github.com/PKU-EMBL/Metagenomics-Skills/tree/e0e19a600d9916736c0713e7659dec3782c51258/skills/lexicmap)
- **License:** `MIT`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info lexicmap
dnakit skills install lexicmap --agent claude-code   # or --agent codex / cursor
```

## `motif-search` {#motif-search}

Find sequence motifs, degenerate IUPAC patterns, and transcription-factor binding sites in DNA/RNA using Biopython and regex, including position weight matrix (PWM/PSSM) scoring. Use when locating regulatory elements, counting overlapping motif occurrences, scanning for binding-site matches above a significance threshold, or reading motif matrices from JASPAR/MEME/TRANSFAC files. For restriction enzyme sites, use restriction-analysis/restriction-sites.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/sequence-manipulation/motif-search)
- **License:** `MIT`
- **Runs with:** Python modules `Bio.motifs`, `dnakit`
- **Also known as:** `bio-motif-search`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info motif-search
dnakit skills install motif-search --agent claude-code   # or --agent codex / cursor
```

## `paired-end-fastq` {#paired-end-fastq}

Handle paired-end FASTQ files (R1/R2) using Biopython while keeping mates synchronized. Use when working with Illumina paired reads, synchronizing pairs, filtering both mates together with orphan routing, interleaving/deinterleaving, or matching mates by read name.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/sequence-io/paired-end-fastq)
- **License:** `MIT`
- **Runs with:** Python modules `Bio.SeqIO`, `dnakit`
- **Also known as:** `bio-paired-end-fastq`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info paired-end-fastq
dnakit skills install paired-end-fastq --agent claude-code   # or --agent codex / cursor
```

## `phylign` {#phylign}

Phylign — BLAST-like alignment of genes, plasmids, or reads against phylogeny-compressed microbial genome collections (e.g. pre-2019 661k bacteria) on ordinary desktops. Uses phylogenetic compression (MiniPhy for building compressed assemblies/COBS indexes; Phylign for search via COBS + Minimap2). Use when searching huge genome warehouses without a cluster. Upstream: https://github.com/karel-brinda/phylign. Project: https://brinda.eu/mof. Route via homology-search; compare lexicmap.

- **Source:** [Copyright (c) 2026 PKU-EMBL @ `e0e19a6`](https://github.com/PKU-EMBL/Metagenomics-Skills/tree/e0e19a600d9916736c0713e7659dec3782c51258/skills/phylign)
- **License:** `MIT`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info phylign
dnakit skills install phylign --agent claude-code   # or --agent codex / cursor
```

## `read-sequences` {#read-sequences}

Read biological sequence files (FASTA, FASTQ, GenBank, EMBL, ABI, SFF) with Biopython Bio.SeqIO, choosing between streaming, in-memory, and on-disk-indexed access. Use when parsing sequence files, iterating multi-record files, randomly accessing records by ID in large files, or maximizing parse throughput.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/sequence-io/read-sequences)
- **License:** `MIT`
- **Runs with:** Python modules `Bio.SeqIO`, `dnakit`
- **Also known as:** `bio-read-sequences`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info read-sequences
dnakit skills install read-sequences --agent claude-code   # or --agent codex / cursor
```

## `reverse-complement` {#reverse-complement}

Generate reverse complements and complements of DNA/RNA sequences using Biopython, including IUPAC ambiguity codes, gapped alignments, and minus-strand features. Use when working with the opposite strand, building reverse primers, normalizing strand orientation before alignment, or extracting a coding sequence from a minus-strand feature.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/sequence-manipulation/reverse-complement)
- **License:** `MIT`
- **Runs with:** Python modules `Bio.Seq`, `dnakit`
- **Also known as:** `bio-reverse-complement`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info reverse-complement
dnakit skills install reverse-complement --agent claude-code   # or --agent codex / cursor
```

## `scikit-bio` {#scikit-bio}

Python library for biology: sequence manipulation (DNA/RNA/protein), pairwise/multiple alignment, phylogenetic trees (NJ, UPGMA), diversity (Shannon, Faith PD, Bray-Curtis, UniFrac), ordination (PCoA, CCA, RDA), stats (PERMANOVA, ANOSIM, Mantel), file I/O (FASTA, FASTQ, Newick, BIOM). Use for microbiome, community ecology, or phylogenetics.

- **Source:** [jaechang-hits; credited K-Dense-AI and Biomni contributors @ `fe505ca`](https://github.com/jaechang-hits/SciAgent-Skills/tree/fe505cae14d20b6c33be2e49666425be98f005bb/skills/genomics-bioinformatics/scikit-bio)
- **License:** `CC-BY-4.0`
- **Also known as:** `alterlab-scikit-bio`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info scikit-bio
dnakit skills install scikit-bio --agent claude-code   # or --agent codex / cursor
```

## `seq-objects` {#seq-objects}

Create and manipulate Seq, MutableSeq, and SeqRecord objects using Biopython. Use when creating sequences from strings, modifying sequence data in-place, building annotated records for file output, or debugging post-1.78 Bio.Alphabet and immutability errors.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/sequence-manipulation/seq-objects)
- **License:** `MIT`
- **Runs with:** Python modules `Bio.Seq`, `dnakit`
- **Also known as:** `bio-seq-objects`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info seq-objects
dnakit skills install seq-objects --agent claude-code   # or --agent codex / cursor
```

## `sequence-analysis` {#sequence-analysis}

Biological sequence analysis — gene/protein sequence retrieval (NCBI, Ensembl, UniProt), nucleotide/protein search, ortholog discovery, and FASTQ QC + alignment workflows (Trimmomatic, BWA, samtools, coverage depth). Use for sequence retrieval, sequence comparison, FASTQ QC analysis, and read alignment pre-processing.

- **Source:** [ToolUniverse contributors, mims-harvard @ `7a0ceb2`](https://github.com/mims-harvard/ToolUniverse/tree/7a0ceb2f3a2ae28594e55ec27344e839df8bc9f7/skills/tooluniverse-sequence-analysis)
- **License:** `Apache-2.0`
- **Also known as:** `biopython`, `biopython-molecular-biology`, `biopython-sequence-analysis`, `tooluniverse-sequence-analysis`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info sequence-analysis
dnakit skills install sequence-analysis --agent claude-code   # or --agent codex / cursor
```

## `sequence-properties` {#sequence-properties}

Calculate nucleotide and protein sequence properties (GC content, GC skew, molecular weight, melting temperature, isoelectric point, instability, hydropathy) with Biopython. Use when analyzing sequence composition, computing primer Tm, estimating DNA or protein mass, or profiling protein biophysical properties.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/sequence-manipulation/sequence-properties)
- **License:** `MIT`
- **Runs with:** Python modules `Bio.SeqUtils`, `dnakit`
- **Also known as:** `bio-sequence-properties`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info sequence-properties
dnakit skills install sequence-properties --agent claude-code   # or --agent codex / cursor
```

## `sequence-slicing` {#sequence-slicing}

Slice, extract, and concatenate biological sequences and annotated records using Biopython. Use when extracting subsequences by position, splicing exons into a transcript, joining sequences, or carrying a sub-region of an annotated record (with quality scores and features) into a new record.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/sequence-manipulation/sequence-slicing)
- **License:** `MIT`
- **Runs with:** Python modules `Bio.Seq`, `dnakit`
- **Also known as:** `bio-sequence-slicing`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info sequence-slicing
dnakit skills install sequence-slicing --agent claude-code   # or --agent codex / cursor
```

## `sequence-statistics` {#sequence-statistics}

Calculate assembly and sequence statistics (N50/L50, auN, NG50/NGA50, length distribution, GC content with ambiguity handling, summary reports) using Biopython. Use when analyzing sequence datasets, generating QC reports, or comparing genome assemblies.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/sequence-io/sequence-statistics)
- **License:** `MIT`
- **Runs with:** Python modules `Bio.SeqIO`, `dnakit`
- **Also known as:** `bio-sequence-statistics`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info sequence-statistics
dnakit skills install sequence-statistics --agent claude-code   # or --agent codex / cursor
```

## `transcription-translation` {#transcription-translation}

Transcribe DNA to RNA and translate to protein using Biopython, with NCBI codon-table selection, CDS validation, and six-frame ORF finding. Use when converting a CDS or ORF to its amino-acid sequence, selecting a non-standard (mitochondrial, bacterial, ciliate) genetic code, validating a coding sequence, or scanning all reading frames.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/sequence-manipulation/transcription-translation)
- **License:** `MIT`
- **Runs with:** Python modules `Bio.Seq`, `dnakit`
- **Also known as:** `bio-transcription-translation`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info transcription-translation
dnakit skills install transcription-translation --agent claude-code   # or --agent codex / cursor
```

## `write-sequences` {#write-sequences}

Write biological sequences to files (FASTA, FASTQ, GenBank, EMBL) using Biopython Bio.SeqIO. Use when saving sequences, creating new sequence files, or outputting modified records.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/sequence-io/write-sequences)
- **License:** `MIT`
- **Runs with:** Python modules `Bio.SeqIO`, `dnakit`
- **Also known as:** `bio-write-sequences`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info write-sequences
dnakit skills install write-sequences --agent claude-code   # or --agent codex / cursor
```
