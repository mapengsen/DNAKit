# Alignment and alignment files

**20** Skills. A Skill is a set of instructions that guides an Agent through an analysis; the software, reference data and APIs it uses are separate dependencies. See [Skills](../skills.md) for searching, installing and using them; the Skills tab of the web demo shows each Skill in full.

| Skill | Source | Summary |
| --- | --- | --- |
| [`alignment-amplicon-clipping`](#alignment-amplicon-clipping) | bioskills | Trim PCR primers from aligned reads in amplicon-panel BAMs using samtools ampliconclip. Use when processing SARS-CoV-2 ARTIC, hereditary cancer panels, ctDNA h… |
| [`alignment-filtering`](#alignment-filtering) | bioskills | Filter alignments by flags, mapping quality, and regions using samtools view and pysam. Use when extracting specific reads, removing low-quality alignments, or… |
| [`alignment-indexing`](#alignment-indexing) | bioskills | Create and use BAI/CSI indices for BAM/CRAM files using samtools and pysam. Use when enabling random access to alignment files or fetching specific genomic reg… |
| [`alignment-io`](#alignment-io) | bioskills | Read, write, and convert multiple sequence alignment files using Biopython Bio.AlignIO. Supports Clustal, PHYLIP, Stockholm, FASTA, Nexus, and other alignment… |
| [`alignment-msa-parsing`](#alignment-msa-parsing) | bioskills | Parse and analyze multiple sequence alignments using Biopython. Extract sequences, identify conserved regions, analyze gaps, work with annotations, and manipul… |
| [`alignment-msa-statistics`](#alignment-msa-statistics) | bioskills | Calculate alignment statistics including sequence identity, conservation scores, substitution matrices, and similarity metrics. Use when comparing alignment qu… |
| [`alignment-multiple`](#alignment-multiple) | bioskills | Perform multiple sequence alignment using MAFFT, MUSCLE5, ClustalOmega, or T-Coffee. Guides tool and algorithm selection based on dataset size, sequence diverg… |
| [`alignment-pairwise`](#alignment-pairwise) | bioskills | Perform pairwise sequence alignment using Biopython Bio.Align.PairwiseAligner. Use when comparing two sequences, finding optimal alignments, scoring similarity… |
| [`alignment-sorting`](#alignment-sorting) | bioskills | Sort alignment files by coordinate or read name using samtools and pysam. Use when preparing BAM files for indexing, variant calling, or paired-end analysis. |
| [`alignment-trimming`](#alignment-trimming) | bioskills | Trim multiple sequence alignments using ClipKIT, trimAl, BMGE, Divvier, or HMMcleaner with mode selection guidance per downstream goal. Use when removing unrel… |
| [`alignment-validation`](#alignment-validation) | bioskills | Validate alignment quality with insert size distribution, proper pairing rates, GC bias, strand balance, and other post-alignment metrics. Use when verifying a… |
| [`bam-statistics`](#bam-statistics) | bioskills | Generate alignment statistics using samtools flagstat, stats, depth, coverage, and mosdepth. Use when assessing alignment quality, calculating coverage, or gen… |
| [`duplicate-handling`](#duplicate-handling) | bioskills | Mark and remove PCR/optical duplicates using samtools fixmate and markdup. Use when preparing alignments for variant calling or when duplicate reads would bias… |
| [`pileup-generation`](#pileup-generation) | bioskills | Generate pileup data for variant calling using samtools mpileup and pysam. Use when preparing data for variant calling, analyzing per-position read data, or ca… |
| [`pysam`](#pysam) | kdense | Python/HTSlib workflows for genomic files. Use when reading, querying, filtering, or writing SAM/BAM/CRAM, VCF/BCF, FASTA/FASTQ, or tabix data with pysam, incl… |
| [`read-alignment-bowtie2-alignment`](#read-alignment-bowtie2-alignment) | bioskills | Aligns DNA short reads to a reference with Bowtie2, choosing end-to-end (whole read must align) vs local (soft-clip read ends) mode and a sensitivity preset; t… |
| [`read-alignment-bwa-alignment`](#read-alignment-bwa-alignment) | bioskills | Aligns DNA short reads (paired- or single-end) to a reference genome with bwa-mem2, the maintained successor to BWA-MEM, for WGS/WES and germline/somatic varia… |
| [`reference-operations`](#reference-operations) | bioskills | Generate consensus sequences and manage reference files using samtools. Use when creating consensus from alignments, indexing references, or creating sequence… |
| [`sam-bam-basics`](#sam-bam-basics) | bioskills | View, convert, and understand SAM/BAM/CRAM alignment files using samtools and pysam. Use when inspecting alignments, converting between formats, or understandi… |
| [`samtools-bam-processing`](#samtools-bam-processing) | sciagent | CLI toolkit for SAM/BAM/CRAM: sort, index, convert, filter, QC alignments. Core commands: view, sort, index, flagstat, stats, depth, markdup, merge. Required b… |

## `alignment-amplicon-clipping` {#alignment-amplicon-clipping}

Trim PCR primers from aligned reads in amplicon-panel BAMs using samtools ampliconclip. Use when processing SARS-CoV-2 ARTIC, hereditary cancer panels, ctDNA hot-spot panels, or any amplicon assay where primer-derived bases would falsely confirm reference at primer footprints.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/alignment-files/alignment-amplicon-clipping)
- **License:** `MIT`
- **Runs with:** commands `samtools`
- **Also known as:** `bio-alignment-amplicon-clipping`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info alignment-amplicon-clipping
dnakit skills install alignment-amplicon-clipping --agent claude-code   # or --agent codex / cursor
```

## `alignment-filtering` {#alignment-filtering}

Filter alignments by flags, mapping quality, and regions using samtools view and pysam. Use when extracting specific reads, removing low-quality alignments, or subsetting to target regions.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/alignment-files/alignment-filtering)
- **License:** `MIT`
- **Runs with:** commands `samtools`
- **Also known as:** `bio-alignment-filtering`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info alignment-filtering
dnakit skills install alignment-filtering --agent claude-code   # or --agent codex / cursor
```

## `alignment-indexing` {#alignment-indexing}

Create and use BAI/CSI indices for BAM/CRAM files using samtools and pysam. Use when enabling random access to alignment files or fetching specific genomic regions.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/alignment-files/alignment-indexing)
- **License:** `MIT`
- **Runs with:** commands `samtools`
- **Also known as:** `bio-alignment-indexing`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info alignment-indexing
dnakit skills install alignment-indexing --agent claude-code   # or --agent codex / cursor
```

## `alignment-io` {#alignment-io}

Read, write, and convert multiple sequence alignment files using Biopython Bio.AlignIO. Supports Clustal, PHYLIP, Stockholm, FASTA, Nexus, and other alignment formats for phylogenetics and conservation analysis. Use when reading, writing, or converting alignment file formats.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/alignment/alignment-io)
- **License:** `MIT`
- **Runs with:** Python modules `Bio.AlignIO`, `dnakit`
- **Also known as:** `bio-alignment-io`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info alignment-io
dnakit skills install alignment-io --agent claude-code   # or --agent codex / cursor
```

## `alignment-msa-parsing` {#alignment-msa-parsing}

Parse and analyze multiple sequence alignments using Biopython. Extract sequences, identify conserved regions, analyze gaps, work with annotations, and manipulate alignment data for downstream analysis. Use when parsing or manipulating multiple sequence alignments.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/alignment/msa-parsing)
- **License:** `MIT`
- **Runs with:** Python modules `Bio.AlignIO`, `dnakit`
- **Also known as:** `bio-alignment-msa-parsing`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info alignment-msa-parsing
dnakit skills install alignment-msa-parsing --agent claude-code   # or --agent codex / cursor
```

## `alignment-msa-statistics` {#alignment-msa-statistics}

Calculate alignment statistics including sequence identity, conservation scores, substitution matrices, and similarity metrics. Use when comparing alignment quality, measuring sequence divergence, and analyzing evolutionary patterns.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/alignment/msa-statistics)
- **License:** `MIT`
- **Runs with:** Python modules `Bio.Align`, `dnakit`
- **Also known as:** `bio-alignment-msa-statistics`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info alignment-msa-statistics
dnakit skills install alignment-msa-statistics --agent claude-code   # or --agent codex / cursor
```

## `alignment-multiple` {#alignment-multiple}

Perform multiple sequence alignment using MAFFT, MUSCLE5, ClustalOmega, or T-Coffee. Guides tool and algorithm selection based on dataset size, sequence divergence, and downstream application. Use when aligning three or more homologous sequences for phylogenetics, conservation analysis, or evolutionary studies.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/alignment/multiple-alignment)
- **License:** `MIT`
- **Also known as:** `bio-alignment-multiple`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info alignment-multiple
dnakit skills install alignment-multiple --agent claude-code   # or --agent codex / cursor
```

## `alignment-pairwise` {#alignment-pairwise}

Perform pairwise sequence alignment using Biopython Bio.Align.PairwiseAligner. Use when comparing two sequences, finding optimal alignments, scoring similarity, and identifying local or global matches between DNA, RNA, or protein sequences.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/alignment/pairwise-alignment)
- **License:** `MIT`
- **Runs with:** Python modules `Bio.Align`, `dnakit`
- **Also known as:** `bio-alignment-pairwise`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info alignment-pairwise
dnakit skills install alignment-pairwise --agent claude-code   # or --agent codex / cursor
```

## `alignment-sorting` {#alignment-sorting}

Sort alignment files by coordinate or read name using samtools and pysam. Use when preparing BAM files for indexing, variant calling, or paired-end analysis.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/alignment-files/alignment-sorting)
- **License:** `MIT`
- **Runs with:** commands `samtools`
- **Also known as:** `bio-alignment-sorting`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info alignment-sorting
dnakit skills install alignment-sorting --agent claude-code   # or --agent codex / cursor
```

## `alignment-trimming` {#alignment-trimming}

Trim multiple sequence alignments using ClipKIT, trimAl, BMGE, Divvier, or HMMcleaner with mode selection guidance per downstream goal. Use when removing unreliable columns or contaminating residues before phylogenetic inference, HMM building, or selection analysis.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/alignment/alignment-trimming)
- **License:** `MIT`
- **Also known as:** `bio-alignment-trimming`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info alignment-trimming
dnakit skills install alignment-trimming --agent claude-code   # or --agent codex / cursor
```

## `alignment-validation` {#alignment-validation}

Validate alignment quality with insert size distribution, proper pairing rates, GC bias, strand balance, and other post-alignment metrics. Use when verifying alignment data quality before variant calling or quantification.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/alignment-files/alignment-validation)
- **License:** `MIT`
- **Also known as:** `bio-alignment-validation`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info alignment-validation
dnakit skills install alignment-validation --agent claude-code   # or --agent codex / cursor
```

## `bam-statistics` {#bam-statistics}

Generate alignment statistics using samtools flagstat, stats, depth, coverage, and mosdepth. Use when assessing alignment quality, calculating coverage, or generating QC reports.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/alignment-files/bam-statistics)
- **License:** `MIT`
- **Runs with:** commands `samtools`
- **Also known as:** `bio-bam-statistics`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info bam-statistics
dnakit skills install bam-statistics --agent claude-code   # or --agent codex / cursor
```

## `duplicate-handling` {#duplicate-handling}

Mark and remove PCR/optical duplicates using samtools fixmate and markdup. Use when preparing alignments for variant calling or when duplicate reads would bias analysis.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/alignment-files/duplicate-handling)
- **License:** `MIT`
- **Runs with:** commands `samtools`
- **Also known as:** `bio-duplicate-handling`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info duplicate-handling
dnakit skills install duplicate-handling --agent claude-code   # or --agent codex / cursor
```

## `pileup-generation` {#pileup-generation}

Generate pileup data for variant calling using samtools mpileup and pysam. Use when preparing data for variant calling, analyzing per-position read data, or calculating allele frequencies.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/alignment-files/pileup-generation)
- **License:** `MIT`
- **Runs with:** commands `samtools`
- **Also known as:** `bio-pileup-generation`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info pileup-generation
dnakit skills install pileup-generation --agent claude-code   # or --agent codex / cursor
```

## `pysam` {#pysam}

Python/HTSlib workflows for genomic files. Use when reading, querying, filtering, or writing SAM/BAM/CRAM, VCF/BCF, FASTA/FASTQ, or tabix data with pysam, including pileup, coverage, indexing, and CRAM references.

- **Source:** [K-Dense Inc.; K-Dense Inc. and contributors @ `c1ed16d`](https://github.com/K-Dense-AI/scientific-agent-skills/tree/c1ed16d97dd61ff50a3bd46dd353e4a55fd77f34/skills/pysam)
- **License:** `MIT`
- **Also known as:** `alterlab-pysam`, `pysam-genomic-files`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info pysam
dnakit skills install pysam --agent claude-code   # or --agent codex / cursor
```

## `read-alignment-bowtie2-alignment` {#read-alignment-bowtie2-alignment}

Aligns DNA short reads to a reference with Bowtie2, choosing end-to-end (whole read must align) vs local (soft-clip read ends) mode and a sensitivity preset; the de-facto aligner for ChIP-seq, ATAC-seq, and CUT&amp;RUN, where fragment-geometry flags (--no-mixed, --no-discordant, --dovetail, -X) and a tool-appropriate MAPQ filter feed the peak caller. Use when aligning ChIP/ATAC/CUT&amp;RUN reads, when read ends are adapter-contaminated and need soft-clipping, or when a tunable sensitivity/speed preset is wanted. DNA variant calling prefers bwa-alignment; RNA spliced alignment is star-alignment/hisat2-alignment; the QC gate and cross-tool MAPQ scale are alignment-files; peak calling is chip-seq/atac-seq; bisulfite uses methylation-analysis/bismark-alignment.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/read-alignment/bowtie2-alignment)
- **License:** `MIT`
- **Runs with:** commands `bowtie2`
- **Also known as:** `bio-read-alignment-bowtie2-alignment`, `bowtie2`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info read-alignment-bowtie2-alignment
dnakit skills install read-alignment-bowtie2-alignment --agent claude-code   # or --agent codex / cursor
```

## `read-alignment-bwa-alignment` {#read-alignment-bwa-alignment}

Aligns DNA short reads (paired- or single-end) to a reference genome with bwa-mem2, the maintained successor to BWA-MEM, for WGS/WES and germline/somatic variant-calling pipelines; covers index build, read-group injection, the collate/fixmate/sort/markdup ordering, soft-clipping for SV split reads, ALT/decoy-aware mapping on GRCh38, -K determinism, and streaming straight to a sorted BAM. Use when mapping DNA short reads to a reference for variant calling, coverage, ChIP/ATAC (alongside bowtie2-alignment), or SV detection. RNA-seq spliced alignment is star-alignment/hisat2-alignment; BAM sort/dedup/stats, the QC gate, and the cross-tool MAPQ scale are alignment-files; read trimming is read-qc; counting reads over features is rna-quantification.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/read-alignment/bwa-alignment)
- **License:** `MIT`
- **Runs with:** commands `bwa-mem2`
- **Also known as:** `alignment-and-mapping`, `bbmap`, `bio-read-alignment-bwa-alignment`, `bwa`, `bwa-mem2`, `bwa-mem2-dna-aligner`, `strobealign`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info read-alignment-bwa-alignment
dnakit skills install read-alignment-bwa-alignment --agent claude-code   # or --agent codex / cursor
```

## `reference-operations` {#reference-operations}

Generate consensus sequences and manage reference files using samtools. Use when creating consensus from alignments, indexing references, or creating sequence dictionaries.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/alignment-files/reference-operations)
- **License:** `MIT`
- **Runs with:** commands `samtools`
- **Also known as:** `bio-reference-operations`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info reference-operations
dnakit skills install reference-operations --agent claude-code   # or --agent codex / cursor
```

## `sam-bam-basics` {#sam-bam-basics}

View, convert, and understand SAM/BAM/CRAM alignment files using samtools and pysam. Use when inspecting alignments, converting between formats, or understanding alignment file structure.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/alignment-files/sam-bam-basics)
- **License:** `MIT`
- **Runs with:** commands `samtools`
- **Also known as:** `bio-sam-bam-basics`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info sam-bam-basics
dnakit skills install sam-bam-basics --agent claude-code   # or --agent codex / cursor
```

## `samtools-bam-processing` {#samtools-bam-processing}

CLI toolkit for SAM/BAM/CRAM: sort, index, convert, filter, QC alignments. Core commands: view, sort, index, flagstat, stats, depth, markdup, merge. Required between alignment and variant/peak calling. Use pysam for Python-native BAM access; deeptools for normalized coverage tracks.

- **Source:** [jaechang-hits; credited K-Dense-AI and Biomni contributors @ `fe505ca`](https://github.com/jaechang-hits/SciAgent-Skills/tree/fe505cae14d20b6c33be2e49666425be98f005bb/skills/genomics-bioinformatics/alignment/samtools-bam-processing)
- **License:** `CC-BY-4.0`
- **Also known as:** `samtools`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info samtools-bam-processing
dnakit skills install samtools-bam-processing --agent claude-code   # or --agent codex / cursor
```
