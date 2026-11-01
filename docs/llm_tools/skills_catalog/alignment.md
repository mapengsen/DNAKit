# 比对与比对文件

本类共 **20** 个 Skills。每个 Skill 是一份让 Agent 按步骤完成分析的说明，分析所需的软件、参考数据和 API 需要另外准备。搜索、安装和使用方法见 [Skills](../skills.md)；网页版的 Skills 标签可以直接阅读每个 Skill 的完整内容。

| Skill | 来源 | 简介 |
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

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/alignment-files/alignment-amplicon-clipping)
- **许可证：**`MIT`
- **运行依赖：**命令 `samtools`
- **别名：**`bio-alignment-amplicon-clipping`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info alignment-amplicon-clipping
dnakit skills install alignment-amplicon-clipping --agent claude-code   # 或 --agent codex / cursor
```

## `alignment-filtering` {#alignment-filtering}

Filter alignments by flags, mapping quality, and regions using samtools view and pysam. Use when extracting specific reads, removing low-quality alignments, or subsetting to target regions.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/alignment-files/alignment-filtering)
- **许可证：**`MIT`
- **运行依赖：**命令 `samtools`
- **别名：**`bio-alignment-filtering`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info alignment-filtering
dnakit skills install alignment-filtering --agent claude-code   # 或 --agent codex / cursor
```

## `alignment-indexing` {#alignment-indexing}

Create and use BAI/CSI indices for BAM/CRAM files using samtools and pysam. Use when enabling random access to alignment files or fetching specific genomic regions.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/alignment-files/alignment-indexing)
- **许可证：**`MIT`
- **运行依赖：**命令 `samtools`
- **别名：**`bio-alignment-indexing`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info alignment-indexing
dnakit skills install alignment-indexing --agent claude-code   # 或 --agent codex / cursor
```

## `alignment-io` {#alignment-io}

Read, write, and convert multiple sequence alignment files using Biopython Bio.AlignIO. Supports Clustal, PHYLIP, Stockholm, FASTA, Nexus, and other alignment formats for phylogenetics and conservation analysis. Use when reading, writing, or converting alignment file formats.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/alignment/alignment-io)
- **许可证：**`MIT`
- **运行依赖：**Python 模块 `Bio.AlignIO`, `dnakit`
- **别名：**`bio-alignment-io`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info alignment-io
dnakit skills install alignment-io --agent claude-code   # 或 --agent codex / cursor
```

## `alignment-msa-parsing` {#alignment-msa-parsing}

Parse and analyze multiple sequence alignments using Biopython. Extract sequences, identify conserved regions, analyze gaps, work with annotations, and manipulate alignment data for downstream analysis. Use when parsing or manipulating multiple sequence alignments.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/alignment/msa-parsing)
- **许可证：**`MIT`
- **运行依赖：**Python 模块 `Bio.AlignIO`, `dnakit`
- **别名：**`bio-alignment-msa-parsing`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info alignment-msa-parsing
dnakit skills install alignment-msa-parsing --agent claude-code   # 或 --agent codex / cursor
```

## `alignment-msa-statistics` {#alignment-msa-statistics}

Calculate alignment statistics including sequence identity, conservation scores, substitution matrices, and similarity metrics. Use when comparing alignment quality, measuring sequence divergence, and analyzing evolutionary patterns.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/alignment/msa-statistics)
- **许可证：**`MIT`
- **运行依赖：**Python 模块 `Bio.Align`, `dnakit`
- **别名：**`bio-alignment-msa-statistics`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info alignment-msa-statistics
dnakit skills install alignment-msa-statistics --agent claude-code   # 或 --agent codex / cursor
```

## `alignment-multiple` {#alignment-multiple}

Perform multiple sequence alignment using MAFFT, MUSCLE5, ClustalOmega, or T-Coffee. Guides tool and algorithm selection based on dataset size, sequence divergence, and downstream application. Use when aligning three or more homologous sequences for phylogenetics, conservation analysis, or evolutionary studies.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/alignment/multiple-alignment)
- **许可证：**`MIT`
- **别名：**`bio-alignment-multiple`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info alignment-multiple
dnakit skills install alignment-multiple --agent claude-code   # 或 --agent codex / cursor
```

## `alignment-pairwise` {#alignment-pairwise}

Perform pairwise sequence alignment using Biopython Bio.Align.PairwiseAligner. Use when comparing two sequences, finding optimal alignments, scoring similarity, and identifying local or global matches between DNA, RNA, or protein sequences.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/alignment/pairwise-alignment)
- **许可证：**`MIT`
- **运行依赖：**Python 模块 `Bio.Align`, `dnakit`
- **别名：**`bio-alignment-pairwise`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info alignment-pairwise
dnakit skills install alignment-pairwise --agent claude-code   # 或 --agent codex / cursor
```

## `alignment-sorting` {#alignment-sorting}

Sort alignment files by coordinate or read name using samtools and pysam. Use when preparing BAM files for indexing, variant calling, or paired-end analysis.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/alignment-files/alignment-sorting)
- **许可证：**`MIT`
- **运行依赖：**命令 `samtools`
- **别名：**`bio-alignment-sorting`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info alignment-sorting
dnakit skills install alignment-sorting --agent claude-code   # 或 --agent codex / cursor
```

## `alignment-trimming` {#alignment-trimming}

Trim multiple sequence alignments using ClipKIT, trimAl, BMGE, Divvier, or HMMcleaner with mode selection guidance per downstream goal. Use when removing unreliable columns or contaminating residues before phylogenetic inference, HMM building, or selection analysis.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/alignment/alignment-trimming)
- **许可证：**`MIT`
- **别名：**`bio-alignment-trimming`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info alignment-trimming
dnakit skills install alignment-trimming --agent claude-code   # 或 --agent codex / cursor
```

## `alignment-validation` {#alignment-validation}

Validate alignment quality with insert size distribution, proper pairing rates, GC bias, strand balance, and other post-alignment metrics. Use when verifying alignment data quality before variant calling or quantification.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/alignment-files/alignment-validation)
- **许可证：**`MIT`
- **别名：**`bio-alignment-validation`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info alignment-validation
dnakit skills install alignment-validation --agent claude-code   # 或 --agent codex / cursor
```

## `bam-statistics` {#bam-statistics}

Generate alignment statistics using samtools flagstat, stats, depth, coverage, and mosdepth. Use when assessing alignment quality, calculating coverage, or generating QC reports.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/alignment-files/bam-statistics)
- **许可证：**`MIT`
- **运行依赖：**命令 `samtools`
- **别名：**`bio-bam-statistics`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info bam-statistics
dnakit skills install bam-statistics --agent claude-code   # 或 --agent codex / cursor
```

## `duplicate-handling` {#duplicate-handling}

Mark and remove PCR/optical duplicates using samtools fixmate and markdup. Use when preparing alignments for variant calling or when duplicate reads would bias analysis.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/alignment-files/duplicate-handling)
- **许可证：**`MIT`
- **运行依赖：**命令 `samtools`
- **别名：**`bio-duplicate-handling`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info duplicate-handling
dnakit skills install duplicate-handling --agent claude-code   # 或 --agent codex / cursor
```

## `pileup-generation` {#pileup-generation}

Generate pileup data for variant calling using samtools mpileup and pysam. Use when preparing data for variant calling, analyzing per-position read data, or calculating allele frequencies.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/alignment-files/pileup-generation)
- **许可证：**`MIT`
- **运行依赖：**命令 `samtools`
- **别名：**`bio-pileup-generation`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info pileup-generation
dnakit skills install pileup-generation --agent claude-code   # 或 --agent codex / cursor
```

## `pysam` {#pysam}

Python/HTSlib workflows for genomic files. Use when reading, querying, filtering, or writing SAM/BAM/CRAM, VCF/BCF, FASTA/FASTQ, or tabix data with pysam, including pileup, coverage, indexing, and CRAM references.

- **来源：**[K-Dense Inc.; K-Dense Inc. and contributors @ `c1ed16d`](https://github.com/K-Dense-AI/scientific-agent-skills/tree/c1ed16d97dd61ff50a3bd46dd353e4a55fd77f34/skills/pysam)
- **许可证：**`MIT`
- **别名：**`alterlab-pysam`, `pysam-genomic-files`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info pysam
dnakit skills install pysam --agent claude-code   # 或 --agent codex / cursor
```

## `read-alignment-bowtie2-alignment` {#read-alignment-bowtie2-alignment}

Aligns DNA short reads to a reference with Bowtie2, choosing end-to-end (whole read must align) vs local (soft-clip read ends) mode and a sensitivity preset; the de-facto aligner for ChIP-seq, ATAC-seq, and CUT&amp;RUN, where fragment-geometry flags (--no-mixed, --no-discordant, --dovetail, -X) and a tool-appropriate MAPQ filter feed the peak caller. Use when aligning ChIP/ATAC/CUT&amp;RUN reads, when read ends are adapter-contaminated and need soft-clipping, or when a tunable sensitivity/speed preset is wanted. DNA variant calling prefers bwa-alignment; RNA spliced alignment is star-alignment/hisat2-alignment; the QC gate and cross-tool MAPQ scale are alignment-files; peak calling is chip-seq/atac-seq; bisulfite uses methylation-analysis/bismark-alignment.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/read-alignment/bowtie2-alignment)
- **许可证：**`MIT`
- **运行依赖：**命令 `bowtie2`
- **别名：**`bio-read-alignment-bowtie2-alignment`, `bowtie2`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info read-alignment-bowtie2-alignment
dnakit skills install read-alignment-bowtie2-alignment --agent claude-code   # 或 --agent codex / cursor
```

## `read-alignment-bwa-alignment` {#read-alignment-bwa-alignment}

Aligns DNA short reads (paired- or single-end) to a reference genome with bwa-mem2, the maintained successor to BWA-MEM, for WGS/WES and germline/somatic variant-calling pipelines; covers index build, read-group injection, the collate/fixmate/sort/markdup ordering, soft-clipping for SV split reads, ALT/decoy-aware mapping on GRCh38, -K determinism, and streaming straight to a sorted BAM. Use when mapping DNA short reads to a reference for variant calling, coverage, ChIP/ATAC (alongside bowtie2-alignment), or SV detection. RNA-seq spliced alignment is star-alignment/hisat2-alignment; BAM sort/dedup/stats, the QC gate, and the cross-tool MAPQ scale are alignment-files; read trimming is read-qc; counting reads over features is rna-quantification.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/read-alignment/bwa-alignment)
- **许可证：**`MIT`
- **运行依赖：**命令 `bwa-mem2`
- **别名：**`alignment-and-mapping`, `bbmap`, `bio-read-alignment-bwa-alignment`, `bwa`, `bwa-mem2`, `bwa-mem2-dna-aligner`, `strobealign`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info read-alignment-bwa-alignment
dnakit skills install read-alignment-bwa-alignment --agent claude-code   # 或 --agent codex / cursor
```

## `reference-operations` {#reference-operations}

Generate consensus sequences and manage reference files using samtools. Use when creating consensus from alignments, indexing references, or creating sequence dictionaries.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/alignment-files/reference-operations)
- **许可证：**`MIT`
- **运行依赖：**命令 `samtools`
- **别名：**`bio-reference-operations`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info reference-operations
dnakit skills install reference-operations --agent claude-code   # 或 --agent codex / cursor
```

## `sam-bam-basics` {#sam-bam-basics}

View, convert, and understand SAM/BAM/CRAM alignment files using samtools and pysam. Use when inspecting alignments, converting between formats, or understanding alignment file structure.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/alignment-files/sam-bam-basics)
- **许可证：**`MIT`
- **运行依赖：**命令 `samtools`
- **别名：**`bio-sam-bam-basics`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info sam-bam-basics
dnakit skills install sam-bam-basics --agent claude-code   # 或 --agent codex / cursor
```

## `samtools-bam-processing` {#samtools-bam-processing}

CLI toolkit for SAM/BAM/CRAM: sort, index, convert, filter, QC alignments. Core commands: view, sort, index, flagstat, stats, depth, markdup, merge. Required between alignment and variant/peak calling. Use pysam for Python-native BAM access; deeptools for normalized coverage tracks.

- **来源：**[jaechang-hits; credited K-Dense-AI and Biomni contributors @ `fe505ca`](https://github.com/jaechang-hits/SciAgent-Skills/tree/fe505cae14d20b6c33be2e49666425be98f005bb/skills/genomics-bioinformatics/alignment/samtools-bam-processing)
- **许可证：**`CC-BY-4.0`
- **别名：**`samtools`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info samtools-bam-processing
dnakit skills install samtools-bam-processing --agent claude-code   # 或 --agent codex / cursor
```
