# Sequencing data analysis

## 1) Summarize BAM and create a regional consensus {#bio-14}

- **Purpose:** Inspect alignment counts, per-base depth and majority consensus.
- **Method:** Summarize BAM flags and count quality-filtered A/C/G/T bases; ties or insufficient depth yield N.
- **API:** `dnakit.bioinformatics.reads.alignment_summary(path)`; `dnakit.bioinformatics.reads.coverage_consensus(path, region, min_base_quality=20, min_mapping_quality=20, min_depth=1)`
- **Input:** Sorted, indexed `aligned.bam`; requires `bio-reads`. Consensus stays in reference-forward coordinates and omits indels.
- **Output:** `alignment_summary` returns alignment counts; `coverage_consensus` returns consensus sequence and coverage information.
- **Example code:**

```python
from dnakit.bioinformatics.reads import Region, alignment_summary, coverage_consensus

summary = alignment_summary("aligned.bam")
result = coverage_consensus("aligned.bam", Region("chr1", 100, 180))
print(summary)
print(result["consensus"])
```

- **Example output:**

`summary` is a count dictionary; `result` contains the region, A/C/G/T counts, depth and consensus sequence. Requires an actual BAM; this example was not executed.

Dependencies and details: [FAQ](../../faq/bioinformatics.md#reads).

## 2) Run a single-sample DNA-seq workflow {#bio-20}

- **Purpose:** Process reads into QC reports, BAM and VCF.
- **Method:** Chain FastQC → BWA → samtools → bcftools and retain each stage result.
- **API:** `dnakit.bioinformatics.ngs.dna_seq_workflow(reference, first, output_directory, second=None, ploidy=2, qc=True, config=None)`
- **Agent/MCP API:** `dnakit.ngs.run_dna_seq_workflow(reference, read1, output_directory, read2=None, ploidy=2, qc=True, ...)`; returns a JSON-serializable stage envelope retaining QC, mapping, BAM and VCF artifacts.
- **Input:** Reference FASTA with `.fai`, FASTQ input, all four programs on PATH, and a new output directory.
- **Output:** Returns a mapping of stage names to `ToolResult` objects, including VCF artifact paths.
- **Example code:**

```python
from dnakit.bioinformatics import ExecutionConfig
from dnakit.bioinformatics.ngs import dna_seq_workflow

stages = dna_seq_workflow("reference.fa", "R1.fq.gz", "new_dna_analysis",
                          second="R2.fq.gz", ploidy=2, config=ExecutionConfig(threads=4))
print(stages["variants"].artifacts["variants.vcf.gz"])
```

- **Example output:**

Returns a dictionary of stage names to `ToolResult`, including `qc`, `mapping`, `bam` and `variants` by default. This external-program example was not executed.

Dependencies and details: [FAQ](../../faq/bioinformatics.md#reads).

## 3) Run a single-sample RNA-seq workflow {#bio-21}

- **Purpose:** Create BAM, transcripts and gene counts from RNA-seq reads.
- **Method:** Run HISAT2, samtools, StringTie and featureCounts.
- **API:** `dnakit.bioinformatics.ngs.rna_seq_workflow(reference, annotation, first, output_directory, second=None, strandedness=0)`
- **Agent/MCP API:** `dnakit.ngs.run_rna_seq_workflow(reference, annotation, read1, output_directory, read2=None, strandedness=0, ...)`; returns a JSON-serializable stage envelope retaining BAM, transcripts, counts and logs.
- **Input:** Reference FASTA, GTF and FASTQ; required programs on PATH; `strandedness` is 0/1/2; a new output directory.
- **Output:** Returns a mapping of stage names to `ToolResult` objects, including expression-count artifact paths.
- **Example code:**

```python
from dnakit.bioinformatics.ngs import rna_seq_workflow

stages = rna_seq_workflow("reference.fa", "genes.gtf", "R1.fq.gz",
                          "new_rna_analysis", second="R2.fq.gz", strandedness=0)
print(stages["rna"].artifacts["counts.tsv"])
```

- **Example output:**

Returns a dictionary with `bam` and `rna`; the RNA stage contains `transcripts.gtf` and `counts.tsv`. This external-program example was not executed.

Dependencies and details: [FAQ](../../faq/bioinformatics.md#reads).

## 4) Call peaks from ChIP-seq BAM {#bio-22}

- **Purpose:** Compare treatment/control BAM files to produce peaks and coverage tracks.
- **Method:** Run MACS3 callpeak.
- **API:** `dnakit.bioinformatics.ngs.chip_seq_workflow(treatment, control, output_directory, genome_size=..., paired_end=False, qvalue=0.05)`
- **Agent/MCP API:** `dnakit.ngs.call_chip_peaks(treatment, control, output_directory, genome_size=..., paired_end=False, qvalue=0.05, ...)`; returns a JSON-serializable stage envelope, runs fixed MACS3 `callpeak` and collects samtools coverage.
- **Citation and license:** The survey records MACS3 as BSD-3-Clause but supplies no DOI; cite the official source for the exact version used in a paper.
- **Input:** Prepared treatment/control BAM files, effective genome size and library layout; MACS3 on PATH; a new output directory.
- **Output:** Returns `ToolResult` with `backend`, artifact paths in `artifacts`, and execution `log`.
- **Example code:**

```python
from dnakit.bioinformatics.ngs import chip_seq_workflow

result = chip_seq_workflow("treatment.bam", "control.bam", "new_chip_analysis",
                           genome_size=2700000000, paired_end=True)
print(sorted(result.artifacts))
```

- **Example output:**

Returns `ToolResult`; `artifacts` contains the peak and coverage files produced by MACS3. This external-program example was not executed.

Dependencies and details: [FAQ](../../faq/bioinformatics.md#reads).

## 5) Build a reference-coordinate Sanger consensus {#bio-56}

- **Purpose:** Combine one or more Sanger reads to inspect a construct sequence.
- **Method:** Locally align trimmed reads in both orientations, filter by quality, identity and coverage, then vote per position.
- **API:** `dnakit.bioinformatics.cloning.sanger_consensus(traces, reference, min_quality=20, min_identity=0.8, min_read_coverage=0.8, min_depth=1)`
- **Input:** Trace dictionaries with `id/sequence/qualities` plus a reference; requires `bio-analysis`. N denotes ties or insufficient coverage; consensus omits indels.
- **Output:** Returns a mapping with a consensus sequence, depth array, and variant information.
- **Example code:**

```python
from dnakit.bioinformatics.cloning import sanger_consensus

trace = {"id": "synthetic", "sequence": "AACCGGTT", "qualities": (30,) * 8}
result = sanger_consensus([trace], "AACCGGTT")
print(result["consensus"])
print(result["depth"].tolist())
```

- **Example output:**

```text
AACCGGTT
[1, 1, 1, 1, 1, 1, 1, 1]
```

Dependencies and details: [FAQ](../../faq/bioinformatics.md#optimization).
