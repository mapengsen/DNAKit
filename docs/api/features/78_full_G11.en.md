# Epigenome functions: methylpy, MACS3 and deepTools

This group lives in `dnakit.epigenome` and wraps the Python APIs of methylpy, MACS3 and deepTools. Each function calls the upstream Python API directly; none runs a command line or downloads data. Optional packages are imported only when a function runs, and a missing package raises an install hint. methylpy and deepTools call `sys.exit` on some bad inputs; those exits become `ConfigurationError`.

This batch wraps 32 functions: 10 from methylpy, 12 from MACS3 and 10 from deepTools. The examples import from the module files (for example `dnakit.epigenome`); once exported from `dnakit.epigenome`, the same functions are callable as `dnakit.epigenome.<name>()`. Functions not wrapped (command-line entry points, plotting, pipelines that need samtools, bowtie2, cutadapt, bgzip, tabix or wigToBigWig, and GC-bias functions that depend on module-global state) are listed in the ledger with a reason each.

## Example data

The example data are small text files. The examples run in order in one directory, and later examples use files that earlier ones create.

Indexed FASTA for reading a chromosome.

**Example file:** `genome.fa`

```text
>chr1
ACGTACGTNNACGT
>chr2
GGCC
```

FASTQ for C-to-T conversion.

**Example file:** `reads.fq`

```text
@r1
ACGTCC
+
IIIIII
@r2
TTCG
+
IIII
```

30 bp reference for the mCH filter.

**Example file:** `ref.fa`

```text
>chr1
ACAACAACAACAACAACAACAAAAAAAAAA
```

Two aligned reads, one methylated and one unmethylated (SAM text, converted to BAM in the example).

**Example file:** `reads.sam`

```text
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:30
r1	0	chr1	1	60	20M	*	0	0	ACAACAACAACAACAACAAC	IIIIIIIIIIIIIIIIIIII
r2	0	chr1	1	60	20M	*	0	0	ATAATAATAATAATAATAAT	IIIIIIIIIIIIIIIIIIII
```

Binomial-test results sorted by column 7 (P value).

**Example file:** `pvalues.tsv`

```text
chr1	10	+	CGN	8	10	0.0001
chr1	40	+	CHH	5	9	0.001
chr1	20	+	CGN	1	10	0.2
chr1	30	+	CHH	0	9	0.6
```

An allc table with contiguous chromosomes.

**Example file:** `allc.tsv`

```text
chr1	10	+	CGN	8	10	1
chr1	20	+	CGN	1	10	1
chr2	30	+	CHH	0	9	1
```

Two BED-format tags.

**Example file:** `tags.bed`

```text
chr1	10	30	r1	0	+
chr1	20	40	r2	0	-
```

Two bedGraph intervals.

**Example file:** `signal.bedgraph`

```text
chr1	0	10	1.5
chr1	10	20	2
```

**Example file (genome.2bit):** written by the example below from base64; the sequence is `ACGTGGCCAANNGGTTACGT`.

## methylpy functions

### 1) `dnakit.epigenome.expand_methylation_context()`

- **Purpose:** Expand an IUPAC cytosine-context pattern such as `CGN` into the concrete context strings an allc table may contain.
- **Method:** Calls methylpy's `expand_nucleotide_code` and sorts the result. methylpy also keeps single-letter codes from the pattern; so does this function.
- **API:** `dnakit.epigenome.expand_methylation_context(mc_type: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** `mc_type`: 1 to 6 IUPAC nucleotide codes, for example `"CGN"` or `"CHH"`.
- **Output:** `contexts` is the sorted list of strings and `count` is its length.
- **Example code:**

```python
from dnakit.epigenome import expand_methylation_context

result = expand_methylation_context("CGN", verbose=False)
print(result["count"])
print("CGN" in result["contexts"], "CHH" in result["contexts"])
```

- **Example output:**

```text
38
True True
```

### 2) `dnakit.epigenome.encode_read_c_positions()`

- **Purpose:** Record where the cytosines of a read are, as printable offset characters, so that a converted read can be restored later.
- **Method:** Calls methylpy's `encode_c_positions`. Each character stores the gap to the next C, so the encoding has one character per C.
- **API:** `dnakit.epigenome.encode_read_c_positions(seq: 'str', *, is_read2: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** `seq`: a sequence of A, C, G, T and N; `is_read2`: whether the read is the second mate.
- **Output:** `encoding` is the offset string; the other fields echo the input.
- **Example code:**

```python
from dnakit.epigenome import encode_read_c_positions

encoded = encode_read_c_positions("ACGTCCGA", verbose=False)
print(encoded["encoding"])
print(len(encoded["encoding"]))
```

- **Example output:**

```text
#$"
3
```

### 3) `dnakit.epigenome.decode_read_c_positions()`

- **Purpose:** Restore the original cytosines of a converted read from the offset encoding made by `encode_read_c_positions`.
- **Method:** Calls methylpy's `decode_c_positions`. On the plus strand the positions become C; on the minus strand they become G, because SAM sequences are stored on the forward strand.
- **API:** `dnakit.epigenome.decode_read_c_positions(seq: 'str', indexes: 'str', strand: 'str', *, is_read2: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** `seq`: the converted sequence; `indexes`: the encoding string (may be empty); `strand`: `"+"` or `"-"`.
- **Output:** `sequence` is the restored sequence.
- **Example code:**

```python
from dnakit.epigenome import decode_read_c_positions, encode_read_c_positions

encoding = encode_read_c_positions("ACGTCCGA", verbose=False)["encoding"]
decoded = decode_read_c_positions("ATGTTTGA", encoding, "+", verbose=False)
print(decoded["sequence"])
```

- **Example output:**

```text
ACGTCCGA
```

### 4) `dnakit.epigenome.count_bisulfite_basecalls()`

- **Purpose:** Count consistent and inconsistent calls of one reference base in pileup characters, plus the converted and unconverted cytosines.
- **Method:** Calls methylpy's `analyze_read_basecalls`. `.` and `,` match the reference, `T` (or `a`) marks a converted C, and insertions and deletions are skipped.
- **API:** `dnakit.epigenome.count_bisulfite_basecalls(ref_base: 'str', read_bases: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** `ref_base`: `A`, `C`, `G`, `T` or `N`; `read_bases`: the pileup base string.
- **Output:** `consistent_calls` and `inconsistent_calls`; for C or G also `unconverted_c` and `converted_c`, otherwise `None`.
- **Example code:**

```python
from dnakit.epigenome import count_bisulfite_basecalls

result = count_bisulfite_basecalls("C", ".,,TTA", verbose=False)
print(result["consistent_calls"], result["inconsistent_calls"])
print(result["unconverted_c"], result["converted_c"])
```

- **Example output:**

```text
3 1
1 2
```

### 5) `dnakit.epigenome.read_fasta_chromosome()`

- **Purpose:** Read the full sequence of one chromosome from an indexed FASTA file.
- **Method:** Checks the length in the `.fai` index first, then calls methylpy's `get_chromosome_sequence`; a mismatch between the read length and the index is an error.
- **API:** `dnakit.epigenome.read_fasta_chromosome(fasta_path: 'str | Path', chrom: 'str', *, max_length: 'int' = 50000000, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** A FASTA file with its `.fai` index (made by `pysam.faidx` or `samtools faidx`); `chrom`; `max_length` as the size limit.
- **Output:** `sequence`, `length`, and `provenance`.
- **Example code:**

```python
import pysam
from dnakit.epigenome import read_fasta_chromosome

pysam.faidx("genome.fa")
result = read_fasta_chromosome("genome.fa", "chr1", verbose=False)
print(result["length"], result["sequence"])
```

- **Example output:**

```text
14 ACGTACGTNNACGT
```

### 6) `dnakit.epigenome.convert_fastq_to_bisulfite()`

- **Purpose:** Convert every C to T in the reads of a FASTQ file, for alignment to a bisulfite-converted reference.
- **Method:** Calls methylpy's `convert_reads`. The original C positions are written into the read name, so `decode_read_c_positions` can restore them. The output file must not exist.
- **API:** `dnakit.epigenome.convert_fastq_to_bisulfite(fastq_path: 'str | Path', output_path: 'str | Path', *, buffer_line_number: 'int' = 100000, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** `fastq_path`: the FASTQ file; `output_path`: a new output file; `buffer_line_number`: the buffer size in lines.
- **Output:** `reads` is the read count; the output file holds the converted FASTQ.
- **Example code:**

```python
from pathlib import Path
from dnakit.epigenome import convert_fastq_to_bisulfite

Path("converted.fq").unlink(missing_ok=True)
result = convert_fastq_to_bisulfite("reads.fq", "converted.fq", verbose=False)
print(result["reads"])
print(Path("converted.fq").read_text().splitlines()[1])
```

- **Example output:**

```text
2
ATGTTT
```

### 7) `dnakit.epigenome.filter_bam_by_mch_level()`

- **Purpose:** Remove aligned reads whose non-CpG methylation (mCH) is high.
- **Method:** Calls methylpy's `bam_quality_mch_filter`. A read is removed when it has at least `min_ch` CH sites and mCH/(mCH+uCH) is at least `max_mch_level`. methylpy accepts `min_mapq` but does not apply it, so this function does not expose it.
- **API:** `dnakit.epigenome.filter_bam_by_mch_level(bam_path: 'str | Path', output_bam: 'str | Path', reference_fasta: 'str | Path', *, min_ch: 'int' = 3, max_mch_level: 'float' = 0.7, buffer_line_number: 'int' = 100000, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** A coordinate-sorted, indexed BAM; a reference FASTA with a `.fai` index; a new output path ending in `.bam`.
- **Output:** `input_reads`, `kept_reads`, and `removed_reads`.
- **Example code:**

```python
import pysam
from pathlib import Path
for name in ("reads.unsorted.bam", "reads.bam", "reads.bam.bai"):
    Path(name).unlink(missing_ok=True)
pysam.view("-b", "-o", "reads.unsorted.bam", "reads.sam", catch_stdout=False)
pysam.sort("-o", "reads.bam", "reads.unsorted.bam")
pysam.index("reads.bam")
import pysam
from pathlib import Path
from dnakit.epigenome import filter_bam_by_mch_level

pysam.faidx("ref.fa")
Path("filtered.bam").unlink(missing_ok=True)
result = filter_bam_by_mch_level("reads.bam", "filtered.bam", "ref.fa", verbose=False)
print(result["input_reads"], result["kept_reads"], result["removed_reads"])
```

- **Example output:**

```text
2 1 1
```

### 8) `dnakit.epigenome.benjamini_hochberg_methylation_cutoffs()`

- **Purpose:** Find the P-value cutoff for each cytosine context at a Benjamini-Hochberg false discovery rate.
- **Method:** Calls methylpy's `benjamini_hochberg_correction_call_methylated_sites`. The input must be sorted by P value (column 7) in ascending order; the function checks the order and the test count of each context first.
- **API:** `dnakit.epigenome.benjamini_hochberg_methylation_cutoffs(binomial_paths: 'Sequence[str | Path]', mc_class_counts: 'Mapping[str, int]', *, sig_cutoff: 'float' = 0.01, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** A list of binomial-test files; `mc_class_counts`: the number of tests per context; `sig_cutoff`: the FDR threshold.
- **Output:** `cutoffs` gives the P-value cutoff per context (0 when no site is significant); `tested_sites` is the number of sites read.
- **Example code:**

```python
from dnakit.epigenome import benjamini_hochberg_methylation_cutoffs

result = benjamini_hochberg_methylation_cutoffs(
    ["pvalues.tsv"], {"CGN": 2, "CHH": 2}, sig_cutoff=0.5, verbose=False
)
print(result["tested_sites"])
print({key: round(value, 6) for key, value in result["cutoffs"].items()})
```

- **Example output:**

```text
4
{'CGN': 0.2, 'CHH': 0.6}
```

### 9) `dnakit.epigenome.index_allc_table()`

- **Purpose:** Build a per-chromosome byte-offset index for an allc table, so that a chromosome can be reached without a full scan.
- **Method:** Calls methylpy's `index_allc_file`. The index `<allc>.idx` is written next to the input, and the first byte offset of each chromosome is read back. Records of one chromosome must be contiguous.
- **API:** `dnakit.epigenome.index_allc_table(allc_path: 'str | Path', *, reindex: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** An allc file sorted by chromosome; `reindex=True` forces the index to be rebuilt.
- **Output:** `index` is the index path; `chromosomes` lists each `chrom` with its `offset`.
- **Example code:**

```python
from pathlib import Path
from dnakit.epigenome import index_allc_table

Path("allc.tsv.idx").unlink(missing_ok=True)
result = index_allc_table("allc.tsv", reindex=True, verbose=False)
for item in result["chromosomes"]:
    print(item["chrom"], item["offset"])
```

- **Example output:**

```text
chr1 0
chr2 42
```

### 10) `dnakit.epigenome.read_allc_index()`

- **Purpose:** Read the chromosome offsets stored in the index of an allc table.
- **Method:** Reads `<allc>.idx` through methylpy's `read_allc_index`. A missing index is reported with a pointer to `index_allc_table()`.
- **API:** `dnakit.epigenome.read_allc_index(allc_path: 'str | Path', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** The allc file path; its `.idx` index must already exist.
- **Output:** `chromosomes` is a list of `chrom` and `offset` pairs sorted by chromosome name.
- **Example code:**

```python
from dnakit.epigenome import read_allc_index

result = read_allc_index("allc.tsv", verbose=False)
print([(item["chrom"], item["offset"]) for item in result["chromosomes"]])
```

- **Example output:**

```text
[('chr1', 0), ('chr2', 42)]
```

## MACS3 functions

### 11) `dnakit.epigenome.build_paired_fragment_pileup()`

- **Purpose:** Compute the coverage (pileup) of paired-end fragments from their left and right coordinates, for fragment-level signals such as ATAC-seq.
- **Method:** Calls MACS3's `pileup_from_LR_as_list`. The output is a segmented coverage: each position is the end of a segment and each value is the coverage over that segment.
- **API:** `dnakit.epigenome.build_paired_fragment_pileup(left_positions: 'Sequence[int]', right_positions: 'Sequence[int]', *, scale_factor: 'float' = 1.0, baseline_value: 'float' = 0.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** `left_positions` and `right_positions`: equal-length integer lists; each right coordinate must be greater than its left one.
- **Output:** `positions` and `values` (equal length) and `fragments`.
- **Example code:**

```python
from dnakit.epigenome import build_paired_fragment_pileup

result = build_paired_fragment_pileup([10, 20, 25], [30, 40, 35], verbose=False)
print(result["positions"])
print(result["values"])
```

- **Example output:**

```text
[10, 20, 25, 30, 35, 40]
[0.0, 1.0, 2.0, 3.0, 2.0, 1.0]
```

### 12) `dnakit.epigenome.build_shifted_read_pileup()`

- **Purpose:** Extend single-end reads into shifted fragments on each strand and compute their coverage, for example Tn5 offsets in ATAC-seq.
- **Method:** Calls MACS3's `pileup_from_PN_shifted`. A plus-strand read covers `[p - five_shift, p + three_shift)` and a minus-strand read covers `[p - three_shift, p + five_shift)`, clipped to the chromosome length.
- **API:** `dnakit.epigenome.build_shifted_read_pileup(plus_positions: 'Sequence[int]', minus_positions: 'Sequence[int]', *, five_shift: 'int', three_shift: 'int', chrom_length: 'int', scale_factor: 'float' = 1.0, baseline_value: 'float' = 0.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** 5' positions of plus and minus reads; `five_shift` and `three_shift`; `chrom_length`.
- **Output:** `positions` and `values`, with the read counts and shifts.
- **Example code:**

```python
from dnakit.epigenome import build_shifted_read_pileup

result = build_shifted_read_pileup(
    [10, 12, 30], [40, 42], five_shift=0, three_shift=10, chrom_length=100, verbose=False
)
print(result["positions"])
print(result["values"])
```

- **Example output:**

```text
[10, 12, 20, 22, 30, 32, 40, 42]
[0.0, 1.0, 2.0, 1.0, 0.0, 2.0, 3.0, 1.0]
```

### 13) `dnakit.epigenome.find_signal_maxima()`

- **Purpose:** Find local maxima in a one-dimensional signal, such as a coverage profile, as candidate peak summits.
- **Method:** Calls MACS3's `SignalProcessing.maxima`, which applies a second-order Savitzky-Golay smoothing before taking local maxima. `window_size` must be odd and no longer than the signal.
- **API:** `dnakit.epigenome.find_signal_maxima(signal: 'Sequence[float]', *, window_size: 'int' = 51, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** A list of at least three numbers and an odd `window_size`.
- **Output:** `maxima` is a list of 0-based indices.
- **Example code:**

```python
from dnakit.epigenome import find_signal_maxima

signal = [0, 1, 2, 3, 2, 1, 0, 1, 4, 6, 4, 1, 0, 0, 0]
print(find_signal_maxima(signal, window_size=3, verbose=False)["maxima"])
```

- **Example output:**

```text
[2, 3, 8, 9]
```

### 14) `dnakit.epigenome.filter_summits_by_peakiness()`

- **Purpose:** Keep the summits whose surrounding peak is wide enough and varied enough, and drop narrow or flat peaks.
- **Method:** Calls MACS3's `SignalProcessing.enforce_peakyness`. A peak must be at least 50 bp wide and have at least six distinct values.
- **API:** `dnakit.epigenome.filter_summits_by_peakiness(signal: 'Sequence[float]', summits: 'Sequence[int]', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** A signal list and a list of summit indices (0-based, ascending, unique).
- **Output:** `summits` lists the kept indices.
- **Example code:**

```python
import math
from dnakit.epigenome import filter_summits_by_peakiness, find_signal_maxima

signal = [round(100 * math.exp(-((i - 60) ** 2) / 450), 4) for i in range(121)]
summits = find_signal_maxima(signal, window_size=21, verbose=False)["maxima"]
print(summits)
print(filter_summits_by_peakiness(signal, summits, verbose=False)["summits"])
```

- **Example output:**

```text
[59]
[59]
```

### 15) `dnakit.epigenome.filter_summits_by_valleys()`

- **Purpose:** Check the valley between neighbouring summits and keep the summits whose valley is low enough.
- **Method:** Calls MACS3's `SignalProcessing.enforce_valleys`. The valley between two summits must be at most `min_valley` times the lower summit.
- **API:** `dnakit.epigenome.filter_summits_by_valleys(signal: 'Sequence[float]', summits: 'Sequence[int]', *, min_valley: 'float' = 0.8, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** A signal list, a summit index list, and `min_valley` between 0 and 1.
- **Output:** `summits` lists the kept indices.
- **Example code:**

```python
from dnakit.epigenome import filter_summits_by_valleys

signal = [0, 5, 10, 5, 1, 4, 8, 4, 0]
print(filter_summits_by_valleys(signal, [2, 6], min_valley=0.8, verbose=False)["summits"])
print(filter_summits_by_valleys(signal, [2, 6], min_valley=0.1, verbose=False)["summits"])
```

- **Example output:**

```text
[2, 6]
[2]
```

### 16) `dnakit.epigenome.refine_peak_summit()`

- **Purpose:** Locate the summit of one called peak from the 5' read positions on each strand.
- **Method:** Calls MACS3's `refinepeak_cmd.find_summit`, which smooths the read ends in a window and takes the densest position. Coordinates are as MACS3 returns them.
- **API:** `dnakit.epigenome.refine_peak_summit(chrom: 'str', plus_positions: 'Sequence[int]', minus_positions: 'Sequence[int]', peak_start: 'int', peak_end: 'int', *, name: 'str' = 'peak', window_size: 'int' = 100, cutoff: 'float' = 5.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** Chromosome name, 5' positions on each strand, the peak start and end, `window_size`, and `cutoff`.
- **Output:** `summit_start`, `summit_end`, `name`, and `score`.
- **Example code:**

```python
from dnakit.epigenome import refine_peak_summit

result = refine_peak_summit(
    "chr1", [100, 110, 120], [105, 115], 90, 130, name="peak1", verbose=False
)
print(result["summit_start"], result["summit_end"], result["name"])
print(round(result["score"], 4))
```

- **Example output:**

```text
101 102 peak1_F
0.8284
```

### 17) `dnakit.epigenome.max_duplicate_tags()`

- **Purpose:** Estimate the largest number of identical tags expected by chance at one position at a given P value, as a duplicate threshold.
- **Method:** Calls MACS3's `callpeak_cmd.cal_max_dup_tags`, which searches the inverse binomial distribution by brute force.
- **API:** `dnakit.epigenome.max_duplicate_tags(genome_size: 'int', tag_count: 'int', *, p_value: 'float' = 1e-05, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** Genome size, total tag count, and `p_value`.
- **Output:** `max_duplicate_tags` is the integer threshold.
- **Example code:**

```python
from dnakit.epigenome import max_duplicate_tags

result = max_duplicate_tags(3_000_000_000, 10_000_000, verbose=False)
print(result["max_duplicate_tags"])
```

- **Example output:**

```text
1
```

### 18) `dnakit.epigenome.guess_tag_file_format()`

- **Purpose:** Identify which MACS3 parser reads a tag file, from its first records.
- **Method:** Calls MACS3's `IO.Parser.guess_parser`, which returns the first parser that accepts the file. An unrecognised file is an error.
- **API:** `dnakit.epigenome.guess_tag_file_format(path: 'str | Path', *, buffer_size: 'int' = 100000, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** A tag file path (BED, SAM, and similar) and the read buffer size `buffer_size`.
- **Output:** `parser` is the parser class name; `format` is the same name without the `Parser` suffix.
- **Example code:**

```python
from dnakit.epigenome import guess_tag_file_format

result = guess_tag_file_format("tags.bed", verbose=False)
print(result["parser"], result["format"])
```

- **Example output:**

```text
BEDParser BED
```

### 19) `dnakit.epigenome.calculate_genotype_quality()`

- **Purpose:** Compute the Phred-scaled genotype quality (GQ) from the natural-log likelihoods of three genotypes.
- **Method:** Calls MACS3's `VariantStat.calculate_GQ`, which computes `GQ = -10 log10((L2 + L3) / (L1 + L2 + L3))` with `L = exp(lnL)`. MACS3 returns an integer.
- **API:** `dnakit.epigenome.calculate_genotype_quality(log_likelihood_1: 'float', log_likelihood_2: 'float', log_likelihood_3: 'float', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** Three natural-log likelihoods: `log_likelihood_1`, `log_likelihood_2`, `log_likelihood_3`.
- **Output:** `genotype_quality` as an integer.
- **Example code:**

```python
from dnakit.epigenome import calculate_genotype_quality

result = calculate_genotype_quality(-10.0, -20.0, -30.0, verbose=False)
print(result["genotype_quality"])
```

- **Example output:**

```text
43
```

### 20) `dnakit.epigenome.calculate_allele_model_likelihood()`

- **Purpose:** Compute the log-likelihood and BIC of a homozygous or heterozygous allele model at one site, from per-read base qualities.
- **Method:** Calls MACS3's `VariantStat.CalModel_Homo` (homozygous) or `CalModel_Heter_noAS` (heterozygous, without allele bias).
- **API:** `dnakit.epigenome.calculate_allele_model_likelihood(top1_bq_t: 'Sequence[int]', top1_bq_c: 'Sequence[int]', top2_bq_t: 'Sequence[int]', top2_bq_c: 'Sequence[int]', *, model: 'str' = 'homozygous', verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** Four equal-length lists of non-negative integers for the T and C qualities of the first and second alleles; `model` is `"homozygous"` or `"heterozygous"`.
- **Output:** `log_likelihood`, `bic`, and `reads`.
- **Example code:**

```python
from dnakit.epigenome import calculate_allele_model_likelihood

result = calculate_allele_model_likelihood(
    [30, 25], [2, 3], [1, 1], [0, 2], model="heterozygous", verbose=False
)
print(result["reads"], round(result["log_likelihood"], 4), round(result["bic"], 4))
```

- **Example output:**

```text
2 -6.5882 15.9489
```

### 21) `dnakit.epigenome.hmmratac_fragment_weights()`

- **Purpose:** Weight each fragment length for the four HMMRATAC signal classes: short, mono, di, and tri.
- **Method:** Calls MACS3's `HMMR_Signal_Processing.generate_weight_mapping`. Each class is a Gaussian with one mean and one standard deviation.
- **API:** `dnakit.epigenome.hmmratac_fragment_weights(fragment_lengths: 'Sequence[int]', means: 'Sequence[float]', stddevs: 'Sequence[float]', *, min_fragment_probability: 'float' = 0.001, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** A list of fragment lengths; four means `means` and four standard deviations `stddevs`, in the order short, mono, di, tri; `min_fragment_probability`.
- **Output:** `weights` has one row per fragment length with the four class weights; `classes` lists the class names.
- **Example code:**

```python
from dnakit.epigenome import hmmratac_fragment_weights

result = hmmratac_fragment_weights(
    [50, 150, 300],
    [50.0, 150.0, 300.0, 450.0],
    [20.0, 30.0, 40.0, 50.0],
    verbose=False,
)
for row in result["weights"]:
    print(row["fragment_length"], round(row["short"], 4), round(row["mono"], 4))
```

- **Example output:**

```text
50 0.9974 0.0026
150 0.0 0.9993
300 0.0 0.0
```

### 22) `dnakit.epigenome.merge_pileup_arrays()`

- **Purpose:** Merge two coverage tracks position by position, taking the maximum, minimum, or mean where they overlap.
- **Method:** Calls MACS3's `Pileup.over_two_pv_array`, which evaluates the reducer at each position of the union of the two `[positions, values]` tracks.
- **API:** `dnakit.epigenome.merge_pileup_arrays(first: 'Mapping[str, Sequence[float]]', second: 'Mapping[str, Sequence[float]]', *, reducer: 'str' = 'max', verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** Two mappings, each with `positions` and `values`; `reducer` is `"max"`, `"min"` or `"mean"`.
- **Output:** The merged `positions` and `values`.
- **Example code:**

```python
from dnakit.epigenome import merge_pileup_arrays

first = {"positions": [10, 20], "values": [1.0, 2.0]}
second = {"positions": [15, 20], "values": [3.0, 0.5]}
result = merge_pileup_arrays(first, second, reducer="max", verbose=False)
print(result["positions"])
print(result["values"])
```

- **Example output:**

```text
[10, 15, 20]
[3.0, 3.0, 2.0]
```

## deeptools functions

### 23) `dnakit.epigenome.gc_fraction_from_twobit()`

- **Purpose:** Compute the GC fraction (or the G+C base count) of one interval in a 2bit genome.
- **Method:** Opens the file with py2bit and calls deepTools' `utilities.getGC_content`. Intervals with more than 5% N bases are rejected.
- **API:** `dnakit.epigenome.gc_fraction_from_twobit(twobit_path: 'str | Path', chrom: 'str', start: 'int', end: 'int', *, fraction: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** A 2bit file, a chromosome name, and a 0-based half-open interval `[start, end)`; `fraction=True` returns the fraction.
- **Output:** `gc` (fraction or count), `unit`, `start`, and `end`.
- **Example code:**

```python
import base64
from pathlib import Path
from dnakit.epigenome import gc_fraction_from_twobit

Path("genome.2bit").write_bytes(base64.b64decode("QydBGgAAAAABAAAAAAAAAARjaHIxGQAAABQAAAABAAAACgAAAAIAAAAAAAAAAAAAAJz1oPCc"))
print(gc_fraction_from_twobit("genome.2bit", "chr1", 0, 8, verbose=False)["gc"])
print(gc_fraction_from_twobit("genome.2bit", "chr1", 0, 8, fraction=False, verbose=False)["gc"])
```

- **Example output:**

```text
0.75
6
```

### 24) `dnakit.epigenome.parse_deeptools_region()`

- **Purpose:** Parse a deepTools-style `chrom:start:end[:tile]` region and correct it to the chromosome lengths.
- **Method:** Calls deepTools' `mapReduce.getUserRegion`. The start is aligned to the tile size and the end is clipped to the chromosome. An unknown chromosome is rejected first.
- **API:** `dnakit.epigenome.parse_deeptools_region(chrom_sizes: 'Mapping[str, int]', region: 'str', *, max_chunk_size: 'int' = 1000000, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** `chrom_sizes` (chromosome to length), the region string, and `max_chunk_size`.
- **Output:** `regions`, `start`, `end`, and `chunk_size`.
- **Example code:**

```python
from dnakit.epigenome import parse_deeptools_region

result = parse_deeptools_region({"chr2": 1000}, "chr2:10:1001", verbose=False)
print(result["regions"])
print(result["start"], result["end"], result["chunk_size"])
```

- **Example output:**

```text
[{'chrom': 'chr2', 'length': 1000}]
10 1000 990
```

### 25) `dnakit.epigenome.estimate_bam_scale_factors()`

- **Purpose:** Estimate scale factors from two indexed BAM files, so that the coverage of two libraries can be compared.
- **Method:** Calls deepTools' `estimateScaleFactor` (the SES method), which samples `number_of_samples` evenly spaced sites across the genome and compares read counts.
- **API:** `dnakit.epigenome.estimate_bam_scale_factors(bam_paths: 'Sequence[str | Path]', *, bin_length: 'int' = 500, number_of_samples: 'int' = 1000, normalization_length: 'int' = 1000, avg_method: 'str' = 'median', verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** Exactly two BAM files, each with a `.bai` index; `bin_length`, `number_of_samples`, `normalization_length`, and `avg_method` (`mean` or `median`).
- **Output:** `size_factors`, `sites_sampled`, and the other estimates in `estimates`.
- **Example code:**

```python
import random
import pysam
from pathlib import Path

random.seed(1)
sequence = "".join(random.choice("ACGT") for _ in range(400))
rows = []
for index in range(60):
    position = 1 + (index * 6) % 300
    bases = sequence[position - 1 : position + 19]
    rows.append(f"r{index}\t0\tchr1\t{position}\t60\t20M\t*\t0\t0\t{bases}\t{'I' * 20}")
sam = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:400\n" + "\n".join(rows) + "\n"
for name in ("scale_a", "scale_b"):
    Path(name + ".sam").write_text(sam)
    pysam.view("-b", "-o", name + ".unsorted.bam", name + ".sam", catch_stdout=False)
    pysam.sort("-o", name + ".bam", name + ".unsorted.bam")
    pysam.index(name + ".bam")
from dnakit.epigenome import estimate_bam_scale_factors

result = estimate_bam_scale_factors(
    ["scale_a.bam", "scale_b.bam"], bin_length=50, number_of_samples=40, verbose=False
)
print(result["size_factors"], result["sites_sampled"])
```

- **Example output:**

```text
[1.0, 1.0] 36
```

### 26) `dnakit.epigenome.bedgraph_to_bigwig()`

- **Purpose:** Write sorted bedGraph files to a bigWig track file.
- **Method:** Checks the chromosome order, that the intervals are sorted and non-overlapping, and that they lie within the chromosome lengths; then calls deepTools' `bedGraphToBigWig` (pyBigWig).
- **API:** `dnakit.epigenome.bedgraph_to_bigwig(chrom_sizes: 'Mapping[str, int]', bedgraph_paths: 'Sequence[str | Path]', bigwig_path: 'str | Path', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** A chromosome length mapping (its order is the bigWig header order), the bedGraph files in the same order, and a new bigWig output path.
- **Output:** `bigwig`, `chromosomes`, `intervals`, and related fields.
- **Example code:**

```python
from pathlib import Path
import pyBigWig
from dnakit.epigenome import bedgraph_to_bigwig

Path("signal.bw").unlink(missing_ok=True)
result = bedgraph_to_bigwig({"chr1": 100}, ["signal.bedgraph"], "signal.bw", verbose=False)
print(result["intervals"])
handle = pyBigWig.open("signal.bw")
print(list(handle.intervals("chr1")))
handle.close()
```

- **Example output:**

```text
2
[(0, 10, 1.5), (10, 20, 2.0)]
```

### 27) `dnakit.epigenome.chop_regions_for_profile()`

- **Purpose:** Split exon intervals into unscaled 5' and 3' flanks and a scaled body, as used for profile plots.
- **Method:** Calls deepTools' `heatmapper.chopRegions`. The input intervals must be sorted and non-overlapping.
- **API:** `dnakit.epigenome.chop_regions_for_profile(exons: 'Sequence[Sequence[int]]', *, left: 'int' = 0, right: 'int' = 0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** A list of `(start, end)` intervals; `left` and `right` give the unscaled flank length at each end.
- **Output:** `left_bins`, `body_bins`, `right_bins`, and `pad_left` and `pad_right` (bases missing from a short flank).
- **Example code:**

```python
from dnakit.epigenome import chop_regions_for_profile

result = chop_regions_for_profile(
    [(0, 200), (300, 400), (800, 900)], left=50, right=50, verbose=False
)
print(result["left_bins"], result["right_bins"])
print(result["body_bins"])
```

- **Example output:**

```text
[[0, 50]] [[850, 900]]
[[50, 200], [300, 400], [800, 850]]
```

### 28) `dnakit.epigenome.match_chromosome_name()`

- **Purpose:** Find the chromosome name in a list that matches a request, allowing the chr prefix to be added or removed.
- **Method:** Calls deepTools' `utilities.mungeChromosome`. The result is `None` when no variant is listed.
- **API:** `dnakit.epigenome.match_chromosome_name(chrom: 'str', chromosome_names: 'Sequence[str]', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** The requested chromosome name and the list of available names.
- **Output:** `matched`: the name found, or `None`.
- **Example code:**

```python
from dnakit.epigenome import match_chromosome_name

print(match_chromosome_name("chr1", ["1", "2"], verbose=False)["matched"])
print(match_chromosome_name("X", ["chr1"], verbose=False)["matched"])
```

- **Example output:**

```text
1
None
```

### 29) `dnakit.epigenome.map_twobit_to_bam_chromosomes()`

- **Purpose:** Map 2bit chromosome names to BAM reference names, allowing a difference in the chr prefix.
- **Method:** Calls deepTools' `utilities.tbitToBamChrName`. `complete` is true when every 2bit name was mapped.
- **API:** `dnakit.epigenome.map_twobit_to_bam_chromosomes(twobit_names: 'Sequence[str]', bam_names: 'Sequence[str]', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** A list of 2bit chromosome names and a list of BAM reference names.
- **Output:** `mapping` (2bit name to BAM name) and `complete`.
- **Example code:**

```python
from dnakit.epigenome import map_twobit_to_bam_chromosomes

result = map_twobit_to_bam_chromosomes(["chr1", "chr2"], ["1", "2"], verbose=False)
print(result["mapping"], result["complete"])
partial = map_twobit_to_bam_chromosomes(["chr1", "chrX"], ["chr1", "chr2"], verbose=False)
print(partial["complete"])
```

- **Example output:**

```text
{'chr1': '1', 'chr2': '2'} True
False
```

### 30) `dnakit.epigenome.fingerprint_jensen_shannon()`

- **Purpose:** Compute the Jensen-Shannon measure between ChIP and input coverage, the JSD that plotFingerprint reports.
- **Method:** Calls deepTools' `plotFingerprint.getJSDcommon`, which compares the cumulative coverage curves of the two. Counts are used in the order given.
- **API:** `dnakit.epigenome.fingerprint_jensen_shannon(chip_counts: 'Sequence[float]', input_counts: 'Sequence[float]', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** Per-bin count lists for ChIP and input (equal length, at least two bins).
- **Output:** `jsd` is the value and `bins` is the bin count.
- **Example code:**

```python
from dnakit.epigenome import fingerprint_jensen_shannon

result = fingerprint_jensen_shannon([0, 1, 2, 3, 10], [1, 1, 1, 1, 1], verbose=False)
print(round(result["jsd"], 6), result["bins"])
print(fingerprint_jensen_shannon([1, 2, 3, 4], [1, 2, 3, 4], verbose=False)["jsd"])
```

- **Example output:**

```text
0.283567 5
0.0
```

### 31) `dnakit.epigenome.estimate_count_size_factors()`

- **Purpose:** Estimate size factors from a bin-by-sample count matrix; the values can be used directly as bamCoverage scale factors.
- **Method:** Calls deepTools' `countReadsPerBin.estimateSizeFactors`, the DESeq2 method on rows without zero counts, and returns its inverse.
- **API:** `dnakit.epigenome.estimate_count_size_factors(counts: 'Sequence[Sequence[float]]', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** An equal-width non-negative matrix with one row per bin and one column per sample (at least two columns).
- **Output:** `size_factors`, one value per sample.
- **Example code:**

```python
from dnakit.epigenome import estimate_count_size_factors

matrix = [[0, 1, 2], [3, 4, 5], [6, 7, 8], [0, 10, 0], [10, 5, 100]]
result = estimate_count_size_factors(matrix, verbose=False)
print([round(value, 4) for value in result["size_factors"]])
```

- **Example output:**

```text
[1.305, 0.9932, 0.783]
```

### 32) `dnakit.epigenome.bigwig_chromosome_sizes()`

- **Purpose:** Read the chromosome names and lengths of bigWig files and list the entries that differ between them.
- **Method:** Calls deepTools' `getScorePerBigWigBin.getChromSizes`, which reads each file's chromosome table through pyBigWig.
- **API:** `dnakit.epigenome.bigwig_chromosome_sizes(bigwig_paths: 'Sequence[str | Path]', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** One to 32 bigWig file paths.
- **Output:** `chromosomes` (those shared by the files) and `non_matching` (entries that differ).
- **Example code:**

```python
import pyBigWig
from pathlib import Path
from dnakit.epigenome import bigwig_chromosome_sizes

Path("track.bw").unlink(missing_ok=True)
handle = pyBigWig.open("track.bw", "w")
handle.addHeader([("chr1", 100)])
handle.addEntries(["chr1"], [0], ends=[10], values=[1.0])
handle.close()
print(bigwig_chromosome_sizes(["track.bw"], verbose=False)["chromosomes"])
```

- **Example output:**

```text
[{'chrom': 'chr1', 'length': 100}]
```
