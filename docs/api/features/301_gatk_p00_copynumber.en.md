# GATK command-line facades: copy number, read counts, coverage and pileups

Typed facades over the `gatk` command line for the somatic and germline copy-number workflow (interval binning, annotation and filtering, read and allelic counts, panel of normals, denoising, segmentation, calling, plots, germline CNV calling), allele-specific read counts, callable loci, depth of coverage, expression counts and the helper tables of Mutect2. Each function builds a fixed `gatk` command line without a shell, writes only to new paths and returns a JSON-serialisable dictionary; GATK itself is licensed BSD-3-Clause (Picard MIT) and is not bundled. The functions live in `dnakit.variants.gatk_p00_copynumber_cli_functions` and are exported from `dnakit.variants`; the examples import from the module so that they run either way.

## Example data

All examples use one small data set that the first example writes into the current folder: a two-contig reference (`toy.fa`, 1,000 bp) with `.fai` and `.dict`, an indexed BAM and CRAM of 38 paired reads with one read group (sample `S1`), a tabix-indexed VCF of 15 heterozygous sites, a BED file and interval list, a haplotype map, a refFlat and a GFF3 annotation, a chain file and a UMI-tagged BAM. The other examples write their outputs into a temporary folder and pass `verbose=False`. The `gatk` launcher is found on PATH (or `DNAKIT_GATK`); the examples that make charts also need `Rscript` on PATH.

**Example code:**

```python
import random
from pathlib import Path

import pysam

rng = random.Random(7)
chr1 = "".join(rng.choice("ACGT") for _ in range(600))
chr2 = "".join(rng.choice("ACGT") for _ in range(400))
contigs = {"chr1": chr1, "chr2": chr2}
with open("toy.fa", "w") as fasta:
    for name, sequence in contigs.items():
        fasta.write(f">{name}\n")
        fasta.writelines(sequence[i : i + 60] + "\n" for i in range(0, len(sequence), 60))
Path("toy.fa.fai").write_text("chr1\t600\t6\t60\t61\nchr2\t400\t622\t60\t61\n")
Path("toy.dict").write_text("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:600\n@SQ\tSN:chr2\tLN:400\n")
Path("toy.bed").write_text("chr1\t100\t300\tt1\nchr2\t50\t200\tt2\n")
Path("toy.interval_list").write_text(
    "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:600\n@SQ\tSN:chr2\tLN:400\n"
    "chr1\t101\t300\t+\tt1\nchr2\t51\t200\t+\tt2\n"
)
header = {
    "HD": {"VN": "1.6", "SO": "coordinate"},
    "SQ": [{"SN": "chr1", "LN": 600}, {"SN": "chr2", "LN": 400}],
    "RG": [{"ID": "rg1", "SM": "S1", "LB": "lib1", "PL": "ILLUMINA", "PU": "unit1"}],
}
complement = str.maketrans("ACGT", "TGCA")
reads = []
for contig, sequence in contigs.items():
    for start in range(0, len(sequence) - 120, 40):
        name = f"r{len(reads) // 2 + 1}"
        reads.append((contig, start, False, sequence[start : start + 50], name))
        mate = sequence[start + 50 : start + 100][::-1].translate(complement)
        reads.append((contig, start + 50, True, mate, name))
reads.sort(key=lambda item: (item[0], item[1]))
with pysam.AlignmentFile("toy.bam", "wb", header=header) as bam:
    for contig, start, reverse, sequence, name in reads:
        read = pysam.AlignedSegment()
        read.query_name = name
        read.query_sequence = sequence
        read.query_qualities = pysam.qualitystring_to_array("".join(chr(33 + rng.randint(25, 40)) for _ in sequence))
        read.reference_id = 0 if contig == "chr1" else 1
        read.reference_start = start
        read.cigarstring = "50M"
        read.mapping_quality = 60
        read.flag = 147 if reverse else 99
        read.next_reference_id = read.reference_id
        read.next_reference_start = start + 50 if not reverse else start - 50
        read.template_length = -100 if reverse else 100
        read.set_tags([("RG", "rg1"), ("NM", 0)])
        bam.write(read)
pysam.index("toy.bam")
pysam.view("-C", "-T", "toy.fa", "-o", "toy.cram", "toy.bam", catch_stdout=False)
vcf = [
    "##fileformat=VCFv4.2",
    "##contig=<ID=chr1,length=600>",
    "##contig=<ID=chr2,length=400>",
    '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
]
sites = []
for contig, sequence in contigs.items():
    for position in range(40, len(sequence) - 40, 60):
        ref = sequence[position - 1]
        alt = next(base for base in "ACGT" if base != ref)
        sites.append((contig, position, ref, alt))
        vcf.append(f"{contig}\t{position}\trs{position}{contig[-1]}\t{ref}\t{alt}\t50\tPASS\tAF=0.4\tGT\t0/1")
Path("toy.vcf").write_text("\n".join(vcf) + "\n")
pysam.tabix_index("toy.vcf", preset="vcf", force=True, keep_original=True)
Path("toy.haplotype_map").write_text(
    "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:600\n@SQ\tSN:chr2\tLN:400\n"
    + "".join(f"{c}\t{p}\trs{p}{c[-1]}\t{r}\t{a}\t0.4\trs{p}{c[-1]}\n" for c, p, r, a in sites)
)
Path("toy.refFlat").write_text(
    "G1\tt1\tchr1\t+\t0\t600\t0\t600\t2\t0,300,\t200,600,\nG2\tt2\tchr2\t+\t0\t400\t0\t400\t1\t0,\t400,\n"
)
Path("toy.gff3").write_text(
    "##gff-version 3\nchr1\tsrc\tgene\t1\t600\t.\t+\t.\tID=G1;Name=G1\n"
    "chr1\tsrc\texon\t1\t200\t.\t+\t.\tParent=G1;Name=G1\nchr1\tsrc\texon\t301\t600\t.\t+\t.\tParent=G1;Name=G1\n"
    "chr2\tsrc\tgene\t1\t400\t.\t+\t.\tID=G2;Name=G2\nchr2\tsrc\texon\t1\t400\t.\t+\t.\tParent=G2;Name=G2\n"
)
Path("shift.chain").write_text("chain 1000 chr1 590 + 0 590 chr1 600 + 10 600 1\n590\n\n")
Path("old.interval_list").write_text("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:590\nchr1\t101\t200\t+\tt1\n")
Path("old.haplotype_map").write_text("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:590\nchr1\t30\trs1\tT\tA\t0.4\trs1\n")
with pysam.AlignmentFile("toy.bam") as source, pysam.AlignmentFile("umi.bam", "wb", template=source) as out:
    for index, read in enumerate(source):
        read.set_tag("RX", "ACGT"[index % 4] * 3 + "-" + "ACGT"[(index // 4) % 4] * 3)
        out.write(read)
pysam.index("umi.bam")
print("toy data written")
```

**Example output:**

```text
toy data written
```

## 1) `dnakit.variants.gatk_preprocess_intervals()`

- **Purpose:** Prepare equal-sized bins (or targets with padding) for coverage collection and write them as an interval list.
- **Method:** Builds the fixed argument list of `gatk PreprocessIntervals` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk PreprocessIntervals --help`.
- **API:** `dnakit.variants.gatk_preprocess_intervals(reference: 'str', output_path: 'str', *, intervals: 'Sequence[str] | None' = None, bin_length: 'int | None' = None, padding: 'int | None' = None, interval_merging_rule: 'str | None' = 'OVERLAPPING_ONLY', options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `reference`, `output_path`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2908.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `interval_count`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants.gatk_p00_copynumber_cli_functions import gatk_preprocess_intervals

work = Path(tempfile.mkdtemp())
result = gatk_preprocess_intervals("toy.fa", work / "bins.interval_list", bin_length=100, padding=0, verbose=False)
print(result["interval_count"])
```

- **Example output:**

```text
10
```

## 2) `dnakit.variants.gatk_annotate_intervals()`

- **Purpose:** Annotate intervals with GC content and, when tracks are given, mappability and segmental-duplication content.
- **Method:** Builds the fixed argument list of `gatk AnnotateIntervals` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk AnnotateIntervals --help`.
- **API:** `dnakit.variants.gatk_annotate_intervals(intervals: 'Sequence[str]', reference: 'str', output_path: 'str', *, mappability_track: 'str | None' = None, segmental_duplication_track: 'str | None' = None, feature_query_lookahead: 'int | None' = None, interval_merging_rule: 'str | None' = 'OVERLAPPING_ONLY', options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `intervals`, `reference`, `output_path`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2813.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `rows`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants.gatk_p00_copynumber_cli_functions import (
    gatk_annotate_intervals,
    gatk_preprocess_intervals,
)

work = Path(tempfile.mkdtemp())
bins = gatk_preprocess_intervals("toy.fa", work / "bins.interval_list", bin_length=100, padding=0, verbose=False)
result = gatk_annotate_intervals([work / "bins.interval_list"], "toy.fa", work / "annotated.tsv", verbose=False)
print(result["rows"][0])
```

- **Example output:**

```text
{'CONTIG': 'chr1', 'START': 1, 'END': 100, 'GC_CONTENT': 0.42}
```

## 3) `dnakit.variants.gatk_collect_read_counts()`

- **Purpose:** Collect read counts at specified intervals (coverage-based copy-number input) into a TSV or HDF5 file.
- **Method:** Builds the fixed argument list of `gatk CollectReadCounts` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk CollectReadCounts --help`.
- **API:** `dnakit.variants.gatk_collect_read_counts(input_paths: 'Sequence[str]', intervals: 'Sequence[str]', output_path: 'str', *, output_format: 'str | None' = None, reference: 'str | None' = None, interval_merging_rule: 'str | None' = 'OVERLAPPING_ONLY', options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_paths`, `intervals`, `output_path`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2833.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `rows`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants.gatk_p00_copynumber_cli_functions import (
    gatk_collect_read_counts,
    gatk_preprocess_intervals,
)

work = Path(tempfile.mkdtemp())
gatk_preprocess_intervals("toy.fa", work / "bins.interval_list", bin_length=100, padding=0, verbose=False)
result = gatk_collect_read_counts(["toy.bam"], [work / "bins.interval_list"], work / "counts.tsv", output_format="TSV", verbose=False)
print([row["COUNT"] for row in result["rows"]])
```

- **Example output:**

```text
[5, 4, 6, 4, 5, 0, 5, 4, 5, 0]
```

## 4) `dnakit.variants.gatk_collect_allelic_counts()`

- **Purpose:** Collect reference and alternate allele counts at the sites of an interval list or VCF.
- **Method:** Builds the fixed argument list of `gatk CollectAllelicCounts` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk CollectAllelicCounts --help`.
- **API:** `dnakit.variants.gatk_collect_allelic_counts(input_paths: 'Sequence[str]', intervals: 'Sequence[str]', reference: 'str', output_path: 'str', *, minimum_base_quality: 'int | None' = None, max_depth_per_sample: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_paths`, `intervals`, `reference`, `output_path`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2831.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `rows`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants.gatk_p00_copynumber_cli_functions import gatk_collect_allelic_counts

work = Path(tempfile.mkdtemp())
result = gatk_collect_allelic_counts(["toy.bam"], ["toy.vcf"], "toy.fa", work / "allelic.tsv", verbose=False)
print(len(result["rows"]), result["rows"][1])
```

- **Example output:**

```text
15 {'CONTIG': 'chr1', 'POSITION': 100, 'REF_COUNT': 1, 'ALT_COUNT': 2, 'REF_NUCLEOTIDE': 'G', 'ALT_NUCLEOTIDE': 'C'}
```

## 5) `dnakit.variants.gatk_filter_intervals()`

- **Purpose:** Filter intervals by annotations (GC content, mappability, segmental duplications) and by count statistics.
- **Method:** Builds the fixed argument list of `gatk FilterIntervals` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk FilterIntervals --help`.
- **API:** `dnakit.variants.gatk_filter_intervals(intervals: 'Sequence[str]', output_path: 'str', *, annotated_intervals: 'str | None' = None, input_paths: 'Sequence[str] | None' = None, minimum_gc_content: 'float | None' = None, maximum_gc_content: 'float | None' = None, minimum_mappability: 'float | None' = None, maximum_mappability: 'float | None' = None, minimum_segmental_duplication_content: 'float | None' = None, maximum_segmental_duplication_content: 'float | None' = None, low_count_filter_count_threshold: 'int | None' = None, low_count_filter_percentage_of_samples: 'float | None' = None, interval_merging_rule: 'str | None' = 'OVERLAPPING_ONLY', options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `intervals`, `output_path`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2819.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `interval_count`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants.gatk_p00_copynumber_cli_functions import (
    gatk_annotate_intervals,
    gatk_filter_intervals,
    gatk_preprocess_intervals,
)

work = Path(tempfile.mkdtemp())
gatk_preprocess_intervals("toy.fa", work / "bins.interval_list", bin_length=100, padding=0, verbose=False)
gatk_annotate_intervals([work / "bins.interval_list"], "toy.fa", work / "annotated.tsv", verbose=False)
result = gatk_filter_intervals([work / "bins.interval_list"], work / "kept.interval_list", annotated_intervals=work / "annotated.tsv", minimum_gc_content=0.45, maximum_gc_content=0.55, verbose=False)
print(result["interval_count"])
```

- **Example output:**

```text
9
```

## 6) `dnakit.variants.gatk_create_read_count_panel_of_normals()`

- **Purpose:** Create a panel of normals (PoN) for read-count denoising from the read counts of normal samples (HDF5).
- **Method:** Builds the fixed argument list of `gatk CreateReadCountPanelOfNormals` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk CreateReadCountPanelOfNormals --help`.
- **API:** `dnakit.variants.gatk_create_read_count_panel_of_normals(input_paths: 'Sequence[str]', output_path: 'str', *, annotated_intervals: 'str | None' = None, number_of_eigensamples: 'int | None' = None, minimum_interval_median_percentile: 'float | None' = None, maximum_zeros_in_sample_percentage: 'float | None' = None, maximum_zeros_in_interval_percentage: 'float | None' = None, do_impute_zeros: 'bool | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_paths`, `output_path`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2816.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`.
- **Note on the example:** the example data has no valid input for this tool (it needs data that cannot be made small: Illumina run folders, GTC files or flow-based reads, a trained or cluster model, or the Python gcnvkernel environment), so the example shows the input check that runs before GATK is started.

- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.variants.gatk_p00_copynumber_cli_functions import gatk_create_read_count_panel_of_normals

work = Path(tempfile.mkdtemp())
try:
    gatk_create_read_count_panel_of_normals(input_paths=["missing.dat"], output_path=work / "o", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] input_paths must be an existing file
```

## 7) `dnakit.variants.gatk_denoise_read_counts()`

- **Purpose:** Denoise read counts with a panel of normals (or annotated intervals) to produce standardized and denoised copy ratios.
- **Method:** Builds the fixed argument list of `gatk DenoiseReadCounts` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk DenoiseReadCounts --help`.
- **API:** `dnakit.variants.gatk_denoise_read_counts(input_path: 'str', standardized_copy_ratios: 'str', denoised_copy_ratios: 'str', *, count_panel_of_normals: 'str | None' = None, annotated_intervals: 'str | None' = None, number_of_eigensamples: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_path`, `standardized_copy_ratios`, `denoised_copy_ratios`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2817.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `rows`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants.gatk_p00_copynumber_cli_functions import (
    gatk_collect_read_counts,
    gatk_denoise_read_counts,
    gatk_preprocess_intervals,
)

work = Path(tempfile.mkdtemp())
gatk_preprocess_intervals("toy.fa", work / "bins.interval_list", bin_length=100, padding=0, verbose=False)
gatk_collect_read_counts(["toy.bam"], [work / "bins.interval_list"], work / "counts.tsv", output_format="TSV", verbose=False)
result = gatk_denoise_read_counts(work / "counts.tsv", work / "standardized.tsv", work / "denoised.tsv", verbose=False)
print(result["rows"][:2])
```

- **Example output:**

```text
[{'CONTIG': 'chr1', 'START': 1, 'END': 100, 'LOG2_COPY_RATIO': 0.160964}, {'CONTIG': 'chr1', 'START': 101, 'END': 200, 'LOG2_COPY_RATIO': -0.160964}]
```

## 8) `dnakit.variants.gatk_model_segments()`

- **Purpose:** Model segmented copy ratios from denoised copy ratios and segmented minor-allele fractions from allelic counts.
- **Method:** Builds the fixed argument list of `gatk ModelSegments` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk ModelSegments --help`.
- **API:** `dnakit.variants.gatk_model_segments(output_dir: 'str', output_prefix: 'str', *, denoised_copy_ratios: 'str | None' = None, allelic_counts: 'str | None' = None, normal_allelic_counts: 'str | None' = None, minimum_total_allele_count_case: 'int | None' = None, minimum_total_allele_count_normal: 'int | None' = None, genotyping_homozygous_log_ratio_threshold: 'float | None' = None, genotyping_base_error_rate: 'float | None' = None, maximum_number_of_segments_per_chromosome: 'int | None' = None, window_sizes: 'Sequence[int] | None' = None, number_of_changepoints_penalty_factor: 'float | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `output_dir`, `output_prefix`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2823.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `segment_files`, `segments`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants.gatk_p00_copynumber_cli_functions import (
    gatk_collect_read_counts,
    gatk_denoise_read_counts,
    gatk_model_segments,
    gatk_preprocess_intervals,
)

work = Path(tempfile.mkdtemp())
gatk_preprocess_intervals("toy.fa", work / "bins.interval_list", bin_length=100, padding=0, verbose=False)
gatk_collect_read_counts(["toy.bam"], [work / "bins.interval_list"], work / "counts.tsv", output_format="TSV", verbose=False)
gatk_denoise_read_counts(work / "counts.tsv", work / "standardized.tsv", work / "denoised.tsv", verbose=False)
result = gatk_model_segments(work / "segments", "S1", denoised_copy_ratios=work / "denoised.tsv", verbose=False)
print(len(result["segments"]), sorted(path.rsplit("/", 1)[1] for path in result["segment_files"]))
```

- **Example output:**

```text
2 ['S1.af.igv.seg', 'S1.cr.igv.seg', 'S1.cr.seg', 'S1.modelBegin.af.param', 'S1.modelBegin.cr.param', 'S1.modelBegin.seg', 'S1.modelFinal.af.param', 'S1.modelFinal.cr.param', 'S1.modelFinal.seg']
```

## 9) `dnakit.variants.gatk_call_copy_ratio_segments()`

- **Purpose:** Call copy-ratio segments as amplified (+), deleted (-) or copy-number neutral (0).
- **Method:** Builds the fixed argument list of `gatk CallCopyRatioSegments` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk CallCopyRatioSegments --help`.
- **API:** `dnakit.variants.gatk_call_copy_ratio_segments(input_path: 'str', output_path: 'str', *, neutral_segment_copy_ratio_lower_bound: 'float | None' = None, neutral_segment_copy_ratio_upper_bound: 'float | None' = None, outlier_neutral_segment_copy_ratio_z_score_threshold: 'float | None' = None, calling_copy_ratio_z_score_threshold: 'float | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_path`, `output_path`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2814.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `segments`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants.gatk_p00_copynumber_cli_functions import (
    gatk_call_copy_ratio_segments,
    gatk_collect_read_counts,
    gatk_denoise_read_counts,
    gatk_model_segments,
    gatk_preprocess_intervals,
)

work = Path(tempfile.mkdtemp())
gatk_preprocess_intervals("toy.fa", work / "bins.interval_list", bin_length=100, padding=0, verbose=False)
gatk_collect_read_counts(["toy.bam"], [work / "bins.interval_list"], work / "counts.tsv", output_format="TSV", verbose=False)
gatk_denoise_read_counts(work / "counts.tsv", work / "standardized.tsv", work / "denoised.tsv", verbose=False)
gatk_model_segments(work / "segments", "S1", denoised_copy_ratios=work / "denoised.tsv", verbose=False)
result = gatk_call_copy_ratio_segments(work / "segments" / "S1.cr.seg", work / "S1.called.seg", verbose=False)
print([row["CALL"] for row in result["segments"]])
```

- **Example output:**

```text
[0, 0]
```

## 10) `dnakit.variants.gatk_combine_segment_breakpoints()`

- **Purpose:** Combine the breakpoints of two segment files and annotate the resulting intervals with chosen columns from each file.
- **Method:** Builds the fixed argument list of `gatk CombineSegmentBreakpoints` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk CombineSegmentBreakpoints --help`.
- **API:** `dnakit.variants.gatk_combine_segment_breakpoints(segments: 'Sequence[str]', columns_of_interest: 'Sequence[str]', output_path: 'str', *, labels: 'Sequence[str] | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `segments`, `columns_of_interest`, `output_path`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2815.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `segments`.
- **Note on the example:** the example data has no valid input for this tool (it needs data that cannot be made small: Illumina run folders, GTC files or flow-based reads, a trained or cluster model, or the Python gcnvkernel environment), so the example shows the input check that runs before GATK is started.

- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.variants.gatk_p00_copynumber_cli_functions import gatk_combine_segment_breakpoints

work = Path(tempfile.mkdtemp())
try:
    gatk_combine_segment_breakpoints(segments=["missing.dat"], columns_of_interest=["x"], output_path=work / "o", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] segments must hold exactly two files.
```

## 11) `dnakit.variants.gatk_merge_annotated_regions()`

- **Purpose:** Merge annotated genomic regions. Without annotations_to_match, touching or overlapping regions are merged (MergeAnnotatedRegions); with annotations_to_match, regions within max_merge_distance whose listed annotation values are identical are merged (MergeAnnotatedRegionsByAnnotation).
- **Method:** Builds the fixed argument list of `gatk MergeAnnotatedRegions / MergeAnnotatedRegionsByAnnotation` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk MergeAnnotatedRegions / MergeAnnotatedRegionsByAnnotation --help`.
- **API:** `dnakit.variants.gatk_merge_annotated_regions(segments: 'str', output_path: 'str', *, reference: 'str | None' = None, annotations_to_match: 'Sequence[str] | None' = None, max_merge_distance: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `segments`, `output_path`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2821, G2822.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `regions`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants.gatk_p00_copynumber_cli_functions import (
    gatk_annotate_intervals,
    gatk_merge_annotated_regions,
    gatk_preprocess_intervals,
)

work = Path(tempfile.mkdtemp())
gatk_preprocess_intervals("toy.fa", work / "bins.interval_list", bin_length=100, padding=0, verbose=False)
gatk_annotate_intervals([work / "bins.interval_list"], "toy.fa", work / "annotated.tsv", verbose=False)
(work / "annotated.seg").write_bytes((work / "annotated.tsv").read_bytes())
plain = gatk_merge_annotated_regions(work / "annotated.seg", work / "plain.seg", reference="toy.fa", verbose=False)
by_annotation = gatk_merge_annotated_regions(work / "annotated.seg", work / "by_gc.seg", reference="toy.fa", annotations_to_match=["GC_CONTENT"], verbose=False)
print(len(plain["regions"]), by_annotation["tool"])
```

- **Example output:**

```text
10 MergeAnnotatedRegionsByAnnotation
```

## 12) `dnakit.variants.gatk_plot_denoised_copy_ratios()`

- **Purpose:** Create plots of standardized and denoised copy ratios and their median absolute deviations.
- **Method:** Builds the fixed argument list of `gatk PlotDenoisedCopyRatios` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk PlotDenoisedCopyRatios --help`. The plots need R (Rscript) with the packages of the GATK conda environment; without it GATK stops with an error.
- **API:** `dnakit.variants.gatk_plot_denoised_copy_ratios(standardized_copy_ratios: 'str', denoised_copy_ratios: 'str', sequence_dictionary: 'str', output_dir: 'str', output_prefix: 'str', *, minimum_contig_length: 'int | None' = None, maximum_copy_ratio: 'float | None' = None, point_size_copy_ratio: 'float | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `standardized_copy_ratios`, `denoised_copy_ratios`, `sequence_dictionary`, `output_dir`, `output_prefix`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2824.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `plot_files`.
- **Note on the example:** the example data has no valid input for this tool (it needs data that cannot be made small: Illumina run folders, GTC files or flow-based reads, a trained or cluster model, or the Python gcnvkernel environment), so the example shows the input check that runs before GATK is started.

- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.variants.gatk_p00_copynumber_cli_functions import gatk_plot_denoised_copy_ratios

work = Path(tempfile.mkdtemp())
try:
    gatk_plot_denoised_copy_ratios(standardized_copy_ratios="missing.dat", denoised_copy_ratios="missing.dat", sequence_dictionary="missing.dat", output_dir=work / "d", output_prefix="x", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] standardized_copy_ratios must be an existing file
```

## 13) `dnakit.variants.gatk_plot_modeled_segments()`

- **Purpose:** Create plots of denoised and segmented copy-ratio and minor-allele-fraction estimates.
- **Method:** Builds the fixed argument list of `gatk PlotModeledSegments` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk PlotModeledSegments --help`. The plots need R (Rscript) with the packages of the GATK conda environment; without it GATK stops with an error.
- **API:** `dnakit.variants.gatk_plot_modeled_segments(segments: 'str', sequence_dictionary: 'str', output_dir: 'str', output_prefix: 'str', *, denoised_copy_ratios: 'str | None' = None, allelic_counts: 'str | None' = None, minimum_contig_length: 'int | None' = None, maximum_copy_ratio: 'float | None' = None, point_size_copy_ratio: 'float | None' = None, point_size_allele_fraction: 'float | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `segments`, `sequence_dictionary`, `output_dir`, `output_prefix`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2825.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `plot_files`.
- **Note on the example:** the example data has no valid input for this tool (it needs data that cannot be made small: Illumina run folders, GTC files or flow-based reads, a trained or cluster model, or the Python gcnvkernel environment), so the example shows the input check that runs before GATK is started.

- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.variants.gatk_p00_copynumber_cli_functions import gatk_plot_modeled_segments

work = Path(tempfile.mkdtemp())
try:
    gatk_plot_modeled_segments(segments="missing.dat", sequence_dictionary="missing.dat", output_dir=work / "d", output_prefix="x", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] segments must be an existing file
```

## 14) `dnakit.variants.gatk_tag_germline_events()`

- **Purpose:** Tag segments of a called tumor segment file that overlap germline events found in the called matched normal.
- **Method:** Builds the fixed argument list of `gatk TagGermlineEvents` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk TagGermlineEvents --help`.
- **API:** `dnakit.variants.gatk_tag_germline_events(segments: 'str', called_matched_normal_seg_file: 'str', reference: 'str', output_path: 'str', *, endpoint_padding: 'int | None' = None, reciprocal_threshold: 'float | None' = None, input_call_header: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `segments`, `called_matched_normal_seg_file`, `reference`, `output_path`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2827.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `segments`.
- **Note on the example:** the example data has no valid input for this tool (it needs data that cannot be made small: Illumina run folders, GTC files or flow-based reads, a trained or cluster model, or the Python gcnvkernel environment), so the example shows the input check that runs before GATK is started.

- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.variants.gatk_p00_copynumber_cli_functions import gatk_tag_germline_events

work = Path(tempfile.mkdtemp())
try:
    gatk_tag_germline_events(segments="missing.dat", called_matched_normal_seg_file="missing.dat", reference="missing.dat", output_path=work / "o", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] segments must be an existing file
```

## 15) `dnakit.variants.gatk_determine_germline_contig_ploidy()`

- **Purpose:** Determine the baseline contig ploidy of germline samples from read counts (cohort mode builds a ploidy model, case mode applies one).
- **Method:** Builds the fixed argument list of `gatk DetermineGermlineContigPloidy` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk DetermineGermlineContigPloidy --help`. This tool runs GATK's Python module gcnvkernel, so the gatk launcher must run inside the GATK conda environment (gatkcondaenv.yml of the release zip).
- **API:** `dnakit.variants.gatk_determine_germline_contig_ploidy(input_paths: 'Sequence[str]', output_dir: 'str', output_prefix: 'str', *, contig_ploidy_priors: 'str | None' = None, model: 'str | None' = None, intervals: 'Sequence[str] | None' = None, mapping_error_rate: 'float | None' = None, interval_merging_rule: 'str | None' = 'OVERLAPPING_ONLY', options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_paths`, `output_dir`, `output_prefix`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2818.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `calls_directories`, `model_directories`.
- **Note on the example:** the example data has no valid input for this tool (it needs data that cannot be made small: Illumina run folders, GTC files or flow-based reads, a trained or cluster model, or the Python gcnvkernel environment), so the example shows the input check that runs before GATK is started.

- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.variants.gatk_p00_copynumber_cli_functions import gatk_determine_germline_contig_ploidy

work = Path(tempfile.mkdtemp())
try:
    gatk_determine_germline_contig_ploidy(input_paths=["missing.dat"], output_dir=work / "d", output_prefix="x", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] input_paths must be an existing file
```

## 16) `dnakit.variants.gatk_germline_cnv_caller()`

- **Purpose:** Call copy-number variants in germline samples from read counts and the contig ploidy calls (COHORT mode learns a denoising model, CASE mode applies one).
- **Method:** Builds the fixed argument list of `gatk GermlineCNVCaller` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk GermlineCNVCaller --help`. This tool runs GATK's Python module gcnvkernel, so the gatk launcher must run inside the GATK conda environment (gatkcondaenv.yml of the release zip).
- **API:** `dnakit.variants.gatk_germline_cnv_caller(run_mode: 'str', input_paths: 'Sequence[str]', contig_ploidy_calls: 'str', output_dir: 'str', output_prefix: 'str', *, model: 'str | None' = None, annotated_intervals: 'str | None' = None, intervals: 'Sequence[str] | None' = None, max_copy_number: 'int | None' = None, p_alt: 'float | None' = None, p_active: 'float | None' = None, interval_merging_rule: 'str | None' = 'OVERLAPPING_ONLY', options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `run_mode`, `input_paths`, `contig_ploidy_calls`, `output_dir`, `output_prefix`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2820.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `calls_directories`, `model_directories`.
- **Note on the example:** the example data has no valid input for this tool (it needs data that cannot be made small: Illumina run folders, GTC files or flow-based reads, a trained or cluster model, or the Python gcnvkernel environment), so the example shows the input check that runs before GATK is started.

- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.variants.gatk_p00_copynumber_cli_functions import gatk_germline_cnv_caller

work = Path(tempfile.mkdtemp())
try:
    gatk_germline_cnv_caller(run_mode="COHORT", input_paths=["missing.dat"], contig_ploidy_calls="missing_dir", output_dir=work / "d", output_prefix="x", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] input_paths must be an existing file
```

## 17) `dnakit.variants.gatk_postprocess_germline_cnv_calls()`

- **Purpose:** Postprocess GermlineCNVCaller output of one sample into an intervals VCF, a segments VCF and denoised copy ratios.
- **Method:** Builds the fixed argument list of `gatk PostprocessGermlineCNVCalls` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk PostprocessGermlineCNVCalls --help`. This tool runs GATK's Python module gcnvkernel, so the gatk launcher must run inside the GATK conda environment (gatkcondaenv.yml of the release zip).
- **API:** `dnakit.variants.gatk_postprocess_germline_cnv_calls(calls_shard_paths: 'Sequence[str]', model_shard_paths: 'Sequence[str]', contig_ploidy_calls: 'str', output_genotyped_intervals: 'str', output_genotyped_segments: 'str', output_denoised_copy_ratios: 'str', *, sample_index: 'int | None' = None, allosomal_contigs: 'Sequence[str] | None' = None, autosomal_ref_copy_number: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `calls_shard_paths`, `model_shard_paths`, `contig_ploidy_calls`, `output_genotyped_intervals`, `output_genotyped_segments`, `output_denoised_copy_ratios`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2826.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `segment_record_count`.
- **Note on the example:** the example data has no valid input for this tool (it needs data that cannot be made small: Illumina run folders, GTC files or flow-based reads, a trained or cluster model, or the Python gcnvkernel environment), so the example shows the input check that runs before GATK is started.

- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.variants.gatk_p00_copynumber_cli_functions import gatk_postprocess_germline_cnv_calls

work = Path(tempfile.mkdtemp())
try:
    gatk_postprocess_germline_cnv_calls(calls_shard_paths=["missing_dir"], model_shard_paths=["missing_dir"], contig_ploidy_calls="missing_dir", output_genotyped_intervals=work / "o", output_genotyped_segments=work / "o", output_denoised_copy_ratios=work / "o", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] calls_shard_paths must be an existing directory
```

## 18) `dnakit.variants.gatk_ase_read_counter()`

- **Purpose:** Generate a table of filtered base counts at heterozygous sites for allele-specific expression.
- **Method:** Builds the fixed argument list of `gatk ASEReadCounter` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk ASEReadCounter --help`.
- **API:** `dnakit.variants.gatk_ase_read_counter(input_paths: 'Sequence[str]', variant_paths: 'Sequence[str]', output_path: 'str', *, reference: 'str | None' = None, intervals: 'Sequence[str] | None' = None, output_format: 'str | None' = None, min_base_quality: 'int | None' = None, min_mapping_quality: 'int | None' = None, min_depth_of_non_filtered_base: 'int | None' = None, max_depth_per_sample: 'int | None' = None, count_overlap_reads_handling: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_paths`, `variant_paths`, `output_path`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2828.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `rows`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants.gatk_p00_copynumber_cli_functions import gatk_ase_read_counter

work = Path(tempfile.mkdtemp())
result = gatk_ase_read_counter(["toy.bam"], ["toy.vcf.gz"], work / "ase.tsv", reference="toy.fa", verbose=False)
print(result["rows"][2]["position"], result["rows"][2]["refCount"], result["rows"][2]["altCount"])
```

- **Example output:**

```text
160 1 1
```

## 19) `dnakit.variants.gatk_analyze_saturation_mutagenesis()`

- **Purpose:** Process reads from a MITESeq or other saturation mutagenesis experiment into per-variant counts of codon and amino-acid changes.
- **Method:** Builds the fixed argument list of `gatk AnalyzeSaturationMutagenesis` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk AnalyzeSaturationMutagenesis --help`.
- **API:** `dnakit.variants.gatk_analyze_saturation_mutagenesis(input_paths: 'Sequence[str]', reference: 'str', orf: 'str', output_file_prefix: 'str', *, min_q: 'int | None' = None, min_mapq: 'int | None' = None, min_length: 'int | None' = None, min_flanking_length: 'int | None' = None, min_variant_obs: 'int | None' = None, paired_mode: 'bool | None' = None, codon_translation: 'str | None' = None, write_rejected_reads: 'bool | None' = None, find_large_deletions: 'bool | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_paths`, `reference`, `orf`, `output_file_prefix`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2829.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `output_files`.
- **Note on the example:** the example data has no valid input for this tool (it needs data that cannot be made small: Illumina run folders, GTC files or flow-based reads, a trained or cluster model, or the Python gcnvkernel environment), so the example shows the input check that runs before GATK is started.

- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.variants.gatk_p00_copynumber_cli_functions import gatk_analyze_saturation_mutagenesis

work = Path(tempfile.mkdtemp())
try:
    gatk_analyze_saturation_mutagenesis(input_paths=["missing.dat"], reference="missing.dat", orf="x", output_file_prefix=work / "o", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] input_paths must be an existing file
```

## 20) `dnakit.variants.gatk_callable_loci()`

- **Purpose:** Determine the callable status of every locus (CALLABLE, LOW_COVERAGE, NO_COVERAGE, POOR_MAPPING_QUALITY, EXCESSIVE_COVERAGE, REF_N).
- **Method:** Builds the fixed argument list of `gatk CallableLoci` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk CallableLoci --help`.
- **API:** `dnakit.variants.gatk_callable_loci(input_paths: 'Sequence[str]', reference: 'str', output_path: 'str', summary_path: 'str', *, intervals: 'Sequence[str] | None' = None, min_depth: 'int | None' = None, max_depth: 'int | None' = None, min_mapping_quality: 'int | None' = None, min_base_quality: 'int | None' = None, min_depth_for_low_mapq: 'int | None' = None, max_low_mapq: 'int | None' = None, max_fraction_of_reads_with_low_mapq: 'float | None' = None, output_format: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_paths`, `reference`, `output_path`, `summary_path`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2830.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `summary`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants.gatk_p00_copynumber_cli_functions import gatk_callable_loci

work = Path(tempfile.mkdtemp())
result = gatk_callable_loci(["toy.bam"], "toy.fa", work / "callable.bed", work / "callable.txt", min_depth=2, verbose=False)
print(result["summary"])
```

- **Example output:**

```text
{'REF_N': 0, 'CALLABLE': 720, 'NO_COVERAGE': 120, 'LOW_COVERAGE': 160, 'EXCESSIVE_COVERAGE': 0, 'POOR_MAPPING_QUALITY': 0}
```

## 21) `dnakit.variants.gatk_collect_f1r2_counts()`

- **Purpose:** Collect F1R2 read counts per sequence context for the Mutect2 orientation-bias mixture model filter.
- **Method:** Builds the fixed argument list of `gatk CollectF1R2Counts` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk CollectF1R2Counts --help`.
- **API:** `dnakit.variants.gatk_collect_f1r2_counts(input_paths: 'Sequence[str]', reference: 'str', output_path: 'str', *, intervals: 'Sequence[str] | None' = None, f1r2_min_bq: 'int | None' = None, f1r2_median_mq: 'int | None' = None, f1r2_max_depth: 'int | None' = None, max_depth_per_sample: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_paths`, `reference`, `output_path`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2832.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants.gatk_p00_copynumber_cli_functions import gatk_collect_f1r2_counts

work = Path(tempfile.mkdtemp())
result = gatk_collect_f1r2_counts(["toy.bam"], "toy.fa", work / "f1r2.tar.gz", verbose=False)
print(result["status"], (work / "f1r2.tar.gz").stat().st_size > 0)
```

- **Example output:**

```text
completed True
```

## 22) `dnakit.variants.gatk_count_reads_or_bases()`

- **Purpose:** Count the reads (unit='reads', CountReads) or the bases (unit='bases', CountBases) of a SAM/BAM/CRAM file after GATK's default read filters.
- **Method:** Builds the fixed argument list of `gatk CountReads / CountBases` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk CountReads / CountBases --help`.
- **API:** `dnakit.variants.gatk_count_reads_or_bases(input_paths: 'Sequence[str]', output_path: 'str', *, unit: 'str | None' = None, intervals: 'Sequence[str] | None' = None, reference: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_paths`, `output_path`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2834, G2835.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `count`, `unit`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants.gatk_p00_copynumber_cli_functions import gatk_count_reads_or_bases

work = Path(tempfile.mkdtemp())
reads = gatk_count_reads_or_bases(["toy.bam"], work / "reads.txt", verbose=False)
bases = gatk_count_reads_or_bases(["toy.bam"], work / "bases.txt", unit="bases", verbose=False)
print(reads["count"], bases["count"])
```

- **Example output:**

```text
38 1900
```

## 23) `dnakit.variants.gatk_depth_of_coverage()`

- **Purpose:** Generate coverage summary information for reads data: depth per locus, per interval, per gene and per sample, with coverage histograms.
- **Method:** Builds the fixed argument list of `gatk DepthOfCoverage` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk DepthOfCoverage --help`.
- **API:** `dnakit.variants.gatk_depth_of_coverage(input_paths: 'Sequence[str]', reference: 'str', intervals: 'Sequence[str]', output_prefix: 'str', *, partition_types: 'Sequence[str] | None' = None, summary_coverage_thresholds: 'Sequence[int] | None' = None, min_base_quality: 'int | None' = None, count_type: 'str | None' = None, output_format: 'str | None' = None, calculate_coverage_over_genes: 'str | None' = None, omit_depth_output_at_each_base: 'bool | None' = None, omit_interval_statistics: 'bool | None' = None, omit_locus_table: 'bool | None' = None, omit_per_sample_statistics: 'bool | None' = None, print_base_counts: 'bool | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_paths`, `reference`, `intervals`, `output_prefix`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2836.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `output_files`, `sample_summary`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants.gatk_p00_copynumber_cli_functions import gatk_depth_of_coverage

work = Path(tempfile.mkdtemp())
result = gatk_depth_of_coverage(["toy.bam"], "toy.fa", ["toy.interval_list"], work / "coverage", verbose=False)
row = result["sample_summary"][0]
print(row["sample_id"], row["total"], row["mean"])
```

- **Example output:**

```text
S1 860 2.46
```

## 24) `dnakit.variants.gatk_gene_expression_evaluation()`

- **Purpose:** Evaluate gene expression from RNA-seq reads aligned to the genome: read counts per gene, exon or other feature of a GFF.
- **Method:** Builds the fixed argument list of `gatk GeneExpressionEvaluation` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk GeneExpressionEvaluation --help`.
- **API:** `dnakit.variants.gatk_gene_expression_evaluation(input_paths: 'Sequence[str]', gff_file: 'str', output_path: 'str', *, reference: 'str | None' = None, intervals: 'Sequence[str] | None' = None, grouping_types: 'Sequence[str] | None' = None, overlap_types: 'Sequence[str] | None' = None, feature_label_key: 'str | None' = None, multi_map_method: 'str | None' = None, multi_overlap_method: 'str | None' = None, read_strands: 'str | None' = None, unspliced: 'bool | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_paths`, `gff_file`, `output_path`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2837.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `rows`.
- **Note on the example:** the example data has no valid input for this tool (it needs data that cannot be made small: Illumina run folders, GTC files or flow-based reads, a trained or cluster model, or the Python gcnvkernel environment), so the example shows the input check that runs before GATK is started.

- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.variants.gatk_p00_copynumber_cli_functions import gatk_gene_expression_evaluation

work = Path(tempfile.mkdtemp())
try:
    gatk_gene_expression_evaluation(input_paths=["missing.dat"], gff_file="missing.dat", output_path=work / "o", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] input_paths must be an existing file
```

## 25) `dnakit.variants.gatk_get_normal_artifact_data()`

- **Purpose:** Collect data for training the normal-artifact filter of Mutect2 (allele counts of normal samples at candidate sites).
- **Method:** Builds the fixed argument list of `gatk GetNormalArtifactData` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk GetNormalArtifactData --help`.
- **API:** `dnakit.variants.gatk_get_normal_artifact_data(input_paths: 'Sequence[str]', normal_samples: 'Sequence[str]', reference: 'str', output_path: 'str', *, intervals: 'Sequence[str] | None' = None, error_prob: 'float | None' = None, max_depth_per_sample: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_paths`, `normal_samples`, `reference`, `output_path`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2838.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants.gatk_p00_copynumber_cli_functions import gatk_get_normal_artifact_data

work = Path(tempfile.mkdtemp())
result = gatk_get_normal_artifact_data(["toy.bam"], ["S1"], "toy.fa", work / "normal_artifact.table", verbose=False)
print(result["status"], (work / "normal_artifact.table").stat().st_size > 0)
```

- **Example output:**

```text
completed True
```

## 26) `dnakit.variants.gatk_get_pileup_summaries()`

- **Purpose:** Tabulate pileup metrics (ref, alt and other-alt counts at common biallelic sites) for inferring cross-sample contamination.
- **Method:** Builds the fixed argument list of `gatk GetPileupSummaries` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk GetPileupSummaries --help`.
- **API:** `dnakit.variants.gatk_get_pileup_summaries(input_paths: 'Sequence[str]', variant: 'str', intervals: 'Sequence[str]', output_path: 'str', *, minimum_population_allele_frequency: 'float | None' = None, maximum_population_allele_frequency: 'float | None' = None, max_depth_per_sample: 'int | None' = None, reference: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_paths`, `variant`, `intervals`, `output_path`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2839.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `rows`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants.gatk_p00_copynumber_cli_functions import gatk_get_pileup_summaries

work = Path(tempfile.mkdtemp())
result = gatk_get_pileup_summaries(["toy.bam"], "toy.vcf.gz", ["toy.vcf.gz"], work / "pileup.table", maximum_population_allele_frequency=0.9, verbose=False)
print(len(result["rows"]), result["rows"][2])
```

- **Example output:**

```text
15 {'contig': 'chr1', 'position': 160, 'ref_count': 1, 'alt_count': 1, 'other_alt_count': 0, 'allele_frequency': 0.4}
```

## 27) `dnakit.variants.gatk_local_assembler()`

- **Purpose:** Assemble small regions locally to discover structural variants (GATK BETA tool) and write the scaffolds or the assembly graph.
- **Method:** Builds the fixed argument list of `gatk LocalAssembler` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk LocalAssembler --help`. GATK marks this tool BETA (work in progress).
- **API:** `dnakit.variants.gatk_local_assembler(input_paths: 'Sequence[str]', intervals: 'Sequence[str]', assembly_name: 'str', *, fasta_file: 'str | None' = None, gfa_file: 'str | None' = None, reference: 'str | None' = None, min_sv_size: 'int | None' = None, no_scaffolding: 'bool | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_paths`, `intervals`, `assembly_name`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2840.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants.gatk_p00_copynumber_cli_functions import gatk_local_assembler

work = Path(tempfile.mkdtemp())
result = gatk_local_assembler(["toy.bam"], ["chr1:1-600"], "demo", fasta_file=work / "scaffolds.fa", gfa_file=work / "graph.gfa", reference="toy.fa", verbose=False)
print(result["status"], sorted(result["outputs"]))
```

- **Example output:**

```text
completed ['fasta_file', 'gfa_file']
```

## 28) `dnakit.variants.gatk_analyze_covariates()`

- **Purpose:** Evaluate and compare base quality score recalibration (BQSR) tables and plot them.
- **Method:** Builds the fixed argument list of `gatk AnalyzeCovariates` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk AnalyzeCovariates --help`. The plots need R (Rscript) with the packages of the GATK conda environment; without it GATK stops with an error.
- **API:** `dnakit.variants.gatk_analyze_covariates(*, before_report_file: 'str | None' = None, after_report_file: 'str | None' = None, bqsr_recal_file: 'str | None' = None, plots_report_file: 'str | None' = None, intermediate_csv_file: 'str | None' = None, ignore_last_modification_times: 'bool | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: . Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2841.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`.
- **Note on the example:** the example data has no valid input for this tool (it needs data that cannot be made small: Illumina run folders, GTC files or flow-based reads, a trained or cluster model, or the Python gcnvkernel environment), so the example shows the input check that runs before GATK is started.

- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.variants.gatk_p00_copynumber_cli_functions import gatk_analyze_covariates

work = Path(tempfile.mkdtemp())
try:
    gatk_analyze_covariates(verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] Give at least one of before_report_file, after_report_file or bqsr_recal_file.
```
