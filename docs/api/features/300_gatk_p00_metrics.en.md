# GATK command-line facades: Illumina runs, read metrics, fingerprints and file checks

Typed facades over the `gatk` command line (GATK 4.6.2.0 with the bundled Picard 3.4.0) for processing Illumina runs (basecalling metrics, barcode extraction, conversion to FASTQ or unmapped BAM, adapter marking), collecting alignment and sequencing metrics (alignment summary, hybrid selection, WGS, RNA-seq, RRBS, artifacts, UMI, variant calling), fingerprinting, and checking, viewing and comparing SAM/BAM/CRAM files. Each function builds a fixed `gatk` command line without a shell, writes only to new paths and returns a JSON-serialisable dictionary; GATK itself is licensed BSD-3-Clause (Picard MIT) and is not bundled. The functions live in `dnakit.sequencing_qc.gatk_p00_metrics_cli_functions` and are exported from `dnakit.sequencing_qc`; the examples import from the module so that they run either way.

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

## 1) `dnakit.sequencing_qc.gatk_collect_illumina_basecalling_metrics()`

- **Purpose:** Collect Illumina basecalling metrics (clusters, PF clusters, per-barcode counts) for one lane of a sequencing run.
- **Method:** Builds the fixed argument list of `gatk CollectIlluminaBasecallingMetrics` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk CollectIlluminaBasecallingMetrics --help`. Needs the BaseCalls directory of an Illumina run (BCL/CBCL, filter, locs or clocs files).
- **API:** `dnakit.sequencing_qc.gatk_collect_illumina_basecalling_metrics(basecalls_dir: 'str', lane: 'int', read_structure: 'str', *, output_path: 'str | None' = None, barcodes_dir: 'str | None' = None, barcode_file: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `basecalls_dir`, `lane`, `read_structure`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2807.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `metrics`, `histogram`.
- **Note on the example:** the example data has no valid input for this tool (it needs data that cannot be made small: Illumina run folders, GTC files or flow-based reads, a trained or cluster model, or the Python gcnvkernel environment), so the example shows the input check that runs before GATK is started.

- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_collect_illumina_basecalling_metrics

work = Path(tempfile.mkdtemp())
try:
    gatk_collect_illumina_basecalling_metrics(basecalls_dir="missing_dir", lane=1, read_structure="x", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] basecalls_dir must be an existing directory
```

## 2) `dnakit.sequencing_qc.gatk_collect_illumina_lane_metrics()`

- **Purpose:** Collect Illumina lane and phasing metrics (cluster density, PF clusters, phasing and prephasing) from a run directory.
- **Method:** Builds the fixed argument list of `gatk CollectIlluminaLaneMetrics` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk CollectIlluminaLaneMetrics --help`. Needs the InterOp files of an Illumina run.
- **API:** `dnakit.sequencing_qc.gatk_collect_illumina_lane_metrics(run_directory: 'str', output_directory: 'str', output_prefix: 'str', *, read_structure: 'str | None' = None, file_extension: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `run_directory`, `output_directory`, `output_prefix`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2808.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `metrics_files`, `outputs`.
- **Note on the example:** the example data has no valid input for this tool (it needs data that cannot be made small: Illumina run folders, GTC files or flow-based reads, a trained or cluster model, or the Python gcnvkernel environment), so the example shows the input check that runs before GATK is started.

- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_collect_illumina_lane_metrics

work = Path(tempfile.mkdtemp())
try:
    gatk_collect_illumina_lane_metrics(run_directory="missing_dir", output_directory="missing_dir", output_prefix="x", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] run_directory must be an existing directory
```

## 3) `dnakit.sequencing_qc.gatk_extract_illumina_barcodes()`

- **Purpose:** Determine the barcode of every read of an Illumina lane and write the *_barcode.txt files and per-barcode metrics.
- **Method:** Builds the fixed argument list of `gatk ExtractIlluminaBarcodes` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk ExtractIlluminaBarcodes --help`. Needs the BaseCalls directory of an Illumina run.
- **API:** `dnakit.sequencing_qc.gatk_extract_illumina_barcodes(basecalls_dir: 'str', lanes: 'Sequence[int]', read_structure: 'str', *, barcodes: 'Sequence[str] | None' = None, barcode_file: 'str | None' = None, metrics_file: 'str | None' = None, output_dir: 'str | None' = None, max_mismatches: 'int | None' = None, min_mismatch_delta: 'int | None' = None, max_no_calls: 'int | None' = None, minimum_base_quality: 'int | None' = None, minimum_quality: 'int | None' = None, distance_mode: 'str | None' = None, compress_outputs: 'bool | None' = None, num_processors: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `basecalls_dir`, `lanes`, `read_structure`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2809.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `metrics`, `histogram`.
- **Note on the example:** the example data has no valid input for this tool (it needs data that cannot be made small: Illumina run folders, GTC files or flow-based reads, a trained or cluster model, or the Python gcnvkernel environment), so the example shows the input check that runs before GATK is started.

- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_extract_illumina_barcodes

work = Path(tempfile.mkdtemp())
try:
    gatk_extract_illumina_barcodes(basecalls_dir="missing_dir", lanes=[1], read_structure="x", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] Give exactly one of barcodes or barcode_file.
```

## 4) `dnakit.sequencing_qc.gatk_illumina_basecalls_to_reads()`

- **Purpose:** Convert Illumina basecalls of one run to FASTQ files (output_format='fastq') or to unmapped SAM/BAM/CRAM files (output_format='sam').
- **Method:** Builds the fixed argument list of `gatk IlluminaBasecallsToFastq / IlluminaBasecallsToSam` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk IlluminaBasecallsToFastq / IlluminaBasecallsToSam --help`. Needs the BaseCalls directory of an Illumina run.
- **API:** `dnakit.sequencing_qc.gatk_illumina_basecalls_to_reads(basecalls_dir: 'str', lanes: 'Sequence[int]', read_structure: 'str', run_barcode: 'str', parameters_file: 'str', output_format: 'str', *, output_prefix: 'str | None' = None, sequencing_center: 'str | None' = None, metrics_file: 'str | None' = None, barcodes_dir: 'str | None' = None, include_non_pf_reads: 'bool | None' = None, max_mismatches: 'int | None' = None, min_mismatch_delta: 'int | None' = None, max_no_calls: 'int | None' = None, minimum_quality: 'int | None' = None, first_tile: 'int | None' = None, tile_limit: 'int | None' = None, num_processors: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `basecalls_dir`, `lanes`, `read_structure`, `run_barcode`, `parameters_file`, `output_format`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2810, G2811.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `output_format`, `metrics`.
- **Note on the example:** the example data has no valid input for this tool (it needs data that cannot be made small: Illumina run folders, GTC files or flow-based reads, a trained or cluster model, or the Python gcnvkernel environment), so the example shows the input check that runs before GATK is started.

- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_illumina_basecalls_to_reads

work = Path(tempfile.mkdtemp())
try:
    gatk_illumina_basecalls_to_reads(basecalls_dir="missing_dir", lanes=[1], read_structure="x", run_barcode="x", parameters_file="missing.dat", output_format="fastq", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] basecalls_dir must be an existing directory
```

## 5) `dnakit.sequencing_qc.gatk_mark_illumina_adapters()`

- **Purpose:** Find Illumina adapter sequences in reads and rewrite the file with XT adapter-trimming tags, plus an adapter-clipping histogram.
- **Method:** Builds the fixed argument list of `gatk MarkIlluminaAdapters` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk MarkIlluminaAdapters --help`.
- **API:** `dnakit.sequencing_qc.gatk_mark_illumina_adapters(input_path: 'str', metrics_file: 'str', *, output_path: 'str | None' = None, adapters: 'Sequence[str] | None' = None, five_prime_adapter: 'str | None' = None, three_prime_adapter: 'str | None' = None, max_error_rate_se: 'float | None' = None, max_error_rate_pe: 'float | None' = None, min_match_bases_se: 'int | None' = None, min_match_bases_pe: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_path`, `metrics_file`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2812.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `metrics`, `histogram`.
- **Note on the example:** the example data has no valid input for this tool (it needs data that cannot be made small: Illumina run folders, GTC files or flow-based reads, a trained or cluster model, or the Python gcnvkernel environment), so the example shows the input check that runs before GATK is started.

- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_mark_illumina_adapters

work = Path(tempfile.mkdtemp())
try:
    gatk_mark_illumina_adapters(input_path="missing.dat", metrics_file=work / "o", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] input_path must be an existing file
```

## 6) `dnakit.sequencing_qc.gatk_collect_alignment_summary_metrics()`

- **Purpose:** Produce a summary of alignment metrics (reads, PF reads, aligned bases, mismatch and indel rates, chimeras, adapters) from a SAM/BAM/CRAM file.
- **Method:** Builds the fixed argument list of `gatk CollectAlignmentSummaryMetrics` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk CollectAlignmentSummaryMetrics --help`.
- **API:** `dnakit.sequencing_qc.gatk_collect_alignment_summary_metrics(input_path: 'str', output_path: 'str', *, reference: 'str | None' = None, metric_accumulation_level: 'Sequence[str] | None' = None, assume_sorted: 'bool | None' = None, is_bisulfite_sequenced: 'bool | None' = None, max_insert_size: 'int | None' = None, expected_pair_orientations: 'Sequence[str] | None' = None, adapter_sequences: 'Sequence[str] | None' = None, histogram_file: 'str | None' = None, stop_after: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_path`, `output_path`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2850.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `metrics`, `histogram`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_collect_alignment_summary_metrics

work = Path(tempfile.mkdtemp())
result = gatk_collect_alignment_summary_metrics("toy.bam", work / "summary.txt", reference="toy.fa", verbose=False)
row = next(r for r in result["metrics"] if r["CATEGORY"] == "PAIR")
print(result["status"], row["TOTAL_READS"], row["PF_READS_ALIGNED"], row["PCT_PF_READS_ALIGNED"])
```

- **Example output:**

```text
completed 38 38 1
```

## 7) `dnakit.sequencing_qc.gatk_collect_base_distribution_by_cycle()`

- **Purpose:** Chart the nucleotide distribution per sequencing cycle of a SAM/BAM/CRAM file and write the per-cycle percentages.
- **Method:** Builds the fixed argument list of `gatk CollectBaseDistributionByCycle` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk CollectBaseDistributionByCycle --help`. The chart output needs R (Rscript) with the packages of the GATK conda environment; without it GATK stops with an error.
- **API:** `dnakit.sequencing_qc.gatk_collect_base_distribution_by_cycle(input_path: 'str', output_path: 'str', chart_output: 'str', *, aligned_reads_only: 'bool | None' = None, pf_reads_only: 'bool | None' = None, reference: 'str | None' = None, assume_sorted: 'bool | None' = None, stop_after: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_path`, `output_path`, `chart_output`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2852.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `metrics`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_collect_base_distribution_by_cycle

work = Path(tempfile.mkdtemp())
result = gatk_collect_base_distribution_by_cycle("toy.bam", work / "base.txt", work / "base.pdf", verbose=False)
print(result["status"], (work / "base.pdf").is_file(), len(result["metrics"]) > 0)
```

- **Example output:**

```text
completed True True
```

## 8) `dnakit.sequencing_qc.gatk_mean_quality_by_cycle()`

- **Purpose:** Collect the mean base quality for every sequencing cycle and chart it.
- **Method:** Builds the fixed argument list of `gatk MeanQualityByCycle` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk MeanQualityByCycle --help`. The chart output needs R (Rscript) with the packages of the GATK conda environment; without it GATK stops with an error.
- **API:** `dnakit.sequencing_qc.gatk_mean_quality_by_cycle(input_path: 'str', output_path: 'str', chart_output: 'str', *, aligned_reads_only: 'bool | None' = None, pf_reads_only: 'bool | None' = None, reference: 'str | None' = None, assume_sorted: 'bool | None' = None, stop_after: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_path`, `output_path`, `chart_output`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2883.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `metrics`, `histogram`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_mean_quality_by_cycle

work = Path(tempfile.mkdtemp())
result = gatk_mean_quality_by_cycle("toy.bam", work / "mean.txt", work / "mean.pdf", verbose=False)
print(result["status"], (work / "mean.pdf").is_file(), len(result["metrics"]) > 0)
```

- **Example output:**

```text
completed True False
```

## 9) `dnakit.sequencing_qc.gatk_quality_score_distribution()`

- **Purpose:** Chart the distribution of base quality scores (number of bases per quality value) of a SAM/BAM/CRAM file.
- **Method:** Builds the fixed argument list of `gatk QualityScoreDistribution` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk QualityScoreDistribution --help`. The chart output needs R (Rscript) with the packages of the GATK conda environment; without it GATK stops with an error.
- **API:** `dnakit.sequencing_qc.gatk_quality_score_distribution(input_path: 'str', output_path: 'str', chart_output: 'str', *, aligned_reads_only: 'bool | None' = None, pf_reads_only: 'bool | None' = None, include_no_calls: 'bool | None' = None, reference: 'str | None' = None, assume_sorted: 'bool | None' = None, stop_after: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_path`, `output_path`, `chart_output`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2884.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `metrics`, `histogram`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_quality_score_distribution

work = Path(tempfile.mkdtemp())
result = gatk_quality_score_distribution("toy.bam", work / "qual.txt", work / "qual.pdf", verbose=False)
print(result["status"], (work / "qual.pdf").is_file(), len(result["metrics"]) > 0)
```

- **Example output:**

```text
completed True False
```

## 10) `dnakit.sequencing_qc.gatk_collect_hiseqx_pf_fail_metrics()`

- **Purpose:** Classify PF-failing reads of a HiSeq X BaseCalls directory into misaligned, empty, polyclonal and unexplained categories.
- **Method:** Builds the fixed argument list of `gatk CollectHiSeqXPfFailMetrics` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk CollectHiSeqXPfFailMetrics --help`. Needs the BaseCalls directory of a HiSeq X run.
- **API:** `dnakit.sequencing_qc.gatk_collect_hiseqx_pf_fail_metrics(basecalls_dir: 'str', lane: 'int', output_path: 'str', *, n_cycles: 'int | None' = None, prob_explicit_reads: 'float | None' = None, num_processors: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `basecalls_dir`, `lane`, `output_path`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2854.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `metrics`.
- **Note on the example:** the example data has no valid input for this tool (it needs data that cannot be made small: Illumina run folders, GTC files or flow-based reads, a trained or cluster model, or the Python gcnvkernel environment), so the example shows the input check that runs before GATK is started.

- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_collect_hiseqx_pf_fail_metrics

work = Path(tempfile.mkdtemp())
try:
    gatk_collect_hiseqx_pf_fail_metrics(basecalls_dir="missing_dir", lane=1, output_path=work / "o", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] basecalls_dir must be an existing directory
```

## 11) `dnakit.sequencing_qc.gatk_collect_hs_metrics()`

- **Purpose:** Collect hybrid-selection (HS) metrics: on-bait and on-target bases, fold enrichment, target coverage and GC dropout.
- **Method:** Builds the fixed argument list of `gatk CollectHsMetrics` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk CollectHsMetrics --help`.
- **API:** `dnakit.sequencing_qc.gatk_collect_hs_metrics(input_path: 'str', output_path: 'str', bait_intervals: 'Sequence[str]', target_intervals: 'Sequence[str]', *, per_target_coverage: 'str | None' = None, per_base_coverage: 'str | None' = None, theoretical_sensitivity_output: 'str | None' = None, reference: 'str | None' = None, metric_accumulation_level: 'Sequence[str] | None' = None, bait_set_name: 'str | None' = None, minimum_mapping_quality: 'int | None' = None, minimum_base_quality: 'int | None' = None, clip_overlapping_reads: 'bool | None' = None, coverage_cap: 'int | None' = None, near_distance: 'int | None' = None, include_indels: 'bool | None' = None, sample_size: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_path`, `output_path`, `bait_intervals`, `target_intervals`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2855.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `metrics`, `histogram`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_collect_hs_metrics

work = Path(tempfile.mkdtemp())
result = gatk_collect_hs_metrics("toy.bam", work / "hs.txt", bait_intervals=["toy.interval_list"], target_intervals=["toy.interval_list"], verbose=False)
row = result["metrics"][0]
print(row["BAIT_TERRITORY"], row["TARGET_TERRITORY"], row["ON_BAIT_BASES"])
```

- **Example output:**

```text
350 350 860
```

## 12) `dnakit.sequencing_qc.gatk_collect_independent_replicate_metrics()`

- **Purpose:** Estimate the fraction of duplicate-marked reads that are independent observations, using heterozygous sites of a VCF (experimental GATK tool).
- **Method:** Builds the fixed argument list of `gatk CollectIndependentReplicateMetrics` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk CollectIndependentReplicateMetrics --help`. GATK marks this tool EXPERIMENTAL (use at your own risk).
- **API:** `dnakit.sequencing_qc.gatk_collect_independent_replicate_metrics(input_path: 'str', vcf: 'str', output_path: 'str', *, matrix_output: 'str | None' = None, sample: 'str | None' = None, barcode_tag: 'str | None' = None, barcode_bq: 'str | None' = None, minimum_gq: 'int | None' = None, minimum_bq: 'int | None' = None, minimum_mq: 'int | None' = None, minimum_barcode_bq: 'int | None' = None, filter_unpaired_reads: 'bool | None' = None, reference: 'str | None' = None, stop_after: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_path`, `vcf`, `output_path`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2856.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `metrics`.
- **Note on the example:** the example data has no valid input for this tool (it needs data that cannot be made small: Illumina run folders, GTC files or flow-based reads, a trained or cluster model, or the Python gcnvkernel environment), so the example shows the input check that runs before GATK is started.

- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_collect_independent_replicate_metrics

work = Path(tempfile.mkdtemp())
try:
    gatk_collect_independent_replicate_metrics(input_path="missing.dat", vcf="missing.dat", output_path=work / "o", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] input_path must be an existing file
```

## 13) `dnakit.sequencing_qc.gatk_collect_jumping_library_metrics()`

- **Purpose:** Collect metrics of jumping (mate-pair) libraries: fragment, jump and non-jump pair counts and their orientations.
- **Method:** Builds the fixed argument list of `gatk CollectJumpingLibraryMetrics` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk CollectJumpingLibraryMetrics --help`.
- **API:** `dnakit.sequencing_qc.gatk_collect_jumping_library_metrics(input_paths: 'Sequence[str]', output_path: 'str', *, chimera_kb_min: 'int | None' = None, minimum_mapping_quality: 'int | None' = None, tail_limit: 'int | None' = None, reference: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_paths`, `output_path`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2857.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `metrics`, `histogram`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_collect_jumping_library_metrics

work = Path(tempfile.mkdtemp())
result = gatk_collect_jumping_library_metrics(["toy.bam"], work / "jump.txt", verbose=False)
row = result["metrics"][0]
print(row["JUMP_PAIRS"], row["NONJUMP_PAIRS"])
```

- **Example output:**

```text
0 19
```

## 14) `dnakit.sequencing_qc.gatk_collect_multiple_metrics()`

- **Purpose:** Run several Picard metric collectors (alignment summary, insert size, quality by cycle, quality distribution, base distribution, GC bias, sequencing artifacts, quality yield) in one pass over a file.
- **Method:** Builds the fixed argument list of `gatk CollectMultipleMetrics` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk CollectMultipleMetrics --help`. All collectors except CollectQualityYieldMetrics and CollectSequencingArtifactMetrics write PDF charts and need R; pass programs=['null', 'CollectQualityYieldMetrics'] to run without R.
- **API:** `dnakit.sequencing_qc.gatk_collect_multiple_metrics(input_path: 'str', output_prefix: 'str', *, reference: 'str | None' = None, programs: 'Sequence[str] | None' = None, extra_arguments: 'Sequence[str] | None' = None, metric_accumulation_level: 'Sequence[str] | None' = None, db_snp: 'str | None' = None, intervals_file: 'str | None' = None, ref_flat: 'str | None' = None, include_unpaired: 'bool | None' = None, assume_sorted: 'bool | None' = None, stop_after: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_path`, `output_prefix`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2858.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `metrics_files`, `collected`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_collect_multiple_metrics

work = Path(tempfile.mkdtemp())
result = gatk_collect_multiple_metrics("toy.bam", work / "multi", reference="toy.fa", programs=["null", "CollectQualityYieldMetrics"], verbose=False)
print(result["collected"])
```

- **Example output:**

```text
['quality_yield_metrics']
```

## 15) `dnakit.sequencing_qc.gatk_collect_oxog_metrics()`

- **Purpose:** Collect metrics that quantify oxidative (8-oxoguanine) artifacts: G>T error rates per sequence context.
- **Method:** Builds the fixed argument list of `gatk CollectOxoGMetrics` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk CollectOxoGMetrics --help`.
- **API:** `dnakit.sequencing_qc.gatk_collect_oxog_metrics(input_path: 'str', output_path: 'str', reference: 'str', *, db_snp: 'str | None' = None, intervals_file: 'str | None' = None, contexts: 'Sequence[str] | None' = None, context_size: 'int | None' = None, minimum_mapping_quality: 'int | None' = None, minimum_quality_score: 'int | None' = None, minimum_insert_size: 'int | None' = None, maximum_insert_size: 'int | None' = None, include_non_pf_reads: 'bool | None' = None, use_oq: 'bool | None' = None, stop_after: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_path`, `output_path`, `reference`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2859.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `metrics`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_collect_oxog_metrics

work = Path(tempfile.mkdtemp())
result = gatk_collect_oxog_metrics("toy.bam", work / "oxog.txt", reference="toy.fa", verbose=False)
print(len(result["metrics"]) > 0, result["metrics"][0]["SAMPLE_ALIAS"])
```

- **Example output:**

```text
True S1
```

## 16) `dnakit.sequencing_qc.gatk_collect_quality_yield_metrics()`

- **Purpose:** Collect yield metrics of reads that pass quality thresholds and vendor filters: variant='standard' (CollectQualityYieldMetrics), 'flow' (flow-based read files, CollectQualityYieldMetricsFlow) or 'snvq' (SNV quality, CollectQualityYieldMetricsSNVQ).
- **Method:** Builds the fixed argument list of `gatk CollectQualityYieldMetrics / CollectQualityYieldMetricsFlow / CollectQualityYieldMetricsSNVQ` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk CollectQualityYieldMetrics / CollectQualityYieldMetricsFlow / CollectQualityYieldMetricsSNVQ --help`.
- **API:** `dnakit.sequencing_qc.gatk_collect_quality_yield_metrics(input_path: 'str', output_path: 'str', variant: 'str', *, use_original_qualities: 'bool | None' = None, include_bq_histogram: 'bool | None' = None, alternate_quality_attribute: 'str | None' = None, include_secondary_alignments: 'bool | None' = None, include_supplemental_alignments: 'bool | None' = None, reference: 'str | None' = None, assume_sorted: 'bool | None' = None, stop_after: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_path`, `output_path`, `variant`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2860, G2861, G2862.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `variant`, `metrics`, `histogram`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_collect_quality_yield_metrics

work = Path(tempfile.mkdtemp())
result = gatk_collect_quality_yield_metrics("toy.bam", work / "yield.txt", variant="standard", verbose=False)
row = result["metrics"][0]
print(result["variant"] if "variant" in result else "standard", row["TOTAL_READS"], row["TOTAL_BASES"], row["Q20_BASES"])
```

- **Example output:**

```text
standard 38 1900 1900
```

## 17) `dnakit.sequencing_qc.gatk_collect_rnaseq_metrics()`

- **Purpose:** Produce RNA-seq alignment metrics: bases in coding, UTR, intronic, intergenic and ribosomal regions, strand specificity and 5'/3' bias.
- **Method:** Builds the fixed argument list of `gatk CollectRnaSeqMetrics` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk CollectRnaSeqMetrics --help`.
- **API:** `dnakit.sequencing_qc.gatk_collect_rnaseq_metrics(input_path: 'str', output_path: 'str', ref_flat: 'str', strand_specificity: 'str', *, ribosomal_intervals: 'str | None' = None, chart_output: 'str | None' = None, end_bias_bases: 'int | None' = None, minimum_length: 'int | None' = None, rrna_fragment_percentage: 'float | None' = None, ignore_sequences: 'Sequence[str] | None' = None, reference: 'str | None' = None, metric_accumulation_level: 'Sequence[str] | None' = None, assume_sorted: 'bool | None' = None, stop_after: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_path`, `output_path`, `ref_flat`, `strand_specificity`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2863.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `metrics`, `histogram`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_collect_rnaseq_metrics

work = Path(tempfile.mkdtemp())
result = gatk_collect_rnaseq_metrics("toy.bam", work / "rna.txt", ref_flat="toy.refFlat", strand_specificity="NONE", verbose=False)
row = result["metrics"][0]
print(row["PF_ALIGNED_BASES"], row["CODING_BASES"], row["INTRONIC_BASES"])
```

- **Example output:**

```text
1900 1640 260
```

## 18) `dnakit.sequencing_qc.gatk_collect_rrbs_metrics()`

- **Purpose:** Collect metrics from reduced-representation bisulfite sequencing (RRBS): CpG conversion, coverage and read-start counts.
- **Method:** Builds the fixed argument list of `gatk CollectRrbsMetrics` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk CollectRrbsMetrics --help`.
- **API:** `dnakit.sequencing_qc.gatk_collect_rrbs_metrics(input_path: 'str', metrics_file_prefix: 'str', reference: 'str', *, c_quality_threshold: 'int | None' = None, next_base_quality_threshold: 'int | None' = None, minimum_read_length: 'int | None' = None, max_mismatch_rate: 'float | None' = None, sequence_names: 'Sequence[str] | None' = None, metric_accumulation_level: 'Sequence[str] | None' = None, assume_sorted: 'bool | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_path`, `metrics_file_prefix`, `reference`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2864.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `metrics_files`, `metrics`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_collect_rrbs_metrics

work = Path(tempfile.mkdtemp())
result = gatk_collect_rrbs_metrics("toy.bam", work / "rrbs", reference="toy.fa", verbose=False)
print(result["status"], len(result["metrics_files"]) >= 2)
```

- **Example output:**

```text
completed True
```

## 19) `dnakit.sequencing_qc.gatk_collect_sequencing_artifact_metrics()`

- **Purpose:** Quantify single-base sequencing artifacts (pre-adapter, bait-bias, error summaries by substitution and context).
- **Method:** Builds the fixed argument list of `gatk CollectSequencingArtifactMetrics` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk CollectSequencingArtifactMetrics --help`.
- **API:** `dnakit.sequencing_qc.gatk_collect_sequencing_artifact_metrics(input_path: 'str', output_prefix: 'str', reference: 'str', *, db_snp: 'str | None' = None, intervals_file: 'str | None' = None, contexts_to_print: 'Sequence[str] | None' = None, context_size: 'int | None' = None, minimum_mapping_quality: 'int | None' = None, minimum_quality_score: 'int | None' = None, minimum_insert_size: 'int | None' = None, maximum_insert_size: 'int | None' = None, include_duplicates: 'bool | None' = None, include_non_pf_reads: 'bool | None' = None, include_unpaired: 'bool | None' = None, tandem_reads: 'bool | None' = None, use_oq: 'bool | None' = None, file_extension: 'str | None' = None, assume_sorted: 'bool | None' = None, stop_after: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_path`, `output_prefix`, `reference`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2866.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `metrics_files`, `metrics`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_collect_sequencing_artifact_metrics

work = Path(tempfile.mkdtemp())
result = gatk_collect_sequencing_artifact_metrics("toy.bam", work / "artifact", reference="toy.fa", verbose=False)
print(sorted(path.rsplit(".", 1)[1] for path in result["metrics_files"]))
```

- **Example output:**

```text
['bait_bias_detail_metrics', 'bait_bias_summary_metrics', 'error_summary_metrics', 'pre_adapter_detail_metrics', 'pre_adapter_summary_metrics']
```

## 20) `dnakit.sequencing_qc.gatk_collect_targeted_pcr_metrics()`

- **Purpose:** Calculate PCR-related metrics from targeted amplicon sequencing: on-amplicon and on-target bases, coverage and sensitivity.
- **Method:** Builds the fixed argument list of `gatk CollectTargetedPcrMetrics` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk CollectTargetedPcrMetrics --help`.
- **API:** `dnakit.sequencing_qc.gatk_collect_targeted_pcr_metrics(input_path: 'str', output_path: 'str', amplicon_intervals: 'str', target_intervals: 'Sequence[str]', *, per_target_coverage: 'str | None' = None, per_base_coverage: 'str | None' = None, theoretical_sensitivity_output: 'str | None' = None, custom_amplicon_set_name: 'str | None' = None, reference: 'str | None' = None, metric_accumulation_level: 'Sequence[str] | None' = None, minimum_mapping_quality: 'int | None' = None, minimum_base_quality: 'int | None' = None, clip_overlapping_reads: 'bool | None' = None, coverage_cap: 'int | None' = None, near_distance: 'int | None' = None, include_indels: 'bool | None' = None, sample_size: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_path`, `output_path`, `amplicon_intervals`, `target_intervals`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2867.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `metrics`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_collect_targeted_pcr_metrics

work = Path(tempfile.mkdtemp())
result = gatk_collect_targeted_pcr_metrics("toy.bam", work / "pcr.txt", amplicon_intervals="toy.interval_list", target_intervals=["toy.interval_list"], verbose=False)
row = result["metrics"][0]
print(row["AMPLICON_TERRITORY"], row["ON_AMPLICON_BASES"])
```

- **Example output:**

```text
350 860
```

## 21) `dnakit.sequencing_qc.gatk_collect_umi_prevalence_metrics()`

- **Purpose:** Tally the counts of UMIs in duplicate sets of a duplicate-marked BAM (how many duplicate sets have 1, 2, 3 ... distinct UMIs).
- **Method:** Builds the fixed argument list of `gatk CollectUmiPrevalenceMetrics` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk CollectUmiPrevalenceMetrics --help`.
- **API:** `dnakit.sequencing_qc.gatk_collect_umi_prevalence_metrics(input_path: 'str', output_path: 'str', *, barcode_tag: 'str | None' = None, barcode_bq: 'str | None' = None, minimum_barcode_bq: 'int | None' = None, minimum_mq: 'int | None' = None, filter_unpaired_reads: 'bool | None' = None, reference: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_path`, `output_path`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2868.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `metrics`.
- **Note on the example:** the example data has no valid input for this tool (it needs data that cannot be made small: Illumina run folders, GTC files or flow-based reads, a trained or cluster model, or the Python gcnvkernel environment), so the example shows the input check that runs before GATK is started.

- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_collect_umi_prevalence_metrics

work = Path(tempfile.mkdtemp())
try:
    gatk_collect_umi_prevalence_metrics(input_path="missing.dat", output_path=work / "o", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] input_path must be an existing file
```

## 22) `dnakit.sequencing_qc.gatk_collect_variant_calling_metrics()`

- **Purpose:** Collect per-sample and aggregate variant-calling metrics (SNP/indel counts, Ti/Tv, het/hom ratio, dbSNP concordance) from a VCF.
- **Method:** Builds the fixed argument list of `gatk CollectVariantCallingMetrics` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk CollectVariantCallingMetrics --help`.
- **API:** `dnakit.sequencing_qc.gatk_collect_variant_calling_metrics(input_path: 'str', dbsnp: 'str', output_prefix: 'str', *, target_intervals: 'str | None' = None, sequence_dictionary: 'str | None' = None, gvcf_input: 'bool | None' = None, thread_count: 'int | None' = None, reference: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_path`, `dbsnp`, `output_prefix`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2869.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `metrics_files`, `summary_metrics`, `detail_metrics`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_collect_variant_calling_metrics

work = Path(tempfile.mkdtemp())
result = gatk_collect_variant_calling_metrics("toy.vcf.gz", "toy.vcf.gz", work / "calling", verbose=False)
row = result["summary_metrics"][0]
print(row["TOTAL_SNPS"], row["NUM_IN_DB_SNP"], row["PCT_DBSNP"])
```

- **Example output:**

```text
15 15 1
```

## 23) `dnakit.sequencing_qc.gatk_collect_wgs_metrics()`

- **Purpose:** Collect coverage and performance metrics of a whole-genome sequencing experiment: mean coverage, fraction of bases at depth thresholds, excluded-base fractions.
- **Method:** Builds the fixed argument list of `gatk CollectWgsMetrics` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk CollectWgsMetrics --help`.
- **API:** `dnakit.sequencing_qc.gatk_collect_wgs_metrics(input_path: 'str', output_path: 'str', reference: 'str', *, intervals_file: 'str | None' = None, minimum_mapping_quality: 'int | None' = None, minimum_base_quality: 'int | None' = None, coverage_cap: 'int | None' = None, read_length: 'int | None' = None, count_unpaired: 'bool | None' = None, include_bq_histogram: 'bool | None' = None, use_fast_algorithm: 'bool | None' = None, sample_size: 'int | None' = None, locus_accumulation_cap: 'int | None' = None, theoretical_sensitivity_output: 'str | None' = None, stop_after: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_path`, `output_path`, `reference`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2870.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `metrics`, `histogram`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_collect_wgs_metrics

work = Path(tempfile.mkdtemp())
result = gatk_collect_wgs_metrics("toy.bam", work / "wgs.txt", reference="toy.fa", minimum_base_quality=1, minimum_mapping_quality=1, verbose=False)
row = result["metrics"][0]
print(row["GENOME_TERRITORY"], row["MEAN_COVERAGE"], row["PCT_EXC_DUPE"])
```

- **Example output:**

```text
1000 1.9 0
```

## 24) `dnakit.sequencing_qc.gatk_calculate_fingerprint_metrics()`

- **Purpose:** Calculate statistics on fingerprints (SAM/BAM/CRAM or VCF) against a haplotype map, checking their viability.
- **Method:** Builds the fixed argument list of `gatk CalculateFingerprintMetrics` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk CalculateFingerprintMetrics --help`.
- **API:** `dnakit.sequencing_qc.gatk_calculate_fingerprint_metrics(haplotype_map: 'str', input_paths: 'Sequence[str]', output_path: 'str', *, calculate_by: 'str | None' = None, genotype_lod_threshold: 'float | None' = None, number_of_sampling: 'int | None' = None, reference: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `haplotype_map`, `input_paths`, `output_path`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2844.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `metrics`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_calculate_fingerprint_metrics

work = Path(tempfile.mkdtemp())
result = gatk_calculate_fingerprint_metrics("toy.haplotype_map", ["toy.bam"], work / "fingerprint_metrics.txt", verbose=False)
row = result["metrics"][0]
print(row["SAMPLE_ALIAS"], row["HAPLOTYPES"], row["HAPLOTYPES_WITH_EVIDENCE"])
```

- **Example output:**

```text
S1 15 13
```

## 25) `dnakit.sequencing_qc.gatk_calculate_read_group_checksum()`

- **Purpose:** Create a hash code over the read groups (RG) of a SAM/BAM/CRAM file; files whose read-group content is the same give the same hash.
- **Method:** Builds the fixed argument list of `gatk CalculateReadGroupChecksum` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk CalculateReadGroupChecksum --help`.
- **API:** `dnakit.sequencing_qc.gatk_calculate_read_group_checksum(input_path: 'str', output_path: 'str', *, reference: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_path`, `output_path`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2845.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `text`, `checksum`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_calculate_read_group_checksum

work = Path(tempfile.mkdtemp())
a = gatk_calculate_read_group_checksum("toy.bam", work / "a.txt", verbose=False)
b = gatk_calculate_read_group_checksum("umi.bam", work / "b.txt", verbose=False)
print(len(a["checksum"]), a["checksum"] == b["checksum"])
```

- **Example output:**

```text
32 True
```

## 26) `dnakit.sequencing_qc.gatk_check_duplicate_marking()`

- **Purpose:** Check that reads of the same name (and mates) carry consistent duplicate markings.
- **Method:** Builds the fixed argument list of `gatk CheckDuplicateMarking` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk CheckDuplicateMarking --help`.
- **API:** `dnakit.sequencing_qc.gatk_check_duplicate_marking(input_path: 'str', *, output_path: 'str | None' = None, mode: 'str | None' = None, reference: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_path`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2846.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `has_findings`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_check_duplicate_marking

work = Path(tempfile.mkdtemp())
result = gatk_check_duplicate_marking("toy.bam", verbose=False)
print(result["status"], result["has_findings"])
```

- **Example output:**

```text
completed False
```

## 27) `dnakit.sequencing_qc.gatk_check_fingerprint()`

- **Purpose:** Compute a fingerprint from a SAM/BAM/CRAM or VCF file and compare it with the provided genotypes (LOD of the same-sample hypothesis).
- **Method:** Builds the fixed argument list of `gatk CheckFingerprint` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk CheckFingerprint --help`.
- **API:** `dnakit.sequencing_qc.gatk_check_fingerprint(input_path: 'str', genotypes: 'str', haplotype_map: 'str', *, output_prefix: 'str | None' = None, summary_output: 'str | None' = None, detail_output: 'str | None' = None, expected_sample_alias: 'str | None' = None, observed_sample_alias: 'str | None' = None, genotype_lod_threshold: 'float | None' = None, ignore_read_groups: 'bool | None' = None, reference: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_path`, `genotypes`, `haplotype_map`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2847.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `summary_metrics`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_check_fingerprint

work = Path(tempfile.mkdtemp())
result = gatk_check_fingerprint("toy.vcf.gz", "toy.vcf.gz", "toy.haplotype_map", output_prefix=work / "check", verbose=False)
row = result["summary_metrics"][0]
print(row["SAMPLE"], row["HAPLOTYPES_WITH_GENOTYPES"], round(row["LOD_EXPECTED_SAMPLE"], 2))
```

- **Example output:**

```text
S1 15 4.71
```

## 28) `dnakit.sequencing_qc.gatk_check_pileup()`

- **Purpose:** Compare GATK's internal pileup with a pileup made by samtools mpileup and report differences.
- **Method:** Builds the fixed argument list of `gatk CheckPileup` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk CheckPileup --help`.
- **API:** `dnakit.sequencing_qc.gatk_check_pileup(input_paths: 'Sequence[str]', pileup: 'str', reference: 'str', *, output_path: 'str | None' = None, intervals: 'Sequence[str] | None' = None, ignore_overlaps: 'bool | None' = None, continue_after_error: 'bool | None' = None, max_depth_per_sample: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_paths`, `pileup`, `reference`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2848.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `text`.
- **Note on the example:** the example data has no valid input for this tool (it needs data that cannot be made small: Illumina run folders, GTC files or flow-based reads, a trained or cluster model, or the Python gcnvkernel environment), so the example shows the input check that runs before GATK is started.

- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_check_pileup

work = Path(tempfile.mkdtemp())
try:
    gatk_check_pileup(input_paths=["missing.dat"], pileup="missing.dat", reference="missing.dat", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] input_paths must be an existing file
```

## 29) `dnakit.sequencing_qc.gatk_check_terminator_block()`

- **Purpose:** Check that the last block of a BGZF-compressed file (for example a BAM) is the well-formed terminator block.
- **Method:** Builds the fixed argument list of `gatk CheckTerminatorBlock` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk CheckTerminatorBlock --help`.
- **API:** `dnakit.sequencing_qc.gatk_check_terminator_block(input_path: 'str', *, reference: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_path`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2849.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `terminator_ok`, `has_findings`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_check_terminator_block

work = Path(tempfile.mkdtemp())
result = gatk_check_terminator_block("toy.bam", verbose=False)
print(result["terminator_ok"])
```

- **Example output:**

```text
True
```

## 30) `dnakit.sequencing_qc.gatk_compare_base_qualities()`

- **Purpose:** Compare the base qualities of two SAM/BAM/CRAM files that hold the same reads in the same order.
- **Method:** Builds the fixed argument list of `gatk CompareBaseQualities` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk CompareBaseQualities --help`.
- **API:** `dnakit.sequencing_qc.gatk_compare_base_qualities(input_paths: 'Sequence[str]', *, output_path: 'str | None' = None, throw_on_diff: 'bool | None' = None, round_down_quantized: 'bool | None' = None, reference: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_paths`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2871.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `text`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_compare_base_qualities

work = Path(tempfile.mkdtemp())
result = gatk_compare_base_qualities(["toy.bam", "umi.bam"], output_path=work / "bq.txt", verbose=False)
print(result["text"].splitlines()[1])
```

- **Example output:**

```text
all 1900 quality scores are the same
```

## 31) `dnakit.sequencing_qc.gatk_compare_duplicates_spark()`

- **Purpose:** Determine whether two potentially identical BAMs have the same duplicate reads (Spark tool, run locally).
- **Method:** Builds the fixed argument list of `gatk CompareDuplicatesSpark` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk CompareDuplicatesSpark --help`.
- **API:** `dnakit.sequencing_qc.gatk_compare_duplicates_spark(input_paths: 'Sequence[str]', input2: 'str', *, output_path: 'str | None' = None, output2_path: 'str | None' = None, throw_on_diff: 'bool | None' = None, print_summary: 'bool | None' = None, intervals: 'Sequence[str] | None' = None, reference: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_paths`, `input2`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2872.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `summary_lines`.
- **Note on the example:** the example data has no valid input for this tool (it needs data that cannot be made small: Illumina run folders, GTC files or flow-based reads, a trained or cluster model, or the Python gcnvkernel environment), so the example shows the input check that runs before GATK is started.

- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_compare_duplicates_spark

work = Path(tempfile.mkdtemp())
try:
    gatk_compare_duplicates_spark(input_paths=["missing.dat"], input2="missing.dat", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] input_paths must be an existing file
```

## 32) `dnakit.sequencing_qc.gatk_compare_sams()`

- **Purpose:** Compare two SAM/BAM/CRAM files (headers, alignments and duplicate marks) strictly or leniently and write the comparison.
- **Method:** Builds the fixed argument list of `gatk CompareSAMs` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk CompareSAMs --help`.
- **API:** `dnakit.sequencing_qc.gatk_compare_sams(input_paths: 'Sequence[str]', *, output_path: 'str | None' = None, lenient_dup: 'bool | None' = None, lenient_header: 'bool | None' = None, lenient_low_mq_alignment: 'bool | None' = None, lenient_unknown_mq_alignment: 'bool | None' = None, low_mq_threshold: 'int | None' = None, compare_mq: 'bool | None' = None, reference: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_paths`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2873.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `files_match`, `has_findings`, `text`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_compare_sams

work = Path(tempfile.mkdtemp())
result = gatk_compare_sams(["toy.bam", "umi.bam"], output_path=work / "cmp.txt", verbose=False)
print(result["files_match"], result["has_findings"])
```

- **Example output:**

```text
True False
```

## 33) `dnakit.sequencing_qc.gatk_convert_haplotype_database_to_vcf()`

- **Purpose:** Convert a haplotype database (haplotype map) file to a VCF.
- **Method:** Builds the fixed argument list of `gatk ConvertHaplotypeDatabaseToVcf` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk ConvertHaplotypeDatabaseToVcf --help`.
- **API:** `dnakit.sequencing_qc.gatk_convert_haplotype_database_to_vcf(input_path: 'str', output_path: 'str', reference: 'str', *, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_path`, `output_path`, `reference`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2874.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `record_count`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_convert_haplotype_database_to_vcf

work = Path(tempfile.mkdtemp())
result = gatk_convert_haplotype_database_to_vcf("toy.haplotype_map", work / "map.vcf", reference="toy.fa", verbose=False)
print(result["record_count"])
```

- **Example output:**

```text
15
```

## 34) `dnakit.sequencing_qc.gatk_crosscheck_fingerprints()`

- **Purpose:** Check that all data in the input files appear to have come from the same individual (or that groups that should match do).
- **Method:** Builds the fixed argument list of `gatk CrosscheckFingerprints` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk CrosscheckFingerprints --help`.
- **API:** `dnakit.sequencing_qc.gatk_crosscheck_fingerprints(haplotype_map: 'str', input_paths: 'Sequence[str]', *, output_path: 'str | None' = None, second_input_paths: 'Sequence[str] | None' = None, crosscheck_by: 'str | None' = None, crosscheck_mode: 'str | None' = None, lod_threshold: 'float | None' = None, expect_all_groups_to_match: 'bool | None' = None, output_errors_only: 'bool | None' = None, matrix_output: 'str | None' = None, num_threads: 'int | None' = None, reference: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `haplotype_map`, `input_paths`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2875.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `metrics`, `has_findings`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_crosscheck_fingerprints

work = Path(tempfile.mkdtemp())
result = gatk_crosscheck_fingerprints("toy.haplotype_map", ["toy.bam", "umi.bam"], output_path=work / "cross.txt", verbose=False)
print(result["has_findings"], sorted({row["RESULT"] for row in result["metrics"]}))
```

- **Example output:**

```text
False ['EXPECTED_MATCH']
```

## 35) `dnakit.sequencing_qc.gatk_dump_tabix_index()`

- **Purpose:** Dump a tabix index file (.tbi) as text: sequence names, bins, chunks and linear index.
- **Method:** Builds the fixed argument list of `gatk DumpTabixIndex` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk DumpTabixIndex --help`.
- **API:** `dnakit.sequencing_qc.gatk_dump_tabix_index(tabix_index: 'str', output_path: 'str', *, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `tabix_index`, `output_path`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2876.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `text`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_dump_tabix_index

work = Path(tempfile.mkdtemp())
result = gatk_dump_tabix_index("toy.vcf.gz.tbi", work / "tbi.txt", verbose=False)
print(result["text"].splitlines()[3])
```

- **Example output:**

```text
chr1 binned index:
```

## 36) `dnakit.sequencing_qc.gatk_extract_fingerprint()`

- **Purpose:** Compute a fingerprint from a SAM/BAM/CRAM file against a haplotype map and write it as a VCF.
- **Method:** Builds the fixed argument list of `gatk ExtractFingerprint` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk ExtractFingerprint --help`.
- **API:** `dnakit.sequencing_qc.gatk_extract_fingerprint(input_path: 'str', haplotype_map: 'str', output_path: 'str', reference: 'str', *, contamination: 'float | None' = None, sample_alias: 'str | None' = None, locus_max_reads: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_path`, `haplotype_map`, `output_path`, `reference`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2878.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `record_count`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_extract_fingerprint

work = Path(tempfile.mkdtemp())
result = gatk_extract_fingerprint("toy.bam", "toy.haplotype_map", work / "fingerprint.vcf", reference="toy.fa", verbose=False)
print(result["record_count"])
```

- **Example output:**

```text
15
```

## 37) `dnakit.sequencing_qc.gatk_get_sample_name()`

- **Purpose:** Emit the single sample name of a SAM/BAM/CRAM file to a file.
- **Method:** Builds the fixed argument list of `gatk GetSampleName` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk GetSampleName --help`.
- **API:** `dnakit.sequencing_qc.gatk_get_sample_name(input_paths: 'Sequence[str]', output_path: 'str', *, use_url_encoding: 'bool | None' = None, reference: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_paths`, `output_path`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2880.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `sample_name`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_get_sample_name

work = Path(tempfile.mkdtemp())
result = gatk_get_sample_name(["toy.bam"], work / "sample.txt", verbose=False)
print(result["sample_name"])
```

- **Example output:**

```text
S1
```

## 38) `dnakit.sequencing_qc.gatk_identify_contaminant()`

- **Purpose:** Compute a fingerprint of the contaminating sample from a SAM/BAM/CRAM file, given a contamination estimate, and write it as a VCF.
- **Method:** Builds the fixed argument list of `gatk IdentifyContaminant` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk IdentifyContaminant --help`.
- **API:** `dnakit.sequencing_qc.gatk_identify_contaminant(input_path: 'str', haplotype_map: 'str', output_path: 'str', reference: 'str', *, contamination: 'float | None' = None, extract_contaminated: 'bool | None' = None, sample_alias: 'str | None' = None, locus_max_reads: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_path`, `haplotype_map`, `output_path`, `reference`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2881.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `record_count`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_identify_contaminant

work = Path(tempfile.mkdtemp())
result = gatk_identify_contaminant("toy.bam", "toy.haplotype_map", work / "contaminant.vcf", reference="toy.fa", contamination=0.1, verbose=False)
print(result["record_count"])
```

- **Example output:**

```text
15
```

## 39) `dnakit.sequencing_qc.gatk_lift_over_haplotype_map()`

- **Purpose:** Lift over a haplotype database (haplotype map) from one reference build to another with a chain file.
- **Method:** Builds the fixed argument list of `gatk LiftOverHaplotypeMap` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk LiftOverHaplotypeMap --help`.
- **API:** `dnakit.sequencing_qc.gatk_lift_over_haplotype_map(input_path: 'str', chain: 'str', sequence_dictionary: 'str', output_path: 'str', *, reference: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_path`, `chain`, `sequence_dictionary`, `output_path`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2882.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_lift_over_haplotype_map

work = Path(tempfile.mkdtemp())
result = gatk_lift_over_haplotype_map("old.haplotype_map", "shift.chain", "toy.dict", work / "lifted.haplotype_map", verbose=False)
print((work / "lifted.haplotype_map").read_text().splitlines()[-1].split("\t")[:3])
```

- **Example output:**

```text
['chr1', '40', 'rs1']
```

## 40) `dnakit.sequencing_qc.gatk_validate_sam_file()`

- **Purpose:** Validate a SAM/BAM/CRAM file (headers, flags, mates, CIGAR, indices) and report errors and warnings.
- **Method:** Builds the fixed argument list of `gatk ValidateSamFile` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk ValidateSamFile --help`.
- **API:** `dnakit.sequencing_qc.gatk_validate_sam_file(input_path: 'str', *, mode: 'str | None' = None, output_path: 'str | None' = None, ignore: 'Sequence[str] | None' = None, ignore_warnings: 'bool | None' = None, max_output: 'int | None' = None, is_bisulfite_sequenced: 'bool | None' = None, skip_mate_validation: 'bool | None' = None, index_validation_stringency: 'str | None' = None, reference: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_path`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2885.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `valid`, `has_findings`, `problems`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_validate_sam_file

work = Path(tempfile.mkdtemp())
result = gatk_validate_sam_file("toy.bam", mode="SUMMARY", verbose=False)
print(result["valid"], result["has_findings"])
```

- **Example output:**

```text
True False
```

## 41) `dnakit.sequencing_qc.gatk_view_sam()`

- **Purpose:** Print the records (and header) of a SAM/BAM/CRAM file, optionally only aligned or unaligned and only PF or non-PF reads.
- **Method:** Builds the fixed argument list of `gatk ViewSam` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk ViewSam --help`.
- **API:** `dnakit.sequencing_qc.gatk_view_sam(input_path: 'str', alignment_status: 'str', pf_status: 'str', *, header_only: 'bool | None' = None, records_only: 'bool | None' = None, interval_list: 'str | None' = None, reference: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_path`, `alignment_status`, `pf_status`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2886.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `header_lines`, `record_count`, `records`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_view_sam

work = Path(tempfile.mkdtemp())
result = gatk_view_sam("toy.bam", alignment_status="All", pf_status="All", records_only=True, verbose=False)
print(result["record_count"], result["records"][0].split("\t")[:6])
```

- **Example output:**

```text
38 ['r1', '99', 'chr1', '1', '60', '50M']
```

## 42) `dnakit.sequencing_qc.gatk_htsget_download()`

- **Purpose:** Download a file from a GA4GH htsget endpoint (a record id on a server URL) to a local file; needs network access.
- **Method:** Builds the fixed argument list of `gatk HtsgetReader` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk HtsgetReader --help`.
- **API:** `dnakit.sequencing_qc.gatk_htsget_download(url: 'str', record_id: 'str', output_path: 'str', *, data_format: 'str | None' = None, data_class: 'str | None' = None, intervals: 'str | None' = None, check_md5: 'bool | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `url`, `record_id`, `output_path`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2887.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`.
- **Note on the example:** the example data has no valid input for this tool (it needs data that cannot be made small: Illumina run folders, GTC files or flow-based reads, a trained or cluster model, or the Python gcnvkernel environment), so the example shows the input check that runs before GATK is started.

- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_htsget_download

work = Path(tempfile.mkdtemp())
try:
    gatk_htsget_download(url="ftp://x", record_id="x", output_path=work / "o", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] url must be an http(s) URL of an htsget endpoint.
```

## 43) `dnakit.sequencing_qc.gatk_split_cram()`

- **Purpose:** Split a CRAM file into smaller CRAM files of a given number of records, without recompressing the containers.
- **Method:** Builds the fixed argument list of `gatk SplitCRAM` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk SplitCRAM --help`.
- **API:** `dnakit.sequencing_qc.gatk_split_cram(input_path: 'str', output_dir: 'str', *, shard_records: 'int | None' = None, shard_max_output_count: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_path`, `output_dir`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2895.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `shards`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_split_cram

work = Path(tempfile.mkdtemp())
result = gatk_split_cram("toy.cram", work / "shards", shard_records=10, verbose=False)
print([path.rsplit("/", 1)[1] for path in result["shards"]])
```

- **Example output:**

```text
['shard_0000.cram', 'shard_0001.cram']
```
