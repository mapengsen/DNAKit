# GATK command-line facades: interval lists, genotyping arrays, flow reads and PathSeq

Typed facades over the `gatk` command line for Picard interval lists (conversion, comparison, set operations, liftover, scattering), Illumina genotyping arrays (GTC to VCF, manifests, ADPC, array metrics), flow-based reads of Ultima-style data, and the PathSeq reference builders. Each function builds a fixed `gatk` command line without a shell, writes only to new paths and returns a JSON-serialisable dictionary; GATK itself is licensed BSD-3-Clause (Picard MIT) and is not bundled. The functions live in `dnakit.annotation.gatk_p00_intervals_cli_functions` and are exported from `dnakit.annotation`; the examples import from the module so that they run either way.

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

## 1) `dnakit.annotation.gatk_convert_interval_list()`

- **Purpose:** Convert between BED and Picard interval-list formats: output_format='interval_list' runs BedToIntervalList, output_format='bed' runs IntervalListToBed.
- **Method:** Builds the fixed argument list of `gatk BedToIntervalList / IntervalListToBed` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk BedToIntervalList / IntervalListToBed --help`.
- **API:** `dnakit.annotation.gatk_convert_interval_list(input_path: 'str', output_path: 'str', output_format: 'str', *, sequence_dictionary: 'str | None' = None, keep_length_zero_intervals: 'bool | None' = None, unique: 'bool | None' = None, score: 'int | None' = None, sort: 'bool | None' = None, reference: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_path`, `output_path`, `output_format`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2903, G2905.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `output_format`, `interval_count`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.annotation.gatk_p00_intervals_cli_functions import gatk_convert_interval_list

work = Path(tempfile.mkdtemp())
to_list = gatk_convert_interval_list("toy.bed", work / "from_bed.interval_list", output_format="interval_list", sequence_dictionary="toy.dict", verbose=False)
to_bed = gatk_convert_interval_list("toy.interval_list", work / "from_list.bed", output_format="bed", verbose=False)
print(to_list["interval_count"], (work / "from_list.bed").read_text().splitlines()[0])
```

- **Example output:**

```text
2 chr1	100	300	t1	500	+
```

## 2) `dnakit.annotation.gatk_compare_interval_lists()`

- **Purpose:** Compare two interval lists for equality and report the intervals that differ.
- **Method:** Builds the fixed argument list of `gatk CompareIntervalLists` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk CompareIntervalLists --help`.
- **API:** `dnakit.annotation.gatk_compare_interval_lists(first: 'str', second: 'str', reference: 'str', *, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `first`, `second`, `reference`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2904.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `identical`, `has_findings`, `report_lines`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.annotation.gatk_p00_intervals_cli_functions import gatk_compare_interval_lists

work = Path(tempfile.mkdtemp())
same = gatk_compare_interval_lists("toy.interval_list", "toy.interval_list", reference="toy.fa", verbose=False)
print(same["identical"], same["report_lines"])
```

- **Example output:**

```text
True ['Intervals are equal']
```

## 3) `dnakit.annotation.gatk_interval_list_tools()`

- **Purpose:** Manipulate Picard interval lists: concatenate, union, intersect, subtract, symmetric difference, pad, invert, merge, break bands and scatter.
- **Method:** Builds the fixed argument list of `gatk IntervalListTools` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk IntervalListTools --help`.
- **API:** `dnakit.annotation.gatk_interval_list_tools(input_paths: 'Sequence[str]', *, output_path: 'str | None' = None, output_directory: 'str | None' = None, action: 'str | None' = None, second_input_paths: 'Sequence[str] | None' = None, padding: 'int | None' = None, invert: 'bool | None' = None, unique: 'bool | None' = None, sort: 'bool | None' = None, dont_merge_abutting: 'bool | None' = None, break_bands_at_multiples_of: 'int | None' = None, scatter_count: 'int | None' = None, scatter_content: 'int | None' = None, subdivision_mode: 'str | None' = None, count_output: 'str | None' = None, output_value: 'str | None' = None, comments: 'Sequence[str] | None' = None, reference: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_paths`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2906.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `interval_count`, `scattered_lists`, `count`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.annotation.gatk_p00_intervals_cli_functions import gatk_interval_list_tools

work = Path(tempfile.mkdtemp())
padded = gatk_interval_list_tools(["toy.interval_list"], output_path=work / "padded.interval_list", padding=10, verbose=False)
print((work / "padded.interval_list").read_text().splitlines()[-2:])
```

- **Example output:**

```text
['chr1\t91\t310\t+\tt1', 'chr2\t41\t210\t+\tt2']
```

## 4) `dnakit.annotation.gatk_lift_over_interval_list()`

- **Purpose:** Lift over an interval list from one reference build to another with a chain file; intervals that cannot be lifted are rejected.
- **Method:** Builds the fixed argument list of `gatk LiftOverIntervalList` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk LiftOverIntervalList --help`.
- **API:** `dnakit.annotation.gatk_lift_over_interval_list(input_path: 'str', chain: 'str', sequence_dictionary: 'str', output_path: 'str', *, reject: 'str | None' = None, min_liftover_pct: 'float | None' = None, reference: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_path`, `chain`, `sequence_dictionary`, `output_path`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2907.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `interval_count`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.annotation.gatk_p00_intervals_cli_functions import gatk_lift_over_interval_list

work = Path(tempfile.mkdtemp())
result = gatk_lift_over_interval_list("old.interval_list", "shift.chain", "toy.dict", work / "lifted.interval_list", reject=work / "rejected.interval_list", verbose=False)
print((work / "lifted.interval_list").read_text().splitlines()[-1])
```

- **Example output:**

```text
chr1	111	210	+	t1
```

## 5) `dnakit.annotation.gatk_split_intervals()`

- **Purpose:** Split an interval list (or a whole reference) into sub-interval files for scatter-gather processing.
- **Method:** Builds the fixed argument list of `gatk SplitIntervals` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk SplitIntervals --help`.
- **API:** `dnakit.annotation.gatk_split_intervals(output_dir: 'str', *, reference: 'str | None' = None, intervals: 'Sequence[str] | None' = None, scatter_count: 'int | None' = None, subdivision_mode: 'str | None' = None, dont_mix_contigs: 'bool | None' = None, min_contig_size: 'int | None' = None, interval_file_prefix: 'str | None' = None, extension: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `output_dir`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2909.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `interval_files`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.annotation.gatk_p00_intervals_cli_functions import gatk_split_intervals

work = Path(tempfile.mkdtemp())
result = gatk_split_intervals(work / "scatter", intervals=["toy.interval_list"], reference="toy.fa", scatter_count=2, verbose=False)
print([path.rsplit("/", 1)[1] for path in result["interval_files"]])
```

- **Example output:**

```text
['0000-scattered.interval_list', '0001-scattered.interval_list']
```

## 6) `dnakit.annotation.gatk_bpm_to_normalization_manifest_csv()`

- **Purpose:** Convert an Illumina bead pool manifest (.bpm) and cluster file into a bpm.csv normalization manifest.
- **Method:** Builds the fixed argument list of `gatk BpmToNormalizationManifestCsv` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk BpmToNormalizationManifestCsv --help`.
- **API:** `dnakit.annotation.gatk_bpm_to_normalization_manifest_csv(input_path: 'str', cluster_file: 'str', output_path: 'str', *, reference: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_path`, `cluster_file`, `output_path`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2896.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `row_count`.
- **Note on the example:** the example data has no valid input for this tool (it needs data that cannot be made small: Illumina run folders, GTC files or flow-based reads, a trained or cluster model, or the Python gcnvkernel environment), so the example shows the input check that runs before GATK is started.

- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.annotation.gatk_p00_intervals_cli_functions import gatk_bpm_to_normalization_manifest_csv

work = Path(tempfile.mkdtemp())
try:
    gatk_bpm_to_normalization_manifest_csv(input_path="missing.dat", cluster_file="missing.dat", output_path=work / "o", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] input_path must be an existing file
```

## 7) `dnakit.annotation.gatk_combine_genotyping_array_vcfs()`

- **Purpose:** Combine several genotyping-array VCF files into one VCF.
- **Method:** Builds the fixed argument list of `gatk CombineGenotypingArrayVcfs` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk CombineGenotypingArrayVcfs --help`.
- **API:** `dnakit.annotation.gatk_combine_genotyping_array_vcfs(input_paths: 'Sequence[str]', output_path: 'str', *, reference: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_paths`, `output_path`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2897.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `record_count`.
- **Note on the example:** the example data has no valid input for this tool (it needs data that cannot be made small: Illumina run folders, GTC files or flow-based reads, a trained or cluster model, or the Python gcnvkernel environment), so the example shows the input check that runs before GATK is started.

- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.annotation.gatk_p00_intervals_cli_functions import gatk_combine_genotyping_array_vcfs

work = Path(tempfile.mkdtemp())
try:
    gatk_combine_genotyping_array_vcfs(input_paths=["missing.dat"], output_path=work / "o", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] input_paths must be an existing file
```

## 8) `dnakit.annotation.gatk_compare_gtc_files()`

- **Purpose:** Compare two Illumina GTC files (genotypes and intensities) against a bead pool manifest.
- **Method:** Builds the fixed argument list of `gatk CompareGtcFiles` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk CompareGtcFiles --help`.
- **API:** `dnakit.annotation.gatk_compare_gtc_files(input_paths: 'Sequence[str]', illumina_bead_pool_manifest_file: 'str', *, reference: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_paths`, `illumina_bead_pool_manifest_file`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2898.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `files_match`, `has_findings`.
- **Note on the example:** the example data has no valid input for this tool (it needs data that cannot be made small: Illumina run folders, GTC files or flow-based reads, a trained or cluster model, or the Python gcnvkernel environment), so the example shows the input check that runs before GATK is started.

- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.annotation.gatk_p00_intervals_cli_functions import gatk_compare_gtc_files

work = Path(tempfile.mkdtemp())
try:
    gatk_compare_gtc_files(input_paths=["missing.dat"], illumina_bead_pool_manifest_file="missing.dat", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] input_paths must hold exactly two GTC files.
```

## 9) `dnakit.annotation.gatk_create_extended_illumina_manifest()`

- **Purpose:** Create the Extended Illumina Manifest (reference-specific, with mapped positions and alleles) that GtcToVcf needs.
- **Method:** Builds the fixed argument list of `gatk CreateExtendedIlluminaManifest` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk CreateExtendedIlluminaManifest --help`.
- **API:** `dnakit.annotation.gatk_create_extended_illumina_manifest(input_path: 'str', output_path: 'str', reference_sequence: 'str', report_file: 'str', *, bad_assays_file: 'str | None' = None, cluster_file: 'str | None' = None, dbsnp_file: 'str | None' = None, target_build: 'str | None' = None, supported_builds: 'Sequence[str] | None' = None, supported_reference_files: 'Sequence[str] | None' = None, supported_chain_files: 'Sequence[str] | None' = None, flag_duplicates: 'bool | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_path`, `output_path`, `reference_sequence`, `report_file`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2899.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `row_count`.
- **Note on the example:** the example data has no valid input for this tool (it needs data that cannot be made small: Illumina run folders, GTC files or flow-based reads, a trained or cluster model, or the Python gcnvkernel environment), so the example shows the input check that runs before GATK is started.

- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.annotation.gatk_p00_intervals_cli_functions import gatk_create_extended_illumina_manifest

work = Path(tempfile.mkdtemp())
try:
    gatk_create_extended_illumina_manifest(input_path="missing.dat", output_path=work / "o", reference_sequence="missing.dat", report_file=work / "o", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] input_path must be an existing file
```

## 10) `dnakit.annotation.gatk_gtc_to_vcf()`

- **Purpose:** Convert an Illumina GTC genotype file to a VCF with an extended manifest, bead pool manifest and cluster file.
- **Method:** Builds the fixed argument list of `gatk GtcToVcf` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk GtcToVcf --help`.
- **API:** `dnakit.annotation.gatk_gtc_to_vcf(input_path: 'str', output_path: 'str', extended_illumina_manifest: 'str', illumina_bead_pool_manifest_file: 'str', cluster_file: 'str', sample_alias: 'str', reference: 'str', *, expected_gender: 'str | None' = None, gender_gtc: 'str | None' = None, fingerprint_genotypes_vcf_file: 'str | None' = None, analysis_version_number: 'int | None' = None, pipeline_version: 'str | None' = None, do_not_allow_calls_on_zeroed_out_assays: 'bool | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_path`, `output_path`, `extended_illumina_manifest`, `illumina_bead_pool_manifest_file`, `cluster_file`, `sample_alias`, `reference`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2900.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `record_count`.
- **Note on the example:** the example data has no valid input for this tool (it needs data that cannot be made small: Illumina run folders, GTC files or flow-based reads, a trained or cluster model, or the Python gcnvkernel environment), so the example shows the input check that runs before GATK is started.

- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.annotation.gatk_p00_intervals_cli_functions import gatk_gtc_to_vcf

work = Path(tempfile.mkdtemp())
try:
    gatk_gtc_to_vcf(input_path="missing.dat", output_path=work / "o", extended_illumina_manifest="missing.dat", illumina_bead_pool_manifest_file="missing.dat", cluster_file="missing.dat", sample_alias="x", reference="missing.dat", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] input_path must be an existing file
```

## 11) `dnakit.annotation.gatk_merge_ped_into_vcf()`

- **Purpose:** Merge a single-sample PED file from zCall into a single-sample VCF.
- **Method:** Builds the fixed argument list of `gatk MergePedIntoVcf` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk MergePedIntoVcf --help`.
- **API:** `dnakit.annotation.gatk_merge_ped_into_vcf(original_vcf: 'str', ped_file: 'str', map_file: 'str', zcall_thresholds_file: 'str', zcall_version: 'str', output_path: 'str', *, reference: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `original_vcf`, `ped_file`, `map_file`, `zcall_thresholds_file`, `zcall_version`, `output_path`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2901.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `record_count`.
- **Note on the example:** the example data has no valid input for this tool (it needs data that cannot be made small: Illumina run folders, GTC files or flow-based reads, a trained or cluster model, or the Python gcnvkernel environment), so the example shows the input check that runs before GATK is started.

- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.annotation.gatk_p00_intervals_cli_functions import gatk_merge_ped_into_vcf

work = Path(tempfile.mkdtemp())
try:
    gatk_merge_ped_into_vcf(original_vcf="missing.dat", ped_file="missing.dat", map_file="missing.dat", zcall_thresholds_file="missing.dat", zcall_version="x", output_path=work / "o", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] original_vcf must be an existing file
```

## 12) `dnakit.annotation.gatk_vcf_to_adpc()`

- **Purpose:** Convert genotyping-array VCF files to an ADPC (adpc.bin) file for verifyIDintensity, with sample and marker-count files.
- **Method:** Builds the fixed argument list of `gatk VcfToAdpc` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk VcfToAdpc --help`.
- **API:** `dnakit.annotation.gatk_vcf_to_adpc(vcf_paths: 'Sequence[str]', output_path: 'str', samples_file: 'str', num_markers_file: 'str', *, reference: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `vcf_paths`, `output_path`, `samples_file`, `num_markers_file`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2902.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `num_markers`, `samples`.
- **Note on the example:** the example data has no valid input for this tool (it needs data that cannot be made small: Illumina run folders, GTC files or flow-based reads, a trained or cluster model, or the Python gcnvkernel environment), so the example shows the input check that runs before GATK is started.

- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.annotation.gatk_p00_intervals_cli_functions import gatk_vcf_to_adpc

work = Path(tempfile.mkdtemp())
try:
    gatk_vcf_to_adpc(vcf_paths=["missing.dat"], output_path=work / "o", samples_file=work / "o", num_markers_file=work / "o", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] vcf_paths must be an existing file
```

## 13) `dnakit.annotation.gatk_collect_arrays_variant_calling_metrics()`

- **Purpose:** Collect per-sample (detail) and aggregate (summary) metrics from a genotyping-array VCF: call rate, hom/het counts, dbSNP concordance.
- **Method:** Builds the fixed argument list of `gatk CollectArraysVariantCallingMetrics` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk CollectArraysVariantCallingMetrics --help`.
- **API:** `dnakit.annotation.gatk_collect_arrays_variant_calling_metrics(input_path: 'str', dbsnp: 'str', output_prefix: 'str', *, call_rate_pf_threshold: 'float | None' = None, sequence_dictionary: 'str | None' = None, num_processors: 'int | None' = None, reference: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_path`, `dbsnp`, `output_prefix`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2851.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `metrics_files`, `summary_metrics`, `detail_metrics`.
- **Note on the example:** the example data has no valid input for this tool (it needs data that cannot be made small: Illumina run folders, GTC files or flow-based reads, a trained or cluster model, or the Python gcnvkernel environment), so the example shows the input check that runs before GATK is started.

- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.annotation.gatk_p00_intervals_cli_functions import gatk_collect_arrays_variant_calling_metrics

work = Path(tempfile.mkdtemp())
try:
    gatk_collect_arrays_variant_calling_metrics(input_path="missing.dat", dbsnp="missing.dat", output_prefix=work / "o", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] input_path must be an existing file
```

## 14) `dnakit.annotation.gatk_add_flow_base_quality()`

- **Purpose:** Add a base-quality attribute to the reads of a flow-based SAM/BAM/CRAM file (Ultima-style flow data).
- **Method:** Builds the fixed argument list of `gatk AddFlowBaseQuality` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk AddFlowBaseQuality --help`. Flow-based tool for Ultima Genomics-style flow data (reads carry flow-order tags); GATK marks the flow tools experimental.
- **API:** `dnakit.annotation.gatk_add_flow_base_quality(input_paths: 'Sequence[str]', output_path: 'str', *, maximal_quality_score: 'int | None' = None, minimal_error_rate: 'float | None' = None, replace_quality_mode: 'bool | None' = None, keep_boundary_flows: 'bool | None' = None, reference: 'str | None' = None, intervals: 'Sequence[str] | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_paths`, `output_path`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2888.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`.
- **Note on the example:** the example data has no valid input for this tool (it needs data that cannot be made small: Illumina run folders, GTC files or flow-based reads, a trained or cluster model, or the Python gcnvkernel environment), so the example shows the input check that runs before GATK is started.

- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.annotation.gatk_p00_intervals_cli_functions import gatk_add_flow_base_quality

work = Path(tempfile.mkdtemp())
try:
    gatk_add_flow_base_quality(input_paths=["missing.dat"], output_path=work / "o", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] input_paths must be an existing file
```

## 15) `dnakit.annotation.gatk_add_flow_snv_quality()`

- **Purpose:** Add an SNV-quality attribute to the reads of a flow-based CRAM/BAM (Ultima-style flow data).
- **Method:** Builds the fixed argument list of `gatk AddFlowSNVQuality` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk AddFlowSNVQuality --help`. Flow-based tool for Ultima Genomics-style flow data (reads carry flow-order tags); GATK marks the flow tools experimental.
- **API:** `dnakit.annotation.gatk_add_flow_snv_quality(input_paths: 'Sequence[str]', output_path: 'str', *, snvq_mode: 'str | None' = None, output_quality_attribute: 'str | None' = None, max_phred_score: 'float | None' = None, keep_supplementary_alignments: 'bool | None' = None, include_qc_failed_read: 'bool | None' = None, reference: 'str | None' = None, intervals: 'Sequence[str] | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_paths`, `output_path`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2889.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`.
- **Note on the example:** the example data has no valid input for this tool (it needs data that cannot be made small: Illumina run folders, GTC files or flow-based reads, a trained or cluster model, or the Python gcnvkernel environment), so the example shows the input check that runs before GATK is started.

- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.annotation.gatk_p00_intervals_cli_functions import gatk_add_flow_snv_quality

work = Path(tempfile.mkdtemp())
try:
    gatk_add_flow_snv_quality(input_paths=["missing.dat"], output_path=work / "o", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] input_paths must be an existing file
```

## 16) `dnakit.annotation.gatk_calculate_average_combined_annotations()`

- **Purpose:** Divide annotations that GenomicsDB summed over samples by the number of samples to give their average.
- **Method:** Builds the fixed argument list of `gatk CalculateAverageCombinedAnnotations` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk CalculateAverageCombinedAnnotations --help`.
- **API:** `dnakit.annotation.gatk_calculate_average_combined_annotations(variant: 'str', summed_annotations_to_divide: 'Sequence[str]', output_path: 'str', *, intervals: 'Sequence[str] | None' = None, reference: 'str | None' = None, variant_output_filtering: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `variant`, `summed_annotations_to_divide`, `output_path`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2890.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `record_count`.
- **Note on the example:** the example data has no valid input for this tool (it needs data that cannot be made small: Illumina run folders, GTC files or flow-based reads, a trained or cluster model, or the Python gcnvkernel environment), so the example shows the input check that runs before GATK is started.

- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.annotation.gatk_p00_intervals_cli_functions import gatk_calculate_average_combined_annotations

work = Path(tempfile.mkdtemp())
try:
    gatk_calculate_average_combined_annotations(variant="missing.dat", summed_annotations_to_divide=["x"], output_path=work / "o", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] variant must be an existing file
```

## 17) `dnakit.annotation.gatk_flow_feature_mapper()`

- **Purpose:** Map features (initially SNVs) in a flow-based BAM and write them as a VCF.
- **Method:** Builds the fixed argument list of `gatk FlowFeatureMapper` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk FlowFeatureMapper --help`. Flow-based tool for Ultima Genomics-style flow data (reads carry flow-order tags); GATK marks the flow tools experimental.
- **API:** `dnakit.annotation.gatk_flow_feature_mapper(input_paths: 'Sequence[str]', output_path: 'str', *, mapping_feature: 'str | None' = None, copy_attrs: 'Sequence[str] | None' = None, copy_attr_prefix: 'str | None' = None, min_score: 'float | None' = None, max_score: 'float | None' = None, limit_score: 'float | None' = None, snv_identical_bases: 'int | None' = None, snv_identical_bases_after: 'int | None' = None, include_dup_reads: 'bool | None' = None, include_qc_failed_reads: 'bool | None' = None, keep_negatives: 'bool | None' = None, exclude_nan_scores: 'bool | None' = None, keep_supplementary_alignments: 'bool | None' = None, emit_ref_confidence: 'str | None' = None, reference: 'str | None' = None, intervals: 'Sequence[str] | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_paths`, `output_path`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2891.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `record_count`.
- **Note on the example:** the example data has no valid input for this tool (it needs data that cannot be made small: Illumina run folders, GTC files or flow-based reads, a trained or cluster model, or the Python gcnvkernel environment), so the example shows the input check that runs before GATK is started.

- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.annotation.gatk_p00_intervals_cli_functions import gatk_flow_feature_mapper

work = Path(tempfile.mkdtemp())
try:
    gatk_flow_feature_mapper(input_paths=["missing.dat"], output_path=work / "o", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] input_paths must be an existing file
```

## 18) `dnakit.annotation.gatk_flow_pair_hmm_align_reads_to_haplotypes()`

- **Purpose:** Produce a read x haplotype matrix with the log-likelihood of every read against every haplotype using the flow-based PairHMM.
- **Method:** Builds the fixed argument list of `gatk FlowPairHMMAlignReadsToHaplotypes` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk FlowPairHMMAlignReadsToHaplotypes --help`. Flow-based tool for Ultima Genomics-style flow data (reads carry flow-order tags); GATK marks the flow tools experimental.
- **API:** `dnakit.annotation.gatk_flow_pair_hmm_align_reads_to_haplotypes(input_paths: 'Sequence[str]', haplotypes: 'str', output_path: 'str', *, aligner: 'str | None' = None, likelihood_calculation_engine: 'str | None' = None, ref_haplotype: 'str | None' = None, concise_output_format: 'bool | None' = None, reference: 'str | None' = None, intervals: 'Sequence[str] | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_paths`, `haplotypes`, `output_path`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2892.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `matrix_written`.
- **Note on the example:** the example data has no valid input for this tool (it needs data that cannot be made small: Illumina run folders, GTC files or flow-based reads, a trained or cluster model, or the Python gcnvkernel environment), so the example shows the input check that runs before GATK is started.

- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.annotation.gatk_p00_intervals_cli_functions import gatk_flow_pair_hmm_align_reads_to_haplotypes

work = Path(tempfile.mkdtemp())
try:
    gatk_flow_pair_hmm_align_reads_to_haplotypes(input_paths=["missing.dat"], haplotypes="missing.dat", output_path=work / "o", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] input_paths must be an existing file
```

## 19) `dnakit.annotation.gatk_ground_truth_reads_builder()`

- **Purpose:** Produce a flexible and robust ground-truth set of reads for base-calling training from reads and the maternal and paternal haplotype references (CSV).
- **Method:** Builds the fixed argument list of `gatk GroundTruthReadsBuilder` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk GroundTruthReadsBuilder --help`. Flow-based tool for Ultima Genomics-style flow data (reads carry flow-order tags); GATK marks the flow tools experimental.
- **API:** `dnakit.annotation.gatk_ground_truth_reads_builder(input_paths: 'Sequence[str]', maternal_ref: 'str', paternal_ref: 'str', ancestral_translators_base_path: 'str', output_csv: 'str', *, subsampling_ratio: 'float | None' = None, max_output_reads: 'int | None' = None, min_haplotype_score: 'float | None' = None, min_mq: 'float | None' = None, include_supp_align: 'bool | None' = None, fill_trimmed_reads: 'bool | None' = None, reference: 'str | None' = None, intervals: 'Sequence[str] | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_paths`, `maternal_ref`, `paternal_ref`, `ancestral_translators_base_path`, `output_csv`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2893.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`.
- **Note on the example:** the example data has no valid input for this tool (it needs data that cannot be made small: Illumina run folders, GTC files or flow-based reads, a trained or cluster model, or the Python gcnvkernel environment), so the example shows the input check that runs before GATK is started.

- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.annotation.gatk_p00_intervals_cli_functions import gatk_ground_truth_reads_builder

work = Path(tempfile.mkdtemp())
try:
    gatk_ground_truth_reads_builder(input_paths=["missing.dat"], maternal_ref="missing.dat", paternal_ref="missing.dat", ancestral_translators_base_path=work / "o", output_csv=work / "o", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] input_paths must be an existing file
```

## 20) `dnakit.annotation.gatk_ground_truth_scorer()`

- **Purpose:** Score flow-based reads against a reference or ground truth and write the per-read CSV and an optional report.
- **Method:** Builds the fixed argument list of `gatk GroundTruthScorer` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk GroundTruthScorer --help`. Flow-based tool for Ultima Genomics-style flow data (reads carry flow-order tags); GATK marks the flow tools experimental.
- **API:** `dnakit.annotation.gatk_ground_truth_scorer(input_paths: 'Sequence[str]', output_csv: 'str', *, report_file: 'str | None' = None, features_file: 'str | None' = None, genome_prior: 'str | None' = None, add_mean_call: 'bool | None' = None, exclude_zero_flows: 'bool | None' = None, omit_zeros_from_report: 'bool | None' = None, quality_percentiles: 'str | None' = None, use_softclipped_bases: 'bool | None' = None, normalized_score_threshold: 'float | None' = None, reference: 'str | None' = None, intervals: 'Sequence[str] | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `input_paths`, `output_csv`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2894.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`.
- **Note on the example:** the example data has no valid input for this tool (it needs data that cannot be made small: Illumina run folders, GTC files or flow-based reads, a trained or cluster model, or the Python gcnvkernel environment), so the example shows the input check that runs before GATK is started.

- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.annotation.gatk_p00_intervals_cli_functions import gatk_ground_truth_scorer

work = Path(tempfile.mkdtemp())
try:
    gatk_ground_truth_scorer(input_paths=["missing.dat"], output_csv=work / "o", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] input_paths must be an existing file
```

## 21) `dnakit.annotation.gatk_pathseq_build_kmers()`

- **Purpose:** Build the set of host reference k-mers (hash set, or Bloom filter when a false-positive probability is given) used by PathSeq to subtract host reads.
- **Method:** Builds the fixed argument list of `gatk PathSeqBuildKmers` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk PathSeqBuildKmers --help`.
- **API:** `dnakit.annotation.gatk_pathseq_build_kmers(reference: 'str', output_path: 'str', *, kmer_size: 'int | None' = None, kmer_mask: 'str | None' = None, kmer_spacing: 'int | None' = None, bloom_false_positive_probability: 'float | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `reference`, `output_path`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2910.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`, `kmer_files`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.annotation.gatk_p00_intervals_cli_functions import gatk_pathseq_build_kmers

work = Path(tempfile.mkdtemp())
result = gatk_pathseq_build_kmers("toy.fa", work / "host", kmer_size=15, verbose=False)
print([path.rsplit("/", 1)[1] for path in result["kmer_files"]])
```

- **Example output:**

```text
['host.hss']
```

## 22) `dnakit.annotation.gatk_pathseq_build_reference_taxonomy()`

- **Purpose:** Build the taxonomy data file (.db) of the microbe reference for PathSeq from the NCBI taxonomy dump and the reference FASTA.
- **Method:** Builds the fixed argument list of `gatk PathSeqBuildReferenceTaxonomy` from the typed parameters (no shell), runs it in the folder of the output with the version recorded, and requires exit code 0 and the expected output files. Other long options of the tool go through `options` and are checked against `gatk PathSeqBuildReferenceTaxonomy --help`.
- **API:** `dnakit.annotation.gatk_pathseq_build_reference_taxonomy(reference: 'str', tax_dump: 'str', output_path: 'str', *, refseq_catalog: 'str | None' = None, genbank_catalog: 'str | None' = None, min_non_virus_contig_length: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** The `gatk` launcher (GATK 4.6.2.0 or later, Java 17) on PATH, in `DNAKIT_GATK`, or given as `gatk_path`; without it `BackendUnavailableError` is raised with the install hint.
- **Input:** Required: `reference`, `tax_dump`, `output_path`. Every optional parameter defaults to None, which keeps the GATK default. Output paths must not exist yet. Covers work-list rows G2911.
- **Output:** The standard result dictionary (`status`, `tool`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance`) plus: `outputs`.
- **Note on the example:** the example data has no valid input for this tool (it needs data that cannot be made small: Illumina run folders, GTC files or flow-based reads, a trained or cluster model, or the Python gcnvkernel environment), so the example shows the input check that runs before GATK is started.

- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.annotation.gatk_p00_intervals_cli_functions import gatk_pathseq_build_reference_taxonomy

work = Path(tempfile.mkdtemp())
try:
    gatk_pathseq_build_reference_taxonomy(reference="missing.dat", tax_dump="missing.dat", output_path=work / "o", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] reference must be an existing file
```
