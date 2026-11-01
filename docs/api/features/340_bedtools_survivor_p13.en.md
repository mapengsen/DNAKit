# BEDTools interval operations and SURVIVOR structural-variant tools

This group adds command-line facades for [BEDTools](https://github.com/arq5x/bedtools2) (`window`, `random`, `shuffle`, `sample`, `multiinter`, `unionbedg`, `pairtopair`, `multicov`, `reldist`, `coverage` histogram and depth, `bamtobed -bedpe`, `split`) in `dnakit.annotation`, and for [SURVIVOR](https://github.com/fritzsedlazeck/SURVIVOR) (simulation, evaluation, filtering, statistics, MUMmer annotation and format conversion) in `dnakit.variants`. Rows of the coverage work list that differ only by a parameter are one function: `bedtools_window(report="count")`, `bedtools_reldist(detail=True)`, `bedtools_coverage_profile(report="depth")`, `bedtools_split(algorithm="simple")` and `survivor_convert_to_sv_vcf(input_format=...)`. SURVIVOR `select_sample` is not wrapped because it segfaults in 1.0.7.

Install the programs separately: `conda install -c bioconda bedtools survivor` (tested with BEDTools 2.27.1 and 2.31.1, SURVIVOR 1.0.7). The functions import without them and raise `BackendUnavailableError` with this hint when a program is missing.

## Example data

**Example file:** `a.bed`

```text
chr1	100	200	a1	5	+
chr1	500	600	a2	6	-
chr2	10	90	a3	7	+
```

**Example file:** `b.bed`

```text
chr1	150	250	b1	1	+
chr1	900	950	b2	2	-
chr2	100	150	b3	3	+
```

**Example file:** `sv.bed`

```text
chr1	1000	2000
chr1	5000	5600
```

## 1) `dnakit.annotation.bedtools_window()`

- **Purpose:** Report the B features that fall in a window around each A interval (``bedtools window``).
- **Method:** Runs the fixed external command or port described in the docstring with a whitelisted argument list (no shell); output files must be new and are never overwritten. A missing program raises `BackendUnavailableError` with the install hint.
- **API:** `dnakit.annotation.bedtools_window(a_path, b_path, output_path, *, window=None, left=None, right=None, strand_aware_window=False, strand=None, report='overlaps', header=False, bedtools_path=None, timeout_seconds=3600.0)`.
- **Input:** Existing input files and the new output path(s) named in the signature; every optional argument is passed to the program only when given, so its own defaults apply.
- **Output:** A JSON-serialisable dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and the parsed summary keys described in the docstring.

## 2) `dnakit.annotation.bedtools_random()`

- **Purpose:** Generate random intervals of a fixed length across a genome (``bedtools random``).
- **Method:** Runs the fixed external command or port described in the docstring with a whitelisted argument list (no shell); output files must be new and are never overwritten. A missing program raises `BackendUnavailableError` with the install hint.
- **API:** `dnakit.annotation.bedtools_random(genome_path, output_path, *, length=100, count=1000000, seed=None, bedtools_path=None, timeout_seconds=3600.0)`.
- **Input:** Existing input files and the new output path(s) named in the signature; every optional argument is passed to the program only when given, so its own defaults apply.
- **Output:** A JSON-serialisable dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and the parsed summary keys described in the docstring.

## 3) `dnakit.annotation.bedtools_shuffle()`

- **Purpose:** Randomly permute the genomic locations of intervals (``bedtools shuffle``).
- **Method:** Runs the fixed external command or port described in the docstring with a whitelisted argument list (no shell); output files must be new and are never overwritten. A missing program raises `BackendUnavailableError` with the install hint.
- **API:** `dnakit.annotation.bedtools_shuffle(intervals_path, genome_path, output_path, *, exclude_path=None, include_path=None, same_chromosome=False, seed=None, max_exclude_overlap=None, choose_chromosome_first=False, bedpe=False, max_tries=None, no_overlapping=False, allow_beyond_chrom_end=False, bedtools_path=None, timeout_seconds=3600.0)`.
- **Input:** Existing input files and the new output path(s) named in the signature; every optional argument is passed to the program only when given, so its own defaults apply.
- **Output:** A JSON-serialisable dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and the parsed summary keys described in the docstring.

## 4) `dnakit.annotation.bedtools_sample()`

- **Purpose:** Take a random sample of records from a BED/GFF/VCF or BAM file by reservoir sampling (``bedtools sample``).
- **Method:** Runs the fixed external command or port described in the docstring with a whitelisted argument list (no shell); output files must be new and are never overwritten. A missing program raises `BackendUnavailableError` with the install hint.
- **API:** `dnakit.annotation.bedtools_sample(input_path, output_path, *, count=1000000, seed=None, strand=None, header=False, bam_as_bed=False, bedtools_path=None, timeout_seconds=3600.0)`.
- **Input:** Existing input files and the new output path(s) named in the signature; every optional argument is passed to the program only when given, so its own defaults apply.
- **Output:** A JSON-serialisable dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and the parsed summary keys described in the docstring.

## 5) `dnakit.annotation.bedtools_multiinter()`

- **Purpose:** Identify the common intervals among several sorted BED/GFF/VCF files (``bedtools multiinter``).
- **Method:** Runs the fixed external command or port described in the docstring with a whitelisted argument list (no shell); output files must be new and are never overwritten. A missing program raises `BackendUnavailableError` with the install hint.
- **API:** `dnakit.annotation.bedtools_multiinter(input_paths, output_path, *, names=None, header=False, cluster=False, genome_path=None, report_empty=False, filler=None, bedtools_path=None, timeout_seconds=3600.0)`.
- **Input:** Existing input files and the new output path(s) named in the signature; every optional argument is passed to the program only when given, so its own defaults apply.
- **Output:** A JSON-serialisable dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and the parsed summary keys described in the docstring.

## 6) `dnakit.annotation.bedtools_unionbedg()`

- **Purpose:** Combine several sorted bedGraph files into one table with a value column per file (``bedtools unionbedg``).
- **Method:** Runs the fixed external command or port described in the docstring with a whitelisted argument list (no shell); output files must be new and are never overwritten. A missing program raises `BackendUnavailableError` with the install hint.
- **API:** `dnakit.annotation.bedtools_unionbedg(input_paths, output_path, *, names=None, header=False, genome_path=None, report_empty=False, filler=None, bedtools_path=None, timeout_seconds=3600.0)`.
- **Input:** Existing input files and the new output path(s) named in the signature; every optional argument is passed to the program only when given, so its own defaults apply.
- **Output:** A JSON-serialisable dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and the parsed summary keys described in the docstring.

## 7) `dnakit.annotation.bedtools_pairtopair()`

- **Purpose:** Report overlaps between two paired-end BED (BEDPE) files (``bedtools pairtopair``).
- **Method:** Runs the fixed external command or port described in the docstring with a whitelisted argument list (no shell); output files must be new and are never overwritten. A missing program raises `BackendUnavailableError` with the install hint.
- **API:** `dnakit.annotation.bedtools_pairtopair(a_path, b_path, output_path, *, min_overlap=None, overlap_type='both', slop=0, strand_aware_slop=False, ignore_strand=False, require_different_names=False, bedtools_path=None, timeout_seconds=3600.0)`.
- **Input:** Existing input files and the new output path(s) named in the signature; every optional argument is passed to the program only when given, so its own defaults apply.
- **Output:** A JSON-serialisable dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and the parsed summary keys described in the docstring.

## 8) `dnakit.annotation.bedtools_multicov()`

- **Purpose:** Count the alignments of several BAM files that overlap each interval (``bedtools multicov``).
- **Method:** Runs the fixed external command or port described in the docstring with a whitelisted argument list (no shell); output files must be new and are never overwritten. A missing program raises `BackendUnavailableError` with the install hint.
- **API:** `dnakit.annotation.bedtools_multicov(bam_paths, intervals_path, output_path, *, split=False, strand=None, min_overlap=None, reciprocal=False, min_mapping_quality=None, include_duplicates=False, include_failed_qc=False, proper_pairs_only=False, bedtools_path=None, timeout_seconds=3600.0)`.
- **Input:** Existing input files and the new output path(s) named in the signature; every optional argument is passed to the program only when given, so its own defaults apply.
- **Output:** A JSON-serialisable dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and the parsed summary keys described in the docstring.

## 9) `dnakit.annotation.bedtools_reldist()`

- **Purpose:** Distribution of the relative distances between two interval sets (``bedtools reldist``).
- **Method:** Runs the fixed external command or port described in the docstring with a whitelisted argument list (no shell); output files must be new and are never overwritten. A missing program raises `BackendUnavailableError` with the install hint.
- **API:** `dnakit.annotation.bedtools_reldist(a_path, b_path, output_path, *, detail=False, bedtools_path=None, timeout_seconds=3600.0)`.
- **Input:** Existing input files and the new output path(s) named in the signature; every optional argument is passed to the program only when given, so its own defaults apply.
- **Output:** A JSON-serialisable dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and the parsed summary keys described in the docstring.

## 10) `dnakit.annotation.bedtools_coverage_profile()`

- **Purpose:** Coverage of A intervals by B features as a histogram or as per-position depth (``bedtools coverage``).
- **Method:** Runs the fixed external command or port described in the docstring with a whitelisted argument list (no shell); output files must be new and are never overwritten. A missing program raises `BackendUnavailableError` with the install hint.
- **API:** `dnakit.annotation.bedtools_coverage_profile(a_path, b_path, output_path, *, report='histogram', strand=None, min_overlap_a=None, min_overlap_b=None, reciprocal=False, either_fraction=False, split=False, sorted_input=False, genome_path=None, header=False, bedtools_path=None, timeout_seconds=3600.0)`.
- **Input:** Existing input files and the new output path(s) named in the signature; every optional argument is passed to the program only when given, so its own defaults apply.
- **Output:** A JSON-serialisable dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and the parsed summary keys described in the docstring.

## 11) `dnakit.annotation.bedtools_bam_to_bedpe()`

- **Purpose:** Convert paired-end BAM alignments to BEDPE (``bedtools bamtobed -bedpe``).
- **Method:** Runs the fixed external command or port described in the docstring with a whitelisted argument list (no shell); output files must be new and are never overwritten. A missing program raises `BackendUnavailableError` with the install hint.
- **API:** `dnakit.annotation.bedtools_bam_to_bedpe(bam_path, output_path, *, mate1_first=False, edit_distance_score=False, bedtools_path=None, timeout_seconds=3600.0)`.
- **Input:** Existing input files and the new output path(s) named in the signature; every optional argument is passed to the program only when given, so its own defaults apply.
- **Output:** A JSON-serialisable dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and the parsed summary keys described in the docstring.

## 12) `dnakit.annotation.bedtools_split()`

- **Purpose:** Split a BED file into ``number_of_files`` files inside a fresh directory (``bedtools split``).
- **Method:** Runs the fixed external command or port described in the docstring with a whitelisted argument list (no shell); output files must be new and are never overwritten. A missing program raises `BackendUnavailableError` with the install hint.
- **API:** `dnakit.annotation.bedtools_split(input_path, output_directory, number_of_files, *, prefix='split', algorithm='size', bedtools_path=None, timeout_seconds=3600.0)`.
- **Input:** Existing input files and the new output path(s) named in the signature; every optional argument is passed to the program only when given, so its own defaults apply.
- **Output:** A JSON-serialisable dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and the parsed summary keys described in the docstring.

## 13) `dnakit.variants.survivor_simulate_sv_genome()`

- **Purpose:** Simulate SVs and SNPs into a reference genome (``SURVIVOR simSV``).
- **Method:** Runs the fixed external command or port described in the docstring with a whitelisted argument list (no shell); output files must be new and are never overwritten. A missing program raises `BackendUnavailableError` with the install hint.
- **API:** `dnakit.variants.survivor_simulate_sv_genome(reference_path, output_prefix, *, duplication=(100, 10000, 3), indel=(20, 500, 1), translocation=(1000, 3000, 2), inversion=(600, 800, 4), inversion_deletion=(600, 800, 2), inversion_duplication=(600, 800, 2), snp_frequency=0.0, coordinates_for='simulated_reads', parameter_file=None, survivor_path=None, timeout_seconds=3600.0)`.
- **Input:** Existing input files and the new output path(s) named in the signature; every optional argument is passed to the program only when given, so its own defaults apply.
- **Output:** A JSON-serialisable dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and the parsed summary keys described in the docstring.

## 14) `dnakit.variants.survivor_scan_read_error_profile()`

- **Purpose:** Extract a per-position error profile from aligned long reads (``SURVIVOR scanreads``).
- **Method:** Runs the fixed external command or port described in the docstring with a whitelisted argument list (no shell); output files must be new and are never overwritten. A missing program raises `BackendUnavailableError` with the install hint.
- **API:** `dnakit.variants.survivor_scan_read_error_profile(alignments_path, output_path, *, min_read_length=1000, reference_path=None, survivor_path=None, samtools_path=None, timeout_seconds=3600.0)`.
- **Input:** Existing input files and the new output path(s) named in the signature; every optional argument is passed to the program only when given, so its own defaults apply.
- **Output:** A JSON-serialisable dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and the parsed summary keys described in the docstring.

## 15) `dnakit.variants.survivor_simulate_long_reads()`

- **Purpose:** Simulate PacBio- or ONT-like long reads from a reference (``SURVIVOR simreads``).
- **Method:** Runs the fixed external command or port described in the docstring with a whitelisted argument list (no shell); output files must be new and are never overwritten. A missing program raises `BackendUnavailableError` with the install hint.
- **API:** `dnakit.variants.survivor_simulate_long_reads(reference_path, error_profile_path, output_path, *, coverage=2.0, survivor_path=None, timeout_seconds=3600.0)`.
- **Input:** Existing input files and the new output path(s) named in the signature; every optional argument is passed to the program only when given, so its own defaults apply.
- **Output:** A JSON-serialisable dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and the parsed summary keys described in the docstring.

## 16) `dnakit.variants.survivor_evaluate_sv_calls()`

- **Purpose:** Evaluate an SV call set against the events written by :func:`survivor_simulate_sv_genome` (``SURVIVOR eval``).
- **Method:** Runs the fixed external command or port described in the docstring with a whitelisted argument list (no shell); output files must be new and are never overwritten. A missing program raises `BackendUnavailableError` with the install hint.
- **API:** `dnakit.variants.survivor_evaluate_sv_calls(vcf_path, simulated_bed_path, output_prefix, *, max_distance=10, survivor_path=None, timeout_seconds=3600.0)`.
- **Input:** Existing input files and the new output path(s) named in the signature; every optional argument is passed to the program only when given, so its own defaults apply.
- **Output:** A JSON-serialisable dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and the parsed summary keys described in the docstring.

## 17) `dnakit.variants.survivor_filter_sv_vcf()`

- **Purpose:** Filter an SV VCF by size, regions to ignore, allele frequency and read support (``SURVIVOR filter``).
- **Method:** Runs the fixed external command or port described in the docstring with a whitelisted argument list (no shell); output files must be new and are never overwritten. A missing program raises `BackendUnavailableError` with the install hint.
- **API:** `dnakit.variants.survivor_filter_sv_vcf(vcf_path, output_path, *, ignore_regions_bed=None, min_size=None, max_size=None, min_allele_frequency=0.0, min_read_support=None, survivor_path=None, timeout_seconds=3600.0)`.
- **Input:** Existing input files and the new output path(s) named in the signature; every optional argument is passed to the program only when given, so its own defaults apply.
- **Output:** A JSON-serialisable dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and the parsed summary keys described in the docstring.

## 18) `dnakit.variants.survivor_sv_vcf_stats()`

- **Purpose:** Summary statistics of an SV VCF: counts by type and by length class (``SURVIVOR stats``).
- **Method:** Runs the fixed external command or port described in the docstring with a whitelisted argument list (no shell); output files must be new and are never overwritten. A missing program raises `BackendUnavailableError` with the install hint.
- **API:** `dnakit.variants.survivor_sv_vcf_stats(vcf_path, summary_path, *, min_size=None, max_size=None, min_read_support=None, survivor_path=None, timeout_seconds=3600.0)`.
- **Input:** Existing input files and the new output path(s) named in the signature; every optional argument is passed to the program only when given, so its own defaults apply.
- **Output:** A JSON-serialisable dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and the parsed summary keys described in the docstring.

## 19) `dnakit.variants.survivor_annotate_sv_with_mummer()`

- **Purpose:** Annotate an SV VCF with the breakpoints found by MUMmer (``SURVIVOR compMUMMer``).
- **Method:** Runs the fixed external command or port described in the docstring with a whitelisted argument list (no shell); output files must be new and are never overwritten. A missing program raises `BackendUnavailableError` with the install hint.
- **API:** `dnakit.variants.survivor_annotate_sv_with_mummer(vcf_path, show_diff_paths, output_path, *, max_distance=100, survivor_path=None, timeout_seconds=3600.0)`.
- **Input:** Existing input files and the new output path(s) named in the signature; every optional argument is passed to the program only when given, so its own defaults apply.
- **Output:** A JSON-serialisable dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and the parsed summary keys described in the docstring.

## 20) `dnakit.variants.survivor_sv_vcf_to_bed()`

- **Purpose:** Convert an SV VCF to a BEDPE-style BED file with an optional size range (``SURVIVOR vcftobed``).
- **Method:** Runs the fixed external command or port described in the docstring with a whitelisted argument list (no shell); output files must be new and are never overwritten. A missing program raises `BackendUnavailableError` with the install hint.
- **API:** `dnakit.variants.survivor_sv_vcf_to_bed(vcf_path, output_path, *, min_size=None, max_size=None, survivor_path=None, timeout_seconds=3600.0)`.
- **Input:** Existing input files and the new output path(s) named in the signature; every optional argument is passed to the program only when given, so its own defaults apply.
- **Output:** A JSON-serialisable dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and the parsed summary keys described in the docstring.

## 21) `dnakit.variants.survivor_convert_to_sv_vcf()`

- **Purpose:** Convert SV calls of another program or format to an SV VCF (``SURVIVOR bedtovcf/smaptovcf/bedpetovcf/convertAssemblytics``).
- **Method:** Runs the fixed external command or port described in the docstring with a whitelisted argument list (no shell); output files must be new and are never overwritten. A missing program raises `BackendUnavailableError` with the install hint.
- **API:** `dnakit.variants.survivor_convert_to_sv_vcf(input_path, output_path, *, input_format='bed', sv_type=None, min_size=0, survivor_path=None, timeout_seconds=3600.0)`.
- **Input:** Existing input files and the new output path(s) named in the signature; every optional argument is passed to the program only when given, so its own defaults apply.
- **Output:** A JSON-serialisable dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and the parsed summary keys described in the docstring.

## 22) `dnakit.variants.survivor_hapcut_to_phased_vcf()`

- **Purpose:** Write a phased VCF from a HapCUT2 final haplotype file and the SNP VCF given to HapCUT2 (``SURVIVOR hapcuttovcf``).
- **Method:** Runs the fixed external command or port described in the docstring with a whitelisted argument list (no shell); output files must be new and are never overwritten. A missing program raises `BackendUnavailableError` with the install hint.
- **API:** `dnakit.variants.survivor_hapcut_to_phased_vcf(snp_vcf_path, hapcut_path, output_path, *, survivor_path=None, timeout_seconds=3600.0)`.
- **Input:** Existing input files and the new output path(s) named in the signature; every optional argument is passed to the program only when given, so its own defaults apply.
- **Output:** A JSON-serialisable dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and the parsed summary keys described in the docstring.

## 23) `dnakit.variants.survivor_low_mapping_quality_regions()`

- **Purpose:** Cluster a MQ0 coverage table into BED regions to filter SVs in poorly mappable sequence (``SURVIVOR bincov``).
- **Method:** Runs the fixed external command or port described in the docstring with a whitelisted argument list (no shell); output files must be new and are never overwritten. A missing program raises `BackendUnavailableError` with the install hint.
- **API:** `dnakit.variants.survivor_low_mapping_quality_regions(coverage_path, output_path, *, border_size=10, min_coverage=2, survivor_path=None, timeout_seconds=3600.0)`.
- **Input:** Existing input files and the new output path(s) named in the signature; every optional argument is passed to the program only when given, so its own defaults apply.
- **Output:** A JSON-serialisable dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and the parsed summary keys described in the docstring.

## 24) `dnakit.variants.survivor_pairwise_overlap_matrix()`

- **Purpose:** Pairwise SV overlap matrix across the samples of a merged multi-sample VCF (``SURVIVOR genComp``).
- **Method:** Runs the fixed external command or port described in the docstring with a whitelisted argument list (no shell); output files must be new and are never overwritten. A missing program raises `BackendUnavailableError` with the install hint.
- **API:** `dnakit.variants.survivor_pairwise_overlap_matrix(merged_vcf_path, output_path, *, normalize=False, survivor_path=None, timeout_seconds=3600.0)`.
- **Input:** Existing input files and the new output path(s) named in the signature; every optional argument is passed to the program only when given, so its own defaults apply.
- **Output:** A JSON-serialisable dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and the parsed summary keys described in the docstring.

## 25) `dnakit.variants.survivor_parental_phasing()`

- **Purpose:** Assign the SNPs of a proband to the parental haplotypes for SV phasing (``SURVIVOR parent_phasing``).
- **Method:** Runs the fixed external command or port described in the docstring with a whitelisted argument list (no shell); output files must be new and are never overwritten. A missing program raises `BackendUnavailableError` with the install hint.
- **API:** `dnakit.variants.survivor_parental_phasing(merged_parental_vcf_path, hapcut_path, gatk_path, proband_snp_vcf_path, output_path, *, survivor_path=None, timeout_seconds=3600.0)`.
- **Input:** Existing input files and the new output path(s) named in the signature; every optional argument is passed to the program only when given, so its own defaults apply.
- **Output:** A JSON-serialisable dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail` and the parsed summary keys described in the docstring.

## Examples

### Example 1: Count B features near A intervals and get per-base depth

- **Example code:**

```python
import tempfile
from pathlib import Path
from dnakit.annotation import bedtools_p13_cli_functions as bt

out = Path(tempfile.mkdtemp())
result = bt.bedtools_window("a.bed", "b.bed", str(out / "w.bed"), window=100, report="count", verbose=False)
print(result["status"], result["record_count"])
print((out / "w.bed").read_text().splitlines()[0])
depth = bt.bedtools_coverage_profile("a.bed", "b.bed", str(out / "d.txt"), report="depth", verbose=False)
print(depth["report"], depth["record_count"])
```

- **Example output:**

```text
completed 3
chr1	100	200	a1	5	+	1
depth 280
```

### Example 2: BED to SV VCF, then filter and summarise

- **Example code:**

```python
import tempfile
from pathlib import Path
from dnakit.variants import survivor_p13_cli_functions as sv

out = Path(tempfile.mkdtemp())
conv = sv.survivor_convert_to_sv_vcf("sv.bed", str(out / "sv.vcf"), input_format="bed", sv_type="DEL", verbose=False)
print(conv["status"], conv["record_count"])
flt = sv.survivor_filter_sv_vcf(str(out / "sv.vcf"), str(out / "big.vcf"), min_size=800, verbose=False)
print(flt["records_in"], flt["records_out"])
stats = sv.survivor_sv_vcf_stats(str(out / "sv.vcf"), str(out / "stats.txt"), verbose=False)
print(stats["totals"])
```

- **Example output:**

```text
completed 2
2 1
{'Tot': 2, 'DEL': 2, 'DUP': 0, 'INS': 0, 'INV': 0, 'TRA': 0}
```
