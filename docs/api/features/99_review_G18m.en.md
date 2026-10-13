# Review batch G18m: other rows of methylpy, MACS3 and deepTools

This page records how the G18m review batch handled the `other` rows of the ledgers for methylpy, MACS3 and deepTools. The three packages had 47 `other` rows in total, and every one now has a final status; no `other` or `cap` row remains. 26 functions are wrapped: 12 from methylpy, 10 from MACS3 and 4 from deepTools. The rest are marked `needs-binary`, `deprecated`, `internal`, `not-dna`, `plot`, `duplicate-of` or `family-variant`. Each package section lists its unwrapped functions with the reason under "Not wrapped".

All wrapped functions live in the `*_review_functions.py` modules of `dnakit.epigenome`. They call the upstream Python APIs directly, run no command line and use no network. The examples import from the module files, and their outputs come from actually running the code.

## Additional methylpy functions (DMR steps, allc splitting, merging and binomial tests)

This section adds 12 methylpy functions in the new module `methylpy_review_functions.py` of `dnakit.epigenome`. The examples import from the module file; once exported from `dnakit.epigenome`, the same functions are callable as `dnakit.epigenome.<name>()`. Every function runs with one process, none uses methylpy's bundled permutation program, and none needs samtools, bowtie2, tabix or bgzip.

The methylpy functions behind them sort input files in place, hand paths to the shell commands `sort`, `wc` and `rm`, leave temporary files behind after errors, and fail with low-level errors on malformed tables. The wrappers check the inputs first (types, ranges, header, field count, contiguous chromosomes, existing output paths), never overwrite an existing file, let methylpy work only on copies or new paths, and turn methylpy exits and low-level errors into `ConfigurationError`. The return value is a dict with `status`, `operation` and `provenance`.

The examples run in order in one directory, and later examples use the `methylpy_rms.tsv.gz` that the first one creates. Output files go to a temporary directory, so the examples can be run repeatedly.

### Example data

RMS permutation-test results for samples A and B (uncompressed; `num_sims=100`, `num_sig_tests=5`). Real tables come from methylpy's `run_rms_tests.out` and are gzip-compressed; the numbers here are hand-made.

**Example file:** `methylpy_rms.tsv`

```text
chr	pos	strand	mc_class	pvalue	mc_A	mc_B	h_A	h_B	frac_A	frac_B	mc_residual_A	mc_residual_B	uc_residual_A	uc_residual_B	num_simulations_sig	num_simulations_run
1	100	+	CGN	0.01	9	1	10	10	0.900000	0.100000	3.577709	-3.577709	-3.577709	3.577709	1	100
1	120	+	CGN	0.01	10	0	10	10	1.000000	0.000000	4.472136	-4.472136	-4.472136	4.472136	1	100
1	140	+	CGN	0.02	8	1	10	10	0.800000	0.100000	3.146266	-3.146266	-3.146266	3.146266	2	100
1	300	+	CGN	0.50	5	6	10	10	0.500000	0.600000	-0.449467	0.449467	0.449467	-0.449467	5	10
1	320	+	CGN	0.42	4	5	10	10	0.400000	0.500000	-0.449467	0.449467	0.449467	-0.449467	5	12
1	340	+	CGN	0.62	3	4	10	10	0.300000	0.400000	-0.468807	0.468807	0.468807	-0.468807	5	8
1	600	+	CGN	0.01	1	9	10	10	0.100000	0.900000	-3.577709	3.577709	3.577709	-3.577709	1	100
1	620	+	CGN	0.03	0	10	10	10	0.000000	1.000000	-4.472136	4.472136	4.472136	-4.472136	3	100
2	50	+	CGN	0.50	5	5	10	10	0.500000	0.500000	0.000000	0.000000	0.000000	0.000000	5	10
2	80	+	CGN	0.01	9	2	10	10	0.900000	0.200000	3.146266	-3.146266	-3.146266	3.146266	1	100
```

The allc table of sample A (chromosomes are contiguous).

**Example file:** `methylpy_allc_A.tsv`

```text
chr1	100	+	CGA	9	10	1
chr1	120	+	CGT	10	10	1
chr1	125	+	CHH	1	4	0
chr1	140	+	CGC	8	10	1
chr1	500	+	CGG	5	10	1
chr2	30	+	CGA	0	9	0
```

The allc table of sample B.

**Example file:** `methylpy_allc_B.tsv`

```text
chr1	100	+	CGA	1	10	0
chr1	120	+	CGT	0	10	0
chr1	140	+	CGC	1	10	0
chr1	500	+	CGG	5	10	1
chr2	30	+	CGA	3	6	1
```

A DMR table with four regions (sorted by start within a chromosome; the second overlaps the first).

**Example file:** `methylpy_dmrs.tsv`

```text
#chr	start	end	number_of_dms	hypermethylated_samples	hypomethylated_samples
chr1	100	140	3	A	B
chr1	120	500	2	A	B
chr1	600	700	1	A	B
chr2	1	40	1	A	B
```

An allc table with three chromosomes holding 3, 1 and 2 records.

**Example file:** `methylpy_allc_three.tsv`

```text
chr1	1	+	CGA	1	2	0
chr1	2	+	CGA	1	2	0
chr1	3	+	CGA	1	2	0
chr2	1	+	CGA	1	2	0
chr3	1	+	CGA	1	2	0
chr3	2	+	CGA	1	2	0
```

An allc table with the unmethylated control chromosome chrL (2 methylated of 200 reads).

**Example file:** `methylpy_allc_S.tsv`

```text
chr1	10	+	CGA	8	10	1
chr1	20	+	CGT	1	10	1
chr1	30	+	CHH	0	9	0
chr1	40	+	CHH	5	9	1
chr1	50	+	CHH	1	50	0
chr1	60	+	CGC	0	0	0
chrL	10	+	CHH	1	100	0
chrL	20	+	CHH	1	100	0
```

The first binomial-test result file (column 7 is the P-value; 2.0 marks an untested site).

**Example file:** `methylpy_binom_1.tsv`

```text
chr1	20	+	CGN	1	10	0.2
chr1	10	+	CGN	8	10	0.0001
chr2	5	+	CHH	0	9	2.0
```

The second binomial-test result file.

**Example file:** `methylpy_binom_2.tsv`

```text
chr1	40	+	CHH	5	9	0.001
chr1	30	+	CHH	0	9	0.6
```

Seven lines of mpileup text.

**Example file:** `methylpy_reads.mpileup`

```text
chr1	1	C	5	.....
chr1	2	G	5	.....
chr1	3	A	5	.....
chr1	4	C	5	.....
chr1	5	C	5	.....
chr1	6	T	5	.....
chr1	7	C	5	.....
```

### 1) `dnakit.epigenome.estimate_dmr_fdr_cutoff()`

- **Purpose:** Find the P-value cutoff of an RMS permutation-test results table at a target FDR, to decide which sites count as differentially methylated sites (DMS).
- **Method:** Calls methylpy's `histogram_correction_DMRfind` (histogram method). First checks the gzip format, the header, the field count of every row, the P-value range, and that the last two columns `num_simulations_sig` and `num_simulations_run` are possible for `num_sims` and `num_sig_tests`. Then it re-estimates the number of null tests m0 until it changes by no more than `convergence_diff`, and takes the cutoff whose estimated FDR is closest to, and not above, `target_fdr`. When the iteration does not converge or no site is significant, methylpy falls back to `1/num_sims` and `outcome` says so. Hand check for the example: 10 tests, m0 about 4.13, cutoff 5/37.
- **API:** `dnakit.epigenome.estimate_dmr_fdr_cutoff(rms_path: 'str | Path', *, num_sims: 'int' = 3000, num_sig_tests: 'int' = 100, target_fdr: 'float' = 0.01, max_iterations: 'int' = 1000, convergence_diff: 'float' = 0.01, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** `rms_path`: a gzip-compressed RMS results table (header chr, pos, strand, mc_class, pvalue, then mc_, h_, frac_, mc_residual_ and uc_residual_ columns per sample, then num_simulations_sig and num_simulations_run); `num_sims` and `num_sig_tests`: the permutation settings the table was made with; `target_fdr`, `max_iterations`, `convergence_diff`: the iteration settings.
- **Output:** `p_value_cutoff`, `total_tests`, `outcome` (`converged`, `not_converged` or `no_significant_result`), `estimated_fdr`, `null_tests_estimate` (m0).
- **Example code:**

```python
import gzip
import shutil

from dnakit.epigenome import estimate_dmr_fdr_cutoff

with open("methylpy_rms.tsv", "rb") as source, gzip.open("methylpy_rms.tsv.gz", "wb") as target:
    shutil.copyfileobj(source, target)
result = estimate_dmr_fdr_cutoff(
    "methylpy_rms.tsv.gz",
    num_sims=100,
    num_sig_tests=5,
    target_fdr=0.1,
    convergence_diff=1,
    verbose=False,
)
print(result["total_tests"], result["outcome"])
print(round(result["p_value_cutoff"], 6), round(result["estimated_fdr"], 4))
```

- **Example output:**

```text
10 converged
0.135135 0.0964
```

### 2) `dnakit.epigenome.dmr_residual_cutoff()`

- **Purpose:** Find the quantile of the residuals of non-significant sites in an RMS results table, used as the threshold for calling a sample high or low at a differential site.
- **Method:** Calls methylpy's `get_resid_cutoff`: it collects all mc and uc residuals of the sites whose P-value is above `pvalue_cutoff` and takes the `100 - resid_cutoff * 100` percentile. Without such sites methylpy returns 2 and `default_used` is true. methylpy also counts the -1 placeholder of samples without data; the function returns the value as methylpy computes it. A P-value starting with `<` makes methylpy fail with a low-level error, so here it raises `ConfigurationError`. Hand check for the example: of the 4 non-significant sites three have residuals of ±0.4495, ±0.4495 and ±0.4688 and one has 0, so the 80th percentile is exactly 0.4495.
- **API:** `dnakit.epigenome.dmr_residual_cutoff(rms_path: 'str | Path', pvalue_cutoff: 'float', *, resid_cutoff: 'float' = 0.01, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** `rms_path`: a gzip-compressed RMS results table; `pvalue_cutoff`: 0 to 1, usually from `estimate_dmr_fdr_cutoff`; `resid_cutoff`: a fraction from 0 to 1, default 0.01.
- **Output:** `residual_cutoff`, `null_sites` (sites in the distribution), `default_used`.
- **Example code:**

```python
from dnakit.epigenome import dmr_residual_cutoff

result = dmr_residual_cutoff("methylpy_rms.tsv.gz", 0.135135, resid_cutoff=0.2, verbose=False)
print(result["null_sites"], result["default_used"])
print(round(result["residual_cutoff"], 4))
```

- **Example output:**

```text
4 False
0.4495
```

### 3) `dnakit.epigenome.collapse_dmr_windows()`

- **Purpose:** Merge significant sites into differentially methylated regions (DMR) by distance and residual direction, and write a DMS and a DMR BED file.
- **Method:** Calls methylpy's `collapse_dmr_windows` (the P-value column is fixed to column 5). A site with a P-value at most `sig_cutoff` joins the current region when it is on the same chromosome, within `max_dist` of the region end and its residual direction agrees; otherwise a new region starts. The direction is only checked when `resid_cutoff` is above 0, and samples are labelled hyper- or hypo-methylated; `collapse_samples` limits the samples considered, `sample_category` and `min_cluster` require a number of same-direction samples per category, and `min_num_dms` drops small regions. The table goes to `<output_path>`, the BED files to `<output_path>.DMS.bed` and `<output_path>.DMR.bed` (chromosomes get a `chr` prefix, starts are 0-based). When no site is significant only the header is written and there are no BED files. Existing files are never overwritten.
- **API:** `dnakit.epigenome.collapse_dmr_windows(rms_path: 'str | Path', output_path: 'str | Path', *, max_dist: 'int' = 100, resid_cutoff: 'float' = 0.0, sig_cutoff: 'float' = 0.01, min_num_dms: 'int' = 0, collapse_samples: 'Sequence[str] | None' = None, sample_category: 'Sequence[int] | None' = None, min_cluster: 'int' = 0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** `rms_path`: a gzip-compressed RMS results table sorted by chromosome and position; `output_path`: a new path; `max_dist`, `resid_cutoff`, `sig_cutoff`, `min_num_dms`, `collapse_samples`, `sample_category`, `min_cluster`: the merging settings.
- **Output:** `dmr_count`, `dms_count`, `regions` (at most 1000, each with `chrom`, `start`, `end`, `dms`, `hypermethylated_samples`, `hypomethylated_samples`), `regions_truncated`, `files` (paths of the table and the two BED files; None for a BED file that was not written).
- **Example code:**

```python
import tempfile

from dnakit.epigenome import collapse_dmr_windows

with tempfile.TemporaryDirectory() as work:
    result = collapse_dmr_windows(
        "methylpy_rms.tsv.gz",
        f"{work}/dmrs.tsv",
        max_dist=50,
        sig_cutoff=0.135135,
        resid_cutoff=2.0,
        verbose=False,
    )
print(result["dmr_count"], result["dms_count"])
for region in result["regions"]:
    print(
        region["chrom"],
        region["start"],
        region["end"],
        region["dms"],
        region["hypermethylated_samples"],
        region["hypomethylated_samples"],
    )
```

- **Example output:**

```text
3 6
1 100 140 3 ['A'] ['B']
1 600 620 2 ['B'] ['A']
2 80 80 1 ['A'] ['B']
```

### 4) `dnakit.epigenome.merge_dms_to_dmrs()`

- **Purpose:** Run the FDR cutoff, the residual cutoff and the region merging in one step, turning an RMS results table into a DMR table.
- **Method:** Calls methylpy's `merge_DMS_to_DMR`, which runs the histogram FDR correction, the residual cutoff and `collapse_dmr_windows` in turn. methylpy does not return the two cutoffs, so the function recomputes them with the same options after the merge and reports them. `num_sims` and `num_sig_tests` must match the settings the table was made with; a row whose simulation counts cannot occur raises an error. The output files are the same as for `collapse_dmr_windows`.
- **API:** `dnakit.epigenome.merge_dms_to_dmrs(rms_path: 'str | Path', output_path: 'str | Path', *, num_sims: 'int' = 3000, num_sig_tests: 'int' = 100, sig_cutoff: 'float' = 0.01, dmr_max_dist: 'int' = 250, min_num_dms: 'int' = 0, resid_cutoff: 'float' = 0.01, collapse_samples: 'Sequence[str] | None' = None, sample_category: 'Sequence[int] | None' = None, min_cluster: 'int' = 0, max_iterations: 'int' = 1000, convergence_diff: 'float' = 1, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** `rms_path`, `output_path`; `num_sims`, `num_sig_tests`; `sig_cutoff`: the target FDR; `dmr_max_dist`; `min_num_dms`; `resid_cutoff`: the residual percentile as a fraction; `collapse_samples`, `sample_category`, `min_cluster`; `max_iterations`, `convergence_diff` (default 1, as in methylpy).
- **Output:** `p_value_cutoff`, `fdr_outcome`, `residual_cutoff`, and the same `dmr_count`, `dms_count`, `regions` and `files` as `collapse_dmr_windows`.
- **Example code:**

```python
import tempfile

from dnakit.epigenome import merge_dms_to_dmrs

with tempfile.TemporaryDirectory() as work:
    result = merge_dms_to_dmrs(
        "methylpy_rms.tsv.gz",
        f"{work}/merged.tsv",
        num_sims=100,
        num_sig_tests=5,
        sig_cutoff=0.1,
        dmr_max_dist=50,
        resid_cutoff=0.2,
        convergence_diff=1,
        verbose=False,
    )
print(round(result["p_value_cutoff"], 6), round(result["residual_cutoff"], 4))
print(result["dmr_count"], result["dms_count"])
print([(r["chrom"], r["start"], r["end"], r["dms"]) for r in result["regions"]])
```

- **Example output:**

```text
0.135135 0.4495
3 6
[('1', 100, 140, 3), ('1', 600, 620, 2), ('2', 80, 80, 1)]
```

### 5) `dnakit.epigenome.add_methylation_levels_to_dmrs()`

- **Purpose:** Add the methylation level of every sample (methylated reads divided by total reads) to each region of a DMR table.
- **Method:** Calls methylpy's `get_methylation_levels_DMRfind` with one process. There is one allc file per sample (all chromosomes, each chromosome contiguous); the `.idx` index is written next to the allc file. Only sites whose context matches `mc_type` and whose coverage is within [`min_cov`, `max_cov`] count; a region without sites gives NA (None). Regions must be sorted by start within a chromosome; a region on a chromosome missing from an allc table raises an error (methylpy would repeat the previous region's value, or fail with a NameError on the first row). The output name must end with `.tsv` and the path may only contain letters, digits and `_ . / + = : , @ % -`, because methylpy removes its temporary files with a shell `rm`; temporary files are cleaned up after an error.
- **API:** `dnakit.epigenome.add_methylation_levels_to_dmrs(dmr_path: 'str | Path', allc_paths: 'Sequence[str | Path]', samples: 'Sequence[str]', output_path: 'str | Path', *, mc_type: 'str | Sequence[str]' = 'CGN', min_cov: 'int' = 0, max_cov: 'int | None' = None, buffer_line_number: 'int' = 100000, input_no_header: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** `dmr_path`: a table whose first three columns are chromosome, start and end (a header line by default); `allc_paths` and `samples`: one allc file per sample, with names of letters, digits, `_`, `.` and `-`; `output_path`; `mc_type`: a context string or list, default CGN; `min_cov`, `max_cov`, `buffer_line_number`; `input_no_header`: true when the DMR table has no header.
- **Output:** `output`, `samples`, `mc_type`, `dmr_count`, `rows` (at most 1000, each with `chrom`, `start`, `end`, `methylation_level`), `rows_truncated`. Hand check for the example: in region chr1:100-140 sample A has 9+10+8 methylated of 30 CGN reads = 0.9 (the CHH site at 125 is excluded).
- **Example code:**

```python
import tempfile

from dnakit.epigenome import add_methylation_levels_to_dmrs

with tempfile.TemporaryDirectory() as work:
    result = add_methylation_levels_to_dmrs(
        "methylpy_dmrs.tsv",
        ["methylpy_allc_A.tsv", "methylpy_allc_B.tsv"],
        ["A", "B"],
        f"{work}/levels.tsv",
        verbose=False,
    )
for row in result["rows"]:
    levels = {k: None if v is None else round(v, 4) for k, v in row["methylation_level"].items()}
    print(row["chrom"], row["start"], row["end"], levels)
```

- **Example output:**

```text
chr1 100 140 {'A': 0.9, 'B': 0.0667}
chr1 120 500 {'A': 0.7667, 'B': 0.2}
chr1 600 700 {'A': None, 'B': None}
chr2 1 40 {'A': 0.0, 'B': 0.5}
```

### 6) `dnakit.epigenome.add_cytosine_info_to_dmrs()`

- **Purpose:** Add the methylation level, the total base calls and the number of covered sites of every sample to each region of a DMR table.
- **Method:** Calls methylpy's `get_c_info_DMRfind`; the rules are those of `add_methylation_levels_to_dmrs`. Three files are written: `<output_path>` (methylation level), `<output_path>.total_basecalls.tsv` and `<output_path>.covered_sites.tsv`. A region without matching sites gives NA (None) in all three.
- **API:** `dnakit.epigenome.add_cytosine_info_to_dmrs(dmr_path: 'str | Path', allc_paths: 'Sequence[str | Path]', samples: 'Sequence[str]', output_path: 'str | Path', *, mc_type: 'str | Sequence[str]' = 'CGN', min_cov: 'int' = 0, max_cov: 'int | None' = None, buffer_line_number: 'int' = 100000, input_no_header: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** The same as `add_methylation_levels_to_dmrs`.
- **Output:** `output`, `total_basecalls_output`, `covered_sites_output`, `samples`, `mc_type`, `dmr_count`, `rows` (each with `chrom`, `start`, `end` and per sample `methylation_level`, `total_basecalls`, `covered_sites`), `rows_truncated`.
- **Example code:**

```python
import tempfile

from dnakit.epigenome import add_cytosine_info_to_dmrs

with tempfile.TemporaryDirectory() as work:
    result = add_cytosine_info_to_dmrs(
        "methylpy_dmrs.tsv",
        ["methylpy_allc_A.tsv", "methylpy_allc_B.tsv"],
        ["A", "B"],
        f"{work}/info.tsv",
        mc_type="CGN",
        verbose=False,
    )
for row in result["rows"]:
    a, b = row["samples"]["A"], row["samples"]["B"]
    print(row["chrom"], row["start"], a["total_basecalls"], a["covered_sites"], b["total_basecalls"])
```

- **Example output:**

```text
chr1 100 30 3 30
chr1 120 30 3 30
chr1 600 None None None
chr2 1 9 1 6
```

### 7) `dnakit.epigenome.split_allc_by_chromosome()`

- **Purpose:** Split an allc table into `allc_<sample>_<chrom>.tsv[.gz]`, one file per chromosome.
- **Method:** Calls methylpy's `do_split_allc_file`. Each chromosome must appear in one block: methylpy overwrites the file of a chromosome that reappears, so this is checked first. Empty files, lines with fewer than 6 columns and chromosome names that cannot be file names are also rejected before the call. methylpy lists the names without `.gz` even when it writes `.tsv.gz`; the paths returned here are the files that exist. The output directory must exist and the target files must not.
- **API:** `dnakit.epigenome.split_allc_by_chromosome(allc_path: 'str | Path', sample: 'str', output_dir: 'str | Path', *, compress_output: 'bool' = True, buffer_line_number: 'int' = 100000, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** `allc_path`; `sample`: 1 to 64 letters, digits, `_`, `.` or `-`; `output_dir`; `compress_output`: true by default, writing `.tsv.gz`; `buffer_line_number`.
- **Output:** `files` (sorted by chromosome, each with `chrom`, `path`, `records`), `records`, `compressed`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import split_allc_by_chromosome

with tempfile.TemporaryDirectory() as work:
    result = split_allc_by_chromosome("methylpy_allc_A.tsv", "A", work, verbose=False)
    for item in result["files"]:
        print(item["chrom"], Path(item["path"]).name, item["records"])
```

- **Example output:**

```text
chr1 allc_A_chr1.tsv.gz 5
chr2 allc_A_chr2.tsv.gz 1
```

### 8) `dnakit.epigenome.split_allc_into_chunks()`

- **Purpose:** Distribute the chromosomes of an allc table over at most `num_chunks` files named `allc_<sample>_<index>.tsv[.gz]`.
- **Method:** Calls methylpy's `do_split_allc_file_chunk`: a chunk holds `ceil(records / num_chunks)` records, chromosomes are packed from the smallest to the largest, a new chunk starts when the current one is full, and a chromosome is never split. With few or very uneven chromosomes fewer files than `num_chunks` are written. Lines keep their input order inside a chunk. An empty file is rejected before the call (methylpy returns 0 and writes nothing).
- **API:** `dnakit.epigenome.split_allc_into_chunks(allc_path: 'str | Path', sample: 'str', num_chunks: 'int', output_dir: 'str | Path', *, compress_output: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** `allc_path`, `sample`, `num_chunks` (1 to 10000), `output_dir`; `compress_output`.
- **Output:** `chunks` (each with `index`, `path`, `records`, `chromosomes`, read back from the files written), `requested_chunks`, `records`, `compressed`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import split_allc_into_chunks

with tempfile.TemporaryDirectory() as work:
    result = split_allc_into_chunks("methylpy_allc_three.tsv", "S", 2, work, verbose=False)
    for chunk in result["chunks"]:
        print(Path(chunk["path"]).name, chunk["records"], chunk["chromosomes"])
```

- **Example output:**

```text
allc_S_0.tsv.gz 3 ['chr2', 'chr3']
allc_S_1.tsv.gz 3 ['chr1']
```

### 9) `dnakit.epigenome.run_allc_binomial_test()`

- **Purpose:** Call methylated sites in an allc table with a binomial test and a Benjamini-Hochberg cutoff per context.
- **Method:** Calls methylpy's `perform_binomial_test` with one process. The non-conversion rate comes from `unmethylated_control`: a number from 0 to 1, or an unmethylated control region such as `chrL`, `chrL:100` or `chrL:100-2000`, for which it is the methylated reads divided by the total reads in the region. The region's chromosome name is looked up verbatim (`perform_binomial_test` always strips a `chr` prefix, so the rate is computed here first). Sites with at least `min_cov` reads and at least one methylated read are tested; other sites with reads get 0, and sites without reads are dropped. Column 7 of the output allc is 1 for significant sites. The allc index is written next to the input; `sort` runs with `LC_ALL=C` and the caller's setting is restored afterwards. The output directory path may only contain letters, digits and `_ . / + = : , @ % -`.
- **API:** `dnakit.epigenome.run_allc_binomial_test(allc_path: 'str | Path', sample: 'str', output_dir: 'str | Path', unmethylated_control: 'float | str', *, min_cov: 'int' = 2, sig_cutoff: 'float' = 0.01, sort_mem: 'str' = '500M', compress_output: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** `allc_path`, `sample`, `output_dir`, `unmethylated_control`; `min_cov`: default 2; `sig_cutoff`: the FDR, default 0.01; `sort_mem`: such as `500M`; `compress_output`.
- **Output:** `output` (`allc_<sample>.tsv[.gz]`), `index`, `non_conversion_rate`, `sites`, `significant_sites`, `contexts` (`sites` and `significant_sites` per context, sorted). Hand check for the example: chrL has 2 of 200 reads, rate 0.01; the sites with 8/10 and 5/9 are significant and those with 1/10, 1/50 and 1/100 are not.
- **Example code:**

```python
import gzip
import tempfile

from dnakit.epigenome import run_allc_binomial_test

with tempfile.TemporaryDirectory() as work:
    result = run_allc_binomial_test("methylpy_allc_S.tsv", "S", work, "chrL", verbose=False)
    with gzip.open(result["output"], "rt") as handle:
        flags = [line.split("\t")[:2] + line.split("\t")[6:7] for line in handle]
print(result["non_conversion_rate"], result["sites"], result["significant_sites"])
print([item["context"] + ":" + str(item["significant_sites"]) for item in result["contexts"]])
print([f"{row[0]}:{row[1]}={row[2].strip()}" for row in flags])
```

- **Example output:**

```text
0.01 7 2
['CGA:1', 'CGT:0', 'CHH:1']
['chr1:10=1', 'chr1:20=0', 'chr1:30=0', 'chr1:40=1', 'chr1:50=0', 'chrL:10=0', 'chrL:20=0']
```

### 10) `dnakit.epigenome.filter_allc_by_pvalue_cutoffs()`

- **Purpose:** Rewrite binomial-test result files as allc records with a significance flag, using a P-value cutoff per context.
- **Method:** Calls methylpy's `filter_files_by_pvalue_combined`. methylpy sorts its input files in place with `sort -o`, so the inputs are copied to a temporary directory first and the originals stay unchanged. Each file is sorted by chromosome and position; a site whose P-value (column 7) is at most the cutoff of its context gets 1 in column 7, any other site 0, and a P-value of 2.0 marks an untested site that is always 0. Several files are written one after another in the order given. methylpy appends `.gz` to an output name that lacks it; the real path is reported. `sort` runs with `LC_ALL=C`.
- **API:** `dnakit.epigenome.filter_allc_by_pvalue_cutoffs(binomial_paths: 'Sequence[str | Path]', output_path: 'str | Path', cutoffs: 'Mapping[str, float]', *, compress_output: 'bool' = True, sort_mem: 'str' = '500M', verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** `binomial_paths`: 1 to 32 result files (at least 7 columns per line); `output_path`; `cutoffs`: a mapping from context to cutoff (0 to 1) that covers every tested context in the files; `compress_output`, `sort_mem`.
- **Output:** `output`, `records`, `significant_sites`, `contexts` (`sites` and `significant_sites` per context, sorted), `inputs`.
- **Example code:**

```python
import gzip
import tempfile

from dnakit.epigenome import filter_allc_by_pvalue_cutoffs

with tempfile.TemporaryDirectory() as work:
    result = filter_allc_by_pvalue_cutoffs(
        ["methylpy_binom_1.tsv", "methylpy_binom_2.tsv"],
        f"{work}/filtered.tsv",
        {"CGN": 0.05, "CHH": 0.01},
        verbose=False,
    )
    with gzip.open(result["output"], "rt") as handle:
        lines = handle.read().splitlines()
print(result["records"], result["significant_sites"])
for line in lines:
    print(line)
```

- **Example output:**

```text
5 2
chr1	10	+	CGN	8	10	1
chr1	20	+	CGN	1	10	0
chr2	5	+	CHH	0	9	0
chr1	30	+	CHH	0	9	0
chr1	40	+	CHH	5	9	1
```

### 11) `dnakit.epigenome.merge_allc_tables()`

- **Purpose:** Add up the methylated and total read counts of several allc tables site by site.
- **Method:** Calls methylpy's `merge_allc_files` with one process (`mini_batch` sets how many files are merged at a time, which is methylpy's minibatch step). When `compress_output` is true and the name lacks `.gz`, methylpy appends it; the real path is reported and that path and its `.idx` must not exist. methylpy writes chromosomes in an order that changes with the process hash seed, so the merged table is rewritten with chromosomes sorted by name and re-indexed, which makes the result repeatable. Inputs are checked first: each chromosome contiguous, counts integer; with `skip_snp_info` false every line needs at least 9 columns and the match and mismatch lists are added position by position. The `.idx` indexes of the inputs are written next to them.
- **API:** `dnakit.epigenome.merge_allc_tables(allc_paths: 'Sequence[str | Path]', output_path: 'str | Path', *, mini_batch: 'int' = 100, compress_output: 'bool' = True, skip_snp_info: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** `allc_paths`: 2 to 1000 allc tables; `output_path`; `mini_batch`: 2 to 1000; `compress_output`; `skip_snp_info`.
- **Output:** `output`, `index`, `records`, `chromosomes` (`chrom` and `records`, sorted by chromosome name), `compressed`, `inputs`.
- **Example code:**

```python
import gzip
import tempfile
from pathlib import Path

from dnakit.epigenome import merge_allc_tables

with tempfile.TemporaryDirectory() as work:
    result = merge_allc_tables(
        ["methylpy_allc_A.tsv", "methylpy_allc_B.tsv"], f"{work}/merged.tsv", verbose=False
    )
    print(Path(result["output"]).name, result["records"])
    with gzip.open(result["output"], "rt") as handle:
        for line in handle:
            print(line.rstrip("\n"))
```

- **Example output:**

```text
merged.tsv.gz 6
chr1	100	+	CGA	10	20	1
chr1	120	+	CGT	10	20	1
chr1	125	+	CHH	1	4	1
chr1	140	+	CGC	9	20	1
chr1	500	+	CGG	10	20	1
chr2	30	+	CGA	3	15	1
```

### 12) `dnakit.epigenome.split_mpileup_chunks()`

- **Purpose:** Cut a plain-text mpileup file into chunks whose neighbours share two lines, so that C contexts at chunk borders stay complete.
- **Method:** Calls methylpy's `split_mpileup_file`: a chunk has `ceil(lines / num_chunks)` lines and every chunk after the first starts with the last two lines of the previous one. methylpy counts lines with `wc -l`, which counts newline bytes of the compressed data for a gzip file, so gzip input is rejected; on an empty file methylpy exits and the exit becomes `ConfigurationError`. Chunk files are `<output_prefix>0`, `<output_prefix>1` and so on, and none may exist. Fewer chunks than `num_chunks` can result.
- **API:** `dnakit.epigenome.split_mpileup_chunks(mpileup_path: 'str | Path', output_prefix: 'str | Path', num_chunks: 'int', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** `mpileup_path`: a plain-text file; `output_prefix`: the path prefix of the chunk files (ends in a file-name prefix, not a directory); `num_chunks`: 1 to 10000.
- **Output:** `records`, `chunk_size`, `chunks` (each with `index`, `path`, `records`, `overlap_records`, read back from the files written).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import split_mpileup_chunks

with tempfile.TemporaryDirectory() as work:
    result = split_mpileup_chunks("methylpy_reads.mpileup", f"{work}/chunk_", 3, verbose=False)
    for chunk in result["chunks"]:
        print(Path(chunk["path"]).name, chunk["records"], chunk["overlap_records"])
```

- **Example output:**

```text
chunk_0 3 0
chunk_1 5 2
chunk_2 3 2
```

### Functions not wrapped

The following 5 functions get no new wrapper; the status and the reason of each follow.

- `methylpy.DMRfind.DMRfind` (`needs-binary`): The end-to-end DMR pipeline needs methylpy's bundled permutation program `run_rms_tests.out` (it links `libgsl.so.0`, which is not installed here) and deletes its temporary allc files; its post-processing steps are wrapped separately.
- `methylpy.DMRfind.run_rms_tests` (`needs-binary`): Launches methylpy's bundled compiled program `run_rms_tests.out` (it links `libgsl.so.0`, which is not installed here).
- `methylpy.call_mc_se.fasta_iter` (`deprecated`): Python 2 code: it calls `.next()` on a `groupby` iterator and raises `AttributeError` on Python 3 (verified on a tiny FASTA); use `read_fasta_chromosome` to read a chromosome.
- `methylpy.utilities.infer_samples` (`internal`): A helper that cuts the sample name out of an `allc_<sample>.tsv` file name; the same object as `DMRfind.infer_samples`.
- `methylpy.utilities.merge_allc_files_minibatch` (`family-variant` → `merge_allc_tables`): The minibatch step inside `merge_allc_files`, reached through the `mini_batch` option of `merge_allc_tables`.

## Additional MACS3 functions (the HMMRATAC pipeline)

This section adds 10 MACS3 functions in the new module `MACS3_review_functions.py` of `dnakit.epigenome`. The examples import from the module file; once exported from `dnakit.epigenome`, the same functions are callable as `dnakit.epigenome.<name>()`. They follow the steps of the HMMRATAC pipeline for open chromatin in ATAC-seq: split fragments into the short, mono-, di- and tri-nucleosomal signal tracks, cut them into bins, train, save, load and use the three-state hidden Markov model, then turn the predicted probabilities into bedGraph tracks, a states path, BED and narrowPeak files. They do not run the `macs3 hmmratac` command line and need no external program.

The MACS3 functions pass fragment tracks and hmmlearn model objects to one another, which a wrapper cannot keep between calls, so plain data is the hand-off: fragments are `(chromosome, left, right)` triples, training data is rows of four numbers per bin, the model is a MACS3 JSON file, and the predicted probabilities are the kind of text file that the MACS3 pipeline keeps in a temporary file (`chromosome,bin end,probability,probability,probability`). On bad input the upstream functions trip assertions, call `sys.exit` or raise low-level errors; the wrappers check the inputs first, never overwrite an existing file, and turn those errors into `ConfigurationError`. The return value is a dict with `status`, `operation` and `provenance`; coordinates are 0-based and right-open.

The example data is made up by hand to show the formats and the calculation; it does not stand for a real experiment. The output files of the examples go to a temporary directory, so the examples can be run repeatedly.

### Example data

Training data for 45 bins, one row per bin with the values of the short, mono, di and tri signal tracks. The rows form 3 groups of 15 bins; each group runs through background, nucleosomal, open and background levels. The values are made up by hand.

**Example file:** `macs3_training_data.tsv`

```text
0.30	0.69	0.27	0.11
0.29	0.58	0.37	0.10
0.30	0.65	0.34	0.10
0.32	0.53	0.29	0.09
0.84	4.91	1.61	0.49
0.98	5.77	2.02	0.42
0.99	6.17	2.18	0.45
7.62	1.52	0.75	0.22
6.64	2.26	0.59	0.33
8.31	1.93	0.84	0.32
9.00	1.94	0.74	0.28
0.88	5.97	1.81	0.56
0.78	5.21	1.77	0.37
0.37	0.43	0.29	0.09
0.36	0.46	0.34	0.09
0.29	0.55	0.32	0.09
0.30	0.63	0.37	0.07
0.35	0.67	0.28	0.10
0.28	0.72	0.31	0.10
0.97	5.86	1.96	0.45
1.25	4.62	1.13	0.49
0.98	6.27	1.95	0.49
8.32	2.23	0.76	0.29
9.86	2.13	0.71	0.38
8.74	1.86	0.69	0.31
7.20	1.75	0.68	0.28
1.13	5.69	1.65	0.54
1.01	6.61	2.29	0.49
0.29	0.60	0.26	0.11
0.35	0.61	0.29	0.10
0.27	0.54	0.29	0.09
0.28	0.49	0.31	0.10
0.26	0.43	0.30	0.11
0.27	0.57	0.28	0.11
0.89	6.71	1.93	0.56
1.00	5.83	1.65	0.46
0.97	6.47	2.06	0.46
8.39	2.24	0.79	0.28
7.62	2.19	0.85	0.27
8.36	1.88	0.73	0.34
8.79	1.83	0.81	0.32
0.92	5.91	2.16	0.39
1.04	6.53	2.12	0.42
0.31	0.54	0.32	0.11
0.31	0.54	0.28	0.11
```

A hand-written MACS3 model file with diagonal covariances: state 0 is background, state 1 is nucleosomal and state 2 is open. Real model files come from `save_hmmratac_model` or `macs3 hmmratac`.

**Example file:** `macs3_model.json`

```text
{
 "startprob": [0.2, 0.3, 0.5],
 "transmat": [[0.8, 0.1, 0.1], [0.1, 0.8, 0.1], [0.1, 0.1, 0.8]],
 "means": [[1.0, 1.0, 1.0, 1.0], [5.0, 5.0, 5.0, 5.0], [10.0, 10.0, 10.0, 10.0]],
 "covars": [[1.0, 1.0, 1.0, 1.0], [4.0, 4.0, 4.0, 4.0], [9.0, 9.0, 9.0, 9.0]],
 "covariance_type": "diag",
 "n_features": 4,
 "i_open_region": 2,
 "i_background_region": 0,
 "i_nucleosomal_region": 1,
 "hmm_binsize": 10,
 "hmm_type": "gaussian"
}
```

A predicted-probability file: each line is `chromosome,bin end,state 0,state 1,state 2` with a bin size of 10. Here state 2 is open, state 0 is nucleosomal and state 1 is background. The values are made up by hand.

**Example file:** `macs3_proba.csv`

```text
chr1,20,0.1,0.8,0.1
chr1,30,0.8,0.1,0.1
chr1,40,0.1,0.1,0.8
chr1,50,0.1,0.1,0.8
chr1,60,0.8,0.1,0.1
chr1,70,0.1,0.8,0.1
chr1,100,0.8,0.1,0.1
chr1,110,0.1,0.1,0.8
chr1,120,0.8,0.1,0.1
chr2,40,0.1,0.8,0.1
chr2,50,0.8,0.1,0.1
chr2,60,0.1,0.1,0.8
chr2,70,0.8,0.1,0.1
```

The bedGraph track of the open-state probability of the file above (the output of `hmmratac_proba_to_bedgraph`), used as the score track of the open regions.

**Example file:** `macs3_open.bdg`

```text
chr1	0	10	0.00000
chr1	10	30	0.10000
chr1	30	50	0.80000
chr1	50	70	0.10000
chr1	70	90	0.00000
chr1	90	100	0.10000
chr1	100	110	0.80000
chr1	110	120	0.10000
chr2	0	30	0.00000
chr2	30	50	0.10000
chr2	50	60	0.80000
chr2	60	70	0.10000
```

### 1) `dnakit.epigenome.digest_hmmratac_fragments()`

- **Purpose:** Split paired-end fragments (for example from ATAC-seq) by length into the short, mono-, di- and tri-nucleosomal signal tracks that HMMRATAC takes as input.
- **Method:** Calls `hmmratac_fragment_weights` for the weight of each fragment length under the four Gaussians, builds a MACS3 `PETrackI` from the fragments, and calls `HMMR_Signal_Processing.generate_digested_signals`: every fragment adds its weight for each class over the bases it covers. A fragment whose density under all four Gaussians is below `min_fragment_probability` gets weight 0. The four bedGraph files can be written in the MACS3 format.
- **API:** `dnakit.epigenome.digest_hmmratac_fragments(fragments: 'Sequence[Sequence[Any]]', means: 'Sequence[float]', stddevs: 'Sequence[float]', *, min_fragment_probability: 'float' = 0.001, output_prefix: 'str | Path | None' = None, include_tracks: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** `fragments` is a list of `(chromosome, left, right)` triples (0-based, right-open, at most 10000 bases long, chromosome names without spaces or commas); `means` and `stddevs` have four values each, in the order short, mono, di, tri; `min_fragment_probability`; with `output_prefix`, the files `<prefix>_short.bdg`, `_mono.bdg`, `_di.bdg` and `_tri.bdg` are written (an existing file is an error, nothing is overwritten); with `include_tracks=False` the tracks are left out of the return value.
- **Output:** `tracks` holds, for each class and chromosome, `positions` (the end of each segment; the first segment starts at 0) and `values` (the summed weights on that segment); `intervals` is the number of segments per class; also `fragments`, `chromosomes` and `classes`; `files` when `output_prefix` is given. In the example the two 50 bp fragments overlap on 120-150, where the short track is 2; the 200 bp fragment belongs to the mono track.
- **Example code:**

```python
from dnakit.epigenome import digest_hmmratac_fragments

fragments = [("chr1", 100, 150), ("chr1", 120, 170), ("chr1", 300, 500)]
result = digest_hmmratac_fragments(
    fragments, [50, 200, 400, 600], [20, 20, 20, 20], verbose=False
)
print(result["fragments"], result["chromosomes"], result["intervals"])
for label in result["classes"]:
    track = result["tracks"][label]["chr1"]
    print(label, track["positions"], [round(value, 3) for value in track["values"]])
```

- **Example output:**

```text
3 ['chr1'] {'short': 6, 'mono': 6, 'di': 2, 'tri': 1}
short [100, 120, 150, 170, 300, 500] [0.0, 1.0, 2.0, 1.0, 0.0, 0.0]
mono [100, 120, 150, 170, 300, 500] [0.0, 0.0, 0.0, 0.0, 0.0, 1.0]
di [300, 500] [0.0, 0.0]
tri [500] [0.0]
```

### 2) `dnakit.epigenome.extract_hmmratac_training_data()`

- **Purpose:** Cut the four signal tracks into equal bins inside regions, giving the data rows for HMM training or decoding.
- **Method:** Builds the four tracks as the previous function does, expands the regions by `flank` on both sides and merges overlaps (as the HMMRATAC pipeline does), then calls `HMMR_Signal_Processing.extract_signals_from_regions`: each bin takes the value of the four tracks there, values below 0.0001 are raised to 0.0001, and `hmm_type="poisson"` truncates to whole numbers. A bin is named by its end and covers `[end - binsize, end)`; no bins are made past the end of the signal; chromosomes come in reverse alphabetical order (MACS3's behaviour). It is an error when no bin falls inside the regions.
- **API:** `dnakit.epigenome.extract_hmmratac_training_data(fragments: 'Sequence[Sequence[Any]]', regions: 'Sequence[Sequence[Any]]', means: 'Sequence[float]', stddevs: 'Sequence[float]', *, binsize: 'int' = 10, flank: 'int' = 0, hmm_type: 'str' = 'gaussian', min_fragment_probability: 'float' = 0.001, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** `fragments` (as above); `regions` is a list of `(chromosome, start, end)` triples; `means` and `stddevs`; `binsize` (default 10); `flank` (default 0); `hmm_type` is `"gaussian"` or `"poisson"`; `min_fragment_probability`.
- **Output:** `bins` (`[chromosome, end]`), `data` (four values per bin), `lengths` (the number of bins of each merged region, used as the HMM sequence lengths), `regions` (the number of merged regions) and `binsize`. In the example each of the two regions has 3 bins.
- **Example code:**

```python
from dnakit.epigenome import extract_hmmratac_training_data

fragments = [("chr1", 100, 150), ("chr1", 120, 170), ("chr1", 300, 500)]
result = extract_hmmratac_training_data(
    fragments,
    [("chr1", 100, 170), ("chr1", 300, 360)],
    [50, 200, 400, 600],
    [20, 20, 20, 20],
    binsize=20,
    verbose=False,
)
print(result["regions"], result["lengths"])
for (chrom, end), row in zip(result["bins"], result["data"]):
    print(chrom, end - 20, end, [round(value, 3) for value in row])
```

- **Example output:**

```text
2 [3, 3]
chr1 100 120 [1.0, 0.0, 0.0, 0.0]
chr1 120 140 [2.0, 0.0, 0.0, 0.0]
chr1 140 160 [1.0, 0.0, 0.0, 0.0]
chr1 300 320 [0.0, 1.0, 0.0, 0.0]
chr1 320 340 [0.0, 1.0, 0.0, 0.0]
chr1 340 360 [0.0, 1.0, 0.0, 0.0]
```

### 3) `dnakit.epigenome.train_hmmratac_model()`

- **Purpose:** Train the three-state hidden Markov model of HMMRATAC on binned four-track data and mark its open, nucleosomal and background states.
- **Method:** Calls `HMMR_HMM.hmm_training` (hmmlearn's `GaussianHMM` or `PoissonHMM`, Baum-Welch training). The random generator is created locally from `random_seed`, so global random state is untouched. States are marked as the HMMRATAC pipeline does: the state whose four means (lambdas for Poisson) have the largest sum is open, the smallest sum is background, and the remaining one is nucleosomal. EM can end in a local optimum and the result depends on the seed, so try another `random_seed` if the states look wrong; two runs with the same seed can differ by about 1e-10 (parallel floating-point summation order). A fitted covariance that is not positive definite is an error.
- **API:** `dnakit.epigenome.train_hmmratac_model(data: 'Sequence[Sequence[float]]', lengths: 'Sequence[int]', *, hmm_type: 'str' = 'gaussian', covariance_type: 'str' = 'full', random_seed: 'int' = 12345, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** `data` has four non-negative numbers per bin (whole numbers for Poisson); `lengths` lists the number of bins of each sequence and must add up to the number of rows; `hmm_type` (default `"gaussian"`); `covariance_type` is `"full"` or `"diag"`; `random_seed` (default 12345).
- **Output:** `open_state`, `nucleosomal_state`, `background_state` (state indexes), `start_probabilities`, `transition_matrix`, `means` and `covariances` for Gaussian models (`lambdas` for Poisson), `converged`, `iterations` and `log_likelihood`. The training data of the example has three clear levels, and the model finds the background, open (high short track) and nucleosomal (high mono track) states.
- **Example code:**

```python
from dnakit.epigenome import train_hmmratac_model

with open("macs3_training_data.tsv", encoding="utf-8") as handle:
    rows = [[float(value) for value in line.split("\t")] for line in handle]
result = train_hmmratac_model(rows, [15, 15, 15], random_seed=10151, verbose=False)
print(result["converged"], result["open_state"], result["nucleosomal_state"], result["background_state"])
for state, means in enumerate(result["means"]):
    print(state, [round(value, 1) for value in means])
print([round(value, 2) for value in result["transition_matrix"][result["open_state"]]])
```

- **Example output:**

```text
True 1 2 0
0 [0.3, 0.6, 0.3, 0.1]
1 [8.2, 2.0, 0.7, 0.3]
2 [1.0, 5.9, 1.9, 0.5]
[0.0, 0.75, 0.25]
```

### 4) `dnakit.epigenome.save_hmmratac_model()`

- **Purpose:** Train the model and write it, with its state indexes and bin size, as a MACS3 JSON model file for prediction.
- **Method:** Training and state marking are those of `train_hmmratac_model`; then `HMMR_HMM.hmm_model_save` writes the file and `hmm_model_init` reads it back as a check. A file that cannot be read back is deleted and reported as an error. An existing file is never overwritten, and the output directory must exist.
- **API:** `dnakit.epigenome.save_hmmratac_model(data: 'Sequence[Sequence[float]]', lengths: 'Sequence[int]', output_path: 'str | Path', *, binsize: 'int' = 10, hmm_type: 'str' = 'gaussian', covariance_type: 'str' = 'full', random_seed: 'int' = 12345, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** `data`, `lengths`, `hmm_type`, `covariance_type` and `random_seed` as in the training function; `output_path` is the JSON file to write; `binsize` records the bin size of the training data, default 10.
- **Output:** `output` (the file path), `binsize`, `converged` and the three state indexes. The file is the MACS3 model JSON: `startprob`, `transmat`, `means` and `covars` (`lambdas` for Poisson), the indexes `i_open_region` and so on, `hmm_binsize` and `hmm_type`.
- **Example code:**

```python
import json
import tempfile

from dnakit.epigenome import save_hmmratac_model

with open("macs3_training_data.tsv", encoding="utf-8") as handle:
    rows = [[float(value) for value in line.split("\t")] for line in handle]
with tempfile.TemporaryDirectory() as work:
    result = save_hmmratac_model(
        rows, [15, 15, 15], f"{work}/model.json", binsize=10, random_seed=10151, verbose=False
    )
    with open(result["output"], encoding="utf-8") as handle:
        model = json.load(handle)
print(result["converged"], result["binsize"])
print(result["open_state"], result["nucleosomal_state"], result["background_state"])
print(sorted(model))
print(model["hmm_type"], model["covariance_type"], model["hmm_binsize"], len(model["means"]))
```

- **Example output:**

```text
True 10
1 2 0
['covariance_type', 'covars', 'hmm_binsize', 'hmm_type', 'i_background_region', 'i_nucleosomal_region', 'i_open_region', 'means', 'n_features', 'startprob', 'transmat']
gaussian full 10 3
```

### 5) `dnakit.epigenome.load_hmmratac_model()`

- **Purpose:** Read a MACS3 HMMRATAC model file and return its start probabilities, transition matrix, emission parameters, state indexes and bin size.
- **Method:** Checks the JSON first (array shapes, `n_features` equal to 4, three distinct state indexes), then calls `HMMR_HMM.hmm_model_init` to rebuild the model from the file; the parameters are read from the rebuilt model. A file without `hmm_type` comes from an older MACS3 and is read as a Gaussian model. A covariance that is not positive definite is an error.
- **API:** `dnakit.epigenome.load_hmmratac_model(model_path: 'str | Path', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** `model_path` is the model JSON file.
- **Output:** `hmm_type`, `binsize`, `open_state`, `nucleosomal_state`, `background_state`, `state_labels` (`open`, `nuc` or `bg` for each state index), `start_probabilities`, `transition_matrix`, `means`, `covariances` and `covariance_type` for Gaussian models (`lambdas` for Poisson), and `features`. In the hand-written model of the example state 2 is open, state 1 is nucleosomal and state 0 is background.
- **Example code:**

```python
from dnakit.epigenome import load_hmmratac_model

result = load_hmmratac_model("macs3_model.json", verbose=False)
print(result["hmm_type"], result["binsize"], result["state_labels"])
print(result["open_state"], result["nucleosomal_state"], result["background_state"])
print(result["start_probabilities"])
print(result["means"][2], result["covariances"][1][0])
```

- **Example output:**

```text
gaussian 10 ['bg', 'nuc', 'open']
2 1 0
[0.2, 0.3, 0.5]
[10.0, 10.0, 10.0, 10.0] [4.0, 0.0, 0.0, 0.0]
```

### 6) `dnakit.epigenome.predict_hmmratac_states()`

- **Purpose:** Compute, with a model file, the posterior probability that each bin is in the open, nucleosomal or background state.
- **Method:** Reads the model with `hmm_model_init`, then calls `HMMR_HMM.hmm_predict` (`predict_proba`, the forward-backward algorithm). The columns of the return value are open, nucleosomal, background, and `labels` takes the most probable state (open wins a tie, then nucleosomal). With both `bins` and `output_path`, it writes the predicted-probability file the MACS3 pipeline uses: `chromosome,end,probability,probability,probability`, with the three probabilities in the model's own state order and six decimals, ready for the next three functions.
- **API:** `dnakit.epigenome.predict_hmmratac_states(model_path: 'str | Path', data: 'Sequence[Sequence[float]]', lengths: 'Sequence[int]', *, bins: 'Sequence[Sequence[Any]] | None' = None, output_path: 'str | Path | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** `model_path`; `data` and `lengths` as in the training function (whole numbers for a Poisson model); `bins` has one `[chromosome, end]` per row; `output_path` is the probability file to write and is never overwritten. `bins` and `output_path` are given together or not at all.
- **Output:** `probabilities`, `columns`, `labels`, `label_counts`, `rows`, `binsize` and the three state indexes; `output` when a file is written. In the example each row is a sequence of one bin, so the posterior is the start probability times the emission density, normalised: the first row is background and the second is open.
- **Example code:**

```python
import tempfile

from dnakit.epigenome import predict_hmmratac_states

data = [[2.0, 2.0, 2.0, 2.0], [10.0, 10.0, 10.0, 10.0]]
with tempfile.TemporaryDirectory() as work:
    result = predict_hmmratac_states(
        "macs3_model.json",
        data,
        [1, 1],
        bins=[["chr1", 10], ["chr1", 20]],
        output_path=f"{work}/proba.csv",
        verbose=False,
    )
    with open(result["output"], encoding="utf-8") as handle:
        print(handle.read().strip())
print(result["columns"])
for row, label in zip(result["probabilities"], result["labels"]):
    print([round(value, 4) for value in row], label)
print(result["label_counts"])
```

- **Example output:**

```text
chr1,10,0.992363,0.007637,0.000000
chr1,20,0.000000,0.000011,0.999989
['open', 'nucleosomal', 'background']
[0.0, 0.0076, 0.9924] bg
[1.0, 0.0, 0.0] open
{'open': 1, 'nuc': 0, 'bg': 1}
```

### 7) `dnakit.epigenome.hmmratac_proba_to_bedgraph()`

- **Purpose:** Turn a predicted-probability file into three bedGraph tracks: the probability that each bin is open, nucleosomal or background.
- **Method:** Checks the file first (five fields per line, non-negative probabilities, one block per chromosome, no overlapping bins), then calls `hmmratac_cmd.save_proba_to_bedGraph`. Gaps between bins and the part before the first bin of each chromosome count as background (open 0, nucleosomal 0, background 1); neighbouring intervals with equal values are merged. bedGraph intervals are 0-based and right-open.
- **API:** `dnakit.epigenome.hmmratac_proba_to_bedgraph(proba_path: 'str | Path', output_prefix: 'str | Path', *, open_state: 'int', nucleosomal_state: 'int', background_state: 'int', binsize: 'int' = 10, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** `proba_path` is the predicted-probability file; `output_prefix` is the output name prefix, which writes `<prefix>_open.bdg`, `<prefix>_nuc.bdg` and `<prefix>_bg.bdg` (an existing file is an error); `open_state`, `nucleosomal_state` and `background_state` are the column indexes of the three states in the file (a permutation of 0, 1, 2, given by the training or loading function); `binsize` (default 10).
- **Output:** `files` (the three paths), `intervals` (the number of lines of each file), `bin_count` and `chromosomes`. In the example state 2 is open, so the open track takes the probability in the third column.
- **Example code:**

```python
import tempfile

from dnakit.epigenome import hmmratac_proba_to_bedgraph

with tempfile.TemporaryDirectory() as work:
    result = hmmratac_proba_to_bedgraph(
        "macs3_proba.csv",
        f"{work}/sample",
        open_state=2,
        nucleosomal_state=0,
        background_state=1,
        binsize=10,
        verbose=False,
    )
    print(result["bin_count"], result["chromosomes"], result["intervals"])
    with open(result["files"]["open"], encoding="utf-8") as handle:
        print(handle.read().strip())
```

- **Example output:**

```text
13 ['chr1', 'chr2'] {'open': 12, 'nucleosomal': 15, 'background': 9}
chr1	0	10	0.00000
chr1	10	30	0.10000
chr1	30	50	0.80000
chr1	50	70	0.10000
chr1	70	90	0.00000
chr1	90	100	0.10000
chr1	100	110	0.80000
chr1	110	120	0.10000
chr2	0	30	0.00000
chr2	30	50	0.10000
chr2	50	60	0.80000
chr2	60	70	0.10000
```

### 8) `dnakit.epigenome.hmmratac_states_path()`

- **Purpose:** Take the most probable state of each bin and merge neighbouring bins of the same state into a states path.
- **Method:** Checks the file and calls `hmmratac_cmd.generate_states_path`. A tie in probability goes to the first of open, nucleosomal, background. Gaps between bins and the part before the first bin of each chromosome are background; background entries on either side of a gap are not merged.
- **API:** `dnakit.epigenome.hmmratac_states_path(proba_path: 'str | Path', *, open_state: 'int', nucleosomal_state: 'int', background_state: 'int', binsize: 'int' = 10, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** `proba_path`, `open_state`, `nucleosomal_state`, `background_state` and `binsize` as in the previous function.
- **Output:** `states` is a list of `{chrom, start, end, state}` with `state` equal to `open`, `nuc` or `bg`; `state_counts` and `state_bases` give the number of entries and bases per state; also `bin_count` and `chromosomes`. It can be passed straight to the next two functions.
- **Example code:**

```python
from dnakit.epigenome import hmmratac_states_path

result = hmmratac_states_path(
    "macs3_proba.csv",
    open_state=2,
    nucleosomal_state=0,
    background_state=1,
    binsize=10,
    verbose=False,
)
print(result["state_counts"], result["state_bases"])
for item in result["states"]:
    print(item["chrom"], item["start"], item["end"], item["state"])
```

- **Example output:**

```text
{'open': 3, 'nuc': 6, 'bg': 6} {'open': 40, 'nuc': 60, 'bg': 90}
chr1 0 10 bg
chr1 10 20 bg
chr1 20 30 nuc
chr1 30 50 open
chr1 50 60 nuc
chr1 60 70 bg
chr1 70 90 bg
chr1 90 100 nuc
chr1 100 110 open
chr1 110 120 nuc
chr2 0 30 bg
chr2 30 40 bg
chr2 40 50 nuc
chr2 50 60 open
chr2 60 70 nuc
```

### 9) `dnakit.epigenome.write_hmmratac_states_bed()`

- **Purpose:** Write the open and nucleosomal entries of a states path as a BED file.
- **Method:** Calls `hmmratac_cmd.save_states_bed`: one line `chromosome start end state` per non-background entry. The input can be the list of dicts returned by `hmmratac_states_path` or lists of `(chromosome, start, end, state)`; entries must be sorted along each chromosome without overlap. An existing file is never overwritten.
- **API:** `dnakit.epigenome.write_hmmratac_states_bed(states: 'Sequence[Mapping[str, Any] | Sequence[Any]]', output_path: 'str | Path', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** `states`; `output_path` is the BED file to write.
- **Output:** `output`, `records` (the number of lines written), `state_counts` and `skipped_background` (the number of background entries left out).
- **Example code:**

```python
import tempfile

from dnakit.epigenome import hmmratac_states_path, write_hmmratac_states_bed

path = hmmratac_states_path(
    "macs3_proba.csv", open_state=2, nucleosomal_state=0, background_state=1, verbose=False
)
with tempfile.TemporaryDirectory() as work:
    result = write_hmmratac_states_bed(path["states"], f"{work}/states.bed", verbose=False)
    with open(result["output"], encoding="utf-8") as handle:
        print(handle.read().strip())
print(result["records"], result["state_counts"], result["skipped_background"])
```

- **Example output:**

```text
chr1	20	30	nuc
chr1	30	50	open
chr1	50	60	nuc
chr1	90	100	nuc
chr1	100	110	open
chr1	110	120	nuc
chr2	40	50	nuc
chr2	50	60	open
chr2	60	70	nuc
9 {'open': 3, 'nuc': 6} 6
```

### 10) `dnakit.epigenome.write_hmmratac_accessible_regions()`

- **Purpose:** Write the open entries that sit between nucleosomal entries in a states path as a narrowPeak file, scored from a bedGraph.
- **Method:** Checks the bedGraph first (sorted, with intervals that touch), then calls `hmmratac_cmd.save_accessible_regions`. An open entry must have a nucleosomal entry right before and right after it, and the three together must span more than `min_length` bases; the score is ten times the highest bedGraph value in the region, truncated to an integer, the summit is the middle of the interval with that value, and fc, pscore and qscore are written as 0. The MACS3 pipeline scores with a fold-change track, and the open-state probability track works as well. Note: in MACS3 3.0.5 `save_accessible_regions` slices its result list with `[:-1]`, so the last qualifying region of the states list is never written (the example has 3 qualifying regions and 2 are written); this function returns the upstream result as it is.
- **API:** `dnakit.epigenome.write_hmmratac_accessible_regions(states: 'Sequence[Mapping[str, Any] | Sequence[Any]]', score_bedgraph: 'str | Path', output_path: 'str | Path', *, min_length: 'int' = 100, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** `states`; `score_bedgraph` is the bedGraph used for scoring; `output_path` is the narrowPeak file to write; `min_length` (default 100, at least 0).
- **Output:** `regions` is a list of `{chrom, start, end, name, score, summit}`, plus `count`, `output`, `min_length` and `score_bedgraph`.
- **Example code:**

```python
import tempfile

from dnakit.epigenome import (
    hmmratac_states_path,
    write_hmmratac_accessible_regions,
)

path = hmmratac_states_path(
    "macs3_proba.csv", open_state=2, nucleosomal_state=0, background_state=1, verbose=False
)
with tempfile.TemporaryDirectory() as work:
    result = write_hmmratac_accessible_regions(
        path["states"], "macs3_open.bdg", f"{work}/regions.narrowPeak", min_length=25, verbose=False
    )
    with open(result["output"], encoding="utf-8") as handle:
        print(handle.read().strip())
print(result["count"])
for region in result["regions"]:
    print(region["name"], region["chrom"], region["start"], region["end"], region["score"], region["summit"])
```

- **Example output:**

```text
chr1	30	50	MACS_peak_1	8	.	0	0	0	10
chr1	100	110	MACS_peak_2	8	.	0	0	0	5
2
MACS_peak_1 chr1 30 50 8 40
MACS_peak_2 chr1 100 110 8 105
```

### Functions not wrapped

The following 10 ledger rows get no wrapper of their own; the status and the reason of each follow. Six of them are names in `hmmratac_cmd` for the very same objects as in `HMMR_HMM` and `HMMR_Signal_Processing`, checked by identity.

- `MACS3.Commands.hmmratac_cmd.generate_digested_signals` (`duplicate-of` → `digest_hmmratac_fragments`): the same object as `HMMR_Signal_Processing.generate_digested_signals` (identity checked), covered by `digest_hmmratac_fragments`.
- `MACS3.Commands.hmmratac_cmd.extract_signals_from_regions` (`duplicate-of` → `extract_hmmratac_training_data`): the same object as `HMMR_Signal_Processing.extract_signals_from_regions` (identity checked), covered by `extract_hmmratac_training_data`.
- `MACS3.Commands.hmmratac_cmd.hmm_training` (`duplicate-of` → `train_hmmratac_model`): the same object as `HMMR_HMM.hmm_training` (identity checked), covered by `train_hmmratac_model`.
- `MACS3.Commands.hmmratac_cmd.hmm_predict` (`duplicate-of` → `predict_hmmratac_states`): the same object as `HMMR_HMM.hmm_predict` (identity checked), covered by `predict_hmmratac_states`.
- `MACS3.Commands.hmmratac_cmd.hmm_model_init` (`duplicate-of` → `load_hmmratac_model`): the same object as `HMMR_HMM.hmm_model_init` (identity checked), covered by `load_hmmratac_model`.
- `MACS3.Commands.hmmratac_cmd.hmm_model_save` (`duplicate-of` → `save_hmmratac_model`): the same object as `HMMR_HMM.hmm_model_save` (identity checked), covered by `save_hmmratac_model`.
- `MACS3.IO.OutputWriter.diag_write` (`internal`): the report writer of the callpeak diagnostic table; no sequence or interval logic.
- `MACS3.Signal.PeakModel.smooth` (`not-dna`): generic window smoothing of an array; not DNA-specific.
- `MACS3.Signal.SignalProcessing.savitzky_golay_order2_deriv1` (`not-dna`): generic Savitzky-Golay derivative filter; not DNA-specific.
- `MACS3.Signal.SignalProcessing.savitzky_golay` (`not-dna`): generic Savitzky-Golay smoothing filter; not DNA-specific.

## Additional deepTools functions (GC-bias tables, synthetic fingerprint JSD and bigWig tile bedGraph)

This section adds 4 deepTools functions in the new module `deeptools_review_functions.py` of `dnakit.epigenome`. The examples import from the module file; once exported from `dnakit.epigenome`, the same functions are callable as `dnakit.epigenome.<name>()`. They call deepTools' Python functions directly, do not run the command line, and run in a single process, so results do not depend on process scheduling.

Two upstream details decide how they are wrapped. First, `computeGCBias.tabulateGCcontent` and `countReadsPerGC` read the module-level dict `global_vars`, which only the command line sets; the wrappers build that dict as the command line does, install it for the call only, and restore the previous state afterwards, also when the call fails. Second, the BAM branch of `writeBedGraph_bam_and_bw.writeBedGraph` calls `getCoverageFromBam` and `getSmoothRange`, which deepTools 4.0.0 does not define, so BAM input raises `NameError`; only the bigWig branch works, so the matching wrapper takes bigWig files and is named `bigwig_bedgraph_tiles`. deepTools calls `sys.exit` on some bad inputs; those exits are converted into `ConfigurationError`. The wrappers check the inputs first and never overwrite an existing output file. The return value is a dict with `status`, `operation` and `provenance`; coordinates are 0-based and right-open.

The example data is made up by hand to show the formats and the calculation; it does not stand for a real experiment. The output files of the examples go to a temporary directory, so the examples can be run repeatedly.

### Example data

A small helper that only makes the input files of the examples; it is not part of DNAKit. `write_twobit` writes a 2bit file from a FASTA file of plain ACGT letters and `sam_to_bam` turns a SAM file into a sorted, indexed BAM file.

**Example file:** `deeptools_fixtures.py`

```python
import os
import struct

import pysam


def write_twobit(fasta_path, twobit_path):
    """Write a one-chromosome 2bit file from a FASTA file of plain ACGT letters."""
    name, *rest = open(fasta_path).read().split()
    sequence = "".join(rest)
    codes = {"T": 0, "C": 1, "A": 2, "G": 3}
    packed = bytearray((len(sequence) + 3) // 4)
    for index, base in enumerate(sequence):
        packed[index // 4] |= codes[base] << (6 - 2 * (index % 4))
    record = struct.pack("<IIII", len(sequence), 0, 0, 0) + bytes(packed)
    raw = name[1:].encode()
    header = struct.pack("<IIII", 0x1A412743, 0, 1, 0) + bytes([len(raw)]) + raw
    header += struct.pack("<I", len(header) + 4)
    with open(twobit_path, "wb") as handle:
        handle.write(header + record)


def sam_to_bam(sam_path, bam_path):
    """Convert a SAM file to a sorted and indexed BAM file."""
    pysam.view("-b", "-o", bam_path + ".unsorted", sam_path, catch_stdout=False)
    pysam.sort("-o", bam_path, bam_path + ".unsorted")
    pysam.index(bam_path)
    os.remove(bam_path + ".unsorted")
```

A 40 bp chromosome made of ten 4-base blocks whose GC counts are 2, 0, 4, 0, 2, 0, 3, 1, 4, 1.

**Example file:** `deeptools_genome.fa`

```text
>chr1
ACGTAAAAGGCCATATGCATTTTTCCCAAGATGCGCACTA
```

17 reads of 5 bp, 14 on the forward strand and 3 on the reverse strand. The positions sampled below at a step of 4 bp are 0, 4, 8, ..., 36; 11 forward reads start on them (2 at position 4, 3 at position 24, 2 at position 32, and 1 each at 0, 8, 16 and 28). The forward reads starting at 2, 5 and 17 are off the sampled positions, and the 3 reverse reads are not counted as forward starts.

**Example file:** `deeptools_reads.sam`

```text
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:40
f0	0	chr1	1	60	5M	*	0	0	AAAAA	IIIII
f11	0	chr1	3	60	5M	*	0	0	AAAAA	IIIII
f1	0	chr1	5	60	5M	*	0	0	AAAAA	IIIII
f2	0	chr1	5	60	5M	*	0	0	AAAAA	IIIII
f12	0	chr1	6	60	5M	*	0	0	AAAAA	IIIII
r0	16	chr1	7	60	5M	*	0	0	AAAAA	IIIII
f3	0	chr1	9	60	5M	*	0	0	AAAAA	IIIII
f4	0	chr1	17	60	5M	*	0	0	AAAAA	IIIII
f13	0	chr1	18	60	5M	*	0	0	AAAAA	IIIII
r1	16	chr1	21	60	5M	*	0	0	AAAAA	IIIII
r2	16	chr1	22	60	5M	*	0	0	AAAAA	IIIII
f5	0	chr1	25	60	5M	*	0	0	AAAAA	IIIII
f6	0	chr1	25	60	5M	*	0	0	AAAAA	IIIII
f7	0	chr1	25	60	5M	*	0	0	AAAAA	IIIII
f8	0	chr1	29	60	5M	*	0	0	AAAAA	IIIII
f10	0	chr1	33	60	5M	*	0	0	AAAAA	IIIII
f9	0	chr1	33	60	5M	*	0	0	AAAAA	IIIII
```

The bedGraph sources of two bigWig files (the example converts them with `bedgraph_to_bigwig`). The first is 2 (0-50) and 4 (50-100) on chr1 and covers only 0-30 on chr2, with value 6; the second is 1 (0-25) and 3 (25-100) on chr1 and 1 over the whole of chr2.

**Example file:** `deeptools_a.bdg`

```text
chr1	0	50	2
chr1	50	100	4
chr2	0	30	6
```

**Example file:** `deeptools_b.bdg`

```text
chr1	0	25	1
chr1	25	100	3
chr2	0	60	1
```

### 1) `dnakit.epigenome.tabulate_gc_fragment_counts()`

- **Purpose:** Take a fragment of `fragment_length` bases every `step_size` positions along the genome, count its G and C bases, and count the forward reads that start at those positions. The result is the frequency table of `computeGCBias` (the input of `correctGCBias`), one row per GC count.
- **Method:** Calls deepTools' `computeGCBias.tabulateGCcontent`. The wrapper first reads the chromosome lengths of the 2bit file and the mapped read count from the BAM index, then builds `global_vars` as the command line does: `reads_per_bp` is the mapped read count divided by `effective_genome_size`; `max_reads` is the upper quantile, at probability `1 / sample_size`, of a Poisson distribution with mean `4 × reads_per_bp × fragment_length`, and a position with at least `max_reads` forward read starts is skipped as a peak (it still counts as a sampled position). The default step is `max(int(2bit total length / sample_size), 1)`. Chromosome names of the 2bit and the BAM are matched with `tbitToBamChrName`; a BAM without the `chr` prefix is matched too, whereas the `computeGCBias` command line fails with `UnboundLocalError` for such a pair. Fragments with more than 5% N bases are not counted, and only forward reads are counted.
- **API:** `dnakit.epigenome.tabulate_gc_fragment_counts(twobit_path: 'str | Path', bam_path: 'str | Path', fragment_length: 'int', effective_genome_size: 'int', *, sample_size: 'int' = 50000000, step_size: 'int | None' = None, region: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** `twobit_path`; `bam_path` (needs a `.bai` index and at least one mapped read); `fragment_length` (1 to 1000000); `effective_genome_size` (the mappable genome size); `sample_size` (default 50000000, at least 2; it sets the default step and `max_reads`); `step_size` (derived from `sample_size` when omitted); `region` (`chromosome:start:end` with the chromosome name of the BAM; only that region is sampled). At most 200000000 positions are sampled.
- **Output:** `gc_table` has one entry `{gc_bases, read_starts, genome_positions, read_ratio}` per GC count: `genome_positions` is the number of sampled positions with that GC count, `read_starts` the forward read starts on them, and `read_ratio` their quotient divided by the same quotient over all sampled positions, a bias relative to the overall density. Also `positions_sampled`, `read_starts_counted`, `step_size`, `mapped_reads`, `reads_per_bp`, `max_reads` and `min_reads`. In the example the 10 sampled positions have 3, 2, 2, 1 and 2 positions with GC counts 0 to 4 and 2, 1, 2, 3 and 3 read starts, so the ratio at GC 3 is (3 / 1) ÷ (11 / 10) = 2.727.
- **Example code:**

```python
from deeptools_fixtures import sam_to_bam, write_twobit
from dnakit.epigenome import tabulate_gc_fragment_counts

write_twobit("deeptools_genome.fa", "deeptools_genome.2bit")
sam_to_bam("deeptools_reads.sam", "deeptools_reads.bam")
result = tabulate_gc_fragment_counts(
    "deeptools_genome.2bit",
    "deeptools_reads.bam",
    4,
    40,
    sample_size=10,
    verbose=False,
)
print(result["step_size"], result["mapped_reads"], result["max_reads"])
print(result["positions_sampled"], result["read_starts_counted"])
for row in result["gc_table"]:
    print(
        row["gc_bases"],
        row["genome_positions"],
        row["read_starts"],
        round(row["read_ratio"], 3),
    )
```

- **Example output:**

```text
4 17 10.0
10 11
0 3 2 0.606
1 2 1 0.455
2 2 2 0.909
3 1 3 2.727
4 2 3 1.364
```

### 2) `dnakit.epigenome.count_reads_per_gc_bins()`

- **Purpose:** Take windows of `region_size` bases along the genome (one every `step_size` positions) and report the GC fraction of each window and the number of reads that overlap it (both strands). This is the data behind the plot of `computeGCBias --biasPlot`.
- **Method:** Calls deepTools' `computeGCBias.countReadsPerGC`; `global_vars` is handled as in the previous function. The default step is `max(int(2bit total length / sample_size), 1) × 10`, the same as the command line uses for the plot. Windows that would run past the chromosome end are not made, and windows with more than 5% N bases are left out. Windows come in genome order; upstream returns only the read count and GC, not the window coordinates.
- **API:** `dnakit.epigenome.count_reads_per_gc_bins(twobit_path: 'str | Path', bam_path: 'str | Path', effective_genome_size: 'int', *, region_size: 'int' = 300, sample_size: 'int' = 50000000, step_size: 'int | None' = None, region: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** `twobit_path`; `bam_path`; `effective_genome_size` (it only sets the chunk size); `region_size` (default 300); `sample_size` (default 50000000); `step_size`; `region` (as above). At most 1000000 windows.
- **Output:** `windows` is a list of `{reads, gc_fraction}`, `window_count` the number of windows and `mean_reads` the mean read count (`None` without windows). In the example the chromosome is 40 bp, the window 8 bp and the step 4, so the starts 0 to 32 give 9 windows. The first window is `ACGTAAAA` with GC 2/8 = 0.25; the forward reads starting at 0, 2, 4, 4 and 5 and the reverse read starting at 6 overlap it, 6 reads in all.
- **Example code:**

```python
from deeptools_fixtures import sam_to_bam, write_twobit
from dnakit.epigenome import count_reads_per_gc_bins

write_twobit("deeptools_genome.fa", "deeptools_genome.2bit")
sam_to_bam("deeptools_reads.sam", "deeptools_reads.bam")
result = count_reads_per_gc_bins(
    "deeptools_genome.2bit",
    "deeptools_reads.bam",
    40,
    region_size=8,
    sample_size=40,
    step_size=4,
    verbose=False,
)
print(result["window_count"], round(result["mean_reads"], 3))
for window in result["windows"]:
    print(window["reads"], window["gc_fraction"])
```

- **Example output:**

```text
9 5.222
6 0.25
7 0.5
5 0.5
3 0.25
4 0.25
7 0.375
6 0.5
6 0.625
3 0.625
```

### 3) `dnakit.epigenome.fingerprint_synthetic_jsd()`

- **Purpose:** Measure how far a vector of per-bin coverage is from a Poisson distribution with the same mean, which is the synthetic JSD that `plotFingerprint` uses when no input sample is given; the larger the value, the more the coverage is concentrated in few bins.
- **Method:** Calls deepTools' `plotFingerprint.getSyntheticJSD`: the bin coverage (integer part; bins with no coverage are dropped; values of `MAXLEN` = 10000000 or more go into the last class) is counted into a histogram, a Poisson distribution with the mean coverage as parameter and the total coverage as size serves as the null, and `getJSDcommon` computes the Jensen-Shannon distance between the two cumulative coverage curves. Nothing is random, so the result is deterministic. One call handles arrays of 10 million entries and takes about 3 to 4 seconds.
- **API:** `dnakit.epigenome.fingerprint_synthetic_jsd(coverage: 'Sequence[float]', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** `coverage` is a list of non-negative numbers (1 to 10000000 values, at least one of them 1 or more). Booleans, negative numbers, NaN and infinity are rejected.
- **Output:** `bins`, `mean_coverage`, `bins_above_maximum` (bins at `MAXLEN` or more) and `jsd`. In the example `even`, close to Poisson, has distance 0.209 and `peaked`, with most of its signal in 2 bins, has 0.292.
- **Example code:**

```python
from dnakit.epigenome import fingerprint_synthetic_jsd

even = [4, 5, 3, 4, 6, 4, 5, 3, 4, 5]
peaked = [0, 0, 1, 12, 15, 1, 0, 0, 0, 1]
for coverage in (even, peaked):
    result = fingerprint_synthetic_jsd(coverage, verbose=False)
    print(result["bins"], result["mean_coverage"], result["jsd"])
```

- **Example output:**

```text
10 4.3 0.209461503453211
10 3.0 0.292041559465403
```

### 4) `dnakit.epigenome.bigwig_bedgraph_tiles()`

- **Purpose:** Average one or more bigWig files over tiles of fixed size, take the mean or the sum of the scaled tile values, and write a bedGraph; useful for bigWig averages, differences or scaling.
- **Method:** Calls deepTools' `writeBedGraph_bam_and_bw.writeBedGraph` (the engine behind `bigwigAverage` and `bigwigCompare`). Each file is averaged over the bases of each tile (with `missing_data_as_zero=True` missing values count as 0, otherwise a tile that includes missing values is NaN and is not written); the tile value of every file is multiplied by its entry of `scale_factors` and then reduced by `reduction`, the mean (`"mean"`) or the sum (`"sum"`); neighbouring tiles with equal values are merged into one line, and `fixed_step=True` writes one line per tile. The BAM branch of upstream raises `NameError` in deepTools 4.0.0, so only bigWig files are accepted and a BAM file gives a `ConfigurationError`. The start of a region is moved down to a tile boundary, and upstream does not write a final tile whose value is 0. Chromosomes with the same name need the same length in all files, and only the chromosomes shared by all files are processed.
- **API:** `dnakit.epigenome.bigwig_bedgraph_tiles(bigwig_paths: 'Sequence[str | Path]', output_path: 'str | Path', *, tile_size: 'int' = 50, reduction: 'str' = 'mean', scale_factors: 'Sequence[float] | None' = None, region: 'str | None' = None, missing_data_as_zero: 'bool' = True, fixed_step: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** `bigwig_paths` (1 to 32 bigWig files); `output_path` (must not exist; nothing is overwritten); `tile_size` (default 50); `reduction` (`"mean"` or `"sum"`); `scale_factors` (one number per file, all 1 by default); `region` (`chromosome:start:end`); `missing_data_as_zero` (default `True`); `fixed_step` (default `False`).
- **Output:** `output` (the bedGraph written), `intervals` (the number of lines), `chromosomes`, `first_intervals` (the first 20 lines as `{chrom, start, end, value}`) and the parameters used. In the example the first tile [0, 25) of chr1 is 2 and 1 in the two files, mean 1.5; the tile [25, 50) of chr2 has the value 6 on only 5 bases in the first file, so its mean is 1.2, and averaged with 1 of the second file it gives 1.1; the second call uses `reduction="sum"` and `scale_factors=[1, -1]` to get the difference of the two files.
- **Example code:**

```python
import tempfile

from dnakit.epigenome import bedgraph_to_bigwig
from dnakit.epigenome import bigwig_bedgraph_tiles

sizes = {"chr1": 100, "chr2": 60}
with tempfile.TemporaryDirectory() as work:
    bedgraph_to_bigwig(sizes, ["deeptools_a.bdg"], f"{work}/a.bw", verbose=False)
    bedgraph_to_bigwig(sizes, ["deeptools_b.bdg"], f"{work}/b.bw", verbose=False)
    tracks = [f"{work}/a.bw", f"{work}/b.bw"]
    mean = bigwig_bedgraph_tiles(
        tracks, f"{work}/mean.bdg", tile_size=25, verbose=False
    )
    print(open(mean["output"]).read().strip())
    print(mean["intervals"], mean["chromosomes"])
    difference = bigwig_bedgraph_tiles(
        tracks,
        f"{work}/difference.bdg",
        tile_size=25,
        reduction="sum",
        scale_factors=[1, -1],
        verbose=False,
    )
    print(open(difference["output"]).read().strip())
```

- **Example output:**

```text
chr1	0	25	1.5
chr1	25	50	2.5
chr1	50	100	3.5
chr2	0	25	3.5
chr2	25	50	1.1
chr2	50	60	0.5
6 ['chr1', 'chr2']
chr1	0	25	1
chr1	25	50	-1
chr1	50	100	1
chr2	0	25	5
chr2	25	50	0.2
chr2	50	60	-1
```

### Functions not wrapped

The following 6 ledger rows get no wrapper of their own; the status and the reason of each follow.

- `deeptools.bamPEFragmentSize.getDensity` (`plot`): a plotting helper that draws the fragment-length histogram with matplotlib's `plt.hist` and `plt.clf`.
- `deeptools.bigwigAverage.getType` (`internal`): guesses the track file type from the file extension.
- `deeptools.bigwigCompare.getType` (`internal`): the same extension-based type guess as `bigwigAverage.getType`.
- `deeptools.computeGCBias.smooth` (`internal`): a window-average smoothing helper that upstream marks as not used.
- `deeptools.heatmapper.trimZones` (`internal`): a bookkeeping helper for the heatmap zones, serving the plot bins.
- `deeptools.heatmapper.computeSilhouetteScore` (`not-dna`): a clustering statistic on a distance matrix; no DNA meaning.
