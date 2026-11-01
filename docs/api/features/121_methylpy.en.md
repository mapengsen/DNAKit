# methylpy: bisulfite alignment, methylation calling and DMR finding

This group of functions in `dnakit.epigenome` is a set of adapters over [methylpy](https://github.com/yupenghe/methylpy), the Python pipeline for whole-genome bisulfite sequencing analysis. They build C-to-T and G-to-A converted references, align single-end and paired-end bisulfite reads, merge and de-duplicate alignments, trim adapters, call methylated cytosines into allc tables, filter and compress those tables, convert them to bigWig tracks, and find differentially methylated regions (DMRs). Each function calls the matching methylpy function and forwards only the arguments you set, so methylpy's own defaults apply to the rest. DNAKit adds three things: the external programs are checked first (a missing program raises `BackendUnavailableError` with an install hint instead of methylpy calling `sys.exit`), methylpy failures are raised as `BackendExecutionError`, and the result is a JSON-serialisable dictionary. Many of these functions wrap alignment and indexing steps that run `samtools`, `bowtie2`, `bgzip` and similar programs, and write several files next to the output prefix.

The package metadata of methylpy 1.4.7 ships the Apache License 2.0, so there is no copyleft caveat for the Python package. The bundled `run_rms_tests.out` binary used by the DMR functions links against the GNU Scientific Library (`libgsl.so.0`).

The dependencies are optional and are imported only when the matching function is called. Install methylpy with `pip install "dnakit[bio-methylpy]"` and `cutadapt` with `pip install "dnakit[bio-cutadapt]"`; `samtools`, `bowtie2` or `bowtie`, `bgzip`, `tabix`, `wigToBigWig`, Picard and GSL are external programs and libraries that you install separately. The examples import from `dnakit.epigenome`.

This environment has methylpy and cutadapt but no `samtools`, `bowtie2`, `bgzip`, `tabix`, `wigToBigWig` or GSL. Five examples therefore run the real computation (`methylpy_filter_allc_files`, the two `find_multi_mappers` functions and the two `quality_trim` functions). The other examples show the function's true behaviour when a program or an input is missing: they point the program path at a folder that does not exist (`path_to_samtools="no_such_dir"`), so the output is the same on every machine, and each of them is followed by an unexecuted "Usage with the dependency installed" block with a realistic call.

## Example data

The examples use three tiny declared files: a five-row allc table, a one-sequence reference FASTA and two FASTQ reads (the first carries a 3' adapter). Other inputs are created by the example itself, and outputs are written to a temporary directory or never created.

**Example file:** `allc.tsv`

```text
chr1	100	+	CGG	10	12	1
chr1	105	-	CHH	2	20	0
chr1	300	+	CAG	3	4	1
chr2	50	+	CGT	0	8	0
chr2	70	-	CCC	5	6	1
```

**Example file:** `ref.fa`

```text
>chr1
ACGTACGTTTGGAACCGGTTACGTACGTTTGGAACCGGTT
```

**Example file:** `reads.fq`

```text
@r1
ACGTACGTTTGGAAGATCGGAAGAGCACACGTCTGAAAA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@r2
TTGGTTGGTTAAGGCCAAGGTTTTAAGGCCTTAAGGAATT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
```

## 1) `dnakit.epigenome.methylpy_bgzip_allc_file()`

- **Purpose:** Compress an allc methylation table with bgzip and index it with tabix.
- **Method:** Check that the file exists, check `bgzip` and `tabix`, then call `methylpy.utilities.bgzip_allc_file`. Upstream writes the block-compressed copy `<allc_file>.bgz` (for an input that already ends in `.gz` it is `<name>.bgz` in place of `.gz`) and indexes it with `tabix -s 1 -b 2 -e 2`, which creates the `.tbi` file. The input file is kept.
- **API:** `dnakit.epigenome.methylpy_bgzip_allc_file(allc_file: 'str', *, path_to_bgzip: 'str | None' = None, path_to_tabix: 'str | None' = None, buffer_line_number: 'int | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs methylpy (extra `dnakit[bio-methylpy]`) and the external programs `bgzip` and `tabix` (htslib) on PATH (or passed through the matching `path_to_...` argument). Without them the call raises `BackendUnavailableError`.
- **Input:** `allc_file`, an existing allc table. `path_to_bgzip` and `path_to_tabix` are folders that contain the programs (methylpy appends the program name). `buffer_line_number` sets how many lines are read at a time.
- **Output:** A dictionary with `status` (`"completed"`), `backend` (`"methylpy"`), the methylpy `version`, `function` (the methylpy function that ran), `outputs` (an empty list: the function declares no output path) and `arguments` (the sorted names of the arguments that were forwarded).
- **Example code:**

```python
from pathlib import Path

from dnakit.epigenome import methylpy_bgzip_allc_file
from dnakit.exceptions import BackendUnavailableError

try:
    methylpy_bgzip_allc_file("allc.tsv", path_to_bgzip="no_such_dir", verbose=False)
except BackendUnavailableError as error:
    print(type(error).__name__, "-", error.message)
```

- **Example output:**

```text
BackendUnavailableError - bgzip was not found at no_such_dir/bgzip.
```

Usage with the dependency installed:

```python
from dnakit.epigenome import methylpy_bgzip_allc_file

result = methylpy_bgzip_allc_file("allc_sample_chr1.tsv", verbose=False)
print(result["status"])  # allc_sample_chr1.tsv.bgz and allc_sample_chr1.tsv.bgz.tbi are written next to the input
```

## 2) `dnakit.epigenome.methylpy_convert_allc_to_bigwig()`

- **Purpose:** Convert an allc table to a binned methylation bigWig track.
- **Method:** Check the input files, check `samtools` and `wigToBigWig`, then call `methylpy.utilities.convert_allc_to_bigwig`, Upstream indexes the FASTA with `samtools faidx` when `<reference_fasta>.fai` is missing, writes temporary `<output_file>.chrom_size` and `<output_file>.wig` files, and runs `wigToBigWig`. In the wig file each bin of `bin_size` bases (default 100) holds the methylation level of the chosen context, the sum of methylated reads divided by the sum of covering reads, and the temporary files are removed afterwards.
- **API:** `dnakit.epigenome.methylpy_convert_allc_to_bigwig(input_allc_file: 'str', output_file: 'str', reference_fasta: 'str', *, mc_type: 'str | None' = None, bin_size: 'int | None' = None, path_to_wigToBigWig: 'str | None' = None, path_to_samtools: 'str | None' = None, min_bin_sites: 'int | None' = None, min_bin_cov: 'int | None' = None, max_site_cov: 'int | None' = None, min_site_cov: 'int | None' = None, remove_chr_prefix: 'bool | None' = None, add_chr_prefix: 'bool | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs methylpy (extra `dnakit[bio-methylpy]`) and the external programs `samtools` and `wigToBigWig` (UCSC) on PATH (or passed through the matching `path_to_...` argument). Without them the call raises `BackendUnavailableError`.
- **Input:** `input_allc_file` and `reference_fasta` must exist; `output_file` is the bigWig path. Optional: `mc_type`, `bin_size`, `min_bin_sites`, `min_bin_cov`, `min_site_cov`, `max_site_cov`, `remove_chr_prefix`, `add_chr_prefix`, `path_to_wigToBigWig`, `path_to_samtools`.
- **Output:** A dictionary with `status` (`"completed"`), `backend` (`"methylpy"`), the methylpy `version`, `function` (the methylpy function that ran), `outputs` (the declared output paths that exist after the call) and `arguments` (the sorted names of the arguments that were forwarded).
- **Example code:**

```python
from pathlib import Path

from dnakit.epigenome import methylpy_convert_allc_to_bigwig
from dnakit.exceptions import BackendUnavailableError

try:
    methylpy_convert_allc_to_bigwig("allc.tsv", "cg.bw", "ref.fa", path_to_samtools="no_such_dir", verbose=False)
except BackendUnavailableError as error:
    print(type(error).__name__, "-", error.message)
```

- **Example output:**

```text
BackendUnavailableError - samtools was not found at no_such_dir/samtools.
```

Usage with the dependency installed:

```python
from dnakit.epigenome import methylpy_convert_allc_to_bigwig

result = methylpy_convert_allc_to_bigwig(
    "allc_sample_chr1.tsv", "sample_cg.bw", "genome.fa", mc_type="CGN", bin_size=100, verbose=False
)
print(result["outputs"])
```

## 3) `dnakit.epigenome.methylpy_dmrfind()`

- **Purpose:** Find differentially methylated regions between samples from several allc files.
- **Method:** Check that the bundled `run_rms_tests.out` binary exists, check that every allc file exists, then call `methylpy.DMRfind.DMRfind`. Upstream finds differentially methylated sites with a root-mean-square test, merges neighbouring sites within `dmr_max_dist`, and runs the permutation tests through the bundled binary. Upstream stops with "Needs more than one allc file as input" when only one file is given, which DNAKit reports as `BackendExecutionError`.
- **API:** `dnakit.epigenome.methylpy_dmrfind(allc_files: 'Sequence[str]', output_prefix: 'str', *, samples: 'Sequence[str] | None' = None, mc_type: 'str | None' = None, chroms: 'Sequence[str] | None' = None, num_procs: 'int | None' = None, min_cov: 'int | None' = None, dmr_max_dist: 'int | None' = None, sig_cutoff: 'float | None' = None, num_sims: 'int | None' = None, num_sig_tests: 'int | None' = None, seed: 'int | None' = None, keep_temp_files: 'bool | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs methylpy (extra `dnakit[bio-methylpy]`); the bundled `run_rms_tests.out` binary needs the GSL shared libraries `libgsl.so.0` and `libgslcblas.so.0` (`apt install libgsl-dev` or `conda install -c conda-forge gsl`).
- **Input:** `allc_files`, a list of at least two existing allc tables; `output_prefix` for the result files. Optional: `samples` (names, one per file), `mc_type`, `chroms`, `num_procs`, `min_cov`, `dmr_max_dist`, `sig_cutoff`, `num_sims`, `num_sig_tests`, `seed`, `keep_temp_files`. `mc_type` is passed to methylpy unchanged (its default is `CGN`).
- **Output:** The result dictionary described above (`outputs` is an empty list; methylpy writes `<output_prefix>_DMRs.tsv` and related files). Missing input files raise `ConfigurationError`.
- **Example code:**

```python
from dnakit.epigenome import methylpy_dmrfind
from dnakit.exceptions import ConfigurationError

try:
    methylpy_dmrfind(["missing_a.tsv", "missing_b.tsv"], "dmr", samples=["a", "b"], verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, "-", error.message.split(":")[0])
```

- **Example output:**

```text
ConfigurationError - allc_files must be an existing file
```

Usage with the dependency installed:

```python
from dnakit.epigenome import methylpy_dmrfind

result = methylpy_dmrfind(
    ["allc_a.tsv.gz", "allc_b.tsv.gz"], "dmr", samples=["a", "b"], mc_type="CGN", num_procs=2, verbose=False
)
print(result["status"])  # dmr_DMRs.tsv is written
```

## 4) `dnakit.epigenome.methylpy_pe_call_methylated_sites()`

- **Purpose:** Call methylated cytosines from a sorted paired-end BAM file.
- **Method:** Check that the BAM and the reference exist, check `samtools` (and `bgzip` and `tabix` when `bgzip=True`), then call `methylpy.call_mc_pe.call_methylated_sites_pe`. Upstream builds a pileup per chromosome, estimates the non-conversion rate from `unmethylated_control` and applies a binomial test with the `sig_cutoff` adjusted p-value, writing one `allc_<sample>_<chrom>.tsv` file per chromosome.
- **API:** `dnakit.epigenome.methylpy_pe_call_methylated_sites(inputf: 'str', sample: 'str', reference_fasta: 'str', *, unmethylated_control: 'str | None' = None, sig_cutoff: 'float | None' = None, num_procs: 'int | None' = None, num_upstr_bases: 'int | None' = None, num_downstr_bases: 'int | None' = None, generate_mpileup_file: 'bool | None' = None, compress_output: 'bool | None' = None, bgzip: 'bool | None' = None, path_to_bgzip: 'str | None' = None, path_to_tabix: 'str | None' = None, buffer_line_number: 'int | None' = None, min_mapq: 'int | None' = None, min_cov: 'int | None' = None, binom_test: 'bool | None' = None, path_to_samtools: 'str | None' = None, remove_chr_prefix: 'bool | None' = None, add_snp_info: 'bool | None' = None, sort_mem: 'str | None' = None, path_to_files: 'str | None' = None, min_base_quality: 'int | None' = None, keep_temp_files: 'bool | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs methylpy (extra `dnakit[bio-methylpy]`) and the external programs `samtools` (plus `bgzip` and `tabix` when `bgzip=True`) on PATH (or passed through the matching `path_to_...` argument). Without them the call raises `BackendUnavailableError`.
- **Input:** `inputf`, a sorted BAM; `sample`, the name used in the output files; `reference_fasta`, an existing FASTA. Optional: `unmethylated_control` (a `chrom:start-end` region or a number), `sig_cutoff`, `num_procs`, `num_upstr_bases`, `num_downstr_bases`, `min_mapq`, `min_cov`, `min_base_quality`, `binom_test`, `generate_mpileup_file`, `compress_output`, `bgzip`, `remove_chr_prefix`, `add_snp_info`, `sort_mem`, `buffer_line_number`, `path_to_files`, `path_to_samtools`, `path_to_bgzip`, `path_to_tabix`, `keep_temp_files`.
- **Output:** A dictionary with `status` (`"completed"`), `backend` (`"methylpy"`), the methylpy `version`, `function` (the methylpy function that ran), `outputs` (an empty list: the function declares no output path) and `arguments` (the sorted names of the arguments that were forwarded).
- **Example code:**

```python
from pathlib import Path

from dnakit.epigenome import methylpy_pe_call_methylated_sites
from dnakit.exceptions import BackendUnavailableError

Path("reads_pe.bam").touch()
try:
    methylpy_pe_call_methylated_sites(
        "reads_pe.bam", "s1", "ref.fa", path_to_samtools="no_such_dir", verbose=False
    )
except BackendUnavailableError as error:
    print(type(error).__name__, "-", error.message)
```

- **Example output:**

```text
BackendUnavailableError - samtools was not found at no_such_dir/samtools.
```

Usage with the dependency installed:

```python
from dnakit.epigenome import methylpy_pe_call_methylated_sites

result = methylpy_pe_call_methylated_sites(
    "s1_processed_reads_no_clonal.bam", "s1", "genome.fa",
    unmethylated_control="chrC:", num_procs=4, verbose=False,
)
print(result["status"])  # allc_s1_<chrom>.tsv.gz files are written
```

## 5) `dnakit.epigenome.methylpy_pe_find_multi_mappers()`

- **Purpose:** Split paired-end alignments into properly paired, uniquely scored reads, restoring the original read sequence.
- **Method:** Call `methylpy.call_mc_pe.find_multi_mappers_pe`. It reads the SAM file written by the aligner, skips header lines, reads that are not properly paired (flag bit 2), reads below `max(3, min_mapq)` and supplementary alignments, and decodes the bisulfite-converted sequence using the positions encoded after the `!` in the read name. The kept lines are written round-robin to `<output>_sorted_<i>` for `i < num_procs`. The input SAM file is deleted unless `keep_temp_files=True`.
- **API:** `dnakit.epigenome.methylpy_pe_find_multi_mappers(inputf: 'str', output: 'str', *, num_procs: 'int | None' = None, min_mapq: 'int | None' = None, keep_temp_files: 'bool | None' = None, append: 'bool | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs methylpy (extra `dnakit[bio-methylpy]`) only; it runs no external program. The input must be a SAM file whose read names end in `!<encoded positions>`, as written by the methylpy alignment steps.
- **Input:** `inputf`, an existing SAM file; `output`, the prefix of the result files. Optional: `num_procs` (number of output files), `min_mapq`, `keep_temp_files`, `append` (append to existing output files instead of overwriting them).
- **Output:** A dictionary with `status` (`"completed"`), `backend` (`"methylpy"`), the methylpy `version`, `function` (the methylpy function that ran), `outputs` (always an empty list here, because the function declares no single output path: methylpy writes files named from the prefix) and `arguments` (the sorted names of the arguments that were forwarded).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import methylpy_pe_find_multi_mappers

folder = Path(tempfile.mkdtemp())
sam = folder / "pairs.sam"
sam.write_text(
    "@HD\tVN:1.0\n"
    "p1!\t99\tchr1\t10\t42\t8M\t=\t30\t28\tACGTACGT\tIIIIIIII\n"
    "p1!\t147\tchr1\t30\t42\t8M\t=\t10\t-28\tTTGGTTAA\tIIIIIIII\n"
    "p2!\t65\tchr1\t10\t42\t8M\t=\t30\t28\tACGTACGT\tIIIIIIII\n"
)
result = methylpy_pe_find_multi_mappers(
    str(sam), str(folder / "pe"), num_procs=2, keep_temp_files=True, verbose=False
)
print(result["status"], result["outputs"])
for index in range(2):
    rows = (folder / f"pe_sorted_{index}").read_text().splitlines()
    print(index, [row.split("\t")[:2] for row in rows])
```

- **Example output:**

```text
completed []
0 [['p1', '99']]
1 [['p1', '147']]
```

Usage with the dependency installed:

```python
from dnakit.epigenome import methylpy_pe_find_multi_mappers

result = methylpy_pe_find_multi_mappers("s1_forward.sam", "s1_forward", num_procs=4, verbose=False)
# s1_forward_sorted_0 ... s1_forward_sorted_3 are written; s1_forward.sam is deleted
```

## 6) `dnakit.epigenome.methylpy_pe_flip_read2_strand()`

- **Purpose:** Flip the strand of read 2 in every mapped read pair of a BAM file.
- **Method:** Check that the input exists and `samtools` is available, then call `methylpy.call_mc_pe.flip_read2_strand`, which rewrites the BAM so that both reads of a pair report the original strand.
- **API:** `dnakit.epigenome.methylpy_pe_flip_read2_strand(input_file: 'str', output_file: 'str', *, path_to_samtools: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs methylpy (extra `dnakit[bio-methylpy]`) and the external program `samtools` on PATH (or passed through the matching `path_to_...` argument). Without them the call raises `BackendUnavailableError`.
- **Input:** `input_file`, an existing paired-end BAM; `output_file`, the BAM to write. `path_to_samtools` is the folder that contains samtools.
- **Output:** A dictionary with `status` (`"completed"`), `backend` (`"methylpy"`), the methylpy `version`, `function` (the methylpy function that ran), `outputs` (the declared output paths that exist after the call) and `arguments` (the sorted names of the arguments that were forwarded).
- **Example code:**

```python
from pathlib import Path

from dnakit.epigenome import methylpy_pe_flip_read2_strand
from dnakit.exceptions import BackendUnavailableError

Path("pairs.bam").touch()
try:
    methylpy_pe_flip_read2_strand("pairs.bam", "flipped.bam", path_to_samtools="no_such_dir", verbose=False)
except BackendUnavailableError as error:
    print(type(error).__name__, "-", error.message)
```

- **Example output:**

```text
BackendUnavailableError - samtools was not found at no_such_dir/samtools.
```

Usage with the dependency installed:

```python
from dnakit.epigenome import methylpy_pe_flip_read2_strand

result = methylpy_pe_flip_read2_strand("s1_pe.bam", "s1_pe_flipped.bam", verbose=False)
print(result["outputs"])
```

## 7) `dnakit.epigenome.methylpy_pe_merge_bam_files()`

- **Purpose:** Merge several paired-end BAM files into one BAM file with a corrected header.
- **Method:** Check that the inputs exist and `samtools` is available, then call `methylpy.call_mc_pe.merge_bam_files`, which writes a temporary `header.sam` in the current directory (the header of the first file plus one `@RG` line per input, then removed) and runs `samtools merge -r -h header.sam`. Every input file name must contain `.bam`.
- **API:** `dnakit.epigenome.methylpy_pe_merge_bam_files(input_files: 'Sequence[str]', output: 'str', *, path_to_samtools: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs methylpy (extra `dnakit[bio-methylpy]`) and the external program `samtools` on PATH (or passed through the matching `path_to_...` argument). Without them the call raises `BackendUnavailableError`.
- **Input:** `input_files`, a non-empty list of existing BAM files; `output`, the merged BAM; `path_to_samtools`.
- **Output:** A dictionary with `status` (`"completed"`), `backend` (`"methylpy"`), the methylpy `version`, `function` (the methylpy function that ran), `outputs` (the declared output paths that exist after the call) and `arguments` (the sorted names of the arguments that were forwarded).
- **Example code:**

```python
from pathlib import Path

from dnakit.epigenome import methylpy_pe_merge_bam_files
from dnakit.exceptions import BackendUnavailableError

Path("lib1.bam").touch()
Path("lib2.bam").touch()
try:
    methylpy_pe_merge_bam_files(
        ["lib1.bam", "lib2.bam"], "merged.bam", path_to_samtools="no_such_dir", verbose=False
    )
except BackendUnavailableError as error:
    print(type(error).__name__, "-", error.message)
```

- **Example output:**

```text
BackendUnavailableError - samtools was not found at no_such_dir/samtools.
```

Usage with the dependency installed:

```python
from dnakit.epigenome import methylpy_pe_merge_bam_files

result = methylpy_pe_merge_bam_files(
    ["s1_lib1_no_clonal.bam", "s1_lib2_no_clonal.bam"], "s1_merged.bam", verbose=False
)
print(result["outputs"])
```

## 8) `dnakit.epigenome.methylpy_pe_merge_sorted_multimap()`

- **Purpose:** Merge the sorted paired-end alignment files of one library and keep the uniquely mapping read pairs.
- **Method:** Check that the files and the reference exist and `samtools` is available, then call `methylpy.call_mc_pe.merge_sorted_multimap_pe`. Upstream indexes the reference with `samtools faidx` when `<reference_fasta>.fai` is missing, writes a SAM header from the index, merges the sorted files read by read and keeps a read pair only when it occurs once across the files, and converts the result to `<prefix>_processed_reads.bam`.
- **API:** `dnakit.epigenome.methylpy_pe_merge_sorted_multimap(current_library: 'str', files: 'Sequence[str]', prefix: 'str', reference_fasta: 'str', *, path_to_samtools: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs methylpy (extra `dnakit[bio-methylpy]`) and the external program `samtools` on PATH (or passed through the matching `path_to_...` argument). Without them the call raises `BackendUnavailableError`.
- **Input:** `current_library`, the library ID; `files`, the `_sorted_<i>` files from `methylpy_pe_find_multi_mappers`; `prefix`, the output prefix; `reference_fasta`, an existing FASTA. `path_to_samtools` is optional.
- **Output:** A dictionary with `status` (`"completed"`), `backend` (`"methylpy"`), the methylpy `version`, `function` (the methylpy function that ran), `outputs` (an empty list: the function declares no output path) and `arguments` (the sorted names of the arguments that were forwarded).
- **Example code:**

```python
from pathlib import Path

from dnakit.epigenome import methylpy_pe_merge_sorted_multimap
from dnakit.exceptions import BackendUnavailableError

Path("lib1_sorted_0").touch()
try:
    methylpy_pe_merge_sorted_multimap(
        "lib1", ["lib1_sorted_0"], "s1_lib1", "ref.fa", path_to_samtools="no_such_dir", verbose=False
    )
except BackendUnavailableError as error:
    print(type(error).__name__, "-", error.message)
```

- **Example output:**

```text
BackendUnavailableError - samtools was not found at no_such_dir/samtools.
```

Usage with the dependency installed:

```python
from dnakit.epigenome import methylpy_pe_merge_sorted_multimap

result = methylpy_pe_merge_sorted_multimap(
    "lib1", ["s1_lib1_sorted_0", "s1_lib1_sorted_1"], "s1_lib1", "genome.fa", verbose=False
)
# s1_lib1_processed_reads.bam is written
```

## 9) `dnakit.epigenome.methylpy_pe_merge_sorted_multimap_max_mapq()`

- **Purpose:** Merge the sorted paired-end alignment files of one library and keep the best-MAPQ pair of each multiply mapped read.
- **Method:** Check that the files and the reference exist and `samtools` is available, then call `methylpy.call_mc_pe.merge_sorted_multimap_pe_max_mapq`. Upstream indexes the reference with `samtools faidx` when `<reference_fasta>.fai` is missing, writes a SAM header from the index, merges the sorted files read by read and keeps, for a read that occurs more than once, the pair with the highest mapping quality, and converts the result to `<prefix>_processed_reads.bam`.
- **API:** `dnakit.epigenome.methylpy_pe_merge_sorted_multimap_max_mapq(current_library: 'str', files: 'Sequence[str]', prefix: 'str', reference_fasta: 'str', *, path_to_samtools: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs methylpy (extra `dnakit[bio-methylpy]`) and the external program `samtools` on PATH (or passed through the matching `path_to_...` argument). Without them the call raises `BackendUnavailableError`.
- **Input:** Same as `methylpy_pe_merge_sorted_multimap`: `current_library`, `files`, `prefix`, an existing `reference_fasta`, and optional `path_to_samtools`.
- **Output:** A dictionary with `status` (`"completed"`), `backend` (`"methylpy"`), the methylpy `version`, `function` (the methylpy function that ran), `outputs` (an empty list: the function declares no output path) and `arguments` (the sorted names of the arguments that were forwarded).
- **Example code:**

```python
from pathlib import Path

from dnakit.epigenome import methylpy_pe_merge_sorted_multimap_max_mapq
from dnakit.exceptions import BackendUnavailableError

Path("lib1_sorted_0").touch()
try:
    methylpy_pe_merge_sorted_multimap_max_mapq(
        "lib1", ["lib1_sorted_0"], "s1_lib1", "ref.fa", path_to_samtools="no_such_dir", verbose=False
    )
except BackendUnavailableError as error:
    print(type(error).__name__, "-", error.message)
```

- **Example output:**

```text
BackendUnavailableError - samtools was not found at no_such_dir/samtools.
```

Usage with the dependency installed:

```python
from dnakit.epigenome import methylpy_pe_merge_sorted_multimap_max_mapq

result = methylpy_pe_merge_sorted_multimap_max_mapq(
    "lib1", ["s1_lib1_sorted_0"], "s1_lib1", "genome.fa", verbose=False
)
# s1_lib1_processed_reads.bam is written
```

## 10) `dnakit.epigenome.methylpy_pe_run_alignment()`

- **Purpose:** Align the paired-end reads of one library on both converted references and merge the two strands.
- **Method:** Check the files and the programs, then call `methylpy.call_mc_pe.run_alignment_pe`. Upstream aligns the reads to the C-to-T converted reference and to the G-to-A converted reference with the chosen aligner (`bowtie2` by default, or `bowtie`), merges the two strands and keeps the reads that map in one place.
- **API:** `dnakit.epigenome.methylpy_pe_run_alignment(current_library: 'str', library_read1_files: 'Sequence[str]', library_read2_files: 'Sequence[str]', sample: 'str', forward_reference: 'str', reverse_reference: 'str', reference_fasta: 'str', *, path_to_output: 'str | None' = None, path_to_samtools: 'str | None' = None, aligner: 'str | None' = None, aligner_options: 'Sequence[str] | None' = None, path_to_aligner: 'str | None' = None, merge_by_max_mapq: 'bool | None' = None, min_mapq: 'int | None' = None, num_procs: 'int | None' = None, keep_temp_files: 'bool | None' = None, sort_mem: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs methylpy (extra `dnakit[bio-methylpy]`) and the external programs `samtools` and `bowtie2` (or `bowtie` with `aligner="bowtie"`) on PATH (or passed through the matching `path_to_...` argument). Without them the call raises `BackendUnavailableError`.
- **Input:** `current_library`, the library ID; `library_read1_files` and `library_read2_files`, lists of existing FASTQ files; `sample`; `forward_reference` and `reverse_reference`, the converted FASTA files built by `methylpy_se_build_ref` (they must exist); `reference_fasta`. Optional: `path_to_output`, `aligner`, `aligner_options`, `path_to_aligner`, `path_to_samtools`, `merge_by_max_mapq`, `min_mapq`, `num_procs`, `keep_temp_files`, `sort_mem`.
- **Output:** A dictionary with `status` (`"completed"`), `backend` (`"methylpy"`), the methylpy `version`, `function` (the methylpy function that ran), `outputs` (an empty list: the function declares no output path) and `arguments` (the sorted names of the arguments that were forwarded).
- **Example code:**

```python
from pathlib import Path

from dnakit.epigenome import methylpy_pe_run_alignment
from dnakit.exceptions import BackendUnavailableError

Path("ref_f.fasta").touch()
Path("ref_r.fasta").touch()
try:
    methylpy_pe_run_alignment(
        "lib1", ["reads.fq"], ["reads.fq"], "s1", "ref_f.fasta", "ref_r.fasta", "ref.fa",
        path_to_samtools="no_such_dir", verbose=False,
    )
except BackendUnavailableError as error:
    print(type(error).__name__, "-", error.message)
```

- **Example output:**

```text
BackendUnavailableError - samtools was not found at no_such_dir/samtools.
```

Usage with the dependency installed:

```python
from dnakit.epigenome import methylpy_pe_run_alignment

result = methylpy_pe_run_alignment(
    "lib1", ["lib1_R1.fastq.gz"], ["lib1_R2.fastq.gz"], "s1",
    "genome_f.fasta", "genome_r.fasta", "genome.fa", num_procs=4, verbose=False,
)
```

## 11) `dnakit.epigenome.methylpy_pe_run_mapping()`

- **Purpose:** Trim, align and merge one library of paired-end bisulfite reads.
- **Method:** Check the files and the programs, then call `methylpy.call_mc_pe.run_mapping_pe`. Upstream aligns the reads to the C-to-T converted reference and to the G-to-A converted reference with the chosen aligner (`bowtie2` by default, or `bowtie`), after optional adapter trimming with cutadapt (`trim_reads`, on by default), merges the two strands and keeps the reads that map in one place.
- **API:** `dnakit.epigenome.methylpy_pe_run_mapping(current_library: 'str', library_read1_files: 'Sequence[str]', library_read2_files: 'Sequence[str]', sample: 'str', forward_reference: 'str', reverse_reference: 'str', reference_fasta: 'str', *, path_to_output: 'str | None' = None, path_to_samtools: 'str | None' = None, path_to_aligner: 'str | None' = None, aligner: 'str | None' = None, aligner_options: 'Sequence[str] | None' = None, merge_by_max_mapq: 'bool | None' = None, min_mapq: 'int | None' = None, num_procs: 'int | None' = None, trim_reads: 'bool | None' = None, path_to_cutadapt: 'str | None' = None, adapter_seq_read1: 'str | None' = None, adapter_seq_read2: 'str | None' = None, min_qual_score: 'int | None' = None, min_read_len: 'int | None' = None, keep_temp_files: 'bool | None' = None, sort_mem: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs methylpy (extra `dnakit[bio-methylpy]`) and the external programs `samtools`, `bowtie2` (or `bowtie`) and `cutadapt` (not needed with `trim_reads=False`) on PATH (or passed through the matching `path_to_...` argument). Without them the call raises `BackendUnavailableError`.
- **Input:** Like `methylpy_pe_run_alignment` (the reference files must exist) plus the trimming options `trim_reads`, `path_to_cutadapt`, `adapter_seq_read1`, `adapter_seq_read2`, `min_qual_score` and `min_read_len`.
- **Output:** A dictionary with `status` (`"completed"`), `backend` (`"methylpy"`), the methylpy `version`, `function` (the methylpy function that ran), `outputs` (an empty list: the function declares no output path) and `arguments` (the sorted names of the arguments that were forwarded).
- **Example code:**

```python
from pathlib import Path

from dnakit.epigenome import methylpy_pe_run_mapping
from dnakit.exceptions import BackendUnavailableError

Path("ref_f.fasta").touch()
Path("ref_r.fasta").touch()
try:
    methylpy_pe_run_mapping(
        "lib1", ["reads.fq"], ["reads.fq"], "s1", "ref_f.fasta", "ref_r.fasta", "ref.fa",
        path_to_samtools="no_such_dir", verbose=False,
    )
except BackendUnavailableError as error:
    print(type(error).__name__, "-", error.message)
```

- **Example output:**

```text
BackendUnavailableError - samtools was not found at no_such_dir/samtools.
```

Usage with the dependency installed:

```python
from dnakit.epigenome import methylpy_pe_run_mapping

result = methylpy_pe_run_mapping(
    "lib1", ["lib1_R1.fastq.gz"], ["lib1_R2.fastq.gz"], "s1",
    "genome_f.fasta", "genome_r.fasta", "genome.fa", num_procs=4, verbose=False,
)
```

## 12) `dnakit.epigenome.methylpy_pe_run_methylation_pipeline()`

- **Purpose:** Run the whole paired-end bisulfite pipeline from FASTQ files to allc tables.
- **Method:** Check the inputs and the programs, then call `methylpy.call_mc_pe.run_methylation_pipeline_pe`. Upstream trims adapters, aligns each library to both converted references, merges the libraries, optionally removes clonal reads with Picard, and calls methylated sites per chromosome. This is the entry point that chains the other paired-end functions of this page.
- **API:** `dnakit.epigenome.methylpy_pe_run_methylation_pipeline(read1_files: 'Sequence[str]', read2_files: 'Sequence[str]', sample: 'str', forward_reference: 'str', reverse_reference: 'str', reference_fasta: 'str', *, libraries: 'Sequence[str] | None' = None, unmethylated_control: 'str | None' = None, path_to_output: 'str | None' = None, sig_cutoff: 'float | None' = None, num_procs: 'int | None' = None, sort_mem: 'str | None' = None, num_upstr_bases: 'int | None' = None, num_downstr_bases: 'int | None' = None, generate_allc_file: 'bool | None' = None, compress_output: 'bool | None' = None, bgzip: 'bool | None' = None, binom_test: 'bool | None' = None, min_cov: 'int | None' = None, trim_reads: 'bool | None' = None, aligner: 'str | None' = None, aligner_options: 'Sequence[str] | None' = None, merge_by_max_mapq: 'bool | None' = None, min_mapq: 'int | None' = None, remove_clonal: 'bool | None' = None, remove_chr_prefix: 'bool | None' = None, pbat: 'bool | None' = None, path_to_samtools: 'str | None' = None, path_to_aligner: 'str | None' = None, path_to_cutadapt: 'str | None' = None, path_to_picard: 'str | None' = None, path_to_bgzip: 'str | None' = None, path_to_tabix: 'str | None' = None, adapter_seq_read1: 'str | None' = None, adapter_seq_read2: 'str | None' = None, min_qual_score: 'int | None' = None, min_read_len: 'int | None' = None, keep_temp_files: 'bool | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs methylpy (extra `dnakit[bio-methylpy]`) and the external programs `samtools`, `bowtie2` (or `bowtie`), `cutadapt` (not needed with `trim_reads=False`), `bgzip` and `tabix` (only with `bgzip=True`); Picard and Java are needed for `remove_clonal=True` and are not checked first on PATH (or passed through the matching `path_to_...` argument). Without them the call raises `BackendUnavailableError`.
- **Input:** `read1_files` and `read2_files`, lists of existing FASTQ files in matching order; `sample`; `forward_reference` and `reverse_reference` (existing converted FASTA files); `reference_fasta`. Optional: `libraries`, `unmethylated_control`, `path_to_output`, `sig_cutoff`, `num_procs`, `min_cov`, `trim_reads`, `aligner`, `aligner_options`, `merge_by_max_mapq`, `min_mapq`, `remove_clonal`, `remove_chr_prefix`, `pbat`, `generate_allc_file`, `compress_output`, `bgzip`, `binom_test`, the `path_to_...` arguments, the adapter and quality options, and `keep_temp_files`.
- **Output:** A dictionary with `status` (`"completed"`), `backend` (`"methylpy"`), the methylpy `version`, `function` (the methylpy function that ran), `outputs` (an empty list: the function declares no output path) and `arguments` (the sorted names of the arguments that were forwarded).
- **Example code:**

```python
from pathlib import Path

from dnakit.epigenome import methylpy_pe_run_methylation_pipeline
from dnakit.exceptions import BackendUnavailableError

Path("ref_f.fasta").touch()
Path("ref_r.fasta").touch()
try:
    methylpy_pe_run_methylation_pipeline(
        ["reads.fq"], ["reads.fq"], "s1", "ref_f.fasta", "ref_r.fasta", "ref.fa",
        path_to_samtools="no_such_dir", verbose=False,
    )
except BackendUnavailableError as error:
    print(type(error).__name__, "-", error.message)
```

- **Example output:**

```text
BackendUnavailableError - samtools was not found at no_such_dir/samtools.
```

Usage with the dependency installed:

```python
from dnakit.epigenome import methylpy_pe_run_methylation_pipeline

result = methylpy_pe_run_methylation_pipeline(
    ["lib1_R1.fastq.gz"], ["lib1_R2.fastq.gz"], "s1",
    "genome_f.fasta", "genome_r.fasta", "genome.fa",
    unmethylated_control="chrC:", num_procs=8, verbose=False,
)
```

## 13) `dnakit.epigenome.methylpy_quality_trim_pe()`

- **Purpose:** Trim adapters and low-quality ends of paired-end bisulfite reads with cutadapt.
- **Method:** Check the inputs and that the outputs do not exist, check `cutadapt`, then call `methylpy.call_mc_pe.quality_trim_pe` with one-element lists. Upstream mishandles plain string paths (it assigns the list conversion to the wrong names), which is why DNAKit always passes lists. Upstream builds a shell command line from the file names, so avoid unusual characters in paths. The default adapters are the Illumina sequences `AGATCGGAAGAGCACACGTCTGAAC` (read 1) and `AGATCGGAAGAGCGTCGTGTAGGGA` (read 2); the defaults for `min_qual_score` and `min_read_len` are 10 and 30. In the example both outputs come from the same input file, so only read 1 loses its adapter (the file carries the read 1 adapter, not the read 2 adapter).
- **API:** `dnakit.epigenome.methylpy_quality_trim_pe(inputf_read1: 'str', outputf_read1: 'str', inputf_read2: 'str', outputf_read2: 'str', *, quality_base: 'int | None' = None, min_qual_score: 'int | None' = None, min_read_len: 'int | None' = None, adapter_seq_read1: 'str | None' = None, adapter_seq_read2: 'str | None' = None, num_procs: 'int | None' = None, input_format: 'str | None' = None, error_rate: 'float | None' = None, max_adapter_removal: 'int | None' = None, overlap_length: 'int | None' = None, zero_cap: 'bool | None' = None, path_to_cutadapt: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs methylpy (extra `dnakit[bio-methylpy]`) and the external program `cutadapt` (extra `dnakit[bio-cutadapt]`). Upstream finds it through `path_to_cutadapt` (a folder) or PATH; the example passes the folder of the running Python, where pip installs it.
- **Input:** `inputf_read1`, `inputf_read2`: existing FASTQ files; `outputf_read1`, `outputf_read2`: new paths that must not exist. Optional: `quality_base`, `min_qual_score`, `min_read_len`, `adapter_seq_read1`, `adapter_seq_read2`, `num_procs`, `input_format`, `error_rate`, `max_adapter_removal`, `overlap_length`, `zero_cap`, `path_to_cutadapt`.
- **Output:** The result dictionary described above, with `outputs` holding the two trimmed FASTQ paths.
- **Example code:**

```python
import sys
import tempfile
from pathlib import Path

from dnakit.epigenome import methylpy_quality_trim_pe

folder = Path(tempfile.mkdtemp())
result = methylpy_quality_trim_pe(
    "reads.fq", str(folder / "r1.fq"), "reads.fq", str(folder / "r2.fq"),
    min_read_len=10,
    path_to_cutadapt=str(Path(sys.executable).parent),
    verbose=False,
)
print(result["status"], len(result["outputs"]))
for name in ("r1.fq", "r2.fq"):
    lines = (folder / name).read_text().splitlines()
    print(name, [lines[i + 1] for i in range(0, len(lines), 4)])
```

- **Example output:**

```text
completed 2
r1.fq ['ACGTACGTTTGGA', 'TTGGTTGGTTAAGGCCAAGGTTTTAAGGCCTTAAGGAATT']
r2.fq ['ACGTACGTTTGGAAGATCGGAAGAGCACACGTCTGAAAA', 'TTGGTTGGTTAAGGCCAAGGTTTTAAGGCCTTAAGGAATT']
```

Usage with the dependency installed:

```python
from dnakit.epigenome import methylpy_quality_trim_pe

result = methylpy_quality_trim_pe(
    "lib1_R1.fastq", "lib1_R1_trimmed.fastq", "lib1_R2.fastq", "lib1_R2_trimmed.fastq",
    min_qual_score=20, min_read_len=30, verbose=False,
)
```

## 14) `dnakit.epigenome.methylpy_quality_trim_se()`

- **Purpose:** Trim adapters and low-quality ends of single-end bisulfite reads with cutadapt.
- **Method:** Check the inputs, that the outputs do not exist, the format and the adapter, find `cutadapt`, and run it once per file as an argument list with no shell (this replaces methylpy's `quality_trim`, which builds a shell command from interpolated strings). The command is `cutadapt --quiet [-z] [-f F] [-e E] [-n N] [-O O] [-m MIN_LEN] [-q MIN_QUAL] [--quality-base=B] -a ADAPTER -o OUT IN`. An option is added only when its value is truthy. Unlike the other trimming function, the defaults of `min_qual_score` and `min_read_len` are not methylpy's: cutadapt's own defaults apply (no quality trimming, no length filter) unless you set them.
- **API:** `dnakit.epigenome.methylpy_quality_trim_se(inputf: 'Sequence[str]', output: 'Sequence[str]', *, quality_base: 'int | None' = None, min_qual_score: 'int | None' = None, min_read_len: 'int | None' = None, adapter_seq: 'str' = 'AGATCGGAAGAGCACACGTCTG', input_format: 'str | None' = None, error_rate: 'float | None' = None, max_adapter_removal: 'int | None' = None, overlap_length: 'int | None' = None, zero_cap: 'bool' = False, path_to_cutadapt: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program `cutadapt` (extra `dnakit[bio-cutadapt]`) on PATH or in the folder of the running Python, or passed as the folder `path_to_cutadapt`. It does not need methylpy itself, except to report its version.
- **Input:** `inputf`, a list of existing FASTQ or FASTA files; `output`, a list of new paths of the same length that must not exist. `adapter_seq` (default `AGATCGGAAGAGCACACGTCTG`, IUPAC letters only), `input_format` (`fasta`, `fastq` or `sra-fastq`), `min_read_len`, `min_qual_score`, `quality_base`, `error_rate`, `max_adapter_removal`, `overlap_length`, `zero_cap`, `path_to_cutadapt`.
- **Output:** A dictionary with `status`, `backend`, `version`, `function`, `outputs` (the trimmed files that exist) and `log_tail` (the last 4000 characters of cutadapt's stderr for each file).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import methylpy_quality_trim_se

output = Path(tempfile.mkdtemp()) / "trimmed.fq"
result = methylpy_quality_trim_se(["reads.fq"], [output], min_read_len=10, verbose=False)
lines = output.read_text().splitlines()
print(result["status"], len(result["outputs"]), len(lines) // 4)
for index in range(0, len(lines), 4):
    print(lines[index], lines[index + 1])
```

- **Example output:**

```text
completed 1 2
@r1 ACGTACGTTTGGA
@r2 TTGGTTGGTTAAGGCCAAGGTTTTAAGGCCTTAAGGAATT
```

Usage with the dependency installed:

```python
from dnakit.epigenome import methylpy_quality_trim_se

result = methylpy_quality_trim_se(
    ["lib1.fastq", "lib2.fastq"], ["lib1_trimmed.fastq", "lib2_trimmed.fastq"],
    min_qual_score=20, min_read_len=30, verbose=False,
)
```

## 15) `dnakit.epigenome.methylpy_remove_clonal_bam()`

- **Purpose:** Remove PCR-duplicate (clonal) reads from a BAM file with Picard MarkDuplicates.
- **Method:** Check that the input exists and the output does not, require `path_to_picard`, check that `<path_to_picard>/picard.jar` exists and that `java` is available, then call `methylpy.call_mc_se.remove_clonal_bam`, which runs `java <java_options> -jar <path_to_picard>/picard.jar MarkDuplicates` with `ASSUME_SORTED=true`, `REMOVE_DUPLICATES=true` and lenient validation, writes the duplicate metrics file `metric` and reads the number of duplicates back from it (DNAKit does not return that number).
- **API:** `dnakit.epigenome.methylpy_remove_clonal_bam(input_bam: 'str', output_bam: 'str', metric: 'str', *, is_pe: 'bool | None' = None, path_to_picard: 'str | None' = None, java_options: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs methylpy (extra `dnakit[bio-methylpy]`), `java` on PATH and a folder that contains `picard.jar`, passed as `path_to_picard`. Without the folder the call raises `ConfigurationError`; with a folder that has no `picard.jar` it raises `BackendUnavailableError`.
- **Input:** `input_bam`, an existing BAM; `output_bam`, a new path; `metric`, the duplicate metrics file; `path_to_picard`, the folder of `picard.jar` (required); optional `is_pe` (set for paired-end data) and `java_options` (upstream default `-Xmx20g`).
- **Output:** A dictionary with `status` (`"completed"`), `backend` (`"methylpy"`), the methylpy `version`, `function` (the methylpy function that ran), `outputs` (the declared output paths that exist after the call) and `arguments` (the sorted names of the arguments that were forwarded).
- **Example code:**

```python
from pathlib import Path

from dnakit.epigenome import methylpy_remove_clonal_bam
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

Path("sample.bam").touch()
for options in ({}, {"path_to_picard": "no_such_dir"}):
    try:
        methylpy_remove_clonal_bam("sample.bam", "dedup.bam", "metrics.txt", verbose=False, **options)
    except (ConfigurationError, BackendUnavailableError) as error:
        print(type(error).__name__, "-", error.message)
```

- **Example output:**

```text
ConfigurationError - path_to_picard must name the folder that contains picard.jar (methylpy reads path_to_picard + '/picard.jar').
BackendUnavailableError - picard.jar was not found at no_such_dir/picard.jar.
```

Usage with the dependency installed:

```python
from dnakit.epigenome import methylpy_remove_clonal_bam

result = methylpy_remove_clonal_bam(
    "s1_lib1.bam", "s1_lib1_no_clonal.bam", "s1_lib1_metrics.txt",
    path_to_picard="/opt/picard", is_pe=True, verbose=False,
)
print(result["outputs"])
```

## 16) `dnakit.epigenome.methylpy_run_rms_tests()`

- **Purpose:** Run the root-mean-square permutation tests of DMRfind on prepared per-chromosome files.
- **Method:** Check that the bundled `run_rms_tests.out` binary and the input files exist, then call `methylpy.DMRfind.run_rms_tests`, which starts the binary on the files and writes the test results to `output`. `DMRfind` calls this step internally; use it directly only when you already hold methylpy's intermediate files.
- **API:** `dnakit.epigenome.methylpy_run_rms_tests(files: 'Sequence[str]', output: 'str', samples: 'Sequence[str]', *, min_cov: 'int | None' = None, num_sims: 'int | str | None' = None, num_sig_tests: 'int | str | None' = None, seed: 'int | None' = None, keep_temp_files: 'bool | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs methylpy (extra `dnakit[bio-methylpy]`); the bundled `run_rms_tests.out` binary needs the GSL shared libraries `libgsl.so.0` and `libgslcblas.so.0` (`apt install libgsl-dev` or `conda install -c conda-forge gsl`).
- **Input:** `files`, a list of existing input files; `output`, the result path; `samples`, the sample names. Optional: `min_cov`, `num_sims` and `num_sig_tests` (integers or strings), `seed`, `keep_temp_files`.
- **Output:** The result dictionary described above (`outputs` is an empty list). Missing input files raise `ConfigurationError`.
- **Example code:**

```python
from dnakit.epigenome import methylpy_run_rms_tests
from dnakit.exceptions import ConfigurationError

try:
    methylpy_run_rms_tests(["missing_input.tsv"], "rms_out.tsv", ["a", "b"], verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, "-", error.message.split(":")[0])
```

- **Example output:**

```text
ConfigurationError - files must be an existing file
```

Usage with the dependency installed:

```python
from dnakit.epigenome import methylpy_run_rms_tests

result = methylpy_run_rms_tests(
    ["dmr_rms_input_chr1.tsv"], "dmr_rms_out.tsv", ["a", "b"], num_sims=1000, verbose=False
)
```

## 17) `dnakit.epigenome.methylpy_se_build_ref()`

- **Purpose:** Build the C-to-T and G-to-A converted reference FASTA files and their aligner indexes.
- **Method:** Check the inputs and the index builder (`bowtie2-build`, or `bowtie-build` when `aligner="bowtie"`), then call `methylpy.call_mc_se.build_ref`. Upstream reads the FASTA files in blocks of `buffsize` bytes and writes `<output>_f.fasta` (every C replaced by T) and `<output>_r.fasta` (every G replaced by A), then builds an index for each. It serves the paired-end functions too.
- **API:** `dnakit.epigenome.methylpy_se_build_ref(input_files: 'Sequence[str]', output: 'str', *, buffsize: 'int | None' = None, aligner: 'str | None' = None, path_to_aligner: 'str | None' = None, num_procs: 'int | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs methylpy (extra `dnakit[bio-methylpy]`) and the external program `bowtie2-build` (from bowtie2; `bowtie-build` with `aligner="bowtie"`) on PATH (or passed through the matching `path_to_...` argument). Without them the call raises `BackendUnavailableError`.
- **Input:** `input_files`, a non-empty list of existing FASTA files; `output`, the prefix of the two converted references. Optional: `buffsize`, `aligner` (`bowtie2` or `bowtie`), `path_to_aligner` (folder), `num_procs`.
- **Output:** A dictionary with `status` (`"completed"`), `backend` (`"methylpy"`), the methylpy `version`, `function` (the methylpy function that ran), `outputs` (an empty list: the function declares no output path) and `arguments` (the sorted names of the arguments that were forwarded).
- **Example code:**

```python
from pathlib import Path

from dnakit.epigenome import methylpy_se_build_ref
from dnakit.exceptions import BackendUnavailableError

try:
    methylpy_se_build_ref(["ref.fa"], "genome", path_to_aligner="no_such_dir", verbose=False)
except BackendUnavailableError as error:
    print(type(error).__name__, "-", error.message)
```

- **Example output:**

```text
BackendUnavailableError - bowtie2-build was not found at no_such_dir/bowtie2-build.
```

Usage with the dependency installed:

```python
from dnakit.epigenome import methylpy_se_build_ref

result = methylpy_se_build_ref(["genome.fa"], "genome", num_procs=4, verbose=False)
# genome_f.fasta, genome_r.fasta and the bowtie2 indexes are written
```

## 18) `dnakit.epigenome.methylpy_se_call_methylated_sites()`

- **Purpose:** Call methylated cytosines from a sorted single-end BAM file.
- **Method:** Check that the BAM and the reference exist, check `samtools` (and `bgzip` and `tabix` when `bgzip=True`), then call `methylpy.call_mc_se.call_methylated_sites`. Upstream builds a pileup per chromosome, estimates the non-conversion rate from `unmethylated_control` and applies a binomial test with the `sig_cutoff` adjusted p-value, writing one `allc_<sample>_<chrom>.tsv` file per chromosome.
- **API:** `dnakit.epigenome.methylpy_se_call_methylated_sites(inputf: 'str', sample: 'str', reference_fasta: 'str', *, unmethylated_control: 'str | None' = None, sig_cutoff: 'float | None' = None, num_procs: 'int | None' = None, num_upstr_bases: 'int | None' = None, num_downstr_bases: 'int | None' = None, generate_mpileup_file: 'bool | None' = None, compress_output: 'bool | None' = None, bgzip: 'bool | None' = None, path_to_bgzip: 'str | None' = None, path_to_tabix: 'str | None' = None, min_mapq: 'int | None' = None, min_cov: 'int | None' = None, binom_test: 'bool | None' = None, path_to_samtools: 'str | None' = None, remove_chr_prefix: 'bool | None' = None, sort_mem: 'str | None' = None, add_snp_info: 'bool | None' = None, min_base_quality: 'int | None' = None, keep_temp_files: 'bool | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs methylpy (extra `dnakit[bio-methylpy]`) and the external programs `samtools` (plus `bgzip` and `tabix` when `bgzip=True`) on PATH (or passed through the matching `path_to_...` argument). Without them the call raises `BackendUnavailableError`.
- **Input:** `inputf`, a sorted single-end BAM; `sample`, the name used in the output files; `reference_fasta`, an existing FASTA. Optional: `unmethylated_control` (a `chrom:start-end` region or a number), `sig_cutoff`, `num_procs`, `num_upstr_bases`, `num_downstr_bases`, `min_mapq`, `min_cov`, `min_base_quality`, `binom_test`, `generate_mpileup_file`, `compress_output`, `bgzip`, `remove_chr_prefix`, `add_snp_info`, `sort_mem`, `path_to_samtools`, `path_to_bgzip`, `path_to_tabix`, `keep_temp_files`.
- **Output:** A dictionary with `status` (`"completed"`), `backend` (`"methylpy"`), the methylpy `version`, `function` (the methylpy function that ran), `outputs` (an empty list: the function declares no output path) and `arguments` (the sorted names of the arguments that were forwarded).
- **Example code:**

```python
from pathlib import Path

from dnakit.epigenome import methylpy_se_call_methylated_sites
from dnakit.exceptions import BackendUnavailableError

Path("reads_se.bam").touch()
try:
    methylpy_se_call_methylated_sites(
        "reads_se.bam", "s1", "ref.fa", path_to_samtools="no_such_dir", verbose=False
    )
except BackendUnavailableError as error:
    print(type(error).__name__, "-", error.message)
```

- **Example output:**

```text
BackendUnavailableError - samtools was not found at no_such_dir/samtools.
```

Usage with the dependency installed:

```python
from dnakit.epigenome import methylpy_se_call_methylated_sites

result = methylpy_se_call_methylated_sites(
    "s1_processed_reads_no_clonal.bam", "s1", "genome.fa",
    unmethylated_control="chrC:", num_procs=4, verbose=False,
)
print(result["status"])  # allc_s1_<chrom>.tsv.gz files are written
```

## 19) `dnakit.epigenome.methylpy_se_call_methylated_sites_with_snp_info()`

- **Purpose:** Call methylated cytosines of a single-end BAM and write the SNP information at each covered site.
- **Method:** Check that the BAM and the reference exist, check `samtools` (and `bgzip` and `tabix` when `bgzip=True`), then call `methylpy.call_mc_se.call_methylated_sites_with_SNP_info`. It works like `methylpy_se_call_methylated_sites`, but the allc rows have nine columns: the usual seven (the last is fixed to 1, with no binomial test) plus, for each position of the context, the base that disagrees with the expected conversion and its coverage.
- **API:** `dnakit.epigenome.methylpy_se_call_methylated_sites_with_snp_info(inputf: 'str', sample: 'str', reference_fasta: 'str', *, unmethylated_control: 'float | None' = None, sig_cutoff: 'float | None' = None, num_procs: 'int | None' = None, num_upstr_bases: 'int | None' = None, num_downstr_bases: 'int | None' = None, generate_mpileup_file: 'bool | None' = None, compress_output: 'bool | None' = None, bgzip: 'bool | None' = None, path_to_bgzip: 'str | None' = None, path_to_tabix: 'str | None' = None, min_mapq: 'int | None' = None, min_cov: 'int | None' = None, binom_test: 'bool | None' = None, path_to_samtools: 'str | None' = None, remove_chr_prefix: 'bool | None' = None, sort_mem: 'str | None' = None, add_snp_info: 'bool | None' = None, min_base_quality: 'int | None' = None, keep_temp_files: 'bool | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs methylpy (extra `dnakit[bio-methylpy]`) and the external programs `samtools` (plus `bgzip` and `tabix` when `bgzip=True`) on PATH (or passed through the matching `path_to_...` argument). Without them the call raises `BackendUnavailableError`.
- **Input:** `inputf`, a sorted single-end BAM; `sample`, the name used in the output files; `reference_fasta`, an existing FASTA. Optional: `unmethylated_control` (typed as a number in this signature), `sig_cutoff`, `num_procs`, `num_upstr_bases`, `num_downstr_bases`, `min_mapq`, `min_cov`, `min_base_quality`, `binom_test`, `generate_mpileup_file`, `compress_output`, `bgzip`, `remove_chr_prefix`, `add_snp_info`, `sort_mem`, `path_to_samtools`, `path_to_bgzip`, `path_to_tabix`, `keep_temp_files`.
- **Output:** A dictionary with `status` (`"completed"`), `backend` (`"methylpy"`), the methylpy `version`, `function` (the methylpy function that ran), `outputs` (an empty list: the function declares no output path) and `arguments` (the sorted names of the arguments that were forwarded).
- **Example code:**

```python
from pathlib import Path

from dnakit.epigenome import methylpy_se_call_methylated_sites_with_snp_info
from dnakit.exceptions import BackendUnavailableError

Path("reads_se.bam").touch()
try:
    methylpy_se_call_methylated_sites_with_snp_info(
        "reads_se.bam", "s1", "ref.fa", path_to_samtools="no_such_dir", verbose=False
    )
except BackendUnavailableError as error:
    print(type(error).__name__, "-", error.message)
```

- **Example output:**

```text
BackendUnavailableError - samtools was not found at no_such_dir/samtools.
```

Usage with the dependency installed:

```python
from dnakit.epigenome import methylpy_se_call_methylated_sites_with_snp_info

result = methylpy_se_call_methylated_sites_with_snp_info(
    "s1_processed_reads_no_clonal.bam", "s1", "genome.fa", num_procs=4, verbose=False
)
```

## 20) `dnakit.epigenome.methylpy_se_find_multi_mappers()`

- **Purpose:** Split single-end alignments into uniquely scored reads, restoring the original read sequence.
- **Method:** Call `methylpy.call_mc_se.find_multi_mappers`. It reads the SAM file written by the aligner, skips header lines, unmapped reads, reads below `max(3, min_mapq)` and supplementary alignments, and decodes the bisulfite-converted sequence using the positions encoded after the `!` in the read name (reverse-strand reads are reversed first). The kept lines are written round-robin to `<output>_sorted_<i>` for `i < num_procs`. The input SAM file is deleted unless `keep_temp_files=True`.
- **API:** `dnakit.epigenome.methylpy_se_find_multi_mappers(inputf: 'str', output: 'str', *, num_procs: 'int | None' = None, min_mapq: 'int | None' = None, keep_temp_files: 'bool | None' = None, append: 'bool | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs methylpy (extra `dnakit[bio-methylpy]`) only; it runs no external program. The input must be a SAM file whose read names end in `!<encoded positions>`, as written by the methylpy alignment steps.
- **Input:** `inputf`, an existing SAM file; `output`, the prefix of the result files. Optional: `num_procs` (number of output files), `min_mapq`, `keep_temp_files`, `append`.
- **Output:** A dictionary with `status` (`"completed"`), `backend` (`"methylpy"`), the methylpy `version`, `function` (the methylpy function that ran), `outputs` (always an empty list here, because the function declares no single output path: methylpy writes files named from the prefix) and `arguments` (the sorted names of the arguments that were forwarded).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import methylpy_se_find_multi_mappers

folder = Path(tempfile.mkdtemp())
sam = folder / "reads.sam"
sam.write_text(
    "@HD\tVN:1.0\n"
    "r1!\t0\tchr1\t10\t42\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\n"
    "r2!\t16\tchr1\t50\t42\t8M\t*\t0\t0\tTTGGTTAA\tIIIIIIII\n"
    "r3!\t0\tchr1\t80\t1\t8M\t*\t0\t0\tCCGGCCGG\tIIIIIIII\n"
    "r4!\t4\t*\t0\t0\t*\t*\t0\t0\tAAAAAAAA\tIIIIIIII\n"
)
result = methylpy_se_find_multi_mappers(str(sam), str(folder / "out"), verbose=False)
rows = (folder / "out_sorted_0").read_text().splitlines()
print(result["status"], result["outputs"], sam.exists())
print([row.split("\t")[:2] for row in rows])
```

- **Example output:**

```text
completed [] False
[['r1', '0'], ['r2', '16']]
```

Usage with the dependency installed:

```python
from dnakit.epigenome import methylpy_se_find_multi_mappers

result = methylpy_se_find_multi_mappers("s1_forward.sam", "s1_forward", num_procs=4, verbose=False)
# s1_forward_sorted_0 ... s1_forward_sorted_3 are written; s1_forward.sam is deleted
```

## 21) `dnakit.epigenome.methylpy_se_merge_bam_files()`

- **Purpose:** Merge several single-end BAM files into one BAM file with a corrected header.
- **Method:** Check that the inputs exist and `samtools` is available, then call `methylpy.call_mc_se.merge_bam_files`, which writes a temporary `header.sam` in the current directory (the header of the first file plus one `@RG` line per input, then removed) and runs `samtools merge -r -h header.sam`. Every input file name must contain `.bam`.
- **API:** `dnakit.epigenome.methylpy_se_merge_bam_files(input_files: 'Sequence[str]', output: 'str', *, path_to_samtools: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs methylpy (extra `dnakit[bio-methylpy]`) and the external program `samtools` on PATH (or passed through the matching `path_to_...` argument). Without them the call raises `BackendUnavailableError`.
- **Input:** `input_files`, a non-empty list of existing BAM files; `output`, the merged BAM; `path_to_samtools`.
- **Output:** A dictionary with `status` (`"completed"`), `backend` (`"methylpy"`), the methylpy `version`, `function` (the methylpy function that ran), `outputs` (the declared output paths that exist after the call) and `arguments` (the sorted names of the arguments that were forwarded).
- **Example code:**

```python
from pathlib import Path

from dnakit.epigenome import methylpy_se_merge_bam_files
from dnakit.exceptions import BackendUnavailableError

Path("lib1.bam").touch()
Path("lib2.bam").touch()
try:
    methylpy_se_merge_bam_files(
        ["lib1.bam", "lib2.bam"], "merged.bam", path_to_samtools="no_such_dir", verbose=False
    )
except BackendUnavailableError as error:
    print(type(error).__name__, "-", error.message)
```

- **Example output:**

```text
BackendUnavailableError - samtools was not found at no_such_dir/samtools.
```

Usage with the dependency installed:

```python
from dnakit.epigenome import methylpy_se_merge_bam_files

result = methylpy_se_merge_bam_files(
    ["s1_lib1_no_clonal.bam", "s1_lib2_no_clonal.bam"], "s1_merged.bam", verbose=False
)
print(result["outputs"])
```

## 22) `dnakit.epigenome.methylpy_se_merge_sorted_multimap()`

- **Purpose:** Merge the sorted single-end alignment files of one library and keep the uniquely mapping reads.
- **Method:** Check that the files and the reference exist and `samtools` is available, then call `methylpy.call_mc_se.merge_sorted_multimap`. Upstream indexes the reference with `samtools faidx` when `<reference_fasta>.fai` is missing, writes a SAM header from the index, merges the sorted files read by read and keeps a read only when it occurs once across the files, and converts the result to `<prefix>_processed_reads.bam`.
- **API:** `dnakit.epigenome.methylpy_se_merge_sorted_multimap(current_library: 'str', files: 'Sequence[str]', prefix: 'str', reference_fasta: 'str', *, path_to_samtools: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs methylpy (extra `dnakit[bio-methylpy]`) and the external program `samtools` on PATH (or passed through the matching `path_to_...` argument). Without them the call raises `BackendUnavailableError`.
- **Input:** `current_library`, the library ID; `files`, the `_sorted_<i>` files from `methylpy_se_find_multi_mappers`; `prefix`, the output prefix; `reference_fasta`, an existing FASTA. `path_to_samtools` is optional.
- **Output:** A dictionary with `status` (`"completed"`), `backend` (`"methylpy"`), the methylpy `version`, `function` (the methylpy function that ran), `outputs` (an empty list: the function declares no output path) and `arguments` (the sorted names of the arguments that were forwarded).
- **Example code:**

```python
from pathlib import Path

from dnakit.epigenome import methylpy_se_merge_sorted_multimap
from dnakit.exceptions import BackendUnavailableError

Path("lib1_sorted_0").touch()
try:
    methylpy_se_merge_sorted_multimap(
        "lib1", ["lib1_sorted_0"], "s1_lib1", "ref.fa", path_to_samtools="no_such_dir", verbose=False
    )
except BackendUnavailableError as error:
    print(type(error).__name__, "-", error.message)
```

- **Example output:**

```text
BackendUnavailableError - samtools was not found at no_such_dir/samtools.
```

Usage with the dependency installed:

```python
from dnakit.epigenome import methylpy_se_merge_sorted_multimap

result = methylpy_se_merge_sorted_multimap(
    "lib1", ["s1_lib1_sorted_0", "s1_lib1_sorted_1"], "s1_lib1", "genome.fa", verbose=False
)
# s1_lib1_processed_reads.bam is written
```

## 23) `dnakit.epigenome.methylpy_se_merge_sorted_multimap_max_mapq()`

- **Purpose:** Merge the sorted single-end alignment files of one library and keep the best-MAPQ alignment of each multiply mapped read.
- **Method:** Check that the files and the reference exist and `samtools` is available, then call `methylpy.call_mc_se.merge_sorted_multimap_max_mapq`. Upstream indexes the reference with `samtools faidx` when `<reference_fasta>.fai` is missing, writes a SAM header from the index, merges the sorted files read by read and keeps, for a read that occurs more than once, the alignment with the highest mapping quality, and converts the result to `<prefix>_processed_reads.bam`.
- **API:** `dnakit.epigenome.methylpy_se_merge_sorted_multimap_max_mapq(current_library: 'str', files: 'Sequence[str]', prefix: 'str', reference_fasta: 'str', *, path_to_samtools: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs methylpy (extra `dnakit[bio-methylpy]`) and the external program `samtools` on PATH (or passed through the matching `path_to_...` argument). Without them the call raises `BackendUnavailableError`.
- **Input:** Same as `methylpy_se_merge_sorted_multimap`: `current_library`, `files`, `prefix`, an existing `reference_fasta`, and optional `path_to_samtools`.
- **Output:** A dictionary with `status` (`"completed"`), `backend` (`"methylpy"`), the methylpy `version`, `function` (the methylpy function that ran), `outputs` (an empty list: the function declares no output path) and `arguments` (the sorted names of the arguments that were forwarded).
- **Example code:**

```python
from pathlib import Path

from dnakit.epigenome import methylpy_se_merge_sorted_multimap_max_mapq
from dnakit.exceptions import BackendUnavailableError

Path("lib1_sorted_0").touch()
try:
    methylpy_se_merge_sorted_multimap_max_mapq(
        "lib1", ["lib1_sorted_0"], "s1_lib1", "ref.fa", path_to_samtools="no_such_dir", verbose=False
    )
except BackendUnavailableError as error:
    print(type(error).__name__, "-", error.message)
```

- **Example output:**

```text
BackendUnavailableError - samtools was not found at no_such_dir/samtools.
```

Usage with the dependency installed:

```python
from dnakit.epigenome import methylpy_se_merge_sorted_multimap_max_mapq

result = methylpy_se_merge_sorted_multimap_max_mapq(
    "lib1", ["s1_lib1_sorted_0"], "s1_lib1", "genome.fa", verbose=False
)
# s1_lib1_processed_reads.bam is written
```

## 24) `dnakit.epigenome.methylpy_se_run_alignment()`

- **Purpose:** Align the single-end reads of one library on both converted references and merge the two strands.
- **Method:** Check the files and the programs, then call `methylpy.call_mc_se.run_alignment`. Upstream aligns the reads to the C-to-T converted reference and to the G-to-A converted reference with the chosen aligner (`bowtie2` by default, or `bowtie`), merges the two strands and keeps the reads that map in one place.
- **API:** `dnakit.epigenome.methylpy_se_run_alignment(current_library: 'str', library_read_files: 'Sequence[str]', sample: 'str', forward_reference: 'str', reverse_reference: 'str', reference_fasta: 'str', *, path_to_output: 'str | None' = None, path_to_samtools: 'str | None' = None, aligner: 'str | None' = None, path_to_aligner: 'str | None' = None, aligner_options: 'Sequence[str] | None' = None, merge_by_max_mapq: 'bool | None' = None, min_mapq: 'int | None' = None, num_procs: 'int | None' = None, keep_temp_files: 'bool | None' = None, sort_mem: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs methylpy (extra `dnakit[bio-methylpy]`) and the external programs `samtools` and `bowtie2` (or `bowtie` with `aligner="bowtie"`) on PATH (or passed through the matching `path_to_...` argument). Without them the call raises `BackendUnavailableError`.
- **Input:** `current_library`, the library ID; `library_read_files`, a list of existing FASTQ files; `sample`; `forward_reference` and `reverse_reference`, the converted FASTA files built by `methylpy_se_build_ref` (they must exist); `reference_fasta`. Optional: `path_to_output`, `aligner`, `aligner_options`, `path_to_aligner`, `path_to_samtools`, `merge_by_max_mapq`, `min_mapq`, `num_procs`, `keep_temp_files`, `sort_mem`.
- **Output:** A dictionary with `status` (`"completed"`), `backend` (`"methylpy"`), the methylpy `version`, `function` (the methylpy function that ran), `outputs` (an empty list: the function declares no output path) and `arguments` (the sorted names of the arguments that were forwarded).
- **Example code:**

```python
from pathlib import Path

from dnakit.epigenome import methylpy_se_run_alignment
from dnakit.exceptions import BackendUnavailableError

Path("ref_f.fasta").touch()
Path("ref_r.fasta").touch()
try:
    methylpy_se_run_alignment(
        "lib1", ["reads.fq"], "s1", "ref_f.fasta", "ref_r.fasta", "ref.fa",
        path_to_samtools="no_such_dir", verbose=False,
    )
except BackendUnavailableError as error:
    print(type(error).__name__, "-", error.message)
```

- **Example output:**

```text
BackendUnavailableError - samtools was not found at no_such_dir/samtools.
```

Usage with the dependency installed:

```python
from dnakit.epigenome import methylpy_se_run_alignment

result = methylpy_se_run_alignment(
    "lib1", ["lib1.fastq.gz"], "s1", "genome_f.fasta", "genome_r.fasta", "genome.fa",
    num_procs=4, verbose=False,
)
```

## 25) `dnakit.epigenome.methylpy_se_run_mapping()`

- **Purpose:** Trim, align and merge one library of single-end bisulfite reads.
- **Method:** Check the files and the programs, then call `methylpy.call_mc_se.run_mapping`. Upstream aligns the reads to the C-to-T converted reference and to the G-to-A converted reference with the chosen aligner (`bowtie2` by default, or `bowtie`), after optional adapter trimming with cutadapt (`trim_reads`, on by default), merges the two strands and keeps the reads that map in one place.
- **API:** `dnakit.epigenome.methylpy_se_run_mapping(current_library: 'str', library_files: 'Sequence[str]', sample: 'str', forward_reference: 'str', reverse_reference: 'str', reference_fasta: 'str', *, path_to_output: 'str | None' = None, path_to_samtools: 'str | None' = None, path_to_aligner: 'str | None' = None, aligner: 'str | None' = None, aligner_options: 'Sequence[str] | None' = None, merge_by_max_mapq: 'bool | None' = None, min_mapq: 'int | None' = None, pbat: 'bool | None' = None, num_procs: 'int | None' = None, trim_reads: 'bool | None' = None, path_to_cutadapt: 'str | None' = None, adapter_seq: 'str | None' = None, min_qual_score: 'int | None' = None, min_read_len: 'int | None' = None, keep_temp_files: 'bool | None' = None, sort_mem: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs methylpy (extra `dnakit[bio-methylpy]`) and the external programs `samtools`, `bowtie2` (or `bowtie`) and `cutadapt` (not needed with `trim_reads=False`) on PATH (or passed through the matching `path_to_...` argument). Without them the call raises `BackendUnavailableError`.
- **Input:** Like `methylpy_se_run_alignment` (with `library_files` instead of `library_read_files`; the reference files must exist) plus `pbat`, `trim_reads`, `path_to_cutadapt`, `adapter_seq`, `min_qual_score` and `min_read_len`.
- **Output:** A dictionary with `status` (`"completed"`), `backend` (`"methylpy"`), the methylpy `version`, `function` (the methylpy function that ran), `outputs` (an empty list: the function declares no output path) and `arguments` (the sorted names of the arguments that were forwarded).
- **Example code:**

```python
from pathlib import Path

from dnakit.epigenome import methylpy_se_run_mapping
from dnakit.exceptions import BackendUnavailableError

Path("ref_f.fasta").touch()
Path("ref_r.fasta").touch()
try:
    methylpy_se_run_mapping(
        "lib1", ["reads.fq"], "s1", "ref_f.fasta", "ref_r.fasta", "ref.fa",
        path_to_samtools="no_such_dir", verbose=False,
    )
except BackendUnavailableError as error:
    print(type(error).__name__, "-", error.message)
```

- **Example output:**

```text
BackendUnavailableError - samtools was not found at no_such_dir/samtools.
```

Usage with the dependency installed:

```python
from dnakit.epigenome import methylpy_se_run_mapping

result = methylpy_se_run_mapping(
    "lib1", ["lib1.fastq.gz"], "s1", "genome_f.fasta", "genome_r.fasta", "genome.fa",
    num_procs=4, verbose=False,
)
```

## 26) `dnakit.epigenome.methylpy_se_run_methylation_pipeline()`

- **Purpose:** Run the whole single-end bisulfite pipeline from FASTQ files to allc tables.
- **Method:** Check the inputs and the programs, then call `methylpy.call_mc_se.run_methylation_pipeline`. Upstream trims adapters, aligns each library to both converted references, merges the libraries, optionally removes clonal reads with Picard, and calls methylated sites per chromosome. This is the entry point that chains the other single-end functions of this page.
- **API:** `dnakit.epigenome.methylpy_se_run_methylation_pipeline(read_files: 'Sequence[str]', sample: 'str', forward_reference: 'str', reverse_reference: 'str', reference_fasta: 'str', *, libraries: 'Sequence[str] | None' = None, unmethylated_control: 'str | None' = None, path_to_output: 'str | None' = None, sig_cutoff: 'float | None' = None, num_procs: 'int | None' = None, sort_mem: 'str | None' = None, num_upstr_bases: 'int | None' = None, num_downstr_bases: 'int | None' = None, generate_allc_file: 'bool | None' = None, compress_output: 'bool | None' = None, bgzip: 'bool | None' = None, binom_test: 'bool | None' = None, min_cov: 'int | None' = None, trim_reads: 'bool | None' = None, aligner: 'str | None' = None, aligner_options: 'Sequence[str] | None' = None, merge_by_max_mapq: 'bool | None' = None, min_mapq: 'int | None' = None, remove_clonal: 'bool | None' = None, remove_chr_prefix: 'bool | None' = None, pbat: 'bool | None' = None, path_to_samtools: 'str | None' = None, path_to_aligner: 'str | None' = None, path_to_cutadapt: 'str | None' = None, path_to_picard: 'str | None' = None, path_to_bgzip: 'str | None' = None, path_to_tabix: 'str | None' = None, adapter_seq: 'str | None' = None, min_qual_score: 'int | None' = None, min_read_len: 'int | None' = None, keep_temp_files: 'bool | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs methylpy (extra `dnakit[bio-methylpy]`) and the external programs `samtools`, `bowtie2` (or `bowtie`), `cutadapt` (not needed with `trim_reads=False`), `bgzip` and `tabix` (only with `bgzip=True`); Picard and Java are needed for `remove_clonal=True` and are not checked first on PATH (or passed through the matching `path_to_...` argument). Without them the call raises `BackendUnavailableError`.
- **Input:** `read_files`, a list of existing FASTQ files; `sample`; `forward_reference` and `reverse_reference` (existing converted FASTA files); `reference_fasta`. Optional: `libraries`, `unmethylated_control`, `path_to_output`, `sig_cutoff`, `num_procs`, `min_cov`, `trim_reads`, `aligner`, `aligner_options`, `merge_by_max_mapq`, `min_mapq`, `remove_clonal`, `remove_chr_prefix`, `pbat`, `generate_allc_file`, `compress_output`, `bgzip`, `binom_test`, the `path_to_...` arguments, `adapter_seq`, `min_qual_score`, `min_read_len` and `keep_temp_files`.
- **Output:** A dictionary with `status` (`"completed"`), `backend` (`"methylpy"`), the methylpy `version`, `function` (the methylpy function that ran), `outputs` (an empty list: the function declares no output path) and `arguments` (the sorted names of the arguments that were forwarded).
- **Example code:**

```python
from pathlib import Path

from dnakit.epigenome import methylpy_se_run_methylation_pipeline
from dnakit.exceptions import BackendUnavailableError

Path("ref_f.fasta").touch()
Path("ref_r.fasta").touch()
try:
    methylpy_se_run_methylation_pipeline(
        ["reads.fq"], "s1", "ref_f.fasta", "ref_r.fasta", "ref.fa",
        path_to_samtools="no_such_dir", verbose=False,
    )
except BackendUnavailableError as error:
    print(type(error).__name__, "-", error.message)
```

- **Example output:**

```text
BackendUnavailableError - samtools was not found at no_such_dir/samtools.
```

Usage with the dependency installed:

```python
from dnakit.epigenome import methylpy_se_run_methylation_pipeline

result = methylpy_se_run_methylation_pipeline(
    ["lib1.fastq.gz"], "s1", "genome_f.fasta", "genome_r.fasta", "genome.fa",
    unmethylated_control="chrC:", num_procs=8, verbose=False,
)
```

## 27) `dnakit.epigenome.methylpy_filter_allc_files()`

- **Purpose:** Filter allc tables by cytosine context, chromosome and read coverage.
- **Method:** Validate the arguments, count the input rows, then call `methylpy.utilities.filter_allc_files` sequentially (`num_procs` is fixed to 1) with its output captured. A row is kept when its chromosome is in `chroms` (if given), its context matches one of the `mc_type` patterns (`H` and `N` are wildcards) and `min_cov <= total reads (column 6) <= max_cov`. A bare string for `mc_type` is refused because methylpy would then match by substring and keep every cytosine. The SNP mismatch filters of methylpy are not offered because they do nothing on Python 3 in methylpy 1.4.7. Existing outputs are never overwritten and partial outputs are removed on failure.
- **API:** `dnakit.epigenome.methylpy_filter_allc_files(allc_files: 'Sequence[str | Path]', output_files: 'Sequence[str | Path]', *, mc_type: 'Sequence[str] | None' = None, chroms: 'Sequence[str] | None' = None, min_cov: 'int' = 0, max_cov: 'int | None' = None, compress_output: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs methylpy (extra `dnakit[bio-methylpy]`) only; it runs no external program.
- **Input:** `allc_files`, a list of 1 to 32 existing allc tables (plain or gzip); `output_files`, one new path per input (a `.gz` name needs `compress_output=True`; with it methylpy appends `.gz` when missing). `mc_type`, a list of 3-letter contexts such as `["CGN"]`; `chroms`, a list of names; `min_cov` (default 0) and `max_cov` (default no limit), both inclusive; `compress_output` (default False).
- **Output:** A dictionary with `status`, `operation`, `inputs`, `outputs` (for each file its `input`, `output`, `records_in` and `records_out`), `filters` (the applied `mc_type`, `chroms`, `min_cov`, `max_cov`), `compressed` and `provenance` (the methylpy version, the rule used and the omitted mismatch filters).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import methylpy_filter_allc_files

output = Path(tempfile.mkdtemp()) / "cg.tsv"
result = methylpy_filter_allc_files(
    ["allc.tsv"], [output], mc_type=["CGN"], min_cov=5, verbose=False
)
row = result["outputs"][0]
print(row["records_in"], row["records_out"], result["filters"]["mc_type"])
for line in output.read_text().splitlines():
    print(line)
```

- **Example output:**

```text
5 2 ['CGN']
chr1	100	+	CGG	10	12	1
chr2	50	+	CGT	0	8	0
```
