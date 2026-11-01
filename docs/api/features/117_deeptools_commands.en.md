# deepTools: coverage, comparison, summaries and read filtering

This group of functions in `dnakit.epigenome` is a set of typed adapters over [deepTools](https://github.com/deeptools/deepTools) (version 4). Twenty-one of them each run one deepTools command-line program (20 different programs, `computeMatrix` is used by two functions: `bamCoverage`, `bamCompare`, `bigwigCompare`, `bigwigAverage`, `computeMatrix`, `multiBamSummary`, `multiBigwigSummary`, `plotHeatmap`, `plotProfile`, `plotCorrelation`, `plotPCA`, `plotFingerprint`, `plotCoverage`, `plotEnrichment`, `estimateReadFiltering`, `estimateScaleFactor`, `bamPEFragmentSize`, `alignmentSieve`, `computeGCBias`, `correctGCBias`) with a fixed, typed list of options and without a shell; `bam_mapping_statistics()` instead imports deepTools and calls `openBam()` in-process. They build coverage and comparison tracks, summarise and plot samples, compute and plot region matrices, estimate scale factors and read filtering, and compute and correct GC bias. They do not expose every option of the tools, and they do not call peaks or interpret the results.

deepTools is distributed under the MIT licence (checked in the `LICENSE.txt` of the installed 4.0.0 package). DNAKit runs the programs as separate processes and copies none of the deepTools code.

Import the functions with `from dnakit.epigenome import ...`. deepTools is an optional dependency (`pip install "dnakit[bio-deeptools]"`) and is used only when a function is called. The command-line facades take the programs from PATH, or from `deeptools_path` (the full path of one executable); if the program is missing they raise `BackendUnavailableError` with an install hint. Every facade checks that the input files exist and that the output paths do not exist yet (`ConfigurationError` otherwise). A tool that exits with an error does not raise: the returned `status` is `failed` and `log_tail` holds the message. Unless stated otherwise, every facade returns the dictionary described in the Output bullet.

## Example data

Both example files are tiny. `reads.sam` is a coordinate-sorted SAM file with 24 read pairs (48 alignments, no sequence) on `chr1` (3000 bases) and `chr2` (2000 bases); the examples convert it to an indexed BAM file with pysam. `genes.bed` holds two regions. Examples that need bigWig files write them with pyBigWig, and the two GC-bias examples build a random 2bit genome and the reads themselves. Each example first puts the folder of the running Python on PATH so that the deepTools commands installed in the same environment are found; the line is not needed when they are already on PATH.

**Example file:** `reads.sam`

```text
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:3000
@SQ	SN:chr2	LN:2000
p1	99	chr1	100	60	50M	=	250	200	*	*
p1	147	chr1	250	60	50M	=	100	-200	*	*
p2	99	chr1	262	60	50M	=	412	200	*	*
p2	147	chr1	412	60	50M	=	262	-200	*	*
p3	99	chr1	425	60	50M	=	575	200	*	*
p3	147	chr1	575	60	50M	=	425	-200	*	*
p4	99	chr1	587	60	50M	=	737	200	*	*
p4	147	chr1	737	60	50M	=	587	-200	*	*
p5	99	chr1	750	60	50M	=	900	200	*	*
p5	147	chr1	900	60	50M	=	750	-200	*	*
p6	99	chr1	912	60	50M	=	1062	200	*	*
p6	147	chr1	1062	60	50M	=	912	-200	*	*
p7	99	chr1	1075	60	50M	=	1225	200	*	*
p7	147	chr1	1225	60	50M	=	1075	-200	*	*
p8	99	chr1	1237	60	50M	=	1387	200	*	*
p8	147	chr1	1387	60	50M	=	1237	-200	*	*
p9	99	chr1	1400	60	50M	=	1550	200	*	*
p9	147	chr1	1550	60	50M	=	1400	-200	*	*
p10	99	chr1	1562	60	50M	=	1712	200	*	*
p10	147	chr1	1712	60	50M	=	1562	-200	*	*
p11	99	chr1	1725	60	50M	=	1875	200	*	*
p11	147	chr1	1875	60	50M	=	1725	-200	*	*
p12	99	chr1	1887	60	50M	=	2037	200	*	*
p12	147	chr1	2037	60	50M	=	1887	-200	*	*
p13	99	chr1	2050	60	50M	=	2200	200	*	*
p13	147	chr1	2200	60	50M	=	2050	-200	*	*
p14	99	chr1	2212	60	50M	=	2362	200	*	*
p14	147	chr1	2362	60	50M	=	2212	-200	*	*
p15	99	chr1	2375	60	50M	=	2525	200	*	*
p15	147	chr1	2525	60	50M	=	2375	-200	*	*
p16	99	chr1	2537	60	50M	=	2687	200	*	*
p16	147	chr1	2687	60	50M	=	2537	-200	*	*
p17	99	chr2	100	60	50M	=	250	200	*	*
p17	147	chr2	250	60	50M	=	100	-200	*	*
p18	99	chr2	300	60	50M	=	450	200	*	*
p18	147	chr2	450	60	50M	=	300	-200	*	*
p19	99	chr2	500	60	50M	=	650	200	*	*
p19	147	chr2	650	60	50M	=	500	-200	*	*
p20	99	chr2	700	60	50M	=	850	200	*	*
p20	147	chr2	850	60	50M	=	700	-200	*	*
p21	99	chr2	900	60	50M	=	1050	200	*	*
p21	147	chr2	1050	60	50M	=	900	-200	*	*
p22	99	chr2	1100	60	50M	=	1250	200	*	*
p22	147	chr2	1250	60	50M	=	1100	-200	*	*
p23	99	chr2	1300	60	50M	=	1450	200	*	*
p23	147	chr2	1450	60	50M	=	1300	-200	*	*
p24	99	chr2	1500	60	50M	=	1650	200	*	*
p24	147	chr2	1650	60	50M	=	1500	-200	*	*
```

**Example file:** `genes.bed`

```text
chr1	1000	2000	g1	0	+
chr1	2200	2800	g2	0	-
```

## 1) `dnakit.epigenome.alignment_sieve()`

- **Purpose:** Filter alignments of a BAM file and optionally shift their ends, writing a new BAM or a BEDPE file.
- **Method:** Check the inputs, build a fixed argument list and run `alignmentSieve` without a shell in the folder of the output file (time limit 300 s). The arguments are `--bam`, `-o` and `-p`, plus `--shift` when `shift` is given and `--BED` when `bed_output` is true. No other filter option is exposed, so with the defaults the tool only rewrites the alignments.
- **API:** `dnakit.epigenome.alignment_sieve(bam: 'str', output: 'str', *, shift: 'Sequence[int] | None' = None, bed_output: 'bool' = False, number_of_processors: 'int' = 1, deeptools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `alignmentSieve` (installed with deepTools, extra `dnakit[bio-deeptools]`) on PATH, or its full path in `deeptools_path`. The example also uses pysam and pyBigWig, which deepTools installs.
- **Input:** `bam` must be an existing file and `output` a new path (it must not exist and its folder must exist). `shift` is `None` or exactly four integers, otherwise a `ConfigurationError` is raised. As in deepTools, the first two are the shifts of the left and right fragment ends for fragments whose read 1 is on the left, and the last two for fragments whose read 1 is on the right (positive moves an end to the right). The error message of the facade describes the four values as forward and reverse starts and ends, which does not match the deepTools help. `bed_output=True` writes BEDPE instead of BAM.
- **Output:** A dictionary with `status` (`completed` when the exit code is 0 and every output file exists, otherwise `failed`), `backend` (`deeptools`), `tool`, `command` (the tool name followed by its arguments), `return_code`, `outputs` (the output files that exist) and `log_tail` (the last 4000 characters of the tool output).
- **Example code:**

```python
import os
import shutil
import sys
import tempfile

import pysam

os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]
work = tempfile.mkdtemp(dir=".")
shutil.copy("reads.sam", work)
os.chdir(work)
pysam.view("-b", "-o", "a.bam", "reads.sam", catch_stdout=False)
pysam.index("a.bam")

from dnakit.epigenome import alignment_sieve

result = alignment_sieve("a.bam", "shifted.bam", shift=[1, -1, 1, -1], verbose=False)
print(result["status"], result["tool"], os.path.basename(result["outputs"][0]))
with pysam.AlignmentFile("a.bam") as before, pysam.AlignmentFile("shifted.bam") as after:
    first_before = [read.reference_start for read in before][:3]
    first_after = [read.reference_start for read in after][:3]
print(first_before, first_after)
```

- **Example output:**

```text
completed alignmentSieve shifted.bam
[99, 249, 261] [100, 249, 262]
```

## 2) `dnakit.epigenome.bam_compare()`

- **Purpose:** Compare two BAM files bin by bin and write the result as a bigWig or bedGraph track.
- **Method:** Check the inputs, build a fixed argument list and run `bamCompare` without a shell in the folder of the output file (time limit 300 s). The arguments are `-b1`, `-b2`, `-o`, `-of`, `--binSize`, `--operation`, `--pseudocount` and `-p`. No scaling option is exposed, so bamCompare uses its default scaling method (`readCount`).
- **API:** `dnakit.epigenome.bam_compare(bam1: 'str', bam2: 'str', output: 'str', *, operation: '_OPERATIONS' = 'log2', output_format: "Literal['bigwig', 'bedgraph']" = 'bigwig', bin_size: 'int' = 50, pseudocount: 'float' = 1.0, number_of_processors: 'int' = 1, deeptools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `bamCompare` (installed with deepTools, extra `dnakit[bio-deeptools]`) on PATH, or its full path in `deeptools_path`. The example also uses pysam and pyBigWig, which deepTools installs.
- **Input:** `bam1` and `bam2` must be existing indexed BAM files and `output` a new path. `operation` is one of `log2`, `ratio`, `subtract`, `add`, `mean`, `reciprocal_ratio`, `first`, `second`; `output_format` is `bigwig` or `bedgraph`; `bin_size` is the bin width in bases; `pseudocount` is added to both samples before a ratio.
- **Output:** A dictionary with `status` (`completed` when the exit code is 0 and every output file exists, otherwise `failed`), `backend` (`deeptools`), `tool`, `command` (the tool name followed by its arguments), `return_code`, `outputs` (the output files that exist) and `log_tail` (the last 4000 characters of the tool output).
- **Example code:**

```python
import os
import shutil
import sys
import tempfile

import pysam

os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]
work = tempfile.mkdtemp(dir=".")
shutil.copy("reads.sam", work)
os.chdir(work)
pysam.view("-b", "-o", "a.bam", "reads.sam", catch_stdout=False)
pysam.index("a.bam")

from dnakit.epigenome import bam_compare

with pysam.AlignmentFile("a.bam") as source, pysam.AlignmentFile("b.bam", "wb", template=source) as half:
    for number, read in enumerate(source):
        if number % 2 == 0:
            half.write(read)
pysam.index("b.bam")
result = bam_compare("a.bam", "b.bam", "diff.bedgraph", operation="subtract", output_format="bedgraph", bin_size=500, verbose=False)
print(result["status"], result["return_code"])
print(open("diff.bedgraph").read().splitlines()[:3])
```

- **Example output:**

```text
completed 0
['chr1\t0\t500\t-0.5', 'chr1\t500\t2500\t0', 'chr1\t2500\t3000\t0.5']
```

## 3) `dnakit.epigenome.bam_coverage()`

- **Purpose:** Build a coverage track (bigWig or bedGraph) from one BAM file.
- **Method:** Check the inputs, build a fixed argument list and run `bamCoverage` without a shell in the folder of the output file (time limit 300 s). The arguments are `-b`, `-o`, `-of`, `--binSize`, `--normalizeUsing` and `-p`, plus `--effectiveGenomeSize`, `--extendReads`, `--scaleFactor` and `--region` when they are given.
- **API:** `dnakit.epigenome.bam_coverage(bam: 'str', output: 'str', *, output_format: "Literal['bigwig', 'bedgraph']" = 'bigwig', bin_size: 'int' = 50, normalize_using: '_NORMALIZATION' = 'None', effective_genome_size: 'int | None' = None, extend_reads: 'int | None' = None, scale_factor: 'float | None' = None, region: 'str | None' = None, number_of_processors: 'int' = 1, deeptools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `bamCoverage` (installed with deepTools, extra `dnakit[bio-deeptools]`) on PATH, or its full path in `deeptools_path`. The example also uses pysam and pyBigWig, which deepTools installs.
- **Input:** `bam` must be an existing indexed BAM file and `output` a new path. `bin_size` must be at least 1. `normalize_using` is `RPGC`, `CPM`, `BPM`, `RPKM` or `None` (the default string `"None"`); `RPGC` needs `effective_genome_size` (deepTools reports the error in the log and the status becomes `failed`). `extend_reads`, `scale_factor` and `region` (for example `chr1:1000:2000`) are optional.
- **Output:** A dictionary with `status` (`completed` when the exit code is 0 and every output file exists, otherwise `failed`), `backend` (`deeptools`), `tool`, `command` (the tool name followed by its arguments), `return_code`, `outputs` (the output files that exist) and `log_tail` (the last 4000 characters of the tool output).
- **Example code:**

```python
import os
import shutil
import sys
import tempfile

import pysam

os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]
work = tempfile.mkdtemp(dir=".")
shutil.copy("reads.sam", work)
os.chdir(work)
pysam.view("-b", "-o", "a.bam", "reads.sam", catch_stdout=False)
pysam.index("a.bam")

from dnakit.epigenome import bam_coverage

result = bam_coverage("a.bam", "cov.bedgraph", output_format="bedgraph", bin_size=500, verbose=False)
print(result["status"], result["tool"], result["command"][0])
print(open("cov.bedgraph").read().splitlines()[:4])
```

- **Example output:**

```text
completed bamCoverage bamCoverage
['chr1\t0\t500\t5', 'chr1\t500\t2500\t6', 'chr1\t2500\t3000\t3', 'chr2\t0\t1000\t5']
```

## 4) `dnakit.epigenome.bam_pe_fragment_size()`

- **Purpose:** Estimate the fragment length distribution of paired-end BAM files and plot it as a histogram.
- **Method:** Check the inputs, build a fixed argument list and run `bamPEFragmentSize` without a shell in the folder of the output file (time limit 300 s). The arguments are `-b`, `-hist` and `-p`, plus `--table` when `table_output` is given and `--maxFragmentLength` when `max_fragment_length` is given. The tool does not reject a single-end BAM file: it completes and only warns in the log that the data may not be paired-end.
- **API:** `dnakit.epigenome.bam_pe_fragment_size(bams: 'Sequence[str]', histogram_output: 'str', *, table_output: 'str | None' = None, max_fragment_length: 'int | None' = None, number_of_processors: 'int' = 1, deeptools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `bamPEFragmentSize` (installed with deepTools, extra `dnakit[bio-deeptools]`) on PATH, or its full path in `deeptools_path`. The example also uses pysam and pyBigWig, which deepTools installs.
- **Input:** `bams` is a non-empty list of existing paired-end BAM files; `histogram_output` is a new image path; `table_output` is an optional new path for the tab-separated summary table; `max_fragment_length` is optional.
- **Output:** A dictionary with `status` (`completed` when the exit code is 0 and every output file exists, otherwise `failed`), `backend` (`deeptools`), `tool`, `command` (the tool name followed by its arguments), `return_code`, `outputs` (the output files that exist) and `log_tail` (the last 4000 characters of the tool output).
- **Example code:**

```python
import os
import shutil
import sys
import tempfile

import pysam

os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]
work = tempfile.mkdtemp(dir=".")
shutil.copy("reads.sam", work)
os.chdir(work)
pysam.view("-b", "-o", "a.bam", "reads.sam", catch_stdout=False)
pysam.index("a.bam")

from dnakit.epigenome import bam_pe_fragment_size

result = bam_pe_fragment_size(["a.bam"], "hist.png", table_output="sizes.tab", verbose=False)
print(result["status"], [os.path.basename(path) for path in result["outputs"]])
print(open("sizes.tab").read().splitlines()[0].split("\t")[1:5])
print(open("sizes.tab").read().splitlines()[1].split("\t")[1:5])
```

- **Example output:**

```text
completed ['hist.png', 'sizes.tab']
['Frag. Sampled', 'Frag. Len. Min.', 'Frag. Len. 1st. Qu.', 'Frag. Len. Mean']
['24', '200', '200.0', '200.0']
```

## 5) `dnakit.epigenome.bigwig_compare()`

- **Purpose:** Compare two bigWig tracks bin by bin and write a bigWig.
- **Method:** Check the inputs, build a fixed argument list and run `bigwigCompare` without a shell in the folder of the output file (time limit 300 s). The arguments are `-b1`, `-b2`, `-o`, `--binSize`, `--operation`, `--pseudocount` and `-p`. The output is always a bigWig.
- **API:** `dnakit.epigenome.bigwig_compare(bigwig1: 'str', bigwig2: 'str', output: 'str', *, operation: '_OPERATIONS' = 'log2', bin_size: 'int' = 50, pseudocount: 'float' = 1.0, number_of_processors: 'int' = 1, deeptools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `bigwigCompare` (installed with deepTools, extra `dnakit[bio-deeptools]`) on PATH, or its full path in `deeptools_path`. The example also uses pysam and pyBigWig, which deepTools installs.
- **Input:** `bigwig1` and `bigwig2` must be existing bigWig files and `output` a new path. `operation` takes the same values as in `bam_compare`; `bin_size` is the bin width in bases; `pseudocount` is added before ratios.
- **Output:** A dictionary with `status` (`completed` when the exit code is 0 and every output file exists, otherwise `failed`), `backend` (`deeptools`), `tool`, `command` (the tool name followed by its arguments), `return_code`, `outputs` (the output files that exist) and `log_tail` (the last 4000 characters of the tool output).
- **Example code:**

```python
import os
import shutil
import sys
import tempfile

os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]
work = tempfile.mkdtemp(dir=".")
os.chdir(work)
import pyBigWig


def write_bigwig(path, values):
    bw = pyBigWig.open(path, "w")
    bw.addHeader([("chr1", 3000), ("chr2", 2000)])
    bw.addEntries(
        ["chr1"] * 3 + ["chr2"] * 2,
        [0, 1000, 2000, 0, 1000],
        ends=[1000, 2000, 3000, 1000, 2000],
        values=values,
    )
    bw.close()

from dnakit.epigenome import bigwig_compare

write_bigwig("one.bw", [4.0, 6.0, 8.0, 2.0, 2.0])
write_bigwig("two.bw", [2.0, 2.0, 2.0, 1.0, 1.0])
result = bigwig_compare("one.bw", "two.bw", "diff.bw", operation="subtract", bin_size=1000, verbose=False)
print(result["status"], result["return_code"])
print(pyBigWig.open("diff.bw").intervals("chr1"))
```

- **Example output:**

```text
completed 0
((0, 1000, 2.0), (1000, 2000, 4.0), (2000, 3000, 6.0))
```

## 6) `dnakit.epigenome.estimate_scale_factor()`

- **Purpose:** Estimate scale factors for two or more BAM files from sampled read counts.
- **Method:** Check the inputs, build the argument list (`-b`, `-l`, `-n`, `-nl`, `-p`) and run `estimateScaleFactor` without a shell in the folder of the first BAM file. Nothing is written by the tool; the factors it prints are returned in the result.
- **API:** `dnakit.epigenome.estimate_scale_factor(bams: 'Sequence[str]', *, sample_window_length: 'int' = 1000, number_of_samples: 'int' = 1000000, normalization_length: 'int' = 1000, number_of_processors: 'int' = 1, deeptools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `estimateScaleFactor` (installed with deepTools, extra `dnakit[bio-deeptools]`) on PATH, or its full path in `deeptools_path`. The example also uses pysam and pyBigWig, which deepTools installs.
- **Input:** `bams` must contain at least two existing BAM files (a `ConfigurationError` is raised otherwise). `sample_window_length` (`-l`), `number_of_samples` (`-n`) and `normalization_length` (`-nl`) control the sampling; `number_of_processors` is `-p`.
- **Output:** `status` (`completed` when the exit code is 0, otherwise `failed`), `backend`, `tool`, `command`, `return_code` and `stdout_tail` (the last 4000 characters of the tool output, which hold the factors and statistics, for example `size_factors_based_on_median`). There is no `outputs` entry.
- **Example code:**

```python
import os
import shutil
import sys
import tempfile

import pysam

os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]
work = tempfile.mkdtemp(dir=".")
shutil.copy("reads.sam", work)
os.chdir(work)
pysam.view("-b", "-o", "a.bam", "reads.sam", catch_stdout=False)
pysam.index("a.bam")

from dnakit.epigenome import estimate_scale_factor

shutil.copy("a.bam", "b.bam")
shutil.copy("a.bam.bai", "b.bam.bai")
result = estimate_scale_factor(["a.bam", "b.bam"], sample_window_length=500, number_of_samples=100, verbose=False)
print(result["status"], result["return_code"])
for line in result["stdout_tail"].splitlines():
    if line.startswith(("size_factors_based_on_median", "sites_sampled")):
        print(line)
```

- **Example output:**

```text
completed 0
size_factors_based_on_median: [1. 1.]
sites_sampled: 82
```

## 7) `dnakit.epigenome.plot_coverage()`

- **Purpose:** Plot the coverage distribution of one or more BAM files.
- **Method:** Check the inputs, build a fixed argument list and run `plotCoverage` without a shell in the folder of the output file (time limit 300 s). The arguments are `-b`, `-o` and `-n`, plus `--labels` when `labels` is given. The tool also prints a coverage statistics table, which is returned in `log_tail`.
- **API:** `dnakit.epigenome.plot_coverage(bams: 'Sequence[str]', plot_output: 'str', *, labels: 'Sequence[str] | None' = None, number_of_samples: 'int' = 1000000, deeptools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `plotCoverage` (installed with deepTools, extra `dnakit[bio-deeptools]`) on PATH, or its full path in `deeptools_path`. The example also uses pysam and pyBigWig, which deepTools installs.
- **Input:** `bams` is a non-empty list of existing indexed BAM files; `plot_output` is a new image path; `labels` is optional and must have one entry per BAM file; `number_of_samples` is the number of genome positions sampled.
- **Output:** A dictionary with `status` (`completed` when the exit code is 0 and every output file exists, otherwise `failed`), `backend` (`deeptools`), `tool`, `command` (the tool name followed by its arguments), `return_code`, `outputs` (the output files that exist) and `log_tail` (the last 4000 characters of the tool output).
- **Example code:**

```python
import os
import shutil
import sys
import tempfile

import pysam

os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]
work = tempfile.mkdtemp(dir=".")
shutil.copy("reads.sam", work)
os.chdir(work)
pysam.view("-b", "-o", "a.bam", "reads.sam", catch_stdout=False)
pysam.index("a.bam")

from dnakit.epigenome import plot_coverage

result = plot_coverage(["a.bam"], "coverage.png", labels=["sample"], number_of_samples=100, verbose=False)
print(result["status"], os.path.isfile("coverage.png"))
print(result["log_tail"].splitlines()[0].split("	"))
```

- **Example output:**

```text
completed True
['sample', 'mean', 'std', 'min', '25%', '50%', '75%', 'max']
```

## 8) `dnakit.epigenome.plot_fingerprint()`

- **Purpose:** Plot the cumulative read-coverage fingerprint of ChIP-like BAM files.
- **Method:** Check the inputs, build a fixed argument list and run `plotFingerprint` without a shell in the folder of the output file (time limit 300 s). The arguments are `-b`, `-plot` and `-p`, plus `--extendReads` and `--minMappingQuality` when they are given.
- **API:** `dnakit.epigenome.plot_fingerprint(bams: 'Sequence[str]', plot_output: 'str', *, extend_reads: 'int | None' = None, min_mapping_quality: 'int | None' = None, number_of_processors: 'int' = 1, deeptools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `plotFingerprint` (installed with deepTools, extra `dnakit[bio-deeptools]`) on PATH, or its full path in `deeptools_path`. The example also uses pysam and pyBigWig, which deepTools installs.
- **Input:** `bams` is a non-empty list of existing indexed BAM files; `plot_output` is a new image path; `extend_reads` and `min_mapping_quality` are optional integers.
- **Output:** A dictionary with `status` (`completed` when the exit code is 0 and every output file exists, otherwise `failed`), `backend` (`deeptools`), `tool`, `command` (the tool name followed by its arguments), `return_code`, `outputs` (the output files that exist) and `log_tail` (the last 4000 characters of the tool output).
- **Example code:**

```python
import os
import shutil
import sys
import tempfile

import pysam

os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]
work = tempfile.mkdtemp(dir=".")
shutil.copy("reads.sam", work)
os.chdir(work)
pysam.view("-b", "-o", "a.bam", "reads.sam", catch_stdout=False)
pysam.index("a.bam")

from dnakit.epigenome import plot_fingerprint

result = plot_fingerprint(["a.bam"], "fingerprint.png", min_mapping_quality=10, verbose=False)
print(result["status"], os.path.isfile("fingerprint.png"))
print([argument for argument in result["command"] if argument.startswith("-")])
```

- **Example output:**

```text
completed True
['-b', '-plot', '-p', '--minMappingQuality']
```

## 9) `dnakit.epigenome.bigwig_average()`

- **Purpose:** Average several bigWig tracks bin by bin into one bigWig.
- **Method:** Check the inputs, build a fixed argument list and run `bigwigAverage` without a shell in the folder of the output file (time limit 300 s). The arguments are `-b`, `-o`, `--binSize` and `-p`.
- **API:** `dnakit.epigenome.bigwig_average(bigwigs: 'Sequence[str]', output: 'str', *, bin_size: 'int' = 50, number_of_processors: 'int' = 1, deeptools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `bigwigAverage` (installed with deepTools, extra `dnakit[bio-deeptools]`) on PATH, or its full path in `deeptools_path`. The example also uses pysam and pyBigWig, which deepTools installs.
- **Input:** `bigwigs` is a non-empty list of existing bigWig files; `output` is a new path; `bin_size` is the bin width in bases.
- **Output:** A dictionary with `status` (`completed` when the exit code is 0 and every output file exists, otherwise `failed`), `backend` (`deeptools`), `tool`, `command` (the tool name followed by its arguments), `return_code`, `outputs` (the output files that exist) and `log_tail` (the last 4000 characters of the tool output).
- **Example code:**

```python
import os
import shutil
import sys
import tempfile

os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]
work = tempfile.mkdtemp(dir=".")
os.chdir(work)
import pyBigWig


def write_bigwig(path, values):
    bw = pyBigWig.open(path, "w")
    bw.addHeader([("chr1", 3000), ("chr2", 2000)])
    bw.addEntries(
        ["chr1"] * 3 + ["chr2"] * 2,
        [0, 1000, 2000, 0, 1000],
        ends=[1000, 2000, 3000, 1000, 2000],
        values=values,
    )
    bw.close()

from dnakit.epigenome import bigwig_average

write_bigwig("one.bw", [1.0, 2.0, 3.0, 4.0, 5.0])
write_bigwig("two.bw", [3.0, 2.0, 1.0, 0.0, 1.0])
result = bigwig_average(["one.bw", "two.bw"], "mean.bw", bin_size=1000, verbose=False)
print(result["status"], result["return_code"])
print(pyBigWig.open("mean.bw").intervals("chr1"))
```

- **Example output:**

```text
completed 0
((0, 3000, 2.0),)
```

## 10) `dnakit.epigenome.compute_matrix_reference_point()`

- **Purpose:** Compute a matrix of signal values around a reference point of each region (TSS, TES or centre).
- **Method:** Check the inputs, build a fixed argument list and run `computeMatrix` without a shell in the folder of the output file (time limit 300 s). The subcommand is `reference-point` with `--referencePoint`, `-S`, `-R`, `-o`, `-b`, `-a`, `--binSize` and `-p`. The output is deepTools' gzipped matrix file, the input of `plot_heatmap` and `plot_profile`.
- **API:** `dnakit.epigenome.compute_matrix_reference_point(score_files: 'Sequence[str]', regions: 'Sequence[str]', output_matrix: 'str', *, reference_point: "Literal['TSS', 'TES', 'center']" = 'TSS', before_region_length: 'int' = 500, after_region_length: 'int' = 500, bin_size: 'int' = 10, number_of_processors: 'int' = 1, deeptools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `computeMatrix` (installed with deepTools, extra `dnakit[bio-deeptools]`) on PATH, or its full path in `deeptools_path`. The example also uses pysam and pyBigWig, which deepTools installs.
- **Input:** `score_files` (bigWig files) and `regions` (BED or GTF files) are non-empty lists of existing files; `output_matrix` is a new path; `reference_point` is `TSS`, `TES` or `center`; `before_region_length` and `after_region_length` are the distances in bases around the point; `bin_size` is the bin width.
- **Output:** A dictionary with `status` (`completed` when the exit code is 0 and every output file exists, otherwise `failed`), `backend` (`deeptools`), `tool`, `command` (the tool name followed by its arguments), `return_code`, `outputs` (the output files that exist) and `log_tail` (the last 4000 characters of the tool output).
- **Example code:**

```python
import os
import shutil
import sys
import tempfile

os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]
work = tempfile.mkdtemp(dir=".")
shutil.copy("genes.bed", work)
os.chdir(work)
import pyBigWig


def write_bigwig(path, values):
    bw = pyBigWig.open(path, "w")
    bw.addHeader([("chr1", 3000), ("chr2", 2000)])
    bw.addEntries(
        ["chr1"] * 3 + ["chr2"] * 2,
        [0, 1000, 2000, 0, 1000],
        ends=[1000, 2000, 3000, 1000, 2000],
        values=values,
    )
    bw.close()

import gzip
import json

from dnakit.epigenome import compute_matrix_reference_point

write_bigwig("one.bw", [1.0, 2.0, 3.0, 4.0, 5.0])
result = compute_matrix_reference_point(
    ["one.bw"], ["genes.bed"], "matrix.gz", reference_point="TSS",
    before_region_length=500, after_region_length=500, bin_size=100, verbose=False,
)
print(result["status"], result["return_code"])
with gzip.open("matrix.gz", "rt") as handle:
    meta = json.loads(handle.readline()[1:])
    rows = handle.read().splitlines()
print(meta["ref point"], meta["sample_boundaries"], meta["group_boundaries"], len(rows))
```

- **Example output:**

```text
completed 0
['TSS'] [0, 10] [0, 2] 2
```

## 11) `dnakit.epigenome.compute_matrix_scale_regions()`

- **Purpose:** Compute a matrix of signal values over regions that are scaled to a common length, with flanks.
- **Method:** Check the inputs, build a fixed argument list and run `computeMatrix` without a shell in the folder of the output file (time limit 300 s). The subcommand is `scale-regions` with `-S`, `-R`, `-o`, `-m`, `-b`, `-a`, `--binSize` and `-p`. Every region is stretched or shrunk to `region_body_length`.
- **API:** `dnakit.epigenome.compute_matrix_scale_regions(score_files: 'Sequence[str]', regions: 'Sequence[str]', output_matrix: 'str', *, region_body_length: 'int' = 1000, before_region_length: 'int' = 500, after_region_length: 'int' = 500, bin_size: 'int' = 10, number_of_processors: 'int' = 1, deeptools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `computeMatrix` (installed with deepTools, extra `dnakit[bio-deeptools]`) on PATH, or its full path in `deeptools_path`. The example also uses pysam and pyBigWig, which deepTools installs.
- **Input:** `score_files` and `regions` are non-empty lists of existing files; `output_matrix` is a new path; `region_body_length` is the length every region is scaled to; `before_region_length` and `after_region_length` are the flank sizes; `bin_size` is the bin width.
- **Output:** A dictionary with `status` (`completed` when the exit code is 0 and every output file exists, otherwise `failed`), `backend` (`deeptools`), `tool`, `command` (the tool name followed by its arguments), `return_code`, `outputs` (the output files that exist) and `log_tail` (the last 4000 characters of the tool output).
- **Example code:**

```python
import os
import shutil
import sys
import tempfile

os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]
work = tempfile.mkdtemp(dir=".")
shutil.copy("genes.bed", work)
os.chdir(work)
import pyBigWig


def write_bigwig(path, values):
    bw = pyBigWig.open(path, "w")
    bw.addHeader([("chr1", 3000), ("chr2", 2000)])
    bw.addEntries(
        ["chr1"] * 3 + ["chr2"] * 2,
        [0, 1000, 2000, 0, 1000],
        ends=[1000, 2000, 3000, 1000, 2000],
        values=values,
    )
    bw.close()

import gzip
import json

from dnakit.epigenome import compute_matrix_scale_regions

write_bigwig("one.bw", [1.0, 2.0, 3.0, 4.0, 5.0])
result = compute_matrix_scale_regions(
    ["one.bw"], ["genes.bed"], "scaled.gz", region_body_length=300,
    before_region_length=100, after_region_length=100, bin_size=50, verbose=False,
)
print(result["status"], result["return_code"])
with gzip.open("scaled.gz", "rt") as handle:
    meta = json.loads(handle.readline()[1:])
    rows = handle.read().splitlines()
print(meta["body"], meta["sample_boundaries"], len(rows))
```

- **Example output:**

```text
completed 0
[300] [0, 10] 2
```

## 12) `dnakit.epigenome.multi_bam_summary()`

- **Purpose:** Count reads of several BAM files in consecutive genomic bins and save the table for correlation or PCA plots.
- **Method:** Check the inputs, build a fixed argument list and run `multiBamSummary` without a shell in the folder of the output file (time limit 300 s). The subcommand is `bins` with `-b`, `-o`, `--binSize`, `-n` (distance between bins) and `-p`, plus `--labels` when `labels` is given. The output is a NumPy `.npz` file.
- **API:** `dnakit.epigenome.multi_bam_summary(bams: 'Sequence[str]', output: 'str', *, bin_size: 'int' = 10000, distance_between_bins: 'int' = 0, labels: 'Sequence[str] | None' = None, number_of_processors: 'int' = 1, deeptools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `multiBamSummary` (installed with deepTools, extra `dnakit[bio-deeptools]`) on PATH, or its full path in `deeptools_path`. The example also uses pysam and pyBigWig, which deepTools installs.
- **Input:** `bams` is a non-empty list of existing indexed BAM files; `output` is a new path; `bin_size` is the bin width; `distance_between_bins` is the gap between sampled bins (0 means consecutive bins); `labels` is optional and must have one entry per BAM file.
- **Output:** A dictionary with `status` (`completed` when the exit code is 0 and every output file exists, otherwise `failed`), `backend` (`deeptools`), `tool`, `command` (the tool name followed by its arguments), `return_code`, `outputs` (the output files that exist) and `log_tail` (the last 4000 characters of the tool output).
- **Example code:**

```python
import os
import shutil
import sys
import tempfile

import pysam

os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]
work = tempfile.mkdtemp(dir=".")
shutil.copy("reads.sam", work)
os.chdir(work)
pysam.view("-b", "-o", "a.bam", "reads.sam", catch_stdout=False)
pysam.index("a.bam")

import numpy as np

from dnakit.epigenome import multi_bam_summary

shutil.copy("a.bam", "b.bam")
shutil.copy("a.bam.bai", "b.bam.bai")
result = multi_bam_summary(["a.bam", "b.bam"], "summary.npz", bin_size=500, labels=["A", "B"], verbose=False)
print(result["status"], result["return_code"])
with np.load("summary.npz") as data:
    print(sorted(data.files), data["matrix"].shape, data["matrix"].sum(axis=0).tolist())
```

- **Example output:**

```text
completed 0
['labels', 'matrix'] (10, 2) [50.0, 50.0]
```

## 13) `dnakit.epigenome.multi_bigwig_summary()`

- **Purpose:** Summarise the values of several bigWig tracks in genomic bins and save the table.
- **Method:** Check the inputs, build a fixed argument list and run `multiBigwigSummary` without a shell in the folder of the output file (time limit 300 s). The subcommand is `bins` with `-b`, `-o`, `--binSize` and `-p`, plus `--labels` when `labels` is given. The output is a NumPy `.npz` file.
- **API:** `dnakit.epigenome.multi_bigwig_summary(bigwigs: 'Sequence[str]', output: 'str', *, bin_size: 'int' = 10000, labels: 'Sequence[str] | None' = None, number_of_processors: 'int' = 1, deeptools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `multiBigwigSummary` (installed with deepTools, extra `dnakit[bio-deeptools]`) on PATH, or its full path in `deeptools_path`. The example also uses pysam and pyBigWig, which deepTools installs.
- **Input:** `bigwigs` is a non-empty list of existing bigWig files; `output` is a new path; `bin_size` is the bin width; `labels` is optional and must have one entry per bigWig file.
- **Output:** A dictionary with `status` (`completed` when the exit code is 0 and every output file exists, otherwise `failed`), `backend` (`deeptools`), `tool`, `command` (the tool name followed by its arguments), `return_code`, `outputs` (the output files that exist) and `log_tail` (the last 4000 characters of the tool output).
- **Example code:**

```python
import os
import shutil
import sys
import tempfile

os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]
work = tempfile.mkdtemp(dir=".")
os.chdir(work)
import pyBigWig


def write_bigwig(path, values):
    bw = pyBigWig.open(path, "w")
    bw.addHeader([("chr1", 3000), ("chr2", 2000)])
    bw.addEntries(
        ["chr1"] * 3 + ["chr2"] * 2,
        [0, 1000, 2000, 0, 1000],
        ends=[1000, 2000, 3000, 1000, 2000],
        values=values,
    )
    bw.close()

import numpy as np

from dnakit.epigenome import multi_bigwig_summary

write_bigwig("one.bw", [1.0, 2.0, 3.0, 4.0, 5.0])
write_bigwig("two.bw", [2.0, 1.0, 4.0, 3.0, 6.0])
result = multi_bigwig_summary(["one.bw", "two.bw"], "summary.npz", bin_size=1000, labels=["one", "two"], verbose=False)
print(result["status"], result["return_code"])
with np.load("summary.npz") as data:
    print(sorted(data.files), [str(label) for label in data["labels"]], data["matrix"].shape)
```

- **Example output:**

```text
completed 0
['labels', 'matrix'] ['one', 'two'] (5, 2)
```

## 14) `dnakit.epigenome.plot_correlation()`

- **Purpose:** Plot the correlation between samples from a multiBamSummary or multiBigwigSummary table.
- **Method:** Check the inputs, build a fixed argument list and run `plotCorrelation` without a shell in the folder of the output file (time limit 300 s). The arguments are `-in`, `-c` (correlation method), `-p` (plot type) and `-o`.
- **API:** `dnakit.epigenome.plot_correlation(summary_matrix: 'str', output_figure: 'str', *, method: "Literal['spearman', 'pearson']" = 'spearman', plot_type: "Literal['heatmap', 'scatterplot']" = 'heatmap', deeptools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `plotCorrelation` (installed with deepTools, extra `dnakit[bio-deeptools]`) on PATH, or its full path in `deeptools_path`. The example also uses pysam and pyBigWig, which deepTools installs.
- **Input:** `summary_matrix` is an existing `.npz` file from `multi_bam_summary` or `multi_bigwig_summary`; `output_figure` is a new image path; `method` is `spearman` or `pearson`; `plot_type` is `heatmap` or `scatterplot`.
- **Output:** A dictionary with `status` (`completed` when the exit code is 0 and every output file exists, otherwise `failed`), `backend` (`deeptools`), `tool`, `command` (the tool name followed by its arguments), `return_code`, `outputs` (the output files that exist) and `log_tail` (the last 4000 characters of the tool output).
- **Example code:**

```python
import os
import shutil
import sys
import tempfile

os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]
work = tempfile.mkdtemp(dir=".")
os.chdir(work)
import pyBigWig


def write_bigwig(path, values):
    bw = pyBigWig.open(path, "w")
    bw.addHeader([("chr1", 3000), ("chr2", 2000)])
    bw.addEntries(
        ["chr1"] * 3 + ["chr2"] * 2,
        [0, 1000, 2000, 0, 1000],
        ends=[1000, 2000, 3000, 1000, 2000],
        values=values,
    )
    bw.close()

from dnakit.epigenome import multi_bigwig_summary, plot_correlation

write_bigwig("one.bw", [1.0, 2.0, 3.0, 4.0, 5.0])
write_bigwig("two.bw", [2.0, 1.0, 4.0, 3.0, 6.0])
write_bigwig("three.bw", [5.0, 4.0, 3.0, 2.0, 1.5])
multi_bigwig_summary(["one.bw", "two.bw", "three.bw"], "summary.npz", bin_size=1000, labels=["one", "two", "three"], verbose=False)
result = plot_correlation("summary.npz", "correlation.png", method="pearson", plot_type="heatmap", verbose=False)
print(result["status"], os.path.isfile("correlation.png"))
print([argument for argument in result["command"] if argument.startswith("-")])
```

- **Example output:**

```text
completed True
['-in', '-c', '-p', '-o']
```

## 15) `dnakit.epigenome.plot_heatmap()`

- **Purpose:** Plot a heatmap from a computeMatrix matrix file.
- **Method:** Check the inputs, build a fixed argument list and run `plotHeatmap` without a shell in the folder of the output file (time limit 300 s). The arguments are `-m`, `-out` and `--dpi`, plus `--kmeans` when `kmeans` is given.
- **API:** `dnakit.epigenome.plot_heatmap(matrix: 'str', output_figure: 'str', *, kmeans: 'int | None' = None, dpi: 'int' = 200, deeptools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `plotHeatmap` (installed with deepTools, extra `dnakit[bio-deeptools]`) on PATH, or its full path in `deeptools_path`. The example also uses pysam and pyBigWig, which deepTools installs.
- **Input:** `matrix` is an existing file written by `compute_matrix_reference_point` or `compute_matrix_scale_regions`; `output_figure` is a new image path; `kmeans` is `None` or a positive integer (otherwise a `ConfigurationError` is raised); `dpi` is the image resolution.
- **Output:** A dictionary with `status` (`completed` when the exit code is 0 and every output file exists, otherwise `failed`), `backend` (`deeptools`), `tool`, `command` (the tool name followed by its arguments), `return_code`, `outputs` (the output files that exist) and `log_tail` (the last 4000 characters of the tool output).
- **Example code:**

```python
import os
import shutil
import sys
import tempfile

os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]
work = tempfile.mkdtemp(dir=".")
shutil.copy("genes.bed", work)
os.chdir(work)
import pyBigWig


def write_bigwig(path, values):
    bw = pyBigWig.open(path, "w")
    bw.addHeader([("chr1", 3000), ("chr2", 2000)])
    bw.addEntries(
        ["chr1"] * 3 + ["chr2"] * 2,
        [0, 1000, 2000, 0, 1000],
        ends=[1000, 2000, 3000, 1000, 2000],
        values=values,
    )
    bw.close()

from dnakit.epigenome import compute_matrix_reference_point, plot_heatmap

write_bigwig("one.bw", [1.0, 2.0, 3.0, 4.0, 5.0])
compute_matrix_reference_point(["one.bw"], ["genes.bed"], "matrix.gz", bin_size=100, verbose=False)
result = plot_heatmap("matrix.gz", "heatmap.png", dpi=100, verbose=False)
print(result["status"], os.path.isfile("heatmap.png"))
print([argument for argument in result["command"] if argument.startswith("-")])
```

- **Example output:**

```text
completed True
['-m', '-out', '--dpi']
```

## 16) `dnakit.epigenome.plot_pca()`

- **Purpose:** Plot a principal component analysis of the samples in a multiBamSummary or multiBigwigSummary table.
- **Method:** Check the inputs, build a fixed argument list and run `plotPCA` without a shell in the folder of the output file (time limit 300 s). The arguments are `-in` and `-o`; no other plotPCA option is exposed.
- **API:** `dnakit.epigenome.plot_pca(summary_matrix: 'str', output_figure: 'str', *, deeptools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `plotPCA` (installed with deepTools, extra `dnakit[bio-deeptools]`) on PATH, or its full path in `deeptools_path`. The example also uses pysam and pyBigWig, which deepTools installs.
- **Input:** `summary_matrix` is an existing `.npz` file from `multi_bam_summary` or `multi_bigwig_summary`; `output_figure` is a new image path.
- **Output:** A dictionary with `status` (`completed` when the exit code is 0 and every output file exists, otherwise `failed`), `backend` (`deeptools`), `tool`, `command` (the tool name followed by its arguments), `return_code`, `outputs` (the output files that exist) and `log_tail` (the last 4000 characters of the tool output).
- **Example code:**

```python
import os
import shutil
import sys
import tempfile

os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]
work = tempfile.mkdtemp(dir=".")
os.chdir(work)
import pyBigWig


def write_bigwig(path, values):
    bw = pyBigWig.open(path, "w")
    bw.addHeader([("chr1", 3000), ("chr2", 2000)])
    bw.addEntries(
        ["chr1"] * 3 + ["chr2"] * 2,
        [0, 1000, 2000, 0, 1000],
        ends=[1000, 2000, 3000, 1000, 2000],
        values=values,
    )
    bw.close()

from dnakit.epigenome import multi_bigwig_summary, plot_pca

write_bigwig("one.bw", [1.0, 2.0, 3.0, 4.0, 5.0])
write_bigwig("two.bw", [2.0, 1.0, 4.0, 3.0, 6.0])
write_bigwig("three.bw", [5.0, 4.0, 3.0, 2.0, 1.5])
multi_bigwig_summary(["one.bw", "two.bw", "three.bw"], "summary.npz", bin_size=1000, labels=["one", "two", "three"], verbose=False)
result = plot_pca("summary.npz", "pca.png", verbose=False)
print(result["status"], os.path.isfile("pca.png"))
print(result["command"][0], result["return_code"])
```

- **Example output:**

```text
completed True
plotPCA 0
```

## 17) `dnakit.epigenome.plot_profile()`

- **Purpose:** Plot average signal profiles from a computeMatrix matrix file.
- **Method:** Check the inputs, build a fixed argument list and run `plotProfile` without a shell in the folder of the output file (time limit 300 s). The arguments are `-m`, `-out` and `--dpi`, plus `--kmeans` when `kmeans` is given.
- **API:** `dnakit.epigenome.plot_profile(matrix: 'str', output_figure: 'str', *, kmeans: 'int | None' = None, dpi: 'int' = 200, deeptools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `plotProfile` (installed with deepTools, extra `dnakit[bio-deeptools]`) on PATH, or its full path in `deeptools_path`. The example also uses pysam and pyBigWig, which deepTools installs.
- **Input:** `matrix` is an existing file written by `compute_matrix_reference_point` or `compute_matrix_scale_regions`; `output_figure` is a new image path; `kmeans` is `None` or a positive integer (otherwise a `ConfigurationError` is raised); `dpi` is the image resolution.
- **Output:** A dictionary with `status` (`completed` when the exit code is 0 and every output file exists, otherwise `failed`), `backend` (`deeptools`), `tool`, `command` (the tool name followed by its arguments), `return_code`, `outputs` (the output files that exist) and `log_tail` (the last 4000 characters of the tool output).
- **Example code:**

```python
import os
import shutil
import sys
import tempfile

os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]
work = tempfile.mkdtemp(dir=".")
shutil.copy("genes.bed", work)
os.chdir(work)
import pyBigWig


def write_bigwig(path, values):
    bw = pyBigWig.open(path, "w")
    bw.addHeader([("chr1", 3000), ("chr2", 2000)])
    bw.addEntries(
        ["chr1"] * 3 + ["chr2"] * 2,
        [0, 1000, 2000, 0, 1000],
        ends=[1000, 2000, 3000, 1000, 2000],
        values=values,
    )
    bw.close()

from dnakit.epigenome import compute_matrix_reference_point, plot_profile

write_bigwig("one.bw", [1.0, 2.0, 3.0, 4.0, 5.0])
compute_matrix_reference_point(["one.bw"], ["genes.bed"], "matrix.gz", bin_size=100, verbose=False)
result = plot_profile("matrix.gz", "profile.png", dpi=100, verbose=False)
print(result["status"], os.path.isfile("profile.png"))
print([argument for argument in result["command"] if argument.startswith("-")])
```

- **Example output:**

```text
completed True
['-m', '-out', '--dpi']
```

## 18) `dnakit.epigenome.estimate_read_filtering()`

- **Purpose:** Estimate how many reads of each BAM file are removed by deepTools' read filters, broken down by reason.
- **Method:** Check the inputs, build a fixed argument list and run `estimateReadFiltering` without a shell in the folder of the output file (time limit 300 s). The arguments are `-b`, `-o`, `--binSize` and `-p`. No filter option is exposed, so the tool's defaults apply. The output is a tab-separated table with one row per BAM file.
- **API:** `dnakit.epigenome.estimate_read_filtering(bams: 'Sequence[str]', output: 'str', *, bin_size: 'int' = 1000, number_of_processors: 'int' = 1, deeptools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `estimateReadFiltering` (installed with deepTools, extra `dnakit[bio-deeptools]`) on PATH, or its full path in `deeptools_path`. The example also uses pysam and pyBigWig, which deepTools installs.
- **Input:** `bams` is a non-empty list of existing indexed BAM files; `output` is a new path; `bin_size` is the sampling bin width.
- **Output:** A dictionary with `status` (`completed` when the exit code is 0 and every output file exists, otherwise `failed`), `backend` (`deeptools`), `tool`, `command` (the tool name followed by its arguments), `return_code`, `outputs` (the output files that exist) and `log_tail` (the last 4000 characters of the tool output).
- **Example code:**

```python
import os
import shutil
import sys
import tempfile

import pysam

os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]
work = tempfile.mkdtemp(dir=".")
shutil.copy("reads.sam", work)
os.chdir(work)
pysam.view("-b", "-o", "a.bam", "reads.sam", catch_stdout=False)
pysam.index("a.bam")

from dnakit.epigenome import estimate_read_filtering

result = estimate_read_filtering(["a.bam"], "filtering.tab", bin_size=500, verbose=False)
print(result["status"], result["return_code"])
rows = open("filtering.tab").read().splitlines()
print(rows[0].split("\t")[1:3])
print(rows[1].split("\t")[1:3])
```

- **Example output:**

```text
completed 0
['Total Reads', 'Mapped Reads']
['48', '48']
```

## 19) `dnakit.epigenome.plot_enrichment()`

- **Purpose:** Plot the fraction of reads that fall into each set of regions.
- **Method:** Check the inputs, build a fixed argument list and run `plotEnrichment` without a shell in the folder of the output file (time limit 300 s). The arguments are `-b`, `--BED` and `-o`, plus `--labels` when `labels` is given.
- **API:** `dnakit.epigenome.plot_enrichment(bams: 'Sequence[str]', regions: 'Sequence[str]', output_figure: 'str', *, labels: 'Sequence[str] | None' = None, deeptools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `plotEnrichment` (installed with deepTools, extra `dnakit[bio-deeptools]`) on PATH, or its full path in `deeptools_path`. The example also uses pysam and pyBigWig, which deepTools installs.
- **Input:** `bams` and `regions` (BED or GTF files) are non-empty lists of existing files; `output_figure` is a new image path; `labels` is optional and must have one entry per BAM file.
- **Output:** A dictionary with `status` (`completed` when the exit code is 0 and every output file exists, otherwise `failed`), `backend` (`deeptools`), `tool`, `command` (the tool name followed by its arguments), `return_code`, `outputs` (the output files that exist) and `log_tail` (the last 4000 characters of the tool output).
- **Example code:**

```python
import os
import shutil
import sys
import tempfile

import pysam

os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]
work = tempfile.mkdtemp(dir=".")
shutil.copy("reads.sam", work)
shutil.copy("genes.bed", work)
os.chdir(work)
pysam.view("-b", "-o", "a.bam", "reads.sam", catch_stdout=False)
pysam.index("a.bam")

from dnakit.epigenome import plot_enrichment

result = plot_enrichment(["a.bam"], ["genes.bed"], "enrichment.png", labels=["sample"], verbose=False)
print(result["status"], os.path.isfile("enrichment.png"))
print([argument for argument in result["command"] if argument.startswith("-")])
```

- **Example output:**

```text
completed True
['-b', '--BED', '-o', '--labels']
```

## 20) `dnakit.epigenome.compute_gc_bias()`

- **Purpose:** Compute the GC-bias frequency table of a BAM file against a 2bit genome.
- **Method:** Check the inputs, build a fixed argument list and run `computeGCBias` without a shell in the folder of the output file (time limit 300 s). The arguments are `-b`, `-g`, `--effectiveGenomeSize`, `-freq` and `-p`, plus `-l` when `fragment_length` is given. The table is the input of `correct_gc_bias`; the plot option of the tool is not exposed.
- **API:** `dnakit.epigenome.compute_gc_bias(bam: 'str', genome_2bit: 'str', effective_genome_size: 'int', frequencies_output: 'str', *, fragment_length: 'int | None' = None, number_of_processors: 'int' = 1, deeptools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `computeGCBias` (installed with deepTools, extra `dnakit[bio-deeptools]`) on PATH, or its full path in `deeptools_path`. The example also uses pysam and pyBigWig, which deepTools installs.
- **Input:** `bam` is an existing indexed BAM file and `genome_2bit` an existing 2bit file; `effective_genome_size` must be a positive integer; `frequencies_output` is a new path; `fragment_length` is `None` or a positive integer (single-end reads need it).
- **Output:** A dictionary with `status` (`completed` when the exit code is 0 and every output file exists, otherwise `failed`), `backend` (`deeptools`), `tool`, `command` (the tool name followed by its arguments), `return_code`, `outputs` (the output files that exist) and `log_tail` (the last 4000 characters of the tool output).
- **Example code:**

```python
import os
import shutil
import sys
import tempfile

os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]
work = tempfile.mkdtemp(dir=".")
os.chdir(work)
import random
import struct

import pysam


def write_twobit(path, name, seq):
    code = {"T": 0, "C": 1, "A": 2, "G": 3}
    packed = bytearray((len(seq) + 3) // 4)
    for i, base in enumerate(seq):
        packed[i // 4] |= code[base] << (6 - 2 * (i % 4))
    record = struct.pack("<IIII", len(seq), 0, 0, 0) + bytes(packed)
    index = struct.pack("B", len(name)) + name.encode() + struct.pack("<I", 16 + 1 + len(name) + 4)
    with open(path, "wb") as handle:
        handle.write(struct.pack("<IIII", 0x1A412743, 0, 1, 0) + index + record)


rng = random.Random(7)
genome = "".join(rng.choice("ACGT") for _ in range(20000))
write_twobit("genome.2bit", "chr1", genome)
header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": len(genome)}]}
starts = sorted(rng.randrange(0, 19800) for _ in range(1500))
with pysam.AlignmentFile("g.bam", "wb", header=header) as out:
    for i, start in enumerate(starts):
        read = pysam.AlignedSegment()
        read.query_name = f"r{i}"
        read.query_sequence = genome[start:start + 50]
        read.flag = 0
        read.reference_id = 0
        read.reference_start = start
        read.mapping_quality = 60
        read.cigarstring = "50M"
        read.query_qualities = pysam.qualitystring_to_array("I" * 50)
        out.write(read)
pysam.index("g.bam")

from dnakit.epigenome import compute_gc_bias

result = compute_gc_bias("g.bam", "genome.2bit", 20000, "freq.txt", fragment_length=200, verbose=False)
print(result["status"], result["tool"], os.path.basename(result["outputs"][0]))
rows = open("freq.txt").read().splitlines()
print(len(rows), len(rows[0].split()))
```

- **Example output:**

```text
completed computeGCBias freq.txt
201 3
```

## 21) `dnakit.epigenome.correct_gc_bias()`

- **Purpose:** Write a GC-corrected BAM file using a frequency table from `compute_gc_bias`.
- **Method:** Check the inputs, build a fixed argument list and run `correctGCBias` without a shell in the folder of the output file (time limit 300 s). The arguments are `-b`, `-g`, `--effectiveGenomeSize`, `-freq`, `-o` and `-p`.
- **API:** `dnakit.epigenome.correct_gc_bias(bam: 'str', genome_2bit: 'str', effective_genome_size: 'int', frequencies: 'str', output_bam: 'str', *, number_of_processors: 'int' = 1, deeptools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `correctGCBias` (installed with deepTools, extra `dnakit[bio-deeptools]`) on PATH, or its full path in `deeptools_path`. The example also uses pysam and pyBigWig, which deepTools installs.
- **Input:** `bam` is an existing indexed BAM file, `genome_2bit` an existing 2bit file and `frequencies` the table from `compute_gc_bias`; `effective_genome_size` must be a positive integer; `output_bam` is a new path.
- **Output:** A dictionary with `status` (`completed` when the exit code is 0 and every output file exists, otherwise `failed`), `backend` (`deeptools`), `tool`, `command` (the tool name followed by its arguments), `return_code`, `outputs` (the output files that exist) and `log_tail` (the last 4000 characters of the tool output).
- **Example code:**

```python
import os
import shutil
import sys
import tempfile

os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]
work = tempfile.mkdtemp(dir=".")
os.chdir(work)
import random
import struct

import pysam


def write_twobit(path, name, seq):
    code = {"T": 0, "C": 1, "A": 2, "G": 3}
    packed = bytearray((len(seq) + 3) // 4)
    for i, base in enumerate(seq):
        packed[i // 4] |= code[base] << (6 - 2 * (i % 4))
    record = struct.pack("<IIII", len(seq), 0, 0, 0) + bytes(packed)
    index = struct.pack("B", len(name)) + name.encode() + struct.pack("<I", 16 + 1 + len(name) + 4)
    with open(path, "wb") as handle:
        handle.write(struct.pack("<IIII", 0x1A412743, 0, 1, 0) + index + record)


rng = random.Random(7)
genome = "".join(rng.choice("ACGT") for _ in range(20000))
write_twobit("genome.2bit", "chr1", genome)
header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": len(genome)}]}
starts = sorted(rng.randrange(0, 19800) for _ in range(1500))
with pysam.AlignmentFile("g.bam", "wb", header=header) as out:
    for i, start in enumerate(starts):
        read = pysam.AlignedSegment()
        read.query_name = f"r{i}"
        read.query_sequence = genome[start:start + 50]
        read.flag = 0
        read.reference_id = 0
        read.reference_start = start
        read.mapping_quality = 60
        read.cigarstring = "50M"
        read.query_qualities = pysam.qualitystring_to_array("I" * 50)
        out.write(read)
pysam.index("g.bam")

from dnakit.epigenome import compute_gc_bias, correct_gc_bias

compute_gc_bias("g.bam", "genome.2bit", 20000, "freq.txt", fragment_length=200, verbose=False)
result = correct_gc_bias("g.bam", "genome.2bit", 20000, "freq.txt", "corrected.bam", verbose=False)
print(result["status"], result["tool"], os.path.basename(result["outputs"][0]))
print(os.path.getsize("corrected.bam") > 0)
```

- **Example output:**

```text
completed correctGCBias corrected.bam
True
```

## 22) `dnakit.epigenome.bam_mapping_statistics()`

- **Purpose:** Count the mapped and unmapped reads of an indexed BAM or CRAM file, with per-reference counts.
- **Method:** Import `deeptools.plotEnrichment` and call its `openBam(returnStats=True)`, which reads the totals from the file index (`mapped`, `unmapped`) and the per-reference counts from `get_index_statistics()`. The alignment handle is closed before the function returns. Because deepTools calls `sys.exit` for a missing file, a missing index or a non-alignment file, those exits are converted to a `ConfigurationError`.
- **API:** `dnakit.epigenome.bam_mapping_statistics(bam_path: 'str | Path', *, threads: 'int' = 1, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package deepTools (extra `dnakit[bio-deeptools]`); it is imported only when the function is called, and no external program is run. The example also uses pysam, which deepTools installs.
- **Input:** `bam_path` is an indexed BAM or CRAM file; `threads` is an integer from 1 to 256 and is used only for CRAM, where deepTools counts the reads itself.
- **Output:** `status`, `operation`, `bam_path`, `format` (`BAM` or `CRAM`), `mapped_reads`, `unmapped_reads`, `references` (a list of `name`, `mapped` and `unmapped` per reference in header order) and `provenance` (package and version, and a description of the method).
- **Example code:**

```python
import os
import shutil
import tempfile

import pysam

work = tempfile.mkdtemp(dir=".")
shutil.copy("reads.sam", work)
os.chdir(work)
pysam.view("-b", "-o", "a.bam", "reads.sam", catch_stdout=False)
pysam.index("a.bam")

from dnakit.epigenome import bam_mapping_statistics

result = bam_mapping_statistics("a.bam", verbose=False)
print(result["status"], result["format"], result["mapped_reads"], result["unmapped_reads"])
print(result["references"])
```

- **Example output:**

```text
completed BAM 48 0
[{'name': 'chr1', 'mapped': 32, 'unmapped': 0}, {'name': 'chr2', 'mapped': 16, 'unmapped': 0}]
```
