# deepTools: matrix operations and plots

This group of functions in `dnakit.epigenome` wraps [deepTools](https://github.com/deeptools/deepTools) (verified here with deepTools 4.0.0). Nine functions are facades over the `computeMatrixOperations` command, which edits the gzipped matrix that `computeMatrix` writes: they show its layout, report its value range, keep or reorder rows, keep or rename groups and samples, and join several matrices. Six functions call deepTools' Python modules directly (`bamPEFragmentSize`, `computeGCBias`, `plotHeatmap`, `plotEnrichment`, `utilities`, `mapReduce`) and the `deeptoolsintervals` package: they compute a fragment-length density, draw GC-bias, matrix-heatmap and enrichment PNG files from numbers or a matrix, read template lengths from alignments, and subtract blacklist intervals from a region. None of them runs `computeMatrix`, counts reads over a BAM and BED file, or aligns anything.

The licence file shipped with deepTools 4.0.0 is a permissive MIT-style licence, and DNAKit copies none of its code. The functions are imported with `from dnakit.epigenome import ...`.

The dependencies are optional and are imported only when a function is called. The matrix operations start the external program `computeMatrixOperations` and raise `BackendUnavailableError` if it is not on PATH. The plot and helper functions import deepTools, matplotlib, numpy, pysam and `deeptoolsintervals` and raise an error that names the extra if one is missing. Install them with `pip install "dnakit[bio-deeptools]"`.

## Example data

Two files are declared. `matrix.txt` is a tiny uncompressed copy of a `computeMatrix` output (an `@` JSON header line and four tab-separated region rows, two groups and two samples with three bins each); the matrix examples compress it to `matrix.gz` first, because deepTools reads gzip. `reads.sam` is three alignments in SAM text, and `template_lengths_from_bam()` accepts SAM as well as BAM. The other examples pass their inputs directly or write small files themselves. The matrix examples set `deeptools_path` to the `computeMatrixOperations` script next to the Python interpreter; with the program on PATH the argument can be left out. All outputs are written to a temporary directory.

**Example file:** `matrix.txt`

```text
@{"upstream":[50,50],"downstream":[50,50],"body":[50,50],"bin size":[50,50],"ref point":[null,null],"verbose":false,"bin avg type":"mean","missing data as zero":false,"min threshold":null,"max threshold":null,"scale":1,"skip zeros":false,"nan after end":false,"proc number":1,"sort regions":"keep","sort using":"mean","startLabel":"TSS","endLabel":"TES","unscaled 5 prime":[0,0],"unscaled 3 prime":[0,0],"group_labels":["up.bed","down.bed"],"group_boundaries":[0,2,4],"sample_labels":["H3K4me3","H3K27ac"],"sample_boundaries":[0,3,6]}
chr1	200	300	geneA	0	+	1.000000	2.000000	3.000000	4.000000	5.000000	6.000000
chr1	500	600	geneB	0	-	6.000000	5.000000	4.000000	3.000000	2.000000	1.000000
chr1	700	800	geneC	0	+	8.000000	7.000000	6.000000	1.000000	2.000000	3.000000
chr1	900	1000	geneD	0	-	2.000000	2.000000	2.000000	9.000000	9.000000	9.000000
```

**Example file:** `reads.sam`

```text
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:1000
r1	99	chr1	101	60	50M	=	301	250	*	*
r2	147	chr1	301	60	50M	=	101	-250	*	*
r3	0	chr1	501	60	30M5D20M	*	0	0	*	*
```

## 1) `dnakit.epigenome.bind_columns()`

- **Purpose:** Join several computeMatrix outputs side by side so that their samples become one matrix.
- **Method:** Check that there are at least two matrices and that each exists, check that the output path is new in an existing folder, and run `computeMatrixOperations cbind -m <matrices> -o <output>`.
- **API:** `bind_columns(matrices: 'Sequence[str]', output: 'str', *, deeptools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program `computeMatrixOperations` on PATH (installed by deepTools, extra `dnakit[bio-deeptools]`), or a path given as `deeptools_path`.
- **Input:** `matrices`: a list of at least two existing gzipped computeMatrix files; `output`: a new file path; `deeptools_path`: optional path of the program. Fewer than two matrices raises `ConfigurationError`.
- **Output:** A dictionary with `status` (`completed` when the program exits with 0 and the output file exists, otherwise `failed`), `backend`, `tool`, the `command` list, `return_code`, `outputs` (the files that exist) and `log_tail` (the last 4000 characters of the program output). Example note: The example splits the demo matrix into one matrix per sample with `subset_matrix()`, joins the two again, and reads the header of the result.
- **Example code:**

```python
import gzip
import json
import shutil
import sys
import tempfile
from pathlib import Path

from dnakit.epigenome import bind_columns, subset_matrix

work = Path(tempfile.mkdtemp())
matrix = work / "matrix.gz"
with open("matrix.txt", "rb") as source, gzip.open(matrix, "wb") as target:
    shutil.copyfileobj(source, target)
# computeMatrixOperations is installed next to the Python interpreter here;
# leave deeptools_path out when the program is on PATH.
tool = str(Path(sys.executable).parent / "computeMatrixOperations")


def summary(path):
    with gzip.open(path, "rt") as handle:
        lines = handle.read().splitlines()
    header = json.loads(lines[0][1:])
    names = [line.split("\t")[3] for line in lines[1:]]
    return header["group_labels"], header["sample_labels"], names

first = work / "h3k4me3.gz"
second = work / "h3k27ac.gz"
subset_matrix(str(matrix), str(first), samples=["H3K4me3"], deeptools_path=tool, verbose=False)
subset_matrix(str(matrix), str(second), samples=["H3K27ac"], deeptools_path=tool, verbose=False)
result = bind_columns(
    [str(first), str(second)], str(work / "joined.gz"), deeptools_path=tool, verbose=False
)
print(result["status"], result["tool"], result["return_code"])
print(summary(work / "joined.gz"))
```

- **Example output:**

```text
completed computeMatrixOperations 0
(['up.bed', 'down.bed'], ['H3K4me3', 'H3K27ac'], ['geneA', 'geneB', 'geneC', 'geneD'])
```

## 2) `dnakit.epigenome.bind_rows()`

- **Purpose:** Stack several computeMatrix outputs one below the other so that their regions become one matrix.
- **Method:** Check that there are at least two matrices that exist and a new output path, and run `computeMatrixOperations rbind -m <matrices> -o <output>`. Rows that carry the same group label end up in the same group of the result.
- **API:** `bind_rows(matrices: 'Sequence[str]', output: 'str', *, deeptools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program `computeMatrixOperations` on PATH (installed by deepTools, extra `dnakit[bio-deeptools]`), or a path given as `deeptools_path`.
- **Input:** `matrices`: a list of at least two existing gzipped computeMatrix files; `output`: a new file path; `deeptools_path`: optional path of the program.
- **Output:** The same dictionary as `bind_columns()`: `status`, `backend`, `tool`, `command`, `return_code`, `outputs` and `log_tail`. Example note: The example stacks the demo matrix on itself, so every region appears twice, grouped by group label.
- **Example code:**

```python
import gzip
import json
import shutil
import sys
import tempfile
from pathlib import Path

from dnakit.epigenome import bind_rows

work = Path(tempfile.mkdtemp())
matrix = work / "matrix.gz"
with open("matrix.txt", "rb") as source, gzip.open(matrix, "wb") as target:
    shutil.copyfileobj(source, target)
# computeMatrixOperations is installed next to the Python interpreter here;
# leave deeptools_path out when the program is on PATH.
tool = str(Path(sys.executable).parent / "computeMatrixOperations")


def summary(path):
    with gzip.open(path, "rt") as handle:
        lines = handle.read().splitlines()
    header = json.loads(lines[0][1:])
    names = [line.split("\t")[3] for line in lines[1:]]
    return header["group_labels"], header["sample_labels"], names

result = bind_rows(
    [str(matrix), str(matrix)], str(work / "stacked.gz"), deeptools_path=tool, verbose=False
)
print(result["status"], result["command"][1])
groups, samples, names = summary(work / "stacked.gz")
print(groups, samples)
print(names)
```

- **Example output:**

```text
completed rbind
['up.bed', 'down.bed'] ['H3K4me3', 'H3K27ac']
['geneA', 'geneB', 'geneA', 'geneB', 'geneC', 'geneD', 'geneC', 'geneD']
```

## 3) `dnakit.epigenome.data_range()`

- **Purpose:** Report the minimum, maximum, median and 10th and 90th percentile of each sample of a computeMatrix output.
- **Method:** Run `computeMatrixOperations dataRange -m <matrix>` in the folder of the matrix and return the text it prints. The text is the combined standard output and error, so an upstream warning can appear before the table; the example prints only the last three lines. A non-zero exit raises `ConfigurationError`.
- **API:** `data_range(matrix: 'str', *, deeptools_path: 'str | None' = None, verbose: bool = True) -> 'str'`.
- **Requirements:** needs the external program `computeMatrixOperations` on PATH (installed by deepTools, extra `dnakit[bio-deeptools]`), or a path given as `deeptools_path`.
- **Input:** `matrix`: an existing gzipped computeMatrix file; `deeptools_path`: optional path of the program.
- **Output:** A string: a tab-separated table with the header `Samples, Min, Max, Median, 10th, 90th` and one row per sample. Example note: The example prints the table with the tabs replaced by vertical bars.
- **Example code:**

```python
import gzip
import json
import shutil
import sys
import tempfile
from pathlib import Path

from dnakit.epigenome import data_range

work = Path(tempfile.mkdtemp())
matrix = work / "matrix.gz"
with open("matrix.txt", "rb") as source, gzip.open(matrix, "wb") as target:
    shutil.copyfileobj(source, target)
# computeMatrixOperations is installed next to the Python interpreter here;
# leave deeptools_path out when the program is on PATH.
tool = str(Path(sys.executable).parent / "computeMatrixOperations")

text = data_range(str(matrix), deeptools_path=tool, verbose=False)
for line in text.splitlines()[-3:]:
    print(" | ".join(line.split("\t")))
```

- **Example output:**

```text
Samples | Min | Max | Median | 10th | 90th
H3K4me3 | 1.0 | 8.0 | 3.5 | 2.0 | 6.9
H3K27ac | 1.0 | 9.0 | 3.5 | 1.1 | 9.0
```

## 4) `dnakit.epigenome.filter_strand()`

- **Purpose:** Keep only the regions on one strand of a computeMatrix output.
- **Method:** Run `computeMatrixOperations filterStrand -m <matrix> -o <output> -s <strand>`. `+`, `-` and `.` select regions with that strand column.
- **API:** `filter_strand(matrix: 'str', output: 'str', *, strand: "Literal['+', '-', '.']", deeptools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program `computeMatrixOperations` on PATH (installed by deepTools, extra `dnakit[bio-deeptools]`), or a path given as `deeptools_path`.
- **Input:** `matrix`: an existing gzipped computeMatrix file; `output`: a new file path; `strand`: keyword-only, `+`, `-` or `.`.
- **Output:** The same dictionary as `bind_columns()`: `status`, `backend`, `tool`, `command`, `return_code`, `outputs` and `log_tail`. Example note: The example keeps the plus-strand regions of the demo matrix (geneA and geneC).
- **Example code:**

```python
import gzip
import json
import shutil
import sys
import tempfile
from pathlib import Path

from dnakit.epigenome import filter_strand

work = Path(tempfile.mkdtemp())
matrix = work / "matrix.gz"
with open("matrix.txt", "rb") as source, gzip.open(matrix, "wb") as target:
    shutil.copyfileobj(source, target)
# computeMatrixOperations is installed next to the Python interpreter here;
# leave deeptools_path out when the program is on PATH.
tool = str(Path(sys.executable).parent / "computeMatrixOperations")


def summary(path):
    with gzip.open(path, "rt") as handle:
        lines = handle.read().splitlines()
    header = json.loads(lines[0][1:])
    names = [line.split("\t")[3] for line in lines[1:]]
    return header["group_labels"], header["sample_labels"], names

result = filter_strand(
    str(matrix), str(work / "plus.gz"), strand="+", deeptools_path=tool, verbose=False
)
print(result["status"], result["command"][-2:])
print(summary(work / "plus.gz")[2])
```

- **Example output:**

```text
completed ['-s', '+']
['geneA', 'geneC']
```

## 5) `dnakit.epigenome.filter_values()`

- **Purpose:** Keep only the rows of a computeMatrix output whose values lie within a range.
- **Method:** Run `computeMatrixOperations filterValues -m <matrix> -o <output> [--min X] [--max Y]`. At least one of `minimum` and `maximum` must be given. The argument check runs after the matrix and output paths are checked and before the program is looked up. With the numpy 2.1.3 installed here, deepTools 4.0.0 stops inside this subcommand with `module 'numpy' has no attribute 'warnings'`, so a real call returns `status` `failed`; the cause is in deepTools, not in the facade.
- **API:** `filter_values(matrix: 'str', output: 'str', *, minimum: 'float | None' = None, maximum: 'float | None' = None, deeptools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program `computeMatrixOperations` on PATH (installed by deepTools, extra `dnakit[bio-deeptools]`), or a path given as `deeptools_path`. The subcommand needs a deepTools build that works with the installed numpy; with deepTools 4.0.0 and numpy 2 it fails.
- **Input:** `matrix`: an existing gzipped computeMatrix file; `output`: a new file path; `minimum` and `maximum`: optional floats, at least one required, otherwise `ConfigurationError`.
- **Output:** The same dictionary as `bind_columns()`: `status`, `backend`, `tool`, `command`, `return_code`, `outputs` and `log_tail`. Example note: The example shows the behaviour for a call without a range, which does not depend on the program. A realistic call follows the output.
- **Example code:**

```python
from pathlib import Path
import gzip
import shutil
import tempfile

from dnakit.epigenome import filter_values
from dnakit.exceptions import ConfigurationError

work = Path(tempfile.mkdtemp())
matrix = work / "matrix.gz"
with open("matrix.txt", "rb") as source, gzip.open(matrix, "wb") as target:
    shutil.copyfileobj(source, target)
try:
    filter_values(str(matrix), str(work / "kept.gz"), verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ")[-1])
```

- **Example output:**

```text
ConfigurationError filter_values needs a minimum or a maximum.
```

Usage with the dependency installed:

```python
from dnakit.epigenome import filter_values

result = filter_values("matrix.gz", "kept.gz", minimum=3.0)
print(result["status"], result["outputs"])
```

## 6) `dnakit.epigenome.matrix_info()`

- **Purpose:** Return the group and sample layout of a computeMatrix output.
- **Method:** Run `computeMatrixOperations info -m <matrix>` in the folder of the matrix and return the text it prints. A non-zero exit raises `ConfigurationError`.
- **API:** `matrix_info(matrix: 'str', *, deeptools_path: 'str | None' = None, verbose: bool = True) -> 'str'`.
- **Requirements:** needs the external program `computeMatrixOperations` on PATH (installed by deepTools, extra `dnakit[bio-deeptools]`), or a path given as `deeptools_path`.
- **Input:** `matrix`: an existing gzipped computeMatrix file; `deeptools_path`: optional path of the program.
- **Output:** A string: `Groups:` followed by one tab-indented group label per line, then `Samples:` followed by one tab-indented sample label per line. Example note: The example splits the text on whitespace.
- **Example code:**

```python
import gzip
import json
import shutil
import sys
import tempfile
from pathlib import Path

from dnakit.epigenome import matrix_info

work = Path(tempfile.mkdtemp())
matrix = work / "matrix.gz"
with open("matrix.txt", "rb") as source, gzip.open(matrix, "wb") as target:
    shutil.copyfileobj(source, target)
# computeMatrixOperations is installed next to the Python interpreter here;
# leave deeptools_path out when the program is on PATH.
tool = str(Path(sys.executable).parent / "computeMatrixOperations")

text = matrix_info(str(matrix), deeptools_path=tool, verbose=False)
print(text.split())
```

- **Example output:**

```text
['Groups:', 'up.bed', 'down.bed', 'Samples:', 'H3K4me3', 'H3K27ac']
```

## 7) `dnakit.epigenome.relabel_matrix()`

- **Purpose:** Rename the groups and/or samples of a computeMatrix output.
- **Method:** Run `computeMatrixOperations relabel -m <matrix> -o <output>` with `--groupLabels` and/or `--sampleLabels`. At least one of the two lists must be given, otherwise `ConfigurationError` is raised before the program is looked up.
- **API:** `relabel_matrix(matrix: 'str', output: 'str', *, group_labels: 'Sequence[str] | None' = None, sample_labels: 'Sequence[str] | None' = None, deeptools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program `computeMatrixOperations` on PATH (installed by deepTools, extra `dnakit[bio-deeptools]`), or a path given as `deeptools_path`.
- **Input:** `matrix`: an existing gzipped computeMatrix file; `output`: a new file path; `group_labels` and `sample_labels`: optional lists of new labels.
- **Output:** The same dictionary as `bind_columns()`: `status`, `backend`, `tool`, `command`, `return_code`, `outputs` and `log_tail`. Example note: The example renames both groups and both samples and reads the new header.
- **Example code:**

```python
import gzip
import json
import shutil
import sys
import tempfile
from pathlib import Path

from dnakit.epigenome import relabel_matrix

work = Path(tempfile.mkdtemp())
matrix = work / "matrix.gz"
with open("matrix.txt", "rb") as source, gzip.open(matrix, "wb") as target:
    shutil.copyfileobj(source, target)
# computeMatrixOperations is installed next to the Python interpreter here;
# leave deeptools_path out when the program is on PATH.
tool = str(Path(sys.executable).parent / "computeMatrixOperations")


def summary(path):
    with gzip.open(path, "rt") as handle:
        lines = handle.read().splitlines()
    header = json.loads(lines[0][1:])
    names = [line.split("\t")[3] for line in lines[1:]]
    return header["group_labels"], header["sample_labels"], names

result = relabel_matrix(
    str(matrix),
    str(work / "renamed.gz"),
    group_labels=["upstream_genes", "downstream_genes"],
    sample_labels=["K4me3", "K27ac"],
    deeptools_path=tool,
    verbose=False,
)
print(result["status"])
groups, samples, names = summary(work / "renamed.gz")
print(groups, samples, names)
```

- **Example output:**

```text
completed
['upstream_genes', 'downstream_genes'] ['K4me3', 'K27ac'] ['geneA', 'geneB', 'geneC', 'geneD']
```

## 8) `dnakit.epigenome.sort_matrix()`

- **Purpose:** Reorder the rows of a computeMatrix output to follow the order of one or more region files.
- **Method:** Run `computeMatrixOperations sort -m <matrix> -o <output> -R <regions>`. An empty list of region files raises `ConfigurationError`, and each region file must exist. deepTools matches each region file to a matrix group by the file name, so the file names must equal the group labels of the matrix, and then orders the rows of every group as the file lists them.
- **API:** `sort_matrix(matrix: 'str', output: 'str', regions: 'Sequence[str]', *, deeptools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program `computeMatrixOperations` on PATH (installed by deepTools, extra `dnakit[bio-deeptools]`), or a path given as `deeptools_path`.
- **Input:** `matrix`: an existing gzipped computeMatrix file; `output`: a new file path; `regions`: a non-empty list of existing BED (or GTF) files named like the groups.
- **Output:** The same dictionary as `bind_columns()`: `status`, `backend`, `tool`, `command`, `return_code`, `outputs` and `log_tail`. If the file names do not match the groups the status is `failed` and `log_tail` holds the upstream message. Example note: The example writes `up.bed` and `down.bed` with the two regions of each group in reverse order and shows the new row order.
- **Example code:**

```python
import gzip
import json
import shutil
import sys
import tempfile
from pathlib import Path

from dnakit.epigenome import sort_matrix

work = Path(tempfile.mkdtemp())
matrix = work / "matrix.gz"
with open("matrix.txt", "rb") as source, gzip.open(matrix, "wb") as target:
    shutil.copyfileobj(source, target)
# computeMatrixOperations is installed next to the Python interpreter here;
# leave deeptools_path out when the program is on PATH.
tool = str(Path(sys.executable).parent / "computeMatrixOperations")


def summary(path):
    with gzip.open(path, "rt") as handle:
        lines = handle.read().splitlines()
    header = json.loads(lines[0][1:])
    names = [line.split("\t")[3] for line in lines[1:]]
    return header["group_labels"], header["sample_labels"], names

# The file names must equal the group labels of the matrix.
(work / "up.bed").write_text("chr1\t500\t600\tgeneB\t0\t-\nchr1\t200\t300\tgeneA\t0\t+\n")
(work / "down.bed").write_text("chr1\t900\t1000\tgeneD\t0\t-\nchr1\t700\t800\tgeneC\t0\t+\n")
result = sort_matrix(
    str(matrix),
    str(work / "sorted.gz"),
    [str(work / "up.bed"), str(work / "down.bed")],
    deeptools_path=tool,
    verbose=False,
)
print(result["status"])
print(summary(work / "sorted.gz")[2])
```

- **Example output:**

```text
completed
['geneB', 'geneA', 'geneD', 'geneC']
```

## 9) `dnakit.epigenome.subset_matrix()`

- **Purpose:** Keep only selected groups and/or samples of a computeMatrix output.
- **Method:** Run `computeMatrixOperations subset -m <matrix> -o <output>` with `--groups` and/or `--samples`. At least one of the two lists must be given, otherwise `ConfigurationError` is raised before the program is looked up.
- **API:** `subset_matrix(matrix: 'str', output: 'str', *, groups: 'Sequence[str] | None' = None, samples: 'Sequence[str] | None' = None, deeptools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program `computeMatrixOperations` on PATH (installed by deepTools, extra `dnakit[bio-deeptools]`), or a path given as `deeptools_path`.
- **Input:** `matrix`: an existing gzipped computeMatrix file; `output`: a new file path; `groups` and `samples`: optional lists of labels to keep.
- **Output:** The same dictionary as `bind_columns()`: `status`, `backend`, `tool`, `command`, `return_code`, `outputs` and `log_tail`. Example note: The example keeps the group `up.bed` and the sample `H3K27ac`; the group label is the file name that `computeMatrix` was given.
- **Example code:**

```python
import gzip
import json
import shutil
import sys
import tempfile
from pathlib import Path

from dnakit.epigenome import subset_matrix

work = Path(tempfile.mkdtemp())
matrix = work / "matrix.gz"
with open("matrix.txt", "rb") as source, gzip.open(matrix, "wb") as target:
    shutil.copyfileobj(source, target)
# computeMatrixOperations is installed next to the Python interpreter here;
# leave deeptools_path out when the program is on PATH.
tool = str(Path(sys.executable).parent / "computeMatrixOperations")


def summary(path):
    with gzip.open(path, "rt") as handle:
        lines = handle.read().splitlines()
    header = json.loads(lines[0][1:])
    names = [line.split("\t")[3] for line in lines[1:]]
    return header["group_labels"], header["sample_labels"], names

result = subset_matrix(
    str(matrix),
    str(work / "subset.gz"),
    groups=["up.bed"],
    samples=["H3K27ac"],
    deeptools_path=tool,
    verbose=False,
)
print(result["status"], result["command"][-4:])
print(summary(work / "subset.gz"))
```

- **Example output:**

```text
completed ['--groups', 'up.bed', '--samples', 'H3K27ac']
(['up.bed'], ['H3K27ac'], ['geneA', 'geneB'])
```

## 10) `dnakit.epigenome.fragment_length_density()`

- **Purpose:** Compute the density of fragment lengths in 100 equal bins over a range.
- **Method:** Call `deeptools.bamPEFragmentSize.getDensity(lengths, min_value, max_value)`, which makes a 100-bin histogram with `density=True` over the range. Lengths outside the range are not counted in the density. The call is made while a scratch matplotlib figure is open, which is closed again afterwards.
- **API:** `fragment_length_density(lengths: 'list[float]', min_value: 'float', max_value: 'float', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package deeptools (extra `dnakit[bio-deeptools]`, which brings matplotlib, numpy and pysam).
- **Input:** `lengths`: a list of 1 to 10 million numbers from 0 to 2e9; `min_value` and `max_value`: numbers with `max_value` greater than `min_value`.
- **Output:** A dictionary with `status`, `operation`, `fragments` (the number of lengths passed in, including those outside the range), `range`, `bins` (100), `bin_edges` (101 values), `density` (100 values that integrate to 1 over the in-range lengths) and `provenance`. Example note: One of the six lengths (900) is outside 100 to 300; `fragments` still counts it, while the densities sum to 1 over the five others (each bin is 2 wide).
- **Example code:**

```python
from dnakit.epigenome import fragment_length_density

result = fragment_length_density([100, 150, 150, 200, 250, 900], 100, 300, verbose=False)
print(result["status"], result["fragments"], result["bins"], result["range"])
print(result["bin_edges"][:3], result["density"][:3])
print(round(sum(density * 2.0 for density in result["density"]), 6))
```

- **Example output:**

```text
completed 6 100 [100.0, 300.0]
[100.0, 102.0, 104.0] [0.1, 0.0, 0.0]
1.0
```

## 11) `dnakit.epigenome.plot_enrichment_counts()`

- **Purpose:** Draw an enrichment bar plot from feature read counts that are already known.
- **Method:** Check the labels, counts and totals, then call `deeptools.plotEnrichment.plotEnrichment` with an argument namespace built from the options. Each bar is the percentage of the sample's total reads that fall in the feature. The plot is saved as a PNG at a new path (its folder is created), the figures that deepTools leaves open are closed, and the percentages are recomputed in the result.
- **API:** `plot_enrichment_counts(features: 'list[str]', sample_labels: 'list[str]', feature_counts: 'list[dict[str, int]]', total_counts: 'list[int]', output_path: 'str', *, per_sample: 'bool' = False, variable_scales: 'bool' = False, plot_title: 'str' = '', colors: 'list[str] | None' = None, alpha: 'float' = 0.9, plot_width_cm: 'float' = 20, plot_height_cm: 'float' = 20, panels_per_row: 'int' = 4, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package deeptools (extra `dnakit[bio-deeptools]`, which brings matplotlib, numpy and pysam).
- **Input:** `features`: feature names; `sample_labels`: one label per sample; `feature_counts`: one dict per sample that has a count for every feature; `total_counts`: one positive total per sample; `output_path`: a new path ending in `.png`; options for grouping (`per_sample`, `variable_scales`), title, colors, alpha, size and panels per row. A missing feature or mismatched list lengths raises `ConfigurationError`.
- **Output:** A dictionary with `path` (the absolute PNG path) and `data`: `operation`, `features`, `per_sample`, `bars_per_panel`, `samples` (per sample its total and the count and percent of each feature) and `provenance`. Example note: The example prints the file name, the PNG signature check and the recomputed percentages.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import plot_enrichment_counts

work = Path(tempfile.mkdtemp())
result = plot_enrichment_counts(
    ["exons", "introns"],
    ["sample1", "sample2"],
    [{"exons": 300, "introns": 100}, {"exons": 50, "introns": 150}],
    [1000, 500],
    str(work / "enrichment.png"),
    verbose=False,
)
print(Path(result["path"]).name, Path(result["path"]).read_bytes()[:4] == b"\x89PNG")
data = result["data"]
print(data["bars_per_panel"])
for sample in data["samples"]:
    print(sample["sample"], sample["total_counts"], sample["features"])
```

- **Example output:**

```text
enrichment.png True
2
sample1 1000 {'exons': {'count': 300, 'percent': 30.0}, 'introns': {'count': 100, 'percent': 10.0}}
sample2 500 {'exons': {'count': 50, 'percent': 10.0}, 'introns': {'count': 150, 'percent': 30.0}}
```

## 12) `dnakit.epigenome.plot_gc_bias()`

- **Purpose:** Draw deepTools' GC-bias plot from a frequency table and per-region read counts.
- **Method:** Check the two tables and call `deeptools.computeGCBias.plotGCbias` to save a PNG, then bin the regions by GC fraction with `computeGCBias.bin_by` (100 bins) and report the bins from 0.2 to 0.7 together with the log2 of the observed/expected ratio. The output path must be a new `.png`.
- **API:** `plot_gc_bias(frequencies: 'list[list[float]]', reads_per_gc: 'list[list[float]]', region_size: 'int', output_path: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package deeptools (extra `dnakit[bio-deeptools]`, which brings matplotlib, numpy and pysam).
- **Input:** `frequencies`: at least three rows `[fragments, reads, observed/expected ratio]`, ratio above 0; `reads_per_gc`: rows `[reads, GC fraction]` with GC from 0 to 1; `region_size`: an integer from 1 to 10 million; `output_path`: a new path ending in `.png`.
- **Output:** A dictionary with `path` (the absolute PNG path) and `data`: `operation`, `region_size`, `gc_bins` (GC fraction, number of regions and median reads of each bin in the window), `gc_fraction`, `log2_observed_expected`, `window` and `provenance`. Example note: The example uses 10 synthetic frequency rows and 40 synthetic regions.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import plot_gc_bias

work = Path(tempfile.mkdtemp())
frequencies = [[10 + i, 10, 1.0 + 0.1 * (i % 5)] for i in range(10)]
reads_per_gc = [[5 + i % 7, 0.2 + 0.05 * (i % 11)] for i in range(40)]
result = plot_gc_bias(frequencies, reads_per_gc, 300, str(work / "gc_bias.png"), verbose=False)
print(Path(result["path"]).name, Path(result["path"]).read_bytes()[:4] == b"\x89PNG")
data = result["data"]
print(data["region_size"], data["window"], len(data["log2_observed_expected"]))
print(data["gc_bins"][0])
print([round(value, 4) for value in data["log2_observed_expected"][:3]])
```

- **Example output:**

```text
gc_bias.png True
300 [0.2, 0.7] 10
{'gc': 0.2, 'regions': 4, 'median_reads': 7.5}
[0.0, 0.1375, 0.263]
```

## 13) `dnakit.epigenome.plot_matrix_heatmap()`

- **Purpose:** Draw a heatmap with a summary profile from a computeMatrix output.
- **Method:** Read the gzipped matrix with `deeptools.heatmapper.heatmapper().read_matrix_file()` and call `deeptools.plotHeatmap.plotMatrix` to save a PNG at a new path. When `z_min` or `z_max` is `None`, deepTools chooses them from percentiles of the matrix. The result reports the layout read from the matrix.
- **API:** `plot_matrix_heatmap(matrix_path: 'str', output_path: 'str', *, plot_title: 'str' = '', colormap: 'str' = 'binary', missing_data_color: 'str' = 'black', alpha: 'float' = 1.0, z_min: 'float | None' = None, z_max: 'float | None' = None, y_min: 'float | None' = None, y_max: 'float | None' = None, average_type: 'str' = 'median', start_label: 'str' = 'TSS', end_label: 'str' = 'TES', heatmap_height_cm: 'float' = 25, heatmap_width_cm: 'float' = 7.5, per_group: 'bool' = False, what_to_show: 'str' = 'plot, heatmap and colorbar', plot_type: 'str' = 'lines', interpolation_method: 'str' = 'auto', dpi: 'int' = 200, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package deeptools (extra `dnakit[bio-deeptools]`, which brings matplotlib, numpy and pysam).
- **Input:** `matrix_path`: an existing computeMatrix file (a file that is not a matrix raises `ConfigurationError`); `output_path`: a new path ending in `.png`; options for title, colormap, missing-data color, alpha, z and y limits (`z_max` must exceed `z_min`), `average_type` (`mean`, `median`, `sum`, `min`, `max`, `std`), start and end labels, size in cm, `per_group`, `what_to_show`, `plot_type` (`lines`, `se`, `std`, `overlapped_lines`, `fill`), `interpolation_method` and `dpi` (50 to 1200).
- **Output:** A dictionary with `path` (the absolute PNG path) and `data`: `operation`, `matrix`, `samples`, `groups`, `regions`, `group_boundaries`, `group_sizes`, `sample_boundaries`, `columns_per_sample`, `parameters` (the header of the matrix), and the colormap, z limits, average type, plot type and what-to-show used. Example note: The example draws the demo matrix at 50 dpi.
- **Example code:**

```python
import gzip
import shutil
import tempfile
from pathlib import Path

from dnakit.epigenome import plot_matrix_heatmap

work = Path(tempfile.mkdtemp())
matrix = work / "matrix.gz"
with open("matrix.txt", "rb") as source, gzip.open(matrix, "wb") as target:
    shutil.copyfileobj(source, target)
result = plot_matrix_heatmap(str(matrix), str(work / "heatmap.png"), dpi=50, verbose=False)
print(Path(result["path"]).name, Path(result["path"]).read_bytes()[:4] == b"\x89PNG")
data = result["data"]
print(data["samples"], data["groups"])
print(data["regions"], data["group_sizes"], data["columns_per_sample"])
print(data["colormap"], data["average_type"], data["z_min"], data["z_max"])
```

- **Example output:**

```text
heatmap.png True
['H3K4me3', 'H3K27ac'] ['up.bed', 'down.bed']
4 [2, 2] [3, 3]
binary median None None
```

## 14) `dnakit.epigenome.subtract_blacklist_intervals()`

- **Purpose:** Return the parts of a region that no blacklist interval covers.
- **Method:** Load the BED file into a `deeptoolsintervals.GTF` tree, find the overlapping blacklist intervals with `findOverlaps`, and compute the remaining pieces with `deeptools.mapReduce.blSubtract`. Coordinates are 0-based and half-open. A file that cannot be parsed as BED raises `ConfigurationError`.
- **API:** `subtract_blacklist_intervals(blacklist_path: 'str', chrom: 'str', start: 'int', end: 'int', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional packages deeptools and deeptoolsintervals (extra `dnakit[bio-deeptools]`).
- **Input:** `blacklist_path`: an existing BED file; `chrom`: a chromosome name; `start` and `end`: integers with `end` greater than `start`.
- **Output:** A dictionary with `status`, `operation`, `blacklist`, `chrom`, `start`, `end`, `blacklisted_overlaps` (the overlapping intervals as stored in the file), `remaining` (the pieces left), `remaining_bases` and `provenance`. Example note: The region 50 to 500 contains two blacklist intervals, so three pieces remain.
- **Example code:**

```python
from dnakit.epigenome import subtract_blacklist_intervals

with open("blacklist.bed", "w") as handle:
    handle.write("chr1\t100\t200\nchr1\t400\t450\n")
result = subtract_blacklist_intervals("blacklist.bed", "chr1", 50, 500, verbose=False)
print(result["blacklisted_overlaps"])
print(result["remaining"], result["remaining_bases"])
```

- **Example output:**

```text
[[100, 200], [400, 450]]
[[50, 100], [200, 400], [450, 500]] 300
```

## 15) `dnakit.epigenome.template_lengths_from_bam()`

- **Purpose:** List the template length of every alignment in a BAM or SAM file.
- **Method:** Open the file with `pysam.AlignmentFile` and call `deeptools.utilities.getTLen` on each alignment. A read with a TLEN field gives its absolute TLEN (the signed value when `not_abs` is true). A read with TLEN 0 gets the sum of its CIGAR M, D, = and X operations. Reading stops after `max_reads` alignments and `truncated` becomes true.
- **API:** `template_lengths_from_bam(bam_path: 'str', *, not_abs: 'bool' = False, max_reads: 'int' = 1000000, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional packages deeptools and pysam (extra `dnakit[bio-deeptools]`).
- **Input:** `bam_path`: an existing BAM or SAM file (an unreadable file raises `ConfigurationError`); `not_abs`: a boolean; `max_reads`: an integer from 1 to 100 million.
- **Output:** A dictionary with `status`, `operation`, `bam`, `not_abs`, `alignments`, `truncated`, `template_lengths` (one integer per alignment), `summary` (`min`, `median` and `max`, or `None` for an empty file) and `provenance`. Example note: The third alignment has TLEN 0 and a `30M5D20M` CIGAR, so its length is 30 + 5 + 20 = 55.
- **Example code:**

```python
from dnakit.epigenome import template_lengths_from_bam

result = template_lengths_from_bam("reads.sam", verbose=False)
print(result["alignments"], result["truncated"], result["template_lengths"], result["summary"])
signed = template_lengths_from_bam("reads.sam", not_abs=True, verbose=False)
print(signed["template_lengths"])
short = template_lengths_from_bam("reads.sam", max_reads=2, verbose=False)
print(short["alignments"], short["truncated"])
```

- **Example output:**

```text
3 False [250, 250, 55] {'min': 55, 'median': 250, 'max': 250}
[250, -250, 55]
2 True
```
