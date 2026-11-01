# BEDOPS interval operations: set operations, bedextract, closest-features, sort-bed

This group of functions in `dnakit.annotation` (module `bedops_p05_cli_functions`) wraps programs of the external toolkit [BEDOPS](https://github.com/bedops/bedops) 2.4.42 that work on sorted BED files: `bedops` (four set operations in one function), `bedextract`, `closest-features` and `sort-bed --unique/--duplicates`. Each function validates its arguments, runs one program with a fixed argument list and no shell, with a timeout and an output size limit, streams the program's standard output to an explicit new file (nothing is overwritten) and returns a dictionary. BEDOPS needs sorted input (`sort-bed` order); the functions do not sort for you, and `error_check=True` makes BEDOPS verify the order.

BEDOPS is an external program: install it (for example `conda install -c bioconda bedops`, or the release archive from https://github.com/bedops/bedops/releases), then put its `bin` folder on PATH or pass `bedops_bin_dir`. If a program cannot be found, `BackendUnavailableError` is raised with the install hint; if it exits with an error, `BackendExecutionError` (code `BEDOPS_EXECUTION_FAILED`) is raised with the tail of its log and no partial output is left behind. BEDOPS is released under the GPL-2.0 licence and is not bundled with DNAKit.

All results share `status` (`completed`), `backend` (`bedops`), `operation`, `program`, `version` (read from the program), `command`, `outputs` (paths of the files written), `log_tail` and `provenance` (tool, version, tested version 2.4.42, executable, citation); the sections below list only the extra fields.

## Example data

The examples read these files from the working directory and write their outputs to a temporary directory.

**Example file:** `ref.bed`

```text
chr1	100	200	r1
chr1	300	400	r2
chr1	1000	1100	r3
chr2	10	60	r4
```

**Example file:** `map.bed`

```text
chr1	150	250	m1	2	+
chr1	180	350	m2	4	-
chr1	320	330	m3	6	+
chr2	0	20	m4	10	+
```

**Example file:** `messy.bed`

```text
chr2	5	9	x
chr1	1	2	y
chr1	1	2	y
chr1	7	8	z
```

## 1) `dnakit.annotation.bedops_set_operation()`

- **Purpose:** Select or combine BED elements with the BEDOPS set operations `--element-of`, `--not-element-of`, `--symmdiff` and `--everything`.
- **Method:** Run `bedops [--chrom C] [--range R] [--ec] [--header] --<operation> [overlap] <files>` and stream standard output to `output_path`. `element-of` keeps the reference (first file) elements that overlap an element of the other files, `not-element-of` keeps those that do not, `symmdiff` returns the symmetric difference of the interval sets and `everything` prints the union of all elements with every column kept. `overlap` is a number of bp (`1`) or a percentage of the reference element (`"50%"`); BEDOPS uses 100% when it is omitted. `bedextract <query> <target>` gives the same rows as `element-of` with `overlap=1` when the query has no fully nested elements. The function covers the work-list rows reference elements overlapping query elements, reference elements not overlapping query elements, symmetric difference of interval sets and union of interval sets keeping all columns.
- **API:** `bedops_set_operation(input_paths: 'Sequence[str]', output_path: 'str', *, operation: 'str' = 'element-of', overlap: 'str | int | None' = None, range: 'str | int | None' = None, chromosome: 'str | None' = None, error_check: 'bool' = False, header: 'bool' = False, bedops_bin_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the BEDOPS program `bedops` (2.4.42 is the tested version) on PATH, or `bedops_bin_dir` pointing to the BEDOPS `bin` folder; no Python package is required. The example output below was produced with BEDOPS 2.4.42 installed. A missing program raises `BackendUnavailableError` with an install hint.
- **Input:** `input_paths` is a list of existing sorted BED (or Starch) files; the first one is the reference. `element-of`, `not-element-of` and `symmdiff` need at least two files, `everything` at least one. `output_path` must not exist yet. `operation` is one of `element-of`, `not-element-of`, `symmdiff`, `everything`. `overlap` (only for `element-of` and `not-element-of`) is an integer or a string such as `"50%"`. `range` pads the non-reference intervals (`"10"` or `"-5:20"`). `chromosome` restricts the run to one chromosome; `error_check` adds `--ec`; `header` adds `--header`.
- **Output:** The common result dictionary plus `line_count` and `preview` (the first 10 lines of the output). `outputs` holds `output`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.annotation.bedops_p05_cli_functions import bedops_set_operation

with tempfile.TemporaryDirectory() as tmp:
    for operation in ("element-of", "not-element-of", "symmdiff", "everything"):
        out = str(Path(tmp) / f"{operation}.bed")
        result = bedops_set_operation(
            ["ref.bed", "map.bed"], out, operation=operation, overlap=1 if "element" in operation else None, verbose=False
        )
        print(operation, result["line_count"], result["preview"][0].replace("\t", " "))
    result = bedops_set_operation(["ref.bed", "map.bed"], str(Path(tmp) / "half.bed"), overlap="50%", verbose=False)
    print("element-of 50%:", [line.split("\t")[3] for line in result["preview"]])
```

- **Example output:**

```text
element-of 3 chr1 100 200 r1
not-element-of 1 chr1 1000 1100 r3
symmdiff 6 chr1 100 150
everything 8 chr1 100 200 r1
element-of 50%: ['r1', 'r2']
```

## 2) `dnakit.annotation.bedextract_list_chromosomes()`

- **Purpose:** List the chromosome names of a sorted BED or Starch file.
- **Method:** Run `bedextract --list-chr <input>` and read its standard output.
- **API:** `bedextract_list_chromosomes(input_path: 'str', *, bedops_bin_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the BEDOPS program `bedextract` (2.4.42 is the tested version) on PATH, or `bedops_bin_dir` pointing to the BEDOPS `bin` folder; no Python package is required. The example output below was produced with BEDOPS 2.4.42 installed. A missing program raises `BackendUnavailableError` with an install hint.
- **Input:** `input_path` is an existing sorted BED or Starch file.
- **Output:** The common result dictionary plus `chromosomes` (list of names in file order) and `chromosome_count`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.annotation.bedops_p05_cli_functions import bedextract_list_chromosomes

result = bedextract_list_chromosomes("ref.bed", verbose=False)
print(result["chromosomes"], result["chromosome_count"])
```

- **Example output:**

```text
['chr1', 'chr2'] 2
```

## 3) `dnakit.annotation.bedextract_chromosome()`

- **Purpose:** Write the rows of one chromosome from a sorted BED or Starch file.
- **Method:** Run `bedextract <chromosome> <input>` and stream the rows to `output_path`; a chromosome that is absent gives an empty file.
- **API:** `bedextract_chromosome(input_path: 'str', chromosome: 'str', output_path: 'str', *, bedops_bin_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the BEDOPS program `bedextract` (2.4.42 is the tested version) on PATH, or `bedops_bin_dir` pointing to the BEDOPS `bin` folder; no Python package is required. The example output below was produced with BEDOPS 2.4.42 installed. A missing program raises `BackendUnavailableError` with an install hint.
- **Input:** `input_path` is an existing sorted BED or Starch file. `chromosome` is a non-empty name that does not start with `-`. `output_path` must not exist yet.
- **Output:** The common result dictionary plus `chromosome`, `line_count` and `preview`. `outputs` holds `output`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.annotation.bedops_p05_cli_functions import bedextract_chromosome

with tempfile.TemporaryDirectory() as tmp:
    result = bedextract_chromosome("ref.bed", "chr2", str(Path(tmp) / "chr2.bed"), verbose=False)
    print(result["chromosome"], result["line_count"], result["preview"][0].replace("\t", " "))
    empty = bedextract_chromosome("ref.bed", "chrZ", str(Path(tmp) / "chrZ.bed"), verbose=False)
    print("absent chromosome rows:", empty["line_count"])
```

- **Example output:**

```text
chr2 1 chr2 10 60 r4
absent chromosome rows: 0
```

## 4) `dnakit.annotation.bedops_closest_features()`

- **Purpose:** For each element of a BED file, report the nearest elements of a query file on its left and right.
- **Method:** Run `closest-features [flags] <input> <query>` with an explicit `--delim` (default `|`) and stream the rows to `output_path`. By default each row is the input element, the nearest left element and the nearest right element; `NA` marks a missing neighbour. An overlapping query element wins over non-overlapping ones. `closest` keeps only the closer neighbour (ties go left), `distances` appends the signed distances (an overlap has distance 0), `no_overlaps` ignores overlapping elements, `no_reference` and `no_query` drop the echoed input or query columns.
- **API:** `bedops_closest_features(input_path: 'str', query_path: 'str', output_path: 'str', *, closest: 'bool' = False, distances: 'bool' = False, no_overlaps: 'bool' = False, no_reference: 'bool' = False, no_query: 'bool' = False, delimiter: 'str' = '|', chromosome: 'str | None' = None, error_check: 'bool' = False, header: 'bool' = False, bedops_bin_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the BEDOPS program `closest-features` (2.4.42 is the tested version) on PATH, or `bedops_bin_dir` pointing to the BEDOPS `bin` folder; no Python package is required. The example output below was produced with BEDOPS 2.4.42 installed. A missing program raises `BackendUnavailableError` with an install hint.
- **Input:** `input_path` and `query_path` are existing sorted BED (or Starch) files. `output_path` must not exist yet. `closest`, `distances`, `no_overlaps`, `no_reference`, `no_query`, `error_check` and `header` are bool. `delimiter` is a non-empty string that does not start with `-` (use a real tab character for a tab). `chromosome` restricts the run to one chromosome.
- **Output:** The common result dictionary plus `line_count` and `preview`. `outputs` holds `output`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.annotation.bedops_p05_cli_functions import bedops_closest_features

with tempfile.TemporaryDirectory() as tmp:
    result = bedops_closest_features("ref.bed", "map.bed", str(Path(tmp) / "closest.txt"), distances=True, verbose=False)
    for line in result["preview"]:
        print(line.replace("\t", " "))
    nearest = bedops_closest_features(
        "ref.bed", "map.bed", str(Path(tmp) / "nearest.txt"), closest=True, no_reference=True, delimiter=";", verbose=False
    )
    print(nearest["preview"][2].replace("\t", " "))
```

- **Example output:**

```text
chr1 100 200 r1|NA|NA|chr1 180 350 m2 4 -|0
chr1 300 400 r2|chr1 180 350 m2 4 -|0|NA|NA
chr1 1000 1100 r3|chr1 180 350 m2 4 -|-651|NA|NA
chr2 10 60 r4|chr2 0 20 m4 10 +|0|NA|NA
chr1 180 350 m2 4 -
```

## 5) `dnakit.annotation.bedops_sort_bed_filter()`

- **Purpose:** Sort BED files and keep only the unique or only the repeated elements.
- **Method:** Run `sort-bed --unique` (like `sort -u`) or `sort-bed --duplicates` (like `uniq -d`) over the input files and stream the sorted result to `output_path`. `max_memory` and `temp_dir` are passed as `--max-mem` and `--tmpdir` for large inputs.
- **API:** `bedops_sort_bed_filter(input_paths: 'Sequence[str]', output_path: 'str', *, mode: 'str' = 'unique', max_memory: 'str | None' = None, temp_dir: 'str | None' = None, bedops_bin_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the BEDOPS program `sort-bed` (2.4.42 is the tested version) on PATH, or `bedops_bin_dir` pointing to the BEDOPS `bin` folder; no Python package is required. The example output below was produced with BEDOPS 2.4.42 installed. A missing program raises `BackendUnavailableError` with an install hint.
- **Input:** `input_paths` is a non-empty list of existing BED files (they need not be sorted). `output_path` must not exist yet. `mode` is `unique` or `duplicates`. `max_memory` looks like `8G`, `8000M` or `8000000000`. `temp_dir` is an existing folder.
- **Output:** The common result dictionary plus `line_count` and `preview`. `outputs` holds `output`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.annotation.bedops_p05_cli_functions import bedops_sort_bed_filter

with tempfile.TemporaryDirectory() as tmp:
    for mode in ("unique", "duplicates"):
        result = bedops_sort_bed_filter(["messy.bed"], str(Path(tmp) / f"{mode}.bed"), mode=mode, verbose=False)
        print(mode, [line.split("\t")[3] for line in result["preview"]])
```

- **Example output:**

```text
unique ['y', 'z', 'x']
duplicates ['y']
```
