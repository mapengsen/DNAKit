# BEDOPS bedmap: statistics and listings of overlapping elements

Two functions in `dnakit.annotation` (module `bedops_p05_cli_functions`) wrap `bedmap` from [BEDOPS](https://github.com/bedops/bedops) 2.4.42, which traverses a reference BED file and, for each element, summarises the elements of a map file that overlap it. `bedmap_score_statistic` returns one number per reference row (score statistics, counts, bases) and `bedmap_echo_elements` prints the reference row or the overlapping elements, IDs, scores, sizes and ranges. The bedmap operation is chosen by a parameter, so one function covers many `bedmap --<operation>` flags. Both input files must be sorted as `sort-bed` sorts them.

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

## 1) `dnakit.annotation.bedmap_score_statistic()`

- **Purpose:** Summarise, for every reference element, the map elements that overlap it with one number: a score statistic, a count or a number of bases.
- **Method:** Run `bedmap [flags] [overlap option] --<statistic> <reference> [<map>]` and stream the one-column result to `output_path` (one row per reference row; `NAN` where nothing overlaps). `statistic` is one of `mean`, `median`, `sum`, `min`, `max`, `variance`, `stdev`, `cv` (stdev / mean), `mad` (median absolute deviation, scaled by `mad_multiplier`), `kth` (the generalized median, the `kth_fraction` quantile), `tmean` (trimmed mean ignoring the bottom `trim_low` and top `trim_high` fractions), `wmean` (mean weighted by the coverage of the reference element), and the overlap measures `count`, `bases`, `bases-uniq`, `bases-uniq-f` and `indicator`. The score statistics read column 5 of the map file. Overlap is defined by `overlap_mode` (`bp-ovr`, `range`, `exact`, `fraction-ref`, `fraction-map`, `fraction-both`, `fraction-either`) and `overlap_value`; BEDOPS defaults to `bp-ovr 1`. The function covers the work-list rows coefficient of variation, generalized median, median absolute deviation, maximum, mean, median, minimum, standard deviation, sum, trimmed mean, variance and coverage-weighted mean of overlapping scores, total overlapping bases and count of overlapping elements.
- **API:** `bedmap_score_statistic(reference_path: 'str', output_path: 'str', *, map_path: 'str | None' = None, statistic: 'str' = 'mean', kth_fraction: 'float | None' = None, mad_multiplier: 'float | None' = None, trim_low: 'float | None' = None, trim_high: 'float | None' = None, overlap_mode: 'str | None' = None, overlap_value: 'float | None' = None, delimiter: 'str | None' = None, precision: 'int | None' = None, scientific: 'bool' = False, skip_unmapped: 'bool' = False, unmapped_value: 'str | None' = None, chromosome: 'str | None' = None, error_check: 'bool' = False, header: 'bool' = False, bedops_bin_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the BEDOPS program `bedmap` (2.4.42 is the tested version) on PATH, or `bedops_bin_dir` pointing to the BEDOPS `bin` folder; no Python package is required. The example output below was produced with BEDOPS 2.4.42 installed. A missing program raises `BackendUnavailableError` with an install hint.
- **Input:** `reference_path` is an existing sorted BED file; `map_path` (optional, default the reference file itself) is an existing sorted BED file with at least 5 columns for score statistics. `output_path` must not exist yet. `statistic` is one of the values above; `kth` needs `kth_fraction` in (0, 1], `tmean` needs `trim_low` and `trim_high` in [0, 1] with sum at most 1, `mad` accepts a positive `mad_multiplier`, and these parameters are refused for other statistics. `overlap_mode` / `overlap_value`: an integer for `bp-ovr` (at least 1) and `range` (at least 0), a number in (0, 1] for the fraction modes, none for `exact`. `delimiter`, `precision`, `scientific`, `skip_unmapped`, `unmapped_value`, `chromosome`, `error_check` and `header` map to `--delim`, `--prec`, `--sci`, `--skip-unmapped`, `--unmapped-val`, `--chrom`, `--ec` and `--header`.
- **Output:** The common result dictionary plus `statistic`, `values` (the first 10 rows parsed to numbers; `NAN` stays text), `line_count` and `preview`. `outputs` holds `output`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.annotation.bedops_p05_cli_functions import bedmap_score_statistic

with tempfile.TemporaryDirectory() as tmp:
    for statistic in ("count", "bases", "sum", "mean", "median", "min", "max"):
        result = bedmap_score_statistic("ref.bed", str(Path(tmp) / f"{statistic}.txt"), map_path="map.bed", statistic=statistic, verbose=False)
        print(statistic, result["values"])
    trimmed = bedmap_score_statistic(
        "ref.bed", str(Path(tmp) / "tmean.txt"), map_path="map.bed", statistic="tmean", trim_low=0.0, trim_high=0.0, verbose=False
    )
    print("tmean", trimmed["values"])
    half = bedmap_score_statistic(
        "ref.bed", str(Path(tmp) / "half.txt"), map_path="map.bed", statistic="count",
        overlap_mode="fraction-ref", overlap_value=0.5, verbose=False,
    )
    print("count with fraction-ref 0.5", half["values"])
```

- **Example output:**

```text
count [2, 2, 0, 1]
bases [70, 60, 0, 10]
sum [6.0, 10.0, 'NAN', 10.0]
mean [3.0, 5.0, 'NAN', 10.0]
median [3.0, 5.0, 'NAN', 10.0]
min [2.0, 4.0, 'NAN', 10.0]
max [4.0, 6.0, 'NAN', 10.0]
tmean [3.0, 5.0, 'NAN', 10.0]
count with fraction-ref 0.5 [1, 1, 0, 0]
```

## 2) `dnakit.annotation.bedmap_echo_elements()`

- **Purpose:** Print, for every reference element, the reference row itself or the map elements, IDs, scores, sizes or ranges that overlap it.
- **Method:** Run `bedmap [flags] [overlap option] --<field> <reference> [<map>]` and stream the result to `output_path`. `field` is one of `echo` (the reference row), `echo-map` (all overlapping map elements), `echo-map-id` and `echo-map-id-uniq` (their IDs from column 4, all or unique), `echo-map-range` (the coordinate range they span), `echo-map-score` (column 5 scores), `echo-map-size` (their full lengths), `echo-overlap-size` (the lengths of the overlaps), `echo-ref-name` (`chrom:start-end`), `echo-ref-row-id` (`id-<line number>`), `echo-ref-size`, and the highest or lowest scoring overlapping element `max-element` and `min-element` (`max-element-rand` and `min-element-rand` choose randomly among ties). Several values in one column are joined with `multi_delimiter` (BEDOPS default `;`). The function covers the work-list rows highest-scoring and lowest-scoring overlapping element, echo reference elements, list overlapping map elements, IDs and unique IDs of overlapping map elements, coordinate range, scores and lengths of overlapping map elements, lengths of overlaps, reference elements as chrom:start-end names and reference row identifiers.
- **API:** `bedmap_echo_elements(reference_path: 'str', output_path: 'str', *, map_path: 'str | None' = None, field: 'str' = 'echo-map', overlap_mode: 'str | None' = None, overlap_value: 'float | None' = None, delimiter: 'str | None' = None, multi_delimiter: 'str | None' = None, precision: 'int | None' = None, scientific: 'bool' = False, skip_unmapped: 'bool' = False, unmapped_value: 'str | None' = None, chromosome: 'str | None' = None, error_check: 'bool' = False, header: 'bool' = False, bedops_bin_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the BEDOPS program `bedmap` (2.4.42 is the tested version) on PATH, or `bedops_bin_dir` pointing to the BEDOPS `bin` folder; no Python package is required. The example output below was produced with BEDOPS 2.4.42 installed. A missing program raises `BackendUnavailableError` with an install hint.
- **Input:** `reference_path`, `map_path` and `output_path` are as in `bedmap_score_statistic` (ID fields need at least 4 map columns, score fields 5). `field` is one of the values above. `overlap_mode`, `overlap_value`, `delimiter`, `multi_delimiter`, `precision`, `scientific`, `skip_unmapped`, `unmapped_value`, `chromosome`, `error_check` and `header` as in `bedmap_score_statistic`.
- **Output:** The common result dictionary plus `field`, `line_count` and `preview`. `outputs` holds `output`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.annotation.bedops_p05_cli_functions import bedmap_echo_elements

with tempfile.TemporaryDirectory() as tmp:
    for field in ("echo-map-id", "echo-map-score", "echo-map-range", "echo-overlap-size", "echo-ref-name", "echo-ref-row-id"):
        result = bedmap_echo_elements("ref.bed", str(Path(tmp) / f"{field}.txt"), map_path="map.bed", field=field, verbose=False)
        print(field, result["preview"])
    best = bedmap_echo_elements(
        "ref.bed", str(Path(tmp) / "best.txt"), map_path="map.bed", field="max-element", skip_unmapped=True, verbose=False
    )
    print("max-element", [line.replace("\t", " ") for line in best["preview"]])
```

- **Example output:**

```text
echo-map-id ['m1;m2', 'm2;m3', '', 'm4']
echo-map-score ['2.000000;4.000000', '4.000000;6.000000', '', '10.000000']
echo-map-range ['chr1\t150\t350', 'chr1\t180\t350', '', 'chr2\t0\t20']
echo-overlap-size ['50;20', '50;10', '', '10']
echo-ref-name ['chr1:100-200', 'chr1:300-400', 'chr1:1000-1100', 'chr2:10-60']
echo-ref-row-id ['id-1', 'id-2', 'id-3', 'id-4']
max-element ['chr1 180 350 m2 4.000000 -', 'chr1 320 330 m3 6.000000 +', 'chr2 0 20 m4 10.000000 +']
```
