# parasail kernels by name parts, and SSW with a reusable profile

Written for: DNAKit users who want any member of parasail's alignment family without a wrapper per member.

This page records the review of the 2,443 `family-variant` ledger rows of parasail. parasail ships about 2,400 alignment functions whose names are built from the same parts: `<algorithm>[_<output>][_<method>][_profile][_<width>]`. The earlier wrappers return a traceback only (`sg_*_trace`, `nw_banded`, `ssw`), so the variants that return a score, match statistics, the whole score matrix or the last row and column had no route. Two functions now cover them: `parasail_dp_variant()` builds the parasail function name from keyword arguments, and `parasail_ssw_profile_align()` builds the SSW query profile once for many targets. Both are exported from `dnakit.alignment`.

## Review result

| Rows | Result |
| --- | --- |
| 2,051 | stay `family-variant`: `parasail_dp_variant(algorithm, output, method, width, profile)` builds exactly the name of the row, and the score equals the direct upstream call. All 2,051 names were rebuilt and called in one scripted pass. |
| 390 | `broken-upstream`: a `stats` kernel with a profile input (`*_stats_*_profile_*`) segfaults in the installed parasail, even with the matching `profile_create_stats_*` profile. `parasail_dp_variant()` refuses `profile=True` with stats outputs. |
| 2 | `wrapped`: `ssw_init` and `ssw_profile` through `parasail_ssw_profile_align()`. parasail does not check `score_size`; a value outside 0 to 2 crashes it, so the wrapper refuses such values. |

## 1) `dnakit.alignment.parasail_dp_variant()`

- **Purpose：** Run any parasail dynamic-programming kernel and get what that kernel returns: score only, match statistics, the score matrix, the last row and column, or the traceback.
- **Method：** The function name is `algorithm` + output suffix + `_method` + `_profile` + `_width`, for example `sw_stats_striped_16`. `algorithm` is `nw` (global), `sw` (local), `sg` or `sg_<ends>` (semi-global with free ends `qb`, `qe`, `qx`, `db`, `de`, `dx`, `qb_db`, `qb_de`, `qe_db`, `qe_de`). `output` is `score`, `stats`, `table`, `rowcol`, `stats_table`, `stats_rowcol` or `trace`. `method` is `reference` (the plain kernel, no width), `scan`, `striped` or `diag`. A name that parasail does not have raises `ConfigurationError`; a saturated narrow lane raises `BackendExecutionError`. `profile=True` builds a query profile first (`scan` or `striped`, with a width). `width=None` selects the kernels without a width suffix, such as `nw_scan`.
- **API：** `dnakit.alignment.parasail_dp_variant(query, target, *, algorithm="nw", output="score", method="striped", width="sat", profile=False, match=5, mismatch=-4, gap_open=10, gap_extend=1, max_cells=25000000, verbose=True)`.
- **Input：** `query` and `target`: strings of ACGTN or DNA objects; `algorithm`, `output`, `method`, `width`, `profile`: pick the kernel; scoring as in the other parasail functions; `max_cells`: cap on DP cells. Table outputs hold every cell and are limited to 250,000 cells.
- **Output：** A dict with `parasail_function`, `score`, `end_query`, `end_target` (zero-based index of the best cell), the lengths and `provenance`. Stats outputs add `matches`, `similar`, `length`; `rowcol` adds `score_row` and `score_col` (and the stat rows and columns for `stats_rowcol`); `table` adds `score_table` (and the stat tables for `stats_table`); `trace` adds `query_aligned`, `target_aligned`, `cigar`, `query_offset` and `target_offset`.
Match statistics of a local alignment.

- **Example code：**

```python
from dnakit.alignment import parasail_dp_variant

stats = parasail_dp_variant(
    "ACGTTACGGATC", "ACGATACGATCC", algorithm="sw", output="stats", method="striped", width=16
)
print(stats["parasail_function"], stats["score"], stats["matches"], stats["similar"], stats["length"])
```

- **Example output：**

```text
sw_stats_striped_16 36 10 10 12
```

The whole score matrix of a global alignment.

- **Example code：**

```python
from dnakit.alignment import parasail_dp_variant

result = parasail_dp_variant("ACGT", "ACGA", algorithm="nw", output="table", method="scan", width=32)
print(result["parasail_function"], result["score"])
for row in result["score_table"]:
    print(row)
```

- **Example output：**

```text
nw_table_scan_32 11
[5, -5, -6, -7]
[-5, 10, 0, -1]
[-6, 0, 15, 5]
[-7, -1, 5, 11]
```

Last row and column of a semi-global alignment with a query profile.

- **Example code：**

```python
from dnakit.alignment import parasail_dp_variant

result = parasail_dp_variant(
    "ACGTTACGGATC",
    "ACGATACGATCC",
    algorithm="sg_qe_de",
    output="rowcol",
    method="striped",
    width=16,
    profile=True,
)
print(result["parasail_function"], result["score"])
print(result["score_row"])
print(result["score_col"])
```

- **Example output：**

```text
sg_qe_de_rowcol_striped_profile_16 36
[-15, -9, -3, -6, 0, 6, 12, 18, 15, 21, 36, 26]
[-15, -9, -3, -6, 0, 6, 12, 18, 15, 14, 20, 26]
```

## 2) `dnakit.alignment.parasail_ssw_profile_align()`

- **Purpose：** Align one query to many targets with Striped Smith-Waterman, building the query profile once.
- **Method：** Calls parasail `ssw_init` once and `ssw_profile` for each target. `score_size` (0 for 8-bit, 1 for 16-bit, 2 for both) is validated first, because parasail crashes on other values. The hit fields match `parasail_ssw_align()`; the tests compare them target by target.
- **API：** `dnakit.alignment.parasail_ssw_profile_align(query, targets, *, score_size=2, match=5, mismatch=-4, gap_open=10, gap_extend=1, max_cells=25000000, verbose=True)`.
- **Input：** `query`: ACGTN string or DNA object; `targets`: list of 1 to 1,000 sequences; `score_size`: 0, 1 or 2; scoring and `max_cells` as above.
- **Output：** A dict with `score_size`, `query_length`, `target_count`, `hits` (in the order of `targets`; each with `target_index`, `score`, `found`, spans, `cigar`, counts and `identity`) and `provenance`.
- **Example code：**

```python
from dnakit.alignment import parasail_ssw_profile_align

result = parasail_ssw_profile_align("ACGTTGCAAGTCCATG", ["ACGTTGCAAGGGTCCATG", "TTTTTTTT"])
for hit in result["hits"]:
    print(hit["target_index"], hit["score"], hit["found"], hit["cigar"])
```

- **Example output：**

```text
0 69 True 9=2D7=
1 10 True 2=
```

## Installing the optional dependency

`dnakit[bio-parasail]`. parasail is imported only when a function runs.

## Verification

`others/tests/unit/test_family_parasail.py` checks every algorithm and output against a direct upstream call, each method, width and profile form, hand-counted statistics and tables, and the refusals above.

## Scope

- Reached through `parasail_dp_variant()`: all `nw`, `sw` and `sg*` kernels except the stats profile kernels.
- Not reached: the 390 stats profile kernels (upstream crash), and the `bindings_v1` and `bindings_v2` mirrors (aliases of the same functions).
