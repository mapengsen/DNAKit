# Three-dimensional genome

This category contains **9** documented candidate capabilities, all of which compute real results natively on Hi-C contact matrices given as dense or sparse text, a path or an in-memory matrix (binary `.hic`/`.cool` files are not read). `annotate`, `calculate`, `compare` and `normalize` need numpy and return `status: "not_implemented"` when it is not installed; none of the functions is unconditionally `not_implemented`. Two results are explicit rules rather than statistical callers: TAD boundaries from `annotate` are insulation-score minima (a window-based rule, no p-values), and the compartment labels are the sign of the leading eigenvector of the observed/expected correlation matrix, whose orientation is arbitrary unless a reference track is supplied. The documented functions are exported through the formal DNAKit namespace; functions with a compatible JSON call plan are also exposed through the Agent/MCP registry.

## 1) aggregate

- **Purpose:** Coarsen a Hi-C contact matrix by summing `factor` x `factor` blocks of bins.

- **Method:** Bin k of the result sums the original bins k*`factor` .. k*`factor`+`factor`-1 in both dimensions (the last bin covers fewer original bins when n is not a multiple of `factor`), so the total number of contacts is conserved. The new bin size is `bin_size` x `factor`. Pure summation: no normalisation is applied.

- **API:** `dnakit.genome3d.aggregate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a square symmetric non-negative contact matrix: a list of rows, `{"matrix": ...}`, dense or sparse text, or a path (`.hic`, `.cool`, `.mcool` and `.scool` files are rejected with a hint to use the cooler wrappers); values must be finite and non-negative, the matrix symmetric up to 1e-9, with at most 3000 bins; for text input `options["format"]` is `"dense"` (default: whitespace/comma/tab-separated rows; a non-numeric header row and/or label column are skipped) or `"sparse"` (lines `bin_i bin_j value`, 0-based, each unordered pair once; optional `n_bins`). `options`: `entity` (only `"contact_matrix"`); `factor` (required integer >= 2); `bin_size` (optional, in bp); `output` (optional dense-TSV path; an existing file needs `overwrite=true`).

- **Output:** The envelope carries `status` (`"ok"`), `module` (`"genome3d"`), `operation`, `provenance` (`implementation`, `backend`), `input_summary` (a generic text-level summary: record count, text length, base counts, SHA-256), `warnings` and `result`. `result` has `factor`, `n_bins`, `original_n_bins`, `bin_size` (only when `bin_size` was given), `matrix` and `output_path` when a file was written. Matrices above 1000 bins are not returned in `result` (it carries `matrix_omitted: true`) and need `output`.

- **Example code:**

```python
import json
from dnakit.genome3d import aggregate

matrix = [[8, 4, 2, 1],
          [4, 8, 4, 2],
          [2, 4, 8, 4],
          [1, 2, 4, 8]]
result = aggregate(matrix, options={"factor": 2, "bin_size": 10000})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "genome3d",
  "operation": "aggregate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 56,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "2fa95e1f8bdd1f9e7d5a022853e4586c4b66a1aab6d3a0088f16470b309031eb"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "factor": 2,
    "n_bins": 2,
    "original_n_bins": 4,
    "bin_size": 20000,
    "matrix": [
      [
        24.0,
        9.0
      ],
      [
        9.0,
        24.0
      ]
    ]
  }
}
```

## 2) annotate

- **Purpose:** Annotate the bins of a contact matrix with TAD boundaries (insulation minima) or A/B compartment labels (`options["entity"]`).

- **Method:** `tad` (default) uses an insulation-minimum rule, not a statistical TAD caller. The insulation score is computed as in `calculate` (`window` w, default 5): S(i) = mean of `M[a][b]` for a in [i-w, i-1] and b in [i, i+w-1] and score = log2(S(i) / mean of S over the defined bins). Bin i is a boundary iff its score is below `threshold` (default 0.0, i.e. insulation below the matrix-wide average) and it is the minimum within `boundary_radius` bins (default = `window`) on both sides, with ties going to the leftmost bin. `tads` are bin intervals [`start_bin`, `end_bin`) (end exclusive). With boundaries b_1 < ... < b_m in a matrix of n bins, the interior domains are [b_k, b_k+1) (`edge: false`); with `include_edge_domains` (default true) the two edge domains are reported as well, flagged `edge: true`: [0, b_1) (omitted when b_1 = 0) and [b_m, n) (the last boundary opens the final domain). So a single boundary gives up to two edge domains, and with no boundary at all the whole matrix [0, n) is one edge domain. Edge domains are cut off by the matrix end, so their size is not comparable to interior domains; `include_edge_domains: false` returns only the interior domains (an empty list with fewer than two boundaries). No p-values or boundary strength are given. `compartment`: A/B labels from the sign of the leading eigenvector of the Pearson correlation matrix of the observed/expected contacts (expected = mean contact per diagonal over non-empty bins); A is eigenvector >= 0 and B < 0. The sign of an eigenvector is arbitrary: without `reference` the entry of largest magnitude is made positive, so which compartment is called A is not meaningful by itself; with a per-bin `reference` track (for example GC content) the sign is flipped so the eigenvector correlates positively with it. At least 3 non-empty bins are needed. Needs numpy; without it the function returns `status: "not_implemented"` with `required_dependency` `pip install numpy`.

- **API:** `dnakit.genome3d.annotate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a square symmetric non-negative contact matrix: a list of rows, `{"matrix": ...}`, dense or sparse text, or a path (`.hic`, `.cool`, `.mcool` and `.scool` files are rejected with a hint to use the cooler wrappers); values must be finite and non-negative, the matrix symmetric up to 1e-9, with at most 3000 bins; for text input `options["format"]` is `"dense"` (default: whitespace/comma/tab-separated rows; a non-numeric header row and/or label column are skipped) or `"sparse"` (lines `bin_i bin_j value`, 0-based, each unordered pair once; optional `n_bins`). `options["entity"]` is `"tad"` (default) or `"compartment"`. For `tad`: `window` (default 5), `boundary_radius` (default = window), `threshold` (default 0.0), `include_edge_domains` (default `true`), `bin_size` (optional bp). For `compartment`: `reference` (optional per-bin numeric list).

- **Output:** The envelope carries `status` (`"ok"`), `module` (`"genome3d"`), `operation`, `provenance` (`implementation`, `backend`), `input_summary` (a generic text-level summary: record count, text length, base counts, SHA-256), `warnings` and `result`. `tad`: `result` has `entity`, `window`, `boundary_radius`, `threshold`, `boundary_bins` (0-based bin indices), `tads` (`start_bin`, `end_bin` exclusive, `edge`: `true` for a domain that touches the start or end of the matrix and has a boundary on one side only, or none) and, with `bin_size`, `boundary_bp` and `bin_size`. `compartment`: `entity`, `eigenvector` (`null` for empty bins), `labels` (`A`/`B`/`null`), `eigenvalue`, `explained_variance` (leading eigenvalue / sum of eigenvalues of the correlation matrix) and `oriented_by_reference`.

- **Example code:**

```python
import json
from dnakit.genome3d import annotate

n = 15
matrix = [[(10.0 / (1 + abs(i - j))) if i // 5 == j // 5 else 1.0 for j in range(n)] for i in range(n)]   # three 5-bin domains
result = annotate(matrix, options={"entity": "tad", "window": 2, "bin_size": 10000})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "genome3d",
  "operation": "annotate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 15,
    "text_length": 1440,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "78895d0d58ea88a971d7a5bed4fda84be9cc4baee71bf82b17abc664229d7c1b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "entity": "tad",
    "window": 2,
    "boundary_radius": 2,
    "threshold": 0.0,
    "boundary_bins": [
      5,
      10
    ],
    "tads": [
      {
        "start_bin": 0,
        "end_bin": 5,
        "edge": true
      },
      {
        "start_bin": 5,
        "end_bin": 10,
        "edge": false
      },
      {
        "start_bin": 10,
        "end_bin": 15,
        "edge": true
      }
    ],
    "boundary_bp": [
      50000,
      100000
    ],
    "bin_size": 10000
  }
}
```

## 3) calculate

- **Purpose:** Calculate distance-decay curves, insulation scores or A/B compartment eigenvectors from a contact matrix (`options["entity"]`).

- **Method:** `distance_decay` (default): P(s) = mean contact over all valid bin pairs on diagonal s = |i - j| (empty bins, whose row sum is zero, are excluded); `slope` is the least-squares slope of log10 P(s) against log10 s for `s_min` <= s <= `s_max`. `insulation`: S(i) = mean of `M[a][b]` for a in [i-w, i-1] and b in [i, i+w-1] (the square on the diagonal whose corner is the start of bin i; defined for w <= i <= n-w, otherwise `null`) and score = log2(S(i) / mean of S over the defined bins); low scores mark insulating positions. `compartment`: leading eigenvector (numpy `eigh`) of the Pearson correlation matrix of the observed/expected matrix (expected = mean per diagonal); the sign gives A (>= 0) or B (< 0), and its orientation is arbitrary unless a per-bin `reference` (for example GC) is supplied, otherwise the entry of largest magnitude is made positive. Needs numpy; without it the function returns `status: "not_implemented"` with `required_dependency` `pip install numpy`.

- **API:** `dnakit.genome3d.calculate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a square symmetric non-negative contact matrix: a list of rows, `{"matrix": ...}`, dense or sparse text, or a path (`.hic`, `.cool`, `.mcool` and `.scool` files are rejected with a hint to use the cooler wrappers); values must be finite and non-negative, the matrix symmetric up to 1e-9, with at most 3000 bins; for text input `options["format"]` is `"dense"` (default: whitespace/comma/tab-separated rows; a non-numeric header row and/or label column are skipped) or `"sparse"` (lines `bin_i bin_j value`, 0-based, each unordered pair once; optional `n_bins`). `options["entity"]` is `"distance_decay"` (default), `"insulation"` or `"compartment"`. `distance_decay`: `bin_size` (optional bp, adds `distance_bp`), `s_min` (default 2), `s_max` (default n/2). `insulation`: `window` (default 5). `compartment`: `reference` (optional per-bin numeric list).

- **Output:** The envelope carries `status` (`"ok"`), `module` (`"genome3d"`), `operation`, `provenance` (`implementation`, `backend`), `input_summary` (a generic text-level summary: record count, text length, base counts, SHA-256), `warnings` and `result`. `distance_decay`: `result` has `entity`, `curve` (per diagonal `distance_bins`, `n_pairs`, `mean_contact`, and `distance_bp` with `bin_size`) and `scaling` (`slope`, `points`). `insulation`: `entity`, `window`, `mean_contact` (S(i)) and `score` per bin. `compartment`: `entity`, `eigenvector`, `labels`, `eigenvalue`, `explained_variance`, `oriented_by_reference`.

- **Example code:**

```python
import json
from dnakit.genome3d import calculate

n = 15
matrix = [[(10.0 / (1 + abs(i - j))) if i // 5 == j // 5 else 1.0 for j in range(n)] for i in range(n)]   # three 5-bin domains
result = calculate(matrix, options={"entity": "distance_decay", "bin_size": 10000})
print(json.dumps(result, indent=2))
```

- **Example output:** Long arrays and strings in the output are shortened with `...`.

```json
{
  "module": "genome3d",
  "operation": "calculate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 15,
    "text_length": 1440,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "78895d0d58ea88a971d7a5bed4fda84be9cc4baee71bf82b17abc664229d7c1b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "entity": "distance_decay",
    "curve": [
      {
        "distance_bins": 0,
        "n_pairs": 15,
        "mean_contact": 10.0,
        "distance_bp": 0
      },
      {
        "distance_bins": 1,
        "n_pairs": 14,
        "mean_contact": 4.428571428571429,
        "distance_bp": 10000
      },
      {
        "distance_bins": 2,
        "n_pairs": 13,
        "mean_contact": 2.6153846153846154,
        "distance_bp": 20000
      },
      {
        "distance_bins": 3,
        "n_pairs": 12,
        "mean_contact": 1.75,
        "distance_bp": 30000
      },
      "..."
    ],
    "scaling": {
      "slope": -0.8199134580718668,
      "points": 6
    }
  }
}
```

## 4) compare

- **Purpose:** Compare two contact matrices by the correlation of their upper-triangle contacts.

- **Method:** Contacts of bin pairs with `min_distance` <= j - i <= `max_distance` are collected from both matrices (pairs where either matrix has an empty bin are skipped) and compared by Pearson correlation, Spearman correlation (Pearson on average ranks) and the mean absolute difference. With `distance_normalize` each entry is first divided by the mean of its diagonal (an observed/expected transform). This is a plain correlation of contacts; it is not stratum-adjusted (HiCRep) and gives no significance. Needs numpy; without it the function returns `status: "not_implemented"` with `required_dependency` `pip install numpy`.

- **API:** `dnakit.genome3d.compare(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: `{"a": matrix, "b": matrix}` or a list of two matrices of the same size, each a square symmetric non-negative contact matrix: a list of rows, `{"matrix": ...}`, dense or sparse text, or a path (`.hic`, `.cool`, `.mcool` and `.scool` files are rejected with a hint to use the cooler wrappers); values must be finite and non-negative, the matrix symmetric up to 1e-9, with at most 3000 bins. `options`: `entity` (only `"contact_matrix"`); `min_distance` (default 1); `max_distance` (optional); `distance_normalize` (default `false`). At least two comparable bin pairs are required.

- **Output:** The envelope carries `status` (`"ok"`), `module` (`"genome3d"`), `operation`, `provenance` (`implementation`, `backend`), `input_summary` (a generic text-level summary: record count, text length, base counts, SHA-256), `warnings` and `result`. `result` has `pairs` (number of compared bin pairs), `pearson`, `spearman` (`null` if one side has zero variance) and `mean_absolute_difference`.

- **Example code:**

```python
import json
from dnakit.genome3d import compare

matrix = [[8, 4, 2, 1],
          [4, 8, 4, 2],
          [2, 4, 8, 4],
          [1, 2, 4, 8]]
other = [[8, 3, 2, 1],
         [3, 9, 5, 2],
         [2, 5, 8, 3],
         [1, 2, 3, 7]]
result = compare({"a": matrix, "b": other}, options={"min_distance": 1})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "genome3d",
  "operation": "compare",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 126,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e1188dc81e4040ecb0c3c3932535a69005cf8552399b86faff6d9b4f1f159e6c"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "pairs": 6,
    "pearson": 0.8443565629465051,
    "spearman": 0.9534625892455922,
    "mean_absolute_difference": 0.5
  }
}
```

## 5) export

- **Purpose:** Serialise a contact matrix as dense TSV or sparse triplet text and return the text, optionally writing it to a file.

- **Method:** `dense` writes one tab-separated row per bin; `sparse` writes `i<TAB>j<TAB>value` for the non-zero entries of the upper triangle including the diagonal (0-based bins), which `read(format="sparse")` accepts. Whole numbers are written without a decimal point. `write` is the same serialisation but requires an output path.

- **API:** `dnakit.genome3d.export(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a square symmetric non-negative contact matrix: a list of rows, `{"matrix": ...}`, dense or sparse text, or a path (`.hic`, `.cool`, `.mcool` and `.scool` files are rejected with a hint to use the cooler wrappers); values must be finite and non-negative, the matrix symmetric up to 1e-9, with at most 3000 bins; for text input `options["format"]` is `"dense"` (default: whitespace/comma/tab-separated rows; a non-numeric header row and/or label column are skipped) or `"sparse"` (lines `bin_i bin_j value`, 0-based, each unordered pair once; optional `n_bins`). `options`: `out_format` (`"dense"` default | `"sparse"`); `output` (optional file path); `overwrite` (default `false`).

- **Output:** The envelope carries `status` (`"ok"`), `module` (`"genome3d"`), `operation`, `provenance` (`implementation`, `backend`), `input_summary` (a generic text-level summary: record count, text length, base counts, SHA-256), `warnings` and `result`. `result` has `format`, `text`, `output_path` (the written file or `null`), `n_bins`, `total_contacts` (sum over the upper triangle including the diagonal) and `nonzero_entries` (non-zero cells of the full matrix).

- **Example code:**

```python
import json
from dnakit.genome3d import export

matrix = [[8, 4, 2, 1],
          [4, 8, 4, 2],
          [2, 4, 8, 4],
          [1, 2, 4, 8]]
result = export(matrix, options={"out_format": "sparse"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "genome3d",
  "operation": "export",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 56,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "2fa95e1f8bdd1f9e7d5a022853e4586c4b66a1aab6d3a0088f16470b309031eb"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "sparse",
    "text": "0\t0\t8\n0\t1\t4\n0\t2\t2\n0\t3\t1\n1\t1\t8\n1\t2\t4\n1\t3\t2\n2\t2\t8\n2\t3\t4\n3\t3\t8\n",
    "output_path": null,
    "n_bins": 4,
    "total_contacts": 49.0,
    "nonzero_entries": 16
  }
}
```

## 6) normalize

- **Purpose:** Normalise (balance) a Hi-C contact matrix by iterative correction (ICE) or vanilla coverage.

- **Method:** `ice` (default; Imakaev et al. 2012): bins with fewer than `min_nonzero` non-zero entries are masked (and `ignore_diagonals` near-diagonals are zeroed first). Each iteration divides the row sums s_i of the valid bins by their mean and divides M_ij by s_i s_j; the product of all s_i is the bias b_i, so the result is raw / (b_i b_j). Iteration stops when max |s_i - 1| < `tol` or after `max_iter` iterations, and the matrix is finally scaled so every valid row sums to `target_row_sum` (as cooler does). `coverage`: M_ij / (w_i w_j) with w_i the row sum divided by the mean row sum; it is not iterated. `converged` tells whether the ICE tolerance was reached, `max_row_sum_deviation` the largest remaining deviation of the row sums from `target_row_sum`. Needs numpy; without it the function returns `status: "not_implemented"` with `required_dependency` `pip install numpy`.

- **API:** `dnakit.genome3d.normalize(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a square symmetric non-negative contact matrix: a list of rows, `{"matrix": ...}`, dense or sparse text, or a path (`.hic`, `.cool`, `.mcool` and `.scool` files are rejected with a hint to use the cooler wrappers); values must be finite and non-negative, the matrix symmetric up to 1e-9, with at most 3000 bins; for text input `options["format"]` is `"dense"` (default: whitespace/comma/tab-separated rows; a non-numeric header row and/or label column are skipped) or `"sparse"` (lines `bin_i bin_j value`, 0-based, each unordered pair once; optional `n_bins`). `options`: `entity` (only `"contact_matrix"`); `method` (`"ice"` default | `"coverage"`); ICE: `max_iter` (default 200), `tol` (default 1e-5), `ignore_diagonals` (default 0; 1 = main diagonal only), `min_nonzero` (default 1), `target_row_sum` (default 1.0); `output` (optional dense-TSV path; an existing file needs `overwrite=true`).

- **Output:** The envelope carries `status` (`"ok"`), `module` (`"genome3d"`), `operation`, `provenance` (`implementation`, `backend`), `input_summary` (a generic text-level summary: record count, text length, base counts, SHA-256), `warnings` and `result`. `result` has `entity`, `method`, `matrix` (balanced, symmetric), `bias` (`null` for masked bins), `masked_bins` and, for `ice`, `iterations`, `converged` and `max_row_sum_deviation`; plus `output_path` when a file was written. Matrices above 1000 bins are not returned in `result` (it carries `matrix_omitted: true`) and need `output`.

- **Example code:**

```python
import json
from dnakit.genome3d import normalize

matrix = [[10, 4, 2],
          [4, 20, 6],
          [2, 6, 30]]
result = normalize(matrix, options={"method": "ice", "max_iter": 200, "tol": 1e-5})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "genome3d",
  "operation": "normalize",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 36,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e5e798d23a5d2fead0ab177eab4c59bd3ac887cf29f21712f11738c75c293e99"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "entity": "contact_matrix",
    "method": "ice",
    "bias": [
      3.719048104659864,
      5.595833606989041,
      6.34110030203879
    ],
    "masked_bins": [],
    "iterations": 26,
    "converged": true,
    "max_row_sum_deviation": 9.753559970127057e-06,
    "matrix": [
      [
        0.7229968552636007,
        0.1922044344481455,
        0.08480736642436058
      ],
      [
        0.1922044344481455,
        0.6387051401499937,
        0.16909152282572365
      ],
      [
        0.08480736642436058,
        0.16909152282572365,
        0.7460913571899457
      ]
    ]
  }
}
```

## 7) plot

- **Purpose:** Draw a contact matrix as an SVG heatmap (white to red on log10(1 + contacts)).

- **Method:** A hand-written SVG (no plotting library). Matrices larger than `max_cells` per side are first coarsened by summing blocks (the factor is the smallest integer that brings the side to at most `max_cells`). Each cell is a `<rect class="cell">` with `data-i`, `data-j` and `data-value`, coloured by log10(1 + value) relative to the maximum of the displayed grid; the root element carries `data-bins` and `data-factor`.

- **API:** `dnakit.genome3d.plot(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a square symmetric non-negative contact matrix: a list of rows, `{"matrix": ...}`, dense or sparse text, or a path (`.hic`, `.cool`, `.mcool` and `.scool` files are rejected with a hint to use the cooler wrappers); values must be finite and non-negative, the matrix symmetric up to 1e-9, with at most 3000 bins; for text input `options["format"]` is `"dense"` (default: whitespace/comma/tab-separated rows; a non-numeric header row and/or label column are skipped) or `"sparse"` (lines `bin_i bin_j value`, 0-based, each unordered pair once; optional `n_bins`). `options`: `entity` (only `"contact_matrix"`); `max_cells` (default 150, range 2-400); `output` (optional `.svg` path).

- **Output:** The envelope carries `status` (`"ok"`), `module` (`"genome3d"`), `operation`, `provenance` (`implementation`, `backend`), `input_summary` (a generic text-level summary: record count, text length, base counts, SHA-256), `warnings` and `result`. `result` has `svg` (the SVG text) and `output_path` (the written file or `null`).

- **Example code:**

```python
import json
from dnakit.genome3d import plot

matrix = [[8, 4, 2, 1],
          [4, 8, 4, 2],
          [2, 4, 8, 4],
          [1, 2, 4, 8]]
result = plot(matrix, options={"max_cells": 150})
print(json.dumps(result, indent=2))
```

- **Example output:** Long arrays and strings in the output are shortened with `...`.

```json
{
  "module": "genome3d",
  "operation": "plot",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 56,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "2fa95e1f8bdd1f9e7d5a022853e4586c4b66a1aab6d3a0088f16470b309031eb"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "svg": "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"40\" height=\"40\" data-bins=\"4\" data-factor=\"1\"><rect class=\"cell\" data-i=\"0\" data-j=...",
    "output_path": null
  }
}
```

## 8) read

- **Purpose:** Read a Hi-C contact matrix from dense or sparse text, a text file path or an in-memory matrix.

- **Method:** Text lines starting with `#` and blank lines are skipped. Dense rows are split at whitespace, commas or tabs; a non-numeric first row and a non-numeric label column are dropped. Sparse text lists `bin_i bin_j value` once per unordered pair (a repeated pair is an error) and is mirrored into a symmetric matrix. Values must be finite and non-negative and the matrix square and symmetric (up to 1e-9); otherwise an input error is raised. Binary `.hic`, `.cool`, `.mcool` and `.scool` files are not read: the call raises an error that points to the cooler wrappers in `dnakit.epigenome`.

- **API:** `dnakit.genome3d.read(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a square symmetric non-negative contact matrix: a list of rows, `{"matrix": ...}`, dense or sparse text, or a path (`.hic`, `.cool`, `.mcool` and `.scool` files are rejected with a hint to use the cooler wrappers); values must be finite and non-negative, the matrix symmetric up to 1e-9, with at most 3000 bins; for text input `options["format"]` is `"dense"` (default: whitespace/comma/tab-separated rows; a non-numeric header row and/or label column are skipped) or `"sparse"` (lines `bin_i bin_j value`, 0-based, each unordered pair once; optional `n_bins`). `options`: `format`, `n_bins` (sparse only), `output` (optional dense-TSV path; an existing file needs `overwrite=true`).

- **Output:** The envelope carries `status` (`"ok"`), `module` (`"genome3d"`), `operation`, `provenance` (`implementation`, `backend`), `input_summary` (a generic text-level summary: record count, text length, base counts, SHA-256), `warnings` and `result`. `result` has `n_bins`, `total_contacts` (sum over the upper triangle including the diagonal), `nonzero_entries` (non-zero cells of the full matrix) and `matrix`; plus `output_path` when a file was written. Matrices above 1000 bins are not returned in `result` (it carries `matrix_omitted: true`) and need `output`.

- **Example code:**

```python
import json
from dnakit.genome3d import read

sparse = "0\t0\t5\n0\t1\t2\n1\t1\t4\n1\t2\t1\n2\t2\t6\n"      # bin_i, bin_j, value (each unordered pair once)
result = read(sparse, options={"format": "sparse"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "genome3d",
  "operation": "read",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 5,
    "text_length": 30,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "d1c502b64c132564d19debb6c077eb023535f38fa72064ccb2321743579d4340"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "n_bins": 3,
    "total_contacts": 18.0,
    "nonzero_entries": 7,
    "matrix": [
      [
        5.0,
        2.0,
        0.0
      ],
      [
        2.0,
        4.0,
        1.0
      ],
      [
        0.0,
        1.0,
        6.0
      ]
    ]
  }
}
```

## 9) write

- **Purpose:** Write a contact matrix to a dense TSV or sparse triplet file.

- **Method:** The same serialisation as `export` (dense: one tab-separated row per bin; sparse: `i<TAB>j<TAB>value` for the non-zero upper triangle including the diagonal, 0-based), but the output path is required and an existing file is an error unless `overwrite` is true.

- **API:** `dnakit.genome3d.write(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a square symmetric non-negative contact matrix: a list of rows, `{"matrix": ...}`, dense or sparse text, or a path (`.hic`, `.cool`, `.mcool` and `.scool` files are rejected with a hint to use the cooler wrappers); values must be finite and non-negative, the matrix symmetric up to 1e-9, with at most 3000 bins; for text input `options["format"]` is `"dense"` (default: whitespace/comma/tab-separated rows; a non-numeric header row and/or label column are skipped) or `"sparse"` (lines `bin_i bin_j value`, 0-based, each unordered pair once; optional `n_bins`). `options`: `output` (required file path); `out_format` (`"dense"` default | `"sparse"`); `overwrite` (default `false`).

- **Output:** The envelope carries `status` (`"ok"`), `module` (`"genome3d"`), `operation`, `provenance` (`implementation`, `backend`), `input_summary` (a generic text-level summary: record count, text length, base counts, SHA-256), `warnings` and `result`. `result` has `format`, `output_path`, `bytes_written`, `n_bins`, `total_contacts` and `nonzero_entries`. The file is written as a side effect.

- **Example code:**

```python
import json
from dnakit.genome3d import write

matrix = [[8, 4, 2, 1],
          [4, 8, 4, 2],
          [2, 4, 8, 4],
          [1, 2, 4, 8]]
result = write(matrix, options={"output": "contacts.tsv", "out_format": "dense"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "genome3d",
  "operation": "write",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 56,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "2fa95e1f8bdd1f9e7d5a022853e4586c4b66a1aab6d3a0088f16470b309031eb"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "dense",
    "output_path": "contacts.tsv",
    "bytes_written": 32,
    "n_bins": 4,
    "total_contacts": 49.0,
    "nonzero_entries": 16
  }
}
```
