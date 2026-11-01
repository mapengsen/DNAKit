# parasail alignment rendering: three-line text view and SAM CIGAR

This page covers the 71 parasail rows of work-list pack P25. 69 of them (Needleman-Wunsch and the twelve semi-global end-gap variants, each as score, statistics, full DP table, last row and column, and the combined statistics outputs) are exactly the members `parasail_dp_variant()` already builds: `algorithm` picks the kernel and `output` picks the result. The other two rows, the SAM CIGAR and the three-line text view of the traceback, are provided by one new function, `parasail_render_alignment()`. parasail 1.3.4 was tested.

## 1) `dnakit.alignment.parasail_render_alignment()`

- **Purpose:** Align two DNA sequences with parasail and return the traceback as parasail's three-line text view, as a SAM-style CIGAR, or both.
- **Method:** Calls `parasail.<algorithm>_trace_striped_sat` with an ACGTN match/mismatch matrix, then `Result.get_traceback()` (query line, comparison line, database line) and `Result.get_cigar()` (operations `=`, `X`, `I`, `D` and the start offsets). `algorithm` is `nw`, `sw`, `sg` or any `sg_<ends>` variant (`qb`, `qe`, `qx`, `db`, `de`, `dx`, `qb_db`, `qb_de`, `qe_db`, `qe_de`); `view` is `text`, `cigar` or `both`. A gap of length k costs `gap_open + (k - 1) * gap_extend`, as in parasail.
- **API:** `dnakit.alignment.parasail_render_alignment(query, target, *, algorithm="sw", view="both", match=5, mismatch=-4, gap_open=10, gap_extend=1, max_cells=25000000, verbose=True)`.
- **Input:** Two ACGTN strings or DNA objects and the scoring values.
- **Output:** A dictionary with `score`, `end_query`, `end_target`, `parasail_function`, and depending on `view` the keys `query_line`, `comparison_line`, `target_line`, `text`, `cigar`, `query_offset`, `target_offset`, plus `provenance`.
- **Example code:**

```python
from dnakit.alignment import parasail_render_alignment

result = parasail_render_alignment("ACGTACGGT", "TTACGTCCGTAA", algorithm="sg_qx", verbose=False)
print(result["text"])
print(result["cigar"], result["score"])
```

- **Example output:**

```text
--ACGTACG---GT
  ||||.||     
TTACGTCCGTAA--
2D4=1X2=3D2I 3
```

## Rows covered by `parasail_dp_variant()`

`parasail_dp_variant(query, target, algorithm=..., output=...)` with `output` in `score`, `stats`, `table`, `rowcol`, `stats_table`, `stats_rowcol`, `trace`. For example `algorithm="sg_qb_de", output="stats_rowcol"` is `parasail.sg_qb_de_stats_rowcol`; the `_scan`, `_scan_64`, `_scan_32`, `_scan_16` and `_scan_8` entry points of each row are `method="scan"` with `width`. See `109_family_parasail.en.md`.
