# parasail 比对结果渲染：三行文本视图与 SAM CIGAR

本页对应工作清单 P25 包中的 71 个 parasail 条目。其中 69 条（Needleman-Wunsch 与十二种半全局末端空位变体，各自的分数、统计、完整 DP 表、最后一行一列以及组合统计输出）正是 `parasail_dp_variant()` 已经能构造的成员：用 `algorithm` 选内核，用 `output` 选结果。另外两条，回溯的 SAM CIGAR 与三行文本视图，由新增的一个函数 `parasail_render_alignment()` 提供。测试版本为 parasail 1.3.4。

## 1) `dnakit.alignment.parasail_render_alignment()`

- **作用：** 用 parasail 比对两条 DNA 序列，并以 parasail 的三行文本视图、SAM 风格 CIGAR 或两者返回回溯结果。
- **处理方法：** 以 ACGTN 的匹配/错配矩阵调用 `parasail.<algorithm>_trace_striped_sat`，再取 `Result.get_traceback()`（查询行、比较行、数据库行）和 `Result.get_cigar()`（操作符 `=`、`X`、`I`、`D` 与起始偏移）。`algorithm` 取 `nw`、`sw`、`sg` 或任意 `sg_<ends>` 变体（`qb`、`qe`、`qx`、`db`、`de`、`dx`、`qb_db`、`qb_de`、`qe_db`、`qe_de`）；`view` 取 `text`、`cigar` 或 `both`。长度为 k 的空位罚分为 `gap_open + (k - 1) * gap_extend`，与 parasail 一致。
- **API：** `dnakit.alignment.parasail_render_alignment(query, target, *, algorithm="sw", view="both", match=5, mismatch=-4, gap_open=10, gap_extend=1, max_cells=25000000, verbose=True)`。
- **输入：** 两条 ACGTN 字符串或 DNA 对象，以及打分参数。
- **输出：** 字典，含 `score`、`end_query`、`end_target`、`parasail_function`，并按 `view` 给出 `query_line`、`comparison_line`、`target_line`、`text`、`cigar`、`query_offset`、`target_offset`，另有 `provenance`。
- **示例代码：**

```python
from dnakit.alignment import parasail_render_alignment

result = parasail_render_alignment("ACGTACGGT", "TTACGTCCGTAA", algorithm="sg_qx", verbose=False)
print(result["text"])
print(result["cigar"], result["score"])
```

- **示例输出：**

```text
--ACGTACG---GT
  ||||.||     
TTACGTCCGTAA--
2D4=1X2=3D2I 3
```

## 由 `parasail_dp_variant()` 覆盖的条目

`parasail_dp_variant(query, target, algorithm=..., output=...)`，`output` 取 `score`、`stats`、`table`、`rowcol`、`stats_table`、`stats_rowcol`、`trace`。例如 `algorithm="sg_qb_de", output="stats_rowcol"` 即 `parasail.sg_qb_de_stats_rowcol`；各行的 `_scan`、`_scan_64`、`_scan_32`、`_scan_16`、`_scan_8` 入口对应 `method="scan"` 加 `width`。详见 `109_family_parasail.md`。
