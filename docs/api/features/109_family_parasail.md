# 按名称部件选择的 parasail 内核，及可复用查询 profile 的 SSW

适用读者：想用到 parasail 比对家族任意成员、又不想为每个成员写一个包装函数的 DNAKit 用户。

本页记录对 parasail 的 2,443 个 `family-variant` 台账行的复核。parasail 约有 2,400 个比对函数，名称由同一套部件拼成：`<algorithm>[_<output>][_<method>][_profile][_<width>]`。此前的包装只返回回溯（`sg_*_trace`、`nw_banded`、`ssw`），因此只返回分数、匹配统计、整张分数矩阵或最后一行一列的变体没有入口。现在由两个函数覆盖：`parasail_dp_variant()` 用关键字参数拼出 parasail 函数名；`parasail_ssw_profile_align()` 只构建一次 SSW 查询 profile，用于多个目标序列。两者都从 `dnakit.alignment` 导出。

## 复核结果

| 行数 | 结果 |
| --- | --- |
| 2,051 | 保持 `family-variant`：`parasail_dp_variant(algorithm, output, method, width, profile)` 拼出的名称与该行完全一致，分数与上游直接调用相等。这 2,051 个名称已在一次脚本运行中全部重建并调用。 |
| 390 | `broken-upstream`：带 profile 输入的 `stats` 内核（`*_stats_*_profile_*`）在已安装的 parasail 中会段错误，即使使用配套的 `profile_create_stats_*` 也一样。`parasail_dp_variant()` 对 stats 输出拒绝 `profile=True`。 |
| 2 | `wrapped`：`ssw_init` 与 `ssw_profile`，通过 `parasail_ssw_profile_align()` 提供。parasail 不检查 `score_size`，0 到 2 以外的值会使其崩溃，所以包装函数拒绝这些值。 |

## 1) `dnakit.alignment.parasail_dp_variant()`

- **作用：** 运行任意 parasail 动态规划内核，并取得该内核的输出：仅分数、匹配统计、分数矩阵、最后一行一列，或回溯。
- **处理方法：** 函数名 = `algorithm` + 输出后缀 + `_method` + `_profile` + `_width`，例如 `sw_stats_striped_16`。`algorithm` 取 `nw`（全局）、`sw`（局部）、`sg` 或 `sg_<ends>`（半全局，自由端为 `qb`、`qe`、`qx`、`db`、`de`、`dx`、`qb_db`、`qb_de`、`qe_db`、`qe_de`）。`output` 取 `score`、`stats`、`table`、`rowcol`、`stats_table`、`stats_rowcol` 或 `trace`。`method` 取 `reference`（普通内核，无宽度）、`scan`、`striped` 或 `diag`。parasail 中不存在的名称抛出 `ConfigurationError`；窄通道饱和抛出 `BackendExecutionError`。`profile=True` 先构建查询 profile（需 `scan` 或 `striped` 并给出宽度）。`width=None` 选择没有宽度后缀的内核，如 `nw_scan`。
- **API：** `dnakit.alignment.parasail_dp_variant(query, target, *, algorithm="nw", output="score", method="striped", width="sat", profile=False, match=5, mismatch=-4, gap_open=10, gap_extend=1, max_cells=25000000, verbose=True)`。
- **输入：** `query`、`target`：ACGTN 字符串或 DNA 对象；`algorithm`、`output`、`method`、`width`、`profile`：选择内核；打分参数同其他 parasail 函数；`max_cells`：动态规划格子数上限。table 输出包含每个格子，限制在 250,000 格以内。
- **输出：** 字典，含 `parasail_function`、`score`、`end_query`、`end_target`（最佳格子的零起点下标）、两条序列长度和 `provenance`。stats 输出增加 `matches`、`similar`、`length`；`rowcol` 增加 `score_row`、`score_col`（`stats_rowcol` 还有各统计量的行和列）；`table` 增加 `score_table`（`stats_table` 还有各统计量的表）；`trace` 增加 `query_aligned`、`target_aligned`、`cigar`、`query_offset`、`target_offset`。
局部比对的匹配统计。

- **示例代码：**

```python
from dnakit.alignment import parasail_dp_variant

stats = parasail_dp_variant(
    "ACGTTACGGATC", "ACGATACGATCC", algorithm="sw", output="stats", method="striped", width=16
)
print(stats["parasail_function"], stats["score"], stats["matches"], stats["similar"], stats["length"])
```

- **示例输出：**

```text
sw_stats_striped_16 36 10 10 12
```

全局比对的整张分数矩阵。

- **示例代码：**

```python
from dnakit.alignment import parasail_dp_variant

result = parasail_dp_variant("ACGT", "ACGA", algorithm="nw", output="table", method="scan", width=32)
print(result["parasail_function"], result["score"])
for row in result["score_table"]:
    print(row)
```

- **示例输出：**

```text
nw_table_scan_32 11
[5, -5, -6, -7]
[-5, 10, 0, -1]
[-6, 0, 15, 5]
[-7, -1, 5, 11]
```

带查询 profile 的半全局比对的最后一行一列。

- **示例代码：**

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

- **示例输出：**

```text
sg_qe_de_rowcol_striped_profile_16 36
[-15, -9, -3, -6, 0, 6, 12, 18, 15, 21, 36, 26]
[-15, -9, -3, -6, 0, 6, 12, 18, 15, 14, 20, 26]
```

## 2) `dnakit.alignment.parasail_ssw_profile_align()`

- **作用：** 用 Striped Smith-Waterman 把一条查询序列与多条目标序列比对，查询 profile 只构建一次。
- **处理方法：** 先调用一次 parasail 的 `ssw_init`，再对每个目标调用 `ssw_profile`。`score_size`（0 为 8 位，1 为 16 位，2 为两者）先经校验，因为其他值会使 parasail 崩溃。每个命中的字段与 `parasail_ssw_align()` 一致，测试逐个目标对比。
- **API：** `dnakit.alignment.parasail_ssw_profile_align(query, targets, *, score_size=2, match=5, mismatch=-4, gap_open=10, gap_extend=1, max_cells=25000000, verbose=True)`。
- **输入：** `query`：ACGTN 字符串或 DNA 对象；`targets`：1 到 1,000 条序列的列表；`score_size`：0、1 或 2；打分与 `max_cells` 同上。
- **输出：** 字典，含 `score_size`、`query_length`、`target_count`、`hits`（顺序与 `targets` 相同；每项含 `target_index`、`score`、`found`、坐标、`cigar`、计数和 `identity`）和 `provenance`。
- **示例代码：**

```python
from dnakit.alignment import parasail_ssw_profile_align

result = parasail_ssw_profile_align("ACGTTGCAAGTCCATG", ["ACGTTGCAAGGGTCCATG", "TTTTTTTT"])
for hit in result["hits"]:
    print(hit["target_index"], hit["score"], hit["found"], hit["cigar"])
```

- **示例输出：**

```text
0 69 True 9=2D7=
1 10 True 2=
```

## 安装可选依赖

`dnakit[bio-parasail]`。只有调用函数时才会导入 parasail。

## 验证方式

`others/tests/unit/test_family_parasail.py` 用上游直接调用核对每种算法与输出、每种 method、宽度与 profile 形式，并用手工计数核对统计量与矩阵，同时检查上述拒绝情形。

## 当前范围

- 通过 `parasail_dp_variant()` 可到达：除 stats profile 内核外的全部 `nw`、`sw`、`sg*` 内核。
- 不可到达：390 个 stats profile 内核（上游崩溃），以及 `bindings_v1`、`bindings_v2` 镜像（同一函数的别名）。
