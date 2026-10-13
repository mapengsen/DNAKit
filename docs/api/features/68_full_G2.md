# 对齐补充封装：parasail、edlib 与 pyhmmer

这一组函数补充 `dnakit.alignment` 中已有的 parasail、edlib 和 pyhmmer 适配器。parasail 的 5479 个清单条目（含 `bindings_v1`、`bindings_v2` 两个绑定镜像）按“算法 × 末端空位模式”划分为 13 个族：每族选一个代表，SIMD 方法和位宽作为参数选择变体。本页封装其中 10 个半全局末端空位族，以及带状全局打分、SSW 局部比对；edlib 的 NICE 比对文本和 pyhmmer 的 DNA 谱图比对各一个。所有函数只处理显式传入的序列或本地文件，不联网。

示例从对应模块导入。

## 示例数据

下面的文件只供 `pyhmmer_align_to_hmm()` 的示例使用，用来生成 DNA 谱图。其余示例直接使用字符串。

**示例文件：** `aligned.fa`

```text
>a
ACGTACGTAC
>b
ACGTTCGTAC
>c
ACGAACGTAC
```

## 与已有函数和上游条目的关系

| 上游包 | 清单条目 | 本页封装 | 其余条目的处理 |
| --- | --- | --- | --- |
| mappy | 3 | 无 | `revcomp` 与 `fastx_read` 重复已有的 `reverse_complement_dna`、`read_fastx`；`verbose` 是日志设置 |
| edlib | 2 | `edlib_nice_alignment` | `align` 已由 `align_dna_edlib` 等封装 |
| parasail | 5479 | 12 个代表（10 个半全局族、`nw_banded`、`ssw`） | `nw`、`sw`、`sg` 的标准族已由 `align_dna_parasail`、`score_dna_parasail` 覆盖；其余为同族变体、绑定镜像或配置函数 |
| pyhmmer | 18 | `pyhmmer_align_to_hmm` | `nhmmer` 已由 `search_dna_hmm` 覆盖；`hmmsearch` 与之重复；`phmmer`、`jackhmmer`、`hmmscan` 是蛋白质搜索；`hmmpress` 只写索引 |

## 1) `dnakit.alignment.parasail_sg_qb_align()`

- **作用：** 半全局比对，查询序列开头可以不参与比对，开头的未比对部分不扣分。
- **处理方法：** 调用 parasail 的 `sg_qb_trace` 族，代表成员是 `sg_qb_trace_scan_sat`。`method`（`scan`、`striped`、`diag`）和 `width`（8、16、32、64 或 `sat`）选择 SIMD 变体，通道饱和时报错。回溯路径包含自由末端的空位列，`cigar` 覆盖整条路径；坐标和一致性统计只计算比对核心区。
- **API：** `dnakit.alignment.parasail_sg_qb_align(query, target, *, match=5, mismatch=-4, gap_open=10, gap_extend=1, method="scan", width="sat", max_cells=25000000, verbose=True)`。
- **输入：** `query`、`target`：只含 ACGTN 的字符串或 DNA 对象；`match`、`mismatch`：碱基打分；`gap_open`、`gap_extend`：仿射空位，长度为 k 的空位扣 `gap_open + (k - 1) * gap_extend`；`method`、`width`：SIMD 变体；`max_cells`：动态规划格子数上限（默认 2500 万）。
- **输出：** 字典，含 `score`、`found`、`free_end_gaps`、`query_start`、`query_end`、`target_start`、`target_end`（零起点半开区间）、`query_aligned`、`target_aligned`、`cigar`、`matches`、`mismatches`、`insertions`、`deletions`、`alignment_columns`、`identity` 和 `provenance`。
- **示例代码：**

```python
from dnakit.alignment import parasail_sg_qb_align

result = parasail_sg_qb_align("TTTTACGTACGGA", "ACGTACGGA")
print(result["score"], result["free_end_gaps"], result["cigar"])
print(result["query_start"], result["query_end"], result["target_start"], result["target_end"])
```

- **示例输出：**

```text
45 ['query_begin'] 4I9=
4 13 0 9
```

## 2) `dnakit.alignment.parasail_sg_qe_align()`

- **作用：** 半全局比对，查询序列末尾可以不参与比对，末尾的未比对部分不扣分。
- **处理方法：** 调用 parasail 的 `sg_qe_trace` 族，代表成员是 `sg_qe_trace_scan_sat`。`method`（`scan`、`striped`、`diag`）和 `width`（8、16、32、64 或 `sat`）选择 SIMD 变体，通道饱和时报错。回溯路径包含自由末端的空位列，`cigar` 覆盖整条路径；坐标和一致性统计只计算比对核心区。
- **API：** `dnakit.alignment.parasail_sg_qe_align(query, target, *, match=5, mismatch=-4, gap_open=10, gap_extend=1, method="scan", width="sat", max_cells=25000000, verbose=True)`。
- **输入：** `query`、`target`：只含 ACGTN 的字符串或 DNA 对象；`match`、`mismatch`：碱基打分；`gap_open`、`gap_extend`：仿射空位，长度为 k 的空位扣 `gap_open + (k - 1) * gap_extend`；`method`、`width`：SIMD 变体；`max_cells`：动态规划格子数上限（默认 2500 万）。
- **输出：** 字典，含 `score`、`found`、`free_end_gaps`、`query_start`、`query_end`、`target_start`、`target_end`（零起点半开区间）、`query_aligned`、`target_aligned`、`cigar`、`matches`、`mismatches`、`insertions`、`deletions`、`alignment_columns`、`identity` 和 `provenance`。
- **示例代码：**

```python
from dnakit.alignment import parasail_sg_qe_align

result = parasail_sg_qe_align("ACGTACGGATTTT", "ACGTACGGA")
print(result["score"], result["free_end_gaps"], result["cigar"])
print(result["query_start"], result["query_end"], result["target_start"], result["target_end"])
```

- **示例输出：**

```text
45 ['query_end'] 9=4I
0 9 0 9
```

## 3) `dnakit.alignment.parasail_sg_qx_align()`

- **作用：** 半全局比对，查询序列两端都可以不参与比对，因此查询可以伸出目标序列的两端。
- **处理方法：** 调用 parasail 的 `sg_qx_trace` 族，代表成员是 `sg_qx_trace_scan_sat`。`method`（`scan`、`striped`、`diag`）和 `width`（8、16、32、64 或 `sat`）选择 SIMD 变体，通道饱和时报错。回溯路径包含自由末端的空位列，`cigar` 覆盖整条路径；坐标和一致性统计只计算比对核心区。
- **API：** `dnakit.alignment.parasail_sg_qx_align(query, target, *, match=5, mismatch=-4, gap_open=10, gap_extend=1, method="scan", width="sat", max_cells=25000000, verbose=True)`。
- **输入：** `query`、`target`：只含 ACGTN 的字符串或 DNA 对象；`match`、`mismatch`：碱基打分；`gap_open`、`gap_extend`：仿射空位，长度为 k 的空位扣 `gap_open + (k - 1) * gap_extend`；`method`、`width`：SIMD 变体；`max_cells`：动态规划格子数上限（默认 2500 万）。
- **输出：** 字典，含 `score`、`found`、`free_end_gaps`、`query_start`、`query_end`、`target_start`、`target_end`（零起点半开区间）、`query_aligned`、`target_aligned`、`cigar`、`matches`、`mismatches`、`insertions`、`deletions`、`alignment_columns`、`identity` 和 `provenance`。
- **示例代码：**

```python
from dnakit.alignment import parasail_sg_qx_align

result = parasail_sg_qx_align("TTTTACGTACGGATTTT", "ACGTACGGA")
print(result["score"], result["free_end_gaps"], result["cigar"])
print(result["query_start"], result["query_end"], result["target_start"], result["target_end"])
```

- **示例输出：**

```text
45 ['query_begin', 'query_end'] 4I9=4I
4 13 0 9
```

## 4) `dnakit.alignment.parasail_sg_db_align()`

- **作用：** 半全局比对，目标序列开头可以不参与比对。
- **处理方法：** 调用 parasail 的 `sg_db_trace` 族，代表成员是 `sg_db_trace_scan_sat`。`method`（`scan`、`striped`、`diag`）和 `width`（8、16、32、64 或 `sat`）选择 SIMD 变体，通道饱和时报错。回溯路径包含自由末端的空位列，`cigar` 覆盖整条路径；坐标和一致性统计只计算比对核心区。
- **API：** `dnakit.alignment.parasail_sg_db_align(query, target, *, match=5, mismatch=-4, gap_open=10, gap_extend=1, method="scan", width="sat", max_cells=25000000, verbose=True)`。
- **输入：** `query`、`target`：只含 ACGTN 的字符串或 DNA 对象；`match`、`mismatch`：碱基打分；`gap_open`、`gap_extend`：仿射空位，长度为 k 的空位扣 `gap_open + (k - 1) * gap_extend`；`method`、`width`：SIMD 变体；`max_cells`：动态规划格子数上限（默认 2500 万）。
- **输出：** 字典，含 `score`、`found`、`free_end_gaps`、`query_start`、`query_end`、`target_start`、`target_end`（零起点半开区间）、`query_aligned`、`target_aligned`、`cigar`、`matches`、`mismatches`、`insertions`、`deletions`、`alignment_columns`、`identity` 和 `provenance`。
- **示例代码：**

```python
from dnakit.alignment import parasail_sg_db_align

result = parasail_sg_db_align("ACGTACGGA", "TTTTACGTACGGA")
print(result["score"], result["free_end_gaps"], result["cigar"])
print(result["query_start"], result["query_end"], result["target_start"], result["target_end"])
```

- **示例输出：**

```text
45 ['database_begin'] 4D9=
0 9 4 13
```

## 5) `dnakit.alignment.parasail_sg_de_align()`

- **作用：** 半全局比对，目标序列末尾可以不参与比对。
- **处理方法：** 调用 parasail 的 `sg_de_trace` 族，代表成员是 `sg_de_trace_scan_sat`。`method`（`scan`、`striped`、`diag`）和 `width`（8、16、32、64 或 `sat`）选择 SIMD 变体，通道饱和时报错。回溯路径包含自由末端的空位列，`cigar` 覆盖整条路径；坐标和一致性统计只计算比对核心区。
- **API：** `dnakit.alignment.parasail_sg_de_align(query, target, *, match=5, mismatch=-4, gap_open=10, gap_extend=1, method="scan", width="sat", max_cells=25000000, verbose=True)`。
- **输入：** `query`、`target`：只含 ACGTN 的字符串或 DNA 对象；`match`、`mismatch`：碱基打分；`gap_open`、`gap_extend`：仿射空位，长度为 k 的空位扣 `gap_open + (k - 1) * gap_extend`；`method`、`width`：SIMD 变体；`max_cells`：动态规划格子数上限（默认 2500 万）。
- **输出：** 字典，含 `score`、`found`、`free_end_gaps`、`query_start`、`query_end`、`target_start`、`target_end`（零起点半开区间）、`query_aligned`、`target_aligned`、`cigar`、`matches`、`mismatches`、`insertions`、`deletions`、`alignment_columns`、`identity` 和 `provenance`。
- **示例代码：**

```python
from dnakit.alignment import parasail_sg_de_align

result = parasail_sg_de_align("ACGTACGGA", "ACGTACGGATTTT")
print(result["score"], result["free_end_gaps"], result["cigar"])
print(result["query_start"], result["query_end"], result["target_start"], result["target_end"])
```

- **示例输出：**

```text
45 ['database_end'] 9=4D
0 9 0 9
```

## 6) `dnakit.alignment.parasail_sg_dx_align()`

- **作用：** 半全局比对，目标序列两端都可以不参与比对，即查询完整落在目标内部（拟合比对）。
- **处理方法：** 调用 parasail 的 `sg_dx_trace` 族，代表成员是 `sg_dx_trace_scan_sat`。`method`（`scan`、`striped`、`diag`）和 `width`（8、16、32、64 或 `sat`）选择 SIMD 变体，通道饱和时报错。回溯路径包含自由末端的空位列，`cigar` 覆盖整条路径；坐标和一致性统计只计算比对核心区。
- **API：** `dnakit.alignment.parasail_sg_dx_align(query, target, *, match=5, mismatch=-4, gap_open=10, gap_extend=1, method="scan", width="sat", max_cells=25000000, verbose=True)`。
- **输入：** `query`、`target`：只含 ACGTN 的字符串或 DNA 对象；`match`、`mismatch`：碱基打分；`gap_open`、`gap_extend`：仿射空位，长度为 k 的空位扣 `gap_open + (k - 1) * gap_extend`；`method`、`width`：SIMD 变体；`max_cells`：动态规划格子数上限（默认 2500 万）。
- **输出：** 字典，含 `score`、`found`、`free_end_gaps`、`query_start`、`query_end`、`target_start`、`target_end`（零起点半开区间）、`query_aligned`、`target_aligned`、`cigar`、`matches`、`mismatches`、`insertions`、`deletions`、`alignment_columns`、`identity` 和 `provenance`。
- **示例代码：**

```python
from dnakit.alignment import parasail_sg_dx_align

result = parasail_sg_dx_align("ACGTACGGA", "TTTTACGTACGGATTTT")
print(result["score"], result["free_end_gaps"], result["cigar"])
print(result["query_start"], result["query_end"], result["target_start"], result["target_end"])
```

- **示例输出：**

```text
45 ['database_begin', 'database_end'] 4D9=4D
0 9 4 13
```

## 7) `dnakit.alignment.parasail_sg_qb_db_align()`

- **作用：** 半全局比对，查询开头或目标开头可以不参与比对，但两者不能同时跳过；两条序列的末端都参与比对。
- **处理方法：** 调用 parasail 的 `sg_qb_db_trace` 族，代表成员是 `sg_qb_db_trace_scan_sat`。`method`（`scan`、`striped`、`diag`）和 `width`（8、16、32、64 或 `sat`）选择 SIMD 变体，通道饱和时报错。回溯路径包含自由末端的空位列，`cigar` 覆盖整条路径；坐标和一致性统计只计算比对核心区。
- **API：** `dnakit.alignment.parasail_sg_qb_db_align(query, target, *, match=5, mismatch=-4, gap_open=10, gap_extend=1, method="scan", width="sat", max_cells=25000000, verbose=True)`。
- **输入：** `query`、`target`：只含 ACGTN 的字符串或 DNA 对象；`match`、`mismatch`：碱基打分；`gap_open`、`gap_extend`：仿射空位，长度为 k 的空位扣 `gap_open + (k - 1) * gap_extend`；`method`、`width`：SIMD 变体；`max_cells`：动态规划格子数上限（默认 2500 万）。
- **输出：** 字典，含 `score`、`found`、`free_end_gaps`、`query_start`、`query_end`、`target_start`、`target_end`（零起点半开区间）、`query_aligned`、`target_aligned`、`cigar`、`matches`、`mismatches`、`insertions`、`deletions`、`alignment_columns`、`identity` 和 `provenance`。
- **示例代码：**

```python
from dnakit.alignment import parasail_sg_qb_db_align

result = parasail_sg_qb_db_align("TTTTACGTACGGA", "ACGTACGGA")
print(result["score"], result["free_end_gaps"], result["cigar"])
print(result["query_start"], result["query_end"], result["target_start"], result["target_end"])
```

- **示例输出：**

```text
45 ['query_begin', 'database_begin'] 4I9=
4 13 0 9
```

## 8) `dnakit.alignment.parasail_sg_qb_de_align()`

- **作用：** 半全局比对，查询开头和目标末尾都可以不参与比对，得到查询后缀与目标前缀的重叠。
- **处理方法：** 调用 parasail 的 `sg_qb_de_trace` 族，代表成员是 `sg_qb_de_trace_scan_sat`。`method`（`scan`、`striped`、`diag`）和 `width`（8、16、32、64 或 `sat`）选择 SIMD 变体，通道饱和时报错。回溯路径包含自由末端的空位列，`cigar` 覆盖整条路径；坐标和一致性统计只计算比对核心区。
- **API：** `dnakit.alignment.parasail_sg_qb_de_align(query, target, *, match=5, mismatch=-4, gap_open=10, gap_extend=1, method="scan", width="sat", max_cells=25000000, verbose=True)`。
- **输入：** `query`、`target`：只含 ACGTN 的字符串或 DNA 对象；`match`、`mismatch`：碱基打分；`gap_open`、`gap_extend`：仿射空位，长度为 k 的空位扣 `gap_open + (k - 1) * gap_extend`；`method`、`width`：SIMD 变体；`max_cells`：动态规划格子数上限（默认 2500 万）。
- **输出：** 字典，含 `score`、`found`、`free_end_gaps`、`query_start`、`query_end`、`target_start`、`target_end`（零起点半开区间）、`query_aligned`、`target_aligned`、`cigar`、`matches`、`mismatches`、`insertions`、`deletions`、`alignment_columns`、`identity` 和 `provenance`。
- **示例代码：**

```python
from dnakit.alignment import parasail_sg_qb_de_align

result = parasail_sg_qb_de_align("TTTTACGTACGGA", "ACGTACGGATTTT")
print(result["score"], result["free_end_gaps"], result["cigar"])
print(result["query_start"], result["query_end"], result["target_start"], result["target_end"])
```

- **示例输出：**

```text
45 ['query_begin', 'database_end'] 4I9=4D
4 13 0 9
```

## 9) `dnakit.alignment.parasail_sg_qe_db_align()`

- **作用：** 半全局比对，查询末尾和目标开头都可以不参与比对，得到查询前缀与目标后缀的重叠。
- **处理方法：** 调用 parasail 的 `sg_qe_db_trace` 族，代表成员是 `sg_qe_db_trace_scan_sat`。`method`（`scan`、`striped`、`diag`）和 `width`（8、16、32、64 或 `sat`）选择 SIMD 变体，通道饱和时报错。回溯路径包含自由末端的空位列，`cigar` 覆盖整条路径；坐标和一致性统计只计算比对核心区。
- **API：** `dnakit.alignment.parasail_sg_qe_db_align(query, target, *, match=5, mismatch=-4, gap_open=10, gap_extend=1, method="scan", width="sat", max_cells=25000000, verbose=True)`。
- **输入：** `query`、`target`：只含 ACGTN 的字符串或 DNA 对象；`match`、`mismatch`：碱基打分；`gap_open`、`gap_extend`：仿射空位，长度为 k 的空位扣 `gap_open + (k - 1) * gap_extend`；`method`、`width`：SIMD 变体；`max_cells`：动态规划格子数上限（默认 2500 万）。
- **输出：** 字典，含 `score`、`found`、`free_end_gaps`、`query_start`、`query_end`、`target_start`、`target_end`（零起点半开区间）、`query_aligned`、`target_aligned`、`cigar`、`matches`、`mismatches`、`insertions`、`deletions`、`alignment_columns`、`identity` 和 `provenance`。
- **示例代码：**

```python
from dnakit.alignment import parasail_sg_qe_db_align

result = parasail_sg_qe_db_align("ACGTACGGATTTT", "TTTTACGTACGGA")
print(result["score"], result["free_end_gaps"], result["cigar"])
print(result["query_start"], result["query_end"], result["target_start"], result["target_end"])
```

- **示例输出：**

```text
45 ['query_end', 'database_begin'] 4D9=4I
0 9 4 13
```

## 10) `dnakit.alignment.parasail_sg_qe_de_align()`

- **作用：** 半全局比对，查询末尾或目标末尾可以不参与比对，但两者不能同时跳过；两条序列的开头都参与比对。
- **处理方法：** 调用 parasail 的 `sg_qe_de_trace` 族，代表成员是 `sg_qe_de_trace_scan_sat`。`method`（`scan`、`striped`、`diag`）和 `width`（8、16、32、64 或 `sat`）选择 SIMD 变体，通道饱和时报错。回溯路径包含自由末端的空位列，`cigar` 覆盖整条路径；坐标和一致性统计只计算比对核心区。
- **API：** `dnakit.alignment.parasail_sg_qe_de_align(query, target, *, match=5, mismatch=-4, gap_open=10, gap_extend=1, method="scan", width="sat", max_cells=25000000, verbose=True)`。
- **输入：** `query`、`target`：只含 ACGTN 的字符串或 DNA 对象；`match`、`mismatch`：碱基打分；`gap_open`、`gap_extend`：仿射空位，长度为 k 的空位扣 `gap_open + (k - 1) * gap_extend`；`method`、`width`：SIMD 变体；`max_cells`：动态规划格子数上限（默认 2500 万）。
- **输出：** 字典，含 `score`、`found`、`free_end_gaps`、`query_start`、`query_end`、`target_start`、`target_end`（零起点半开区间）、`query_aligned`、`target_aligned`、`cigar`、`matches`、`mismatches`、`insertions`、`deletions`、`alignment_columns`、`identity` 和 `provenance`。
- **示例代码：**

```python
from dnakit.alignment import parasail_sg_qe_de_align

result = parasail_sg_qe_de_align("ACGTACGGATTTT", "ACGTACGGA")
print(result["score"], result["free_end_gaps"], result["cigar"])
print(result["query_start"], result["query_end"], result["target_start"], result["target_end"])
```

- **示例输出：**

```text
45 ['query_end', 'database_end'] 9=4I
0 9 0 9
```

## 11) `dnakit.alignment.parasail_nw_banded_score()`

- **作用：** 全局比对打分，只允许路径落在主对角线两侧 `band` 格以内；路径偏离过远时分数会低于完整的全局比对。
- **处理方法：** 调用 parasail 的 `nw_banded`，只返回分数，不回溯。若两条序列长度不同，parasail 会放宽带宽，保证终点可达。
- **API：** `dnakit.alignment.parasail_nw_banded_score(query, target, *, band=10, match=5, mismatch=-4, gap_open=10, gap_extend=1, max_cells=25000000, verbose=True)`。
- **输入：** `query`、`target`：ACGTN 序列；`band`：带宽半径（1 到 100000）；打分参数同上；`max_cells`：带宽内格子数上限（序列长度乘以 2 × band + 1）。
- **输出：** 字典，含 `score`、`band`、`query_length`、`target_length` 和 `provenance`。
- **示例代码：**

```python
from dnakit.alignment import parasail_nw_banded_score

query = "ACGTTGCAAGTCCATGGATCCAGTTACGGTAACGTTAGCCA"
target = query[:12] + "TTTTTT" + query[12:30] + query[36:]
print(parasail_nw_banded_score(query, target, band=1)["score"])
print(parasail_nw_banded_score(query, target, band=6)["score"])
```

- **示例输出：**

```text
81
145
```

## 12) `dnakit.alignment.parasail_ssw_align()`

- **作用：** 用 Striped Smith-Waterman（SSW）求最优局部比对，返回分数、坐标、CIGAR 和比对文本，适合短读段的局部定位。
- **处理方法：** 调用 parasail 的 `ssw`。SSW 使用 16 位分数，分数接近 32767 时 parasail 不返回结果，本函数报错。CIGAR 由 SSW 的操作码还原为 `=`、`X`、`I`、`D`，比对文本据此生成。
- **API：** `dnakit.alignment.parasail_ssw_align(query, target, *, match=5, mismatch=-4, gap_open=10, gap_extend=1, max_cells=25000000, verbose=True)`。
- **输入：** `query`、`target`：ACGTN 序列；打分与空位参数同上；`max_cells`：动态规划格子数上限。
- **输出：** 字典，含 `score`、`found`、`query_start`、`query_end`、`target_start`、`target_end`（零起点半开区间）、`query_aligned`、`target_aligned`、`cigar`、各类计数、`identity` 和 `provenance`。未找到比对时 `found` 为 False。
- **示例代码：**

```python
from dnakit.alignment import parasail_ssw_align

result = parasail_ssw_align("ACGTTGCAAGTCCATG", "ACGTTGCAAGGGTCCATG")
print(result["score"], result["cigar"])
print(result["query_aligned"])
print(result["target_aligned"])
```

- **示例输出：**

```text
69 9=2D7=
ACGTTGCAA--GTCCATG
ACGTTGCAAGGGTCCATG
```

## 13) `dnakit.alignment.edlib_nice_alignment()`

- **作用：** 用 Edlib 比对两条 DNA 序列，并返回 NICE 三行文本：查询序列、匹配线和目标序列。
- **处理方法：** 调用 `edlib.align(task="path")`，再调用 `edlib.getNiceAlignment()` 生成三行文本。编辑距离和扩展 CIGAR 与 `align_dna_edlib()` 一致。`mode` 为 `NW`（全局）或 `HW`（查询全局、目标局部）。
- **API：** `dnakit.alignment.edlib_nice_alignment(query, target, *, mode="NW", k=-1, iupac=False, gap_symbol="-", max_length=100000, verbose=True)`。
- **输入：** `query`、`target`：DNA 序列；`mode`：`NW` 或 `HW`；`k`：编辑距离上限（-1 表示不限）；`iupac`：是否让 IUPAC 歧义码等于其覆盖的碱基；`gap_symbol`：空位字符（单个非字母、非空格的可打印字符）；`max_length`：单条序列长度上限。
- **输出：** 字典，含 `found`、`edit_distance`、`query_aligned`、`matched_aligned`、`target_aligned`、`cigar`、`target_start`、`target_end`、`identity` 和 `provenance`。编辑距离超过 `k` 时 `found` 为 False。
- **示例代码：**

```python
from dnakit.alignment import edlib_nice_alignment

result = edlib_nice_alignment("ACGTTGCA", "GGACGTAGCAGG")
print(result["edit_distance"], result["cigar"])
print(result["query_aligned"])
print(result["matched_aligned"])
print(result["target_aligned"])
```

- **示例输出：**

```text
5 2D4=1X3=2D
--ACGTTGCA--
--||||.|||--
GGACGTAGCAGG
```

## 14) `dnakit.alignment.pyhmmer_align_to_hmm()`

- **作用：** 把命名的 DNA 序列对齐到一个 DNA 谱图 HMM，输出每条序列在模型列上的比对行。
- **处理方法：** 读取 HMMER3 文件中唯一的一个 DNA 谱图，用其字母表数字化序列，调用 `pyhmmer.hmmalign`。大写字母是匹配列，小写字母是模型外的插入残基，`.` 表示空位；`trim=True` 只保留模型列。蛋白质谱图会报错。
- **API：** `dnakit.alignment.pyhmmer_align_to_hmm(profile, sequences, *, trim=False, all_consensus_cols=True, max_sequences=10000, max_length=100000, verbose=True)`。
- **输入：** `profile`：HMMER3 文件路径（一个 DNA 谱图）；`sequences`：DNA 对象或 `名称 -> 序列` 映射；`trim`：是否只保留模型列；`all_consensus_cols`：是否保留全部共识列；`max_sequences`、`max_length`：数量与长度上限。
- **输出：** 字典，含 `profile`、`profile_length`、`consensus`、`sequence_count`、`alignment_columns`、`rows`（每项为 `name` 与 `aligned`）和 `provenance`。
- **示例代码：**

```python
from dnakit.alignment.backends import build_dna_hmm
from dnakit.alignment import pyhmmer_align_to_hmm

build_dna_hmm("aligned.fa", name="probe", output_path="probe.hmm", overwrite=True)
result = pyhmmer_align_to_hmm("probe.hmm", {"s1": "GGACGTACGTACGG", "s2": "ACGTTCGTAC"}, trim=True)
print(result["profile"], result["profile_length"], result["consensus"])
for row in result["rows"]:
    print(row["name"], row["aligned"])
```

- **示例输出：**

```text
probe 10 acgtacgtac
s1 ACGTACGTAC
s2 ACGTTCGTAC
```

## 安装可选依赖

`parasail`、`edlib` 和 `pyhmmer` 都是可选依赖，只有调用对应函数时才会导入。缺少时会报错并给出安装命令：
`dnakit[bio-parasail]`、`dnakit[bio-edlib]`、`dnakit[bio-pyhmmer]`。

## 验证方式

`others/tests/unit/test_full_G2.py` 用上游函数的直接调用和手工计算核对每个封装；网络相关函数不存在。

## 当前范围

- 半全局族只封装 `trace` 输出；`table`、`rowcol` 形式的动态规划矩阵属于调试输出，不在本页封装。
- `mappy.Aligner` 的映射与索引接口不在清单中，已由 `backends.py` 中的 `map_query_sequences_paf()` 等封装。
- 半全局族中，同一比对的两个起点（或两个终点）不能同时跳过：动态规划从第一行或第一列开始，在最后一行或最后一列结束，因此 `qb_db` 与 `qe_de` 每次只跳过一侧。
- SIMD 变体通过 `method` 与 `width` 选择；profile 输入形式（`*_profile_*`）不在本页暴露。
