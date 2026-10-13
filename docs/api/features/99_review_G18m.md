# 审查批次 G18m：methylpy、MACS3 与 deepTools 的 other 行

本页记录 G18m 审查批次对账本中 `other` 行的处理结果。三个包共 47 个 `other` 行，全部已重新标注状态，没有剩余的 `other` 或 `cap` 行。其中封装了 26 个函数：methylpy 12 个、MACS3 10 个、deepTools 4 个；其余函数按原因标记为 `needs-binary`、`deprecated`、`internal`、`not-dna`、`plot`、`duplicate-of` 或 `family-variant`，每个未封装函数的原因列在对应包的"未封装的函数"小节中。

所有封装函数都在 `dnakit.epigenome` 的 `*_review_functions.py` 模块中，调用上游 Python API，不运行命令行，不联网。示例从模块文件导入，示例输出由实际运行得到。

## methylpy 补充函数（DMR 步骤、allc 拆分、合并与二项检验）

这一节补充 methylpy 的 12 个函数，位于 `dnakit.epigenome` 的新模块 `methylpy_review_functions.py`，示例从该模块文件导入；导出到 `dnakit.epigenome` 后同样可以用 `dnakit.epigenome.<名称>()` 调用。所有函数都单进程运行，不使用 methylpy 自带的置换检验程序，也不需要 samtools、bowtie2、tabix 或 bgzip。

methylpy 的这些函数会就地排序输入文件、通过 shell 的 `sort`、`wc`、`rm` 处理路径、在出错时留下临时文件，并在表格格式不对时抛出底层错误。包装函数先检查输入（类型、范围、表头、字段数、染色体是否连续、输出路径是否已存在），不覆盖任何已有文件，让 methylpy 只处理副本或新路径，并把 methylpy 的退出和底层错误转为 `ConfigurationError`。返回值是带 `status`、`operation` 和 `provenance` 的字典。

示例按顺序在同一目录中运行，后面的示例会使用前面生成的 `methylpy_rms.tsv.gz`；输出文件写入临时目录，因此示例可以重复运行。

### 示例数据

两个样本 A、B 的 RMS 置换检验结果（未压缩；`num_sims=100`，`num_sig_tests=5`）。真实的结果表由 methylpy 的 `run_rms_tests.out` 生成并经 gzip 压缩，这里的数值是手工构造的。

**示例文件：** `methylpy_rms.tsv`

```text
chr	pos	strand	mc_class	pvalue	mc_A	mc_B	h_A	h_B	frac_A	frac_B	mc_residual_A	mc_residual_B	uc_residual_A	uc_residual_B	num_simulations_sig	num_simulations_run
1	100	+	CGN	0.01	9	1	10	10	0.900000	0.100000	3.577709	-3.577709	-3.577709	3.577709	1	100
1	120	+	CGN	0.01	10	0	10	10	1.000000	0.000000	4.472136	-4.472136	-4.472136	4.472136	1	100
1	140	+	CGN	0.02	8	1	10	10	0.800000	0.100000	3.146266	-3.146266	-3.146266	3.146266	2	100
1	300	+	CGN	0.50	5	6	10	10	0.500000	0.600000	-0.449467	0.449467	0.449467	-0.449467	5	10
1	320	+	CGN	0.42	4	5	10	10	0.400000	0.500000	-0.449467	0.449467	0.449467	-0.449467	5	12
1	340	+	CGN	0.62	3	4	10	10	0.300000	0.400000	-0.468807	0.468807	0.468807	-0.468807	5	8
1	600	+	CGN	0.01	1	9	10	10	0.100000	0.900000	-3.577709	3.577709	3.577709	-3.577709	1	100
1	620	+	CGN	0.03	0	10	10	10	0.000000	1.000000	-4.472136	4.472136	4.472136	-4.472136	3	100
2	50	+	CGN	0.50	5	5	10	10	0.500000	0.500000	0.000000	0.000000	0.000000	0.000000	5	10
2	80	+	CGN	0.01	9	2	10	10	0.900000	0.200000	3.146266	-3.146266	-3.146266	3.146266	1	100
```

样本 A 的 allc 表（染色体连续排列）。

**示例文件：** `methylpy_allc_A.tsv`

```text
chr1	100	+	CGA	9	10	1
chr1	120	+	CGT	10	10	1
chr1	125	+	CHH	1	4	0
chr1	140	+	CGC	8	10	1
chr1	500	+	CGG	5	10	1
chr2	30	+	CGA	0	9	0
```

样本 B 的 allc 表。

**示例文件：** `methylpy_allc_B.tsv`

```text
chr1	100	+	CGA	1	10	0
chr1	120	+	CGT	0	10	0
chr1	140	+	CGC	1	10	0
chr1	500	+	CGG	5	10	1
chr2	30	+	CGA	3	6	1
```

四个区域的 DMR 表（按染色体内起点排序；第二个区域与第一个重叠）。

**示例文件：** `methylpy_dmrs.tsv`

```text
#chr	start	end	number_of_dms	hypermethylated_samples	hypomethylated_samples
chr1	100	140	3	A	B
chr1	120	500	2	A	B
chr1	600	700	1	A	B
chr2	1	40	1	A	B
```

三条染色体的 allc 表，记录数分别为 3、1、2。

**示例文件：** `methylpy_allc_three.tsv`

```text
chr1	1	+	CGA	1	2	0
chr1	2	+	CGA	1	2	0
chr1	3	+	CGA	1	2	0
chr2	1	+	CGA	1	2	0
chr3	1	+	CGA	1	2	0
chr3	2	+	CGA	1	2	0
```

含未甲基化对照染色体 chrL 的 allc 表（chrL 上 200 条 reads 中有 2 条甲基化）。

**示例文件：** `methylpy_allc_S.tsv`

```text
chr1	10	+	CGA	8	10	1
chr1	20	+	CGT	1	10	1
chr1	30	+	CHH	0	9	0
chr1	40	+	CHH	5	9	1
chr1	50	+	CHH	1	50	0
chr1	60	+	CGC	0	0	0
chrL	10	+	CHH	1	100	0
chrL	20	+	CHH	1	100	0
```

第一个二项检验结果文件（第 7 列为 P 值，2.0 表示未检验的位点）。

**示例文件：** `methylpy_binom_1.tsv`

```text
chr1	20	+	CGN	1	10	0.2
chr1	10	+	CGN	8	10	0.0001
chr2	5	+	CHH	0	9	2.0
```

第二个二项检验结果文件。

**示例文件：** `methylpy_binom_2.tsv`

```text
chr1	40	+	CHH	5	9	0.001
chr1	30	+	CHH	0	9	0.6
```

七行 mpileup 文本。

**示例文件：** `methylpy_reads.mpileup`

```text
chr1	1	C	5	.....
chr1	2	G	5	.....
chr1	3	A	5	.....
chr1	4	C	5	.....
chr1	5	C	5	.....
chr1	6	T	5	.....
chr1	7	C	5	.....
```

### 1) `dnakit.epigenome.estimate_dmr_fdr_cutoff()`

- **作用：** 在目标 FDR 下，从 RMS 置换检验结果表求出 P 值阈值，用来判定哪些位点是差异甲基化位点（DMS）。
- **处理方法：** 调用 methylpy 的 `histogram_correction_DMRfind`（直方图法）。先检查 gzip 压缩、表头、每行字段数、P 值范围，以及最后两列 `num_simulations_sig`、`num_simulations_run` 是否与 `num_sims`、`num_sig_tests` 相容；再迭代估计零假设检验数 m0，直到变化不超过 `convergence_diff`，取估计 FDR 不超过 `target_fdr` 且最接近它的阈值。不收敛或没有显著位点时 methylpy 退回 `1/num_sims`，`outcome` 会标出。示例手工核对：10 个检验、m0 约 4.13、阈值 5/37。
- **API：** `dnakit.epigenome.estimate_dmr_fdr_cutoff(rms_path: 'str | Path', *, num_sims: 'int' = 3000, num_sig_tests: 'int' = 100, target_fdr: 'float' = 0.01, max_iterations: 'int' = 1000, convergence_diff: 'float' = 0.01, verbose: bool = True) -> 'dict[str, Any]'`.
- **输入：** `rms_path`：gzip 压缩的 RMS 结果表（表头为 chr、pos、strand、mc_class、pvalue，随后是每个样本的 mc_、h_、frac_、mc_residual_、uc_residual_ 列，最后是 num_simulations_sig、num_simulations_run）；`num_sims`、`num_sig_tests`：生成该表时的置换设置；`target_fdr`、`max_iterations`、`convergence_diff`：迭代设置。
- **输出：** `p_value_cutoff`、`total_tests`、`outcome`（`converged`、`not_converged` 或 `no_significant_result`）、`estimated_fdr`、`null_tests_estimate`（m0）。
- **示例代码：**

```python
import gzip
import shutil

from dnakit.epigenome import estimate_dmr_fdr_cutoff

with open("methylpy_rms.tsv", "rb") as source, gzip.open("methylpy_rms.tsv.gz", "wb") as target:
    shutil.copyfileobj(source, target)
result = estimate_dmr_fdr_cutoff(
    "methylpy_rms.tsv.gz",
    num_sims=100,
    num_sig_tests=5,
    target_fdr=0.1,
    convergence_diff=1,
    verbose=False,
)
print(result["total_tests"], result["outcome"])
print(round(result["p_value_cutoff"], 6), round(result["estimated_fdr"], 4))
```

- **示例输出：**

```text
10 converged
0.135135 0.0964
```

### 2) `dnakit.epigenome.dmr_residual_cutoff()`

- **作用：** 求 RMS 结果表中非显著位点的残差分布分位数，作为判断样本在差异位点上偏高或偏低的残差阈值。
- **处理方法：** 调用 methylpy 的 `get_resid_cutoff`：取 P 值大于 `pvalue_cutoff` 的位点，收集全部 mc 与 uc 残差，求第 `100 - resid_cutoff * 100` 百分位。没有这样的位点时 methylpy 返回 2，`default_used` 为真。methylpy 把没有数据的样本占位值 -1 也算进分布，函数照此返回。P 值以 `<` 开头的行会使 methylpy 抛出底层错误，这里改为 `ConfigurationError`。示例手工核对：4 个非显著位点中三个的残差为 ±0.4495、±0.4495、±0.4688，另一个为 0，第 80 百分位恰为 0.4495。
- **API：** `dnakit.epigenome.dmr_residual_cutoff(rms_path: 'str | Path', pvalue_cutoff: 'float', *, resid_cutoff: 'float' = 0.01, verbose: bool = True) -> 'dict[str, Any]'`.
- **输入：** `rms_path`：gzip 压缩的 RMS 结果表；`pvalue_cutoff`：0 到 1，通常取自 `estimate_dmr_fdr_cutoff`；`resid_cutoff`：0 到 1 的比例，默认 0.01。
- **输出：** `residual_cutoff`、`null_sites`（参与分布的位点数）、`default_used`。
- **示例代码：**

```python
from dnakit.epigenome import dmr_residual_cutoff

result = dmr_residual_cutoff("methylpy_rms.tsv.gz", 0.135135, resid_cutoff=0.2, verbose=False)
print(result["null_sites"], result["default_used"])
print(round(result["residual_cutoff"], 4))
```

- **示例输出：**

```text
4 False
0.4495
```

### 3) `dnakit.epigenome.collapse_dmr_windows()`

- **作用：** 把显著位点按距离和残差方向合并为差异甲基化区域（DMR），并写出 DMS 与 DMR 两个 BED 文件。
- **处理方法：** 调用 methylpy 的 `collapse_dmr_windows`（P 值列固定为第 5 列）。P 值不超过 `sig_cutoff` 的位点，若与当前区域同染色体、相距不超过 `max_dist`、残差方向一致，就并入该区域，否则开始新区域。`resid_cutoff` 大于 0 时才检查方向，并把样本标为高甲基化或低甲基化；`collapse_samples` 限定参与判断的样本，`sample_category` 与 `min_cluster` 要求各类别的同向样本数，`min_num_dms` 过滤过小的区域。输出表为 `<output_path>`，BED 为 `<output_path>.DMS.bed` 与 `<output_path>.DMR.bed`（染色体加 chr 前缀，起点从 0 计）。没有显著位点时只写表头，不写 BED。文件已存在时报错，不覆盖。
- **API：** `dnakit.epigenome.collapse_dmr_windows(rms_path: 'str | Path', output_path: 'str | Path', *, max_dist: 'int' = 100, resid_cutoff: 'float' = 0.0, sig_cutoff: 'float' = 0.01, min_num_dms: 'int' = 0, collapse_samples: 'Sequence[str] | None' = None, sample_category: 'Sequence[int] | None' = None, min_cluster: 'int' = 0, verbose: bool = True) -> 'dict[str, Any]'`.
- **输入：** `rms_path`：gzip 压缩的 RMS 结果表，按染色体和位置排序；`output_path`：新输出路径；`max_dist`、`resid_cutoff`、`sig_cutoff`、`min_num_dms`、`collapse_samples`、`sample_category`、`min_cluster`：合并设置。
- **输出：** `dmr_count`、`dms_count`、`regions`（至多 1000 个，含 `chrom`、`start`、`end`、`dms`、`hypermethylated_samples`、`hypomethylated_samples`）、`regions_truncated`、`files`（表与两个 BED 的路径，BED 未写出时为 None）。
- **示例代码：**

```python
import tempfile

from dnakit.epigenome import collapse_dmr_windows

with tempfile.TemporaryDirectory() as work:
    result = collapse_dmr_windows(
        "methylpy_rms.tsv.gz",
        f"{work}/dmrs.tsv",
        max_dist=50,
        sig_cutoff=0.135135,
        resid_cutoff=2.0,
        verbose=False,
    )
print(result["dmr_count"], result["dms_count"])
for region in result["regions"]:
    print(
        region["chrom"],
        region["start"],
        region["end"],
        region["dms"],
        region["hypermethylated_samples"],
        region["hypomethylated_samples"],
    )
```

- **示例输出：**

```text
3 6
1 100 140 3 ['A'] ['B']
1 600 620 2 ['B'] ['A']
2 80 80 1 ['A'] ['B']
```

### 4) `dnakit.epigenome.merge_dms_to_dmrs()`

- **作用：** 一步完成 FDR 阈值、残差阈值和区域合并，把 RMS 结果表变成 DMR 表。
- **处理方法：** 调用 methylpy 的 `merge_DMS_to_DMR`，它依次运行直方图 FDR 校正、残差阈值计算和 `collapse_dmr_windows`。上游不返回两个阈值，函数在合并后用相同参数重新计算它们并报告。`num_sims`、`num_sig_tests` 必须与生成 RMS 表时一致，否则（某行的模拟次数不可能出现时）报错。输出文件与 `collapse_dmr_windows` 相同。
- **API：** `dnakit.epigenome.merge_dms_to_dmrs(rms_path: 'str | Path', output_path: 'str | Path', *, num_sims: 'int' = 3000, num_sig_tests: 'int' = 100, sig_cutoff: 'float' = 0.01, dmr_max_dist: 'int' = 250, min_num_dms: 'int' = 0, resid_cutoff: 'float' = 0.01, collapse_samples: 'Sequence[str] | None' = None, sample_category: 'Sequence[int] | None' = None, min_cluster: 'int' = 0, max_iterations: 'int' = 1000, convergence_diff: 'float' = 1, verbose: bool = True) -> 'dict[str, Any]'`.
- **输入：** `rms_path`、`output_path`；`num_sims`、`num_sig_tests`；`sig_cutoff`：目标 FDR；`dmr_max_dist`；`min_num_dms`；`resid_cutoff`：残差百分位的比例；`collapse_samples`、`sample_category`、`min_cluster`；`max_iterations`、`convergence_diff`（默认 1，与上游一致）。
- **输出：** `p_value_cutoff`、`fdr_outcome`、`residual_cutoff`，以及与 `collapse_dmr_windows` 相同的 `dmr_count`、`dms_count`、`regions`、`files`。
- **示例代码：**

```python
import tempfile

from dnakit.epigenome import merge_dms_to_dmrs

with tempfile.TemporaryDirectory() as work:
    result = merge_dms_to_dmrs(
        "methylpy_rms.tsv.gz",
        f"{work}/merged.tsv",
        num_sims=100,
        num_sig_tests=5,
        sig_cutoff=0.1,
        dmr_max_dist=50,
        resid_cutoff=0.2,
        convergence_diff=1,
        verbose=False,
    )
print(round(result["p_value_cutoff"], 6), round(result["residual_cutoff"], 4))
print(result["dmr_count"], result["dms_count"])
print([(r["chrom"], r["start"], r["end"], r["dms"]) for r in result["regions"]])
```

- **示例输出：**

```text
0.135135 0.4495
3 6
[('1', 100, 140, 3), ('1', 600, 620, 2), ('2', 80, 80, 1)]
```

### 5) `dnakit.epigenome.add_methylation_levels_to_dmrs()`

- **作用：** 为 DMR 表的每个区域加上每个样本的甲基化水平（甲基化 reads 数除以总 reads 数）。
- **处理方法：** 调用 methylpy 的 `get_methylation_levels_DMRfind`，单进程运行。每个样本一个 allc 文件（全部染色体、同一染色体连续），索引 `.idx` 写在 allc 旁边。只统计上下文匹配 `mc_type` 且覆盖度在 [`min_cov`, `max_cov`] 内的位点；区域内没有位点时为 NA（返回 None）。区域须在同一染色体内按起点排序；区域所在染色体缺少于某个 allc 时报错（上游会沿用上一区域的数值，或在第一行抛出 NameError）。输出文件名须以 `.tsv` 结尾，路径只含字母、数字和 `_ . / + = : , @ % -`，因为 methylpy 通过 shell 的 `rm` 删除临时文件；出错时清理临时文件。
- **API：** `dnakit.epigenome.add_methylation_levels_to_dmrs(dmr_path: 'str | Path', allc_paths: 'Sequence[str | Path]', samples: 'Sequence[str]', output_path: 'str | Path', *, mc_type: 'str | Sequence[str]' = 'CGN', min_cov: 'int' = 0, max_cov: 'int | None' = None, buffer_line_number: 'int' = 100000, input_no_header: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **输入：** `dmr_path`：前三列为染色体、起点、终点的表（默认带一行表头）；`allc_paths`、`samples`：每个样本一个 allc 文件，名称为字母、数字、`_`、`.`、`-`；`output_path`；`mc_type`：上下文字符串或列表，默认 CGN；`min_cov`、`max_cov`、`buffer_line_number`；`input_no_header`：DMR 表无表头时为真。
- **输出：** `output`、`samples`、`mc_type`、`dmr_count`、`rows`（至多 1000 行，含 `chrom`、`start`、`end`、`methylation_level`）、`rows_truncated`。示例手工核对：区域 chr1:100-140 内 A 的 CGN 位点为 9+10+8 个甲基化 / 30 条 reads = 0.9（125 位的 CHH 位点被排除）。
- **示例代码：**

```python
import tempfile

from dnakit.epigenome import add_methylation_levels_to_dmrs

with tempfile.TemporaryDirectory() as work:
    result = add_methylation_levels_to_dmrs(
        "methylpy_dmrs.tsv",
        ["methylpy_allc_A.tsv", "methylpy_allc_B.tsv"],
        ["A", "B"],
        f"{work}/levels.tsv",
        verbose=False,
    )
for row in result["rows"]:
    levels = {k: None if v is None else round(v, 4) for k, v in row["methylation_level"].items()}
    print(row["chrom"], row["start"], row["end"], levels)
```

- **示例输出：**

```text
chr1 100 140 {'A': 0.9, 'B': 0.0667}
chr1 120 500 {'A': 0.7667, 'B': 0.2}
chr1 600 700 {'A': None, 'B': None}
chr2 1 40 {'A': 0.0, 'B': 0.5}
```

### 6) `dnakit.epigenome.add_cytosine_info_to_dmrs()`

- **作用：** 为 DMR 表的每个区域加上每个样本的甲基化水平、总 base call 数和覆盖位点数。
- **处理方法：** 调用 methylpy 的 `get_c_info_DMRfind`，规则与 `add_methylation_levels_to_dmrs` 相同。输出三个文件：`<output_path>`（甲基化水平）、`<output_path>.total_basecalls.tsv` 与 `<output_path>.covered_sites.tsv`。区域内没有符合条件的位点时三项均为 NA（None）。
- **API：** `dnakit.epigenome.add_cytosine_info_to_dmrs(dmr_path: 'str | Path', allc_paths: 'Sequence[str | Path]', samples: 'Sequence[str]', output_path: 'str | Path', *, mc_type: 'str | Sequence[str]' = 'CGN', min_cov: 'int' = 0, max_cov: 'int | None' = None, buffer_line_number: 'int' = 100000, input_no_header: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **输入：** 与 `add_methylation_levels_to_dmrs` 相同。
- **输出：** `output`、`total_basecalls_output`、`covered_sites_output`、`samples`、`mc_type`、`dmr_count`、`rows`（每行含 `chrom`、`start`、`end` 与按样本的 `methylation_level`、`total_basecalls`、`covered_sites`）、`rows_truncated`。
- **示例代码：**

```python
import tempfile

from dnakit.epigenome import add_cytosine_info_to_dmrs

with tempfile.TemporaryDirectory() as work:
    result = add_cytosine_info_to_dmrs(
        "methylpy_dmrs.tsv",
        ["methylpy_allc_A.tsv", "methylpy_allc_B.tsv"],
        ["A", "B"],
        f"{work}/info.tsv",
        mc_type="CGN",
        verbose=False,
    )
for row in result["rows"]:
    a, b = row["samples"]["A"], row["samples"]["B"]
    print(row["chrom"], row["start"], a["total_basecalls"], a["covered_sites"], b["total_basecalls"])
```

- **示例输出：**

```text
chr1 100 30 3 30
chr1 120 30 3 30
chr1 600 None None None
chr2 1 9 1 6
```

### 7) `dnakit.epigenome.split_allc_by_chromosome()`

- **作用：** 把 allc 表按染色体拆成 `allc_<sample>_<chrom>.tsv[.gz]`。
- **处理方法：** 调用 methylpy 的 `do_split_allc_file`。同一染色体必须连续出现，否则 methylpy 会覆盖前面写好的文件，所以先检查并报错；空文件、少于 6 列的行和不能用作文件名的染色体名也在调用前报错。methylpy 返回的文件名不带 `.gz`，即使写出的是 `.tsv.gz`，这里报告实际存在的文件。输出目录须已存在，目标文件不能已存在。
- **API：** `dnakit.epigenome.split_allc_by_chromosome(allc_path: 'str | Path', sample: 'str', output_dir: 'str | Path', *, compress_output: 'bool' = True, buffer_line_number: 'int' = 100000, verbose: bool = True) -> 'dict[str, Any]'`.
- **输入：** `allc_path`；`sample`：1 到 64 个字母、数字、`_`、`.`、`-`；`output_dir`；`compress_output`：默认真，写 `.tsv.gz`；`buffer_line_number`。
- **输出：** `files`（按染色体排序，含 `chrom`、`path`、`records`）、`records`、`compressed`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import split_allc_by_chromosome

with tempfile.TemporaryDirectory() as work:
    result = split_allc_by_chromosome("methylpy_allc_A.tsv", "A", work, verbose=False)
    for item in result["files"]:
        print(item["chrom"], Path(item["path"]).name, item["records"])
```

- **示例输出：**

```text
chr1 allc_A_chr1.tsv.gz 5
chr2 allc_A_chr2.tsv.gz 1
```

### 8) `dnakit.epigenome.split_allc_into_chunks()`

- **作用：** 把 allc 表的染色体分配到至多 `num_chunks` 个文件，名称为 `allc_<sample>_<index>.tsv[.gz]`。
- **处理方法：** 调用 methylpy 的 `do_split_allc_file_chunk`：每块容量为 `ceil(总记录数 / num_chunks)`，染色体按记录数从小到大装入，当前块满后开新块，染色体不会被切开。染色体较少或大小悬殊时，写出的文件比 `num_chunks` 少。块内保持输入中的行序。空文件在调用前报错；methylpy 对空文件返回 0，不写文件。
- **API：** `dnakit.epigenome.split_allc_into_chunks(allc_path: 'str | Path', sample: 'str', num_chunks: 'int', output_dir: 'str | Path', *, compress_output: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`.
- **输入：** `allc_path`、`sample`、`num_chunks`（1 到 10000）、`output_dir`；`compress_output`。
- **输出：** `chunks`（含 `index`、`path`、`records`、`chromosomes`，由写出的文件读回）、`requested_chunks`、`records`、`compressed`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import split_allc_into_chunks

with tempfile.TemporaryDirectory() as work:
    result = split_allc_into_chunks("methylpy_allc_three.tsv", "S", 2, work, verbose=False)
    for chunk in result["chunks"]:
        print(Path(chunk["path"]).name, chunk["records"], chunk["chromosomes"])
```

- **示例输出：**

```text
allc_S_0.tsv.gz 3 ['chr2', 'chr3']
allc_S_1.tsv.gz 3 ['chr1']
```

### 9) `dnakit.epigenome.run_allc_binomial_test()`

- **作用：** 用二项检验和每个上下文的 Benjamini-Hochberg 阈值从 allc 表中调用甲基化位点。
- **处理方法：** 调用 methylpy 的 `perform_binomial_test`（单进程）。未转化率来自 `unmethylated_control`：0 到 1 的数，或未甲基化对照区域（如 `chrL`、`chrL:100`、`chrL:100-2000`），后者取区域内甲基化 reads 数除以总 reads 数。区域的染色体名按原样查找（`perform_binomial_test` 会去掉 `chr` 前缀，所以这里先自行求出比率）。总 reads 数不小于 `min_cov` 且有甲基化 reads 的位点参与检验，其余有 reads 的位点记为 0，没有 reads 的位点丢弃。输出 allc 的第 7 列为 1 表示显著。allc 索引写在输入旁边；`sort` 在 `LC_ALL=C` 下运行，结束后恢复原环境变量。输出目录路径只能含字母、数字和 `_ . / + = : , @ % -`。
- **API：** `dnakit.epigenome.run_allc_binomial_test(allc_path: 'str | Path', sample: 'str', output_dir: 'str | Path', unmethylated_control: 'float | str', *, min_cov: 'int' = 2, sig_cutoff: 'float' = 0.01, sort_mem: 'str' = '500M', compress_output: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`.
- **输入：** `allc_path`、`sample`、`output_dir`、`unmethylated_control`；`min_cov`：默认 2；`sig_cutoff`：FDR，默认 0.01；`sort_mem`：如 `500M`；`compress_output`。
- **输出：** `output`（`allc_<sample>.tsv[.gz]`）、`index`、`non_conversion_rate`、`sites`、`significant_sites`、`contexts`（按上下文排序的 `sites` 与 `significant_sites`）。示例手工核对：chrL 为 2/200，比率 0.01；8/10 与 5/9 的位点显著，1/10、1/50、1/100 不显著。
- **示例代码：**

```python
import gzip
import tempfile

from dnakit.epigenome import run_allc_binomial_test

with tempfile.TemporaryDirectory() as work:
    result = run_allc_binomial_test("methylpy_allc_S.tsv", "S", work, "chrL", verbose=False)
    with gzip.open(result["output"], "rt") as handle:
        flags = [line.split("\t")[:2] + line.split("\t")[6:7] for line in handle]
print(result["non_conversion_rate"], result["sites"], result["significant_sites"])
print([item["context"] + ":" + str(item["significant_sites"]) for item in result["contexts"]])
print([f"{row[0]}:{row[1]}={row[2].strip()}" for row in flags])
```

- **示例输出：**

```text
0.01 7 2
['CGA:1', 'CGT:0', 'CHH:1']
['chr1:10=1', 'chr1:20=0', 'chr1:30=0', 'chr1:40=1', 'chr1:50=0', 'chrL:10=0', 'chrL:20=0']
```

### 10) `dnakit.epigenome.filter_allc_by_pvalue_cutoffs()`

- **作用：** 按每个上下文的 P 值阈值，把二项检验结果文件改写成带显著标记的 allc 记录。
- **处理方法：** 调用 methylpy 的 `filter_files_by_pvalue_combined`。上游会用 `sort -o` 就地排序输入文件，所以先把输入复制到临时目录再调用，原文件不变。每个文件按染色体和位置排序，P 值（第 7 列）不超过所属上下文阈值的位点第 7 列写 1，其余写 0；P 值为 2.0 表示未检验，永远是 0。多个文件按给定顺序依次写出。上游在输出名不带 `.gz` 时会补上，这里报告真实路径。`sort` 在 `LC_ALL=C` 下运行。
- **API：** `dnakit.epigenome.filter_allc_by_pvalue_cutoffs(binomial_paths: 'Sequence[str | Path]', output_path: 'str | Path', cutoffs: 'Mapping[str, float]', *, compress_output: 'bool' = True, sort_mem: 'str' = '500M', verbose: bool = True) -> 'dict[str, Any]'`.
- **输入：** `binomial_paths`：1 到 32 个结果文件（每行至少 7 列）；`output_path`；`cutoffs`：上下文到阈值（0 到 1）的映射，必须覆盖文件中所有已检验的上下文；`compress_output`、`sort_mem`。
- **输出：** `output`、`records`、`significant_sites`、`contexts`（按上下文排序的 `sites` 与 `significant_sites`）、`inputs`。
- **示例代码：**

```python
import gzip
import tempfile

from dnakit.epigenome import filter_allc_by_pvalue_cutoffs

with tempfile.TemporaryDirectory() as work:
    result = filter_allc_by_pvalue_cutoffs(
        ["methylpy_binom_1.tsv", "methylpy_binom_2.tsv"],
        f"{work}/filtered.tsv",
        {"CGN": 0.05, "CHH": 0.01},
        verbose=False,
    )
    with gzip.open(result["output"], "rt") as handle:
        lines = handle.read().splitlines()
print(result["records"], result["significant_sites"])
for line in lines:
    print(line)
```

- **示例输出：**

```text
5 2
chr1	10	+	CGN	8	10	1
chr1	20	+	CGN	1	10	0
chr2	5	+	CHH	0	9	0
chr1	30	+	CHH	0	9	0
chr1	40	+	CHH	5	9	1
```

### 11) `dnakit.epigenome.merge_allc_tables()`

- **作用：** 按位点累加多个 allc 表的甲基化 reads 数和总 reads 数。
- **处理方法：** 调用 methylpy 的 `merge_allc_files`（单进程；`mini_batch` 控制每次合并的文件数，即上游的 minibatch 步骤）。上游在 `compress_output` 为真且输出名不带 `.gz` 时会补上 `.gz`，这里报告真实路径并要求该路径及其 `.idx` 不存在。上游的染色体写出顺序随进程的哈希种子变化，函数在合并后按染色体名排序重写并重建索引，使结果可重复。输入先检查：同一染色体连续、计数为整数；`skip_snp_info` 为假时每行至少 9 列，匹配与错配列按位置相加。输入的 `.idx` 写在输入旁边。
- **API：** `dnakit.epigenome.merge_allc_tables(allc_paths: 'Sequence[str | Path]', output_path: 'str | Path', *, mini_batch: 'int' = 100, compress_output: 'bool' = True, skip_snp_info: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`.
- **输入：** `allc_paths`：2 到 1000 个 allc 表；`output_path`；`mini_batch`：2 到 1000；`compress_output`；`skip_snp_info`。
- **输出：** `output`、`index`、`records`、`chromosomes`（按染色体名排序的 `chrom` 与 `records`）、`compressed`、`inputs`。
- **示例代码：**

```python
import gzip
import tempfile
from pathlib import Path

from dnakit.epigenome import merge_allc_tables

with tempfile.TemporaryDirectory() as work:
    result = merge_allc_tables(
        ["methylpy_allc_A.tsv", "methylpy_allc_B.tsv"], f"{work}/merged.tsv", verbose=False
    )
    print(Path(result["output"]).name, result["records"])
    with gzip.open(result["output"], "rt") as handle:
        for line in handle:
            print(line.rstrip("\n"))
```

- **示例输出：**

```text
merged.tsv.gz 6
chr1	100	+	CGA	10	20	1
chr1	120	+	CGT	10	20	1
chr1	125	+	CHH	1	4	1
chr1	140	+	CGC	9	20	1
chr1	500	+	CGG	10	20	1
chr2	30	+	CGA	3	15	1
```

### 12) `dnakit.epigenome.split_mpileup_chunks()`

- **作用：** 把纯文本 mpileup 文件切成若干块，相邻块共享两行，使块边界处的 C 上下文保持完整。
- **处理方法：** 调用 methylpy 的 `split_mpileup_file`：每块 `ceil(行数 / num_chunks)` 行，第二块起开头重复上一块的最后两行。上游用 `wc -l` 数行，对 gzip 文件会数到压缩字节里的换行，所以 gzip 输入被拒绝；空文件时上游退出，转为 `ConfigurationError`。块文件名为 `<output_prefix>0`、`<output_prefix>1`……，调用前检查它们都不存在。实际块数可以少于 `num_chunks`。
- **API：** `dnakit.epigenome.split_mpileup_chunks(mpileup_path: 'str | Path', output_prefix: 'str | Path', num_chunks: 'int', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **输入：** `mpileup_path`：纯文本文件；`output_prefix`：块文件的路径前缀（以文件名前缀结尾，不是目录）；`num_chunks`：1 到 10000。
- **输出：** `records`、`chunk_size`、`chunks`（含 `index`、`path`、`records`、`overlap_records`，由写出的文件读回）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import split_mpileup_chunks

with tempfile.TemporaryDirectory() as work:
    result = split_mpileup_chunks("methylpy_reads.mpileup", f"{work}/chunk_", 3, verbose=False)
    for chunk in result["chunks"]:
        print(Path(chunk["path"]).name, chunk["records"], chunk["overlap_records"])
```

- **示例输出：**

```text
chunk_0 3 0
chunk_1 5 2
chunk_2 3 2
```

### 未封装的函数

以下 5 个函数没有新增包装函数，状态与理由如下（台账中的理由为英文）。

- `methylpy.DMRfind.DMRfind`（`needs-binary`）：端到端 DMR 流程，需要 methylpy 自带的置换检验程序 `run_rms_tests.out`（依赖 `libgsl.so.0`，本环境没有），并会删除临时 allc 文件；它的后处理步骤已分别封装。
- `methylpy.DMRfind.run_rms_tests`（`needs-binary`）：调用 methylpy 自带的编译程序 `run_rms_tests.out`（依赖 `libgsl.so.0`，本环境没有）。
- `methylpy.call_mc_se.fasta_iter`（`deprecated`）：Python 2 代码：对 `groupby` 迭代器调用 `.next()`，在 Python 3 上抛出 `AttributeError`（已在小型 FASTA 上验证）；读取染色体序列请用 `read_fasta_chromosome`。
- `methylpy.utilities.infer_samples`（`internal`）：从 `allc_<sample>.tsv` 文件名截取样本名的辅助函数，与 `DMRfind.infer_samples` 是同一对象。
- `methylpy.utilities.merge_allc_files_minibatch`（`family-variant` → `merge_allc_tables`）：`merge_allc_files` 内部的分批合并步骤，通过 `merge_allc_tables` 的 `mini_batch` 选项使用。

## MACS3 补充函数（HMMRATAC 流程）

这一节补充 MACS3 的 10 个函数，位于 `dnakit.epigenome` 的新模块 `MACS3_review_functions.py`，示例从该模块文件导入；导出到 `dnakit.epigenome` 后同样可以用 `dnakit.epigenome.<名称>()` 调用。它们对应 HMMRATAC（ATAC-seq 开放染色质区域）流程的各个步骤：把片段拆成 short、mono、di、tri 四条信号轨道，切成 bin，训练、保存、读取和使用三状态隐马尔可夫模型，再把预测概率转成 bedGraph、状态路径、BED 和 narrowPeak。不运行 `macs3 hmmratac` 命令行，也不需要外部程序。

MACS3 的这些函数之间传递的是片段轨道和 hmmlearn 模型对象，包装函数无法跨调用保存它们，因此用普通数据交接：片段是 `(染色体, 左端, 右端)` 三元组，训练数据是每个 bin 四个数的行，模型是 MACS3 的 JSON 文件，预测概率是 MACS3 流程里放在临时文件中的那种文本文件（`染色体,bin 终点,概率,概率,概率`）。上游函数遇到非法输入会触发断言、`sys.exit` 或底层异常，包装函数先检查输入，不覆盖任何已有文件，并把这些错误转为 `ConfigurationError`。返回值是带 `status`、`operation` 和 `provenance` 的字典，坐标都是 0 起点、右开。

示例数据是手工构造的，用来说明格式和计算，不代表真实实验的结果。示例的输出文件写入临时目录，因此可以重复运行。

### 示例数据

45 个 bin 的训练数据，每行是 short、mono、di、tri 四条信号轨道的取值；按 15 个 bin 一组分为 3 组，每组依次经过背景、核小体、开放和背景状态。数值是手工构造的。

**示例文件：** `macs3_training_data.tsv`

```text
0.30	0.69	0.27	0.11
0.29	0.58	0.37	0.10
0.30	0.65	0.34	0.10
0.32	0.53	0.29	0.09
0.84	4.91	1.61	0.49
0.98	5.77	2.02	0.42
0.99	6.17	2.18	0.45
7.62	1.52	0.75	0.22
6.64	2.26	0.59	0.33
8.31	1.93	0.84	0.32
9.00	1.94	0.74	0.28
0.88	5.97	1.81	0.56
0.78	5.21	1.77	0.37
0.37	0.43	0.29	0.09
0.36	0.46	0.34	0.09
0.29	0.55	0.32	0.09
0.30	0.63	0.37	0.07
0.35	0.67	0.28	0.10
0.28	0.72	0.31	0.10
0.97	5.86	1.96	0.45
1.25	4.62	1.13	0.49
0.98	6.27	1.95	0.49
8.32	2.23	0.76	0.29
9.86	2.13	0.71	0.38
8.74	1.86	0.69	0.31
7.20	1.75	0.68	0.28
1.13	5.69	1.65	0.54
1.01	6.61	2.29	0.49
0.29	0.60	0.26	0.11
0.35	0.61	0.29	0.10
0.27	0.54	0.29	0.09
0.28	0.49	0.31	0.10
0.26	0.43	0.30	0.11
0.27	0.57	0.28	0.11
0.89	6.71	1.93	0.56
1.00	5.83	1.65	0.46
0.97	6.47	2.06	0.46
8.39	2.24	0.79	0.28
7.62	2.19	0.85	0.27
8.36	1.88	0.73	0.34
8.79	1.83	0.81	0.32
0.92	5.91	2.16	0.39
1.04	6.53	2.12	0.42
0.31	0.54	0.32	0.11
0.31	0.54	0.28	0.11
```

手写的 MACS3 模型文件（对角协方差）：状态 0 是背景，状态 1 是核小体，状态 2 是开放。真实的模型文件由 `save_hmmratac_model` 或 `macs3 hmmratac` 生成。

**示例文件：** `macs3_model.json`

```text
{
 "startprob": [0.2, 0.3, 0.5],
 "transmat": [[0.8, 0.1, 0.1], [0.1, 0.8, 0.1], [0.1, 0.1, 0.8]],
 "means": [[1.0, 1.0, 1.0, 1.0], [5.0, 5.0, 5.0, 5.0], [10.0, 10.0, 10.0, 10.0]],
 "covars": [[1.0, 1.0, 1.0, 1.0], [4.0, 4.0, 4.0, 4.0], [9.0, 9.0, 9.0, 9.0]],
 "covariance_type": "diag",
 "n_features": 4,
 "i_open_region": 2,
 "i_background_region": 0,
 "i_nucleosomal_region": 1,
 "hmm_binsize": 10,
 "hmm_type": "gaussian"
}
```

预测概率文件：每行为 `染色体,bin 终点,状态0概率,状态1概率,状态2概率`，bin 大小为 10。这里状态 2 是开放，状态 0 是核小体，状态 1 是背景。数值是手工构造的。

**示例文件：** `macs3_proba.csv`

```text
chr1,20,0.1,0.8,0.1
chr1,30,0.8,0.1,0.1
chr1,40,0.1,0.1,0.8
chr1,50,0.1,0.1,0.8
chr1,60,0.8,0.1,0.1
chr1,70,0.1,0.8,0.1
chr1,100,0.8,0.1,0.1
chr1,110,0.1,0.1,0.8
chr1,120,0.8,0.1,0.1
chr2,40,0.1,0.8,0.1
chr2,50,0.8,0.1,0.1
chr2,60,0.1,0.1,0.8
chr2,70,0.8,0.1,0.1
```

上面的预测概率文件中开放状态概率的 bedGraph 轨道（`hmmratac_proba_to_bedgraph` 的输出），用作开放区域的评分轨道。

**示例文件：** `macs3_open.bdg`

```text
chr1	0	10	0.00000
chr1	10	30	0.10000
chr1	30	50	0.80000
chr1	50	70	0.10000
chr1	70	90	0.00000
chr1	90	100	0.10000
chr1	100	110	0.80000
chr1	110	120	0.10000
chr2	0	30	0.00000
chr2	30	50	0.10000
chr2	50	60	0.80000
chr2	60	70	0.10000
```

### 1) `dnakit.epigenome.digest_hmmratac_fragments()`

- **作用：** 把双端测序（如 ATAC-seq）的片段按长度拆成 short、mono、di、tri 四条信号轨道，这是 HMMRATAC 的输入信号。
- **处理方法：** 先调用 `hmmratac_fragment_weights` 算出每个片段长度在四个高斯分布下的权重，再用片段建立 MACS3 的 `PETrackI` 轨道，调用 `HMMR_Signal_Processing.generate_digested_signals`：每个片段在其覆盖范围内给四条轨道各加上对应权重。四个高斯的密度都低于 `min_fragment_probability` 的片段权重为 0。可选地按 MACS3 的格式写出四个 bedGraph 文件。
- **API：** `dnakit.epigenome.digest_hmmratac_fragments(fragments: 'Sequence[Sequence[Any]]', means: 'Sequence[float]', stddevs: 'Sequence[float]', *, min_fragment_probability: 'float' = 0.001, output_prefix: 'str | Path | None' = None, include_tracks: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`.
- **输入：** `fragments` 为 `(染色体, 左端, 右端)` 三元组列表（0 起点、右开，长度不超过 10000，染色体名不含空格和逗号）；`means` 与 `stddevs` 各四个值，顺序为 short、mono、di、tri；`min_fragment_probability`；`output_prefix` 不为空时写出 `<前缀>_short.bdg`、`_mono.bdg`、`_di.bdg`、`_tri.bdg`（已存在则报错，不覆盖）；`include_tracks=False` 时返回值里不放轨道。
- **输出：** `tracks` 为每类轨道、每条染色体的 `positions`（每一段的终点，第一段从 0 开始）与 `values`（该段的权重和）；`intervals` 为每类轨道的段数；还有 `fragments`、`chromosomes`、`classes`；有 `output_prefix` 时还有 `files`。示例中两个 50 bp 片段重叠在 120–150，short 轨道在那里为 2；200 bp 的片段属于 mono 轨道。
- **示例代码：**

```python
from dnakit.epigenome import digest_hmmratac_fragments

fragments = [("chr1", 100, 150), ("chr1", 120, 170), ("chr1", 300, 500)]
result = digest_hmmratac_fragments(
    fragments, [50, 200, 400, 600], [20, 20, 20, 20], verbose=False
)
print(result["fragments"], result["chromosomes"], result["intervals"])
for label in result["classes"]:
    track = result["tracks"][label]["chr1"]
    print(label, track["positions"], [round(value, 3) for value in track["values"]])
```

- **示例输出：**

```text
3 ['chr1'] {'short': 6, 'mono': 6, 'di': 2, 'tri': 1}
short [100, 120, 150, 170, 300, 500] [0.0, 1.0, 2.0, 1.0, 0.0, 0.0]
mono [100, 120, 150, 170, 300, 500] [0.0, 0.0, 0.0, 0.0, 0.0, 1.0]
di [300, 500] [0.0, 0.0]
tri [500] [0.0]
```

### 2) `dnakit.epigenome.extract_hmmratac_training_data()`

- **作用：** 在给定区域内把四条信号轨道切成等宽的 bin，得到 HMM 训练或解码用的数据行。
- **处理方法：** 先按上一个函数生成四条轨道；把区域向两侧扩展 `flank`、合并重叠（与 HMMRATAC 流程一致），再调用 `HMMR_Signal_Processing.extract_signals_from_regions`：每个 bin 取四条轨道在该处的值，低于 0.0001 的值抬到 0.0001；`hmm_type="poisson"` 时取整。bin 以终点表示，覆盖 `[终点 - binsize, 终点)`；超出信号末端的 bin 不会产生；染色体按字母倒序排列（MACS3 的行为）。区域内没有任何 bin 时报错。
- **API：** `dnakit.epigenome.extract_hmmratac_training_data(fragments: 'Sequence[Sequence[Any]]', regions: 'Sequence[Sequence[Any]]', means: 'Sequence[float]', stddevs: 'Sequence[float]', *, binsize: 'int' = 10, flank: 'int' = 0, hmm_type: 'str' = 'gaussian', min_fragment_probability: 'float' = 0.001, verbose: bool = True) -> 'dict[str, Any]'`.
- **输入：** `fragments`（同上）；`regions` 为 `(染色体, 起点, 终点)` 三元组列表；`means`、`stddevs`；`binsize`（默认 10）；`flank`（默认 0）；`hmm_type` 为 `"gaussian"` 或 `"poisson"`；`min_fragment_probability`。
- **输出：** `bins`（`[染色体, 终点]`）、`data`（每个 bin 四个值）、`lengths`（每个合并后区域的 bin 数，用作 HMM 的序列长度）、`regions`（合并后的区域数）、`binsize`。示例中两个区域各有 3 个 bin。
- **示例代码：**

```python
from dnakit.epigenome import extract_hmmratac_training_data

fragments = [("chr1", 100, 150), ("chr1", 120, 170), ("chr1", 300, 500)]
result = extract_hmmratac_training_data(
    fragments,
    [("chr1", 100, 170), ("chr1", 300, 360)],
    [50, 200, 400, 600],
    [20, 20, 20, 20],
    binsize=20,
    verbose=False,
)
print(result["regions"], result["lengths"])
for (chrom, end), row in zip(result["bins"], result["data"]):
    print(chrom, end - 20, end, [round(value, 3) for value in row])
```

- **示例输出：**

```text
2 [3, 3]
chr1 100 120 [1.0, 0.0, 0.0, 0.0]
chr1 120 140 [2.0, 0.0, 0.0, 0.0]
chr1 140 160 [1.0, 0.0, 0.0, 0.0]
chr1 300 320 [0.0, 1.0, 0.0, 0.0]
chr1 320 340 [0.0, 1.0, 0.0, 0.0]
chr1 340 360 [0.0, 1.0, 0.0, 0.0]
```

### 3) `dnakit.epigenome.train_hmmratac_model()`

- **作用：** 在分 bin 的四轨道数据上训练 HMMRATAC 的三状态隐马尔可夫模型，并标出开放、核小体和背景状态。
- **处理方法：** 调用 `HMMR_HMM.hmm_training`（hmmlearn 的 `GaussianHMM` 或 `PoissonHMM`，Baum-Welch 训练）。随机数发生器由 `random_seed` 在局部创建，不改变全局随机状态。状态按 HMMRATAC 流程标记：四个均值（Poisson 为 lambda）之和最大的是开放态，最小的是背景态，另一个是核小体态。EM 可能落入局部最优，结果依赖种子，可换 `random_seed` 重试；同一种子的两次训练数值可能相差约 1e-10（并行浮点求和顺序）。拟合出的协方差不正定时报错。
- **API：** `dnakit.epigenome.train_hmmratac_model(data: 'Sequence[Sequence[float]]', lengths: 'Sequence[int]', *, hmm_type: 'str' = 'gaussian', covariance_type: 'str' = 'full', random_seed: 'int' = 12345, verbose: bool = True) -> 'dict[str, Any]'`.
- **输入：** `data` 为每个 bin 四个非负数（Poisson 需为整数）；`lengths` 为各序列的 bin 数，总和等于行数；`hmm_type`（默认 `"gaussian"`）；`covariance_type` 为 `"full"` 或 `"diag"`；`random_seed`（默认 12345）。
- **输出：** `open_state`、`nucleosomal_state`、`background_state`（状态编号）、`start_probabilities`、`transition_matrix`、高斯模型的 `means` 与 `covariances`（Poisson 为 `lambdas`）、`converged`、`iterations`、`log_likelihood`。示例的训练数据有三个清晰的水平，模型找回了背景、开放（short 轨道高）和核小体（mono 轨道高）三种状态。
- **示例代码：**

```python
from dnakit.epigenome import train_hmmratac_model

with open("macs3_training_data.tsv", encoding="utf-8") as handle:
    rows = [[float(value) for value in line.split("\t")] for line in handle]
result = train_hmmratac_model(rows, [15, 15, 15], random_seed=10151, verbose=False)
print(result["converged"], result["open_state"], result["nucleosomal_state"], result["background_state"])
for state, means in enumerate(result["means"]):
    print(state, [round(value, 1) for value in means])
print([round(value, 2) for value in result["transition_matrix"][result["open_state"]]])
```

- **示例输出：**

```text
True 1 2 0
0 [0.3, 0.6, 0.3, 0.1]
1 [8.2, 2.0, 0.7, 0.3]
2 [1.0, 5.9, 1.9, 0.5]
[0.0, 0.75, 0.25]
```

### 4) `dnakit.epigenome.save_hmmratac_model()`

- **作用：** 训练模型并写成 MACS3 的 JSON 模型文件（含状态编号和 bin 大小），供预测使用。
- **处理方法：** 训练和状态标记同 `train_hmmratac_model`，再调用 `HMMR_HMM.hmm_model_save` 写文件，并立刻用 `hmm_model_init` 读回检验。写出的文件无法读回时删除该文件并报错。不覆盖已有文件，输出目录必须已存在。
- **API：** `dnakit.epigenome.save_hmmratac_model(data: 'Sequence[Sequence[float]]', lengths: 'Sequence[int]', output_path: 'str | Path', *, binsize: 'int' = 10, hmm_type: 'str' = 'gaussian', covariance_type: 'str' = 'full', random_seed: 'int' = 12345, verbose: bool = True) -> 'dict[str, Any]'`.
- **输入：** `data`、`lengths`、`hmm_type`、`covariance_type`、`random_seed` 同训练函数；`output_path` 为要写入的 JSON 文件；`binsize` 记录训练数据的 bin 大小，默认 10。
- **输出：** `output`（文件路径）、`binsize`、`converged`、三个状态编号。文件内容为 MACS3 的模型 JSON：`startprob`、`transmat`、`means` 与 `covars`（Poisson 为 `lambdas`）、`i_open_region` 等编号、`hmm_binsize`、`hmm_type`。
- **示例代码：**

```python
import json
import tempfile

from dnakit.epigenome import save_hmmratac_model

with open("macs3_training_data.tsv", encoding="utf-8") as handle:
    rows = [[float(value) for value in line.split("\t")] for line in handle]
with tempfile.TemporaryDirectory() as work:
    result = save_hmmratac_model(
        rows, [15, 15, 15], f"{work}/model.json", binsize=10, random_seed=10151, verbose=False
    )
    with open(result["output"], encoding="utf-8") as handle:
        model = json.load(handle)
print(result["converged"], result["binsize"])
print(result["open_state"], result["nucleosomal_state"], result["background_state"])
print(sorted(model))
print(model["hmm_type"], model["covariance_type"], model["hmm_binsize"], len(model["means"]))
```

- **示例输出：**

```text
True 10
1 2 0
['covariance_type', 'covars', 'hmm_binsize', 'hmm_type', 'i_background_region', 'i_nucleosomal_region', 'i_open_region', 'means', 'n_features', 'startprob', 'transmat']
gaussian full 10 3
```

### 5) `dnakit.epigenome.load_hmmratac_model()`

- **作用：** 读取 MACS3 的 HMMRATAC 模型文件，返回起始概率、转移矩阵、发射参数、状态编号和 bin 大小。
- **处理方法：** 先检查 JSON（数组形状、`n_features` 为 4、三个状态编号互不相同），再调用 `HMMR_HMM.hmm_model_init` 由文件重建模型，参数取自重建后的模型。没有 `hmm_type` 的旧版文件按高斯模型读取。协方差不正定时报错。
- **API：** `dnakit.epigenome.load_hmmratac_model(model_path: 'str | Path', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **输入：** `model_path` 为模型 JSON 文件。
- **输出：** `hmm_type`、`binsize`、`open_state`、`nucleosomal_state`、`background_state`、`state_labels`（按状态编号给出 `open`、`nuc`、`bg`）、`start_probabilities`、`transition_matrix`、高斯模型的 `means`、`covariances`、`covariance_type`（Poisson 为 `lambdas`）、`features`。示例的手写模型中状态 2 为开放，状态 1 为核小体，状态 0 为背景。
- **示例代码：**

```python
from dnakit.epigenome import load_hmmratac_model

result = load_hmmratac_model("macs3_model.json", verbose=False)
print(result["hmm_type"], result["binsize"], result["state_labels"])
print(result["open_state"], result["nucleosomal_state"], result["background_state"])
print(result["start_probabilities"])
print(result["means"][2], result["covariances"][1][0])
```

- **示例输出：**

```text
gaussian 10 ['bg', 'nuc', 'open']
2 1 0
[0.2, 0.3, 0.5]
[10.0, 10.0, 10.0, 10.0] [4.0, 0.0, 0.0, 0.0]
```

### 6) `dnakit.epigenome.predict_hmmratac_states()`

- **作用：** 用模型文件计算每个 bin 属于开放、核小体、背景状态的后验概率。
- **处理方法：** 先用 `hmm_model_init` 读取模型，再调用 `HMMR_HMM.hmm_predict`（前向-后向算法的 `predict_proba`）。返回值的列按开放、核小体、背景排列，`labels` 取概率最大的状态（并列时开放优先，其次核小体）。同时给出 `bins` 与 `output_path` 时，写出 MACS3 流程使用的预测概率文件：`染色体,终点,概率,概率,概率`，三个概率按模型自己的状态顺序、保留 6 位小数，可直接交给后面三个函数。
- **API：** `dnakit.epigenome.predict_hmmratac_states(model_path: 'str | Path', data: 'Sequence[Sequence[float]]', lengths: 'Sequence[int]', *, bins: 'Sequence[Sequence[Any]] | None' = None, output_path: 'str | Path | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **输入：** `model_path`；`data` 与 `lengths` 同训练函数（Poisson 模型需要整数）；`bins` 为每行一个 `[染色体, 终点]`；`output_path` 为要写入的预测概率文件，不覆盖已有文件。`bins` 与 `output_path` 要么同时给出，要么都不给。
- **输出：** `probabilities`、`columns`、`labels`、`label_counts`、`rows`、`binsize`、三个状态编号；写文件时还有 `output`。示例中两行数据以单个 bin 作序列，后验就是起始概率乘以发射密度再归一化，第一行属背景，第二行属开放。
- **示例代码：**

```python
import tempfile

from dnakit.epigenome import predict_hmmratac_states

data = [[2.0, 2.0, 2.0, 2.0], [10.0, 10.0, 10.0, 10.0]]
with tempfile.TemporaryDirectory() as work:
    result = predict_hmmratac_states(
        "macs3_model.json",
        data,
        [1, 1],
        bins=[["chr1", 10], ["chr1", 20]],
        output_path=f"{work}/proba.csv",
        verbose=False,
    )
    with open(result["output"], encoding="utf-8") as handle:
        print(handle.read().strip())
print(result["columns"])
for row, label in zip(result["probabilities"], result["labels"]):
    print([round(value, 4) for value in row], label)
print(result["label_counts"])
```

- **示例输出：**

```text
chr1,10,0.992363,0.007637,0.000000
chr1,20,0.000000,0.000011,0.999989
['open', 'nucleosomal', 'background']
[0.0, 0.0076, 0.9924] bg
[1.0, 0.0, 0.0] open
{'open': 1, 'nuc': 0, 'bg': 1}
```

### 7) `dnakit.epigenome.hmmratac_proba_to_bedgraph()`

- **作用：** 把预测概率文件转成三条 bedGraph 轨道：每个 bin 属于开放、核小体、背景状态的概率。
- **处理方法：** 先检查文件（每行 5 个字段、概率非负、同一染色体连续且 bin 不重叠），再调用 `hmmratac_cmd.save_proba_to_bedGraph`。bin 之间的缺口以及每条染色体第一个 bin 之前的部分记为背景（开放 0、核小体 0、背景 1）；相邻且值相同的区间合并。bedGraph 区间是 0 起点、右开。
- **API：** `dnakit.epigenome.hmmratac_proba_to_bedgraph(proba_path: 'str | Path', output_prefix: 'str | Path', *, open_state: 'int', nucleosomal_state: 'int', background_state: 'int', binsize: 'int' = 10, verbose: bool = True) -> 'dict[str, Any]'`.
- **输入：** `proba_path` 为预测概率文件；`output_prefix` 为输出文件名前缀，写出 `<前缀>_open.bdg`、`<前缀>_nuc.bdg`、`<前缀>_bg.bdg`（已存在则报错）；`open_state`、`nucleosomal_state`、`background_state` 为三个状态在文件中的列编号（0、1、2 的一个排列，由训练或加载函数给出）；`binsize`（默认 10）。
- **输出：** `files`（三个文件路径）、`intervals`（每个文件的行数）、`bin_count`、`chromosomes`。示例中状态 2 是开放，所以开放轨道取第三列的概率。
- **示例代码：**

```python
import tempfile

from dnakit.epigenome import hmmratac_proba_to_bedgraph

with tempfile.TemporaryDirectory() as work:
    result = hmmratac_proba_to_bedgraph(
        "macs3_proba.csv",
        f"{work}/sample",
        open_state=2,
        nucleosomal_state=0,
        background_state=1,
        binsize=10,
        verbose=False,
    )
    print(result["bin_count"], result["chromosomes"], result["intervals"])
    with open(result["files"]["open"], encoding="utf-8") as handle:
        print(handle.read().strip())
```

- **示例输出：**

```text
13 ['chr1', 'chr2'] {'open': 12, 'nucleosomal': 15, 'background': 9}
chr1	0	10	0.00000
chr1	10	30	0.10000
chr1	30	50	0.80000
chr1	50	70	0.10000
chr1	70	90	0.00000
chr1	90	100	0.10000
chr1	100	110	0.80000
chr1	110	120	0.10000
chr2	0	30	0.00000
chr2	30	50	0.10000
chr2	50	60	0.80000
chr2	60	70	0.10000
```

### 8) `dnakit.epigenome.hmmratac_states_path()`

- **作用：** 取每个 bin 概率最大的状态，把相邻同状态的 bin 合并，得到状态路径。
- **处理方法：** 检查文件后调用 `hmmratac_cmd.generate_states_path`。概率并列时按开放、核小体、背景的顺序取第一个。bin 之间的缺口和每条染色体第一个 bin 之前的部分记为背景；缺口前后的背景片段不会合并。
- **API：** `dnakit.epigenome.hmmratac_states_path(proba_path: 'str | Path', *, open_state: 'int', nucleosomal_state: 'int', background_state: 'int', binsize: 'int' = 10, verbose: bool = True) -> 'dict[str, Any]'`.
- **输入：** `proba_path`、`open_state`、`nucleosomal_state`、`background_state`、`binsize` 同上一个函数。
- **输出：** `states` 为 `{chrom, start, end, state}` 的列表，`state` 是 `open`、`nuc` 或 `bg`；`state_counts` 与 `state_bases` 按状态给出片段数与碱基数；`bin_count`、`chromosomes`。可直接交给下面两个函数。
- **示例代码：**

```python
from dnakit.epigenome import hmmratac_states_path

result = hmmratac_states_path(
    "macs3_proba.csv",
    open_state=2,
    nucleosomal_state=0,
    background_state=1,
    binsize=10,
    verbose=False,
)
print(result["state_counts"], result["state_bases"])
for item in result["states"]:
    print(item["chrom"], item["start"], item["end"], item["state"])
```

- **示例输出：**

```text
{'open': 3, 'nuc': 6, 'bg': 6} {'open': 40, 'nuc': 60, 'bg': 90}
chr1 0 10 bg
chr1 10 20 bg
chr1 20 30 nuc
chr1 30 50 open
chr1 50 60 nuc
chr1 60 70 bg
chr1 70 90 bg
chr1 90 100 nuc
chr1 100 110 open
chr1 110 120 nuc
chr2 0 30 bg
chr2 30 40 bg
chr2 40 50 nuc
chr2 50 60 open
chr2 60 70 nuc
```

### 9) `dnakit.epigenome.write_hmmratac_states_bed()`

- **作用：** 把状态路径里的开放和核小体片段写成 BED 文件。
- **处理方法：** 调用 `hmmratac_cmd.save_states_bed`：每个非背景片段写一行 `染色体 起点 终点 状态`。输入可以是 `hmmratac_states_path` 返回的字典列表，也可以是 `(染色体, 起点, 终点, 状态)` 四项列表；片段需沿染色体排序且不重叠。不覆盖已有文件。
- **API：** `dnakit.epigenome.write_hmmratac_states_bed(states: 'Sequence[Mapping[str, Any] | Sequence[Any]]', output_path: 'str | Path', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **输入：** `states`；`output_path` 为要写入的 BED 文件。
- **输出：** `output`、`records`（写出的行数）、`state_counts`、`skipped_background`（未写出的背景片段数）。
- **示例代码：**

```python
import tempfile

from dnakit.epigenome import hmmratac_states_path, write_hmmratac_states_bed

path = hmmratac_states_path(
    "macs3_proba.csv", open_state=2, nucleosomal_state=0, background_state=1, verbose=False
)
with tempfile.TemporaryDirectory() as work:
    result = write_hmmratac_states_bed(path["states"], f"{work}/states.bed", verbose=False)
    with open(result["output"], encoding="utf-8") as handle:
        print(handle.read().strip())
print(result["records"], result["state_counts"], result["skipped_background"])
```

- **示例输出：**

```text
chr1	20	30	nuc
chr1	30	50	open
chr1	50	60	nuc
chr1	90	100	nuc
chr1	100	110	open
chr1	110	120	nuc
chr2	40	50	nuc
chr2	50	60	open
chr2	60	70	nuc
9 {'open': 3, 'nuc': 6} 6
```

### 10) `dnakit.epigenome.write_hmmratac_accessible_regions()`

- **作用：** 把状态路径里“核小体-开放-核小体”结构中的开放片段写成 narrowPeak 文件，并用一条 bedGraph 给区域打分。
- **处理方法：** 先检查 bedGraph（区间排序且首尾相接），再调用 `hmmratac_cmd.save_accessible_regions`。开放片段两侧必须紧邻核小体片段，三段合起来的长度要大于 `min_length`；区域内 bedGraph 的最高值乘以 10 取整作为分数，峰顶取最高值所在区间的中点，fc、pscore、qscore 写 0。MACS3 流程用折叠倍数轨道打分，开放态概率轨道同样可用。注意：MACS3 3.0.5 的 `save_accessible_regions` 把结果列表切成 `[:-1]`，状态列表中最后一个符合条件的区域不会写出（示例中有 3 个符合条件的区域，只写出 2 个），本函数如实返回上游的结果。
- **API：** `dnakit.epigenome.write_hmmratac_accessible_regions(states: 'Sequence[Mapping[str, Any] | Sequence[Any]]', score_bedgraph: 'str | Path', output_path: 'str | Path', *, min_length: 'int' = 100, verbose: bool = True) -> 'dict[str, Any]'`.
- **输入：** `states`；`score_bedgraph` 为评分用 bedGraph；`output_path` 为要写入的 narrowPeak 文件；`min_length`（默认 100，不小于 0）。
- **输出：** `regions` 为 `{chrom, start, end, name, score, summit}` 的列表，`count`、`output`、`min_length`、`score_bedgraph`。
- **示例代码：**

```python
import tempfile

from dnakit.epigenome import (
    hmmratac_states_path,
    write_hmmratac_accessible_regions,
)

path = hmmratac_states_path(
    "macs3_proba.csv", open_state=2, nucleosomal_state=0, background_state=1, verbose=False
)
with tempfile.TemporaryDirectory() as work:
    result = write_hmmratac_accessible_regions(
        path["states"], "macs3_open.bdg", f"{work}/regions.narrowPeak", min_length=25, verbose=False
    )
    with open(result["output"], encoding="utf-8") as handle:
        print(handle.read().strip())
print(result["count"])
for region in result["regions"]:
    print(region["name"], region["chrom"], region["start"], region["end"], region["score"], region["summit"])
```

- **示例输出：**

```text
chr1	30	50	MACS_peak_1	8	.	0	0	0	10
chr1	100	110	MACS_peak_2	8	.	0	0	0	5
2
MACS_peak_1 chr1 30 50 8 40
MACS_peak_2 chr1 100 110 8 105
```

### 未封装的函数

以下 10 个台账条目没有新增包装函数，状态与理由逐条如下（台账中的理由为英文）。其中 6 个是 `hmmratac_cmd` 中与 `HMMR_HMM`、`HMMR_Signal_Processing` 里同一对象的重复名称，已用恒等检查核对。

- `MACS3.Commands.hmmratac_cmd.generate_digested_signals`（`duplicate-of` → `digest_hmmratac_fragments`）：与 `HMMR_Signal_Processing.generate_digested_signals` 是同一个对象（已核对），由 `digest_hmmratac_fragments` 覆盖。
- `MACS3.Commands.hmmratac_cmd.extract_signals_from_regions`（`duplicate-of` → `extract_hmmratac_training_data`）：与 `HMMR_Signal_Processing.extract_signals_from_regions` 是同一个对象（已核对），由 `extract_hmmratac_training_data` 覆盖。
- `MACS3.Commands.hmmratac_cmd.hmm_training`（`duplicate-of` → `train_hmmratac_model`）：与 `HMMR_HMM.hmm_training` 是同一个对象（已核对），由 `train_hmmratac_model` 覆盖。
- `MACS3.Commands.hmmratac_cmd.hmm_predict`（`duplicate-of` → `predict_hmmratac_states`）：与 `HMMR_HMM.hmm_predict` 是同一个对象（已核对），由 `predict_hmmratac_states` 覆盖。
- `MACS3.Commands.hmmratac_cmd.hmm_model_init`（`duplicate-of` → `load_hmmratac_model`）：与 `HMMR_HMM.hmm_model_init` 是同一个对象（已核对），由 `load_hmmratac_model` 覆盖。
- `MACS3.Commands.hmmratac_cmd.hmm_model_save`（`duplicate-of` → `save_hmmratac_model`）：与 `HMMR_HMM.hmm_model_save` 是同一个对象（已核对），由 `save_hmmratac_model` 覆盖。
- `MACS3.IO.OutputWriter.diag_write`（`internal`）：callpeak 诊断表的报告写出函数，不含序列或区间逻辑。
- `MACS3.Signal.PeakModel.smooth`（`not-dna`）：对任意数组做窗口平滑，与 DNA 无关。
- `MACS3.Signal.SignalProcessing.savitzky_golay_order2_deriv1`（`not-dna`）：通用的 Savitzky-Golay 导数滤波，与 DNA 无关。
- `MACS3.Signal.SignalProcessing.savitzky_golay`（`not-dna`）：通用的 Savitzky-Golay 平滑滤波，与 DNA 无关。

## deepTools 补充函数（GC 偏差表、合成指纹 JSD 与 bigWig 瓦片 bedGraph）

这一节补充 deepTools 的 4 个函数，位于 `dnakit.epigenome` 的新模块 `deeptools_review_functions.py`，示例从该模块文件导入；导出到 `dnakit.epigenome` 后同样可以用 `dnakit.epigenome.<名称>()` 调用。它们直接调用 deepTools 的 Python 函数，不运行命令行，全部单进程运行，因此结果不依赖进程调度。

有两处上游细节决定了包装方式。其一，`computeGCBias.tabulateGCcontent` 与 `countReadsPerGC` 读取模块级字典 `global_vars`，只有命令行才会设置它；包装函数按命令行的方式构造这个字典，只在调用期间安装，调用结束（包括出错时）后恢复原来的状态。其二，`writeBedGraph_bam_and_bw.writeBedGraph` 的 BAM 分支调用 `getCoverageFromBam` 和 `getSmoothRange`，deepTools 4.0.0 没有定义这两个函数，BAM 输入会触发 `NameError`；只有 bigWig 分支能工作，所以对应的包装函数接受 bigWig 文件，名为 `bigwig_bedgraph_tiles`。deepTools 在输入不合法时会调用 `sys.exit`，这些退出被转换为 `ConfigurationError`；包装函数先检查输入，不覆盖已有的输出文件。返回值是带 `status`、`operation` 和 `provenance` 的字典，坐标是 0 起点、右开。

示例数据是手工构造的，用来说明格式和计算，不代表真实实验。示例的输出文件写入临时目录，因此可以重复运行。

### 示例数据

演示用的小工具，只用于生成示例的输入文件，不属于 DNAKit：`write_twobit` 从只含 ACGT 的 FASTA 写出 2bit 文件，`sam_to_bam` 把 SAM 转成排序并建索引的 BAM。

**示例文件：** `deeptools_fixtures.py`

```python
import os
import struct

import pysam


def write_twobit(fasta_path, twobit_path):
    """Write a one-chromosome 2bit file from a FASTA file of plain ACGT letters."""
    name, *rest = open(fasta_path).read().split()
    sequence = "".join(rest)
    codes = {"T": 0, "C": 1, "A": 2, "G": 3}
    packed = bytearray((len(sequence) + 3) // 4)
    for index, base in enumerate(sequence):
        packed[index // 4] |= codes[base] << (6 - 2 * (index % 4))
    record = struct.pack("<IIII", len(sequence), 0, 0, 0) + bytes(packed)
    raw = name[1:].encode()
    header = struct.pack("<IIII", 0x1A412743, 0, 1, 0) + bytes([len(raw)]) + raw
    header += struct.pack("<I", len(header) + 4)
    with open(twobit_path, "wb") as handle:
        handle.write(header + record)


def sam_to_bam(sam_path, bam_path):
    """Convert a SAM file to a sorted and indexed BAM file."""
    pysam.view("-b", "-o", bam_path + ".unsorted", sam_path, catch_stdout=False)
    pysam.sort("-o", bam_path, bam_path + ".unsorted")
    pysam.index(bam_path)
    os.remove(bam_path + ".unsorted")
```

40 bp 的染色体，由 10 个 4 碱基块拼成，各块的 GC 数依次为 2、0、4、0、2、0、3、1、4、1。

**示例文件：** `deeptools_genome.fa`

```text
>chr1
ACGTAAAAGGCCATATGCATTTTTCCCAAGATGCGCACTA
```

17 条 5 bp 的 read，14 条正链、3 条反链。后面按 4 bp 取样的位置是 0、4、8、……、36；起点落在这些位置上的正链 read 有 11 条（位置 4 有 2 条、位置 24 有 3 条、位置 32 有 2 条，位置 0、8、16、28 各 1 条）。起点为 2、5、17 的正链 read 不在取样位置上，3 条反链 read 不计入正链起点。

**示例文件：** `deeptools_reads.sam`

```text
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:40
f0	0	chr1	1	60	5M	*	0	0	AAAAA	IIIII
f11	0	chr1	3	60	5M	*	0	0	AAAAA	IIIII
f1	0	chr1	5	60	5M	*	0	0	AAAAA	IIIII
f2	0	chr1	5	60	5M	*	0	0	AAAAA	IIIII
f12	0	chr1	6	60	5M	*	0	0	AAAAA	IIIII
r0	16	chr1	7	60	5M	*	0	0	AAAAA	IIIII
f3	0	chr1	9	60	5M	*	0	0	AAAAA	IIIII
f4	0	chr1	17	60	5M	*	0	0	AAAAA	IIIII
f13	0	chr1	18	60	5M	*	0	0	AAAAA	IIIII
r1	16	chr1	21	60	5M	*	0	0	AAAAA	IIIII
r2	16	chr1	22	60	5M	*	0	0	AAAAA	IIIII
f5	0	chr1	25	60	5M	*	0	0	AAAAA	IIIII
f6	0	chr1	25	60	5M	*	0	0	AAAAA	IIIII
f7	0	chr1	25	60	5M	*	0	0	AAAAA	IIIII
f8	0	chr1	29	60	5M	*	0	0	AAAAA	IIIII
f10	0	chr1	33	60	5M	*	0	0	AAAAA	IIIII
f9	0	chr1	33	60	5M	*	0	0	AAAAA	IIIII
```

两个 bigWig 文件的来源 bedGraph（示例中用 `bedgraph_to_bigwig` 转成 bigWig）。第一个在 chr1 上是 2（0–50）和 4（50–100），在 chr2 上只覆盖 0–30，值为 6；第二个在 chr1 上是 1（0–25）和 3（25–100），在 chr2 上整条为 1。

**示例文件：** `deeptools_a.bdg`

```text
chr1	0	50	2
chr1	50	100	4
chr2	0	30	6
```

**示例文件：** `deeptools_b.bdg`

```text
chr1	0	25	1
chr1	25	100	3
chr2	0	60	1
```

### 1) `dnakit.epigenome.tabulate_gc_fragment_counts()`

- **作用：** 沿基因组每隔 `step_size` 个位置取一个长度为 `fragment_length` 的片段，统计它的 GC 碱基数，同时统计起点落在这些位置上的正链 read 数，得到 `computeGCBias` 的频率表（`correctGCBias` 的输入）：每个 GC 数一行。
- **处理方法：** 调用 deepTools 的 `computeGCBias.tabulateGCcontent`。包装函数先读 2bit 的染色体长度和 BAM 索引中的比对 read 数，像命令行一样构造 `global_vars`：`reads_per_bp` 是比对 read 数除以 `effective_genome_size`；`max_reads` 是均值为 `4 × reads_per_bp × fragment_length` 的泊松分布在概率 `1 / sample_size` 处的上分位数，起点处正链 read 数达到 `max_reads` 的位置被当作峰跳过（仍计入取样位置数）。默认步长是 `max(int(2bit 总长 / sample_size), 1)`。2bit 与 BAM 的染色体名用 `tbitToBamChrName` 对应；BAM 缺少 `chr` 前缀时也能对上，而 `computeGCBias` 命令行在这种情况下会因 `UnboundLocalError` 失败。含 N 超过 5% 的片段不计入，只统计正链 read。
- **API：** `dnakit.epigenome.tabulate_gc_fragment_counts(twobit_path: 'str | Path', bam_path: 'str | Path', fragment_length: 'int', effective_genome_size: 'int', *, sample_size: 'int' = 50000000, step_size: 'int | None' = None, region: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **输入：** `twobit_path`；`bam_path`（需要 `.bai` 索引、至少有一条比对 read）；`fragment_length`（1 到 1000000）；`effective_genome_size`（可比对的基因组大小）；`sample_size`（默认 50000000，至少 2，决定默认步长和 `max_reads`）；`step_size`（默认由 `sample_size` 算出）；`region`（`染色体:起点:终点`，染色体名用 BAM 中的名字，只取样该区域）。取样位置总数不超过 200000000。
- **输出：** `gc_table` 每个 GC 数一项 `{gc_bases, read_starts, genome_positions, read_ratio}`：`genome_positions` 是该 GC 数的取样位置数，`read_starts` 是这些位置上的正链 read 起点数，`read_ratio` 是二者之比再除以全部取样位置上的同一个比值，即相对于整体密度的偏差；还有 `positions_sampled`、`read_starts_counted`、`step_size`、`mapped_reads`、`reads_per_bp`、`max_reads` 和 `min_reads`。示例中 10 个取样位置的 GC 数 0 到 4 分别有 3、2、2、1、2 个，正链 read 起点数为 2、1、2、3、3，所以 GC 数为 3 的比值是 (3 / 1) ÷ (11 / 10) = 2.727。
- **示例代码：**

```python
from deeptools_fixtures import sam_to_bam, write_twobit
from dnakit.epigenome import tabulate_gc_fragment_counts

write_twobit("deeptools_genome.fa", "deeptools_genome.2bit")
sam_to_bam("deeptools_reads.sam", "deeptools_reads.bam")
result = tabulate_gc_fragment_counts(
    "deeptools_genome.2bit",
    "deeptools_reads.bam",
    4,
    40,
    sample_size=10,
    verbose=False,
)
print(result["step_size"], result["mapped_reads"], result["max_reads"])
print(result["positions_sampled"], result["read_starts_counted"])
for row in result["gc_table"]:
    print(
        row["gc_bases"],
        row["genome_positions"],
        row["read_starts"],
        round(row["read_ratio"], 3),
    )
```

- **示例输出：**

```text
4 17 10.0
10 11
0 3 2 0.606
1 2 1 0.455
2 2 2 0.909
3 1 3 2.727
4 2 3 1.364
```

### 2) `dnakit.epigenome.count_reads_per_gc_bins()`

- **作用：** 沿基因组取长度为 `region_size` 的窗口（每 `step_size` 个位置一个），给出每个窗口的 GC 比例和与它重叠的 read 数（两条链都计），这是 `computeGCBias --biasPlot` 画图用的数据。
- **处理方法：** 调用 deepTools 的 `computeGCBias.countReadsPerGC`，`global_vars` 的处理与上一个函数相同。默认步长是 `max(int(2bit 总长 / sample_size), 1) × 10`，与命令行画图时一致。超出染色体末端的窗口不产生，含 N 超过 5% 的窗口被略去。窗口按基因组顺序排列；上游只返回 read 数和 GC，不返回窗口坐标。
- **API：** `dnakit.epigenome.count_reads_per_gc_bins(twobit_path: 'str | Path', bam_path: 'str | Path', effective_genome_size: 'int', *, region_size: 'int' = 300, sample_size: 'int' = 50000000, step_size: 'int | None' = None, region: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **输入：** `twobit_path`；`bam_path`；`effective_genome_size`（只影响分块大小）；`region_size`（默认 300）；`sample_size`（默认 50000000）；`step_size`；`region`（同上）。窗口总数不超过 1000000。
- **输出：** `windows` 是 `{reads, gc_fraction}` 的列表，`window_count` 是窗口数，`mean_reads` 是平均 read 数（没有窗口时为 `None`）。示例的染色体 40 bp、窗口 8 bp、步长 4，起点 0 到 32 共 9 个窗口。第一个窗口是 `ACGTAAAA`，GC 为 2/8 = 0.25；起点为 0、2、4、4、5 的正链 read 和起点为 6 的反链 read 与它重叠，共 6 条。
- **示例代码：**

```python
from deeptools_fixtures import sam_to_bam, write_twobit
from dnakit.epigenome import count_reads_per_gc_bins

write_twobit("deeptools_genome.fa", "deeptools_genome.2bit")
sam_to_bam("deeptools_reads.sam", "deeptools_reads.bam")
result = count_reads_per_gc_bins(
    "deeptools_genome.2bit",
    "deeptools_reads.bam",
    40,
    region_size=8,
    sample_size=40,
    step_size=4,
    verbose=False,
)
print(result["window_count"], round(result["mean_reads"], 3))
for window in result["windows"]:
    print(window["reads"], window["gc_fraction"])
```

- **示例输出：**

```text
9 5.222
6 0.25
7 0.5
5 0.5
3 0.25
4 0.25
7 0.375
6 0.5
6 0.625
3 0.625
```

### 3) `dnakit.epigenome.fingerprint_synthetic_jsd()`

- **作用：** 衡量一份每个 bin 的覆盖度向量偏离同均值泊松分布的程度，即 `plotFingerprint` 在没有对照样本时使用的合成 JSD；值越大，覆盖越集中在少数 bin。
- **处理方法：** 调用 deepTools 的 `plotFingerprint.getSyntheticJSD`：把 bin 覆盖度（取整数部分，覆盖为 0 的 bin 舍去，不小于 `MAXLEN` = 10000000 的值计入最后一格）统计成直方图；以覆盖度均值为参数、总覆盖度为规模生成泊松分布作为对照；再用 `getJSDcommon` 计算两条累积覆盖曲线之间的 Jensen-Shannon 距离。没有随机成分，结果是确定的。一次调用要处理 1000 万长的数组，约需 3 到 4 秒。
- **API：** `dnakit.epigenome.fingerprint_synthetic_jsd(coverage: 'Sequence[float]', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **输入：** `coverage` 是非负数的列表（1 到 10000000 个值，至少一个不小于 1）。布尔值、负数、NaN 和无穷大会被拒绝。
- **输出：** `bins`、`mean_coverage`、`bins_above_maximum`（不小于 `MAXLEN` 的 bin 数）和 `jsd`。示例中接近泊松的 `even` 距离为 0.209，集中在 2 个 bin 里的 `peaked` 距离为 0.292。
- **示例代码：**

```python
from dnakit.epigenome import fingerprint_synthetic_jsd

even = [4, 5, 3, 4, 6, 4, 5, 3, 4, 5]
peaked = [0, 0, 1, 12, 15, 1, 0, 0, 0, 1]
for coverage in (even, peaked):
    result = fingerprint_synthetic_jsd(coverage, verbose=False)
    print(result["bins"], result["mean_coverage"], result["jsd"])
```

- **示例输出：**

```text
10 4.3 0.209461503453211
10 3.0 0.292041559465403
```

### 4) `dnakit.epigenome.bigwig_bedgraph_tiles()`

- **作用：** 把一个或多个 bigWig 文件按固定大小的瓦片取平均，按缩放系数求平均或求和，写出 bedGraph；可用来做 bigWig 的平均、差值或缩放。
- **处理方法：** 调用 deepTools 的 `writeBedGraph_bam_and_bw.writeBedGraph`（`bigwigAverage` 与 `bigwigCompare` 使用的引擎）。每个文件在每个瓦片内取碱基值的平均（`missing_data_as_zero=True` 时缺失值当作 0，否则含缺失值的瓦片为 NaN，不写出）；各文件的瓦片值乘以对应的 `scale_factors`，再按 `reduction` 取平均（`"mean"`）或求和（`"sum"`）；相邻且值相等的瓦片合并成一行，`fixed_step=True` 时每个瓦片一行。上游的 BAM 分支在 deepTools 4.0.0 中会触发 `NameError`，所以只接受 bigWig 文件，BAM 输入会报 `ConfigurationError`。区域的起点向下取整到瓦片边界；染色体末尾值为 0 的瓦片上游不写出。多个文件的同名染色体需要等长，只处理共有的染色体。
- **API：** `dnakit.epigenome.bigwig_bedgraph_tiles(bigwig_paths: 'Sequence[str | Path]', output_path: 'str | Path', *, tile_size: 'int' = 50, reduction: 'str' = 'mean', scale_factors: 'Sequence[float] | None' = None, region: 'str | None' = None, missing_data_as_zero: 'bool' = True, fixed_step: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **输入：** `bigwig_paths`（1 到 32 个 bigWig 文件）；`output_path`（必须不存在，不会覆盖）；`tile_size`（默认 50）；`reduction`（`"mean"` 或 `"sum"`）；`scale_factors`（每个文件一个数，默认全为 1）；`region`（`染色体:起点:终点`）；`missing_data_as_zero`（默认 `True`）；`fixed_step`（默认 `False`）。
- **输出：** `output`（写出的 bedGraph 路径）、`intervals`（行数）、`chromosomes`、`first_intervals`（前 20 行的 `{chrom, start, end, value}`）和所用参数。示例中 chr1 的第一个瓦片 [0, 25) 在两个文件中是 2 和 1，平均值 1.5；chr2 的瓦片 [25, 50) 在第一个文件里只有 5 个碱基有值 6，平均 1.2，与第二个文件的 1 平均得 1.1；第二次调用用 `reduction="sum"` 和 `scale_factors=[1, -1]` 得到两个文件之差。
- **示例代码：**

```python
import tempfile

from dnakit.epigenome import bedgraph_to_bigwig
from dnakit.epigenome import bigwig_bedgraph_tiles

sizes = {"chr1": 100, "chr2": 60}
with tempfile.TemporaryDirectory() as work:
    bedgraph_to_bigwig(sizes, ["deeptools_a.bdg"], f"{work}/a.bw", verbose=False)
    bedgraph_to_bigwig(sizes, ["deeptools_b.bdg"], f"{work}/b.bw", verbose=False)
    tracks = [f"{work}/a.bw", f"{work}/b.bw"]
    mean = bigwig_bedgraph_tiles(
        tracks, f"{work}/mean.bdg", tile_size=25, verbose=False
    )
    print(open(mean["output"]).read().strip())
    print(mean["intervals"], mean["chromosomes"])
    difference = bigwig_bedgraph_tiles(
        tracks,
        f"{work}/difference.bdg",
        tile_size=25,
        reduction="sum",
        scale_factors=[1, -1],
        verbose=False,
    )
    print(open(difference["output"]).read().strip())
```

- **示例输出：**

```text
chr1	0	25	1.5
chr1	25	50	2.5
chr1	50	100	3.5
chr2	0	25	3.5
chr2	25	50	1.1
chr2	50	60	0.5
6 ['chr1', 'chr2']
chr1	0	25	1
chr1	25	50	-1
chr1	50	100	1
chr2	0	25	5
chr2	25	50	0.2
chr2	50	60	-1
```

### 未封装的函数

以下 6 个台账条目没有新增包装函数，状态与理由逐条如下（台账中的理由为英文）。

- `deeptools.bamPEFragmentSize.getDensity`（`plot`）：用 matplotlib 的 `plt.hist` 和 `plt.clf` 画片段长度直方图的绘图辅助函数。
- `deeptools.bigwigAverage.getType`（`internal`）：按文件扩展名猜测轨道文件类型。
- `deeptools.bigwigCompare.getType`（`internal`）：与 `bigwigAverage.getType` 相同的扩展名判断。
- `deeptools.computeGCBias.smooth`（`internal`）：窗口平均的平滑辅助函数，上游注明目前未使用。
- `deeptools.heatmapper.trimZones`（`internal`）：热图分区的记账辅助函数，服务于绘图分箱。
- `deeptools.heatmapper.computeSilhouetteScore`（`not-dna`）：对距离矩阵的聚类统计，与 DNA 无关。
