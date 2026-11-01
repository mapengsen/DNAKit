# BEDOPS bedmap：重叠元素的统计与列表

`dnakit.annotation`（模块 `bedops_p05_cli_functions`）中的两个函数封装 [BEDOPS](https://github.com/bedops/bedops) 2.4.42 的 `bedmap`：它遍历参考 BED 文件，对每个元素概括映射文件中与之重叠的元素。`bedmap_score_statistic` 为每个参考行返回一个数值（得分统计量、数量、碱基数），`bedmap_echo_elements` 输出参考行或重叠元素、ID、得分、长度和范围。bedmap 的操作由参数选择，因此一个函数覆盖多个 `bedmap --<operation>` 开关。两个输入文件都必须按 `sort-bed` 顺序排序。

BEDOPS 是外部程序：先安装（例如 `conda install -c bioconda bedops`，或从 https://github.com/bedops/bedops/releases 下载发布包），再把其 `bin` 目录放入 PATH 或传入 `bedops_bin_dir`。找不到程序时抛出带安装提示的 `BackendUnavailableError`；程序出错退出时抛出 `BackendExecutionError`（代码 `BEDOPS_EXECUTION_FAILED`）并附日志末尾，且不会留下不完整的输出。BEDOPS 采用 GPL-2.0 许可，DNAKit 不内置它。

所有结果共有这些字段：`status`（`completed`）、`backend`（`bedops`）、`operation`、`program`、`version`（从程序读取）、`command`、`outputs`（写出文件的路径）、`log_tail` 和 `provenance`（工具、版本、测试版本 2.4.42、可执行文件、引用文献）；下面各函数小节只列出额外字段。

## 示例数据

示例从当前目录读取下面的文件，并把输出写到临时目录。

**示例文件：** `ref.bed`

```text
chr1	100	200	r1
chr1	300	400	r2
chr1	1000	1100	r3
chr2	10	60	r4
```

**示例文件：** `map.bed`

```text
chr1	150	250	m1	2	+
chr1	180	350	m2	4	-
chr1	320	330	m3	6	+
chr2	0	20	m4	10	+
```

## 1) `dnakit.annotation.bedmap_score_statistic()`

- **作用：** 对每个参考元素，用一个数值概括与之重叠的映射元素：得分统计量、数量或碱基数。
- **处理方法：** 运行 `bedmap [flags] [overlap option] --<statistic> <reference> [<map>]`，把单列结果写入 `output_path`（每个参考行一行，没有重叠时为 `NAN`）。`statistic` 取 `mean`、`median`、`sum`、`min`、`max`、`variance`、`stdev`、`cv`（标准差/均值）、`mad`（中位数绝对偏差，按 `mad_multiplier` 缩放）、`kth`（广义中位数，即 `kth_fraction` 分位数）、`tmean`（截尾均值，忽略最低 `trim_low` 和最高 `trim_high` 比例）、`wmean`（按参考元素被覆盖比例加权的均值），以及重叠度量 `count`、`bases`、`bases-uniq`、`bases-uniq-f` 和 `indicator`。得分统计量读取映射文件的第 5 列。重叠由 `overlap_mode`（`bp-ovr`、`range`、`exact`、`fraction-ref`、`fraction-map`、`fraction-both`、`fraction-either`）和 `overlap_value` 决定，BEDOPS 默认 `bp-ovr 1`。该函数覆盖工作清单中的：重叠得分的变异系数、广义中位数、中位数绝对偏差、最大值、均值、中位数、最小值、标准差、总和、截尾均值、方差和覆盖加权均值，以及重叠碱基总数和重叠元素数。
- **API：** `bedmap_score_statistic(reference_path: 'str', output_path: 'str', *, map_path: 'str | None' = None, statistic: 'str' = 'mean', kth_fraction: 'float | None' = None, mad_multiplier: 'float | None' = None, trim_low: 'float | None' = None, trim_high: 'float | None' = None, overlap_mode: 'str | None' = None, overlap_value: 'float | None' = None, delimiter: 'str | None' = None, precision: 'int | None' = None, scientific: 'bool' = False, skip_unmapped: 'bool' = False, unmapped_value: 'str | None' = None, chromosome: 'str | None' = None, error_check: 'bool' = False, header: 'bool' = False, bedops_bin_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要 BEDOPS 程序 `bedmap`（测试版本 2.4.42）位于 PATH 中，或用 `bedops_bin_dir` 指向 BEDOPS 的 `bin` 目录；不需要任何 Python 包。下面的示例输出是在已安装 BEDOPS 2.4.42 的环境中生成的。找不到程序时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** `reference_path` 是已存在的已排序 BED 文件；`map_path`（可选，默认用参考文件本身）是已存在的已排序 BED 文件，得分统计量需要至少 5 列。`output_path` 必须尚不存在。`statistic` 取上述值；`kth` 需要 (0, 1] 内的 `kth_fraction`，`tmean` 需要 `trim_low` 和 `trim_high`（均在 [0, 1] 内且和不超过 1），`mad` 可带正数 `mad_multiplier`，其他统计量传这些参数会被拒绝。`overlap_mode` / `overlap_value`：`bp-ovr`（至少 1）和 `range`（至少 0）用整数，分数模式用 (0, 1] 的数，`exact` 不带值。`delimiter`、`precision`、`scientific`、`skip_unmapped`、`unmapped_value`、`chromosome`、`error_check`、`header` 分别对应 `--delim`、`--prec`、`--sci`、`--skip-unmapped`、`--unmapped-val`、`--chrom`、`--ec`、`--header`。
- **输出：** 通用结果字典，另有 `statistic`、`values`（前 10 行解析为数值，`NAN` 保持文本）、`line_count` 和 `preview`。`outputs` 含 `output`。
- **示例代码：**

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

- **示例输出：**

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

- **作用：** 对每个参考元素，输出参考行本身，或与之重叠的映射元素、ID、得分、长度、范围。
- **处理方法：** 运行 `bedmap [flags] [overlap option] --<field> <reference> [<map>]`，把结果写入 `output_path`。`field` 取 `echo`（参考行）、`echo-map`（所有重叠的映射元素）、`echo-map-id` 和 `echo-map-id-uniq`（第 4 列的 ID，全部或去重）、`echo-map-range`（它们覆盖的坐标范围）、`echo-map-score`（第 5 列得分）、`echo-map-size`（它们的完整长度）、`echo-overlap-size`（重叠部分的长度）、`echo-ref-name`（`chrom:start-end`）、`echo-ref-row-id`（`id-<行号>`）、`echo-ref-size`，以及得分最高或最低的重叠元素 `max-element` 和 `min-element`（`max-element-rand` 和 `min-element-rand` 在并列时随机选择）。同一列中的多个值用 `multi_delimiter` 连接（BEDOPS 默认 `;`）。该函数覆盖工作清单中的：得分最高和最低的重叠元素、回显参考元素、列出重叠映射元素、重叠映射元素的 ID 与去重 ID、坐标范围、得分和长度、重叠长度、`chrom:start-end` 形式的参考元素名称以及参考行标识。
- **API：** `bedmap_echo_elements(reference_path: 'str', output_path: 'str', *, map_path: 'str | None' = None, field: 'str' = 'echo-map', overlap_mode: 'str | None' = None, overlap_value: 'float | None' = None, delimiter: 'str | None' = None, multi_delimiter: 'str | None' = None, precision: 'int | None' = None, scientific: 'bool' = False, skip_unmapped: 'bool' = False, unmapped_value: 'str | None' = None, chromosome: 'str | None' = None, error_check: 'bool' = False, header: 'bool' = False, bedops_bin_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要 BEDOPS 程序 `bedmap`（测试版本 2.4.42）位于 PATH 中，或用 `bedops_bin_dir` 指向 BEDOPS 的 `bin` 目录；不需要任何 Python 包。下面的示例输出是在已安装 BEDOPS 2.4.42 的环境中生成的。找不到程序时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** `reference_path`、`map_path` 和 `output_path` 与 `bedmap_score_statistic` 相同（ID 字段需要映射文件至少 4 列，得分字段需要 5 列）。`field` 取上述值。其余参数与 `bedmap_score_statistic` 相同。
- **输出：** 通用结果字典，另有 `field`、`line_count` 和 `preview`。`outputs` 含 `output`。
- **示例代码：**

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

- **示例输出：**

```text
echo-map-id ['m1;m2', 'm2;m3', '', 'm4']
echo-map-score ['2.000000;4.000000', '4.000000;6.000000', '', '10.000000']
echo-map-range ['chr1\t150\t350', 'chr1\t180\t350', '', 'chr2\t0\t20']
echo-overlap-size ['50;20', '50;10', '', '10']
echo-ref-name ['chr1:100-200', 'chr1:300-400', 'chr1:1000-1100', 'chr2:10-60']
echo-ref-row-id ['id-1', 'id-2', 'id-3', 'id-4']
max-element ['chr1 180 350 m2 4.000000 -', 'chr1 320 330 m3 6.000000 +', 'chr2 0 20 m4 10.000000 +']
```
