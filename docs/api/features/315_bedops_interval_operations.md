# BEDOPS 区间操作：集合运算、bedextract、closest-features、sort-bed

这一组 `dnakit.annotation`（模块 `bedops_p05_cli_functions`）中的函数封装外部工具集 [BEDOPS](https://github.com/bedops/bedops) 2.4.42 中处理已排序 BED 文件的程序：`bedops`（四种集合运算合并为一个函数）、`bedextract`、`closest-features` 以及 `sort-bed --unique/--duplicates`。每个函数先校验参数，用固定参数列表且不经过 shell 运行一个程序，带超时和输出大小限制，把程序的标准输出写入明确指定的新文件（绝不覆盖）并返回字典。BEDOPS 要求输入按 `sort-bed` 顺序排序；这些函数不会替你排序，`error_check=True` 可让 BEDOPS 检查顺序。

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

**示例文件：** `messy.bed`

```text
chr2	5	9	x
chr1	1	2	y
chr1	1	2	y
chr1	7	8	z
```

## 1) `dnakit.annotation.bedops_set_operation()`

- **作用：** 用 BEDOPS 的集合运算 `--element-of`、`--not-element-of`、`--symmdiff` 和 `--everything` 选取或合并 BED 元素。
- **处理方法：** 运行 `bedops [--chrom C] [--range R] [--ec] [--header] --<operation> [overlap] <文件>`，并把标准输出写入 `output_path`。`element-of` 保留与其他文件元素重叠的参考（第一个文件）元素，`not-element-of` 保留不重叠的元素，`symmdiff` 返回区间集合的对称差，`everything` 输出所有元素的并集并保留全部列。`overlap` 是 bp 数（`1`）或占参考元素的百分比（`"50%"`）；省略时 BEDOPS 使用 100%。当查询文件没有完全嵌套的元素时，`bedextract <query> <target>` 与 `overlap=1` 的 `element-of` 结果相同。该函数覆盖工作清单中的：与查询元素重叠的参考元素、与查询元素不重叠的参考元素、区间集合的对称差、保留全部列的区间集合并集。
- **API：** `bedops_set_operation(input_paths: 'Sequence[str]', output_path: 'str', *, operation: 'str' = 'element-of', overlap: 'str | int | None' = None, range: 'str | int | None' = None, chromosome: 'str | None' = None, error_check: 'bool' = False, header: 'bool' = False, bedops_bin_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要 BEDOPS 程序 `bedops`（测试版本 2.4.42）位于 PATH 中，或用 `bedops_bin_dir` 指向 BEDOPS 的 `bin` 目录；不需要任何 Python 包。下面的示例输出是在已安装 BEDOPS 2.4.42 的环境中生成的。找不到程序时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** `input_paths` 是已存在的、已排序的 BED（或 Starch）文件列表，第一个是参考文件。`element-of`、`not-element-of` 和 `symmdiff` 至少需要两个文件，`everything` 至少一个。`output_path` 必须尚不存在。`operation` 取 `element-of`、`not-element-of`、`symmdiff`、`everything` 之一。`overlap`（仅用于 `element-of` 和 `not-element-of`）是整数或 `"50%"` 这样的字符串。`range` 对非参考文件的区间做扩展（`"10"` 或 `"-5:20"`）。`chromosome` 限定只处理一条染色体；`error_check` 加 `--ec`；`header` 加 `--header`。
- **输出：** 通用结果字典，另有 `line_count` 和 `preview`（输出的前 10 行）。`outputs` 含 `output`。
- **示例代码：**

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

- **示例输出：**

```text
element-of 3 chr1 100 200 r1
not-element-of 1 chr1 1000 1100 r3
symmdiff 6 chr1 100 150
everything 8 chr1 100 200 r1
element-of 50%: ['r1', 'r2']
```

## 2) `dnakit.annotation.bedextract_list_chromosomes()`

- **作用：** 列出已排序 BED 或 Starch 文件中的染色体名称。
- **处理方法：** 运行 `bedextract --list-chr <input>` 并读取其标准输出。
- **API：** `bedextract_list_chromosomes(input_path: 'str', *, bedops_bin_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要 BEDOPS 程序 `bedextract`（测试版本 2.4.42）位于 PATH 中，或用 `bedops_bin_dir` 指向 BEDOPS 的 `bin` 目录；不需要任何 Python 包。下面的示例输出是在已安装 BEDOPS 2.4.42 的环境中生成的。找不到程序时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** `input_path` 是已存在的已排序 BED 或 Starch 文件。
- **输出：** 通用结果字典，另有 `chromosomes`（按文件顺序的名称列表）和 `chromosome_count`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.annotation.bedops_p05_cli_functions import bedextract_list_chromosomes

result = bedextract_list_chromosomes("ref.bed", verbose=False)
print(result["chromosomes"], result["chromosome_count"])
```

- **示例输出：**

```text
['chr1', 'chr2'] 2
```

## 3) `dnakit.annotation.bedextract_chromosome()`

- **作用：** 从已排序的 BED 或 Starch 文件中取出一条染色体的全部行。
- **处理方法：** 运行 `bedextract <chromosome> <input>` 并把结果行写入 `output_path`；染色体不存在时得到空文件。
- **API：** `bedextract_chromosome(input_path: 'str', chromosome: 'str', output_path: 'str', *, bedops_bin_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要 BEDOPS 程序 `bedextract`（测试版本 2.4.42）位于 PATH 中，或用 `bedops_bin_dir` 指向 BEDOPS 的 `bin` 目录；不需要任何 Python 包。下面的示例输出是在已安装 BEDOPS 2.4.42 的环境中生成的。找不到程序时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** `input_path` 是已存在的已排序 BED 或 Starch 文件。`chromosome` 是不以 `-` 开头的非空名称。`output_path` 必须尚不存在。
- **输出：** 通用结果字典，另有 `chromosome`、`line_count` 和 `preview`。`outputs` 含 `output`。
- **示例代码：**

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

- **示例输出：**

```text
chr2 1 chr2 10 60 r4
absent chromosome rows: 0
```

## 4) `dnakit.annotation.bedops_closest_features()`

- **作用：** 对 BED 文件的每个元素，报告查询文件中位于其左侧和右侧最近的元素。
- **处理方法：** 运行 `closest-features [flags] <input> <query>`，显式传入 `--delim`（默认 `|`），把结果行写入 `output_path`。默认每行是输入元素、左侧最近元素和右侧最近元素；缺少邻居时输出 `NA`。重叠的查询元素优先于不重叠的元素。`closest` 只保留较近的一侧（距离相同取左侧），`distances` 追加带符号的距离（重叠距离为 0），`no_overlaps` 忽略重叠元素，`no_reference` 和 `no_query` 去掉回显的输入或查询列。
- **API：** `bedops_closest_features(input_path: 'str', query_path: 'str', output_path: 'str', *, closest: 'bool' = False, distances: 'bool' = False, no_overlaps: 'bool' = False, no_reference: 'bool' = False, no_query: 'bool' = False, delimiter: 'str' = '|', chromosome: 'str | None' = None, error_check: 'bool' = False, header: 'bool' = False, bedops_bin_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要 BEDOPS 程序 `closest-features`（测试版本 2.4.42）位于 PATH 中，或用 `bedops_bin_dir` 指向 BEDOPS 的 `bin` 目录；不需要任何 Python 包。下面的示例输出是在已安装 BEDOPS 2.4.42 的环境中生成的。找不到程序时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** `input_path` 和 `query_path` 是已存在的已排序 BED（或 Starch）文件。`output_path` 必须尚不存在。`closest`、`distances`、`no_overlaps`、`no_reference`、`no_query`、`error_check`、`header` 是布尔值。`delimiter` 是不以 `-` 开头的非空字符串（制表符请传真正的制表符）。`chromosome` 限定只处理一条染色体。
- **输出：** 通用结果字典，另有 `line_count` 和 `preview`。`outputs` 含 `output`。
- **示例代码：**

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

- **示例输出：**

```text
chr1 100 200 r1|NA|NA|chr1 180 350 m2 4 -|0
chr1 300 400 r2|chr1 180 350 m2 4 -|0|NA|NA
chr1 1000 1100 r3|chr1 180 350 m2 4 -|-651|NA|NA
chr2 10 60 r4|chr2 0 20 m4 10 +|0|NA|NA
chr1 180 350 m2 4 -
```

## 5) `dnakit.annotation.bedops_sort_bed_filter()`

- **作用：** 对 BED 文件排序，并只保留唯一元素或只保留重复元素。
- **处理方法：** 对输入文件运行 `sort-bed --unique`（类似 `sort -u`）或 `sort-bed --duplicates`（类似 `uniq -d`），把排序后的结果写入 `output_path`。大文件可用 `max_memory` 和 `temp_dir`，分别对应 `--max-mem` 和 `--tmpdir`。
- **API：** `bedops_sort_bed_filter(input_paths: 'Sequence[str]', output_path: 'str', *, mode: 'str' = 'unique', max_memory: 'str | None' = None, temp_dir: 'str | None' = None, bedops_bin_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要 BEDOPS 程序 `sort-bed`（测试版本 2.4.42）位于 PATH 中，或用 `bedops_bin_dir` 指向 BEDOPS 的 `bin` 目录；不需要任何 Python 包。下面的示例输出是在已安装 BEDOPS 2.4.42 的环境中生成的。找不到程序时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** `input_paths` 是非空的已存在 BED 文件列表（不必已排序）。`output_path` 必须尚不存在。`mode` 取 `unique` 或 `duplicates`。`max_memory` 形如 `8G`、`8000M` 或 `8000000000`。`temp_dir` 是已存在的文件夹。
- **输出：** 通用结果字典，另有 `line_count` 和 `preview`。`outputs` 含 `output`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.annotation.bedops_p05_cli_functions import bedops_sort_bed_filter

with tempfile.TemporaryDirectory() as tmp:
    for mode in ("unique", "duplicates"):
        result = bedops_sort_bed_filter(["messy.bed"], str(Path(tmp) / f"{mode}.bed"), mode=mode, verbose=False)
        print(mode, [line.split("\t")[3] for line in result["preview"]])
```

- **示例输出：**

```text
unique ['y', 'z', 'x']
duplicates ['y']
```
