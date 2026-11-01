# deepTools：矩阵操作与绘图

这一组位于 `dnakit.epigenome` 的函数封装 [deepTools](https://github.com/deeptools/deepTools)（此处用 deepTools 4.0.0 验证）。其中 9 个函数是 `computeMatrixOperations` 命令的外壳，用于编辑 `computeMatrix` 写出的 gzip 矩阵：查看布局、报告取值范围、保留或重排行、保留或重命名分组与样本，以及合并多个矩阵。另外 6 个函数直接调用 deepTools 的 Python 模块（`bamPEFragmentSize`、`computeGCBias`、`plotHeatmap`、`plotEnrichment`、`utilities`、`mapReduce`）和 `deeptoolsintervals` 包：根据数值或矩阵计算片段长度密度、绘制 GC 偏好、矩阵热图和富集 PNG 图，读取比对的模板长度，以及从区域中扣除黑名单区间。它们都不会运行 `computeMatrix`，也不会在 BAM 与 BED 文件上计数 reads，更不会做比对。

deepTools 4.0.0 自带的许可文件是宽松的 MIT 类许可，DNAKit 不复制其任何代码。函数通过 `from dnakit.epigenome import ...` 导入。

依赖都是可选的，只有调用对应函数时才会导入。矩阵操作会启动外部程序 `computeMatrixOperations`，若它不在 PATH 上则抛出 `BackendUnavailableError`。绘图与辅助函数会导入 deepTools、matplotlib、numpy、pysam 和 `deeptoolsintervals`，缺少任何一个时会报错并指出对应的 extra。用 `pip install "dnakit[bio-deeptools]"` 安装。

## 示例数据

共声明两个文件。`matrix.txt` 是一个很小的、未压缩的 `computeMatrix` 输出副本（一行以 `@` 开头的 JSON 表头，加四行以制表符分隔的区域，两个分组、两个样本，每个样本三个 bin）；矩阵示例会先把它压缩成 `matrix.gz`，因为 deepTools 读取 gzip。`reads.sam` 是三条 SAM 文本比对，`template_lengths_from_bam()` 除 BAM 外也接受 SAM。其余示例直接传入输入，或自己写出小文件。矩阵示例把 `deeptools_path` 指向 Python 解释器旁边的 `computeMatrixOperations` 脚本；程序在 PATH 上时可以省略该参数。所有输出都写入临时目录。

**示例文件：** `matrix.txt`

```text
@{"upstream":[50,50],"downstream":[50,50],"body":[50,50],"bin size":[50,50],"ref point":[null,null],"verbose":false,"bin avg type":"mean","missing data as zero":false,"min threshold":null,"max threshold":null,"scale":1,"skip zeros":false,"nan after end":false,"proc number":1,"sort regions":"keep","sort using":"mean","startLabel":"TSS","endLabel":"TES","unscaled 5 prime":[0,0],"unscaled 3 prime":[0,0],"group_labels":["up.bed","down.bed"],"group_boundaries":[0,2,4],"sample_labels":["H3K4me3","H3K27ac"],"sample_boundaries":[0,3,6]}
chr1	200	300	geneA	0	+	1.000000	2.000000	3.000000	4.000000	5.000000	6.000000
chr1	500	600	geneB	0	-	6.000000	5.000000	4.000000	3.000000	2.000000	1.000000
chr1	700	800	geneC	0	+	8.000000	7.000000	6.000000	1.000000	2.000000	3.000000
chr1	900	1000	geneD	0	-	2.000000	2.000000	2.000000	9.000000	9.000000	9.000000
```

**示例文件：** `reads.sam`

```text
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:1000
r1	99	chr1	101	60	50M	=	301	250	*	*
r2	147	chr1	301	60	50M	=	101	-250	*	*
r3	0	chr1	501	60	30M5D20M	*	0	0	*	*
```

## 1) `dnakit.epigenome.bind_columns()`

- **作用：** 把多个 computeMatrix 输出并排拼接，使它们的样本合成一个矩阵。
- **处理方法：** 检查至少有两个矩阵且都存在，检查输出路径是新文件且所在文件夹存在，然后运行 `computeMatrixOperations cbind -m <矩阵> -o <输出>`。
- **API：** `bind_columns(matrices: 'Sequence[str]', output: 'str', *, deeptools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要 PATH 上的外部程序 `computeMatrixOperations`（随 deepTools 安装，extra `dnakit[bio-deeptools]`），或用 `deeptools_path` 指定其路径。
- **输入：** `matrices`：至少两个已存在的 gzip computeMatrix 文件的列表；`output`：新文件路径；`deeptools_path`：可选的程序路径。少于两个矩阵会抛出 `ConfigurationError`。
- **输出：** 一个字典，含 `status`（程序退出码为 0 且输出文件存在时为 `completed`，否则为 `failed`）、`backend`、`tool`、`command` 列表、`return_code`、`outputs`（存在的文件）和 `log_tail`（程序输出的最后 4000 个字符）。 示例说明：示例先用 `subset_matrix()` 把演示矩阵按样本拆成两个矩阵，再合并回去，并读取结果的表头。
- **示例代码：**

```python
import gzip
import json
import shutil
import sys
import tempfile
from pathlib import Path

from dnakit.epigenome import bind_columns, subset_matrix

work = Path(tempfile.mkdtemp())
matrix = work / "matrix.gz"
with open("matrix.txt", "rb") as source, gzip.open(matrix, "wb") as target:
    shutil.copyfileobj(source, target)
# computeMatrixOperations is installed next to the Python interpreter here;
# leave deeptools_path out when the program is on PATH.
tool = str(Path(sys.executable).parent / "computeMatrixOperations")


def summary(path):
    with gzip.open(path, "rt") as handle:
        lines = handle.read().splitlines()
    header = json.loads(lines[0][1:])
    names = [line.split("\t")[3] for line in lines[1:]]
    return header["group_labels"], header["sample_labels"], names

first = work / "h3k4me3.gz"
second = work / "h3k27ac.gz"
subset_matrix(str(matrix), str(first), samples=["H3K4me3"], deeptools_path=tool, verbose=False)
subset_matrix(str(matrix), str(second), samples=["H3K27ac"], deeptools_path=tool, verbose=False)
result = bind_columns(
    [str(first), str(second)], str(work / "joined.gz"), deeptools_path=tool, verbose=False
)
print(result["status"], result["tool"], result["return_code"])
print(summary(work / "joined.gz"))
```

- **示例输出：**

```text
completed computeMatrixOperations 0
(['up.bed', 'down.bed'], ['H3K4me3', 'H3K27ac'], ['geneA', 'geneB', 'geneC', 'geneD'])
```

## 2) `dnakit.epigenome.bind_rows()`

- **作用：** 把多个 computeMatrix 输出上下堆叠，使它们的区域合成一个矩阵。
- **处理方法：** 检查至少有两个存在的矩阵和一个新的输出路径，然后运行 `computeMatrixOperations rbind -m <矩阵> -o <输出>`。分组标签相同的行在结果中属于同一个分组。
- **API：** `bind_rows(matrices: 'Sequence[str]', output: 'str', *, deeptools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要 PATH 上的外部程序 `computeMatrixOperations`（随 deepTools 安装，extra `dnakit[bio-deeptools]`），或用 `deeptools_path` 指定其路径。
- **输入：** `matrices`：至少两个已存在的 gzip computeMatrix 文件的列表；`output`：新文件路径；`deeptools_path`：可选的程序路径。
- **输出：** 与 `bind_columns()` 相同的字典：`status`、`backend`、`tool`、`command`、`return_code`、`outputs` 和 `log_tail`。 示例说明：示例把演示矩阵与自身堆叠，因此每个区域出现两次，并按分组标签归类。
- **示例代码：**

```python
import gzip
import json
import shutil
import sys
import tempfile
from pathlib import Path

from dnakit.epigenome import bind_rows

work = Path(tempfile.mkdtemp())
matrix = work / "matrix.gz"
with open("matrix.txt", "rb") as source, gzip.open(matrix, "wb") as target:
    shutil.copyfileobj(source, target)
# computeMatrixOperations is installed next to the Python interpreter here;
# leave deeptools_path out when the program is on PATH.
tool = str(Path(sys.executable).parent / "computeMatrixOperations")


def summary(path):
    with gzip.open(path, "rt") as handle:
        lines = handle.read().splitlines()
    header = json.loads(lines[0][1:])
    names = [line.split("\t")[3] for line in lines[1:]]
    return header["group_labels"], header["sample_labels"], names

result = bind_rows(
    [str(matrix), str(matrix)], str(work / "stacked.gz"), deeptools_path=tool, verbose=False
)
print(result["status"], result["command"][1])
groups, samples, names = summary(work / "stacked.gz")
print(groups, samples)
print(names)
```

- **示例输出：**

```text
completed rbind
['up.bed', 'down.bed'] ['H3K4me3', 'H3K27ac']
['geneA', 'geneB', 'geneA', 'geneB', 'geneC', 'geneD', 'geneC', 'geneD']
```

## 3) `dnakit.epigenome.data_range()`

- **作用：** 报告 computeMatrix 输出中每个样本的最小值、最大值、中位数以及 10% 和 90% 分位数。
- **处理方法：** 在矩阵所在文件夹运行 `computeMatrixOperations dataRange -m <矩阵>`，返回它打印的文本。该文本是标准输出与标准错误的合并，所以表格前面可能出现上游的警告；示例只打印最后三行。退出码非零时抛出 `ConfigurationError`。
- **API：** `data_range(matrix: 'str', *, deeptools_path: 'str | None' = None, verbose: bool = True) -> 'str'`。
- **运行要求：** 需要 PATH 上的外部程序 `computeMatrixOperations`（随 deepTools 安装，extra `dnakit[bio-deeptools]`），或用 `deeptools_path` 指定其路径。
- **输入：** `matrix`：已存在的 gzip computeMatrix 文件；`deeptools_path`：可选的程序路径。
- **输出：** 一个字符串：以制表符分隔的表格，表头为 `Samples, Min, Max, Median, 10th, 90th`，每个样本一行。 示例说明：示例把制表符替换成竖线后打印表格。
- **示例代码：**

```python
import gzip
import json
import shutil
import sys
import tempfile
from pathlib import Path

from dnakit.epigenome import data_range

work = Path(tempfile.mkdtemp())
matrix = work / "matrix.gz"
with open("matrix.txt", "rb") as source, gzip.open(matrix, "wb") as target:
    shutil.copyfileobj(source, target)
# computeMatrixOperations is installed next to the Python interpreter here;
# leave deeptools_path out when the program is on PATH.
tool = str(Path(sys.executable).parent / "computeMatrixOperations")

text = data_range(str(matrix), deeptools_path=tool, verbose=False)
for line in text.splitlines()[-3:]:
    print(" | ".join(line.split("\t")))
```

- **示例输出：**

```text
Samples | Min | Max | Median | 10th | 90th
H3K4me3 | 1.0 | 8.0 | 3.5 | 2.0 | 6.9
H3K27ac | 1.0 | 9.0 | 3.5 | 1.1 | 9.0
```

## 4) `dnakit.epigenome.filter_strand()`

- **作用：** 只保留 computeMatrix 输出中位于某一条链上的区域。
- **处理方法：** 运行 `computeMatrixOperations filterStrand -m <矩阵> -o <输出> -s <链>`。`+`、`-` 和 `.` 选择链列等于该值的区域。
- **API：** `filter_strand(matrix: 'str', output: 'str', *, strand: "Literal['+', '-', '.']", deeptools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要 PATH 上的外部程序 `computeMatrixOperations`（随 deepTools 安装，extra `dnakit[bio-deeptools]`），或用 `deeptools_path` 指定其路径。
- **输入：** `matrix`：已存在的 gzip computeMatrix 文件；`output`：新文件路径；`strand`：仅限关键字参数，取 `+`、`-` 或 `.`。
- **输出：** 与 `bind_columns()` 相同的字典：`status`、`backend`、`tool`、`command`、`return_code`、`outputs` 和 `log_tail`。 示例说明：示例保留演示矩阵中位于正链的区域（geneA 和 geneC）。
- **示例代码：**

```python
import gzip
import json
import shutil
import sys
import tempfile
from pathlib import Path

from dnakit.epigenome import filter_strand

work = Path(tempfile.mkdtemp())
matrix = work / "matrix.gz"
with open("matrix.txt", "rb") as source, gzip.open(matrix, "wb") as target:
    shutil.copyfileobj(source, target)
# computeMatrixOperations is installed next to the Python interpreter here;
# leave deeptools_path out when the program is on PATH.
tool = str(Path(sys.executable).parent / "computeMatrixOperations")


def summary(path):
    with gzip.open(path, "rt") as handle:
        lines = handle.read().splitlines()
    header = json.loads(lines[0][1:])
    names = [line.split("\t")[3] for line in lines[1:]]
    return header["group_labels"], header["sample_labels"], names

result = filter_strand(
    str(matrix), str(work / "plus.gz"), strand="+", deeptools_path=tool, verbose=False
)
print(result["status"], result["command"][-2:])
print(summary(work / "plus.gz")[2])
```

- **示例输出：**

```text
completed ['-s', '+']
['geneA', 'geneC']
```

## 5) `dnakit.epigenome.filter_values()`

- **作用：** 只保留 computeMatrix 输出中数值落在某个范围内的行。
- **处理方法：** 运行 `computeMatrixOperations filterValues -m <矩阵> -o <输出> [--min X] [--max Y]`。`minimum` 与 `maximum` 至少要给一个。参数检查在矩阵与输出路径检查之后、查找程序之前进行。在本机的 numpy 2.1.3 下，deepTools 4.0.0 会在这个子命令内部以 `module 'numpy' has no attribute 'warnings'` 终止，所以真实调用返回的 `status` 为 `failed`；原因在 deepTools，不在这个外壳。
- **API：** `filter_values(matrix: 'str', output: 'str', *, minimum: 'float | None' = None, maximum: 'float | None' = None, deeptools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要 PATH 上的外部程序 `computeMatrixOperations`（随 deepTools 安装，extra `dnakit[bio-deeptools]`），或用 `deeptools_path` 指定其路径。该子命令还需要与已安装 numpy 兼容的 deepTools；deepTools 4.0.0 搭配 numpy 2 会失败。
- **输入：** `matrix`：已存在的 gzip computeMatrix 文件；`output`：新文件路径；`minimum` 与 `maximum`：可选浮点数，至少给一个，否则抛出 `ConfigurationError`。
- **输出：** 与 `bind_columns()` 相同的字典：`status`、`backend`、`tool`、`command`、`return_code`、`outputs` 和 `log_tail`。 示例说明：示例展示不给范围时的行为，它不依赖外部程序。输出之后给出一个实际调用的写法。
- **示例代码：**

```python
from pathlib import Path
import gzip
import shutil
import tempfile

from dnakit.epigenome import filter_values
from dnakit.exceptions import ConfigurationError

work = Path(tempfile.mkdtemp())
matrix = work / "matrix.gz"
with open("matrix.txt", "rb") as source, gzip.open(matrix, "wb") as target:
    shutil.copyfileobj(source, target)
try:
    filter_values(str(matrix), str(work / "kept.gz"), verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ")[-1])
```

- **示例输出：**

```text
ConfigurationError filter_values needs a minimum or a maximum.
```

依赖可用时的用法（不执行）：

```python
from dnakit.epigenome import filter_values

result = filter_values("matrix.gz", "kept.gz", minimum=3.0)
print(result["status"], result["outputs"])
```

## 6) `dnakit.epigenome.matrix_info()`

- **作用：** 返回 computeMatrix 输出的分组与样本布局。
- **处理方法：** 在矩阵所在文件夹运行 `computeMatrixOperations info -m <矩阵>`，返回它打印的文本。退出码非零时抛出 `ConfigurationError`。
- **API：** `matrix_info(matrix: 'str', *, deeptools_path: 'str | None' = None, verbose: bool = True) -> 'str'`。
- **运行要求：** 需要 PATH 上的外部程序 `computeMatrixOperations`（随 deepTools 安装，extra `dnakit[bio-deeptools]`），或用 `deeptools_path` 指定其路径。
- **输入：** `matrix`：已存在的 gzip computeMatrix 文件；`deeptools_path`：可选的程序路径。
- **输出：** 一个字符串：先是 `Groups:`，后面每行一个缩进的分组标签，再是 `Samples:`，后面每行一个缩进的样本标签。 示例说明：示例按空白字符切分该文本。
- **示例代码：**

```python
import gzip
import json
import shutil
import sys
import tempfile
from pathlib import Path

from dnakit.epigenome import matrix_info

work = Path(tempfile.mkdtemp())
matrix = work / "matrix.gz"
with open("matrix.txt", "rb") as source, gzip.open(matrix, "wb") as target:
    shutil.copyfileobj(source, target)
# computeMatrixOperations is installed next to the Python interpreter here;
# leave deeptools_path out when the program is on PATH.
tool = str(Path(sys.executable).parent / "computeMatrixOperations")

text = matrix_info(str(matrix), deeptools_path=tool, verbose=False)
print(text.split())
```

- **示例输出：**

```text
['Groups:', 'up.bed', 'down.bed', 'Samples:', 'H3K4me3', 'H3K27ac']
```

## 7) `dnakit.epigenome.relabel_matrix()`

- **作用：** 重命名 computeMatrix 输出的分组和/或样本。
- **处理方法：** 运行 `computeMatrixOperations relabel -m <矩阵> -o <输出>`，并带上 `--groupLabels` 和/或 `--sampleLabels`。两个列表至少要给一个，否则在查找程序之前就抛出 `ConfigurationError`。
- **API：** `relabel_matrix(matrix: 'str', output: 'str', *, group_labels: 'Sequence[str] | None' = None, sample_labels: 'Sequence[str] | None' = None, deeptools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要 PATH 上的外部程序 `computeMatrixOperations`（随 deepTools 安装，extra `dnakit[bio-deeptools]`），或用 `deeptools_path` 指定其路径。
- **输入：** `matrix`：已存在的 gzip computeMatrix 文件；`output`：新文件路径；`group_labels` 与 `sample_labels`：可选的新标签列表。
- **输出：** 与 `bind_columns()` 相同的字典：`status`、`backend`、`tool`、`command`、`return_code`、`outputs` 和 `log_tail`。 示例说明：示例重命名两个分组和两个样本，并读取新的表头。
- **示例代码：**

```python
import gzip
import json
import shutil
import sys
import tempfile
from pathlib import Path

from dnakit.epigenome import relabel_matrix

work = Path(tempfile.mkdtemp())
matrix = work / "matrix.gz"
with open("matrix.txt", "rb") as source, gzip.open(matrix, "wb") as target:
    shutil.copyfileobj(source, target)
# computeMatrixOperations is installed next to the Python interpreter here;
# leave deeptools_path out when the program is on PATH.
tool = str(Path(sys.executable).parent / "computeMatrixOperations")


def summary(path):
    with gzip.open(path, "rt") as handle:
        lines = handle.read().splitlines()
    header = json.loads(lines[0][1:])
    names = [line.split("\t")[3] for line in lines[1:]]
    return header["group_labels"], header["sample_labels"], names

result = relabel_matrix(
    str(matrix),
    str(work / "renamed.gz"),
    group_labels=["upstream_genes", "downstream_genes"],
    sample_labels=["K4me3", "K27ac"],
    deeptools_path=tool,
    verbose=False,
)
print(result["status"])
groups, samples, names = summary(work / "renamed.gz")
print(groups, samples, names)
```

- **示例输出：**

```text
completed
['upstream_genes', 'downstream_genes'] ['K4me3', 'K27ac'] ['geneA', 'geneB', 'geneC', 'geneD']
```

## 8) `dnakit.epigenome.sort_matrix()`

- **作用：** 让 computeMatrix 输出的行顺序与一个或多个区域文件的顺序一致。
- **处理方法：** 运行 `computeMatrixOperations sort -m <矩阵> -o <输出> -R <区域文件>`。区域文件列表为空会抛出 `ConfigurationError`，每个区域文件都必须存在。deepTools 按文件名把每个区域文件对应到矩阵的分组，所以文件名必须等于矩阵的分组标签，然后让每个分组的行按文件中的顺序排列。
- **API：** `sort_matrix(matrix: 'str', output: 'str', regions: 'Sequence[str]', *, deeptools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要 PATH 上的外部程序 `computeMatrixOperations`（随 deepTools 安装，extra `dnakit[bio-deeptools]`），或用 `deeptools_path` 指定其路径。
- **输入：** `matrix`：已存在的 gzip computeMatrix 文件；`output`：新文件路径；`regions`：非空的、已存在的 BED（或 GTF）文件列表，文件名要与分组名一致。
- **输出：** 与 `bind_columns()` 相同的字典：`status`、`backend`、`tool`、`command`、`return_code`、`outputs` 和 `log_tail`。若文件名与分组不匹配，则 `status` 为 `failed`，`log_tail` 中是上游的提示。 示例说明：示例写出 `up.bed` 和 `down.bed`，每个分组的两个区域都按相反顺序排列，并显示新的行顺序。
- **示例代码：**

```python
import gzip
import json
import shutil
import sys
import tempfile
from pathlib import Path

from dnakit.epigenome import sort_matrix

work = Path(tempfile.mkdtemp())
matrix = work / "matrix.gz"
with open("matrix.txt", "rb") as source, gzip.open(matrix, "wb") as target:
    shutil.copyfileobj(source, target)
# computeMatrixOperations is installed next to the Python interpreter here;
# leave deeptools_path out when the program is on PATH.
tool = str(Path(sys.executable).parent / "computeMatrixOperations")


def summary(path):
    with gzip.open(path, "rt") as handle:
        lines = handle.read().splitlines()
    header = json.loads(lines[0][1:])
    names = [line.split("\t")[3] for line in lines[1:]]
    return header["group_labels"], header["sample_labels"], names

# The file names must equal the group labels of the matrix.
(work / "up.bed").write_text("chr1\t500\t600\tgeneB\t0\t-\nchr1\t200\t300\tgeneA\t0\t+\n")
(work / "down.bed").write_text("chr1\t900\t1000\tgeneD\t0\t-\nchr1\t700\t800\tgeneC\t0\t+\n")
result = sort_matrix(
    str(matrix),
    str(work / "sorted.gz"),
    [str(work / "up.bed"), str(work / "down.bed")],
    deeptools_path=tool,
    verbose=False,
)
print(result["status"])
print(summary(work / "sorted.gz")[2])
```

- **示例输出：**

```text
completed
['geneB', 'geneA', 'geneD', 'geneC']
```

## 9) `dnakit.epigenome.subset_matrix()`

- **作用：** 只保留 computeMatrix 输出中选定的分组和/或样本。
- **处理方法：** 运行 `computeMatrixOperations subset -m <矩阵> -o <输出>`，并带上 `--groups` 和/或 `--samples`。两个列表至少要给一个，否则在查找程序之前就抛出 `ConfigurationError`。
- **API：** `subset_matrix(matrix: 'str', output: 'str', *, groups: 'Sequence[str] | None' = None, samples: 'Sequence[str] | None' = None, deeptools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要 PATH 上的外部程序 `computeMatrixOperations`（随 deepTools 安装，extra `dnakit[bio-deeptools]`），或用 `deeptools_path` 指定其路径。
- **输入：** `matrix`：已存在的 gzip computeMatrix 文件；`output`：新文件路径；`groups` 与 `samples`：要保留的标签的可选列表。
- **输出：** 与 `bind_columns()` 相同的字典：`status`、`backend`、`tool`、`command`、`return_code`、`outputs` 和 `log_tail`。 示例说明：示例保留分组 `up.bed` 和样本 `H3K27ac`；分组标签就是给 `computeMatrix` 的文件名。
- **示例代码：**

```python
import gzip
import json
import shutil
import sys
import tempfile
from pathlib import Path

from dnakit.epigenome import subset_matrix

work = Path(tempfile.mkdtemp())
matrix = work / "matrix.gz"
with open("matrix.txt", "rb") as source, gzip.open(matrix, "wb") as target:
    shutil.copyfileobj(source, target)
# computeMatrixOperations is installed next to the Python interpreter here;
# leave deeptools_path out when the program is on PATH.
tool = str(Path(sys.executable).parent / "computeMatrixOperations")


def summary(path):
    with gzip.open(path, "rt") as handle:
        lines = handle.read().splitlines()
    header = json.loads(lines[0][1:])
    names = [line.split("\t")[3] for line in lines[1:]]
    return header["group_labels"], header["sample_labels"], names

result = subset_matrix(
    str(matrix),
    str(work / "subset.gz"),
    groups=["up.bed"],
    samples=["H3K27ac"],
    deeptools_path=tool,
    verbose=False,
)
print(result["status"], result["command"][-4:])
print(summary(work / "subset.gz"))
```

- **示例输出：**

```text
completed ['--groups', 'up.bed', '--samples', 'H3K27ac']
(['up.bed'], ['H3K27ac'], ['geneA', 'geneB'])
```

## 10) `dnakit.epigenome.fragment_length_density()`

- **作用：** 在一个范围内把片段长度分成 100 个等宽 bin，计算其密度。
- **处理方法：** 调用 `deeptools.bamPEFragmentSize.getDensity(lengths, min_value, max_value)`，它在该范围上以 `density=True` 做 100 个 bin 的直方图。范围之外的长度不计入密度。调用期间会打开一个临时的 matplotlib 图，结束后再关闭。
- **API：** `fragment_length_density(lengths: 'list[float]', min_value: 'float', max_value: 'float', *, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 deeptools（extra `dnakit[bio-deeptools]`，会一并装上 matplotlib、numpy 和 pysam）。
- **输入：** `lengths`：1 到 1000 万个、取值 0 到 2e9 的数字列表；`min_value` 与 `max_value`：数字，且 `max_value` 大于 `min_value`。
- **输出：** 一个字典，含 `status`、`operation`、`fragments`（传入的长度个数，包含范围外的）、`range`、`bins`（100）、`bin_edges`（101 个值）、`density`（100 个值，对范围内的长度积分为 1）和 `provenance`。 示例说明：六个长度中有一个（900）在 100 到 300 之外；`fragments` 仍把它计入，而密度在其余五个上求和为 1（每个 bin 宽 2）。
- **示例代码：**

```python
from dnakit.epigenome import fragment_length_density

result = fragment_length_density([100, 150, 150, 200, 250, 900], 100, 300, verbose=False)
print(result["status"], result["fragments"], result["bins"], result["range"])
print(result["bin_edges"][:3], result["density"][:3])
print(round(sum(density * 2.0 for density in result["density"]), 6))
```

- **示例输出：**

```text
completed 6 100 [100.0, 300.0]
[100.0, 102.0, 104.0] [0.1, 0.0, 0.0]
1.0
```

## 11) `dnakit.epigenome.plot_enrichment_counts()`

- **作用：** 根据已知的特征 reads 计数绘制富集柱状图。
- **处理方法：** 检查标签、计数和总数，然后用由各选项构造的参数命名空间调用 `deeptools.plotEnrichment.plotEnrichment`。每根柱子是该特征占样本总 reads 的百分比。图以 PNG 保存到新路径（会创建其文件夹），deepTools 遗留的图会被关闭，结果中重新计算百分比。
- **API：** `plot_enrichment_counts(features: 'list[str]', sample_labels: 'list[str]', feature_counts: 'list[dict[str, int]]', total_counts: 'list[int]', output_path: 'str', *, per_sample: 'bool' = False, variable_scales: 'bool' = False, plot_title: 'str' = '', colors: 'list[str] | None' = None, alpha: 'float' = 0.9, plot_width_cm: 'float' = 20, plot_height_cm: 'float' = 20, panels_per_row: 'int' = 4, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 deeptools（extra `dnakit[bio-deeptools]`，会一并装上 matplotlib、numpy 和 pysam）。
- **输入：** `features`：特征名；`sample_labels`：每个样本一个标签；`feature_counts`：每个样本一个字典，且包含每个特征的计数；`total_counts`：每个样本一个正的总数；`output_path`：以 `.png` 结尾的新路径；另有分组选项（`per_sample`、`variable_scales`）、标题、颜色、透明度、尺寸和每行面板数。缺少特征或列表长度不一致会抛出 `ConfigurationError`。
- **输出：** 一个字典，含 `path`（PNG 的绝对路径）和 `data`：`operation`、`features`、`per_sample`、`bars_per_panel`、`samples`（每个样本的总数及各特征的计数与百分比）和 `provenance`。 示例说明：示例打印文件名、PNG 文件头检查和重新计算的百分比。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import plot_enrichment_counts

work = Path(tempfile.mkdtemp())
result = plot_enrichment_counts(
    ["exons", "introns"],
    ["sample1", "sample2"],
    [{"exons": 300, "introns": 100}, {"exons": 50, "introns": 150}],
    [1000, 500],
    str(work / "enrichment.png"),
    verbose=False,
)
print(Path(result["path"]).name, Path(result["path"]).read_bytes()[:4] == b"\x89PNG")
data = result["data"]
print(data["bars_per_panel"])
for sample in data["samples"]:
    print(sample["sample"], sample["total_counts"], sample["features"])
```

- **示例输出：**

```text
enrichment.png True
2
sample1 1000 {'exons': {'count': 300, 'percent': 30.0}, 'introns': {'count': 100, 'percent': 10.0}}
sample2 500 {'exons': {'count': 50, 'percent': 10.0}, 'introns': {'count': 150, 'percent': 30.0}}
```

## 12) `dnakit.epigenome.plot_gc_bias()`

- **作用：** 根据频率表和每个区域的 reads 计数绘制 deepTools 的 GC 偏好图。
- **处理方法：** 检查两个表，调用 `deeptools.computeGCBias.plotGCbias` 保存 PNG，再用 `computeGCBias.bin_by` 把区域按 GC 比例分成 100 个 bin，报告 0.2 到 0.7 的 bin，以及观测/期望比值的 log2。输出路径必须是新的 `.png`。
- **API：** `plot_gc_bias(frequencies: 'list[list[float]]', reads_per_gc: 'list[list[float]]', region_size: 'int', output_path: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 deeptools（extra `dnakit[bio-deeptools]`，会一并装上 matplotlib、numpy 和 pysam）。
- **输入：** `frequencies`：至少三行 `[片段数, reads 数, 观测/期望比值]`，比值大于 0；`reads_per_gc`：若干行 `[reads 数, GC 比例]`，GC 在 0 到 1 之间；`region_size`：1 到 1000 万的整数；`output_path`：以 `.png` 结尾的新路径。
- **输出：** 一个字典，含 `path`（PNG 的绝对路径）和 `data`：`operation`、`region_size`、`gc_bins`（窗口内每个 bin 的 GC 比例、区域数和 reads 中位数）、`gc_fraction`、`log2_observed_expected`、`window` 和 `provenance`。 示例说明：示例使用 10 行合成的频率表和 40 个合成区域。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import plot_gc_bias

work = Path(tempfile.mkdtemp())
frequencies = [[10 + i, 10, 1.0 + 0.1 * (i % 5)] for i in range(10)]
reads_per_gc = [[5 + i % 7, 0.2 + 0.05 * (i % 11)] for i in range(40)]
result = plot_gc_bias(frequencies, reads_per_gc, 300, str(work / "gc_bias.png"), verbose=False)
print(Path(result["path"]).name, Path(result["path"]).read_bytes()[:4] == b"\x89PNG")
data = result["data"]
print(data["region_size"], data["window"], len(data["log2_observed_expected"]))
print(data["gc_bins"][0])
print([round(value, 4) for value in data["log2_observed_expected"][:3]])
```

- **示例输出：**

```text
gc_bias.png True
300 [0.2, 0.7] 10
{'gc': 0.2, 'regions': 4, 'median_reads': 7.5}
[0.0, 0.1375, 0.263]
```

## 13) `dnakit.epigenome.plot_matrix_heatmap()`

- **作用：** 根据 computeMatrix 输出绘制带汇总曲线的热图。
- **处理方法：** 用 `deeptools.heatmapper.heatmapper().read_matrix_file()` 读取 gzip 矩阵，再调用 `deeptools.plotHeatmap.plotMatrix` 把 PNG 保存到新路径。`z_min` 或 `z_max` 为 `None` 时，由 deepTools 根据矩阵的百分位数选取。结果报告从矩阵中读到的布局。
- **API：** `plot_matrix_heatmap(matrix_path: 'str', output_path: 'str', *, plot_title: 'str' = '', colormap: 'str' = 'binary', missing_data_color: 'str' = 'black', alpha: 'float' = 1.0, z_min: 'float | None' = None, z_max: 'float | None' = None, y_min: 'float | None' = None, y_max: 'float | None' = None, average_type: 'str' = 'median', start_label: 'str' = 'TSS', end_label: 'str' = 'TES', heatmap_height_cm: 'float' = 25, heatmap_width_cm: 'float' = 7.5, per_group: 'bool' = False, what_to_show: 'str' = 'plot, heatmap and colorbar', plot_type: 'str' = 'lines', interpolation_method: 'str' = 'auto', dpi: 'int' = 200, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 deeptools（extra `dnakit[bio-deeptools]`，会一并装上 matplotlib、numpy 和 pysam）。
- **输入：** `matrix_path`：已存在的 computeMatrix 文件（不是矩阵的文件会抛出 `ConfigurationError`）；`output_path`：以 `.png` 结尾的新路径；另有标题、颜色表、缺失值颜色、透明度、z 与 y 范围（`z_max` 必须大于 `z_min`）、`average_type`（`mean`、`median`、`sum`、`min`、`max`、`std`）、起止标签、以厘米计的尺寸、`per_group`、`what_to_show`、`plot_type`（`lines`、`se`、`std`、`overlapped_lines`、`fill`）、`interpolation_method` 和 `dpi`（50 到 1200）等选项。
- **输出：** 一个字典，含 `path`（PNG 的绝对路径）和 `data`：`operation`、`matrix`、`samples`、`groups`、`regions`、`group_boundaries`、`group_sizes`、`sample_boundaries`、`columns_per_sample`、`parameters`（矩阵的表头），以及实际使用的颜色表、z 范围、平均方式、绘图类型和显示内容。 示例说明：示例以 50 dpi 绘制演示矩阵。
- **示例代码：**

```python
import gzip
import shutil
import tempfile
from pathlib import Path

from dnakit.epigenome import plot_matrix_heatmap

work = Path(tempfile.mkdtemp())
matrix = work / "matrix.gz"
with open("matrix.txt", "rb") as source, gzip.open(matrix, "wb") as target:
    shutil.copyfileobj(source, target)
result = plot_matrix_heatmap(str(matrix), str(work / "heatmap.png"), dpi=50, verbose=False)
print(Path(result["path"]).name, Path(result["path"]).read_bytes()[:4] == b"\x89PNG")
data = result["data"]
print(data["samples"], data["groups"])
print(data["regions"], data["group_sizes"], data["columns_per_sample"])
print(data["colormap"], data["average_type"], data["z_min"], data["z_max"])
```

- **示例输出：**

```text
heatmap.png True
['H3K4me3', 'H3K27ac'] ['up.bed', 'down.bed']
4 [2, 2] [3, 3]
binary median None None
```

## 14) `dnakit.epigenome.subtract_blacklist_intervals()`

- **作用：** 返回一个区域中不被任何黑名单区间覆盖的部分。
- **处理方法：** 把 BED 文件载入 `deeptoolsintervals.GTF` 树，用 `findOverlaps` 找出重叠的黑名单区间，再用 `deeptools.mapReduce.blSubtract` 计算剩余的片段。坐标为 0 起始、半开区间。无法按 BED 解析的文件会抛出 `ConfigurationError`。
- **API：** `subtract_blacklist_intervals(blacklist_path: 'str', chrom: 'str', start: 'int', end: 'int', *, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 deeptools 和 deeptoolsintervals（extra `dnakit[bio-deeptools]`）。
- **输入：** `blacklist_path`：已存在的 BED 文件；`chrom`：染色体名；`start` 与 `end`：整数，且 `end` 大于 `start`。
- **输出：** 一个字典，含 `status`、`operation`、`blacklist`、`chrom`、`start`、`end`、`blacklisted_overlaps`（重叠的区间，按文件中存储的样子给出）、`remaining`（剩余片段）、`remaining_bases` 和 `provenance`。 示例说明：区域 50 到 500 内有两个黑名单区间，所以剩下三段。
- **示例代码：**

```python
from dnakit.epigenome import subtract_blacklist_intervals

with open("blacklist.bed", "w") as handle:
    handle.write("chr1\t100\t200\nchr1\t400\t450\n")
result = subtract_blacklist_intervals("blacklist.bed", "chr1", 50, 500, verbose=False)
print(result["blacklisted_overlaps"])
print(result["remaining"], result["remaining_bases"])
```

- **示例输出：**

```text
[[100, 200], [400, 450]]
[[50, 100], [200, 400], [450, 500]] 300
```

## 15) `dnakit.epigenome.template_lengths_from_bam()`

- **作用：** 列出 BAM 或 SAM 文件中每条比对的模板长度。
- **处理方法：** 用 `pysam.AlignmentFile` 打开文件，对每条比对调用 `deeptools.utilities.getTLen`。有 TLEN 字段的 read 给出 TLEN 的绝对值（`not_abs` 为真时给出带符号的值）。TLEN 为 0 的 read 取其 CIGAR 中 M、D、= 和 X 操作的长度之和。读到 `max_reads` 条比对后停止，此时 `truncated` 为真。
- **API：** `template_lengths_from_bam(bam_path: 'str', *, not_abs: 'bool' = False, max_reads: 'int' = 1000000, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 deeptools 和 pysam（extra `dnakit[bio-deeptools]`）。
- **输入：** `bam_path`：已存在的 BAM 或 SAM 文件（无法读取的文件会抛出 `ConfigurationError`）；`not_abs`：布尔值；`max_reads`：1 到 1 亿的整数。
- **输出：** 一个字典，含 `status`、`operation`、`bam`、`not_abs`、`alignments`、`truncated`、`template_lengths`（每条比对一个整数）、`summary`（`min`、`median` 和 `max`，文件为空时为 `None`）和 `provenance`。 示例说明：第三条比对的 TLEN 为 0，CIGAR 为 `30M5D20M`，所以长度是 30 + 5 + 20 = 55。
- **示例代码：**

```python
from dnakit.epigenome import template_lengths_from_bam

result = template_lengths_from_bam("reads.sam", verbose=False)
print(result["alignments"], result["truncated"], result["template_lengths"], result["summary"])
signed = template_lengths_from_bam("reads.sam", not_abs=True, verbose=False)
print(signed["template_lengths"])
short = template_lengths_from_bam("reads.sam", max_reads=2, verbose=False)
print(short["alignments"], short["truncated"])
```

- **示例输出：**

```text
3 False [250, 250, 55] {'min': 55, 'median': 250, 'max': 250}
[250, -250, 55]
2 True
```
