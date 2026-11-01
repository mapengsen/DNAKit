# cooler 与 cooltools：Hi-C 接触矩阵

这一组函数是 `dnakit.epigenome` 中基于可选依赖的适配器，封装 [cooler](https://github.com/open2c/cooler)（BSD-3-Clause，Hi-C 接触矩阵的 HDF5 格式与命令行工具）和 [cooltools](https://github.com/open2c/cooltools)（MIT，建立在 cooler 之上的分析工具包）。它们可以由读段对、像素表或 hiclib 文件构建 cooler，对接触列表分箱、排序和建索引，对矩阵做平衡、粗化、合并、多分辨率缩放和导出，查看并整理 cooler 文件，把接触矩阵区域画成图像文件，列出和下载 cooltools 示例数据，以及绘制 saddle 图。它们不做 loop、结构域或区室的调用。

大多数函数只运行一条 `cooler` 子命令，参数是带类型的关键字列表，对应固定的命令行选项；不使用 shell，也不能传入其他选项。`cooler` 程序的查找顺序是显式的 `cooler_path`、`PATH`、当前 Python 解释器所在目录。这些函数返回一个字典，包含 `status`（`"completed"` 或 `"failed"`）、`backend`、`backend_version`、`subcommand`、`command`、`return_code`、`outputs`（已存在的输出文件）和 `log_tail`（cooler 日志的末尾），另有各节列出的额外字段。参数无效时抛出 `ConfigurationError`；缺少程序或软件包时抛出 `BackendUnavailableError`。少数函数在进程内运行（`cooler_csort_flip_pairs`、`cooler_pairs_header`、`cooler_show_view`、`cooler_resolve_bins` 以及 cooltools 的几个函数）。除非某节另有说明，输出文件必须是新路径。

这些函数通过 `from dnakit.epigenome import ...` 导入。依赖都是可选的，只有调用对应函数时才会导入；缺少包时会报错并给出安装命令。cooler 的函数需要 `dnakit[bio-cooler]`，cooltools 的函数需要 `dnakit[bio-cooltools]`，画图的函数还需要 `dnakit[bio-plots]`。

## 示例数据

所有示例都使用总长 160 bp 的两条染色体，按 20 bp 分箱（共 8 个箱）。每个示例用前两个文件在临时文件夹里自己构建 cooler，因此不会在示例文件旁边写入任何内容。

**示例文件：** `chrom.sizes`

```text
chr1	100
chr2	60
```

**示例文件：** `pairs.txt`

```text
## pairs format v1.0
#columns: readID chrom1 pos1 chrom2 pos2
r1	chr1	5	chr1	45
r2	chr1	10	chr1	30
r3	chr1	25	chr1	85
r4	chr1	50	chr2	10
r5	chr2	5	chr2	50
r6	chr1	95	chr1	3
r7	chr1	12	chr2	45
r8	chr1	5	chr3	7
```

最后一条读段对位于 `chr3`，不在 `chrom.sizes` 中，所以分箱时被丢弃。第六条读段对的两端顺序反了，会被翻转。

**示例文件：** `pixels.coo`

```text
0	1	3
0	2	1
3	3	4
```

**示例文件：** `genome.fa`

```text
>chr1
AAAAAAAAAAAAAAAAAAAAAAGCTTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCAAGCTTGGGGGGGGGGGGGGGGGGGGGGGGGGGG
>chr2
TTTTTTTTTTTTTTTTTTTTTTTTTTTTAAGCTTCCCCCCCCCCCCCCCCCCCCCCCCCC
```

**示例文件：** `datasets.tsv`

```text
# key	comment	link	filename	checksum
tiny	A tiny text file	https://example.org/tiny.txt	tiny.txt	b1946ac92492d2347c6235b4d2611184
small	Another file	https://example.org/small.txt	small.txt	00000000000000000000000000000000
```

## 1) `dnakit.epigenome.cooler_attrs()`

- **作用：** 以 YAML 文本和字典两种形式返回 cooler 文件或组的 HDF5 属性。
- **处理方法：** 运行 `cooler attrs`，可带 `--level`。cooler 打印的文本（`cooler.fileops.pprint_attr_tree`）放在 `attributes_yaml` 中；若安装了 PyYAML，则用 `yaml.safe_load` 解析后放在 `attributes` 中。
- **API：** `cooler_attrs(uri: 'str', *, level: 'int | None' = None, cooler_path: 'str | None' = None, timeout_seconds: 'float' = 600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 cooler（额外依赖 `dnakit[bio-cooler]`）；PyYAML 是 cooler 的依赖。
- **输入：** `uri` 是已存在的 cooler 文件或 `file::/group`。`level` 是最大深度（不小于 0 的整数）。
- **输出：** 通用字段，加上 `uri`、`attributes_yaml` 和 `attributes`（运行失败时为 `None`）。根组的属性在键 `"@attrs"` 下。
- **示例代码：**

```python
import tempfile

from dnakit.epigenome import cooler_attrs, cooler_cload_pairs

with tempfile.TemporaryDirectory() as folder:
    cool = f"{folder}/a.cool"
    cooler_cload_pairs("chrom.sizes:20", "pairs.txt", cool, chrom1=2, pos1=3, chrom2=4, pos2=5, verbose=False)
    result = cooler_attrs(cool, verbose=False)
    root = result["attributes"]["@attrs"]
    print(result["status"], root["bin-size"], root["nbins"], root["sum"])
    print(sorted(key for key in result["attributes"] if key != "@attrs"))
```

- **示例输出：**

```text
completed 20 8 7
['bins', 'chroms', 'indexes', 'pixels']
```

## 2) `dnakit.epigenome.cooler_balance()`

- **作用：** 用迭代校正法平衡 cooler，并把权重存入其箱表。
- **处理方法：** 对文件运行 `cooler balance`，文件会被就地修改（在 `bins` 中添加 `weight` 列，或 `name` 指定的列）。运行完成后，DNAKit 调用 `cooler_balance_check()`，把结果作为 `weights_stored` 返回。保持为 `None` 的选项沿用 cooler 的默认值（ignore_diags 2、mad_max 5、min_nnz 10、min_count 0、tol 1e-5、max_iters 200、convergence_policy store_final、8 个进程）。不提供 `cooler balance --stdout`：请用 `cooler_dump(table="bins", columns=["weight"])` 读取权重。
- **API：** `cooler_balance(cool_uri: 'str', *, cis_only: 'bool' = False, trans_only: 'bool' = False, ignore_diags: 'int | None' = None, ignore_dist: 'int | None' = None, mad_max: 'int | None' = None, min_nnz: 'int | None' = None, min_count: 'int | None' = None, blacklist: 'str | None' = None, nproc: 'int | None' = None, chunksize: 'int | None' = None, tol: 'float | None' = None, max_iters: 'int | None' = None, name: 'str | None' = None, force: 'bool' = False, convergence_policy: "Literal['store_final', 'store_nan', 'discard', 'error'] | None" = None, cooler_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 cooler（额外依赖 `dnakit[bio-cooler]`）。
- **输入：** `cool_uri` 是已存在的 cooler。`cis_only` 与 `trans_only` 不能同时为真。整数选项不小于 0（`nproc`、`chunksize`、`max_iters` 不小于 1）。`blacklist` 是至少两行的已存在 BED 文件（cooler 0.10.4 会把只有一行的文件当作表头而失败）。`force` 会重新计算已有的权重。
- **输出：** 通用字段，加上 `cool_uri`、`weights_column` 和 `weights_stored`。在默认的 `store_final` 策略下，未收敛的运行同样会保存结果，原因记录在 `log_tail` 中。
- **示例代码：**

```python
import tempfile

from dnakit.epigenome import cooler_balance, cooler_cload_pairs

with tempfile.TemporaryDirectory() as folder:
    cool = f"{folder}/a.cool"
    cooler_cload_pairs("chrom.sizes:20", "pairs.txt", cool, chrom1=2, pos1=3, chrom2=4, pos2=5, verbose=False)
    result = cooler_balance(cool, ignore_diags=0, min_nnz=0, mad_max=0, cis_only=True, nproc=1, verbose=False)
    print(result["status"], result["weights_column"], result["weights_stored"])
```

- **示例输出：**

```text
completed weight True
```

## 3) `dnakit.epigenome.cooler_balance_check()`

- **作用：** 报告 cooler 是否已经有平衡权重列。
- **处理方法：** 运行 `cooler balance --check`。cooler 对没有权重的文件和无法读取的文件都以退出码 1 结束，所以 DNAKit 读取它的消息：含 “is balanced” 则为 `True`，含 “column found” 则为 `False`，其他情况按失败处理，`balanced` 为 `None`。
- **API：** `cooler_balance_check(cool_uri: 'str', *, name: 'str | None' = None, cooler_path: 'str | None' = None, timeout_seconds: 'float' = 600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 cooler（额外依赖 `dnakit[bio-cooler]`）。
- **输入：** `cool_uri` 是已存在的 cooler；`name` 是要查找的权重列（默认 `weight`）。
- **输出：** 通用字段，加上 `cool_uri`、`weights_column` 和 `balanced`。
- **示例代码：**

```python
import tempfile

from dnakit.epigenome import cooler_balance, cooler_balance_check, cooler_cload_pairs

with tempfile.TemporaryDirectory() as folder:
    cool = f"{folder}/a.cool"
    cooler_cload_pairs("chrom.sizes:20", "pairs.txt", cool, chrom1=2, pos1=3, chrom2=4, pos2=5, verbose=False)
    print(cooler_balance_check(cool, verbose=False)["balanced"])
    cooler_balance(cool, verbose=False)
    print(cooler_balance_check(cool, verbose=False)["balanced"])
```

- **示例输出：**

```text
False
True
```

## 4) `dnakit.epigenome.cooler_cload_hiclib()`

- **作用：** 把 hiclib 的 HDF5 接触列表分箱为 cooler。
- **处理方法：** 运行 `cooler cload hiclib`。HDF5 文件包含数据集 `chrms1`、`cuts1`、`chrms2` 和 `cuts2`（按箱顺序的 0 起始染色体编号和 0 起始位置）。接触对必须按第一条染色体排序，并且对 `bins` 的染色体顺序是上三角。
- **API：** `cooler_cload_hiclib(bins: 'str', pairs_path: 'str', output_path: 'str', *, chunksize: 'int | None' = None, metadata_path: 'str | None' = None, assembly: 'str | None' = None, cooler_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 cooler（额外依赖 `dnakit[bio-cooler]`）；h5py 是 cooler 的依赖。
- **输入：** `bins` 是 `<chromsizes 文件>:<箱大小 bp>` 或定义箱的 BED 文件路径。`pairs_path` 是已存在的 HDF5 文件；`output_path` 必须不存在。`chunksize`（cooler 默认 100,000,000）是正整数。`metadata_path` 是已存在的 JSON 文件，`assembly` 是 `hg38` 之类的基因组名。
- **输出：** 通用字段（`subcommand` 为 `"cload hiclib"`），加上 `output_uri`。
- **示例代码：**

```python
import tempfile

import h5py
import numpy as np

from dnakit.epigenome import cooler_cload_hiclib, cooler_info

with tempfile.TemporaryDirectory() as folder:
    with h5py.File(f"{folder}/contacts.h5", "w") as handle:
        handle["chrms1"] = np.array([0, 0, 0, 1], dtype="int32")
        handle["cuts1"] = np.array([5, 10, 50, 5], dtype="int64")
        handle["chrms2"] = np.array([0, 0, 1, 1], dtype="int32")
        handle["cuts2"] = np.array([45, 30, 10, 50], dtype="int64")
    result = cooler_cload_hiclib("chrom.sizes:20", f"{folder}/contacts.h5", f"{folder}/h.cool", verbose=False)
    info = cooler_info(f"{folder}/h.cool", verbose=False)["info"]
    print(result["status"], result["subcommand"], info["nbins"], info["sum"])
```

- **示例输出：**

```text
completed cload hiclib 8 4
```

## 5) `dnakit.epigenome.cooler_cload_pairix()`

- **作用：** 把 bgzip 压缩并用 pairix 建了索引的接触对文件分箱为 cooler。
- **处理方法：** 先检查文件旁是否有 `.px2` 索引，再运行 `cooler cload pairix`。cooler 的这个模式需要 Python 包 `pypairix`；没有它时，运行会被报告为失败，导入错误记录在 `log_tail` 中。
- **API：** `cooler_cload_pairix(bins: 'str', pairs_path: 'str', output_path: 'str', *, zero_based: 'bool' = False, nproc: 'int | None' = None, max_split: 'int | None' = None, block_char: 'str | None' = None, metadata_path: 'str | None' = None, assembly: 'str | None' = None, cooler_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 cooler（额外依赖 `dnakit[bio-cooler]`）和 Python 包 `pypairix`（`python -m pip install pypairix`）；建立索引需要 pairix 程序（`conda install -c bioconda pairix`）。示例使用占位文件：它通过索引检查后，展示未安装 `pypairix` 时 cooler 的报告；安装该包后，请使用输出下方给出的调用。
- **输入：** `bins` 同 `cooler_cload_hiclib()`。`pairs_path` 是已存在的 bgzip 文件，旁边有 `.px2` 索引（否则抛出 `ConfigurationError`）。`nproc`（cooler 默认 8）和 `max_split`（默认 2）是正整数；`block_char`（默认 `|`）是不含引号和反斜杠的短字符串。
- **输出：** 通用字段（`subcommand` 为 `"cload pairix"`），加上 `output_uri`。
- **示例代码：**

```python
import tempfile

from dnakit.epigenome import cooler_cload_pairix

with tempfile.TemporaryDirectory() as folder:
    pairs = f"{folder}/pairs.gz"
    for path in (pairs, pairs + ".px2"):
        with open(path, "wb") as handle:
            handle.write(b"x")
    result = cooler_cload_pairix("chrom.sizes:20", pairs, f"{folder}/p.cool", verbose=False)
    print(result["status"], result["subcommand"], "pypairix is required" in result["log_tail"])
```

- **示例输出：**

```text
failed cload pairix True
```

安装了依赖后的用法：

```python
from dnakit.epigenome import cooler_csort, cooler_cload_pairix

cooler_csort("pairs.txt", "chrom.sizes", "pairs.sorted.gz", chrom1=2, pos1=3, chrom2=4, pos2=5, index="pairix")
cooler_cload_pairix("chrom.sizes:20", "pairs.sorted.gz", "pairix.cool")
```

## 6) `dnakit.epigenome.cooler_cload_pairs()`

- **作用：** 把未排序的读段对文本文件分箱为 cooler。
- **处理方法：** 按四个字段编号运行 `cooler cload pairs`。cooler 自己完成排序和合并，把每个接触转成上三角，并丢弃不在箱中的染色体上的接触对。
- **API：** `cooler_cload_pairs(bins: 'str', pairs_path: 'str', output_path: 'str', *, chrom1: 'int', pos1: 'int', chrom2: 'int', pos2: 'int', zero_based: 'bool' = False, comment_char: 'str | None' = None, no_symmetric_upper: 'bool' = False, input_copy_status: "Literal['unique', 'duplex'] | None" = None, fields: 'Sequence[str] | None' = None, chunksize: 'int | None' = None, mergebuf: 'int | None' = None, max_merge: 'int | None' = None, temp_dir: 'str | None' = None, no_delete_temp: 'bool' = False, storage_options: 'str | None' = None, metadata_path: 'str | None' = None, assembly: 'str | None' = None, append: 'bool' = False, cooler_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 cooler（额外依赖 `dnakit[bio-cooler]`）。
- **输入：** `bins` 是 `<chromsizes 文件>:<箱大小 bp>` 或 BED 文件。`pairs_path` 是已存在的文本文件，可以是 gzip 压缩的；不支持标准输入。`chrom1`、`pos1`、`chrom2`、`pos2` 是四个互不相同的 1 起始列号。位置默认 1 起始，`zero_based` 为真时为 0 起始。注释行以 `comment_char`（默认 `#`）开头。`fields` 的条目形如 `score=7:dtype=float,agg=mean`。`storage_options` 形如 `compression=gzip,compression_opts=6`。`append=True` 时把 cooler 添加到已有文件，`output_path` 可以是 `file.mcool::resolutions/10000`；否则 `output_path` 必须不存在。
- **输出：** 通用字段（`subcommand` 为 `"cload pairs"`），加上 `output_uri`。
- **示例代码：**

```python
import tempfile

from dnakit.epigenome import cooler_cload_pairs, cooler_info

with tempfile.TemporaryDirectory() as folder:
    cool = f"{folder}/a.cool"
    result = cooler_cload_pairs("chrom.sizes:20", "pairs.txt", cool, chrom1=2, pos1=3, chrom2=4, pos2=5, verbose=False)
    info = cooler_info(cool, verbose=False)["info"]
    print(result["status"], result["subcommand"], result["return_code"])
    print(info["bin-size"], info["nbins"], info["nnz"], info["sum"])
```

- **示例输出：**

```text
completed cload pairs 0
20 8 7 7
```

## 7) `dnakit.epigenome.cooler_cload_tabix()`

- **作用：** 把 bgzip 压缩并用 tabix 建了索引的接触对文件分箱为 cooler。
- **处理方法：** 先检查文件旁是否有 `.tbi` 或 `.csi` 索引，再运行 `cooler cload tabix`。索引必须建在 chrom1 和 pos1 上，`bins` 的染色体顺序必须与文件中的顺序一致。
- **API：** `cooler_cload_tabix(bins: 'str', pairs_path: 'str', output_path: 'str', *, chrom2: 'int | None' = None, pos2: 'int | None' = None, zero_based: 'bool' = False, nproc: 'int | None' = None, max_split: 'int | None' = None, metadata_path: 'str | None' = None, assembly: 'str | None' = None, cooler_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 cooler（额外依赖 `dnakit[bio-cooler]`）以及 cooler 用来读取文件的 pysam。示例用 pysam 建立索引；`tabix` 和 `bgzip` 程序（htslib）只有 `cooler_csort()` 才需要。
- **输入：** `bins` 同上。`pairs_path` 是已存在的 bgzip 文件，旁边有索引（否则抛出 `ConfigurationError`）。`chrom2`（cooler 默认 4）和 `pos2`（默认 5）是 1 起始列号；`nproc`（默认 8）和 `max_split`（默认 2）是正整数。
- **输出：** 通用字段（`subcommand` 为 `"cload tabix"`），加上 `output_uri`。
- **示例代码：**

```python
import tempfile

import pysam

from dnakit.epigenome import cooler_cload_tabix, cooler_csort_flip_pairs, cooler_csort_sort_pairs, cooler_info

with tempfile.TemporaryDirectory() as folder:
    fields = dict(chrom1=2, pos1=3, chrom2=4, pos2=5)
    cooler_csort_flip_pairs("pairs.txt", "chrom.sizes", f"{folder}/flipped.tsv", verbose=False, **fields)
    cooler_csort_sort_pairs(f"{folder}/flipped.tsv", f"{folder}/sorted.tsv", index="tabix", verbose=False, **fields)
    pysam.tabix_compress(f"{folder}/sorted.tsv", f"{folder}/sorted.tsv.gz")
    pysam.tabix_index(f"{folder}/sorted.tsv.gz", seq_col=1, start_col=2, end_col=2, zerobased=False, force=True)
    result = cooler_cload_tabix("chrom.sizes:20", f"{folder}/sorted.tsv.gz", f"{folder}/t.cool", chrom2=4, pos2=5, verbose=False)
    print(result["status"], result["subcommand"], cooler_info(f"{folder}/t.cool", verbose=False)["info"]["sum"])
```

- **示例输出：**

```text
completed cload tabix 7
```

## 8) `dnakit.epigenome.cooler_coarsen()`

- **作用：** 把 k×k 像素块合并，将 cooler 粗化为更低的分辨率。
- **处理方法：** 运行 `cooler coarsen --out`。`fields` 中列出的数值列会被聚合（默认是 `count` 列求和）。
- **API：** `cooler_coarsen(cool_uri: 'str', output_path: 'str', *, factor: 'int | None' = None, nproc: 'int | None' = None, chunksize: 'int | None' = None, fields: 'Sequence[str] | None' = None, append: 'bool' = False, cooler_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 cooler（额外依赖 `dnakit[bio-cooler]`）。
- **输入：** `cool_uri` 是已存在的 cooler。`factor`（cooler 默认 2）是不小于 2 的整数；`nproc`（默认 1）和 `chunksize`（默认 10,000,000）是正整数。`fields` 的条目形如 `score:dtype=float,agg=mean`。`append=True` 时结果添加到已有文件，`output_path` 可带 `::group`；否则它必须不存在。
- **输出：** 通用字段，加上 `input_uri` 和 `output_uri`。
- **示例代码：**

```python
import tempfile

from dnakit.epigenome import cooler_cload_pairs, cooler_coarsen, cooler_info

with tempfile.TemporaryDirectory() as folder:
    cool = f"{folder}/a.cool"
    cooler_cload_pairs("chrom.sizes:20", "pairs.txt", cool, chrom1=2, pos1=3, chrom2=4, pos2=5, verbose=False)
    result = cooler_coarsen(cool, f"{folder}/c.cool", factor=2, verbose=False)
    info = cooler_info(f"{folder}/c.cool", verbose=False)["info"]
    print(result["status"], info["bin-size"], info["nbins"], info["nnz"], info["sum"])
```

- **示例输出：**

```text
completed 40 5 6 7
```

## 9) `dnakit.epigenome.cooler_cp()`

- **作用：** 把 cooler 复制到另一个文件或文件内的另一个组。
- **处理方法：** 运行 `cooler cp`。`overwrite=True` 时 DNAKit 传入 `-w`，它会截断并替换已存在的整个目标文件；当源和目标在同一个文件内时，cooler 0.10.4 无法执行，会以 “unable to truncate a file which is already open” 失败。
- **API：** `cooler_cp(source_uri: 'str', destination_uri: 'str', *, overwrite: 'bool' = False, cooler_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 cooler（额外依赖 `dnakit[bio-cooler]`）。
- **输入：** `source_uri` 是已存在的 cooler（`file` 或 `file::/group`）；`destination_uri` 是位于已存在文件夹中的 `file` 或 `file::/group`。
- **输出：** 通用字段，加上 `source_uri`、`destination_uri` 和 `source_file`；目标文件存在时，`outputs` 列出它。
- **示例代码：**

```python
import tempfile

from dnakit.epigenome import cooler_cload_pairs, cooler_cp, cooler_info

with tempfile.TemporaryDirectory() as folder:
    cool = f"{folder}/a.cool"
    cooler_cload_pairs("chrom.sizes:20", "pairs.txt", cool, chrom1=2, pos1=3, chrom2=4, pos2=5, verbose=False)
    result = cooler_cp(cool, f"{folder}/b.cool", verbose=False)
    print(result["status"], result["return_code"], cooler_info(f"{folder}/b.cool", field="sum", verbose=False)["value"])
```

- **示例输出：**

```text
completed 0 7
```

## 10) `dnakit.epigenome.cooler_csort()`

- **作用：** 对接触列表定向、排序、bgzip 压缩并建索引，使其可用于 `cooler_cload_tabix()` 或 `cooler_cload_pairix()`。
- **处理方法：** 运行 `cooler csort`。两端会被排成对 `chromosomes_path` 的染色体顺序是上三角（其他染色体上的接触对被丢弃），各行排序（tabix 用位置排序，pairix 用块排序），用 bgzip 压缩并建索引。运行前，DNAKit 检查 `sort`、`bgzip` 以及 `tabix` 或 `pairix` 是否在 `PATH`（或 cooler 旁边），否则抛出 `BackendUnavailableError`。结果列出排序后的文件及其索引（`.tbi` 或 `.px2`）。
- **API：** `cooler_csort(pairs_path: 'str', chromosomes_path: 'str', output_path: 'str', *, chrom1: 'int', pos1: 'int', chrom2: 'int', pos2: 'int', index: "Literal['tabix', 'pairix']" = 'pairix', nproc: 'int | None' = None, zero_based: 'bool' = False, sep: 'str | None' = None, comment_char: 'str | None' = None, cooler_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 cooler（额外依赖 `dnakit[bio-cooler]`），以及 `PATH` 上的外部程序 `sort`、`bgzip` 和 `tabix`（htslib，`conda install -c bioconda htslib`）或 `pairix`。不支持读取标准输入和自由的 `--sort-options` 文本；排序调优见 `cooler_csort_sort_pairs()`。生成本页的环境没有安装 `bgzip` 和 `tabix`，所以示例展示这种情况下抛出的错误。
- **输入：** `pairs_path` 是已存在的文本文件；`chromosomes_path` 按期望的顺序列出染色体（chromsizes 文件可用）；`output_path` 必须不存在，且不同于 `pairs_path`。`chrom1`、`pos1`、`chrom2`、`pos2` 是四个互不相同的 1 起始列号。`index` 为 `"tabix"` 或 `"pairix"`（默认）。`sep` 默认为制表符，`comment_char` 默认为 `#`。
- **输出：** 通用字段（`outputs` 包含排序后的文件及其索引），加上 `pairs_path`、`index` 和 `index_file`。
- **示例代码：**

```python
import tempfile

from dnakit.epigenome import cooler_csort
from dnakit.exceptions import BackendUnavailableError

with tempfile.TemporaryDirectory() as folder:
    try:
        cooler_csort("pairs.txt", "chrom.sizes", f"{folder}/pairs.sorted.gz", chrom1=2, pos1=3, chrom2=4, pos2=5, index="tabix", verbose=False)
    except BackendUnavailableError as error:
        print(type(error).__name__, "cooler csort needs bgzip" in str(error))
```

- **示例输出：**

```text
BackendUnavailableError True
```

安装了依赖后的用法：

```python
from dnakit.epigenome import cooler_csort

result = cooler_csort("pairs.txt", "chrom.sizes", "pairs.sorted.gz", chrom1=2, pos1=3, chrom2=4, pos2=5, index="tabix")
print(result["status"], result["index_file"])
```

## 11) `dnakit.epigenome.cooler_csort_flip_pairs()`

- **作用：** 按给定的染色体顺序，把每个接触转成上三角（csort 的翻转步骤）。
- **处理方法：** 在进程内运行，算法与 `cooler csort` 作为 `python -c` 管道阶段运行的相同。未列出的染色体上的接触被丢弃，以 `comment_char` 开头的行被跳过；当 chrom1 排在 chrom2 之后，或在同一染色体上 pos1 大于 pos2 时，两端互换。`.gz` 输入通过 gzip 读取。输出列用制表符连接。不需要 cooler 或其他程序。
- **API：** `cooler_csort_flip_pairs(pairs_path: 'str', chromosomes_path: 'str', output_path: 'str', *, chrom1: 'int', pos1: 'int', chrom2: 'int', pos2: 'int', sep: 'str' = '\t', comment_char: 'str' = '#', verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 安装后即可使用（普通 pip install）。
- **输入：** `pairs_path` 和 `chromosomes_path` 是已存在的文件（`chromosomes_path` 的第一列（制表符分隔）给出染色体顺序）。`output_path` 必须不存在。四个字段编号是互不相同的 1 起始整数；列数不足或位置不是整数的行会抛出 `ConfigurationError`，并删除未写完的输出。
- **输出：** `status`、`backend`（`"python"`）、`outputs`、`chromosomes`，以及计数 `lines_read`、`comment_lines`、`dropped`、`flipped` 和 `written`。
- **示例代码：**

```python
import tempfile

from dnakit.epigenome import cooler_csort_flip_pairs

with tempfile.TemporaryDirectory() as folder:
    out = f"{folder}/flipped.tsv"
    result = cooler_csort_flip_pairs("pairs.txt", "chrom.sizes", out, chrom1=2, pos1=3, chrom2=4, pos2=5, verbose=False)
    print({key: value for key, value in result.items() if key not in ("outputs", "backend")})
    lines = open(out).read().splitlines()
    print(" ".join(lines[5].split()))
```

- **示例输出：**

```text
{'status': 'completed', 'chromosomes': 2, 'lines_read': 10, 'comment_lines': 2, 'dropped': 1, 'flipped': 1, 'written': 7}
r6 chr1 3 chr1 95
```

## 12) `dnakit.epigenome.cooler_csort_index_pairs()`

- **作用：** 用 tabix 或 pairix 为已排序、bgzip 压缩的接触列表建索引（csort 的索引步骤）。
- **处理方法：** 检查文件是否以 gzip 魔数开头，查找 `tabix` 或 `pairix` 程序，并用与 `cooler csort` 相同的参数运行它。tabix 对 chrom1 和 pos1 建索引并写出 `<pairs_path>.tbi`；pairix 对 chrom1、chrom2、pos1 和 pos2 建索引并写出 `<pairs_path>.px2`。
- **API：** `cooler_csort_index_pairs(pairs_path: 'str', *, chrom1: 'int', pos1: 'int', chrom2: 'int', pos2: 'int', index: "Literal['tabix', 'pairix']" = 'pairix', zero_based: 'bool' = False, index_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要 `PATH` 上、Python 解释器旁或通过 `index_path` 给出的外部程序 `tabix`（htslib，`conda install -c bioconda htslib`）或 `pairix`（`conda install -c bioconda pairix`）。示例给出一个不存在的 `index_path`，使结果不依赖机器。
- **输入：** `pairs_path` 是已存在的 bgzip 文件，其排序方式与 `cooler_csort_sort_pairs()` 对同一 `index` 写出的相同（否则抛出 `ConfigurationError`）。四个字段编号是互不相同的 1 起始整数。`index_path` 是要使用的可执行文件。
- **输出：** `status`、`backend`、`command`、`return_code`、`outputs`、`index_file` 和 `log_tail`。
- **示例代码：**

```python
import gzip
import tempfile

from dnakit.epigenome import cooler_csort_index_pairs
from dnakit.exceptions import BackendUnavailableError

with tempfile.TemporaryDirectory() as folder:
    pairs = f"{folder}/pairs.sorted.gz"
    with gzip.open(pairs, "wt") as handle:
        handle.write("r1\tchr1\t5\tchr1\t45\n")
    try:
        cooler_csort_index_pairs(pairs, chrom1=2, pos1=3, chrom2=4, pos2=5, index="tabix", index_path="/nonexistent/tabix", verbose=False)
    except BackendUnavailableError as error:
        print(type(error).__name__, "tabix_path does not exist" in str(error))
```

- **示例输出：**

```text
BackendUnavailableError True
```

安装了依赖后的用法：

```python
from dnakit.epigenome import cooler_csort_index_pairs

result = cooler_csort_index_pairs("pairs.sorted.gz", chrom1=2, pos1=3, chrom2=4, pos2=5, index="tabix")
print(result["status"], result["index_file"])
```

## 13) `dnakit.epigenome.cooler_csort_read_pairs()`

- **作用：** 把接触列表写成纯文本，必要时解压 `.gz` 输入（csort 的读取步骤）。
- **处理方法：** 与 `cooler csort` 一致：以 `.gz` 结尾的文件在有 pigz 时用 `pigz -dc` 读取，否则用 `gzip -dc`；程序的输出直接写入 `output_path`。其他文件原样复制（cooler 让它经过 `cat`），此时 `command` 为 `None`。
- **API：** `cooler_csort_read_pairs(pairs_path: 'str', output_path: 'str', *, pigz_path: 'str | None' = None, gzip_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 普通文件不需要任何东西。`.gz` 文件需要 `PATH` 上的外部程序 `pigz`（`conda install -c conda-forge pigz`）或 `gzip`；两者都没有时抛出 `BackendUnavailableError`。
- **输入：** `pairs_path` 是已存在的文件；`output_path` 必须不存在。`pigz_path` 和 `gzip_path` 是可选的显式程序路径，给出时必须存在。
- **输出：** `status`、`backend`（`"copy"`、`"pigz"` 或 `"gzip"`）、`command`、`return_code`、`outputs`、`bytes_written` 和 `log_tail`。
- **示例代码：**

```python
import gzip
import os
import shutil
import tempfile

from dnakit.epigenome import cooler_csort_read_pairs

with tempfile.TemporaryDirectory() as folder:
    with open("pairs.txt", "rb") as source, gzip.open(f"{folder}/pairs.txt.gz", "wb") as target:
        shutil.copyfileobj(source, target)
    result = cooler_csort_read_pairs(f"{folder}/pairs.txt.gz", f"{folder}/plain.txt", verbose=False)
    print(result["status"], result["return_code"], result["bytes_written"] == os.path.getsize("pairs.txt"))
    copied = cooler_csort_read_pairs("pairs.txt", f"{folder}/copy.txt", verbose=False)
    print(copied["backend"], copied["command"], copied["bytes_written"])
```

- **示例输出：**

```text
completed 0 True
copy None 209
```

## 14) `dnakit.epigenome.cooler_csort_sort_pairs()`

- **作用：** 按 cooler csort 的排序键对接触列表排序（csort 的排序步骤）。
- **处理方法：** 运行系统的 `sort`，设置 `LC_ALL=C`（只对子进程生效），并用 `sort -o` 写出结果。`index="tabix"` 按 chrom1、pos1、chrom2、pos2 排序（位置排序）；`index="pairix"` 按 chrom1、chrom2、pos1、pos2 排序（块排序）。位置按数值排序，名称按字节排序。`nproc` 添加 `--parallel`，`buffer_size` 添加 `--buffer-size`，`temp_dir` 添加 `-T`。
- **API：** `cooler_csort_sort_pairs(pairs_path: 'str', output_path: 'str', *, chrom1: 'int', pos1: 'int', chrom2: 'int', pos2: 'int', index: "Literal['tabix', 'pairix']" = 'pairix', nproc: 'int | None' = None, buffer_size: 'str | None' = None, temp_dir: 'str | None' = None, sort_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要外部程序 `sort`（GNU coreutils，Linux 和 macOS 自带）。
- **输入：** `pairs_path` 是已存在的、以制表符分隔的文件，应当已经翻转过（见 `cooler_csort_flip_pairs()`）；`output_path` 必须不存在。四个字段编号是互不相同的 1 起始整数。`buffer_size` 形如 `50%` 或 `2G`；`temp_dir` 是已存在的文件夹；`sort_path` 是显式的程序路径。
- **输出：** `status`、`backend`（`"sort"`）、`command`、`return_code`、`outputs` 和 `log_tail`。
- **示例代码：**

```python
import tempfile

from dnakit.epigenome import cooler_csort_flip_pairs, cooler_csort_sort_pairs

with tempfile.TemporaryDirectory() as folder:
    fields = dict(chrom1=2, pos1=3, chrom2=4, pos2=5)
    cooler_csort_flip_pairs("pairs.txt", "chrom.sizes", f"{folder}/flipped.tsv", verbose=False, **fields)
    result = cooler_csort_sort_pairs(f"{folder}/flipped.tsv", f"{folder}/sorted.tsv", index="tabix", verbose=False, **fields)
    print(result["status"], result["return_code"])
    for line in open(f"{folder}/sorted.tsv").read().splitlines()[:3]:
        print(" ".join(line.split()))
```

- **示例输出：**

```text
completed 0
r6 chr1 3 chr1 95
r1 chr1 5 chr1 45
r2 chr1 10 chr1 30
```

## 15) `dnakit.epigenome.cooler_digest()`

- **作用：** 把基因组的限制性酶切片段箱写成 BED 文件。
- **处理方法：** 运行 `cooler digest --out`。cooler 在每条染色体中查找 Biopython 限制性内切酶的切点，每个片段写成一个箱，染色体顺序与 chromsizes 文件一致。结果另外包含行数和 BED 文件的前几行。
- **API：** `cooler_digest(chromsizes_path: 'str', fasta_path: 'str', enzyme: 'str', output_path: 'str', *, header: 'bool' = False, rel_ids: 'Literal[0, 1] | None' = None, cooler_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 cooler（额外依赖 `dnakit[bio-cooler]`），此子命令还需要 Biopython 和 pyfaidx；pyfaidx 可能会在 FASTA 文件旁创建 `.fai` 索引。示例打印 cooler 对 `chr1` 的两个 `AAGCTT` 位点报告的片段边界。
- **输入：** `chromsizes_path` 是已存在的 chromsizes 文件。`fasta_path` 是未压缩的基因组 FASTA 文件，或由 `*.fa` / `*.fasta` 文件组成的文件夹。`enzyme` 是 Biopython 的酶名，如 `HindIII` 或 `DpnII`。`header` 添加一行表头；`rel_ids`（0 或 1）添加每条染色体内的片段编号列。`output_path` 必须不存在。
- **输出：** 通用字段，加上 `line_count` 和 `preview`（最多五行，字符串列表）。
- **示例代码：**

```python
import tempfile

from dnakit.epigenome import cooler_digest

with tempfile.TemporaryDirectory() as folder:
    result = cooler_digest("chrom.sizes", "genome.fa", "HindIII", f"{folder}/fragments.bed", verbose=False)
    print(result["status"], result["line_count"])
    for line in result["preview"]:
        print(" ".join(line.split()))
```

- **示例输出：**

```text
completed 5
chr1 0 23
chr1 23 69
chr1 69 100
chr2 0 31
chr2 31 60
```

## 16) `dnakit.epigenome.cooler_dump()`

- **作用：** 把 cooler 的表（chroms、bins 或 pixels）导出为制表符分隔的文本。
- **处理方法：** 运行 `cooler dump --out`。对于 pixels，`region` 和 `region2` 选择 UCSC 风格的查询框，`balanced` 添加 `balanced` 列（weight1 * weight2 * count），`join` 把箱编号替换为染色体坐标，`annotate` 连接额外的箱表列（后缀为 1 和 2），`fill_lower` 为对称上三角存储的 cooler 生成隐含的下三角。`one_based_ids` 和 `one_based_starts` 只有与 `join`、`balanced` 或 `annotate` 一起使用时，cooler 才会应用。结果另外包含行数和文件的前几行。
- **API：** `cooler_dump(cool_uri: 'str', output_path: 'str', *, table: "Literal['chroms', 'bins', 'pixels']" = 'pixels', columns: 'Sequence[str] | None' = None, header: 'bool' = False, na_rep: 'str | None' = None, float_format: 'str | None' = None, region: 'str | None' = None, region2: 'str | None' = None, fill_lower: 'bool' = False, balanced: 'bool' = False, join: 'bool' = False, annotate: 'Sequence[str] | None' = None, one_based_ids: 'bool' = False, one_based_starts: 'bool' = False, chunksize: 'int | None' = None, cooler_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 cooler（额外依赖 `dnakit[bio-cooler]`）。
- **输入：** `cool_uri` 是已存在的 cooler；`output_path` 必须不存在。`table` 选择要导出的表。`columns`（普通列名）适用于 chroms 和 bins。`region2` 需要同时给出 `region`。`float_format` 形如 `g`、`.12g` 或 `03.2f`（默认 `g`）。给出 `na_rep` 时不能为空（cooler 的默认值为空）。`chunksize` 是正整数。
- **输出：** 通用字段，加上 `input_uri`、`table`、`line_count` 和 `preview`。
- **示例代码：**

```python
import tempfile

from dnakit.epigenome import cooler_cload_pairs, cooler_dump

with tempfile.TemporaryDirectory() as folder:
    cool = f"{folder}/a.cool"
    cooler_cload_pairs("chrom.sizes:20", "pairs.txt", cool, chrom1=2, pos1=3, chrom2=4, pos2=5, verbose=False)
    result = cooler_dump(cool, f"{folder}/pixels.tsv", join=True, header=True, verbose=False)
    print(result["status"], result["table"], result["line_count"])
    for line in result["preview"][:3]:
        print(" ".join(line.split()))
    bins = cooler_dump(cool, f"{folder}/bins.tsv", table="bins", columns=["chrom", "end"], verbose=False)
    print(bins["line_count"], bins["preview"])
```

- **示例输出：**

```text
completed pixels 8
chrom1 start1 end1 chrom2 start2 end2 count
chr1 0 20 chr1 20 40 1
chr1 0 20 chr1 40 60 1
8 ['chr1\t20', 'chr1\t40', 'chr1\t60', 'chr1\t80', 'chr1\t100']
```

## 17) `dnakit.epigenome.cooler_info()`

- **作用：** 返回 cooler 的信息块：箱大小、箱数和染色体数、像素总和。
- **处理方法：** 运行 `cooler info --out` 写入临时文件并解析 JSON。`field` 返回单个值（尽可能按 JSON 解析，否则按文本）；`metadata=True` 则改为返回用户元数据。
- **API：** `cooler_info(cool_uri: 'str', *, field: 'str | None' = None, metadata: 'bool' = False, cooler_path: 'str | None' = None, timeout_seconds: 'float' = 600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 cooler（额外依赖 `dnakit[bio-cooler]`）。
- **输入：** `cool_uri` 是已存在的 cooler 或 `file::/group`。`field` 是 `nbins` 之类的信息键；不能与 `metadata` 同时使用。
- **输出：** 通用字段，加上 `cool_uri`、`info`（所有字段）、`value`（请求的单个字段）和 `metadata`；未请求的为 `None`。
- **示例代码：**

```python
import tempfile

from dnakit.epigenome import cooler_cload_pairs, cooler_info

with tempfile.TemporaryDirectory() as folder:
    cool = f"{folder}/a.cool"
    cooler_cload_pairs("chrom.sizes:20", "pairs.txt", cool, chrom1=2, pos1=3, chrom2=4, pos2=5, verbose=False)
    info = cooler_info(cool, verbose=False)["info"]
    print({key: info[key] for key in ("bin-size", "bin-type", "nbins", "nchroms", "nnz", "storage-mode", "sum")})
    print(cooler_info(cool, field="nbins", verbose=False)["value"], cooler_info(cool, metadata=True, verbose=False)["metadata"])
```

- **示例输出：**

```text
{'bin-size': 20, 'bin-type': 'fixed', 'nbins': 8, 'nchroms': 2, 'nnz': 7, 'storage-mode': 'symmetric-upper', 'sum': 7}
8 {}
```

## 18) `dnakit.epigenome.cooler_ln()`

- **作用：** 创建指向 cooler 的链接，而不是复制它。
- **处理方法：** 运行 `cooler ln`：默认是同一文件内的硬链接；`soft=True` 时在同一文件内建软链接，目标是另一个文件时建外部链接。`overwrite=True` 传入 `-w`，它会截断目标文件（在同一文件内会失败，见 `cooler_cp()`）。不要把根组 `file::/` 链接到它自己的子组中：这会形成环，之后 cooler 列出该文件时会因 `RecursionError` 失败。
- **API：** `cooler_ln(source_uri: 'str', destination_uri: 'str', *, overwrite: 'bool' = False, soft: 'bool' = False, cooler_path: 'str | None' = None, timeout_seconds: 'float' = 600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 cooler（额外依赖 `dnakit[bio-cooler]`）。
- **输入：** `source_uri` 是已存在的 cooler（`file` 或 `file::/group`）；`destination_uri` 是位于已存在文件夹中的 `file` 或 `file::/group`。
- **输出：** 通用字段，加上 `source_uri`、`destination_uri` 和 `source_file`。
- **示例代码：**

```python
import tempfile

from dnakit.epigenome import cooler_cload_pairs, cooler_cp, cooler_ln, cooler_ls

with tempfile.TemporaryDirectory() as folder:
    cool = f"{folder}/a.cool"
    multi = f"{folder}/m.mcool"
    cooler_cload_pairs("chrom.sizes:20", "pairs.txt", cool, chrom1=2, pos1=3, chrom2=4, pos2=5, verbose=False)
    cooler_cp(cool, f"{multi}::resolutions/20", verbose=False)
    result = cooler_ln(f"{multi}::resolutions/20", f"{multi}::alias", verbose=False)
    print(result["status"], [uri.split("::")[1] for uri in cooler_ls(multi, verbose=False)["coolers"]])
```

- **示例输出：**

```text
completed ['/alias', '/resolutions/20']
```

## 19) `dnakit.epigenome.cooler_load()`

- **作用：** 由 COO 或 BG2 文本格式的预分箱矩阵创建 cooler。
- **处理方法：** 运行 `cooler load --format`。`coo` 的行是 `bin1_id`、`bin2_id`、`count`；`bg2` 的行是 `chrom1`、`start1`、`end1`、`chrom2`、`start2`、`end2`、`count`。一旦给出任何 `fields` 条目，只有把 count 列也列出来（`count=7` 或 `count:dtype=float`）才会加载它。
- **API：** `cooler_load(bins: 'str', pixels_path: 'str', output_path: 'str', *, input_format: "Literal['coo', 'bg2']", metadata_path: 'str | None' = None, assembly: 'str | None' = None, fields: 'Sequence[str] | None' = None, count_as_float: 'bool' = False, one_based: 'bool' = False, comment_char: 'str | None' = None, no_symmetric_upper: 'bool' = False, input_copy_status: "Literal['unique', 'duplex'] | None" = None, chunksize: 'int | None' = None, mergebuf: 'int | None' = None, temp_dir: 'str | None' = None, max_merge: 'int | None' = None, no_delete_temp: 'bool' = False, storage_options: 'str | None' = None, append: 'bool' = False, cooler_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 cooler（额外依赖 `dnakit[bio-cooler]`）。
- **输入：** `bins` 是 `<chromsizes 文件>:<箱大小 bp>` 或 BED 文件。`pixels_path` 是已存在的文本文件，可以是 gzip 压缩的（不支持标准输入）。`input_format` 为 `"coo"` 或 `"bg2"`。`fields` 的条目形如 `score=4:dtype=float`。`count_as_float` 把计数存成浮点数；`one_based` 表示 COO 文件中的箱编号是 1 起始。除非 `append=True`，`output_path` 必须不存在。
- **输出：** 通用字段，加上 `output_uri`。
- **示例代码：**

```python
import tempfile

from dnakit.epigenome import cooler_info, cooler_load

with tempfile.TemporaryDirectory() as folder:
    result = cooler_load("chrom.sizes:20", "pixels.coo", f"{folder}/l.cool", input_format="coo", verbose=False)
    info = cooler_info(f"{folder}/l.cool", verbose=False)["info"]
    print(result["status"], info["nbins"], info["nnz"], info["sum"])
```

- **示例输出：**

```text
completed 8 3 8
```

## 20) `dnakit.epigenome.cooler_ls()`

- **作用：** 列出 HDF5 文件中的所有 cooler，例如 `.mcool` 文件的各个分辨率。
- **处理方法：** 运行 `cooler ls`（需要时带 `--long`）并解析每一行。`long=True` 时报告每个 cooler 的箱大小，可变箱大小时为 `None`。
- **API：** `cooler_ls(cool_path: 'str', *, long: 'bool' = False, cooler_path: 'str | None' = None, timeout_seconds: 'float' = 600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 cooler（额外依赖 `dnakit[bio-cooler]`）。
- **输入：** `cool_path` 是已存在的 cooler 或多分辨率文件（普通路径，不是 `::group` URI）。
- **输出：** 通用字段，加上 `cool_path`、`coolers`（`file::/group` URI 的列表），以及 `long=True` 时的 `bin_sizes`（URI 到箱大小的映射）。
- **示例代码：**

```python
import tempfile

from dnakit.epigenome import cooler_cload_pairs, cooler_coarsen, cooler_cp, cooler_ls

with tempfile.TemporaryDirectory() as folder:
    cool = f"{folder}/a.cool"
    multi = f"{folder}/m.mcool"
    cooler_cload_pairs("chrom.sizes:20", "pairs.txt", cool, chrom1=2, pos1=3, chrom2=4, pos2=5, verbose=False)
    cooler_cp(cool, f"{multi}::resolutions/20", verbose=False)
    cooler_coarsen(cool, f"{multi}::resolutions/40", append=True, verbose=False)
    result = cooler_ls(multi, long=True, verbose=False)
    print({uri.split("::")[1]: size for uri, size in result["bin_sizes"].items()})
```

- **示例输出：**

```text
{'/resolutions/20': 20, '/resolutions/40': 40}
```

## 21) `dnakit.epigenome.cooler_makebins()`

- **作用：** 由 chromsizes 文件把固定宽度的基因组箱写成 BED 文件。
- **处理方法：** 运行 `cooler makebins --out`。`rel_ids` 添加每条染色体内的箱编号列，0 起始（0）或 1 起始（1）。结果另外包含行数和文件的前几行。
- **API：** `cooler_makebins(chromsizes_path: 'str', bin_size: 'int', output_path: 'str', *, header: 'bool' = False, rel_ids: 'Literal[0, 1] | None' = None, cooler_path: 'str | None' = None, timeout_seconds: 'float' = 600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 cooler（额外依赖 `dnakit[bio-cooler]`）。
- **输入：** `chromsizes_path` 是已存在的 chromsizes 文件；`bin_size` 是以 bp 为单位的正整数；`output_path` 必须不存在；`header` 添加一行表头。
- **输出：** 通用字段，加上 `line_count` 和 `preview`（最多五行）。
- **示例代码：**

```python
import tempfile

from dnakit.epigenome import cooler_makebins

with tempfile.TemporaryDirectory() as folder:
    result = cooler_makebins("chrom.sizes", 40, f"{folder}/bins.bed", rel_ids=1, verbose=False)
    print(result["status"], result["line_count"])
    for line in result["preview"]:
        print(" ".join(line.split()))
```

- **示例输出：**

```text
completed 5
chr1 0 40 1
chr1 40 80 2
chr1 80 100 3
chr2 0 40 1
chr2 40 60 2
```

## 22) `dnakit.epigenome.cooler_merge()`

- **作用：** 把箱完全相同的多个 cooler 合并为一个，像素相加。
- **处理方法：** 运行 `cooler merge`，输出在前，输入在后。`fields` 同 `cooler_coarsen()`。
- **API：** `cooler_merge(output_path: 'str', input_paths: 'Sequence[str]', *, chunksize: 'int | None' = None, fields: 'Sequence[str] | None' = None, append: 'bool' = False, cooler_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 cooler（额外依赖 `dnakit[bio-cooler]`）。
- **输入：** `input_paths` 是至少两个已存在的 cooler（文件或 URI）的列表，染色体和箱必须完全相同。`chunksize`（合并缓冲区，单位为像素行，cooler 默认 20,000,000）是正整数。除非 `append=True`，`output_path` 必须不存在。
- **输出：** 通用字段，加上 `output_uri` 和 `input_uris`。
- **示例代码：**

```python
import tempfile

from dnakit.epigenome import cooler_cload_pairs, cooler_cp, cooler_info, cooler_merge

with tempfile.TemporaryDirectory() as folder:
    cool = f"{folder}/a.cool"
    cooler_cload_pairs("chrom.sizes:20", "pairs.txt", cool, chrom1=2, pos1=3, chrom2=4, pos2=5, verbose=False)
    cooler_cp(cool, f"{folder}/b.cool", verbose=False)
    result = cooler_merge(f"{folder}/merged.cool", [cool, f"{folder}/b.cool"], verbose=False)
    info = cooler_info(f"{folder}/merged.cool", verbose=False)["info"]
    print(result["status"], len(result["input_uris"]), info["nnz"], info["sum"])
```

- **示例输出：**

```text
completed 2 7 14
```

## 23) `dnakit.epigenome.cooler_mv()`

- **作用：** 在同一个 HDF5 文件内重命名 cooler 组。
- **处理方法：** 运行 `cooler mv`。`overwrite=True` 传入 `-w`，cooler 0.10.4 无法在同一文件内执行它（见 `cooler_cp()`）；保留此参数只是为了完整。
- **API：** `cooler_mv(source_uri: 'str', destination_uri: 'str', *, overwrite: 'bool' = False, cooler_path: 'str | None' = None, timeout_seconds: 'float' = 600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 cooler（额外依赖 `dnakit[bio-cooler]`）。
- **输入：** `source_uri` 是已存在的 `file::/group`；`destination_uri` 是同一文件内的 `file::/group`。
- **输出：** 通用字段，加上 `source_uri`、`destination_uri` 和 `source_file`。
- **示例代码：**

```python
import tempfile

from dnakit.epigenome import cooler_cload_pairs, cooler_cp, cooler_ls, cooler_mv

with tempfile.TemporaryDirectory() as folder:
    cool = f"{folder}/a.cool"
    multi = f"{folder}/m.mcool"
    cooler_cload_pairs("chrom.sizes:20", "pairs.txt", cool, chrom1=2, pos1=3, chrom2=4, pos2=5, verbose=False)
    cooler_cp(cool, f"{multi}::resolutions/20", verbose=False)
    result = cooler_mv(f"{multi}::resolutions/20", f"{multi}::resolutions/renamed", verbose=False)
    print(result["status"], [uri.split("::")[1] for uri in cooler_ls(multi, verbose=False)["coolers"]])
```

- **示例输出：**

```text
completed ['/resolutions/renamed']
```

## 24) `dnakit.epigenome.cooler_pairs_header()`

- **作用：** 返回接触列表开头的注释行，纯文本或 gzip 压缩均可。
- **处理方法：** 在进程内运行，与 `cooler cload` 读取表头的方式相同。逐行读取，直到第一行不以 `comment_char` 开头为止；各行去掉行尾空格和换行，保留注释字符。`.gz` 文件用 gzip 打开。不需要 cooler 或其他程序。
- **API：** `cooler_pairs_header(pairs_path: 'str', *, comment_char: 'str' = '#', verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 安装后即可使用（普通 pip install）。
- **输入：** `pairs_path` 是已存在的文本或 `.gz` 文件；`comment_char` 是不含引号和反斜杠的短字符串。无法读取或不是 UTF-8 的文件会抛出 `ConfigurationError`。
- **输出：** `status`、`pairs_path`、`header`（行列表）、`header_line_count` 和 `has_data_lines`。
- **示例代码：**

```python
from dnakit.epigenome import cooler_pairs_header

result = cooler_pairs_header("pairs.txt", verbose=False)
print(result["header_line_count"], result["has_data_lines"])
for line in result["header"]:
    print(line)
```

- **示例输出：**

```text
2 True
## pairs format v1.0
#columns: readID chrom1 pos1 chrom2 pos2
```

## 25) `dnakit.epigenome.cooler_show()`

- **作用：** 不需要显示器，把 cooler 的一个接触矩阵区域画成图像文件。
- **处理方法：** 在 `MPLBACKEND=Agg` 下运行 `cooler show --out`。如果能导入 cooler 包，DNAKit 随后用 `Cooler.matrix(...).fetch()` 读取同一区域，并返回图像背后的数值：一份摘要，以及在不超过 40,000 个单元格时的矩阵本身（非有限值变为 `None`）。超过 4,000,000 个单元格时省略数值。颜色刻度只作用于图像，返回的矩阵不做变换。
- **API：** `cooler_show(cool_uri: 'str', region: 'str', output_path: 'str', *, region2: 'str | None' = None, balanced: 'bool' = False, dpi: 'int | None' = None, scale: "Literal['linear', 'log2', 'log10'] | None" = None, force: 'bool' = False, zmin: 'float | None' = None, zmax: 'float | None' = None, cmap: 'str | None' = None, field: 'str | None' = None, cooler_path: 'str | None' = None, timeout_seconds: 'float' = 1800.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 cooler（额外依赖 `dnakit[bio-cooler]`）和 matplotlib（额外依赖 `dnakit[bio-plots]`）。
- **输入：** `region` 是不含空白的 UCSC 风格区域（`chr1:10,000,000-11,000,000`）；`region2` 给出列范围，用于跨染色体或非对称视图。`output_path` 必须不存在，并以 `.png`、`.jpg`、`.jpeg`、`.svg`、`.pdf`、`.ps` 或 `.eps` 结尾。`scale` 为 `linear`、`log2` 或 `log10`（cooler 默认）；`dpi` 不小于 10；`cmap` 是 matplotlib 色图名（cooler 默认 `YlOrRd`）；`field` 是像素列（默认 `count`）；`zmin` 和 `zmax` 是有限数。
- **输出：** 通用字段，加上 `image_path`、`image_bytes`、`cool_uri`、`region`、`region2`、`balanced`、`scale`、`field` 和 `data`（`binsize`、`summary`（含 `shape`、`finite_cells`、`nan_cells`、`min`、`max`、`sum`、`mean`）以及 `matrix`）。
- **示例代码：**

```python
import tempfile

from dnakit.epigenome import cooler_cload_pairs, cooler_show

with tempfile.TemporaryDirectory() as folder:
    cool = f"{folder}/a.cool"
    cooler_cload_pairs("chrom.sizes:20", "pairs.txt", cool, chrom1=2, pos1=3, chrom2=4, pos2=5, verbose=False)
    result = cooler_show(cool, "chr1", f"{folder}/chr1.png", verbose=False)
    print(result["status"], result["image_bytes"] > 0, result["scale"], result["data"]["summary"]["shape"])
    print(result["data"]["matrix"][0])
    trans = cooler_show(cool, "chr1", f"{folder}/trans.png", region2="chr2", scale="linear", verbose=False)
    print(trans["data"]["summary"]["shape"], trans["data"]["summary"]["sum"])
```

- **示例输出：**

```text
completed True log10 [5, 5]
[0.0, 1.0, 1.0, 0.0, 1.0]
[5, 3] 2.0
```

## 26) `dnakit.epigenome.cooler_show_view()`

- **作用：** 不需要显示器，画出交互式 `cooler show` 窗口在缩放或平移后显示的热图。
- **处理方法：** 在进程内使用 cooler 和 matplotlib 包运行（不启动 `cooler` 程序）。与交互式窗口一样，它把窗口向外取整到整箱，裁剪到染色体范围内（最后一个不完整的箱被排除），并在窗口包含 `max_matrix_size` 个或更多单元格时画出占位文字而不是矩阵。不给 `view` 时画出整个 `region`。
- **API：** `cooler_show_view(cool_uri: 'str', region: 'str', output_path: 'str', *, region2: 'str | None' = None, view: 'Sequence[float] | None' = None, balanced: 'bool' = False, field: 'str' = 'count', scale: "Literal['linear', 'log2', 'log10']" = 'log10', cmap: 'str' = 'YlOrRd', zmin: 'float | None' = None, zmax: 'float | None' = None, dpi: 'int | None' = None, max_matrix_size: 'int' = 10000000, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 cooler（`dnakit[bio-cooler]`）、numpy（`dnakit[bio-analysis]`）和 matplotlib（`dnakit[bio-plots]`）。cooler 必须是固定箱大小。
- **输入：** `region` 和 `region2` 同 `cooler_show()`。`view` 是 `[列起点, 列终点, 行起点, 行终点]`，单位为 bp，分别位于 `region2`（默认 `region`）和 `region` 的染色体上，起点必须小于终点。`balanced=True` 需要 `weight` 列（先运行 `cooler_balance()`）。`output_path` 必须不存在，并带图像后缀。`max_matrix_size` 是正整数。
- **输出：** `status`、`backend`、`image_path`、`image_bytes`、`cool_uri`、`row_region` 和 `column_region`（取整后的 `[chrom, start, end]`）、`view_applied`、`placeholder`、`binsize`、`balanced`、`scale`、`field` 和 `data`（`summary`、`matrix`、`matrix_omitted`）。
- **示例代码：**

```python
import tempfile

from dnakit.epigenome import cooler_cload_pairs, cooler_show_view

with tempfile.TemporaryDirectory() as folder:
    cool = f"{folder}/a.cool"
    cooler_cload_pairs("chrom.sizes:20", "pairs.txt", cool, chrom1=2, pos1=3, chrom2=4, pos2=5, verbose=False)
    result = cooler_show_view(cool, "chr1", f"{folder}/view.png", view=[10, 70, 10, 70], verbose=False)
    print(result["row_region"], result["column_region"], result["view_applied"], result["placeholder"])
    print(result["data"]["matrix"])
```

- **示例输出：**

```text
['chr1', 0, 80] ['chr1', 0, 80] True False
[[0.0, 1.0, 1.0, 0.0], [1.0, 0.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0]]
```

## 27) `dnakit.epigenome.cooler_tree()`

- **作用：** 以文本形式返回 cooler 文件或组的 HDF5 数据层级。
- **处理方法：** 运行 `cooler tree`，可带 `--level`。cooler 打印的文本（`cooler.fileops.pprint_data_tree`）放在 `tree` 中。
- **API：** `cooler_tree(uri: 'str', *, level: 'int | None' = None, cooler_path: 'str | None' = None, timeout_seconds: 'float' = 600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 cooler（额外依赖 `dnakit[bio-cooler]`）。
- **输入：** `uri` 是已存在的 cooler 文件或 `file::/group`。`level` 是最大深度（不小于 0 的整数）。
- **输出：** 通用字段，加上 `uri` 和 `tree`（运行失败时为 `None`）。
- **示例代码：**

```python
import tempfile

from dnakit.epigenome import cooler_cload_pairs, cooler_tree

with tempfile.TemporaryDirectory() as folder:
    cool = f"{folder}/a.cool"
    cooler_cload_pairs("chrom.sizes:20", "pairs.txt", cool, chrom1=2, pos1=3, chrom2=4, pos2=5, verbose=False)
    print(cooler_tree(cool, level=1, verbose=False)["tree"])
```

- **示例输出：**

```text
/
 ├── bins
 ├── chroms
 ├── indexes
 └── pixels
```

## 28) `dnakit.epigenome.cooler_zoomify()`

- **作用：** 通过递归粗化构建多分辨率 cooler（`.mcool`）。
- **处理方法：** 运行 `cooler zoomify --out`。`resolutions` 是逗号分隔的列表，如 `1000,5000,25000`，或一种递进方式：`B`（二进制，每步翻倍）、`N`（“nice”，1-2-5 步进）、`<起点>N` 或 `4DN`；默认为 `B`。递进方式 `<起点>B`（如 `1000B`）在 cooler 文档中有说明，但在 0.10.4 中会失败，所以请用 `B` 或显式列表。`balance=True` 会平衡每一级缩放，`balance_options` 传入 cooler balance 的设置。
- **API：** `cooler_zoomify(cool_uri: 'str', output_path: 'str', *, resolutions: 'str | None' = None, nproc: 'int | None' = None, chunksize: 'int | None' = None, balance: 'bool' = False, balance_options: 'Mapping[str, Any] | None' = None, base_uris: 'Sequence[str] | None' = None, fields: 'Sequence[str] | None' = None, legacy: 'bool' = False, cooler_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 cooler（额外依赖 `dnakit[bio-cooler]`）。
- **输入：** `cool_uri` 是已存在的 cooler；`output_path` 必须不存在。`balance_options` 需要 `balance=True`，可包含 `cis_only`、`trans_only`、`ignore_diags`、`ignore_dist`、`mad_max`、`min_nnz`、`min_count`、`nproc`、`max_iters` 和 `tol`。`base_uris` 是用于聚合的额外基础 cooler。`legacy=True` 写出旧布局中以整数标记的缩放级别。`nproc` 和 `chunksize` 是正整数。
- **输出：** 通用字段，加上 `input_uri`、`output_path`、`legacy` 和 `balanced`。
- **示例代码：**

```python
import tempfile

from dnakit.epigenome import cooler_cload_pairs, cooler_ls, cooler_zoomify

with tempfile.TemporaryDirectory() as folder:
    cool = f"{folder}/a.cool"
    cooler_cload_pairs("chrom.sizes:20", "pairs.txt", cool, chrom1=2, pos1=3, chrom2=4, pos2=5, verbose=False)
    result = cooler_zoomify(cool, f"{folder}/z.mcool", resolutions="20,40,80", verbose=False)
    listing = cooler_ls(f"{folder}/z.mcool", long=True, verbose=False)
    print(result["status"], result["legacy"], result["balanced"])
    print({uri.split("::")[1]: size for uri, size in listing["bin_sizes"].items()})
```

- **示例输出：**

```text
completed False False
{'/resolutions/20': 20, '/resolutions/40': 40, '/resolutions/80': 80}
```

## 29) `dnakit.epigenome.cooler_resolve_bins()`

- **作用：** 解析 cooler 的 `BINS` 参数，并把箱表写入新的 BED 文件。
- **处理方法：** 遵循 `cooler.cli._util.parse_bins`（cooler 0.10.4）。含冒号的字符串是 `<chromsizes 文件>:<箱大小 bp>`，由 `cooler.util.binnify` 生成固定大小的箱；否则它是定义箱的三列 BED 文件（chrom、start、end）的路径。无法解析的 BED 文件会抛出 `ConfigurationError`（在 pandas 2 下上游此处抛出 `AttributeError`）。表写出时不带表头。
- **API：** `cooler_resolve_bins(bins: 'str', output_bed: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 cooler（额外依赖 `dnakit[bio-cooler]`）；pandas 是 cooler 的依赖。
- **输入：** `bins` 是上述形式的非空字符串，chromsizes 或 BED 文件必须存在，箱大小为正整数。`output_bed` 必须不存在。
- **输出：** `operation`、`status`、`backend`、`package`、`source`（`"chromsizes"` 或 `"bed"`）、`binsize`（BED 来源时为 `None`）、`bin_count`、`chromosomes`（按首次出现的顺序）、`output_bed` 和 `coordinates`（`"0-based half-open BED"`）。
- **示例代码：**

```python
import tempfile

from dnakit.epigenome import cooler_resolve_bins

with tempfile.TemporaryDirectory() as folder:
    result = cooler_resolve_bins("chrom.sizes:40", f"{folder}/bins.bed", verbose=False)
    print(result["source"], result["binsize"], result["bin_count"], result["chromosomes"])
    print(open(f"{folder}/bins.bed").read().split("\n")[0].split("\t"))
    again = cooler_resolve_bins(f"{folder}/bins.bed", f"{folder}/bins2.bed", verbose=False)
    print(again["source"], again["binsize"], again["bin_count"])
```

- **示例输出：**

```text
chromsizes 40 5 ['chr1', 'chr2']
['chr1', '0', '40']
bed None 5
```

## 30) `dnakit.epigenome.download_cooltools_example_data()`

- **作用：** 把一个 cooltools 示例数据集或全部数据集下载到指定文件夹。
- **处理方法：** 调用 `cooltools.lib.io.download_data()`。`cache=True` 时，`data_dir` 中已存在的文件会被直接复用，不再下载，但会先把它的 md5 校验和与数据集表中的比较（除非 `ignore_checksum` 为真）。否则通过 HTTP(S) 获取文件。
- **API：** `download_cooltools_example_data(name: 'str' = 'all', *, data_dir: 'str', cache: 'bool' = True, ignore_checksum: 'bool' = False, info_url: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 cooltools（额外依赖 `dnakit[bio-cooltools]`），未缓存的文件还需要联网；示例使用本地表和已缓存的文件，所以不需要网络。
- **输入：** `name` 是数据集键或 `"all"`；`data_dir` 是文件夹（不存在则创建）。`cache` 和 `ignore_checksum` 是布尔值。`info_url` 是 http(s) URL 或本地数据集表的路径（默认使用 cooltools 发布的表）。表以 `# ` 开头的一行为表头，其后是以制表符分隔的行：键、说明、链接、文件名和 md5 校验和。
- **输出：** `status`（`"complete"`）、`operation`、`service`、`backend`、`version`、`request`、`retrieved_utc`，以及 `data`，其中有 `files`（路径）、`file_count` 和 `data_dir`。未知的数据集或校验和不符会抛出 `BackendExecutionError`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import download_cooltools_example_data
from dnakit.exceptions import BackendExecutionError

with tempfile.TemporaryDirectory() as folder:
    Path(folder, "tiny.txt").write_text("hello\n")
    result = download_cooltools_example_data("tiny", data_dir=folder, info_url="datasets.tsv", verbose=False)
    print(result["status"], result["data"]["file_count"], [Path(path).name for path in result["data"]["files"]])
    try:
        download_cooltools_example_data("missing", data_dir=folder, info_url="datasets.tsv", verbose=False)
    except BackendExecutionError as error:
        print(type(error).__name__, "is not available" in str(error))
```

- **示例输出：**

```text
complete 1 ['tiny.txt']
BackendExecutionError True
```

联网时的用法：

```python
from dnakit.epigenome import download_cooltools_example_data

result = download_cooltools_example_data("HFF_MicroC", data_dir="cooltools_data")
print(result["data"]["files"])
```

## 31) `dnakit.epigenome.download_cooltools_file()`

- **作用：** 把一个 http(s) URL 上的文件下载到新的本地路径。
- **处理方法：** 调用 `cooltools.lib.io.download_file()`（以流的方式 HTTP GET 到文件），然后报告文件的大小和 md5 校验和。没有写入任何数据时抛出 `BackendExecutionError`。
- **API：** `download_cooltools_file(url: 'str', output_path: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 cooltools（额外依赖 `dnakit[bio-cooltools]`）和网络。示例不需要网络：其他协议的 URL 会在下载之前就被拒绝。
- **输入：** `url` 是 http 或 https URL（不超过 2,048 个字符）；`output_path` 必须不存在（其文件夹不存在时会被创建）。
- **输出：** 与上一节相同的外层字段（`status`、`operation`、`service`、`backend`、`version`、`request`、`retrieved_utc`），以及 `data`，其中有 `output_path`、`bytes` 和 `md5`。
- **示例代码：**

```python
from dnakit.epigenome import download_cooltools_file
from dnakit.exceptions import ConfigurationError

try:
    download_cooltools_file("ftp://example.org/data.cool", "data.cool", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, "url must be an http or https URL" in str(error))
```

- **示例输出：**

```text
ConfigurationError True
```

联网时的用法：

```python
from dnakit.epigenome import download_cooltools_file

result = download_cooltools_file("https://example.org/data.txt", "data.txt")
print(result["data"]["bytes"], result["data"]["md5"])
```

## 32) `dnakit.epigenome.list_cooltools_example_datasets()`

- **作用：** 返回 `print_available_datasets` 所打印的 cooltools 示例数据集表。
- **处理方法：** 用 cooltools 自己的解析器（`cooltools.lib.io._get_datasets_info`）从 `url` 读取表；`url` 为 `None` 时读取 cooltools 发布的表。每一行包含数据集键、说明、来源链接、本地文件名和 md5 校验和，另有行号 `index`。
- **API：** `list_cooltools_example_datasets(*, url: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 cooltools（额外依赖 `dnakit[bio-cooltools]`）。不给 `url` 时需要联网；使用本地表时可离线运行。
- **输入：** `url` 是 http(s) URL，或已存在的本地表的路径，表有以 `# ` 开头的表头行和以制表符分隔的行。
- **输出：** `status`（`"complete"`）、`operation`、`service`、`backend`、`version`、`request`、`retrieved_utc`，以及 `data`，其中有 `datasets`（行字典的列表）和 `count`。
- **示例代码：**

```python
from dnakit.epigenome import list_cooltools_example_datasets

result = list_cooltools_example_datasets(url="datasets.tsv", verbose=False)
print(result["status"], result["operation"], result["data"]["count"])
for row in result["data"]["datasets"]:
    print(row["index"], row["key"], row["filename"], row["checksum"][:8])
```

- **示例输出：**

```text
complete list_cooltools_example_datasets 2
1 tiny tiny.txt b1946ac9
2 small small.txt 00000000
```

## 33) `dnakit.epigenome.plot_cooltools_saddle()`

- **作用：** 把轨道相对于 saddle 矩阵的 saddle 图画进新的 PNG 文件。
- **处理方法：** 用 `cooltools.api.saddle.digitize` 对轨道离散化，调用 `cooltools.api.saddle.saddleplot`（它绘制在 pyplot 的当前图上），以 150 dpi 保存该图，并关闭本次调用打开的所有图。返回的计数来自离散化后的轨道；热图数值是去掉两侧各一行一列“外侧”后的 saddle 矩阵。
- **API：** `plot_cooltools_saddle(track_rows: 'Sequence[Sequence[object]]', saddle_matrix: 'Sequence[Sequence[object]]', n_bins: 'int', output_png: 'str', *, vrange: 'Sequence[float] | None' = None, qrange: 'Sequence[float] | None' = (0.0, 1.0), scale: 'str' = 'log', vmin: 'float' = 0.5, vmax: 'float' = 2.0, cmap: 'str' = 'coolwarm', title: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 cooltools（额外依赖 `dnakit[bio-cooltools]`）和 matplotlib（额外依赖 `dnakit[bio-plots]`）。
- **输入：** `track_rows` 是非空的 `[chrom, start, end, value]` 行列表（位置为非负整数，end 大于 start；缺失值用 `None`）。`saddle_matrix` 是 `n_bins + 2` 行的方阵，包含 `cooltools.api.saddle` 返回的两个外侧行列。`vrange`（数值范围）和 `qrange`（分位数范围，默认 `(0.0, 1.0)`）必须恰好给出一个，所以使用 `vrange` 时要传入 `qrange=None`。`scale` 为 `"log"` 或 `"linear"`；`vmin` 必须小于 `vmax`（`log` 时还必须大于 0）。`output_png` 必须不存在，并以 `.png` 结尾。
- **输出：** `status`、`operation`、`service`、`backend`、`version`、`image_path`、`data` 和 `upstream_warnings`。`data` 包含 `n_bins`、`vrange`、`qrange`、`scale`、`vmin`、`vmax`、`cmap`、`track_rows`、`bin_edges`、`bin_counts`、`left_outlier_count`、`right_outlier_count`、`missing_count`、`heatmap_shape` 和 `heatmap_values`。
- **示例代码：**

```python
import tempfile

from dnakit.epigenome import plot_cooltools_saddle

track = [["chr1", 0, 20, 0.1], ["chr1", 20, 40, 0.5], ["chr1", 40, 60, -0.3],
         ["chr1", 60, 80, 0.9], ["chr1", 80, 100, -0.8], ["chr2", 0, 20, 0.2]]
saddle = [[0] * 5, [0, 1.2, 0.8, 0.5, 0], [0, 0.8, 1.0, 0.9, 0], [0, 0.5, 0.9, 1.5, 0], [0] * 5]
with tempfile.TemporaryDirectory() as folder:
    result = plot_cooltools_saddle(track, saddle, 3, f"{folder}/saddle.png", vrange=[-1, 1], qrange=None, verbose=False)
    data = result["data"]
    print(result["status"], data["bin_counts"], data["left_outlier_count"], data["right_outlier_count"])
    print([round(edge, 3) for edge in data["bin_edges"]], data["heatmap_shape"])
    print(data["heatmap_values"][0])
```

- **示例输出：**

```text
completed [1, 3, 2] 0 0
[-1.0, -0.333, 0.333, 1.0] [3, 3]
[1.2, 0.8, 0.5]
```
