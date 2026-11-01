# 变异与注释 I/O：bed-reader、bgen-reader、hgvs、polars-bio 和 pyensembl

这一组函数位于 `dnakit.variants` 与 `dnakit.annotation`，是五个可选软件包的适配器：[bed-reader](https://github.com/fastlmm/bed-reader)（PLINK 1 的 `.bed` 文件，也支持 URL）、[bgen-reader](https://github.com/limix/bgen-reader-py)（BGEN 示例文件）、[hgvs](https://github.com/biocommons/hgvs)（HGVS 变异字符串与序列上下文文本）、[polars-bio](https://github.com/biodatageeks/polars-bio)（BGEN、PLINK 2 PGEN、VCF Zarr 和 CRAM 表，以及区间图）和 [pyensembl](https://github.com/openvax/pyensembl)（Ensembl 发布版缓存与注释索引）。它们把基因型文件读成普通字典、格式化 HGVS 文本、把一个变异映射为 g.、c.、n.、p. 四种写法、列出 Ensembl 注释日期，并构建或检查 pyensembl 注释缓存。它们不做变异检出，也不预测变异效应。

函数用 `from dnakit.variants import ...`（第 1 至 10 节）和 `from dnakit.annotation import ...`（第 11 至 33 节）导入。按已安装包的元数据，bed-reader、hgvs、polars-bio 和 pyensembl 使用 Apache-2.0 许可，bgen-reader 使用 MIT 许可，都不是 copyleft 许可。依赖都是可选的，只有调用对应函数时才会导入；缺少包时会抛出带安装提示的 `BackendUnavailableError`。

有几个函数需要的不只是普通 Python 包：bed-reader 的两个示例数据函数和 pyensembl 的三个列表函数会访问网络，`pyensembl_install_ensembl_release()` 会从 Ensembl FTP 服务器下载，两个 HGVS 映射函数需要建立在 UTA 数据库上的 hgvs 映射器。对这些函数，示例展示的是函数对无效输入的真实行为，并另附一段未执行的代码，展示依赖齐全时的调用。

## 示例数据

BGEN 和 PGEN 文件是二进制的，所以这里以十六进制文本声明，由示例代码写入磁盘。`cohort.bgen.hex` 是一个 150 字节的 BGEN v1.2 文件（layout 2，8 位概率，无压缩），含三个样本（S1、S2、S3）和两个变异（`rs1` 位于 1 号染色体 100 位，`rs2` 位于 250 位）。`cohort.pgen.hex` 是 PLINK 1 模式的 PGEN，样本和变异相同：三个魔数字节 `6c1b01` 加上每个变异一个打包字节。`cohort.pvar` 和 `cohort.psam` 是它的配套文件。`make_vcz.py` 用 `zarr` 包写出一个含两个变异、两个样本的 VCF Zarr 存储。`a.bed` 和 `b.bed` 是两个区间有重叠的小 BED 文件，`mini.gtf` 是一个含三个外显子的 GTF 注释。

**示例文件：** `cohort.bgen.hex`

```text
280000001400000002000000030000006267656e08000080140000000300
000002005331020053320200533304007661723003007273310100316400
000002000100000041010000004713000000030000000200020202020200
08ff0000ff00000400766172310300727332010031fa0000000200010000
0043010000005413000000030000000200020202020200080000ff0000ff
```

**示例文件：** `cohort.pgen.hex`

```text
6c1b010b1e
```

**示例文件：** `cohort.pvar`

```text
#CHROM	POS	ID	REF	ALT
1	100	rs1	A	G
1	250	rs2	C	T
```

**示例文件：** `cohort.psam`

```text
#IID	SEX
S1	NA
S2	NA
S3	NA
```

**示例文件：** `make_vcz.py`

```text
import numpy as np
import zarr
from zarr.core.dtype import VariableLengthUTF8


def make_vcz(path):
    """Write a 2-variant, 2-sample VCF Zarr store (Zarr v2 layout)."""
    root = zarr.open_group(path, mode="w", zarr_format=2)
    root.attrs["vcf_zarr_version"] = "0.4"
    root.attrs["vcf_header"] = (
        "##fileformat=VCFv4.3\n"
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
    )
    root.attrs["contigs"] = [{"ID": "chr1", "length": 1000}]

    def put(name, data, dtype=None):
        data = np.asarray(data, dtype=dtype)
        zdtype = VariableLengthUTF8() if data.dtype == object else data.dtype
        array = root.create_array(name, shape=data.shape, dtype=zdtype, chunks=data.shape)
        array[...] = data

    put("contig_id", ["chr1"], object)
    put("sample_id", ["S1", "S2"], object)
    put("variant_contig", [0, 0], np.int32)
    put("variant_position", [100, 250], np.int32)
    put("variant_id", ["rs1", "rs2"], object)
    put("variant_allele", [["A", "G"], ["C", "T"]], object)
    put("variant_quality", [50.0, 60.0], np.float32)
    put("variant_filter", np.zeros((2, 1), dtype=bool))
    put("filter_id", ["PASS"], object)
    put("variant_DP", [14, 20], np.int32)
    genotype = np.array([[[0, 0], [0, 1]], [[1, 1], [0, 0]]], dtype=np.int8)
    put("call_genotype", genotype)
    put("call_genotype_mask", genotype < 0)
    put("call_genotype_phased", np.zeros((2, 2), dtype=bool))
```

**示例文件：** `a.bed`

```text
chr1	10	20	a
chr1	40	60	b
```

**示例文件：** `b.bed`

```text
chr1	15	30	c
chr1	50	55	d
```

**示例文件：** `mini.gtf`

```text
chr1	toy	exon	100	300	.	+	.	gene_id "G1"; transcript_id "T1"; exon_number "1"; gene_name "GENEA";
chr1	toy	exon	600	900	.	+	.	gene_id "G1"; transcript_id "T1"; exon_number "2"; gene_name "GENEA";
chr1	toy	exon	2000	2500	.	-	.	gene_id "G2"; transcript_id "T2"; exon_number "1"; gene_name "GENEB";
```

## 1) `dnakit.variants.check_cloud_location()`

- **作用：** 报告 bed-reader 能否读取某个 URL 上的对象。
- **处理方法：** 先校验 URL，再调用 bed-reader 的原生函数 `check_file_cloud(url, options)`。原生检查对缺失或不可读对象抛出的 `ValueError` 会被报告为 `False`，其他错误则继续抛出。
- **API：** `check_cloud_location(location: 'str', *, cloud_options: 'dict[str, str] | None' = None, verbose: bool = True) -> 'bool'`。
- **运行要求：** 需要可选包 bed-reader（extra `dnakit[bio-bed-reader]`）。`file://` URL 读取本地文件；`https://` 和 `s3://` URL 需要联网，私有存储桶还需要在 `cloud_options` 中提供凭据。
- **输入：** `location` 必须是带协议（`https://`、`s3://`、`file://`）、长度不超过 4096 的字符串；普通本地路径会被 `ConfigurationError` 拒绝。`cloud_options` 是可选的字符串字典，传给存储客户端。
- **输出：** 一个布尔值，对象可以打开读取时为 `True`。
- **示例代码：**

```python
from pathlib import Path

import numpy as np
from bed_reader import to_bed

from dnakit.variants import check_cloud_location

to_bed("tiny.bed", np.array([[1, 0, np.nan], [2, 0, 2]], dtype=np.float32))
url = Path("tiny.bed").resolve().as_uri()
print(check_cloud_location(url, verbose=False))
print(check_cloud_location(url.replace("tiny", "absent"), verbose=False))
```

- **示例输出：**

```text
True
False
```

## 2) `dnakit.variants.read_cloud_genotypes()`

- **作用：** 读取以 URL 给出的 PLINK 1 `.bed` 文件的基因型矩阵。
- **处理方法：** 用 bed-reader 的 `open_bed(url, cloud_options=...)` 打开 URL，再调用 `read(index=(iid, sid), dtype=dtype)`，底层使用原生的 `read_cloud_i8`、`read_cloud_f32` 或 `read_cloud_f64`。缺失的基因型（`int8` 中的 -127、浮点中的 NaN）变成 `None`。读取之前先对所选单元格数与 `max_cells` 做检查。
- **API：** `read_cloud_genotypes(location: 'str', *, dtype: 'str' = 'float32', iid_index: 'list[int] | None' = None, sid_index: 'list[int] | None' = None, cloud_options: 'dict[str, str] | None' = None, max_cells: 'int' = 5000000, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 bed-reader 和 NumPy（extra `dnakit[bio-bed-reader]`）。`.fam` 和 `.bim` 文件必须与 `.bed` 对象放在同一位置。远程 URL 需要联网。
- **输入：** `location` 同上。`dtype` 为 `"int8"`、`"float32"` 或 `"float64"`。`iid_index`（样本）和 `sid_index`（变异）是非空的、范围内的 0 起始位置列表，`None` 表示全部。`max_cells` 上限为 100,000,000；所选单元格超过该值时抛出 `ConfigurationError`。
- **输出：** 字典，含 `status`、`operation`、`location`、`dtype`、`shape`（样本数乘变异数）、`missing_calls`、`values`（嵌套列表，缺失为 `None`）和 `provenance`。
- **示例代码：**

```python
from pathlib import Path

import numpy as np
from bed_reader import to_bed

from dnakit.variants import read_cloud_genotypes

to_bed("tiny.bed", np.array([[1, 0, np.nan], [2, 0, 2]], dtype=np.float32))
url = Path("tiny.bed").resolve().as_uri()
result = read_cloud_genotypes(url, verbose=False)
print(result["shape"], result["missing_calls"], result["values"])
subset = read_cloud_genotypes(url, dtype="int8", iid_index=[1], sid_index=[2, 0], verbose=False)
print(subset["dtype"], subset["shape"], subset["values"])
```

- **示例输出：**

```text
[2, 3] 1 [[1.0, 0.0, None], [2.0, 0.0, 2.0]]
int8 [1, 2] [[2, 2]]
```

## 3) `dnakit.variants.sample_bed_file()`

- **作用：** 下载 bed-reader 自带的某个示例 `.bed` 文件（未缓存时）并返回其本地路径。
- **处理方法：** 先校验名称，再调用 `bed_reader.sample_file(name)`，由 pooch 把文件取到 bed-reader 的缓存中。名称以 `.bed` 结尾时，若对应的 `.fam` 和 `.bim` 文件存在，会把它们的路径一并返回。缺少 pooch 时抛出 `BackendUnavailableError`，并给出 `bed-reader[samples]` 的安装提示。
- **API：** `sample_bed_file(filepath: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 bed-reader 和 pooch（extra `dnakit[bio-bed-reader]`），首次获取示例文件时需要联网。下面的示例展示的是在下载之前执行的名称检查；依赖和网络齐备时，请使用用法代码块中的调用。
- **输入：** `filepath` 是 bed-reader 示例文件的名称，例如 `small.bed`：非空、不超过 255 个字符、不含控制字符的字符串。
- **输出：** `service`、`source`（名称）和 `result`，其中含 `bed`、`fam`、`bim`（绝对路径，`fam` 和 `bim` 可能为 `None`）。
- **示例代码：**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.variants import sample_bed_file

try:
    sample_bed_file("", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error.code)
    print(error.message)
```

- **示例输出：**

```text
ConfigurationError CONFIGURATION_ERROR
filepath must be the name of a bed-reader sample file, such as small.bed.
```

安装依赖后的用法：

```python
from dnakit.variants import sample_bed_file

result = sample_bed_file("small.bed")
print(result["result"]["bed"], result["result"]["fam"], result["result"]["bim"])
```

## 4) `dnakit.variants.sample_bed_url()`

- **作用：** 返回 bed-reader 自带的某个示例 `.bed` 文件的 URL。
- **处理方法：** 先校验名称，再调用 `bed_reader.sample_url(name)`；文件未缓存时先由 pooch 获取，然后返回 `file://` URL。缺少 pooch 时抛出 `BackendUnavailableError`。
- **API：** `sample_bed_url(filepath: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 bed-reader 和 pooch（extra `dnakit[bio-bed-reader]`），首次获取示例文件时需要联网。示例展示名称检查，用法代码块展示联网时的调用。
- **输入：** `filepath`，校验规则与 `sample_bed_file()` 相同。
- **输出：** `service`、`source` 和只含一个键 `url` 的 `result`。该 URL 可以传给 `read_cloud_genotypes()`。
- **示例代码：**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.variants import sample_bed_url

try:
    sample_bed_url("bad\nname.bed", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error.code)
    print(error.message)
```

- **示例输出：**

```text
ConfigurationError CONFIGURATION_ERROR
filepath must not contain control characters.
```

安装依赖后的用法：

```python
from dnakit.variants import read_cloud_genotypes, sample_bed_url

url = sample_bed_url("small.bed")["result"]["url"]
print(read_cloud_genotypes(url)["shape"])
```

## 5) `dnakit.variants.fetch_bgen_example()`

- **作用：** 返回 bgen-reader 自带的某个示例 BGEN 文件的本地路径，首次使用时下载。
- **处理方法：** 先用 `bgen_reader._example` 中保存的名单检查名称，再调用 `bgen_reader.example_filepath(name)`，它会从 GitHub 上的 `fastlmm/bgen-sample-files` 仓库把文件下载到 bgen-reader 的缓存并校验哈希。结果中另外加上文件大小和来源 URL。
- **API：** `fetch_bgen_example(filename: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 bgen-reader（extra `dnakit[bio-bgen-reader]`，同时会安装 `requests`），首次获取文件时需要联网。示例展示对未知名称的拒绝，这一步发生在任何下载之前。
- **输入：** `filename` 必须是 bgen-reader 认识的名称，例如 `example.bgen` 或 `complex.bgen`；否则抛出 `ConfigurationError`，并列出已知名称。
- **输出：** `status`、`backend`、`name`、`path`（本地绝对路径）、`size_bytes` 和 `source`（下载 URL）。
- **示例代码：**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.variants import fetch_bgen_example

try:
    fetch_bgen_example("nope.bgen", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__)
    print(error.message.split(" Known names")[0])
```

- **示例输出：**

```text
ConfigurationError
Unknown bgen_reader example file 'nope.bgen'.
```

安装依赖后的用法：

```python
from dnakit.variants import fetch_bgen_example

result = fetch_bgen_example("example.bgen")
print(result["path"], result["size_bytes"])
```

## 6) `dnakit.variants.format_sequence_blocks()`

- **作用：** 把 DNA 序列按碱基分组排版，并附上坐标标尺。
- **处理方法：** 检查参数后调用 hgvs 的 `hgvs.utils.context.format_sequence(sequence, start, end, group_size=...)`，把它逐行产出的内容连成一个字符串。每行宽 100 列：一个坐标列、一行标注每组最后一个碱基位置的标尺，以及碱基本身。显示的坐标从 1 开始。
- **API：** `format_sequence_blocks(sequence: 'str', *, start: 'int | None' = None, end: 'int | None' = None, group_size: 'int' = 3, verbose: bool = True) -> 'str'`。
- **运行要求：** 需要可选包 hgvs（extra `dnakit[bio-hgvs]`）。只使用 hgvs 的排版代码，不需要数据库。
- **输入：** `sequence` 是非空字符串。`start` 和 `end` 是 0 起始的切片偏移，需满足 `0 <= start < end <= len(sequence)`（默认取整条序列）；显示的编号从 1 开始。`group_size` 是 1 到 88 的整数。
- **输出：** 一个多行字符串，行块之间以空行分隔。
- **示例代码：**

```python
from dnakit.variants import format_sequence_blocks

text = format_sequence_blocks("ACGTACGTACGTACGTAC", verbose=False)
part = format_sequence_blocks("ACGTACGTACGTACGTAC", start=3, end=12, group_size=4, verbose=False)
for block in (text, part):
    for line in block.splitlines():
        if line.strip():
            print("|" + line.rstrip())
```

- **示例输出：**

```text
|          :   3   6   9  12  15  18  21  24  27  30  33  36  39  42  45  48  51  54  57  60
|        1 : ACG TAC GTA CGT ACG TAC
|          :    7   11   15   19   23   27   31   35   39   43   47   51   55   59   63   67   71
|        4 : TACG TACG T
```

## 7) `dnakit.variants.hgvs_full_house()`

- **作用：** 把一个变异表示为 g.、c.、n.、p. 四种 HGVS 字符串。
- **处理方法：** 用 hgvs 解析器解析变异，再调用 `hgvs.utils.context.full_house(mapper, variant, tx_ac=transcript)`。对未给出 `transcript` 的基因组变异，它向映射器请求 `relevant_transcripts`，并要求恰好有一个；随后调用 `g_to_n`、`n_to_c` 和 `c_to_p`（输入为 n. 或 c. 时按需调用 `n_to_g`、`c_to_g`、`c_to_n`）。hgvs 抛出的 `RuntimeError` 会被转成 `ConfigurationError`。
- **API：** `hgvs_full_house(mapper: 'Any', variant: 'str', *, transcript: 'str | None' = None, verbose: bool = True) -> 'dict[str, str]'`。
- **运行要求：** 需要可选包 hgvs（extra `dnakit[bio-hgvs]`），以及建立在 UTA 数据库连接和序列库上的 hgvs 映射器；hgvs 两者都不自带，UTA 数据库要通过网络或本地副本访问。docstring 把映射器称为 `AlignmentMapper`，但实际用到的方法（`relevant_transcripts`、`g_to_n`、`n_to_c`、`c_to_p`）属于 hgvs 的 `AssemblyMapper`。示例用一个列出两个转录本的简单替身对象来展示拒绝行为；用法代码块展示真实映射器。
- **输入：** `mapper`（必填，不能为 `None`）；`variant` 是 HGVS 字符串，例如 `NC_000001.11:g.100A>G`；`transcript` 在映射器对某个基因组变异报告出多个重叠转录本时必须给出。
- **输出：** 一个字典，含键 `g`、`c`、`n`、`p`，值均为字符串。
- **示例代码：**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.variants import hgvs_full_house


class TwoTranscripts:
    """Stand-in mapper: two transcripts overlap every variant."""

    def relevant_transcripts(self, variant):
        return ["NM_000001.1", "NM_000002.1"]


try:
    hgvs_full_house(TwoTranscripts(), "NC_000001.11:g.100A>G", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__)
    print(error.message)
```

- **示例输出：**

```text
ConfigurationError
2 relevant transcripts for NC_000001.11; you need to pick one
```

安装依赖后的用法：

```python
import hgvs.assemblymapper
import hgvs.dataproviders.uta

from dnakit.variants import hgvs_full_house

hdp = hgvs.dataproviders.uta.connect()  # needs the UTA database and a sequence store
mapper = hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name="GRCh37")
print(hgvs_full_house(mapper, "NC_000007.13:g.117199646_117199648del", transcript="NM_000492.3"))
```

## 8) `dnakit.variants.hgvs_pointer_line()`

- **作用：** 画出在所显示序列区间内指向变异的标记行。
- **处理方法：** 解析变异并调用 `hgvs.utils.context.pointer_line(variant, (start, end))`。插入用 `><` 标记，其他编辑每个碱基一个 `*`。行末是变异文本。
- **API：** `hgvs_pointer_line(variant: 'str', span: 'Sequence[int]', *, verbose: bool = True) -> 'str'`。
- **运行要求：** 需要可选包 hgvs（extra `dnakit[bio-hgvs]`）。只使用 hgvs 解析器，不需要数据库。
- **输入：** `variant` 是带位置的基因组、编码或非编码 HGVS 变异（g.、c. 或 n.，其他类型会抛出 `ConfigurationError`）。`span` 是整数对，需满足 `0 <= start <= end`。
- **输出：** 一行文本。
- **示例代码：**

```python
from dnakit.variants import hgvs_pointer_line

for variant in ("NM_000399.3:c.5A>G", "NM_000399.3:c.5_6insA"):
    line = hgvs_pointer_line(variant, (0, 10), verbose=False)
    print("|" + line.rstrip())
```

- **示例输出：**

```text
|                                    *                                         NM_000399.3:c.5A>G
|                                    ><                                        NM_000399.3:c.5_6insA
```

## 9) `dnakit.variants.hgvs_sequence_line()`

- **作用：** 格式化变异上下文显示中的一行序列。
- **处理方法：** 解析变异并调用 `hgvs.utils.context.seq_line_fmt(var=..., span=..., content=..., dir=...)`。这一行左侧是登录号和变异类型，其后依次是起始坐标、链方向标记、序列文本、再次出现的链方向标记、结束坐标，以及作为行尾注释的变异字符串。
- **API：** `hgvs_sequence_line(variant: 'str', span: 'Sequence[int]', content: 'str', *, direction: 'str' = '', verbose: bool = True) -> 'str'`。
- **运行要求：** 需要可选包 hgvs（extra `dnakit[bio-hgvs]`）。只使用 hgvs 解析器，不需要数据库。
- **输入：** `variant` 是任何能被解析的 HGVS 变异字符串。`span` 是整数对 `(start, end)`，需满足 `0 <= start <= end`。`content` 是要打印的序列文本。`direction` 为 `""`、`">"` 或 `"<"`。
- **输出：** 一行文本。
- **示例代码：**

```python
from dnakit.variants import hgvs_sequence_line

line = hgvs_sequence_line("NM_000399.3:c.5A>G", (1, 10), "ACGTACGTA", direction=">", verbose=False)
print("|" + line.rstrip())
```

- **示例输出：**

```text
|   NM_000399.3  c          1 > ACGTACGTA                                     >       10 NM_000399.3:c.5A>G
```

## 10) `dnakit.variants.hgvs_variant_context_with_alignment()`

- **作用：** 把一个编码区变异放在基因组与转录本的比对上显示，并附带基因组和 cDNA 上下文。
- **处理方法：** 先检查 `uta_align` 包能否导入，再解析变异，然后调用 hgvs 的实验性函数 `variant_context_w_alignment(mapper, variant, margin=..., tx_ac=...)`。它从映射器的数据提供者取出基因组与转录本序列，用 `uta_align` 做全局比对，返回一幅多行文本图。该函数只适用于编码区变异。
- **API：** `hgvs_variant_context_with_alignment(mapper: 'Any', variant: 'str', *, margin: 'int' = 20, transcript: 'str | None' = None, verbose: bool = True) -> 'str'`。
- **运行要求：** 需要可选包 hgvs 和 uta-align（extra `dnakit[bio-hgvs-alignment]`），以及建立在 UTA 数据库连接和序列库上的 hgvs 映射器。缺少 `uta_align` 时会抛出带安装提示的 `BackendUnavailableError`，示例展示的正是在没有该包的机器上的这一行为。
- **输入：** `mapper`、`variant` 和 `transcript` 与 `hgvs_full_house()` 相同。`margin` 是两侧各增加的碱基数，不能为负。
- **输出：** 一个多行字符串。
- **示例代码：**

```python
from dnakit.exceptions import BackendUnavailableError
from dnakit.variants import hgvs_variant_context_with_alignment

try:
    hgvs_variant_context_with_alignment(None, "NM_000399.3:c.5A>G", verbose=False)
except BackendUnavailableError as error:
    print(type(error).__name__, error.code)
    print("uta-align" in error.hint)
```

- **示例输出：**

```text
BackendUnavailableError BACKEND_UNAVAILABLE
True
```

安装依赖后的用法：

```python
import hgvs.assemblymapper
import hgvs.dataproviders.uta

from dnakit.variants import hgvs_variant_context_with_alignment

hdp = hgvs.dataproviders.uta.connect()  # needs the UTA database and a sequence store
mapper = hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name="GRCh37")
print(hgvs_variant_context_with_alignment(mapper, "NM_000492.3:c.1521_1523del", margin=10))
```

## 11) `dnakit.annotation.polars_bio_describe_bgen()`

- **作用：** 列出本地 BGEN 文件产生的列，及其类型和 layout。
- **处理方法：** 调用 polars-bio 的 `describe_bgen(path, sample_path=..., bgi_path=...)`，把结果的每一行转成普通字典。
- **API：** `polars_bio_describe_bgen(path: 'str', *, sample_path: 'str | None' = None, bgi_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 polars-bio（extra `dnakit[bio-polarsbio]`）。
- **输入：** `path` 是以 `.bgen` 结尾的本地文件；URL（含 `://`）会被拒绝。`sample_path` 和 `bgi_path` 是可选的本地配套文件。
- **输出：** `status`、`operation`、`field_count`、`fields`（每个输出列一个字典，含 `name`、`type`、`layout`、`index` 和坐标系字段）以及 `provenance`。
- **示例代码：**

```python
from pathlib import Path

from dnakit.annotation import polars_bio_describe_bgen

Path("cohort.bgen").write_bytes(bytes.fromhex(Path("cohort.bgen.hex").read_text()))
result = polars_bio_describe_bgen("cohort.bgen", verbose=False)
print(result["field_count"])
print([field["name"] for field in result["fields"]])
print(result["fields"][0]["layout"], result["fields"][0]["coordinate_system_zero_based"])
```

- **示例输出：**

```text
9
['chrom', 'start', 'end', 'id', 'rsid', 'alleles', 'phased', 'bits', 'genotypes']
2 false
```

## 12) `dnakit.annotation.polars_bio_describe_pgen()`

- **作用：** 列出本地 PLINK 2 PGEN 文件集产生的列。
- **处理方法：** 调用 polars-bio 的 `describe_pgen(path, pvar_path=..., psam_path=..., pgi_path=...)`，把每一行转成普通字典。
- **API：** `polars_bio_describe_pgen(path: 'str', *, pvar_path: 'str | None' = None, psam_path: 'str | None' = None, pgi_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 polars-bio（extra `dnakit[bio-polarsbio]`）。
- **输入：** `path` 是以 `.pgen` 结尾的本地文件。`.pvar` 和 `.psam` 文件默认在它旁边查找，也可以用 `pvar_path`、`psam_path` 或 `pgi_path` 指定。
- **输出：** `status`、`operation`、`field_count`、`fields`（每个输出列一个字典，含 `name`、`type`、`storage_mode`、`index` 和坐标系字段）以及 `provenance`。
- **示例代码：**

```python
from pathlib import Path

from dnakit.annotation import polars_bio_describe_pgen

Path("cohort.pgen").write_bytes(bytes.fromhex(Path("cohort.pgen.hex").read_text()))
result = polars_bio_describe_pgen("cohort.pgen", verbose=False)
print(result["field_count"], result["fields"][0]["storage_mode"])
print([field["name"] for field in result["fields"]])
```

- **示例输出：**

```text
7 0x01
['chrom', 'start', 'end', 'id', 'ref', 'alt', 'genotypes']
```

## 13) `dnakit.annotation.polars_bio_describe_vcf_zarr()`

- **作用：** 列出本地 VCF Zarr 存储的 INFO 和 FORMAT 字段。
- **处理方法：** 调用 polars-bio 的 `describe_vcf_zarr(path)`，把每一行转成普通字典。
- **API：** `polars_bio_describe_vcf_zarr(path: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 polars-bio（extra `dnakit[bio-polarsbio]`）。构建示例存储还需要 `zarr`。
- **输入：** `path` 是已存在的本地目录（Zarr 存储）；URL 会被拒绝。
- **输出：** `status`、`operation`、`field_count`、`fields`（每项含 `name`、取值 `INFO` 或 `FORMAT` 的 `field_type`、`data_type` 和 `description`）以及 `provenance`。
- **示例代码：**

```python
from dnakit.annotation import polars_bio_describe_vcf_zarr
from make_vcz import make_vcz

make_vcz("cohort.vcz")
result = polars_bio_describe_vcf_zarr("cohort.vcz", verbose=False)
print(result["field_count"])
for field in result["fields"]:
    print(field["name"], field["field_type"], field["data_type"])
```

- **示例输出：**

```text
2
DP INFO Integer
genotypes FORMAT Struct
```

## 14) `dnakit.annotation.polars_bio_read_bgen()`

- **作用：** 把本地 BGEN 文件读成每个变异一条记录，含概率或剂量。
- **处理方法：** 调用 polars-bio 的 `read_bgen(path, genotype_output=..., samples=..., sample_path=..., bgi_path=...)`，把 polars 数据框转成字典，并拒绝返回超过 `max_rows` 的行数。
- **API：** `polars_bio_read_bgen(path: 'str', *, genotype_output: 'str' = 'probability', samples: 'list[str] | None' = None, sample_path: 'str | None' = None, bgi_path: 'str | None' = None, max_rows: 'int' = 100000, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 polars-bio（extra `dnakit[bio-polarsbio]`）。
- **输入：** `path` 是本地 `.bgen` 文件。`genotype_output` 为 `"probability"`（字段 `GP`）或 `"dosage"`（字段 `DS`）。`samples` 是可选的非空样本名列表，用于只保留这些样本。`max_rows` 是 1 到 1,000,000 的整数；变异数超过该值的文件会抛出 `ConfigurationError`。
- **输出：** `status`、`operation`、`genotype_output`、`variant_count`、`variants`（每项含 `chrom`、`start`、`end`、`id`、`rsid`、`alleles`、`phased`、`bits` 和 `genotypes` 映射）以及 `provenance`。坐标从 1 开始。
- **示例代码：**

```python
from pathlib import Path

from dnakit.annotation import polars_bio_read_bgen

Path("cohort.bgen").write_bytes(bytes.fromhex(Path("cohort.bgen.hex").read_text()))
result = polars_bio_read_bgen("cohort.bgen", verbose=False)
print(result["variant_count"], result["genotype_output"])
first = result["variants"][0]
print(first["rsid"], first["start"], first["alleles"], first["genotypes"]["GP"])
dosage = polars_bio_read_bgen("cohort.bgen", genotype_output="dosage", verbose=False)
print(dosage["variants"][1]["genotypes"]["DS"])
```

- **示例输出：**

```text
2 probability
rs1 100 ['A', 'G'] [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
[2.0, 0.0, 1.0]
```

## 15) `dnakit.annotation.polars_bio_read_bgen_matrix()`

- **作用：** 把本地 BGEN 文件的 ALT 剂量读成稠密的“变异乘样本”矩阵。
- **处理方法：** 调用 polars-bio 的 `read_bgen_matrix(path, samples=..., missing=..., ...)`，把矩阵、位置和样本名转成列表。矩阵大小限制为 5,000,000 个单元格。
- **API：** `polars_bio_read_bgen_matrix(path: 'str', *, samples: 'list[str] | None' = None, missing: 'float | None' = None, sample_path: 'str | None' = None, bgi_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 polars-bio（extra `dnakit[bio-polarsbio]`）。
- **输入：** `path` 是本地 `.bgen` 文件。`samples` 按名称选择并排序样本。`missing` 是用于缺失剂量的数值，`None` 表示使用 polars-bio 的默认值。
- **输出：** `status`、`operation`、`variant_count`、`sample_count`、`positions`、`sample_names`、`values`（行列表，每个变异一行）以及 `provenance`。
- **示例代码：**

```python
from pathlib import Path

from dnakit.annotation import polars_bio_read_bgen_matrix

Path("cohort.bgen").write_bytes(bytes.fromhex(Path("cohort.bgen.hex").read_text()))
result = polars_bio_read_bgen_matrix("cohort.bgen", verbose=False)
print(result["variant_count"], result["sample_count"], result["positions"])
print(result["sample_names"], result["values"])
subset = polars_bio_read_bgen_matrix("cohort.bgen", samples=["S3", "S1"], verbose=False)
print(subset["sample_names"], subset["values"])
```

- **示例输出：**

```text
2 3 [100, 250]
['S1', 'S2', 'S3'] [[0.0, 1.0, 2.0], [2.0, 0.0, 1.0]]
['S3', 'S1'] [[2.0, 0.0], [1.0, 2.0]]
```

## 16) `dnakit.annotation.polars_bio_read_pgen()`

- **作用：** 把本地 PLINK 2 PGEN 文件集读成每个变异一条记录。
- **处理方法：** 调用 polars-bio 的 `read_pgen(path, genotype_fields=..., samples=..., psam_id_mode=..., ...)`，把数据框转成字典，并拒绝超过 `max_rows` 的行数。字段列表默认为 `["GT"]`。
- **API：** `polars_bio_read_pgen(path: 'str', *, genotype_fields: 'list[str] | None' = None, samples: 'list[str] | None' = None, psam_id_mode: 'str' = 'iid', pvar_path: 'str | None' = None, psam_path: 'str | None' = None, pgi_path: 'str | None' = None, max_rows: 'int' = 100000, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 polars-bio（extra `dnakit[bio-polarsbio]`）。
- **输入：** `path` 是本地 `.pgen` 文件，`.pvar` 和 `.psam` 放在它旁边。`genotype_fields` 是非空的字段名列表（默认 `GT`）。`psam_id_mode` 为 `"iid"`、`"fid_iid"` 或 `"fid_iid_sid"`。`samples` 保留指定的样本。`max_rows` 是 1 到 1,000,000 的整数。
- **输出：** `status`、`operation`、`variant_count`、`variants`（每项含 `chrom`、`start`、`end`、`id`、`ref`、`alt` 和 `genotypes` 映射，缺失的基因型为 `None`）以及 `provenance`。坐标从 1 开始。
- **示例代码：**

```python
from pathlib import Path

from dnakit.annotation import polars_bio_read_pgen

Path("cohort.pgen").write_bytes(bytes.fromhex(Path("cohort.pgen.hex").read_text()))
result = polars_bio_read_pgen("cohort.pgen", verbose=False)
print(result["variant_count"])
for variant in result["variants"]:
    print(variant["id"], variant["start"], variant["ref"], variant["alt"], variant["genotypes"]["GT"])
```

- **示例输出：**

```text
2
rs1 100 A ['G'] [[0, 0], [0, 1], [1, 1]]
rs2 250 C ['T'] [[0, 1], [0, 0], None]
```

## 17) `dnakit.annotation.polars_bio_read_pgen_matrix()`

- **作用：** 把 PGEN 文件集的某一个基因型字段读成稠密的“变异乘样本”矩阵。
- **处理方法：** 调用 polars-bio 的 `read_pgen_matrix(path, field=..., samples=..., missing=..., psam_id_mode=..., ...)`，把矩阵、位置和样本名转成列表。矩阵大小限制为 5,000,000 个单元格。
- **API：** `polars_bio_read_pgen_matrix(path: 'str', *, field: 'str' = 'ALT_COUNT', samples: 'list[str] | None' = None, missing: 'float | None' = None, psam_id_mode: 'str' = 'iid', pvar_path: 'str | None' = None, psam_path: 'str | None' = None, pgi_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 polars-bio（extra `dnakit[bio-polarsbio]`）。
- **输入：** `path` 是本地 `.pgen` 文件。`field` 为 `"ALT_COUNT"` 或 `"DS"`。`missing` 是用于缺失基因型的数值，`None` 表示使用 polars-bio 的默认值（`ALT_COUNT` 为 -9）。其余参数与 `polars_bio_read_pgen()` 相同。
- **输出：** `status`、`operation`、`variant_count`、`sample_count`、`positions`、`sample_names`、`values` 以及 `provenance`。
- **示例代码：**

```python
from pathlib import Path

from dnakit.annotation import polars_bio_read_pgen_matrix

Path("cohort.pgen").write_bytes(bytes.fromhex(Path("cohort.pgen.hex").read_text()))
default = polars_bio_read_pgen_matrix("cohort.pgen", verbose=False)
print(default["sample_names"], default["values"])
marked = polars_bio_read_pgen_matrix("cohort.pgen", missing=-1, verbose=False)
print(marked["values"])
```

- **示例输出：**

```text
['S1', 'S2', 'S3'] [[0, 1, 2], [1, 0, -9]]
[[0, 1, 2], [1, 0, -1]]
```

## 18) `dnakit.annotation.polars_bio_read_vcf_zarr()`

- **作用：** 把本地 VCF Zarr 存储读成每个变异一条记录。
- **处理方法：** 调用 polars-bio 的 `read_vcf_zarr(path, info_fields=..., format_fields=..., samples=...)`，把数据框转成字典，并拒绝超过 `max_rows` 的行数。
- **API：** `polars_bio_read_vcf_zarr(path: 'str', *, info_fields: 'list[str] | None' = None, format_fields: 'list[str] | None' = None, samples: 'list[str] | None' = None, max_rows: 'int' = 100000, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 polars-bio（extra `dnakit[bio-polarsbio]`）。构建示例存储还需要 `zarr`。
- **输入：** `path` 是已存在的本地目录。`info_fields`、`format_fields` 和 `samples` 是可选的非空名称列表。`max_rows` 是 1 到 1,000,000 的整数。
- **输出：** `status`、`operation`、`variant_count`、`variants`（每项含 `chrom`、`start`、`end`、`id`、`ref`、`alt`、`qual`、`filter`、所请求的 INFO 字段和 `genotypes` 映射）以及 `provenance`。坐标从 1 开始。
- **示例代码：**

```python
from dnakit.annotation import polars_bio_read_vcf_zarr
from make_vcz import make_vcz

make_vcz("cohort.vcz")
result = polars_bio_read_vcf_zarr("cohort.vcz", info_fields=["DP"], verbose=False)
print(result["variant_count"])
for variant in result["variants"]:
    print(variant["id"], variant["start"], variant["ref"], variant["alt"], variant["DP"], variant["genotypes"]["GT"])
```

- **示例输出：**

```text
2
rs1 100 A G 14 [[0, 0], [0, 1]]
rs2 250 C T 20 [[1, 1], [0, 0]]
```

## 19) `dnakit.annotation.polars_bio_register_bgen()`

- **作用：** 把本地 BGEN 文件加载为临时 SQL 表，并返回它的全部行。
- **处理方法：** 用 `register_bgen(path, table_name, ...)` 把文件注册到 polars-bio 会话中，执行 `SELECT * FROM <table> LIMIT max_rows+1`，转换各行，并在最后无论如何都注销该表。会话中已有同名表时抛出 `ConfigurationError`。
- **API：** `polars_bio_register_bgen(path: 'str', table_name: 'str', *, genotype_output: 'str' = 'probability', samples: 'list[str] | None' = None, sample_path: 'str | None' = None, bgi_path: 'str | None' = None, max_rows: 'int' = 100000, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 polars-bio（extra `dnakit[bio-polarsbio]`）。
- **输入：** `path` 是本地 `.bgen` 文件。`table_name` 由小写字母、数字和下划线组成，以字母或下划线开头，最长 63 个字符。其余参数与 `polars_bio_read_bgen()` 相同。
- **输出：** `status`、`operation`、`table_name`、`columns`、`row_count`、`rows` 以及 `provenance`。调用结束后该表不再存在。
- **示例代码：**

```python
from pathlib import Path

from dnakit.annotation import polars_bio_register_bgen

Path("cohort.bgen").write_bytes(bytes.fromhex(Path("cohort.bgen.hex").read_text()))
result = polars_bio_register_bgen("cohort.bgen", "cohort_bgen", verbose=False)
print(result["table_name"], result["row_count"])
print(result["columns"])
print([row["rsid"] for row in result["rows"]])
```

- **示例输出：**

```text
cohort_bgen 2
['chrom', 'start', 'end', 'id', 'rsid', 'alleles', 'phased', 'bits', 'genotypes']
['rs1', 'rs2']
```

## 20) `dnakit.annotation.polars_bio_register_cram()`

- **作用：** 把本地 CRAM 文件加载为临时 SQL 表，并返回它的全部行。
- **处理方法：** 用 polars-bio 的 `register_cram(path, table_name, tag_fields=...)` 注册文件，执行 `SELECT * FROM <table> LIMIT max_rows+1`，转换各行并注销该表。这里 polars-bio 不接受参考 FASTA，所以只有内嵌序列的 CRAM 文件才能解码。
- **API：** `polars_bio_register_cram(path: 'str', table_name: 'str', *, tag_fields: 'list[str] | None' = None, max_rows: 'int' = 100000, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 polars-bio（extra `dnakit[bio-polarsbio]`）。示例用 `pysam` 写出 CRAM 文件。
- **输入：** `path` 是本地 CRAM 文件。`table_name` 的规则与 `polars_bio_register_bgen()` 相同。`tag_fields` 是可选的非空比对标签列表，这些标签会作为列加入。`max_rows` 是 1 到 1,000,000 的整数。
- **输出：** `status`、`operation`、`table_name`、`columns`、`row_count`、`rows` 以及 `provenance`。
- **示例代码：**

```python
import pysam

from dnakit.annotation import polars_bio_register_cram

header = {"HD": {"VN": "1.6", "SO": "unsorted"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
with pysam.AlignmentFile("unmapped.cram", "wc", header=header) as out:
    for index, sequence in enumerate(["ACGTACGTAC", "GGGGCCCCAA"]):
        read = pysam.AlignedSegment()
        read.query_name = f"r{index}"
        read.query_sequence = sequence
        read.query_qualities = pysam.qualitystring_to_array("I" * len(sequence))
        read.flag = 4
        read.reference_id = -1
        read.reference_start = -1
        read.set_tag("NM", 0)
        out.write(read)

result = polars_bio_register_cram("unmapped.cram", "reads_cram", tag_fields=["NM"], verbose=False)
print(result["table_name"], result["row_count"])
print(result["columns"])
print([(row["name"], row["sequence"], row["NM"]) for row in result["rows"]])
```

- **示例输出：**

```text
reads_cram 2
['name', 'chrom', 'start', 'end', 'flags', 'cigar', 'mapping_quality', 'mate_chrom', 'mate_start', 'sequence', 'quality_scores', 'template_length', 'NM']
[('r0', 'ACGTACGTAC', 0), ('r1', 'GGGGCCCCAA', 0)]
```

## 21) `dnakit.annotation.polars_bio_register_pgen()`

- **作用：** 把本地 PGEN 文件集加载为临时 SQL 表，并返回它的全部行。
- **处理方法：** 用 polars-bio 的 `register_pgen(path, table_name, ...)` 注册文件集，执行 `SELECT * FROM <table> LIMIT max_rows+1`，转换各行并注销该表。
- **API：** `polars_bio_register_pgen(path: 'str', table_name: 'str', *, genotype_fields: 'list[str] | None' = None, samples: 'list[str] | None' = None, psam_id_mode: 'str' = 'iid', pvar_path: 'str | None' = None, psam_path: 'str | None' = None, pgi_path: 'str | None' = None, max_rows: 'int' = 100000, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 polars-bio（extra `dnakit[bio-polarsbio]`）。
- **输入：** `path` 是本地 `.pgen` 文件；`table_name` 的规则与 `polars_bio_register_bgen()` 相同；其余参数与 `polars_bio_read_pgen()` 相同。
- **输出：** `status`、`operation`、`table_name`、`columns`、`row_count`、`rows` 以及 `provenance`。
- **示例代码：**

```python
from pathlib import Path

from dnakit.annotation import polars_bio_register_pgen

Path("cohort.pgen").write_bytes(bytes.fromhex(Path("cohort.pgen.hex").read_text()))
result = polars_bio_register_pgen("cohort.pgen", "cohort_pgen", verbose=False)
print(result["table_name"], result["row_count"])
print(result["columns"])
print([row["genotypes"]["GT"] for row in result["rows"]])
```

- **示例输出：**

```text
cohort_pgen 2
['chrom', 'start', 'end', 'id', 'ref', 'alt', 'genotypes']
[[[0, 0], [0, 1], [1, 1]], [[0, 1], [0, 0], None]]
```

## 22) `dnakit.annotation.polars_bio_register_vcf_zarr()`

- **作用：** 把本地 VCF Zarr 存储加载为临时 SQL 表，并返回它的全部行。
- **处理方法：** 用 polars-bio 的 `register_vcf_zarr(path, table_name, ...)` 注册存储，执行 `SELECT * FROM <table> LIMIT max_rows+1`，转换各行并注销该表。
- **API：** `polars_bio_register_vcf_zarr(path: 'str', table_name: 'str', *, info_fields: 'list[str] | None' = None, format_fields: 'list[str] | None' = None, samples: 'list[str] | None' = None, max_rows: 'int' = 100000, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 polars-bio（extra `dnakit[bio-polarsbio]`）。构建示例存储还需要 `zarr`。
- **输入：** `path` 是已存在的本地目录；`table_name` 的规则与 `polars_bio_register_bgen()` 相同；其余参数与 `polars_bio_read_vcf_zarr()` 相同。
- **输出：** `status`、`operation`、`table_name`、`columns`、`row_count`、`rows` 以及 `provenance`。
- **示例代码：**

```python
from dnakit.annotation import polars_bio_register_vcf_zarr
from make_vcz import make_vcz

make_vcz("cohort.vcz")
result = polars_bio_register_vcf_zarr("cohort.vcz", "cohort_vcz", info_fields=["DP"], verbose=False)
print(result["table_name"], result["row_count"])
print(result["columns"])
print([(row["id"], row["DP"]) for row in result["rows"]])
```

- **示例输出：**

```text
cohort_vcz 2
['chrom', 'start', 'end', 'id', 'ref', 'alt', 'qual', 'filter', 'DP', 'genotypes']
[('rs1', 14), ('rs2', 20)]
```

## 23) `dnakit.annotation.polars_bio_scan_bgen()`

- **作用：** 惰性扫描本地 BGEN 文件，统计全部变异数，并最多返回 `max_rows` 个变异。
- **处理方法：** 用 polars-bio 的 `scan_bgen(path, ...)` 打开文件，用惰性的 `pl.len()` 计数，然后只收集前 `max_rows` 行。结果会报告列结构以及输出是否被截断。
- **API：** `polars_bio_scan_bgen(path: 'str', *, genotype_output: 'str' = 'probability', samples: 'list[str] | None' = None, sample_path: 'str | None' = None, bgi_path: 'str | None' = None, max_rows: 'int' = 100000, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 polars-bio（extra `dnakit[bio-polarsbio]`）。
- **输入：** 与 `polars_bio_read_bgen()` 相同。与读取函数不同的是，变异数超过 `max_rows` 的文件会被截断而不是被拒绝。
- **输出：** `status`、`operation`、`columns`（名称与 dtype 的配对）、`variant_count`（全部变异数）、`returned_count`、`truncated`、`variants` 以及 `provenance`。
- **示例代码：**

```python
from pathlib import Path

from dnakit.annotation import polars_bio_scan_bgen

Path("cohort.bgen").write_bytes(bytes.fromhex(Path("cohort.bgen.hex").read_text()))
result = polars_bio_scan_bgen("cohort.bgen", max_rows=1, verbose=False)
print(result["variant_count"], result["returned_count"], result["truncated"])
print([column["name"] for column in result["columns"]])
print(result["variants"][0]["rsid"])
```

- **示例输出：**

```text
2 1 True
['chrom', 'start', 'end', 'id', 'rsid', 'alleles', 'phased', 'bits', 'genotypes']
rs1
```

## 24) `dnakit.annotation.polars_bio_scan_pgen()`

- **作用：** 惰性扫描本地 PGEN 文件集，统计全部变异数，并最多返回 `max_rows` 个变异。
- **处理方法：** 用 polars-bio 的 `scan_pgen(path, ...)` 打开文件集，惰性计数，然后只收集前 `max_rows` 行。
- **API：** `polars_bio_scan_pgen(path: 'str', *, genotype_fields: 'list[str] | None' = None, samples: 'list[str] | None' = None, psam_id_mode: 'str' = 'iid', pvar_path: 'str | None' = None, psam_path: 'str | None' = None, pgi_path: 'str | None' = None, max_rows: 'int' = 100000, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 polars-bio（extra `dnakit[bio-polarsbio]`）。
- **输入：** 与 `polars_bio_read_pgen()` 相同。
- **输出：** `status`、`operation`、`columns`、`variant_count`、`returned_count`、`truncated`、`variants` 以及 `provenance`。
- **示例代码：**

```python
from pathlib import Path

from dnakit.annotation import polars_bio_scan_pgen

Path("cohort.pgen").write_bytes(bytes.fromhex(Path("cohort.pgen.hex").read_text()))
result = polars_bio_scan_pgen("cohort.pgen", max_rows=1, verbose=False)
print(result["variant_count"], result["returned_count"], result["truncated"])
print([column["name"] for column in result["columns"]])
print(result["variants"][0]["id"])
```

- **示例输出：**

```text
2 1 True
['chrom', 'start', 'end', 'id', 'ref', 'alt', 'genotypes']
rs1
```

## 25) `dnakit.annotation.polars_bio_scan_vcf_zarr()`

- **作用：** 惰性扫描本地 VCF Zarr 存储，统计全部变异数，并最多返回 `max_rows` 个变异。
- **处理方法：** 用 polars-bio 的 `scan_vcf_zarr(path, ...)` 打开存储，惰性计数，然后只收集前 `max_rows` 行。
- **API：** `polars_bio_scan_vcf_zarr(path: 'str', *, info_fields: 'list[str] | None' = None, format_fields: 'list[str] | None' = None, samples: 'list[str] | None' = None, max_rows: 'int' = 100000, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 polars-bio（extra `dnakit[bio-polarsbio]`）。构建示例存储还需要 `zarr`。
- **输入：** 与 `polars_bio_read_vcf_zarr()` 相同。
- **输出：** `status`、`operation`、`columns`、`variant_count`、`returned_count`、`truncated`、`variants` 以及 `provenance`。
- **示例代码：**

```python
from dnakit.annotation import polars_bio_scan_vcf_zarr
from make_vcz import make_vcz

make_vcz("cohort.vcz")
result = polars_bio_scan_vcf_zarr("cohort.vcz", info_fields=["DP"], max_rows=1, verbose=False)
print(result["variant_count"], result["returned_count"], result["truncated"])
print([column["name"] for column in result["columns"]])
print(result["variants"][0]["id"], result["variants"][0]["DP"])
```

- **示例输出：**

```text
2 1 True
['chrom', 'start', 'end', 'id', 'ref', 'alt', 'qual', 'filter', 'DP', 'genotypes']
rs1 14
```

## 26) `dnakit.annotation.polars_bio_visualize_intervals()`

- **作用：** 把两个 BED 文件中相互重叠的区间对画进一张 PNG，每对一个子图。
- **处理方法：** 用 polars-bio 的 `read_bed(..., use_zero_based=True)` 读取两个文件，计算 `overlap(...)`，对前 `max_pairs` 个区间对逐个调用 polars-bio 的 `visualize_intervals`，它会自己打开一个 matplotlib 图。这些图被截取并纵向拼成一张图，写入 `output_path`。两个文件没有重叠时抛出 `ConfigurationError`。
- **API：** `polars_bio_visualize_intervals(left_bed: 'str', right_bed: 'str', output_path: 'str', *, label: 'str' = 'overlapping pair', max_pairs: 'int' = 20, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 polars-bio、matplotlib（extra `dnakit[bio-plots]`）和 bioframe（由 `polars-bio[viz]` 安装；缺少 bioframe 时抛出 `BackendUnavailableError`）。
- **输入：** `left_bed` 和 `right_bed` 是本地 BED 文件。`output_path` 是新的 `.png` 路径（已有文件不会被覆盖，缺失的上级目录会被创建）。`label` 是非空且不超过 128 个字符的字符串。`max_pairs` 是 1 到 50 的整数。
- **输出：** `path`（PNG 文件）和 `data`，后者含 `label`、`pair_count`（全部重叠数）、`plotted_count`、`pairs`（已画出的区间对，BED 列带 `_1` 和 `_2` 后缀）以及 `provenance`。区间是 0 起始的左闭右开区间。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.annotation import polars_bio_visualize_intervals

target = Path(tempfile.mkdtemp()) / "pairs.png"
result = polars_bio_visualize_intervals("a.bed", "b.bed", str(target), verbose=False)
print(result["data"]["pair_count"], result["data"]["plotted_count"])
for pair in result["data"]["pairs"]:
    print(pair["name_1"], pair["start_1"], pair["end_1"], pair["name_2"], pair["start_2"], pair["end_2"])
print(target.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n")
```

- **示例输出：**

```text
2 2
a 10 20 c 15 30
b 40 60 d 50 55
True
```

## 27) `dnakit.annotation.cached_ensembl_release()`

- **作用：** 描述已安装在本地 pyensembl 缓存中的某个 Ensembl 发布版。
- **处理方法：** 调用 `pyensembl.ensembl_release.cached_release(release, species=species)`，再调用只读缓存的 `inspect_data()`。发布版已安装时，结果列出来自注释索引的基因数和 contig 数以及缓存文件；未安装时抛出 `BackendUnavailableError`，其提示是 `pyensembl install` 命令。不会下载任何内容。
- **API：** `cached_ensembl_release(release: 'int | str', *, species: 'str' = 'human', verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 pyensembl（extra `dnakit[bio-pyensembl]`）和一个已安装的发布版（例如通过 `pyensembl install --release 110 --species human` 或 `pyensembl_install_ensembl_release()` 安装）。缓存目录是 `$PYENSEMBL_CACHE_DIR` 或 `~/.cache/pyensembl`。示例把缓存指向一个空目录并展示未安装时的错误；用法代码块展示在已安装发布版上的调用。
- **输入：** `release` 是正整数，或形如 `2026_04` 的注释日期字符串。`species` 是 pyensembl 认识的非空物种名；未知名称会抛出 `ConfigurationError`。
- **输出：** `status`、`operation`、`service`、`backend`、`version`、`request`、`retrieved_utc`，以及 `data`，其中含 `species`、`reference_name`、`release`、`annotation_name`、`annotation_version`、`cache_directory`、`installed`、`annotation_status`、`reference_dna_status`、`gene_count`、`contig_count` 和 `files`。
- **示例代码：**

```python
import os
from pathlib import Path

os.environ["PYENSEMBL_CACHE_DIR"] = str(Path("empty_cache").resolve())

from dnakit.annotation import cached_ensembl_release
from dnakit.exceptions import BackendUnavailableError

try:
    cached_ensembl_release(110, verbose=False)
except BackendUnavailableError as error:
    print(type(error).__name__, error.code)
    print(error.hint.split("; ")[0])
```

- **示例输出：**

```text
BackendUnavailableError BACKEND_UNAVAILABLE
Install the annotation with the pyensembl command line first, for example: pyensembl install --release 110 --species homo_sapiens
```

安装依赖后的用法：

```python
from dnakit.annotation import cached_ensembl_release

result = cached_ensembl_release(110)  # after: pyensembl install --release 110 --species human
print(result["data"]["gene_count"], result["data"]["contig_count"])
```

## 28) `dnakit.annotation.ensembl_release_for_reference_name()`

- **作用：** 描述 pyensembl 为 GRCh38 这类参考基因组名称所选中的、已安装的 Ensembl 发布版。
- **处理方法：** 调用 `pyensembl.reference_name.genome_for_reference_name(name, allow_older_downloaded_release=...)`。标志为真时选择已安装的最新发布版，否则选择已下载文件的最新发布版。随后按 `cached_ensembl_release()` 的方式描述所选发布版；未安装时抛出 `BackendUnavailableError`。不会下载任何内容。
- **API：** `ensembl_release_for_reference_name(reference_name: 'str', *, allow_older_downloaded_release: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 pyensembl（extra `dnakit[bio-pyensembl]`）和该参考基因组已安装的发布版。示例展示对未知参考名称的拒绝；用法代码块展示已安装的情形。
- **输入：** `reference_name` 是非空字符串，例如 `GRCh38`；未知名称会抛出 `ConfigurationError`。`allow_older_downloaded_release` 必须是布尔值。
- **输出：** 与 `cached_ensembl_release()` 相同的字典，其中 `request` 含 `reference_name` 和 `allow_older_downloaded_release`。
- **示例代码：**

```python
import os
from pathlib import Path

os.environ["PYENSEMBL_CACHE_DIR"] = str(Path("empty_cache").resolve())

from dnakit.annotation import ensembl_release_for_reference_name
from dnakit.exceptions import ConfigurationError

try:
    ensembl_release_for_reference_name("GRCh99", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error.code)
    print(error.message)
```

- **示例输出：**

```text
ConfigurationError CONFIGURATION_ERROR
pyensembl genome_for_reference_name rejected the input: Reference genome 'GRCh99' not found
```

安装依赖后的用法：

```python
from dnakit.annotation import ensembl_release_for_reference_name

result = ensembl_release_for_reference_name("GRCh38")  # needs an installed GRCh38 release
print(result["data"]["release"], result["data"]["gene_count"])
```

## 29) `dnakit.annotation.fetch_all_ensembl_dated_releases()`

- **作用：** 收集所有带日期发布版的物种的注释日期。
- **处理方法：** 调用 `pyensembl.dated_releases.fetch_all_dated_releases(expire_after=..., timeout=...)`。比 `expire_after_days` 新的缓存列表原样使用；先获取人类的列表，若读取失败，则所有物种都回退到各自缓存的日期。日期未知的物种映射为 `None`。如果什么都取不到、也没有缓存，则抛出 `BackendExecutionError`。
- **API：** `fetch_all_ensembl_dated_releases(*, expire_after_days: 'float' = 1.0, timeout_seconds: 'int' = 10, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 pyensembl（extra `dnakit[bio-pyensembl]`），对缺失或过期的列表还需要联网。示例展示参数检查，它发生在任何网络访问之前。
- **输入：** `expire_after_days` 是大于等于零的有限数（负值会抛出 `ConfigurationError`）。`timeout_seconds` 是正整数。
- **输出：** 标准的列表信封（`status`、`operation`、`service`、`backend`、`version`、`request`、`retrieved_utc`）以及 `data`，其中含 `dates_by_species`、`species_count`、`species_with_dates` 和 `species_without_dates`。
- **示例代码：**

```python
from dnakit.annotation import fetch_all_ensembl_dated_releases
from dnakit.exceptions import ConfigurationError

try:
    fetch_all_ensembl_dated_releases(expire_after_days=-1, verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error.code)
    print(error.message)
```

- **示例输出：**

```text
ConfigurationError CONFIGURATION_ERROR
expire_after_days must be a finite number of zero or more.
```

安装依赖后的用法：

```python
from dnakit.annotation import fetch_all_ensembl_dated_releases

result = fetch_all_ensembl_dated_releases(expire_after_days=7)
print(result["data"]["species_count"], result["data"]["dates_by_species"]["mus_musculus"])
```

## 30) `dnakit.annotation.list_ensembl_annotation_dates()`

- **作用：** 列出 Ensembl 为某个基因组组装登录号发布的注释日期。
- **处理方法：** 校验登录号和提供者，然后调用 `pyensembl.dated_releases.available_annotation_dates(accession, provider, refresh=...)`。列表第一次从 Ensembl 读取，之后从缓存复用，除非 `refresh` 为真。日期按从旧到新的顺序以 `YYYY_MM` 形式返回。
- **API：** `list_ensembl_annotation_dates(assembly_accession: 'str', *, provider: 'str' = 'ensembl', refresh: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 pyensembl（extra `dnakit[bio-pyensembl]`），对未缓存的列表还需要联网。示例展示登录号检查，它发生在任何网络访问之前。
- **输入：** `assembly_accession` 是带版本号的 GenBank 或 RefSeq 登录号，格式为 `GCA_` 或 `GCF_` 加九位数字、一个点和版本号，例如 `GCA_000001405.29`。`provider` 是小写的目录名（字母和下划线，最长 32 个字符）。`refresh` 是布尔值。
- **输出：** 标准的列表信封以及 `data`，其中含 `dates`（字符串）和 `count`。读不到的列表会抛出 `BackendExecutionError`；非法登录号会抛出 `ConfigurationError`。
- **示例代码：**

```python
from dnakit.annotation import list_ensembl_annotation_dates
from dnakit.exceptions import ConfigurationError

try:
    list_ensembl_annotation_dates("GRCh38", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error.code)
    print(error.message)
```

- **示例输出：**

```text
ConfigurationError CONFIGURATION_ERROR
assembly_accession must be a versioned GCA or GCF accession, for example GCA_000001405.29.
```

安装依赖后的用法：

```python
from dnakit.annotation import list_ensembl_annotation_dates

result = list_ensembl_annotation_dates("GCA_000001405.29")
print(result["data"]["count"], result["data"]["dates"][:3])
```

## 31) `dnakit.annotation.list_ensembl_dated_releases()`

- **作用：** 列出某个物种当前组装版本已发布的注释日期。
- **处理方法：** 用 pyensembl 的物种注册表（`check_species_object`）解析物种；没有带日期发布版的物种会在任何网络访问之前被拒绝。否则调用 `pyensembl.dated_releases.available_dated_releases(species, refresh=...)`，它读取该物种组装登录号和提供者对应的 Ensembl 列表，或复用缓存的列表。
- **API：** `list_ensembl_dated_releases(species: 'str' = 'human', *, refresh: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 pyensembl（extra `dnakit[bio-pyensembl]`），除非存在缓存列表且 `refresh` 为假，否则需要联网。示例展示对没有带日期发布版的物种的拒绝，这一步无需联网。
- **输入：** `species` 是 pyensembl 认识的物种名（未知名称会抛出 `ConfigurationError`）。`refresh` 是布尔值。
- **输出：** 标准的列表信封以及 `data`，其中含 `latin_name`、`assembly_accession`、`provider`、`dates` 和 `count`。
- **示例代码：**

```python
from dnakit.annotation import list_ensembl_dated_releases
from dnakit.exceptions import ConfigurationError

try:
    list_ensembl_dated_releases("toxoplasma_gondii", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error.code)
    print(error.message)
```

- **示例输出：**

```text
ConfigurationError CONFIGURATION_ERROR
toxoplasma_gondii has no dated Ensembl releases.
```

安装依赖后的用法：

```python
from dnakit.annotation import list_ensembl_dated_releases

result = list_ensembl_dated_releases("mouse")
print(result["data"]["assembly_accession"], result["data"]["dates"])
```

## 32) `dnakit.annotation.pyensembl_install_local_annotation()`

- **作用：** 用本地 GTF 文件（及可选的本地 FASTA 文件）构建 pyensembl 索引，使其可以像 Ensembl 发布版一样被查询。
- **处理方法：** 用本地路径创建 `pyensembl.Genome`，然后调用它的 `download()`（会让本地文件保持原位）和 `index()`，与带 `--gtf` 的 `pyensembl install` 动作一致。全程不使用网络。结果统计已建立索引的基因数和转录本数。
- **API：** `pyensembl_install_local_annotation(reference_name: 'str', annotation_name: 'str', gtf_path: 'str | Path', *, annotation_version: 'int | None' = None, transcript_fasta_paths: 'Sequence[str | Path]' = (), protein_fasta_paths: 'Sequence[str | Path]' = (), genome_fasta_path: 'str | Path | None' = None, cache_directory_path: 'str | Path | None' = None, overwrite: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 pyensembl（extra `dnakit[bio-pyensembl]`）。索引会写入缓存目录。
- **输入：** `reference_name` 和 `annotation_name` 是 pyensembl 保存该基因组所用的名称（最长 128 个字符）。`gtf_path` 是已存在的本地文件。`annotation_version` 是正整数或 `None`。FASTA 参数是已存在本地文件的列表（最多 10 个文件），`genome_fasta_path` 则是单个文件。`cache_directory_path` 必须是已存在的目录，否则使用 pyensembl 的默认缓存。`overwrite=True` 会重建已有索引。
- **输出：** `status`（`installed` 或 `not_installed`）、`operation`、`reference_name`、`annotation_name`、`annotation_version`、`gtf`（文件名）、`cache_directory`、`genes`、`transcripts` 以及 `provenance`。
- **示例代码：**

```python
from pathlib import Path

from dnakit.annotation import pyensembl_install_local_annotation

Path("toy_cache").mkdir(exist_ok=True)
result = pyensembl_install_local_annotation(
    "toyref", "toyann", "mini.gtf", cache_directory_path="toy_cache", verbose=False
)
print(result["status"], result["gtf"], result["genes"], result["transcripts"])
print(sorted(path.name for path in Path("toy_cache").iterdir()))
```

- **示例输出：**

```text
installed mini.gtf 2 2
['mini.db']
```

## 33) `dnakit.annotation.pyensembl_install_ensembl_release()`

- **作用：** 下载某个物种的一个 Ensembl 发布版，并在 pyensembl 缓存中建立索引。
- **处理方法：** 创建 `pyensembl.EnsemblRelease(release=..., species=...)`，调用 `download()` 和 `index()`，与 `pyensembl install --release` 一致。pyensembl 从 Ensembl FTP 服务器获取 GTF 注释以及转录本和蛋白质 FASTA 文件，并建立基因索引。结果记录发布版、物种、GTF 的 URL，以及基因数和转录本数。
- **API：** `pyensembl_install_ensembl_release(release: 'int | str', species: 'str' = 'human', *, overwrite: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 pyensembl（extra `dnakit[bio-pyensembl]`）和访问 Ensembl FTP 服务器的网络；发布版注释和基因组 FASTA 的下载量较大。示例展示在任何下载之前执行的参数检查。
- **输入：** `release` 是 1 到 1000 的整数，或形如 `2026_04` 的字符串。`species` 是 pyensembl 的物种名，例如 `human` 或 `caenorhabditis_elegans`。`overwrite` 是布尔值。pyensembl 拒绝的物种或发布版会抛出 `ConfigurationError`。
- **输出：** `status`（`installed` 或 `not_installed`）、`operation`、`ensembl_release`、`species`、`gtf_url`、`cache_directory`、`genes`、`transcripts` 以及 `provenance`。
- **示例代码：**

```python
from dnakit.annotation import pyensembl_install_ensembl_release
from dnakit.exceptions import ConfigurationError

try:
    pyensembl_install_ensembl_release(0, verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error.code)
    print(error.message)
```

- **示例输出：**

```text
ConfigurationError CONFIGURATION_ERROR
release must be between 1 and 1000.
```

安装依赖后的用法：

```python
from dnakit.annotation import pyensembl_install_ensembl_release

result = pyensembl_install_ensembl_release(110, species="caenorhabditis_elegans")
print(result["status"], result["genes"], result["transcripts"])
```
