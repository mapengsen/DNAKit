# Biopython：NCBI、KEGG 与 TogoWS 服务、系统发育树与绘图

这一组函数是 `dnakit.sequence` 中对 [Biopython](https://biopython.org/)（包名 `Bio`，[源码](https://github.com/biopython/biopython)）的封装。其中 19 个函数查询远程服务，并把返回内容整理成便于 JSON 序列化的字典：通过 NCBI QBLAST 做 BLAST 检索、NCBI Entrez E-utilities（`einfo`、`esearch`、`efetch`、`esummary`、`elink`、`epost`）、KEGG REST 服务和 TogoWS 服务。另外 3 个函数在本地运行：用 `Bio.Phylo` 把 Newick 树画成 ASCII 字符图或 PNG，以及用 `Bio.SeqUtils` 计算并绘制 DNA 序列的 GC 偏斜。本模块只是包装这些上游函数，不镜像任何数据库，也不解释返回内容的生物学含义。

每个联网函数先检查参数，调用时才会调用一个上游函数，在此之前不会向远程服务发送任何内容。数量有上限（每次请求最多 200 个 NCBI 标识符，`kegg_get_entries()` 最多 10 个 KEGG 条目，每次转换或链接请求最多 100 个 KEGG 标识符，每次 TogoWS 迭代最多 1000 个标识符）。NCBI 邮箱地址是可选的：它只在调用期间设置到 Biopython 模块上，不会写入结果。输出文件必须是新路径，已存在的文件会被拒绝。远程服务产生的错误（未知标识符、HTTP 错误、NCBI 错误回复）会原样向上抛出。发送请求前请查看 NCBI、KEGG 和 TogoWS 的使用条款。

Biopython 使用 Biopython License Agreement 许可（部分文件同时以 BSD 3-Clause 许可双重授权），已核对已安装包的元数据及其 `LICENSE.rst`。这是宽松许可，不是 GPL。依赖都是可选的，只有调用对应函数时才会导入；缺少包时会报错并给出安装命令。示例从 `dnakit.sequence` 导入。

## 示例数据

不需要输入文件。三个离线示例直接传入输入，并把文件写到临时目录。联网函数的示例用会被参数校验拒绝的参数调用函数，因此无需联网即可运行并打印真实的错误；这些小节还给出一段不执行的、联网时的实际调用方式。

## 1) `dnakit.sequence.blast_qblast_search()`

- **作用：** 在 NCBI 服务器上运行 BLAST 检索并汇总命中结果。
- **处理方法：** 先校验参数，再以 `format_type="XML"` 调用 `Bio.Blast.qblast()`（较新的 `Bio.Blast` 接口中的模块级函数），完整读取 XML，可选保存到 `output_xml`，并用旧式的 `Bio.Blast.NCBIXML.read()` 解析成简单的汇总字典。
- **API：** `blast_qblast_search(program: 'str', database: 'str', sequence: 'str', *, hitlist_size: 'int' = 50, expect: 'float' = 10.0, entrez_query: 'str' = '(none)', output_xml: 'str | None' = None, email: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 Biopython（extra `dnakit[bio-formats]`）以及访问 NCBI BLAST（QBLAST）的网络。不需要密钥，邮箱地址可选。检索会在 NCBI 服务器上排队，可能需要几分钟。下面的示例传入的参数会在发出任何请求之前被拒绝，因此可以离线运行。
- **输入：** `program` 取 `blastn`、`blastp`、`blastx`、`tblastn`、`tblastx` 之一。`database` 是 `nt` 这样的短名称。`sequence` 是不超过 100000 个字符的非空字符串（不检查是否为 DNA 字母，尽管文档字符串写的是 DNA 查询）。`hitlist_size` 是 1 到 500 的整数，`expect` 是正数，`entrez_query` 不超过 512 个字符，`output_xml` 是可选的新文件路径，`email` 是可选的邮箱地址。
- **输出：** 所有联网函数都返回一个字典，包含 `status`（`completed`）、`service`、`upstream`（实际调用的 Biopython 函数）、`parameters`（校验后的请求参数）、`retrieved_utc`、`biopython_version`、`output_path`（路径或 `None`）和 `data`。`parameters` 记录的是查询长度而不是查询序列本身。`data` 含 `query_id`、`query`、`query_letters`、`database`、`alignment_count` 和 `alignments`；每个比对含 `hit_id`、`hit_def`、`length` 和 `hsps`（`score`、`bits`、`expect`、`identities`、`align_length`、`query_start`、`query_end`、`subject_start`、`subject_end`）。给出 `output_xml` 时，原始 XML 也会写入该文件。
- **示例代码：**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.sequence import blast_qblast_search

try:
    blast_qblast_search("blastz", "nt", "ACGT", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error)
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] program must be one of blastn, blastp, blastx, tblastn, tblastx.
```

联网并安装依赖后的用法：

```python
from dnakit.sequence import blast_qblast_search

result = blast_qblast_search(
    "blastn", "nt", "ATGGCGGCGCTGAGCGGTGGCGGTGGTGGCGGCGCGGAGCCGGGCCAGGCTCTGTTCAACGGGGACATGGAGCCCGAGGCCGGC",
    hitlist_size=5, email="you@example.org",
)
print(result["data"]["alignment_count"])
```

## 2) `dnakit.sequence.gc_skew_window_plot()`

- **作用：** 按固定窗口计算 DNA 序列的 GC 偏斜，并把偏斜及其累积和的折线图保存为新的 PNG 文件。
- **处理方法：** 检查序列只含 IUPAC DNA 字母，用 `Bio.SeqUtils.GC_skew()` 按 (G - C) / (G + C) 计算每个窗口的值，用 `itertools.accumulate` 求累积和，并用 matplotlib 对象接口（`Figure` 和 `FigureCanvasAgg`，不使用 pyplot 全局状态）画出两条线。它不调用需要 Tk 显示的 `Bio.SeqUtils.xGC_skew()`。`GC_skew()` 会丢弃末尾不足一个窗口的部分。
- **API：** `gc_skew_window_plot(sequence: 'str', image_path: 'str', *, window: 'int' = 1000, dpi: 'int' = 150, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 Biopython（extra `dnakit[bio-formats]`）和 matplotlib（extra `dnakit[bio-plots]`）。离线运行，不需要显示器。
- **输入：** `sequence` 是由 IUPAC DNA 字母组成的非空字符串，最长 10000000 个字符。`image_path` 必须是已有文件夹中以 `.png` 结尾的新路径。`window` 是 1 到 10000000 的整数（默认 1000），`dpi` 是 50 到 600 的整数。
- **输出：** 返回字典，含 `status`（`completed`）、`image_path`、`format`（`png`）和 `data`。`data` 含 `sequence_length`、`window`、`window_count`、`window_starts`、`gc_skew` 和 `accumulated_gc_skew`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import gc_skew_window_plot

image = Path(tempfile.mkdtemp()) / "skew.png"
result = gc_skew_window_plot("GGGGGGCCAAGGGCCCCCCTT", image, window=7, verbose=False)
data = result["data"]
print(result["status"], result["format"], image.read_bytes()[:4])
print(data["sequence_length"], data["window_count"], data["window_starts"])
print([round(value, 3) for value in data["gc_skew"]])
print([round(value, 3) for value in data["accumulated_gc_skew"]])
```

- **示例输出：**

```text
completed png b'\x89PNG'
21 3 [0, 7, 14]
[0.714, 0.2, -1.0]
[0.714, 0.914, -0.086]
```

## 3) `dnakit.sequence.kegg_convert_identifiers()`

- **作用：** 在 KEGG 标识符与外部数据库标识符之间互相转换。
- **处理方法：** 校验参数后调用 `Bio.KEGG.REST.kegg_conv()`，把以制表符分隔的返回文本拆成行。
- **API：** `kegg_convert_identifiers(target_db: 'str', source: 'str | Sequence[str]', *, option: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 Biopython（extra `dnakit[bio-formats]`）以及访问 KEGG REST 服务（rest.kegg.jp）的网络。不需要密钥。下面的示例传入的参数会在发出任何请求之前被拒绝，因此可以离线运行；要真正运行该函数需要联网。
- **输入：** `target_db` 是 `ncbi-geneid` 或 `hsa` 这样的短名称。`source` 是数据库名或单个标识符字符串，或包含 1 到 100 个标识符的列表（只允许字母、数字和 `_ . : | + -`）。`option` 取 `None`、`turtle` 或 `n-triple`。
- **输出：** 所有联网函数都返回一个字典，包含 `status`（`completed`）、`service`、`upstream`（实际调用的 Biopython 函数）、`parameters`（校验后的请求参数）、`retrieved_utc`、`biopython_version`、`output_path`（路径或 `None`）和 `data`。`data` 含 `rows`（字符串列表的列表）和 `row_count`。
- **示例代码：**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.sequence import kegg_convert_identifiers

try:
    kegg_convert_identifiers("ncbi-geneid", "hsa:10458", option="xml", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error)
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] option must be one of None, turtle, n-triple.
```

联网并安装依赖后的用法：

```python
from dnakit.sequence import kegg_convert_identifiers

result = kegg_convert_identifiers("ncbi-geneid", ["hsa:10458", "hsa:10459"])
print(result["data"]["rows"])
```

## 4) `dnakit.sequence.kegg_database_info()`

- **作用：** 描述一个 KEGG 数据库，包括版本和统计信息。
- **处理方法：** 校验名称后调用 `Bio.KEGG.REST.kegg_info()`，原样返回得到的文本。
- **API：** `kegg_database_info(database: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 Biopython（extra `dnakit[bio-formats]`）以及访问 KEGG REST 服务（rest.kegg.jp）的网络。不需要密钥。下面的示例传入的参数会在发出任何请求之前被拒绝，因此可以离线运行；要真正运行该函数需要联网。
- **输入：** `database` 是 `pathway` 或 `hsa` 这样的 KEGG 数据库或物种标识符：只允许字母、数字和 `_ . : | + -`，最长 256 个字符。
- **输出：** 所有联网函数都返回一个字典，包含 `status`（`completed`）、`service`、`upstream`（实际调用的 Biopython 函数）、`parameters`（校验后的请求参数）、`retrieved_utc`、`biopython_version`、`output_path`（路径或 `None`）和 `data`。`data` 含 `text` 和 `line_count`。
- **示例代码：**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.sequence import kegg_database_info

try:
    kegg_database_info("bad name", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error)
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] database must be an identifier without spaces, commas or slashes.
```

联网并安装依赖后的用法：

```python
from dnakit.sequence import kegg_database_info

result = kegg_database_info("pathway")
print(result["data"]["text"].splitlines()[0])
```

## 5) `dnakit.sequence.kegg_find_entries()`

- **作用：** 按关键词、分子式、精确质量或分子量检索 KEGG 数据库。
- **处理方法：** 校验参数后调用 `Bio.KEGG.REST.kegg_find()`，把以制表符分隔的返回文本拆成行。
- **API：** `kegg_find_entries(database: 'str', query: 'str', *, option: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 Biopython（extra `dnakit[bio-formats]`）以及访问 KEGG REST 服务（rest.kegg.jp）的网络。不需要密钥。下面的示例传入的参数会在发出任何请求之前被拒绝，因此可以离线运行；要真正运行该函数需要联网。
- **输入：** `database` 是 `compound` 或 `genes` 这样的短名称。`query` 是不超过 512 个字符的非空字符串。`option` 取 `None`、`formula`、`exact_mass` 或 `mol_weight`。
- **输出：** 所有联网函数都返回一个字典，包含 `status`（`completed`）、`service`、`upstream`（实际调用的 Biopython 函数）、`parameters`（校验后的请求参数）、`retrieved_utc`、`biopython_version`、`output_path`（路径或 `None`）和 `data`。`data` 含 `rows` 和 `row_count`。
- **示例代码：**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.sequence import kegg_find_entries

try:
    kegg_find_entries("compound", "C7H10O5", option="weight", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error)
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] option must be one of None, formula, exact_mass, mol_weight.
```

联网并安装依赖后的用法：

```python
from dnakit.sequence import kegg_find_entries

result = kegg_find_entries("compound", "C7H10O5", option="formula")
print(result["data"]["row_count"])
```

## 6) `dnakit.sequence.kegg_get_entries()`

- **作用：** 以 `option` 指定的格式取回 KEGG 条目文本。
- **处理方法：** 校验参数后调用 `Bio.KEGG.REST.kegg_get()`，原样返回文本。
- **API：** `kegg_get_entries(entries: 'str | Sequence[str]', *, option: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 Biopython（extra `dnakit[bio-formats]`）以及访问 KEGG REST 服务（rest.kegg.jp）的网络。不需要密钥。下面的示例传入的参数会在发出任何请求之前被拒绝，因此可以离线运行；要真正运行该函数需要联网。
- **输入：** `entries` 是单个标识符字符串，或包含 1 到 10 个标识符的列表。`option` 取 `None`、`aaseq`、`ntseq`、`mol`、`kcf`、`kgml` 或 `json`；`image` 会被拒绝。
- **输出：** 所有联网函数都返回一个字典，包含 `status`（`completed`）、`service`、`upstream`（实际调用的 Biopython 函数）、`parameters`（校验后的请求参数）、`retrieved_utc`、`biopython_version`、`output_path`（路径或 `None`）和 `data`。`data` 含 `text` 和 `line_count`。
- **示例代码：**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.sequence import kegg_get_entries

try:
    kegg_get_entries("hsa:10458", option="image", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error)
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] option must be one of None, aaseq, ntseq, mol, kcf, kgml, json; image is not supported.
```

联网并安装依赖后的用法：

```python
from dnakit.sequence import kegg_get_entries

result = kegg_get_entries(["hsa:10458"], option="aaseq")
print(result["data"]["text"].splitlines()[0])
```

## 7) `dnakit.sequence.kegg_link_entries()`

- **作用：** 查找与某个数据库或给定条目相关联的 KEGG 条目。
- **处理方法：** 校验参数后调用 `Bio.KEGG.REST.kegg_link()`，把以制表符分隔的返回文本拆成行。可用的 `option` 与 `kegg_convert_identifiers()` 相同。
- **API：** `kegg_link_entries(target_db: 'str', source: 'str | Sequence[str]', *, option: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 Biopython（extra `dnakit[bio-formats]`）以及访问 KEGG REST 服务（rest.kegg.jp）的网络。不需要密钥。下面的示例传入的参数会在发出任何请求之前被拒绝，因此可以离线运行；要真正运行该函数需要联网。
- **输入：** `target_db` 是 `pathway` 这样的短名称。`source` 是数据库名或单个标识符字符串，或包含 1 到 100 个标识符的列表。`option` 取 `None`、`turtle` 或 `n-triple`。
- **输出：** 所有联网函数都返回一个字典，包含 `status`（`completed`）、`service`、`upstream`（实际调用的 Biopython 函数）、`parameters`（校验后的请求参数）、`retrieved_utc`、`biopython_version`、`output_path`（路径或 `None`）和 `data`。`data` 含 `rows` 和 `row_count`。
- **示例代码：**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.sequence import kegg_link_entries

try:
    kegg_link_entries("pathway", "hsa:10458", option="rdf", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error)
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] option must be one of None, turtle, n-triple.
```

联网并安装依赖后的用法：

```python
from dnakit.sequence import kegg_link_entries

result = kegg_link_entries("pathway", ["hsa:10458"])
print(result["data"]["rows"])
```

## 8) `dnakit.sequence.kegg_list_entries()`

- **作用：** 列出某个 KEGG 数据库的条目，或某个物种的通路。
- **处理方法：** 校验参数后用数据库名和可选的物种调用 `Bio.KEGG.REST.kegg_list()`，把以制表符分隔的返回文本拆成行。
- **API：** `kegg_list_entries(database: 'str | Sequence[str]', *, organism: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 Biopython（extra `dnakit[bio-formats]`）以及访问 KEGG REST 服务（rest.kegg.jp）的网络。不需要密钥。下面的示例传入的参数会在发出任何请求之前被拒绝，因此可以离线运行；要真正运行该函数需要联网。
- **输入：** `database` 是数据库名或单个条目字符串，或包含 1 到 100 个标识符的列表。`organism` 为 `None` 或 `hsa` 这样的短名称，作为 list 请求的第二个参数传给 KEGG。
- **输出：** 所有联网函数都返回一个字典，包含 `status`（`completed`）、`service`、`upstream`（实际调用的 Biopython 函数）、`parameters`（校验后的请求参数）、`retrieved_utc`、`biopython_version`、`output_path`（路径或 `None`）和 `data`。`data` 含 `rows` 和 `row_count`。
- **示例代码：**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.sequence import kegg_list_entries

try:
    kegg_list_entries("pathway", organism="hsa!", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error)
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] organism must be a short name such as 'pubmed', 'nt' or 'hsa'.
```

联网并安装依赖后的用法：

```python
from dnakit.sequence import kegg_list_entries

result = kegg_list_entries("pathway", organism="hsa")
print(result["data"]["row_count"])
```

## 9) `dnakit.sequence.ncbi_efetch()`

- **作用：** 以文本形式取回 NCBI Entrez 记录，可选写入新文件。
- **处理方法：** 校验参数后用逗号连接标识符调用 `Bio.Entrez.efetch()`，把返回内容读成文本；给出 `output_path` 时写入该新文件。
- **API：** `ncbi_efetch(db: 'str', ids: 'Sequence[str]', *, rettype: 'str' = 'fasta', retmode: 'str' = 'text', output_path: 'str | None' = None, email: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 Biopython（extra `dnakit[bio-formats]`）以及访问 NCBI Entrez E-utilities 的网络。不需要密钥，邮箱地址可选。下面的示例传入的参数会在发出任何请求之前被拒绝，因此可以离线运行；要真正运行该函数需要联网。
- **输入：** `db` 是 `nuccore` 或 `pubmed` 这样的短名称。`ids` 是包含 1 到 200 个标识符字符串的列表。`rettype`（默认 `fasta`）和 `retmode`（默认 `text`）是短名称。`output_path` 为 `None` 或尚不存在的路径。`email` 为 `None` 或邮箱地址。
- **输出：** 所有联网函数都返回一个字典，包含 `status`（`completed`）、`service`、`upstream`（实际调用的 Biopython 函数）、`parameters`（校验后的请求参数）、`retrieved_utc`、`biopython_version`、`output_path`（路径或 `None`）和 `data`。`data` 含 `text` 和 `characters`；`output_path` 是写出的文件或 `None`。
- **示例代码：**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.sequence import ncbi_efetch

try:
    ncbi_efetch("nuccore", [], verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error)
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] ids must hold 1 to 200 identifiers.
```

联网并安装依赖后的用法：

```python
from dnakit.sequence import ncbi_efetch

result = ncbi_efetch("nuccore", ["NM_000546.6"], rettype="fasta", output_path="tp53.fa")
print(result["data"]["characters"])
```

## 10) `dnakit.sequence.ncbi_einfo()`

- **作用：** 列出 NCBI Entrez 数据库，或描述其中一个数据库。
- **处理方法：** 校验参数后带或不带 `db` 调用 `Bio.Entrez.einfo()`，用 `Bio.Entrez.read()` 解析 XML 返回并转换为普通 Python 类型。
- **API：** `ncbi_einfo(*, db: 'str | None' = None, email: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 Biopython（extra `dnakit[bio-formats]`）以及访问 NCBI Entrez E-utilities 的网络。不需要密钥，邮箱地址可选。下面的示例传入的参数会在发出任何请求之前被拒绝，因此可以离线运行；要真正运行该函数需要联网。
- **输入：** `db` 为 `None`（列出所有数据库）或 `pubmed` 这样的短名称。`email` 为 `None` 或邮箱地址。
- **输出：** 所有联网函数都返回一个字典，包含 `status`（`completed`）、`service`、`upstream`（实际调用的 Biopython 函数）、`parameters`（校验后的请求参数）、`retrieved_utc`、`biopython_version`、`output_path`（路径或 `None`）和 `data`。`data` 是解析后的返回内容：不给 `db` 时含 `DbList`，给出 `db` 时含 `DbInfo`。
- **示例代码：**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.sequence import ncbi_einfo

try:
    ncbi_einfo(email="not-an-email", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error)
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] email must be an e-mail address when it is given.
```

联网并安装依赖后的用法：

```python
from dnakit.sequence import ncbi_einfo

result = ncbi_einfo()
print(len(result["data"]["DbList"]))
```

## 11) `dnakit.sequence.ncbi_elink()`

- **作用：** 查找 NCBI Entrez 中与另一数据库标识符相关联的记录。
- **处理方法：** 校验参数后调用 `Bio.Entrez.elink()`，用 `Bio.Entrez.read()` 解析 XML 返回并转换为普通 Python 类型。
- **API：** `ncbi_elink(dbfrom: 'str', db: 'str', ids: 'Sequence[str]', *, linkname: 'str | None' = None, email: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 Biopython（extra `dnakit[bio-formats]`）以及访问 NCBI Entrez E-utilities 的网络。不需要密钥，邮箱地址可选。下面的示例传入的参数会在发出任何请求之前被拒绝，因此可以离线运行；要真正运行该函数需要联网。
- **输入：** `dbfrom` 和 `db` 是短名称。`ids` 是包含 1 到 200 个标识符的列表，每个标识符不得含空格、逗号或斜杠。`linkname` 为 `None` 或短名称。`email` 为 `None` 或邮箱地址。
- **输出：** 所有联网函数都返回一个字典，包含 `status`（`completed`）、`service`、`upstream`（实际调用的 Biopython 函数）、`parameters`（校验后的请求参数）、`retrieved_utc`、`biopython_version`、`output_path`（路径或 `None`）和 `data`。`data` 是解析后的返回内容，为列表，每个请求一项，找到链接时含 `LinkSetDb`。
- **示例代码：**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.sequence import ncbi_elink

try:
    ncbi_elink("gene", "nuccore", ["1,2"], verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error)
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] ids must be an identifier without spaces, commas or slashes.
```

联网并安装依赖后的用法：

```python
from dnakit.sequence import ncbi_elink

result = ncbi_elink("gene", "nuccore", ["7157"], linkname="gene_nuccore_refseqrna")
print(len(result["data"]))
```

## 12) `dnakit.sequence.ncbi_epost()`

- **作用：** 把标识符提交到 NCBI 历史服务器并返回 WebEnv 键。
- **处理方法：** 校验参数后调用 `Bio.Entrez.epost()`，用 `Bio.Entrez.read()` 解析 XML 返回并转换为普通 Python 类型。
- **API：** `ncbi_epost(db: 'str', ids: 'Sequence[str]', *, email: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 Biopython（extra `dnakit[bio-formats]`）以及访问 NCBI Entrez E-utilities 的网络。不需要密钥，邮箱地址可选。下面的示例传入的参数会在发出任何请求之前被拒绝，因此可以离线运行；要真正运行该函数需要联网。
- **输入：** `db` 是短名称。`ids` 是包含 1 到 200 个标识符字符串的列表。`email` 为 `None` 或邮箱地址。
- **输出：** 所有联网函数都返回一个字典，包含 `status`（`completed`）、`service`、`upstream`（实际调用的 Biopython 函数）、`parameters`（校验后的请求参数）、`retrieved_utc`、`biopython_version`、`output_path`（路径或 `None`）和 `data`。`data` 是解析后的返回内容，含 `WebEnv` 和 `QueryKey`。
- **示例代码：**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.sequence import ncbi_epost

try:
    ncbi_epost("nuccore", [str(number) for number in range(201)], verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error)
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] ids must hold 1 to 200 identifiers.
```

联网并安装依赖后的用法：

```python
from dnakit.sequence import ncbi_epost

result = ncbi_epost("nuccore", ["NM_000546.6", "NM_007294.4"])
print(sorted(result["data"]))
```

## 13) `dnakit.sequence.ncbi_esearch()`

- **作用：** 检索 NCBI Entrez 数据库，返回总命中数和匹配的 UID。
- **处理方法：** 校验参数后调用 `Bio.Entrez.esearch()`，用 `Bio.Entrez.read()` 解析 XML 返回，转换为普通 Python 类型，并补充 `count` 和 `ids`。
- **API：** `ncbi_esearch(db: 'str', term: 'str', *, retmax: 'int' = 20, email: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 Biopython（extra `dnakit[bio-formats]`）以及访问 NCBI Entrez E-utilities 的网络。不需要密钥，邮箱地址可选。下面的示例传入的参数会在发出任何请求之前被拒绝，因此可以离线运行；要真正运行该函数需要联网。
- **输入：** `db` 是 `pubmed` 这样的短名称。`term` 是不超过 4096 个字符的非空字符串。`retmax` 是 0 到 10000 的整数（默认 20）。`email` 为 `None` 或邮箱地址。
- **输出：** 所有联网函数都返回一个字典，包含 `status`（`completed`）、`service`、`upstream`（实际调用的 Biopython 函数）、`parameters`（校验后的请求参数）、`retrieved_utc`、`biopython_version`、`output_path`（路径或 `None`）和 `data`。`data` 是解析后的返回内容（`Count`、`RetMax`、`IdList` 等），另加 `count`（整数）和 `ids`（字符串列表）。
- **示例代码：**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.sequence import ncbi_esearch

try:
    ncbi_esearch("pubmed", "BRCA1", retmax=10001, verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error)
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] retmax must be an integer from 0 to 10000.
```

联网并安装依赖后的用法：

```python
from dnakit.sequence import ncbi_esearch

result = ncbi_esearch("pubmed", "BRCA1[Title] AND 2020[PDAT]", retmax=5)
print(result["data"]["count"], result["data"]["ids"])
```

## 14) `dnakit.sequence.ncbi_esummary()`

- **作用：** 取回一组标识符的 NCBI Entrez 文档摘要。
- **处理方法：** 校验参数后调用 `Bio.Entrez.esummary()`，用 `Bio.Entrez.read()` 解析 XML 返回并转换为普通 Python 类型。
- **API：** `ncbi_esummary(db: 'str', ids: 'Sequence[str]', *, email: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 Biopython（extra `dnakit[bio-formats]`）以及访问 NCBI Entrez E-utilities 的网络。不需要密钥，邮箱地址可选。下面的示例传入的参数会在发出任何请求之前被拒绝，因此可以离线运行；要真正运行该函数需要联网。
- **输入：** `db` 是短名称。`ids` 是包含 1 到 200 个标识符字符串的列表。`email` 为 `None` 或邮箱地址。
- **输出：** 所有联网函数都返回一个字典，包含 `status`（`completed`）、`service`、`upstream`（实际调用的 Biopython 函数）、`parameters`（校验后的请求参数）、`retrieved_utc`、`biopython_version`、`output_path`（路径或 `None`）和 `data`。`data` 是解析后的返回内容，结构取决于数据库（多数数据库为摘要记录列表）。
- **示例代码：**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.sequence import ncbi_esummary

try:
    ncbi_esummary("pub med", ["1"], verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error)
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] db must be a short name such as 'pubmed', 'nt' or 'hsa'.
```

联网并安装依赖后的用法：

```python
from dnakit.sequence import ncbi_esummary

result = ncbi_esummary("pubmed", ["31452104"])
print(result["data"][0]["Title"])
```

## 15) `dnakit.sequence.ncbiwww_qblast_search()`

- **作用：** 通过较早的 `NCBIWWW` 模块在 NCBI 服务器上运行 BLAST 检索并汇总命中结果。
- **处理方法：** 与 `blast_qblast_search()` 相同，区别只在于实际调用的是 `Bio.Blast.NCBIWWW.qblast()`（`format_type="XML"`）。XML 会被完整读取，可选保存到 `output_xml`，并用 `Bio.Blast.NCBIXML.read()` 解析成简单的汇总字典。
- **API：** `ncbiwww_qblast_search(program: 'str', database: 'str', sequence: 'str', *, hitlist_size: 'int' = 50, expect: 'float' = 10.0, entrez_query: 'str' = '(none)', output_xml: 'str | None' = None, email: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 Biopython（extra `dnakit[bio-formats]`）以及访问 NCBI BLAST（QBLAST）的网络。不需要密钥，邮箱地址可选。检索会在 NCBI 服务器上排队，可能需要几分钟。下面的示例传入的参数会在发出任何请求之前被拒绝，因此可以离线运行。
- **输入：** 参数和限制与 `blast_qblast_search()` 相同。`program` 取 `blastn`、`blastp`、`blastx`、`tblastn`、`tblastx` 之一。`database` 是 `nt` 这样的短名称。`sequence` 是不超过 100000 个字符的非空字符串（不检查是否为 DNA 字母，尽管文档字符串写的是 DNA 查询）。`hitlist_size` 是 1 到 500 的整数，`expect` 是正数，`entrez_query` 不超过 512 个字符，`output_xml` 是可选的新文件路径，`email` 是可选的邮箱地址。
- **输出：** 所有联网函数都返回一个字典，包含 `status`（`completed`）、`service`、`upstream`（实际调用的 Biopython 函数）、`parameters`（校验后的请求参数）、`retrieved_utc`、`biopython_version`、`output_path`（路径或 `None`）和 `data`。`parameters` 记录的是查询长度而不是查询序列本身。`data` 含 `query_id`、`query`、`query_letters`、`database`、`alignment_count` 和 `alignments`；每个比对含 `hit_id`、`hit_def`、`length` 和 `hsps`（`score`、`bits`、`expect`、`identities`、`align_length`、`query_start`、`query_end`、`subject_start`、`subject_end`）。给出 `output_xml` 时，原始 XML 也会写入该文件。
- **示例代码：**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.sequence import ncbiwww_qblast_search

try:
    ncbiwww_qblast_search("blastn", "nt", "   ", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error)
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] sequence must be a non-empty string.
```

联网并安装依赖后的用法：

```python
from dnakit.sequence import ncbiwww_qblast_search

result = ncbiwww_qblast_search(
    "blastn", "nt", "ATGGCGGCGCTGAGCGGTGGCGGTGGTGGCGGCGCGGAGCCGGGCCAGGCTCTGTTCAACGGGGACATGGAGCCCGAGGCCGGC",
    hitlist_size=5, expect=1e-5, output_xml="blast_result.xml",
)
print(result["data"]["alignment_count"])
```

## 16) `dnakit.sequence.phylo_draw_ascii_text()`

- **作用：** 把 Newick 树画成 ASCII 字符图，并保存为新的文本文件。
- **处理方法：** 用 `Bio.Phylo.read()`（`newick` 格式）解析树，用 `Bio.Phylo.draw_ascii()` 按指定宽度画入字符串缓冲区，把文本写入 `output_path` 并返回。
- **API：** `phylo_draw_ascii_text(newick: 'str', output_path: 'str', *, column_width: 'int' = 80, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 Biopython（extra `dnakit[bio-formats]`）。离线运行，不需要 matplotlib。
- **输入：** `newick` 是不超过 5000000 个字符的非空 Newick 字符串。`output_path` 必须是已有文件夹中的新文件路径。`column_width` 是 20 到 1000 的整数（默认 80）。
- **输出：** 返回字典，含 `status`（`completed`）、`output_path` 和 `data`。`data` 含 `text`（字符图）、`line_count`、`tip_names` 和 `column_width`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import phylo_draw_ascii_text

target = Path(tempfile.mkdtemp()) / "tree.txt"
result = phylo_draw_ascii_text("((A:1,B:2)0.9:1,C:3);", target, column_width=40, verbose=False)
data = result["data"]
for line in data["text"].splitlines():
    if line.strip():
        print("|" + line.rstrip())
print(data["tip_names"], data["column_width"])
print(target.read_text() == data["text"])
```

- **示例输出：**

```text
|              ___________ A
|  ___________|
|_|           |_______________________ B
| |
| |___________________________________ C
['A', 'B', 'C'] 40
True
```

## 17) `dnakit.sequence.phylo_draw_newick_png()`

- **作用：** 把 Newick 树画成新的 PNG 图像，并返回树的数据。
- **处理方法：** 用 `Bio.Phylo.read()`（`newick` 格式）解析树，用 `Bio.Phylo.draw()` 在使用 Agg 画布的 matplotlib `Figure`（6 x 4 英寸，不使用 pyplot 状态，不弹窗）上绘制并保存 PNG。同时用 `Bio.Phylo.write()` 把树写回，得到规范化的 Newick 字符串。
- **API：** `phylo_draw_newick_png(newick: 'str', image_path: 'str', *, show_confidence: 'bool' = True, dpi: 'int' = 150, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 Biopython（extra `dnakit[bio-formats]`）和 matplotlib（extra `dnakit[bio-plots]`）。离线运行，不需要显示器。
- **输入：** `newick` 是不超过 5000000 个字符的非空 Newick 字符串。`image_path` 必须是已有文件夹中以 `.png` 结尾的新路径。`show_confidence`（默认 `True`）为带置信度的内部节点加标签。`dpi` 是 50 到 600 的整数（默认 150）。
- **输出：** 返回字典，含 `status`（`completed`）、`image_path`、`format`（`png`）和 `data`。`data` 含 `newick`（规范化后）、`tip_names`、`terminal_count`、`clade_count` 和 `total_branch_length`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import phylo_draw_newick_png

image = Path(tempfile.mkdtemp()) / "tree.png"
result = phylo_draw_newick_png("((A:1,B:2)0.9:1,C:3);", image, verbose=False)
data = result["data"]
print(result["status"], result["format"], image.read_bytes()[:4])
print(data["newick"])
print(data["tip_names"], data["terminal_count"], data["clade_count"], data["total_branch_length"])
```

- **示例输出：**

```text
completed png b'\x89PNG'
((A:1,B:2)0.90:1,C:3):0;
['A', 'B', 'C'] 3 5 7.0
```

## 18) `dnakit.sequence.togows_convert_data()`

- **作用：** 借助 TogoWS 服务在文件格式之间转换记录。
- **处理方法：** 校验参数后调用 `Bio.TogoWS.convert()`，把返回内容读成文本。
- **API：** `togows_convert_data(data: 'str', in_format: 'str', out_format: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 Biopython（extra `dnakit[bio-formats]`）以及访问 TogoWS（togows.dbcls.jp）的网络。不需要密钥。下面的示例传入的参数会在发出任何请求之前被拒绝，因此可以离线运行；要真正运行该函数需要联网。
- **输入：** `data` 是不超过 10000000 个字符的非空字符串。`in_format` 和 `out_format` 是 TogoWS 支持的短格式名（例如 `genbank` 和 `fasta`）。
- **输出：** 所有联网函数都返回一个字典，包含 `status`（`completed`）、`service`、`upstream`（实际调用的 Biopython 函数）、`parameters`（校验后的请求参数）、`retrieved_utc`、`biopython_version`、`output_path`（路径或 `None`）和 `data`。`parameters` 记录两种格式和输入长度。`data` 含 `text` 和 `characters`。
- **示例代码：**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.sequence import togows_convert_data

try:
    togows_convert_data(">x\nACGT", "1bad", "fasta", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error)
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] in_format must be a short name such as 'pubmed', 'nt' or 'hsa'.
```

联网并安装依赖后的用法：

```python
from dnakit.sequence import togows_convert_data

genbank_text = open("record.gb").read()
result = togows_convert_data(genbank_text, "genbank", "fasta")
print(result["data"]["text"].splitlines()[0])
```

## 19) `dnakit.sequence.togows_count_matches()`

- **作用：** 统计某个查询在 TogoWS 数据库中的匹配数。
- **处理方法：** 校验参数后调用 `Bio.TogoWS.search_count()`，以整数返回计数。
- **API：** `togows_count_matches(database: 'str', query: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 Biopython（extra `dnakit[bio-formats]`）以及访问 TogoWS（togows.dbcls.jp）的网络。不需要密钥。下面的示例传入的参数会在发出任何请求之前被拒绝，因此可以离线运行；要真正运行该函数需要联网。
- **输入：** `database` 是 `ncbi-nucleotide` 或 `pubmed` 这样的短名称。`query` 是不超过 1024 个字符的非空字符串。
- **输出：** 所有联网函数都返回一个字典，包含 `status`（`completed`）、`service`、`upstream`（实际调用的 Biopython 函数）、`parameters`（校验后的请求参数）、`retrieved_utc`、`biopython_version`、`output_path`（路径或 `None`）和 `data`。`data` 含 `count`。
- **示例代码：**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.sequence import togows_count_matches

try:
    togows_count_matches("ncbi-nucleotide", "  ", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error)
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] query must be a non-empty string.
```

联网并安装依赖后的用法：

```python
from dnakit.sequence import togows_count_matches

result = togows_count_matches("pubmed", "lung and cancer")
print(result["data"]["count"])
```

## 20) `dnakit.sequence.togows_fetch_entry()`

- **作用：** 取回一条 TogoWS 记录，可指定格式或只取一个字段。
- **处理方法：** 校验参数后调用 `Bio.TogoWS.entry()`，把返回内容读成文本。
- **API：** `togows_fetch_entry(database: 'str', identifier: 'str', *, output_format: 'str | None' = None, field: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 Biopython（extra `dnakit[bio-formats]`）以及访问 TogoWS（togows.dbcls.jp）的网络。不需要密钥。下面的示例传入的参数会在发出任何请求之前被拒绝，因此可以离线运行；要真正运行该函数需要联网。
- **输入：** `database` 是短名称。`identifier` 是单个标识符，只允许字母、数字和 `_ . : | + -`（不含斜杠或逗号）。`output_format` 和 `field` 为 `None` 或短名称。
- **输出：** 所有联网函数都返回一个字典，包含 `status`（`completed`）、`service`、`upstream`（实际调用的 Biopython 函数）、`parameters`（校验后的请求参数）、`retrieved_utc`、`biopython_version`、`output_path`（路径或 `None`）和 `data`。`data` 含 `text` 和 `characters`。
- **示例代码：**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.sequence import togows_fetch_entry

try:
    togows_fetch_entry("ncbi-nucleotide", "a/b", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error)
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] identifier must be an identifier without spaces, commas or slashes.
```

联网并安装依赖后的用法：

```python
from dnakit.sequence import togows_fetch_entry

result = togows_fetch_entry("ncbi-nucleotide", "X52280", output_format="fasta")
print(result["data"]["text"].splitlines()[0])
```

## 21) `dnakit.sequence.togows_iterate_ids()`

- **作用：** 为一个 TogoWS 查询收集数量有上限的标识符。
- **处理方法：** 校验参数后调用 `Bio.TogoWS.search_iter()`，并用 `itertools.islice` 从迭代器中最多取 `limit` 个标识符。
- **API：** `togows_iterate_ids(database: 'str', query: 'str', *, limit: 'int' = 100, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 Biopython（extra `dnakit[bio-formats]`）以及访问 TogoWS（togows.dbcls.jp）的网络。不需要密钥。下面的示例传入的参数会在发出任何请求之前被拒绝，因此可以离线运行；要真正运行该函数需要联网。
- **输入：** `database` 是短名称。`query` 是不超过 1024 个字符的非空字符串。`limit` 是 1 到 1000 的整数（默认 100）。
- **输出：** 所有联网函数都返回一个字典，包含 `status`（`completed`）、`service`、`upstream`（实际调用的 Biopython 函数）、`parameters`（校验后的请求参数）、`retrieved_utc`、`biopython_version`、`output_path`（路径或 `None`）和 `data`。`data` 含 `ids` 和 `returned_count`。
- **示例代码：**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.sequence import togows_iterate_ids

try:
    togows_iterate_ids("ncbi-nucleotide", "x", limit=1001, verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error)
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] limit must be an integer from 1 to 1000.
```

联网并安装依赖后的用法：

```python
from dnakit.sequence import togows_iterate_ids

result = togows_iterate_ids("pubmed", "lung and cancer", limit=25)
print(result["data"]["returned_count"])
```

## 22) `dnakit.sequence.togows_search_ids()`

- **作用：** 检索 TogoWS 数据库，返回一页标识符，或指定格式的原始返回文本。
- **处理方法：** 校验参数后调用 `Bio.TogoWS.search()`。不给 `output_format` 时按空白拆分返回内容得到标识符；给出时返回原始文本。
- **API：** `togows_search_ids(database: 'str', query: 'str', *, offset: 'int | None' = None, limit: 'int | None' = None, output_format: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 Biopython（extra `dnakit[bio-formats]`）以及访问 TogoWS（togows.dbcls.jp）的网络。不需要密钥。下面的示例传入的参数会在发出任何请求之前被拒绝，因此可以离线运行；要真正运行该函数需要联网。
- **输入：** `database` 是短名称。`query` 是不超过 1024 个字符的非空字符串。`offset` 为 `None` 或 1 到 1000000000 的整数（TogoWS 从 1 开始计数）。`limit` 为 `None` 或 1 到 1000000 的整数。`output_format` 为 `None` 或 `json` 这样的短名称。
- **输出：** 所有联网函数都返回一个字典，包含 `status`（`completed`）、`service`、`upstream`（实际调用的 Biopython 函数）、`parameters`（校验后的请求参数）、`retrieved_utc`、`biopython_version`、`output_path`（路径或 `None`）和 `data`。`data` 含 `ids` 和 `count`；给出 `output_format` 时则含 `text`。
- **示例代码：**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.sequence import togows_search_ids

try:
    togows_search_ids("ncbi-nucleotide", "x", offset=0, verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error)
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] offset must be an integer from 1 to 1000000000.
```

联网并安装依赖后的用法：

```python
from dnakit.sequence import togows_search_ids

result = togows_search_ids("pubmed", "lung and cancer", offset=1, limit=10)
print(result["data"]["count"], result["data"]["ids"][:3])
```
