# sourmash：生成、比较、搜索和索引 DNA 签名

这一组函数是 `dnakit.microbiome` 中对 [sourmash](https://github.com/sourmash-bio/sourmash) 命令行的带类型封装。sourmash 把 DNA 压缩成很小的 MinHash 签名，并对签名做比较、搜索、索引和分类。这一组函数涵盖：对 FASTA 文件生成草图（DNA 以及六框翻译后的蛋白）、两两比较与绘图、对数据库的 search、gather、prefetch 和 categorize、签名的集合运算（merge、intersect、subtract、downsample、flatten、filter）、日常管理（describe、manifest、split、extract、rename、picklist 检查）、Mash 的导入与导出，以及最低共同祖先（LCA）分类数据库。每个函数都用固定的参数构造 `sourmash` 命令行，不经过 shell，在输出所在文件夹中运行，并返回可 JSON 序列化的字典；输出文件必须是新路径。它们都不会下载数据或访问服务器。

sourmash 采用 BSD-3-Clause 许可（已核对包元数据），没有 copyleft 限制。请从 `dnakit.microbiome` 导入函数，例如 `from dnakit.microbiome import sketch_dna`。该依赖是可选的，只有调用函数时才会查找程序：如果 `sourmash` 不在 PATH 中且没有给出 `sourmash_path`，函数会抛出 `BackendUnavailableError` 并给出安装命令。所有示例都用 `shutil.which` 找到的路径作为 `sourmash_path`，找不到时退回到 Python 解释器旁边的 `sourmash` 脚本，在临时文件夹中工作，并传入 `verbose=False`。

## 示例数据

示例使用两个单条记录的 FASTA 文件，长度分别为 60 和 65 个碱基。前 35 个碱基相同，所以在 k-mer 大小为 21 时，两个草图（用 `scaled=1` 生成，保留全部哈希）在各自的 40 个和 45 个 k-mer 中共有 15 个。示例会自己把这些文件生成草图并放入临时文件夹。

**示例文件：** `a.fa`

```text
>seqA
GCTAAAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCG
```

**示例文件：** `b.fa`

```text
>seqB
GCTAAAGACAATTACATAACATACACGTCAGCACGCTTAAGGGTTAAGTAAGTGTGATGCATACG
```

## 1) `dnakit.microbiome.check_signatures_against_picklist()`

- **作用：** 检查挑选列表（picklist）中的哪些行存在于签名文件中，并列出缺失的行。
- **处理方法：** 运行 `sourmash sig check <signatures> --picklist <picklist> -o <缺失行 csv> --no-require-manifest`，给出 `output_manifest` 时再加 `-m <manifest>`。退出码 0 和 1 都算完成。函数还会读取缺失行 CSV，只有退出码为 0 且 CSV 中没有数据行时才把 `all_picklist_rows_found` 设为真。
- **API：** `dnakit.microbiome.check_signatures_against_picklist(signatures: 'Sequence[str]', picklist: 'str', output_missing_csv: 'str', *, output_manifest: 'str | None' = None, sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `sourmash` 在 PATH 中（或通过 `sourmash_path` 指定其位置）。它由可选包 sourmash 4.9 及以上版本提供（extra `dnakit[bio-sourmash]`）；缺少时抛出 `BackendUnavailableError`。
- **输入：** `signatures` 是非空的已有签名文件列表。`picklist` 原样传给 sourmash，格式为 `file.csv:column:coltype`（例如 `names.csv:name:name`）；DNAKit 不检查它，相对路径相对于 `output_missing_csv` 所在目录解析，所以最好用绝对路径。`output_missing_csv`（以及给出的 `output_manifest`）必须尚不存在。
- **输出：** 标准结果字典（此处退出码为 0 或 1 时 `status` 都是 `completed`，`outputs` 列出缺失行 CSV 和可选的 manifest），另有布尔值 `all_picklist_rows_found`。
- **示例代码：**

```python
import shutil, sys, tempfile
from pathlib import Path

from dnakit.microbiome import sketch_dna, check_signatures_against_picklist

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())
sketch_dna("a.fa", work / "a.sig", scaled=1, name="sigA first", sourmash_path=tool, verbose=False)
sketch_dna("b.fa", work / "b.sig", scaled=1, name="sigB second", sourmash_path=tool, verbose=False)

(work / "wanted.csv").write_text("name\nsigA first\nnot a signature\n")
result = check_signatures_against_picklist(
    [work / "a.sig", work / "b.sig"], f"{work / 'wanted.csv'}:name:name", work / "missing.csv",
    sourmash_path=tool, verbose=False,
)
print(result["status"], result["return_code"], result["all_picklist_rows_found"])
print((work / "missing.csv").read_text().splitlines())
```

- **示例输出：**

```text
completed 0 False
['name', 'not a signature']
```

## 2) `dnakit.microbiome.classify_with_lca()`

- **作用：** 用最低共同祖先（LCA）数据库给每个查询签名指定分类谱系。
- **处理方法：** 运行 `sourmash lca classify --query <queries> --db <databases> --threshold <n> -o <csv>`，`majority` 为真时再加 `--majority`。每个查询对应一行，报告数据库中命中的哈希所支持的谱系。
- **API：** `dnakit.microbiome.classify_with_lca(queries: 'Sequence[str]', databases: 'Sequence[str]', output_csv: 'str', *, threshold: 'int' = 5, majority: 'bool' = False, sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `sourmash` 在 PATH 中（或通过 `sourmash_path` 指定其位置）。它由可选包 sourmash 4.9 及以上版本提供（extra `dnakit[bio-sourmash]`）；缺少时抛出 `BackendUnavailableError`。
- **输入：** `queries` 和 `databases` 是非空的已有文件列表，数据库是 LCA JSON 文件，例如 `index_lca_database()` 写出的文件。`threshold` 是支持一条谱系所需的最少哈希数（默认 5，DNAKit 不检查）。`majority` 使用多数投票分类，而不是最低共同祖先。`output_csv` 必须尚不存在。
- **输出：** 标准结果字典：`status`（sourmash 退出码为 0 且所有预期输出文件都存在时为 `completed`，否则为 `failed`）、`backend`（`sourmash`）、`command`、`return_code`、`outputs`（实际存在的输出文件）和 `log_tail`（sourmash 输出的最后 4000 个字符，其中含终端控制码）。 CSV 的列为 `ID`、`status` 以及每个分类等级各一列。
- **示例代码：**

```python
import shutil, sys, tempfile
from pathlib import Path
import csv

from dnakit.microbiome import sketch_dna, classify_with_lca, index_lca_database

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())
sketch_dna("a.fa", work / "a.sig", scaled=1, name="sigA first", sourmash_path=tool, verbose=False)
sketch_dna("b.fa", work / "b.sig", scaled=1, name="sigB second", sourmash_path=tool, verbose=False)

(work / "lineages.csv").write_text(
    "accession,superkingdom,phylum,class,order,family,genus,species\n"
    "sigA,Bacteria,Firmicutes,Bacilli,Order1,Family1,Genusone,Genusone alpha\n"
    "sigB,Bacteria,Firmicutes,Bacilli,Order1,Family1,Genusone,Genusone beta\n"
)
index_lca_database([work / "a.sig", work / "b.sig"], work / "lineages.csv", work / "lca.lca.json", scaled=1, sourmash_path=tool, verbose=False)
result = classify_with_lca([work / "a.sig"], [work / "lca.lca.json"], work / "lca.csv", threshold=1, sourmash_path=tool, verbose=False)
rows = list(csv.DictReader(open(work / "lca.csv")))
print(result["status"], len(rows), rows[0]["status"], rows[0]["genus"], rows[0]["species"])
```

- **示例输出：**

```text
completed 1 found Genusone Genusone alpha
```

## 3) `dnakit.microbiome.collect_signature_manifests()`

- **作用：** 把多个签名文件的 manifest 收集成一个独立的 manifest 数据库。
- **处理方法：** 运行 `sourmash sig collect <signatures> -o <manifest> --no-require-manifest`。输出是指向原始文件的 SQLite manifest（`.sqlmf`），签名本身不会被复制。
- **API：** `dnakit.microbiome.collect_signature_manifests(signatures: 'Sequence[str]', output_manifest: 'str', *, sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `sourmash` 在 PATH 中（或通过 `sourmash_path` 指定其位置）。它由可选包 sourmash 4.9 及以上版本提供（extra `dnakit[bio-sourmash]`）；缺少时抛出 `BackendUnavailableError`。
- **输入：** `signatures` 是非空的已有签名文件列表。`output_manifest` 必须尚不存在；使用 `.sqlmf` 扩展名，sourmash 才会写出 SQLite。
- **输出：** 标准结果字典：`status`（sourmash 退出码为 0 且所有预期输出文件都存在时为 `completed`，否则为 `failed`）、`backend`（`sourmash`）、`command`、`return_code`、`outputs`（实际存在的输出文件）和 `log_tail`（sourmash 输出的最后 4000 个字符，其中含终端控制码）。
- **示例代码：**

```python
import shutil, sys, tempfile
from pathlib import Path
import sqlite3

from dnakit.microbiome import sketch_dna, collect_signature_manifests

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())
sketch_dna("a.fa", work / "a.sig", scaled=1, name="sigA first", sourmash_path=tool, verbose=False)
sketch_dna("b.fa", work / "b.sig", scaled=1, name="sigB second", sourmash_path=tool, verbose=False)

result = collect_signature_manifests([work / "a.sig", work / "b.sig"], work / "collection.sqlmf", sourmash_path=tool, verbose=False)
connection = sqlite3.connect(work / "collection.sqlmf")
print(result["status"], connection.execute("select count(*) from sourmash_sketches").fetchone()[0])
connection.close()
```

- **示例输出：**

```text
completed 2
```

## 4) `dnakit.microbiome.index_lca_database()`

- **作用：** 用签名和谱系 CSV 构建最低共同祖先（LCA）分类数据库。
- **处理方法：** 运行 `sourmash lca index <谱系 csv> <database> <signatures> -k <ksize> --scaled <scaled> --dna`，`split_identifiers` 为真时再加 `--split-identifiers`。启用拆分后，只用签名名称的第一个词（到第一个空格或点为止）去匹配谱系标识符。如果数据库名称不以 `.lca.json` 或 `.lca.json.gz` 结尾，sourmash 会自动追加 `.lca.json`，而 DNAKit 只查找你给出的原路径，所以像 `db.json` 这样的名称会被报告为 `failed`（返回码 0，没有 `outputs`），尽管 `db.json.lca.json` 已经写出。请使用以 `.lca.json` 结尾的名称。
- **API：** `dnakit.microbiome.index_lca_database(signatures: 'Sequence[str]', lineages_csv: 'str', output_database: 'str', *, ksize: 'int' = 21, scaled: 'int' = 1000, split_identifiers: 'bool' = True, sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `sourmash` 在 PATH 中（或通过 `sourmash_path` 指定其位置）。它由可选包 sourmash 4.9 及以上版本提供（extra `dnakit[bio-sourmash]`）；缺少时抛出 `BackendUnavailableError`。
- **输入：** `signatures` 是非空的已有签名文件列表，`lineages_csv` 是已有的 CSV，第一列为标识符列（`accession`、`ident` 或 `identifiers`），之后是各分类等级列（`superkingdom,phylum,class,order,family,genus,species`，可选 `strain`）。`ksize` 和 `scaled` 必须与签名一致（默认 21 和 1000）。`output_database` 必须尚不存在，且应以 `.lca.json` 结尾。
- **输出：** 标准结果字典：`status`（sourmash 退出码为 0 且所有预期输出文件都存在时为 `completed`，否则为 `failed`）、`backend`（`sourmash`）、`command`、`return_code`、`outputs`（实际存在的输出文件）和 `log_tail`（sourmash 输出的最后 4000 个字符，其中含终端控制码）。
- **示例代码：**

```python
import shutil, sys, tempfile
from pathlib import Path
import json

from dnakit.microbiome import sketch_dna, index_lca_database

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())
sketch_dna("a.fa", work / "a.sig", scaled=1, name="sigA first", sourmash_path=tool, verbose=False)
sketch_dna("b.fa", work / "b.sig", scaled=1, name="sigB second", sourmash_path=tool, verbose=False)

(work / "lineages.csv").write_text(
    "accession,superkingdom,phylum,class,order,family,genus,species\n"
    "sigA,Bacteria,Firmicutes,Bacilli,Order1,Family1,Genusone,Genusone alpha\n"
    "sigB,Bacteria,Firmicutes,Bacilli,Order1,Family1,Genusone,Genusone beta\n"
)
result = index_lca_database([work / "a.sig", work / "b.sig"], work / "lineages.csv", work / "lca.lca.json", scaled=1, sourmash_path=tool, verbose=False)
data = json.loads((work / "lca.lca.json").read_text())
print(result["status"], data["ksize"], data["scaled"], len(data["lid_to_lineage"]), len(data["hashval_to_idx"]))
```

- **示例输出：**

```text
completed 21 1 2 70
```

## 5) `dnakit.microbiome.ingest_mash_csv()`

- **作用：** 把 Mash CSV 哈希转储转换成 sourmash 签名文件。
- **处理方法：** 运行 `sourmash sig ingest --csv <files> -o <signature>`。每一行输入变成一个 `num` 草图。
- **API：** `dnakit.microbiome.ingest_mash_csv(csv_files: 'Sequence[str]', output_signature: 'str', *, sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `sourmash` 在 PATH 中（或通过 `sourmash_path` 指定其位置）。它由可选包 sourmash 4.9 及以上版本提供（extra `dnakit[bio-sourmash]`）；缺少时抛出 `BackendUnavailableError`。
- **输入：** `csv_files` 是非空的已有文件列表。在 sourmash 4.9 中，每行必须正好有五个以逗号分隔的字段：`murmur64`、`42`、k-mer 大小、名称，以及用单个空格分隔的全部整数哈希。`output_signature` 必须尚不存在。
- **输出：** 标准结果字典：`status`（sourmash 退出码为 0 且所有预期输出文件都存在时为 `completed`，否则为 `failed`）、`backend`（`sourmash`）、`command`、`return_code`、`outputs`（实际存在的输出文件）和 `log_tail`（sourmash 输出的最后 4000 个字符，其中含终端控制码）。
- **示例代码：**

```python
import shutil, sys, tempfile
from pathlib import Path
from sourmash import load_file_as_signatures

from dnakit.microbiome import sketch_dna, ingest_mash_csv

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())

(work / "mash.csv").write_text("murmur64,42,21,mash sample,101 202 303 404 505\n")
result = ingest_mash_csv([work / "mash.csv"], work / "mash.sig", sourmash_path=tool, verbose=False)
(sig,) = load_file_as_signatures(str(work / "mash.sig"))
print(result["status"], sig.minhash.ksize, sig.minhash.num, sorted(sig.minhash.hashes))
```

- **示例输出：**

```text
completed 21 5 [101, 202, 303, 404, 505]
```

## 6) `dnakit.microbiome.categorize_signatures()`

- **作用：** 把每个查询签名归入数据库中最匹配的签名。
- **处理方法：** 运行 `sourmash categorize <database> <queries> -k <ksize> --dna --threshold <t> --csv <csv>`。CSV 没有表头，各列依次是查询文件、查询名称、最佳匹配名称和相似度。
- **API：** `dnakit.microbiome.categorize_signatures(database: 'str', queries: 'Sequence[str]', output_csv: 'str', *, ksize: 'int' = 21, threshold: 'float' = 0.0, sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `sourmash` 在 PATH 中（或通过 `sourmash_path` 指定其位置）。它由可选包 sourmash 4.9 及以上版本提供（extra `dnakit[bio-sourmash]`）；缺少时抛出 `BackendUnavailableError`。
- **输入：** `database` 是已有的签名集合，`queries` 是非空的已有签名文件列表。`threshold` 是 0 到 1 之间的相似度（默认 0.0），`ksize` 默认 21。`output_csv` 必须尚不存在。
- **输出：** 标准结果字典：`status`（sourmash 退出码为 0 且所有预期输出文件都存在时为 `completed`，否则为 `failed`）、`backend`（`sourmash`）、`command`、`return_code`、`outputs`（实际存在的输出文件）和 `log_tail`（sourmash 输出的最后 4000 个字符，其中含终端控制码）。
- **示例代码：**

```python
import shutil, sys, tempfile
from pathlib import Path
import csv

from dnakit.microbiome import sketch_dna, categorize_signatures

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())
sketch_dna("a.fa", work / "a.sig", scaled=1, name="sigA first", sourmash_path=tool, verbose=False)
sketch_dna("b.fa", work / "b.sig", scaled=1, name="sigB second", sourmash_path=tool, verbose=False)

result = categorize_signatures(work / "b.sig", [work / "a.sig"], work / "categories.csv", threshold=0.1, sourmash_path=tool, verbose=False)
rows = list(csv.reader(open(work / "categories.csv")))
print(result["status"], len(rows), rows[0][1], rows[0][2], round(float(rows[0][3]), 3))
```

- **示例输出：**

```text
completed 1 sigA first sigB second 0.214
```

## 7) `dnakit.microbiome.export_signature_to_mash()`

- **作用：** 把签名导出为 Mash JSON 转储格式。
- **处理方法：** 运行 `sourmash sig export <signature> -k <ksize> -o <json>`。JSON 含有 `kmer`、`sketchSize`、`hashType`、`hashBits`、`hashSeed`，以及带哈希值的 `sketches` 列表。
- **API：** `dnakit.microbiome.export_signature_to_mash(signature: 'str', output_json: 'str', *, ksize: 'int' = 21, sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `sourmash` 在 PATH 中（或通过 `sourmash_path` 指定其位置）。它由可选包 sourmash 4.9 及以上版本提供（extra `dnakit[bio-sourmash]`）；缺少时抛出 `BackendUnavailableError`。
- **输入：** `signature` 是已有的签名文件，其中含有所请求 `ksize`（默认 21）的草图。`output_json` 必须尚不存在。
- **输出：** 标准结果字典：`status`（sourmash 退出码为 0 且所有预期输出文件都存在时为 `completed`，否则为 `failed`）、`backend`（`sourmash`）、`command`、`return_code`、`outputs`（实际存在的输出文件）和 `log_tail`（sourmash 输出的最后 4000 个字符，其中含终端控制码）。
- **示例代码：**

```python
import shutil, sys, tempfile
from pathlib import Path
import json

from dnakit.microbiome import sketch_dna, export_signature_to_mash

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())
sketch_dna("a.fa", work / "a.sig", scaled=1, name="sigA first", sourmash_path=tool, verbose=False)
sketch_dna("b.fa", work / "b.sig", scaled=1, name="sigB second", sourmash_path=tool, verbose=False)

result = export_signature_to_mash(work / "a.sig", work / "a_mash.json", sourmash_path=tool, verbose=False)
data = json.loads((work / "a_mash.json").read_text())
print(result["status"], data["kmer"], data["sketchSize"], len(data["sketches"][0]["hashes"]))
```

- **示例输出：**

```text
completed 21 40 40
```

## 8) `dnakit.microbiome.extract_signature_kmers()`

- **作用：** 写出包含签名中所哈希 k-mer 的序列。
- **处理方法：** 运行 `sourmash sig kmers --signatures <signature> --sequences <files> -k <ksize> --save-sequences <output>`。输出包含那些含有哈希在签名中的 k-mer 的序列；`scaled` 较大的签名只含全部哈希的一小部分，所以示例中使用 `scaled=1`。
- **API：** `dnakit.microbiome.extract_signature_kmers(signature: 'str', sequences: 'Sequence[str]', output_sequences: 'str', *, ksize: 'int' = 21, sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `sourmash` 在 PATH 中（或通过 `sourmash_path` 指定其位置）。它由可选包 sourmash 4.9 及以上版本提供（extra `dnakit[bio-sourmash]`）；缺少时抛出 `BackendUnavailableError`。
- **输入：** `signature` 是已有的签名文件，`sequences` 是非空的已有 FASTA 或 FASTQ 文件列表，`ksize` 默认 21。`output_sequences` 必须尚不存在。
- **输出：** 标准结果字典：`status`（sourmash 退出码为 0 且所有预期输出文件都存在时为 `completed`，否则为 `failed`）、`backend`（`sourmash`）、`command`、`return_code`、`outputs`（实际存在的输出文件）和 `log_tail`（sourmash 输出的最后 4000 个字符，其中含终端控制码）。
- **示例代码：**

```python
import shutil, sys, tempfile
from pathlib import Path

from dnakit.microbiome import sketch_dna, extract_signature_kmers

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())
sketch_dna("a.fa", work / "a.sig", scaled=1, name="sigA first", sourmash_path=tool, verbose=False)
sketch_dna("b.fa", work / "b.sig", scaled=1, name="sigB second", sourmash_path=tool, verbose=False)

result = extract_signature_kmers(work / "b.sig", ["a.fa"], work / "matches.fa", sourmash_path=tool, verbose=False)
print(result["status"], (work / "matches.fa").read_text().splitlines()[0])
```

- **示例输出：**

```text
completed >seqA
```

## 9) `dnakit.microbiome.multigather_signatures()`

- **作用：** 对每个查询签名分别对数据库运行 gather 搜索。
- **处理方法：** 运行 `sourmash multigather --query <queries> --db <databases> -k <ksize> --dna --threshold-bp <bp> --output-dir <dir>`。若输出文件夹不存在，DNAKit 会创建它，并在运行前后比较其内容，把新增文件作为 `outputs`；状态只取决于退出码。
- **API：** `dnakit.microbiome.multigather_signatures(queries: 'Sequence[str]', databases: 'Sequence[str]', output_dir: 'str', *, ksize: 'int' = 21, threshold_bp: 'int' = 0, sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `sourmash` 在 PATH 中（或通过 `sourmash_path` 指定其位置）。它由可选包 sourmash 4.9 及以上版本提供（extra `dnakit[bio-sourmash]`）；缺少时抛出 `BackendUnavailableError`。
- **输入：** `queries` 和 `databases` 是非空的已有签名文件列表。`output_dir` 必要时会被创建，也可以已含文件。`ksize` 默认 21，`threshold_bp`（最小重叠碱基数）默认 0。
- **输出：** 标准结果字典，退出码为 0 时 `status` 为 `completed`，`outputs` 列出 `output_dir` 中新出现的文件（每个查询以其源文件命名，各有 `.csv`、`.matches.sig` 和 `.unassigned.sig` 三个文件）。
- **示例代码：**

```python
import shutil, sys, tempfile
from pathlib import Path

from dnakit.microbiome import sketch_dna, multigather_signatures

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())
sketch_dna("a.fa", work / "a.sig", scaled=1, name="sigA first", sourmash_path=tool, verbose=False)
sketch_dna("b.fa", work / "b.sig", scaled=1, name="sigB second", sourmash_path=tool, verbose=False)

result = multigather_signatures([work / "a.sig"], [work / "b.sig"], work / "gathers", sourmash_path=tool, verbose=False)
print(result["status"], [Path(name).suffix for name in result["outputs"]])
```

- **示例输出：**

```text
completed ['.csv', '.sig', '.sig']
```

## 10) `dnakit.microbiome.overlap_signatures()`

- **作用：** 报告两个签名的 Jaccard 相似度、包含度以及共有和独有的哈希数。
- **处理方法：** 运行 `sourmash sig overlap <a> <b>`，并把每一行形如 `名称: 数值` 的输出解析成浮点数字典。每个文件必须只含一个签名，否则 sourmash 失败，`status` 为 `failed`。注意 `size` 和 `sum hashes` 两行每个签名各打印一次，所以字典里保留的是第二个签名的值。
- **API：** `dnakit.microbiome.overlap_signatures(signature_a: 'str', signature_b: 'str', *, sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `sourmash` 在 PATH 中（或通过 `sourmash_path` 指定其位置）。它由可选包 sourmash 4.9 及以上版本提供（extra `dnakit[bio-sourmash]`）；缺少时抛出 `BackendUnavailableError`。
- **输入：** `signature_a` 和 `signature_b` 是已有的签名文件，各含一个签名，k-mer 大小和 scaled 值相同。
- **输出：** `status`、`backend`、`return_code`、`values`（以 sourmash 的标签为键的浮点数，例如 `jaccard similarity`、`number of hashes in common`、`only in first`、`only in second` 和 `total (union)`）以及 `report`（原始输出的最后 4000 个字符）。
- **示例代码：**

```python
import shutil, sys, tempfile
from pathlib import Path

from dnakit.microbiome import sketch_dna, overlap_signatures

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())
sketch_dna("a.fa", work / "a.sig", scaled=1, name="sigA first", sourmash_path=tool, verbose=False)
sketch_dna("b.fa", work / "b.sig", scaled=1, name="sigB second", sourmash_path=tool, verbose=False)

result = overlap_signatures(work / "a.sig", work / "b.sig", sourmash_path=tool, verbose=False)
print(result["status"])
values = result["values"]
for key in ("jaccard similarity", "number of hashes in common", "only in first", "only in second", "total (union)"):
    print(key, values[key])
```

- **示例输出：**

```text
completed
jaccard similarity 0.21429
number of hashes in common 15.0
only in first 25.0
only in second 30.0
total (union) 70.0
```

## 11) `dnakit.microbiome.plot_comparison_matrix()`

- **作用：** 把 sourmash compare 矩阵绘制成聚类相似度热图（PNG）。
- **处理方法：** 运行 `sourmash plot <matrix> --output-dir <dir>`。sourmash 读取 `compare_signatures()` 写在矩阵旁边的 `.labels.txt` 文件，对样本聚类，并用 matplotlib 和 scipy（随 sourmash 一起安装）画出热图。若文件夹不存在，DNAKit 会创建它并报告新增文件。
- **API：** `dnakit.microbiome.plot_comparison_matrix(comparison_matrix: 'str', output_dir: 'str', *, sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `sourmash` 在 PATH 中（或通过 `sourmash_path` 指定其位置）。它由可选包 sourmash 4.9 及以上版本提供（extra `dnakit[bio-sourmash]`）；缺少时抛出 `BackendUnavailableError`；绘图还需要 matplotlib 和 scipy，它们随 sourmash 一起安装。
- **输入：** `comparison_matrix` 是 `compare_signatures()` 写出的已有矩阵文件（旁边必须有它的 `.labels.txt` 文件）。`output_dir` 必要时会被创建。
- **输出：** 标准结果字典，`outputs` 列出所有新文件，只有退出码为 0 且写出了 PNG 时 `status` 才是 `completed`；另有 `image`，即第一个新 PNG 的路径（没有则为 `None`）。
- **示例代码：**

```python
import shutil, sys, tempfile
from pathlib import Path

from dnakit.microbiome import sketch_dna, compare_signatures, plot_comparison_matrix

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())
sketch_dna("a.fa", work / "a.sig", scaled=1, name="sigA first", sourmash_path=tool, verbose=False)
sketch_dna("b.fa", work / "b.sig", scaled=1, name="sigB second", sourmash_path=tool, verbose=False)

compare_signatures([work / "a.sig", work / "b.sig"], work / "cmp", sourmash_path=tool, verbose=False)
result = plot_comparison_matrix(work / "cmp", work / "plots", sourmash_path=tool, verbose=False)
print(result["status"], Path(result["image"]).suffix, Path(result["image"]).stat().st_size > 1000)
```

- **示例输出：**

```text
completed .png True
```

## 12) `dnakit.microbiome.prefetch_signatures()`

- **作用：** 选出与查询签名有重叠的数据库签名。
- **处理方法：** 运行 `sourmash prefetch <query> <databases> -k <ksize> --dna --threshold-bp <bp> -o <csv>`。prefetch 是快速的第一轮筛选，列出所有与查询共享哈希的数据库签名，不做 `gather` 那样的贪心分解。
- **API：** `dnakit.microbiome.prefetch_signatures(query: 'str', databases: 'Sequence[str]', output_csv: 'str', *, ksize: 'int' = 21, threshold_bp: 'int' = 0, sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `sourmash` 在 PATH 中（或通过 `sourmash_path` 指定其位置）。它由可选包 sourmash 4.9 及以上版本提供（extra `dnakit[bio-sourmash]`）；缺少时抛出 `BackendUnavailableError`。
- **输入：** `query` 是已有的签名文件，`databases` 是非空的已有签名集合列表。`ksize` 默认 21，`threshold_bp` 默认 0。`output_csv` 必须尚不存在。
- **输出：** 标准结果字典：`status`（sourmash 退出码为 0 且所有预期输出文件都存在时为 `completed`，否则为 `failed`）、`backend`（`sourmash`）、`command`、`return_code`、`outputs`（实际存在的输出文件）和 `log_tail`（sourmash 输出的最后 4000 个字符，其中含终端控制码）。 CSV 中每个有重叠的数据库签名占一行。
- **示例代码：**

```python
import shutil, sys, tempfile
from pathlib import Path
import csv

from dnakit.microbiome import sketch_dna, prefetch_signatures

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())
sketch_dna("a.fa", work / "a.sig", scaled=1, name="sigA first", sourmash_path=tool, verbose=False)
sketch_dna("b.fa", work / "b.sig", scaled=1, name="sigB second", sourmash_path=tool, verbose=False)

result = prefetch_signatures(work / "a.sig", [work / "b.sig"], work / "prefetch.csv", sourmash_path=tool, verbose=False)
rows = list(csv.DictReader(open(work / "prefetch.csv")))
print(result["status"], len(rows), rows[0]["match_name"], rows[0]["intersect_bp"])
```

- **示例输出：**

```text
completed 1 sigB second 15
```

## 13) `dnakit.microbiome.sketch_translated_dna()`

- **作用：** 把 DNA 按六个阅读框翻译，并对蛋白 k-mer 生成草图。
- **处理方法：** 运行 `sourmash sketch translate <fasta> -p k=<ksize>,scaled=<scaled> -o <signature>`。签名中保存的是蛋白草图；`k=5` 时 sourmash 对该蛋白草图报告的 ksize 为 5。
- **API：** `dnakit.microbiome.sketch_translated_dna(fasta: 'str', output_signature: 'str', *, ksize: 'int' = 10, scaled: 'int' = 200, sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `sourmash` 在 PATH 中（或通过 `sourmash_path` 指定其位置）。它由可选包 sourmash 4.9 及以上版本提供（extra `dnakit[bio-sourmash]`）；缺少时抛出 `BackendUnavailableError`。
- **输入：** `fasta` 是已有的 FASTA 文件，`ksize`（默认 10）和 `scaled`（默认 200）是正整数，否则抛出 `ConfigurationError`。`output_signature` 必须尚不存在。
- **输出：** 标准结果字典：`status`（sourmash 退出码为 0 且所有预期输出文件都存在时为 `completed`，否则为 `failed`）、`backend`（`sourmash`）、`command`、`return_code`、`outputs`（实际存在的输出文件）和 `log_tail`（sourmash 输出的最后 4000 个字符，其中含终端控制码）。
- **示例代码：**

```python
import shutil, sys, tempfile
from pathlib import Path
from sourmash import load_file_as_signatures

from dnakit.microbiome import sketch_dna, sketch_translated_dna

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())

result = sketch_translated_dna("a.fa", work / "protein.sig", ksize=5, scaled=1, sourmash_path=tool, verbose=False)
(sig,) = load_file_as_signatures(str(work / "protein.sig"))
print(result["status"], sig.minhash.moltype, sig.minhash.ksize, len(sig.minhash))
```

- **示例输出：**

```text
completed protein 5 92
```

## 14) `dnakit.microbiome.compare_signatures()`

- **作用：** 两两比较签名，并写出相似度（或包含度）矩阵。
- **处理方法：** 运行 `sourmash compare <signatures> -k <ksize> -o <matrix>`，`containment` 为真时再加 `--containment`。sourmash 把矩阵写成 NumPy 文件（可用 `numpy.load` 读取），并把样本顺序写在旁边的 `<matrix>.labels.txt` 中。
- **API：** `dnakit.microbiome.compare_signatures(signatures: 'Sequence[str]', output_matrix: 'str', *, ksize: 'int' = 21, containment: 'bool' = False, sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `sourmash` 在 PATH 中（或通过 `sourmash_path` 指定其位置）。它由可选包 sourmash 4.9 及以上版本提供（extra `dnakit[bio-sourmash]`）；缺少时抛出 `BackendUnavailableError`。
- **输入：** `signatures` 是至少两个已有签名文件的列表，其中含有 `ksize`（默认 21）的草图；少于两个会抛出 `ConfigurationError`。`output_matrix` 必须尚不存在。
- **输出：** 标准结果字典：`status`（sourmash 退出码为 0 且所有预期输出文件都存在时为 `completed`，否则为 `failed`）、`backend`（`sourmash`）、`command`、`return_code`、`outputs`（实际存在的输出文件）和 `log_tail`（sourmash 输出的最后 4000 个字符，其中含终端控制码）。 `outputs` 只含矩阵文件，不列出 labels 文件。
- **示例代码：**

```python
import shutil, sys, tempfile
from pathlib import Path
import numpy as np

from dnakit.microbiome import sketch_dna, compare_signatures

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())
sketch_dna("a.fa", work / "a.sig", scaled=1, name="sigA first", sourmash_path=tool, verbose=False)
sketch_dna("b.fa", work / "b.sig", scaled=1, name="sigB second", sourmash_path=tool, verbose=False)

result = compare_signatures([work / "a.sig", work / "b.sig"], work / "cmp", sourmash_path=tool, verbose=False)
matrix = np.load(work / "cmp")
print(result["status"], matrix.shape, matrix.round(3).tolist())
```

- **示例输出：**

```text
completed (2, 2) [[1.0, 0.214], [0.214, 1.0]]
```

## 15) `dnakit.microbiome.concatenate_signatures()`

- **作用：** 把签名连接成一个文件，但不合并它们。
- **处理方法：** 运行 `sourmash sig cat <signatures> -o <signature>`。每个签名在输出中仍是独立条目。
- **API：** `dnakit.microbiome.concatenate_signatures(signatures: 'Sequence[str]', output_signature: 'str', *, sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `sourmash` 在 PATH 中（或通过 `sourmash_path` 指定其位置）。它由可选包 sourmash 4.9 及以上版本提供（extra `dnakit[bio-sourmash]`）；缺少时抛出 `BackendUnavailableError`。
- **输入：** `signatures` 是非空的已有签名文件列表。`output_signature` 必须尚不存在。
- **输出：** 标准结果字典：`status`（sourmash 退出码为 0 且所有预期输出文件都存在时为 `completed`，否则为 `failed`）、`backend`（`sourmash`）、`command`、`return_code`、`outputs`（实际存在的输出文件）和 `log_tail`（sourmash 输出的最后 4000 个字符，其中含终端控制码）。
- **示例代码：**

```python
import shutil, sys, tempfile
from pathlib import Path
from sourmash import load_file_as_signatures

from dnakit.microbiome import sketch_dna, concatenate_signatures

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())
sketch_dna("a.fa", work / "a.sig", scaled=1, name="sigA first", sourmash_path=tool, verbose=False)
sketch_dna("b.fa", work / "b.sig", scaled=1, name="sigB second", sourmash_path=tool, verbose=False)

result = concatenate_signatures([work / "a.sig", work / "b.sig"], work / "both.sig", sourmash_path=tool, verbose=False)
sigs = load_file_as_signatures(str(work / "both.sig"))
print(result["status"], [sig.name for sig in sigs])
```

- **示例输出：**

```text
completed ['sigA first', 'sigB second']
```

## 16) `dnakit.microbiome.extract_signatures()`

- **作用：** 从集合中提取某一 k-mer 大小的签名。
- **处理方法：** 运行 `sourmash sig extract <collection> -k <ksize> -o <signature>`，`dna` 为假时用 `--protein`，否则用 `--dna`（默认）。
- **API：** `dnakit.microbiome.extract_signatures(collection: 'Sequence[str]', output_signature: 'str', *, ksize: 'int', dna: 'bool' = True, sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `sourmash` 在 PATH 中（或通过 `sourmash_path` 指定其位置）。它由可选包 sourmash 4.9 及以上版本提供（extra `dnakit[bio-sourmash]`）；缺少时抛出 `BackendUnavailableError`。
- **输入：** `collection` 是非空的已有签名文件列表，`ksize` 必填（仅限关键字参数），`dna` 选择核苷酸或蛋白草图。`output_signature` 必须尚不存在。
- **输出：** 标准结果字典：`status`（sourmash 退出码为 0 且所有预期输出文件都存在时为 `completed`，否则为 `failed`）、`backend`（`sourmash`）、`command`、`return_code`、`outputs`（实际存在的输出文件）和 `log_tail`（sourmash 输出的最后 4000 个字符，其中含终端控制码）。
- **示例代码：**

```python
import shutil, sys, tempfile
from pathlib import Path
from sourmash import load_file_as_signatures

from dnakit.microbiome import sketch_dna, extract_signatures

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())
sketch_dna("a.fa", work / "a.sig", scaled=1, name="sigA first", sourmash_path=tool, verbose=False)
sketch_dna("b.fa", work / "b.sig", scaled=1, name="sigB second", sourmash_path=tool, verbose=False)

sketch_dna("b.fa", work / "multi.sig", ksizes=[21, 31], scaled=1, sourmash_path=tool, verbose=False)
result = extract_signatures([work / "multi.sig"], work / "k31.sig", ksize=31, sourmash_path=tool, verbose=False)
(sig,) = load_file_as_signatures(str(work / "k31.sig"))
print(result["status"], sig.minhash.ksize, len(sig.minhash))
```

- **示例输出：**

```text
completed 31 35
```

## 17) `dnakit.microbiome.filter_signature_abundance()`

- **作用：** 保留丰度落在闭区间内的哈希。
- **处理方法：** 运行 `sourmash sig filter <signature> -o <output>`，并按需加 `-m <min_abundance>` 和 `-M <max_abundance>`。它只对记录了丰度的签名有意义；`sketch_dna()` 生成的签名不含丰度，所以示例用 sourmash 的 Python API 构建带丰度的签名。
- **API：** `dnakit.microbiome.filter_signature_abundance(signature: 'str', output_signature: 'str', *, min_abundance: 'int | None' = None, max_abundance: 'int | None' = None, sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `sourmash` 在 PATH 中（或通过 `sourmash_path` 指定其位置）。它由可选包 sourmash 4.9 及以上版本提供（extra `dnakit[bio-sourmash]`）；缺少时抛出 `BackendUnavailableError`。
- **输入：** `signature` 是已有的签名文件；`min_abundance` 和 `max_abundance` 至少要给一个，否则抛出 `ConfigurationError`。`output_signature` 必须尚不存在。
- **输出：** 标准结果字典：`status`（sourmash 退出码为 0 且所有预期输出文件都存在时为 `completed`，否则为 `failed`）、`backend`（`sourmash`）、`command`、`return_code`、`outputs`（实际存在的输出文件）和 `log_tail`（sourmash 输出的最后 4000 个字符，其中含终端控制码）。
- **示例代码：**

```python
import shutil, sys, tempfile
from pathlib import Path
from sourmash import MinHash, SourmashSignature, load_file_as_signatures
from sourmash.sourmash_args import SaveSignaturesToLocation

from dnakit.microbiome import sketch_dna, filter_signature_abundance

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())

minhash = MinHash(n=0, ksize=21, scaled=1, track_abundance=True)
minhash.add_sequence("GCTAAAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCG")
minhash.add_sequence("GCTAAAGACAATTACATAACATACACGTCAGCACGA")
with SaveSignaturesToLocation(str(work / "abund.sig")) as saver:
    saver.add(SourmashSignature(minhash, name="abund"))
result = filter_signature_abundance(work / "abund.sig", work / "kept.sig", min_abundance=2, sourmash_path=tool, verbose=False)
(sig,) = load_file_as_signatures(str(work / "kept.sig"))
print(result["status"], len(minhash), len(sig.minhash), sorted(set(sig.minhash.hashes.values())))
```

- **示例输出：**

```text
completed 40 16 [2]
```

## 18) `dnakit.microbiome.manifest_signature()`

- **作用：** 写出描述文件中每个签名的 CSV manifest。
- **处理方法：** 运行 `sourmash sig manifest <signature> -o <csv>`。CSV 以 `# SOURMASH-MANIFEST-VERSION` 注释行开头，随后是表头，列包括 `md5`、`ksize`、`moltype`、`num`、`scaled`、`n_hashes`、`with_abundance` 和 `name` 等。
- **API：** `dnakit.microbiome.manifest_signature(signature: 'str', output_csv: 'str', *, sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `sourmash` 在 PATH 中（或通过 `sourmash_path` 指定其位置）。它由可选包 sourmash 4.9 及以上版本提供（extra `dnakit[bio-sourmash]`）；缺少时抛出 `BackendUnavailableError`。
- **输入：** `signature` 是已有的签名文件。`output_csv` 必须尚不存在。
- **输出：** 标准结果字典：`status`（sourmash 退出码为 0 且所有预期输出文件都存在时为 `completed`，否则为 `failed`）、`backend`（`sourmash`）、`command`、`return_code`、`outputs`（实际存在的输出文件）和 `log_tail`（sourmash 输出的最后 4000 个字符，其中含终端控制码）。
- **示例代码：**

```python
import shutil, sys, tempfile
from pathlib import Path
import csv

from dnakit.microbiome import sketch_dna, manifest_signature

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())
sketch_dna("a.fa", work / "a.sig", scaled=1, name="sigA first", sourmash_path=tool, verbose=False)
sketch_dna("b.fa", work / "b.sig", scaled=1, name="sigB second", sourmash_path=tool, verbose=False)

result = manifest_signature(work / "a.sig", work / "manifest.csv", sourmash_path=tool, verbose=False)
lines = (work / "manifest.csv").read_text().splitlines()
rows = list(csv.DictReader(lines[1:]))
print(result["status"], len(rows), rows[0]["name"], rows[0]["ksize"], rows[0]["n_hashes"])
```

- **示例输出：**

```text
completed 1 sigA first 21 40
```

## 19) `dnakit.microbiome.split_signatures()`

- **作用：** 把签名文件拆成每个签名一个文件，并使用标准化的文件名。
- **处理方法：** 运行 `sourmash sig split <signatures> --outdir <dir>`。必要时 DNAKit 会创建该文件夹（必须是新建的或为空），`outputs` 是其中排序后的 `*.sig` 文件列表；只有退出码为 0 且至少有一个 `.sig` 文件时 `status` 才是 `completed`。
- **API：** `dnakit.microbiome.split_signatures(signatures: 'Sequence[str]', output_directory: 'str', *, sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `sourmash` 在 PATH 中（或通过 `sourmash_path` 指定其位置）。它由可选包 sourmash 4.9 及以上版本提供（extra `dnakit[bio-sourmash]`）；缺少时抛出 `BackendUnavailableError`。
- **输入：** `signatures` 是非空的已有签名文件列表。`output_directory` 必须不存在或为空，否则抛出 `ConfigurationError`。
- **输出：** 标准结果字典：`status`（sourmash 退出码为 0 且所有预期输出文件都存在时为 `completed`，否则为 `failed`）、`backend`（`sourmash`）、`command`、`return_code`、`outputs`（实际存在的输出文件）和 `log_tail`（sourmash 输出的最后 4000 个字符，其中含终端控制码）。 此处 `outputs` 列出写入该文件夹的每个 `.sig` 文件。
- **示例代码：**

```python
import shutil, sys, tempfile
from pathlib import Path

from dnakit.microbiome import sketch_dna, split_signatures

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())
sketch_dna("a.fa", work / "a.sig", scaled=1, name="sigA first", sourmash_path=tool, verbose=False)
sketch_dna("b.fa", work / "b.sig", scaled=1, name="sigB second", sourmash_path=tool, verbose=False)

result = split_signatures([work / "a.sig", work / "b.sig"], work / "parts", sourmash_path=tool, verbose=False)
print(result["status"], len(result["outputs"]), all(name.endswith(".sig") for name in result["outputs"]))
```

- **示例输出：**

```text
completed 2 True
```

## 20) `dnakit.microbiome.downsample_signature()`

- **作用：** 把签名的分辨率降低到更大的 scaled 值。
- **处理方法：** 运行 `sourmash sig downsample <signature> --scaled <n> -o <output>`。只保留小于 `2^64 / n` 的哈希，因此签名只能降采样，不能变得更精细。
- **API：** `dnakit.microbiome.downsample_signature(signature: 'str', output_signature: 'str', *, scaled: 'int', sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `sourmash` 在 PATH 中（或通过 `sourmash_path` 指定其位置）。它由可选包 sourmash 4.9 及以上版本提供（extra `dnakit[bio-sourmash]`）；缺少时抛出 `BackendUnavailableError`。
- **输入：** `signature` 是已有的签名文件，`scaled`（必填，仅限关键字参数）是比签名自身值更大的正整数。`output_signature` 必须尚不存在。
- **输出：** 标准结果字典：`status`（sourmash 退出码为 0 且所有预期输出文件都存在时为 `completed`，否则为 `failed`）、`backend`（`sourmash`）、`command`、`return_code`、`outputs`（实际存在的输出文件）和 `log_tail`（sourmash 输出的最后 4000 个字符，其中含终端控制码）。
- **示例代码：**

```python
import shutil, sys, tempfile
from pathlib import Path
from sourmash import load_file_as_signatures

from dnakit.microbiome import sketch_dna, downsample_signature

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())
sketch_dna("a.fa", work / "a.sig", scaled=1, name="sigA first", sourmash_path=tool, verbose=False)
sketch_dna("b.fa", work / "b.sig", scaled=1, name="sigB second", sourmash_path=tool, verbose=False)

result = downsample_signature(work / "a.sig", work / "a_scaled4.sig", scaled=4, sourmash_path=tool, verbose=False)
(sig,) = load_file_as_signatures(str(work / "a_scaled4.sig"))
print(result["status"], sig.minhash.scaled, len(sig.minhash))
```

- **示例输出：**

```text
completed 4 10
```

## 21) `dnakit.microbiome.flatten_signature()`

- **作用：** 去除签名中的丰度信息。
- **处理方法：** 运行 `sourmash sig flatten <signature> -o <output>`。哈希保留，计数被丢弃。
- **API：** `dnakit.microbiome.flatten_signature(signature: 'str', output_signature: 'str', *, sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `sourmash` 在 PATH 中（或通过 `sourmash_path` 指定其位置）。它由可选包 sourmash 4.9 及以上版本提供（extra `dnakit[bio-sourmash]`）；缺少时抛出 `BackendUnavailableError`。
- **输入：** `signature` 是已有的签名文件。`output_signature` 必须尚不存在。
- **输出：** 标准结果字典：`status`（sourmash 退出码为 0 且所有预期输出文件都存在时为 `completed`，否则为 `failed`）、`backend`（`sourmash`）、`command`、`return_code`、`outputs`（实际存在的输出文件）和 `log_tail`（sourmash 输出的最后 4000 个字符，其中含终端控制码）。
- **示例代码：**

```python
import shutil, sys, tempfile
from pathlib import Path
from sourmash import MinHash, SourmashSignature, load_file_as_signatures
from sourmash.sourmash_args import SaveSignaturesToLocation

from dnakit.microbiome import sketch_dna, flatten_signature

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())

minhash = MinHash(n=0, ksize=21, scaled=1, track_abundance=True)
minhash.add_sequence("GCTAAAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCG")
minhash.add_sequence("GCTAAAGACAATTACATAACATACACGTCAGCACGA")
with SaveSignaturesToLocation(str(work / "abund.sig")) as saver:
    saver.add(SourmashSignature(minhash, name="abund"))
result = flatten_signature(work / "abund.sig", work / "flat.sig", sourmash_path=tool, verbose=False)
(sig,) = load_file_as_signatures(str(work / "flat.sig"))
print(result["status"], minhash.track_abundance, sig.minhash.track_abundance, len(sig.minhash))
```

- **示例输出：**

```text
completed True False 40
```

## 22) `dnakit.microbiome.rename_signature()`

- **作用：** 给签名指定新名称。
- **处理方法：** 运行 `sourmash sig rename <signature> <新名称> -o <output>`。
- **API：** `dnakit.microbiome.rename_signature(signature: 'str', new_name: 'str', output_signature: 'str', *, sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `sourmash` 在 PATH 中（或通过 `sourmash_path` 指定其位置）。它由可选包 sourmash 4.9 及以上版本提供（extra `dnakit[bio-sourmash]`）；缺少时抛出 `BackendUnavailableError`。
- **输入：** `signature` 是已有的签名文件，`new_name` 是非空字符串（空白会抛出 `ConfigurationError`）。`output_signature` 必须尚不存在。
- **输出：** 标准结果字典：`status`（sourmash 退出码为 0 且所有预期输出文件都存在时为 `completed`，否则为 `failed`）、`backend`（`sourmash`）、`command`、`return_code`、`outputs`（实际存在的输出文件）和 `log_tail`（sourmash 输出的最后 4000 个字符，其中含终端控制码）。
- **示例代码：**

```python
import shutil, sys, tempfile
from pathlib import Path
from sourmash import load_file_as_signatures

from dnakit.microbiome import sketch_dna, rename_signature

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())
sketch_dna("a.fa", work / "a.sig", scaled=1, name="sigA first", sourmash_path=tool, verbose=False)
sketch_dna("b.fa", work / "b.sig", scaled=1, name="sigB second", sourmash_path=tool, verbose=False)

result = rename_signature(work / "a.sig", "renamed sample", work / "renamed.sig", sourmash_path=tool, verbose=False)
(sig,) = load_file_as_signatures(str(work / "renamed.sig"))
print(result["status"], sig.name)
```

- **示例输出：**

```text
completed renamed sample
```

## 23) `dnakit.microbiome.signature_file_info()`

- **作用：** 报告签名文件的文件级信息。
- **处理方法：** 运行 `sourmash sig fileinfo <signature>` 并返回其文本：文件类型、是否有 manifest、签名数、哈希总数，以及按 k-mer 大小和分子类型给出的汇总。
- **API：** `dnakit.microbiome.signature_file_info(signature: 'str', *, sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `sourmash` 在 PATH 中（或通过 `sourmash_path` 指定其位置）。它由可选包 sourmash 4.9 及以上版本提供（extra `dnakit[bio-sourmash]`）；缺少时抛出 `BackendUnavailableError`。
- **输入：** `signature` 是已有的签名文件。
- **输出：** `status`（退出码为 0 时为 `completed`）、`backend`、`return_code` 和 `description`，后者是 sourmash 文本输出的最后 4000 个字符（可能含终端控制码）。
- **示例代码：**

```python
import shutil, sys, tempfile
from pathlib import Path

from dnakit.microbiome import sketch_dna, signature_file_info

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())
sketch_dna("a.fa", work / "a.sig", scaled=1, name="sigA first", sourmash_path=tool, verbose=False)
sketch_dna("b.fa", work / "b.sig", scaled=1, name="sigB second", sourmash_path=tool, verbose=False)

result = signature_file_info(work / "a.sig", sourmash_path=tool, verbose=False)
print(result["status"])
print([line.strip() for line in result["description"].splitlines() if "total hashes" in line or "sketches with" in line])
```

- **示例输出：**

```text
completed
['total hashes: 40', '1 sketches with DNA, k=21, scaled=1                40 total hashes']
```

## 24) `dnakit.microbiome.intersect_signatures()`

- **作用：** 只保留所有签名共有的 k-mer。
- **处理方法：** 运行 `sourmash sig intersect <signatures> -o <output>`。
- **API：** `dnakit.microbiome.intersect_signatures(signatures: 'Sequence[str]', output_signature: 'str', *, sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `sourmash` 在 PATH 中（或通过 `sourmash_path` 指定其位置）。它由可选包 sourmash 4.9 及以上版本提供（extra `dnakit[bio-sourmash]`）；缺少时抛出 `BackendUnavailableError`。
- **输入：** `signatures` 是至少两个已有签名文件的列表，k-mer 大小和 scaled 值相同；少于两个会抛出 `ConfigurationError`。`output_signature` 必须尚不存在。
- **输出：** 标准结果字典：`status`（sourmash 退出码为 0 且所有预期输出文件都存在时为 `completed`，否则为 `failed`）、`backend`（`sourmash`）、`command`、`return_code`、`outputs`（实际存在的输出文件）和 `log_tail`（sourmash 输出的最后 4000 个字符，其中含终端控制码）。
- **示例代码：**

```python
import shutil, sys, tempfile
from pathlib import Path
from sourmash import load_file_as_signatures

from dnakit.microbiome import sketch_dna, intersect_signatures

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())
sketch_dna("a.fa", work / "a.sig", scaled=1, name="sigA first", sourmash_path=tool, verbose=False)
sketch_dna("b.fa", work / "b.sig", scaled=1, name="sigB second", sourmash_path=tool, verbose=False)

result = intersect_signatures([work / "a.sig", work / "b.sig"], work / "shared.sig", sourmash_path=tool, verbose=False)
(sig,) = load_file_as_signatures(str(work / "shared.sig"))
print(result["status"], len(sig.minhash))
```

- **示例输出：**

```text
completed 15
```

## 25) `dnakit.microbiome.merge_signatures()`

- **作用：** 把同一 k-mer 大小的多个签名合并成一个。
- **处理方法：** 运行 `sourmash sig merge <signatures> -o <output>`，`flatten` 为真时再加 `--flatten`（去掉所有丰度，这样带丰度和不带丰度的签名可以混合合并）。结果是哈希的并集。
- **API：** `dnakit.microbiome.merge_signatures(signatures: 'Sequence[str]', output_signature: 'str', *, flatten: 'bool' = False, sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `sourmash` 在 PATH 中（或通过 `sourmash_path` 指定其位置）。它由可选包 sourmash 4.9 及以上版本提供（extra `dnakit[bio-sourmash]`）；缺少时抛出 `BackendUnavailableError`。
- **输入：** `signatures` 是至少两个已有签名文件的列表，k-mer 大小相同；少于两个会抛出 `ConfigurationError`。`output_signature` 必须尚不存在。
- **输出：** 标准结果字典：`status`（sourmash 退出码为 0 且所有预期输出文件都存在时为 `completed`，否则为 `failed`）、`backend`（`sourmash`）、`command`、`return_code`、`outputs`（实际存在的输出文件）和 `log_tail`（sourmash 输出的最后 4000 个字符，其中含终端控制码）。
- **示例代码：**

```python
import shutil, sys, tempfile
from pathlib import Path
from sourmash import load_file_as_signatures

from dnakit.microbiome import sketch_dna, merge_signatures

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())
sketch_dna("a.fa", work / "a.sig", scaled=1, name="sigA first", sourmash_path=tool, verbose=False)
sketch_dna("b.fa", work / "b.sig", scaled=1, name="sigB second", sourmash_path=tool, verbose=False)

result = merge_signatures([work / "a.sig", work / "b.sig"], work / "merged.sig", sourmash_path=tool, verbose=False)
(sig,) = load_file_as_signatures(str(work / "merged.sig"))
print(result["status"], len(sig.minhash))
```

- **示例输出：**

```text
completed 70
```

## 26) `dnakit.microbiome.subtract_signatures()`

- **作用：** 从一个签名中去掉其他签名的 k-mer。
- **处理方法：** 运行 `sourmash sig subtract <signature> <others> -o <output>`。
- **API：** `dnakit.microbiome.subtract_signatures(signature: 'str', others: 'Sequence[str]', output_signature: 'str', *, sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `sourmash` 在 PATH 中（或通过 `sourmash_path` 指定其位置）。它由可选包 sourmash 4.9 及以上版本提供（extra `dnakit[bio-sourmash]`）；缺少时抛出 `BackendUnavailableError`。
- **输入：** `signature` 是已有的签名文件，`others` 是非空的已有签名文件列表（为空会抛出 `ConfigurationError`）。`output_signature` 必须尚不存在。
- **输出：** 标准结果字典：`status`（sourmash 退出码为 0 且所有预期输出文件都存在时为 `completed`，否则为 `failed`）、`backend`（`sourmash`）、`command`、`return_code`、`outputs`（实际存在的输出文件）和 `log_tail`（sourmash 输出的最后 4000 个字符，其中含终端控制码）。
- **示例代码：**

```python
import shutil, sys, tempfile
from pathlib import Path
from sourmash import load_file_as_signatures

from dnakit.microbiome import sketch_dna, subtract_signatures

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())
sketch_dna("a.fa", work / "a.sig", scaled=1, name="sigA first", sourmash_path=tool, verbose=False)
sketch_dna("b.fa", work / "b.sig", scaled=1, name="sigB second", sourmash_path=tool, verbose=False)

result = subtract_signatures(work / "a.sig", [work / "b.sig"], work / "only_a.sig", sourmash_path=tool, verbose=False)
(sig,) = load_file_as_signatures(str(work / "only_a.sig"))
print(result["status"], len(sig.minhash))
```

- **示例输出：**

```text
completed 25
```

## 27) `dnakit.microbiome.describe_signature()`

- **作用：** 描述签名文件的内容：名称、k-mer 大小、scaled 值和大小。
- **处理方法：** 运行 `sourmash sig describe <signature>` 并返回其文本。每个签名会显示名称、源文件、md5、一行含 `k=`、分子类型、`num`、`scaled`、种子和丰度标志的信息，以及哈希数量。
- **API：** `dnakit.microbiome.describe_signature(signature: 'str', *, sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `sourmash` 在 PATH 中（或通过 `sourmash_path` 指定其位置）。它由可选包 sourmash 4.9 及以上版本提供（extra `dnakit[bio-sourmash]`）；缺少时抛出 `BackendUnavailableError`。
- **输入：** `signature` 是已有的签名文件。
- **输出：** `status`（退出码为 0 时为 `completed`）、`backend`、`return_code` 和 `description`，后者是 sourmash 文本输出的最后 4000 个字符（可能含终端控制码）。
- **示例代码：**

```python
import shutil, sys, tempfile
from pathlib import Path

from dnakit.microbiome import sketch_dna, describe_signature

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())
sketch_dna("a.fa", work / "a.sig", scaled=1, name="sigA first", sourmash_path=tool, verbose=False)
sketch_dna("b.fa", work / "b.sig", scaled=1, name="sigB second", sourmash_path=tool, verbose=False)

result = describe_signature(work / "a.sig", sourmash_path=tool, verbose=False)
print(result["status"])
print([line for line in result["description"].splitlines() if line.startswith(("signature:", "k=", "size:"))])
```

- **示例输出：**

```text
completed
['signature: sigA first', 'k=21 molecule=DNA num=0 scaled=1 seed=42 track_abundance=0', 'size: 40']
```

## 28) `dnakit.microbiome.gather_signatures()`

- **作用：** 把查询签名分解为数据库中最佳匹配的签名。
- **处理方法：** 运行 `sourmash gather <query> <database> -k <ksize> -o <csv> --threshold-bp <bp>`。gather 反复选取剩余重叠最大的数据库签名，并从查询中去掉其哈希，因此各行是互不重叠的匹配。`threshold_bp` 是报告的最小重叠；默认值 50000 是 sourmash 的默认值；对很小的测试签名，没有任何匹配能达到该值，sourmash 不会写出 CSV，`status` 为 `failed`。因此示例中设为 0。
- **API：** `dnakit.microbiome.gather_signatures(query: 'str', database: 'str', output_csv: 'str', *, ksize: 'int' = 21, threshold_bp: 'int' = 50000, sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `sourmash` 在 PATH 中（或通过 `sourmash_path` 指定其位置）。它由可选包 sourmash 4.9 及以上版本提供（extra `dnakit[bio-sourmash]`）；缺少时抛出 `BackendUnavailableError`。
- **输入：** `query` 是已有的签名文件，`database` 是已有的签名集合，`ksize` 默认 21，`threshold_bp` 不能为负（否则抛出 `ConfigurationError`）。`output_csv` 必须尚不存在。
- **输出：** 标准结果字典：`status`（sourmash 退出码为 0 且所有预期输出文件都存在时为 `completed`，否则为 `failed`）、`backend`（`sourmash`）、`command`、`return_code`、`outputs`（实际存在的输出文件）和 `log_tail`（sourmash 输出的最后 4000 个字符，其中含终端控制码）。 CSV 中每个 gather 匹配占一行，列包括 `intersect_bp`、`f_orig_query`、`f_match` 和 `name` 等。
- **示例代码：**

```python
import shutil, sys, tempfile
from pathlib import Path
import csv

from dnakit.microbiome import sketch_dna, gather_signatures

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())
sketch_dna("a.fa", work / "a.sig", scaled=1, name="sigA first", sourmash_path=tool, verbose=False)
sketch_dna("b.fa", work / "b.sig", scaled=1, name="sigB second", sourmash_path=tool, verbose=False)

result = gather_signatures(work / "a.sig", work / "b.sig", work / "gather.csv", threshold_bp=0, sourmash_path=tool, verbose=False)
rows = list(csv.DictReader(open(work / "gather.csv")))
print(result["status"], len(rows), rows[0]["name"], rows[0]["intersect_bp"], round(float(rows[0]["f_orig_query"]), 3))
```

- **示例输出：**

```text
completed 1 sigB second 15 0.375
```

## 29) `dnakit.microbiome.index_signatures()`

- **作用：** 用签名文件构建可搜索的签名数据库（SBT）。
- **处理方法：** 运行 `sourmash index <database> <signatures> -k <ksize>`。sourmash 会在名称后追加 `.sbt.zip`，因此 DNAKit 检查并报告 `<output_database>.sbt.zip`。
- **API：** `dnakit.microbiome.index_signatures(signatures: 'Sequence[str]', output_database: 'str', *, ksize: 'int' = 21, sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `sourmash` 在 PATH 中（或通过 `sourmash_path` 指定其位置）。它由可选包 sourmash 4.9 及以上版本提供（extra `dnakit[bio-sourmash]`）；缺少时抛出 `BackendUnavailableError`。
- **输入：** `signatures` 是非空的已有签名文件列表，`ksize` 默认 21。`output_database` 和 `<output_database>.sbt.zip` 都不能已存在。
- **输出：** 标准结果字典：`status`（sourmash 退出码为 0 且所有预期输出文件都存在时为 `completed`，否则为 `failed`）、`backend`（`sourmash`）、`command`、`return_code`、`outputs`（实际存在的输出文件）和 `log_tail`（sourmash 输出的最后 4000 个字符，其中含终端控制码）。 另有 `database`，即 `.sbt.zip` 文件的路径。
- **示例代码：**

```python
import shutil, sys, tempfile
from pathlib import Path

from dnakit.microbiome import sketch_dna, index_signatures

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())
sketch_dna("a.fa", work / "a.sig", scaled=1, name="sigA first", sourmash_path=tool, verbose=False)
sketch_dna("b.fa", work / "b.sig", scaled=1, name="sigB second", sourmash_path=tool, verbose=False)

result = index_signatures([work / "a.sig", work / "b.sig"], work / "db", sourmash_path=tool, verbose=False)
print(result["status"], Path(result["database"]).name, len(result["outputs"]))
```

- **示例输出：**

```text
completed db.sbt.zip 1
```

## 30) `dnakit.microbiome.search_signatures()`

- **作用：** 在签名数据库中搜索与查询签名匹配的签名。
- **处理方法：** 运行 `sourmash search <query> <database> -k <ksize> --threshold <t> -o <csv>`，`containment` 为真时再加 `--containment`。匹配按 Jaccard 相似度排序，若要求则按查询的包含度排序。
- **API：** `dnakit.microbiome.search_signatures(query: 'str', database: 'str', output_csv: 'str', *, ksize: 'int' = 21, threshold: 'float' = 0.08, containment: 'bool' = False, sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `sourmash` 在 PATH 中（或通过 `sourmash_path` 指定其位置）。它由可选包 sourmash 4.9 及以上版本提供（extra `dnakit[bio-sourmash]`）；缺少时抛出 `BackendUnavailableError`。
- **输入：** `query` 是已有的签名文件，`database` 是已有的签名集合，`ksize` 默认 21，`threshold` 在 0 到 1 之间（默认 0.08，否则抛出 `ConfigurationError`）。`output_csv` 必须尚不存在。
- **输出：** 标准结果字典：`status`（sourmash 退出码为 0 且所有预期输出文件都存在时为 `completed`，否则为 `failed`）、`backend`（`sourmash`）、`command`、`return_code`、`outputs`（实际存在的输出文件）和 `log_tail`（sourmash 输出的最后 4000 个字符，其中含终端控制码）。 CSV 中每个匹配占一行，列包括 `similarity`、`md5`、`name` 和 `filename` 等。
- **示例代码：**

```python
import shutil, sys, tempfile
from pathlib import Path
import csv

from dnakit.microbiome import sketch_dna, search_signatures

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())
sketch_dna("a.fa", work / "a.sig", scaled=1, name="sigA first", sourmash_path=tool, verbose=False)
sketch_dna("b.fa", work / "b.sig", scaled=1, name="sigB second", sourmash_path=tool, verbose=False)

result = search_signatures(work / "a.sig", work / "b.sig", work / "search.csv", threshold=0.1, sourmash_path=tool, verbose=False)
rows = list(csv.DictReader(open(work / "search.csv")))
print(result["status"], len(rows), rows[0]["name"], round(float(rows[0]["similarity"]), 3))
```

- **示例输出：**

```text
completed 1 sigB second 0.214
```

## 31) `dnakit.microbiome.sketch_dna()`

- **作用：** 对 FASTA 文件计算带 scaled 的 DNA MinHash 签名。
- **处理方法：** 运行 `sourmash sketch dna <fasta> -o <signature> -p k=<k1>,k=<k2>,...,scaled=<n>`，给出 `name` 时再加 `--name <name>`。文件中所有记录对每个 k-mer 大小汇入一个签名，不记录丰度。
- **API：** `dnakit.microbiome.sketch_dna(fasta: 'str', output_signature: 'str', *, ksizes: 'Sequence[int]' = (21,), scaled: 'int' = 1000, name: 'str | None' = None, sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `sourmash` 在 PATH 中（或通过 `sourmash_path` 指定其位置）。它由可选包 sourmash 4.9 及以上版本提供（extra `dnakit[bio-sourmash]`）；缺少时抛出 `BackendUnavailableError`。
- **输入：** `fasta` 是已有的 FASTA 文件；`ksizes` 是非空的正整数序列（默认 `(21,)`），`scaled` 是正整数（默认 1000），否则抛出 `ConfigurationError`。`scaled` 低于 100 时 sourmash 会警告但继续运行，所以示例中使用 1。`output_signature` 必须尚不存在。
- **输出：** 标准结果字典：`status`（sourmash 退出码为 0 且所有预期输出文件都存在时为 `completed`，否则为 `failed`）、`backend`（`sourmash`）、`command`、`return_code`、`outputs`（实际存在的输出文件）和 `log_tail`（sourmash 输出的最后 4000 个字符，其中含终端控制码）。
- **示例代码：**

```python
import shutil, sys, tempfile
from pathlib import Path
from sourmash import load_file_as_signatures

from dnakit.microbiome import sketch_dna

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())

result = sketch_dna("a.fa", work / "a.sig", scaled=1, name="sigA first", sourmash_path=tool, verbose=False)
(sig,) = load_file_as_signatures(str(work / "a.sig"))
print(result["status"], result["return_code"], len(sig.minhash), sig.name)
```

- **示例输出：**

```text
completed 0 40 sigA first
```
