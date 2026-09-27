# 文件读写

读取使用 `dnakit.read()`，写出使用 `dnakit.write()`。

## 1) 序列格式 {#1-io-001}

- **作用：** 读取和写出 FASTA、FASTQ、GenBank 文件，将序列及可用的 ID、质量值和注释转换为 DNAKit 记录，便于继续分析或格式转换。
- **计算方法：** 按文件格式逐条解析记录，保留顺序；mode="stream" 惰性读取，写出时按目标格式序列化。
- **API：** `dnakit.read(source[必须], format[可选], config[可选], mode[可选])`、`dnakit.write(records[必须], target[必须], format[可选], config[可选])`；普通读取设 `mode="dna"`，大文件设 `mode="stream"`。
- **输入：** `mode="dna"` 无论文件含一条还是多条记录都返回 `DNA`；`mode="stream"` 返回单次消费的数据流。写出可直接传 `DNA`，也可直接传流式记录或其他记录迭代器，不需要“大文件写出”专用方法。
- **输出：** `read` 返回 `DNA` 或流式 `RecordSource`；`write` 返回含写入记录数的 `WriteResult`。
- **示例代码：**

```python
from io import StringIO

from dnakit import read, write

source = StringIO(">seq1 demo\nACGTN\n")
records = read(source, format="fasta", mode="dna")

output = StringIO()
written = write(records, output, format="fasta")
print(type(records).__name__, records.ids, records[0].symbols)
print(written.record_count)
```

- **示例结果：**

```text
DNA ('seq1',) ACGTN
1
```

<span id="5-io-005"></span>[FAQ：补充说明](../../faq/tutorial_basics.md#large-files)

## 2) 注释格式 {#2-io-002}

- **作用：** 读取和写出 GFF3、BED、AGP 注释文件，并统一不同格式的坐标规则，生成可用于区域查询、序列切分和注释转换的结构化数据。
- **计算方法：** 按 GFF3/BED/AGP 格式规则解析字段，并将坐标转换为内部区间。
- **API：** `dnakit.io.read_gff3(source[必须], max_records[可选], max_line_length[可选], max_header_lines[可选])`、`dnakit.io.write_gff3(document[必须], target[必须], overwrite[可选], create_parents[可选])`、`dnakit.io.read_bed(source[必须], max_records[可选], max_line_length[可选], max_header_lines[可选])`、`dnakit.io.write_bed(document[必须], target[必须], overwrite[可选], create_parents[可选])`、`dnakit.io.read_agp(source[必须], max_records[可选], max_line_length[可选], max_header_lines[可选])`、`dnakit.io.write_agp(document[必须], target[必须], overwrite[可选], create_parents[可选])`。
- **输入：** 读取时必填注释文件路径或文本流，可选记录数、行长与 header 上限；写出时传入对应 document/entry 和目标路径或文本流。
- **输出：** `read_gff3` 返回 `AnnotationDocument`；`write_gff3` 返回写出的特征数量。
- **示例代码：**

```python
from io import StringIO

from dnakit.io import read_gff3, write_gff3

source = StringIO(
    "##gff-version 3\n"
    "chr1\ttest\tgene\t2\t4\t.\t+\t.\tID=g1\n"
)
document = read_gff3(source)

output = StringIO()
count = write_gff3(document, output)
print(document.entries[0].feature.location)
print(count)
```

- **示例结果：**

```text
Interval(start=1, end=4)
1
```

## 3) 表格格式 {#3-io-003}

- **作用：** 按指定表结构读写 CSV、TSV、JSON 和 Parquet，在 DNA 记录、metadata 与通用表格之间进行可复现转换。
- **计算方法：** 按显式 schema 校验列、类型和缺失值，再解析或序列化表格。
- **API：** `dnakit.io.TableSchema(columns[必须], schema_version[可选], column_types[可选], nullable[可选])`、`dnakit.io.export_table(rows[必须], target[必须], format[必须], schema[必须], overwrite[可选], max_rows[可选], max_columns[可选], max_cell_characters[可选], max_output_bytes[可选], null_value[可选], parquet_compression[可选])`、`dnakit.io.read_table(source[必须], format[必须], schema[必须], missing_values[可选], max_rows[可选], max_columns[可选], max_cell_characters[可选], max_file_bytes[可选], max_decoded_bytes[可选])`、`dnakit.io.export_result(result[必须], target[必须], format[可选], overwrite[可选])`、`dnakit.io.parquet_backend_status()`。
- **输入：** `export_table()`/`read_table()` 使用二维行数据或目标路径、显式 `format` 和 `TableSchema`，可选列类型、nullable 字段、缺失值标记、压缩方式和资源上限；`export_result()` 接收带 `to_dict()` 的单个结果对象和目标路径，默认 JSON 并自动推导单行 schema。
- **输出：** `export_table` 返回写出行数；`read_table` 返回表格行和 schema 信息。
- **示例代码：**

```python
from pathlib import Path
from tempfile import TemporaryDirectory

from dnakit.io import TableSchema, export_table, read_table

schema = TableSchema(
    ("id", "score"),
    column_types={"id": "string", "score": "number"},
    nullable=(),
)
with TemporaryDirectory() as directory:
    path = Path(directory) / "scores.json"
    exported = export_table(
        [{"id": "a", "score": 0.75}],
        path,
        format="json",
        schema=schema,
    )
    loaded = read_table(path, format="json", schema=schema)
    print(exported.row_count, loaded.rows[0]["score"])
```

- **示例结果：**

```text
1 0.75
```

## 4) 压缩文件 {#4-io-004}

- **作用：** 直接读写 gzip 压缩的序列或表格文件，无需调用方先解压，并保留与普通文件相同的解析结果。
- **计算方法：** 按配置识别 gzip 压缩，在数据流外层解压或压缩后复用普通格式读写。
- **API：** `dnakit.read(source[必须], format[可选], config[可选], mode[可选])`、`dnakit.write(records[必须], target[必须], format[可选], config[可选])`；旧的 `read_one()` 仅为兼容入口。
- **输入：** 必填压缩文件路径或支持二进制的流；用 `compression="auto"|"none"|"gzip"` 控制压缩，可选写出压缩级别。
- **输出：** 读取返回 `DNA` 或 `RecordSource`；写入返回文件写入结果，压缩格式由路径和参数决定。
- **示例代码：**

```python
from pathlib import Path
from tempfile import TemporaryDirectory

from dnakit import DNA, ReadConfig, WriteConfig, read, write

with TemporaryDirectory() as directory:
    path = Path(directory) / "seq.fa.gz"
    record = DNA("ACGT", id="seq1")
    write(
        record,
        path,
        format="fasta",
        config=WriteConfig(compression="gzip"),
    )
    loaded = read(
        path,
        format="fasta",
        config=ReadConfig(compression="auto"),
        mode="dna",
    )
    print(loaded.id, loaded.symbols)
```

- **示例结果：**

```text
seq1 ACGT
```

[FAQ：补充说明](../../faq/tutorial_basics.md#compatibility-02_io_data)

## 5) 读写多序列比对文件 {#bio-2}

- **作用：** 在 FASTA、Clustal、Stockholm 等格式之间转换比对。
- **计算方法：** 使用 Biopython 解码与编码单个比对，保留序列 ID。
- **API：** `dnakit.bioinformatics.alignment.read_alignment(path, format='fasta')`; `dnakit.bioinformatics.alignment.write_alignment(alignment, path, format='fasta')`
- **输入：** 比对对象或已有比对文件；需要 `bio-formats`。示例自行创建输入。
- **输出：** 读取返回 `MultipleAlignment`，写入生成所选格式的比对文件。
- **示例代码：**

```python
from tempfile import TemporaryDirectory
from pathlib import Path
from dnakit.bioinformatics.alignment import MultipleAlignment, read_alignment, write_alignment

with TemporaryDirectory() as folder:
    msa = MultipleAlignment.from_sequences({"a": "ACGT", "b": "AC-T"})
    path = write_alignment(msa, Path(folder) / "aligned.fa")
    restored = read_alignment(path)
    print(restored.ids, restored.sequences)
```

- **示例结果：**

```text
('a', 'b') ('ACGT', 'AC-T')
```

[依赖与实现说明](../../faq/bioinformatics.md#alignment)。

## 6) 读写进化树文件 {#bio-7-io}

- **作用：** 将进化树保存为 Newick 等文件，并读取树结构。
- **计算方法：** 使用 Biopython 编解码树，保留枝长与叶节点。
- **API：** `dnakit.bioinformatics.phylogeny.read_tree(path, format='newick')`; `dnakit.bioinformatics.phylogeny.write_tree(tree, path, format='newick')`
- **输入：** 树对象或已有树文件；需要 `bio-formats`。示例自行创建树。
- **输出：** 读取返回进化树 `TreeNode`，写入生成所选格式的树文件。
- **示例代码：**

```python
from tempfile import TemporaryDirectory
from pathlib import Path
from dnakit.bioinformatics.phylogeny import TreeNode, read_tree, write_tree

tree = TreeNode(children=(TreeNode("a", 1), TreeNode("b", 1)))
with TemporaryDirectory() as folder:
    path = write_tree(tree, Path(folder) / "tree.nwk")
    restored = read_tree(path)
    print(restored.tips)
```

- **示例结果：**

```text
('a', 'b')
```

[依赖与实现说明](../../faq/bioinformatics.md#alignment)。

## 7) 读写 FASTA/FASTQ 并转换质量编码 {#bio-9}

- **作用：** 读取测序记录、保存序列并显式转换 Phred+64/33。
- **计算方法：** 流式解析 FASTA/FASTQ（含 gzip），质量以整数 Phred 保存。
- **API：** `dnakit.bioinformatics.reads.Read(id, sequence, qualities=None)`; `dnakit.bioinformatics.reads.read_fastx(path, format='fastq', phred_offset=33)`; `dnakit.bioinformatics.reads.write_fastx(reads, path, format='fastq', phred_offset=33)`
- **输入：** 文件路径或 `Read` 记录；FASTA 转 FASTQ 需显式提供 `default_quality`。
- **输出：** 读取返回 `Read` 记录迭代器，转换返回更新质量分数的记录。
- **示例代码：**

```python
from tempfile import TemporaryDirectory
from pathlib import Path
from dnakit.bioinformatics.reads import Read, read_fastx, write_fastx

with TemporaryDirectory() as folder:
    source = Path(folder) / "phred64.fq"
    target = Path(folder) / "phred33.fq"
    write_fastx([Read("r1", "ACGT", (30, 31, 32, 33))], source, phred_offset=64)
    count = write_fastx(read_fastx(source, phred_offset=64), target)
    print(count, next(read_fastx(target)).qualities)
```

- **示例结果：**

```text
1 (30, 31, 32, 33)
```

[依赖与实现说明](../../faq/bioinformatics.md#reads)。

## 8) 读取 JASPAR、MEME 或 TRANSFAC motif {#bio-50}

- **作用：** 将已有 motif 文件转成可扫描的原生模型。
- **计算方法：** 使用 Biopython 解码 counts，再生成 `MotifModel`。
- **API：** `dnakit.bioinformatics.motifs.read_motifs(path, format='jaspar', pseudocount=0.5)`
- **输入：** motif 文件；需要 `bio-formats` 和 `bio-analysis`。示例自行创建小型 JASPAR 文件。
- **输出：** 返回 `MotifModel` 元组，每项包含 motif 名称和概率矩阵。
- **示例代码：**

```python
from tempfile import TemporaryDirectory
from pathlib import Path
from dnakit.bioinformatics.motifs import read_motifs

with TemporaryDirectory() as folder:
    path = Path(folder) / "motif.jaspar"
    path.write_text(">MA0000.1 example\nA [ 3 0 0 ]\nC [ 0 3 0 ]\nG [ 0 0 3 ]\nT [ 0 0 0 ]\n")
    models = read_motifs(path)
    print(len(models), models[0].consensus)
```

- **示例结果：**

```text
1 ACG
```

[依赖与实现说明](../../faq/bioinformatics.md#optimization)。

## 9) 读取 Sanger ABI 碱基、质量与信号 {#bio-55-io}

- **作用：** 读取测序信号和碱基，并得到质量裁剪后的序列。
- **计算方法：** 使用 Biopython 解码 ABI，再通过 Mott 算法裁剪低质量两端。
- **API：** `dnakit.bioinformatics.cloning.read_sanger(path, error_probability_cutoff=0.05)`
- **输入：** 已有 `read.ab1` 文件；需要 `bio-formats`。
- **输出：** 返回字典，包含修剪区间、修剪后的序列、质量值和信号。
- **示例代码：**

```python
from dnakit.bioinformatics.cloning import read_sanger

trace = read_sanger("read.ab1")
print(trace["trim_interval"], trace["trimmed_sequence"])
```

- **示例结果：**

`trace` 含 `sequence`、`qualities`、`trim_interval`、`trimmed_sequence`、`channels` 和 `peak_positions`；打印裁剪区间与裁剪后序列。需要真实 ABI，此处未执行。

[依赖与实现说明](../../faq/bioinformatics.md#optimization)。

## 10) 读取 SnapGene 或 GenBank 构建文件 {#bio-57}

- **作用：** 读取构建序列及其 feature 注释。
- **计算方法：** 默认严格使用 Biopython 解码；显式 `permissive=True` 时使用 pydna。
- **API：** `dnakit.bioinformatics.cloning.read_constructs(path, format='snapgene', permissive=False)`
- **输入：** 已有 `.dna` 或 GenBank 文件；严格模式需 `bio-formats`，宽容模式需 `bio-assembly`。
- **输出：** 返回构建记录元组，包含序列和注释。
- **示例代码：**

```python
from dnakit.bioinformatics.cloning import read_constructs

records = read_constructs("construct.gb", format="genbank")
for record in records:
    print(record.id, len(record.seq), len(record.features))
```

- **示例结果：**

返回记录元组；循环输出构建 ID、碱基长度和 feature 数量。需要实际构建文件，此处未执行。

[依赖与实现说明](../../faq/bioinformatics.md#optimization)。

## 便捷行数据与结果保存 {#streaming-row-artifacts}

普通结果文件可以直接使用 `dnakit.io` 的以下接口，无需在脚本中编写 CSV writer、JSON 编码或文件打开逻辑。原有 `read_table`／`export_table` 继续用于需要显式类型、schema 和资源上限的表格。

| API | 用途 |
| --- | --- |
| `iter_rows(source, batch_size=4096)` | 按行读取 CSV、TSV、JSONL、Parquet；JSON 数组会整体加载 |
| `write_rows(rows, target, columns=None, overwrite=False, batch_size=4096)` | 将行迭代器写到一个路径或多个路径，返回行数 |
| `RowWriter(target, columns=None, overwrite=False, batch_size=4096)` | 在 `with` 中通过 `.write(row)` 逐条写入，适合按组分流 |
| `save_json(value, target, overwrite=False)`／`load_json(source)` | 保存、读取嵌套字典或数组；保存也接受有 `to_dict()` 的结果对象 |
| `save_arrays(arrays, target, compressed=True, overwrite=False)`／`load_arrays(source)` | 保存、读取命名数组组成的 NPZ；不使用 pickle，禁止 object 数组 |
| `file_digest(source, algorithm="sha256")` | 分块计算文件校验值，默认 SHA256 |

写出格式由目标后缀决定：CSV、TSV、JSON 数组、JSONL、Parquet 或 FASTA；文本格式支持 `.gz`。多目标写出只消费一次输入迭代器。列顺序默认取第一行，也可通过 `columns` 指定；后续行必须具有相同的键。空输入须提供列名。CSV／TSV 的 `None` 写为空单元格，读入时保持字符串；JSON、JSONL 和 Parquet 保留格式本身的类型。FASTA 仅使用 `id` 和 `sequence`，其他字段需同时保存到表格中。FASTA 读取仍使用 `dnakit.read()`。

```python
from dnakit.io import iter_rows, write_rows, save_json

rows = ({"id": f"seq{i}", "sequence": sequence, "label": i}
        for i, sequence in enumerate(["ACGT", "NNRY"]))
count = write_rows(rows, ["results/sequences.csv", "results/sequences.fasta"], overwrite=True)
save_json({"count": count, "parameters": {"source": "demo"}}, "results/run.json", overwrite=True)
print(list(iter_rows("results/sequences.csv")))
```

写出接口自动创建父目录，默认拒绝覆盖；只有成功完成后才替换目标文件。多个文件各自原子替换，不构成跨文件事务。CSV／TSV／JSONL／FASTA 写出不积累完整数据集，Parquet 每批最多缓存 `batch_size` 行并从第一批确定类型；需严格指定类型时使用 `TableSchema` 和 `export_table`。Parquet 需要 PyArrow，数组接口需要 NumPy。`load_arrays` 会将数组加载到内存并关闭 NPZ 文件。
