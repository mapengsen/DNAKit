# 文件读写

读取和写入常见 DNA 序列、注释、表格及压缩文件。普通用户读取序列文件时只需记住 `dnakit.read()`，写出时只需记住 `dnakit.write()`；文件较大时仍使用相同方法，仅将读取参数设为 `mode="stream"`。

大文件不是一种新的文件格式，也不需要新的普通读写入口。原 `IO-005` 编号继续用于需求追踪，其流式、分块和索引能力统一归入 `IO-001` 的大文件用法。

读取接口默认保留输入顺序。

## 1) IO-001 序列格式

- **作用：** 读取和写出 FASTA、FASTQ、GenBank 文件，将序列及可用的 ID、质量值和注释转换为 DNAKit 记录，便于继续分析或格式转换。
- **统一 API：** `dnakit.read(source[必须], format[可选], config[可选], mode[可选])`、`dnakit.write(records[必须], target[必须], format[可选], config[可选])`；普通读取设 `mode="dna"`，大文件设 `mode="stream"`。
- **兼容 API：** `read_one()`、`read_set()` 仍可用，但新代码不需要记住这两个名称。
- **输入：** `mode="dna"` 无论文件含一条还是多条记录都返回 `DNA`；`mode="stream"` 返回单次消费的数据流。写出可直接传 `DNA`，也可直接传流式记录或其他记录迭代器，不需要“大文件写出”专用方法。
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

### 大文件参数与高级索引（需求追踪：IO-005） {#5-io-005}

- **普通用户：** 使用同一个 `dnakit.read(..., mode="stream")` 逐条读取，再把返回的数据流直接交给 `dnakit.write()`；只有需要固定分块时才使用 `iter_chunks()`。

```python
import dnakit

with dnakit.read("large.fa", mode="stream") as records:
    dnakit.write(records, "copy.fa")
```

- **高级用户：** 需要按记录 ID 或坐标随机访问大型 FASTA/FASTQ 时，才使用索引接口；这些接口不是普通文件读取的必学入口。
- **高级 API：** `dnakit.io.iter_chunks(values[必须], chunk_size[必须])`、`dnakit.io.build_fasta_index(source_path[必须], index_path[可选], overwrite[可选], max_records[可选], max_line_length[可选])`、`dnakit.io.load_fasta_index(index_path[必须], source_path[可选], verify_checksum[可选], max_index_bytes[可选])`、`dnakit.io.build_fastq_index(source_path[必须], index_path[可选], overwrite[可选], phred_offset[可选], max_records[可选], max_line_length[可选], max_record_bytes[可选], max_source_bytes[可选])`、`dnakit.io.load_fastq_index(index_path[必须], source_path[可选], verify_checksum[可选], max_index_bytes[可选], max_entries[可选])`、`dnakit.io.FastaIndex.fetch(record_id[必须], start[可选], end[可选], strand[可选], max_record_bytes[可选])`、`dnakit.io.FastqIndex.fetch(record_id[必须], start[可选], end[可选], strand[可选], max_record_bytes[可选])`。

## 2) IO-002 注释格式

- **作用：** 读取和写出 GFF3、BED、AGP 注释文件，并统一不同格式的坐标规则，生成可用于区域查询、序列切分和注释转换的结构化数据。
- **API：** `dnakit.io.read_gff3(source[必须], max_records[可选], max_line_length[可选], max_header_lines[可选])`、`dnakit.io.write_gff3(document[必须], target[必须], overwrite[可选], create_parents[可选])`、`dnakit.io.read_bed(source[必须], max_records[可选], max_line_length[可选], max_header_lines[可选])`、`dnakit.io.write_bed(document[必须], target[必须], overwrite[可选], create_parents[可选])`、`dnakit.io.read_agp(source[必须], max_records[可选], max_line_length[可选], max_header_lines[可选])`、`dnakit.io.write_agp(document[必须], target[必须], overwrite[可选], create_parents[可选])`。
- **输入：** 读取时必填注释文件路径或文本流，可选记录数、行长与 header 上限；写出时传入对应 document/entry 和目标路径或文本流。
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

## 3) IO-003 表格格式

- **作用：** 按指定表结构读写 CSV、TSV、JSON 和 Parquet，在 DNA 记录、metadata 与通用表格之间进行可复现转换。
- **API：** `dnakit.io.TableSchema(columns[必须], schema_version[可选], column_types[可选], nullable[可选])`、`dnakit.io.export_table(rows[必须], target[必须], format[必须], schema[必须], overwrite[可选], max_rows[可选], max_columns[可选], max_cell_characters[可选], max_output_bytes[可选], null_value[可选], parquet_compression[可选])`、`dnakit.io.read_table(source[必须], format[必须], schema[必须], missing_values[可选], max_rows[可选], max_columns[可选], max_cell_characters[可选], max_file_bytes[可选], max_decoded_bytes[可选])`、`dnakit.io.export_result(result[必须], target[必须], format[可选], overwrite[可选])`、`dnakit.io.parquet_backend_status()`。
- **输入：** `export_table()`/`read_table()` 使用二维行数据或目标路径、显式 `format` 和 `TableSchema`，可选列类型、nullable 字段、缺失值标记、压缩方式和资源上限；`export_result()` 接收带 `to_dict()` 的单个结果对象和目标路径，默认 JSON 并自动推导单行 schema。
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

## 4) IO-004 压缩文件

- **作用：** 直接读写 gzip 压缩的序列或表格文件，无需调用方先解压，并保留与普通文件相同的解析结果。
- **API：** `dnakit.read(source[必须], format[可选], config[可选], mode[可选])`、`dnakit.write(records[必须], target[必须], format[可选], config[可选])`；旧的 `read_one()` 仅为兼容入口。
- **输入：** 必填压缩文件路径或支持二进制的流；用 `compression="auto"|"none"|"gzip"` 控制压缩，可选写出压缩级别。
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
