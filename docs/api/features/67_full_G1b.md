# G1b：pyfaidx、dnaio 与 screed 离线封装

本组从 pyfastx、pyfaidx、dnaio、screed 和 py2bit 的公共函数清单中，选出与 DNA 序列、FASTA/FASTQ 读取、区间解析直接相关的函数，封装为 `dnakit.sequence` 下的离线函数。新增封装共 8 个：pyfaidx 3 个、dnaio 2 个、screed 3 个。pyfastx 与 py2bit 没有需要新增的封装。其余函数均已记录原因：与已有 dnakit 函数重复、已由 `sequence/backends.py` 封装、命令行入口、内部辅助、弃用别名，或上游测试代码。所有函数只读取显式给出的本地文件或字符串，不访问网络。

## 覆盖结果

下表由台账统计，状态含义见本组的台账记录。

| 包 | 新增封装 | 已有封装 | 与 dnakit 重复 | 同族变体 | 跳过 |
| --- | --- | --- | --- | --- | --- |
| pyfastx | 0 | 0 | 1 | 0 | 2 |
| pyfaidx | 3 | 0 | 2 | 0 | 15 |
| dnaio | 2 | 0 | 1 | 0 | 6 |
| screed | 3 | 3 | 16 | 5 | 68 |
| py2bit | 0 | 1 | 0 | 0 | 0 |

## 示例数据

下面的示例直接从模块导入函数。文件只用于演示接口，配对文件的读段名遵循 `/1` 与 `/2` 约定。数据库由 `build_screed_database` 在临时目录中生成。

**示例文件：** `reads.fa`

```text
>r1 first read
ACGTNACGT
>r2 second read
ggcc
```

**示例文件：** `reads.fq`

```text
@q1 desc
ACGTAC
+
IIII#!
@q2
GGCCAA
+
ABCDEF
```

**示例文件：** `mate_1.fq`

```text
@pair7/1 lane1
ACGTAC
+
IIIIII
@pair8/1
GGCC
+
IIII
```

**示例文件：** `mate_2.fq`

```text
@pair7/2 lane1
TTGCAA
+
IIIIII
@pair8/2
CCAA
+
IIII
```

## 1) 按显式对照表转换染色体名称

- **作用：** 按显式给出的旧名到新名对照表转换染色体名称。
- **处理方法：** 调用 pyfaidx.translate_chr_name 建立映射。两个列表长度必须相等，旧名不能重复。没有对照的名称进入 `unmapped`，不做猜测。
- **API：** `dnakit.sequence.translate_chromosome_names(names[必须], *, from_names[必须], to_names[必须], verbose[可选])`
- **输入：** 待转换的名称列表，以及两个等长的旧名、新名列表。
- **输出：** `translations` 按输入顺序给出 `name` 与 `translated`；`unmapped` 为排序后的未匹配名称；另有 `requested_count` 与 `translated_count`。
- **示例代码：**

```python
from dnakit.sequence import translate_chromosome_names

result = translate_chromosome_names(
    ["chr1", "chrM", "chr2", "chrX"],
    from_names=["chr1", "chr2", "chrM"],
    to_names=["1", "2", "MT"],
    verbose=False,
)
print(result["translations"])
print(result["unmapped"])
```

**示例输出：**

```text
[{'name': 'chr1', 'translated': '1'}, {'name': 'chrM', 'translated': 'MT'}, {'name': 'chr2', 'translated': '2'}]
['chrX']
```

## 2) 解析一行 BED 区间

- **作用：** 解析一行 BED 文本，得到染色体、0 起点和半开终点。
- **处理方法：** 调用 pyfaidx.bed_split。以 `#` 开头的行标记为注释。起点为负或终点早于起点时报错，空行和格式错误的行也报错。
- **API：** `dnakit.sequence.split_bed_entry(line[必须], verbose[可选])`
- **输入：** 一行至少包含三列（染色体、起点、终点）的 BED 文本。
- **输出：** `is_comment`、`chromosome`、`start`、`end` 与 `length`。注释行的坐标字段为 `None`。
- **示例代码：**

```python
from dnakit.sequence import split_bed_entry

entry = split_bed_entry("chr1\t10\t20\tgeneA", verbose=False)
print(entry["chromosome"], entry["start"], entry["end"], entry["length"])
print(split_bed_entry("# header", verbose=False)["is_comment"])
```

**示例输出：**

```text
chr1 10 20 10
True
```

## 3) 解析 UCSC 区域

- **作用：** 把 `chr:start-end` 形式的 UCSC 区域转换为 0 起点、半开区间；只给染色体名时表示整条序列。
- **处理方法：** 调用 pyfaidx.ucsc_split，并把 1 起点闭区间的起点减一。千分位逗号等无法解析的坐标报错，不会被当作整条序列。名称本身含冒号时必须带区间。
- **API：** `dnakit.sequence.split_ucsc_region(region[必须], verbose[可选])`
- **输入：** 一个区域字符串，例如 `chr1:101-200` 或 `chr7`。
- **输出：** `chromosome`、`start`、`end`、`length` 与 `is_whole_chromosome`。整条序列时坐标为 `None`。
- **示例代码：**

```python
from dnakit.sequence import split_ucsc_region

region = split_ucsc_region("chr1:101-200", verbose=False)
print(region["chromosome"], region["start"], region["end"], region["length"])
print(split_ucsc_region("chr7", verbose=False)["is_whole_chromosome"])
```

**示例输出：**

```text
chr1 100 200 100
True
```

## 4) 读取 FASTA/FASTQ 与配对读段

- **作用：** 用 dnaio 读取 FASTA 或 FASTQ 文件；也可以读取两个配对 FASTQ 文件。
- **处理方法：** 调用 `dnaio.open`。`max_records` 限制读取的条数，配对读取时按对计数。两个文件的读段名不匹配时，dnaio 会报错，封装转为 ConfigurationError。`interleaved=True` 用于单个交错排列的文件。
- **API：** `dnakit.sequence.read_fastx_records(source[必须], *, mate_source[可选], interleaved[可选], max_records[可选], verbose[可选])`
- **输入：** FASTA 或 FASTQ 文件路径；可选 `mate_source` 指定第二个配对文件，或用 `interleaved=True` 读取交错文件。
- **输出：** `layout`、`format`、`entry_unit`、`entry_count` 与 `records`。单端记录含 `name`、`id`、`comment`、`sequence`、`length`，FASTQ 另含 `qualities`；配对记录含 `read1` 与 `read2`。
- **示例代码：**

```python
from dnakit.sequence import read_fastx_records

result = read_fastx_records("reads.fq", verbose=False)
print(result["format"], result["entry_count"])
for record in result["records"]:
    print(record["id"], record["sequence"], record["qualities"])
paired = read_fastx_records("mate_1.fq", mate_source="mate_2.fq", verbose=False)
print(paired["layout"], paired["entry_count"])
for pair in paired["records"]:
    print(pair["read1"]["id"], pair["read2"]["id"], pair["read1"]["length"])
```

**示例输出：**

```text
fastq 2
q1 ACGTAC IIII#!
q2 GGCCAA ABCDEF
paired-end 2
pair7/1 pair7/2 6
pair8/1 pair8/2 4
```

## 5) 检查配对读段名

- **作用：** 判断一组读段名是否互为配对读段。
- **处理方法：** 只比较名称：为每个名称构造无序列、无质量的 dnaio.SequenceRecord，再调用 `dnaio.records_are_mates`。名称为空或含换行时报错。
- **API：** `dnakit.sequence.check_read_mates(names[必须], *, verbose[可选])`
- **输入：** 至少两个读段名的列表。
- **输出：** `mates` 布尔值，以及 `read_count` 与原样返回的 `names`。
- **示例代码：**

```python
from dnakit.sequence import check_read_mates

print(check_read_mates(["pair7/1", "pair7/2"], verbose=False)["mates"])
print(check_read_mates(["pair7/1", "pair8/2"], verbose=False)["mates"])
```

**示例输出：**

```text
True
False
```

## 6) 导出 screed 数据库

- **作用：** 把 screed 数据库（SQLite 文件）中的全部记录导出为 FASTA 或 FASTQ 文件。
- **处理方法：** 先统计记录数：空库或超过 `max_records` 时报错。随后调用 `screed.ToFasta` 或 `screed.ToFastq`。输出文件已存在时拒绝覆盖；写入失败或记录数不一致时删除部分输出。screed 会把序列按 80 列换行。
- **API：** `dnakit.sequence.export_screed_database(database[必须], output[必须], *, output_format[可选], max_records[可选], verbose[可选])`
- **输入：** `build_screed_database` 生成的数据库路径，以及一个尚不存在的输出路径。`output_format` 可为 `fasta`（默认）或 `fastq`。
- **输出：** `record_count`、`output_format`、`database` 与 `output`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import build_screed_database
from dnakit.sequence import export_screed_database

work = Path(tempfile.mkdtemp())
built = build_screed_database("reads.fa", database=work / "reads", verbose=False)
exported = export_screed_database(built["database"], work / "out.fa", verbose=False)
print(exported["record_count"], exported["output_format"])
print((work / "out.fa").read_text(encoding="utf-8"), end="")
```

**示例输出：**

```text
2 fasta
>r1 first read
ACGTNACGT
>r2 second read
ggcc
```

## 7) 解析 FASTA/FASTQ（screed 迭代器）

- **作用：** 用 screed 自带的 FASTA 或 FASTQ 迭代器解析未压缩文件。
- **处理方法：** `file_format` 选择 `screed.fasta.fasta_iter` 或 `screed.fastq.fastq_iter`。screed 把 FASTQ 标题存为 `annotations`，封装统一输出为 `description`。格式不符或超过 `max_records` 时报错。
- **API：** `dnakit.sequence.read_screed_fastx_records(source[必须], *, file_format[可选], parse_description[可选], max_records[可选], verbose[可选])`
- **输入：** 未压缩的 FASTA 或 FASTQ 文件，并通过 `file_format` 声明格式。
- **输出：** `format`、`record_count` 与 `records`。每条记录含 `name`、`description`、`length` 和 `sequence`，FASTQ 另含 `quality`。
- **示例代码：**

```python
from dnakit.sequence import read_screed_fastx_records

result = read_screed_fastx_records("reads.fa", verbose=False)
print(result["format"], result["record_count"])
for record in result["records"]:
    print(record["name"], record["description"], record["length"])
fastq = read_screed_fastx_records("reads.fq", file_format="fastq", verbose=False)
print(fastq["records"][0]["quality"])
```

**示例输出：**

```text
fasta 2
r1 first read 9
r2 second read 4
IIII#!
```

## 8) 将记录格式化为 FASTA/FASTQ 文本

- **作用：** 把显式给出的记录格式化为 FASTA 或 FASTQ 文本，不写入文件。
- **处理方法：** 每条记录通过 `screed.screedRecord.write_fastx` 写入内存缓冲区。名称不能含空白；序列只允许 ASCII 字母、`-` 和 `*`；FASTQ 要求质量串与序列等长且为可打印 ASCII。序列不做换行。
- **API：** `dnakit.sequence.format_fastx_records(records[必须], *, output_format[可选], max_records[可选], verbose[可选])`
- **输入：** 字典列表。每个字典含 `name` 与 `sequence`，可选 `description`；选择 `fastq` 时还需 `quality`。
- **输出：** `text`、`record_count` 与 `output_format`。
- **示例代码：**

```python
from dnakit.sequence import format_fastx_records

fasta = format_fastx_records(
    [
        {"name": "r1", "description": "first read", "sequence": "ACGTN"},
        {"name": "r2", "sequence": "ggcc"},
    ],
    verbose=False,
)
print(fasta["text"], end="")
fastq = format_fastx_records(
    [{"name": "q1", "description": "desc", "sequence": "ACGTAC", "quality": "IIII#!"}],
    output_format="fastq",
    verbose=False,
)
print(fastq["text"], end="")
```

**示例输出：**

```text
>r1 first read
ACGTN
>r2
ggcc
@q1 desc
ACGTAC
+
IIII#!
```

## 注意事项

- 封装函数从 `dnakit.sequence` 导出。
- `screed.screedRecord.write_fastx_pair` 的上游实现调用了未定义的 `write_record`，会抛出 `NameError`，因此没有封装。
- screed 对空数据库执行 `ToFasta` 或 `ToFastq` 会抛出 `TypeError`。`export_screed_database` 在导出前统计记录数，并给出明确错误。
- `pyfaidx.ucsc_split` 对 `chr1:0-5` 返回负起点；`split_ucsc_region` 将其视为错误。
- `pyfastx.reverse_complement` 不校验输入字符，而 `reverse_complement_dna` 会拒绝非 IUPAC 字符。
