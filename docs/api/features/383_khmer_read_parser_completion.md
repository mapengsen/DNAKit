# khmer 有界读段解析

## khmer_parse_reads

**作用：** Read FASTA/FASTQ records, preserving names, sequences and FASTQ qualities.

**处理方法：** 调用真实 khmer.ReadParser；配对使用 iter_read_pairs。达到数量或碱基上限后停止收集并返回 truncated。

**API：** `dnakit.microbiome.khmer_p27_read_parser_functions.khmer_parse_reads(path: 'str | Path', *, max_records: 'int' = 1000, max_bases: 'int' = 1000000, verbose: bool = True) -> 'dict[str, Any]'`

**输入：** `path`, `max_records=1000`, `max_bases=1000000`, `verbose=True`.

**输出：** `records`, `record_count`, `base_count`, `truncated`, `provenance`.

**示例文件：** `pairs.fa`

```text
>read/1
ACGT
>read/2
TGCA
```

**示例代码：**

```python
from dnakit.microbiome.khmer_p27_read_parser_functions import khmer_parse_reads
result = khmer_parse_reads("pairs.fa", verbose=False)
print(result["record_count"], result["base_count"], result["truncated"])
```

**示例输出：**

```text
2 8 False
```

## khmer_parse_read_pairs

**作用：** Read interleaved FASTA/FASTQ pairs with khmer's name-pairing validation.

**处理方法：** 调用真实 khmer.ReadParser；配对使用 iter_read_pairs。达到数量或碱基上限后停止收集并返回 truncated。

**API：** `dnakit.microbiome.khmer_p27_read_parser_functions.khmer_parse_read_pairs(path: 'str | Path', *, max_records: 'int' = 1000, max_bases: 'int' = 1000000, verbose: bool = True) -> 'dict[str, Any]'`

**输入：** `path`, `max_records=1000`, `max_bases=1000000`, `verbose=True`.

**输出：** `records`, `record_count`, `base_count`, `truncated`, `provenance`.

**示例文件：** `pairs.fa`

```text
>read/1
ACGT
>read/2
TGCA
```

**示例代码：**

```python
from dnakit.microbiome.khmer_p27_read_parser_functions import khmer_parse_read_pairs
result = khmer_parse_read_pairs("pairs.fa", verbose=False)
print(result["record_count"], result["base_count"], result["truncated"])
```

**示例输出：**

```text
1 8 False
```

