# OPS-001 序列补全

公共处理接口支持[统一的文件与文本输入](00_processing_inputs.md)，具体入口及输出规则见该说明。

对 DNA 序列进行逆序、互补或反向互补，从而得到方向转换后的新序列。

## 1) OPS-001.1 逆序

- **作用：** 只把 DNA 序列的字符顺序倒过来而不替换碱基，用于检查或构造相反的位置顺序。
- **API：** `dnakit.ops.reverse(sequence[必须])`。
- **输入：** 必填一条 `DNASequence`；无额外必填参数。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.ops import reverse

seq = DNASequence("AATCGC")
print(reverse(seq).symbols)
```

- **示例结果：**

```text
CGCTAA
```

## 2) OPS-001.2 互补

- **作用：** 按碱基配对规则把每个 DNA 碱基替换为互补碱基，但保持原有位置顺序，用于构造对应互补链。
- **API：** `dnakit.ops.complement(sequence[必须])`。
- **输入：** 必填一条 `DNASequence`；无额外必填参数。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.ops import complement

seq = DNASequence("AATCGC")
print(complement(seq).symbols)
```

- **示例结果：**

```text
TTAGCG
```

## 3) OPS-001.3 反向互补

大型未压缩 FASTA/FNA（至少 16 MiB）在默认读取配置下会使用文件区间处理：
各工作进程直接读取自己的区间、校验并转换碱基，按原记录顺序写入临时输出的独立区间。
进程间只传递路径和偏移，不传递整条染色体；校验成功后才安装输出文件。
`workers=1` 使用相同的分块算法串行执行，`workers=-1` 在可用 CPU 和内存范围内并行。
使用多进程的脚本应保留 `if __name__ == "__main__":` 入口保护。
不规则换行、压缩文件、其他格式、显式 `read_config` 或输出字节限制继续使用通用处理路径。
读取配置和错误检查不会因启用多核而被放宽。

- **作用：** 同时反转序列顺序并替换为互补碱基，得到另一条链的 5′→3′ 表示，可用于反向链搜索、方向统一和引物分析。
- **API：** `dnakit.ops.reverse_complement(value, output_path=None, *, feature_policy="preserve", read_config=None, write_config=None, sequence_column=None, show_progress=True)`。
- **输入：** DNA 文本、文件路径，或既有的单记录 `DNA`、`DNASequence`、`DNARecord`。
  文本返回反向互补字符串；既有对象保留原返回类型，并同步 feature、strand 和逐碱基注释。
  文件按后缀识别 `.fna/.fa/.fasta`、`.txt`、`.csv/.tsv` 等 DNAKit I/O 格式，也支持 `.gz`。
  字符串中包含点号或路径分隔符时作为路径；无后缀文件请显式传 `Path`。
- **示例代码：**

```python
from dnakit import DNA
from dnakit.ops import reverse_complement

dna = DNA("AATCGC", id="seq-1")
print(reverse_complement(dna).symbols)
```

- **示例结果：**

```text
GCGATT
```

也可以直接处理文本或文件，不必自行循环读取、转换和写出：

```python
from dnakit.ops import reverse_complement

print(reverse_complement("AATCGC"))  # GCGATT；接受小写和空白，输出大写 IUPAC DNA
result = reverse_complement("hg38.fna")  # 自动生成 hg38_reverse_complement.fna
result = reverse_complement("sequences.txt", "rc.txt")
result = reverse_complement("sequences.csv", "rc.csv", sequence_column="sequence")
```

文件处理逐条读取记录，默认显示开始提示、已完成记录数和完成路径，可用 `show_progress=False` 关闭。
返回 `WriteResult`，其中 `record_count` 是记录数，`target_artifact.relative_path` 是相对于当前工作目录的输出路径。
未指定读取配置时，未压缩文件按实际文件大小设置读取上限，完整 hg38 无需手写三个上限参数；
gzip 默认限制解压后读取量为 20 GB，可传 `read_config=ReadConfig(...)` 显式调整。
原始输入不能作为输出，已有输出默认不覆盖；重新生成时传 `write_config=WriteConfig(overwrite=True)`。

TXT 默认每个非空行是一条 DNA 序列，输出每条序列一行，不包含记录名称。
若整个 TXT 是一条换行排版的序列，传 `read_config=ReadConfig(text_mode="sequence")`。
CSV 默认读取 `sequence` 列，可通过 `sequence_column` 指定其他列名。CSV 输出保留原列名、顺序及标签等附加列；记录注释按对应的反向互补规则映射。
