# cutadapt 与 umi-tools 的离线序列和 UMI 函数

这组函数位于 `dnakit.sequencing_qc`，包装 [cutadapt](https://github.com/marcelm/cutadapt)（5.x）与 [umi-tools](https://github.com/CGATOxford/UMI-tools)（1.1.x）的 Python 接口。cutadapt 部分计算邻近序列集合、接头 k-mer 搜索集合、poly-A 与 NextSeq 修剪位置，以及接头展开和接头文件读取；umi-tools 部分计算 UMI 距离与邻接关系、条码提取与白名单、CIGAR 剪接位置、BAM 标签和基因级 UMI 计数。所有函数只处理显式给出的序列、文件或计数表，不修剪 reads，也不运行 cutadapt 或 umi_tools 命令行，更不做 UMI 去重。

依赖都是可选的，只有调用对应函数时才会导入。缺少包时会报错并给出安装命令。

## 安装可选依赖

```bash
python -m pip install "dnakit[bio-cutadapt]"    # cutadapt 与 dnaio
python -m pip install "dnakit[bio-umi-tools]"   # umi-tools 及其运行时依赖
python -m pip install "dnakit[bio-reads]"       # pysam：BAM 读取，以及 umi-tools 的导入
```

umi-tools 导入时还会加载 pysam、matplotlib、scipy 和 pybktree；`detect_bam_tags()` 还需要 pandas。目前 `bio-umi-tools` 只列出 umi-tools 与 regex，因此需要另外安装这些包（见“当前范围”）。

各包的项目地址：[cutadapt](https://github.com/marcelm/cutadapt)、[umi-tools](https://github.com/CGATOxford/UMI-tools)、[pysam](https://github.com/pysam-developers/pysam)。

下面每个示例都由真实运行生成输出。示例直接从模块 `dnakit.sequencing_qc` 与 `dnakit.sequencing_qc` 导入；API 行给出的包级名称与之指向同一函数。

## 1) `dnakit.sequencing_qc.hamming_neighbours()`

- **作用：** 列出与给定 DNA 序列的 Hamming 距离恰好为（`exact=True`）或不超过（默认）`max_distance` 的全部序列，常用于 barcode 纠错候选集合。
- **处理方法：** `exact=True` 调用 `cutadapt.align.hamming_sphere`；默认调用 `cutadapt.align.hamming_environment`，返回距离 0 到 `max_distance` 的序列。每条结果给出距离和匹配碱基数。结果按距离、再按序列排序。数量超过 `max_neighbours` 时在枚举过程中报错。字母表只接受 A、C、G、T。
- **API：** `dnakit.sequencing_qc.hamming_neighbours(sequence, max_distance, *, exact=False, max_neighbours=10000)`。
- **输入：** DNA 字符串；`max_distance` 为 0 到 3 的整数且不超过序列长度；`exact` 为布尔值；`max_neighbours` 为显式上限。
- **输出：** 含 `count` 和 `neighbours` 的字典。`neighbours` 中每项有 `sequence`、`distance` 和 `matches`；`provenance` 记录 cutadapt 版本与所用函数。
- **示例代码：**

```python
from dnakit.sequencing_qc import hamming_neighbours

exact = hamming_neighbours("ACGT", 1, exact=True)
print(exact["count"], [row["sequence"] for row in exact["neighbours"]][:4])
within = hamming_neighbours("ACGT", 1)
print(within["count"], within["neighbours"][0])
```

- **示例输出：**

```text
12 ['AAGT', 'ACAT', 'ACCT', 'ACGA']
13 {'sequence': 'ACGT', 'distance': 0, 'matches': 4}
```

## 2) `dnakit.sequencing_qc.edit_neighbours()`

- **作用：** 列出与给定 DNA 序列的编辑距离（允许插入、删除和替换）不超过 `max_distance` 的全部序列。
- **处理方法：** 调用 `cutadapt.align.edit_environment`（Cython 实现）。它为每个序列给出编辑距离和最优比对中的匹配数。同一序列只保留最小距离。结果按距离、再按序列排序。
- **API：** `dnakit.sequencing_qc.edit_neighbours(sequence, max_distance, *, max_neighbours=10000)`。
- **输入：** DNA 字符串；`max_distance` 为 0 到 3 的整数；`max_neighbours` 为显式上限。
- **输出：** `count`、`neighbours`（每项含 `sequence`、`distance`、`matches`）和 `provenance`。
- **示例代码：**

```python
from dnakit.sequencing_qc import edit_neighbours

result = edit_neighbours("AC", 1)
print(result["count"], [row["sequence"] for row in result["neighbours"]][:5])
```

- **示例输出：**

```text
19 ['AC', 'A', 'AA', 'AAC', 'ACA']
```

## 3) `dnakit.sequencing_qc.adapter_kmer_positions()`

- **作用：** 为接头序列生成 cutadapt 的 k-mer 搜索集合。对每个搜索偏移，给出必须出现在 read 中的 k-mer，用于快速定位接头。
- **处理方法：** 调用 `cutadapt.kmer_heuristic.create_positions_and_kmers`。`start` 与 `stop` 是 cutadapt 的偏移量，负数表示从 read 末端往回数；`stop` 为 `None` 表示延伸到 read 末端。条目内的 k-mer 排序，条目按 `(start, stop)` 排序，保证结果确定。
- **API：** `dnakit.sequencing_qc.adapter_kmer_positions(adapter, *, min_overlap=3, max_error_rate=0.1, back_adapter=True, front_adapter=False, internal=True)`。
- **输入：** 接头 DNA 字符串；`min_overlap` 为不超过接头长度的正整数；`max_error_rate` 在 [0, 1) 内；三个布尔开关分别对应 3' 接头、5' 接头和内部匹配。
- **输出：** `entries`（每项含 `start`、`stop`、`kmers`）、`entry_count`、`kmer_count` 和 `provenance`。
- **示例代码：**

```python
from dnakit.sequencing_qc import adapter_kmer_positions

result = adapter_kmer_positions("ACGTACGT", min_overlap=3)
for entry in result["entries"]:
    print(entry["start"], entry["stop"], entry["kmers"])
```

- **示例输出：**

```text
-8 None ['ACGTA']
-4 None ['ACGT']
-3 None ['ACG']
0 None ['ACGTACGT']
```

## 4) `dnakit.sequencing_qc.adapter_kmer_table()`

- **作用：** 估计接头 k-mer 搜索集合在随机序列中被误命中的概率，用于判断某个接头的 k-mer 是否过于常见。
- **处理方法：** 先调用 `create_positions_and_kmers`，再调用 `cutadapt.kmer_heuristic.kmer_probability_analysis`，并以 `default_length` 作为读长。上游返回制表符分隔的文本表，本函数把每一行解析为一条记录，并提取整体命中概率。原始文本保存在 `table_text` 中。
- **API：** `dnakit.sequencing_qc.adapter_kmer_table(adapter, *, min_overlap=3, max_error_rate=0.1, back_adapter=True, front_adapter=False, internal=True, default_length=150)`。
- **输入：** 与 `adapter_kmer_positions()` 相同，另加 `default_length`（1 到 1000000 之间的整数）。
- **输出：** `rows`（每项含 `kmer`、`start`、`stop`、`considered_sites`、`hit_chance_percent`）、`profile_hit_chance_percent`、`table_text` 和 `provenance`。
- **示例代码：**

```python
from dnakit.sequencing_qc import adapter_kmer_table

result = adapter_kmer_table("ACGTACGT")
for row in result["rows"]:
    print(row["kmer"], row["start"], row["considered_sites"], row["hit_chance_percent"])
print(result["profile_hit_chance_percent"])
```

- **示例输出：**

```text
ACG -3 1 1.56
ACGT -4 1 0.39
ACGTA -8 4 0.39
ACGTACGT 0 143 0.22
2.54
```

## 5) `dnakit.sequencing_qc.sequence_kmer_chunks()`

- **作用：** 把序列切成长度几乎相等的若干段，返回去重后的段序列。cutadapt 用它生成 k-mer 种子。
- **处理方法：** 调用 `cutadapt.kmer_heuristic.kmer_chunks`。序列长度除以段数的余数决定有多少段比其余段长一个碱基。结果排序后返回。输入允许 IUPAC 字母。
- **API：** `dnakit.sequencing_qc.sequence_kmer_chunks(sequence, chunks)`。
- **输入：** 序列字符串；`chunks` 为 1 到序列长度之间的整数。
- **输出：** `kmers`（排序列表）、`count`、`sequence`、`chunks` 和 `provenance`。
- **示例代码：**

```python
from dnakit.sequencing_qc import sequence_kmer_chunks

result = sequence_kmer_chunks("ACGTACGTAC", 3)
print(result["kmers"])
```

- **示例输出：**

```text
['ACG', 'ACGT', 'TAC']
```

## 6) `dnakit.sequencing_qc.poly_a_tail_start()`

- **作用：** 找到序列 3' 端的 poly-A 尾；若 `reverse_complement=True`，则找 5' 端的 poly-T 头。
- **处理方法：** 调用 `cutadapt.qualtrim.poly_a_trim_index`。`index` 是尾部起点（或头部终点），`length` 是尾部或头部的碱基数。cutadapt 忽略短于 3 个碱基的尾；没有尾时 `index` 等于序列长度。
- **API：** `dnakit.sequencing_qc.poly_a_tail_start(sequence, *, reverse_complement=False)`。
- **输入：** 序列字符串（允许 IUPAC 字母）；布尔开关 `reverse_complement`。
- **输出：** `index`、`length`、`sequence_length`、`reverse_complement` 和 `provenance`。
- **示例代码：**

```python
from dnakit.sequencing_qc import poly_a_tail_start

tail = poly_a_tail_start("ACGTAAAAAAAAA")
print(tail["index"], tail["length"])
head = poly_a_tail_start("TTTTTACGT", reverse_complement=True)
print(head["index"], head["length"])
```

- **示例输出：**

```text
4 9
5 5
```

## 7) `dnakit.sequencing_qc.nextseq_trim_bounds()`

- **作用：** 计算 NextSeq 双色化学的 3' 端修剪位置。3' 端连续的 G 碱基（无颜色信号，常出现在片段末端）会被去掉。
- **处理方法：** 把序列和质量串封装为 `dnaio` 记录后调用 `cutadapt.qualtrim.nextseq_trim_index`。算法与 `find_quality_trim_bounds()` 相同（BWA 的累加和规则），区别是 G 碱基的质量按 `cutoff - 1` 计算。只修剪 3' 端，因此 `index` 即保留长度。
- **API：** `dnakit.sequencing_qc.nextseq_trim_bounds(sequence, qualities, *, cutoff=20, phred_offset=33)`。
- **输入：** 序列；与之等长的 Phred 质量串；`cutoff` 为 0 到 93；`phred_offset` 为 33 或 64。
- **输出：** `index`、`retained_length`、`length`、`cutoff`、`phred_offset` 和 `provenance`。
- **示例代码：**

```python
from dnakit.sequencing_qc import nextseq_trim_bounds

result = nextseq_trim_bounds("ACGTGGGGGG", "I" * 10, cutoff=20)
print(result["index"], result["retained_length"])
```

- **示例输出：**

```text
4 4
```

## 8) `dnakit.sequencing_qc.expand_adapter_braces()`

- **作用：** 展开 cutadapt 的 `x{n}` 重复写法，便于书写含重复碱基的接头。
- **处理方法：** 调用 `cutadapt.parser.expand_braces`，把每个 `x{n}` 替换为 n 个 x。展开前先按重复次数估算长度，估算值超过 `max_length` 时直接拒绝。格式错误（缺少数字、括号不配对）会转为 `ConfigurationError`。
- **API：** `dnakit.sequencing_qc.expand_adapter_braces(sequence, *, max_length=1000000)`。
- **输入：** 含 `x{n}` 写法的字符串；`max_length` 为显式上限。
- **输出：** `expanded`、`input_length`、`expanded_length` 和 `provenance`。
- **示例代码：**

```python
from dnakit.sequencing_qc import expand_adapter_braces

print(expand_adapter_braces("TGA{5}CT")["expanded"])
print(expand_adapter_braces("AC{3}G{2}")["expanded"])
```

- **示例输出：**

```text
TGAAAAACT
ACCCGG
```

**示例文件：** `adapters.fa`

```text
>ad1 description
ACGTACGT
>ad2
TTGGAACC
```

## 9) `dnakit.sequencing_qc.read_adapter_fasta()`

- **作用：** 从本地 FASTA 文件读取接头名称和序列。
- **处理方法：** 调用 `cutadapt.parser.read_adapters_fasta`。名称取表头的第一个词，空表头得到 `None`；序列保持文件中的写法；结果按文件顺序返回。不是 FASTA 的文件，或记录数超过 `max_records` 的文件，都会报错。
- **API：** `dnakit.sequencing_qc.read_adapter_fasta(source, *, max_records=10000)`。
- **输入：** FASTA 文件路径；`max_records` 为显式上限。
- **输出：** `records`（每项含 `name` 和 `sequence`）、`count`、`source` 和 `provenance`。
- **示例代码：**

```python
from dnakit.sequencing_qc import read_adapter_fasta

result = read_adapter_fasta("adapters.fa")
for record in result["records"]:
    print(record["name"], record["sequence"])
```

- **示例输出：**

```text
ad1 ACGTACGT
ad2 TTGGAACC
```

## 10) `dnakit.sequencing_qc.parse_adapter_spec()`

- **作用：** 解析一条 cutadapt 接头规格，返回接头对象及其设置。
- **处理方法：** 调用 `cutadapt.parser.make_adapters_from_one_specification`。规格可以是 `ACGT`、`name=ACGT`、连接接头 `ACGT...TTTT`（不能用于 `anywhere`），或 `file:<路径>`（读取本地 FASTA 中的全部接头）。`file:` 后可用分号附加参数，如 `;max_errors=0.2`。`adapter_type` 取 `front`（-g）、`back`（-a）或 `anywhere`（-b）。连接接头返回前后两段序列。
- **API：** `dnakit.sequencing_qc.parse_adapter_spec(spec, *, adapter_type="back", max_errors=0.1, min_overlap=3, indels=True, max_adapters=1000)`。
- **输入：** 规格字符串；`adapter_type`；`max_errors` 在 [0, 1) 内；`min_overlap`；`indels` 为布尔值；`max_adapters` 为显式上限。
- **输出：** `adapters`（每项含 `name`、`adapter_class`、`sequence`、`max_error_rate`、`min_overlap`、`indels` 等；连接接头为 `front_sequence` 和 `back_sequence`）、`count` 和 `provenance`。
- **示例代码：**

```python
from dnakit.sequencing_qc import parse_adapter_spec

result = parse_adapter_spec("file:adapters.fa", adapter_type="back", max_errors=0.2)
for row in result["adapters"]:
    print(row["name"], row["adapter_class"], row["sequence"], row["max_error_rate"])
linked = parse_adapter_spec("ACGT...TTTT", adapter_type="back")
print(linked["adapters"][0]["front_sequence"], linked["adapters"][0]["back_sequence"])
```

- **示例输出：**

```text
ad1 BackAdapter ACGTACGT 0.2
ad2 BackAdapter TTGGAACC 0.2
ACGT TTTT
```

## 11) `dnakit.sequencing_qc.average_umi_distance()`

- **作用：** 计算一组等长 UMI 两两之间的平均 Hamming 距离，用于评估 UMI 多样性。
- **处理方法：** 把 UMI 转为字节后调用 `umi_tools.umi_methods.get_average_umi_distance`（上游要求字节输入）。每一对都计入，重复 UMI 的距离为 0。只有一个 UMI 时上游返回哨兵值 -1，本函数改为 `None`。
- **API：** `dnakit.sequencing_qc.average_umi_distance(umis, *, max_umis=2000)`。
- **输入：** UMI 字符串列表（A、C、G、T、N），长度必须相同；`max_umis` 为显式上限。
- **输出：** `umi_count`、`umi_length`、`pair_count`、`average_distance` 和 `provenance`。
- **示例代码：**

```python
from dnakit.sequencing_qc import average_umi_distance

result = average_umi_distance(["AAAA", "AAAT", "CCCC"])
print(result["average_distance"], result["pair_count"])
```

- **示例输出：**

```text
3.0 3
```

## 12) `dnakit.sequencing_qc.find_umi_neighbours()`

- **作用：** 找出 Hamming 距离不超过 `max_distance` 的 UMI 对，用于 UMI 聚类前的候选筛选。
- **处理方法：** 先去重并排序。`umi_tools.network.build_substr_idx` 把每个 UMI 切成 `max_distance + 1` 段，距离内的两个 UMI 至少有一段相同；`iter_nearest_neighbours` 生成候选对；再用上游 `umi_tools.umi_methods.edit_distance` 核对每一对。结果中 `umi_a` 小于 `umi_b`，按字典序排列。
- **API：** `dnakit.sequencing_qc.find_umi_neighbours(umis, *, max_distance=1, max_umis=2000)`。
- **输入：** 等长的 UMI 列表；`max_distance` 为 0 到 3，且小于 UMI 长度。
- **输出：** `unique_umis`、`pairs`（每项含 `umi_a`、`umi_b`、`distance`）、`pair_count` 和 `provenance`。
- **示例代码：**

```python
from dnakit.sequencing_qc import find_umi_neighbours

result = find_umi_neighbours(["AAAA", "AAAT", "CCCC", "AAAG"], max_distance=1)
for pair in result["pairs"]:
    print(pair["umi_a"], pair["umi_b"], pair["distance"])
```

- **示例输出：**

```text
AAAA AAAG 1
AAAA AAAT 1
AAAG AAAT 1
```

## 13) `dnakit.sequencing_qc.umi_connected_component()`

- **作用：** 返回无向 UMI 图中包含某个节点的连通分量，用于把通过错配相连的 UMI 归为一组。
- **处理方法：** 调用 `umi_tools.network.breadth_first_search`。邻接表中缺少的邻居会被补为没有边的节点。结果排序后返回。
- **API：** `dnakit.sequencing_qc.umi_connected_component(adjacency, node, *, max_nodes=100000)`。
- **输入：** `adjacency`，即 UMI 到邻居列表的映射；`node` 为起点 UMI；`max_nodes` 为显式上限。
- **输出：** `component`（排序列表）、`size`、`node` 和 `provenance`。
- **示例代码：**

```python
from dnakit.sequencing_qc import umi_connected_component

adjacency = {"A": ["B"], "B": ["A", "C"], "D": ["E"]}
result = umi_connected_component(adjacency, "A")
print(result["component"], result["size"])
```

- **示例输出：**

```text
['A', 'B', 'C'] 3
```

## 14) `dnakit.sequencing_qc.mask_low_quality_umi()`

- **作用：** 把质量低于阈值的 UMI 碱基替换为 N。
- **处理方法：** 调用 `umi_tools.extract_methods.mask_umi`，并用 `get_below_threshold` 找出被替换的位置。`quality_encoding` 可取 `phred33`、`phred64` 或 `solexa`。只有质量严格小于阈值的碱基才会被屏蔽。
- **API：** `dnakit.sequencing_qc.mask_low_quality_umi(umi, qualities, *, quality_encoding="phred33", threshold=10)`。
- **输入：** UMI；与之等长的质量串；质量编码；阈值（0 到 93）。
- **输出：** `masked_umi`、`masked_positions`（从 0 开始）、`masked_count` 和 `provenance`。
- **示例代码：**

```python
from dnakit.sequencing_qc import mask_low_quality_umi

result = mask_low_quality_umi("ACGT", "II!I", threshold=10)
print(result["masked_umi"], result["masked_positions"])
```

- **示例输出：**

```text
ACNT [2]
```

## 15) `dnakit.sequencing_qc.extract_barcode_positions()`

- **作用：** 按 0 起始的位置，从 read 中去掉 UMI、细胞条码和需丢弃的碱基，并给出条码序列与质量。
- **处理方法：** 调用 `umi_tools.extract_methods.extractSeqAndQuals`。上游只返回保留序列的质量，以及 UMI 和细胞条码的质量，不返回条码碱基；本函数按位置从输入中取出条码序列（按 read 中的顺序）。三组位置不能重叠。
- **API：** `dnakit.sequencing_qc.extract_barcode_positions(sequence, qualities, *, umi_positions=(), cell_positions=(), discard_positions=(), retain_umi=False)`。
- **输入：** read 序列；Phred+33 质量串；三组 0 起始的位置列表；`retain_umi` 为布尔值。
- **输出：** `read_sequence`、`read_qualities`、`umi_sequence`、`umi_qualities`、`cell_sequence`、`cell_qualities`、`positions` 和 `provenance`。注意：`retain_umi=True` 时上游会把所有保留碱基的质量复制到 `umi_qualities`，本函数原样返回。
- **示例代码：**

```python
from dnakit.sequencing_qc import extract_barcode_positions

result = extract_barcode_positions(
    "ACGTACGT", "IIIIHHHH",
    umi_positions=[0, 1], cell_positions=[2], discard_positions=[7],
)
print(result["read_sequence"], result["read_qualities"])
print(result["umi_sequence"], result["umi_qualities"], result["cell_sequence"])
```

- **示例输出：**

```text
TACG IHHH
AC II G
```

## 16) `dnakit.sequencing_qc.tag_read_name_with_umi()`

- **作用：** 把 UMI（以及可选的细胞条码）加到 FASTQ 标题的第一个词之后。
- **处理方法：** 调用 `umi_tools.extract_methods.addBarcodesToIdentifier`。上游按空格切分标题，在第一个词后追加 `分隔符 + [细胞条码 + 分隔符] + UMI`，标题其余部分保持不变。
- **API：** `dnakit.sequencing_qc.tag_read_name_with_umi(read_name, umi, *, cell_barcode=None, separator="_")`。
- **输入：** 标题字符串；UMI；可选的细胞条码；分隔符（1 到 4 个非空白字符）。
- **输出：** `tagged_name`、`umi`、`cell_barcode`、`separator` 和 `provenance`。
- **示例代码：**

```python
from dnakit.sequencing_qc import tag_read_name_with_umi

result = tag_read_name_with_umi("read1 extra words", "ACGT", cell_barcode="TTGG")
print(result["tagged_name"])
```

- **示例输出：**

```text
read1_TTGG_ACGT extra words
```

## 17) `dnakit.sequencing_qc.cigar_splice_offset()`

- **作用：** 返回 CIGAR 中第一个剪接位置（N）或末端软剪切（S）相对 read 起点的偏移，用于 RNA-seq 的剪接判断。
- **处理方法：** 把 CIGAR 解析为 `(操作码, 长度)`，调用 `umi_tools.sam_methods.find_splice`。偏移只计入消耗参考序列的操作（M、D、=、X），跳过插入、硬剪切和填充。开头的软剪切计入偏移，但不算剪接。`reverse=True` 从 3' 端扫描，与 umi-tools 处理反向链 read 的方式一致。没有剪接时 `splice_offset` 为 `None`。
- **API：** `dnakit.sequencing_qc.cigar_splice_offset(cigar, *, reverse=False)`。
- **输入：** CIGAR 字符串，如 `10M100N20M`；布尔值 `reverse`。
- **输出：** `splice_offset`、`has_splice`、`reverse`、`cigar` 和 `provenance`。
- **示例代码：**

```python
from dnakit.sequencing_qc import cigar_splice_offset

print(cigar_splice_offset("10M100N20M")["splice_offset"])
print(cigar_splice_offset("30M")["splice_offset"])
print(cigar_splice_offset("10M100N20M", reverse=True)["splice_offset"])
```

- **示例输出：**

```text
10
None
20
```

## 18) `dnakit.sequencing_qc.parse_umi_read_name()`

- **作用：** 从 read 名称末尾取出 UMI，必要时同时取出细胞条码。
- **处理方法：** 调用 `umi_tools.sam_methods.get_umi_read_string`（只取 UMI）或 `get_cell_umi_read_string`（先取细胞条码，再取 UMI）。两者都按分隔符切分名称，取最后一个字段（以及它前面的字段）。字段不足时报错。
- **API：** `dnakit.sequencing_qc.parse_umi_read_name(read_name, *, separator="_", cell_barcode=False)`。
- **输入：** read 名称；分隔符；是否同时取细胞条码。
- **输出：** `umi`、`cell_barcode`（未取时为 `None`）、`mode` 和 `provenance`。
- **示例代码：**

```python
from dnakit.sequencing_qc import parse_umi_read_name

result = parse_umi_read_name("read1_CELL1_UMI1", cell_barcode=True)
print(result["cell_barcode"], result["umi"])
print(parse_umi_read_name("read1_UMI1")["umi"])
```

- **示例输出：**

```text
CELL1 UMI1
UMI1
```

## 19) `dnakit.sequencing_qc.barcode_indel_matches()`

- **作用：** 在白名单中寻找与条码相差若干替换、插入或缺失的条码，用于条码纠错。
- **处理方法：** 调用 `umi_tools.whitelist_methods.checkError`，它使用模糊正则匹配。上游跳过与条码本身相同的白名单条目，并在找到两个匹配后立即停止，因此 `matches` 最多包含两个条码；`ambiguous` 为真表示至少找到两个。白名单中有与条码完全相同的条目时，`exact_in_whitelist` 为真，该条目不计入匹配。
- **API：** `dnakit.sequencing_qc.barcode_indel_matches(barcode, whitelist, *, max_errors=1)`。
- **输入：** 条码（A、C、G、T）；白名单列表；`max_errors` 为 0 到 3 的整数。
- **输出：** `matches`、`match_count`、`ambiguous`、`exact_in_whitelist`、`whitelist_size` 和 `provenance`。
- **示例代码：**

```python
from dnakit.sequencing_qc import barcode_indel_matches

result = barcode_indel_matches("ACGT", ["ACGA", "ACGTT", "TTTT"])
print(result["matches"], result["ambiguous"])
```

- **示例输出：**

```text
['ACGA', 'ACGTT'] True
```

**示例文件：** `whitelist.tsv`

```text
# barcode<TAB>corrections
AAAA	CCCC,CCCG
GGGG	TTTT
```

## 20) `dnakit.sequencing_qc.read_barcode_whitelist()`

- **作用：** 读取用户提供的细胞条码白名单 TSV，并可选地给出纠错映射。
- **处理方法：** 调用 `umi_tools.whitelist_methods.getUserDefinedBarcodes`。文件第一列是白名单条码，以 `#` 开头的行是注释，空行会报错。`error_correction` 取值为：`none`（只读白名单）；`file`（第二列列出可纠正为第一列条码的逗号分隔条码）；`derived`（自动推导：白名单条码的每个碱基可替换为 A、C、G、T、N 中的其他字符，`threshold` 为最多替换的碱基数，取 1 或 2）。`threshold=2` 时，同时与两个白名单条码都在范围内的条码不进入映射。
- **API：** `dnakit.sequencing_qc.read_barcode_whitelist(source, *, error_correction="none", threshold=1, max_barcodes=100000)`。
- **输入：** TSV 文件路径；`error_correction`；`threshold`（1 或 2）；`max_barcodes` 为显式上限。
- **输出：** `barcodes`、`barcode_count`、`error_map`（不纠错时为 `None`）、`error_correction`、`threshold` 和 `provenance`。`error_map` 把被纠正的条码映射到白名单条码。
- **示例代码：**

```python
from dnakit.sequencing_qc import read_barcode_whitelist

result = read_barcode_whitelist("whitelist.tsv", error_correction="file")
print(result["barcodes"])
print(result["error_map"])
```

- **示例输出：**

```text
['AAAA', 'GGGG']
{'CCCC': 'AAAA', 'CCCG': 'AAAA', 'TTTT': 'GGGG'}
```

## 21) `dnakit.sequencing_qc.call_cell_whitelist()`

- **作用：** 根据各条码的 read 计数，用拐点法判断真正的细胞条码。
- **处理方法：** 调用 `umi_tools.whitelist_methods.getCellWhitelist`。`knee_method="distance"` 取累计曲线上离首尾连线最远的点（`getKneeEstimateDistance`）；`"density"` 使用高斯密度（`getKneeEstimateDensity`）。`expect_cells` 仅用于 density；`cell_number` 固定接受的条码数。`error_correct_threshold` 大于 0 时，返回每个被接受条码在该 Hamming 距离内的非白名单条码（`error_map`）。上游要求输入为 `collections.Counter`，本函数会自动构造。
- **API：** `dnakit.sequencing_qc.call_cell_whitelist(cell_barcode_counts, *, knee_method="distance", expect_cells=None, cell_number=None, error_correct_threshold=0, max_barcodes=100000)`。
- **输入：** 条码到计数的映射（计数为正整数）；`knee_method`；可选的 `expect_cells` 与 `cell_number`；`error_correct_threshold`（0 到 2）。
- **输出：** `whitelist`、`cell_count`、`input_barcodes`、`knee_method`、`error_map`（或 `None`）和 `provenance`。
- **示例代码：**

```python
from dnakit.sequencing_qc import call_cell_whitelist

counts = {"AAAA": 100, "AAAT": 90, "CCCC": 10, "GGGG": 9, "TTTT": 1, "ACGT": 50}
result = call_cell_whitelist(counts, knee_method="distance")
print(result["whitelist"], result["cell_count"])
```

- **示例输出：**

```text
['AAAA', 'AAAT', 'ACGT'] 3
```

**示例文件：** `reads.sam`

```text
@HD	VN:1.6	SO:unsorted
@SQ	SN:chr1	LN:100
r1	0	chr1	1	60	4M	*	0	0	ACGT	IIII	NH:i:1
r2	0	chr1	9	60	4M	*	0	0	ACGT	IIII	NH:i:1	XT:A:U
```

## 22) `dnakit.sequencing_qc.detect_bam_tags()`

- **作用：** 检查 SAM 或 BAM 文件前若干条已比对 read 是否带有 NH、X0 和 XT 标签。
- **处理方法：** 调用 `umi_tools.dedup.detect_bam_features`，通过 pysam 读取文件。只检查已比对的 read；某个标签只有在所有被检查的 read 都带有它时才记为 `True`。上游检查 `n_entries + 1` 条记录。
- **API：** `dnakit.sequencing_qc.detect_bam_tags(source, *, n_entries=1000)`。
- **输入：** `.bam` 或 `.sam` 文件路径；`n_entries` 为 1 到 1000000 之间的整数。
- **输出：** `tags`（NH、X0、XT 的布尔值）、`n_entries`、`source` 和 `provenance`。
- **示例代码：**

```python
from dnakit.sequencing_qc import detect_bam_tags

result = detect_bam_tags("reads.sam")
print(result["tags"])
```

- **示例输出：**

```text
{'NH': True, 'X0': False, 'XT': False}
```

**示例文件：** `assignments.tsv`

```text
r1_CELL1_UMI1	GENE_A
r2_CELL1_UMI1	GENE_A
r3_CELL1_UMI2	GENE_A
r4_CELL2_UMI1	GENE_B
r5_CELL1_UMI1	GENE_B
```

## 23) `dnakit.sequencing_qc.count_umis_per_gene()`

- **作用：** 按基因统计 read 数，以及不同 UMI 的数量（或不同细胞条码与 UMI 组合的数量）。
- **处理方法：** 输入是两列、无表头的 TSV：read 名称和分配的基因。调用 `umi_tools.sam_methods.get_gene_count_tab`，从 read 名称中取出 UMI（`cell_barcode=True` 时同时取出细胞条码）。上游假定同一基因的行连续，本函数会按基因合并，不依赖输入顺序。上游遇到空输入会产生一个 `None` 基因，本函数忽略它。
- **API：** `dnakit.sequencing_qc.count_umis_per_gene(source, *, cell_barcode=False, separator="_", max_rows=1000000)`。
- **输入：** TSV 文件（每行两列，均不能为空）；read 名称中的分隔符；`max_rows` 为显式上限。
- **输出：** `genes`（按基因排序，每项含 `gene`、`reads`、`unique_umis`）、`gene_count`、`row_count`、`cell_barcode` 和 `provenance`。
- **示例代码：**

```python
from dnakit.sequencing_qc import count_umis_per_gene

result = count_umis_per_gene("assignments.tsv", cell_barcode=True)
for gene in result["genes"]:
    print(gene["gene"], gene["reads"], gene["unique_umis"])
```

- **示例输出：**

```text
GENE_A 3 2
GENE_B 2 2
```

## 验证方式

测试文件是 `others/tests/unit/test_full_G9.py`，共 30 项。主要检查：

- 每个包装函数至少有一项测试。能直接与上游比较的，比较上游的返回值，例如 `cutadapt.align.hamming_sphere`、`cutadapt.qualtrim.poly_a_trim_index`、`umi_tools.network.breadth_first_search` 和 `umi_tools.whitelist_methods.getCellWhitelist`。其余的用手工计算的结果核对：例如与 `ACGT` 的 Hamming 距离恰好为 1 的序列有 12 个；与 `AC` 的编辑距离不超过 1 的序列有 19 个（这一数字也由上游 `edit_distance` 的穷举结果核对）。
- UMI 配对搜索用随机生成的 UMI 与上游 Hamming 距离做穷举比较，检查 Hamming 距离内的所有配对都被找到。
- 文件类函数（接头 FASTA、白名单、BAM、基因计数表）使用 `tmp_path` 中的小文件，不依赖仓库数据。
- 错误输入（非法碱基、越界位置、重叠位置、超过上限、格式错误的接头括号等）会抛出 `ConfigurationError`。
- 缺少 `cutadapt` 或 `umi_tools` 时，通过屏蔽 `sys.modules` 验证会抛出 `BackendUnavailableError` 并提示 `bio-cutadapt` 或 `bio-umi-tools`。
- 分组结果（k-mer 搜索集合、UMI 图、白名单映射、基因计数）在 `PYTHONHASHSEED` 为 1 与 2 时输出完全一致。
- 没有函数访问网络，因此没有需要模拟的网络函数。

## 当前范围

- 本批包装 23 个函数：cutadapt 10 个，umi-tools 13 个。两个上游包的其余条目按下表分类（数字来自本批的账本记录）。所有函数都接受仅关键字参数 `verbose=True`，用于控制进度输出，未列在 API 行中。

| 状态 | 含义 | cutadapt | umi-tools |
| --- | --- | ---: | ---: |
| `wrapped`（包装） | 本批包装的函数 | 10 | 13 |
| `already-wrapped`（已有封装） | `dnakit.sequencing_qc` 中已有对应封装 | 4 | 2 |
| `duplicate-of`（与现有函数重复） | dnakit 其他函数已经提供相同计算 | 3 | 5 |
| `family-variant`（族内变体） | 由包装函数的参数或别名覆盖 | 13 | 7 |
| `cli`（命令行入口） | 命令行与脚本入口，不重新包装 | 22 | 13 |
| `internal`（内部辅助） | 日志、临时文件、资源限制或内部步骤 | 27 | 16 |
| `other`（其他） | 需要 cutadapt 报告对象、pysam 记录或 pandas 数据 | 11 | 11 |

- 与已有封装的关系：`find_quality_trim_bounds()`、`calculate_expected_errors()`（cutadapt 的 `qualtrim`）、`extract_read_barcodes()` 与 `map_barcodes_to_whitelist()`（umi-tools）已在 `dnakit.sequencing_qc` 中，本批不重复封装。`find_adapter_match()` 与 `trim_fastq_adapters()` 包装的是 cutadapt 的 `BackAdapter`、`AdapterCutter` 等类，清单中没有对应的函数条目。
- pyslow5：清单中没有函数条目。pyslow5 1.5.0 是单个扩展模块，没有子模块，公开接口是类 `pyslow5.Open`（方法有 `seq_reads`、`get_read_list`、`write_record` 等）。已有的 `read_slow5_signal()`、`summarize_slow5_signals()` 与 `write_slow5_signals()` 位于 `dnakit.sequencing_qc`。**注意：pyslow5 1.5.0 遇到损坏的 SLOW5/BLOW5 文件时可能使解释器崩溃**，已有封装会先检查文件头的魔数再打开文件。
- 依赖：umi-tools 的模块在导入时加载 pysam、matplotlib、scipy 与 pybktree，`umi_tools.dedup` 还加载 pandas，而 `bio-umi-tools` 只列出 umi-tools 与 regex。本批函数在导入上游模块前先检查这几个包，缺少时给出的安装提示是 `bio-umi-tools`，因此需要同步补全该额外依赖。
- 上游问题（只记录，未修复）：
  1. `umi_tools.umi_methods.edit_distance` 对长度不同的输入抛出 `AttributeError`（上游使用了 NumPy 2 已移除的 `np.Inf`）。本批只在等长 UMI 上调用它。
  2. `extractSeqAndQuals` 在 `retain_umi=True` 时把所有保留碱基的质量写入 `umi_qualities`，而不只是 UMI 碱基。本函数原样返回。
  3. 白名单与拐点函数要求 `collections.Counter`，普通字典会报 `AttributeError`。本函数会自动构造。
  4. `detect_bam_features` 通过 pysam 打开的文件不会显式关闭，只能等待垃圾回收。
  5. `get_gene_count_tab` 要求同一基因的行连续，空输入会产生一个 `None` 基因。本函数会合并分组，并忽略 `None`。
  6. `checkError` 找到两个匹配后就停止，所以 `matches` 不是完整列表。本函数在结果中给出 `ambiguous` 标记。
