# 测序质控适配器、UMI 与 SLOW5 信号后端

这一组函数是 `dnakit.sequencing_qc` 中基于可选依赖的适配器。适配器比对、FASTQ 适配器修剪和质量修剪使用 [cutadapt](https://github.com/marcelm/cutadapt) 的 Python API，FASTQ 读写使用它的依赖 [dnaio](https://github.com/marcelm/dnaio)；UMI 聚类、条码提取和白名单校正使用 [umi-tools](https://github.com/CGATOxford/UMI-tools)；原始信号的读写使用 [pyslow5](https://github.com/hasindu2008/slow5lib)。所有函数只处理显式给出的序列、文件和记录，不下载参考数据，也不调用任何命令行。它们计算适配器位置、修剪结果、质量修剪区间、期望错误数、UMI 簇、提取出的条码和白名单映射，以及 SLOW5/BLOW5 信号的统计、读取和写入；不做碱基识别或比对。

依赖都是可选的，只有调用对应函数时才会导入。缺少包时会报错并给出安装命令。

## 示例数据

下面的 FASTQ 只有 3 条 read，用于演示适配器修剪。`r1` 的 3' 端含有适配器 `AGATCGGAAGAGC`，`r3` 整条都是适配器，`r2` 不含适配器。SLOW5 信号数据由示例代码自己生成。

**示例文件：** `reads.fastq`

```text
@r1
ACGTTTACGTAGATCGGAAGAGCTTT
+
IIIIIIIIIIIIIIIIIIIIIIIIII
@r2
GGGGGGGGGGGGGGGG
+
IIIIIIIIIIIIIIII
@r3
AGATCGGAAGAGCACGT
+
IIIIIIIIIIIIIIIII
```

## 1) `dnakit.sequencing_qc.find_adapter_match()`

- **作用：** 在一条 DNA 序列中定位一个适配器，并返回修剪后的序列。
- **处理方法：** 使用 cutadapt 的适配器比对器，支持 3' 适配器（`"3prime"`）、5' 适配器（`"5prime"`）和任意位置（`"anywhere"`）。`max_errors` 是允许的错误率，`min_overlap` 是最短重叠长度。坐标为 0-based 半开区间。
- **API：** `dnakit.sequencing_qc.find_adapter_match(sequence, adapter, *, position="3prime", max_errors=0.1, min_overlap=3)`。
- **输入：** DNA 序列和适配器序列，只允许 ACGT 与 IUPAC 简并碱基；`max_errors` 须在 [0, 1) 内。
- **输出：** 含 `matched`、`start`、`stop`、`errors`、`score`、`trimmed_sequence` 和 `provenance` 的字典。未匹配时坐标为 `None`，修剪序列与输入相同。
- **示例代码：**

```python
from dnakit.sequencing_qc import find_adapter_match

result = find_adapter_match("ACGTTTACGTAGATCGGAAGAGCTTT", "AGATCGGAAGAGC")
print(result["matched"], result["start"], result["stop"], result["errors"])
print(result["trimmed_sequence"])
```

- **示例输出：**

```text
True 10 23 0
ACGTTTACGT
```

## 2) `dnakit.sequencing_qc.trim_fastq_adapters()`

- **作用：** 对 FASTQ 文件中的每条 read 修剪同一个适配器，并把保留的 read 写入新文件。
- **处理方法：** 用 cutadapt 的 `AdapterCutter` 逐条修剪，用 dnaio 读写 FASTQ。修剪后长度小于 `minimum_length` 的 read 被丢弃并计数。
- **API：** `dnakit.sequencing_qc.trim_fastq_adapters(source, output, adapter, *, position="3prime", max_errors=0.1, min_overlap=3, minimum_length=1, max_reads=1000000)`。
- **输入：** FASTQ 路径 `source`；`output` 必须是尚不存在的文件；`max_reads` 是显式上限，超过时报错并删除未完成的输出。
- **输出：** `reads`、`reads_with_adapter`、`reads_written`、`reads_too_short`、`bases_removed`（被删除的碱基总数）和 `provenance`。
- **示例代码：**

```python
import tempfile
from pathlib import Path
from dnakit.sequencing_qc import trim_fastq_adapters

with tempfile.TemporaryDirectory() as directory:
    result = trim_fastq_adapters("reads.fastq", Path(directory) / "trimmed.fastq", "AGATCGGAAGAGC")
    print({key: result[key] for key in ("reads", "reads_with_adapter", "reads_written", "reads_too_short", "bases_removed")})
```

- **示例输出：**

```text
{'reads': 3, 'reads_with_adapter': 2, 'reads_written': 2, 'reads_too_short': 1, 'bases_removed': 33}
```

## 3) `dnakit.sequencing_qc.find_quality_trim_bounds()`

- **作用：** 返回质量修剪后保留的区间。
- **处理方法：** 使用 cutadapt 的 `quality_trim_index`（BWA 算法）。从两端累计 `cutoff - Q`，取累计和最大的位置作为修剪点。
- **API：** `dnakit.sequencing_qc.find_quality_trim_bounds(qualities, *, cutoff_front=0, cutoff_back=20, phred_offset=33)`。
- **输入：** Phred 质量字符串；`cutoff_front` 和 `cutoff_back` 在 [0, 93] 内；`phred_offset` 为 33 或 64。
- **输出：** `length`、`start`、`stop`（0-based 半开区间）、`retained_length` 和 `provenance`。
- **示例代码：**

```python
from dnakit.sequencing_qc import find_quality_trim_bounds

result = find_quality_trim_bounds("IIII!!!!", cutoff_back=20)
print(result["start"], result["stop"], result["retained_length"])
```

- **示例输出：**

```text
0 4 4
```

## 4) `dnakit.sequencing_qc.calculate_expected_errors()`

- **作用：** 计算一条 read 的期望错误数。
- **处理方法：** 每个碱基贡献 10^(−Q/10)，求和由 cutadapt 的 `expected_errors` 完成。
- **API：** `dnakit.sequencing_qc.calculate_expected_errors(qualities, *, phred_offset=33)`。
- **输入：** Phred 质量字符串，`phred_offset` 为 33 或 64。
- **输出：** `length`、`expected_errors` 和 `provenance`。
- **示例代码：**

```python
from dnakit.sequencing_qc import calculate_expected_errors

result = calculate_expected_errors("++II")
print(round(result["expected_errors"], 4))
```

- **示例输出：**

```text
0.2002
```

## 5) `dnakit.sequencing_qc.write_slow5_signals()`

- **作用：** 把原始信号记录写入新的 SLOW5（`.slow5`）或 BLOW5（`.blow5`）文件。
- **处理方法：** 用 pyslow5 的写入接口写入单个 read group。所有记录必须共享同一个采样率，它被写入 `sample_frequency` 头字段。文件格式由输出文件的后缀决定。
- **API：** `dnakit.sequencing_qc.write_slow5_signals(records, output, *, max_reads=100000)`。
- **输入：** 记录列表。每条记录含 `read_id`、`signal`（int16 范围内的整数列表）、`digitisation`、`offset`、`range` 和 `sampling_rate`。`output` 必须尚不存在。
- **输出：** `reads_written`、`signal_samples_written`、`sampling_rate`、`output` 和 `provenance`。
- **示例代码：**

```python
import tempfile
from pathlib import Path
from dnakit.sequencing_qc import write_slow5_signals

common = {"digitisation": 8192.0, "offset": 10.0, "range": 1443.0, "sampling_rate": 4000.0}
records = [
    {"read_id": "r1", "signal": [500, 510, 520, 530], **common},
    {"read_id": "r2", "signal": [100, 200, 300], **common},
]
with tempfile.TemporaryDirectory() as directory:
    result = write_slow5_signals(records, Path(directory) / "signals.slow5")
    print(result["reads_written"], result["signal_samples_written"])
```

- **示例输出：**

```text
2 7
```

## 6) `dnakit.sequencing_qc.summarize_slow5_signals()`

- **作用：** 统计 SLOW5 或 BLOW5 文件中的 read 数、信号点数和采样率。
- **处理方法：** 顺序扫描一遍文件（`seq_reads`），不在输入文件旁边生成索引。read 数超过 `max_reads` 时报错。
- **API：** `dnakit.sequencing_qc.summarize_slow5_signals(source, *, max_reads=100000)`。
- **输入：** `.slow5` 或 `.blow5` 文件路径。后缀不符或文件头魔数不符时直接拒绝。
- **输出：** `reads`、`read_groups`、`signal_samples`、`min_signal_samples`、`mean_signal_samples`、`max_signal_samples`、`sampling_rates`、`header_sample_frequency` 和 `provenance`。
- **示例代码：**

```python
import tempfile
from pathlib import Path
from dnakit.sequencing_qc import summarize_slow5_signals, write_slow5_signals

common = {"digitisation": 8192.0, "offset": 10.0, "range": 1443.0, "sampling_rate": 4000.0}
records = [
    {"read_id": "r1", "signal": [500, 510, 520, 530], **common},
    {"read_id": "r2", "signal": [100, 200, 300], **common},
]
with tempfile.TemporaryDirectory() as directory:
    path = Path(directory) / "signals.slow5"
    write_slow5_signals(records, path)
    summary = summarize_slow5_signals(path)
    print(summary["reads"], summary["signal_samples"], summary["mean_signal_samples"], summary["sampling_rates"])
```

- **示例输出：**

```text
2 7 3.5 [4000.0]
```

## 7) `dnakit.sequencing_qc.read_slow5_signal()`

- **作用：** 读取一条 read 的原始信号，可选转换为 pA。
- **处理方法：** 顺序扫描直到找到 `read_id`，不建立索引。`pA=True` 时由 pyslow5 按 `(raw + offset) × range / digitisation` 转换，结果保留 7 位有效数字。
- **API：** `dnakit.sequencing_qc.read_slow5_signal(source, read_id, *, pA=False, max_reads=100000, max_samples=100000)`。
- **输入：** 文件路径和不含空白的 `read_id`。`max_samples` 限制单条 read 的信号长度，超过时报错。
- **输出：** `read_id`、`read_group`、`units`（`"raw DAC"` 或 `"pA"`）、`length`、`digitisation`、`offset`、`range`、`sampling_rate`、`signal` 和 `provenance`。
- **示例代码：**

```python
import tempfile
from pathlib import Path
from dnakit.sequencing_qc import read_slow5_signal, write_slow5_signals

common = {"digitisation": 8192.0, "offset": 10.0, "range": 1443.0, "sampling_rate": 4000.0}
records = [
    {"read_id": "r1", "signal": [500, 510, 520, 530], **common},
    {"read_id": "r2", "signal": [100, 200, 300], **common},
]
with tempfile.TemporaryDirectory() as directory:
    path = Path(directory) / "signals.slow5"
    write_slow5_signals(records, path)
    raw = read_slow5_signal(path, "r2")
    picoamps = read_slow5_signal(path, "r1", pA=True)
    print(raw["signal"], picoamps["units"], picoamps["signal"])
```

- **示例输出：**

```text
[100, 200, 300] pA [89.83521, 91.59668, 93.35815, 95.11963]
```

## 8) `dnakit.sequencing_qc.cluster_umi_counts()`

- **作用：** 按 Hamming 距离把 UMI 聚类，每个簇的代表是计数最高的 UMI。
- **处理方法：** 调用 umi-tools 的 `UMIClusterer`。`directional` 方法要求父 UMI 的计数不小于子 UMI 计数的 2 倍减 1 才会合并；`adjacency` 只看距离；`cluster` 合并连通分量；`unique` 不合并。距离按 umi-tools 的 Hamming 距离计算。
- **API：** `dnakit.sequencing_qc.cluster_umi_counts(umi_counts, *, method="directional", max_distance=1, max_umis=10000)`。
- **输入：** `umi_counts` 是 `{UMI: 正整数计数}`；所有 UMI 等长且只含 A、C、G、T、N；`method` 取 `"directional"`、`"adjacency"`、`"cluster"` 或 `"unique"`；`max_distance` 在 0 到 5 之间。
- **输出：** `clusters`（每个簇是 UMI 列表，代表在第一位）、`representatives`、`cluster_counts`（每个簇的计数之和）、`umis`、`umi_length` 和 `provenance`。
- **示例代码：**

```python
from dnakit.sequencing_qc import cluster_umi_counts

result = cluster_umi_counts({"AAAA": 10, "AAAT": 3, "CCCC": 5, "AAAC": 1})
print(result["clusters"])
print(result["cluster_counts"])
```

- **示例输出：**

```text
[['AAAA', 'AAAT', 'AAAC'], ['CCCC']]
[14, 5]
```

## 9) `dnakit.sequencing_qc.extract_read_barcodes()`

- **作用：** 用正则表达式从一条 read 中提取细胞条码和 UMI，并返回去除条码后的序列。
- **处理方法：** 调用 umi-tools 的 `ExtractBarcodes`。命名组 `cell_<n>` 是细胞条码，`umi_<n>` 是 UMI，`discard_<n>` 是要丢弃的碱基（`discard=True` 时）。正则从 read 起点匹配，每个条码组都必须参与匹配。umi-tools 不返回细胞条码的质量值，因此本函数只返回 UMI 的质量值。
- **API：** `dnakit.sequencing_qc.extract_read_barcodes(sequence, qualities, pattern, *, discard=True)`。
- **输入：** DNA 序列、长度相同的 Phred+33 质量字符串，以及不超过 1000 个字符的正则表达式字符串。
- **输出：** `matched`、`cell_barcode`、`umi`、`umi_qualities`、`read_sequence`、`read_qualities` 和 `provenance`。未匹配时两个条码为空字符串，read 保持不变。
- **示例代码：**

```python
from dnakit.sequencing_qc import extract_read_barcodes

read = "ACGTTTAGCCATTGGGAAA"
pattern = r"^(?P<cell_1>.{4})(?P<umi_1>.{3})(?P<discard_1>.{2})"
result = extract_read_barcodes(read, "I" * len(read), pattern)
print(result["cell_barcode"], result["umi"], result["umi_qualities"], result["read_sequence"])
```

- **示例输出：**

```text
ACGT TTA III CATTGGGAAA
```

## 10) `dnakit.sequencing_qc.map_barcodes_to_whitelist()`

- **作用：** 把观察到的细胞条码校正到白名单上，只接受唯一的近似匹配。
- **处理方法：** 调用 umi-tools 的 `getErrorCorrectMapping`，按 Hamming 距离匹配，阈值为 1 到 3。一个非白名单条码只有在恰好一个白名单条码落在阈值内时才会被分配；没有候选或有多个候选的条码进入 `unassigned_barcodes`。已经在白名单中的条码计入 `exact_matches`。
- **API：** `dnakit.sequencing_qc.map_barcodes_to_whitelist(observed, whitelist, *, threshold=1, max_barcodes=100000)`。
- **输入：** 观察到的 DNA 条码列表和白名单列表。所有条码等长且只含 A、C、G、T。
- **输出：** `corrected`（白名单条码到被校正的观察条码列表的映射）、`corrected_barcodes`、`exact_matches`、`observed_unique`、`unassigned_barcodes` 和 `provenance`。
- **示例代码：**

```python
from dnakit.sequencing_qc import map_barcodes_to_whitelist

result = map_barcodes_to_whitelist(
    ["AAAT", "AAAG", "CCCA", "GGGG", "TTTT"], ["AAAA", "CCCC", "GGGG"]
)
print(result["corrected"])
print(result["unassigned_barcodes"])
```

- **示例输出：**

```text
{'AAAA': ['AAAG', 'AAAT'], 'CCCC': ['CCCA']}
['TTTT']
```

## 安装可选依赖

```bash
python -m pip install "dnakit[bio-cutadapt]"    # cutadapt 与 dnaio，适配器与 FASTQ 修剪、质量修剪
python -m pip install "dnakit[bio-pyslow5]"     # pyslow5 与 numpy，SLOW5/BLOW5 函数
python -m pip install "dnakit[bio-umi-tools]"   # umi-tools 与 regex，UMI 聚类、条码提取、白名单校正
```

各包的项目地址：[cutadapt](https://github.com/marcelm/cutadapt)、[dnaio](https://github.com/marcelm/dnaio)、[umi-tools](https://github.com/CGATOxford/UMI-tools)、[pyslow5](https://github.com/hasindu2008/slow5lib)。

## 验证方式

测试文件是 `others/tests/unit/test_sequencing_backends.py`，共 32 项，主要检查：

- 适配器坐标、错误数和修剪结果与手工计算一致；3' 与 5' 适配器、`anywhere` 位置和错误率的作用都经过检查。
- 适配器坐标与直接调用 cutadapt 的 `BackAdapter.match_to` 的结果一致。
- FASTQ 修剪的读数、丢弃数和删除的碱基数与手工计算一致，输出文件内容逐行核对；超过 `max_reads` 时报错并删除未完成的输出，已存在的输出文件不会被覆盖。
- 质量修剪区间、期望错误数与手工公式一致（`10^(-Q/10)`）。
- SLOW5 的读数、信号点数和 pA 转换与手工公式一致；统计结果与直接调用 pyslow5 顺序扫描的结果一致；`.blow5` 输出可以往返读取，并且不生成索引文件。
- 缺少 cutadapt、pyslow5 或 umi-tools 时给出安装提示；越界的信号值、重复的 `read_id`、不一致的采样率、非 SLOW5 文件头都会被拒绝。
- UMI 聚类的结果与手工按 `directional` 规则计算的结果一致：计数 10 的 `AAAA` 合并计数 3 的 `AAAT` 和计数 1 的 `AAAC`。相同的两个 UMI（计数均为 3）在 `directional` 下不合并，在 `adjacency` 下合并。结果与直接调用 umi-tools 的 `UMIClusterer` 一致。
- 条码提取的结果与手工切分一致，并与直接调用 umi-tools 的 `ExtractBarcodes` 一致；未匹配的 read 保持不变，缺少条码组和非法正则都会被拒绝。
- 白名单校正的结果与手工计算一致，并与直接调用 `getErrorCorrectMapping` 一致；有多个候选的条码不分配。

## 当前范围

这一批覆盖 cutadapt 的适配器比对、FASTQ 适配器修剪、质量修剪区间与期望错误数，umi-tools 的 UMI 聚类、正则条码提取和白名单校正，以及 pyslow5 的 SLOW5/BLOW5 写入、统计和单条读取。

以下内容未包含在本批中：

- **cutadapt 的其他功能：** poly-A 修剪、双端修剪和命令行调用未包装。`run_cutadapt` 是 `dnakit.sequencing_qc` 中已有的命令行适配器，与本批的 Python API 函数互不替代。

**pyslow5 的崩溃风险：** pyslow5 1.5.0 读取损坏的 SLOW5/BLOW5 文件时，可能直接让 Python 解释器段错误退出。包装器会先检查后缀和文件头魔数，拒绝空文件、文本文件和后缀不符的文件；但魔数正确的文件仍可能因文件头或内容损坏而崩溃，实际测试中只要文件头行损坏就会发生。损坏的记录也可能被静默跳过，统计结果会偏小。只打开可信来源的 SLOW5/BLOW5 文件。

pyslow5 的随机访问接口（按 `read_id` 直接查找）会在输入文件旁边生成 `.slow5.idx` 索引文件。本包装器全部使用顺序扫描，因此不会产生这个文件。
