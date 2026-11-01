# GimmeMotifs FASTA 工具、统计、阈值、匹配计数与基序数据库

这一组函数位于 `dnakit.regulation`，封装 [GimmeMotifs](https://github.com/vanheeringen-lab/gimmemotifs) 0.18.4 的 FASTA、峰输入、统计、阈值、匹配计数和基序数据库例程。每个函数对本地文件或显式列表调用 GimmeMotifs 例程，返回可 JSON 序列化的字典，含 `status`、结果和 `provenance`。基因组参数是 FASTA 路径（不通过 genomepy 下载）；输出文件必须不存在。两个基序数据库函数需要联网。

GimmeMotifs 是可选依赖（`pip install "dnakit[bio-gimmemotifs]"`；BED 与 narrowPeak 输入准备还需要 `genomepy`）。示例从 `dnakit.regulation.gimmemotifs_p28_data_functions` 导入；区域包接入后，同名函数也可从 `dnakit.regulation` 导入。

覆盖的例程有 `Fasta.hardmask`、`Fasta.get_random`、`utils.as_seqdict`、`stats.rank_motifs`、`stats.write_stats`、`Scanner.set_threshold`（及 `gimme threshold` 命令）、`Scanner.count` 与 `Scanner.total_count`、`utils.narrowpeak_to_bed`、`motif.denovo.prepare_denovo_input_fa`、`_bed`、`_narrowpeak`，以及 JASPAR、HOMER、HOCOMOCO、ENCODE、SwissRegulon 和 IMAGE 基序数据库的 `download` 与 `annotate_factors` 方法。已有函数覆盖的行有 `rocmetrics`（`motif_score_metric`）、`preprocessing.peakfile2bedfile`（`normalize_peak_widths`）和 `utils.get_jaspar_motif_info`（`gimmemotifs_jaspar_motif_info`）。

## 示例数据

第一个示例读取含两条记录的软屏蔽 FASTA 文件 `soft.fa`。其余示例用一个小的确定性生成器（`make_fasta`）把输入写入临时目录，或直接传入。

**示例文件：** `soft.fa`

```text
>a
ACGTacgtnNAC
>b
ggggCCCC
```

## 1) `dnakit.regulation.gimmemotifs_hardmask_fasta()`

- **作用：** 把 FASTA 文件中的小写（软屏蔽）碱基替换成 N。
- **处理方法：** 用 `gimmemotifs.fasta.Fasta` 读取文件，调用 `Fasta.hardmask`（小写的 a、c、g、t、n 变成 `N`，大写字母保持不变），再用 `Fasta.writefasta` 写出记录，不折行。输出路径必须是新的。
- **API：** `dnakit.regulation.gimmemotifs_hardmask_fasta(input_fasta, output_fasta)`。
- **输入：** 以 `>` 开头的 FASTA 文件和一个新的输出路径。
- **输出：** `masked_bases`、`n_sequences`，以及 `outputs` 中写出的文件。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.regulation.gimmemotifs_p28_data_functions import gimmemotifs_hardmask_fasta

out = Path(tempfile.mkdtemp()) / "hard.fa"
result = gimmemotifs_hardmask_fasta("soft.fa", out)
print(result["masked_bases"], result["n_sequences"])
print(out.read_text().strip())
```

- **示例输出：**

```text
9 2
>a
ACGTNNNNNNAC
>b
NNNNCCCC
```

## 2) `dnakit.regulation.gimmemotifs_random_fasta_subset()`

- **作用：** 从 FASTA 文件中随机抽取整条记录或定长窗口。
- **处理方法：** 调用 `Fasta.get_random`。不给 `length` 时无放回地抽取 `n_sequences` 条记录并保留名称。给出 `length` 时按随机顺序遍历记录，从每条足够长的记录中取一个长度为 `length` 的随机窗口，命名为 `random1`、`random2`……；足够长的记录少于 `n_sequences` 条时失败。`seed` 为本次调用设置 Python 随机数种子（之后恢复原状态）。子集可写入新的 FASTA 文件。
- **API：** `dnakit.regulation.gimmemotifs_random_fasta_subset(input_fasta, n_sequences, *, length=None, seed=None, output_fasta=None)`。
- **输入：** FASTA 文件、序列数、可选的窗口 `length`、`seed` 与 `output_fasta`。
- **输出：** `sequences`（名称到序列）、`n_sequences`，以及 `outputs` 中可选的文件。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.regulation.gimmemotifs_p28_data_functions import gimmemotifs_random_fasta_subset

def make_fasta(path, count, length, seed):
    state, rows = seed, []
    for index in range(count):
        bases = []
        for _ in range(length):
            state = (state * 1103515245 + 12345) % 2**31
            bases.append("ACGT"[(state >> 16) % 4])
        rows.append(f">r{index}\n{''.join(bases)}\n")
    path.write_text("".join(rows))
    return path

fasta = make_fasta(Path(tempfile.mkdtemp()) / "r.fa", 20, 60, 7)
records = gimmemotifs_random_fasta_subset(fasta, 4, seed=3)
print(len(records["sequences"]), all(len(s) == 60 for s in records["sequences"].values()))
again = gimmemotifs_random_fasta_subset(fasta, 4, seed=3)
print(records["sequences"] == again["sequences"])
windows = gimmemotifs_random_fasta_subset(fasta, 3, length=25, seed=3)
print(list(windows["sequences"]), [len(s) for s in windows["sequences"].values()])
```

- **示例输出：**

```text
4 True
True
['random1', 'random2', 'random3'] [25, 25, 25]
```

## 3) `dnakit.regulation.gimmemotifs_sequence_dict()`

- **作用：** 把序列读入以序列名为键的字典。
- **处理方法：** 调用 `gimmemotifs.utils.as_seqdict`。FASTA 文件以完整标题行作键；DNA 字符串列表的键为 `1`、`2`……；BED 文件或 `chrom:start-end` 区域列表从 `genome_fasta`（按 FASTA 路径作为 genomepy 基因组读取）中截取。给出 `minsize` 时每条序列至少要这么长，否则调用失败。
- **API：** `dnakit.regulation.gimmemotifs_sequence_dict(source, *, genome_fasta=None, minsize=None)`。
- **输入：** FASTA、BED 或区域文件的路径，或 DNA 字符串/区域列表；可选的 `genome_fasta` 与 `minsize`。
- **输出：** `sequences`（名称到序列）与 `n_sequences`。
- **示例代码：**

```python
from dnakit.regulation.gimmemotifs_p28_data_functions import gimmemotifs_sequence_dict

print(gimmemotifs_sequence_dict("soft.fa")["sequences"])
print(gimmemotifs_sequence_dict(["ACGT", "GGAA"])["sequences"])
try:
    gimmemotifs_sequence_dict("soft.fa", minsize=10)
except Exception as error:
    print(type(error).__name__, "shorter than 10" in str(error))
```

- **示例输出：**

```text
{'a': 'ACGTacgtnNAC', 'b': 'ggggCCCC'}
{'1': 'ACGT', '2': 'GGAA'}
BackendExecutionError True
```

## 4) `dnakit.regulation.gimmemotifs_rank_motif_stats()`

- **作用：** 按多个富集指标的平均名次给基序排名。
- **处理方法：** 对形如 `{motif: {background: {metric: value}}}` 的表（基序富集统计的输出形式）调用 `gimmemotifs.stats.rank_motifs`。对每个指标，先把各背景上的基序值取平均，再用 `scipy.stats.rankdata` 排名（1 为最低值，并列取平均名次）；结果是这些名次的平均值，数值越高基序越好。
- **API：** `dnakit.regulation.gimmemotifs_rank_motif_stats(stats, *, metrics=("roc_auc", "recall_at_fdr"))`。
- **输入：** 嵌套的 `stats` 映射（每个基序的背景和指标相同）以及要使用的 `metrics` 列表。
- **输出：** 每个基序的 `mean_rank`、`order`（最好的在前）、`metrics` 与 `backgrounds`。
- **示例代码：**

```python
from dnakit.regulation.gimmemotifs_p28_data_functions import gimmemotifs_rank_motif_stats

stats = {
    "m_1": {"bg": {"roc_auc": 0.9, "recall_at_fdr": 0.5}},
    "m_2": {"bg": {"roc_auc": 0.6, "recall_at_fdr": 0.7}},
    "m_3": {"bg": {"roc_auc": 0.7, "recall_at_fdr": 0.1}},
}
result = gimmemotifs_rank_motif_stats(stats)
print(result["mean_rank"], result["order"])
print(gimmemotifs_rank_motif_stats(stats, metrics=["recall_at_fdr"])["order"])
```

- **示例输出：**

```text
{'m_1': 2.5, 'm_2': 2.0, 'm_3': 1.5} ['m_1', 'm_2', 'm_3']
['m_2', 'm_1', 'm_3']
```

## 5) `dnakit.regulation.gimmemotifs_write_motif_stats()`

- **作用：** 写出基序统计表，每个背景一个制表符分隔的文件。
- **处理方法：** 调用 `gimmemotifs.stats.write_stats`。`output_template` 是含一个 `{}` 的路径，`{}` 被替换为各背景名。每个文件含可选的 `header` 文本、`Motif` 列和按字母顺序排列的各指标列。与上游一致，基序名写出时去掉最后一个 `_` 之后的部分（`m_1` 变成 `m`）。不会覆盖已有文件。
- **API：** `dnakit.regulation.gimmemotifs_write_motif_stats(stats, output_template, *, header=None)`。
- **输入：** 与 `gimmemotifs_rank_motif_stats` 相同的 `stats` 映射、含一个 `{}` 的模板和可选的表头。
- **输出：** `outputs`（背景到文件路径）、`backgrounds` 与 `n_motifs`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.regulation.gimmemotifs_p28_data_functions import gimmemotifs_write_motif_stats

stats = {
    "m_1": {"bg": {"roc_auc": 0.9, "recall_at_fdr": 0.5}},
    "m_2": {"bg": {"roc_auc": 0.6, "recall_at_fdr": 0.7}},
}
folder = Path(tempfile.mkdtemp())
result = gimmemotifs_write_motif_stats(stats, str(folder / "stats.{}.txt"), header="# demo\n")
print(list(result["outputs"]))
print((folder / "stats.bg.txt").read_text().strip())
```

- **示例输出：**

```text
['bg']
# demo
Motif	recall_at_fdr	roc_auc
m	0.5	0.9
m	0.7	0.6
```

## 6) `dnakit.regulation.gimmemotifs_motif_thresholds()`

- **作用：** 在背景序列上按目标假阳性率给出基序分数阈值。
- **处理方法：** 用 GimmeMotifs 扫描器在背景 FASTA 上扫描 PFM 文件中的每个基序（每条序列取最佳命中，含两条链），对每个基序取这些最佳分数的 `1 - fpr` 分位数作为阈值（`Scanner.get_gc_thresholds`，即 `Scanner.set_threshold` 与 `gimme threshold` 背后的步骤）。GimmeMotifs 会对背景分数重采样，因此需要可重复结果时请传 `seed`。不做 GC 校正。阈值为 `None` 表示分位数已达基序最大分数，没有序列能通过。阈值同时以其占基序最小到最大分数之间距离的比例给出。
- **API：** `dnakit.regulation.gimmemotifs_motif_thresholds(motifs_pfm, background_fasta, *, fpr=0.01, seed=None)`。
- **输入：** PFM 文件、背景 FASTA 文件、严格介于 0 与 1 之间的 `fpr` 和可选的整数 `seed`。
- **输出：** 每个基序的 `thresholds`（`threshold`、`fraction`、`max_score`）和 `n_background_sequences`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.regulation.gimmemotifs_p28_data_functions import gimmemotifs_motif_thresholds

def make_fasta(path, count, length, seed):
    state, rows = seed, []
    for index in range(count):
        bases = []
        for _ in range(length):
            state = (state * 1103515245 + 12345) % 2**31
            bases.append("ACGT"[(state >> 16) % 4])
        rows.append(f">r{index}\n{''.join(bases)}\n")
    path.write_text("".join(rows))
    return path

folder = Path(tempfile.mkdtemp())
background = make_fasta(folder / "bg.fa", 60, 100, 11)
(folder / "m.pfm").write_text((
    ">weak\n6 2 1 1\n1 6 2 1\n1 1 6 2\n2 1 1 6\n3 3 2 2\n"
    ">acgt\n10 0 0 0\n0 10 0 0\n0 0 10 0\n0 0 0 10\n"
))
result = gimmemotifs_motif_thresholds(folder / "m.pfm", background, fpr=0.5, seed=1)
print(result["n_background_sequences"], list(result["thresholds"]))
weak = result["thresholds"]["weak"]
print(round(weak["threshold"], 4), round(weak["fraction"], 4), round(weak["max_score"], 4))
print(result["thresholds"]["acgt"]["threshold"])
```

- **示例输出：**

```text
60 ['weak', 'acgt']
2.6183 0.8544 3.7085
-0.40047363884467657
```

## 7) `dnakit.regulation.gimmemotifs_count_motif_matches()`

- **作用：** 按序列和总体统计基序匹配数。
- **处理方法：** 使用 GimmeMotifs 的 `Scanner`：`Scanner.count` 给出每个基序在每条序列中的匹配数，`Scanner.total_count` 给出所有序列的总和。匹配指分数超过基序截断值的命中，`scan_rc` 为真时包含两条链，每个基序在每条序列中最多计 `nreport` 个命中。默认截断值是 `threshold`，即占基序最小到最大分数之间距离的比例（0.95 是无基因组时 GimmeMotifs 的默认值）。给出 `fpr` 与 `background_fasta` 时，截断值改为占背景序列比例 `fpr` 所能达到的分数（见 `gimmemotifs_motif_thresholds`）。
- **API：** `dnakit.regulation.gimmemotifs_count_motif_matches(sequences_fasta, motifs_pfm, *, threshold=0.95, background_fasta=None, fpr=None, nreport=100, scan_rc=True, seed=None)`。
- **输入：** FASTA 文件、PFM 文件，以及 `threshold` 或成对的 `fpr` 与 `background_fasta`。
- **输出：** `counts`（序列到基序到计数）、`total`（基序到计数）、`motifs` 与截断模式。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.regulation.gimmemotifs_p28_data_functions import gimmemotifs_count_motif_matches

def make_fasta(path, count, length, seed):
    state, rows = seed, []
    for index in range(count):
        bases = []
        for _ in range(length):
            state = (state * 1103515245 + 12345) % 2**31
            bases.append("ACGT"[(state >> 16) % 4])
        rows.append(f">r{index}\n{''.join(bases)}\n")
    path.write_text("".join(rows))
    return path

folder = Path(tempfile.mkdtemp())
sequences = make_fasta(folder / "seq.fa", 10, 100, 21)
(folder / "m.pfm").write_text((
    ">weak\n6 2 1 1\n1 6 2 1\n1 1 6 2\n2 1 1 6\n3 3 2 2\n"
    ">acgt\n10 0 0 0\n0 10 0 0\n0 0 10 0\n0 0 0 10\n"
))
result = gimmemotifs_count_motif_matches(sequences, folder / "m.pfm", threshold=0.8)
print(result["motifs"], result["total"])
print(result["counts"]["r0"])
print(sum(sum(row.values()) for row in result["counts"].values()) == sum(result["total"].values()))
```

- **示例输出：**

```text
['weak', 'acgt'] {'weak': 41, 'acgt': 14}
{'weak': 4, 'acgt': 0}
True
```

## 8) `dnakit.regulation.gimmemotifs_narrowpeak_to_bed()`

- **作用：** 把 narrowPeak 文件转换成四列 BED 文件，可选地转换成以峰顶为中心的定长窗口。
- **处理方法：** 调用 `gimmemotifs.utils.narrowpeak_to_bed`。`size=0` 时每个峰保持原坐标。`size` 为正时区间长 `size` 个碱基，起点为 `峰起点 + 峰顶 - size // 2`；峰顶取 narrowPeak 第 10 列，-1（无峰顶）时退回 `(end - start) // 2`。BED 的列为染色体、起点、终点和信号值（narrowPeak 第 7 列）。注释行、`track` 行和 `browser` 行被跳过。输出文件必须不存在。
- **API：** `dnakit.regulation.gimmemotifs_narrowpeak_to_bed(input_file, output_bed, *, size=0)`。
- **输入：** 10 列的 narrowPeak 文件、新的 BED 路径和 `size`。
- **输出：** `intervals`（`[染色体, 起点, 终点, 信号]`）、`n_intervals`，以及 `outputs` 中的文件。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.regulation.gimmemotifs_p28_data_functions import gimmemotifs_narrowpeak_to_bed

folder = Path(tempfile.mkdtemp())
peaks = folder / "p.narrowPeak"
peaks.write_text(
    "chr1\t100\t300\tp0\t0\t.\t12.5\t-1\t-1\t80\n"
    "chr2\t500\t900\tp1\t0\t.\t7.0\t-1\t-1\t-1\n"
)
print(gimmemotifs_narrowpeak_to_bed(peaks, folder / "plain.bed")["intervals"])
print(gimmemotifs_narrowpeak_to_bed(peaks, folder / "summit.bed", size=60)["intervals"])
```

- **示例输出：**

```text
[['chr1', 100, 300, '12.5'], ['chr2', 500, 900, '7.0']]
[['chr1', 150, 210, '12.5'], ['chr2', 670, 730, '7.0']]
```

## 9) `dnakit.regulation.gimmemotifs_prepare_denovo_input()`

- **作用：** 准备 de novo 基序预测的输入文件：预测集、验证集和定位集。
- **处理方法：** `input_type="fasta"` 调用 `prepare_denovo_input_fa`：记录被分成 `prediction.fa`（占比 `fraction`，至多 `abs_max` 条）和 `validation.fa`，`localization.fa` 是验证集的副本。`"bed"` 调用 `prepare_denovo_input_bed`：区域以中点为中心调整为 `size` 个碱基（`input.bed`），拆成 `prediction.bed` 与 `validation.bed`，并从 `genome_fasta` 截取成对应的 FASTA；定位用 FASTA 扩展到 `lsize` 个碱基。`"narrowpeak"` 调用 `prepare_denovo_input_narrowpeak`，先把峰转换成以峰顶为中心、长 `size` 的 BED（`<名称>_to.bed`），再按 BED 处理。未设置的 `size`、`fraction`、`abs_max`、`lsize` 取 GimmeMotifs 默认值（200、0.2、1000、500）。`seed` 为随机拆分设置种子。
- **API：** `dnakit.regulation.gimmemotifs_prepare_denovo_input(input_file, output_dir, *, input_type="fasta", genome_fasta=None, size=None, fraction=None, abs_max=None, lsize=None, use_strand=False, seed=None)`。
- **输入：** 输入文件、新建或为空的 `output_dir`、`input_type`，BED 与 narrowPeak 还需 `genome_fasta`。
- **输出：** `outputs`（文件名到路径）、`records_per_file` 与所用 `parameters`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.regulation.gimmemotifs_p28_data_functions import gimmemotifs_prepare_denovo_input

def make_fasta(path, count, length, seed):
    state, rows = seed, []
    for index in range(count):
        bases = []
        for _ in range(length):
            state = (state * 1103515245 + 12345) % 2**31
            bases.append("ACGT"[(state >> 16) % 4])
        rows.append(f">r{index}\n{''.join(bases)}\n")
    path.write_text("".join(rows))
    return path

folder = Path(tempfile.mkdtemp())
genome = make_fasta(folder / "genome.fa", 1, 5000, 3)
genome.write_text(genome.read_text().replace(">r0", ">chr1"))
peaks = folder / "peaks.bed"
peaks.write_text("".join(f"chr1\t{100 + 150 * i}\t{300 + 150 * i}\tp{i}\t{100 - i}\n" for i in range(10)))
fa = make_fasta(folder / "in.fa", 20, 60, 5)

result = gimmemotifs_prepare_denovo_input(fa, folder / "from_fasta", fraction=0.5, seed=1)
print(result["records_per_file"])
result = gimmemotifs_prepare_denovo_input(
    peaks, folder / "from_bed", input_type="bed", genome_fasta=genome, size=100, fraction=0.5, lsize=140, seed=1
)
print(result["records_per_file"])
print(len((folder / "from_bed" / "prediction.fa").read_text().splitlines()[1]),
      len((folder / "from_bed" / "localization.fa").read_text().splitlines()[1]))
```

- **示例输出：**

```text
{'localization.fa': 10, 'prediction.fa': 10, 'validation.fa': 10}
{'input.bed': 10, 'localization.fa': 5, 'prediction.bed': 5, 'prediction.fa': 5, 'validation.bed': 5, 'validation.fa': 5}
100 140
```

## 10) `dnakit.regulation.gimmemotifs_download_motif_database()`

- **作用：** 下载公共基序数据库，写成 GimmeMotifs PFM 格式并附上因子注释（联网）。
- **处理方法：** 创建所选数据库的 GimmeMotifs `MotifDb` 并调用其 `download` 方法。`database` 取 `jaspar`（`jaspar_version` 的全部 CORE 分组）、`homer`、`hocomoco`（v11 人和小鼠）、`encode`、`swissregulon` 或 `image`。每个 PFM 文件（含 `# Retrieved from` 与 `# Date` 头部行）与由 `annotate_factors` 生成的 `.motif2factors.txt` 文件一起写出。目录需要时会创建，文件与上游一样会被覆盖。Factorbook、CIS-BP 和 RSAT 的上游下载未提供，因为它们的源地址返回 404（2026-10-20 核查）。示例用本地文本代替网络请求，以便离线运行。
- **API：** `dnakit.regulation.gimmemotifs_download_motif_database(database, output_dir, *, jaspar_version="2022")`。
- **输入：** 数据库名、输出目录，JASPAR 另需四位数的发布年份。
- **输出：** `outputs`（所写文件的文件名到路径）、`date` 与 `database`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

import gimmemotifs.db as db

from dnakit.regulation.gimmemotifs_p28_data_functions import gimmemotifs_download_motif_database

HOMER = (
    ">ATGACTCATC_PU.1(ETS)/ThioMac-PU.1-ChIP-Seq(GSE21512)/Homer\t6.04\t-1782\t0\t9805\n"
    "0.419\t0.275\t0.277\t0.028\n0.001\t0.001\t0.001\t0.997\n0.010\t0.002\t0.965\t0.023\n"
    ">GGAAGTGGGG_Gfi1b(Zf)/Bone-Gfi1b-ChIP-Seq(GSE23725)/Homer\t6.1\t-12\t0\t9\n"
    "0.1\t0.1\t0.7\t0.1\n0.1\t0.7\t0.1\t0.1\n0.7\t0.1\t0.1\t0.1\n"
)


class Reply(list):
    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


# a local text stands in for the HOMER web server so that the example runs offline;
# remove the next line to download the real database
db.urlopen = lambda url: Reply(HOMER.encode().splitlines(keepends=True))

folder = Path(tempfile.mkdtemp())
result = gimmemotifs_download_motif_database("homer", folder)
print(sorted(result["outputs"]))
print((folder / "HOMER.pfm").read_text().count(">"))
print((folder / "HOMER.motif2factors.txt").read_text().strip())
```

- **示例输出：**

```text
['HOMER.motif2factors.txt', 'HOMER.pfm']
2
Motif	Factor	Evidence	Curated
ATGACTCATC_PU.1(ETS)/ThioMac-PU.1-ChIP-Seq(GSE21512)/Homer_6.04	PU.1	ChIP-seq	Y
GGAAGTGGGG_Gfi1b(Zf)/Bone-Gfi1b-ChIP-Seq(GSE23725)/Homer_6.1	Gfi1b	ChIP-seq	Y
```

## 11) `dnakit.regulation.gimmemotifs_annotate_motif_factors()`

- **作用：** 用基序数据库中的转录因子注释 PFM 文件里的基序（联网）。
- **处理方法：** 读取基序并调用所选 GimmeMotifs `MotifDb` 的 `annotate_factors`。`database` 取 `jaspar`（到 JASPAR 逐个查询 `MA` 矩阵）、`homer`、`hocomoco`（`hocomoco_group` 取 `HUMAN` 或 `MOUSE`，决定注释表）、`encode`、`swissregulon` 或 `image`。基序 ID 必须遵循该数据库自己的命名，例如 `gimmemotifs_download_motif_database` 写出的 PFM 文件；没有匹配的基序不出现在结果中。给出 `output_dir` 时，表还会写成 `<pfm 名>.motif2factors.txt`（列为 Motif、Factor、Evidence、Curated）。RSAT、CIS-BP 与 Factorbook 的注释未提供：前者在上游是空实现，后两者读取的地址返回 404（2026-10-20 核查）。示例用本地文本代替网络请求。
- **API：** `dnakit.regulation.gimmemotifs_annotate_motif_factors(database, motifs_pfm, *, hocomoco_group="HUMAN", output_dir=None)`。
- **输入：** 数据库名、PFM 文件、`hocomoco_group` 与可选的输出目录。
- **输出：** `annotation`（基序到 `[factor, evidence, curated]` 列表）、`n_motifs`、`n_annotated`，以及 `outputs` 中可选的文件。
- **示例代码：**

```python
import tempfile
from pathlib import Path

import gimmemotifs.db as db

from dnakit.regulation.gimmemotifs_p28_data_functions import gimmemotifs_annotate_motif_factors

HOMER = (
    ">ATGACTCATC_PU.1(ETS)/ThioMac-PU.1-ChIP-Seq(GSE21512)/Homer\t6.04\t-1782\t0\t9805\n"
    "0.419\t0.275\t0.277\t0.028\n0.001\t0.001\t0.001\t0.997\n0.010\t0.002\t0.965\t0.023\n"
    ">GGAAGTGGGG_Gfi1b(Zf)/Bone-Gfi1b-ChIP-Seq(GSE23725)/Homer\t6.1\t-12\t0\t9\n"
    "0.1\t0.1\t0.7\t0.1\n0.1\t0.7\t0.1\t0.1\n0.7\t0.1\t0.1\t0.1\n"
)


class Reply(list):
    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


# a local text stands in for the HOMER web server so that the example runs offline;
# remove the next line to download the real database
db.urlopen = lambda url: Reply(HOMER.encode().splitlines(keepends=True))

folder = Path(tempfile.mkdtemp())
(folder / "homer.pfm").write_text(
    ">ATGACTCATC_PU.1(ETS)/ThioMac-PU.1-ChIP-Seq(GSE21512)/Homer\n0.4 0.3 0.2 0.1\n0.1 0.1 0.1 0.7\n"
    ">NNNNNNNN_AP-1(bZIP)/ThioMac-PU.1-ChIP-Seq(GSE21512)/Homer\n0.4 0.3 0.2 0.1\n0.1 0.1 0.1 0.7\n"
)
result = gimmemotifs_annotate_motif_factors("homer", folder / "homer.pfm", output_dir=folder / "out")
print(result["n_motifs"], result["n_annotated"])
print(result["annotation"])
print(sorted(p.name for p in (folder / "out").iterdir()))
```

- **示例输出：**

```text
2 1
{'ATGACTCATC_PU.1(ETS)/ThioMac-PU.1-ChIP-Seq(GSE21512)/Homer': [['PU.1', 'ChIP-seq', 'Y']]}
['homer.motif2factors.txt']
```
