# GimmeMotifs 基序比较、矩阵转换与扫描

这一组函数位于 `dnakit.regulation`，封装 [GimmeMotifs](https://github.com/vanheeringen-lab/gimmemotifs) 0.18.4 的基序层面例程：基序比较、基序矩阵的转换与导出，以及对 FASTA 序列的基序扫描。基序用 PFM 表示，即每个位置一行、每行为 A、C、G、T 的计数的列表（也接受概率行，与 GimmeMotifs 一致）。每个函数调用同名的 GimmeMotifs 例程，返回可 JSON 序列化的字典，含 `status`、结果和 `provenance`（给出例程名与 GimmeMotifs 版本）。不下载任何数据。

GimmeMotifs 是可选依赖（`pip install "dnakit[bio-gimmemotifs]"`），只有调用函数时才会导入。示例从模块 `dnakit.regulation.gimmemotifs_p28_motif_functions` 导入；区域包接入后，同名函数也可从 `dnakit.regulation` 导入。

覆盖的 GimmeMotifs 例程有 `MotifComparer.compare_motifs`、`MotifComparer.pvalue`、`MotifComparer.generate_score_dist`、`comparison.akl` 与 `ssd`（已由 `compare_motif_positions` 提供）、`Motif.pcc`、`max_pcc`、`ic`、`matrix_ic`、`pwm_min_score`、`pwm_max_score`、`rc`、`shuffle`、`pfm_to_ppm`、`to_motevo`、`to_meme`、`to_ppm`、`to_consensusv2`、`average_motifs`、`sample`、`pwm_scan_score`、`pwm_scan_to_gff`、`consensus_scan` 以及 `utils.ppm_pseudocount`。

## 示例数据

扫描类示例读取同一个三条记录的 FASTA 文件；其余示例直接传入基序矩阵 `A`（4 个位置，共有序列 ACGT）和 `B`（5 个位置）。

**示例文件：** `seqs.fa`

```text
>s1
TTACGTTTACGTAA
>s2
GGGGACGTCCCC
>s3
CCCCCCCCCC
```

**示例文件：** `set.pfm`

```text
>a
10 0 0 0
0 10 0 0
0 0 10 0
0 0 0 10
>b
0 10 0 0
0 0 10 0
0 0 0 10
10 0 0 0
>c
5 5 0 0
0 5 5 0
0 0 5 5
5 0 0 5
>d
10 0 0 0
10 0 0 0
0 0 10 0
0 0 10 0
0 10 0 0
```

## 1) `dnakit.regulation.gimmemotifs_compare_motif_pair()`

- **作用：** 在所有对齐位置和两条链上比较两个基序，返回最佳相似度分数或其 p 值。
- **处理方法：** 构建两个 `Motif` 对象并调用 `MotifComparer.compare_motifs`。`metric` 可取 `wic`（加权信息量，默认）、`pcc`、`ed`、`distance`、`chisq`、`akl`、`ssd` 或 `seqcor`。`match` 为 `total`（整个匹配，含未对齐位置）、`subtotal` 或 `partial`（只用重叠位置）；`combine` 为 `mean` 或 `sum`。`pvalue=True` 时由 `MotifComparer.pvalue` 把分数换算成按基序长度对拟合的正态分布的上尾概率。GimmeMotifs 只自带 `wic`、`total`、`mean` 的分布；其他设置请传入 `score_distribution`，即 `gimmemotifs_motif_score_distribution` 写出的文件。
- **API：** `dnakit.regulation.gimmemotifs_compare_motif_pair(pfm_a, pfm_b, *, metric="wic", match="total", combine="mean", pvalue=False, score_distribution=None)`。
- **输入：** 两个 PFM，以及上述选项；`score_distribution` 只与 `pvalue=True` 一起使用。
- **输出：** `value`（分数或 p 值，见 `value_kind`）、`position`（第二个基序相对第一个的偏移）、`strand`（1 或 -1）、两个基序的长度和 `provenance`。
- **示例代码：**

```python
from dnakit.regulation.gimmemotifs_p28_motif_functions import gimmemotifs_compare_motif_pair

A = [[10, 0, 0, 0], [0, 10, 0, 0], [0, 0, 8, 2], [1, 1, 1, 7]]
B = [[0, 10, 0, 0], [0, 0, 10, 0], [1, 1, 1, 7], [5, 5, 0, 0], [0, 0, 0, 10]]

for kwargs in ({}, {"metric": "pcc"}, {"metric": "akl", "match": "partial"}, {"pvalue": True}):
    result = gimmemotifs_compare_motif_pair(A, B, **kwargs)
    print(kwargs, result["value_kind"], round(result["value"], 4), result["position"], result["strand"])
```

- **示例输出：**

```text
{} score -1.6806 1 1
{'metric': 'pcc'} score 0.4948 1 1
{'metric': 'akl', 'match': 'partial'} score 9.4234 -2 -1
{'pvalue': True} pvalue 0.1521 1 1
```

## 2) `dnakit.regulation.gimmemotifs_motif_score_distribution()`

- **作用：** 按长度计算所有基序对的分数分布，即基序相似性 p 值所依据的表。
- **处理方法：** 读取 PFM 文件中的基序，串行运行 `MotifComparer.generate_score_dist`：每个基序与每个基序比较，对每一对基序长度写出分数的均值和标准差（`长度1<TAB>长度2<TAB>均值<TAB>sd`）。上游把文件写进 GimmeMotifs 安装目录；这里写入 `output_file`（必须不存在），不改动已安装的表。`metric` 取 `wic`、`pcc`、`ed` 或 `distance`，`match` 取 `total` 或 `subtotal`。
- **API：** `dnakit.regulation.gimmemotifs_motif_score_distribution(motifs_pfm, output_file, *, metric="wic", match="total", combine="mean")`。
- **输入：** 含 1 到 1000 个基序的 PFM 文件，以及一个新的输出路径。
- **输出：** `distribution` 各行（`length_a`、`length_b`、`mean`、`sd`）、`n_motifs`，以及 `outputs` 中写出的文件。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.regulation.gimmemotifs_p28_motif_functions import gimmemotifs_compare_motif_pair, gimmemotifs_motif_score_distribution

A = [[10, 0, 0, 0], [0, 10, 0, 0], [0, 0, 8, 2], [1, 1, 1, 7]]
B = [[0, 10, 0, 0], [0, 0, 10, 0], [1, 1, 1, 7], [5, 5, 0, 0], [0, 0, 0, 10]]

out = Path(tempfile.mkdtemp()) / "pcc_dist.txt"
dist = gimmemotifs_motif_score_distribution("set.pfm", out, metric="pcc")
print(dist["n_motifs"], [(r["length_a"], r["length_b"]) for r in dist["distribution"]])
scored = gimmemotifs_compare_motif_pair(A, B, metric="pcc", pvalue=True, score_distribution=out)
print(scored["value_kind"], 0 <= scored["value"] <= 1)
```

- **示例输出：**

```text
4 [(4, 4), (4, 5), (5, 4), (5, 5)]
pvalue True
```

## 3) `dnakit.regulation.gimmemotifs_motif_pcc()`

- **作用：** 两个基序各位置概率的 Pearson 相关，可取某个固定对齐，也可取所有对齐中最好的一个。
- **处理方法：** 给出 `position`（第二个基序相对第一个的偏移，-2 表示向左两格）时调用 `Motif.pcc`：对齐两个概率矩阵、截去未对齐的两端，再把各行四个值的 Pearson 相关相加。不给 `position` 时调用 `Motif.max_pcc`：尝试所有偏移（`revcomp` 为真时还包括第二个基序的反向互补），取最高的和。
- **API：** `dnakit.regulation.gimmemotifs_motif_pcc(pfm_a, pfm_b, *, position=None, revcomp=True)`。
- **输入：** 两个 PFM、可选的整数偏移和 `revcomp`。
- **输出：** `mode`（`pcc` 或 `max_pcc`）、`score`、`position` 与 `strand`（1 或 -1）。
- **示例代码：**

```python
from dnakit.regulation.gimmemotifs_p28_motif_functions import gimmemotifs_motif_pcc

A = [[10, 0, 0, 0], [0, 10, 0, 0], [0, 0, 8, 2], [1, 1, 1, 7]]
B = [[0, 10, 0, 0], [0, 0, 10, 0], [1, 1, 1, 7], [5, 5, 0, 0], [0, 0, 0, 10]]

best = gimmemotifs_motif_pcc(A, B)
print(best["mode"], round(best["score"], 4), best["position"], best["strand"])
fixed = gimmemotifs_motif_pcc(A, B, position=-1)
print(fixed["mode"], round(fixed["score"], 4), fixed["position"])
```

- **示例输出：**

```text
max_pcc 2.9685 1 1
pcc -0.4291 -1
```

## 4) `dnakit.regulation.gimmemotifs_motif_ic_similarity()`

- **作用：** 两个基序的信息量相似度，同时给出每个基序自身的信息量。
- **处理方法：** `mode="matrix"` 调用 `Motif.matrix_ic`：对两条链和所有偏移，用两个概率行相对均匀背景共有的信息量打分，返回最佳分数、偏移和链。`mode="aligned"` 调用 `Motif.ic`，只评分 `position` 处的一个对齐；没有重叠的位置与背景比较，并乘以 `bg_factor`。同时报告两个基序的 `Motif.information_content`（单位 bit，完全保守的位置为 2）。当第二个基序比第一个长时，GimmeMotifs 的 `matrix_ic` 对某些组合会抛出广播错误（上游缺陷）；函数会报告该错误并提示交换两个基序或改用 `aligned`。
- **API：** `dnakit.regulation.gimmemotifs_motif_ic_similarity(pfm_a, pfm_b, *, mode="matrix", position=0, bg_factor=1.0)`。
- **输入：** 两个 PFM、`mode`、整数 `position` 和非负的 `bg_factor`。
- **输出：** `score`、`position`、`strand`、`information_content_a` 与 `information_content_b`。
- **示例代码：**

```python
from dnakit.regulation.gimmemotifs_p28_motif_functions import gimmemotifs_motif_ic_similarity

A = [[10, 0, 0, 0], [0, 10, 0, 0], [0, 0, 8, 2], [1, 1, 1, 7]]
B = [[0, 10, 0, 0], [0, 0, 10, 0], [1, 1, 1, 7], [5, 5, 0, 0], [0, 0, 0, 10]]

best = gimmemotifs_motif_ic_similarity(A, B)
print(best["mode"], round(best["score"], 4), best["position"], best["strand"])
print(round(best["information_content_a"], 4), round(best["information_content_b"], 4))
one = gimmemotifs_motif_ic_similarity(A, B, mode="aligned", position=1)
print(one["mode"], round(one["score"], 4))
```

- **示例输出：**

```text
matrix 3.5532 1 1
5.909 7.6267
aligned -14.7849
```

## 5) `dnakit.regulation.gimmemotifs_motif_score_range()`

- **作用：** 基序可能达到的最小和最大 PWM（对数几率）分数。
- **处理方法：** 构建 `Motif` 并读取 `min_score` 与 `max_score`，即 `Motif.pwm_min_score` 与 `pwm_max_score` 背后的值（上游已弃用这两个方法，改用属性）。单元格的对数几率为 `log(p / 0.25 + 0.01)`；分数是各位置最小值之和与最大值之和。
- **API：** `dnakit.regulation.gimmemotifs_motif_score_range(pfm)`。
- **输入：** 一个 PFM。
- **输出：** `min_score`、`max_score` 与 `length`。
- **示例代码：**

```python
from dnakit.regulation.gimmemotifs_p28_motif_functions import gimmemotifs_motif_score_range

A = [[10, 0, 0, 0], [0, 10, 0, 0], [0, 0, 8, 2], [1, 1, 1, 7]]
B = [[0, 10, 0, 0], [0, 0, 10, 0], [1, 1, 1, 7], [5, 5, 0, 0], [0, 0, 0, 10]]

result = gimmemotifs_motif_score_range(A)
print(result["length"], round(result["min_score"], 4), round(result["max_score"], 4))
```

- **示例输出：**

```text
4 -14.5894 4.9759
```

## 6) `dnakit.regulation.gimmemotifs_motif_consensus()`

- **作用：** 基序的 IUPAC 共有序列。
- **处理方法：** 调用 `Motif.to_consensusv2`，返回 GimmeMotifs 的共有序列：某碱基占优势的位置写大写碱基，两到三个碱基合起来占优势时写大写 IUPAC 码，较弱的简并写小写码，信息量不足的位置写 `n`。
- **API：** `dnakit.regulation.gimmemotifs_motif_consensus(pfm)`。
- **输入：** 一个 PFM。
- **输出：** `consensus` 与 `length`。
- **示例代码：**

```python
from dnakit.regulation.gimmemotifs_p28_motif_functions import gimmemotifs_motif_consensus

A = [[10, 0, 0, 0], [0, 10, 0, 0], [0, 0, 8, 2], [1, 1, 1, 7]]
B = [[0, 10, 0, 0], [0, 0, 10, 0], [1, 1, 1, 7], [5, 5, 0, 0], [0, 0, 0, 10]]

print(gimmemotifs_motif_consensus(A)["consensus"])
print(gimmemotifs_motif_consensus(B)["consensus"])
print(gimmemotifs_motif_consensus([[3, 1, 1, 1], [2, 2, 2, 2]])["consensus"])
```

- **示例输出：**

```text
ACGT
CGTmT
nn
```

## 7) `dnakit.regulation.gimmemotifs_motif_reverse_complement()`

- **作用：** 对基序矩阵取反向互补。
- **处理方法：** 调用 `Motif.rc`：反转位置顺序并反转 A、C、G、T 列，使 A 与 T、C 与 G 配对，再由计数重建基序。
- **API：** `dnakit.regulation.gimmemotifs_motif_reverse_complement(pfm)`。
- **输入：** 一个 PFM。
- **输出：** 反向互补基序的 `pfm`（计数）、`ppm`（概率）和 `consensus`。
- **示例代码：**

```python
from dnakit.regulation.gimmemotifs_p28_motif_functions import gimmemotifs_motif_reverse_complement

A = [[10, 0, 0, 0], [0, 10, 0, 0], [0, 0, 8, 2], [1, 1, 1, 7]]
B = [[0, 10, 0, 0], [0, 0, 10, 0], [1, 1, 1, 7], [5, 5, 0, 0], [0, 0, 0, 10]]

result = gimmemotifs_motif_reverse_complement(A)
print(result["pfm"])
print(result["consensus"])
```

- **示例输出：**

```text
[[7.0, 1.0, 1.0, 1.0], [2.0, 8.0, 0.0, 0.0], [0.0, 0.0, 10.0, 0.0], [0.0, 0.0, 0.0, 10.0]]
ACGT
```

## 8) `dnakit.regulation.gimmemotifs_motif_shuffle()`

- **作用：** 打乱基序各位置的顺序，得到组成相同的常用对照基序。
- **处理方法：** 调用 `Motif.shuffle`，用 numpy 对矩阵的行做置换，A、C、G、T 列顺序不变。上游使用 numpy 的全局随机状态；`seed` 在调用前为其设定种子（调用后恢复原状态），使置换可重复。
- **API：** `dnakit.regulation.gimmemotifs_motif_shuffle(pfm, *, seed=None)`。
- **输入：** 一个 PFM 和可选的整数 `seed`。
- **输出：** `pfm`、`ppm` 与 `seed`。
- **示例代码：**

```python
from dnakit.regulation.gimmemotifs_p28_motif_functions import gimmemotifs_motif_shuffle

A = [[10, 0, 0, 0], [0, 10, 0, 0], [0, 0, 8, 2], [1, 1, 1, 7]]
B = [[0, 10, 0, 0], [0, 0, 10, 0], [1, 1, 1, 7], [5, 5, 0, 0], [0, 0, 0, 10]]

first = gimmemotifs_motif_shuffle(A, seed=1)
again = gimmemotifs_motif_shuffle(A, seed=1)
print(first["pfm"])
print(first["pfm"] == again["pfm"])
```

- **示例输出：**

```text
[[1.0, 1.0, 1.0, 7.0], [0.0, 0.0, 8.0, 2.0], [10.0, 0.0, 0.0, 0.0], [0.0, 10.0, 0.0, 0.0]]
True
```

## 9) `dnakit.regulation.gimmemotifs_pfm_to_ppm()`

- **作用：** 用伪计数把基序计数转换成频率。
- **处理方法：** 调用 `Motif.pfm_to_ppm`：每个单元格变为 `(计数 + pseudo) / (行和 + 4 * pseudo)`，每行之和为 1，且没有零值。
- **API：** `dnakit.regulation.gimmemotifs_pfm_to_ppm(pfm, *, pseudo=0.001)`。
- **输入：** 一个计数 PFM 和非负的 `pseudo`。
- **输出：** `ppm`，每个位置一行四个概率。
- **示例代码：**

```python
from dnakit.regulation.gimmemotifs_p28_motif_functions import gimmemotifs_pfm_to_ppm

A = [[10, 0, 0, 0], [0, 10, 0, 0], [0, 0, 8, 2], [1, 1, 1, 7]]
B = [[0, 10, 0, 0], [0, 0, 10, 0], [1, 1, 1, 7], [5, 5, 0, 0], [0, 0, 0, 10]]

ppm = gimmemotifs_pfm_to_ppm(A)["ppm"]
print([round(v, 5) for v in ppm[0]])
print([round(v, 5) for v in ppm[2]])
print(round(sum(ppm[3]), 6))
```

- **示例输出：**

```text
[0.9997, 0.0001, 0.0001, 0.0001]
[0.0001, 0.0001, 0.79978, 0.20002]
1.0
```

## 10) `dnakit.regulation.gimmemotifs_ppm_pseudocount()`

- **作用：** 给概率矩阵加伪计数并对每行重新归一化。
- **处理方法：** 调用 `gimmemotifs.utils.ppm_pseudocount`：每行计算 `(p + pseudo) / sum(p + pseudo)`。默认 `1e-6` 是 GimmeMotifs 比较基序前加入的值，使概率不会恰好为零。
- **API：** `dnakit.regulation.gimmemotifs_ppm_pseudocount(ppm, *, pseudo=1e-06)`。
- **输入：** 概率矩阵（每行四个值）和非负的 `pseudo`。
- **输出：** 加了伪计数后的 `ppm`。
- **示例代码：**

```python
from dnakit.regulation.gimmemotifs_p28_motif_functions import gimmemotifs_ppm_pseudocount

result = gimmemotifs_ppm_pseudocount([[1, 0, 0, 0], [0.25, 0.25, 0.25, 0.25]])
print([f"{v:.9f}" for v in result["ppm"][0]])
print(result["ppm"][1])
```

- **示例输出：**

```text
['0.999997000', '0.000001000', '0.000001000', '0.000001000']
[0.25, 0.25, 0.25, 0.25]
```

## 11) `dnakit.regulation.gimmemotifs_format_motif()`

- **作用：** 把基序写成 GimmeMotifs 导出格式的文本。
- **处理方法：** `output_format` 选择 GimmeMotifs 的方法：`ppm`（`Motif.to_ppm`，`>名称` 行加用制表符分隔、保留 `precision` 位小数的概率）、`pfm`（`Motif.to_pfm`，计数）、`meme`（`Motif.to_meme`，MEME 的 `letter-probability matrix` 块）、`motevo`（`Motif.to_motevo`，MotEvo 计数格式）或 `transfac`（`Motif.to_transfac`）。文本直接返回，不写文件。
- **API：** `dnakit.regulation.gimmemotifs_format_motif(pfm, *, output_format="ppm", motif_id="motif", precision=4)`。
- **输入：** 一个 PFM、格式、单行的 `motif_id` 和 0 到 12 的 `precision`。
- **输出：** `text`、`output_format` 与 `motif_id`。
- **示例代码：**

```python
from dnakit.regulation.gimmemotifs_p28_motif_functions import gimmemotifs_format_motif

A = [[10, 0, 0, 0], [0, 10, 0, 0], [0, 0, 8, 2], [1, 1, 1, 7]]
B = [[0, 10, 0, 0], [0, 0, 10, 0], [1, 1, 1, 7], [5, 5, 0, 0], [0, 0, 0, 10]]

for fmt in ("ppm", "motevo", "meme"):
    print(gimmemotifs_format_motif(A, output_format=fmt, motif_id="mA", precision=2)["text"])
```

- **示例输出：**

```text
>mA
1.00	0.00	0.00	0.00
0.00	1.00	0.00	0.00
0.00	0.00	0.80	0.20
0.10	0.10	0.10	0.70
//
NA mA
P0	A	C	G	T
0	10	0	0	0
1	0	10	0	0
2	0	0	8	2
3	1	1	1	7
//
MOTIF mA
BL   MOTIF mA width=0 seqs=0
letter-probability matrix: alength= 4 w= 4 nsites= 10 E= 0
0.9997	0.0001	0.0001	0.0001
0.0001	0.9997	0.0001	0.0001
0.0001	0.0001	0.7998	0.2
0.1	0.1001	0.1001	0.6998
```

## 12) `dnakit.regulation.gimmemotifs_average_motifs()`

- **作用：** 按给定对齐把两个基序平均成一个新基序。
- **处理方法：** 调用 `Motif.average_motifs`。`position` 是第二个基序相对第一个的偏移（-2 表示向左两格），`orientation` 为 1 表示按原样使用第二个基序，-1 表示使用其反向互补。只被一个基序覆盖的位置用另一个基序的频率补全（按第一个基序的计数缩放），`include_bg` 为真时则用均匀背景补全；然后把对齐后的矩阵计数相加。
- **API：** `dnakit.regulation.gimmemotifs_average_motifs(pfm_a, pfm_b, *, position=0, orientation=1, include_bg=False)`。
- **输入：** 两个 PFM、整数 `position`、取 1 或 -1 的 `orientation` 和 `include_bg`。
- **输出：** 平均基序的 `pfm`、`ppm` 与 `consensus`。
- **示例代码：**

```python
from dnakit.regulation.gimmemotifs_p28_motif_functions import gimmemotifs_average_motifs

A = [[10, 0, 0, 0], [0, 10, 0, 0], [0, 0, 8, 2], [1, 1, 1, 7]]
B = [[0, 10, 0, 0], [0, 0, 10, 0], [1, 1, 1, 7], [5, 5, 0, 0], [0, 0, 0, 10]]

result = gimmemotifs_average_motifs(A, B, position=-1, orientation=-1)
print(result["pfm"])
print(result["consensus"])
```

- **示例输出：**

```text
[[20.0, 0.0, 0.0, 0.0], [10.0, 0.0, 5.0, 5.0], [7.0, 11.0, 1.0, 1.0], [0.0, 10.0, 8.0, 2.0], [1.0, 1.0, 11.0, 7.0]]
Awmsk
```

## 13) `dnakit.regulation.gimmemotifs_sample_motif_sequences()`

- **作用：** 按基序的概率矩阵随机抽样序列。
- **处理方法：** 用 numpy `default_rng(seed)` 调用 `Motif.sample`：每条序列的每个位置独立地按该位置的概率抽取碱基。某碱基概率为 1 的位置总是得到该碱基。
- **API：** `dnakit.regulation.gimmemotifs_sample_motif_sequences(pfm, n_sequences, *, seed=None)`。
- **输入：** 一个 PFM、序列数（1 到 1,000,000）和可选的整数 `seed`。
- **输出：** `sequences`、`n_sequences` 与 `seed`。
- **示例代码：**

```python
from dnakit.regulation.gimmemotifs_p28_motif_functions import gimmemotifs_sample_motif_sequences

A = [[10, 0, 0, 0], [0, 10, 0, 0], [0, 0, 8, 2], [1, 1, 1, 7]]
B = [[0, 10, 0, 0], [0, 0, 10, 0], [1, 1, 1, 7], [5, 5, 0, 0], [0, 0, 0, 10]]

result = gimmemotifs_sample_motif_sequences(A, 6, seed=4)
print(result["sequences"])
print(all(s.startswith("AC") for s in result["sequences"]))
```

- **示例输出：**

```text
['ACGC', 'ACGT', 'ACTT', 'ACGT', 'ACTT', 'ACTC']
True
```

## 14) `dnakit.regulation.gimmemotifs_pwm_scan()`

- **作用：** 用一个基序扫描 FASTA 序列，输出分数、位置、完整命中或 GFF 文件。
- **处理方法：** 截断值是从基序最小 PWM 分数到最大 PWM 分数之间的比例（0 到 1）。`output_format="scores"` 调用 `Motif.pwm_scan_score`（每个报告匹配的分数；默认 `cutoff=0`、`nreport=1`，即每条序列的最佳匹配），`"positions"` 调用 `Motif.scan`（匹配位置），`"hits"` 调用 `Motif.scan_all`（每个匹配的位置、分数和链），`"gff"` 调用 `Motif.pwm_scan_to_gff`，并把命中写入新的 GFF 文件（`gff_path`），基序名为 `motif_id`。后三者默认 `cutoff=0.9`、`nreport=50`。`scan_rc` 表示包含反向链。序列名是完整的 FASTA 标题行。
- **API：** `dnakit.regulation.gimmemotifs_pwm_scan(sequences_fasta, pfm, *, output_format="scores", cutoff=None, nreport=None, scan_rc=True, gff_path=None, motif_id="motif")`。
- **输入：** 至多 100,000 条序列的 FASTA 文件、一个 PFM 和上述选项。
- **输出：** 每条序列的 `matches`、所用参数；`gff` 时在 `outputs` 中给出文件。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.regulation.gimmemotifs_p28_motif_functions import gimmemotifs_pwm_scan

A = [[10, 0, 0, 0], [0, 10, 0, 0], [0, 0, 8, 2], [1, 1, 1, 7]]
B = [[0, 10, 0, 0], [0, 0, 10, 0], [1, 1, 1, 7], [5, 5, 0, 0], [0, 0, 0, 10]]

best = gimmemotifs_pwm_scan("seqs.fa", A)
print({name: [round(v, 3) for v in values] for name, values in best["matches"].items()})
hits = gimmemotifs_pwm_scan("seqs.fa", A, output_format="hits", cutoff=0.95)
print(hits["matches"]["s1"])
gff = Path(tempfile.mkdtemp()) / "hits.gff"
gimmemotifs_pwm_scan("seqs.fa", A, output_format="gff", cutoff=0.95, gff_path=gff, motif_id="mA")
print(gff.read_text().splitlines()[0])
```

- **示例输出：**

```text
{'s1': [4.976], 's2': [4.976], 's3': [-8.634]}
[{'position': 2, 'score': 4.975905319526848, 'strand': 1}, {'position': 8, 'score': 4.975905319526848, 'strand': 1}, {'position': 2, 'score': 4.975905319526848, 'strand': -1}, {'position': 8, 'score': 4.975905319526848, 'strand': -1}]
s1	pfmscan	misc_feature	2	6	4.976	+	.	motif_name "mA" ; motif_instance "ACGT"
```

## 15) `dnakit.regulation.gimmemotifs_consensus_scan()`

- **作用：** 把基序的 IUPAC 共有序列作为精确模式在 FASTA 序列中搜索。
- **处理方法：** 调用 `Motif.consensus_scan`：把共有序列（`Motif.to_consensusv2`）的每个符号展开为字符类，再在每条序列的正向链上搜索得到的正则表达式。报告每个匹配的中间坐标（0 起始，因此可能以 .5 结尾）。
- **API：** `dnakit.regulation.gimmemotifs_consensus_scan(sequences_fasta, pfm)`。
- **输入：** 至多 100,000 条序列的 FASTA 文件和一个 PFM。
- **输出：** `matches`（每条序列的匹配中点）和所搜索的 `consensus`。
- **示例代码：**

```python
from dnakit.regulation.gimmemotifs_p28_motif_functions import gimmemotifs_consensus_scan

A = [[10, 0, 0, 0], [0, 10, 0, 0], [0, 0, 8, 2], [1, 1, 1, 7]]
B = [[0, 10, 0, 0], [0, 0, 10, 0], [1, 1, 1, 7], [5, 5, 0, 0], [0, 0, 0, 10]]

result = gimmemotifs_consensus_scan("seqs.fa", A)
print(result["consensus"], result["matches"])
```

- **示例输出：**

```text
ACGT {'s1': [4.0, 10.0], 's2': [6.0], 's3': []}
```
