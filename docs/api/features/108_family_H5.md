# 108 复核 family-variant（H5）：微生物组、调控、表观基因组、变异、预测

本页记录 H5 组对台账中 `family-variant` 行的复核：共 122 行，涉及 17 个包。每一行都被检查是否真是已有封装的别名、重导出，或可通过已有封装的某个参数到达的分支。核实后：40 行改为 `duplicate-of`（同一对象或同一代码），51 行保留 `family-variant`（已确认参数存在且结果与上游直接调用一致），31 行改为 `wrapped`，由 26 个新函数封装。新函数全部在本地计算，不联网；输入为草图、数字列表或本地文件。可选依赖：`bio-sourmash`、`bio-gimmemotifs`、`bio-tangermeme`、`bio-macs3`、`bio-cooltools`、`bio-pairtools`、`bio-deeptools`、`bio-methylpy`、`bio-varcode`、`bio-reads`（pysam）。模块为 `dnakit.microbiome`、`dnakit.regulation`、`dnakit.regulation`、`dnakit.epigenome`、`dnakit.epigenome`、`dnakit.epigenome`、`dnakit.epigenome`、`dnakit.epigenome`、`dnakit.variants` 与 `dnakit.variants`，接入包导出由总负责人完成，示例因此直接从模块导入。

## 复核结果

| 包 | 复核前 | duplicate-of | family-variant（已核实） | wrapped |
| --- | --- | --- | --- | --- |
| sourmash | 36 | 19 | 5 | 12 |
| gimmemotifs | 16 | 0 | 13 | 3 |
| MACS3 | 14 | 7 | 2 | 5 |
| grelu | 10 | 2 | 8 | 0 |
| cooltools | 8 | 1 | 5 | 2 |
| alphagenome | 7 | 0 | 7 | 0 |
| pairtools | 5 | 2 | 1 | 2 |
| tangermeme | 5 | 2 | 1 | 2 |
| deeptools | 5 | 4 | 0 | 1 |
| varcode | 4 | 1 | 1 | 2 |
| bed_reader | 3 | 1 | 2 | 0 |
| methylpy | 2 | 0 | 1 | 1 |
| hgvs | 2 | 0 | 2 | 0 |
| vcfpy | 2 | 0 | 2 | 0 |
| pysam | 1 | 0 | 0 | 1 |
| dendropy | 1 | 1 | 0 | 0 |
| enformer_pytorch | 1 | 0 | 1 | 0 |
| **合计** | 122 | 40 | 51 | 31 |

## 1) `dnakit.microbiome.write_sourmash_signature_archive()`

- **作用：** 把 MinHash 草图写成 sourmash 的 zip 压缩集合或 SQLite 索引。
- **处理方法：** 按输出文件名选择格式（.zip 与 .sig.zip 为 zip，.sqldb 为 SQLite），用 sourmash.save_load.SaveSignaturesToLocation 写入；草图带 abundances 时保留丰度。已有的 write_sourmash_signatures 只写 .sig 和 .sig.gz。
- **API：** `dnakit.microbiome.write_sourmash_signature_archive(sketches, output, *, verbose=True)`。
- **输入：** sketches 为草图字典列表（ksize、scaled、hashes，可选 id 与 abundances）；output 为尚不存在的路径。
- **输出：** format（zip 或 sqlite）、signatures（数量）、output 与 provenance；文件可用 read_sourmash_signatures 读回。
- **示例代码：**

```python
a = {"id": "a", "ksize": 21, "scaled": 1, "hashes": [1, 2, 3, 4, 5, 6]}
b = {"id": "b", "ksize": 21, "scaled": 1, "hashes": [4, 5, 6, 7, 8, 9]}
c = {"id": "c", "ksize": 21, "scaled": 1, "hashes": [20, 21, 22]}
import tempfile
from pathlib import Path

from dnakit.microbiome import read_sourmash_signatures
from dnakit.microbiome import write_sourmash_signature_archive

with tempfile.TemporaryDirectory() as folder:
    for name in ("sketches.sig.zip", "sketches.sqldb"):
        result = write_sourmash_signature_archive([a, b], Path(folder) / name, verbose=False)
        back = read_sourmash_signatures(result["output"], verbose=False)
        print(result["format"], result["signatures"], [s["id"] for s in back["sketches"]])
        print(back["sketches"][1]["hashes"])
```

- **示例输出：**

```text
zip 2 ['a', 'b']
[4, 5, 6, 7, 8, 9]
sqlite 2 ['a', 'b']
[4, 5, 6, 7, 8, 9]
```

## 2) `dnakit.microbiome.parse_sourmash_signature_json()`

- **作用：** 把 sourmash 签名的 JSON 文本解析为草图，不读文件。
- **处理方法：** require_single=False 时调用 sourmash.signature.load_signatures_from_json；True 时调用 load_one_signature_from_json，要求恰好一个签名匹配。可用 ksize 选择 k-mer 大小。输出按 id、ksize、md5sum 排序。
- **API：** `dnakit.microbiome.parse_sourmash_signature_json(json_text, *, ksize=None, require_single=False, max_signatures=1000, verbose=True)`。
- **输入：** json_text 为 JSON 字符串；ksize 可选；require_single 为布尔值；max_signatures 为上限。
- **输出：** signatures（数量）、sketches（id、ksize、scaled、hash_count、hashes、md5sum，有丰度时含 abundances）与 provenance。
- **示例代码：**

```python
a = {"id": "a", "ksize": 21, "scaled": 1, "hashes": [1, 2, 3, 4, 5, 6]}
b = {"id": "b", "ksize": 21, "scaled": 1, "hashes": [4, 5, 6, 7, 8, 9]}
c = {"id": "c", "ksize": 21, "scaled": 1, "hashes": [20, 21, 22]}
from dnakit.microbiome import write_sourmash_signatures
from dnakit.microbiome import parse_sourmash_signature_json
import tempfile
from pathlib import Path

k31 = {"id": "k31", "ksize": 31, "scaled": 1, "hashes": [1, 2, 3]}
with tempfile.TemporaryDirectory() as folder:
    path = Path(folder) / "all.sig"
    write_sourmash_signatures([a, b, k31], path, verbose=False)
    text = path.read_text()
every = parse_sourmash_signature_json(text, verbose=False)
print(every["signatures"], [s["id"] for s in every["sketches"]])
one = parse_sourmash_signature_json(text, ksize=31, require_single=True, verbose=False)
print(one["sketches"][0]["id"], one["sketches"][0]["ksize"], one["sketches"][0]["hashes"])
```

- **示例输出：**

```text
3 ['a', 'b', 'k31']
k31 31 [1, 2, 3]
```

## 3) `dnakit.microbiome.select_sourmash_signatures()`

- **作用：** 从 sourmash 文件按 k-mer 大小、名称和 md5 选取 DNA 签名。
- **处理方法：** expect='many' 调用 sourmash_args.load_many_signatures；'one' 对每个来源调用 load_one_signature，要求恰好一个；'query' 调用 load_query_signature，要求单一来源且筛选后只剩一个。name_pattern 为不区分大小写的正则，在名称、文件名和 md5 中搜索。
- **API：** `dnakit.microbiome.select_sourmash_signatures(sources, *, ksize=None, name_pattern=None, md5_prefix=None, expect='many', max_signatures=1000, verbose=True)`。
- **输入：** sources 为文件路径列表；ksize、name_pattern、md5_prefix 为可选筛选；expect 为 many、one 或 query。
- **输出：** signatures（数量）、sketches（含 source 与 md5sum，按来源、id、md5sum 排序）与 provenance。
- **示例代码：**

```python
a = {"id": "a", "ksize": 21, "scaled": 1, "hashes": [1, 2, 3, 4, 5, 6]}
b = {"id": "b", "ksize": 21, "scaled": 1, "hashes": [4, 5, 6, 7, 8, 9]}
c = {"id": "c", "ksize": 21, "scaled": 1, "hashes": [20, 21, 22]}
import tempfile
from pathlib import Path

from dnakit.microbiome import write_sourmash_signatures
from dnakit.microbiome import select_sourmash_signatures

k31 = {"id": "k31", "ksize": 31, "scaled": 1, "hashes": [1, 2, 3]}
with tempfile.TemporaryDirectory() as folder:
    path = Path(folder) / "all.sig"
    write_sourmash_signatures([a, b, k31], path, verbose=False)
    def ids(**options):
        result = select_sourmash_signatures([path], verbose=False, **options)
        return [item["id"] for item in result["sketches"]]
    print(ids())
    print(ids(ksize=21))
    print(ids(name_pattern="^B"))
    print(ids(expect="one", ksize=31))
    print(ids(expect="query", ksize=31))
```

- **示例输出：**

```text
['a', 'b', 'k31']
['a', 'b']
['b']
['k31']
['k31']
```

## 4) `dnakit.microbiome.compare_minhash_abundance_matrix()`

- **作用：** 两两比较草图，可使用丰度，也可多进程运行。
- **处理方法：** 调用 sourmash.compare.compare_all_pairs：n_jobs=1 走 compare_serial，否则走 compare_parallel。草图带 abundances 时给出角相似度（余弦），除非 ignore_abundance 为真；否则给出 Jaccard。return_ani=True 返回估计的 ANI。
- **API：** `dnakit.microbiome.compare_minhash_abundance_matrix(sketches, *, ignore_abundance=False, return_ani=False, n_jobs=1, verbose=True)`。
- **输入：** sketches（至少 2 个，ksize 与 scaled 相同；要么都带 abundances 要么都不带）；ignore_abundance、return_ani、n_jobs（1 到 16）。
- **输出：** measure（jaccard、angular_similarity 或 ani）、abundance_weighted、ids、matrix 与 provenance。
- **示例代码：**

```python
from dnakit.microbiome import compare_minhash_abundance_matrix

hashes = [1, 2, 3, 4]
a = {"id": "a", "ksize": 21, "scaled": 1, "hashes": hashes, "abundances": [1, 1, 1, 1]}
b = {"id": "b", "ksize": 21, "scaled": 1, "hashes": hashes, "abundances": [1, 1, 2, 2]}
result = compare_minhash_abundance_matrix([a, b], verbose=False)
print(result["measure"], round(result["matrix"][0][1], 4))
flat = compare_minhash_abundance_matrix([a, b], ignore_abundance=True, verbose=False)
print(flat["measure"], flat["matrix"][0][1])
parallel = compare_minhash_abundance_matrix([a, b], n_jobs=2, verbose=False)
print(parallel["matrix"] == result["matrix"])
```

- **示例输出：**

```text
angular_similarity 0.7952
jaccard 1.0
True
```

## 5) `dnakit.microbiome.search_minhash_database_ranked()`

- **作用：** 按与查询草图的相似度对数据库草图排序，如 sourmash search 命令。
- **处理方法：** mode 为 jaccard、containment、max_containment 时调用 sourmash.search.search_databases_with_flat_query；mode='abundance' 时调用 search_databases_with_abund_query，需要查询和数据库草图都带 abundances，分数为角相似度。低于 threshold 的匹配被丢弃，其余按分数、id 排序。
- **API：** `dnakit.microbiome.search_minhash_database_ranked(query, database, *, mode='jaccard', threshold=0.08, verbose=True)`。
- **输入：** query 与 database（草图字典）；mode；threshold（0 到 1，默认 0.08）。
- **输出：** query_id、mode、match_count、matches（id、md5、score）与 provenance。
- **示例代码：**

```python
a = {"id": "a", "ksize": 21, "scaled": 1, "hashes": list(range(1, 2001))}
b = {"id": "b", "ksize": 21, "scaled": 1, "hashes": list(range(1001, 3001))}
c = {"id": "c", "ksize": 21, "scaled": 1, "hashes": list(range(5001, 5501))}
from dnakit.microbiome import search_minhash_database_ranked

for mode in ("jaccard", "containment", "max_containment"):
    result = search_minhash_database_ranked(a, [a, b, c], mode=mode, verbose=False)
    print(mode, [(m["id"], round(m["score"], 4)) for m in result["matches"]])
```

- **示例输出：**

```text
jaccard [('a', 1.0), ('b', 0.3333)]
containment [('a', 1.0), ('b', 0.5)]
max_containment [('a', 1.0), ('b', 0.5)]
```

## 6) `dnakit.microbiome.search_sourmash_sbt_index()`

- **作用：** 用查询草图搜索已保存的 sourmash 序列布隆树（SBT）索引。
- **处理方法：** 用 sourmash.sbtmh.load_sbt_index 读取索引，再调用 sourmash.sbtmh.search_sbt_index。分数为 Jaccard 指数，因此匹配与 search_minhash_database_ranked(mode='jaccard') 在同一组草图上相同。
- **API：** `dnakit.microbiome.search_sourmash_sbt_index(index, query, *, threshold=0.08, verbose=True)`。
- **输入：** index 为 .sbt.zip 或 .sbt.json 文件（由 sourmash index 生成）；query 为草图，ksize 与 scaled 与索引一致；threshold。
- **输出：** query_id、match_count、matches（id、md5、jaccard，按 jaccard、id 排序）与 provenance。
- **示例代码：**

```python
a = {"id": "a", "ksize": 21, "scaled": 1, "hashes": list(range(1, 2001))}
b = {"id": "b", "ksize": 21, "scaled": 1, "hashes": list(range(1001, 3001))}
c = {"id": "c", "ksize": 21, "scaled": 1, "hashes": list(range(5001, 5501))}
import tempfile
from pathlib import Path

import sourmash
from sourmash.sbtmh import create_sbt_index

from dnakit.microbiome import search_sourmash_sbt_index

tree = create_sbt_index()
for sketch in (a, b, c):
    minhash = sourmash.MinHash(n=0, ksize=21, scaled=1)
    minhash.add_many(sketch["hashes"])
    tree.insert(sourmash.SourmashSignature(minhash, name=sketch["id"]))
with tempfile.TemporaryDirectory() as folder:
    path = Path(folder) / "index.sbt.zip"
    tree.save(str(path))
    result = search_sourmash_sbt_index(path, a, threshold=0.1, verbose=False)
print([(m["id"], round(m["jaccard"], 4)) for m in result["matches"]])
```

- **示例输出：**

```text
[('a', 1.0), ('b', 0.3333)]
```

## 7) `dnakit.microbiome.downsample_minhash_num_sketch()`

- **作用：** 把 num 型 MinHash 草图降到更少的哈希数。
- **处理方法：** 调用 sourmash.minhash.flatten_and_downsample_num，保留最小的 num 个哈希。新的 num 不能超过草图自身的 num。已有的 downsample_minhash_sketch 只处理 scaled 草图。
- **API：** `dnakit.microbiome.downsample_minhash_num_sketch(sketch, *, num, verbose=True)`。
- **输入：** sketch 含 ksize、num、hashes（可选 id）；num 为新的哈希数。
- **输出：** input_num、num、hash_count、hashes 与 provenance。
- **示例代码：**

```python
from dnakit.microbiome import downsample_minhash_num_sketch

sketch = {"id": "n", "ksize": 21, "num": 6, "hashes": [34, 3, 21, 5, 13, 8]}
result = downsample_minhash_num_sketch(sketch, num=3, verbose=False)
print(result["input_num"], result["num"], result["hashes"])
```

- **示例输出：**

```text
6 3 [3, 5, 8]
```

## 8) `dnakit.regulation.normalize_peak_widths()`

- **作用：** 把峰值缩放到围绕中点的固定窗口，按链方向处理。
- **处理方法：** 调用 gimmemotifs.preprocessing.peakfile2bedfile。中点为 (start + end) // 2；+ 链（缺省）为中点 ± window // 2；其他链的起止坐标被上游写反，start 大于 end，此类行带 start_after_end=True。起点小于 0 的行被丢弃而不是截断。缺失的名称、分数、链分别记为空串、0、+。
- **API：** `dnakit.regulation.normalize_peak_widths(path, *, window=200, verbose=True)`。
- **输入：** path 为制表符分隔的 BED 类文件（至少 3 列）；window 为窗口宽度。
- **输出：** intervals（chrom、start、end、name、score、strand、start_after_end）、input_rows、dropped_rows 与 provenance。

**示例文件：** `peaks.bed`

```text
chr1	100	300	p1	5	+
chr1	1000	1100	p2	7	-
chr2	10	50
chr3	500	700	p4
```

- **示例代码：**

```python
from dnakit.regulation import normalize_peak_widths

result = normalize_peak_widths("peaks.bed", window=100, verbose=False)
for item in result["intervals"]:
    print(item["chrom"], item["start"], item["end"], item["name"], item["strand"], item["start_after_end"])
print(result["input_rows"], result["dropped_rows"])
```

- **示例输出：**

```text
chr1 150 250 p1 + False
chr1 1100 1000 p2 - True
chr3 550 650 p4 + False
4 1
```

## 9) `dnakit.regulation.count_motif_matches_at_fpr()`

- **作用：** 统计在给定假阳性率下超过得分阈值的前景和背景序列数。
- **处理方法：** 调用 gimmemotifs.rocmetrics.matches_at_fpr：阈值为背景得分的 (100 - fpr * 100) 百分位，结果为前景得分和背景得分中不低于阈值的个数。motif_score_metric 只返回单个数字，无法返回这两个计数。
- **API：** `dnakit.regulation.count_motif_matches_at_fpr(foreground, background, *, fpr=0.01, verbose=True)`。
- **输入：** foreground 与 background 为数字列表；fpr 在 (0, 1]。
- **输出：** foreground_matches、background_matches、两者占比、各自总数、fpr 与 provenance。
- **示例代码：**

```python
from dnakit.regulation import count_motif_matches_at_fpr

background = [float(value) for value in range(100)]
foreground = [95.0, 99.0, 50.0, 80.0, 98.5]
result = count_motif_matches_at_fpr(foreground, background, fpr=0.1, verbose=False)
print(result["foreground_matches"], result["background_matches"], result["foreground_fraction"])
```

- **示例输出：**

```text
3 10 0.6
```

## 10) `dnakit.regulation.calculate_partial_motif_roc_auc()`

- **作用：** 计算得分 ROC 曲线在假阳性率上限以内的面积。
- **处理方法：** 调用 gimmemotifs.rocmetrics.roc_auc_xlim，传入前景和背景得分列表（上游参数虽名为 x_bla、y_bla，实为得分）。面积是真阳性率对 0 到 fpr_limit 假阳性率的积分，不超过 fpr_limit。motif_score_metric 把 ROC 曲线点传给该函数，会使 gimmemotifs 报错。
- **API：** `dnakit.regulation.calculate_partial_motif_roc_auc(foreground, background, *, fpr_limit=0.1, verbose=True)`。
- **输入：** foreground 与 background 为有限数字列表；fpr_limit 在 (0, 1]。
- **输出：** area、area_fraction_of_limit、各自数量、fpr_limit 与 provenance。
- **示例代码：**

```python
from dnakit.regulation import calculate_partial_motif_roc_auc

foreground = [3.0, 4.0, 5.0]
background = [1.0, 2.0, 3.5, 4.5]
result = calculate_partial_motif_roc_auc(foreground, background, fpr_limit=0.5, verbose=False)
print(result["area"], result["area_fraction_of_limit"])
```

- **示例输出：**

```text
0.25 0.5
```

## 11) `dnakit.regulation.design_beam_substitution()`

- **作用：** 用束搜索把基序替换进序列，使模型输出接近目标。
- **处理方法：** 调用 tangermeme.design.beam_substitution。每轮把每个基序（reverse_complement=True 时含反向互补）铺到束中每条序列的每个位置，按与 target 的均方误差给完整序列打分，保留 beam_size 条互不相同的最佳序列。最佳损失的改进不超过 tol 或达到 max_iter 时停止。beam_size=1 等于 design_greedy_substitution。
- **API：** `dnakit.regulation.design_beam_substitution(model, sequence, target, *, motifs=None, reverse_complement=True, beam_size=4, n_best=1, max_iter=10, tol=0.001, verbose=True)`。
- **输入：** model 为 eval 模式的 PyTorch 模型；sequence 为单条 DNA；target 为每个输出一个数；motifs 可选；beam_size、n_best（不超过 beam_size）、max_iter、tol。
- **输出：** designed_sequences（损失从低到高）、designed_sequence（最佳）、input_sequence、beam_size 与 provenance。
- **示例代码：**

```python
import torch

from dnakit.regulation import design_greedy_substitution
from dnakit.regulation import design_beam_substitution

length = 20
weights = torch.tensor(
    [[((base * 7 + pos * 3) % 11) / 4.0 - 1.0 for pos in range(length)] for base in range(4)]
)


class Linear(torch.nn.Module):
    def __init__(self):
        super().__init__()
        self.register_buffer("weights", weights)

    def forward(self, x, *args):
        return (x * self.weights).sum(dim=(1, 2)).unsqueeze(-1)


model = Linear().eval()
sequence = "ACGTTGATAAGCCGTAGGAT"
beam = design_beam_substitution(
    model, sequence, [5.0], motifs=["GGG", "AC"], beam_size=4, n_best=2, max_iter=3, verbose=False
)
print(beam["designed_sequences"])
single = design_beam_substitution(model, sequence, [5.0], motifs=["GGG"], beam_size=1, max_iter=2, verbose=False)
greedy = design_greedy_substitution(model, sequence, [5.0], motifs=["GGG"], max_iter=2, verbose=False)
print(single["designed_sequence"] == greedy["designed_sequence"])
```

- **示例输出：**

```text
['ACCCCGATAAGCCGTACCCT', 'ACCCCGATACCCCGTACCCT']
True
```

## 12) `dnakit.regulation.insert_motif()`

- **作用：** 把基序插入每条序列，并保留变长后的结果。
- **处理方法：** 调用 tangermeme.ersatz.insert。基序插在第 start 个碱基之前（0 起，缺省为序列中心），每条输出比输入长 len(motif)，不丢失碱基。place_motif(mode='insert') 会从右侧截回原长度，丢掉末端碱基。
- **API：** `dnakit.regulation.insert_motif(sequences, motif, *, start=None, verbose=True)`。
- **输入：** sequences 为 DNA 字符串或名称到字符串的映射（等长）；motif 为 DNA；start 可选。
- **输出：** sequences（name、sequence）、input_length、output_length、start 与 provenance。
- **示例代码：**

```python
from dnakit.regulation import place_motif
from dnakit.regulation import insert_motif

sequence = "ACGTTGATAAGCCGTAGGAT"
longer = insert_motif({"s": sequence}, "GG", start=3, verbose=False)
print(longer["sequences"][0]["sequence"], longer["input_length"], longer["output_length"])
trimmed = place_motif({"s": sequence}, "GG", start=3, mode="insert", verbose=False)
print(trimmed["sequences"][0]["sequence"])
```

- **示例输出：**

```text
ACGGGTTGATAAGCCGTAGGAT 20 22
ACGGGTTGATAAGCCGTAGG
```

## 13) `dnakit.epigenome.build_counted_fragment_pileup()`

- **作用：** 对每个代表多条相同读段的双端片段做计数加权堆叠。
- **处理方法：** 调用 MACS3.Signal.PileupV2.pileup_from_LRC_as_list：第 i 个片段覆盖 [left, right)，每个碱基加 counts[i]。left_sorted=True 声明片段已按左端排序（会检查），让 MACS3 跳过排序。结果布局与 build_paired_fragment_pileup 相同：positions 为段末端，values 为该段的值。
- **API：** `dnakit.epigenome.build_counted_fragment_pileup(left_positions, right_positions, counts, *, scale_factor=1.0, baseline_value=0.0, left_sorted=False, verbose=True)`。
- **输入：** left_positions、right_positions、counts（0 到 65535）为等长整数列表；scale_factor、baseline_value、left_sorted。
- **输出：** fragments、total_count、positions、values 与 provenance。
- **示例代码：**

```python
from dnakit.epigenome import build_counted_fragment_pileup

left = [10, 20, 40, 100]
right = [50, 60, 90, 150]
counts = [2, 1, 3, 1]
result = build_counted_fragment_pileup(left, right, counts, verbose=False)
print(result["positions"])
print(result["values"])
```

- **示例输出：**

```text
[10, 20, 40, 50, 60, 90, 100, 150]
[0.0, 2.0, 3.0, 6.0, 4.0, 3.0, 0.0, 1.0]
```

## 14) `dnakit.epigenome.build_fragment_end_pileup()`

- **作用：** 堆叠以计数加权片段两端为中心的窗口。
- **处理方法：** 调用 MACS3.Signal.PileupV2.pileup_from_LRC_centers_as_list。每个片段贡献两个宽 window 的窗口，权重 counts[i]：一个从 left - window // 2 开始，另一个从 right - window // 2 开始。片段中间不被覆盖，因此图形显示片段起止位置。
- **API：** `dnakit.epigenome.build_fragment_end_pileup(left_positions, right_positions, counts, *, window, scale_factor=1.0, baseline_value=0.0, verbose=True)`。
- **输入：** left_positions、right_positions、counts 为等长整数列表；window 为窗口宽度；scale_factor、baseline_value。
- **输出：** fragments、window、total_count、positions、values 与 provenance。
- **示例代码：**

```python
from dnakit.epigenome import build_fragment_end_pileup

result = build_fragment_end_pileup([100], [200], [2], window=20, verbose=False)
print(result["positions"])
print(result["values"])
```

- **示例输出：**

```text
[90, 110, 190, 210]
[0.0, 2.0, 0.0, 2.0]
```

## 15) `dnakit.epigenome.build_hmmratac_fragment_pileup()`

- **作用：** 按 HMMRATAC 的方式，用片段长度权重堆叠双端片段。
- **处理方法：** 调用 MACS3.Signal.PileupV2.pileup_from_LR_hmmratac。length_weights 把片段长度（right - left，写成十进制字符串）映射到该片段每个碱基加的权重；出现的每个长度都必须有权重，可来自 hmmratac_fragment_weights。
- **API：** `dnakit.epigenome.build_hmmratac_fragment_pileup(left_positions, right_positions, length_weights, *, verbose=True)`。
- **输入：** left_positions、right_positions 为等长整数列表；length_weights 为长度字符串到权重的映射。
- **输出：** fragments、positions、values 与 provenance。
- **示例代码：**

```python
from dnakit.epigenome import build_hmmratac_fragment_pileup

left = [10, 20, 40, 100]
right = [50, 60, 90, 150]
result = build_hmmratac_fragment_pileup(left, right, {"40": 0.5, "50": 2.0}, verbose=False)
print(result["positions"])
print(result["values"])
```

- **示例输出：**

```text
[10, 20, 40, 50, 60, 90, 100, 150]
[0.0, 0.5, 1.0, 3.0, 2.5, 2.0, 0.0, 2.0]
```

## 16) `dnakit.epigenome.calculate_allele_specific_model_likelihood()`

- **作用：** 计算带等位基因比例的杂合模型的对数似然和 BIC。
- **处理方法：** 调用 MACS3.Signal.VariantStat.CalModel_Heter_AS。处理组读段（_t）用上限为 max_allowed_ratio 的自由等位比例拟合；对照组（_c）不设等位偏好。每个参数列出支持第一等位（top1）或第二等位（top2）的读段碱基质量。calculate_allele_model_likelihood 只覆盖没有自由比例的模型。
- **API：** `dnakit.epigenome.calculate_allele_specific_model_likelihood(top1_bq_t, top1_bq_c, top2_bq_t, top2_bq_c, *, max_allowed_ratio=0.99, verbose=True)`。
- **输入：** top1_bq_t、top1_bq_c、top2_bq_t、top2_bq_c 为等长整数列表；max_allowed_ratio 在 [0.5, 1]。
- **输出：** reads、max_allowed_ratio、log_likelihood、bic 与 provenance。
- **示例代码：**

```python
from dnakit.epigenome import calculate_allele_specific_model_likelihood

reads = ([30, 32, 35, 28, 30, 33], [2, 3, 1, 2, 3, 2], [0] * 6, [0] * 6)
result = calculate_allele_specific_model_likelihood(*reads, verbose=False)
print(round(result["log_likelihood"], 4), round(result["bic"], 4))
```

- **示例输出：**

```text
-11.2048 29.8643
```

## 17) `dnakit.epigenome.calculate_allele_specific_genotype_quality()`

- **作用：** 计算等位特异模型相对第二个模型的 Phred 式得分。
- **处理方法：** 调用 MACS3.Signal.VariantStat.calculate_GQ_heterASsig。取 L1 = 1，L2 = min(1, exp(lnL2 - lnL1))（不小于 1e-110），得分为 int(-10 * log10(L2 / (L1 + L2)))，比值不超过 1e-110 时为 255。两个对数似然相等时为 3，模型 1 拟合越好得分越高。三模型的 calculate_genotype_quality 使用不同公式。
- **API：** `dnakit.epigenome.calculate_allele_specific_genotype_quality(log_likelihood_1, log_likelihood_2, *, verbose=True)`。
- **输入：** log_likelihood_1 与 log_likelihood_2 为自然对数似然。
- **输出：** score（整数）与 provenance。
- **示例代码：**

```python
from dnakit.epigenome import calculate_allele_specific_genotype_quality

for first, second in ((-10.0, -15.0), (-10.0, -10.0), (-15.0, -10.0), (-10.0, -300.0)):
    score = calculate_allele_specific_genotype_quality(first, second, verbose=False)["score"]
    print(first, second, score)
```

- **示例输出：**

```text
-10.0 -15.0 21
-10.0 -10.0 3
-15.0 -10.0 3
-10.0 -300.0 255
```

## 18) `dnakit.epigenome.pairwise_region_sums()`

- **作用：** 汇总每对区域的对角线计数，并为每对不同染色体增加一行总计。
- **处理方法：** 调用 cooltools.sandbox.obs_over_exp_cooler.sum_pairwise，使用 cooler 的染色体视图。region1 等于 region2 且 dist 不小于 0 的行与 pairwise_diagonal_sums 相同；每对不同区域另有一行 dist = -1，count.sum 为该块总计数，n_valid 为有效像素数，n_total 为空。
- **API：** `dnakit.epigenome.pairwise_region_sums(cool_path, *, weight='weight', ignore_diags=2, verbose=True)`。
- **输入：** cool_path 为已平衡的 cooler；weight 为权重列名；ignore_diags 为忽略的对角线数。
- **输出：** rows（region1、region2、dist、n_total、n_valid、count.sum）、row_count、block_row_count 与 provenance。
- **示例代码：**

```python
import tempfile
from pathlib import Path

import cooler
import numpy as np
import pandas as pd

from dnakit.epigenome import pairwise_region_sums

rng = np.random.default_rng(0)
bins = cooler.util.binnify(pd.Series({"chr1": 1500, "chr2": 2000}, name="length"), 100)
rows = []
for i in range(len(bins)):
    for j in range(i, len(bins)):
        same = bins.chrom[i] == bins.chrom[j]
        count = rng.poisson(6.0 / (1 + abs(i - j)) if same else 0.5) + (3 if i == j else 0)
        if count > 0:
            rows.append((i, j, int(count)))
pixels = pd.DataFrame(rows, columns=["bin1_id", "bin2_id", "count"])
with tempfile.TemporaryDirectory() as folder:
    path = Path(folder) / "toy.cool"
    cooler.create_cooler(str(path), bins, pixels, dtypes={"count": "int"}, ensure_sorted=True)
    cooler.balance_cooler(cooler.Cooler(str(path)), store=True)
    result = pairwise_region_sums(path, verbose=False)
block = [row for row in result["rows"] if row["dist"] == -1][0]
print(result["row_count"], result["block_row_count"])
print(block["region1"], block["region2"], block["count.sum"], block["n_valid"])
```

- **示例输出：**

```text
36 1
chr1 chr2 144.0 300.0
```

## 19) `dnakit.epigenome.compute_region_contact_scaling()`

- **作用：** 对一对区域，给出对数间隔距离区间内的原始接触数和可能碱基对面积。
- **处理方法：** 调用 cooltools.sandbox.pairs_scaling_functions.compute_scaling。接触 i 连接同一染色体上的 positions1[i] 与 positions2[i]；只统计 positions1 在 region1、positions2 在 region2（均为 [start, end)，region2 缺省同 region1）内的接触，按 positions2 - positions1 落入的距离区间计数，距离从 min_distance 到 max_distance（不含）。区间面积是可能位置对数的积分，frequency = observed / area。与 compute_contact_scaling 不同，输入是位置而非 pairs 文件，也不需要染色体大小。
- **API：** `dnakit.epigenome.compute_region_contact_scaling(positions1, positions2, *, region1, region2=None, min_distance=10, max_distance=10000000, n_bins=50, verbose=True)`。
- **输入：** positions1、positions2 为等长整数列表；region1（及可选 region2）为 [start, end]；min_distance、max_distance、n_bins。
- **输出：** bins（min_dist、max_dist、observed、area、frequency）、contacts_counted 与 provenance。
- **示例代码：**

```python
from dnakit.epigenome import compute_region_contact_scaling

result = compute_region_contact_scaling(
    [100, 150, 200, 1000, 5000, 200],
    [130, 400, 1000, 1500, 9000, 260],
    region1=[0, 10000],
    min_distance=10,
    max_distance=10000,
    n_bins=5,
    verbose=False,
)
for item in result["bins"]:
    print(item["min_dist"], item["max_dist"], item["observed"], round(item["area"]))
print(result["contacts_counted"])
```

- **示例输出：**

```text
10 56 1 458482
56 316 2 2551640
316 1778 2 13089286
1778 10000 1 33800642
6
```

## 20) `dnakit.epigenome.deduplicate_pairs_online()`

- **作用：** 用 pairtools 的 Cython 在线算法标记 PCR 重复的 pair。
- **处理方法：** 调用 pairtools.lib.dedup.streaming_dedup_cython。若一个 pair 的两侧都与更早保留的 pair 相差不超过 max_mismatch bp，则为重复（method='max' 取两侧位置差的较大值，'sum' 取和）。重复不会串联：deduplicate_pairs 会把靠近重复的 pair 并入同一簇，因此对间隔 2 bp 的一串 pair 它保留得更少。mark_duplicates=True 时写入 duplicates_path 的重复 pair 类型改为 DD。
- **API：** `dnakit.epigenome.deduplicate_pairs_online(pairs_path, *, output_path=None, duplicates_path=None, unmapped_path=None, max_mismatch=3, method='max', mark_duplicates=False, verbose=True)`。
- **输入：** pairs_path 为已排序的 .pairs 文件；output_path、duplicates_path、unmapped_path 为尚不存在的可选路径；max_mismatch、method、mark_duplicates。
- **输出：** total、total_mapped、total_nodups、total_dups、duplicate_fraction、unique_read_ids、duplicate_read_ids（各至多 10000）与 provenance。

**示例文件：** `dup.pairs`

```text
## pairs format v1.0.0
#shape: upper triangle
#chromsize: chr1 10000
#chromsize: chr2 10000
#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type
r1	chr1	100	chr2	500	+	-	UU
r5	chr1	100	chr2	502	+	-	UU
r2	chr1	102	chr2	500	+	-	UU
r3	chr1	104	chr2	500	+	-	UU
r4	chr1	106	chr2	500	+	-	UU
r6	chr1	3000	chr1	4000	+	-	UU
r8	chr1	3000	chr1	4000	-	+	UU
r7	chr1	3001	chr1	4002	+	-	UU
```

- **示例代码：**

```python
from dnakit.epigenome import deduplicate_pairs
from dnakit.epigenome import deduplicate_pairs_online

online = deduplicate_pairs_online("dup.pairs", max_mismatch=2, verbose=False)
print(online["total_nodups"], online["total_dups"])
print(online["unique_read_ids"])
print(sorted(online["duplicate_read_ids"]))
chained = deduplicate_pairs("dup.pairs", max_mismatch=2, verbose=False)
print(chained["total_nodups"], chained["total_dups"])
```

- **示例输出：**

```text
4 4
['r1', 'r3', 'r6', 'r8']
['r2', 'r4', 'r5', 'r7']
3 5
```

## 21) `dnakit.epigenome.normalize_scaling_curve_unity()`

- **作用：** 缩放接触频率曲线，使曲线下面积为 1。
- **处理方法：** 调用 pairtools.lib.scaling.unity_norm_scaling。每个频率除以 sum(bin_length * frequency)，求和只含几何中点 sqrt(low * high) 严格落在 norm_range 内的有限区间；norm_range=None 使用所有有限区间。normalize_scaling_curve 则除以窗口内的平均频率。
- **API：** `dnakit.epigenome.normalize_scaling_curve_unity(bin_edges, frequencies, *, norm_range=(10000, 1000000000), verbose=True)`。
- **输入：** bin_edges 为 n + 1 个升序整数边界；frequencies 为 n 个非负频率；norm_range 为 [low, high] 或 None。
- **输出：** bins（min_dist、max_dist、midpoint、frequency）、norm_range 与 provenance。
- **示例代码：**

```python
from dnakit.epigenome import normalize_scaling_curve_unity

edges = [10, 100, 1000, 10000]
frequencies = [0.5, 0.1, 0.01]
result = normalize_scaling_curve_unity(edges, frequencies, norm_range=None, verbose=False)
print([round(item["frequency"], 6) for item in result["bins"]])
area = sum((high - low) * f for low, high, f in zip(edges[:-1], edges[1:], frequencies))
print(area)
```

- **示例输出：**

```text
[0.002222, 0.000444, 4.4e-05]
225.0
```

## 22) `dnakit.epigenome.chop_regions_from_middle()`

- **作用：** 在一组外显子的中心两侧各切出固定数量的碱基。
- **处理方法：** 调用 deeptools.heatmapper.chopRegionsFromMiddle。中心为外显子总长 // 2 个外显子碱基之后的位置，跨过中心的外显子被拆开。返回中心前最后 left 个碱基和中心后前 right 个碱基，均为 [start, end] 区间列表；某侧不足时 pad_left 或 pad_right 给出缺少的碱基数。与 chop_regions_for_profile 不同，区间从中心切出，也没有主体部分。
- **API：** `dnakit.epigenome.chop_regions_from_middle(exons, *, left=0, right=0, verbose=True)`。
- **输入：** exons 为有序、不重叠的 [start, end] 列表；left、right 为碱基数。
- **输出：** left_bins、right_bins、pad_left、pad_right 与 provenance。
- **示例代码：**

```python
from dnakit.epigenome import chop_regions_from_middle

exons = [[0, 200], [300, 400], [800, 900]]
result = chop_regions_from_middle(exons, left=150, right=80, verbose=False)
print(result["left_bins"], result["right_bins"], result["pad_left"], result["pad_right"])
short = chop_regions_from_middle([[0, 50]], left=100, right=10, verbose=False)
print(short["left_bins"], short["right_bins"], short["pad_left"], short["pad_right"])
```

- **示例输出：**

```text
[[50, 200]] [[300, 380]] 0 0
[[0, 25]] [[25, 35]] 75 0
```

## 23) `dnakit.epigenome.convert_paired_fastq_to_bisulfite()`

- **作用：** 为亚硫酸氢盐比对转换双端文库一个配对文件的读段。
- **处理方法：** 调用 methylpy.call_mc_pe.convert_reads_pe。读段 1（is_read2=False）中每个 C 变为 T，与 convert_fastq_to_bisulfite 相同；读段 2（is_read2=True）中每个 G 变为 A。两种情形下原始 C（读段 1）或 G（读段 2）的位置都编码在读段名的 ! 之后，名称中的空格也换成 !。输出文件须尚不存在。
- **API：** `dnakit.epigenome.convert_paired_fastq_to_bisulfite(fastq_path, output_path, *, is_read2=False, buffer_line_number=100000, verbose=True)`。
- **输入：** fastq_path 为 FASTQ 文件；output_path 为新文件；is_read2 为布尔值；buffer_line_number。
- **输出：** reads、converted_base、is_read2、input、output 与 provenance。

**示例文件：** `reads.fq`

```text
@r1 x
ACCGTTCA
+
IIIIIIII
@r2
GGATCC
+
IIIIII
```

- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import convert_paired_fastq_to_bisulfite

with tempfile.TemporaryDirectory() as folder:
    for mate in (False, True):
        output = Path(folder) / f"mate_{int(mate) + 1}.fq"
        result = convert_paired_fastq_to_bisulfite("reads.fq", output, is_read2=mate, verbose=False)
        lines = output.read_text().splitlines()
        print(result["converted_base"], result["reads"], lines[1], lines[5])
```

- **示例输出：**

```text
C to T 2 ATTGTTTA GGATTT
G to A 2 ACCATTCA AAATCC
```

## 24) `dnakit.variants.load_maf_table()`

- **作用：** 把 TCGA 风格的 MAF 文件读成含全部列的表。
- **处理方法：** 调用 varcode.maf.load_maf_dataframe。文件须含 19 个标准 MAF 列且顺序一致（列名忽略大小写并改为标准写法），可有更多列。load_maf_variants 只保留 contig、start、ref、alt；本函数还返回基因符号、变异类别、样本条码等所有列。空单元格为 null。max_rows 限制读取行数。
- **API：** `dnakit.variants.load_maf_table(path, *, max_rows=None, verbose=True)`。
- **输入：** path 为本地 MAF 文件；max_rows 可选。
- **输出：** row_count、columns、rows（每行一个字典）与 provenance。

**示例文件：** `sample.maf`

```text
Hugo_Symbol	Entrez_Gene_Id	Center	NCBI_Build	Chromosome	Start_Position	End_Position	Strand	Variant_Classification	Variant_Type	Reference_Allele	Tumor_Seq_Allele1	Tumor_Seq_Allele2	dbSNP_RS	dbSNP_Val_Status	Tumor_Sample_Barcode	Matched_Norm_Sample_Barcode	Match_Norm_Seq_Allele1	Match_Norm_Seq_Allele2
TP53			GRCh38	17	7676154	7676154	+	Missense_Mutation	SNP	C	C	T			TCGA-AA-0001-01A			
```

- **示例代码：**

```python
from dnakit.variants import load_maf_table

result = load_maf_table("sample.maf", verbose=False)
row = result["rows"][0]
print(result["row_count"], len(result["columns"]))
print(row["Hugo_Symbol"], row["Variant_Classification"], row["Start_Position"], row["Entrez_Gene_Id"])
```

- **示例输出：**

```text
1 19
TP53 Missense_Mutation 7676154 None
```

## 25) `dnakit.variants.read_vcf_table()`

- **作用：** 用 varcode 把 VCF 读成原始文本表。
- **处理方法：** 调用 varcode.vcf.read_vcf_into_dataframe。每列都是 VCF 字段的文本：ALT 保持逗号分隔，INFO（仅 include_info=True 时）不解析。给出 sample_names（需要 include_info=True）时加入 FORMAT 和这些样本列；varcode 只按文件顺序给前几个样本列命名，所以名称必须是 VCF 表头中最前面的样本。样本单元格是完整的冒号分隔字符串，保留 GQ、DP 等，这是 read_vcf 所没有的。
- **API：** `dnakit.variants.read_vcf_table(path, *, include_info=False, sample_names=None, verbose=True)`。
- **输入：** path 为 VCF（可 .gz）；include_info 为布尔值；sample_names 可选。
- **输出：** row_count、columns、rows 与 provenance。

**示例文件：** `calls.vcf`

```text
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
1	100	.	A	G,T	50	PASS	DP=10;AF=0.1,0.2	GT:GQ	0/1:35	1/2:40
1	200	rs1	CT	C	60	PASS	DP=20	GT:GQ	0/0:99	0/1:12
```

- **示例代码：**

```python
from dnakit.variants import read_vcf_table

result = read_vcf_table("calls.vcf", include_info=True, sample_names=["S1", "S2"], verbose=False)
print(result["columns"])
first = result["rows"][0]
print(first["ALT"], first["INFO"], first["FORMAT"], first["S1"], first["S2"])
```

- **示例输出：**

```text
['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'S1', 'S2']
G,T DP=10;AF=0.1,0.2 GT:GQ 0/1:35 1/2:40
```

## 26) `dnakit.variants.export_vcf_to_phasing_format()`

- **作用：** 用 bcftools convert 把 VCF 基因型导出为 IMPUTE2 或 SHAPEIT 文本格式。
- **处理方法：** fmt='gensample' 写 <prefix>.gen.gz（每个样本每个位点三个基因型概率）和 <prefix>.samples；'hapsample' 写 <prefix>.hap.gz（每个样本两个单倍型等位，应为定相基因型）和 <prefix>.samples；'haplegendsample' 写 .hap.gz、.legend.gz 和 .samples。位点名为 chrom:pos_ref_alt。以 catch_stdout=False 调用 pysam 自带的 bcftools，否则它不写文件。输出文件须尚不存在。
- **API：** `dnakit.variants.export_vcf_to_phasing_format(source, output_prefix, *, fmt='gensample', max_rows=100000, verbose=True)`。
- **输入：** source 为 VCF 路径；output_prefix 为输出前缀（其目录须存在）；fmt；max_rows 限制返回行数。
- **输出：** files、sample_ids、site_count、rows（主数据文件的行）、rows_truncated 与 provenance。

**示例文件：** `ph.vcf`

```text
##fileformat=VCFv4.2
##contig=<ID=1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
1	100	rs1	A	G	.	PASS	.	GT	0|1	1|1
1	200	rs2	C	T	.	PASS	.	GT	0|0	1|0
```

- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants import export_vcf_to_phasing_format

with tempfile.TemporaryDirectory() as folder:
    for fmt in ("gensample", "hapsample", "haplegendsample"):
        result = export_vcf_to_phasing_format("ph.vcf", Path(folder) / fmt, fmt=fmt, verbose=False)
        print(fmt, [Path(name).name for name in result["files"]], result["sample_ids"])
        print(result["rows"][0])
```

- **示例输出：**

```text
gensample ['gensample.gen.gz', 'gensample.samples'] ['S1', 'S2']
1:100_A_G 1:100_A_G 100 A G 0 1 0 0 0 1
hapsample ['hapsample.hap.gz', 'hapsample.samples'] ['S1', 'S2']
1 1:100_A_G 100 A G 0 1 1 1
haplegendsample ['haplegendsample.hap.gz', 'haplegendsample.legend.gz', 'haplegendsample.samples'] ['S1', 'S2']
0 1 1 1
```

## 主要的错误结论

- gimmemotifs：`motif_score_metric(metric='matches_at_fpr')` 会因 `float()` 无法转换列表而抛出 TypeError（上游返回两个计数）；`metric='roc_auc_xlim'` 把 ROC 曲线点传给实际要求得分列表的函数，抛出 BackendExecutionError。这两行改为 `wrapped`。`peakfile2bedfile` 并不等于 `equal_size_intervals(keep_columns=False)`：它按链处理（非 + 链 start 大于 end）、丢弃 start 小于 0 的行而不是截断为 0，并总是写六列。
- sourmash：`nodegraph.extract_nodegraph_info` 与 khmer 的版本代码完全相同，两者都把表大小取整到百位，台账中“不取整”的说法有误；`nodegraph.calc_expected_collisions` 在同一 FASTA 上与 `expected_kmer_table_collisions` 给出相同的哈希大小、占用桶数和假阳性率（改为 `duplicate-of`）。`compare_all_pairs` 带来丰度加权角相似度和 ANI 输出，`search_databases_with_flat_query` 和 `abund_query` 按分数排序，现有封装都不覆盖，已新封装。
- pairtools：`streaming_dedup_cython` 不能通过 `deduplicate_pairs` 到达（后者固定 scipy 后端），而且结果不同：8 个 pair、max_mismatch 为 2 时在线算法保留 4 个，聚类算法保留 3 个。`unity_norm_scaling` 除以曲线下面积，不等于按窗口平均值归一化的 `norm_window`。
- cooltools：`sandbox.sum_pairwise` 不是“另一组对角线”，而是在对角线行之外为每对不同染色体多返回一行总计；`compute_scaling` 以位置列表为输入并返回原始计数和面积，不是 `compute_contact_scaling` 的原始计数版本。
- tangermeme：`place_motif(mode='insert')` 会把插入后的序列截回原长度，丢掉末端碱基；`beam_substitution` 的束宽和 n_best 参数在 `design_greedy_substitution` 中不存在。
- MACS3：`pileup_from_PN` 的结果只比 `build_shifted_read_pileup(five_shift=0, three_shift=extsize)` 多出一个位于负坐标的前导零点；`se_all_in_one_pileup` 与 `pileup_from_PN_shifted` 输出相同。`CalModel_Heter_AS` 与 `calculate_GQ_heterASsig` 无法通过现有参数到达。
- varcode：`load_maf_dataframe` 与 `read_vcf_into_dataframe` 的表保留全部列和完整的 FORMAT 字符串，现有封装只留 contig、start、ref、alt 或解析后的基因型；上游的 `sample_names` 只按顺序给前几个样本列命名，传入后面的样本名会得到错误的数据，封装因此检查表头。
- pysam：台账称 `bcftools convert` 导出 GEN/HAP 不写文件。实际是默认的 `catch_stdout=True` 吞掉了输出；`catch_stdout=False` 时三种格式都正常写出。
- methylpy：`convert_reads_pe` 对读段 2 做 G 到 A 转换，`convert_fastq_to_bisulfite` 无法到达。
- deeptools：`chopRegionsFromMiddle` 不是 `chopRegions` 的变体，它从外显子中心向两侧切，没有主体部分。

## 验证方式

- 测试文件 `others/tests/unit/test_family_H5.py`，共 55 个测试，离线运行。每个新函数至少有一个手算值或与上游直接调用逐项相等的测试；保留的 51 个 `family-variant` 行中，对每个包都有测试证明所述参数到达上游函数且结果一致，别名行用 `is` 检查为同一对象。
- 所有 26 个新函数都通过了代理工具注册表的参数分析（`analyze_callable`），映射参数的键为字符串。
