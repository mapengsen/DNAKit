# gtars 区间集合、基因注释、信号矩阵、分词器与 LOLA 数据库

本页对应工作清单 P24 包中 gtars 处理 BED 区间集合（`gtars.models`）、基因注释（GTF）、信号矩阵、区域分词器（`gtars.tokenizers`）和 LOLA 区域数据库（`gtars.lola`）的条目。坐标为 0 起点、半开区间。仅参数不同的条目合并：区间集合统计（宽度、平均宽度、核苷酸长度、最大终点、分染色体统计）为一个函数；两个邻接距离条目为带 `kind` 的一个函数；基因模型、注释和分区列表条目为一个概览函数；词元与编号互转为带 `direction` 的一个函数。已被现有函数覆盖的条目（`region_gc_content`、`partition_region_counts`、`sequence_digest_record`）列在文末。测试版本为 gtars 0.10.0。

## 示例文件

- **示例文件：** `a.bed`

```text
chr1	100	200
chr1	150	300
chr1	500	600
chr2	10	50
```

- **示例文件：** `b.bed`

```text
chr1	180	520
chr2	40	90
chr3	5	15
```

- **示例文件：** `u.bed`

```text
chr1	0	1000
chr1	1000	2000
chr2	0	500
```

- **示例文件：** `q.bed`

```text
chr1	900	1000
chr1	20000	20100
```

- **示例文件：** `sm.tsv`

```text
label	c1	c2
chr1_100_500	1.5	2
chr1_1000_2000	0.5	3
```

- **示例文件：** `g.gtf`

```text
chr1	demo	gene	1001	4000	.	+	.	gene_id "G0"; transcript_id "T0"; gene_biotype "protein_coding"; gene_name "N0";
chr1	demo	transcript	1001	4000	.	+	.	gene_id "G0"; transcript_id "T0"; gene_biotype "protein_coding"; gene_name "N0";
chr1	demo	exon	1001	1500	.	+	.	gene_id "G0"; transcript_id "T0"; gene_biotype "protein_coding"; gene_name "N0";
chr1	demo	exon	2501	4000	.	+	.	gene_id "G0"; transcript_id "T0"; gene_biotype "protein_coding"; gene_name "N0";
chr1	demo	gene	6001	9000	.	-	.	gene_id "G1"; transcript_id "T1"; gene_biotype "protein_coding"; gene_name "N1";
chr1	demo	transcript	6001	9000	.	-	.	gene_id "G1"; transcript_id "T1"; gene_biotype "protein_coding"; gene_name "N1";
chr1	demo	exon	6001	6500	.	-	.	gene_id "G1"; transcript_id "T1"; gene_biotype "protein_coding"; gene_name "N1";
chr1	demo	exon	7501	9000	.	-	.	gene_id "G1"; transcript_id "T1"; gene_biotype "protein_coding"; gene_name "N1";
chr1	demo	gene	11001	14000	.	+	.	gene_id "G2"; transcript_id "T2"; gene_biotype "protein_coding"; gene_name "N2";
chr1	demo	transcript	11001	14000	.	+	.	gene_id "G2"; transcript_id "T2"; gene_biotype "protein_coding"; gene_name "N2";
chr1	demo	exon	11001	11500	.	+	.	gene_id "G2"; transcript_id "T2"; gene_biotype "protein_coding"; gene_name "N2";
chr1	demo	exon	12501	14000	.	+	.	gene_id "G2"; transcript_id "T2"; gene_biotype "protein_coding"; gene_name "N2";
```

- **示例文件：** `c.toml`

```text
universe = "u.bed"
tokenizer_type = "ailist"
```

## 1) `dnakit.annotation.region_set_summary()`

- **作用：** 概览 BED 区间集合：区间数、宽度、平均宽度、总核苷酸数、每条染色体的最大终点以及分染色体统计。
- **处理方法：** 用 `gtars.models.RegionSet` 读取文件，调用 `widths`、`mean_region_width`、`get_nucleotide_length`、`get_max_end_per_chr` 与 `chromosome_statistics`。坐标为 0 起点、半开区间。
- **API：** `dnakit.annotation.region_set_summary(regions, *, include_widths=False, max_rows=100000, verbose=True)`。
- **输入：** BED 路径、`include_widths`、`max_rows`。
- **输出：** `region_count`、`mean_region_width`、`total_nucleotides`、`max_end_per_chromosome`、`chromosomes`，可选 `widths`。
- **示例代码：**

```python
from dnakit.annotation.gtars_p24_regions_functions import region_set_summary

r = region_set_summary('a.bed', include_widths=True, verbose=False)
print(r['widths'], r['mean_region_width'], r['total_nucleotides'])
print(r['max_end_per_chromosome'])
print(r['chromosomes'][0]['number_of_regions'])
```

- **示例输出：**

```text
[100, 150, 100, 40] 97.5 390
{'chr1': 600, 'chr2': 50}
3
```

## 2) `dnakit.annotation.region_neighbor_distances()`

- **作用：** 每条染色体上相邻区间的距离：连续区间之间的有符号间隔，或每个区间到最近邻的距离。
- **处理方法：** `kind='signed_gap'` 调用 `RegionSet.neighbor_distances`，`kind='nearest'` 调用 `RegionSet.nearest_neighbors`。区间过少的染色体被跳过，因此结果与输入不一一对应；测量间隔前 gtars 会合并重叠区间。
- **API：** `dnakit.annotation.region_neighbor_distances(regions, *, kind='signed_gap', max_rows=100000, verbose=True)`。
- **输入：** BED 路径、`kind`、`max_rows`。
- **输出：** `distances`、`distance_count`。
- **示例代码：**

```python
from dnakit.annotation.gtars_p24_regions_functions import region_neighbor_distances

print(region_neighbor_distances('a.bed', verbose=False)['distances'])
print(region_neighbor_distances('a.bed', kind='nearest', verbose=False)['distances'])
```

- **示例输出：**

```text
[200]
[0, 0, 200]
```

## 3) `dnakit.annotation.region_distribution_bins()`

- **作用：** 按等宽基因组分箱统计区间中点数。
- **处理方法：** 调用 `RegionSet.distribution(n_bins, chrom_sizes)`。不给 `chrom_sizes` 时箱宽由文件中最大终点决定；给出后每条染色体各切成 `n_bins` 个箱，超出所给长度的区间被忽略。
- **API：** `dnakit.annotation.region_distribution_bins(regions, *, n_bins=250, chrom_sizes=None, verbose=True)`。
- **输入：** BED 路径、`n_bins`、可选 `chrom_sizes`。
- **输出：** `bins`（chrom、start、end、count、bin）。
- **示例代码：**

```python
from dnakit.annotation.gtars_p24_regions_functions import region_distribution_bins

r = region_distribution_bins('a.bed', n_bins=2, chrom_sizes={'chr1': 1000, 'chr2': 100}, verbose=False)
for b in r['bins']:
    print(b)
```

- **示例输出：**

```text
{'chrom': 'chr1', 'start': 0, 'end': 500, 'count': 2, 'bin': 0}
{'chrom': 'chr1', 'start': 500, 'end': 1000, 'count': 1, 'bin': 1}
{'chrom': 'chr2', 'start': 0, 'end': 100, 'count': 1, 'bin': 0}
```

## 4) `dnakit.annotation.region_promoters()`

- **作用：** 围绕每个区间起点生成启动子窗口 `[start - upstream, start + downstream)`。
- **处理方法：** 调用 `RegionSet.promoters(upstream, downstream)`。gtars 从区间起点计算且不考虑链方向，窗口下限为 0。可用 `output` 同时写出新的 BED 文件。
- **API：** `dnakit.annotation.region_promoters(regions, upstream, downstream, *, output=None, max_rows=100000, verbose=True)`。
- **输入：** BED 路径、`upstream`、`downstream`、可选 `output`、`max_rows`。
- **输出：** `regions`、`region_count`、`truncated`、`output`。
- **示例代码：**

```python
from dnakit.annotation.gtars_p24_regions_functions import region_promoters

r = region_promoters('a.bed', 50, 10, verbose=False)
print([(x['chrom'], x['start'], x['end']) for x in r['regions']])
```

- **示例输出：**

```text
[('chr1', 50, 110), ('chr1', 100, 160), ('chr1', 450, 510), ('chr2', 0, 20)]
```

## 5) `dnakit.annotation.region_pairwise_operation()`

- **作用：** 对两个 BED 区间集合做按位置求交、拼接或全对全求交。
- **处理方法：** `operation` 选择 `RegionSet.pintersect`（两集合按位置一一求交）、`concat`（全部区间，不合并）或 `intersect_all`（所有重叠对；gtars 会先合并同一集合内的重叠区间）。
- **API：** `dnakit.annotation.region_pairwise_operation(regions_a, regions_b, *, operation='pintersect', output=None, max_rows=100000, verbose=True)`。
- **输入：** 两个 BED 路径、`operation`、可选 `output`、`max_rows`。
- **输出：** `regions`、`region_count`、`truncated`、`output`。
- **示例代码：**

```python
from dnakit.annotation.gtars_p24_regions_functions import region_pairwise_operation

for op in ('pintersect', 'concat', 'intersect_all'):
    r = region_pairwise_operation('a.bed', 'b.bed', operation=op, verbose=False)
    print(op, r['region_count'], [(x['chrom'], x['start'], x['end']) for x in r['regions']][:3])
```

- **示例输出：**

```text
pintersect 3 [('chr1', 180, 200), ('chr1', 150, 150), ('chr1', 500, 500)]
concat 7 [('chr1', 100, 200), ('chr1', 150, 300), ('chr1', 500, 600)]
intersect_all 3 [('chr1', 180, 300), ('chr1', 500, 520), ('chr2', 40, 50)]
```

## 6) `dnakit.annotation.region_overlap_coefficient()`

- **作用：** 核苷酸水平的重叠系数：共有碱基数除以较小集合的碱基数。
- **处理方法：** 调用 `RegionSet.overlap_coefficient`。
- **API：** `dnakit.annotation.region_overlap_coefficient(regions_a, regions_b, *, verbose=True)`。
- **输入：** 两个 BED 路径。
- **输出：** `overlap_coefficient`，取值 [0, 1]。
- **示例代码：**

```python
from dnakit.annotation.gtars_p24_regions_functions import region_overlap_coefficient

print(region_overlap_coefficient('a.bed', 'b.bed', verbose=False)['overlap_coefficient'])
```

- **示例输出：**

```text
0.4411764705882353
```

## 7) `dnakit.annotation.region_subset_by_overlaps()`

- **作用：** 保留第一个集合中与第二个集合至少一个区间重叠的区间。
- **处理方法：** 调用 `RegionSet.subset_by_overlaps`（为第二个集合建 AIList 索引，用第一个集合的每个区间查询）。保留的区间不做修改。
- **API：** `dnakit.annotation.region_subset_by_overlaps(regions_a, regions_b, *, output=None, max_rows=100000, verbose=True)`。
- **输入：** 两个 BED 路径、可选 `output`、`max_rows`。
- **输出：** `regions`、`region_count`、`input_region_count`。
- **示例代码：**

```python
from dnakit.annotation.gtars_p24_regions_functions import region_subset_by_overlaps

r = region_subset_by_overlaps('a.bed', 'b.bed', verbose=False)
print(r['input_region_count'], r['region_count'])
```

- **示例输出：**

```text
4 4
```

## 8) `dnakit.annotation.region_overlap_indices()`

- **作用：** 对第一个集合的每个区间，给出第二个集合中与之重叠区间的 0 起点索引。
- **处理方法：** 调用 `RegionSet.find_overlaps`。
- **API：** `dnakit.annotation.region_overlap_indices(regions_a, regions_b, *, max_rows=100000, verbose=True)`。
- **输入：** 两个 BED 路径、`max_rows`。
- **输出：** `overlaps`（索引列表的列表）、`overlapping_region_count`。
- **示例代码：**

```python
from dnakit.annotation.gtars_p24_regions_functions import region_overlap_indices

r = region_overlap_indices('a.bed', 'b.bed', verbose=False)
print(r['overlaps'], r['overlapping_region_count'])
```

- **示例输出：**

```text
[[0], [0], [0], [1]] 4
```

## 9) `dnakit.annotation.region_gaps()`

- **作用：** 返回每条染色体上区间之间的空隙，范围从位置 0 到染色体长度。
- **处理方法：** 调用 `RegionSet.gaps(chrom_sizes)`：包括起始空隙、区间间空隙、末尾空隙，`chrom_sizes` 中没有区间的染色体给出一个整条染色体的空隙。其他染色体上的区间被忽略。
- **API：** `dnakit.annotation.region_gaps(regions, chrom_sizes, *, output=None, max_rows=100000, verbose=True)`。
- **输入：** BED 路径、`chrom_sizes`、可选 `output`、`max_rows`。
- **输出：** `regions`、`region_count`。
- **示例代码：**

```python
from dnakit.annotation.gtars_p24_regions_functions import region_gaps

r = region_gaps('a.bed', {'chr1': 1000, 'chr2': 100}, verbose=False)
print([(x['chrom'], x['start'], x['end']) for x in r['regions']])
```

- **示例输出：**

```text
[('chr1', 0, 100), ('chr1', 300, 500), ('chr1', 600, 1000), ('chr2', 0, 10), ('chr2', 50, 100)]
```

## 10) `dnakit.annotation.region_disjoin()`

- **作用：** 把所有区间拆成互不重叠的片段：每个区间边界都成为片段边界。
- **处理方法：** 调用 `RegionSet.disjoin`；片段恰好铺满被覆盖的位置。
- **API：** `dnakit.annotation.region_disjoin(regions, *, output=None, max_rows=100000, verbose=True)`。
- **输入：** BED 路径、可选 `output`、`max_rows`。
- **输出：** `regions`、`region_count`。
- **示例代码：**

```python
from dnakit.annotation.gtars_p24_regions_functions import region_disjoin

r = region_disjoin('a.bed', verbose=False)
print([(x['chrom'], x['start'], x['end']) for x in r['regions']])
```

- **示例输出：**

```text
[('chr1', 100, 150), ('chr1', 150, 200), ('chr1', 200, 300), ('chr1', 500, 600), ('chr2', 10, 50)]
```

## 11) `dnakit.annotation.gene_annotation_summary()`

- **作用：** 由 GTF 文件（或已保存的注释二进制文件）构建 gtars 基因模型、分区列表和 TSS 索引，并报告其规模。
- **处理方法：** 先 `GenomicDistAnnotation.from_gtf`（或 `load_bin`），再 `gene_model`、`partition_list(core_prom, prox_prom, chrom_sizes)` 与 `tss_index`。默认值与 gtars 一致：仅蛋白编码基因、Ensembl 名转 UCSC、核心启动子 100 bp、近端启动子 2000 bp。
- **API：** `dnakit.annotation.gene_annotation_summary(*, gtf=None, annotation_binary=None, core_promoter_bp=100, proximal_promoter_bp=2000, protein_coding_only=True, convert_ensembl_ucsc=True, chrom_sizes=None, verbose=True)`。
- **输入：** `gtf` 或 `annotation_binary`、启动子长度、`protein_coding_only`、`convert_ensembl_ucsc`、可选 `chrom_sizes`。
- **输出：** `n_genes`、`n_exons`、`partition_names`、`tss_count`。
- **示例代码：**

```python
from dnakit.annotation.gtars_p24_regions_functions import gene_annotation_summary

r = gene_annotation_summary(gtf='g.gtf', verbose=False)
print(r['n_genes'], r['n_exons'], r['tss_count'])
print(r['partition_names'])
```

- **示例输出：**

```text
3 6 3
['promoterCore', 'promoterProx', 'exon', 'intron']
```

## 12) `dnakit.annotation.tss_distances()`

- **作用：** 用 gtars `TssIndex` 计算每个查询区间到最近 TSS（或特征）的距离。
- **处理方法：** 索引可由 TSS 位置的 BED 文件（`TssIndex.from_regionset`）、GTF 文件（考虑链方向的 TSS）或已保存的注释二进制文件构建。`signed=False` 调用 `calc_tss_distances`，`signed=True` 调用 `feature_distances`（无特征的染色体为 null）。距离从区间中点算起。
- **API：** `dnakit.annotation.tss_distances(regions, *, tss_regions=None, gtf=None, annotation_binary=None, signed=False, protein_coding_only=True, max_rows=100000, verbose=True)`。
- **输入：** 查询 BED 路径，以及 `tss_regions`、`gtf`、`annotation_binary` 三者之一；`signed`。
- **输出：** `distances`（区间及 `distance`）、`tss_count`。
- **示例代码：**

```python
from dnakit.annotation.gtars_p24_regions_functions import tss_distances

for signed in (False, True):
    r = tss_distances('q.bed', gtf='g.gtf', signed=signed, verbose=False)
    print(signed, [x['distance'] for x in r['distances']])
```

- **示例输出：**

```text
False [50, 9050]
True [50.0, -9050.0]
```

## 13) `dnakit.annotation.signal_matrix_info()`

- **作用：** 加载 gtars 信号矩阵并报告区间数、条件数和条件名称。
- **处理方法：** `SignalMatrix.from_tsv`（首列为 `chrom_start_end` 标签，其余每列一个条件）；`binary=True` 时用 `SignalMatrix.load_bin`。
- **API：** `dnakit.annotation.signal_matrix_info(path, *, binary=False, verbose=True)`。
- **输入：** 矩阵路径、`binary`。
- **输出：** `n_regions`、`n_conditions`、`condition_names`。
- **示例代码：**

```python
from dnakit.annotation.gtars_p24_regions_functions import signal_matrix_info

print(signal_matrix_info('sm.tsv', verbose=False))
```

- **示例输出：**

```text
{'status': 'completed', 'operation': 'signal_matrix_info', 'format': 'tsv', 'n_regions': 2, 'n_conditions': 2, 'condition_names': ['c1', 'c2'], 'provenance': {'package': 'gtars 0.10.0', 'convention': 'TSV first column holds labels chrom_start_end'}}
```

## 14) `dnakit.annotation.region_tokenizer_vocabulary()`

- **作用：** 构建 gtars 区域分词器并返回其词表与特殊符号。
- **处理方法：** `source_type='bed'` 用 `Tokenizer.from_bed`（BED 全集），`'config'` 用 `from_config`（含 `universe` 与 `tokenizer_type = "ailist"` 或 `"bits"` 的 TOML 文件），`'pretrained'` 用 `from_pretrained`（Hugging Face 模型名，需联网）。词元形如 `chr1:0-1000`。
- **API：** `dnakit.annotation.region_tokenizer_vocabulary(source, *, source_type='bed', include_vocabulary=True, verbose=True)`。
- **输入：** `source`、`source_type`、`include_vocabulary`。
- **输出：** `vocab_size`、`special_tokens`、`vocabulary`。
- **示例代码：**

```python
from dnakit.annotation.gtars_p24_regions_functions import region_tokenizer_vocabulary

r = region_tokenizer_vocabulary('u.bed', verbose=False)
print(r['vocab_size'])
print(r['vocabulary'])
print(r['special_tokens']['cls_token'])
c = region_tokenizer_vocabulary('c.toml', source_type='config', include_vocabulary=False, verbose=False)
print(c['vocab_size'])
```

- **示例输出：**

```text
10
{'chr1:0-1000': 0, 'chr1:1000-2000': 1, 'chr2:0-500': 2, '<unk>': 3, '<pad>': 4, '<mask>': 5, '<cls>': 6, '<eos>': 7, '<bos>': 8, '<sep>': 9}
<cls>
10
```

## 15) `dnakit.annotation.region_token_ids()`

- **作用：** 用 gtars 分词器在区域词元与整数编号之间互相转换。
- **处理方法：** `direction='to_ids'` 调用 `convert_tokens_to_ids`（不在词表中的词元映射为未知符编号）；`direction='to_tokens'` 调用 `convert_ids_to_tokens`。
- **API：** `dnakit.annotation.region_token_ids(source, values, *, source_type='bed', direction='to_ids', verbose=True)`。
- **输入：** `source`、`values`、`source_type`、`direction`。
- **输出：** `outputs`，与输入一一对应。
- **示例代码：**

```python
from dnakit.annotation.gtars_p24_regions_functions import region_token_ids

print(region_token_ids('u.bed', ['chr1:1000-2000', 'unknown'], verbose=False)['outputs'])
print(region_token_ids('u.bed', [0, 6], direction='to_tokens', verbose=False)['outputs'])
```

- **示例输出：**

```text
[1, 3]
['chr1:0-1000', '<cls>']
```

## 16) `dnakit.annotation.region_special_tokens_mask()`

- **作用：** 标出给定词元中哪些是 `<cls>`、`<pad>` 之类的特殊符号。
- **处理方法：** 对词元字符串列表调用 `Tokenizer.get_special_tokens_mask`。
- **API：** `dnakit.annotation.region_special_tokens_mask(source, tokens, *, source_type='bed', verbose=True)`。
- **输入：** `source`、`tokens`、`source_type`。
- **输出：** `mask`（true 表示特殊符号）。
- **示例代码：**

```python
from dnakit.annotation.gtars_p24_regions_functions import region_special_tokens_mask

print(region_special_tokens_mask('u.bed', ['<cls>', 'chr1:0-1000', '<pad>'], verbose=False)['mask'])
```

- **示例输出：**

```text
[True, False, True]
```

## 17) `dnakit.annotation.region_database_list_sets()`

- **作用：** 由 BED 文件构建 LOLA 区域数据库（或读取 LOLA 文件夹）并列出其中的区域集。
- **处理方法：** 先 `RegionDB.from_bed_files(files, names)` 或 `RegionDB.from_folder(folder)`，再 `list_region_sets`。名称默认取文件名。
- **API：** `dnakit.annotation.region_database_list_sets(*, bed_files=None, folder=None, names=None, verbose=True)`。
- **输入：** `bed_files` 或 `folder`、可选 `names`。
- **输出：** `set_count`、`sets`（index、name）。
- **示例代码：**

```python
from dnakit.annotation.gtars_p24_regions_functions import region_database_list_sets

r = region_database_list_sets(bed_files=['a.bed', 'b.bed'], verbose=False)
print(r['set_count'], r['sets'])
```

- **示例输出：**

```text
2 [{'index': 0, 'name': 'a.bed'}, {'index': 1, 'name': 'b.bed'}]
```

## 18) `dnakit.annotation.region_database_get_sets()`

- **作用：** 返回 LOLA 区域数据库中所选区域集（按 0 起点索引，默认全部）的区间。
- **处理方法：** 调用 `RegionDB.get_region_sets(indices)`。
- **API：** `dnakit.annotation.region_database_get_sets(*, bed_files=None, folder=None, names=None, indices=None, max_rows=100000, verbose=True)`。
- **输入：** `bed_files` 或 `folder`、`names`、`indices`、`max_rows`。
- **输出：** `sets`（index、name、region_count、regions）。
- **示例代码：**

```python
from dnakit.annotation.gtars_p24_regions_functions import region_database_get_sets

r = region_database_get_sets(bed_files=['a.bed', 'b.bed'], indices=[1], verbose=False)
s = r['sets'][0]
print(s['name'], s['region_count'], [(x['chrom'], x['start'], x['end']) for x in s['regions']])
```

- **示例输出：**

```text
b.bed 3 [('chr1', 180, 520), ('chr2', 40, 90), ('chr3', 5, 15)]
```

## 已被现有函数覆盖的条目

- `gtars.refget.digest_sequence`（sequence_record_from_raw_sequence）：`sequence_digest_record()`。
- `gtars.genomic_distributions.calc_gc_content` 与 `GenomeAssembly.__init__`：`region_gc_content()`（FASTA 基因组）。
- `gtars.genomic_distributions.calc_expected_partitions`：带 `chrom_sizes` 的 `partition_region_counts()`。

未实现：`BinaryGenomeAssembly`（.fab 文件只能由 gtars 命令行 `gtars prep` 生成，Python 接口无法生成）；`GenomicDistAnnotation.load_bin` 的输入文件同样来自 gtars 命令行。若用户已有此类文件，函数可以接收（`annotation_binary`、`binary=True`），但测试不构造该类文件。
