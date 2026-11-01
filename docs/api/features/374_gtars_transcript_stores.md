# gtars 转录本存储：坐标转换与查询

本页对应工作清单 P24 包中 gtars `reftx` 的条目。转录本存储是由转录本记录、cdot JSON 文件和 MANE 汇总表构建的单个二进制文件。三个坐标转换条目（c. 转 g.、n. 转 g.、按基因的 c. 转 g.）合并为带 `coordinate_type` 或 `gene` 的一个函数；查询、染色体访问号、链方向条目合并为一个查询函数；添加转录本、染色体映射、读取 cdot 与读取 MANE 是同一个构建函数的参数。测试版本为 gtars 0.10.0。

## 示例文件

- **示例文件：** `c.json`

```text
{"transcripts": {"NM_9.1": {"id": "NM_9.1", "gene_name": "GX", "contig": "NC_000001.11", "strand": 1, "cds_start": 150, "cds_end": 950, "exons": [[100, 400], [600, 1000]]}}}
```

- **示例文件：** `mane.tsv`

```text
#NCBI_GeneID	Ensembl_Gene	HGNC_ID	symbol	name	RefSeq_nuc	RefSeq_prot	Ensembl_nuc	Ensembl_prot	MANE_status	GRCh38_chr	chr_start	chr_end	chr_strand
1	ENSG1	HGNC:1	GX	gene x	NM_9.1	NP_9.1	ENST1	ENSP1	MANE Select	NC_000001.11	100	1000	+
```

## 1) `dnakit.annotation.transcript_store_build()`

- **作用：** 由转录本记录、cdot JSON 文件和 MANE 汇总表构建 gtars 转录本存储文件。
- **处理方法：** 逐条 `TxStoreBuilder.add_transcript`（染色体为 refget 摘要，链为 `+` 或 `-`）、用 `add_chrom_mapping` 登记 contig 名、用 `load_cdot` 读取 cdot 文件（链为 1 或 -1，自动识别 gzip）、用 `load_mane_summary` 读取 NCBI MANE 汇总 TSV，最后 `build(path)`。输出文件必须不存在。
- **API：** `dnakit.annotation.transcript_store_build(output, *, transcripts=None, cdot_files=None, mane_summary=None, chrom_mappings=None, verbose=True)`。
- **输入：** `output`、`transcripts`、`cdot_files`、`mane_summary`、`chrom_mappings`。
- **输出：** `transcript_count`、`transcripts_added`、`cdot_transcripts_loaded`、`mane_rows_loaded`、`output`。
- **示例代码：**

```python
from dnakit.annotation.gtars_p24_reftx_functions import transcript_store_build

import os
CH = 'ixSNGKM4iP6C8F3mkCr_2N0I8RB4mD5U'
def tx(acc, gene, strand, mane):
    return {'accession': acc, 'gene': gene, 'chrom': CH, 'strand': strand, 'cds_start': 150, 'cds_end': 950,
            'exons': [{'start': 100, 'end': 400}, {'start': 600, 'end': 1000}], 'mane': mane}
if os.path.exists('tx.reftx'):
    os.remove('tx.reftx')
r = transcript_store_build('tx.reftx', transcripts=[tx('NM_1.1', 'G1', '+', {'select': True, 'plus_clinical': False}),
                                                   tx('NM_2.1', 'G2', '-', None)], verbose=False)
print(r['transcripts_added'], r['transcript_count'])
if os.path.exists('cdot.reftx'):
    os.remove('cdot.reftx')
c = transcript_store_build('cdot.reftx', cdot_files=['c.json'], mane_summary='mane.tsv', chrom_mappings={'NC_000001.11': CH}, verbose=False)
print(c['cdot_transcripts_loaded'], c['mane_rows_loaded'], c['transcript_count'])
```

- **示例输出：**

```text
2 2
1 2 1
```

## 2) `dnakit.annotation.transcript_record_lookup()`

- **作用：** 在转录本存储中查询转录本（基因、染色体、链、CDS 范围、外显子、MANE 标记）。
- **处理方法：** 先 `ReadonlyTxStore.open`，再 `lookup`、`get_chrom_accession` 与 `get_strand`。访问号不在存储中时 `found` 为 false。
- **API：** `dnakit.annotation.transcript_record_lookup(transcript_store, accession, *, verbose=True)`。
- **输入：** 存储路径、带版本号的访问号。
- **输出：** `found`、`transcript`、`chromosome_accession`、`strand`。
- **示例代码：**

```python
from dnakit.annotation.gtars_p24_reftx_functions import transcript_record_lookup

r = transcript_record_lookup('tx.reftx', 'NM_2.1', verbose=False)
print(r['found'], r['strand'], r['chromosome_accession'])
print(r['transcript']['exons'], r['transcript']['mane'])
print(transcript_record_lookup('tx.reftx', 'NM_404.1', verbose=False)['found'])
```

- **示例输出：**

```text
True - SQ.ixSNGKM4iP6C8F3mkCr_2N0I8RB4mD5U
[{'start': 100, 'end': 400}, {'start': 600, 'end': 1000}] None
False
```

## 3) `dnakit.annotation.transcript_coordinate_to_genomic()`

- **作用：** 把转录本坐标转为基因组坐标：编码（c.）坐标、非编码（n.）坐标，或某基因 MANE Select 转录本的 c. 坐标。
- **处理方法：** `coordinate_type='c'` 调用 `CoordinateMapper.c_to_g_full`，`'n'` 调用 `n_to_g_full`，给出 `gene` 时调用 `c_to_g_by_gene`（存储需含 MANE 标记）。`datum` 为 gtars 对 c. 坐标的可选 datum 参数。
- **API：** `dnakit.annotation.transcript_coordinate_to_genomic(transcript_store, position, *, accession=None, gene=None, coordinate_type='c', datum=None, verbose=True)`。
- **输入：** 存储路径、`position`、`accession` 或 `gene`、`coordinate_type`、`datum`。
- **输出：** `genomic_pos`、`chrom_accession`、`strand`。
- **示例代码：**

```python
from dnakit.annotation.gtars_p24_reftx_functions import transcript_coordinate_to_genomic

a = transcript_coordinate_to_genomic('tx.reftx', 10, accession='NM_1.1', verbose=False)
print(a['genomic_pos'], a['strand'])
b = transcript_coordinate_to_genomic('tx.reftx', 5, accession='NM_2.1', coordinate_type='n', verbose=False)
print(b['genomic_pos'], b['strand'])
g = transcript_coordinate_to_genomic('tx.reftx', 10, gene='G1', verbose=False)
print(g['genomic_pos'], g['accession'])
```

- **示例输出：**

```text
159 +
995 -
159 NM_1.1
```
