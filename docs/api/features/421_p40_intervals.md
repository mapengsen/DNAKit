# 基因组区间与组装

## 作用

为 P39/P40 整合批次提供明确的 DNA 操作接口，返回可序列化结果。

## 处理方法

函数调用其对应的上游操作；显式指定本地文件，新输出路径不得已存在。UCSC/黑名单下载仅在主动调用时联网。验证环境缺少 NUPACK 与 BLAST+，相关科学计算流程尚未实测。Enformer 使用与架构和适配器配置一致的本地权重。

## API

`dnakit.annotation.genomicranges_p40_functions`

| Function | Operation |
| --- | --- |
| `genomicranges_binned_average` | Average the BED score (column 5) of ``source`` over each interval of ``bins`` (GenomicRanges.binned_average). |
| `genomicranges_disjoint_bins` | Assign each interval to a bin so that intervals within one bin do not overlap (GenomicRanges.disjoint_bins). |
| `genomicranges_interval_coverage` | Per-position coverage depth of the intervals on each chromosome (GenomicRanges.coverage). |
| `genomicranges_is_disjoint` | Test whether no two intervals overlap on the same chromosome (and strand) (GenomicRanges.is_disjoint). |
| `genomicranges_match_intervals` | Find, for every query interval, the identical interval in ``source`` (GenomicRanges.match). |
| `genomicranges_narrow_intervals` | Narrow every interval by relative ``start``, ``width`` or ``end`` (GenomicRanges.narrow). |
| `genomicranges_order_intervals` | Indices that sort the intervals by chromosome, strand and position (GenomicRanges.order). |
| `genomicranges_out_of_bound_intervals` | Find intervals that extend past the chromosome length (GenomicRanges.get_out_of_bound_index). |
| `genomicranges_rank_intervals` | Rank of each interval in coordinate order (GenomicRanges.rank). |
| `genomicranges_restrict_intervals` | Restrict intervals to a coordinate window (genomicranges GenomicRanges.restrict). |
| `genomicranges_span_intervals` | Span (minimum start to maximum end) of the intervals per chromosome and strand (GenomicRanges.range). |
| `genomicranges_tile_intervals` | Tile every interval into ``n`` equal sub-intervals or into pieces of ``width`` (GenomicRanges.tile). |

`dnakit.annotation.polars_bio_p40_functions`

| Function | Operation |
| --- | --- |
| `read_bigwig_intervals` | Read the intervals of a BigWig signal track with polars-bio (read_bigwig). |
| `read_msa_records` | Read an A2M or A3M alignment file with polars-bio (read_a2m, read_a3m). |

`dnakit.bioinformatics.bioframe_p40_functions`

| Function | Operation |
| --- | --- |
| `bioframe_fetch_ucsc_table` | Fetch a UCSC assembly table through bioframe's UCSCClient. |
| `bioframe_is_complete_ucsc_string` | Test whether a UCSC-style region string names a chromosome with a known end (bioframe.is_complete_ucsc_string). |
| `bioframe_load_fasta` | Open an indexed FASTA lazily with bioframe.load_fasta and report its records. |

`dnakit.assembly.spades_p40_cli_functions`

| Function | Operation |
| --- | --- |
| `spades_assemble` | Assemble reads with SPAdes in one of its modes (``spades.py`` and the ``*spades.py`` wrappers). |

## 输入

DNA 字符串、JSON 字典/列表、本地基因组文件或模型权重。区间坐标、热力学参数见各函数 docstring；NUPACK 默认材料为 DNA。全部接口支持关键字参数 `verbose`。

## 输出

含完成状态、后端信息和操作结果的字典；文件操作返回输出路径，绘图返回图片路径及对应数据。缺少可选依赖时抛出 `BackendUnavailableError`。各批 manifest 单独记录未验收流程。

**示例代码：**

```python
from dnakit.bioinformatics.bioframe_p40_functions import bioframe_is_complete_ucsc_string
print(bioframe_is_complete_ucsc_string('chr1:10-20', verbose=False)['is_complete'])
```

**示例输出：**

```text
True
```
