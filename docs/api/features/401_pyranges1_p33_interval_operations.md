# PyRanges 1.5 区间操作

本页记录 P33 相关的真实 DNAKit 接口。函数按所需可选依赖延迟导入；坐标、默认值与参数选择以各接口说明为准。

验证范围说明：Requires pyranges1 >=1.5. The default validation environment has 1.4.5; examples and supplementary parity tests were run with a 1.5.0 wheel extracted in scratch without changing the environment. RLE and bigWig coverage additionally require runlength.

## `dnakit.annotation.pyranges1_p33_functions.read_pyranges1_intervals`

**作用：** Read WIG (``wig``, pyranges1.read_wig), ENCODE broadPeak (``broadPeak``), gappedPeak

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def read_pyranges1_intervals(path: str | Path, file_format: str, *, max_rows: int=_c.DEFAULT_ROWS, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Read WIG (``wig``, pyranges1.read_wig), ENCODE broadPeak (``broadPeak``), gappedPeak
(``gappedPeak``), Parquet interval tables (``parquet``) or UCSC chain files (``chain``,
pyranges1.read_chain) into interval records.

Rows carry ``Chromosome``, ``Start`` and ``End`` plus the format's own columns (for example
``Value`` for WIG, the peak statistics for the peak formats, chain header fields for chain).
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.annotation.pyranges1_p33_functions.has_pyranges1_tag`

**作用：** Test, for every row of a GTF annotation, whether it carries a tag such as

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def has_pyranges1_tag(annotation_path: str | Path, tag: str, *, tag_column: str='tag', max_rows: int=_c.DEFAULT_ROWS, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Test, for every row of a GTF annotation, whether it carries a tag such as
``Ensembl_canonical`` or ``MANE_Select`` (pyranges1.genes.has_tag).

The GTF is read with ``duplicate_attr=True`` so that all tags of a row are kept. Whole tags are
matched. Returns the flag per row (as ``has_tag``) and the number of rows that carry the tag.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.annotation.pyranges1_p33_functions.select_pyranges1_transcripts`

**作用：** Keep one representative transcript per gene from a GTF annotation

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def select_pyranges1_transcripts(annotation_path: str | Path, *, how: str='canonical', gene_column: str='gene_id', transcript_column: str='transcript_id', tag_column: str='tag', output_path: str | Path | None=None, max_rows: int=_c.DEFAULT_ROWS, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Keep one representative transcript per gene from a GTF annotation
(pyranges1.genes.select_transcripts).

``how``: ``canonical`` (the Ensembl_canonical transcript), ``mane`` (MANE_Select, else
Ensembl_canonical) or ``longest`` (longest span). Genes without the tag fall back to the next
rule. Rows of the chosen transcripts and rows without a transcript (gene lines) are kept.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.annotation.pyranges1_p33_functions.label_pyranges1_utrs`

**作用：** Relabel unsplit ``UTR`` rows of a GTF as ``five_prime_utr`` or ``three_prime_utr``

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def label_pyranges1_utrs(annotation_path: str | Path, *, transcript_column: str='transcript_id', feature_column: str='Feature', output_path: str | Path | None=None, max_rows: int=_c.DEFAULT_ROWS, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Relabel unsplit ``UTR`` rows of a GTF as ``five_prime_utr`` or ``three_prime_utr``
from the transcript's coding extent (pyranges1.genes.label_utrs).

A UTR row wholly before the first CDS base (in transcription order) becomes the 5' UTR, wholly
after the last CDS base the 3' UTR; UTR rows of transcripts without CDS keep ``UTR``.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.annotation.pyranges1_p33_functions.outer_pyranges1_ranges`

**作用：** Outer boundaries of interval groups: the smallest interval that spans each group

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def outer_pyranges1_ranges(intervals_path: str | Path, *, input_format: str='bed', group_by: str | Sequence[str] | None=None, use_strand: str='auto', output_path: str | Path | None=None, max_rows: int=_c.DEFAULT_ROWS, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Outer boundaries of interval groups: the smallest interval that spans each group
(PyRanges.outer_ranges; ``boundaries`` in pyranges 0.1).

Groups are chromosome (and strand when ``use_strand`` is ``auto`` and strands exist) plus the
columns in ``group_by`` (for example a transcript id).
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.annotation.pyranges1_p33_functions.resize_pyranges1_ranges`

**作用：** Resize intervals to a fixed width, keeping the start, end or centre fixed

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def resize_pyranges1_ranges(intervals_path: str | Path, width: int, *, fix: str='center', input_format: str='bed', use_strand: str='auto', output_path: str | Path | None=None, max_rows: int=_c.DEFAULT_ROWS, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Resize intervals to a fixed width, keeping the start, end or centre fixed
(PyRanges.resize_ranges).

``fix`` is ``start``, ``end`` or ``center``; with ``use_strand="auto"`` and strands present,
``start`` and ``end`` mean the 5' and 3' end of each interval.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.annotation.pyranges1_p33_functions.end_positions_pyranges1`

**作用：** 5' end (``end="five"``, PyRanges.five_end) or 3' end (``end="three"``, PyRanges.three_end)

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def end_positions_pyranges1(intervals_path: str | Path, *, end: str='five', input_format: str='bed', group_by: str | Sequence[str] | None=None, extension: int=0, output_path: str | Path | None=None, max_rows: int=_c.DEFAULT_ROWS, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
5' end (``end="five"``, PyRanges.five_end) or 3' end (``end="three"``, PyRanges.three_end)
of intervals, optionally of whole groups such as transcripts.

Without ``group_by`` every interval is reduced to its end; with it, one end per group.
``extension`` widens the reported end by that many bases.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.annotation.pyranges1_p33_functions.join_pyranges1_overlaps`

**作用：** Join intervals of A with the B intervals they overlap, one output row per pair

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def join_pyranges1_overlaps(a_path: str | Path, b_path: str | Path, *, a_format: str='bed', b_format: str='bed', multiple: str='all', strand_behavior: str='auto', join_type: str='inner', match_by: str | Sequence[str] | None=None, contained_intervals_only: bool=False, slack: int=0, suffix: str='_b', overlap_column: str | None=None, output_path: str | Path | None=None, max_rows: int=_c.DEFAULT_ROWS, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Join intervals of A with the B intervals they overlap, one output row per pair
(PyRanges.join_overlaps; ``join`` in pyranges 0.1).

``join_type``: inner, left, right or outer. ``multiple``: all, first or last partner.
``slack`` widens A before testing. B columns get ``suffix``; ``overlap_column`` names an
optional column with the overlap length.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.annotation.pyranges1_p33_functions.interval_lengths_pyranges1`

**作用：** Length (End minus Start) of every interval (PyRanges.lengths).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def interval_lengths_pyranges1(intervals_path: str | Path, *, input_format: str='bed', max_rows: int=_c.DEFAULT_ROWS, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Length (End minus Start) of every interval (PyRanges.lengths).

Returns the first ``max_rows`` lengths in input order plus count, minimum, maximum, mean and sum.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.annotation.pyranges1_p33_functions.map_pyranges1_coordinates`

**作用：** Map intervals between a local reference frame (for example transcript coordinates) and

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def map_pyranges1_coordinates(query_path: str | Path, reference_path: str | Path, reference_column: str, *, direction: str='to_global', query_format: str='bed', reference_format: str='gtf', local_column: str='Chromosome', match_by: str | Sequence[str] | None=None, keep_id: bool=False, keep_chrom: bool=False, keep_loc: bool=False, pep_to_cds: bool=False, split_at_junctions: bool=True, output_path: str | Path | None=None, max_rows: int=_c.DEFAULT_ROWS, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Map intervals between a local reference frame (for example transcript coordinates) and
global (genome) coordinates.

``direction="to_global"`` (PyRanges.map_to_global) takes query intervals whose
``local_column`` (default the chromosome field, holding the reference id) names rows of the
reference and returns genome coordinates; ``keep_id``, ``keep_loc`` and ``pep_to_cds`` apply.
``direction="to_local"`` (PyRanges.map_to_local) takes genome intervals and returns coordinates
local to the reference rows matched by ``match_by``; ``keep_chrom``, ``keep_loc`` and
``split_at_junctions`` apply. ``reference_column`` is the reference id column (``global_on``
or ``ref_on``), for example ``transcript_id`` of a GTF.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.annotation.pyranges1_p33_functions.max_disjoint_pyranges1_overlaps`

**作用：** Largest set of intervals that do not overlap each other (PyRanges.max_disjoint_overlaps).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def max_disjoint_pyranges1_overlaps(intervals_path: str | Path, *, input_format: str='bed', use_strand: str='auto', slack: int=0, match_by: str | Sequence[str] | None=None, output_path: str | Path | None=None, max_rows: int=_c.DEFAULT_ROWS, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Largest set of intervals that do not overlap each other (PyRanges.max_disjoint_overlaps).

``slack`` also separates intervals closer than that many bases; ``match_by`` restricts the
disjointness to rows with equal values in the given columns.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.annotation.pyranges1_p33_functions.slice_pyranges1_ranges`

**作用：** Cut intervals at start/end offsets, genomic or spliced (PyRanges.slice_ranges;

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def slice_pyranges1_ranges(intervals_path: str | Path, *, start: int=0, end: int | None=None, input_format: str='bed', group_by: str | Sequence[str] | None=None, use_strand: str='auto', count_introns: bool=False, output_path: str | Path | None=None, max_rows: int=_c.DEFAULT_ROWS, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Cut intervals at start/end offsets, genomic or spliced (PyRanges.slice_ranges;
``spliced_subsequence`` in pyranges 0.1).

Without ``group_by`` each interval is sliced alone (offsets are within the interval). With
``group_by`` (for example ``transcript_id``) the rows of a group are one spliced sequence:
offsets count transcript bases, introns are skipped unless ``count_introns``. Negative offsets
count from the end; minus-strand groups are counted from their 5' end when ``use_strand`` is
``auto``.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.annotation.pyranges1_p33_functions.make_pyranges1_strand_valid`

**作用：** Replace invalid strand values (anything but ``+`` and ``-``, such as ``.``) by ``+``

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def make_pyranges1_strand_valid(intervals_path: str | Path, *, input_format: str='bed', output_path: str | Path | None=None, max_rows: int=_c.DEFAULT_ROWS, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Replace invalid strand values (anything but ``+`` and ``-``, such as ``.``) by ``+``
(PyRanges.make_strand_valid).

Returns the repaired intervals and how many strand values changed.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.annotation.pyranges1_p33_functions.summarize_pyranges1_intervals`

**作用：** Summary statistics of interval lengths and of the merged coverage (PyRanges.summary).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def summarize_pyranges1_intervals(intervals_path: str | Path, *, input_format: str='bed', verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Summary statistics of interval lengths and of the merged coverage (PyRanges.summary).

Returns the table pyranges1 prints: the number of intervals, total and merged (non-overlapping)
bases, and the distribution of interval lengths.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.annotation.pyranges1_p33_functions.write_pyranges1_intervals`

**作用：** Write intervals as a Parquet table (``parquet``, PyRanges.to_parquet), ENCODE ``narrowPeak``,

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def write_pyranges1_intervals(intervals_path: str | Path, output_path: str | Path, output_format: str, *, input_format: str='bed', value_column: str='Value', verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Write intervals as a Parquet table (``parquet``, PyRanges.to_parquet), ENCODE ``narrowPeak``,
``broadPeak`` or ``gappedPeak`` (PyRanges.to_narrowPeak, to_broadPeak, to_gappedPeak), a Hi-C
``pairs`` file (PyRanges.to_pairs) or a ``bedgraph`` file (PyRanges.to_bedgraph).

The input must carry the columns the output format needs (for example the peak statistics for
the peak formats, which a narrowPeak, broadPeak or gappedPeak input provides). For ``bedgraph``
``value_column`` names the signal column (``Value`` by default, as WIG and bigWig readers
produce it). The new output file is written in the requested format.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.annotation.pyranges1_p33_functions.bigwig_stats_pyranges1`

**作用：** Summarise bigWig signal over each interval (PyRanges.bigwig_stats).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def bigwig_stats_pyranges1(intervals_path: str | Path, bigwig_path: str | Path, *, input_format: str='bed', stats: str | Sequence[str]='mean', missing: float | None=None, pad: int=0, bins: int | None=None, use_strand: str='auto', prefix: str='', output_path: str | Path | None=None, max_rows: int=_c.DEFAULT_ROWS, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Summarise bigWig signal over each interval (PyRanges.bigwig_stats).

``stats`` is any of mean, max, min, sum, std (one name or a list); one column per statistic is
added. ``missing`` is the value for bases without data (None leaves them out; 0.0 gives the
mean over the whole interval), ``pad`` widens each interval, and ``bins`` splits each interval
into equal bins, one column per bin and statistic (a regions x bins matrix). Needs pyBigWig,
runlength and pyarrow.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.annotation.pyranges1_p33_functions.coverage_to_bigwig_pyranges1`

**作用：** Compute coverage of intervals and write it as a bigWig file (PyRanges.to_bigwig).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def coverage_to_bigwig_pyranges1(intervals_path: str | Path, output_path: str | Path, *, input_format: str='bed', chromsizes_path: str | Path | None=None, value_column: str | None=None, divide: bool=False, rpm: bool=True, precomputed: bool=False, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Compute coverage of intervals and write it as a bigWig file (PyRanges.to_bigwig).

The score at each position is the number of intervals covering it (or the sum of
``value_column``); ``rpm`` scales to reads per million, ``divide`` (with ``value_column``)
divides by plain coverage and takes log2, and ``precomputed`` writes the input ranges with their
``value_column`` as they are (they must not overlap). ``chromsizes_path`` is a two-column size
file; without it the largest end per chromosome is used. Needs pyBigWig.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.annotation.pyranges1_p33_functions.coverage_runs_pyranges1`

**作用：** Coverage (or the sum of ``value_column``) as run-length encoded sequences per chromosome and

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def coverage_runs_pyranges1(intervals_path: str | Path, *, input_format: str='bed', value_column: str | None=None, use_strand: str='auto', rpm: bool=False, max_rows: int=_c.DEFAULT_ROWS, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Coverage (or the sum of ``value_column``) as run-length encoded sequences per chromosome and
strand (PyRanges.to_rle).

Returns one row per run (``Chromosome``, ``Strand``, ``start``, ``end``, ``value``; ``rpm``
scales to reads per million). Needs the runlength package.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.annotation.pyranges1_p33_functions.flank_pyranges1_ranges`

**作用：** Regions of ``length`` bases upstream (``side="upstream"``, PyRanges.upstream) or downstream

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def flank_pyranges1_ranges(intervals_path: str | Path, length: int, *, side: str='upstream', gap: int=0, input_format: str='bed', group_by: str | Sequence[str] | None=None, use_strand: str='auto', output_path: str | Path | None=None, max_rows: int=_c.DEFAULT_ROWS, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Regions of ``length`` bases upstream (``side="upstream"``, PyRanges.upstream) or downstream
(``side="downstream"``, PyRanges.downstream) of each interval or interval group.

``gap`` leaves that many bases between the interval and the flank; with strands present
(``use_strand="auto"``) upstream means 5' of the interval.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.annotation.pyranges1_p33_functions.window_pyranges1_ranges`

**作用：** Split intervals into windows of ``window_size`` anchored at the 5' end (PyRanges.window_ranges;

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def window_pyranges1_ranges(intervals_path: str | Path, window_size: int, *, input_format: str='bed', use_strand: str='auto', group_by: str | Sequence[str] | None=None, add_window_id: bool=False, output_path: str | Path | None=None, max_rows: int=_c.DEFAULT_ROWS, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Split intervals into windows of ``window_size`` anchored at the 5' end (PyRanges.window_ranges;
``window`` in pyranges 0.1).

The last window of an interval may be shorter. ``add_window_id`` adds a window number column.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.annotation.pyranges1_p33_functions.group_cumsum_pyranges1`

**作用：** Cumulative interval length within each group, in transcription order for minus-strand groups

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def group_cumsum_pyranges1(intervals_path: str | Path, *, input_format: str='bed', group_by: str | Sequence[str] | None=None, use_strand: str='auto', start_column: str | None=None, end_column: str | None=None, keep_order: bool=True, output_path: str | Path | None=None, max_rows: int=_c.DEFAULT_ROWS, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Cumulative interval length within each group, in transcription order for minus-strand groups
(PyRanges.group_cumsum).

Adds two columns (default ``CumsumStart`` and ``CumsumEnd``, or ``start_column`` and
``end_column``) giving each interval's start and end offset within the concatenated group, for
example exon offsets within a transcript when ``group_by="transcript_id"``.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.annotation.pyranges1_p33_functions.interval_metrics_pyranges1`

**作用：** Relationship metrics for every overlapping pair of A and B intervals (join_overlaps followed

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def interval_metrics_pyranges1(a_path: str | Path, b_path: str | Path, *, a_format: str='bed', b_format: str='bed', metrics: str | Sequence[str]='fraction', denominator: str='first', strand_behavior: str='auto', output_path: str | Path | None=None, max_rows: int=_c.DEFAULT_ROWS, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Relationship metrics for every overlapping pair of A and B intervals (join_overlaps followed
by PyRanges.compute_interval_metrics).

``metrics`` is one name or a list from pyranges1 (for example fraction, jaccard, distance); each
becomes a column. ``denominator`` (first, second, union or min/max as pyranges1 documents) sets
the reference length for ``fraction``.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.annotation.pyranges1_p33_functions.combine_pyranges1_interval_columns`

**作用：** Combine the start/end pairs of overlapping A and B intervals into a new start and end

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def combine_pyranges1_interval_columns(a_path: str | Path, b_path: str | Path, *, function: str='intersect', a_format: str='bed', b_format: str='bed', drop_old_columns: bool=True, strand_behavior: str='auto', output_path: str | Path | None=None, max_rows: int=_c.DEFAULT_ROWS, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Combine the start/end pairs of overlapping A and B intervals into a new start and end
(join_overlaps followed by PyRanges.combine_interval_columns).

``function``: ``intersect`` (the overlap), ``union`` (the span covering both) or ``swap``
(exchange the A and B coordinates). ``drop_old_columns`` removes the B coordinate columns.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.annotation.pyranges1_p33_functions.liftover_pyranges1_ranges`

**作用：** Convert intervals between assemblies with a UCSC chain file (PyRanges.liftover_ranges).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def liftover_pyranges1_ranges(intervals_path: str | Path, chain_path: str | Path, *, input_format: str='bed', min_match: float=0.95, return_unmapped: bool=False, output_path: str | Path | None=None, max_rows: int=_c.DEFAULT_ROWS, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Convert intervals between assemblies with a UCSC chain file (PyRanges.liftover_ranges).

``min_match`` is the minimum fraction of an interval that must map. With ``return_unmapped``
the intervals that failed to map are returned as well (``unmapped_rows``).
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.annotation.pyranges1_p33_functions.drop_pyranges1_duplicate_positions`

**作用：** Drop rows that repeat a genomic position, comparing Chromosome, Start, End and (when

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def drop_pyranges1_duplicate_positions(intervals_path: str | Path, *, input_format: str='bed', keep: str='first', use_strand: str='auto', output_path: str | Path | None=None, max_rows: int=_c.DEFAULT_ROWS, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Drop rows that repeat a genomic position, comparing Chromosome, Start, End and (when
present and ``use_strand`` is ``auto`` or ``yes``) Strand (pyranges 0.1
``drop_duplicate_positions``, here PyRanges.drop_duplicates on the position columns).

``keep``: ``first`` or ``last`` row of each duplicate set, or ``none`` to drop all rows of
duplicated positions.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## 可执行示例

**示例代码：**

```python
from pathlib import Path
from dnakit.annotation.pyranges1_p33_functions import interval_lengths_pyranges1
Path('input.bed').write_text('chr1\t10\t20\nchr1\t30\t35\n')
r = interval_lengths_pyranges1('input.bed', verbose=False)
print(r['lengths'], r['sum'])
```

**示例输出：**

```text
[10, 5] 15
```
