# cooler 与 cooltools 接触图操作

本页记录 P38 相关的真实 DNAKit 接口。函数按所需可选依赖延迟导入；坐标、默认值与参数选择以各接口说明为准。

## `dnakit.epigenome.cooler_cooltools_p38_functions.cooler_region_bin_ids`

**作用：** Map a genomic region to cooler bin IDs, using ``cooler.Cooler.offset`` and ``Cooler.extent``.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def cooler_region_bin_ids(cool_path: str, region: str, *, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Map a genomic region to cooler bin IDs, using ``cooler.Cooler.offset`` and ``Cooler.extent``.

``region`` is a UCSC-style string such as ``chr1`` or ``chr1:1,000,000-2,000,000``. Returns the ID of
the first bin that overlaps the region (``offset``) and the half-open range of bin IDs ``[lo, hi)``
(``extent``), plus the number of bins.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.epigenome.cooler_cooltools_p38_functions.query_cooler_bins`

**作用：** Query the bin table of a cooler (chrom, start, end and balancing weights), using ``Cooler.bins``.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def query_cooler_bins(cool_path: str, *, region: str | None=None, columns: Sequence[str] | None=None, limit: int=1000, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Query the bin table of a cooler (chrom, start, end and balancing weights), using ``Cooler.bins``.

``region`` restricts the table to a UCSC-style region (all bins when omitted); ``columns`` selects columns
such as ``["chrom", "start", "end", "weight"]`` (all stored columns when omitted); ``limit`` caps the
number of returned rows (the total is reported).
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.epigenome.cooler_cooltools_p38_functions.sanitize_cooler_pixels`

**作用：** Sanitise pre-binned pixel records (bin-ID base, upper-triangle orientation), using

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def sanitize_cooler_pixels(pixels_path: str, chromsizes: Mapping[str, int], binsize: int, *, output_path: str | None=None, is_one_based: bool=False, tril_action: str='reflect', sort: bool=True, limit: int=100, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Sanitise pre-binned pixel records (bin-ID base, upper-triangle orientation), using
``cooler.create.sanitize_pixels``.

The input is a tab-separated file with a header containing ``bin1_id``, ``bin2_id`` and optional value
columns such as ``count``. Bin IDs refer to fixed bins of ``binsize`` over ``chromsizes``.
``is_one_based`` shifts the IDs to zero-based, ``tril_action`` is ``reflect`` (flip lower-triangle
records), ``drop`` or ``raise``. The sanitised table is written to ``output_path`` when given.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.epigenome.cooler_cooltools_p38_functions.sanitize_cooler_records`

**作用：** Sanitise genomic contact records and assign bin IDs, using ``cooler.create.sanitize_records``.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def sanitize_cooler_records(records_path: str, chromsizes: Mapping[str, int], binsize: int, *, schema: str='bg2', output_path: str | None=None, is_one_based: bool=False, tril_action: str='reflect', sort: bool=True, limit: int=100, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Sanitise genomic contact records and assign bin IDs, using ``cooler.create.sanitize_records``.

``schema`` selects the column layout: ``bg2`` (chrom1, start1, end1, chrom2, start2, end2, count) or
``pairs`` (chrom1, pos1, chrom2, pos2 ...). The tab-separated input must have a header with those names.
Records are binned with fixed ``binsize`` bins over ``chromsizes``; ``tril_action`` is ``reflect``,
``drop`` or ``raise``. The result has ``bin1_id`` and ``bin2_id`` added and is written to ``output_path``.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.epigenome.cooler_cooltools_p38_functions.detect_cooler_file_type`

**作用：** Detect whether a file or URI is a single-resolution cooler, a multi-resolution cooler or a

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def detect_cooler_file_type(path: str, *, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Detect whether a file or URI is a single-resolution cooler, a multi-resolution cooler or a
single-cell cooler, using ``cooler.fileops.is_cooler``, ``is_multires_file`` and ``is_scool_file``.

``path`` is a file path, optionally with a ``::/group`` suffix. ``file_type`` is ``multires``, ``scool``,
``cooler`` or ``unknown``.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.epigenome.cooler_cooltools_p38_functions.build_dot_calling_kernels`

**作用：** Build the convolution kernels of the dot caller (donut, vertical, horizontal, lowleft, upright).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def build_dot_calling_kernels(*, binsize: int | None=None, kernel_half_width: int | None=None, donut_inner_half_width: int | None=None, kernel_types: Sequence[str]=('donut', 'vertical', 'horizontal', 'lowleft'), verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Build the convolution kernels of the dot caller (donut, vertical, horizontal, lowleft, upright).

With ``binsize`` the cooltools recommendation (``dotfinder.recommend_kernels``) is returned for 4 to 28 kb
resolutions (donut w, p = 7, 4 near 5 kb; 5, 2 near 10 kb; 3, 1 near 25 kb) using the four HiCCUPS kernel
types. With ``kernel_half_width`` (w) and ``donut_inner_half_width`` (p) the kernels of the requested
``kernel_types`` are built directly (``cooltools.lib.numutils.get_kernel``).
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.epigenome.cooler_cooltools_p38_functions.validate_dot_calling_kernels`

**作用：** Validate user-supplied dot-calling kernels, using ``cooltools.api.dotfinder.is_compatible_kernels``.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def validate_dot_calling_kernels(kernels: Mapping[str, Sequence[Sequence[int]]], binsize: int, *, max_nans_tolerated: int=1, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Validate user-supplied dot-calling kernels, using ``cooltools.api.dotfinder.is_compatible_kernels``.

``kernels`` maps names to square integer matrices of the same odd size greater than 1; ``max_nans_tolerated``
must not exceed the kernel width. Returns ``valid`` and, when invalid, the upstream ``reason``.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.epigenome.cooler_cooltools_p38_functions.calculate_dot_calling_stage`

**作用：** Run the lambda-binning dot-calling procedure of cooltools up to one stage and write that stage's table.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def calculate_dot_calling_stage(cool_path: str, expected_path: str, output_path: str, *, stage: str='pixels', weight: str='weight', expected_value_col: str='balanced.avg', kernel_half_width: int | None=None, donut_inner_half_width: int | None=None, max_loci_separation: int=10000000, tile_size: int=5000000, n_lambda_bins: int=40, lambda_bin_fdr: float=0.1, max_nans_tolerated: int=1, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Run the lambda-binning dot-calling procedure of cooltools up to one stage and write that stage's table.

``stage`` selects the output: ``scored`` (locally adjusted expected per kernel and Poisson-ready counts for
every pixel of every tile: ``score_tile`` with ``get_adjusted_expected_tile_some_nans``), ``histogram``
(lambda-bin by observed-count histogram per kernel: ``scoring_and_histogramming_step``), ``thresholds`` (FDR
threshold per lambda bin and kernel: ``determine_thresholds``), ``pixels`` (pixels passing the thresholds for
every kernel: ``scoring_and_extraction_step`` / ``extract_scored_pixels``) or ``qvalues`` (those pixels
annotated with q-values: ``annotate_pixels_with_qvalues``, then genomic coordinates chrom1 ... end2 as cooltools ``dots`` does). ``cool_path`` must be balanced; ``expected_path``
is a tab-separated cis expected table (region1, region2, dist, balanced.avg). Kernels are the cooltools
recommendation for the resolution unless ``kernel_half_width`` (and ``donut_inner_half_width``) are given.
The table goes to ``output_path`` (histogram and thresholds in long format with a ``kernel`` column).
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.epigenome.cooler_cooltools_p38_functions.cluster_dot_pixels`

**作用：** Cluster enriched pixels into dot centroids with Birch clustering, using

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def cluster_dot_pixels(pixels_path: str, output_path: str, *, clustering_radius: int=20000, count_column: str='count', verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Cluster enriched pixels into dot centroids with Birch clustering, using
``cooltools.api.dotfinder.clustering_step`` (which calls ``clust_2D_pixels``).

``pixels_path`` is a tab-separated table of enriched pixels annotated with chrom1, start1, chrom2, start2 and
the observed count column (the ``qvalues`` output of ``calculate_dot_calling_stage`` joined to genomic
coordinates). Clustering is done per chromosome (or per ``region`` column when present) on the start
coordinates with ``clustering_radius`` in bp; the brightest pixel of every cluster is kept with its
centroid (``cstart1``, ``cstart2``), ``c_label`` and ``c_size``.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.epigenome.cooler_cooltools_p38_functions.filter_dot_centroids`

**作用：** Filter dot centroids by the HiCCUPS enrichment criteria, using

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def filter_dot_centroids(centroids_path: str, output_path: str, *, count_column: str='count', enrichment_factor_vh: float=1.5, enrichment_factor_d_and_ll: float=1.75, enrichment_factor_d_or_ll: float=2.0, fdr_orphan_threshold: float=0.02, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Filter dot centroids by the HiCCUPS enrichment criteria, using
``cooltools.api.dotfinder.cluster_filtering_hiccups``.

The tab-separated input is the output of ``cluster_dot_pixels`` scored with the four HiCCUPS kernels (donut,
vertical, horizontal, lowleft; columns ``la_exp.<kernel>.value`` and ``la_exp.<kernel>.qval``). Defaults are
those of cooltools (1.5, 1.75, 2.0 and 0.02).
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.epigenome.cooler_cooltools_p38_functions.count_block_pixels`

**作用：** Count total and masked (bad) pixels in a rectangular block of a contact map, using

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def count_block_pixels(width: int, height: int, *, bad_bins_x: int=0, bad_bins_y: int=0, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Count total and masked (bad) pixels in a rectangular block of a contact map, using
``cooltools.api.expected.count_all_pixels_per_block`` and ``count_bad_pixels_per_block``.

``width`` and ``height`` are the block sides in bins; ``bad_bins_x`` and ``bad_bins_y`` count masked bins on
each side. Returns total, bad and valid pixel counts.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.epigenome.cooler_cooltools_p38_functions.combine_logbinned_expected`

**作用：** Combine by-region log-binned expected into a genome-wide P(s) curve and slope, using

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def combine_logbinned_expected(binned_expected_path: str, output_path: str, *, binned_slope_path: str | None=None, pc_name: str='balanced.avg', spread_funcs: str='logstd', spread_funcs_slope: str='std', minmax_drop_bins: int=2, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Combine by-region log-binned expected into a genome-wide P(s) curve and slope, using
``cooltools.api.expected.combine_binned_expected``.

``binned_expected_path`` is the tab-separated output of ``logbin_expected`` (columns region1, region2,
dist_bin_start, dist_bin_end, n_valid, balanced.avg ...). With ``binned_slope_path`` (the slope table of
``logbin_expected`` as a tab-separated file without index) the spread of the slopes is estimated too. ``spread_funcs`` is ``logstd``, ``std`` or
``minmax``. The combined P(s) table is written to ``output_path``.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.epigenome.cooler_cooltools_p38_functions.interpolate_logbinned_expected`

**作用：** Interpolate log-binned expected back to per-diagonal values, using

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def interpolate_logbinned_expected(expected_path: str, binned_expected_path: str, output_path: str, *, columns: Sequence[str]=('balanced.avg',), kind: str='quadratic', by_region: bool=True, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Interpolate log-binned expected back to per-diagonal values, using
``cooltools.api.expected.interpolate_expected``.

``expected_path`` is the per-diagonal expected (region1, region2, dist ...), ``binned_expected_path`` the
log-binned or combined expected (the tab-separated ``logbin_expected`` or ``combine_logbinned_expected``
table with ``dist.avg``). ``kind`` is a scipy ``interp1d`` kind. With
``by_region`` false one combined curve is used for every region.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.epigenome.cooler_cooltools_p38_functions.snip_expected_around_features`

**作用：** Cut expected-matrix snippets around features, using ``cooltools.api.snipping.ExpectedSnipper``

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def snip_expected_around_features(cool_path: str, expected_path: str, features_path: str, *, flank: int=100000, expected_value_col: str='balanced.avg', min_diag: int=2, limit: int=20, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Cut expected-matrix snippets around features, using ``cooltools.api.snipping.ExpectedSnipper``
(``select`` and ``snip``).

``features_path`` is a BED-like tab-separated file with a header (chrom, start, end) of features to centre
on; every snippet spans the feature midpoint plus and minus ``flank`` bp, aligned to bins. Diagonals below
``min_diag`` are NaN (reported as null); windows leaving the chromosome are all null. Features on different
chromosomes use cis expected of that chromosome (the whole chromosome is the region).
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.epigenome.cooler_cooltools_p38_functions.iterative_correction_rectangular`

**作用：** Iterative correction (matrix balancing) of an asymmetric, possibly rectangular contact matrix, using

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def iterative_correction_rectangular(matrix: Sequence[Sequence[float]], *, max_iter: int=1000, tol: float=1e-05, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Iterative correction (matrix balancing) of an asymmetric, possibly rectangular contact matrix, using
``cooltools.lib.numutils.iterative_correction_asymmetric``.

The matrix is rows by columns (for example two chromosomes in trans). Returns the corrected matrix, the bias
vector of the rows (``row_bias``, length = rows) and of the columns (``column_bias``) and the convergence
report. Note: upstream names its two returned vectors in the order column, row.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## 可执行示例

**示例代码：**

```python
from dnakit.epigenome.cooler_cooltools_p38_functions import build_dot_calling_kernels
r = build_dot_calling_kernels(binsize=10000, verbose=False)
print(r['width'], sorted(r['kernels']))
```

**示例输出：**

```text
11 ['donut', 'horizontal', 'lowleft', 'vertical']
```
