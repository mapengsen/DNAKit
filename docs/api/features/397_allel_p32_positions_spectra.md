# scikit-allel 坐标、频谱与状态区块

本页记录 P32 相关的真实 DNAKit 接口。函数按所需可选依赖延迟导入；坐标、默认值与参数选择以各接口说明为准。

## `dnakit.variants.allel_p32_functions.count_allele_count_sites`

**作用：** Number of sites of one kind: ``variant``, ``non_variant``, ``segregating``, ``non_segregating``, ``singleton``, ``doubleton``.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def count_allele_count_sites(allele_counts: Sequence[Sequence[int]], *, kind: str='variant', allele: int=1, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Number of sites of one kind: ``variant``, ``non_variant``, ``segregating``, ``non_segregating``, ``singleton``, ``doubleton``.

``allele`` applies to ``singleton`` and ``doubleton`` (count_singleton, count_doubleton).
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.variants.allel_p32_functions.feature_coverage_mask`

**作用：** Boolean mask along a contig of the positions covered by features (``FeatureTable.to_mask``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def feature_coverage_mask(features: Sequence[Sequence[int]], size: int, *, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Boolean mask along a contig of the positions covered by features (``FeatureTable.to_mask``).

``features`` is a list of 1-based inclusive ``[start, end]`` pairs and ``size`` the contig length; the mask has
``size + 1`` entries (index 0 unused, as in scikit-allel).
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.variants.allel_p32_functions.locate_sorted_positions`

**作用：** Locate one position (``locate_key``) or several (``locate_keys``) in a sorted position index.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def locate_sorted_positions(positions: Sequence[int], keys: int | Sequence[int], *, strict: bool=True, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Locate one position (``locate_key``) or several (``locate_keys``) in a sorted position index.

A single integer key returns its index; a list returns a boolean mask over ``positions``. With ``strict=False`` missing keys
are ignored instead of raising. Positions must be unique and sorted.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.variants.allel_p32_functions.intersect_sorted_positions`

**作用：** Intersect two sets of sorted genomic positions (``locate_intersection`` masks and ``intersect`` values).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def intersect_sorted_positions(positions: Sequence[int], other: Sequence[int], *, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Intersect two sets of sorted genomic positions (``locate_intersection`` masks and ``intersect`` values).
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.variants.allel_p32_functions.select_positions_in_ranges`

**作用：** Positions that fall inside one or several 1-based inclusive ranges (``locate_range``, ``intersect_range``, ``locate_ranges``, ``intersect_ranges``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def select_positions_in_ranges(positions: Sequence[int], ranges: Sequence[Sequence[int]], *, strict: bool=False, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Positions that fall inside one or several 1-based inclusive ranges (``locate_range``, ``intersect_range``, ``locate_ranges``, ``intersect_ranges``).

``ranges`` is a list of ``[start, stop]`` pairs; one pair equals the single-range methods. Returns the boolean mask over
``positions`` and the selected positions. With ``strict=True`` a range without any position raises.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.variants.allel_p32_functions.overlap_positions_with_ranges`

**作用：** Mask of positions inside any range and mask of ranges that contain at least one position (``locate_intersection_ranges``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def overlap_positions_with_ranges(positions: Sequence[int], ranges: Sequence[Sequence[int]], *, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Mask of positions inside any range and mask of ranges that contain at least one position (``locate_intersection_ranges``).
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.variants.allel_p32_functions.create_allele_index_mapping`

**作用：** Mapping from VCF allele indices (REF, ALT...) to indices in a per-variant allele list such as ancestral/derived (``create_allele_mapping``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def create_allele_index_mapping(ref: Sequence[str], alt: Sequence[Sequence[str]], alleles: Sequence[Sequence[str]], *, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Mapping from VCF allele indices (REF, ALT...) to indices in a per-variant allele list such as ancestral/derived (``create_allele_mapping``).

``ref`` has one base string per variant, ``alt`` one list of alternate alleles per variant (padded
with ``''``), ``alleles`` the target alleles per variant. Alleles absent from the target list map to -1.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.variants.allel_p32_functions.principal_coordinates_analysis`

**作用：** Principal coordinate analysis of a square genetic distance matrix (``allel.pcoa``); returns coordinates and explained variance ratios.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def principal_coordinates_analysis(distance_matrix: Sequence[Sequence[float]], *, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Principal coordinate analysis of a square genetic distance matrix (``allel.pcoa``); returns coordinates and explained variance ratios.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.variants.allel_p32_functions.tabulate_state_transitions`

**作用：** Table of transitions between the listed states along a series (``allel.tabulate_state_transitions``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def tabulate_state_transitions(x: Sequence[int], states: Sequence[int], *, pos: Sequence[int] | None=None, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Table of transitions between the listed states along a series (``allel.tabulate_state_transitions``).

``x`` is the state series, ``states`` the states of interest, ``pos`` optional sorted 1-based positions.
Rows hold lstate, rstate, lidx, ridx and, with ``pos``, lpos and rpos.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.variants.allel_p32_functions.tabulate_state_blocks`

**作用：** Table of contiguous blocks of the listed states along a series (``allel.tabulate_state_blocks``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def tabulate_state_blocks(x: Sequence[int], states: Sequence[int], *, pos: Sequence[int] | None=None, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Table of contiguous blocks of the listed states along a series (``allel.tabulate_state_blocks``).

Rows hold the state, support, start/stop indices, positions with ``pos``, and the minimum and maximum
block length. -1 marks an undefined boundary.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.variants.allel_p32_functions.scale_genotype_dosages`

**作用：** Scale a ``variants x samples`` alternate-allele dosage matrix per variant before PCA.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def scale_genotype_dosages(dosages: Sequence[Sequence[float]], *, method: str='patterson', ploidy: int=2, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Scale a ``variants x samples`` alternate-allele dosage matrix per variant before PCA.

``method='center'`` subtracts the variant mean (``CenterScaler``); ``'patterson'`` also divides by
sqrt(p(1-p)) with p the allele frequency (``PattersonScaler`` with ``ploidy``). Monomorphic variants give undefined values, returned as null.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.variants.allel_p32_functions.fold_site_frequency_spectrum`

**作用：** Fold a one-population site frequency spectrum onto minor-allele counts (``allel.fold_sfs``); ``sample_size`` is the chromosome count n and a shorter spectrum is zero-padded to n + 1 entries. As in scikit-allel, bin k is s[k] + s[n - k] for k < (n + 1) // 2, so use an odd n for a symmetric fold.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def fold_site_frequency_spectrum(spectrum: Sequence[int], sample_size: int, *, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Fold a one-population site frequency spectrum onto minor-allele counts (``allel.fold_sfs``); ``sample_size`` is the chromosome count n and a shorter spectrum is zero-padded to n + 1 entries. As in scikit-allel, bin k is s[k] + s[n - k] for k < (n + 1) // 2, so use an odd n for a symmetric fold.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## 可执行示例

**示例代码：**

```python
from dnakit.variants.allel_p32_functions import intersect_sorted_positions
r = intersect_sorted_positions([1, 3, 5], [3, 4, 5], verbose=False)
print(r['positions'])
```

**示例输出：**

```text
[3, 5]
```
