# scikit-allel 基因型与等位基因操作

本页记录 P32 相关的真实 DNAKit 接口。函数按所需可选依赖延迟导入；坐标、默认值与参数选择以各接口说明为准。

## `dnakit.variants.allel_p32_functions.write_sites_vcf`

**作用：** Write a sites-only VCF (no sample columns) from variant records with ``allel.write_vcf``.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def write_sites_vcf(output_vcf: str, variants: Sequence[Mapping[str, Any]], *, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Write a sites-only VCF (no sample columns) from variant records with ``allel.write_vcf``.

Each record has ``chrom``, 1-based ``pos``, ``ref`` and ``alt`` (a string or a list of strings) and
optionally ``id``, ``qual`` (missing is written as ``.``) and ``filters`` (a list of failed filter names; none means PASS). Covers ``write_vcf`` for a variant table
(``VariantTable.to_vcf`` and ``VariantChunkedTable.to_vcf`` build the same callset). An existing file is refused.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.variants.allel_p32_functions.haploidify_genotypes`

**作用：** Make pseudo-haplotypes by picking one random allele per sample and variant (``GenotypeArray.haploidify_samples``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def haploidify_genotypes(genotypes: Sequence[Sequence[Sequence[int]]], *, seed: int | None=None, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Make pseudo-haplotypes by picking one random allele per sample and variant (``GenotypeArray.haploidify_samples``).

``seed`` makes the random choice reproducible by seeding NumPy's legacy generator before the call.
Returns a ``variants x samples`` list.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.variants.allel_p32_functions.fill_masked_genotypes`

**作用：** Replace the calls flagged by a ``variants x samples`` boolean mask with ``fill_value`` (``Genotypes.fill_masked``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def fill_masked_genotypes(genotypes: Sequence[Sequence[Sequence[int]]], mask: Sequence[Sequence[bool]], *, fill_value: int=-1, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Replace the calls flagged by a ``variants x samples`` boolean mask with ``fill_value`` (``Genotypes.fill_masked``).
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.variants.allel_p32_functions.locate_genotype_calls`

**作用：** Boolean ``variants x samples`` mask of genotype calls of one kind.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def locate_genotype_calls(genotypes: Sequence[Sequence[Sequence[int]]], *, kind: str='called', allele: int | None=None, call: Sequence[int] | None=None, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Boolean ``variants x samples`` mask of genotype calls of one kind.

``kind`` selects the scikit-allel method: ``called`` (is_called), ``missing`` (is_missing), ``het`` (is_het),
``hom`` (is_hom; with ``allele`` only homozygous for that allele), ``hom_ref`` (is_hom_ref), ``hom_alt``
(is_hom_alt) or ``call`` (is_call; ``call`` is the genotype such as ``[0, 1]`` to match).
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.variants.allel_p32_functions.locate_haplotype_calls`

**作用：** Boolean ``variants x haplotypes`` mask of haplotype calls of one kind.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def locate_haplotype_calls(haplotypes: Sequence[Sequence[int]], *, kind: str='called', allele: int | None=None, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Boolean ``variants x haplotypes`` mask of haplotype calls of one kind.

``kind``: ``called`` (is_called), ``missing`` (is_missing), ``ref`` (is_ref), ``alt`` (is_alt) or ``call``
(is_call; ``allele`` is the allele index to match).
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.variants.allel_p32_functions.count_genotype_calls`

**作用：** Count genotype calls of one kind in total (``axis`` None), per variant (1) or per sample (0).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def count_genotype_calls(genotypes: Sequence[Sequence[Sequence[int]]], *, kind: str='hom', axis: int | None=None, allele: int | None=None, call: Sequence[int] | None=None, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Count genotype calls of one kind in total (``axis`` None), per variant (1) or per sample (0).

``kind``: ``hom`` (count_hom; with ``allele``), ``hom_ref``, ``hom_alt``, ``het`` or ``call``
(count_call; ``call`` is the genotype to match). Note: scikit-allel axis 0 sums over variants, giving one
value per sample; axis 1 gives one value per variant.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.variants.allel_p32_functions.count_haplotype_calls`

**作用：** Count haplotype calls of one kind (``ref`` count_ref, ``alt`` count_alt, ``call`` count_call with ``allele``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def count_haplotype_calls(haplotypes: Sequence[Sequence[int]], *, kind: str='ref', axis: int | None=None, allele: int | None=None, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Count haplotype calls of one kind (``ref`` count_ref, ``alt`` count_alt, ``call`` count_call with ``allele``).

``axis`` None gives the total, 0 one value per haplotype, 1 one value per variant.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.variants.allel_p32_functions.genotype_dosage`

**作用：** Genotypes to the number of reference (``allele='ref'``, to_n_ref) or alternate (``'alt'``, to_n_alt) alleles.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def genotype_dosage(genotypes: Sequence[Sequence[Sequence[int]]], *, allele: str='alt', fill: int=0, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Genotypes to the number of reference (``allele='ref'``, to_n_ref) or alternate (``'alt'``, to_n_alt) alleles.

``fill`` is the value for missing calls.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.variants.allel_p32_functions.genotype_allele_counts_per_call`

**作用：** Per-call allele counts: for each call, how many times each allele ``0..max_allele`` occurs (``to_allele_counts``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def genotype_allele_counts_per_call(genotypes: Sequence[Sequence[Sequence[int]]], *, max_allele: int | None=None, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Per-call allele counts: for each call, how many times each allele ``0..max_allele`` occurs (``to_allele_counts``).
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.variants.allel_p32_functions.genotypes_to_gt_strings`

**作用：** Format genotype calls as VCF GT strings such as ``0/1`` and ``./.`` (``Genotypes.to_gt``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def genotypes_to_gt_strings(genotypes: Sequence[Sequence[Sequence[int]]], *, max_allele: int | None=None, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Format genotype calls as VCF GT strings such as ``0/1`` and ``./.`` (``Genotypes.to_gt``).
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.variants.allel_p32_functions.remap_alleles`

**作用：** Recode allele indices through a per-variant ``variants x old_alleles`` mapping (``map_alleles``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def remap_alleles(values: Sequence[Any], mapping: Sequence[Sequence[int]], *, kind: str='genotypes', max_allele: int | None=None, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Recode allele indices through a per-variant ``variants x old_alleles`` mapping (``map_alleles``).

``kind``: ``genotypes`` (variants x samples x ploidy), ``haplotypes`` (variants x haplotypes) or
``allele_counts`` (variants x alleles; ``max_allele`` limits the new width). A mapped value of -1 drops the allele.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.variants.allel_p32_functions.summarize_alleles_per_variant`

**作用：** Per-variant integer summary of allele counts.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def summarize_alleles_per_variant(allele_counts: Sequence[Sequence[int]], *, kind: str='allelism', verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Per-variant integer summary of allele counts.

``kind='allelism'`` is the number of distinct alleles observed (``allelism``); ``kind='max_allele'`` is the
highest allele index observed (``max_allele``, -1 when nothing is observed).
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.variants.allel_p32_functions.locate_allele_count_sites`

**作用：** Boolean per-variant mask of sites of one kind from allele counts.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def locate_allele_count_sites(allele_counts: Sequence[Sequence[int]], *, kind: str='variant', allele: int=1, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Boolean per-variant mask of sites of one kind from allele counts.

``kind``: ``variant`` (is_variant), ``non_variant``, ``segregating``, ``non_segregating``, ``singleton`` and
``doubleton`` (with ``allele``), ``biallelic`` or ``biallelic_01`` (biallelic for alleles 0 and 1).
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## 可执行示例

**示例代码：**

```python
from dnakit.variants.allel_p32_functions import genotype_dosage
r = genotype_dosage([[[0, 0], [0, 1], [-1, -1]]], verbose=False)
print(r['dosage'])
```

**示例输出：**

```text
[[0, 1, 0]]
```
