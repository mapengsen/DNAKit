# scikit-allel genotype and allele calls

These DNAKit interfaces call their documented upstream implementations. Optional dependencies are imported only when needed. Coordinate conventions, defaults and parameter choices are specified per function.

## `dnakit.variants.allel_p32_functions.write_sites_vcf`

**Purpose:** Write a sites-only VCF (no sample columns) from variant records with ``allel.write_vcf``.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def write_sites_vcf(output_vcf: str, variants: Sequence[Mapping[str, Any]], *, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Write a sites-only VCF (no sample columns) from variant records with ``allel.write_vcf``.

Each record has ``chrom``, 1-based ``pos``, ``ref`` and ``alt`` (a string or a list of strings) and
optionally ``id``, ``qual`` (missing is written as ``.``) and ``filters`` (a list of failed filter names; none means PASS). Covers ``write_vcf`` for a variant table
(``VariantTable.to_vcf`` and ``VariantChunkedTable.to_vcf`` build the same callset). An existing file is refused.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.variants.allel_p32_functions.haploidify_genotypes`

**Purpose:** Make pseudo-haplotypes by picking one random allele per sample and variant (``GenotypeArray.haploidify_samples``).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def haploidify_genotypes(genotypes: Sequence[Sequence[Sequence[int]]], *, seed: int | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Make pseudo-haplotypes by picking one random allele per sample and variant (``GenotypeArray.haploidify_samples``).

``seed`` makes the random choice reproducible by seeding NumPy's legacy generator before the call.
Returns a ``variants x samples`` list.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.variants.allel_p32_functions.fill_masked_genotypes`

**Purpose:** Replace the calls flagged by a ``variants x samples`` boolean mask with ``fill_value`` (``Genotypes.fill_masked``).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def fill_masked_genotypes(genotypes: Sequence[Sequence[Sequence[int]]], mask: Sequence[Sequence[bool]], *, fill_value: int=-1, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Replace the calls flagged by a ``variants x samples`` boolean mask with ``fill_value`` (``Genotypes.fill_masked``).
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.variants.allel_p32_functions.locate_genotype_calls`

**Purpose:** Boolean ``variants x samples`` mask of genotype calls of one kind.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def locate_genotype_calls(genotypes: Sequence[Sequence[Sequence[int]]], *, kind: str='called', allele: int | None=None, call: Sequence[int] | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Boolean ``variants x samples`` mask of genotype calls of one kind.

``kind`` selects the scikit-allel method: ``called`` (is_called), ``missing`` (is_missing), ``het`` (is_het),
``hom`` (is_hom; with ``allele`` only homozygous for that allele), ``hom_ref`` (is_hom_ref), ``hom_alt``
(is_hom_alt) or ``call`` (is_call; ``call`` is the genotype such as ``[0, 1]`` to match).
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.variants.allel_p32_functions.locate_haplotype_calls`

**Purpose:** Boolean ``variants x haplotypes`` mask of haplotype calls of one kind.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def locate_haplotype_calls(haplotypes: Sequence[Sequence[int]], *, kind: str='called', allele: int | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Boolean ``variants x haplotypes`` mask of haplotype calls of one kind.

``kind``: ``called`` (is_called), ``missing`` (is_missing), ``ref`` (is_ref), ``alt`` (is_alt) or ``call``
(is_call; ``allele`` is the allele index to match).
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.variants.allel_p32_functions.count_genotype_calls`

**Purpose:** Count genotype calls of one kind in total (``axis`` None), per variant (1) or per sample (0).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def count_genotype_calls(genotypes: Sequence[Sequence[Sequence[int]]], *, kind: str='hom', axis: int | None=None, allele: int | None=None, call: Sequence[int] | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Count genotype calls of one kind in total (``axis`` None), per variant (1) or per sample (0).

``kind``: ``hom`` (count_hom; with ``allele``), ``hom_ref``, ``hom_alt``, ``het`` or ``call``
(count_call; ``call`` is the genotype to match). Note: scikit-allel axis 0 sums over variants, giving one
value per sample; axis 1 gives one value per variant.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.variants.allel_p32_functions.count_haplotype_calls`

**Purpose:** Count haplotype calls of one kind (``ref`` count_ref, ``alt`` count_alt, ``call`` count_call with ``allele``).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def count_haplotype_calls(haplotypes: Sequence[Sequence[int]], *, kind: str='ref', axis: int | None=None, allele: int | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Count haplotype calls of one kind (``ref`` count_ref, ``alt`` count_alt, ``call`` count_call with ``allele``).

``axis`` None gives the total, 0 one value per haplotype, 1 one value per variant.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.variants.allel_p32_functions.genotype_dosage`

**Purpose:** Genotypes to the number of reference (``allele='ref'``, to_n_ref) or alternate (``'alt'``, to_n_alt) alleles.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def genotype_dosage(genotypes: Sequence[Sequence[Sequence[int]]], *, allele: str='alt', fill: int=0, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Genotypes to the number of reference (``allele='ref'``, to_n_ref) or alternate (``'alt'``, to_n_alt) alleles.

``fill`` is the value for missing calls.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.variants.allel_p32_functions.genotype_allele_counts_per_call`

**Purpose:** Per-call allele counts: for each call, how many times each allele ``0..max_allele`` occurs (``to_allele_counts``).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def genotype_allele_counts_per_call(genotypes: Sequence[Sequence[Sequence[int]]], *, max_allele: int | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Per-call allele counts: for each call, how many times each allele ``0..max_allele`` occurs (``to_allele_counts``).
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.variants.allel_p32_functions.genotypes_to_gt_strings`

**Purpose:** Format genotype calls as VCF GT strings such as ``0/1`` and ``./.`` (``Genotypes.to_gt``).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def genotypes_to_gt_strings(genotypes: Sequence[Sequence[Sequence[int]]], *, max_allele: int | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Format genotype calls as VCF GT strings such as ``0/1`` and ``./.`` (``Genotypes.to_gt``).
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.variants.allel_p32_functions.remap_alleles`

**Purpose:** Recode allele indices through a per-variant ``variants x old_alleles`` mapping (``map_alleles``).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def remap_alleles(values: Sequence[Any], mapping: Sequence[Sequence[int]], *, kind: str='genotypes', max_allele: int | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Recode allele indices through a per-variant ``variants x old_alleles`` mapping (``map_alleles``).

``kind``: ``genotypes`` (variants x samples x ploidy), ``haplotypes`` (variants x haplotypes) or
``allele_counts`` (variants x alleles; ``max_allele`` limits the new width). A mapped value of -1 drops the allele.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.variants.allel_p32_functions.summarize_alleles_per_variant`

**Purpose:** Per-variant integer summary of allele counts.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def summarize_alleles_per_variant(allele_counts: Sequence[Sequence[int]], *, kind: str='allelism', verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Per-variant integer summary of allele counts.

``kind='allelism'`` is the number of distinct alleles observed (``allelism``); ``kind='max_allele'`` is the
highest allele index observed (``max_allele``, -1 when nothing is observed).
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.variants.allel_p32_functions.locate_allele_count_sites`

**Purpose:** Boolean per-variant mask of sites of one kind from allele counts.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def locate_allele_count_sites(allele_counts: Sequence[Sequence[int]], *, kind: str='variant', allele: int=1, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Boolean per-variant mask of sites of one kind from allele counts.

``kind``: ``variant`` (is_variant), ``non_variant``, ``segregating``, ``non_segregating``, ``singleton`` and
``doubleton`` (with ``allele``), ``biallelic`` or ``biallelic_01`` (biallelic for alleles 0 and 1).
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## Executed example

**Example code:**

```python
from dnakit.variants.allel_p32_functions import genotype_dosage
r = genotype_dosage([[[0, 0], [0, 1], [-1, -1]]], verbose=False)
print(r['dosage'])
```

**Example output:**

```text
[[0, 1, 0]]
```
