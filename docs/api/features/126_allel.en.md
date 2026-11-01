# scikit-allel: population-genetic plots and genotype stores

This group of functions in `dnakit.population` wraps [scikit-allel](https://github.com/cggh/scikit-allel) for population-genetic figures, runs of homozygosity and VCF-to-Zarr conversion. Fifteen functions draw site frequency spectra, haplotype paintings and frequencies, pairwise distance and linkage-disequilibrium matrices and a variant locator with scikit-allel's own plotting functions. They write a new PNG file (150 dpi, Agg backend, an existing file is never overwritten) and return its path together with the numbers that were drawn, never a plotting object. `roh_poissonhmm()` calls scikit-allel's Poisson hidden Markov model for runs of homozygosity, and `allel_vcf_to_zarr()` copies the arrays of a VCF file into a Zarr store with [zarr](https://github.com/zarr-developers/zarr-python). The plot functions draw what you pass in: they do not compute the spectra, LD or distances themselves.

scikit-allel is released under the MIT licence (checked in the installed package metadata). The dependencies are optional and are imported only when a function is called; a missing package raises an error with the install command. The examples import from `dnakit.population`.

## Example data

Only the last function reads an input file: a three-variant, three-sample VCF. All other examples pass small lists directly and write PNG files to a temporary directory.

**Example file:** `mini.vcf`

```text
##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3
chr1	100	v1	A	G	50	PASS	.	GT	0/0	0/1	1/1
chr1	250	v2	C	T	60	PASS	.	GT	0/1	0/0	0/1
chr1	400	v3	G	A	70	PASS	.	GT	1/1	0/1	./.
```

## 1) `dnakit.population.fig_voight_painting()`

- **Purpose:** Draw scikit-allel's two-flank Voight painting figure, with EHH decay curves, around one variant of a haplotype matrix.
- **Method:** Validate the matrix, then call `allel.fig_voight_painting()` and save the figure as a PNG. The painting of each flank is recomputed with `allel.voight_painting()`: the left flank on `matrix[:index + 1][::-1]` and the right flank on `matrix[index:]`, so row 0 is the variant at `index` and later rows move away from it.
- **API:** `fig_voight_painting(haplotypes: 'list[list[int]]', output_path: 'str', *, index: 'int | None' = None, palette: 'str' = 'colorblind', height_factor: 'float' = 0.01, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package scikit-allel (extra `dnakit[bio-scikit-allel]`, which also installs scipy and seaborn) and matplotlib (extra `dnakit[bio-plots]`). Works offline.
- **Input:** `haplotypes` is a non-empty variants x haplotypes matrix of non-negative whole numbers (biallelic data is 0 and 1). `output_path` must be a new `.png` path. `index` is an integer from 0 to the number of variants minus 1; the default is the middle variant (`n_variants // 2`). `palette` is a seaborn palette name and `height_factor` scales the painting height. It also needs seaborn.
- **Output:** `path` (the PNG) and `data` with `index`, `n_variants`, `n_haplotypes`, `palette`, `left_painting` and `right_painting` (variants from the index outward x haplotypes; 0 marks haplotypes that are no longer shared with another one).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.population import fig_voight_painting

out = Path(tempfile.mkdtemp())
haplotypes = [
    [0, 0, 1, 1, 0, 1],
    [0, 0, 1, 1, 1, 0],
    [0, 0, 1, 1, 0, 0],
    [0, 0, 0, 0, 1, 1],
    [0, 0, 0, 0, 1, 1],
]
result = fig_voight_painting(haplotypes, str(out / "voight.png"), index=2, verbose=False)
data = result["data"]
print(data["index"], data["n_variants"], data["n_haplotypes"])
print("left ", data["left_painting"])
print("right", data["right_painting"])
```

- **Example output:**

```text
2 5 6
left  [[1, 1, 1, 1, 2, 2], [1, 1, 1, 0, 2, 2], [1, 1, 0, 0, 2, 2]]
right [[1, 1, 1, 1, 2, 2], [1, 1, 3, 3, 2, 2], [1, 1, 3, 3, 2, 2]]
```

## 2) `dnakit.population.plot_haplotype_frequencies()`

- **Purpose:** Draw the frequencies of the distinct haplotypes of a region as one bar of coloured spans.
- **Method:** Wrap the matrix in `allel.HaplotypeArray`, count identical haplotypes with `distinct_counts()` and draw them with `allel.plot_haplotype_frequencies()`. Haplotypes seen more than once get a palette colour, singletons get `singleton_color`.
- **API:** `plot_haplotype_frequencies(haplotypes: 'list[list[int]]', output_path: 'str', *, palette: 'str' = 'Paired', singleton_color: 'str' = 'w', verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package scikit-allel (extra `dnakit[bio-scikit-allel]`, which also installs scipy and seaborn) and matplotlib (extra `dnakit[bio-plots]`). Works offline.
- **Input:** `haplotypes` is a variants x haplotypes matrix of non-negative whole numbers; two columns are the same haplotype when all their alleles match. `output_path` must be a new `.png` path. `palette` (default `Paired`) and `singleton_color` (default `w`, white) are matplotlib/seaborn colour settings. It also needs seaborn.
- **Output:** `path` and `data` with `n_variants`, `n_haplotypes`, `haplotype_counts` (the count of each distinct haplotype, in scikit-allel's `distinct_counts` order), `singleton_count`, `palette` and `singleton_color`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.population import plot_haplotype_frequencies

out = Path(tempfile.mkdtemp())
haplotypes = [
    [0, 0, 1, 1, 0, 1],
    [0, 0, 1, 1, 1, 0],
    [0, 0, 1, 1, 0, 0],
    [0, 0, 0, 0, 1, 1],
    [0, 0, 0, 0, 1, 1],
]
result = plot_haplotype_frequencies(haplotypes, str(out / "freq.png"), verbose=False)
data = result["data"]
print(data["haplotype_counts"], data["singleton_count"])
print(data["n_variants"], data["n_haplotypes"])
```

- **Example output:**

```text
[2, 2, 1, 1] 2
5 6
```

## 3) `dnakit.population.plot_joint_sfs()`

- **Purpose:** Draw a joint site frequency spectrum of two populations as a log-scaled image.
- **Method:** Call `allel.plot_joint_sfs()` on the matrix and save the figure as a PNG. The default colour scale is logarithmic (`LogNorm`), so zero cells are left blank. The drawn matrix is the input transposed, so population 1 is on the x axis.
- **API:** `plot_joint_sfs(spectrum: 'list[list[float]]', output_path: 'str', *, imshow_kwargs: 'dict[str, Any] | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages scikit-allel (extra `dnakit[bio-scikit-allel]`) and matplotlib (extra `dnakit[bio-plots]`). Works offline.
- **Input:** `spectrum` is a non-empty 2-D list of non-negative numbers: rows are the derived allele count in population 1 and columns the count in population 2. `output_path` must be a new `.png` path. `imshow_kwargs` is passed to matplotlib's `imshow`.
- **Output:** `path` and `data` with `spectrum`, `shape`, `drawn_matrix` (the transposed matrix as drawn) and the axis labels `xlabel` and `ylabel` (derived allele count of each population).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.population import plot_joint_sfs

out = Path(tempfile.mkdtemp())
spectrum = [[0, 3, 1], [4, 2, 1], [1, 1, 0]]
result = plot_joint_sfs(spectrum, str(out / "joint.png"), verbose=False)
data = result["data"]
print(data["shape"], data["xlabel"], "|", data["ylabel"])
print(data["drawn_matrix"])
print(result["path"].endswith("joint.png"))
```

- **Example output:**

```text
[3, 3] derived allele count (population 1) | derived allele count (population 2)
[[0.0, 4.0, 1.0], [3.0, 2.0, 1.0], [1.0, 1.0, 0.0]]
True
```

## 4) `dnakit.population.plot_joint_sfs_folded()`

- **Purpose:** Draw a folded joint site frequency spectrum (minor allele counts) of two populations.
- **Method:** Call `allel.plot_joint_sfs_folded()` on the matrix and save the figure as a PNG. The colour scale is logarithmic, as for the unfolded plot. The drawn matrix is the input transposed, so population 1 is on the x axis.
- **API:** `plot_joint_sfs_folded(spectrum: 'list[list[float]]', output_path: 'str', *, imshow_kwargs: 'dict[str, Any] | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages scikit-allel (extra `dnakit[bio-scikit-allel]`) and matplotlib (extra `dnakit[bio-plots]`). Works offline.
- **Input:** `spectrum` is a non-empty 2-D list of non-negative numbers: rows are the derived allele count in population 1 and columns the count in population 2. `output_path` must be a new `.png` path. `imshow_kwargs` is passed to matplotlib's `imshow`.
- **Output:** `path` and `data` with `spectrum`, `shape`, `drawn_matrix` (the transposed matrix as drawn) and the axis labels `xlabel` and `ylabel` (minor allele count of each population).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.population import plot_joint_sfs_folded

out = Path(tempfile.mkdtemp())
spectrum = [[0, 3, 1], [4, 2, 1], [1, 1, 0]]
result = plot_joint_sfs_folded(spectrum, str(out / "joint.png"), verbose=False)
data = result["data"]
print(data["shape"], data["xlabel"], "|", data["ylabel"])
print(data["drawn_matrix"])
print(result["path"].endswith("joint.png"))
```

- **Example output:**

```text
[3, 3] minor allele count (population 1) | minor allele count (population 2)
[[0.0, 4.0, 1.0], [3.0, 2.0, 1.0], [1.0, 1.0, 0.0]]
True
```

## 5) `dnakit.population.plot_joint_sfs_folded_scaled()`

- **Purpose:** Draw a folded, scaled joint site frequency spectrum with scikit-allel.
- **Method:** Call `allel.plot_joint_sfs_folded_scaled()` on the matrix and save the figure as a PNG. scikit-allel does not rescale the values: they are drawn as given, on a linear colour scale; pass the output of scikit-allel's `joint_sfs_folded_scaled`. The drawn matrix is the input transposed, so population 1 is on the x axis.
- **API:** `plot_joint_sfs_folded_scaled(spectrum: 'list[list[float]]', output_path: 'str', *, imshow_kwargs: 'dict[str, Any] | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages scikit-allel (extra `dnakit[bio-scikit-allel]`) and matplotlib (extra `dnakit[bio-plots]`). Works offline.
- **Input:** `spectrum` is a non-empty 2-D list of non-negative numbers: rows are the derived allele count in population 1 and columns the count in population 2. `output_path` must be a new `.png` path. `imshow_kwargs` is passed to matplotlib's `imshow`.
- **Output:** `path` and `data` with `spectrum`, `shape`, `drawn_matrix` (the transposed matrix as drawn) and the axis labels `xlabel` and `ylabel` (minor allele count of each population).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.population import plot_joint_sfs_folded_scaled

out = Path(tempfile.mkdtemp())
spectrum = [[0, 3, 1], [4, 2, 1], [1, 1, 0]]
result = plot_joint_sfs_folded_scaled(spectrum, str(out / "joint.png"), verbose=False)
data = result["data"]
print(data["shape"], data["xlabel"], "|", data["ylabel"])
print(data["drawn_matrix"])
print(result["path"].endswith("joint.png"))
```

- **Example output:**

```text
[3, 3] minor allele count (population 1) | minor allele count (population 2)
[[0.0, 4.0, 1.0], [3.0, 2.0, 1.0], [1.0, 1.0, 0.0]]
True
```

## 6) `dnakit.population.plot_joint_sfs_scaled()`

- **Purpose:** Draw a scaled joint site frequency spectrum with scikit-allel's `plot_joint_sfs_scaled`.
- **Method:** Call `allel.plot_joint_sfs_scaled()` on the matrix and save the figure as a PNG. scikit-allel does not rescale the values: they are drawn as given, on a linear colour scale, with the labels of the unscaled plot; pass the output of scikit-allel's `joint_sfs_scaled`. The drawn matrix is the input transposed, so population 1 is on the x axis.
- **API:** `plot_joint_sfs_scaled(spectrum: 'list[list[float]]', output_path: 'str', *, imshow_kwargs: 'dict[str, Any] | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages scikit-allel (extra `dnakit[bio-scikit-allel]`) and matplotlib (extra `dnakit[bio-plots]`). Works offline.
- **Input:** `spectrum` is a non-empty 2-D list of non-negative numbers: rows are the derived allele count in population 1 and columns the count in population 2. `output_path` must be a new `.png` path. `imshow_kwargs` is passed to matplotlib's `imshow`.
- **Output:** `path` and `data` with `spectrum`, `shape`, `drawn_matrix` (the transposed matrix as drawn) and the axis labels `xlabel` and `ylabel` (derived allele count of each population).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.population import plot_joint_sfs_scaled

out = Path(tempfile.mkdtemp())
spectrum = [[0, 3, 1], [4, 2, 1], [1, 1, 0]]
result = plot_joint_sfs_scaled(spectrum, str(out / "joint.png"), verbose=False)
data = result["data"]
print(data["shape"], data["xlabel"], "|", data["ylabel"])
print(data["drawn_matrix"])
print(result["path"].endswith("joint.png"))
```

- **Example output:**

```text
[3, 3] derived allele count (population 1) | derived allele count (population 2)
[[0.0, 4.0, 1.0], [3.0, 2.0, 1.0], [1.0, 1.0, 0.0]]
True
```

## 7) `dnakit.population.plot_moving_haplotype_frequencies()`

- **Purpose:** Draw how haplotype frequencies change along consecutive windows of variants.
- **Method:** Wrap the matrix in `allel.HaplotypeArray` and call `allel.plot_moving_haplotype_frequencies()`. The first and last position of each window are computed with `allel.moving_statistic()`, and the rank matrix is read back from the drawn mesh. Each plot column is a window of `window_size` variants and each row a frequency rank of the haplotypes in it (rank 1 is the most frequent non-singleton, 0 is a singleton).
- **API:** `plot_moving_haplotype_frequencies(positions: 'list[int]', haplotypes: 'list[list[int]]', window_size: 'int', output_path: 'str', *, start: 'int' = 0, stop: 'int | None' = None, n: 'int | None' = None, palette: 'str' = 'Paired', singleton_color: 'str' = 'w', verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package scikit-allel (extra `dnakit[bio-scikit-allel]`, which also installs scipy and seaborn) and matplotlib (extra `dnakit[bio-plots]`). Works offline.
- **Input:** `positions` is a list of non-negative integers, one per row of `haplotypes`. `haplotypes` is a variants x haplotypes matrix of non-negative whole numbers. `window_size` is an integer of at least 1. `output_path` must be a new `.png` path. `start` and `stop` are variant indices, `n` (at least 1) keeps only the n most frequent ranks. It also needs seaborn.
- **Output:** `path` and `data` with `window_size`, `start`, `stop`, `n`, `window_bounds` (first and last position of each window), `rank_matrix` (windows x haplotypes, as drawn), `palette` and `singleton_color`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.population import plot_moving_haplotype_frequencies

out = Path(tempfile.mkdtemp())
positions = [10, 25, 40, 70, 90, 120, 150, 170]
haplotypes = [
    [0, 0, 1, 1, 0, 1],
    [0, 0, 1, 1, 1, 0],
    [0, 0, 1, 1, 0, 0],
    [0, 0, 0, 0, 1, 1],
    [0, 0, 0, 0, 1, 1],
    [0, 0, 1, 1, 0, 0],
    [0, 0, 1, 1, 0, 0],
    [0, 0, 1, 1, 1, 1],
]
result = plot_moving_haplotype_frequencies(
    positions, haplotypes, 4, str(out / "moving.png"), verbose=False
)
data = result["data"]
print(data["window_bounds"])
print(data["rank_matrix"])
top = plot_moving_haplotype_frequencies(
    positions, haplotypes, 4, str(out / "top1.png"), n=1, verbose=False
)
print(top["data"]["rank_matrix"])
```

- **Example output:**

```text
[[10, 70], [90, 170]]
[[1, 1, 2, 2, 0, 0], [1, 1, 2, 2, 3, 3]]
[[1, 1, 0, 0, 0, 0], [1, 1, 0, 0, 0, 0]]
```

## 8) `dnakit.population.plot_pairwise_distance()`

- **Purpose:** Draw a pairwise distance matrix between samples as an image with optional sample labels.
- **Method:** Check that the input is a square matrix or a condensed vector of n*(n-1)/2 values, then call `allel.plot_pairwise_distance()`, which turns a condensed vector into a square matrix with a zero diagonal. The matrix is read back from the drawn image.
- **API:** `plot_pairwise_distance(distances: 'list[float]', output_path: 'str', *, labels: 'list[str] | None' = None, colorbar: 'bool' = True, imshow_kwargs: 'dict[str, Any] | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages scikit-allel (extra `dnakit[bio-scikit-allel]`) and matplotlib (extra `dnakit[bio-plots]`). Works offline.
- **Input:** `distances` is a list of finite numbers: either the condensed upper triangle, row by row, or a square n x n list of lists (any other length is a `ConfigurationError`). `labels` is a list of n strings. `colorbar` adds a colour bar. `output_path` must be a new `.png` path.
- **Output:** `path` and `data` with `n_samples`, `labels`, `colorbar` and `square_matrix` (the matrix as drawn).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.population import plot_pairwise_distance

out = Path(tempfile.mkdtemp())
condensed = [0.1, 0.2, 0.3]  # (a,b), (a,c), (b,c)
result = plot_pairwise_distance(
    condensed, str(out / "distance.png"), labels=["a", "b", "c"], verbose=False
)
data = result["data"]
print(data["n_samples"], data["labels"])
print(data["square_matrix"])
```

- **Example output:**

```text
3 ['a', 'b', 'c']
[[0.0, 0.1, 0.2], [0.1, 0.0, 0.3], [0.2, 0.3, 0.0]]
```

## 9) `dnakit.population.plot_pairwise_ld()`

- **Purpose:** Draw linkage disequilibrium (r squared) between all pairs of variants as a triangle image.
- **Method:** Check that the input is a square matrix or a condensed vector of n*(n-1)/2 values and call `allel.plot_pairwise_ld()`. scikit-allel keeps the lower triangle of the square matrix and flips it top to bottom; the matrix is read back from the drawn image. The function only draws values you pass in; it does not compute LD.
- **API:** `plot_pairwise_ld(ld_values: 'list[float]', output_path: 'str', *, colorbar: 'bool' = True, imshow_kwargs: 'dict[str, Any] | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages scikit-allel (extra `dnakit[bio-scikit-allel]`) and matplotlib (extra `dnakit[bio-plots]`). Works offline.
- **Input:** `ld_values` is a list of finite numbers: the condensed upper triangle, row by row, or a square matrix. `colorbar` adds a colour bar and `imshow_kwargs` goes to matplotlib's `imshow`. `output_path` must be a new `.png` path.
- **Output:** `path` and `data` with `n_variants`, `ld_values` (the input), `colorbar` and `drawn_matrix` (the flipped lower triangle as drawn).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.population import plot_pairwise_ld

out = Path(tempfile.mkdtemp())
condensed = [0.5, 0.1, 0.9]  # r squared of pairs (0,1), (0,2), (1,2)
result = plot_pairwise_ld(condensed, str(out / "ld.png"), verbose=False)
data = result["data"]
print(data["n_variants"])
print(data["drawn_matrix"])
```

- **Example output:**

```text
3
[[0.1, 0.9, 0.0], [0.5, 0.0, 0.0], [0.0, 0.0, 0.0]]
```

## 10) `dnakit.population.plot_sfs()`

- **Purpose:** Draw a one-dimensional site frequency spectrum.
- **Method:** Call `allel.plot_sfs()` on the spectrum and save the figure as a PNG; the points are read from the drawn line. 
- **API:** `plot_sfs(spectrum: 'list[float]', output_path: 'str', *, yscale: 'str' = 'log', bins: 'int | list[float] | None' = None, n: 'int | None' = None, clip_endpoints: 'bool' = True, label: 'str | None' = None, plot_kwargs: 'dict[str, Any] | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages scikit-allel (extra `dnakit[bio-scikit-allel]`) and matplotlib (extra `dnakit[bio-plots]`). Works offline.
- **Input:** `spectrum` is a non-empty 1-D list of non-negative site counts, indexed by allele count. `output_path` must be a new `.png` path. `yscale` is a matplotlib scale name (default `log`). `bins` is a whole number or a list of bin edges. `n` is the number of chromosomes sampled; when set, the x axis becomes allele frequency. `clip_endpoints` (default true) drops the first and last entries. `label` and `plot_kwargs` go to the line.
- **Output:** `path` and `data` with `spectrum`, the drawn points `x` and `y`, `xlabel` (`derived allele count`, or `derived allele frequency` when `n` is given), `ylabel` (`site frequency`), `yscale`, `n`, `bins`, `clip_endpoints` and `label`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.population import plot_sfs

out = Path(tempfile.mkdtemp())
spectrum = [0, 40, 20, 10, 5, 2, 0]
result = plot_sfs(spectrum, str(out / "sfs.png"), verbose=False)
data = result["data"]
print(data["x"], data["y"])
print(data["xlabel"], "|", data["ylabel"], "|", data["yscale"])
with_n = plot_sfs(spectrum, str(out / "sfs_n.png"), n=6, yscale="linear", verbose=False)
print(with_n["data"]["x"][:2], with_n["data"]["xlabel"])
```

- **Example output:**

```text
[1, 2, 3, 4, 5] [40.0, 20.0, 10.0, 5.0, 2.0]
derived allele count | site frequency | log
[0.16666666666666666, 0.3333333333333333] derived allele frequency
```

## 11) `dnakit.population.plot_sfs_folded()`

- **Purpose:** Draw a folded one-dimensional site frequency spectrum (minor allele counts).
- **Method:** Call `allel.plot_sfs_folded()` on the spectrum and save the figure as a PNG; the points are read from the drawn line. It takes the same arguments as `plot_sfs`; only the x label differs.
- **API:** `plot_sfs_folded(spectrum: 'list[float]', output_path: 'str', *, yscale: 'str' = 'log', bins: 'int | list[float] | None' = None, n: 'int | None' = None, clip_endpoints: 'bool' = True, label: 'str | None' = None, plot_kwargs: 'dict[str, Any] | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages scikit-allel (extra `dnakit[bio-scikit-allel]`) and matplotlib (extra `dnakit[bio-plots]`). Works offline.
- **Input:** `spectrum` is a non-empty 1-D list of non-negative site counts, indexed by allele count. `output_path` must be a new `.png` path. `yscale` is a matplotlib scale name (default `log`). `bins` is a whole number or a list of bin edges. `n` is the number of chromosomes sampled; when set, the x axis becomes allele frequency. `clip_endpoints` (default true) drops the first and last entries. `label` and `plot_kwargs` go to the line.
- **Output:** `path` and `data` with `spectrum`, the drawn points `x` and `y`, `xlabel` (`minor allele count`), `ylabel` (`site frequency`), `yscale`, `n`, `bins`, `clip_endpoints` and `label`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.population import plot_sfs_folded

out = Path(tempfile.mkdtemp())
spectrum = [0, 40, 20, 10, 5, 2, 0]
result = plot_sfs_folded(spectrum, str(out / "sfs.png"), verbose=False)
data = result["data"]
print(data["x"], data["y"])
print(data["xlabel"], "|", data["ylabel"], "|", data["yscale"])
```

- **Example output:**

```text
[1, 2, 3, 4, 5] [40.0, 20.0, 10.0, 5.0, 2.0]
minor allele count | site frequency | log
```

## 12) `dnakit.population.plot_sfs_folded_scaled()`

- **Purpose:** Draw a folded, scaled one-dimensional site frequency spectrum.
- **Method:** Call `allel.plot_sfs_folded_scaled()` on the spectrum and save the figure as a PNG; the points are read from the drawn line. scikit-allel draws the values as given and only changes the labels; pass the output of its `sfs_folded_scaled` to plot scaled values.
- **API:** `plot_sfs_folded_scaled(spectrum: 'list[float]', output_path: 'str', *, yscale: 'str' = 'log', bins: 'int | list[float] | None' = None, n: 'int | None' = None, clip_endpoints: 'bool' = True, label: 'str | None' = None, plot_kwargs: 'dict[str, Any] | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages scikit-allel (extra `dnakit[bio-scikit-allel]`) and matplotlib (extra `dnakit[bio-plots]`). Works offline.
- **Input:** `spectrum` is a non-empty 1-D list of non-negative site counts, indexed by allele count. `output_path` must be a new `.png` path. `yscale` is a matplotlib scale name (default `log`). `bins` is a whole number or a list of bin edges. `n` is the number of chromosomes sampled; when set, the x axis becomes allele frequency. `clip_endpoints` (default true) drops the first and last entries. `label` and `plot_kwargs` go to the line.
- **Output:** `path` and `data` with `spectrum`, the drawn points `x` and `y`, `xlabel` (`minor allele count`), `ylabel` (`scaled site frequency`), `yscale`, `n`, `bins`, `clip_endpoints` and `label`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.population import plot_sfs_folded_scaled

out = Path(tempfile.mkdtemp())
spectrum = [0, 40, 20, 10, 5, 2, 0]
result = plot_sfs_folded_scaled(spectrum, str(out / "sfs.png"), verbose=False)
data = result["data"]
print(data["x"], data["y"])
print(data["xlabel"], "|", data["ylabel"], "|", data["yscale"])
```

- **Example output:**

```text
[1, 2, 3, 4, 5] [40.0, 20.0, 10.0, 5.0, 2.0]
minor allele count | scaled site frequency | log
```

## 13) `dnakit.population.plot_sfs_scaled()`

- **Purpose:** Draw a scaled one-dimensional site frequency spectrum.
- **Method:** Call `allel.plot_sfs_scaled()` on the spectrum and save the figure as a PNG; the points are read from the drawn line. scikit-allel draws the values as given and only changes the y label; pass the output of its `sfs_scaled` to plot scaled values.
- **API:** `plot_sfs_scaled(spectrum: 'list[float]', output_path: 'str', *, yscale: 'str' = 'log', bins: 'int | list[float] | None' = None, n: 'int | None' = None, clip_endpoints: 'bool' = True, label: 'str | None' = None, plot_kwargs: 'dict[str, Any] | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages scikit-allel (extra `dnakit[bio-scikit-allel]`) and matplotlib (extra `dnakit[bio-plots]`). Works offline.
- **Input:** `spectrum` is a non-empty 1-D list of non-negative site counts, indexed by allele count. `output_path` must be a new `.png` path. `yscale` is a matplotlib scale name (default `log`). `bins` is a whole number or a list of bin edges. `n` is the number of chromosomes sampled; when set, the x axis becomes allele frequency. `clip_endpoints` (default true) drops the first and last entries. `label` and `plot_kwargs` go to the line.
- **Output:** `path` and `data` with `spectrum`, the drawn points `x` and `y`, `xlabel` (`derived allele count`), `ylabel` (`scaled site frequency`), `yscale`, `n`, `bins`, `clip_endpoints` and `label`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.population import plot_sfs_scaled

out = Path(tempfile.mkdtemp())
spectrum = [0, 40, 20, 10, 5, 2, 0]
result = plot_sfs_scaled(spectrum, str(out / "sfs.png"), verbose=False)
data = result["data"]
print(data["x"], data["y"])
print(data["xlabel"], "|", data["ylabel"], "|", data["yscale"])
```

- **Example output:**

```text
[1, 2, 3, 4, 5] [40.0, 20.0, 10.0, 5.0, 2.0]
derived allele count | scaled site frequency | log
```

## 14) `dnakit.population.plot_variant_locator()`

- **Purpose:** Draw lines from the genome position of each variant to its place in the variant index.
- **Method:** Check that the positions are sorted and call `allel.plot_variant_locator()`. When `step` is not given, DNAKit uses `max(1, n_variants // 100)` (scikit-allel itself would fail with a zero step for fewer than 100 variants). The x limits and the number of drawn lines are read from the axes.
- **API:** `plot_variant_locator(positions: 'list[int]', output_path: 'str', *, step: 'int | None' = None, start: 'int | None' = None, stop: 'int | None' = None, flip: 'bool' = False, line_kwargs: 'dict[str, Any] | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages scikit-allel (extra `dnakit[bio-scikit-allel]`) and matplotlib (extra `dnakit[bio-plots]`). Works offline.
- **Input:** `positions` is a list of non-negative integers for one contig, in ascending order (unsorted input is a `ConfigurationError`). `step` is an integer of at least 1. `start` and `stop` bound the genome region. `flip` flips the plot and `line_kwargs` goes to the lines. `output_path` must be a new `.png` path.
- **Output:** `path` and `data` with `n_variants`, `step`, `start` and `stop` (the x limits as drawn), `flip` and `lines_drawn`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.population import plot_variant_locator

out = Path(tempfile.mkdtemp())
positions = [10, 25, 40, 70, 90, 120, 150, 170]
result = plot_variant_locator(
    positions, str(out / "locator.png"), step=2, start=0, stop=200, verbose=False
)
data = result["data"]
print(data["n_variants"], data["step"], data["start"], data["stop"], data["lines_drawn"])
default = plot_variant_locator(positions, str(out / "locator2.png"), verbose=False)
print(default["data"]["step"], default["data"]["lines_drawn"])
```

- **Example output:**

```text
8 2 0.0 200.0 4
1 8
```

## 15) `dnakit.population.plot_voight_painting()`

- **Purpose:** Draw one Voight painting, in which each shared haplotype prefix is a coloured block.
- **Method:** Validate the painting array and call `allel.plot_voight_painting()`; `flank='left'` reverses the variant axis. The number of colours is the largest group id in the array.
- **API:** `plot_voight_painting(painting: 'list[list[int]]', output_path: 'str', *, palette: 'str' = 'colorblind', flank: 'str' = 'right', height_factor: 'float' = 0.01, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package scikit-allel (extra `dnakit[bio-scikit-allel]`, which also installs scipy and seaborn) and matplotlib (extra `dnakit[bio-plots]`). Works offline.
- **Input:** `painting` is a non-empty 2-D list of non-negative whole numbers, normally the first result of scikit-allel's `voight_painting` (variants x haplotypes; 0 marks singletons). `flank` is `right` or `left`, anything else is a `ConfigurationError`. `palette` and `height_factor` are as in `fig_voight_painting`. `output_path` must be a new `.png` path. It also needs seaborn.
- **Output:** `path` and `data` with `painting`, `flank`, `n_colors`, `palette` and `height_factor`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.population import plot_voight_painting

out = Path(tempfile.mkdtemp())
painting = [
    [1, 1, 2, 2, 0],
    [1, 1, 2, 3, 0],
    [1, 1, 0, 0, 0],
]
result = plot_voight_painting(painting, str(out / "painting.png"), flank="left", verbose=False)
data = result["data"]
print(data["flank"], data["n_colors"])
try:
    plot_voight_painting(painting, str(out / "bad.png"), flank="up", verbose=False)
except Exception as error:
    print(type(error).__name__, error.message)
```

- **Example output:**

```text
left 3
ConfigurationError flank must be 'right' or 'left'.
```

## 16) `dnakit.population.roh_poissonhmm()`

- **Purpose:** Call runs of homozygosity (ROH) in one individual with scikit-allel's Poisson hidden Markov model.
- **Method:** Validate the genotypes and positions, check the window and contig settings, require the `protopunica` package, then call `allel.roh_poissonhmm()` with the heterozygote rate inside a ROH (`phet_roh`) and the two rates outside it (`phet_nonroh_low`, `phet_nonroh_high`). The ROH table is returned as records together with FROH, the fraction of the contig covered by ROH.
- **API:** `roh_poissonhmm(genotypes: 'list[list[int]]', positions: 'list[int]', *, phet_roh: 'float' = 0.001, phet_nonroh_low: 'float' = 0.0025, phet_nonroh_high: 'float' = 0.01, transition: 'float' = 0.001, window_size: 'int' = 1000, min_roh: 'int' = 0, is_accessible: 'list[bool] | None' = None, contig_size: 'int | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package scikit-allel (extra `dnakit[bio-scikit-allel]`) and the package `protopunica`. `protopunica` is not part of any DNAKit extra and must be installed separately with pip; without it a `BackendUnavailableError` with an install hint is raised. It is not installed in the environment used for this page, so the example below shows the argument check and the real missing-dependency behaviour.
- **Input:** `genotypes` is a variants x ploidy list of non-negative allele indices (two columns for a diploid) and `positions` the matching 1-based positions. Give either `is_accessible` (one true/false per base of the contig) or `contig_size`, otherwise a `ConfigurationError` is raised. `window_size` is at least 1, `min_roh` at least 0.
- **Output:** A dictionary with `status`, `roh_blocks` (a list of records from scikit-allel's ROH table), `block_count`, `froh`, `window_size`, `phet_roh`, `phet_nonroh`, `transition` and `min_roh`.
- **Example code:**

```python
from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.population import roh_poissonhmm

genotypes = [[0, 0], [0, 1], [1, 1]]
positions = [100, 200, 300]
try:
    roh_poissonhmm(genotypes, positions, verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error.message)
try:
    roh_poissonhmm(genotypes, positions, contig_size=1000, verbose=False)
except BackendUnavailableError as error:
    print(type(error).__name__, "protopunica could not be imported" in error.message)
```

- **Example output:**

```text
ConfigurationError Provide is_accessible (one flag per base) or contig_size.
BackendUnavailableError True
```

Usage with the dependency installed:

```python
from dnakit.population import roh_poissonhmm

# genotypes: one row per variant, two allele columns for a diploid individual
result = roh_poissonhmm(
    genotypes,
    positions,
    contig_size=5_000_000,
    window_size=1000,
    verbose=False,
)
print(result["block_count"], result["froh"])
```

## 17) `dnakit.population.allel_vcf_to_zarr()`

- **Purpose:** Read a VCF file with scikit-allel and write its variant and genotype arrays to a new Zarr store.
- **Method:** Read the file with `allel.read_vcf()` and write each array with zarr 3's `create_array()` into the `vcf_to_zarr` layout (`variants/<FIELD>`, `calldata/<FIELD>`, `samples`). It does not call scikit-allel's own `vcf_to_zarr`, which writes with the zarr 2 API that zarr 3 lacks. String arrays are stored as strings. If the VCF has no variants, no store is written.
- **API:** `allel_vcf_to_zarr(vcf_path: 'str | Path', zarr_path: 'str | Path', *, fields: 'Sequence[str] | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages scikit-allel (extra `dnakit[bio-scikit-allel]`) and zarr 3 (extra `dnakit[io-zarr]`). Works offline. Note that if zarr is missing the error hint names the extra `bio-scikit-allel`, which does not contain zarr.
- **Input:** `vcf_path` must be an existing file. `zarr_path` must not exist yet and its parent directory must exist, otherwise a `ConfigurationError` is raised. `fields` is a non-empty list of up to 200 scikit-allel field names such as `variants/CHROM` or `calldata/GT`, or `None` for scikit-allel's default fields. The whole VCF is read into memory.
- **Output:** A dictionary with `status`, `operation`, `output_path` (absolute path of the store, or `None` for a VCF without variants), `variants`, `samples`, `arrays` (name, shape and dtype of each array) and `provenance` (package versions and layout). A VCF without variants returns `variants` 0, an empty `arrays` list and no `samples` key.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.population import allel_vcf_to_zarr

out = Path(tempfile.mkdtemp())
result = allel_vcf_to_zarr("mini.vcf", out / "mini.zarr", verbose=False)
print(result["status"], result["variants"], result["samples"])
for array in result["arrays"]:
    print(array["name"], array["shape"])
only_gt = allel_vcf_to_zarr(
    "mini.vcf", out / "gt.zarr", fields=["calldata/GT"], verbose=False
)
print([array["name"] for array in only_gt["arrays"]])
try:
    allel_vcf_to_zarr("mini.vcf", out / "mini.zarr", verbose=False)
except Exception as error:
    print(type(error).__name__, "already exists" in error.message)
```

- **Example output:**

```text
completed 3 3
calldata/GT [3, 3, 2]
samples [3]
variants/ALT [3, 3]
variants/CHROM [3]
variants/FILTER_PASS [3]
variants/ID [3]
variants/POS [3]
variants/QUAL [3]
variants/REF [3]
['calldata/GT']
ConfigurationError True
```
