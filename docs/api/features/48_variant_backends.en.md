# Optional VCF backends and population-genetics statistics

This group of functions in `dnakit.variants` is a set of adapters over optional packages. VCF reading can use `cyvcf2`, `pysam`, or `vcfpy`; the three return identical record structures. Population statistics use `scikit-allel` and read those same records. Every function works only on an explicit VCF file or record dictionary. None of them queries a reference database or calls a variant caller.

The dependencies are optional and are imported only when the matching function is called. If a package is missing, the function raises an error that includes the install command.

## Example data

All examples below use this VCF. It has 5 sites and 4 samples. It is only for demonstrating the interface, and its statistics do not describe a real population.

**Example file:** `sample.vcf`

```text
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="AF">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4
chr1	10	rs1	A	G	50	PASS	DP=20;AF=0.5;DB	GT	0/1	1/1	0/0	./.
chr1	20	.	C	T,A	.	q10	DP=7;AF=0.2,0.1	GT	0|1	0/0	0/2	1/2
chr1	30	rs3	G	C	99	.	DP=12	GT	0/0	0/1	1/1	0/1
chr1	40	rs4	T	A	60	PASS	DP=9;AF=0.4	GT	1/1	0/1	0/0	0/0
chr1	50	rs5	G	A	70	PASS	DP=11;AF=0.3	GT	0/1	0/1	./.	0/0
```

## 1) `dnakit.variants.read_vcf()`

- **Purpose:** Read a VCF through one backend and return normalized variant records.
- **Method:** Parse each site and normalize genotypes to `"0/1"`, `"0|1"`, and `"./."` forms. The filter field keeps only failing filter names, so `PASS` becomes an empty list. INFO fields with a single value are returned as scalars.
- **API:** `dnakit.variants.read_vcf(source, *, backend="cyvcf2", max_records=100000)`.
- **Input:** A VCF path. `backend` is `"cyvcf2"`, `"pysam"`, or `"vcfpy"`. `max_records` is an explicit record limit; exceeding it raises an error.
- **Output:** A dictionary with `samples`, `variants`, and `provenance`. Each variant has `chrom`, `pos`, `id`, `ref`, `alt`, `qual`, `filter`, `info`, and `genotypes`.
- **Example code:**

```python
from dnakit.variants import read_vcf

result = read_vcf("sample.vcf", backend="cyvcf2")
print(result["samples"])
print(result["variants"][0])
```

- **Example output:**

```text
['S1', 'S2', 'S3', 'S4']
{'chrom': 'chr1', 'pos': 10, 'id': 'rs1', 'ref': 'A', 'alt': ['G'], 'qual': 50.0, 'filter': [], 'info': {'DP': 20, 'AF': 0.5, 'DB': True}, 'genotypes': ['0/1', '1/1', '0/0', './.']}
```

## 2) `dnakit.variants.fetch_vcf_region()`

- **Purpose:** Read only the variants in one region, without scanning the whole file.
- **Method:** Call tabix-style region queries in `pysam` or `cyvcf2`. Region coordinates are 1-based and inclusive; they are converted internally to 0-based half-open intervals.
- **API:** `dnakit.variants.fetch_vcf_region(source, region, *, backend="pysam", max_records=100000)`.
- **Input:** A bgzip-compressed VCF with a tabix or CSI index. `region` looks like `"chr1:1000-2000"`, `"chr1:1000"`, or `"chr1"`.
- **Output:** The same structure as `read_vcf()`, with the requested region recorded in the `region` field.
- **Example code:**

```python
import pysam
from dnakit.variants import fetch_vcf_region

pysam.tabix_compress("sample.vcf", "sample.vcf.gz", force=True)
pysam.tabix_index("sample.vcf.gz", preset="vcf", force=True)
result = fetch_vcf_region("sample.vcf.gz", "chr1:15-35", backend="pysam")
print([variant["pos"] for variant in result["variants"]])
```

- **Example output:**

```text
[20, 30]
```

## 3) `dnakit.variants.calculate_observed_heterozygosity()`

- **Purpose:** Compute mean observed heterozygosity for each population.
- **Method:** Use scikit-allel's `heterozygosity_observed` to get the heterozygous fraction at each site, where the denominator is the number of called genotypes. Average across sites, excluding sites where every call is missing.
- **API:** `dnakit.variants.calculate_observed_heterozygosity(data, *, populations=None)`.
- **Input:** The output of `read_vcf()` or a dictionary with the same structure. `populations` is `{population_name: [sample_name, ...]}`. If omitted, all samples form one population.
- **Output:** For each population: sample count, number of usable sites, and mean heterozygosity.
- **Example code:**

```python
from dnakit.variants import read_vcf, calculate_observed_heterozygosity

data = read_vcf("sample.vcf", backend="cyvcf2")
result = calculate_observed_heterozygosity(
    data, populations={"A": ["S1", "S2"], "B": ["S3", "S4"]}
)
print(result["populations"])
```

- **Example output:**

```text
{'A': {'samples': 2, 'variants': 5, 'mean_observed_heterozygosity': 0.6}, 'B': {'samples': 2, 'variants': 5, 'mean_observed_heterozygosity': 0.3}}
```

## 4) `dnakit.variants.calculate_fst()`

- **Purpose:** Compute genome-wide Fst between two populations.
- **Method:** Use biallelic sites only. `method="hudson"` uses the Hudson estimator (numerator and denominator are summed separately, then divided). `method="weir_cockerham"` uses the Weir–Cockerham estimator. The standard error comes from a block jackknife with about ten blocks. Weir–Cockerham requires at least two samples in each population.
- **API:** `dnakit.variants.calculate_fst(data, *, populations, method="hudson")`.
- **Input:** A data dictionary. `populations` must contain exactly two populations. `method` is `"hudson"` or `"weir_cockerham"`.
- **Output:** Fst, standard error, number of sites, block size, and counts of skipped sites.
- **Example code:**

```python
from dnakit.variants import read_vcf, calculate_fst

data = read_vcf("sample.vcf", backend="cyvcf2")
result = calculate_fst(data, populations={"A": ["S1", "S2"], "B": ["S3", "S4"]})
print(result["fst"], result["standard_error"])
```

- **Example output:**

```text
0.4920634920634921 0.12077896535034918
```

## 5) `dnakit.variants.calculate_pca()`

- **Purpose:** Run principal component analysis on sample genotypes.
- **Method:** Use biallelic, diploid sites with complete calls. Remove monomorphic sites, then run scikit-allel's PCA with Patterson scaling.
- **API:** `dnakit.variants.calculate_pca(data, *, n_components=10, samples=None)`.
- **Input:** A data dictionary. `n_components` is a positive integer and cannot exceed the number of samples or polymorphic sites. `samples` optionally selects a subset.
- **Output:** Sample coordinates, explained-variance ratio for each component, and the number of sites used.
- **Example code:**

```python
from dnakit.variants import read_vcf, calculate_pca

data = read_vcf("sample.vcf", backend="cyvcf2")
result = calculate_pca(data, n_components=2)
print(result["explained_variance_ratio"], result["polymorphic_variants"])
```

- **Example output:**

```text
[0.9292512536048889, 0.07074874639511108] 2
```

## 6) `dnakit.variants.calculate_tajimas_d()`

- **Purpose:** Compute Tajima's D for one set of samples on one chromosome.
- **Method:** Use biallelic sites and call scikit-allel's `tajima_d`. When fewer than three sites segregate, scikit-allel returns NaN, and this function returns `None`.
- **API:** `dnakit.variants.calculate_tajimas_d(data, *, samples=None)`.
- **Input:** A data dictionary. All sites must be on the same chromosome. `samples` is optional.
- **Output:** Tajima's D, sample count, and site count.
- **Example code:**

```python
from dnakit.variants import read_vcf, calculate_tajimas_d

data = read_vcf("sample.vcf", backend="cyvcf2")
print(calculate_tajimas_d(data)["tajimas_d"])
```

- **Example output:**

```text
1.973246148487123
```

## 7) `dnakit.variants.calculate_watterson_theta()`

- **Purpose:** Compute Watterson's theta per base.
- **Method:** Computed by scikit-allel's `watterson_theta` as S / a_n / L. Here S is the number of segregating sites, a_n depends on the number of chromosomes, and L is `sequence_length`.
- **API:** `dnakit.variants.calculate_watterson_theta(data, *, sequence_length, samples=None)`.
- **Input:** A data dictionary. `sequence_length` is a positive integer and cannot be smaller than the largest site position. `samples` is optional.
- **Output:** Theta per base, sample count, and site count.
- **Example code:**

```python
from dnakit.variants import read_vcf, calculate_watterson_theta

data = read_vcf("sample.vcf", backend="cyvcf2")
print(calculate_watterson_theta(data, sequence_length=100)["theta_per_base"])
```

- **Example output:**

```text
0.015426997245179066
```

## 8) `dnakit.variants.calculate_nucleotide_diversity()`

- **Purpose:** Compute nucleotide diversity π per base.
- **Method:** For each site, compute ac × (an − ac) / C(an, 2), sum the values, and divide by `sequence_length`. The calculation is done by scikit-allel's `sequence_diversity`.
- **API:** `dnakit.variants.calculate_nucleotide_diversity(data, *, sequence_length, samples=None)`.
- **Input:** Same as Watterson's theta.
- **Output:** π per base, sample count, and site count.
- **Example code:**

```python
from dnakit.variants import read_vcf, calculate_nucleotide_diversity

data = read_vcf("sample.vcf", backend="cyvcf2")
print(calculate_nucleotide_diversity(data, sequence_length=100)["pi_per_base"])
```

- **Example output:**

```text
0.0224047619047619
```

## 9) `dnakit.variants.calculate_site_frequency_spectrum()`

- **Purpose:** Compute the unfolded site frequency spectrum (SFS).
- **Method:** Count the alternate-allele total at each biallelic site. With `n` diploid samples, the spectrum has length `2n + 1`.
- **API:** `dnakit.variants.calculate_site_frequency_spectrum(data, *, samples=None)`.
- **Input:** A data dictionary. `samples` is optional.
- **Output:** A `spectrum` list where item i is the number of sites with exactly i alternate alleles, plus a `chromosomes` field giving the total chromosome count.
- **Example code:**

```python
from dnakit.variants import read_vcf, calculate_site_frequency_spectrum

data = read_vcf("sample.vcf", backend="cyvcf2")
print(calculate_site_frequency_spectrum(data)["spectrum"])
```

- **Example output:**

```text
[0, 0, 1, 2, 1, 0, 0, 0, 0]
```

## 10) `dnakit.variants.calculate_windowed_ld()`

- **Purpose:** Compute the median linkage-disequilibrium r² in fixed-size windows.
- **Method:** Use scikit-allel's `windowed_r_squared` (Rogers and Huff method). Only biallelic, diploid sites with complete calls are used. Windows start at position 1 and advance by `window_size`. The last window ends at the largest retained site.
- **API:** `dnakit.variants.calculate_windowed_ld(data, *, window_size=10000, samples=None)`.
- **Input:** A data dictionary for one chromosome. `window_size` is in bp. `samples` is optional.
- **Output:** For each window: start and end coordinates (1-based, inclusive), number of sites, and median r².
- **Example code:**

```python
from dnakit.variants import read_vcf, calculate_windowed_ld

data = read_vcf("sample.vcf", backend="cyvcf2")
print(calculate_windowed_ld(data, window_size=100)["windows"])
```

- **Example output:**

```text
[{'start': 1, 'end': 40, 'variants': 2, 'median_r_squared': 0.7272727489471436}]
```

## Installing the optional dependencies

```bash
python -m pip install "dnakit[bio-cyvcf2]"          # cyvcf2: read_vcf / fetch_vcf_region
python -m pip install "dnakit[bio-reads]"           # pysam: read_vcf / fetch_vcf_region
python -m pip install "dnakit[bio-vcfpy]"           # vcfpy: read_vcf
python -m pip install "dnakit[bio-scikit-allel]"    # scikit-allel and scipy: calculate_* population statistics
```

Project pages: [cyvcf2](https://github.com/brentp/cyvcf2), [pysam](https://github.com/pysam-developers/pysam), [vcfpy](https://github.com/dahlo/vcfpy), [scikit-allel](https://github.com/cggh/scikit-allel).

## How these functions are verified

The test file is `others/tests/unit/test_variant_backends.py`, with 25 tests. They check the following:

- The three VCF backends return identical output, including filters, missing genotypes, and phasing.
- Watterson's theta, π, heterozygosity, and the SFS match values calculated by hand from the formulas.
- Tajima's D, Hudson Fst, Weir–Cockerham Fst, and PCA match direct calls to scikit-allel.
- A missing optional package raises an install hint. Unknown backends, exceeding `max_records`, malformed genotypes, and unknown sample names are all rejected.

## Current scope

This batch covers VCF reading and scikit-allel population statistics. The remaining DNA and gene functions will be added in later domain batches, with their documentation updated alongside them.
