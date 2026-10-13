# Optional Hail backends: MatrixTable import, QC and population analysis

This group of functions in `dnakit.variants` is a set of optional adapters over [hail](https://hail.is/). They import a VCF into a Hail MatrixTable, compute site and sample QC metrics, the exact Hardy-Weinberg test, linkage-disequilibrium pruning, genotype PCA, and intercept-only linear or logistic association tests. Every function works only on an explicit VCF, MatrixTable directory, or phenotype file. None of them accesses the network, calls a variant caller, or reads BGEN, PLINK, or gVCF input.

The dependencies are optional and are imported only when a function is called. The first call starts a local Spark session, which needs Java. The examples on this page were run with Java 17; Hail prints a Java 11 compatibility notice, which goes to stderr.

## Example data

All examples below use these two files. The VCF has 4 samples, 3 biallelic sites, and 1 multi-allelic site; the multi-allelic site is skipped and counted by analyses that require biallelic sites. The data only demonstrates the interface, and its statistics do not describe a real population.

**Example file:** `sample.vcf`

```text
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4
chr1	10	rs1	A	G	50	PASS	.	GT	0/1	1/1	0/0	./.
chr1	20	rs2	C	T	.	PASS	.	GT	0|1	0/0	0/1	1/1
chr1	30	rs3	G	C	99	PASS	.	GT	0/0	0/1	1/1	0/1
chr2	5	rs4	A	C,T	60	PASS	.	GT	0/1	1/2	./.	0/0
```

**Example file:** `phenotype.tsv`

```text
sample	y	case
S1	1.0	0
S2	2.5	1
S3	3.0	0
S4	2.0	1
```

## 1) `dnakit.variants.import_vcf_matrix_table()`

- **Purpose:** Import a VCF into a Hail MatrixTable and write it to a new directory.
- **Method:** Call `hl.import_vcf` with the reference genome named by `reference_genome` (`GRCh37` or `GRCh38`). The output path must not exist, and the write never overwrites existing content.
- **API:** `dnakit.variants.import_vcf_matrix_table(source, output_path, *, reference_genome="GRCh38")`.
- **Input:** A VCF path (plain or bgzip-compressed); `output_path` is a path that does not exist yet.
- **Output:** `n_variants`, `n_samples`, the written `path`, and `provenance`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants import import_vcf_matrix_table

target = Path(tempfile.mkdtemp()) / "sample.mt"
result = import_vcf_matrix_table("sample.vcf", target)
print(result["n_variants"], result["n_samples"])
```

- **Example output:**

```text
4 4
```

## 2) `dnakit.variants.summarize_vcf_matrix_table()`

- **Purpose:** Summarize the variant count, biallelic and SNV counts, contig distribution, and sample names of a VCF or MatrixTable directory.
- **Method:** Aggregate with `hl.agg.count_where` and `hl.agg.counter`. Contigs are sorted by name and sample names by string order. The call raises an error if the contig count exceeds `max_contigs` or the sample count exceeds `max_samples`.
- **API:** `dnakit.variants.summarize_vcf_matrix_table(source, *, reference_genome="GRCh38", max_contigs=100, max_samples=100000)`.
- **Input:** A VCF path, or a MatrixTable directory written by `import_vcf_matrix_table()`.
- **Output:** `n_variants`, `n_samples`, `biallelic_variants`, `multiallelic_variants`, `snv_variants`, `contigs`, and `samples`.
- **Example code:**

```python
from dnakit.variants import summarize_vcf_matrix_table

result = summarize_vcf_matrix_table("sample.vcf")
print(result["biallelic_variants"], result["multiallelic_variants"], result["contigs"])
```

- **Example output:**

```text
3 1 [{'contig': 'chr1', 'variants': 3}, {'contig': 'chr2', 'variants': 1}]
```

## 3) `dnakit.variants.calculate_variant_qc()`

- **Purpose:** Compute the call rate, allele frequency, and heterozygote and homozygote counts for each biallelic site.
- **Method:** Remove multi-allelic sites, then call `hl.variant_qc`. The allele frequency is `AC / AN`, and missing genotypes are not counted in `AN`. The Hardy-Weinberg test is computed separately by `calculate_hardy_weinberg_pvalues()`.
- **API:** `dnakit.variants.calculate_variant_qc(source, *, reference_genome="GRCh38", max_variants=100000)`.
- **Input:** A VCF path or MatrixTable directory; `max_variants` is the upper limit on the number of biallelic sites.
- **Output:** A `variants` list (each item has `chrom`, `pos`, `ref`, `alt`, `n_called`, `call_rate`, `allele_count`, `allele_number`, `allele_frequency`, `n_het`, and `n_hom_alt`) and `skipped_multiallelic`.
- **Example code:**

```python
from dnakit.variants import calculate_variant_qc

result = calculate_variant_qc("sample.vcf")
first = result["variants"][0]
print(first["pos"], first["call_rate"], first["allele_frequency"], first["n_het"])
```

- **Example output:**

```text
10 0.75 0.5 1
```

## 4) `dnakit.variants.calculate_sample_qc()`

- **Purpose:** Compute the call rate, genotype class counts, and transition-to-transversion ratio for each sample.
- **Method:** Call `hl.sample_qc` over all sites. Transitions and transversions are counted per alternate-allele copy, so a homozygous alternate genotype counts twice. `r_ti_tv` is `None` when there are no SNVs.
- **API:** `dnakit.variants.calculate_sample_qc(source, *, reference_genome="GRCh38", max_samples=100000)`.
- **Input:** A VCF path or MatrixTable directory.
- **Output:** A `samples` list sorted by sample name (each item has `sample`, `call_rate`, `n_called`, `n_hom_ref`, `n_het`, `n_hom_var`, `r_het_hom_var`, `n_transition`, `n_transversion`, and `r_ti_tv`).
- **Example code:**

```python
from dnakit.variants import calculate_sample_qc

result = calculate_sample_qc("sample.vcf")
print([(row["sample"], row["call_rate"]) for row in result["samples"]])
```

- **Example output:**

```text
[('S1', 1.0), ('S2', 1.0), ('S3', 0.75), ('S4', 0.75)]
```

## 5) `dnakit.variants.calculate_hardy_weinberg_pvalues()`

- **Purpose:** Test each biallelic site for Hardy-Weinberg equilibrium.
- **Method:** Count the homozygous reference, heterozygous, and homozygous alternate genotypes, then call `hl.hardy_weinberg_test`. The default is a two-sided exact test with mid-p correction (Levene–Haldane distribution). With `one_sided=True`, the test looks for heterozygote excess.
- **API:** `dnakit.variants.calculate_hardy_weinberg_pvalues(source, *, reference_genome="GRCh38", one_sided=False, max_variants=100000)`.
- **Input:** A VCF path or MatrixTable directory; `one_sided` is a boolean.
- **Output:** A `variants` list (each item has the genotype counts, `het_freq_hwe`, and `p_value`), `alternative` (`two_sided` or `excess_heterozygosity`), and `skipped_multiallelic`.
- **Example code:**

```python
from dnakit.variants import calculate_hardy_weinberg_pvalues

result = calculate_hardy_weinberg_pvalues("sample.vcf")
print([(row["pos"], round(row["p_value"], 4)) for row in result["variants"]])
```

- **Example output:**

```text
[(10, 0.7), (20, 0.6571), (30, 0.6571)]
```

## 6) `dnakit.variants.prune_linkage_disequilibrium()`

- **Purpose:** Remove correlated sites within linkage-disequilibrium windows and return the retained set.
- **Method:** Call `hl.ld_prune` with diploid genotypes from biallelic sites, ignoring phase. A site is kept when its pairwise r² with every other kept site in the window is strictly below `r2`. The result is a maximal independent set, not necessarily a maximum one.
- **API:** `dnakit.variants.prune_linkage_disequilibrium(source, *, reference_genome="GRCh38", r2=0.2, bp_window_size=1000000, max_variants=100000)`.
- **Input:** A VCF path or MatrixTable directory; `r2` is in (0, 1]; `bp_window_size` is a positive integer in bp.
- **Output:** `n_input_variants`, `n_retained`, `skipped_multiallelic`, and `retained_variants`, a list sorted by coordinate.
- **Example code:**

```python
from dnakit.variants import prune_linkage_disequilibrium

result = prune_linkage_disequilibrium("sample.vcf", r2=0.2)
print(result["n_input_variants"], [row["pos"] for row in result["retained_variants"]])
```

- **Example output:**

```text
3 [10]
```

## 7) `dnakit.variants.calculate_hwe_normalized_pca()`

- **Purpose:** Run PCA on Hardy-Weinberg-normalized genotypes and return the eigenvalues and sample scores.
- **Method:** Call `hl.hwe_normalized_pca`. Hail removes sites that are all homozygous reference or all homozygous alternate. The sign of each score is set by the SVD, so compare absolute values or relative positions between runs rather than raw signs.
- **API:** `dnakit.variants.calculate_hwe_normalized_pca(source, *, reference_genome="GRCh38", n_components=10, max_variants=100000, max_samples=100000)`.
- **Input:** A VCF path or MatrixTable directory; `n_components` is a positive integer that cannot exceed the number of samples.
- **Output:** `eigenvalues`, `samples` (each item has `sample` and `scores`), `n_variants`, `n_samples`, and `skipped_multiallelic`.
- **Example code:**

```python
from dnakit.variants import calculate_hwe_normalized_pca

result = calculate_hwe_normalized_pca("sample.vcf", n_components=2)
print([round(value, 4) for value in result["eigenvalues"]])
first = result["samples"][0]
print(first["sample"], [round(value, 4) for value in first["scores"]])
```

- **Example output:**

```text
[2.2761, 1.3333]
S1 [-0.4082, 0.5774]
```

## 8) `dnakit.variants.calculate_association_tests()`

- **Purpose:** Test the association of one phenotype with allele dosage at each biallelic site. The model has an intercept only, with no covariates.
- **Method:** The phenotype table is read in Python and aligned to samples through `hl.literal`. A numeric phenotype uses linear regression (`hl.linear_regression_rows`, t test); a 0/1 phenotype uses logistic regression (`hl.logistic_regression_rows`, Wald test). Hail mean-imputes missing dosages. Samples without a phenotype are excluded.
- **API:** `dnakit.variants.calculate_association_tests(source, phenotype_path, *, phenotype_column, model="linear", sample_column="sample", reference_genome="GRCh38", max_samples=100000, max_variants=100000)`.
- **Input:** A VCF path or MatrixTable directory. The phenotype file is tab-separated text and must contain `sample_column` and `phenotype_column`. `model` is `"linear"` or `"logistic"`.
- **Output:** A `results` list (each item has `beta`, `standard_error`, `statistic`, and `p_value`; logistic results also have `converged`), `n_samples_with_phenotype`, `statistic` (`t` or `z`), `covariates` (always an empty list), and `skipped_multiallelic`.
- **Example code:**

```python
from dnakit.variants import calculate_association_tests

result = calculate_association_tests("sample.vcf", "phenotype.tsv", phenotype_column="y")
for row in result["results"]:
    print(row["pos"], row["beta"], round(row["p_value"], 4))
```

- **Example output:**

```text
10 -0.25 0.761
20 -0.25 0.761
30 1.0 0.0438
```

## How these functions are verified

The test file is `others/tests/unit/test_hail_backends.py`, with 14 tests. In an environment with hail installed, all of them run. Without hail, the 7 argument-validation and missing-package tests run and the other 7 are skipped with a reason.

- Each function that needs Hail has at least one test that runs on real Hail.
- Site and sample call rates, allele frequencies, and heterozygote and homozygote counts are checked against counts made genotype by genotype. Sample S2's transition and transversion counts are checked per alternate-allele copy.
- The Hardy-Weinberg exact test matches a mid-p value from the Levene–Haldane distribution computed with exact fractions, in both the two-sided and one-sided forms.
- The linear regression's beta, standard error, and p value match the closed-form formulas (four samples and two parameters, so the t distribution has 2 degrees of freedom).
- The logistic results for sites 20 and 30 match the hand values (beta 0, standard error sqrt(2), p 1). Site 10 has quasi-complete separation and no finite maximum-likelihood estimate, so it returns `converged: false` and `beta: None`.
- The LD pruning result has pairwise r² below the threshold among kept sites, every dropped site has r² at or above the threshold with a kept site (computed by hand), and the result matches a direct `hl.ld_prune` call.
- The PCA eigenvalues and sample scores match a direct `hl.hwe_normalized_pca` call.
- A missing hail package reports the `bio-hail` install hint. Invalid arguments (booleans, out-of-range values, unsupported reference genomes, duplicate phenotype samples, non-numeric phenotypes, non-0/1 logistic phenotypes, and existing output paths) are rejected before Spark starts.

## Install optional dependency

```bash
python -m pip install "dnakit[bio-hail]"   # hail==0.2.139; requires Java
```

Project link: [hail](https://hail.is/).

## Current scope

- Wrapped: VCF import and MatrixTable write, MatrixTable summary, site and sample QC, the Hardy-Weinberg exact test (two-sided and one-sided for heterozygote excess), LD pruning, HWE-normalized PCA, and intercept-only linear and logistic regression (Wald test).
- Biallelic sites only. Multi-allelic sites are skipped and counted in `skipped_multiallelic`; the summary also reports `multiallelic_variants`.
- Not wrapped: regression with covariates; relatedness and GRM tools (for example `pc_relate`); gVCF and VDS combination; BGEN and PLINK input; imputation workflows; and arbitrary Hail expression queries. These need separate parameter design and validation on large data, which this batch does not cover.
- `max_variants` and `max_samples` limit the size of results returned to Python and the number of phenotype rows; they do not limit Spark's computation. Hail still scans the complete VCF.
- Different from existing code: `calculate_hwe_normalized_pca()` uses Hail's HWE normalization. It is not `dnakit.genomics.genotype_pca()` (PCA on a dense dosage matrix) or `dnakit.variants.calculate_pca()` (scikit-allel with Patterson scaling).
- On the first call, Hail starts a local Spark session with 1 GB of driver memory, and it writes its log to `dnakit-hail.log` in the system temporary directory.
