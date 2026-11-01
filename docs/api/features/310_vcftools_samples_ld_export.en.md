# VCFtools: individual statistics, linkage disequilibrium, PCA, export and comparison

This group of functions in `dnakit.variants` continues the typed `vcftools` facades of the previous document (`309_vcftools_filters_site_reports`; the common behaviour, result dictionary, `filters` mapping and `output_prefix` rules are described there). They compute per-individual heterozygosity and inbreeding, variant burden, long runs of homozygosity and relatedness; pairwise linkage disequilibrium in six forms; the number of unique haplotypes per BED region; principal components and SNP loadings; export genotypes to 0/1/2, IMPUTE, LDhat, LDhelmet, BEAGLE, PLINK PED/MAP and TPED/TFAM files; and compare two VCF files (individuals, genotype discordance, phase switch errors).

VCFtools is distributed under the LGPL-3.0; DNAKit runs `vcftools` as a separate process and copies none of its code. Install it yourself (`conda install -c bioconda vcftools`, `apt install vcftools` or a source build) and put it on PATH or pass `vcftools_path`; `vcftools_pca()` needs a build configured with `./configure --enable-pca` (LAPACK) and raises `BackendUnavailableError` with that hint otherwise.

Three upstream behaviours of VCFtools 0.1.16 are reproduced as they are. `--hapcount` treats the first line of the BED file as a header and skips it, and it reads freed memory when the last record of the VCF lies inside the last region of the last chromosome, which can zero the last row; `vcftools_haplotype_count()` adds `warnings` for both. The long-run, LDhat, LDhelmet and BEAGLE outputs require one chromosome.

The examples import from `dnakit.variants` and write to a fresh temporary directory.

## Example data

`toy.vcf` is the same six-sample VCF as in the previous document: 20 phased sites on `chr1` (14) and `chr2` (6), three indels, two `LowQual` sites, one missing genotype and the FORMAT fields GT, DP, GQ and PL. `win.bed` has a header line and three regions; `pos.txt` lists two sites.

**Example file:** `toy.vcf`

```text
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=LowQual,Description="Low quality">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled likelihoods">
##contig=<ID=chr1,length=20000>
##contig=<ID=chr2,length=20000>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4	S5	S6
chr1	100	v1	G	T	38	PASS	DP=91;AF=0.1	GT:DP:GQ:PL	1|0:9:45:30,0,30	0|1:31:60:30,0,30	0|1:9:30:30,0,30	0|1:20:99:30,0,30	1|0:8:12:30,0,30	0|0:14:99:0,30,300
chr1	800	v2	G	A	55	PASS	DP=120;AF=0.25	GT:DP:GQ:PL	1|0:25:99:30,0,30	1|1:21:60:300,30,0	0|0:17:45:0,30,300	0|0:21:60:0,30,300	0|0:11:45:0,30,300	0|0:25:99:0,30,300
chr1	1500	v3	A	G	67	PASS	DP=99;AF=0.5	GT:DP:GQ:PL	1|1:10:45:300,30,0	0|1:11:99:30,0,30	0|1:23:60:30,0,30	0|0:36:45:0,30,300	0|1:10:60:30,0,30	1|1:9:45:300,30,0
chr1	2200	v4	A	AAC	38	PASS	DP=166;AF=0.5	GT:DP:GQ:PL	1|0:33:30:30,0,30	1|1:23:12:300,30,0	1|0:32:30:30,0,30	1|0:28:30:30,0,30	1|1:30:30:300,30,0	1|1:20:30:300,30,0
chr1	2900	v5	A	G	94	PASS	DP=172;AF=0.25	GT:DP:GQ:PL	0|1:31:60:30,0,30	0|0:31:12:0,30,300	0|0:31:30:0,30,300	0|0:31:12:0,30,300	0|0:12:30:0,30,300	1|0:36:60:30,0,30
chr1	3600	v6	C	T	22	LowQual	DP=180;AF=0.1	GT:DP:GQ:PL	0|0:37:30:0,30,300	0|0:35:12:0,30,300	0|1:36:45:30,0,30	0|0:36:45:0,30,300	0|0:25:60:0,30,300	0|0:11:30:0,30,300
chr1	4300	v7	A	G	81	PASS	DP=137;AF=0.4	GT:DP:GQ:PL	0|0:27:30:0,30,300	1|1:20:99:300,30,0	0|0:18:60:0,30,300	0|0:21:45:0,30,300	0|0:31:12:0,30,300	0|0:20:12:0,30,300
chr1	5000	v8	G	A	81	PASS	DP=124;AF=0.25	GT:DP:GQ:PL	0|0:6:30:0,30,300	0|0:36:60:0,30,300	0|0:28:30:0,30,300	0|1:11:60:30,0,30	1|1:13:45:300,30,0	1|0:30:12:30,0,30
chr1	5700	v9	TGT	T	55	PASS	DP=93;AF=0.5	GT:DP:GQ:PL	1|0:15:30:30,0,30	0|1:14:60:30,0,30	1|1:7:30:300,30,0	0|0:6:30:0,30,300	1|0:12:12:30,0,30	0|0:39:45:0,30,300
chr1	6400	v10	C	T	38	PASS	DP=173;AF=0.25	GT:DP:GQ:PL	0|0:40:30:0,30,300	0|0:15:99:0,30,300	1|0:39:12:30,0,30	0|0:38:30:0,30,300	./.:.:.:.	0|0:34:30:0,30,300
chr1	7100	v11	C	T	22	LowQual	DP=177;AF=0.1	GT:DP:GQ:PL	0|0:34:60:0,30,300	0|0:26:99:0,30,300	0|0:38:99:0,30,300	0|0:38:60:0,30,300	0|0:18:99:0,30,300	0|0:23:30:0,30,300
chr1	7800	v12	G	GAC	94	PASS	DP=135;AF=0.25	GT:DP:GQ:PL	0|1:29:12:30,0,30	1|0:15:60:30,0,30	1|1:22:60:300,30,0	1|0:14:30:30,0,30	0|0:35:30:0,30,300	1|0:20:30:30,0,30
chr1	8500	v13	T	G	38	PASS	DP=93;AF=0.4	GT:DP:GQ:PL	0|1:12:45:30,0,30	1|1:11:30:300,30,0	1|0:22:60:30,0,30	0|1:23:45:30,0,30	0|1:8:60:30,0,30	0|1:17:30:30,0,30
chr1	9200	v14	T	A	55	PASS	DP=137;AF=0.4	GT:DP:GQ:PL	1|0:32:12:30,0,30	0|1:23:30:30,0,30	0|0:14:45:0,30,300	0|1:8:12:30,0,30	0|1:39:30:30,0,30	0|0:21:30:0,30,300
chr2	300	v15	G	C	94	PASS	DP=204;AF=0.4	GT:DP:GQ:PL	0|0:33:30:0,30,300	1|0:37:30:30,0,30	0|1:40:45:30,0,30	1|0:31:30:30,0,30	0|0:38:30:0,30,300	1|0:25:60:30,0,30
chr2	1200	v16	G	C	38	PASS	DP=141;AF=0.1	GT:DP:GQ:PL	1|0:23:99:30,0,30	0|0:34:45:0,30,300	1|0:6:30:30,0,30	0|0:22:12:0,30,300	0|0:29:45:0,30,300	0|0:27:30:0,30,300
chr2	2100	v17	G	A	38	PASS	DP=158;AF=0.4	GT:DP:GQ:PL	1|0:11:30:30,0,30	0|1:39:45:30,0,30	0|1:15:99:30,0,30	1|1:30:30:300,30,0	1|1:26:12:300,30,0	1|1:37:99:300,30,0
chr2	3000	v18	T	A	94	PASS	DP=133;AF=0.1	GT:DP:GQ:PL	0|0:7:60:0,30,300	0|0:40:12:0,30,300	0|0:21:99:0,30,300	1|0:37:99:30,0,30	0|0:22:12:0,30,300	0|0:6:99:0,30,300
chr2	3900	v19	A	T	94	PASS	DP=117;AF=0.4	GT:DP:GQ:PL	0|0:18:99:0,30,300	1|0:10:99:30,0,30	1|0:15:30:30,0,30	0|0:27:12:0,30,300	1|0:22:60:30,0,30	1|1:25:12:300,30,0
chr2	4800	v20	T	C	67	PASS	DP=154;AF=0.1	GT:DP:GQ:PL	0|0:10:30:0,30,300	0|0:38:12:0,30,300	0|0:34:99:0,30,300	0|0:23:12:0,30,300	0|0:30:30:0,30,300	0|1:19:99:30,0,30
```

**Example file:** `win.bed`

```text
#CHROM	START	END
chr1	0	5000
chr1	5000	12000
chr2	0	4000
```

**Example file:** `pos.txt`

```text
chr1	100
chr1	800
```

## 1) `dnakit.variants.vcftools_heterozygosity()`

- **Purpose:** Compute the homozygosity and the inbreeding coefficient F of each individual.
- **Method:** Run `vcftools --het`; the table is `<output_prefix>.het`. F is the method-of-moments estimate (observed minus expected homozygous sites over the number of sites minus the expectation), as in PLINK `--het`.
- **API:** `dnakit.variants.vcftools_heterozygosity(vcf_path: 'str', output_prefix: 'str', *, filters: 'Mapping[str, object] | None' = None, max_rows: 'int' = 1000, vcftools_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the `vcftools` program (VCFtools 0.1.16 or later; `conda install -c bioconda vcftools`). DNAKit looks for it at `vcftools_path`, then on PATH, then next to the running Python interpreter; if it is missing, `BackendUnavailableError` is raised before anything is written.
- **Input:** `vcf_path`, `output_prefix`, `filters`, `max_rows`.
- **Output:** The common result dictionary with `table` (columns `INDV, O(HOM), E(HOM), N_SITES, F`).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants import vcftools_heterozygosity

result = vcftools_heterozygosity("toy.vcf", str(Path(tempfile.mkdtemp()) / "het"))
print(result["table"]["columns"])
for row in result["table"]["rows"][:3]:
    print(row)
```

- **Example output:**

```text
['INDV', 'O(HOM)', 'E(HOM)', 'N_SITES', 'F']
['S1', 9, 11.5, 19, -0.33549]
['S2', 11, 11.5, 19, -0.06839]
['S3', 9, 11.5, 19, -0.33549]
```

## 2) `dnakit.variants.vcftools_individual_burden()`

- **Purpose:** Count the variant burden of each individual by genotype class or by allele-frequency class.
- **Method:** Run `vcftools --indv-burden` (`burden="genotype_classes"`, `.iburden`: homozygous-reference, heterozygous, homozygous-alternative and missing genotypes), `--indv-freq-burden` (`"frequency_classes"`, `.ifreqburden`: non-reference alleles counted by the allele count of their site, both alleles of a homozygous genotype counted) or `--indv-freq-burden2` (`"frequency_classes_hom_once"`, same file name, a homozygous genotype counted once).
- **API:** `dnakit.variants.vcftools_individual_burden(vcf_path: 'str', output_prefix: 'str', *, burden: 'str' = 'genotype_classes', filters: 'Mapping[str, object] | None' = None, max_rows: 'int' = 1000, vcftools_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the `vcftools` program (VCFtools 0.1.16 or later; `conda install -c bioconda vcftools`). DNAKit looks for it at `vcftools_path`, then on PATH, then next to the running Python interpreter; if it is missing, `BackendUnavailableError` is raised before anything is written.
- **Input:** `vcf_path`, `output_prefix`, `burden`, `filters`, `max_rows`.
- **Output:** The common result dictionary with `table` (`INDV, N_HOM_REF, N_HET, N_HOM_ALT, N_MISS`, or `INDV` and one column per allele count).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants import vcftools_individual_burden

out = Path(tempfile.mkdtemp())
by_genotype = vcftools_individual_burden("toy.vcf", str(out / "g"))
print(by_genotype["table"]["columns"])
print(by_genotype["table"]["rows"][:2])
both = vcftools_individual_burden("toy.vcf", str(out / "a"), burden="frequency_classes")
once = vcftools_individual_burden("toy.vcf", str(out / "b"), burden="frequency_classes_hom_once")
print(sum(sum(r[1:]) for r in both["table"]["rows"]), sum(sum(r[1:]) for r in once["table"]["rows"]))
```

- **Example output:**

```text
['INDV', 'N_HOM_REF', 'N_HET', 'N_HOM_ALT', 'N_MISS']
[['S1', 9, 10, 1, 0], ['S2', 8, 8, 4, 0]]
78 63
```

## 3) `dnakit.variants.vcftools_long_runs_of_homozygosity()`

- **Purpose:** Find long runs of homozygosity on one chromosome.
- **Method:** Run `vcftools --LROH --chr <chromosome>`; the table is `<output_prefix>.LROH`. VCFtools requires the chromosome.
- **API:** `dnakit.variants.vcftools_long_runs_of_homozygosity(vcf_path: 'str', output_prefix: 'str', chromosome: 'str', *, filters: 'Mapping[str, object] | None' = None, max_rows: 'int' = 1000, vcftools_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the `vcftools` program (VCFtools 0.1.16 or later; `conda install -c bioconda vcftools`). DNAKit looks for it at `vcftools_path`, then on PATH, then next to the running Python interpreter; if it is missing, `BackendUnavailableError` is raised before anything is written.
- **Input:** `vcf_path`, `output_prefix`, `chromosome` (third positional argument; `filters` must not contain `chromosomes`), `filters`, `max_rows`.
- **Output:** The common result dictionary with `table` (VCFtools reports runs with their start, end and length; the table is empty when there is none).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants import vcftools_long_runs_of_homozygosity

result = vcftools_long_runs_of_homozygosity("toy.vcf", str(Path(tempfile.mkdtemp()) / "roh"), "chr1")
print(result["status"], result["sites_kept"], Path(result["outputs"][0]).suffix)
print(result["table"]["columns"], result["table"]["row_count"])
```

- **Example output:**

```text
completed 14 .LROH
['CHROM', 'AUTO_START', 'AUTO_END', 'MIN_START', 'MAX_END', 'N_VARIANTS_BETWEEN_MAX_BOUNDARIES', 'N_MISMATCHES', 'INDV'] 0
```

## 4) `dnakit.variants.vcftools_relatedness()`

- **Purpose:** Estimate the pairwise relatedness between individuals.
- **Method:** Run `vcftools --relatedness` (`method="yang"`, Yang et al. 2010, `.relatedness`, statistic `RELATEDNESS_AJK`) or `--relatedness2` (`"manichaikul"`, Manichaikul et al. 2010, `.relatedness2`, statistic `RELATEDNESS_PHI`).
- **API:** `dnakit.variants.vcftools_relatedness(vcf_path: 'str', output_prefix: 'str', *, method: 'str' = 'yang', filters: 'Mapping[str, object] | None' = None, max_rows: 'int' = 1000, vcftools_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the `vcftools` program (VCFtools 0.1.16 or later; `conda install -c bioconda vcftools`). DNAKit looks for it at `vcftools_path`, then on PATH, then next to the running Python interpreter; if it is missing, `BackendUnavailableError` is raised before anything is written.
- **Input:** `vcf_path`, `output_prefix`, `method`, `filters`, `max_rows`.
- **Output:** The common result dictionary with `table` (all pairs including self pairs).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants import vcftools_relatedness

out = Path(tempfile.mkdtemp())
yang = vcftools_relatedness("toy.vcf", str(out / "y"), max_rows=2)
print(yang["table"]["columns"], yang["table"]["row_count"], yang["table"]["rows"])
phi = vcftools_relatedness("toy.vcf", str(out / "m"), method="manichaikul", max_rows=2)
print(phi["table"]["columns"])
print(phi["table"]["rows"])
```

- **Example output:**

```text
['INDV1', 'INDV2', 'RELATEDNESS_AJK'] 21 [['S1', 'S1', 0.632232], ['S1', 'S2', -0.0434305]]
['INDV1', 'INDV2', 'N_AaAa', 'N_AAaa', 'N1_Aa', 'N2_Aa', 'RELATEDNESS_PHI']
[['S1', 'S1', 10, 0, 10, 10, 0.5], ['S1', 'S2', 5, 1, 10, 8, 0.166667]]
```

## 5) `dnakit.variants.vcftools_linkage_disequilibrium()`

- **Purpose:** Compute pairwise linkage disequilibrium between sites.
- **Method:** `statistic="geno_r2"` runs `--geno-r2` (r2 from the correlation of 0/1/2 genotypes, `.geno.ld`), `"geno_chisq"` runs `--geno-chisq` (chi-squared test of genotype independence, `.geno.chisq`), `"hap_r2"` runs `--hap-r2` (r2 from phased haplotypes, `.hap.ld`). With `positions_file` only the listed sites are compared with all other sites (`--geno-r2-positions` `.list.geno.ld`, `--hap-r2-positions` `.list.hap.ld`). With `interchromosomal=True` sites on different chromosomes are compared (`--interchrom-geno-r2` `.interchrom.geno.ld`, `--interchrom-hap-r2` `.interchrom.hap.ld`). `ld_window_bp`, `ld_window_snps`, `ld_window_bp_min`, `ld_window_snps_min` and `min_r2` map to `--ld-window-bp`, `--ld-window`, `--ld-window-bp-min`, `--ld-window-min` and `--min-r2`.
- **API:** `dnakit.variants.vcftools_linkage_disequilibrium(vcf_path: 'str', output_prefix: 'str', *, statistic: 'str' = 'geno_r2', positions_file: 'str | None' = None, interchromosomal: 'bool' = False, ld_window_bp: 'int | None' = None, ld_window_snps: 'int | None' = None, ld_window_bp_min: 'int | None' = None, ld_window_snps_min: 'int | None' = None, min_r2: 'float | None' = None, filters: 'Mapping[str, object] | None' = None, max_rows: 'int' = 1000, vcftools_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the `vcftools` program (VCFtools 0.1.16 or later; `conda install -c bioconda vcftools`). DNAKit looks for it at `vcftools_path`, then on PATH, then next to the running Python interpreter; if it is missing, `BackendUnavailableError` is raised before anything is written.
- **Input:** `vcf_path`, `output_prefix`, `statistic`, `positions_file` (tab-separated CHROM, POS), `interchromosomal`, the window options, `min_r2` (0 to 1), `filters`, `max_rows`. The chi-squared test has no list or inter-chromosomal form; `positions_file` and `interchromosomal` cannot be combined.
- **Output:** The common result dictionary with `table` (`CHR, POS1, POS2, N_INDV, R^2`; `N_CHR, R^2, D, Dprime` for haplotypes; chi-squared columns for `geno_chisq`).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants import vcftools_linkage_disequilibrium

out = Path(tempfile.mkdtemp())
pairs = vcftools_linkage_disequilibrium("toy.vcf", str(out / "g"), max_rows=2)
print(pairs["table"]["columns"], pairs["table"]["row_count"])
print(pairs["table"]["rows"])
near = vcftools_linkage_disequilibrium("toy.vcf", str(out / "w"), ld_window_bp=1500, min_r2=0.1)
print(near["table"]["row_count"])
hap = vcftools_linkage_disequilibrium("toy.vcf", str(out / "h"), statistic="hap_r2", max_rows=1)
print(hap["table"]["columns"], hap["table"]["rows"])
listed = vcftools_linkage_disequilibrium("toy.vcf", str(out / "l"), positions_file="pos.txt")
print(Path(listed["outputs"][0]).name, listed["table"]["row_count"])
inter = vcftools_linkage_disequilibrium("toy.vcf", str(out / "i"), interchromosomal=True)
print(Path(inter["outputs"][0]).name, inter["table"]["row_count"])
```

- **Example output:**

```text
['CHR', 'POS1', 'POS2', 'N_INDV', 'R^2'] 106
[['chr1', 100, 800, 6, 0.0857143], ['chr1', 100, 1500, 6, 0.294118]]
16
['CHR', 'POS1', 'POS2', 'N_CHR', 'R^2', 'D', 'Dprime'] [['chr1', 100, 800, 12, 0.0857143, 0.0625, 0.428571]]
l.list.geno.ld 13
i.interchrom.geno.ld 84
```

## 6) `dnakit.variants.vcftools_haplotype_count()`

- **Purpose:** Count the unique haplotypes in each BED region of a phased VCF.
- **Method:** Run `vcftools --hapcount <regions_bed>`; the table is `<output_prefix>.hapcount`. The regions must not overlap. VCFtools skips the first BED line (a header is expected) and can zero the last row if the last VCF record lies in the last region; `warnings` says so.
- **API:** `dnakit.variants.vcftools_haplotype_count(vcf_path: 'str', output_prefix: 'str', regions_bed: 'str', *, filters: 'Mapping[str, object] | None' = None, max_rows: 'int' = 1000, vcftools_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the `vcftools` program (VCFtools 0.1.16 or later; `conda install -c bioconda vcftools`). DNAKit looks for it at `vcftools_path`, then on PATH, then next to the running Python interpreter; if it is missing, `BackendUnavailableError` is raised before anything is written.
- **Input:** `vcf_path` (phased), `output_prefix`, `regions_bed`, `filters`, `max_rows`.
- **Output:** The common result dictionary with `table` (`CHROM, BIN_START, BIN_END, N_SNP, N_UNIQ_HAPS, N_GROUPS` and the multiplicity:frequency pairs) and `warnings`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants import vcftools_haplotype_count

result = vcftools_haplotype_count("toy.vcf", str(Path(tempfile.mkdtemp()) / "hap"), "win.bed")
for row in result["table"]["rows"]:
    print(row[:5])
print(len(result["warnings"]))
```

- **Example output:**

```text
['chr1', 0, 5000, 8, 11]
['chr1', 5000, 12000, 5, 9]
['chr2', 0, 4000, 5, 9]
1
```

## 7) `dnakit.variants.vcftools_pca()`

- **Purpose:** Run a principal component analysis of the genotypes, or write the SNP loadings.
- **Method:** Run `vcftools --pca` (`output="components"`, `.pca`; `normalise=False` runs `--pca-no-norm`) or `vcftools --pca-snp-loadings <n>` (`output="snp_loadings"`, `.pca.loadings`, only with normalisation), following Patterson, Price and Reich 2006. VCFtools needs bi-allelic diploid sites and more sites than individuals, and a build with PCA support (`--enable-pca`).
- **API:** `dnakit.variants.vcftools_pca(vcf_path: 'str', output_prefix: 'str', *, output: 'str' = 'components', normalise: 'bool' = True, snp_loadings: 'int | None' = None, filters: 'Mapping[str, object] | None' = None, max_rows: 'int' = 1000, vcftools_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the `vcftools` program (VCFtools 0.1.16 or later; `conda install -c bioconda vcftools`). DNAKit looks for it at `vcftools_path`, then on PATH, then next to the running Python interpreter; if it is missing, `BackendUnavailableError` is raised before anything is written.
- **Input:** `vcf_path`, `output_prefix`, `output`, `normalise`, `snp_loadings` (number of components, required for loadings), `filters`, `max_rows`.
- **Output:** The common result dictionary with `table`: `INDV, EIG_0, ...` whose first row, `EIGENVALUE`, holds the eigenvalues and which has one further row per individual, or `CHROM, POS, GAMMA_0, ...`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants import vcftools_pca

out = Path(tempfile.mkdtemp())
pcs = vcftools_pca("toy.vcf", str(out / "pca"))
print(pcs["status"], pcs["table"]["columns"][:3], pcs["table"]["row_count"])
print(pcs["table"]["rows"][0][:3])
loadings = vcftools_pca("toy.vcf", str(out / "load"), output="snp_loadings", snp_loadings=2)
print(loadings["table"]["columns"], loadings["table"]["row_count"], loadings["table"]["rows"][0])
```

- **Example output:**

```text
completed ['INDV', 'EIG_0', 'EIG_1'] 7
['EIGENVALUE', 3.03327, 2.69104]
['CHROM', 'POS', 'GAMMA_0', 'GAMMA_1'] 19 ['chr1', 100, 0.493707, -0.370979]
```

## 8) `dnakit.variants.vcftools_export_genotypes()`

- **Purpose:** Export the genotypes of a VCF to the input format of another program.
- **Method:** `export_format` selects the VCFtools option: `"012"` (`--012`; `.012`, `.012.indv`, `.012.pos`), `"impute"` (`--IMPUTE`; `.impute.hap`, `.impute.hap.indv`, `.impute.legend`; completely phased bi-allelic SNPs), `"ldhat"` (`--ldhat`; `.ldhat.sites`, `.ldhat.locs`), `"ldhelmet"` (`--ldhelmet`; `.ldhelmet.snps`, `.ldhelmet.pos`), `"ldhat_geno"` (`--ldhat-geno`; unphased LDhat), `"beagle_gl"` and `"beagle_pl"` (`--BEAGLE-GL`, `--BEAGLE-PL`; `.BEAGLE.GL`, `.BEAGLE.PL`; need a GL or PL FORMAT field), `"plink"` (`--plink`; `.ped`, `.map`) and `"plink_tped"` (`--plink-tped`; `.tped`, `.tfam`). LDhat, LDhelmet and BEAGLE need one `chromosome`.
- **API:** `dnakit.variants.vcftools_export_genotypes(vcf_path: 'str', output_prefix: 'str', *, export_format: 'str' = '012', chromosome: 'str | None' = None, filters: 'Mapping[str, object] | None' = None, vcftools_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the `vcftools` program (VCFtools 0.1.16 or later; `conda install -c bioconda vcftools`). DNAKit looks for it at `vcftools_path`, then on PATH, then next to the running Python interpreter; if it is missing, `BackendUnavailableError` is raised before anything is written.
- **Input:** `vcf_path`, `output_prefix`, `export_format`, `chromosome`, `filters`.
- **Output:** The common result dictionary; `outputs` lists every written file.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants import vcftools_export_genotypes

out = Path(tempfile.mkdtemp())
for fmt in ("012", "impute", "plink", "plink_tped"):
    result = vcftools_export_genotypes("toy.vcf", str(out / fmt), export_format=fmt)
    print(fmt, result["status"], sorted(Path(p).name.split(".", 1)[1] for p in result["outputs"]))
beagle = vcftools_export_genotypes("toy.vcf", str(out / "bgl"), export_format="beagle_pl", chromosome="chr1")
print(beagle["status"], beagle["sites_kept"])
print(Path(out / "plink.ped").read_text().splitlines()[0][:60])
print(Path(out / "012.012").read_text().splitlines()[0])
```

- **Example output:**

```text
012 completed ['012', '012.indv', '012.pos']
impute completed ['impute.hap', 'impute.hap.indv', 'impute.legend']
plink completed ['map', 'ped']
plink_tped completed ['tfam', 'tped']
completed 14
S1	S1	0	0	0	0	T	G	A	G	G	G	AAC	A	A	G	C	C	A	A	G	G	T	TGT	C	C	C	
0	1	1	2	1	1	0	0	0	1	0	0	1	1	1	0	1	1	0	0	0
```

## 9) `dnakit.variants.vcftools_compare_vcfs()`

- **Purpose:** Compare two VCF files: shared individuals, genotype discordance per site, or phase switch errors.
- **Method:** Run `vcftools --vcf <first> --diff <second>` with `--diff-indv` (`mode="individuals"`, `.diff.indv_in_files`, FILES is `1`, `2` or `B` for both), `--diff-site-discordance` (`"site_discordance"`, `.diff.sites`) or `--diff-switch-error` (`"switch_error"`, `.diff.switch` and `.diff.indv.switch`). The second file may be VCF, VCF.gz or BCF. `individual_map_file` renames individuals of the second file (`--diff-indv-map`). `filters` apply to the first file.
- **API:** `dnakit.variants.vcftools_compare_vcfs(vcf_path: 'str', second_vcf_path: 'str', output_prefix: 'str', *, mode: 'str' = 'individuals', individual_map_file: 'str | None' = None, filters: 'Mapping[str, object] | None' = None, max_rows: 'int' = 1000, vcftools_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the `vcftools` program (VCFtools 0.1.16 or later; `conda install -c bioconda vcftools`). DNAKit looks for it at `vcftools_path`, then on PATH, then next to the running Python interpreter; if it is missing, `BackendUnavailableError` is raised before anything is written.
- **Input:** `vcf_path`, `second_vcf_path`, `output_prefix`, `mode`, `individual_map_file`, `filters`, `max_rows`.
- **Output:** The common result dictionary with `table` of the first output file.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants import vcftools_compare_vcfs

lines = Path("toy.vcf").read_text().splitlines()
changed = [x.replace("1|1:", "0|0:", 1) if x.startswith("chr1\t800") else x for x in lines]
Path("toy_changed.vcf").write_text("\n".join(changed) + "\n")
out = Path(tempfile.mkdtemp())
people = vcftools_compare_vcfs("toy.vcf", "toy_changed.vcf", str(out / "i"))
print(people["table"]["columns"], people["table"]["rows"][:2])
sites = vcftools_compare_vcfs("toy.vcf", "toy_changed.vcf", str(out / "s"), mode="site_discordance")
print(sites["table"]["columns"])
print([row for row in sites["table"]["rows"] if row[5] > 0])
switch = vcftools_compare_vcfs("toy.vcf", "toy_changed.vcf", str(out / "w"), mode="switch_error")
print(switch["status"], sorted(Path(p).name.split(".", 1)[1] for p in switch["outputs"]))
```

- **Example output:**

```text
['INDV', 'FILES'] [['S1', 'B'], ['S2', 'B']]
['CHROM', 'POS', 'FILES', 'MATCHING_ALLELES', 'N_COMMON_CALLED', 'N_DISCORD', 'DISCORDANCE']
[['chr1', 800, 'B', 1, 6, 1, 0.166667]]
completed ['diff.indv.switch', 'diff.switch']
```
