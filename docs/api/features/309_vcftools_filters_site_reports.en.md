# VCFtools: filtering and site, depth and quality reports

This group of functions in `dnakit.variants` is a set of typed command-line facades over the `vcftools` program of [VCFtools](https://github.com/vcftools/vcftools) (tested with 0.1.16). Each function runs one fixed `vcftools` command with a typed argument list and without a shell. They filter sites, individuals and genotypes, list the sites that pass a filter, and write depth, quality, missingness, singleton, indel-length, SNP-density, Ts/Tv, per-site diversity and INFO/FORMAT-field tables. Only the options of the VCFtools operations they cover are exposed.

VCFtools is distributed under the LGPL-3.0. DNAKit runs `vcftools` as a separate process and copies none of its code. The program is not a Python dependency: install it yourself (`conda install -c bioconda vcftools`, `apt install vcftools`, or build it from source) and put it on PATH or pass `vcftools_path`. If it is missing the function raises `BackendUnavailableError` with an install hint.

Common behaviour. `vcf_path` is a `.vcf`, a `.vcf.gz` (read with `--gzvcf`) or a `.bcf`. `output_prefix` is a file-name prefix, exactly as `--out` in VCFtools: the result files are `<output_prefix>.<suffix>` and the suffix is named in each section. A result file or the prefix folder that does not fit raises `ConfigurationError`, and an existing output is never overwritten. The report functions take `filters`, a mapping that uses the keyword names of `vcftools_filter_vcf()` (for example `{"min_maf": 0.05, "remove_filtered_all": True}`), so every report can be computed on a filtered subset. Defaults are those of VCFtools: an option is only passed when it is given. A tool that exits with an error does not raise: `status` is `failed` and `log_tail` holds the message. Every result is a JSON-serialisable dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (the files that exist), `log_tail`, `output_prefix` and, when VCFtools reports them, `sites_kept`, `sites_total`, `individuals_kept` and `individuals_total`. The report functions also return `table`: `columns`, the first `max_rows` rows (`rows`, numbers converted, `nan` as `None`), `row_count` and `truncated`.

The examples import from `dnakit.variants`, call every function with the default verbosity and write to a fresh temporary directory.

## Example data

The only declared input is `toy.vcf`: six phased samples (S1 to S6), 20 sites on `chr1` and `chr2`, three indels, two `LowQual` sites, four singletons, one doubleton, one missing genotype and the FORMAT fields GT, DP, GQ and PL. `chr1` has 14 sites and `chr2` has 6.

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

## 1) `dnakit.variants.vcftools_filter_vcf()`

- **Purpose:** Filter the sites, individuals and genotypes of a VCF or BCF with VCFtools and write the result as VCF or BCF.
- **Method:** Build one `vcftools` command from the given options and add `--recode` (`output_format="vcf"`) or `--recode-bcf` (`"bcf"`). All given filters are applied together, as in VCFtools. Site filters: `chromosomes`, `exclude_chromosomes` (`--chr`, `--not-chr`), `from_bp`, `to_bp`, `positions_file`, `exclude_positions_file`, `positions_overlap_file`, `exclude_positions_overlap_file`, `bed_file`, `exclude_bed_file`, `thin` (minimum spacing in bp), `mask_file`, `invert_mask_file`, `mask_min`, `snp_ids`, `snp_ids_file`, `exclude_ids_file`, `only_indels`, `remove_indels`, `remove_filtered_all`, `keep_filtered`, `remove_filtered`, `keep_info`, `remove_info`, `min_maf`, `max_maf`, `min_mac`, `max_mac`, `min_alleles`, `max_alleles`, `min_mean_depth`, `max_mean_depth`, `hwe_pvalue`, `phased_only`, `min_quality`, `max_missing`. Individual filters: `keep_individuals`, `keep_individuals_files`, `remove_individuals`, `remove_individuals_files`, `max_individuals` (random thinning). Genotype filters: `remove_filtered_genotypes`, `remove_filtered_genotype_flags`, `min_genotype_quality`, `min_genotype_depth`, `max_genotype_depth`. Like VCFtools, INFO is dropped from the output unless `recode_info` is `True` (all INFO) or a list of INFO keys.
- **API:** `dnakit.variants.vcftools_filter_vcf(vcf_path: 'str', output_prefix: 'str', *, output_format: 'str' = 'vcf', recode_info: 'bool | Sequence[str]' = False, chromosomes: 'Sequence[str] | None' = None, exclude_chromosomes: 'Sequence[str] | None' = None, from_bp: 'int | None' = None, to_bp: 'int | None' = None, positions_file: 'str | None' = None, exclude_positions_file: 'str | None' = None, positions_overlap_file: 'str | None' = None, exclude_positions_overlap_file: 'str | None' = None, bed_file: 'str | None' = None, exclude_bed_file: 'str | None' = None, thin: 'int | None' = None, mask_file: 'str | None' = None, invert_mask_file: 'str | None' = None, mask_min: 'int | None' = None, snp_ids: 'Sequence[str] | None' = None, snp_ids_file: 'str | None' = None, exclude_ids_file: 'str | None' = None, only_indels: 'bool' = False, remove_indels: 'bool' = False, remove_filtered_all: 'bool' = False, keep_filtered: 'Sequence[str] | None' = None, remove_filtered: 'Sequence[str] | None' = None, keep_info: 'Sequence[str] | None' = None, remove_info: 'Sequence[str] | None' = None, min_maf: 'float | None' = None, max_maf: 'float | None' = None, min_mac: 'int | None' = None, max_mac: 'int | None' = None, min_alleles: 'int | None' = None, max_alleles: 'int | None' = None, min_mean_depth: 'float | None' = None, max_mean_depth: 'float | None' = None, hwe_pvalue: 'float | None' = None, phased_only: 'bool' = False, min_quality: 'float | None' = None, max_missing: 'float | None' = None, keep_individuals: 'Sequence[str] | None' = None, keep_individuals_files: 'Sequence[str] | None' = None, remove_individuals: 'Sequence[str] | None' = None, remove_individuals_files: 'Sequence[str] | None' = None, max_individuals: 'int | None' = None, remove_filtered_genotypes: 'bool' = False, remove_filtered_genotype_flags: 'Sequence[str] | None' = None, min_genotype_quality: 'float | None' = None, min_genotype_depth: 'int | None' = None, max_genotype_depth: 'int | None' = None, vcftools_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the `vcftools` program (VCFtools 0.1.16 or later; `conda install -c bioconda vcftools`). DNAKit looks for it at `vcftools_path`, then on PATH, then next to the running Python interpreter; if it is missing, `BackendUnavailableError` is raised before anything is written.
- **Input:** `vcf_path` (VCF, VCF.gz or BCF), `output_prefix`, `output_format` (`"vcf"` or `"bcf"`), `recode_info`, and the filter options above. Lists are `Sequence[str]`; file options are existing files (tab-separated CHROM, POS for position lists; CHROM, START, END for BED). Numbers are checked before the program starts.
- **Output:** The common result dictionary. `outputs` holds `<output_prefix>.recode.vcf` or `<output_prefix>.recode.bcf`; `sites_kept`/`sites_total` and `individuals_kept`/`individuals_total` come from the VCFtools log.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants import vcftools_filter_vcf

out = Path(tempfile.mkdtemp())
kept = vcftools_filter_vcf(
    "toy.vcf", str(out / "pass"), min_maf=0.2, remove_filtered_all=True, min_genotype_quality=20
)
print(kept["status"], kept["sites_kept"], kept["sites_total"], kept["individuals_kept"])
subset = vcftools_filter_vcf(
    "toy.vcf", str(out / "sub"), keep_individuals=["S1", "S2"], remove_indels=True, recode_info=["DP"]
)
print(subset["sites_kept"], subset["individuals_kept"], Path(subset["outputs"][0]).name)
print(Path(subset["outputs"][0]).read_text().splitlines()[-1])
bcf = vcftools_filter_vcf("toy.vcf", str(out / "b"), output_format="bcf", chromosomes=["chr2"])
print(Path(bcf["outputs"][0]).name, bcf["sites_kept"])
```

- **Example output:**

```text
completed 14 20 6
17 2 sub.recode.vcf
chr2	4800	v20	T	C	67	PASS	DP=154	GT:DP:GQ:PL	0|0:10:30:0,30,300	0|0:38:12:0,30,300
b.recode.bcf 6
```

## 2) `dnakit.variants.vcftools_list_sites()`

- **Purpose:** List the sites that pass, or that fail, a set of VCFtools filters.
- **Method:** Run `vcftools` with the filters of `filters` and `--kept-sites` (`which="kept"`, file `<output_prefix>.kept.sites`) or `--removed-sites` (`which="removed"`, file `<output_prefix>.removed.sites`), then read the two-column table.
- **API:** `dnakit.variants.vcftools_list_sites(vcf_path: 'str', output_prefix: 'str', *, which: 'str' = 'kept', filters: 'Mapping[str, object] | None' = None, max_rows: 'int' = 1000, vcftools_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the `vcftools` program (VCFtools 0.1.16 or later; `conda install -c bioconda vcftools`). DNAKit looks for it at `vcftools_path`, then on PATH, then next to the running Python interpreter; if it is missing, `BackendUnavailableError` is raised before anything is written.
- **Input:** `vcf_path`, `output_prefix`, `which` (`"kept"` or `"removed"`), `filters` (keys of `vcftools_filter_vcf()`), `max_rows`.
- **Output:** The common result dictionary with `table` (columns `CHROM`, `POS`).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants import vcftools_list_sites

out = Path(tempfile.mkdtemp())
flt = {"min_maf": 0.3, "remove_filtered_all": True}
kept = vcftools_list_sites("toy.vcf", str(out / "k"), filters=flt)
removed = vcftools_list_sites("toy.vcf", str(out / "r"), which="removed", filters=flt)
print(kept["table"]["columns"], kept["table"]["row_count"], removed["table"]["row_count"])
print(kept["table"]["rows"][:3])
print(Path(removed["outputs"][0]).name)
```

- **Example output:**

```text
['CHROM', 'POS'] 9 11
[['chr1', 100], ['chr1', 1500], ['chr1', 5000]]
r.removed.sites
```

## 3) `dnakit.variants.vcftools_depth()`

- **Purpose:** Report read depth per individual, per site or per genotype from the DP genotype field.
- **Method:** Run `vcftools` with `--depth` (`level="individual_mean"`, `.idepth`), `--site-depth` (`"site_sum"`, `.ldepth`), `--site-mean-depth` (`"site_mean"`, `.ldepth.mean`) or `--geno-depth` (`"genotype"`, `.gdepth`).
- **API:** `dnakit.variants.vcftools_depth(vcf_path: 'str', output_prefix: 'str', *, level: 'str' = 'individual_mean', filters: 'Mapping[str, object] | None' = None, max_rows: 'int' = 1000, vcftools_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the `vcftools` program (VCFtools 0.1.16 or later; `conda install -c bioconda vcftools`). DNAKit looks for it at `vcftools_path`, then on PATH, then next to the running Python interpreter; if it is missing, `BackendUnavailableError` is raised before anything is written.
- **Input:** `vcf_path`, `output_prefix`, `level`, `filters`, `max_rows`. The VCF must carry FORMAT/DP.
- **Output:** The common result dictionary with `table`; columns depend on `level` (`INDV, N_SITES, MEAN_DEPTH`; `CHROM, POS, SUM_DEPTH, SUMSQ_DEPTH`; `CHROM, POS, MEAN_DEPTH, VAR_DEPTH`; or one column per individual).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants import vcftools_depth

out = Path(tempfile.mkdtemp())
per_sample = vcftools_depth("toy.vcf", str(out / "i"))
print(per_sample["table"]["columns"], per_sample["table"]["rows"][:2])
per_site = vcftools_depth("toy.vcf", str(out / "m"), level="site_mean", max_rows=2)
print(per_site["table"]["columns"], per_site["table"]["rows"], per_site["table"]["truncated"])
genotype = vcftools_depth("toy.vcf", str(out / "g"), level="genotype", max_rows=1)
print(genotype["table"]["columns"], genotype["table"]["rows"][0])
```

- **Example output:**

```text
['INDV', 'N_SITES', 'MEAN_DEPTH'] [['S1', 20, 22.1], ['S2', 20, 25.5]]
['CHROM', 'POS', 'MEAN_DEPTH', 'VAR_DEPTH'] [['chr1', 100, 15.1667, 80.5667], ['chr1', 800, 20, 28.4]] True
['CHROM', 'POS', 'S1', 'S2', 'S3', 'S4', 'S5', 'S6'] ['chr1', 100, 9, 31, 9, 20, 8, 14]
```

## 4) `dnakit.variants.vcftools_site_quality()`

- **Purpose:** Extract the QUAL value of every site.
- **Method:** Run `vcftools --site-quality`; the table is `<output_prefix>.lqual`.
- **API:** `dnakit.variants.vcftools_site_quality(vcf_path: 'str', output_prefix: 'str', *, filters: 'Mapping[str, object] | None' = None, max_rows: 'int' = 1000, vcftools_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the `vcftools` program (VCFtools 0.1.16 or later; `conda install -c bioconda vcftools`). DNAKit looks for it at `vcftools_path`, then on PATH, then next to the running Python interpreter; if it is missing, `BackendUnavailableError` is raised before anything is written.
- **Input:** `vcf_path`, `output_prefix`, `filters`, `max_rows`.
- **Output:** The common result dictionary with `table` (columns `CHROM`, `POS`, `QUAL`).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants import vcftools_site_quality

result = vcftools_site_quality("toy.vcf", str(Path(tempfile.mkdtemp()) / "q"), max_rows=3)
print(result["table"]["columns"])
print(result["table"]["rows"])
```

- **Example output:**

```text
['CHROM', 'POS', 'QUAL']
[['chr1', 100, 38], ['chr1', 800, 55], ['chr1', 1500, 67]]
```

## 5) `dnakit.variants.vcftools_missingness()`

- **Purpose:** Report the fraction of missing genotypes per individual or per site.
- **Method:** Run `vcftools --missing-indv` (`by="individual"`, `.imiss`) or `--missing-site` (`by="site"`, `.lmiss`).
- **API:** `dnakit.variants.vcftools_missingness(vcf_path: 'str', output_prefix: 'str', *, by: 'str' = 'individual', filters: 'Mapping[str, object] | None' = None, max_rows: 'int' = 1000, vcftools_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the `vcftools` program (VCFtools 0.1.16 or later; `conda install -c bioconda vcftools`). DNAKit looks for it at `vcftools_path`, then on PATH, then next to the running Python interpreter; if it is missing, `BackendUnavailableError` is raised before anything is written.
- **Input:** `vcf_path`, `output_prefix`, `by` (`"individual"` or `"site"`), `filters`, `max_rows`.
- **Output:** The common result dictionary with `table` (`INDV, N_DATA, N_GENOTYPES_FILTERED, N_MISS, F_MISS` or `CHR, POS, N_DATA, N_GENOTYPE_FILTERED, N_MISS, F_MISS`).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants import vcftools_missingness

out = Path(tempfile.mkdtemp())
sites = vcftools_missingness("toy.vcf", str(out / "s"), by="site")
print([row for row in sites["table"]["rows"] if row[4] > 0])
people = vcftools_missingness("toy.vcf", str(out / "i"))
print(people["table"]["columns"])
print(people["table"]["rows"][4])
```

- **Example output:**

```text
[['chr1', 6400, 12, 0, 2, 0.166667]]
['INDV', 'N_DATA', 'N_GENOTYPES_FILTERED', 'N_MISS', 'F_MISS']
['S5', 20, 0, 1, 0.05]
```

## 6) `dnakit.variants.vcftools_singletons()`

- **Purpose:** Locate singleton and doubleton alleles.
- **Method:** Run `vcftools --singletons`; the table is `<output_prefix>.singletons`. A singleton allele is seen once in the sample, a doubleton twice in one individual.
- **API:** `dnakit.variants.vcftools_singletons(vcf_path: 'str', output_prefix: 'str', *, filters: 'Mapping[str, object] | None' = None, max_rows: 'int' = 1000, vcftools_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the `vcftools` program (VCFtools 0.1.16 or later; `conda install -c bioconda vcftools`). DNAKit looks for it at `vcftools_path`, then on PATH, then next to the running Python interpreter; if it is missing, `BackendUnavailableError` is raised before anything is written.
- **Input:** `vcf_path`, `output_prefix`, `filters`, `max_rows`.
- **Output:** The common result dictionary with `table` (columns `CHROM, POS, SINGLETON/DOUBLETON, ALLELE, INDV`).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants import vcftools_singletons

result = vcftools_singletons("toy.vcf", str(Path(tempfile.mkdtemp()) / "s"))
for row in result["table"]["rows"]:
    print(row)
```

- **Example output:**

```text
['chr1', 3600, 'S', 'T', 'S3']
['chr1', 4300, 'D', 'G', 'S2']
['chr1', 6400, 'S', 'T', 'S3']
['chr2', 3000, 'S', 'A', 'S4']
['chr2', 4800, 'S', 'C', 'S6']
```

## 7) `dnakit.variants.vcftools_indel_length_histogram()`

- **Purpose:** Count the lengths of the indels in a VCF.
- **Method:** Run `vcftools --hist-indel-len`; the table is `<output_prefix>.indel.hist`. Negative lengths are deletions, positive lengths insertions and 0 are SNPs.
- **API:** `dnakit.variants.vcftools_indel_length_histogram(vcf_path: 'str', output_prefix: 'str', *, filters: 'Mapping[str, object] | None' = None, max_rows: 'int' = 1000, vcftools_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the `vcftools` program (VCFtools 0.1.16 or later; `conda install -c bioconda vcftools`). DNAKit looks for it at `vcftools_path`, then on PATH, then next to the running Python interpreter; if it is missing, `BackendUnavailableError` is raised before anything is written.
- **Input:** `vcf_path`, `output_prefix`, `filters`, `max_rows`.
- **Output:** The common result dictionary with `table` (columns `LENGTH, COUNT, PRCT`).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants import vcftools_indel_length_histogram

result = vcftools_indel_length_histogram("toy.vcf", str(Path(tempfile.mkdtemp()) / "h"))
print(result["table"]["columns"])
print(result["table"]["rows"])
```

- **Example output:**

```text
['LENGTH', 'COUNT', 'PRCT']
[[-2, 1, 5], [0, 17, 85], [2, 2, 10]]
```

## 8) `dnakit.variants.vcftools_snp_density()`

- **Purpose:** Count the SNPs in non-overlapping windows along each chromosome.
- **Method:** Run `vcftools --SNPdensity <bin_size>`; the table is `<output_prefix>.snpden`.
- **API:** `dnakit.variants.vcftools_snp_density(vcf_path: 'str', output_prefix: 'str', bin_size: 'int', *, filters: 'Mapping[str, object] | None' = None, max_rows: 'int' = 1000, vcftools_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the `vcftools` program (VCFtools 0.1.16 or later; `conda install -c bioconda vcftools`). DNAKit looks for it at `vcftools_path`, then on PATH, then next to the running Python interpreter; if it is missing, `BackendUnavailableError` is raised before anything is written.
- **Input:** `vcf_path`, `output_prefix`, `bin_size` (bp, at least 1), `filters`, `max_rows`.
- **Output:** The common result dictionary with `table` (columns `CHROM, BIN_START, SNP_COUNT, VARIANTS/KB`).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants import vcftools_snp_density

result = vcftools_snp_density("toy.vcf", str(Path(tempfile.mkdtemp()) / "d"), 5000)
print(result["table"]["columns"])
for row in result["table"]["rows"]:
    print(row)
```

- **Example output:**

```text
['CHROM', 'BIN_START', 'SNP_COUNT', 'VARIANTS/KB']
['chr1', 0, 7, 1.4]
['chr1', 5000, 7, 1.4]
['chr2', 0, 6, 1.2]
```

## 9) `dnakit.variants.vcftools_transition_transversion()`

- **Purpose:** Count transitions and transversions of bi-allelic SNPs and report the Ts/Tv ratio in four groupings.
- **Method:** Run `vcftools --TsTv <bin_size>` (`mode="distance_bins"`, `.TsTv`), `--TsTv-by-count` (`"alt_allele_count"`, `.TsTv.count`), `--TsTv-by-qual` (`"qual_threshold"`, `.TsTv.qual`) or `--FILTER-summary` (`"filter_category"`, `.FILTER.summary`). VCFtools 0.1.16 prints uninitialised memory in the three `GT` columns of the highest QUAL threshold of `--TsTv-by-qual`, where no site lies above the threshold; those cells are not meaningful.
- **API:** `dnakit.variants.vcftools_transition_transversion(vcf_path: 'str', output_prefix: 'str', *, mode: 'str' = 'distance_bins', bin_size: 'int | None' = None, filters: 'Mapping[str, object] | None' = None, max_rows: 'int' = 1000, vcftools_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the `vcftools` program (VCFtools 0.1.16 or later; `conda install -c bioconda vcftools`). DNAKit looks for it at `vcftools_path`, then on PATH, then next to the running Python interpreter; if it is missing, `BackendUnavailableError` is raised before anything is written.
- **Input:** `vcf_path`, `output_prefix`, `mode`, `bin_size` (required for `distance_bins`), `filters`, `max_rows`.
- **Output:** The common result dictionary with `table`; the columns follow the VCFtools file (`CHROM, BinStart, SNP_count, Ts/Tv`; `ALT_ALLELE_COUNT, N_Ts, N_Tv, Ts/Tv`; `QUAL_THRESHOLD, N_Ts_LT_QUAL_THRESHOLD, ...`; or `FILTER, N_VARIANTS, N_Ts, N_Tv, Ts/Tv`).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants import vcftools_transition_transversion

out = Path(tempfile.mkdtemp())
bins = vcftools_transition_transversion("toy.vcf", str(out / "b"), bin_size=5000)
print(bins["table"]["columns"], bins["table"]["rows"])
by_filter = vcftools_transition_transversion("toy.vcf", str(out / "f"), mode="filter_category")
print(by_filter["table"]["rows"])
```

- **Example output:**

```text
['CHROM', 'BinStart', 'SNP_count', 'Ts/Tv'] [['chr1', 0, 6, 5], ['chr1', 5000, 5, 1.5], ['chr2', 0, 6, 0.5]]
[['PASS', 18, 8, 7, 1.14286], ['LowQual', 2, 2, 0, None]]
```

## 10) `dnakit.variants.vcftools_site_nucleotide_diversity()`

- **Purpose:** Compute the nucleotide diversity (pi) of every site.
- **Method:** Run `vcftools --site-pi`; the table is `<output_prefix>.sites.pi`.
- **API:** `dnakit.variants.vcftools_site_nucleotide_diversity(vcf_path: 'str', output_prefix: 'str', *, filters: 'Mapping[str, object] | None' = None, max_rows: 'int' = 1000, vcftools_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the `vcftools` program (VCFtools 0.1.16 or later; `conda install -c bioconda vcftools`). DNAKit looks for it at `vcftools_path`, then on PATH, then next to the running Python interpreter; if it is missing, `BackendUnavailableError` is raised before anything is written.
- **Input:** `vcf_path`, `output_prefix`, `filters`, `max_rows`.
- **Output:** The common result dictionary with `table` (columns `CHROM, POS, PI`).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants import vcftools_site_nucleotide_diversity

result = vcftools_site_nucleotide_diversity("toy.vcf", str(Path(tempfile.mkdtemp()) / "pi"), max_rows=3)
print(result["table"]["columns"])
print(result["table"]["rows"])
```

- **Example output:**

```text
['CHROM', 'POS', 'PI']
[['chr1', 100, 0.530303], ['chr1', 800, 0.409091], ['chr1', 1500, 0.530303]]
```

## 11) `dnakit.variants.vcftools_extract_fields()`

- **Purpose:** Extract INFO keys per site, or one FORMAT key per genotype, into a table.
- **Method:** Run `vcftools --get-INFO <key>` once per key (`source="info"`, `.INFO`) or `--extract-FORMAT-info <key>` (`source="format"`, `.<key>.FORMAT`). VCFtools extracts one FORMAT key per run, so `fields` has one entry for `source="format"`.
- **API:** `dnakit.variants.vcftools_extract_fields(vcf_path: 'str', output_prefix: 'str', fields: 'Sequence[str]', *, source: 'str' = 'info', filters: 'Mapping[str, object] | None' = None, max_rows: 'int' = 1000, vcftools_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the `vcftools` program (VCFtools 0.1.16 or later; `conda install -c bioconda vcftools`). DNAKit looks for it at `vcftools_path`, then on PATH, then next to the running Python interpreter; if it is missing, `BackendUnavailableError` is raised before anything is written.
- **Input:** `vcf_path`, `output_prefix`, `fields` (names made of letters, digits, `_` and `.`), `source`, `filters`, `max_rows`.
- **Output:** The common result dictionary with `table`: `CHROM, POS, REF, ALT` and the INFO keys, or `CHROM, POS` and one column per individual.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants import vcftools_extract_fields

out = Path(tempfile.mkdtemp())
info = vcftools_extract_fields("toy.vcf", str(out / "i"), ["DP", "AF"], max_rows=2)
print(info["table"]["columns"], info["table"]["rows"])
fmt = vcftools_extract_fields("toy.vcf", str(out / "f"), ["GQ"], source="format", max_rows=1)
print(Path(fmt["outputs"][0]).name, fmt["table"]["columns"], fmt["table"]["rows"])
```

- **Example output:**

```text
['CHROM', 'POS', 'REF', 'ALT', 'DP', 'AF'] [['chr1', 100, 'G', 'T', 91, 0.1], ['chr1', 800, 'G', 'A', 120, 0.25]]
f.GQ.FORMAT ['CHROM', 'POS', 'S1', 'S2', 'S3', 'S4', 'S5', 'S6'] [['chr1', 100, 45, 60, 30, 99, 12, 99]]
```
