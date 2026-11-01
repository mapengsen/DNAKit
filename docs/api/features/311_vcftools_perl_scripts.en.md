# VCFtools Perl scripts: vcf-annotate, vcf-subset, vcf-compare and the other vcf-* tools

This group of functions in `dnakit.variants` is a set of typed facades over the Perl scripts that ship with [VCFtools](https://github.com/vcftools/vcftools) (tested with 0.1.16): `vcf-annotate`, `vcf-subset`, `vcf-convert`, `vcf-to-tab`, `vcf-fix-ploidy`, `vcf-shuffle-cols`, `vcf-contrast`, `vcf-compare`, `vcf-indel-stats`, `vcf-phased-join`, `fill-aa`, `fill-fs` and `fill-ref-md5`. Each function runs one script with a fixed argument list built from typed arguments and without a shell; standard input and output are connected to files, so every function reads a file and writes a new file. The sibling `vcftools` program is covered by the two previous documents (`309_vcftools_filters_site_reports`, `310_vcftools_samples_ld_export`).

VCFtools is distributed under the LGPL-3.0; DNAKit runs the scripts as separate processes and copies none of their code. Install VCFtools yourself (`conda install -c bioconda vcftools`, `apt install vcftools` or a source build) and put the scripts on PATH or pass `script_path` (the full path of the one script). The scripts need their `Vcf.pm` library: DNAKit sets `PERL5LIB` to the folder next to the installation (`../share/perl*`) when it is not already found. Scripts that read compressed files need `tabix` (htslib) on PATH, and `fill-aa`, `fill-fs` and `fill-ref-md5` need `samtools`. If a script is missing, the function raises `BackendUnavailableError` with an install hint.

Common behaviour. `output_path` must not exist yet, its folder must exist, and a failed run leaves no partial file (`status` is `failed` and `log_tail` holds the Perl message). Results are JSON-serialisable dictionaries with `status`, `backend`, `version` (the version of `Vcf.pm`), `command`, `return_code`, `outputs`, `log_tail`, `output_path` and, for VCF outputs, `records` (the number of data lines). Several scripts iterate Perl hashes, so the order of INFO keys, header lines and per-sample rows is not stable between runs.

The examples import from `dnakit.variants`, use pysam to bgzip-compress and tabix-index the files that need it, and write to a fresh temporary directory.

## Example data

The only declared input is `toy.vcf`: six phased samples (S1 to S6), 20 sites on `chr1` (14) and `chr2` (6), three indels, two `LowQual` sites, four singletons, one doubleton, one missing genotype and the FORMAT fields GT, DP, GQ and PL. The examples that need a reference FASTA write one whose bases agree with the REF alleles of `toy.vcf`.

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

## 1) `dnakit.variants.vcf_annotate()`

- **Purpose:** Add FILTER flags, fill AC, AN, HWE and TYPE, normalise alleles, drop tags or hard-filter sites with `vcf-annotate`.
- **Method:** Run `vcf-annotate` with `-f <filters>` (the built-in site filters, for example `+/Qual=40/MinDP=5`: `+` switches on every filter with its default, `-X` excludes filter X; user-defined filter files are not supported), `-H` (`hard_filter`: drop sites whose FILTER is not PASS or "."), `-n` (`normalize_alleles`: compact REF and ALT, for example TA,TAA to T,TA), `-r <tags>` (`remove_tags`: `ID`, `QUAL`, `FILTER`, `INFO/DP`, `FORMAT/DP`, `FILTER/MinDP` or a bare tag), `--fill-AC-AN` (`fill_ac_an`, the same as the `fill-an-ac` script), `--fill-HWE` (`fill_hwe`: also adds the HWE p-value) and `--fill-type` (`fill_type`: INFO/TYPE as snp, del, ins, mnp or complex, the classification of `Vcf::event_type`). At least one operation must be chosen.
- **API:** `dnakit.variants.vcf_annotate(vcf_path: 'str', output_path: 'str', *, filters: 'str | None' = None, hard_filter: 'bool' = False, normalize_alleles: 'bool' = False, remove_tags: 'Sequence[str] | None' = None, fill_ac_an: 'bool' = False, fill_hwe: 'bool' = False, fill_type: 'bool' = False, script_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the VCFtools Perl script of this section (VCFtools 0.1.16 or later), perl, and `tabix`/`samtools` where the text says so. DNAKit looks for the script at `script_path`, then on PATH, then next to the running Python interpreter; if it is missing, `BackendUnavailableError` is raised before anything is written.
- **Input:** `vcf_path` (VCF or VCF.gz), `output_path`, `filters`, `hard_filter`, `normalize_alleles`, `remove_tags`, `fill_ac_an`, `fill_hwe`, `fill_type`, `script_path`, `timeout_seconds`.
- **Output:** The common result dictionary with `records` and the annotated VCF at `output_path`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants import vcf_annotate

out = Path(tempfile.mkdtemp())
flagged = vcf_annotate("toy.vcf", str(out / "flag.vcf"), filters="+/Qual=40")
rows = [line.split("\t") for line in Path(flagged["output_path"]).read_text().splitlines() if line[0] != "#"]
print(flagged["status"], [row[6] for row in rows][:6])
hard = vcf_annotate("toy.vcf", str(out / "hard.vcf"), filters="+/Qual=40", hard_filter=True)
print(hard["records"])
filled = vcf_annotate("toy.vcf", str(out / "fill.vcf"), fill_ac_an=True, fill_type=True)
print(sorted(Path(filled["output_path"]).read_text().splitlines()[-1].split("\t")[7].split(";")))
bare = vcf_annotate("toy.vcf", str(out / "bare.vcf"), remove_tags=["INFO/AF", "QUAL"])
print(Path(bare["output_path"]).read_text().splitlines()[-1].split("\t")[5:8])
```

- **Example output:**

```text
completed ['Qual', 'PASS', 'PASS', 'Qual', 'PASS', 'LowQual;Qual']
12
['AC=1', 'AF=0.1', 'AN=12', 'DP=154', 'TYPE=snp']
['.', 'PASS', 'DP=154']
```

## 2) `dnakit.variants.vcf_subset()`

- **Purpose:** Select samples and variant classes of a VCF with `vcf-subset`.
- **Method:** Run `vcf-subset` with `-c` (`samples`, or `samples_file` with one name per line), `-e` (`exclude_ref_only`: drop sites without an ALT allele in the selection), `-p` (`private_only`: keep sites where only the selected samples carry an alternate allele), `-a` (`trim_alt_alleles`: drop ALT alleles the selection does not use), `-t` (`variant_types`, among `ref`, `SNPs`, `indels`, `MNPs`, `other`), `-r` (`replace_with_ref`: excluded types become the reference allele instead of "."), `-u` (`keep_uncalled`) and `-f` (`force`: continue when a sample is absent).
- **API:** `dnakit.variants.vcf_subset(vcf_path: 'str', output_path: 'str', *, samples: 'Sequence[str] | None' = None, samples_file: 'str | None' = None, exclude_ref_only: 'bool' = False, private_only: 'bool' = False, trim_alt_alleles: 'bool' = False, variant_types: 'Sequence[str] | None' = None, replace_with_ref: 'bool' = False, keep_uncalled: 'bool' = False, force: 'bool' = False, script_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the VCFtools Perl script of this section (VCFtools 0.1.16 or later), perl, and `tabix`/`samtools` where the text says so. DNAKit looks for the script at `script_path`, then on PATH, then next to the running Python interpreter; if it is missing, `BackendUnavailableError` is raised before anything is written.
- **Input:** `vcf_path`, `output_path`, `samples` or `samples_file` (not both), `exclude_ref_only`, `private_only`, `trim_alt_alleles`, `variant_types`, `replace_with_ref`, `keep_uncalled`, `force`, `script_path`, `timeout_seconds`.
- **Output:** The common result dictionary with `records`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants import vcf_subset

out = Path(tempfile.mkdtemp())
two = vcf_subset("toy.vcf", str(out / "two.vcf"), samples=["S1", "S2"])
header = [x for x in Path(two["output_path"]).read_text().splitlines() if x.startswith("#CHROM")][0]
print(two["records"], header.split("\t")[9:])
indels = vcf_subset("toy.vcf", str(out / "indel.vcf"), samples=["S1", "S2"], variant_types=["indels"], exclude_ref_only=True)
print(indels["records"])
private = vcf_subset("toy.vcf", str(out / "private.vcf"), samples=["S3"], private_only=True)
print(private["records"])
```

- **Example output:**

```text
20 ['S1', 'S2']
3
2
```

## 3) `dnakit.variants.vcf_convert()`

- **Purpose:** Convert a VCF between format versions 4.0, 4.1 and 4.2 with `vcf-convert`.
- **Method:** Run `vcf-convert -v <target_version>` with the VCF on standard input and, when given, `-r <reference_fasta>`. The reference (a FASTA indexed with `samtools faidx`) is required only when the file has indels. Downgrading is marked experimental by the script.
- **API:** `dnakit.variants.vcf_convert(vcf_path: 'str', output_path: 'str', *, target_version: 'str' = '4.2', reference_fasta: 'str | None' = None, script_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the VCFtools Perl script of this section (VCFtools 0.1.16 or later), perl, and `tabix`/`samtools` where the text says so. DNAKit looks for the script at `script_path`, then on PATH, then next to the running Python interpreter; if it is missing, `BackendUnavailableError` is raised before anything is written.
- **Input:** `vcf_path`, `output_path`, `target_version` (`"4.0"`, `"4.1"` or `"4.2"`, default `"4.2"`), `reference_fasta`, `script_path`, `timeout_seconds`.
- **Output:** The common result dictionary with `records`; a file with indels and no reference gives `status` `failed` and the script message in `log_tail`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants import vcf_convert

lines = Path("toy.vcf").read_text().splitlines()
snps = [x for x in lines if x[0] == "#" or (len(x.split("\t")[3]) == 1 and len(x.split("\t")[4]) == 1)]
Path("snps.vcf").write_text("\n".join(snps) + "\n")
out = Path(tempfile.mkdtemp())
converted = vcf_convert("snps.vcf", str(out / "v40.vcf"), target_version="4.0")
print(converted["status"], converted["records"], Path(converted["output_path"]).read_text().splitlines()[0])
refused = vcf_convert("toy.vcf", str(out / "indels.vcf"), target_version="4.0")
print(refused["status"], "Indels present" in refused["log_tail"])
```

- **Example output:**

```text
completed 17 ##fileformat=VCFv4.0
failed True
```

## 4) `dnakit.variants.vcf_to_tab()`

- **Purpose:** Convert a VCF to a tab-delimited table of genotypes written as bases, with `vcf-to-tab`.
- **Method:** Run `vcf-to-tab` with the VCF on standard input and, for `iupac=True`, `-i` (heterozygous genotypes as one-letter IUPAC codes; the script stops on indels). The columns are CHROM, POS, REF and one column of bases per sample. The script itself says that `bcftools query` should be used instead in the future.
- **API:** `dnakit.variants.vcf_to_tab(vcf_path: 'str', output_path: 'str', *, iupac: 'bool' = False, script_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the VCFtools Perl script of this section (VCFtools 0.1.16 or later), perl, and `tabix`/`samtools` where the text says so. DNAKit looks for the script at `script_path`, then on PATH, then next to the running Python interpreter; if it is missing, `BackendUnavailableError` is raised before anything is written.
- **Input:** `vcf_path`, `output_path`, `iupac`, `script_path`, `timeout_seconds`.
- **Output:** The common result dictionary with `records` (data lines).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants import vcf_to_tab

out = Path(tempfile.mkdtemp())
table = vcf_to_tab("toy.vcf", str(out / "t.tab"))
lines = Path(table["output_path"]).read_text().splitlines()
print(table["status"], table["records"])
print(lines[0])
print(lines[1])
iupac = vcf_to_tab("toy.vcf", str(out / "t2.tab"), iupac=True)
print(iupac["status"], "IUPAC" in iupac["log_tail"])
```

- **Example output:**

```text
completed 20
#CHROM	POS	REF	S1	S2	S3	S4	S5	S6
chr1	100	G	T/G	G/T	G/T	G/T	T/G	G/G
failed True
```

## 5) `dnakit.variants.vcf_fix_ploidy()`

- **Purpose:** Set the genotype ploidy of chromosomes X, Y and MT by the sex of the samples, with `vcf-fix-ploidy`.
- **Method:** Run `vcf-fix-ploidy -s <samples_file>` with the VCF on standard input, plus `-a <M|F>` (`assumed_sex`, for samples missing from the list), `-l` (`fix_likelihoods`: add or remove heterozygous likelihoods) and `-p <ploidy_file>` (a ploidy definition replacing the default human one). The default makes male X outside the pseudo-autosomal regions haploid, Y absent in females and haploid in males, and MT haploid.
- **API:** `dnakit.variants.vcf_fix_ploidy(vcf_path: 'str', output_path: 'str', samples_file: 'str', *, assumed_sex: 'str | None' = None, fix_likelihoods: 'bool' = False, ploidy_file: 'str | None' = None, script_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the VCFtools Perl script of this section (VCFtools 0.1.16 or later), perl, and `tabix`/`samtools` where the text says so. DNAKit looks for the script at `script_path`, then on PATH, then next to the running Python interpreter; if it is missing, `BackendUnavailableError` is raised before anything is written.
- **Input:** `vcf_path`, `output_path`, `samples_file` (`sample_name M|F` per line), `assumed_sex`, `fix_likelihoods`, `ploidy_file`, `script_path`, `timeout_seconds`.
- **Output:** The common result dictionary with `records`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants import vcf_fix_ploidy

lines = Path("toy.vcf").read_text().splitlines()
x_lines = [x.replace("ID=chr1", "ID=X") if x[0] == "#" else x.replace("chr1", "X", 1) for x in lines]
Path("x.vcf").write_text("\n".join(x_lines) + "\n")
Path("sex.txt").write_text("S1 M\nS2 F\nS3 M\nS4 F\nS5 M\nS6 F\n")
fixed = vcf_fix_ploidy("x.vcf", str(Path(tempfile.mkdtemp()) / "fixed.vcf"), "sex.txt")
row = Path(fixed["output_path"]).read_text().splitlines()[-20].split("\t")
print(fixed["status"], fixed["records"])
print([cell.split(":")[0] for cell in row[9:15]])
```

- **Example output:**

```text
completed 20
['1', '0|1', '0', '0|1', '1', '0|0']
```

## 6) `dnakit.variants.vcf_shuffle_columns()`

- **Purpose:** Reorder the sample columns of a VCF to the order of a template VCF, with `vcf-shuffle-cols`.
- **Method:** Run `vcf-shuffle-cols -t <template_vcf_path> <vcf_path>`. Both files hold the same samples and are normally bgzip-compressed and tabix-indexed.
- **API:** `dnakit.variants.vcf_shuffle_columns(vcf_path: 'str', template_vcf_path: 'str', output_path: 'str', *, script_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the VCFtools Perl script of this section (VCFtools 0.1.16 or later), perl, and `tabix`/`samtools` where the text says so. DNAKit looks for the script at `script_path`, then on PATH, then next to the running Python interpreter; if it is missing, `BackendUnavailableError` is raised before anything is written.
- **Input:** `vcf_path`, `template_vcf_path`, `output_path`, `script_path`, `timeout_seconds`.
- **Output:** The common result dictionary with `records`.
- **Example code:**

```python
import shutil
import tempfile
from pathlib import Path

import pysam
from dnakit.variants import vcf_shuffle_columns

lines = Path("toy.vcf").read_text().splitlines()
columns = [x for x in lines if x.startswith("#CHROM")][0].split("\t")
order = [*range(9), *range(len(columns) - 1, 8, -1)]
template = [
    "\t".join(x.split("\t")[i] for i in order) if x[:2] != "##" else x for x in lines
]
Path("template.vcf").write_text("\n".join(template) + "\n")
shutil.copy("toy.vcf", "toy_c.vcf")
for name in ("template.vcf", "toy_c.vcf"):
    pysam.tabix_index(name, preset="vcf", force=True)
result = vcf_shuffle_columns("toy_c.vcf.gz", "template.vcf.gz", str(Path(tempfile.mkdtemp()) / "shuffled.vcf"))
header = [x for x in Path(result["output_path"]).read_text().splitlines() if x.startswith("#CHROM")][0]
print(result["status"], result["records"], header.split("\t")[9:])
```

- **Example output:**

```text
completed 20 ['S6', 'S5', 'S4', 'S3', 'S2', 'S1']
```

## 7) `dnakit.variants.vcf_contrast()`

- **Purpose:** Find genotypes and alleles that occur only in some samples, with `vcf-contrast`.
- **Method:** Run `vcf-contrast +<variant samples> -<background samples> [-d <min_depth>] [-f] [-n] <vcf>`. NOVELGT, NOVELAL and NOVELTY are added to INFO for sites where a variant sample has a genotype or allele absent from all background samples. `min_depth` is the minimum depth across the background samples (`-d`), `apply_filters` skips sites whose FILTER is not PASS or "." (`-f`), `novel_sites_only` prints only records with novel genotypes (`-n`). Haploid genotypes are treated as homozygous diploid.
- **API:** `dnakit.variants.vcf_contrast(vcf_path: 'str', output_path: 'str', variant_samples: 'Sequence[str]', background_samples: 'Sequence[str]', *, min_depth: 'int | None' = None, apply_filters: 'bool' = False, novel_sites_only: 'bool' = False, script_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the VCFtools Perl script of this section (VCFtools 0.1.16 or later), perl, and `tabix`/`samtools` where the text says so. DNAKit looks for the script at `script_path`, then on PATH, then next to the running Python interpreter; if it is missing, `BackendUnavailableError` is raised before anything is written.
- **Input:** `vcf_path` (VCF or bgzip-compressed VCF), `output_path`, `variant_samples`, `background_samples`, `min_depth`, `apply_filters`, `novel_sites_only`, `script_path`, `timeout_seconds`.
- **Output:** The common result dictionary with `records`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants import vcf_contrast

result = vcf_contrast(
    "toy.vcf", str(Path(tempfile.mkdtemp()) / "novel.vcf"), ["S1"], ["S2", "S3"], novel_sites_only=True
)
print(result["status"], result["records"])
for line in Path(result["output_path"]).read_text().splitlines():
    if line[0] != "#":
        fields = line.split("\t")
        print(fields[0], fields[1], sorted(fields[7].split(";")))
```

- **Example output:**

```text
completed 5
chr1 800 ['AF=0.25', 'DP=120', 'NOVELGT=S1', 'NOVELTY=30']
chr1 1500 ['AF=0.5', 'DP=99', 'NOVELGT=S1', 'NOVELTY=30']
chr1 2900 ['AF=0.25', 'DP=172', 'NOVELAL=S1', 'NOVELTY=30']
chr2 300 ['AF=0.4', 'DP=204', 'NOVELGT=S1', 'NOVELTY=30']
chr2 3900 ['AF=0.4', 'DP=117', 'NOVELGT=S1', 'NOVELTY=30']
```

## 8) `dnakit.variants.vcf_compare()`

- **Purpose:** Compare the sites and genotypes of two or more VCF files, with `vcf-compare`.
- **Method:** Run `vcf-compare` with `-g` (`compare_genotypes`), `--ignore-indels` (`ignore_indels`), `-a` (`apply_filters`: ignore sites whose FILTER is not PASS or "."), `-r` (`regions`), `-s` (`samples`), `-w` (`window`: treat records this many bp apart as matching) and `-R` (`reference_fasta`: compare the sequence, not only the positions). The files must be bgzip-compressed and tabix-indexed. The report is written to `output_path` and parsed into `venn` (sites per file combination), `summary` (the `SN` lines) and `genotype_classes` (the `GS` lines).
- **API:** `dnakit.variants.vcf_compare(vcf_paths: 'Sequence[str]', output_path: 'str', *, compare_genotypes: 'bool' = False, apply_filters: 'bool' = False, ignore_indels: 'bool' = False, regions: 'Sequence[str] | None' = None, samples: 'Sequence[str] | None' = None, window: 'int | None' = None, reference_fasta: 'str | None' = None, script_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the VCFtools Perl script of this section (VCFtools 0.1.16 or later), perl, and `tabix`/`samtools` where the text says so. DNAKit looks for the script at `script_path`, then on PATH, then next to the running Python interpreter; if it is missing, `BackendUnavailableError` is raised before anything is written.
- **Input:** `vcf_paths` (at least two), `output_path`, `compare_genotypes`, `apply_filters`, `ignore_indels`, `regions`, `samples`, `window`, `reference_fasta`, `script_path`, `timeout_seconds`.
- **Output:** The common result dictionary with `venn`, `summary` and `genotype_classes`; the full report is at `output_path`.
- **Example code:**

```python
import tempfile
from pathlib import Path

import pysam
from dnakit.variants import vcf_compare

lines = Path("toy.vcf").read_text().splitlines()
head = [x for x in lines if x[0] == "#"]
body = [x for x in lines if x[0] != "#"]
Path("all.vcf").write_text("\n".join(head + body) + "\n")
Path("half.vcf").write_text("\n".join(head + body[::2]) + "\n")
for name in ("all.vcf", "half.vcf"):
    pysam.tabix_index(name, preset="vcf", force=True)
result = vcf_compare(["all.vcf.gz", "half.vcf.gz"], str(Path(tempfile.mkdtemp()) / "cmp.txt"), compare_genotypes=True)
print(result["status"])
for row in sorted(result["venn"], key=lambda r: len(r["files"])):
    print(row["sites"], [Path(name.split(" ")[0]).name for name in row["files"]])
print(result["summary"]["Number of REF matches"], result["summary"]["Number of samples in GT comparison"])
```

- **Example output:**

```text
completed
10 ['all.vcf.gz']
10 ['all.vcf.gz', 'half.vcf.gz']
10 6
```

## 9) `dnakit.variants.vcf_indel_stats()`

- **Purpose:** Count the in-frame and frameshift indels of a VCF given exon coordinates, with `vcf-indel-stats`.
- **Method:** Run `vcf-indel-stats -e <exons_file>` with the VCF on standard input. The script classifies the indels that fall in an exon by whether their length is a multiple of three and reports the in-frame ratio.
- **API:** `dnakit.variants.vcf_indel_stats(vcf_path: 'str', output_path: 'str', exons_file: 'str', *, script_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the VCFtools Perl script of this section (VCFtools 0.1.16 or later), perl, and `tabix`/`samtools` where the text says so. DNAKit looks for the script at `script_path`, then on PATH, then next to the running Python interpreter; if it is missing, `BackendUnavailableError` is raised before anything is written.
- **Input:** `vcf_path`, `output_path`, `exons_file` (tab-separated CHROM, FROM, TO; 1-based, inclusive), `script_path`, `timeout_seconds`.
- **Output:** The common result dictionary with `statistics` (the lines of the report as a dictionary).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants import vcf_indel_stats

Path("exons.txt").write_text("chr1\t2000\t3000\nchr1\t5600\t6000\n")
result = vcf_indel_stats("toy.vcf", str(Path(tempfile.mkdtemp()) / "stats.txt"), "exons.txt")
print(result["status"], result["statistics"])
```

- **Example output:**

```text
completed {'total': 3, 'in-frame': 0, 'frameshift': 2, 'ratio': 1.0}
```

## 10) `dnakit.variants.vcf_phased_join()`

- **Purpose:** Join overlapping pre-phased VCF chunks into one VCF, with `vcf-phased-join`.
- **Method:** Run `vcf-phased-join -o <output_path> [-j <min_join_quality>] [-q <min_phase_quality>] <chunks>`. Heterozygous calls in the overlaps decide the phase of each chunk. `min_join_quality` is the quality for gluing blocks together (script default 10) and `min_phase_quality` breaks pre-phased segments whose PQ value is lower (script default 0.6). The chunks are given in genomic order and are normally bgzip-compressed and tabix-indexed.
- **API:** `dnakit.variants.vcf_phased_join(vcf_paths: 'Sequence[str]', output_path: 'str', *, min_join_quality: 'float | None' = None, min_phase_quality: 'float | None' = None, script_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the VCFtools Perl script of this section (VCFtools 0.1.16 or later), perl, and `tabix`/`samtools` where the text says so. DNAKit looks for the script at `script_path`, then on PATH, then next to the running Python interpreter; if it is missing, `BackendUnavailableError` is raised before anything is written.
- **Input:** `vcf_paths` (at least two chunks), `output_path`, `min_join_quality`, `min_phase_quality`, `script_path`, `timeout_seconds`.
- **Output:** The common result dictionary with `records` of the joined VCF.
- **Example code:**

```python
import tempfile
from pathlib import Path

import pysam
from dnakit.variants import vcf_phased_join

lines = Path("toy.vcf").read_text().splitlines()
head = [x for x in lines if x[0] == "#"]
body = [x for x in lines if x[0] != "#"]
Path("chunk1.vcf").write_text("\n".join(head + body[:12]) + "\n")
Path("chunk2.vcf").write_text("\n".join(head + body[7:14]) + "\n")
for name in ("chunk1.vcf", "chunk2.vcf"):
    pysam.tabix_index(name, preset="vcf", force=True)
result = vcf_phased_join(["chunk1.vcf.gz", "chunk2.vcf.gz"], str(Path(tempfile.mkdtemp()) / "joined.vcf"))
print(result["status"], result["records"])
```

- **Example output:**

```text
completed 14
```

## 11) `dnakit.variants.vcf_fill_ancestral_allele()`

- **Purpose:** Add the ancestral allele (INFO/AA) from per-chromosome ancestral FASTA files, with `fill-aa`.
- **Method:** Run `fill-aa -a <ancestral_prefix> -t <variant_types>` with the VCF on standard input. The files are named `<ancestral_prefix><chromosome>.fa.gz` and must be compressed and indexed with `samtools faidx`. `variant_types` selects among `all`, `indel`, `ref` and `snp`.
- **API:** `dnakit.variants.vcf_fill_ancestral_allele(vcf_path: 'str', output_path: 'str', ancestral_prefix: 'str', *, variant_types: 'Sequence[str]' = ('all',), script_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the VCFtools Perl script of this section (VCFtools 0.1.16 or later), perl, and `tabix`/`samtools` where the text says so. DNAKit looks for the script at `script_path`, then on PATH, then next to the running Python interpreter; if it is missing, `BackendUnavailableError` is raised before anything is written.
- **Input:** `vcf_path`, `output_path`, `ancestral_prefix`, `variant_types` (default `("all",)`), `script_path`, `timeout_seconds`. The script needs `samtools` on PATH.
- **Output:** The common result dictionary with `records`; `log_tail` holds the script's summary of filled sites.
- **Example code:**

```python
import tempfile
from pathlib import Path

import pysam
from dnakit.variants import vcf_fill_ancestral_allele

sequences = {"chr1": ["A"] * 20000, "chr2": ["A"] * 20000}
for line in Path("toy.vcf").read_text().splitlines():
    if not line.startswith("#"):
        fields = line.split("\t")
        for offset, base in enumerate(fields[3]):
            sequences[fields[0]][int(fields[1]) - 1 + offset] = base
for name, bases in sequences.items():
    Path(f"anc_{name}.fa").write_text(f">{name}\n" + "".join(bases) + "\n")
    pysam.tabix_compress(f"anc_{name}.fa", f"anc_{name}.fa.gz", force=True)
    pysam.faidx(f"anc_{name}.fa.gz")
result = vcf_fill_ancestral_allele("toy.vcf", str(Path(tempfile.mkdtemp()) / "aa.vcf"), "anc_")
print(result["status"], result["records"])
print(Path(result["output_path"]).read_text().splitlines()[-1].split("\t")[7])
```

- **Example output:**

```text
completed 20
DP=154;AF=0.1;AA=T
```

## 12) `dnakit.variants.vcf_fill_flanking_sequence()`

- **Purpose:** Add the flanking reference sequence (INFO/FS) around every variant, with `fill-fs`.
- **Method:** Run `fill-fs -r <reference_fasta> -l <length> [-m <mask_character> -b <bed mask>] [-m <mask_character> -v <vcf mask>] [-m <mask_character> -c <cluster_distance>] <vcf>`. The reference is a FASTA indexed with `samtools faidx`. Flanks are `length` bases on each side (script default 100); `bed_mask_files` and `vcf_mask_files` (tabix-indexed) mask regions and known variants in the flank with `mask_character` (one character or `lc` for lower case, default `N`), and `cluster_distance` self-masks clustered variants. The script requires `-m` to precede each mask option, so DNAKit repeats it before every mask; it has no effect without masks.
- **API:** `dnakit.variants.vcf_fill_flanking_sequence(vcf_path: 'str', output_path: 'str', reference_fasta: 'str', *, length: 'int' = 100, mask_character: 'str' = 'N', bed_mask_files: 'Sequence[str] | None' = None, vcf_mask_files: 'Sequence[str] | None' = None, cluster_distance: 'int | None' = None, script_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the VCFtools Perl script of this section (VCFtools 0.1.16 or later), perl, and `tabix`/`samtools` where the text says so. DNAKit looks for the script at `script_path`, then on PATH, then next to the running Python interpreter; if it is missing, `BackendUnavailableError` is raised before anything is written.
- **Input:** `vcf_path`, `output_path`, `reference_fasta`, `length`, `mask_character`, `bed_mask_files`, `vcf_mask_files`, `cluster_distance`, `script_path`, `timeout_seconds`.
- **Output:** The common result dictionary with `records`; INFO/FS has the form `LEFT[REF/ALT]RIGHT`.
- **Example code:**

```python
import tempfile
from pathlib import Path

import pysam
from dnakit.variants import vcf_fill_flanking_sequence

sequences = {"chr1": ["A"] * 20000, "chr2": ["A"] * 20000}
for line in Path("toy.vcf").read_text().splitlines():
    if not line.startswith("#"):
        fields = line.split("\t")
        for offset, base in enumerate(fields[3]):
            sequences[fields[0]][int(fields[1]) - 1 + offset] = base
Path("ref.fa").write_text("".join(f">{n}\n{''.join(b)}\n" for n, b in sequences.items()))
pysam.faidx("ref.fa")
result = vcf_fill_flanking_sequence("toy.vcf", str(Path(tempfile.mkdtemp()) / "fs.vcf"), "ref.fa", length=5)
rows = [x.split("\t") for x in Path(result["output_path"]).read_text().splitlines() if x[0] != "#"]
print(result["status"], result["records"])
print(rows[0][7])
print(rows[1][7])
```

- **Example output:**

```text
completed 20
DP=91;AF=0.1;FS=AAAAA[G/T]AAAAA
DP=120;AF=0.25;FS=AAAAA[G/A]AAAAA
```

## 13) `dnakit.variants.vcf_fill_reference_md5()`

- **Purpose:** Insert the MD5 sums of the reference contigs into the header of a VCF, with `fill-ref-md5`.
- **Method:** Run `fill-ref-md5 -r <reference_fasta> [-d <dictionary_path>] [-i AS:..,SP:..,TX:..] <vcf.gz> <out.vcf.gz>`. The script adds `reference` and `contig` header lines as recommended by VCF 4.1, using `tabix -r` to replace the header. The input must be bgzip-compressed and tabix-indexed. `dictionary_path` caches computed MD5 sums and is opened in append mode, so it may already exist. The script ignores a failing `tabix`, so DNAKit reports an empty output as `failed`.
- **API:** `dnakit.variants.vcf_fill_reference_md5(vcf_path: 'str', output_path: 'str', reference_fasta: 'str', *, dictionary_path: 'str | None' = None, assembly: 'str | None' = None, species: 'str | None' = None, taxonomy: 'str | None' = None, script_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the VCFtools Perl script of this section (VCFtools 0.1.16 or later), perl, and `tabix`/`samtools` where the text says so. DNAKit looks for the script at `script_path`, then on PATH, then next to the running Python interpreter; if it is missing, `BackendUnavailableError` is raised before anything is written.
- **Input:** `vcf_path` (`.vcf.gz`), `output_path` (a new `.vcf.gz`), `reference_fasta`, `dictionary_path`, `assembly`, `species`, `taxonomy` (the AS, SP and TX header values; no commas), `script_path`, `timeout_seconds`. Needs `tabix` and `samtools` on PATH.
- **Output:** The common result dictionary; the new compressed VCF is at `output_path`.
- **Example code:**

```python
import gzip
import shutil
import tempfile
from pathlib import Path

import pysam
from dnakit.variants import vcf_fill_reference_md5

sequences = {"chr1": ["A"] * 20000, "chr2": ["A"] * 20000}
for line in Path("toy.vcf").read_text().splitlines():
    if not line.startswith("#"):
        fields = line.split("\t")
        for offset, base in enumerate(fields[3]):
            sequences[fields[0]][int(fields[1]) - 1 + offset] = base
Path("ref.fa").write_text("".join(f">{n}\n{''.join(b)}\n" for n, b in sequences.items()))
pysam.faidx("ref.fa")
shutil.copy("toy.vcf", "toy_md5.vcf")
pysam.tabix_index("toy_md5.vcf", preset="vcf", force=True)
result = vcf_fill_reference_md5("toy_md5.vcf.gz", str(Path(tempfile.mkdtemp()) / "md5.vcf.gz"), "ref.fa", assembly="toy")
print(result["status"])
for line in gzip.open(result["output_path"], "rt").read().splitlines():
    if line.startswith("##contig"):
        print(line)
```

- **Example output:**

```text
completed
##contig=<ID=chr1,assembly=toy,length=20000,md5=93e08f94707a244b4a8ccfd891707c3d>
##contig=<ID=chr2,assembly=toy,length=20000,md5=6a1cea00ce86defc019f8e80d5492f86>
```
