# bcftools plugin facades: VCF editing, splitting and trio analysis

This group of functions in `dnakit.variants` (modules `bcftools_plugins_p06_cli_functions` and `bcftools_trio_p06_cli_functions`) wraps plugins of the external program [bcftools](https://github.com/samtools/bcftools) (`bcftools +plugin`, 1.24 is the tested version) that change or split VCF/BCF files, plus the trio plugins: Mendelian consistency, de novo mutation screening, transmission statistics, phase switch rate and the parental origin of a copy-number variant. Each function validates its arguments, runs one plugin with a fixed argument list and no shell, captures standard output and standard error separately, and returns a JSON-serialisable dictionary. Record-changing plugins write a VCF/BCF to `output_path` (`output_type` `v`, `z`, `b` or `u`, default from the file name); plugins that write many files use a new folder; report plugins save the report and parse it into `report`. An existing output is never overwritten and a failed run leaves nothing behind.

bcftools is an external program with its plugins: install it (for example `conda install -c bioconda bcftools`, or a release from https://www.htslib.org/download/), then put it on PATH or pass `bcftools_path`. Plugins are found through `BCFTOOLS_PLUGINS` or the `libexec/bcftools` folder next to the installation. If the program cannot be found, `BackendUnavailableError` is raised with the install hint; if it exits with an error, `BackendExecutionError` (code `HTSLIB_EXECUTION_FAILED`) is raised with the tail of its log. A trio is given either as a PED file (`ped_path`) or as a list of sample names (`trio`). `bcftools +trio-dnm3` model `dmm` needs bcftools 1.24 or later. The common options `regions` (`-r`), `include` (`-i`) and `exclude` (`-e`) are available where bcftools accepts them. Results share `status`, `backend`, `operation`, `version`, `command`, `outputs`, `log_tail` and `provenance`; the sections below list only the extra fields.

## Example data

The examples read these files from the working directory and write their outputs to a temporary directory. `trio.vcf` has a child C, a father F, a mother M and an unrelated sample U (fourteen sites cycling through seven genotype patterns, two of which break Mendelian inheritance); `trio_phased.vcf` is the same with phased genotypes; `trio.ped` names the trio; `csq.vcf` carries a VEP-style `CSQ` tag.

**Example file:** `genome.fa`

```text
>chr1
GCTAAAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCG
CTTAAGGGTTAAGTAAGTGTGATGCATACGCCTTTACTTGCTGTGTCCACCCCATCGGAC
TGGCATTTTTATTACACTCAGAAACAGAACTCGGGTAATTTTGACAGGTCACGCAGAGGC
GCGCCCTCCTGAAGTGCGTGGACACTCGCTATGAATCTCTGATTTACCCACTCTGCCAAA
CTCCAGCGCGGTCAGTTCCATCACCCTAAGTAACCGAATAATGCGTTCGCTCTATTGACT
ACGACGCGCTCATTCCCTTGTCGGAGAGTTATGGAACAAGGACGCTGTCTGAGACTAGAA
GACAGATAGTGCACACGACCGGCGTCGGAGAAACTCTATTTGCCGCCTGACAAGTCAATG
CGATCCGTAGGGGCAGCGCAGTATGCCAAGACTATAGGCACTGTCGCATCACAAACGATT
AACTGATAAATGAGCCCTTTATGACACGGGCATATGACTGGTTTACGATAGTATGTCCAA
CGGCGAGCTTTACATTTGCTGTGAGAGGTACAGGGATTAGTGAGAAGCCGTGCGTATCAA
>chr2
TTCGTACCTTGGGGGTCGTTACCACTCTGTTCCCACGAGCGGCATTTCTGGATGGCCAGC
TTTTGACATTTAATTTCACCCATAAACCAGCGTAAAGCTGCAAGTGGCTCCATGAACTTA
GCTGCTAGTGTCAGACTCGCCTCGGATCCTTACTACACTAACTTGAACGCCTAGTGGTCA
AAGAGTACTGGTAATCGTCGGTATCTATATAAGCAGGGGAGGGGAAACATTTGTTCTCAG
CCGGTGACTCCTAATGCTAAGACATTTCCCTTCAGGGGGGGCTCCCCCGCGATGCCATAA
ATCTGAGCAACCAGCTGAAGCAGGCACGACAGTGCGACATTATATCACTGTGGTAGGTTA
GCTTCATCTAATGTCCAACTAGCCGGCCAATTCGCATGAT
```

**Example file:** `variants.vcf`

```text
##fileformat=VCFv4.2
##contig=<ID=chr1,length=600>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4	S5	S6
chr1	60	var1	G	T	40	PASS	DP=50;AF=0.417	GT	0/1	0/0	1/1	0/1	0/0	0/1
chr1	120	var2	C	T	41	PASS	DP=55;AF=0.333	GT	0/0	0/1	0/1	1/1	0/0	0/0
chr1	180	var3	C	A	42	PASS	DP=60;AF=0.5	GT	1/1	0/1	0/0	0/1	0/1	./.
chr1	240	var4	ACT	A	43	PASS	DP=65;AF=0.417	GT	0/1	0/1	0/0	0/0	1/1	0/1
chr1	300	var5	T	C	44	PASS	DP=70;AF=0.417	GT	0/0	0/0	0/1	0/1	0/1	1/1
chr1	360	var6	A	T	45	PASS	DP=75;AF=0.417	GT	0/1	1/1	0/1	0/0	0/0	0/1
chr1	420	var7	G	GCG	46	PASS	DP=80;AF=0.417	GT	0/0	0/1	0/0	0/1	1/1	0/1
chr1	480	var8	T	C	47	PASS	DP=85;AF=0.5	GT	0/1	0/0	1/1	1/1	0/1	0/0
chr1	520	var9	G	C	48	PASS	DP=90;AF=0.25	GT	0/0	0/0	0/1	0/0	0/1	0/1
chr1	570	var10	A	C	49	PASS	DP=95;AF=0.417	GT	1/1	0/1	0/1	0/0	0/0	0/1
```

**Example file:** `trio.vcf`

```text
##fileformat=VCFv4.2
##contig=<ID=chr1,length=600>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred likelihoods">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depth">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=QS,Number=R,Type=Integer,Description="Phred-scaled quality sum">
##FORMAT=<ID=SP,Number=1,Type=Integer,Description="Strand bias">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	C	F	M	U
chr1	40	.	T	A	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	0/1:60,0,60:15,15:30:60:570,570:0	0/1:60,0,60:15,15:30:60:570,570:0	0/0:0,60,255:30,0:30:60:1140,0:0	0/0:0,60,255:30,0:30:60:1140,0:0
chr1	80	.	T	A	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	1/1:255,60,0:0,30:30:60:0,1140:0	0/1:60,0,60:15,15:30:60:570,570:0	0/1:60,0,60:15,15:30:60:570,570:0	0/1:60,0,60:15,15:30:60:570,570:0
chr1	120	.	C	G	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	0/1:60,0,60:15,15:30:60:570,570:0	1/1:255,60,0:0,30:30:60:0,1140:0	0/0:0,60,255:30,0:30:60:1140,0:0	0/0:0,60,255:30,0:30:60:1140,0:0
chr1	160	.	T	A	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	0/0:0,60,255:30,0:30:60:1140,0:0	0/0:0,60,255:30,0:30:60:1140,0:0	0/0:0,60,255:30,0:30:60:1140,0:0	0/1:60,0,60:15,15:30:60:570,570:0
chr1	200	.	G	T	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	0/1:60,0,60:15,15:30:60:570,570:0	0/0:0,60,255:30,0:30:60:1140,0:0	0/0:0,60,255:30,0:30:60:1140,0:0	0/0:0,60,255:30,0:30:60:1140,0:0
chr1	240	.	A	C	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	0/0:0,60,255:30,0:30:60:1140,0:0	0/1:60,0,60:15,15:30:60:570,570:0	0/1:60,0,60:15,15:30:60:570,570:0	1/1:255,60,0:0,30:30:60:0,1140:0
chr1	280	.	A	C	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	1/1:255,60,0:0,30:30:60:0,1140:0	0/0:0,60,255:30,0:30:60:1140,0:0	0/1:60,0,60:15,15:30:60:570,570:0	0/0:0,60,255:30,0:30:60:1140,0:0
chr1	320	.	G	T	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	0/1:60,0,60:15,15:30:60:570,570:0	0/1:60,0,60:15,15:30:60:570,570:0	0/0:0,60,255:30,0:30:60:1140,0:0	0/0:0,60,255:30,0:30:60:1140,0:0
chr1	360	.	A	C	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	1/1:255,60,0:0,30:30:60:0,1140:0	0/1:60,0,60:15,15:30:60:570,570:0	0/1:60,0,60:15,15:30:60:570,570:0	0/1:60,0,60:15,15:30:60:570,570:0
chr1	400	.	T	A	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	0/1:60,0,60:15,15:30:60:570,570:0	1/1:255,60,0:0,30:30:60:0,1140:0	0/0:0,60,255:30,0:30:60:1140,0:0	0/0:0,60,255:30,0:30:60:1140,0:0
chr1	440	.	A	C	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	0/0:0,60,255:30,0:30:60:1140,0:0	0/0:0,60,255:30,0:30:60:1140,0:0	0/0:0,60,255:30,0:30:60:1140,0:0	0/1:60,0,60:15,15:30:60:570,570:0
chr1	480	.	T	A	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	0/1:60,0,60:15,15:30:60:570,570:0	0/0:0,60,255:30,0:30:60:1140,0:0	0/0:0,60,255:30,0:30:60:1140,0:0	0/0:0,60,255:30,0:30:60:1140,0:0
chr1	520	.	G	T	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	0/0:0,60,255:30,0:30:60:1140,0:0	0/1:60,0,60:15,15:30:60:570,570:0	0/1:60,0,60:15,15:30:60:570,570:0	1/1:255,60,0:0,30:30:60:0,1140:0
chr1	560	.	T	A	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	1/1:255,60,0:0,30:30:60:0,1140:0	0/0:0,60,255:30,0:30:60:1140,0:0	0/1:60,0,60:15,15:30:60:570,570:0	0/0:0,60,255:30,0:30:60:1140,0:0
```

**Example file:** `trio_phased.vcf`

```text
##fileformat=VCFv4.2
##contig=<ID=chr1,length=600>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred likelihoods">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depth">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=QS,Number=R,Type=Integer,Description="Phred-scaled quality sum">
##FORMAT=<ID=SP,Number=1,Type=Integer,Description="Strand bias">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	C	F	M	U
chr1	40	.	T	A	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	0|1:60,0,60:15,15:30:60:570,570:0	0|1:60,0,60:15,15:30:60:570,570:0	0|0:0,60,255:30,0:30:60:1140,0:0	0|0:0,60,255:30,0:30:60:1140,0:0
chr1	80	.	T	A	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	1|1:255,60,0:0,30:30:60:0,1140:0	0|1:60,0,60:15,15:30:60:570,570:0	0|1:60,0,60:15,15:30:60:570,570:0	0|1:60,0,60:15,15:30:60:570,570:0
chr1	120	.	C	G	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	0|1:60,0,60:15,15:30:60:570,570:0	1|1:255,60,0:0,30:30:60:0,1140:0	0|0:0,60,255:30,0:30:60:1140,0:0	0|0:0,60,255:30,0:30:60:1140,0:0
chr1	160	.	T	A	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	0|0:0,60,255:30,0:30:60:1140,0:0	0|0:0,60,255:30,0:30:60:1140,0:0	0|0:0,60,255:30,0:30:60:1140,0:0	0|1:60,0,60:15,15:30:60:570,570:0
chr1	200	.	G	T	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	0|1:60,0,60:15,15:30:60:570,570:0	0|0:0,60,255:30,0:30:60:1140,0:0	0|0:0,60,255:30,0:30:60:1140,0:0	0|0:0,60,255:30,0:30:60:1140,0:0
chr1	240	.	A	C	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	0|0:0,60,255:30,0:30:60:1140,0:0	0|1:60,0,60:15,15:30:60:570,570:0	0|1:60,0,60:15,15:30:60:570,570:0	1|1:255,60,0:0,30:30:60:0,1140:0
chr1	280	.	A	C	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	1|1:255,60,0:0,30:30:60:0,1140:0	0|0:0,60,255:30,0:30:60:1140,0:0	0|1:60,0,60:15,15:30:60:570,570:0	0|0:0,60,255:30,0:30:60:1140,0:0
chr1	320	.	G	T	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	0|1:60,0,60:15,15:30:60:570,570:0	0|1:60,0,60:15,15:30:60:570,570:0	0|0:0,60,255:30,0:30:60:1140,0:0	0|0:0,60,255:30,0:30:60:1140,0:0
chr1	360	.	A	C	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	1|1:255,60,0:0,30:30:60:0,1140:0	0|1:60,0,60:15,15:30:60:570,570:0	0|1:60,0,60:15,15:30:60:570,570:0	0|1:60,0,60:15,15:30:60:570,570:0
chr1	400	.	T	A	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	0|1:60,0,60:15,15:30:60:570,570:0	1|1:255,60,0:0,30:30:60:0,1140:0	0|0:0,60,255:30,0:30:60:1140,0:0	0|0:0,60,255:30,0:30:60:1140,0:0
chr1	440	.	A	C	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	0|0:0,60,255:30,0:30:60:1140,0:0	0|0:0,60,255:30,0:30:60:1140,0:0	0|0:0,60,255:30,0:30:60:1140,0:0	0|1:60,0,60:15,15:30:60:570,570:0
chr1	480	.	T	A	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	0|1:60,0,60:15,15:30:60:570,570:0	0|0:0,60,255:30,0:30:60:1140,0:0	0|0:0,60,255:30,0:30:60:1140,0:0	0|0:0,60,255:30,0:30:60:1140,0:0
chr1	520	.	G	T	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	0|0:0,60,255:30,0:30:60:1140,0:0	0|1:60,0,60:15,15:30:60:570,570:0	0|1:60,0,60:15,15:30:60:570,570:0	1|1:255,60,0:0,30:30:60:0,1140:0
chr1	560	.	T	A	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	1|1:255,60,0:0,30:30:60:0,1140:0	0|0:0,60,255:30,0:30:60:1140,0:0	0|1:60,0,60:15,15:30:60:570,570:0	0|0:0,60,255:30,0:30:60:1140,0:0
```

**Example file:** `csq.vcf`

```text
##fileformat=VCFv4.2
##contig=<ID=chr1,length=600>
##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Feature">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	10	.	A	G	.	.	CSQ=G|missense_variant|MODERATE|GENE1|T1,G|synonymous_variant|LOW|GENE1|T2
chr1	20	.	C	T	.	.	CSQ=T|stop_gained|HIGH|GENE2|T3
chr1	30	.	G	A	.	.	.
```

**Example file:** `trio.ped`

```text
fam1	C	F	M	1	0
```

## 1) `dnakit.variants.bcftools_prune_sites()`

- **Purpose:** Prune or annotate sites by linkage disequilibrium (r2, Lewontin D, Ragsdale D) or by site density in windows (`bcftools +prune`).
- **Method:** `statistic` is `r2`, `LD`, `RD` or `count`; `max_value` removes sites above it (`-m`) or, with `set_filter` (`-f`), marks them with that FILTER. `annotate` adds the biggest `r2`/`LD`/`RD` value and its position, or the `count` (`-a`). `window` is `100kb`, `1000` (sites) and similar (`-w`). `sites_per_window` keeps at most N sites per window (`-n`), chosen by `sites_per_window_mode` `maxAF`, `1st` or `rand`. `af_tag`, `keep_sites`, `random_seed` and `randomize_missing` are the matching options. One of `max_value`, `sites_per_window` or `annotate` is required.
- **API:** `dnakit.variants.bcftools_prune_sites(vcf_path, output_path, *, statistic='r2', max_value=None, annotate=None, window='100kb', sites_per_window=None, sites_per_window_mode='maxAF', af_tag=None, set_filter=None, keep_sites=False, random_seed=None, randomize_missing=False, regions=None, include=None, exclude=None, output_type=None, bcftools_path=None)`.
- **Input:** A VCF/BCF with `FORMAT/GT` and a new output path.
- **Output:** The common fields, `output_type` and `records` (plain VCF output only).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_plugins_p06_cli_functions import bcftools_prune_sites

work = Path(tempfile.mkdtemp())
pruned = bcftools_prune_sites("variants.vcf", work / "pruned.vcf", max_value=0.3, window="1000")
thinned = bcftools_prune_sites("variants.vcf", work / "thin.vcf", sites_per_window=1, window="150bp", sites_per_window_mode="1st")
annotated = bcftools_prune_sites("variants.vcf", work / "r2.vcf", annotate=["r2"], window="1", set_filter=".")
print(pruned["records"], thinned["records"], annotated["records"])
print(annotated["outputs"]["output"].endswith("r2.vcf"), "R2=" in (work / "r2.vcf").read_text())
```

- **Example output:**

```text
4 4 10
True True
```

## 2) `dnakit.variants.bcftools_remove_overlaps()`

- **Purpose:** Remove, list or mark overlapping and duplicate variants (`bcftools +remove-overlaps`).
- **Method:** `mark` chooses the sites (`-m`): `overlap`, `dup` (all duplicate sites) or `min(QUAL)` (the lowest-QUAL sites until the overlaps are resolved; `missing_value` `0` or `DP`). The chosen sites are removed, or with `mark_tag` kept and flagged `INFO/TAG` (`-M`). `reverse` applies the reverse logic. `output_type` `t` or `tz` writes a plain (compressed) list of chromosome and position instead of a VCF.
- **API:** `dnakit.variants.bcftools_remove_overlaps(vcf_path, output_path, *, mark='overlap', mark_tag=None, missing_value='0', reverse=False, output_type=None, regions=None, bcftools_path=None)`.
- **Input:** A VCF/BCF and a new output path.
- **Output:** The common fields, `output_type`, `records` (plain VCF output only) and `messages`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_plugins_p06_cli_functions import bcftools_remove_overlaps

work = Path(tempfile.mkdtemp())
dup = Path("variants.vcf").read_text().splitlines()
row = dup[-1].split("\t")
row[2] = "copy"
(work / "dup.vcf").write_text("\n".join(dup + ["\t".join(row)]) + "\n")
removed = bcftools_remove_overlaps(work / "dup.vcf", work / "removed.vcf", mark="dup")
marked = bcftools_remove_overlaps(work / "dup.vcf", work / "marked.vcf", mark="dup", mark_tag="DUP")
listed = bcftools_remove_overlaps(work / "dup.vcf", work / "sites.txt", mark="dup", output_type="t")
print(removed["records"], marked["records"], (work / "sites.txt").read_text().split())
```

- **Example output:**

```text
9 11 ['chr1', '60', 'chr1', '120', 'chr1', '180', 'chr1', '240', 'chr1', '300', 'chr1', '360', 'chr1', '420', 'chr1', '480', 'chr1', '520']
```

## 3) `dnakit.variants.bcftools_scatter_vcf()`

- **Purpose:** Split a VCF/BCF into chunks of N sites or into region windows (`bcftools +scatter`).
- **Method:** Exactly one of `sites_per_chunk` (`-n`), `regions_to_scatter` (`-s`, one output per region) and `scatter_file` (`-S`, regions with an optional second column naming the output) is required. `extra_name` writes the records outside the listed regions to a separate file (`-x`), `prefix` is prepended to the output names (`-p`). The files go into the new folder `output_dir`.
- **API:** `dnakit.variants.bcftools_scatter_vcf(vcf_path, output_dir, *, sites_per_chunk=None, regions_to_scatter=None, scatter_file=None, extra_name=None, prefix=None, output_type='v', bcftools_path=None)`.
- **Input:** A VCF/BCF (indexed for region scattering) and a new output folder.
- **Output:** The common fields, `output_type` and `files`, the names of the files written.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_plugins_p06_cli_functions import bcftools_scatter_vcf

result = bcftools_scatter_vcf("variants.vcf", Path(tempfile.mkdtemp()) / "chunks", sites_per_chunk=4, prefix="part-")
print(result["files"])
```

- **Example output:**

```text
['part-0.vcf', 'part-1.vcf', 'part-2.vcf']
```

## 4) `dnakit.variants.bcftools_set_genotypes()`

- **Purpose:** Set genotypes by rule: missing to reference, major or minor allele, phase or unphase, custom genotype (`bcftools +setGT` and `bcftools +missing2ref`).
- **Method:** `target_genotypes` (`-t`): `./.`, `./x`, `.`, `a`, `q` (genotypes for which `genotype_include` / `genotype_exclude` is true, the plugin's `-i` / `-e`), `b` (heterozygotes failing a binomial test, `b:AD<1e-3`) or `r:FLOAT` (a random proportion). `new_genotype` (`-n`): `.`, `0`, `c:GT`, `i`, `m`, `M`, `X`, `p`, `u`, optionally `0p`/`Mp`. `+missing2ref` is `target_genotypes='.'` with `new_genotype` `0` (`0p` for its `-p`) or `M` (`Mp` for `-m -p`); the outputs are identical.
- **API:** `dnakit.variants.bcftools_set_genotypes(vcf_path, output_path, *, target_genotypes='.', new_genotype, genotype_include=None, genotype_exclude=None, random_seed=0, regions=None, include=None, exclude=None, output_type=None, bcftools_path=None)`.
- **Input:** A VCF/BCF with `FORMAT/GT` and a new output path.
- **Output:** The common fields, `output_type` and `records` (plain VCF output only).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_plugins_p06_cli_functions import bcftools_set_genotypes

work = Path(tempfile.mkdtemp())
bcftools_set_genotypes("variants.vcf", work / "ref.vcf", target_genotypes=".", new_genotype="0p")
bcftools_set_genotypes("trio.vcf", work / "low_depth.vcf", target_genotypes="q", new_genotype=".", genotype_include="FMT/DP<31")
print([line.split("\t")[11] for line in (work / "ref.vcf").read_text().splitlines()[-8:-6]])
print(sum("./." in line for line in (work / "low_depth.vcf").read_text().splitlines() if not line.startswith("#")))
```

- **Example output:**

```text
['0/0', '0/0']
14
```

## 5) `dnakit.variants.bcftools_split_vcf_by_sample()`

- **Purpose:** Split a VCF/BCF by sample or by sample group (`bcftools +split`).
- **Method:** By default one single-sample file per sample is written, named after the sample. `samples_file` (`-S`) lists the samples of each output file (optional second column renames, optional third names the file); `groups_file` (`-G`) assigns each sample to one or more outputs. `keep_tags` limits the tags kept, for example `['INFO', 'FMT/GT']` (`-k`). The files go into the new folder `output_dir`.
- **API:** `dnakit.variants.bcftools_split_vcf_by_sample(vcf_path, output_dir, *, samples_file=None, groups_file=None, keep_tags=None, regions=None, include=None, exclude=None, output_type='v', bcftools_path=None)`.
- **Input:** A multi-sample VCF/BCF and a new output folder.
- **Output:** The common fields, `output_type` and `files`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_plugins_p06_cli_functions import bcftools_split_vcf_by_sample

work = Path(tempfile.mkdtemp())
(work / "groups.txt").write_text("S1,S2\t-\tpair\nS3\n")
print(bcftools_split_vcf_by_sample("variants.vcf", work / "each", keep_tags=["INFO", "FMT/GT"])["files"])
print(sorted(bcftools_split_vcf_by_sample("variants.vcf", work / "grouped", samples_file=work / "groups.txt")["files"]))
```

- **Example output:**

```text
['S1.vcf', 'S2.vcf', 'S3.vcf', 'S4.vcf', 'S5.vcf', 'S6.vcf']
['S3.vcf', 'pair.vcf']
```

## 6) `dnakit.variants.bcftools_split_vep()`

- **Purpose:** Query or extract structured consequence annotations: VEP `CSQ`, bcftools `BCSQ` and SnpEff `ANN` (`bcftools +split-vep`).
- **Method:** `annotation` names the INFO tag (auto-detected). `list_fields` lists the annotation fields (`-l`). `columns` extracts fields by name or 0-based index into new INFO tags (`Consequence,IMPACT`, `1-3`, `-` for all; `:Integer`/`:Float` set types). `select` is `TR:CSQ[:PRN]` such as `worst` or `primary:missense+`. With `text_format` (a `bcftools query` format, `-f`) the output is a text table instead of a VCF. `annotation_prefix`, `all_fields_delimiter`, `severity_file`, `gene_list`, `duplicate`, `drop_sites`, `allow_undefined_tags` and `print_header` are `-p -A -S -g -d -x/-X -u -H`.
- **API:** `dnakit.variants.bcftools_split_vep(vcf_path, output_path, *, annotation=None, columns=None, select=None, text_format=None, list_fields=False, annotation_prefix=None, all_fields_delimiter=None, severity_file=None, gene_list=None, duplicate=False, drop_sites=None, allow_undefined_tags=False, print_header=False, regions=None, include=None, exclude=None, output_type=None, bcftools_path=None)`.
- **Input:** A VCF/BCF with a consequence tag and a new output path.
- **Output:** The common fields; `records` for VCF output, `report` for text output and, with `list_fields`, `fields` (index and name).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_plugins_p06_cli_functions import bcftools_split_vep

work = Path(tempfile.mkdtemp())
fields = bcftools_split_vep("csq.vcf", work / "fields.txt", list_fields=True)
print([field["name"] for field in fields["fields"]])
table = bcftools_split_vep("csq.vcf", work / "table.txt", select="worst", text_format="%CHROM\t%POS\t%Consequence\t%IMPACT\n")
print(table["report"]["rows"])
vcf = bcftools_split_vep("csq.vcf", work / "vep.vcf", columns="Consequence,IMPACT", select="worst")
print(vcf["records"], (work / "vep.vcf").read_text().splitlines()[-3].split("\t")[7])
```

- **Example output:**

```text
['Allele', 'Consequence', 'IMPACT', 'SYMBOL', 'Feature']
[['chr1', 10, 'missense_variant', 'MODERATE'], ['chr1', 20, 'stop_gained', 'HIGH']]
3 CSQ=G|missense_variant|MODERATE|GENE1|T1,G|synonymous_variant|LOW|GENE1|T2;Consequence=missense_variant;IMPACT=MODERATE
```

## 7) `dnakit.variants.bcftools_convert_format_tags()`

- **Purpose:** Convert between equivalent FORMAT tags: GL, PL, GP and GT, QR/QA to QS, and localized alleles (`bcftools +tag2tag`).
- **Method:** `conversion` is `GL-to-PL`, `GL-to-GP`, `GL-to-GT`, `PL-to-GL`, `PL-to-GP`, `PL-to-GT`, `GP-to-GL`, `GP-to-PL`, `GP-to-GT`, `QR-QA-to-QS`, `LXX-to-XX`, `LPL-to-PL` or `LAD-to-AD`. `replace_source` drops the source tag (`-r`), `gp_threshold` is the largest error for a GP to GT hard call (`-t`), `skip_nalt` leaves sites with more ALT alleles unchanged in the localized conversions (`-s`), `defaults` fills missing localized fields (`-d`). The opposite direction `--XX-to-LXX` is advertised by bcftools 1.24 but prints "todo", so it is not offered.
- **API:** `dnakit.variants.bcftools_convert_format_tags(vcf_path, output_path, *, conversion, replace_source=False, gp_threshold=0.1, skip_nalt=0, defaults=None, regions=None, include=None, exclude=None, output_type=None, bcftools_path=None)`.
- **Input:** A VCF/BCF containing the source tags and a new output path.
- **Output:** The common fields, `output_type` and `records` (plain VCF output only).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_plugins_p06_cli_functions import bcftools_convert_format_tags

output = Path(tempfile.mkdtemp()) / "gl.vcf"
result = bcftools_convert_format_tags("trio.vcf", output, conversion="PL-to-GL")
row = output.read_text().splitlines()[-1].split("\t")
print(result["records"], row[8])
print(row[9])
```

- **Example output:**

```text
14 GT:PL:AD:DP:GQ:QS:SP:GL
1/1:255,60,0:0,30:30:60:0,1140:0:-25.5,-6,-0
```

## 8) `dnakit.variants.bcftools_variant_distance()`

- **Purpose:** Annotate sites with the distance to the nearest variant (`bcftools +variant-distance`).
- **Method:** `direction` is `nearest` (default), `fwd`, `rev` or `both` (two values: distance to the previous and to the next variant); `tag_name` is the INFO tag written (`-n`).
- **API:** `dnakit.variants.bcftools_variant_distance(vcf_path, output_path, *, direction='nearest', tag_name='DIST', regions=None, output_type=None, bcftools_path=None)`.
- **Input:** A VCF/BCF and a new output path.
- **Output:** The common fields, `output_type` and `records` (plain VCF output only).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_plugins_p06_cli_functions import bcftools_variant_distance

output = Path(tempfile.mkdtemp()) / "distance.vcf"
bcftools_variant_distance("variants.vcf", output, direction="both")
for line in output.read_text().splitlines()[-9:-6]:
    print(line.split("\t")[1], line.split("\t")[7])
```

- **Example output:**

```text
120 DP=55;AF=0.333;DIST=60,60
180 DP=60;AF=0.5;DIST=60,60
240 DP=65;AF=0.417;DIST=60,60
```

## 9) `dnakit.variants.bcftools_fixref()`

- **Purpose:** Report or repair REF/ALT strand orientation and swaps against a reference (`bcftools +fixref`; the work-list rows for `-m flip|flip-all|id|ref-alt|swap|top` and `-m stats`).
- **Method:** `mode="stats"` (default) saves the strand and mismatch statistics (`SC`, `ST`, `NS` rows, parsed into `sections`). The conversion modes write a VCF/BCF: `flip`, `flip-all` (also A/T and C/G sites), `id` (uses the ID column against `dbsnp_vcf`), `ref-alt`, `swap` and `top` (Illumina TOP strand to forward). `discard_unresolved` drops sites that cannot be repaired (`-d`); `tag_name` renames `INFO/FIXREF` (`-t`). Do not convert without understanding the strand convention of the data.
- **API:** `dnakit.variants.bcftools_fixref(vcf_path, output_path, *, reference_fasta, mode='stats', dbsnp_vcf=None, discard_unresolved=False, tag_name=None, regions=None, output_type=None, bcftools_path=None)`.
- **Input:** A VCF/BCF, an indexed reference FASTA and a new output path.
- **Output:** The common fields and `mode`; `report` and `sections` in `stats` mode, `records` for a conversion.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_plugins_p06_cli_functions import bcftools_fixref

import subprocess

subprocess.run(["samtools", "faidx", "genome.fa"], check=True)
work = Path(tempfile.mkdtemp())
stats = bcftools_fixref("variants.vcf", work / "stats.txt", reference_fasta="genome.fa")
print({str(row[0]).strip(): row[1] for row in stats["sections"]["NS"]})
fixed = bcftools_fixref("variants.vcf", work / "flip.vcf", reference_fasta="genome.fa", mode="flip")
print(fixed["mode"], fixed["records"])
```

- **Example output:**

```text
{'total': 10, 'ref match': 8, 'ref mismatch': 0, 'errors': 0, 'skipped': 2, 'non-ACGT': 0, 'non-SNP': 2, 'non-biallelic': 0}
flip 10
```

## 10) `dnakit.variants.bcftools_color_chromosomes()`

- **Purpose:** Find chromosomal segments shared in a trio or between two unrelated samples from phased genotypes (`bcftools +color-chrs`).
- **Method:** `trio` is `[mother, father, child]` (`-t`) and `unrelated` two sample names (`-u`); give exactly one. The segments are written to `<prefix>.dat` in the new folder `output_dir`, the format that bcftools' `color-chrs.pl` draws.
- **API:** `dnakit.variants.bcftools_color_chromosomes(vcf_path, output_dir, *, trio=None, unrelated=None, prefix='color-chrs', regions=None, bcftools_path=None)`.
- **Input:** A VCF/BCF with phased genotypes and a new output folder.
- **Output:** The common fields, `files`, `segments` (`SG` rows: chromosome, start, end, haplotype labels) and `switches` (`SW` rows).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_plugins_p06_cli_functions import bcftools_color_chromosomes

import subprocess

work = Path(tempfile.mkdtemp())
subprocess.run(["bcftools", "view", "-Oz", "-o", str(work / "phased.vcf.gz"), "trio_phased.vcf"], check=True)
subprocess.run(["bcftools", "index", str(work / "phased.vcf.gz")], check=True)
result = bcftools_color_chromosomes(work / "phased.vcf.gz", work / "colors", trio=["M", "F", "C"])
print(result["files"], result["segments"], result["switches"])
```

- **Example output:**

```text
['color-chrs.dat'] [['chr1', 1, 399, 'M:2', 'F:2'], ['chr1', 400, 439, 'M:1', 'F:2'], ['chr1', 440, 479, 'M:1', 'F:1']] [['M', 'chr1', 4, 1, 0.333333], ['F', 'chr1', 6, 1, 0.2]]
```

## 11) `dnakit.variants.bcftools_parental_origin()`

- **Purpose:** Determine the parental origin of a copy-number variant in a trio (`bcftools +parental-origin`).
- **Method:** `trio` is `[proband, father, mother]` (`-p`), `cnv_type` `del` or `dup` (`-t`) and `region` the CNV as `chr:start-end` (`-r`). Parental heterozygous sites with a skewed ALT fraction are excluded below `min_binomial_probability` (`-b`); `greedy` also uses ambiguous sites (`-g`); `list_informative_sites` lists them (`-d`). The VCF needs an index, `FORMAT/PL` and `FORMAT/AD`.
- **API:** `dnakit.variants.bcftools_parental_origin(vcf_path, output_path, *, trio, cnv_type, region, min_binomial_probability=0.01, greedy=False, list_informative_sites=False, include=None, exclude=None, bcftools_path=None)`.
- **Input:** An indexed VCF/BCF and a new output path.
- **Output:** The common fields, `report`, `predicted_origin` (`paternal`, `maternal` or `uncertain`), `quality` and `markers`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_trio_p06_cli_functions import bcftools_parental_origin

import subprocess

work = Path(tempfile.mkdtemp())
lines = [
    "##fileformat=VCFv4.2", "##contig=<ID=chr1,length=100000>",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depth">',
    '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred likelihoods">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tP\tF\tM",
]
for index in range(30):  # proband carries an extra ALT-bearing copy that the father's heterozygous sites explain
    lines.append(f"chr1\t{1000 + 100 * index}\t.\tA\tG\t50\tPASS\t.\tGT:AD:PL\t0/1:10,20:60,0,60\t0/1:15,15:60,0,60\t0/0:30,0:0,60,255")
(work / "cnv.vcf").write_text("\n".join(lines) + "\n")
subprocess.run(["bcftools", "view", "-Oz", "-o", str(work / "cnv.vcf.gz"), str(work / "cnv.vcf")], check=True)
subprocess.run(["bcftools", "index", str(work / "cnv.vcf.gz")], check=True)
result = bcftools_parental_origin(work / "cnv.vcf.gz", work / "origin.txt", trio=["P", "F", "M"], cnv_type="dup", region="chr1:1-5000")
print(result["predicted_origin"], result["markers"])
```

- **Example output:**

```text
paternal 30
```

## 12) `dnakit.variants.bcftools_trio_stats()`

- **Purpose:** Transmission rate and Mendelian statistics in trio children (`bcftools +trio-stats`).
- **Method:** The trios come from `ped_path` (`-p`) or `trio=[proband, father, mother]` (`-P`). `include` / `exclude` may hold curly-bracket ranges such as `GQ>{10,20}` that are scanned together, one `FLT` row per filter and trio. `max_alt_trios` limits the transmission rate to sites with at most that many alternate trios (`-a`); `debug` is a list from `mendel-errors` and `transmitted` (`-d`).
- **API:** `dnakit.variants.bcftools_trio_stats(vcf_path, output_path, *, ped_path=None, trio=None, max_alt_trios=0, debug=None, regions=None, include=None, exclude=None, bcftools_path=None)`.
- **Input:** A VCF/BCF, a PED file or trio names, and a new output path.
- **Output:** The common fields, `report`, `trio_columns` and `trios`: one dictionary per filter and trio with valid genotypes, non-reference genotypes, DNMs or Mendelian errors, singleton counts, ts, tv, ts/tv, homozygous and recurrent DNMs.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_trio_p06_cli_functions import bcftools_trio_stats

result = bcftools_trio_stats("trio.vcf", Path(tempfile.mkdtemp()) / "stats.txt", ped_path="trio.ped", include="GQ>{10,20}")
for row in result["trios"]:
    print(row["filter"], row["child"], row["valid_gts"], row["dnm_or_mendel_errors"], row["transmitted_singletons"])
```

- **Example output:**

```text
FLT0 C 14 4 2
FLT1 C 14 4 2
```

## 13) `dnakit.variants.bcftools_trio_switch_rate()`

- **Purpose:** Phase switch rate in trio children from phased genotypes (`bcftools +trio-switch-rate`).
- **Method:** `ped_path` lists the trios; an optional 7th column groups the results by population.
- **API:** `dnakit.variants.bcftools_trio_switch_rate(vcf_path, output_path, *, ped_path, regions=None, bcftools_path=None)`.
- **Input:** A VCF/BCF with phased genotypes, a PED file and a new output path.
- **Output:** The common fields, `report` and `trios`: father, mother, child, genotypes tested, Mendelian errors, switches and switch percentage.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_trio_p06_cli_functions import bcftools_trio_switch_rate

result = bcftools_trio_switch_rate("trio_phased.vcf", Path(tempfile.mkdtemp()) / "switch.txt", ped_path="trio.ped")
print(result["trios"])
```

- **Example output:**

```text
[{'father': 'F', 'mother': 'M', 'child': 'C', 'tested': 4, 'mendelian_errors': 2, 'switches': 0, 'switch_percent': 0.0}]
```

## 14) `dnakit.variants.bcftools_mendelian_consistency()`

- **Purpose:** Count, annotate or filter genotypes by Mendelian consistency in trios (`bcftools +mendelian2`; the work-list rows for `-m c` and for `-m a|d|e|E|g|m|M|S`).
- **Method:** The trios come from `ped_path` (`-P`) or `trio=[child, father, mother]` (`-p`; a `1X:`/`2X:` prefix selects the chrX pattern). `mode="c"` (default) saves a text summary with the number of good, erroneous and missing genotypes per trio. Otherwise `mode` is a string of VCF-output letters: `a` add `INFO/MERR` and related counts, `d` delete genotypes of inconsistent trios, `e` keep only sites with an error, `E` drop them, `g` keep sites with a good trio, `m` keep sites with missing genotypes, `M` drop them, `S` drop skipped sites. `rules` is a predefined inheritance rule set and `rules_file` a custom one.
- **API:** `dnakit.variants.bcftools_mendelian_consistency(vcf_path, output_path, *, ped_path=None, trio=None, mode='c', rules='GRCh37', rules_file=None, regions=None, include=None, exclude=None, output_type=None, bcftools_path=None)`.
- **Input:** A VCF/BCF with `FORMAT/GT`, trios and a new output path.
- **Output:** The common fields and `mode`; `report` for `c`, `records` for the VCF modes.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_trio_p06_cli_functions import bcftools_mendelian_consistency

work = Path(tempfile.mkdtemp())
counts = bcftools_mendelian_consistency("trio.vcf", work / "counts.txt", ped_path="trio.ped")
print({row[0]: row[1] for row in counts["report"]["rows"] if row[0] in ("sites_merr", "sites_good")})
for mode in ("E", "e"):
    print(mode, bcftools_mendelian_consistency("trio.vcf", work / f"{mode}.vcf", ped_path="trio.ped", mode=mode)["records"])
```

- **Example output:**

```text
{'sites_merr': 4, 'sites_good': 10}
E 10
e 4
```

## 15) `dnakit.variants.bcftools_trio_denovo()`

- **Purpose:** Screen variants for de novo mutations in trios and write `FORMAT/DNM` (`bcftools +trio-dnm3`; the four model rows `--use-DMM`, `--use-ALM`, `--use-DNG`, `--use-NAIVE`).
- **Method:** `model` is `dmm` (Dirichlet-multinomial with site-noise awareness, the default of bcftools 1.24), `alm` (allele-likelihood model, default up to 1.23), `dng` (original DeNovoGear) or `naive` (genotypes only, `DNM:flag`). The trios come from `ped_path` (`-P`) or `trio=[child, father, mother]` (`-p`). `dnm_tag` is the output tag and score type (`DNM:log` by default; `log`, `flag`, `phred`, `prob`), `min_score` omits scores below it (`-m`), `chrx_regions` the chrX regions or preset (`-X`), `mutation_rate` the prior (`--mrate`), `strictly_novel` `-n`. Model options: `noise_prior` and `phi` (dmm, alm), `min_vaf` and `strand_bias` (dmm, alm; `strand_bias` also dng), `max_qm` (dmm; negative runs without `FORMAT/QM`), `allele_dropout` (alm), `use_parental_ad`, `use_child_ad` (dmm), `force_ad`. An option that does not belong to the chosen model is refused.
- **API:** `dnakit.variants.bcftools_trio_denovo(vcf_path, output_path, *, ped_path=None, trio=None, model='dmm', min_score=None, dnm_tag=None, chrx_regions='GRCh37', mutation_rate=1e-08, strictly_novel=False, noise_prior=None, phi=None, min_vaf=None, max_qm=None, strand_bias=None, allele_dropout=None, use_parental_ad=False, use_child_ad=False, force_ad=False, regions=None, include=None, exclude=None, output_type=None, bcftools_path=None)`.
- **Input:** A VCF/BCF with `GT`, `PL`, `AD` (and `QS`), trios and a new output path.
- **Output:** The common fields, `model`, `output_type` and `records` (plain VCF output only).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_trio_p06_cli_functions import bcftools_trio_denovo

work = Path(tempfile.mkdtemp())
naive = bcftools_trio_denovo("trio.vcf", work / "naive.vcf", trio=["C", "F", "M"], model="naive")
flagged = [
    line.split("\t")[1]
    for line in (work / "naive.vcf").read_text().splitlines()
    if not line.startswith("#") and "DNM" in line.split("\t")[8].split(":")
]
alm = bcftools_trio_denovo("trio.vcf", work / "alm.vcf", trio=["C", "F", "M"], model="alm")
scores = {}
for line in (work / "alm.vcf").read_text().splitlines():
    if not line.startswith("#"):
        fields = line.split("\t")
        scores[int(fields[1])] = float(fields[9].split(":")[fields[8].split(":").index("DNM")])
print(naive["model"], flagged)
print(alm["records"], max(scores, key=scores.get), round(max(scores.values()), 2))
```

- **Example output:**

```text
naive ['200', '280', '480', '560']
14 280 -3.93
```
