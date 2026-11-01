# bcftools plugin facades: statistics, reports and annotation

This group of functions in `dnakit.variants` (module `bcftools_plugins_p06_cli_functions`) wraps plugins of the external program [bcftools](https://github.com/samtools/bcftools) (`bcftools +plugin`, 1.24 is the tested version): allele and genotype statistics, annotation, genotype comparison and subsetting. Each function validates its arguments, runs one plugin with a fixed argument list and no shell, captures standard output and standard error separately, and returns a JSON-serialisable dictionary. Plugins that print a report have it saved to `output_path` and parsed into `report` (comment lines, column names and numeric rows). Plugins that change records write a VCF/BCF to `output_path` (`output_type` `v`, `z`, `b` or `u`, default from the file name) and return the record count for plain VCF output. An existing output is never overwritten and a failed run leaves nothing behind.

bcftools is an external program with its plugins: install it (for example `conda install -c bioconda bcftools`, or a release from https://www.htslib.org/download/), then put it on PATH or pass `bcftools_path`. Plugins are found through the `BCFTOOLS_PLUGINS` environment variable or, when it is not set, the `libexec/bcftools` folder next to the installation. If the program cannot be found, `BackendUnavailableError` is raised with the install hint; if it exits with an error, `BackendExecutionError` (code `HTSLIB_EXECUTION_FAILED`) is raised with the tail of its log.

The common bcftools options are `regions` (`-r`, needs an index), `include` (`-i`) and `exclude` (`-e`). Not offered: `bcftools +frameshifts`, which in 1.24 writes `OOF=-1` for every indel (its type test `var[i].type != VCF_INDEL` never matches since htslib marks indels `VCF_INS`/`VCF_DEL`), and the indel-stats option `--nvaf`. Results share `status`, `backend`, `operation`, `version`, `command`, `outputs`, `log_tail` and `provenance`; the sections below list only the extra fields.

## Example data

The examples read these files from the working directory and write their outputs to a temporary directory. `variants.vcf` has six samples and ten sites (a deletion and an insertion among them); `trio.vcf` has a child C, a father F, a mother M and an unrelated sample U with `GT`, `PL`, `AD`, `DP`, `GQ`, `QS` and `SP`; `trio_indel.vcf` is the same with indels at every third site.

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

**Example file:** `trio_indel.vcf`

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
chr1	40	.	T	TCG	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	0/1:60,0,60:15,15:30:60:570,570:0	0/1:60,0,60:15,15:30:60:570,570:0	0/0:0,60,255:30,0:30:60:1140,0:0	0/0:0,60,255:30,0:30:60:1140,0:0
chr1	80	.	T	A	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	1/1:255,60,0:0,30:30:60:0,1140:0	0/1:60,0,60:15,15:30:60:570,570:0	0/1:60,0,60:15,15:30:60:570,570:0	0/1:60,0,60:15,15:30:60:570,570:0
chr1	120	.	C	G	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	0/1:60,0,60:15,15:30:60:570,570:0	1/1:255,60,0:0,30:30:60:0,1140:0	0/0:0,60,255:30,0:30:60:1140,0:0	0/0:0,60,255:30,0:30:60:1140,0:0
chr1	160	.	TTT	T	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	0/0:0,60,255:30,0:30:60:1140,0:0	0/0:0,60,255:30,0:30:60:1140,0:0	0/0:0,60,255:30,0:30:60:1140,0:0	0/1:60,0,60:15,15:30:60:570,570:0
chr1	200	.	G	T	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	0/1:60,0,60:15,15:30:60:570,570:0	0/0:0,60,255:30,0:30:60:1140,0:0	0/0:0,60,255:30,0:30:60:1140,0:0	0/0:0,60,255:30,0:30:60:1140,0:0
chr1	240	.	A	C	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	0/0:0,60,255:30,0:30:60:1140,0:0	0/1:60,0,60:15,15:30:60:570,570:0	0/1:60,0,60:15,15:30:60:570,570:0	1/1:255,60,0:0,30:30:60:0,1140:0
chr1	280	.	A	ACG	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	1/1:255,60,0:0,30:30:60:0,1140:0	0/0:0,60,255:30,0:30:60:1140,0:0	0/1:60,0,60:15,15:30:60:570,570:0	0/0:0,60,255:30,0:30:60:1140,0:0
chr1	320	.	G	T	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	0/1:60,0,60:15,15:30:60:570,570:0	0/1:60,0,60:15,15:30:60:570,570:0	0/0:0,60,255:30,0:30:60:1140,0:0	0/0:0,60,255:30,0:30:60:1140,0:0
chr1	360	.	A	C	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	1/1:255,60,0:0,30:30:60:0,1140:0	0/1:60,0,60:15,15:30:60:570,570:0	0/1:60,0,60:15,15:30:60:570,570:0	0/1:60,0,60:15,15:30:60:570,570:0
chr1	400	.	TTG	T	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	0/1:60,0,60:15,15:30:60:570,570:0	1/1:255,60,0:0,30:30:60:0,1140:0	0/0:0,60,255:30,0:30:60:1140,0:0	0/0:0,60,255:30,0:30:60:1140,0:0
chr1	440	.	A	C	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	0/0:0,60,255:30,0:30:60:1140,0:0	0/0:0,60,255:30,0:30:60:1140,0:0	0/0:0,60,255:30,0:30:60:1140,0:0	0/1:60,0,60:15,15:30:60:570,570:0
chr1	480	.	T	A	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	0/1:60,0,60:15,15:30:60:570,570:0	0/0:0,60,255:30,0:30:60:1140,0:0	0/0:0,60,255:30,0:30:60:1140,0:0	0/0:0,60,255:30,0:30:60:1140,0:0
chr1	520	.	G	GCG	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	0/0:0,60,255:30,0:30:60:1140,0:0	0/1:60,0,60:15,15:30:60:570,570:0	0/1:60,0,60:15,15:30:60:570,570:0	1/1:255,60,0:0,30:30:60:0,1140:0
chr1	560	.	T	A	50	PASS	DP=120	GT:PL:AD:DP:GQ:QS:SP	1/1:255,60,0:0,30:30:60:0,1140:0	0/0:0,60,255:30,0:30:60:1140,0:0	0/1:60,0,60:15,15:30:60:570,570:0	0/0:0,60,255:30,0:30:60:1140,0:0
```

**Example file:** `reads.sam`

```text
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:600
@SQ	SN:chr2	LN:400
r1	0	chr1	1	60	50M	*	0	0	GCTAAAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCA	??II:I:I::II::5F:I?5D:IF5?:5?DID:??:DFFDI?DI?I5?FI
r2	0	chr1	41	60	50M	*	0	0	TGTTGGCCCAGTGTGAATCGCTTAAGGGTTAAGTAAGTGTGATGCATACG	?IDFFI:DIIDF:?::F???D?F?D:FDF:F::I?:5DIIDD5ID5::DI
r3	0	chr1	81	60	50M	*	0	0	GATGCATACGCCTTTACTTGCTGTGTCCACCCCATCGGACTGGCATTTTT	DI?:IDI??:II?FII:IF:ID:F5FIF5?FIF?:?::IF5::5DIFFFI
r4	0	chr1	121	60	50M	*	0	0	TGGCATTTTTATTACACTCAGAAACAGAACTCGGGTAATTTTGACAGGTC	5D?I5F5:?5D??FFI5DF?55??:DFFFD:I:?5DFI:I:?DI:?DF:?
r5	0	chr1	161	60	50M	*	0	0	TTGACAGGTCACGCAGAGGCGCGCCCTCCTGAAGTGCGTGGACACTCGCT	:F?D:I?D5?DD555F5I?:DDFD?5?FD::55?::??FD?5?:5?DDFF
r6	0	chr1	201	60	50M	*	0	0	GACACTCGCTATGAATCTCTGATTTACCCACTCTGCCAAACTCCAGCGCG	D5?F?:5?:FI:5:I:5D?5??:5I:I:FD:FDF5F55F:F?FD?:F??5
r7	0	chr1	241	60	50M	*	0	0	CTCCAGCGCGGTCAGTTCCATCACCCTAAGTAACCGAATAATGCGTTCGC	?55:I55FI?5?ID:FIII?555?DDDD?::D?5F5F???::?D:DD?:D
r8	0	chr1	281	60	50M	*	0	0	ATGCGTTCGCTCTATTGACTACGACGCGCTCATTCCCTTGTCGGAGAGTT	FFD:DID:5F5DF::55F:?5DFFI5:IIIID5F5I:?F:I5IFD?:D5F
r9	0	chr2	1	60	50M	*	0	0	TTCGTACCTTGGGGGTCGTTACCACTCTGTTCCCACGAGCGGCATTTCTG	5:?FF5FD?I:F5FID:FI55?:D?:I?I:I?F5?I::I:F55FI:FIF5
r10	0	chr2	61	60	50M	*	0	0	TTTTGACATTTAATTTCACCCATAAACCAGCGTAAAGCTGCAAGTGGCTC	D?DFFDFFD5DF??DI5:I?DF5DD?5IDDDDFI55?:5:F??F:DF?5:
r11	0	chr2	121	60	50M	*	0	0	GCTGCTAGTGTCAGACTCGCCTCGGATCCTTACTACACTAACTTGAACGC	I5I5I:5FDI:?:5:I:55:DF:IIDI5DF:55?F5FFI:5ID?D?:??F
r12	0	chr2	181	60	50M	*	0	0	AAGAGTACTGGTAATCGTCGGTATCTATATAAGCAGGGGAGGGGAAACAT	IFD??II:5I5:55?::F?55DDF55IFFD:ID:D??55ID5I5?D5D:I
```

**Example file:** `sites.txt`

```text
chr1	60	G	T
chr1	120	C	G
```

## 1) `dnakit.variants.bcftools_ad_bias()`

- **Purpose:** Find sites where the ALT allele fraction (`FORMAT/AD`) differs between a sample and its control, with a Fisher exact test (`bcftools +ad-bias`).
- **Method:** Write the `[sample, control]` pairs to a temporary pair file and run the plugin. Sites with a p-value below `threshold` (`-t`) give one `FT` row with the allele depths of both samples. `min_alt_depth` and `min_depth` are `-a` and `-d`, `variant_type` is `snp` or `indel` (`-v`), `extra_columns` a `bcftools query` format appended to each row (`-f`). With `clean_vcf=True` (`-c`) the output is a VCF without the sites and ALT alleles that fail the threshold.
- **API:** `dnakit.variants.bcftools_ad_bias(vcf_path, output_path, *, sample_pairs, threshold=0.001, min_alt_depth=1, min_depth=0, variant_type=None, extra_columns=None, clean_vcf=False, regions=None, include=None, exclude=None, bcftools_path=None)`.
- **Input:** A VCF/BCF with `FORMAT/AD`, a new output path and the sample pairs.
- **Output:** The common fields and `report` (comment lines, column names, numeric `rows`); in `clean_vcf` mode `records` instead.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_plugins_p06_cli_functions import bcftools_ad_bias

output = Path(tempfile.mkdtemp()) / "adbias.txt"
result = bcftools_ad_bias("trio.vcf", output, sample_pairs=[["F", "M"]])
for row in result["report"]["rows"][:3]:
    print(row[:7], f"{row[-1]:.2e}")
```

- **Example output:**

```text
['FT', 'F', 'M', 'chr1', 40, 'T', 'A'] 5.83e-06
['FT', 'F', 'M', 'chr1', 120, 'C', 'G'] 1.69e-17
['FT', 'F', 'M', 'chr1', 280, 'A', 'C'] 5.83e-06
```

## 2) `dnakit.variants.bcftools_af_distribution()`

- **Purpose:** Collect the distribution of allele-frequency deviations and of genotype probabilities under Hardy-Weinberg equilibrium (`bcftools +af-dist`).
- **Method:** Only non-reference genotypes at sites with a known allele frequency are used; their probability is `2*AF*(1-AF)` for heterozygous and `AF**2` for homozygous alternate. `deviation_bins` and `probability_bins` are the bin edges (`-d`, `-p`; default 0, 0.1, ..., 1), `af_tag` the INFO tag with the frequency (`-t`), `per_sample` adds the geometric mean probability per sample (`-s`). `regions`, `include`, `exclude` are the bcftools site selection options `-r`, `-i`, `-e` (`regions` needs an indexed file).
- **API:** `dnakit.variants.bcftools_af_distribution(vcf_path, output_path, *, deviation_bins=None, probability_bins=None, per_sample=False, af_tag='AF', regions=None, include=None, exclude=None, bcftools_path=None)`.
- **Input:** A VCF/BCF with an INFO allele-frequency tag and `FORMAT/GT`, and a new output path.
- **Output:** The common fields and `report`; `PROB_DIST` and `DEV_DIST` rows are bin start, bin end and count.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_plugins_p06_cli_functions import bcftools_af_distribution

output = Path(tempfile.mkdtemp()) / "afdist.txt"
result = bcftools_af_distribution("variants.vcf", output, deviation_bins=[0, 0.25, 0.5, 1], probability_bins=[0, 0.5, 1])
for row in result["report"]["rows"]:
    print(row)
```

- **Example output:**

```text
['PROB_DIST', 0.0, 0.5, 33]
['PROB_DIST', 0.5, 1.0, 5]
['DEV_DIST', 0.0, 0.25, 10]
['DEV_DIST', 0.25, 0.5, 0]
['DEV_DIST', 0.5, 1.0, 0]
```

## 3) `dnakit.variants.bcftools_allele_length_counts()`

- **Purpose:** Count how often each REF and ALT allele length occurs (`bcftools +allele-length`).
- **Method:** Run the plugin and parse its table. `regions`, `include`, `exclude` are the bcftools site selection options `-r`, `-i`, `-e` (`regions` needs an indexed file).
- **API:** `dnakit.variants.bcftools_allele_length_counts(vcf_path, output_path, *, regions=None, include=None, exclude=None, bcftools_path=None)`.
- **Input:** A VCF/BCF and a new output path.
- **Output:** The common fields and `report` with the columns `LENGTH`, `REF`, `ALT`, `REF+ALT` and `REF+ALT WITH NON-BASE NUCLEOTIDES`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_plugins_p06_cli_functions import bcftools_allele_length_counts

output = Path(tempfile.mkdtemp()) / "lengths.txt"
result = bcftools_allele_length_counts("variants.vcf", output)
print(result["report"]["columns"][:3])
print([row[:3] for row in result["report"]["rows"] if any(row[1:])])
```

- **Example output:**

```text
['LENGTH', 'REF', 'ALT']
[[1, 9, 9], [2, 0, 0], [3, 1, 1], [4, 0, 0], ['', '', '']]
```

## 4) `dnakit.variants.bcftools_check_ploidy()`

- **Purpose:** Report the ploidy of every sample along each chromosome and whether it is consistent across sites (`bcftools +check-ploidy`).
- **Method:** Run the plugin. Missing and half-missing genotypes are ignored unless `use_missing` is set (`-m`).
- **API:** `dnakit.variants.bcftools_check_ploidy(vcf_path, output_path, *, use_missing=False, regions=None, bcftools_path=None)`.
- **Input:** A VCF/BCF with `FORMAT/GT` and a new output path.
- **Output:** The common fields and `report`: one row per sample and region with sample, chromosome, region start, region end and ploidy.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_plugins_p06_cli_functions import bcftools_check_ploidy

output = Path(tempfile.mkdtemp()) / "ploidy.txt"
result = bcftools_check_ploidy("variants.vcf", output)
print(result["report"]["columns"])
print(result["report"]["rows"][:2], len(result["report"]["rows"]))
```

- **Example output:**

```text
['Sample', 'Chromosome', 'Region Start', 'Region End', 'Ploidy']
[['S1', 'chr1', 60, 570, 2], ['S2', 'chr1', 60, 570, 2]] 6
```

## 5) `dnakit.variants.bcftools_check_sparsity()`

- **Purpose:** List the samples that have fewer than a minimum number of genotypes in a region or chromosome (`bcftools +check-sparsity`).
- **Method:** Run the plugin with `-n min_markers`; without `regions` or `regions_file` (`-r`, `-R`) every chromosome is checked. An empty report means that every sample has enough markers.
- **API:** `dnakit.variants.bcftools_check_sparsity(vcf_path, output_path, *, min_markers=1, regions=None, regions_file=None, bcftools_path=None)`.
- **Input:** A VCF/BCF with `FORMAT/GT` and a new output path.
- **Output:** The common fields and the sample list saved to `output_path`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_plugins_p06_cli_functions import bcftools_check_sparsity

work = Path(tempfile.mkdtemp())
bcftools_check_sparsity("variants.vcf", work / "all.txt", min_markers=1)
bcftools_check_sparsity("variants.vcf", work / "dense.txt", min_markers=10)
print(repr((work / "all.txt").read_text()))
print((work / "dense.txt").read_text().split())
```

- **Example output:**

```text
''
['chr1', 'S6']
```

## 6) `dnakit.variants.bcftools_genotype_dosage()`

- **Purpose:** Print the genotype dosage (0 to 2, -1 for missing) per site and sample, from PL, GL or GT (`bcftools +dosage`).
- **Method:** Run the plugin; `tags` is the priority list of FORMAT tags the dosage is taken from (`-t`).
- **API:** `dnakit.variants.bcftools_genotype_dosage(vcf_path, output_path, *, tags=('PL', 'GL', 'GT'), regions=None, include=None, exclude=None, bcftools_path=None)`.
- **Input:** A VCF/BCF and a new output path.
- **Output:** The common fields and `report`: one row per site with CHROM, POS, REF, ALT and one dosage per sample.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_plugins_p06_cli_functions import bcftools_genotype_dosage

output = Path(tempfile.mkdtemp()) / "dosage.txt"
result = bcftools_genotype_dosage("variants.vcf", output, tags=["GT"])
for row in result["report"]["rows"][:3]:
    print(row)
```

- **Example output:**

```text
['chr1', 60, 'G', 'T', 1.0, 0.0, 2.0, 1.0, 0.0, 1.0]
['chr1', 120, 'C', 'T', 0.0, 1.0, 1.0, 2.0, 0.0, 0.0]
['chr1', 180, 'C', 'A', 2.0, 1.0, 0.0, 1.0, 1.0, -1.0]
```

## 7) `dnakit.variants.bcftools_genotype_intersections()`

- **Purpose:** Count genotype intersections across all subsets of samples (`bcftools +GTisec`).
- **Method:** Run the plugin. `include_missing` adds the number of missing genotypes per sample (`-m`), `annotate_subsets` names the samples of each subset (`-v`), `human_readable` sorts the output by sample (`-H`, implies `annotate_subsets`).
- **API:** `dnakit.variants.bcftools_genotype_intersections(vcf_path, output_path, *, include_missing=False, annotate_subsets=False, human_readable=False, regions=None, include=None, exclude=None, bcftools_path=None)`.
- **Input:** A multi-sample VCF/BCF and a new output path.
- **Output:** The common fields and `report`: one row per sample subset with the number of sites at which exactly these samples share a genotype.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_plugins_p06_cli_functions import bcftools_genotype_intersections

output = Path(tempfile.mkdtemp()) / "gtisec.txt"
result = bcftools_genotype_intersections("trio.vcf", output, annotate_subsets=True)
for row in result["report"]["rows"][:5]:
    print(row)
```

- **Example output:**

```text
['@SMPS U M F C']
[4, 'U']
[2, 'M']
[2, 'F']
[10, 'C']
```

## 8) `dnakit.variants.bcftools_indel_stats()`

- **Purpose:** Per-sample indel statistics, or de novo indel statistics for the trios of a PED file (`bcftools +indel-stats`).
- **Method:** Run the plugin. `max_length` is the longest indel counted (`--max-len`), `csq_tag` the INFO tag with VEP/BCSQ consequences used to separate frameshift from inframe indels (`-c`). With `ped_path` (`-p`) only de novo indels are counted; `alt_to_ref_dnm` also counts GT errors such as 0/1 + 1/1 -> 0/0 as de novo. `include` may hold curly-bracket ranges such as `GQ>{10,20}` that are scanned together. The plugin needs `FORMAT/AD`. The documented `--nvaf` option of bcftools 1.24 is rejected by the plugin (it parses a fraction instead of a bin count), so it is not offered.
- **API:** `dnakit.variants.bcftools_indel_stats(vcf_path, output_path, *, ped_path=None, csq_tag='CSQ', max_length=20, alt_to_ref_dnm=False, regions=None, include=None, exclude=None, bcftools_path=None)`.
- **Input:** A VCF/BCF with `FORMAT/GT` and `FORMAT/AD`, and a new output path.
- **Output:** The common fields and `report` with the `SN*`, `FRM`, `VAF`, `NFRAC` and related rows.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_plugins_p06_cli_functions import bcftools_indel_stats

output = Path(tempfile.mkdtemp()) / "indels.txt"
result = bcftools_indel_stats("trio_indel.vcf", output, include="GQ>{10,20}")
for row in result["report"]["rows"]:
    if str(row[0]).startswith("SN"):
        print(row[:6])
```

- **Example output:**

```text
['SN0', 4, 5, 5, 10, 3]
['SN1', 4, 5, 5, 10, 3]
```

## 9) `dnakit.variants.bcftools_vcf_to_table()`

- **Purpose:** Print variants as text tables: a variant table, an INFO table, genotype tables and consequence tables (`bcftools +vcf2table`).
- **Method:** Run the plugin and save its text. `hide` lists the tables to leave out: `HOM_REF`/`RR`, `HET`/`AR`, `NO_CALL`/`MISSING`, `CSQ`/`VEP`, `SPLICEAI`, `ANN`/`SNPEFF`, `LOF`, `VC`, `INFO`, `GT`, `GTTYPES`, `URL`.
- **API:** `dnakit.variants.bcftools_vcf_to_table(vcf_path, output_path, *, hide=None, regions=None, include=None, exclude=None, bcftools_path=None)`.
- **Input:** A VCF/BCF and a new output path.
- **Output:** The common fields, `report` (`variants`: the number of variant blocks) and `text_preview`, the first 2000 characters of the text.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_plugins_p06_cli_functions import bcftools_vcf_to_table

output = Path(tempfile.mkdtemp()) / "table.txt"
result = bcftools_vcf_to_table("variants.vcf", output, hide=["HOM_REF", "NO_CALL", "INFO"], include="POS<=120")
print(result["report"]["variants"])
print(output.read_text().splitlines()[0])
print(output.read_text().splitlines()[3])
```

- **Example output:**

```text
2
<<< chr1:60:G (n. 1)
+--------+-------+
```

## 10) `dnakit.variants.bcftools_variant_read_frequency_score()`

- **Purpose:** Assess how noisy each site is (variant read frequency score) from many BAM/CRAM files of unaffected samples (`bcftools +vrfs`).
- **Method:** Scoring mode: `alignment_list` is a text file with one BAM/CRAM path per line (`-a`), `reference_fasta` the reference (`-f`) and `sites_file` a tab-separated file of chromosome, position, ref, alt (`-s`); `min_depth` and `vaf_bins` are `-d` and `-n`; `batch` is `I/N` to run the I-th of N batches (`-b`). Merge mode: `merge_batch_files` lists batch outputs to merge (`-M`) and the scoring arguments must be left out.
- **API:** `dnakit.variants.bcftools_variant_read_frequency_score(output_path, *, alignment_list=None, reference_fasta=None, sites_file=None, min_depth=10, vaf_bins=20, batch=None, merge_batch_files=None, bcftools_path=None)`.
- **Input:** Alignment files listed in a text file, an indexed reference FASTA and a sites file; or batch outputs to merge.
- **Output:** The common fields and `report` with the score table saved to `output_path`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_plugins_p06_cli_functions import bcftools_variant_read_frequency_score

import subprocess

work = Path(tempfile.mkdtemp())
subprocess.run(["samtools", "view", "-b", "-o", str(work / "reads.bam"), "reads.sam"], check=True)
(work / "bams.txt").write_text(f"{work / 'reads.bam'}\n")
result = bcftools_variant_read_frequency_score(
    work / "scores.txt", alignment_list=work / "bams.txt", reference_fasta="genome.fa", sites_file="sites.txt", min_depth=1
)
print(result["status"], len(result["report"]["rows"]))
print([row for row in result["report"]["rows"] if row and row[0] == "SC"][:1])
```

- **Example output:**

```text
completed 4
[]
```

## 11) `dnakit.variants.bcftools_add_variantkey()`

- **Purpose:** Add the INFO fields `VKX` (hexadecimal 64-bit VariantKey) and `RSX` (hexadecimal rsID number) (`bcftools +add-variantkey`).
- **Method:** Run the plugin and write the annotated VCF/BCF; `output_type` is `v`, `z`, `b` or `u` (default: from the file name).
- **API:** `dnakit.variants.bcftools_add_variantkey(vcf_path, output_path, *, regions=None, include=None, exclude=None, output_type=None, bcftools_path=None)`.
- **Input:** A VCF/BCF and a new output path.
- **Output:** The common fields, `output_type` and `records` (plain VCF output only).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_plugins_p06_cli_functions import bcftools_add_variantkey

output = Path(tempfile.mkdtemp()) / "keys.vcf"
result = bcftools_add_variantkey("variants.vcf", output)
print(result["records"])
print(output.read_text().splitlines()[-1].split("\t")[7])
```

- **Example output:**

```text
10
DP=95;AF=0.417;VKX=0800011c88880000;RSX=00000000
```

## 12) `dnakit.variants.bcftools_variantkey_tables()`

- **Purpose:** Generate unsorted VariantKey lookup tables in hexadecimal format (`bcftools +variantkey-hex`).
- **Method:** Run the plugin inside the new folder `output_dir`, where it writes `vkrs.unsorted.hex` (VariantKey to rsID), `rsvk.unsorted.hex` (rsID to VariantKey) and `nrvk.unsorted.tsv` (non-reversible VariantKeys).
- **API:** `dnakit.variants.bcftools_variantkey_tables(vcf_path, output_dir, *, bcftools_path=None)`.
- **Input:** A VCF/BCF and a new output folder.
- **Output:** The common fields, `files` and `counters` (the numbers bcftools prints: VariantKeys and Non-reversible VariantKeys).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_plugins_p06_cli_functions import bcftools_variantkey_tables

result = bcftools_variantkey_tables("variants.vcf", Path(tempfile.mkdtemp()) / "keys")
print(result["files"], result["counters"])
```

- **Example output:**

```text
['nrvk.unsorted.tsv', 'rsvk.unsorted.hex', 'vkrs.unsorted.hex'] {'VariantKeys': 10, 'Non-reversible VariantKeys': 0}
```

## 13) `dnakit.variants.bcftools_contrast_groups()`

- **Purpose:** Run a case/control association test and detect novel alleles or genotypes between two sample groups (`bcftools +contrast`).
- **Method:** `control_samples` is the background group (`-0`), `case_samples` the group where novel alleles or genotypes are expected (`-1`). `annotations` chooses the INFO tags written from `PASSOC` (Fisher test, REF vs non-REF), `FASSOC` (non-REF fraction in controls and cases), `NASSOC` (allele counts), `NOVELAL` and `NOVELGT` (default PASSOC, FASSOC, NOVELAL). `max_allele_freq` (`-f`) tests the enrichment of rare alleles: a value in (0, 1) is an ALT frequency, an integer an ALT count.
- **API:** `dnakit.variants.bcftools_contrast_groups(vcf_path, output_path, *, control_samples, case_samples, annotations=None, max_allele_freq=None, force_samples=False, regions=None, include=None, exclude=None, output_type=None, bcftools_path=None)`.
- **Input:** A multi-sample VCF/BCF, two lists of samples and a new output path.
- **Output:** The common fields, `output_type` and `records` (plain VCF output only).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_plugins_p06_cli_functions import bcftools_contrast_groups

output = Path(tempfile.mkdtemp()) / "contrast.vcf"
result = bcftools_contrast_groups("trio.vcf", output, control_samples=["F", "M"], case_samples=["C"])
print(result["records"])
for line in output.read_text().splitlines()[-40:-37]:
    print(line.split("\t")[1], line.split("\t")[7])
```

- **Example output:**

```text
14
```

## 14) `dnakit.variants.bcftools_fill_from_fasta()`

- **Purpose:** Fill the REF column or an INFO tag from a FASTA file (`bcftools +fill-from-fasta`).
- **Method:** `column` is `REF` (repair REF alleles that are N or wrong) or an INFO tag such as `AA` for an ancestral allele. The FASTA must be indexed (`samtools faidx`). `header_lines_path` holds the `##INFO` line for a new tag (`-h`). `include` / `exclude` (plugin options `-i`, `-e`) limit the annotated records, `replace_non_acgtn` replaces other characters with N (`-N`).
- **API:** `dnakit.variants.bcftools_fill_from_fasta(vcf_path, output_path, *, fasta_path, column, header_lines_path=None, include=None, exclude=None, replace_non_acgtn=False, output_type=None, bcftools_path=None)`.
- **Input:** A VCF/BCF, an indexed FASTA, a column name and a new output path.
- **Output:** The common fields, `output_type` and `records` (plain VCF output only).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_plugins_p06_cli_functions import bcftools_fill_from_fasta

import subprocess

subprocess.run(["samtools", "faidx", "genome.fa"], check=True)
Path("aa.hdr").write_text('##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral allele">\n')
output = Path(tempfile.mkdtemp()) / "aa.vcf"
result = bcftools_fill_from_fasta(
    "variants.vcf", output, fasta_path="genome.fa", column="AA", header_lines_path="aa.hdr", include='TYPE="snp"'
)
print(result["records"])
for line in output.read_text().splitlines()[-10:-7]:
    print(line.split("\t")[1], line.split("\t")[3], line.split("\t")[7])
```

- **Example output:**

```text
10
60 G DP=50;AF=0.417;AA=G
120 C DP=55;AF=0.333;AA=C
180 C DP=60;AF=0.5;AA=C
```

## 15) `dnakit.variants.bcftools_fill_tags()`

- **Purpose:** Compute INFO and FORMAT tags from the genotypes (AF, AC, AN, MAF, HWE, ExcHet, NS, F_MISSING, VAF, and custom expressions) (`bcftools +fill-tags`).
- **Method:** `tags` (`-t`) lists what is filled: `AF`, `AC`, `AC_Hom`, `AC_Het`, `AC_Hemi`, `AN`, `ExcHet`, `HWE`, `MAF`, `NS`, `F_MISSING`, `END`, `TYPE`, `VAF`, `VAF1`, `all`, or custom expressions such as `DP:1=int(sum(FORMAT/DP))`. `populations_file` (`-S`) lists samples with comma-separated population names, to compute the tags per population. `drop_missing` does not count half-missing genotypes as hemizygous (`-d`).
- **API:** `dnakit.variants.bcftools_fill_tags(vcf_path, output_path, *, tags, populations_file=None, drop_missing=False, regions=None, include=None, exclude=None, output_type=None, bcftools_path=None)`.
- **Input:** A VCF/BCF, a list of tags and a new output path.
- **Output:** The common fields, `output_type` and `records` (plain VCF output only).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_plugins_p06_cli_functions import bcftools_fill_tags

output = Path(tempfile.mkdtemp()) / "tags.vcf"
result = bcftools_fill_tags("variants.vcf", output, tags=["AN", "AC", "NS", "F_MISSING"])
print(result["records"])
for line in output.read_text().splitlines()[-10:-8]:
    print(line.split("\t")[1], line.split("\t")[7])
```

- **Example output:**

```text
10
60 DP=50;AF=0.417;F_MISSING=0;NS=6;AN=12;AC=5
120 DP=55;AF=0.333;F_MISSING=0;NS=6;AN=12;AC=4
```

## 16) `dnakit.variants.bcftools_exclusive_genotype_sites()`

- **Purpose:** Output only the sites where the requested samples all share one genotype and no other sample has it (`bcftools +GTsubset`).
- **Method:** Run the plugin with the comma-joined `samples` (`-s`).
- **API:** `dnakit.variants.bcftools_exclusive_genotype_sites(vcf_path, output_path, *, samples, regions=None, output_type=None, bcftools_path=None)`.
- **Input:** A multi-sample VCF/BCF, a list of samples and a new output path.
- **Output:** The common fields, `output_type` and `records` (plain VCF output only).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_plugins_p06_cli_functions import bcftools_exclusive_genotype_sites

output = Path(tempfile.mkdtemp()) / "exclusive.vcf"
result = bcftools_exclusive_genotype_sites("variants.vcf", output, samples=["S1", "S2"])
print(result["records"])
for line in output.read_text().splitlines():
    if not line.startswith("#"):
        print(line.split("\t")[1], line.split("\t")[9:15])
```

- **Example output:**

```text
1
300 ['0/0', '0/0', '0/1', '0/1', '0/1', '1/1']
```

## 17) `dnakit.variants.bcftools_isec_genotypes()`

- **Purpose:** Compare two files and set the genotypes of the first file that differ from the second to missing (`bcftools +isecGT`).
- **Method:** Run the plugin on the two files; both must be bgzip-compressed or BCF and indexed, and samples and sites are matched by name and position.
- **API:** `dnakit.variants.bcftools_isec_genotypes(vcf_path, other_vcf_path, output_path, *, regions=None, output_type=None, bcftools_path=None)`.
- **Input:** Two indexed VCF/BCF files and a new output path.
- **Output:** The common fields, `output_type` and `records` (plain VCF output only).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_plugins_p06_cli_functions import bcftools_isec_genotypes

import subprocess

work = Path(tempfile.mkdtemp())
text = Path("trio.vcf").read_text().replace("0/1:60,0,60:15,15", "1/1:60,0,60:15,15", 3)
(work / "other.vcf").write_text(text)
for name, source in (("first", "trio.vcf"), ("other", str(work / "other.vcf"))):
    subprocess.run(["bcftools", "view", "-Oz", "-o", str(work / f"{name}.vcf.gz"), source], check=True)
    subprocess.run(["bcftools", "index", str(work / f"{name}.vcf.gz")], check=True)
output = work / "isec.vcf"
result = bcftools_isec_genotypes(work / "first.vcf.gz", work / "other.vcf.gz", output)
rows = [line.split("\t") for line in output.read_text().splitlines() if not line.startswith("#")]
print(result["records"], sum(call.startswith("./.") for row in rows for call in row[9:]))
```

- **Example output:**

```text
14 3
```
