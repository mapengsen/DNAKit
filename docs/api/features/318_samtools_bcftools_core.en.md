# samtools and bcftools command-line facades: FASTQ regions, coverage histograms, checksums, convert, index, norm, polysomy

This group of functions in `dnakit.alignment` (module `samtools_p06_cli_functions`) and `dnakit.variants` (module `bcftools_p06_cli_functions`) wraps commands of the external programs [samtools](https://github.com/samtools/samtools) and [bcftools](https://github.com/samtools/bcftools) (1.24 is the tested version). Each function validates its arguments, runs one command with a fixed argument list and no shell, captures standard output and standard error separately, and returns a JSON-serialisable dictionary. Outputs go to explicit new paths; an existing file is never overwritten, and a failed run leaves no partial output behind.

samtools and bcftools are external programs: install them (for example `conda install -c bioconda samtools bcftools`, or a release from https://www.htslib.org/download/), then put them on PATH or pass `samtools_path` / `bcftools_path`. If a program cannot be found, `BackendUnavailableError` is raised with the install hint; if it exits with an error, `BackendExecutionError` (code `HTSLIB_EXECUTION_FAILED`) is raised with the tail of its log. Both programs are MIT licensed and are not bundled with DNAKit. `bcftools polysomy` is a GPL feature that exists only in builds made with GSL.

All results share `status` (`completed`), `backend`, `operation`, `version` (read from the program), `command`, `outputs` (the files written), `log_tail` (the tail of standard error) and `provenance` (tool, version, tested version, executable, citation); the sections below list only the extra fields. The pysam-bundled wrappers for `samtools consensus`, `flagstat`, `stats`, `ampliconstats`, `phase`, `reference -e` and for `bcftools call -c` are the existing functions in `dnakit.bioinformatics.pysam_gap_functions`, and the plugins are in the two documents that follow.

## Example data

The examples read these files from the working directory and write their outputs to a temporary directory. Each run also creates the index next to the input or in the temporary directory.

**Example file:** `reads.fastq`

```text
@read1
TTTCCTCATGCAATTCAAAACCATGTCCGTAATGTAGGCG
+
55F5D5?D55I5::5DDIDD5FI:I??5??5D5::I555D
@read2
TCCTCCTTATTCAGGACCTAACCTGAGGTAAACCAGGTCT
+
F:D:I:?:F:I::IIF:IDDF5D5555F?:IID?DFD?F:
@read3
ACCTAGCCAAGTTCAACGGCAGCTGCAATGGAAATAGGCA
+
5DF?I5III:??5ID5D5DF5FIDDF5F555I5?DI?DII
@read4
TTTAAGATACATTGAGGCCCGTTCGTGCTCCTCGCCCTGA
+
I5?:5DD?:DDIFDI?IFD?55?F5II?F??IF5I:DD:5
@read5
GCCAATAGACCTGCATACCGGCTCATTCTTCATGTGCAAC
+
:D5?F??I5?F55D?FD?D5:5D5:F?I:D:FIFIIDDFI
@read6
ACTGCGGTCGCGTCTAATAATATACATTTGCTTCGTTGAC
+
DF5?:I55:::5I??I?:FD5DIF5FF?I:IFD5DIDFF:
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

**Example file:** `raw.tsv`

```text
rs1	chr1	60	GT
rs2	chr1	120	CC
rs3	chr1	180	--
```

**Example file:** `sample.g.vcf`

```text
##fileformat=VCFv4.2
##INFO=<ID=END,Number=1,Type=Integer,Description="End">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=chr1,length=600>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
chr1	1	.	G	<*>	.	.	END=5	GT	0/0
chr1	6	.	A	C,<*>	50	.	.	GT	1/1
```

**Example file:** `multi.vcf`

```text
##fileformat=VCFv4.2
##contig=<ID=chr1,length=600>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
chr1	10	.	A	C,G	.	.	.	GT	1/2
chr1	20	.	CA	GC	.	.	.	GT	0/1
chr1	30	.	A	C	.	.	.	GT	0/1
chr1	30	.	A	C	.	.	.	GT	0/1
```

**Example file:** `variants.vcf`

```text
##fileformat=VCFv4.2
##contig=<ID=chr1,length=600>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4	S5	S6
chr1	60	var1	G	T	40	PASS	DP=50	GT	0/1	0/0	1/1	0/1	0/0	0/1
chr1	120	var2	C	T	41	PASS	DP=55	GT	0/0	0/1	0/1	1/1	0/0	0/0
chr1	180	var3	C	A	42	PASS	DP=60	GT	1/1	0/1	0/0	0/1	0/1	./.
chr1	240	var4	ACT	A	43	PASS	DP=65	GT	0/1	0/1	0/0	0/0	1/1	0/1
chr1	300	var5	T	C	44	PASS	DP=70	GT	0/0	0/0	0/1	0/1	0/1	1/1
chr1	360	var6	A	T	45	PASS	DP=75	GT	0/1	1/1	0/1	0/0	0/0	0/1
chr1	420	var7	G	GCG	46	PASS	DP=80	GT	0/0	0/1	0/0	0/1	1/1	0/1
chr1	480	var8	T	C	47	PASS	DP=85	GT	0/1	0/0	1/1	1/1	0/1	0/0
chr1	520	var9	G	C	48	PASS	DP=90	GT	0/0	0/0	0/1	0/0	0/1	0/1
chr1	570	var10	A	C	49	PASS	DP=95	GT	1/1	0/1	0/1	0/0	0/0	0/1
```

## 1) `dnakit.alignment.samtools_extract_fastq_regions()`

- **Purpose:** Extract whole records or sub-regions from a FASTQ file through its index (the work-list rows `samtools faidx -f REGION` and `samtools fqidx REGION`, which run the same code).
- **Method:** Run `samtools fqidx` with a fixed argument list. The FASTQ index (`.fai`) is created next to the input file, or at `index_path` (use it when the folder is read-only). Regions are `name` or `name:start-end` (1-based, inclusive); `reverse_complement` reverse-complements the sequence and reverses the quality string, `mark_strand` is `rc`, `no` or `sign`, and `continue_on_missing` keeps going after a region that cannot be retrieved. The output file must not exist.
- **API:** `dnakit.alignment.samtools_extract_fastq_regions(fastq_path, output_path, *, regions=(), region_file=None, line_length=None, reverse_complement=False, mark_strand=None, continue_on_missing=False, index_path=None, samtools_path=None)`.
- **Input:** A FASTQ file, a new output path, and `regions` and/or a `region_file` with one region per line.
- **Output:** `status`, `backend`, `version`, `command`, `outputs` (`output`), `log_tail`, `provenance`, and `records`: the name and sequence length of every FASTQ record written.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.alignment.samtools_p06_cli_functions import samtools_extract_fastq_regions

output = Path(tempfile.mkdtemp()) / "regions.fq"
result = samtools_extract_fastq_regions(
    "reads.fastq", output, regions=["read2:5-14", "read3"], reverse_complement=True, mark_strand="sign"
)
print(result["status"], [(row["name"], row["length"]) for row in result["records"]])
print(output.read_text().splitlines()[1])
```

- **Example output:**

```text
completed [('read2:5-14(-)', 10), ('read3(-)', 40)]
CTGAATAAGG
```

## 2) `dnakit.alignment.samtools_coverage_histogram()`

- **Purpose:** Draw the coverage of each reference sequence as a text histogram (`samtools coverage -m`).
- **Method:** Run `samtools coverage -m` and save the panels. `bins` is the number of histogram columns (`-w`); samtools itself uses the terminal width minus 40, which depends on the output when no terminal is attached, so the facade fixes it to 40 (an 80-column terminal). `ascii_only` replaces the block characters with ASCII (`-A`). `region`, the read and quality filters and the flag masks (`--rf`, `--ff`, a number or flag names such as `UNMAP,SECONDARY`) are passed through unchanged.
- **API:** `dnakit.alignment.samtools_coverage_histogram(alignment_paths, output_path, *, region=None, bins=40, ascii_only=False, min_read_length=0, min_mapping_quality=0, min_base_quality=0, required_flags=None, excluded_flags=None, max_depth=1000000, min_depth=1, reference_fasta=None, samtools_path=None)`.
- **Input:** One or more indexed BAM/CRAM files (the region option needs the index) and a new output path.
- **Output:** The common fields and `contigs`: for every reference sequence its name, length text and `statistics` (number of reads, covered bases, percent covered, mean coverage, mean base and mapping quality, histogram bin width and maximum bin).
- **Example code:**

```python
import subprocess
import tempfile
from pathlib import Path

from dnakit.alignment.samtools_p06_cli_functions import samtools_coverage_histogram

work = Path(tempfile.mkdtemp())
subprocess.run(["samtools", "view", "-b", "-o", str(work / "reads.bam"), "reads.sam"], check=True)
result = samtools_coverage_histogram([work / "reads.bam"], work / "histogram.txt", ascii_only=True, bins=20)
for contig in result["contigs"]:
    stats = contig["statistics"]
    print(contig["name"], contig["length"], stats["Number of reads"], stats["Percent covered"])
print((work / "histogram.txt").read_text().splitlines()[1])
```

- **Example output:**

```text
chr1 600bp 8 55
chr2 400bp 4 50
>  90.00% |:::::::::::         | Number of reads: 8
```

## 3) `dnakit.alignment.samtools_checksum()`

- **Purpose:** Order-agnostic checksums of the sequence, quality, name and tag content of alignment files, and the merge of checksum reports (the work-list rows `samtools checksum` and `samtools checksum -m`).
- **Method:** Run `samtools checksum`. The checksum does not depend on record order unless `in_order` is set (`-O`); `check_position`, `check_cigar` and `check_mate` add CHR/POS, MAPQ/CIGAR and PNEXT/RNEXT/TLEN (`-P -C -M`). `exclude_flags`, `require_flags` and `flag_mask` are `-F`, `-f` and `-b`; `tags` is `-t`; `max_records` is `-N`. With `merge=True` (`-m`) the inputs are checksum reports written earlier by this function and the merged report is written; the report options must match those used for the inputs.
- **API:** `dnakit.alignment.samtools_checksum(input_paths, output_path, *, merge=False, exclude_flags=None, require_flags=None, flag_mask=None, no_reverse_complement=False, tags=None, in_order=False, check_position=False, check_cigar=False, check_mate=False, max_records=0, show_qc=False, tab_separated=False, bamseqchksum_format=False, samtools_path=None)`.
- **Input:** SAM/BAM/CRAM files (or, with `merge=True`, report files) and a new output path.
- **Output:** The common fields, `merged`, and `groups`: one row per group with `Group`, `QC`, `count` and the checksums `flag+seq`, `+name`, `+qual`, `+aux` and `combined`.
- **Example code:**

```python
import subprocess
import tempfile
from pathlib import Path

from dnakit.alignment.samtools_p06_cli_functions import samtools_checksum

work = Path(tempfile.mkdtemp())
bam = work / "reads.bam"
by_name = work / "by_name.bam"
subprocess.run(["samtools", "view", "-b", "-o", str(bam), "reads.sam"], check=True)
subprocess.run(["samtools", "sort", "-n", "-o", str(by_name), str(bam)], check=True)
first = samtools_checksum([bam], work / "a.chk")
second = samtools_checksum([by_name], work / "b.chk")
print(first["groups"][0]["count"], first["groups"][0]["combined"] == second["groups"][0]["combined"])
merged = samtools_checksum([work / "a.chk", work / "b.chk"], work / "merged.chk", merge=True)
print(merged["merged"], merged["groups"][0]["count"])
```

- **Example output:**

```text
12 True
True 24
```

## 4) `dnakit.variants.bcftools_convert_to_vcf()`

- **Purpose:** Convert genotype text formats to VCF/BCF: gen/sample, hap/sample, hap/legend/sample, TSV genotype tables such as 23andMe exports, and gVCF reference blocks (`bcftools convert -G`, `--hapsample2vcf`, `-H`, `--tsv2vcf`, `--gvcf2vcf`).
- **Method:** `input_format` selects the source and the files `input_files` must list in this order: `gensample` (gen, sample), `hapsample` (hap, sample), `haplegendsample` (hap, legend, sample), `tsv` (one table; needs `reference_fasta` and `sample_names`; `columns` as in `-c`, default `ID,CHROM,POS,AA`) and `gvcf` (one gVCF whose reference blocks are expanded to one record per position; needs `reference_fasta`). `haploid_to_diploid` turns haploid hap genotypes into diploid homozygotes. `output_type` is `v`, `z`, `b` or `u` (default: from the file name).
- **API:** `dnakit.variants.bcftools_convert_to_vcf(input_format, input_files, output_path, *, reference_fasta=None, columns=None, sample_names=None, haploid_to_diploid=False, regions=None, include=None, exclude=None, output_type=None, bcftools_path=None)`.
- **Input:** `input_format`, the input files, a new output path and, for `tsv` and `gvcf`, a reference FASTA.
- **Output:** The common fields, `input_format`, `output_type`, `records` (plain VCF output only) and `summary`: the counters bcftools prints (for TSV input rows total, rows skipped, sites written, missing GTs, hom RR, het RA, hom AA, het AA).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_p06_cli_functions import bcftools_convert_to_vcf

output = Path(tempfile.mkdtemp()) / "sample.vcf"
result = bcftools_convert_to_vcf(
    "tsv", ["raw.tsv"], output, reference_fasta="genome.fa", sample_names=["SAMP1"]
)
print(result["records"], result["summary"]["Missing GTs"])
print([line.split("\t")[9] for line in output.read_text().splitlines() if not line.startswith("#")])
```

- **Example output:**

```text
3 1
['0/1', '0/0', './.']
```

- **Example code:**

```python
import subprocess
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_p06_cli_functions import bcftools_convert_to_vcf

work = Path(tempfile.mkdtemp())
subprocess.run(["bcftools", "convert", "-g", str(work / "panel"), "variants.vcf"], check=True, capture_output=True)
result = bcftools_convert_to_vcf(
    "gensample", [work / "panel.gen.gz", work / "panel.samples"], work / "back.vcf"
)
print(result["records"], result["summary"])
expanded = bcftools_convert_to_vcf("gvcf", ["sample.g.vcf"], work / "expanded.vcf", reference_fasta="genome.fa")
print(expanded["records"])
```

- **Example output:**

```text
10 {'Number of processed rows': 10}
6
```

## 5) `dnakit.variants.bcftools_index_records()`

- **Purpose:** Read the number of records per contig, or the total, from the index of a VCF/BCF (`bcftools index -s` and `bcftools index -n`).
- **Method:** Run `bcftools index -s` (`per_contig=True`, with `-a` when `include_empty_contigs`) or `bcftools index -n`. The file needs an existing CSI or TBI index; nothing is rewritten.
- **API:** `dnakit.variants.bcftools_index_records(vcf_path, *, per_contig=True, include_empty_contigs=False, bcftools_path=None)`.
- **Input:** A bgzip-compressed VCF or BCF with an index.
- **Output:** The common fields and `total_records`; with `per_contig` also `contigs`, one row per contig with `contig`, `length` (None when the header gives none) and `records`.
- **Example code:**

```python
import subprocess
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_p06_cli_functions import bcftools_index_records

work = Path(tempfile.mkdtemp())
compressed = work / "variants.vcf.gz"
subprocess.run(["bcftools", "view", "-Oz", "-o", str(compressed), "variants.vcf"], check=True)
subprocess.run(["bcftools", "index", str(compressed)], check=True)
result = bcftools_index_records(compressed)
print(result["contigs"], result["total_records"])
print(bcftools_index_records(compressed, per_contig=False)["total_records"])
```

- **Example output:**

```text
[{'contig': 'chr1', 'length': 600, 'records': 10}] 10
10
```

## 6) `dnakit.variants.bcftools_norm_variants()`

- **Purpose:** Normalise a VCF/BCF: left-align indels, split or join multiallelic records, atomize complex variants, remove duplicates and check REF alleles (`bcftools norm`). It completes `bcftools_normalize_vcf`, which only offers `-f` and `-m`.
- **Method:** One function for all the modes of `bcftools norm`. `reference_fasta` (`-f`) enables left-alignment and the REF check; `check_ref` is `e`, `w`, `x` or `s` (`-c`); `multiallelics` is `-snps`, `-indels`, `-both`, `-any` to split or `+snps`, `+indels`, `+both`, `+any` to join (`-m`); `atomize` is `-a` with `atom_overlaps` `*` or `.`; `remove_duplicates` is `snps`, `indels`, `both`, `all` or `exact` (`-d`); `no_realign` is `-N`; `gff_annotation` right-aligns in transcripts (`-g`). At least one action must be requested.
- **API:** `dnakit.variants.bcftools_norm_variants(vcf_path, output_path, *, reference_fasta=None, multiallelics=None, atomize=False, atom_overlaps='*', check_ref=None, remove_duplicates=None, no_realign=False, gff_annotation=None, old_record_tag=None, strict_filter=False, site_window=1000, regions=None, include=None, exclude=None, output_type=None, bcftools_path=None)`.
- **Input:** A VCF/BCF (bgzip and an index only for `regions`), a new output path and the options above.
- **Output:** The common fields, `output_type`, `records` (plain VCF output only) and `summary` (counters printed on standard error).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_p06_cli_functions import bcftools_norm_variants

work = Path(tempfile.mkdtemp())
result = bcftools_norm_variants("multi.vcf", work / "split.vcf", multiallelics="-any")
print(result["records"])
atoms = bcftools_norm_variants("multi.vcf", work / "atoms.vcf", atomize=True, remove_duplicates="exact")
print(atoms["records"])
for line in (work / "atoms.vcf").read_text().splitlines():
    if not line.startswith("#"):
        print(*line.split("\t")[1:5])
```

- **Example output:**

```text
5
5
10 . A C,*
10 . A G,*
20 . C G
21 . A C
30 . A C
```

## 7) `dnakit.variants.bcftools_polysomy()`

- **Purpose:** Detect whole-chromosome copy-number aberrations from B-allele frequencies (`bcftools polysomy`).
- **Method:** Run `bcftools polysomy`. The VCF/BCF needs `FORMAT/BAF`. For each chromosome the peaks of the BAF distribution are fitted for 2, 3 and 4 copies and the best fit is the estimated copy number. `peak_size`, `cn_penalty`, `fit_threshold`, `min_fraction`, `peak_symmetry` and `include_aa` are `-b -c -f -m -p -i`. The command exists only in bcftools built with GSL (`--enable-libgsl`; the conda builds are); otherwise `BackendUnavailableError` is raised. The plot script `dist.py` is written but not run.
- **API:** `dnakit.variants.bcftools_polysomy(vcf_path, output_dir, *, sample=None, regions=None, targets=None, peak_size=0.1, cn_penalty=0.7, fit_threshold=3.3, min_fraction=0.1, peak_symmetry=0.5, include_aa=False, bcftools_path=None)`.
- **Input:** A VCF/BCF with `FORMAT/BAF` and a new output folder.
- **Output:** The common fields, `files` (`dist.dat` with the BAF distributions, fitted functions and `CN` lines, and the plotting script `dist.py`) and `copy_numbers`: chromosome, estimated copy number and absolute fit deviation.
- **Example code:**

```python
import random
import subprocess
import tempfile
from pathlib import Path

from dnakit.variants.bcftools_p06_cli_functions import bcftools_polysomy

rng = random.Random(3)
lines = [
    "##fileformat=VCFv4.2", "##contig=<ID=1,length=1000000>", "##contig=<ID=2,length=1000000>",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=BAF,Number=1,Type=Float,Description="B allele frequency">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
]
for chrom, centers in (("1", [0.0, 0.5, 1.0]), ("2", [0.0, 1 / 3, 2 / 3, 1.0])):  # two and three copies
    for index in range(2500):
        center = rng.choice(centers)
        value = abs(rng.gauss(0, 0.02)) if center == 0 else 1 - abs(rng.gauss(0, 0.02)) if center == 1 else rng.gauss(center, 0.025)
        genotype = "0/0" if center == 0 else "1/1" if center == 1 else "0/1"
        lines.append(f"{chrom}\t{1000 + index * 100}\t.\tA\tG\t.\t.\t.\tGT:BAF\t{genotype}:{min(1.0, max(0.0, value)):.4f}")
work = Path(tempfile.mkdtemp())
(work / "baf.vcf").write_text("\n".join(lines) + "\n")
subprocess.run(["bcftools", "view", "-Oz", "-o", str(work / "baf.vcf.gz"), str(work / "baf.vcf")], check=True)
subprocess.run(["bcftools", "index", str(work / "baf.vcf.gz")], check=True)
result = bcftools_polysomy(work / "baf.vcf.gz", work / "polysomy", sample="S1")
print([(row["chromosome"], row["copy_number"]) for row in result["copy_numbers"]], sorted(result["files"]))
```

- **Example output:**

```text
[('1', 2.0), ('2', 3.0)] ['dist.dat', 'dist.py']
```
