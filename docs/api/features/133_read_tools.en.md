# Read tools: pysam, pyfaidx, cutadapt and UMI-tools

This group of functions wraps four packages that work on sequencing reads and FASTA files: [pysam](https://github.com/pysam-developers/pysam) (with the samtools and bcftools code that it bundles), [pyfaidx](https://github.com/mdshw5/pyfaidx), [cutadapt](https://github.com/marcelm/cutadapt) and [UMI-tools](https://github.com/CGATOxford/UMI-tools). The pysam functions convert quality strings, read tabix-style files, run samtools and bcftools commands (depadding, FASTQ export, collating, CRAM sizes, copy-number calling) and render a text alignment view. The pyfaidx functions parse regions and BED lines, and fetch, mask, tabulate and split FASTA sequences. The cutadapt functions trim adapters, quality ends, fixed numbers of bases and N bases from FASTQ reads, with several report formats. The UMI-tools functions move UMIs into read names, estimate cell whitelists, group and deduplicate BAM files by UMI, count UMIs per gene and prepare BAM files for RSEM. Import them from `dnakit.alignment`, `dnakit.sequence`, `dnakit.fasta_fastq` and `dnakit.sequencing_qc`, for example `from dnakit.alignment import qualitystring_to_array`.

The licences of the installed packages (from their package metadata) are MIT for pysam, cutadapt and UMI-tools and BSD-3-Clause for pyfaidx; none of them is GPL.

The dependencies are optional and are imported only when a function is called. If a package or program is missing, the function raises an error that includes the install hint. The examples call every function with `verbose=False`, write their outputs into a new temporary folder, and read the inputs from the example files below.

## Example data

The examples use small files that are declared below: a SAM file (`reads.sam`) and its reference (`ref.fa`), a BED file, a padded alignment with its padded reference, a VCF with `BAF` and `LRR` values, a FASTA file (`genes.fa`), two FASTQ files, and UMI SAM files. Examples that need a BAM or CRAM file make it from the SAM file with `pysam.view()` inside the example.

**Example file:** `reads.sam`

```text
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:40
r1	0	chr1	5	60	10M	*	0	0	ACGTACGTTG	IIIIIHHHHH
r2	16	chr1	12	60	8M	*	0	0	TTGACCGA	FFFFFFFF
r3	0	chr1	20	60	6M	*	0	0	GGGCCC	IIIIII
```

**Example file:** `ref.fa`

```text
>chr1
ACGTACGTACGTTGACCGAGGGCCCAAATTTGGGCCCAAA
```

**Example file:** `regions.bed`

```text
g1	2	10
g2	0	8
```

**Example file:** `padded.fa`

```text
>chr1
ACGT**ACGTAC
```

**Example file:** `pad.sam`

```text
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:12
r1	0	chr1	1	60	4M2D6M	*	0	0	ACGTACGTAC	IIIIIIIIII
r2	0	chr1	3	60	2M2D4M	*	0	0	GTCGTA	IIIIII
```

**Example file:** `cnv.vcf`

```text
##fileformat=VCFv4.2
##contig=<ID=chr1,length=100000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=BAF,Number=1,Type=Float,Description="B-allele frequency">
##FORMAT=<ID=LRR,Number=1,Type=Float,Description="Log R ratio">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
chr1	1000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.48:0.02
chr1	2000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.52:-0.02
chr1	3000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.48:0.02
chr1	4000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.52:-0.02
chr1	5000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.48:0.02
chr1	6000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.52:-0.02
chr1	7000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.48:0.02
chr1	8000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.52:-0.02
chr1	9000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.48:0.02
chr1	10000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.52:-0.02
chr1	11000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.48:0.02
chr1	12000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.52:-0.02
chr1	13000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.33:0.4
chr1	14000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.67:0.4
chr1	15000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.33:0.4
chr1	16000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.67:0.4
chr1	17000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.33:0.4
chr1	18000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.67:0.4
chr1	19000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.33:0.4
chr1	20000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.67:0.4
chr1	21000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.33:0.4
chr1	22000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.67:0.4
chr1	23000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.33:0.4
chr1	24000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.67:0.4
chr1	25000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.48:0.02
chr1	26000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.52:-0.02
chr1	27000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.48:0.02
chr1	28000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.52:-0.02
chr1	29000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.48:0.02
chr1	30000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.52:-0.02
chr1	31000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.48:0.02
chr1	32000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.52:-0.02
chr1	33000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.48:0.02
chr1	34000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.52:-0.02
chr1	35000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.48:0.02
chr1	36000	.	A	G	.	.	.	GT:BAF:LRR	0/1:0.52:-0.02
```

**Example file:** `genes.fa`

```text
>g1
ACGTACGTACGTTGACCGAGGGCCC
>g2
TTTTGGGGCCCCAAAANNACGTACGT
>g3
GATTACA
```

**Example file:** `reads.fq`

```text
@r1
ACGTTGCAAGATCGGAAGAGCACACG
+
IIIIIIIIIIIIIIIIIIIIIIIIII
@r2
TTGACCGATCAGATCGGAAGAGC
+
IIIIIIIIIIIIIIIIIIIIIII
@r3
GGGCCCAAATTTGGGAAA
+
IIIIIIIIIIIIIIIIII
@r4
NNACGTACGTTGAACGNN
+
##IIIIIIIIIIIIII##
```

**Example file:** `umi_reads.fq`

```text
@r1
ACGTAAGGGCCCAAATTTGGG
+
IIIIIIIIIIIIIIIIIIIII
@r2
TTGCCAGGGCCCAAATTTGGG
+
IIIIIIIIIIIIIIIIIIIII
@r3
ACGTAAGGGCCCAAATTTGGG
+
IIIIIIIIIIIIIIIIIIIII
@r4
ACGTTTGGGCCCAAATTTGGG
+
IIIIIIIIIIIIIIIIIIIII
```

**Example file:** `umi.sam`

```text
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:100
a1_AAAA	0	chr1	10	60	10M	*	0	0	ACGTACGTAC	IIIIIIIIII	XT:Z:geneA
a2_AAAA	0	chr1	10	60	10M	*	0	0	ACGTACGTAC	IIIIIIIIII	XT:Z:geneA
a3_AAAT	0	chr1	10	60	10M	*	0	0	ACGTACGTAC	IIIIIIIIII	XT:Z:geneA
a4_CCCC	0	chr1	10	60	10M	*	0	0	ACGTACGTAC	IIIIIIIIII	XT:Z:geneA
a5_GGGG	0	chr1	30	60	10M	*	0	0	ACGTACGTAC	IIIIIIIIII	XT:Z:geneB
a6_GGGG	0	chr1	30	60	10M	*	0	0	ACGTACGTAC	IIIIIIIIII	XT:Z:geneB
```

**Example file:** `pair.sam`

```text
@HD	VN:1.6	SO:queryname
@SQ	SN:tx1	LN:100
p1_AAAA	99	tx1	5	60	10M	=	40	45	ACGTACGTAC	IIIIIIIIII	RX:Z:AAAA
p1_AAAA	147	tx1	40	60	10M	=	5	-45	ACGTACGTAC	IIIIIIIIII
p2_CCCC	99	tx1	10	60	10M	=	50	45	ACGTACGTAC	IIIIIIIIII	RX:Z:CCCC
p2_CCCC	147	tx1	50	60	10M	=	10	-45	ACGTACGTAC	IIIIIIIIII
```

## 1) `dnakit.alignment.array_to_qualitystring()`

- **Purpose:** Convert integer Phred quality values to a quality string.
- **Method:** Wraps the integers in an unsigned-byte array and calls `pysam.array_to_qualitystring()`, which adds `offset` to every value and returns the characters.
- **API:** `array_to_qualitystring(qualities: 'Iterable[int]', *, offset: 'int' = 33, verbose: bool = True) -> 'str'`.
- **Requirements:** Needs the optional package pysam (extra `dnakit[bio-reads]`). The samtools code is bundled with pysam, so no external program is needed.
- **Input:** `qualities` is an iterable of integers (each converted with `int()` and stored as an unsigned byte, so 0 to 255). `offset` is 33 for Sanger/Illumina 1.8 and 64 for the old Illumina encoding.
- **Output:** A string with one character per quality value.
- **Example code:**

```python
from dnakit.alignment import array_to_qualitystring

print(array_to_qualitystring([40, 40, 20, 2], verbose=False))
print(array_to_qualitystring([40, 40, 20, 2], offset=64, verbose=False))
```

- **Example output:**

```text
II5#
hhTB
```

## 2) `dnakit.alignment.bam_to_fastq()`

- **Purpose:** Convert alignments to FASTQ text with the samtools `fastq` command.
- **Method:** Calls `pysam.bam2fq()` with your arguments plus `-0 <output_path>`. pysam would send the output to the process standard output, so the function has samtools write it to `output_path` and then reads that file back.
- **API:** `bam_to_fastq(arguments: 'Sequence[str]', output_path: 'str', *, verbose: bool = True) -> 'str'`.
- **Requirements:** Needs the optional package pysam (extra `dnakit[bio-reads]`). The samtools code is bundled with pysam, so no external program is needed.
- **Input:** `arguments` is a list of strings passed to `samtools fastq` (a string instead of a list, or non-string items, raise `ConfigurationError`); give at least the input BAM, SAM or CRAM path. `output_path` must not exist yet. In samtools, `-0` receives the reads that have neither the READ1 nor the READ2 flag, such as the single-end reads of the example.
- **Output:** The FASTQ text as a string. The same text is in the new file `output_path`.
- **Example code:**

```python
import tempfile
from pathlib import Path

import pysam
from dnakit.alignment import bam_to_fastq

work = Path(tempfile.mkdtemp())
bam = work / "reads.bam"
pysam.view("-b", "-o", str(bam), "reads.sam", catch_stdout=False)
text = bam_to_fastq(["-n", str(bam)], str(work / "reads.fq"), verbose=False)
print(text)
```

- **Example output:**

```text
@r1
ACGTACGTTG
+
IIIIIHHHHH
@r2
TCGGTCAA
+
FFFFFFFF
@r3
GGGCCC
+
IIIIII
```

## 3) `dnakit.alignment.cram_container_sizes()`

- **Purpose:** Report the container and block sizes of a CRAM file.
- **Method:** Calls `pysam.cram_size()` (the samtools `cram_size` command) with your arguments and returns its text.
- **API:** `cram_container_sizes(arguments: 'Sequence[str]', *, verbose: bool = True) -> 'str | list[str]'`.
- **Requirements:** Needs the optional package pysam (extra `dnakit[bio-reads]`). The samtools code is bundled with pysam, so no external program is needed.
- **Input:** `arguments` is a list of strings for `samtools cram_size`, normally the CRAM path and optionally options such as `-v`. A non-CRAM file raises pysam's `SamtoolsError` ("Input is not a CRAM file").
- **Output:** The report as text: a block table (content id, uncompressed and compressed size, ratio, method, data series) followed by the lines `Number of containers`, `Number of slices`, `Number of sequences`, `Number of bases`, `Total file size` and `Format overhead size`.
- **Example code:**

```python
import tempfile
from pathlib import Path

import pysam
from dnakit.alignment import cram_container_sizes

work = Path(tempfile.mkdtemp())
cram = work / "reads.cram"
pysam.view("-C", "-T", "ref.fa", "-o", str(cram), "reads.sam", catch_stdout=False)
report = cram_container_sizes([str(cram)], verbose=False)
for line in report.splitlines():
    if line.startswith("Number of"):
        print(line)
```

- **Example output:**

```text
Number of containers                   1
Number of slices                       1
Number of sequences                    3
Number of bases                       24
```

## 4) `dnakit.alignment.pad_to_unpadded()`

- **Purpose:** Convert a padded alignment to unpadded coordinates with the legacy samtools name `pad2unpad` and your own arguments.
- **Method:** Calls `pysam.pad2unpad()` with your arguments. This pysam wrapper does not capture standard output: without `-o` the alignments go to the process standard output and the function returns an empty string; with `-o <file>` samtools writes the file and the function also returns an empty string.
- **API:** `pad_to_unpadded(arguments: 'Sequence[str]', *, verbose: bool = True) -> 'str | list[str]'`.
- **Requirements:** Needs the optional package pysam (extra `dnakit[bio-reads]`). The samtools code is bundled with pysam, so no external program is needed.
- **Input:** `arguments` is a list of strings for `samtools pad2unpad`: the padded reference is given with `-T <fasta>`, `-s` asks for SAM, `-o <file>` names the output, and the input alignment comes last. Reads must not use P or I CIGAR operations.
- **Output:** An empty string in the usual case. Read the result from the `-o` file. For a call that returns the data and for a structured result, use `remove_padding()` or `pad_to_unpad_alignments()`.
- **Example code:**

```python
import tempfile
from pathlib import Path

import pysam
from dnakit.alignment import pad_to_unpadded

work = Path(tempfile.mkdtemp())
bam = work / "pad.bam"
pysam.view("-b", "-o", str(bam), "pad.sam", catch_stdout=False)
out = work / "unpadded.sam"
returned = pad_to_unpadded(["-T", "padded.fa", "-s", "-o", str(out), str(bam)], verbose=False)
print(repr(returned))
for line in out.read_text().splitlines():
    if not line.startswith("@"):
        print(line)
```

- **Example output:**

```text
''
r1	0	chr1	1	60	10M	*	0	0	ACGTACGTAC	IIIIIIIIII
r2	0	chr1	3	60	6M	*	0	0	GTCGTA	IIIIII
```

## 5) `dnakit.alignment.qualities_to_qualitystring()`

- **Purpose:** Convert a list of quality scores to a quality string.
- **Method:** Wraps the integers in an unsigned-byte array and calls `pysam.qualities_to_qualitystring()`. For integer input it gives the same string as `array_to_qualitystring()`.
- **API:** `qualities_to_qualitystring(qualities: 'Iterable[int]', *, offset: 'int' = 33, verbose: bool = True) -> 'str'`.
- **Requirements:** Needs the optional package pysam (extra `dnakit[bio-reads]`). The samtools code is bundled with pysam, so no external program is needed.
- **Input:** `qualities` is an iterable of integers (any iterable, for example a list or a `range`; each is converted with `int()` and stored as an unsigned byte). `offset` defaults to 33.
- **Output:** A string with one character per quality value.
- **Example code:**

```python
from dnakit.alignment import qualities_to_qualitystring

print(qualities_to_qualitystring([40, 40, 20, 2], verbose=False))
print(qualities_to_qualitystring(range(0, 41, 10), offset=64, verbose=False))
```

- **Example output:**

```text
II5#
@JT^h
```

## 6) `dnakit.alignment.qualitystring_to_array()`

- **Purpose:** Convert a Phred quality string to a list of integer quality values.
- **Method:** Calls `pysam.qualitystring_to_array()`, which subtracts `offset` from the code of every character, and converts the result to a list of Python integers.
- **API:** `qualitystring_to_array(quality: 'str', *, offset: 'int' = 33, verbose: bool = True) -> 'list[int]'`.
- **Requirements:** Needs the optional package pysam (extra `dnakit[bio-reads]`). The samtools code is bundled with pysam, so no external program is needed.
- **Input:** `quality` must be a string (otherwise `ConfigurationError`). `offset` defaults to 33; use 64 for the old Illumina encoding.
- **Output:** A list of integers, one per character.
- **Example code:**

```python
from dnakit.alignment import qualitystring_to_array

print(qualitystring_to_array("II5#", verbose=False))
print(qualitystring_to_array("hhTB", offset=64, verbose=False))
```

- **Example output:**

```text
[40, 40, 20, 2]
[40, 40, 20, 2]
```

## 7) `dnakit.alignment.remove_padding()`

- **Purpose:** Remove alignment padding from a padded alignment, with your own samtools arguments.
- **Method:** Calls `pysam.depad()` (the samtools `depad` command) with your arguments and returns what pysam captured from standard output. pysam adds its own temporary `-o` option, so the returned text is the whole SAM, including a `@PG` header line that contains a temporary file name.
- **API:** `remove_padding(arguments: 'Sequence[str]', *, verbose: bool = True) -> 'str | list[str]'`.
- **Requirements:** Needs the optional package pysam (extra `dnakit[bio-reads]`). The samtools code is bundled with pysam, so no external program is needed.
- **Input:** `arguments` is a list of strings for `samtools depad`: `-T <fasta>` gives the padded reference (pad columns written as `*`), `-s` asks for SAM text, and the input alignment comes last. Do not add `-o`; the captured text is the result. Reads must not use P or I CIGAR operations.
- **Output:** The SAM text (with `-s`) as a string. Without `-s` samtools writes BAM, which is not text.
- **Example code:**

```python
import tempfile
from pathlib import Path

import pysam
from dnakit.alignment import remove_padding

work = Path(tempfile.mkdtemp())
bam = work / "pad.bam"
pysam.view("-b", "-o", str(bam), "pad.sam", catch_stdout=False)
text = remove_padding(["-T", "padded.fa", "-s", str(bam)], verbose=False)
for line in text.splitlines():
    if not line.startswith("@"):
        print(line)
```

- **Example output:**

```text
r1	0	chr1	1	60	10M	*	0	0	ACGTACGTAC	IIIIIIIIII
r2	0	chr1	3	60	6M	*	0	0	GTCGTA	IIIIII
```

## 8) `dnakit.alignment.shuffle_bam()`

- **Purpose:** Run the samtools `bamshuf` command with your own arguments.
- **Method:** Calls `pysam.bamshuf()` with your arguments and returns what pysam captured from standard output. In the bundled samtools 1.21 `bamshuf` is the `collate` command (the usage text it prints is that of `samtools collate`), which groups reads with the same name together.
- **API:** `shuffle_bam(arguments: 'Sequence[str]', *, verbose: bool = True) -> 'str | list[str]'`.
- **Requirements:** Needs the optional package pysam (extra `dnakit[bio-reads]`). The samtools code is bundled with pysam, so no external program is needed.
- **Input:** `arguments` is a list of strings: options such as `-o <file>` (output file), `-n <int>` (number of temporary files) and `-T <prefix>` (temporary file prefix), then the input BAM and, without `-o`, an output prefix. Give `-o` or a prefix, so the result is written to a file; `-O` sends BAM to the process standard output, where it is not returned.
- **Output:** The captured standard output: an empty string when the result goes to a file.
- **Example code:**

```python
import tempfile
from pathlib import Path

import pysam
from dnakit.alignment import shuffle_bam

work = Path(tempfile.mkdtemp())
bam = work / "reads.bam"
pysam.view("-b", "-o", str(bam), "reads.sam", catch_stdout=False)
out = work / "shuffled.bam"
returned = shuffle_bam(
    ["-o", str(out), "-n", "2", "-T", str(work / "tmp"), str(bam)], verbose=False
)
print(repr(returned))
print(sorted(read.query_name for read in pysam.AlignmentFile(str(out))))
```

- **Example output:**

```text
''
['r1', 'r2', 'r3']
```

## 9) `dnakit.alignment.tabix_records()`

- **Purpose:** Read every record of a text file with a pysam tabix parser.
- **Method:** Opens the file in binary mode and returns `list(pysam.tabix_iterator(handle, parser))`.
- **API:** `tabix_records(path: 'str', parser: 'Any', *, verbose: bool = True) -> 'list[Any]'`.
- **Requirements:** Needs the optional package pysam (extra `dnakit[bio-reads]`). The samtools code is bundled with pysam, so no external program is needed.
- **Input:** `path` is an existing file (it is opened directly; the example uses a plain BED file and no index). `parser` is a pysam parser object such as `pysam.asBed()`, `pysam.asGTF()`, `pysam.asVCF()` or `pysam.asTuple()`.
- **Output:** A list of parsed records. With `pysam.asBed()` each record has `contig`, `start` (0-based), `end` and `name`.
- **Example code:**

```python
import pysam
from dnakit.alignment import tabix_records

records = tabix_records("regions.bed", pysam.asBed(), verbose=False)
print(len(records))
for record in records:
    print(record.contig, record.start, record.end)
```

- **Example output:**

```text
2
g1 2 10
g2 0 8
```

## 10) `dnakit.alignment.target_cut()`

- **Purpose:** Run the samtools `targetcut` command with your own arguments.
- **Method:** Calls `pysam.targetcut()` with your arguments and returns what pysam captured from standard output. The command reads a BAM file and prints SAM lines; in the example it prints one pseudo-read for the covered region, built from the base calls and qualities of the reads (positions with no coverage are `N` with quality `!`).
- **API:** `target_cut(arguments: 'Sequence[str]', *, verbose: bool = True) -> 'str | list[str]'`.
- **Requirements:** Needs the optional package pysam (extra `dnakit[bio-reads]`). The samtools code is bundled with pysam, so no external program is needed.
- **Input:** `arguments` is a list of strings for `samtools targetcut`: options such as `-Q <int>` (minimum base quality) followed by the BAM path. If the arguments do not name an input, samtools waits on standard input, so always give the BAM.
- **Output:** The SAM text as a string, one line per cut region.
- **Example code:**

```python
import tempfile
from pathlib import Path

import pysam
from dnakit.alignment import target_cut

work = Path(tempfile.mkdtemp())
bam = work / "reads.bam"
pysam.view("-b", "-o", str(bam), "reads.sam", catch_stdout=False)
text = target_cut([str(bam)], verbose=False)
for line in text.splitlines():
    fields = line.split("\t")
    print(fields[0], fields[1], fields[3], fields[5], fields[9])
```

- **Example output:**

```text
chr1:2-25 0 2 24M NNNACGTACGTTGACCGAGGGCCC
```

## 11) `dnakit.alignment.call_copy_number_states()`

- **Purpose:** Call copy-number states along a chromosome from B-allele frequency and log R ratio values in a VCF.
- **Method:** Runs `bcftools cnv` through `pysam.bcftools.cnv` with `-s <query_sample> -o <output_dir>` (plus `-c`, `-f` and `-r` when given) and `catch_stdout=False`. It then reads `cn.<sample>.tab` (one row per site) and `summary.<sample>.tab` (the `RG` rows are the segments).
- **API:** `call_copy_number_states(vcf: 'str', output_dir: 'str', *, query_sample: 'str', control_sample: 'str | None' = None, allele_frequency_file: 'str | None' = None, regions: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package pysam (extra `dnakit[bio-reads]`). The samtools code is bundled with pysam, so no external program is needed.
- **Input:** `vcf` is an existing VCF with the FORMAT fields `BAF` and `LRR` for the query sample. `output_dir` must not exist yet and its parent folder must exist; the function creates it. `query_sample` is a non-empty name that does not start with `-`. `control_sample` (`-c`), `allele_frequency_file` (`-f`, an existing file) and `regions` (`-r`, a region string that does not start with `-`) are optional.
- **Output:** A dictionary with `status`, `backend`, `command`, `sample`, `outputs` (paths of `cn`, `dat`, `summary` and `plot_script` files that exist), `segments` (each with `chromosome`, `start`, `end`, `copy_number`, `quality`, `n_sites`, `n_hets`) and `sites` (each with `chromosome`, `position`, `copy_number` and four `probabilities`). The `quality` of a segment can be `nan`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.alignment import call_copy_number_states

out_dir = Path(tempfile.mkdtemp()) / "cnv"
result = call_copy_number_states("cnv.vcf", str(out_dir), query_sample="S1", verbose=False)
print(result["status"], result["sample"], sorted(result["outputs"]))
for seg in result["segments"]:
    print(seg["chromosome"], seg["start"], seg["end"], seg["copy_number"], seg["n_sites"])
print(len(result["sites"]), result["sites"][0]["position"], result["sites"][0]["copy_number"])
```

- **Example output:**

```text
completed S1 ['cn', 'dat', 'plot_script', 'summary']
chr1 1000 13999 2 14
chr1 14000 25999 3 12
chr1 26000 36000 2 10
36 1000 2
```

## 12) `dnakit.alignment.depad_alignments()`

- **Purpose:** Convert a padded alignment to unpadded coordinates and return the records.
- **Method:** Runs `pysam.depad()` with `-T <padded_reference> -o <output_alignment> <alignment>` (and `-s` for SAM) and `catch_stdout=False`, then reads the output file back with `pysam.AlignmentFile`.
- **API:** `depad_alignments(alignment: 'str', padded_reference: 'str', output_alignment: 'str', *, output_sam: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package pysam (extra `dnakit[bio-reads]`). The samtools code is bundled with pysam, so no external program is needed.
- **Input:** `alignment` and `padded_reference` (a FASTA whose pad columns are `*`) must be existing files. `output_alignment` must not exist yet, and its folder must exist. `output_sam=True` writes SAM, otherwise BAM. Reads must not use P or I CIGAR operations.
- **Output:** A dictionary with `status`, `backend`, `command`, `output_format` (`BAM` or `SAM`), `outputs`, `record_count` and `records` (at most the first 200, each with `query_name`, `flag`, `reference_name`, 1-based `reference_start`, `cigar`, `sequence`). A samtools failure raises `BackendExecutionError`.
- **Example code:**

```python
import tempfile
from pathlib import Path

import pysam
from dnakit.alignment import depad_alignments

work = Path(tempfile.mkdtemp())
bam = work / "pad.bam"
pysam.view("-b", "-o", str(bam), "pad.sam", catch_stdout=False)
result = depad_alignments(str(bam), "padded.fa", str(work / "unpadded.bam"), verbose=False)
print(result["status"], result["output_format"], result["record_count"])
for read in result["records"]:
    print(read["query_name"], read["reference_start"], read["cigar"], read["sequence"])
```

- **Example output:**

```text
completed BAM 2
r1 1 10M ACGTACGTAC
r2 3 6M GTCGTA
```

## 13) `dnakit.alignment.pad_to_unpad_alignments()`

- **Purpose:** Convert a padded alignment to unpadded coordinates through the legacy samtools command name `pad2unpad` and return the records.
- **Method:** Same code path as `depad_alignments()`, but it calls `pysam.pad2unpad()`. With the bundled samtools the two commands gave the same records for this page's data.
- **API:** `pad_to_unpad_alignments(alignment: 'str', padded_reference: 'str', output_alignment: 'str', *, output_sam: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package pysam (extra `dnakit[bio-reads]`). The samtools code is bundled with pysam, so no external program is needed.
- **Input:** The same as `depad_alignments()`: existing `alignment` and `padded_reference`, a new `output_alignment`, optional `output_sam`.
- **Output:** The same dictionary as `depad_alignments()`; `command` shows `pad2unpad`.
- **Example code:**

```python
import tempfile
from pathlib import Path

import pysam
from dnakit.alignment import pad_to_unpad_alignments

work = Path(tempfile.mkdtemp())
bam = work / "pad.bam"
pysam.view("-b", "-o", str(bam), "pad.sam", catch_stdout=False)
result = pad_to_unpad_alignments(
    str(bam), "padded.fa", str(work / "unpadded.sam"), output_sam=True, verbose=False
)
print(result["status"], result["output_format"], result["record_count"])
print(result["command"][:2])
for read in result["records"]:
    print(read["query_name"], read["reference_start"], read["cigar"])
```

- **Example output:**

```text
completed SAM 2
['samtools', 'pad2unpad']
r1 1 10M
r2 3 6M
```

## 14) `dnakit.alignment.pysam_tview_text()`

- **Purpose:** Render one screen of a BAM or CRAM alignment as text with `samtools tview`.
- **Method:** Runs the samtools found on PATH (or `samtools=`) as `samtools tview -d T -p <chrom:pos> -w <width> <alignment> [reference]`, writes its standard output to the new file `output`, and returns the lines. pysam's bundled samtools has no `tview`, so the function does not use `pysam.tview`. The version comes from `samtools --version`.
- **API:** `pysam_tview_text(alignment: 'str | Path', position: 'str', *, output: 'str | Path', reference: 'str | Path | None' = None, width: 'int' = 80, samtools: 'str | Path' = 'samtools', verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program samtools (the install hint in the source says it was tested with 1.24) on PATH or given with `samtools=`; pysam's bundled copy does not work. Without it the function raises `BackendUnavailableError` with an install hint. This environment has no samtools on PATH, so the example shows that error; the call that runs for real is shown after it.
- **Input:** `alignment` is an existing BAM or CRAM and `position` is `chrom:pos` with a 1-based position. `output` must be a new file whose folder exists, and `<output>.stderr` must not exist. `width` is an integer from 20 to 1000. `reference` is an optional FASTA. The checks run before samtools is looked for, so a missing alignment gives `ConfigurationError` first.
- **Output:** A dictionary with `status`, `backend`, `program`, `version`, `alignment`, `reference`, `position`, `width`, `output`, `line_count`, `lines` (the text view) and `log`.
- **Example code:**

```python
import tempfile
from pathlib import Path

import pysam
from dnakit.alignment import pysam_tview_text
from dnakit.exceptions import BackendUnavailableError

work = Path(tempfile.mkdtemp())
bam = work / "reads.bam"
pysam.view("-b", "-o", str(bam), "reads.sam", catch_stdout=False)
try:
    pysam_tview_text(
        bam, "chr1:12", output=work / "view.txt", samtools="samtools_not_installed", verbose=False
    )
except BackendUnavailableError as error:
    print(type(error).__name__, "samtools executable 'samtools_not_installed' was not found" in str(error))
print((work / "view.txt").exists())
```

- **Example output:**

```text
BackendUnavailableError True
False
```

Usage with the dependency installed:

```python
from dnakit.alignment import pysam_tview_text

result = pysam_tview_text("reads.bam", "chr1:12", output="view.txt", reference="ref.fa", width=60)
print(result["version"], result["line_count"])
print("\n".join(result["lines"][:5]))
```

## 15) `dnakit.sequence.faidx_check_default_seq()`

- **Purpose:** Check the single-character default base used for missing positions and masking.
- **Method:** Returns `None` for `None` and the value itself when it is exactly one character (not NUL). Anything else raises `ConfigurationError`. It follows the rule of `pyfaidx.cli.check_seq_length` for `--default-seq`; the function does not call pyfaidx.
- **API:** `faidx_check_default_seq(value: 'str | None', *, verbose: bool = True) -> 'str | None'`.
- **Requirements:** Works after plain pip install (no optional package is imported).
- **Input:** `value` is `None` or a string. A dash is accepted as a value.
- **Output:** `None` or the one-character string.
- **Example code:**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.sequence import faidx_check_default_seq

print(faidx_check_default_seq(None, verbose=False))
print(faidx_check_default_seq("-", verbose=False))
try:
    faidx_check_default_seq("NN", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, "default_seq must be a single character" in str(error))
```

- **Example output:**

```text
None
-
ConfigurationError True
```

## 16) `dnakit.sequence.faidx_fetch_sequences()`

- **Purpose:** Extract sequences for regions (or all sequences) into a new FASTA file.
- **Method:** Runs `faidx` with the flags built from the arguments (`-c` complement, `-r` reverse, `-y` auto strand, `-n` no names, `-f` long names, `-t` no coordinates, `-l` lazy, `-s` default base, `-d` delimiter, `--no-rebuild`, and `-b` for a BED file) and streams its standard output into `output`. It does not use `faidx -o`, because pyfaidx 0.8.2 can lose the tail of a larger output file written that way.
- **API:** `faidx_fetch_sequences(fasta: 'str | Path', output: 'str | Path', *, regions: 'Sequence[str]' = (), bed: 'str | Path | None' = None, complement: 'bool' = False, reverse: 'bool' = False, auto_strand: 'bool' = False, no_names: 'bool' = False, long_names: 'bool' = False, no_coords: 'bool' = False, lazy: 'bool' = False, default_seq: 'str | None' = None, delimiter: 'str | None' = None, no_rebuild: 'bool' = False, faidx_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package pyfaidx (extra `dnakit[bio-pyfaidx]`), which installs the `faidx` program. DNAKit looks for `faidx` in `faidx_path`, then on PATH, then in the folder of the running Python. faidx writes a `.fai` index next to the FASTA when there is none.
- **Input:** `fasta` is an existing file and `output` a new file in an existing folder. Give `regions` (`name:start-end`, 1-based and inclusive) or `bed`, not both; with neither, every sequence is written. `no_names` and `long_names` cannot both be true. `default_seq` is one character.
- **Output:** A dictionary with `status` (`completed` for exit code 0, otherwise `failed`), `backend`, `version`, `command`, `return_code`, `outputs` (the new file) and `log_tail` (faidx's standard error). A region that is not in the FASTA is reported in `log_tail` and gives no record.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import faidx_fetch_sequences

out = Path(tempfile.mkdtemp()) / "regions.fa"
result = faidx_fetch_sequences(
    "genes.fa", str(out), regions=["g1:13-20", "g2:9-16"], complement=True, reverse=True, verbose=False
)
print(result["status"], result["backend"], result["return_code"])
print(out.read_text(), end="")
```

- **Example output:**

```text
completed faidx 0
>g1:20-13 (complement)
CTCGGTCA
>g2:16-9 (complement)
TTTTGGGG
```

## 17) `dnakit.sequence.faidx_mask_sequences()`

- **Purpose:** Replace the bases of given regions with a mask character in a new copy of a FASTA.
- **Method:** Copies the FASTA to `output`, then runs `faidx -m -s <default_seq>` on the copy for the regions. The input is never changed. Lowercase masking (`-M`) is not offered because pyfaidx 0.8.2 fails on it.
- **API:** `faidx_mask_sequences(fasta: 'str | Path', output: 'str | Path', *, regions: 'Sequence[str]' = (), bed: 'str | Path | None' = None, default_seq: 'str' = 'N', delimiter: 'str | None' = None, no_rebuild: 'bool' = False, faidx_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package pyfaidx (extra `dnakit[bio-pyfaidx]`), which installs the `faidx` program. DNAKit looks for `faidx` in `faidx_path`, then on PATH, then in the folder of the running Python.
- **Input:** `fasta` is an existing file and `output` a new file. At least one region (`name:start-end`, 1-based and inclusive) or a `bed` file is required (`ConfigurationError` otherwise). `default_seq` is the one-character replacement (default `N`).
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (the masked copy) and `log_tail`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import faidx_mask_sequences

out = Path(tempfile.mkdtemp()) / "masked.fa"
result = faidx_mask_sequences("genes.fa", str(out), regions=["g1:1-4", "g3:5-7"], verbose=False)
print(result["status"])
print(out.read_text(), end="")
print(Path("genes.fa").read_text().splitlines()[1])
```

- **Example output:**

```text
completed
>g1
NNNNACGTACGTTGACCGAGGGCCC
>g2
TTTTGGGGCCCCAAAANNACGTACGT
>g3
GATTNNN
ACGTACGTACGTTGACCGAGGGCCC
```

## 18) `dnakit.sequence.faidx_parse_bed_entry()`

- **Purpose:** Parse one BED line into a chromosome name and 0-based start and end.
- **Method:** Calls `pyfaidx.bed_split()` in the Python process. A comment line (starting with `#`) gives `None` for all three fields. A malformed line raises `ConfigurationError`.
- **API:** `faidx_parse_bed_entry(line: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package pyfaidx (extra `dnakit[bio-pyfaidx]`); no external program is run.
- **Input:** `line` is a string (tab separated; the first three columns are used).
- **Output:** A dictionary with `name`, `start` and `end` (integers, BED coordinates, 0-based half-open), or all `None` for a comment.
- **Example code:**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.sequence import faidx_parse_bed_entry

print(faidx_parse_bed_entry("g1\t2\t10", verbose=False))
print(faidx_parse_bed_entry("# a comment", verbose=False))
try:
    faidx_parse_bed_entry("g1\tx\t10", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, "Malformed BED entry" in str(error))
```

- **Example output:**

```text
{'name': 'g1', 'start': 2, 'end': 10}
{'name': None, 'start': None, 'end': None}
ConfigurationError True
```

## 19) `dnakit.sequence.faidx_parse_ucsc_region()`

- **Purpose:** Parse a UCSC region string into a name and 0-based start and exclusive end.
- **Method:** Calls `pyfaidx.ucsc_split()`. The input is 1-based and inclusive (`g1:3-10`); the result is 0-based and end-exclusive (start 2, end 10). A bare name gives `None` for start and end, which means the whole sequence.
- **API:** `faidx_parse_ucsc_region(region: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package pyfaidx (extra `dnakit[bio-pyfaidx]`); no external program is run.
- **Input:** `region` is a non-empty string that does not start with `-`.
- **Output:** A dictionary with `name`, `start` and `end`.
- **Example code:**

```python
from dnakit.sequence import faidx_parse_ucsc_region

print(faidx_parse_ucsc_region("g1:3-10", verbose=False))
print(faidx_parse_ucsc_region("g1", verbose=False))
```

- **Example output:**

```text
{'name': 'g1', 'start': 2, 'end': 10}
{'name': 'g1', 'start': None, 'end': None}
```

## 20) `dnakit.sequence.faidx_split_regions()`

- **Purpose:** Turn region strings or a BED file into name, start and end records.
- **Method:** For region strings it calls `faidx_parse_ucsc_region()` on each; for a BED file it calls `faidx_parse_bed_entry()` on each line and drops the comment lines. The order of the input is kept.
- **API:** `faidx_split_regions(regions: 'Sequence[str]' = (), *, bed: 'str | Path | None' = None, verbose: bool = True) -> 'list[dict[str, Any]]'`.
- **Requirements:** Needs the optional package pyfaidx (extra `dnakit[bio-pyfaidx]`); no external program is run.
- **Input:** `regions` is a sequence of `name:start-end` strings; `bed` is an existing BED file. Giving both raises `ConfigurationError`. With neither, the result is an empty list.
- **Output:** A list of dictionaries with `name`, `start` and `end`.
- **Example code:**

```python
from dnakit.sequence import faidx_split_regions

print(faidx_split_regions(["g1:3-10", "g2"], verbose=False))
print(faidx_split_regions(bed="regions.bed", verbose=False))
print(faidx_split_regions(verbose=False))
```

- **Example output:**

```text
[{'name': 'g1', 'start': 2, 'end': 10}, {'name': 'g2', 'start': None, 'end': None}]
[{'name': 'g1', 'start': 2, 'end': 10}, {'name': 'g2', 'start': 0, 'end': 8}]
[]
```

## 21) `dnakit.sequence.faidx_transform_sequences()`

- **Purpose:** Describe the sequences or regions of a FASTA as a table instead of FASTA.
- **Method:** Runs `faidx -i <transform>` and streams its standard output into `output`. `bed` writes name, 0-based start and end; `chromsizes` writes name and length; `nucleotide` writes a header line and the counts of A, T, C, G, N and other characters; `transposed` writes name, start, end and the sequence.
- **API:** `faidx_transform_sequences(fasta: 'str | Path', output: 'str | Path', transform: 'str', *, regions: 'Sequence[str]' = (), bed: 'str | Path | None' = None, complement: 'bool' = False, reverse: 'bool' = False, delimiter: 'str | None' = None, no_rebuild: 'bool' = False, faidx_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package pyfaidx (extra `dnakit[bio-pyfaidx]`), which installs the `faidx` program. DNAKit looks for `faidx` in `faidx_path`, then on PATH, then in the folder of the running Python.
- **Input:** `transform` is one of `bed`, `chromsizes`, `nucleotide`, `transposed`. `fasta` is an existing file and `output` a new file. `regions` (or `bed`), `complement`, `reverse`, `delimiter` and `no_rebuild` are optional.
- **Output:** A dictionary with `status` (`completed` for exit code 0, otherwise `failed`), `backend`, `version`, `command`, `return_code`, `outputs` (the new file) and `log_tail` (faidx's standard error).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import faidx_transform_sequences

out = Path(tempfile.mkdtemp()) / "composition.tsv"
result = faidx_transform_sequences("genes.fa", str(out), "nucleotide", verbose=False)
print(result["status"])
for line in out.read_text().splitlines():
    print(line.rstrip("	"))
```

- **Example output:**

```text
completed
name	start	end	A	T	C	G	N	others
g1	1	25	5	4	8	8	0
g2	1	26	6	6	6	6	2
g3	1	7	3	2	1	1	0
```

## 22) `dnakit.sequence.faidx_wrap_sequence()`

- **Purpose:** Wrap a sequence into lines of a fixed width.
- **Method:** Calls `pyfaidx.wrap_sequence()` and returns its lines as a list. Every line, including the last one, ends with a newline.
- **API:** `faidx_wrap_sequence(sequence: 'str', width: 'int', *, verbose: bool = True) -> 'list[str]'`.
- **Requirements:** Needs the optional package pyfaidx (extra `dnakit[bio-pyfaidx]`); no external program is run.
- **Input:** `sequence` is a string; `width` is an integer of at least 1.
- **Output:** A list of strings, each at most `width` characters plus the newline.
- **Example code:**

```python
from dnakit.sequence import faidx_wrap_sequence

print(faidx_wrap_sequence("ACGTACGTAC", 4, verbose=False))
```

- **Example output:**

```text
['ACGT\n', 'ACGT\n', 'AC\n']
```

## 23) `dnakit.sequence.faidx_write_sequences()`

- **Purpose:** Write FASTA sequences selected by a name regular expression and a length range, to one file or one file per sequence.
- **Method:** Runs `faidx -g <regex>` (with `-v` to invert, `-a low,high` for the size range, `-x` for split files, and the same output flags as `faidx_fetch_sequences()`). With `output` the standard output is streamed into the new file; with `split_directory` faidx writes one file per selected sequence into the new directory, named after the sequence and its coordinates, and the function lists them.
- **API:** `faidx_write_sequences(fasta: 'str | Path', *, output: 'str | Path | None' = None, split_directory: 'str | Path | None' = None, regions: 'Sequence[str]' = (), bed: 'str | Path | None' = None, regex: 'str' = '.*', invert_match: 'bool' = False, size_range: 'Sequence[int] | None' = None, complement: 'bool' = False, reverse: 'bool' = False, auto_strand: 'bool' = False, no_names: 'bool' = False, long_names: 'bool' = False, no_coords: 'bool' = False, lazy: 'bool' = False, default_seq: 'str | None' = None, delimiter: 'str | None' = None, no_rebuild: 'bool' = False, faidx_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package pyfaidx (extra `dnakit[bio-pyfaidx]`), which installs the `faidx` program. DNAKit looks for `faidx` in `faidx_path`, then on PATH, then in the folder of the running Python.
- **Input:** Give exactly one of `output` (a new file) and `split_directory` (a new directory) or `ConfigurationError` is raised. `regex` selects sequence names and is not anchored, so use `^...$` for an exact name. `size_range=(low, high)` is inclusive and keeps sequences (or regions) whose length is in the range. `regions`/`bed` limit the regions; the other flags are as in `faidx_fetch_sequences()`.
- **Output:** A dictionary with `status` (`completed` for exit code 0, otherwise `failed`), `backend`, `version`, `command`, `return_code`, `outputs` (the new file) and `log_tail` (faidx's standard error). With `split_directory`, `outputs` lists the files written there.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import faidx_write_sequences

work = Path(tempfile.mkdtemp())
result = faidx_write_sequences(
    "genes.fa", output=str(work / "long.fa"), size_range=(10, 100), verbose=False
)
print(result["status"], [Path(path).name for path in result["outputs"]])
print((work / "long.fa").read_text(), end="")
split = faidx_write_sequences("genes.fa", split_directory=str(work / "split"), regex="^g3$", verbose=False)
print(split["status"], [Path(path).name for path in split["outputs"]])
```

- **Example output:**

```text
completed ['long.fa']
>g1
ACGTACGTACGTTGACCGAGGGCCC
>g2
TTTTGGGGCCCCAAAANNACGTACGT
completed ['g3.fa']
```

## 24) `dnakit.fasta_fastq.cutadapt_full_report()`

- **Purpose:** Trim adapters and low-quality ends of a FASTQ file and return cutadapt's full text report.
- **Method:** The same code as `cutadapt_trim()` with `--report full`: it runs `cutadapt -e <error_rate> -j <cores> [-a/-g/-b adapter] [-q] [-m] [-M] [--discard-untrimmed] -o <output> <input>` without a shell. The report has the summary, a read-fate breakdown, one section per adapter and the length histograms.
- **API:** `cutadapt_full_report(input_path: 'str', output_path: 'str', *, adapter_3p: 'str | None' = None, adapter_5p: 'str | None' = None, anywhere_adapter: 'str | None' = None, error_rate: 'float' = 0.1, quality_cutoff: 'str | None' = None, minimum_length: 'int | None' = None, maximum_length: 'int | None' = None, discard_untrimmed: 'bool' = False, cores: 'int' = 1, cutadapt_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package cutadapt (extra `dnakit[bio-cutadapt]`), which installs the `cutadapt` program. DNAKit looks for it in `cutadapt_path`, then on PATH, then in the folder of the running Python.
- **Input:** `input_path` is an existing FASTQ and `output_path` a new file. Adapters (`adapter_3p` for `-a`, `adapter_5p` for `-g`, `anywhere_adapter` for `-b`) are IUPAC DNA letters and are upper-cased. `error_rate` is between 0 and 1; `quality_cutoff` looks like `20` or `10,20`; `minimum_length` and `maximum_length` are non-negative (maximum not below minimum); `cores` is at least 1.
- **Output:** A dictionary with `status` (`completed` when the exit code is 0 and the output file exists, otherwise `failed`), `backend`, `version`, `command`, `return_code`, `outputs` (the files that exist), `report` (cutadapt's standard output) and `log_tail` (the tail of its standard error). The report contains the command line, which includes the file paths.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.fasta_fastq import cutadapt_full_report

out = Path(tempfile.mkdtemp()) / "trimmed.fq"
result = cutadapt_full_report("reads.fq", str(out), adapter_3p="AGATCGGAAGAGC", verbose=False)
print(result["status"], result["return_code"])
for line in result["report"].splitlines():
    if line.startswith(("Total reads", "Reads with adapters", "Reads written", "Sequence:")):
        print(line)
```

- **Example output:**

```text
completed 0
Total reads processed:                       4
Reads with adapters:                         2 (50.0%)
Reads written (passing filters):             4 (100.0%)
Sequence: AGATCGGAAGAGC; Type: regular 3'; Length: 13; Trimmed: 2 times
```

## 25) `dnakit.fasta_fastq.cutadapt_json_report()`

- **Purpose:** Trim adapters and low-quality ends of a FASTQ file and also write cutadapt's JSON report.
- **Method:** The same code as `cutadapt_trim()` with `--report full --json <json_path>`. After the run the JSON file is parsed and returned under `json`.
- **API:** `cutadapt_json_report(input_path: 'str', output_path: 'str', json_path: 'str', *, adapter_3p: 'str | None' = None, adapter_5p: 'str | None' = None, anywhere_adapter: 'str | None' = None, error_rate: 'float' = 0.1, quality_cutoff: 'str | None' = None, minimum_length: 'int | None' = None, maximum_length: 'int | None' = None, discard_untrimmed: 'bool' = False, cores: 'int' = 1, cutadapt_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package cutadapt (extra `dnakit[bio-cutadapt]`), which installs the `cutadapt` program. DNAKit looks for it in `cutadapt_path`, then on PATH, then in the folder of the running Python.
- **Input:** As `cutadapt_full_report()`, plus `json_path`, a new file path (it must not exist).
- **Output:** A dictionary with `status` (`completed` when the exit code is 0 and the output file exists, otherwise `failed`), `backend`, `version`, `command`, `return_code`, `outputs` (the files that exist), `report` (cutadapt's standard output) and `log_tail` (the tail of its standard error). Also `json`, the parsed JSON report (keys such as `read_counts`, `basepair_counts`, `adapters_read1`, `command_line_arguments`); `outputs` then holds the trimmed FASTQ and the JSON file.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.fasta_fastq import cutadapt_json_report

work = Path(tempfile.mkdtemp())
result = cutadapt_json_report(
    "reads.fq", str(work / "trimmed.fq"), str(work / "report.json"),
    adapter_3p="AGATCGGAAGAGC", verbose=False,
)
print(result["status"], len(result["outputs"]))
counts = result["json"]["read_counts"]
print(counts["input"], counts["output"], counts["read1_with_adapter"])
```

- **Example output:**

```text
completed 2
4 4 2
```

## 26) `dnakit.fasta_fastq.cutadapt_minimal_report()`

- **Purpose:** Trim adapters and low-quality ends of a FASTQ file and return cutadapt's two-line tab-separated minimal report.
- **Method:** The same code as `cutadapt_trim()` with `--report minimal`.
- **API:** `cutadapt_minimal_report(input_path: 'str', output_path: 'str', *, adapter_3p: 'str | None' = None, adapter_5p: 'str | None' = None, anywhere_adapter: 'str | None' = None, error_rate: 'float' = 0.1, quality_cutoff: 'str | None' = None, minimum_length: 'int | None' = None, maximum_length: 'int | None' = None, discard_untrimmed: 'bool' = False, cores: 'int' = 1, cutadapt_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package cutadapt (extra `dnakit[bio-cutadapt]`), which installs the `cutadapt` program. DNAKit looks for it in `cutadapt_path`, then on PATH, then in the folder of the running Python.
- **Input:** As `cutadapt_full_report()`.
- **Output:** A dictionary with `status` (`completed` when the exit code is 0 and the output file exists, otherwise `failed`), `backend`, `version`, `command`, `return_code`, `outputs` (the files that exist), `report` (cutadapt's standard output) and `log_tail` (the tail of its standard error). In `report`, the first line holds the column names (`status`, `in_reads`, `in_bp`, `too_short`, `too_long`, `too_many_n`, `out_reads`, `w/adapters`, `qualtrim_bp`, `out_bp`) and the second line the values.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.fasta_fastq import cutadapt_minimal_report

out = Path(tempfile.mkdtemp()) / "trimmed.fq"
result = cutadapt_minimal_report("reads.fq", str(out), adapter_3p="AGATCGGAAGAGC", verbose=False)
header, values = result["report"].strip().splitlines()
for name, value in zip(header.split("\t"), values.split("\t")):
    print(name, value)
```

- **Example output:**

```text
status OK
in_reads 4
in_bp 85
too_short 0
too_long 0
too_many_n 0
out_reads 4
w/adapters 2
qualtrim_bp 0
out_bp 54
```

## 27) `dnakit.fasta_fastq.cutadapt_trim()`

- **Purpose:** Trim 3' and 5' adapters and low-quality ends of one FASTQ file.
- **Method:** Runs `cutadapt --report <full|minimal> -e <error_rate> -j <cores>` with `-a`, `-g`, `-b`, `-q`, `-m`, `-M` and `--discard-untrimmed` as requested, writing the reads to `output_path`. The other three report functions are this function with a fixed report type.
- **API:** `cutadapt_trim(input_path: 'str', output_path: 'str', *, adapter_3p: 'str | None' = None, adapter_5p: 'str | None' = None, anywhere_adapter: 'str | None' = None, error_rate: 'float' = 0.1, quality_cutoff: 'str | None' = None, minimum_length: 'int | None' = None, maximum_length: 'int | None' = None, discard_untrimmed: 'bool' = False, cores: 'int' = 1, report: 'str' = 'full', cutadapt_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package cutadapt (extra `dnakit[bio-cutadapt]`), which installs the `cutadapt` program. DNAKit looks for it in `cutadapt_path`, then on PATH, then in the folder of the running Python.
- **Input:** As `cutadapt_full_report()`, plus `report`, `"full"` (default) or `"minimal"`. Invalid values raise `ConfigurationError` before cutadapt starts. `discard_untrimmed=True` drops reads without an adapter match.
- **Output:** A dictionary with `status` (`completed` when the exit code is 0 and the output file exists, otherwise `failed`), `backend`, `version`, `command`, `return_code`, `outputs` (the files that exist), `report` (cutadapt's standard output) and `log_tail` (the tail of its standard error).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.fasta_fastq import cutadapt_trim

out = Path(tempfile.mkdtemp()) / "trimmed.fq"
result = cutadapt_trim(
    "reads.fq", str(out), adapter_3p="AGATCGGAAGAGC", quality_cutoff="10",
    minimum_length=10, report="minimal", verbose=False,
)
print(result["status"], result["return_code"])
lines = out.read_text().splitlines()
for name, sequence in zip(lines[0::4], lines[1::4]):
    print(name, sequence)
```

- **Example output:**

```text
completed 0
@r2 TTGACCGATC
@r3 GGGCCCAAATTTGGGAAA
@r4 NNACGTACGTTGAACG
```

## 28) `dnakit.fasta_fastq.cutadapt_adapter_action()`

- **Purpose:** Find adapters in every read and apply one action to each match: trim, retain, mask with N, lowercase, crop or none.
- **Method:** Runs `cutadapt --report minimal -e <error_rate> --action <action> -n <times> [--rc]` with the adapters given by `-a`, `-g` and `-b`. With `--rc` cutadapt also searches the reverse complement of each read and writes such a read in its reverse-complemented form, adding ` rc` to its name.
- **API:** `cutadapt_adapter_action(input_path: 'str', output_path: 'str', *, adapter_3p: 'str | None' = None, adapter_5p: 'str | None' = None, anywhere_adapter: 'str | None' = None, action: 'str' = 'trim', times: 'int' = 1, reverse_complement: 'bool' = False, error_rate: 'float' = 0.1, cutadapt_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package cutadapt (extra `dnakit[bio-cutadapt]`), which installs the `cutadapt` program. DNAKit looks for it in `cutadapt_path`, then on PATH, then in the folder of the running Python.
- **Input:** `input_path` is an existing FASTQ and `output_path` a new file. At least one adapter is required. `action` is `trim` (default), `retain`, `mask`, `lowercase`, `crop` or `none`. `times` (at least 1) is how often adapters are searched for in a read; `reverse_complement` sets `--rc`; `error_rate` is between 0 and 1.
- **Output:** A dictionary with `status` (`completed` when the exit code is 0 and the output file exists, otherwise `failed`), `backend`, `version`, `command`, `return_code`, `outputs` (the files that exist), `report` (cutadapt's standard output) and `log_tail` (the tail of its standard error). Here `report` is the one-line minimal report.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.fasta_fastq import cutadapt_adapter_action

work = Path(tempfile.mkdtemp())
result = cutadapt_adapter_action(
    "reads.fq", str(work / "masked.fq"), adapter_3p="AGATCGGAAGAGC", action="mask", verbose=False
)
print(result["status"])
lines = (work / "masked.fq").read_text().splitlines()
print(lines[1::4])
rc = cutadapt_adapter_action(
    "reads.fq", str(work / "rc.fq"), adapter_5p="TTTCCCAAATTTGGG",
    reverse_complement=True, verbose=False,
)
lines = (work / "rc.fq").read_text().splitlines()
print(lines[8], lines[9])
```

- **Example output:**

```text
completed
['ACGTTGCANNNNNNNNNNNNNNNNNN', 'TTGACCGATCNNNNNNNNNNNNN', 'GGGCCCAAATTTGGGAAA', 'NNACGTACGTTGAACGNN']
@r3 rc CCC
```

## 29) `dnakit.fasta_fastq.cutadapt_shorten()`

- **Purpose:** Shorten every read to a fixed length.
- **Method:** Runs `cutadapt --report minimal -l <length>`. A positive length keeps the first bases of each read, a negative one keeps the last bases. Reads that are already shorter are not changed.
- **API:** `cutadapt_shorten(input_path: 'str', output_path: 'str', *, length: 'int', cutadapt_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package cutadapt (extra `dnakit[bio-cutadapt]`), which installs the `cutadapt` program. DNAKit looks for it in `cutadapt_path`, then on PATH, then in the folder of the running Python.
- **Input:** `input_path` is an existing FASTQ and `output_path` a new file. `length` is a non-zero integer (`ConfigurationError` for 0 or a non-integer).
- **Output:** A dictionary with `status` (`completed` when the exit code is 0 and the output file exists, otherwise `failed`), `backend`, `version`, `command`, `return_code`, `outputs` (the files that exist), `report` (the one-line minimal report) and `log_tail` (the tail of its standard error).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.fasta_fastq import cutadapt_shorten

out = Path(tempfile.mkdtemp()) / "short.fq"
result = cutadapt_shorten("reads.fq", str(out), length=-6, verbose=False)
print(result["status"])
print(out.read_text().splitlines()[1::4])
```

- **Example output:**

```text
completed
['CACACG', 'AAGAGC', 'GGGAAA', 'AACGNN']
```

## 30) `dnakit.fasta_fastq.cutadapt_trim_n_ends()`

- **Purpose:** Remove N bases from both ends of every read.
- **Method:** Runs `cutadapt --report minimal --trim-n`.
- **API:** `cutadapt_trim_n_ends(input_path: 'str', output_path: 'str', *, cutadapt_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package cutadapt (extra `dnakit[bio-cutadapt]`), which installs the `cutadapt` program. DNAKit looks for it in `cutadapt_path`, then on PATH, then in the folder of the running Python.
- **Input:** `input_path` is an existing FASTQ and `output_path` a new file.
- **Output:** A dictionary with `status` (`completed` when the exit code is 0 and the output file exists, otherwise `failed`), `backend`, `version`, `command`, `return_code`, `outputs` (the files that exist), `report` (the one-line minimal report) and `log_tail` (the tail of its standard error).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.fasta_fastq import cutadapt_trim_n_ends

out = Path(tempfile.mkdtemp()) / "no_n.fq"
result = cutadapt_trim_n_ends("reads.fq", str(out), verbose=False)
print(result["status"])
print(out.read_text().splitlines()[1::4][-1])
```

- **Example output:**

```text
completed
ACGTACGTTGAACG
```

## 31) `dnakit.fasta_fastq.cutadapt_unconditional_cut()`

- **Purpose:** Remove a fixed number of bases from the 5' and/or 3' end of every read.
- **Method:** Runs `cutadapt --report minimal -u <cut_5p>` for the 5' end and `-u -<cut_3p>` for the 3' end, whichever are non-zero.
- **API:** `cutadapt_unconditional_cut(input_path: 'str', output_path: 'str', *, cut_5p: 'int' = 0, cut_3p: 'int' = 0, cutadapt_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package cutadapt (extra `dnakit[bio-cutadapt]`), which installs the `cutadapt` program. DNAKit looks for it in `cutadapt_path`, then on PATH, then in the folder of the running Python.
- **Input:** `input_path` is an existing FASTQ and `output_path` a new file. `cut_5p` and `cut_3p` are non-negative integers; at least one must be positive (`ConfigurationError` otherwise).
- **Output:** A dictionary with `status` (`completed` when the exit code is 0 and the output file exists, otherwise `failed`), `backend`, `version`, `command`, `return_code`, `outputs` (the files that exist), `report` (the one-line minimal report) and `log_tail` (the tail of its standard error).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.fasta_fastq import cutadapt_unconditional_cut

out = Path(tempfile.mkdtemp()) / "cut.fq"
result = cutadapt_unconditional_cut("reads.fq", str(out), cut_5p=2, cut_3p=3, verbose=False)
print(result["status"])
print(out.read_text().splitlines()[1::4])
```

- **Example output:**

```text
completed
['GTTGCAAGATCGGAAGAGCAC', 'GACCGATCAGATCGGAAG', 'GCCCAAATTTGGG', 'ACGTACGTTGAAC']
```

## 32) `dnakit.sequencing_qc.umi_count_per_gene()`

- **Purpose:** Count UMI groups per gene in a BAM file and write a gene count table.
- **Method:** Runs `umi_tools count` with the UMI options (`--extract-umi-method`, `--umi-separator` or `--umi-tag`, `--method`, `--edit-distance-threshold`, `--paired`), `--gene-tag <tag>` or `--per-contig`, and `--mapping-quality`. Reads are grouped by UMI and position as in `umi_group_bam()` and the groups of each gene are counted.
- **API:** `umi_count_per_gene(bam: 'str', output: 'str', *, gene_tag: 'str | None' = None, per_contig: 'bool' = False, extract_umi_method: 'str' = 'read_id', umi_tag: 'str | None' = None, umi_separator: 'str' = '_', method: 'str' = 'directional', edit_distance_threshold: 'int' = 1, mapping_quality: 'int' = 0, paired: 'bool' = False, umi_tools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package umi_tools (extra `dnakit[bio-umi-tools]`), which installs the `umi_tools` program. DNAKit looks for it in `umi_tools_path`, then on PATH, then in the folder of the running Python.
- **Input:** `bam` is an existing coordinate-sorted BAM with an index and `output` a new file. Give either `gene_tag` (a BAM tag naming each read's gene) or `per_contig=True`, not both and not neither (`ConfigurationError`). `extract_umi_method` is `read_id` (UMI after `umi_separator` in the read name), `tag` (`umi_tag`, default `RX`) or `umis`. `method` is `adjacency`, `directional`, `percentile`, `unique` or `cluster`. `edit_distance_threshold` and `mapping_quality` are non-negative.
- **Output:** A dictionary with `status` (`completed` when the exit code is 0 and every expected output file exists, otherwise `failed`), `backend`, `version`, `command`, `return_code`, `outputs` (the files that exist) and `log_tail` (the tail of umi_tools' log). The results are in the output files. The output TSV has the columns `gene` and `count`.
- **Example code:**

```python
import tempfile
from pathlib import Path

import pysam
from dnakit.sequencing_qc import umi_count_per_gene

work = Path(tempfile.mkdtemp())
bam = str(work / "umi.bam")
pysam.view("-b", "-o", bam, "umi.sam", catch_stdout=False)
pysam.index(bam)
result = umi_count_per_gene(bam, str(work / "counts.tsv"), gene_tag="XT", verbose=False)
print(result["status"])
print((work / "counts.tsv").read_text(), end="")
```

- **Example output:**

```text
completed
gene	count
geneA	2
geneB	1
```

## 33) `dnakit.sequencing_qc.umi_count_table()`

- **Purpose:** Count UMI groups per gene from a two-column table of read ids and genes.
- **Method:** Runs `umi_tools count_tab` with `--barcode-separator`, `--method`, `--edit-distance-threshold` and, for `per_cell`, `--per-cell`. The UMIs of each gene are clustered and the groups are counted.
- **API:** `umi_count_table(table: 'str', output: 'str', *, per_cell: 'bool' = False, barcode_separator: 'str' = '_', method: 'str' = 'directional', edit_distance_threshold: 'int' = 1, umi_tools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package umi_tools (extra `dnakit[bio-umi-tools]`), which installs the `umi_tools` program. DNAKit looks for it in `umi_tools_path`, then on PATH, then in the folder of the running Python.
- **Input:** `table` is an existing tab-separated file with the read id in column 1 and the gene in column 2; the UMI is the last field of the read id after `barcode_separator` (default `_`). With `per_cell=True` the read id must end in `cell_barcode<sep>UMI`: the code takes the field before the last one as the cell (the program's help text says the reverse). `output` is a new file.
- **Output:** A dictionary with `status` (`completed` when the exit code is 0 and every expected output file exists, otherwise `failed`), `backend`, `version`, `command`, `return_code`, `outputs` (the files that exist) and `log_tail` (the tail of umi_tools' log). The results are in the output files. The output TSV has the columns `gene` and `count`, or `cell`, `gene` and `count` with `per_cell`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc import umi_count_table

work = Path(tempfile.mkdtemp())
table = work / "read_gene.tsv"
table.write_text("r1_CELL1_AAAA\tgeneA\nr2_CELL1_AAAA\tgeneA\nr3_CELL2_AAAA\tgeneA\nr4_CELL2_CCCC\tgeneB\n")
result = umi_count_table(str(table), str(work / "counts.tsv"), per_cell=True, verbose=False)
print(result["status"])
print((work / "counts.tsv").read_text(), end="")
```

- **Example output:**

```text
completed
cell	gene	count
CELL1	geneA	1
CELL2	geneA	1
CELL2	geneB	1
```

## 34) `dnakit.sequencing_qc.umi_dedup_bam()`

- **Purpose:** Remove PCR duplicates from a coordinate-sorted BAM file, keeping one read per UMI group and position.
- **Method:** Runs `umi_tools dedup` with the UMI options (as in `umi_group_bam()`), `--ignore-umi`, `--per-contig` or `--per-gene --gene-tag <tag>` as requested, and writes the reads to `output_bam`.
- **API:** `umi_dedup_bam(bam: 'str', output_bam: 'str', *, extract_umi_method: 'str' = 'read_id', umi_tag: 'str | None' = None, umi_separator: 'str' = '_', method: 'str' = 'directional', edit_distance_threshold: 'int' = 1, paired: 'bool' = False, ignore_umi: 'bool' = False, per_contig: 'bool' = False, gene_tag: 'str | None' = None, umi_tools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package umi_tools (extra `dnakit[bio-umi-tools]`), which installs the `umi_tools` program. DNAKit looks for it in `umi_tools_path`, then on PATH, then in the folder of the running Python.
- **Input:** `bam` is an existing coordinate-sorted BAM with an index and `output_bam` a new file. `extract_umi_method`, `umi_tag`, `umi_separator`, `method`, `edit_distance_threshold` and `paired` are as in `umi_count_per_gene()`. `ignore_umi=True` deduplicates on position only; `per_contig=True` or `gene_tag` deduplicates within contigs or genes.
- **Output:** A dictionary with `status` (`completed` when the exit code is 0 and every expected output file exists, otherwise `failed`), `backend`, `version`, `command`, `return_code`, `outputs` (the files that exist) and `log_tail` (the tail of umi_tools' log). The results are in the output files. The output is a BAM with one read for each UMI group at each position.
- **Example code:**

```python
import tempfile
from pathlib import Path

import pysam
from dnakit.sequencing_qc import umi_dedup_bam

work = Path(tempfile.mkdtemp())
bam = str(work / "umi.bam")
pysam.view("-b", "-o", bam, "umi.sam", catch_stdout=False)
pysam.index(bam)
result = umi_dedup_bam(bam, str(work / "dedup.bam"), verbose=False)
print(result["status"])
print(sorted(read.query_name for read in pysam.AlignmentFile(str(work / "dedup.bam"))))
```

- **Example output:**

```text
completed
['a2_AAAA', 'a4_CCCC', 'a6_GGGG']
```

## 35) `dnakit.sequencing_qc.umi_extract_to_names()`

- **Purpose:** Move the UMI and cell barcode from the start of each read into the read name.
- **Method:** Runs `umi_tools extract --extract-method=<string|regex> --bc-pattern=<pattern> --umi-separator=<sep>` (plus `--3prime`, and for paired reads `--read2-in`, `--read2-out`, `--bc-pattern2`) with the FASTQ on standard input and the new FASTQ on standard output. The barcode bases are removed from the read and appended to its name: the cell barcode and the UMI, separated by `_`.
- **API:** `umi_extract_to_names(fastq: 'str', output: 'str', *, bc_pattern: 'str', extract_method: 'str' = 'string', umi_separator: 'str' = '_', three_prime: 'bool' = False, read2_in: 'str | None' = None, read2_out: 'str | None' = None, bc_pattern2: 'str | None' = None, umi_tools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package umi_tools (extra `dnakit[bio-umi-tools]`), which installs the `umi_tools` program. DNAKit looks for it in `umi_tools_path`, then on PATH, then in the folder of the running Python.
- **Input:** `fastq` is an existing file and `output` a new file. With `extract_method="string"`, `bc_pattern` uses `C` (cell barcode), `N` (UMI) and `X` (discard) and must contain an `N`; with `"regex"` it is a regular expression with named groups. For paired reads give `read2_in` and `read2_out` (and `bc_pattern2` if read 2 has a barcode); `read2_out` or `bc_pattern2` without `read2_in` raises `ConfigurationError`. `three_prime=True` means the barcode is at the 3' end.
- **Output:** A dictionary with `status` (`completed` when the exit code is 0 and every expected output file exists, otherwise `failed`), `backend`, `version`, `command`, `return_code`, `outputs` (the files that exist) and `log_tail` (the tail of umi_tools' log). The results are in the output files. The output FASTQ has the barcodes in the read names.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc import umi_extract_to_names

out = Path(tempfile.mkdtemp()) / "extracted.fq"
result = umi_extract_to_names("umi_reads.fq", str(out), bc_pattern="CCCCCNNN", verbose=False)
print(result["status"], result["return_code"])
lines = out.read_text().splitlines()
for name, sequence in zip(lines[0::4], lines[1::4]):
    print(name, sequence)
```

- **Example output:**

```text
completed 0
@r1_ACGTA_AGG GCCCAAATTTGGG
@r2_TTGCC_AGG GCCCAAATTTGGG
@r3_ACGTA_AGG GCCCAAATTTGGG
@r4_ACGTT_TGG GCCCAAATTTGGG
```

## 36) `dnakit.sequencing_qc.umi_group_bam()`

- **Purpose:** Group the reads of a BAM file by UMI and mapping position.
- **Method:** Runs `umi_tools group --output-bam` with the UMI options, `--umi-group-tag=<tag>` and, if requested, `--group-out=<group_tsv>`. Each read of the output BAM gets the group in the tag `umi_group_tag` (`BX` by default).
- **API:** `umi_group_bam(bam: 'str', output_bam: 'str', *, group_tsv: 'str | None' = None, extract_umi_method: 'str' = 'read_id', umi_tag: 'str | None' = None, umi_separator: 'str' = '_', umi_group_tag: 'str' = 'BX', method: 'str' = 'directional', edit_distance_threshold: 'int' = 1, paired: 'bool' = False, umi_tools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package umi_tools (extra `dnakit[bio-umi-tools]`), which installs the `umi_tools` program. DNAKit looks for it in `umi_tools_path`, then on PATH, then in the folder of the running Python.
- **Input:** `bam` is an existing coordinate-sorted BAM with an index and `output_bam` a new file; `group_tsv` is an optional new file. `extract_umi_method`, `umi_tag`, `umi_separator`, `method`, `edit_distance_threshold` and `paired` are as in `umi_count_per_gene()`.
- **Output:** A dictionary with `status` (`completed` when the exit code is 0 and every expected output file exists, otherwise `failed`), `backend`, `version`, `command`, `return_code`, `outputs` (the files that exist) and `log_tail` (the tail of umi_tools' log). The results are in the output files. `group_tsv` has the columns `read_id`, `contig`, `position`, `gene`, `umi`, `umi_count`, `final_umi`, `final_umi_count` and `unique_id`.
- **Example code:**

```python
import tempfile
from pathlib import Path

import pysam
from dnakit.sequencing_qc import umi_group_bam

work = Path(tempfile.mkdtemp())
bam = str(work / "umi.bam")
pysam.view("-b", "-o", bam, "umi.sam", catch_stdout=False)
pysam.index(bam)
result = umi_group_bam(bam, str(work / "grouped.bam"), group_tsv=str(work / "groups.tsv"), verbose=False)
print(result["status"], len(result["outputs"]))
for read in pysam.AlignmentFile(str(work / "grouped.bam")):
    print(read.query_name, read.get_tag("BX"))
rows = (work / "groups.tsv").read_text().splitlines()
print(rows[0].split("\t")[-4:])
```

- **Example output:**

```text
completed 2
a1_AAAA AAAA
a2_AAAA AAAA
a3_AAAT AAAA
a4_CCCC CCCC
a5_GGGG GGGG
a6_GGGG GGGG
['umi_count', 'final_umi', 'final_umi_count', 'unique_id']
```

## 37) `dnakit.sequencing_qc.umi_prepare_bam_for_rsem()`

- **Purpose:** Prepare a paired transcriptome-aligned BAM for RSEM.
- **Method:** Runs `umi_tools prepare-for-rsem [--tags=<tags>] [--sam]`. The program pairs the alignments of each read name, copies the comma-separated tags from read 1 to read 2 and writes the BAM in the form RSEM expects.
- **API:** `umi_prepare_bam_for_rsem(bam: 'str', output_bam: 'str', *, tags: 'str | None' = None, sam: 'bool' = False, umi_tools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package umi_tools (extra `dnakit[bio-umi-tools]`), which installs the `umi_tools` program. DNAKit looks for it in `umi_tools_path`, then on PATH, then in the folder of the running Python.
- **Input:** `bam` is an existing paired-end file (the example uses a name-sorted one) and `output_bam` a new file. `tags` is a comma-separated list of tag names (for example `RX`). `sam=True` reads and writes SAM text. Single-end input produced a BAM with no reads in the authors' check, so use paired data.
- **Output:** A dictionary with `status` (`completed` when the exit code is 0 and every expected output file exists, otherwise `failed`), `backend`, `version`, `command`, `return_code`, `outputs` (the files that exist) and `log_tail` (the tail of umi_tools' log). The results are in the output files.
- **Example code:**

```python
import tempfile
from pathlib import Path

import pysam
from dnakit.sequencing_qc import umi_prepare_bam_for_rsem

work = Path(tempfile.mkdtemp())
bam = str(work / "pair.bam")
pysam.view("-b", "-o", bam, "pair.sam", catch_stdout=False)
result = umi_prepare_bam_for_rsem(bam, str(work / "rsem.bam"), tags="RX", verbose=False)
print(result["status"])
for read in pysam.AlignmentFile(str(work / "rsem.bam")):
    print(read.query_name, read.flag, read.get_tag("RX"))
```

- **Example output:**

```text
completed
p1_AAAA 99 AAAA
p1_AAAA 147 AAAA
p2_CCCC 99 CCCC
p2_CCCC 147 CCCC
```

## 38) `dnakit.sequencing_qc.umi_whitelist_cells()`

- **Purpose:** Estimate which cell barcodes are present in a FASTQ file and write them as a whitelist.
- **Method:** Runs `umi_tools whitelist --extract-method=<...> --bc-pattern=<...> --method=<reads|umis> --knee-method=<distance|density>` (plus `--expect-cells` and `--set-cell-number` when given). The plot options of umi_tools are not used.
- **API:** `umi_whitelist_cells(fastq: 'str', output: 'str', *, bc_pattern: 'str', extract_method: 'str' = 'string', method: 'str' = 'reads', knee_method: 'str' = 'distance', expect_cells: 'int | None' = None, set_cell_number: 'int | None' = None, umi_tools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package umi_tools (extra `dnakit[bio-umi-tools]`), which installs the `umi_tools` program. DNAKit looks for it in `umi_tools_path`, then on PATH, then in the folder of the running Python.
- **Input:** `fastq` is an existing file and `output` a new file. `bc_pattern` and `extract_method` are as in `umi_extract_to_names()`, but with the `string` method the pattern must contain a `C`. `method` counts reads (`reads`) or unique UMIs (`umis`) per barcode. `knee_method="density"` needs `expect_cells`. `expect_cells` and `set_cell_number` are integers of at least 1.
- **Output:** A dictionary with `status` (`completed` when the exit code is 0 and every expected output file exists, otherwise `failed`), `backend`, `version`, `command`, `return_code`, `outputs` (the files that exist) and `log_tail` (the tail of umi_tools' log). The results are in the output files. The whitelist file has one barcode per line with four tab-separated columns: the barcode, the error-corrected barcodes that map to it, its count and their counts.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc import umi_whitelist_cells

out = Path(tempfile.mkdtemp()) / "whitelist.txt"
result = umi_whitelist_cells("umi_reads.fq", str(out), bc_pattern="CCCCCNNN", verbose=False)
print(result["status"])
for line in out.read_text().splitlines():
    print(line.split("\t"))
```

- **Example output:**

```text
completed
['ACGTA', '', '2', '']
['ACGTT', '', '1', '']
['TTGCC', '', '1', '']
```
