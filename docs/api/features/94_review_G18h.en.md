# Review G18h: pysam and hgvs variant helpers

This page covers offline wrappers for the pysam subcommands (samtools and bcftools) and the hgvs functions that were not yet covered. Every function reads local inputs only and writes to new paths. The remaining bed-reader, pandas-plink and bgen-reader rows are internal helpers or test entry points with no public genotype operation, so they are not wrapped here. The hgvs reference functions need seqrepo and PostgreSQL and are marked needs-data, so they are not on this page.

Modules: `dnakit.variants` and `dnakit.variants`.

## Example files

These files are written to the working directory before the examples run; each example writes its results to a new temporary directory.

**Example file:**
`reads.sam`

```text
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:120
@SQ	SN:chr2	LN:80
@RG	ID:rg1	SM:sampleA
r1	0	chr1	1	60	10M	*	0	0	ACGTACGTTA	IIIIIIIIII	RG:Z:rg1
r2	0	chr1	5	60	10M	*	0	0	ACGTTTAGCC	IIIIIIIIII	RG:Z:rg1
r3	16	chr1	20	30	8M2I	*	0	0	GATCGATCGG	IIIIIIIIII	RG:Z:rg1
r4	0	chr2	3	60	10M	*	0	0	GGCCTTAAGG	IIIIIIIIII	RG:Z:rg1
r5	4	*	0	0	*	*	0	0	ACGTACGTAC	IIIIIIIIII	RG:Z:rg1
```

**Example file:**
`ref.fa`

```text
>chr1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
>chr2
GGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCC
```

**Example file:**
`primers.bed`

```text
chr1	0	10	amp1_LEFT	60	+
chr1	100	110	amp1_RIGHT	60	-
```

**Example file:**
`intervals.bed`

```text
chr1	0	10
chr1	0	120
```

**Example file:**
`genes.gff3`

```text
##gff-version 3
chr1	test	gene	1	120	.	+	.	ID=gene1;Name=G1
chr1	test	mRNA	1	120	.	+	.	ID=tx1;Parent=gene1;Name=T1
chr1	test	CDS	1	120	.	+	0	ID=cds1;Parent=tx1
```

**Example file:**
`names.txt`

```text
A1
A2
```

**Example file:**
`header.sam`

```text
@HD	VN:1.6	SO:unsorted
@SQ	SN:chr1	LN:120
@SQ	SN:chr2	LN:80
@RG	ID:new	SM:sampleB
```

**Example file:**
`reads.fq`

```text
@r1
ACGTACGTAC
+
IIIIIIIIII
@r2
GGGGCCCCAA
+
IIIIIIIIII
```

**Example file:**
`calls.vcf`

```text
##fileformat=VCFv4.2
##contig=<ID=chr1,length=120>
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##FILTER=<ID=LowQual,Description="Low quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
chr1	5	.	A	G	50	PASS	AF=0.5	GT:DP	0/1:20	1/1:18
chr1	30	.	C	T	60	PASS	AF=0.5	GT:DP	0/0:20	0/1:18
chr1	60	.	T	C	40	LowQual	AF=0.5	GT:DP	0/1:20	0/0:18
```

**Example file:**
`calls2.vcf`

```text
##fileformat=VCFv4.2
##contig=<ID=chr1,length=120>
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##FILTER=<ID=LowQual,Description="Low quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S3	S4
chr1	5	.	A	G	50	PASS	AF=0.5	GT:DP	0/1:20	1/1:18
chr1	30	.	C	T	60	PASS	AF=0.5	GT:DP	0/0:20	0/1:18
chr1	60	.	T	C	40	LowQual	AF=0.5	GT:DP	0/1:20	0/0:18
```

**Example file:**
`phased.vcf`

```text
##fileformat=VCFv4.2
##contig=<ID=chr1,length=120>
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##FILTER=<ID=LowQual,Description="Low quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
chr1	5	.	A	G	50	PASS	AF=0.5	GT:DP	0|1:20	1|1:18
chr1	30	.	C	T	60	PASS	AF=0.5	GT:DP	0|0:20	0|1:18
chr1	60	.	T	C	40	LowQual	AF=0.5	GT:DP	0|1:20	0|0:18
```

## samtools alignment-file functions

## 1) Add a read group to every read

- **Purpose:** Write an @RG header line and tag every read with its ID.
- **Method:** Calls `samtools addreplacerg`; `mode` is overwrite_all or orphan_only.
- **API:** `dnakit.variants.add_bam_read_group(source, output_path, read_group, *, mode='overwrite_all', force=False, verbose=True)`
- **Input:** Alignment file, @RG line and mode; the output must be a new file.
- **Output:** Returns read_group_id, mode, records and the output path.

**Example code:**

```python
from dnakit.variants import add_bam_read_group

import tempfile
from pathlib import Path
from dnakit.variants import select_alignment_records

work = Path(tempfile.mkdtemp())
bam = select_alignment_records("reads.sam", work / "reads.bam")["output"]
result = add_bam_read_group(bam, work / "rg.bam", "@RG\tID:rgX\tSM:sampleB")
print(result["read_group_id"], result["records"])
```

**Example output:**

```text
rgX 5
```

## 2) Soft-clip amplicon primers

- **Purpose:** Soft-clip the bases inside primer intervals; read count is kept.
- **Method:** Calls `samtools ampliconclip -b BED`; BED intervals are 0-based.
- **API:** `dnakit.variants.clip_amplicon_primers(source, primers_bed, output_path, *, force=False, verbose=True)`
- **Input:** Alignment file and primer BED; the output is a new file.
- **Output:** Returns records and the output path.

**Example code:**

```python
import shutil
import tempfile
from pathlib import Path
from dnakit.variants import clip_amplicon_primers, select_alignment_records

work = Path(tempfile.mkdtemp())
shutil.copy("primers.bed", work / "primers.bed")
bam = select_alignment_records("reads.sam", work / "reads.bam")["output"]
result = clip_amplicon_primers(bam, work / "primers.bed", work / "clipped.bam")
print(result["records"])
```

**Example output:**

```text
5
```

## 3) Amplicon layout and statistics

- **Purpose:** Read amplicon positions and per-reference read statistics.
- **Method:** Calls `samtools ampliconstats` and parses the AMPLICON and FSS rows.
- **API:** `dnakit.variants.summarize_amplicon_stats(source, primers_bed, *, verbose=True)`
- **Input:** Alignment file and primer BED.
- **Output:** Returns amplicons, summary and count.

**Example code:**

```python
import shutil
import tempfile
from pathlib import Path
from dnakit.variants import select_alignment_records, summarize_amplicon_stats

work = Path(tempfile.mkdtemp())
shutil.copy("primers.bed", work / "primers.bed")
bam = select_alignment_records("reads.sam", work / "reads.bam")["output"]
result = summarize_amplicon_stats(bam, work / "primers.bed")
print(result["count"], result["amplicons"][0]["left"], result["amplicons"][0]["right"])
```

**Example output:**

```text
1 1-10 101-110
```

## 4) Export FASTQ from a BAM

- **Purpose:** Write BAM reads to four FASTQ files by read role; bases and qualities are unchanged.
- **Method:** Calls `samtools fastq -0 -1 -2 -s` with all four outputs named, so nothing goes to stdout.
- **API:** `dnakit.variants.extract_fastq_from_bam(source, output_directory, *, verbose=True)`
- **Input:** Alignment file and a new directory.
- **Output:** Returns per-role record counts (records), files and directory.

**Example code:**

```python
import tempfile
from pathlib import Path
from dnakit.variants import select_alignment_records

work = Path(tempfile.mkdtemp())
bam = select_alignment_records("reads.sam", work / "reads.bam")["output"]
from dnakit.variants import extract_fastq_from_bam

result = extract_fastq_from_bam(bam, work / "fastq")
print(result["records"])
```

**Example output:**

```text
{'other': 5, 'read1': 0, 'read2': 0, 'singletons': 0}
```

## 5) Export FASTA from a BAM

- **Purpose:** Write BAM sequences to four FASTA files by read role.
- **Method:** Calls `samtools fasta -0 -1 -2 -s` with every output named.
- **API:** `dnakit.variants.extract_fasta_from_bam(source, output_directory, *, verbose=True)`
- **Input:** Alignment file and a new directory.
- **Output:** Returns per-role record counts, files and directory.

**Example code:**

```python
import tempfile
from pathlib import Path
from dnakit.variants import select_alignment_records

work = Path(tempfile.mkdtemp())
bam = select_alignment_records("reads.sam", work / "reads.bam")["output"]
from dnakit.variants import extract_fasta_from_bam

result = extract_fasta_from_bam(bam, work / "fasta")
print(result["records"])
```

**Example output:**

```text
{'other': 5, 'read1': 0, 'read2': 0, 'singletons': 0}
```

## 6) Group reads by name

- **Purpose:** Keep reads with the same name adjacent, for later paired-read steps.
- **Method:** Calls `samtools collate -o`.
- **API:** `dnakit.variants.collate_bam_reads(source, output_path, *, force=False, verbose=True)`
- **Input:** Alignment file and a new output path.
- **Output:** Returns records and the output path.

**Example code:**

```python
import tempfile
from pathlib import Path
from dnakit.variants import select_alignment_records

work = Path(tempfile.mkdtemp())
bam = select_alignment_records("reads.sam", work / "reads.bam")["output"]
from dnakit.variants import collate_bam_reads

result = collate_bam_reads(bam, work / "collated.bam")
print(result["records"])
```

**Example output:**

```text
5
```

## 7) Concatenate BAM files

- **Purpose:** Concatenate BAM files with the same header, in input order.
- **Method:** Calls `samtools cat -o`.
- **API:** `dnakit.variants.concatenate_bam_files(sources, output_path, *, force=False, verbose=True)`
- **Input:** Two or more alignment files; the output is a new file.
- **Output:** Returns records, inputs and the output path.

**Example code:**

```python
import tempfile
from pathlib import Path
from dnakit.variants import select_alignment_records

work = Path(tempfile.mkdtemp())
bam = select_alignment_records("reads.sam", work / "reads.bam")["output"]
from dnakit.variants import concatenate_bam_files

result = concatenate_bam_files([bam, bam], work / "cat.bam")
print(result["records"])
```

**Example output:**

```text
10
```

## 8) Covered bases per BED interval

- **Purpose:** For each BED interval, sum the aligned bases of all reads that overlap it.
- **Method:** Calls `samtools bedcov` and needs a BAM index. A base covered by several reads is counted once per read.
- **API:** `dnakit.variants.count_bam_reads_per_interval(intervals_bed, source, *, verbose=True)`
- **Input:** BED intervals and an alignment file.
- **Output:** Returns chrom, start, end and covered_bases for each interval.

**Example code:**

```python
import tempfile
from pathlib import Path
from dnakit.variants import select_alignment_records

work = Path(tempfile.mkdtemp())
bam = select_alignment_records("reads.sam", work / "reads.bam")["output"]
import shutil
from dnakit.variants import count_bam_reads_per_interval, index_alignment_file

index_alignment_file(bam)
shutil.copy("intervals.bed", work / "intervals.bed")
result = count_bam_reads_per_interval(work / "intervals.bed", bam)
print([row["covered_bases"] for row in result["intervals"]])
```

**Example output:**

```text
[16, 28]
```

## 9) Recompute NM and MD tags

- **Purpose:** Recompute the NM and MD tags of each alignment against a reference FASTA.
- **Method:** Calls `samtools calmd -b`; the reference needs a .fai index first.
- **API:** `dnakit.variants.recompute_nm_md_tags(source, reference, output_path, *, force=False, verbose=True)`
- **Input:** Alignment file, reference FASTA and a new output path.
- **Output:** Returns records and the output path.

**Example code:**

```python
import tempfile
from pathlib import Path
from dnakit.variants import select_alignment_records

work = Path(tempfile.mkdtemp())
bam = select_alignment_records("reads.sam", work / "reads.bam")["output"]
import shutil
from dnakit.bioinformatics.reads import index_reference
from dnakit.variants import recompute_nm_md_tags

shutil.copy("ref.fa", work / "ref.fa")
index_reference(work / "ref.fa")
result = recompute_nm_md_tags(bam, work / "ref.fa", work / "calmd.bam")
print(result["records"])
```

**Example output:**

```text
5
```

## 10) Coverage per contig

- **Purpose:** Per-contig read count, covered bases and mean depth.
- **Method:** Calls `samtools coverage` and parses the table rows.
- **API:** `dnakit.variants.summarize_bam_coverage(source, *, verbose=True)`
- **Input:** Alignment file.
- **Output:** Returns the contigs list and count.

**Example code:**

```python
import tempfile
from pathlib import Path
from dnakit.variants import select_alignment_records

work = Path(tempfile.mkdtemp())
bam = select_alignment_records("reads.sam", work / "reads.bam")["output"]
from dnakit.variants import summarize_bam_coverage

result = summarize_bam_coverage(bam)
print({row["rname"]: row["covbases"] for row in result["contigs"]})
```

**Example output:**

```text
{'chr1': 22, 'chr2': 10}
```

## 11) Per-base read depth

- **Purpose:** Depth at each position, optionally limited to a region.
- **Method:** Calls `samtools depth -a`; a region query needs an index.
- **API:** `dnakit.variants.calculate_bam_read_depth(source, *, region=None, all_positions=True, max_positions=100000, verbose=True)`
- **Input:** Alignment file, optional region and a row limit.
- **Output:** Returns positions (contig, position, depth) and count.

**Example code:**

```python
import tempfile
from pathlib import Path
from dnakit.variants import select_alignment_records

work = Path(tempfile.mkdtemp())
bam = select_alignment_records("reads.sam", work / "reads.bam")["output"]
from dnakit.variants import calculate_bam_read_depth, index_alignment_file

index_alignment_file(bam)
result = calculate_bam_read_depth(bam, region="chr1:1-10")
print([row["depth"] for row in result["positions"]])
```

**Example output:**

```text
[1, 1, 1, 1, 2, 2, 2, 2, 2, 2]
```

## 12) Write a FASTA sequence dictionary

- **Purpose:** Write the reference .dict file and list each contig's name, length and MD5.
- **Method:** Calls `samtools dict`.
- **API:** `dnakit.variants.write_fasta_sequence_dictionary(reference, output_path, *, force=False, verbose=True)`
- **Input:** Reference FASTA and a new output path ending in .dict.
- **Output:** Returns contigs (name, length, md5) and the output path.

**Example code:**

```python
import shutil
import tempfile
from pathlib import Path
from dnakit.variants import write_fasta_sequence_dictionary

work = Path(tempfile.mkdtemp())
shutil.copy("ref.fa", work / "ref.fa")
result = write_fasta_sequence_dictionary(work / "ref.fa", work / "ref.dict")
print([(item["name"], item["length"]) for item in result["contigs"]])
```

**Example output:**

```text
[('chr1', 120), ('chr2', 80)]
```

## 13) Decode SAM FLAG values

- **Purpose:** Name the bits of each FLAG integer, such as 99 = PAIRED, PROPER_PAIR, ...
- **Method:** Calls `samtools flags`; names are sorted.
- **API:** `dnakit.variants.decode_sam_flags(flags, *, verbose=True)`
- **Input:** A list of integers from 0 to 65535.
- **Output:** Returns flag, hex and names for each FLAG.

**Example code:**

```python
from dnakit.variants import decode_sam_flags

result = decode_sam_flags([99, 4])
print([(item["flag"], item["names"]) for item in result["flags"]])
```

**Example output:**

```text
[(99, ['MREVERSE', 'PAIRED', 'PROPER_PAIR', 'READ1']), (4, ['UNMAP'])]
```

## 14) Read the header and first records

- **Purpose:** Return the header lines and the first N records of an alignment file.
- **Method:** Calls `samtools head -n`.
- **API:** `dnakit.variants.read_bam_header_and_records(source, *, records=5, verbose=True)`
- **Input:** Alignment file and record count (0 to 1000).
- **Output:** Returns header and records (SAM columns).

**Example code:**

```python
import tempfile
from pathlib import Path
from dnakit.variants import select_alignment_records

work = Path(tempfile.mkdtemp())
bam = select_alignment_records("reads.sam", work / "reads.bam")["output"]
from dnakit.variants import read_bam_header_and_records

result = read_bam_header_and_records(bam, records=2)
print([row[0] for row in result["records"]], result["header"][0])
```

**Example output:**

```text
['r1', 'r2'] @HD	VN:1.6	SO:coordinate
```

## 15) Index a BAM or CRAM

- **Purpose:** Write a BAI (or CSI) index for region queries.
- **Method:** Calls `samtools index`; the index is written beside the input.
- **API:** `dnakit.variants.index_alignment_file(source, *, csi=False, force=False, verbose=True)`
- **Input:** A coordinate-sorted alignment file.
- **Output:** Returns the index path and format.

**Example code:**

```python
import tempfile
from pathlib import Path
from dnakit.variants import select_alignment_records

work = Path(tempfile.mkdtemp())
bam = select_alignment_records("reads.sam", work / "reads.bam")["output"]
from dnakit.variants import index_alignment_file

print(index_alignment_file(bam)["format"], index_alignment_file(bam, csi=True)["format"])
```

**Example output:**

```text
BAI CSI
```

## 16) Mapped and unmapped counts per contig

- **Purpose:** Read mapped and unmapped counts per contig from the index.
- **Method:** Calls `samtools idxstats`; the index must exist.
- **API:** `dnakit.variants.summarize_bam_index_stats(source, *, verbose=True)`
- **Input:** An indexed alignment file.
- **Output:** Returns contigs, unplaced_unmapped and provenance.

**Example code:**

```python
import tempfile
from pathlib import Path
from dnakit.variants import select_alignment_records

work = Path(tempfile.mkdtemp())
bam = select_alignment_records("reads.sam", work / "reads.bam")["output"]
from dnakit.variants import index_alignment_file, summarize_bam_index_stats

index_alignment_file(bam)
result = summarize_bam_index_stats(bam)
print([(row["contig"], row["mapped"]) for row in result["contigs"]], result["unplaced_unmapped"])
```

**Example output:**

```text
[('chr1', 3), ('chr2', 1)] 1
```

## 17) Mark duplicate reads

- **Purpose:** Set flag 0x400 on duplicate reads without removing them.
- **Method:** Calls `samtools markdup`; the input must be coordinate-sorted.
- **API:** `dnakit.variants.mark_bam_duplicates(source, output_path, *, force=False, verbose=True)`
- **Input:** A coordinate-sorted alignment file and a new output path.
- **Output:** Returns records, duplicates and the output path.

**Example code:**

```python
import tempfile
from pathlib import Path
from dnakit.variants import select_alignment_records

work = Path(tempfile.mkdtemp())
bam = select_alignment_records("reads.sam", work / "reads.bam")["output"]
from dnakit.variants import mark_bam_duplicates

result = mark_bam_duplicates(bam, work / "marked.bam")
print(result["records"], result["duplicates"])
```

**Example output:**

```text
5 0
```

## 18) Merge BAM files

- **Purpose:** Merge coordinate-sorted BAM files; the output stays sorted.
- **Method:** Calls `samtools merge`.
- **API:** `dnakit.variants.merge_bam_files(sources, output_path, *, force=False, verbose=True)`
- **Input:** Two or more alignment files and a new output path.
- **Output:** Returns records, inputs and the output path.

**Example code:**

```python
import tempfile
from pathlib import Path
from dnakit.variants import select_alignment_records

work = Path(tempfile.mkdtemp())
bam = select_alignment_records("reads.sam", work / "reads.bam")["output"]
from dnakit.variants import merge_bam_files

result = merge_bam_files([bam, bam], work / "merged.bam")
print(result["records"])
```

**Example output:**

```text
10
```

## 19) Per-position base pileup text

- **Purpose:** Write the base pileup text for each covered position, as input to variant calling.
- **Method:** Calls `samtools mpileup`; a reference needs a .fai index.
- **API:** `dnakit.variants.pileup_alignment_positions(source, output_path, *, reference=None, force=False, verbose=True)`
- **Input:** Alignment file, optional reference FASTA and a new output path.
- **Output:** Returns positions (covered positions) and the output path.

**Example code:**

```python
import shutil
import tempfile
from pathlib import Path
from dnakit.bioinformatics.reads import index_reference
from dnakit.variants import select_alignment_records
from dnakit.variants import pileup_alignment_positions

work = Path(tempfile.mkdtemp())
bam = select_alignment_records("reads.sam", work / "reads.bam")["output"]
shutil.copy("ref.fa", work / "ref.fa")
index_reference(work / "ref.fa")
result = pileup_alignment_positions(bam, work / "pileup.txt", reference=work / "ref.fa")
print(result["positions"])
```

**Example output:**

```text
32
```

## 20) Fill paired-read mate fields

- **Purpose:** Fill mate position, insert size and mate tags of name-sorted paired reads.
- **Method:** Calls `samtools fixmate -m`; the input must be name-sorted or an error is raised.
- **API:** `dnakit.variants.fix_bam_mate_fields(source, output_path, *, force=False, verbose=True)`
- **Input:** A name-sorted alignment file and a new output path.
- **Output:** Returns records and the output path.

**Example code:**

```python
import tempfile
from pathlib import Path
from dnakit.variants import select_alignment_records

work = Path(tempfile.mkdtemp())
bam = select_alignment_records("reads.sam", work / "reads.bam")["output"]
from dnakit.variants import fix_bam_mate_fields, sort_alignment_file

named = sort_alignment_file(bam, work / "byname.bam", by_name=True)["output"]
print(fix_bam_mate_fields(named, work / "fixed.bam")["records"])
```

**Example output:**

```text
5
```

## 21) Check alignment files are readable

- **Purpose:** Check that each BAM or CRAM opens and has a valid header.
- **Method:** Calls `samtools quickcheck` on each file separately.
- **API:** `dnakit.variants.check_alignment_file_integrity(sources, *, verbose=True)`
- **Input:** A list of local files.
- **Output:** Returns the valid state of each file and the valid and invalid counts.

**Example code:**

```python
import tempfile
from pathlib import Path
from dnakit.variants import select_alignment_records

work = Path(tempfile.mkdtemp())
bam = select_alignment_records("reads.sam", work / "reads.bam")["output"]
from pathlib import Path
from dnakit.variants import check_alignment_file_integrity

(work / "junk.bam").write_text("not a bam", encoding="utf-8")
result = check_alignment_file_integrity([bam, work / "junk.bam"])
print([(Path(item["path"]).name, item["valid"]) for item in result["files"]])
```

**Example output:**

```text
[('reads.bam', True), ('junk.bam', False)]
```

## 22) Replace the BAM header

- **Purpose:** Replace only the header; records are copied unchanged.
- **Method:** Calls `samtools reheader`; the result is written to a new file.
- **API:** `dnakit.variants.replace_bam_header(header_sam, source, output_path, *, force=False, verbose=True)`
- **Input:** A SAM header file starting with @ lines, an alignment file and a new output path.
- **Output:** Returns records and the output path.

**Example code:**

```python
import tempfile
from pathlib import Path
from dnakit.variants import select_alignment_records

work = Path(tempfile.mkdtemp())
bam = select_alignment_records("reads.sam", work / "reads.bam")["output"]
import shutil
from dnakit.variants import replace_bam_header

shutil.copy("header.sam", work / "header.sam")
print(replace_bam_header(work / "header.sam", bam, work / "rh.bam")["records"])
```

**Example output:**

```text
5
```

## 23) Clear alignment fields

- **Purpose:** Clear position, CIGAR and mapping flags so reads become unmapped; sequences, qualities and read-group tags are kept.
- **Method:** Calls `samtools reset -o`.
- **API:** `dnakit.variants.reset_alignments_to_unmapped(source, output_path, *, force=False, verbose=True)`
- **Input:** Alignment file and a new output path.
- **Output:** Returns records and the output path.

**Example code:**

```python
import tempfile
from pathlib import Path
from dnakit.variants import select_alignment_records

work = Path(tempfile.mkdtemp())
bam = select_alignment_records("reads.sam", work / "reads.bam")["output"]
from dnakit.variants import reset_alignments_to_unmapped

print(reset_alignments_to_unmapped(bam, work / "reset.bam")["records"])
```

**Example output:**

```text
5
```

## 24) List read-group sample names

- **Purpose:** List the SM sample names from the @RG lines of each alignment file.
- **Method:** Calls `samtools samples`.
- **API:** `dnakit.variants.list_read_group_samples(sources, *, verbose=True)`
- **Input:** One or more alignment files.
- **Output:** Returns the samples list (sample and source).

**Example code:**

```python
import tempfile
from pathlib import Path
from dnakit.variants import select_alignment_records

work = Path(tempfile.mkdtemp())
bam = select_alignment_records("reads.sam", work / "reads.bam")["output"]
from dnakit.variants import list_read_group_samples

print([item["sample"] for item in list_read_group_samples([bam])["samples"]])
```

**Example output:**

```text
['sampleA']
```

## 25) Sort an alignment file

- **Purpose:** Sort by coordinate, or by read name with by_name=True.
- **Method:** Calls `samtools sort -o`, with `-n` when by_name is True.
- **API:** `dnakit.variants.sort_alignment_file(source, output_path, *, by_name=False, force=False, verbose=True)`
- **Input:** Alignment file and a new output path.
- **Output:** Returns order, records and the output path.

**Example code:**

```python
import tempfile
from pathlib import Path
from dnakit.variants import select_alignment_records

work = Path(tempfile.mkdtemp())
bam = select_alignment_records("reads.sam", work / "reads.bam")["output"]
from dnakit.variants import sort_alignment_file

print(sort_alignment_file(bam, work / "coord.bam")["order"], sort_alignment_file(bam, work / "name.bam", by_name=True)["order"])
```

**Example output:**

```text
coordinate queryname
```

## 26) Split a BAM by read group

- **Purpose:** Write one BAM per read group; every read must carry an RG tag.
- **Method:** Calls `samtools split -f`; file names are the read-group IDs.
- **API:** `dnakit.variants.split_bam_by_read_group(source, output_directory, *, verbose=True)`
- **Input:** Alignment file and a new directory.
- **Output:** Returns each file name and record count.

**Example code:**

```python
import tempfile
from pathlib import Path
from dnakit.variants import select_alignment_records

work = Path(tempfile.mkdtemp())
bam = select_alignment_records("reads.sam", work / "reads.bam")["output"]
from dnakit.variants import split_bam_by_read_group

result = split_bam_by_read_group(bam, work / "split")
print([(item["file"], item["records"]) for item in result["files"]])
```

**Example output:**

```text
[('rg1.bam', 5)]
```

## 27) Alignment summary statistics

- **Purpose:** Take read, mapping and length totals from the SN lines of samtools stats.
- **Method:** Calls `samtools stats` and parses the SN lines.
- **API:** `dnakit.variants.summarize_alignment_stats(source, *, verbose=True)`
- **Input:** Alignment file.
- **Output:** Returns the summary dictionary.

**Example code:**

```python
import tempfile
from pathlib import Path
from dnakit.variants import select_alignment_records

work = Path(tempfile.mkdtemp())
bam = select_alignment_records("reads.sam", work / "reads.bam")["output"]
from dnakit.variants import summarize_alignment_stats

print(summarize_alignment_stats(bam)["summary"]["raw total sequences"])
```

**Example output:**

```text
5
```

## 28) Select alignment records

- **Purpose:** Keep reads by region, FLAG masks and minimum mapping quality.
- **Method:** Calls `samtools view`; a region needs an index. A .sam output is text, anything else is BAM.
- **API:** `dnakit.variants.select_alignment_records(source, output_path, *, region=None, required_flags=0, excluded_flags=0, min_mapping_quality=0, force=False, verbose=True)`
- **Input:** Alignment file, optional region, FLAG masks and a new output path.
- **Output:** Returns records and the output path.

**Example code:**

```python
import tempfile
from pathlib import Path
from dnakit.variants import select_alignment_records

work = Path(tempfile.mkdtemp())
bam = select_alignment_records("reads.sam", work / "reads.bam")["output"]
from dnakit.variants import index_alignment_file, select_alignment_records

index_alignment_file(bam)
print(select_alignment_records(bam, work / "sel.bam", region="chr1:1-10", min_mapping_quality=30)["records"])
```

**Example output:**

```text
2
```

## 29) Rebuild reference bases from alignments

- **Purpose:** Rebuild reference bases from aligned reads with MD tags; uncovered positions are N.
- **Method:** Calls `samtools reference`; without MD tags the result is almost all N.
- **API:** `dnakit.variants.reconstruct_reference_from_alignments(source, output_path, *, force=False, verbose=True)`
- **Input:** Alignment file with MD tags and a new FASTA path.
- **Output:** Returns each contig's length and number of determined bases.

**Example code:**

```python
import shutil
import tempfile
from pathlib import Path
from dnakit.bioinformatics.reads import index_reference
from dnakit.variants import recompute_nm_md_tags, reconstruct_reference_from_alignments, select_alignment_records

work = Path(tempfile.mkdtemp())
bam = select_alignment_records("reads.sam", work / "reads.bam")["output"]
shutil.copy("ref.fa", work / "ref.fa")
index_reference(work / "ref.fa")
tagged = recompute_nm_md_tags(bam, work / "ref.fa", work / "calmd.bam")["output"]
result = reconstruct_reference_from_alignments(tagged, work / "rebuilt.fa")
print({item["name"]: item["determined_bases"] for item in result["contigs"]})
```

**Example output:**

```text
{'chr1': 22, 'chr2': 10}
```

## 30) Split reads by haplotype

- **Purpose:** Tag reads with HP haplotypes and write them to separate files.
- **Method:** Calls `samtools phase -k -b`; writes prefix.0.bam, prefix.1.bam and prefix.chimera.bam.
- **API:** `dnakit.variants.phase_bam_haplotypes(source, output_prefix, *, block_length=13, verbose=True)`
- **Input:** Alignment file and an output prefix; existing files under the prefix are refused.
- **Output:** Returns each output file name and record count.

**Example code:**

```python
import tempfile
from pathlib import Path
from dnakit.variants import select_alignment_records

work = Path(tempfile.mkdtemp())
bam = select_alignment_records("reads.sam", work / "reads.bam")["output"]
from dnakit.variants import phase_bam_haplotypes

result = phase_bam_haplotypes(bam, work / "phased")
print(sorted(item["file"].rsplit("/", 1)[-1] for item in result["files"]))
```

**Example output:**

```text
['phased.0.bam', 'phased.1.bam', 'phased.chimera.bam']
```

## 31) Index a FASTQ file

- **Purpose:** Write a .fai index with one row per record.
- **Method:** Calls `samtools fqidx`; the index is written beside the FASTQ.
- **API:** `dnakit.variants.index_fastq_file(source, *, force=False, verbose=True)`
- **Input:** A FASTQ file.
- **Output:** Returns the index path and entry count.

**Example code:**

```python
import shutil
import tempfile
from pathlib import Path
from dnakit.variants import index_fastq_file

work = Path(tempfile.mkdtemp())
shutil.copy("reads.fq", work / "reads.fq")
print(index_fastq_file(work / "reads.fq")["entries"])
```

**Example output:**

```text
2
```

## 32) Import FASTQ as an unaligned BAM

- **Purpose:** Convert one FASTQ (single-end) or two FASTQ files (paired) to an unaligned BAM.
- **Method:** Calls `samtools import -o`.
- **API:** `dnakit.variants.import_fastq_to_bam(sources, output_path, *, force=False, verbose=True)`
- **Input:** One or two FASTQ files and a new output path.
- **Output:** Returns records, layout and the output path.

**Example code:**

```python
import shutil
import tempfile
from pathlib import Path
from dnakit.variants import import_fastq_to_bam

work = Path(tempfile.mkdtemp())
shutil.copy("reads.fq", work / "reads.fq")
result = import_fastq_to_bam([work / "reads.fq"], work / "imported.bam")
print(result["records"], result["layout"])
```

**Example output:**

```text
2 single
```

## bcftools variant-file functions

## 33) Filter VCF sites by expression

- **Purpose:** Keep (include) or drop (exclude) sites that match a bcftools expression.
- **Method:** Calls `bcftools filter -i` or `-e`; the expression is one argument, not run through a shell.
- **API:** `dnakit.variants.filter_vcf_by_expression(source, expression, output_path, *, mode='include', force=False, verbose=True)`
- **Input:** A bgzipped, indexed VCF, an expression and a new output path.
- **Output:** Returns records (sites kept) and the output path.

**Example code:**

```python
import shutil
import tempfile
from pathlib import Path
from dnakit.variants import index_tabix_file

work = Path(tempfile.mkdtemp())
shutil.copy("calls.vcf", work / "calls.vcf")
gz = index_tabix_file(work / "calls.vcf")["compressed"]
from dnakit.variants import filter_vcf_by_expression

print(filter_vcf_by_expression(gz, "QUAL>45", work / "kept.vcf")["records"])
```

**Example output:**

```text
2
```

## 34) Annotate VCF protein consequences

- **Purpose:** Add BCSQ consequence tags from a GFF3 gene model and a reference FASTA.
- **Method:** Calls `bcftools csq -g -f`. The input must be phased; bcftools stops on unphased heterozygous sites.
- **API:** `dnakit.variants.annotate_vcf_consequence(source, gene_annotation, reference, output_path, *, force=False, verbose=True)`
- **Input:** A phased VCF, a GFF3 file, a reference FASTA and a new output path.
- **Output:** Returns records, annotated_records and the output path.

**Example code:**

```python
import shutil
import tempfile
from pathlib import Path
from dnakit.variants import annotate_vcf_consequence

work = Path(tempfile.mkdtemp())
for name in ("phased.vcf", "genes.gff3", "ref.fa"):
    shutil.copy(name, work / name)
result = annotate_vcf_consequence(work / "phased.vcf", work / "genes.gff3", work / "ref.fa", work / "csq.vcf")
print(result["annotated_records"])
```

**Example output:**

```text
3
```

## 35) Check genotype concordance between samples

- **Purpose:** Compare genotype discordance for each sample pair; by default the VCF is compared with itself.
- **Method:** Calls `bcftools gtcheck -g`.
- **API:** `dnakit.variants.check_vcf_sample_concordance(source, *, genotype_source=None, verbose=True)`
- **Input:** A VCF and an optional second genotype source.
- **Output:** Returns pairs (query, genotyped, discordance, sites, matching).

**Example code:**

```python
import shutil
import tempfile
from pathlib import Path
from dnakit.variants import index_tabix_file

work = Path(tempfile.mkdtemp())
shutil.copy("calls.vcf", work / "calls.vcf")
gz = index_tabix_file(work / "calls.vcf")["compressed"]
from dnakit.variants import check_vcf_sample_concordance

result = check_vcf_sample_concordance(gz)
print([(item["query"], item["genotyped"], item["discordance"]) for item in result["pairs"]])
```

**Example output:**

```text
[('S1', 'S1', 0.0), ('S1', 'S2', 27.63102), ('S2', 'S1', 27.63102), ('S2', 'S2', 0.0)]
```

## 36) Read a VCF header

- **Purpose:** Return the header lines and sample names of a VCF or BCF.
- **Method:** Calls `bcftools head`.
- **API:** `dnakit.variants.read_vcf_header(source, *, verbose=True)`
- **Input:** A VCF or BCF file.
- **Output:** Returns header_lines and samples.

**Example code:**

```python
import shutil
import tempfile
from pathlib import Path
from dnakit.variants import read_vcf_header

work = Path(tempfile.mkdtemp())
shutil.copy("calls.vcf", work / "calls.vcf")
print(read_vcf_header(work / "calls.vcf")["samples"])
```

**Example output:**

```text
['S1', 'S2']
```

## 37) Index a VCF with CSI

- **Purpose:** Write a CSI index for a bgzipped VCF or a BCF.
- **Method:** Calls `bcftools index -c`.
- **API:** `dnakit.variants.index_vcf_csi(source, *, force=False, verbose=True)`
- **Input:** A file ending in .vcf.gz or .bcf.
- **Output:** Returns the index path.

**Example code:**

```python
import shutil
import tempfile
from pathlib import Path
from dnakit.variants import index_tabix_file

work = Path(tempfile.mkdtemp())
shutil.copy("calls.vcf", work / "calls.vcf")
gz = index_tabix_file(work / "calls.vcf")["compressed"]
from dnakit.variants import index_vcf_csi

print(Path(index_vcf_csi(gz)["index"]).name)
```

**Example output:**

```text
calls.vcf.gz.csi
```

## 38) Merge samples of indexed VCFs

- **Purpose:** Join VCFs with different samples into one multi-sample VCF.
- **Method:** Calls `bcftools merge`; each input needs a .tbi or .csi index.
- **API:** `dnakit.variants.merge_indexed_vcfs(sources, output_path, *, force=False, verbose=True)`
- **Input:** Two or more bgzipped VCFs and a new output path.
- **Output:** Returns samples, records and the output path.

**Example code:**

```python
import shutil
import tempfile
from pathlib import Path
from dnakit.variants import index_tabix_file
from dnakit.variants import merge_indexed_vcfs

work = Path(tempfile.mkdtemp())
shutil.copy("calls.vcf", work / "calls.vcf")
shutil.copy("calls2.vcf", work / "calls2.vcf")
first = index_tabix_file(work / "calls.vcf")["compressed"]
second = index_tabix_file(work / "calls2.vcf")["compressed"]
print(merge_indexed_vcfs([first, second], work / "merged.vcf")["samples"])
```

**Example output:**

```text
['S1', 'S2', 'S3', 'S4']
```

## 39) BAM to VCF likelihoods

- **Purpose:** Write genotype likelihoods (PL) of a BAM at each site, an intermediate before calling.
- **Method:** Calls `bcftools mpileup -f -Ov`; the reference needs a .fai index.
- **API:** `dnakit.variants.pileup_bam_to_vcf(source, reference, output_path, *, force=False, verbose=True)`
- **Input:** Alignment file, reference FASTA and a new output path.
- **Output:** Returns records and the output path.

**Example code:**

```python
import shutil
import tempfile
from pathlib import Path
from dnakit.bioinformatics.reads import index_reference
from dnakit.variants import pileup_bam_to_vcf, select_alignment_records

work = Path(tempfile.mkdtemp())
bam = select_alignment_records("reads.sam", work / "reads.bam")["output"]
shutil.copy("ref.fa", work / "ref.fa")
index_reference(work / "ref.fa")
print(pileup_bam_to_vcf(bam, work / "ref.fa", work / "pileup.vcf")["records"])
```

**Example output:**

```text
32
```

## 40) Rename VCF samples

- **Purpose:** Replace the sample names with the names in a file; everything else is unchanged.
- **Method:** Calls `bcftools reheader -s`.
- **API:** `dnakit.variants.rename_vcf_samples(source, sample_names, output_path, *, force=False, verbose=True)`
- **Input:** A VCF, a file with one name per line (names must be unique) and a new output path.
- **Output:** Returns samples and the output path.

**Example code:**

```python
import shutil
import tempfile
from pathlib import Path
from dnakit.variants import rename_vcf_samples

work = Path(tempfile.mkdtemp())
shutil.copy("calls.vcf", work / "calls.vcf")
shutil.copy("names.txt", work / "names.txt")
print(rename_vcf_samples(work / "calls.vcf", work / "names.txt", work / "renamed.vcf")["samples"])
```

**Example output:**

```text
['A1', 'A2']
```

## 41) Detect runs of homozygosity with bcftools

- **Purpose:** Detect runs of homozygosity (ROH) for one sample with bcftools' HMM.
- **Method:** Calls `bcftools roh -G`; only GT fields are needed. Small inputs may report no runs.
- **API:** `dnakit.variants.detect_bcftools_roh(source, *, sample, genotype_quality=30.0, verbose=True)`
- **Input:** A bgzipped VCF, a sample name and the PL stand-in value (0 to 1000).
- **Output:** Returns runs (chromosome, start, end, length, markers, quality) and count.

**Example code:**

```python
import shutil
import tempfile
from pathlib import Path
from dnakit.variants import index_tabix_file
from dnakit.variants import detect_bcftools_roh

work = Path(tempfile.mkdtemp())
shutil.copy("calls.vcf", work / "calls.vcf")
gz = index_tabix_file(work / "calls.vcf")["compressed"]
result = detect_bcftools_roh(gz, sample="S1")
print(result["sample"], result["count"])
```

**Example output:**

```text
S1 0
```

## hgvs helper functions

## 42) Trim the common suffix of alleles

- **Purpose:** Remove the suffix shared by allele strings and report how many bases were removed.
- **Method:** Calls hgvs `trim_common_suffixes`.
- **API:** `dnakit.variants.trim_common_suffix_strings(strings, *, min_length=0, verbose=True)`
- **Input:** At least two allele strings made of A, C, G, T or N.
- **Output:** Returns trimmed_length and strings.

**Example code:**

```python
from dnakit.variants import trim_common_suffix_strings

result = trim_common_suffix_strings(["GGCTA", "AGCTA"])
print(result["trimmed_length"], result["strings"])
```

**Example output:**

```text
4 ['G', 'A']
```

## 43) 0-based half-open span of an HGVS variant

- **Purpose:** Parse one HGVS string (no reference data) and return its 0-based, half-open start and end.
- **Method:** Parses with the hgvs parser, then calls `get_start_end_interbase`.
- **API:** `dnakit.variants.hgvs_interbase_span(variant, *, outer_confidence=True, verbose=True)`
- **Input:** One HGVS string such as NC_000001.11:g.100_102del.
- **Output:** Returns start and end (0-based, half-open).

**Example code:**

```python
from dnakit.variants import hgvs_interbase_span

result = hgvs_interbase_span("NC_000001.11:g.100_102del")
print(result["start"], result["end"])
```

**Example output:**

```text
99 102
```

## 44) Repeat notation for alleles

- **Purpose:** Format HGVS-style repeat notation, such as AC[3], for one allele from its repeat blocks.
- **Method:** Calls hgvs `get_repeat_str`; the result is None when no notation applies.
- **API:** `dnakit.variants.repeat_notation_for_alleles(sequence, other_sequence, blocks, other_blocks=(), *, reverse=False, verbose=True)`
- **Input:** Two sequences and their repeat blocks (the output of detect_repeat_blocks).
- **Output:** Returns the formatted string, or None.

**Example code:**

```python
from dnakit.variants import detect_repeat_blocks
from dnakit.variants import repeat_notation_for_alleles

blocks = detect_repeat_blocks("ACACAC")["blocks"]
print(repeat_notation_for_alleles("ACACAC", "ACAC", blocks)["formatted"])
```

**Example output:**

```text
AC[3]
```
