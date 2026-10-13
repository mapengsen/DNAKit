# 70 Full G4: polars-bio, gtars, ncls, and genomicranges annotation functions

This group wraps the DNA- and gene-level public functions of polars-bio, gtars, ncls, and genomicranges that run offline on local files. Each function reads only the local files it is given and returns a JSON-serializable dictionary with `status`, `operation`, result fields, and `provenance`. Coordinates use the BED convention (0-based, half-open) unless a section says otherwise.

## Notes

**Import path:** the functions live in `dnakit.annotation`, `dnakit.annotation`, and `dnakit.annotation`; the examples import from these modules.
**Dependencies:** a missing package raises an installation hint when a function is called; polars-bio, gtars, genomicranges, and ncls are imported lazily. The `read_cooler_matrix` example also needs `cooler` and `pandas`.
**gtars:** the upstream functions live in the `gtars.genomic_distributions`, `gtars.lola`, `gtars.refget`, `gtars.tokenizers`, `gtars.vrs`, and `gtars.vrs.hgvs` submodules, so the top-level namespace exports none and a top-level scan produced an empty inventory.
**ncls:** the two index constructors in the inventory are the integer index (already wrapped in `dnakit.annotation.backends`) and the float index (its coordinates are not base positions, so it is not wrapped); this group adds no ncls function.
**genomicranges:** the upstream GTF reader stores the GTF end minus one, which `read_gtf_ranges` restores; `compute_up_down` returns widths rather than ends, which `upstream_downstream_regions` converts.

## Example data

The examples share these files. Files that example code writes (such as `copy.vcf` and `contacts.cool`) are removed before each run.

**Example file:** `a.bed`

```text
chr1	100	200	a
chr1	150	300	b
chr1	400	500	c
chr2	0	50	d
```

**Example file:** `b.bed`

```text
chr1	180	220	x
chr1	450	460	y
chr2	10	20	z
```

**Example file:** `genes.gtf`

```text
chr1	ens	gene	1001	2000	.	+	.	gene_id "G1"; gene_name "ABC"; gene_type "protein_coding";
chr1	ens	exon	1001	1200	.	+	.	gene_id "G1"; transcript_id "T1"; gene_type "protein_coding";
chr1	ens	exon	1801	2000	.	+	.	gene_id "G1"; transcript_id "T1"; gene_type "protein_coding";
chr1	ens	gene	5001	6000	.	-	.	gene_id "G2"; gene_name "DEF"; gene_type "protein_coding";
chr1	ens	exon	5001	6000	.	-	.	gene_id "G2"; transcript_id "T2"; gene_type "protein_coding";
```

**Example file:** `genes.gff3`

```text
##gff-version 3
chr1	src	gene	100	500	.	+	.	ID=gene1;Name=GENE1
chr1	src	mRNA	100	500	.	+	.	ID=tx1;Parent=gene1
chr1	src	exon	100	200	.	+	.	ID=exon1;Parent=tx1
chr1	src	gene	800	900	.	-	.	ID=gene2;Name=GENE2
```

**Example file:** `variants.vcf`

```text
##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s1	s2
chr1	10	rs1	A	G	50	PASS	DP=12	GT	0/1	1/1
chr1	25	.	C	T	40	PASS	DP=8	GT	0/0	0/1
```

**Example file:** `genome.fa`

```text
>chr1 first
ACGTACGTAC
GGTT
>chr2
TTTTAAAACCCC
```

**Example file:** `reads.fq`

```text
@r1
ACGTN
+
IIII#
@r2
GGCCAA
+
IIIIII
```

**Example file:** `aln.sam`

```text
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:1000
r1	0	chr1	11	60	5M	*	0	0	ACGTA	IIIII	NM:i:0
r2	0	chr1	13	60	4M	*	0	0	GTAC	IIII	NM:i:1
```

**Example file:** `contacts.pairs`

```text
## pairs format v1.0
#columns: readID chr1 pos1 chr2 pos2 strand1 strand2
r1	chr1	100	chr1	500	+	-
r2	chr1	200	chr2	40	+	+
```

**Example file:** `align.sto`

```text
# STOCKHOLM 1.0
#=GF ID test_family
seq1  ACGT-ACGT
seq2  ACGTTACGA
//
```

**Example file:** `gc.bed`

```text
chr1	0	10
chr1	10	14
chr2	0	4
chr2	8	12
```

**Example file:** `regions.bed`

```text
chr1	950	1010	p1
chr1	1500	1600	p2
chr1	5950	6100	p3
chr1	9000	9100	p4
```

**Example file:** `signal.tsv`

```text
region	condA	condB
chr1_0_100	1.0	2.0
chr2_0_100	3.0	4.0
```

**Example file:** `frags.tsv`

```text
chr1	10	20	AAA	1
chr1	150	160	AAA	2
chr2	5	9	BBB	3
```

**Example file:** `universe.bed`

```text
chr1	0	100
chr1	100	200
chr2	0	50
```

**Example file:** `user.bed`

```text
chr1	0	10
chr1	20	30
```

**Example file:** `db1.bed`

```text
chr1	0	10
```

**Example file:** `db2.bed`

```text
chr1	20	30
chr2	0	50
```

**Example file:** `set1.bed`

```text
chr1	100	200
chr1	300	400
```

**Example file:** `set2.bed`

```text
chr1	100	200
chr1	500	600
```

## polars-bio: intervals, gene annotation, variants, sequences, and alignments

### 1) `dnakit.annotation.overlap_intervals()`

- **Purpose:** List every pair of overlapping intervals from two BED files.
- **Method:** Read both files with polars-bio (0-based, half-open), call `overlap`, and sort the rows by coordinate.
- **API:** `overlap_intervals(left, right, *, max_rows=100000, verbose=True)`
- **Input:** Two BED file paths; `max_rows` limits the number of output rows.
- **Output:** `pairs`, one record per overlapping pair with `chrom_*`, `start_*`, `end_*`, and `name_*` for each side, plus `pair_count`.
- **Example code:**

```python
from dnakit.annotation import overlap_intervals

result = overlap_intervals("a.bed", "b.bed")
print(result["pair_count"])
print([(row["start_1"], row["end_1"], row["start_2"], row["end_2"]) for row in result["pairs"]])
```

- **Example output:**

```text
4
[(100, 200, 180, 220), (150, 300, 180, 220), (400, 500, 450, 460), (0, 50, 10, 20)]
```

### 2) `dnakit.annotation.count_interval_overlaps()`

- **Purpose:** Count how many intervals of `other` overlap each interval of `source`.
- **Method:** Call polars-bio `count_overlaps`; each source interval keeps one row, and an interval with no overlap has count 0.
- **API:** `count_interval_overlaps(source, other, *, max_rows=100000, verbose=True)`
- **Input:** Source and comparison BED paths; `max_rows` limits output rows.
- **Output:** `intervals` with a `count` per source interval, plus `total_overlaps` and `interval_count`.
- **Example code:**

```python
from dnakit.annotation import count_interval_overlaps

result = count_interval_overlaps("a.bed", "b.bed")
print(result["total_overlaps"])
print([(row["chrom"], row["start"], row["count"]) for row in result["intervals"]])
```

- **Example output:**

```text
4
[('chr1', 100, 1), ('chr1', 150, 1), ('chr1', 400, 1), ('chr2', 0, 1)]
```

### 3) `dnakit.annotation.cluster_intervals()`

- **Purpose:** Group intervals that overlap or lie within `min_distance` bases into clusters.
- **Method:** Call polars-bio `cluster`, then renumber clusters in coordinate order so the result does not depend on processing order.
- **API:** `cluster_intervals(source, *, min_distance=0, max_rows=100000, verbose=True)`
- **Input:** A BED file; `min_distance` is a non-negative base count.
- **Output:** `intervals` with `cluster`, `cluster_start`, and `cluster_end`, plus `cluster_count`.
- **Example code:**

```python
from dnakit.annotation import cluster_intervals

result = cluster_intervals("a.bed", min_distance=0)
print(result["cluster_count"])
print([(row["name"], row["cluster"]) for row in result["intervals"]])
```

- **Example output:**

```text
3
[('a', 0), ('b', 0), ('c', 1), ('d', 2)]
```

### 4) `dnakit.annotation.complement_intervals()`

- **Purpose:** Return the gaps that no interval covers on each chromosome.
- **Method:** Build view intervals [0, length) from `chrom_sizes` and call polars-bio `complement` to take the difference.
- **API:** `complement_intervals(source, *, chrom_sizes, max_rows=100000, verbose=True)`
- **Input:** A BED file; `chrom_sizes` maps chromosome names to lengths.
- **Output:** `gaps` with `chrom`, `start`, and `end`, plus `gap_count`.
- **Example code:**

```python
from dnakit.annotation import complement_intervals

result = complement_intervals("a.bed", chrom_sizes={"chr1": 1000, "chr2": 100})
print(result["gap_count"])
print([(row["chrom"], row["start"], row["end"]) for row in result["gaps"]])
```

- **Example output:**

```text
4
[('chr1', 0, 100), ('chr1', 300, 400), ('chr1', 500, 1000), ('chr2', 50, 100)]
```

### 5) `dnakit.annotation.subtract_intervals()`

- **Purpose:** Remove the bases covered by `other` from the intervals of `source`.
- **Method:** Call polars-bio `subtract`; partially covered intervals are split into the remaining pieces.
- **API:** `subtract_intervals(source, other, *, max_rows=100000, verbose=True)`
- **Input:** Source and subtrahend BED paths; `max_rows` limits output rows.
- **Output:** `pieces` keeping the source columns, plus `piece_count`.
- **Example code:**

```python
from dnakit.annotation import subtract_intervals

result = subtract_intervals("a.bed", "b.bed")
print(result["piece_count"])
print([(row["chrom"], row["start"], row["end"]) for row in result["pieces"]][:3])
```

- **Example output:**

```text
7
[('chr1', 100, 180), ('chr1', 150, 180), ('chr1', 220, 300)]
```

### 6) `dnakit.annotation.read_gene_annotation()`

- **Purpose:** Read a GTF or GFF3 gene annotation file with one record per feature.
- **Method:** Call polars-bio `read_gtf` or `read_gff` according to `file_format`, with `attributes` selecting attribute fields.
- **API:** `read_gene_annotation(path, *, file_format="gtf", attributes=None, max_rows=100000, verbose=True)`
- **Input:** A GTF or GFF3 file; `file_format` is `gtf` or `gff`.
- **Output:** `features` with `chrom`, `start`, `end`, `type`, `strand`, and attribute columns, plus `feature_count`.
- **Example code:**

```python
from dnakit.annotation import read_gene_annotation

result = read_gene_annotation("genes.gtf", attributes=["gene_id", "gene_name"])
print(result["feature_count"])
print([(row["type"], row["gene_id"], row["start"], row["end"]) for row in result["features"]])
gff = read_gene_annotation("genes.gff3", file_format="gff")
print(gff["feature_count"], sorted({row["type"] for row in gff["features"]}))
```

- **Example output:**

```text
5
[('exon', 'G1', 1000, 1200), ('gene', 'G1', 1000, 2000), ('exon', 'G1', 1800, 2000), ('exon', 'G2', 5000, 6000), ('gene', 'G2', 5000, 6000)]
4 ['exon', 'gene', 'mRNA']
```

### 7) `dnakit.annotation.read_variant_calls()`

- **Purpose:** Read a VCF or BCF variant file with one record per variant, including INFO and FORMAT fields.
- **Method:** Call polars-bio `read_vcf` or `read_bcf` according to `file_format`; coordinates are 0-based half-open.
- **API:** `read_variant_calls(path, *, file_format="vcf", info_fields=None, format_fields=None, max_rows=100000, verbose=True)`
- **Input:** A VCF or BCF file; `info_fields` and `format_fields` select the fields to read.
- **Output:** `variants` with `chrom`, `start`, `end`, `id`, `ref`, `alt`, INFO fields, and `genotypes`, plus `variant_count`.
- **Example code:**

```python
from dnakit.annotation import read_variant_calls

result = read_variant_calls("variants.vcf", info_fields=["DP"], format_fields=["GT"])
print(result["variant_count"])
print([(row["start"], row["end"], row["id"], row["DP"]) for row in result["variants"]])
```

- **Example output:**

```text
2
[(9, 10, 'rs1', 12), (24, 25, '', 8)]
```

### 8) `dnakit.annotation.write_vcf_variants()`

- **Purpose:** Read a VCF file and write its records to a new VCF file.
- **Method:** Read the source with polars-bio and call `write_vcf`; an existing output file is never overwritten.
- **API:** `write_vcf_variants(source, destination, *, info_fields=None, format_fields=None, verbose=True)`
- **Input:** Source VCF path and destination path; `info_fields` and `format_fields` select the fields to keep.
- **Output:** `records_written`, `destination`, and `status`.
- **Example code:**

```python
import pathlib

from dnakit.annotation import write_vcf_variants

pathlib.Path("copy.vcf").unlink(missing_ok=True)
result = write_vcf_variants("variants.vcf", "copy.vcf", info_fields=["DP"], format_fields=["GT"])
print(result["records_written"], pathlib.Path(result["destination"]).name)
```

- **Example output:**

```text
2 copy.vcf
```

### 9) `dnakit.annotation.read_fasta_sequences()`

- **Purpose:** Read every sequence in a FASTA file with its length.
- **Method:** Call polars-bio `read_fasta`; sequences keep their original case, and rows are sorted by name.
- **API:** `read_fasta_sequences(path, *, max_rows=100000, verbose=True)`
- **Input:** A FASTA file path.
- **Output:** `sequences` with `name`, `description`, `sequence`, and `length`, plus `sequence_count`.
- **Example code:**

```python
from dnakit.annotation import read_fasta_sequences

result = read_fasta_sequences("genome.fa")
print(result["sequence_count"], [(row["name"], row["length"]) for row in result["sequences"]])
```

- **Example output:**

```text
2 [('chr1', 14), ('chr2', 12)]
```

### 10) `dnakit.annotation.read_fastq_reads()`

- **Purpose:** Read every read in a FASTQ file with its quality string.
- **Method:** Call polars-bio `read_fastq`; Phred+33 quality strings are kept as stored.
- **API:** `read_fastq_reads(path, *, max_rows=100000, verbose=True)`
- **Input:** A FASTQ file path.
- **Output:** `reads` with `name`, `sequence`, `quality_scores`, and `length`, plus `read_count`.
- **Example code:**

```python
from dnakit.annotation import read_fastq_reads

result = read_fastq_reads("reads.fq")
print(result["read_count"], [(row["name"], row["sequence"], row["quality_scores"]) for row in result["reads"]])
```

- **Example output:**

```text
2 [('r1', 'ACGTN', 'IIII#'), ('r2', 'GGCCAA', 'IIIIII')]
```

### 11) `dnakit.annotation.write_fasta_sequences()`

- **Purpose:** Write DNA sequences, given as `DNA` records or an ID-to-sequence mapping, to a FASTA file.
- **Method:** Sequences are checked against the DNA alphabet, uppercased, and written with polars-bio `write_fasta`; an existing output file is never overwritten.
- **API:** `write_fasta_sequences(sequences, destination, *, verbose=True)`
- **Input:** A sequence mapping or `DNA` object, and the output path.
- **Output:** `records_written`, `destination`, and `status`.
- **Example code:**

```python
import pathlib

from dnakit.annotation import write_fasta_sequences

pathlib.Path("written.fa").unlink(missing_ok=True)
result = write_fasta_sequences({"s1": "acgtn", "s2": "GGCC"}, "written.fa")
print(result["records_written"])
print(pathlib.Path("written.fa").read_text(encoding="utf-8"))
```

- **Example output:**

```text
2
>s1
ACGTN
>s2
GGCC
```

### 12) `dnakit.annotation.read_alignment_records()`

- **Purpose:** Read SAM or BAM alignments with one record per read, including the CIGAR string and sequence.
- **Method:** Call polars-bio `read_sam` or `read_bam` according to `file_format`; coordinates are 0-based half-open.
- **API:** `read_alignment_records(path, *, file_format="sam", tag_fields=None, max_rows=100000, verbose=True)`
- **Input:** A SAM or BAM file; `tag_fields` selects SAM tags to read.
- **Output:** `alignments` with `name`, `chrom`, `start`, `end`, `cigar`, `sequence`, and related fields, plus `record_count`.
- **Example code:**

```python
from dnakit.annotation import read_alignment_records

result = read_alignment_records("aln.sam")
print(result["record_count"])
print([(row["name"], row["start"], row["end"], row["cigar"]) for row in result["alignments"]])
```

- **Example output:**

```text
2
[('r1', 10, 15, '5M'), ('r2', 12, 16, '4M')]
```

### 13) `dnakit.annotation.alignment_depth()`

- **Purpose:** Return the runs of constant read depth in an alignment file.
- **Method:** Call polars-bio `depth` with `filter_flag` and `min_mapping_quality` (the default flag drops unmapped, duplicate, and QC-failed reads).
- **API:** `alignment_depth(path, *, min_mapping_quality=0, filter_flag=1796, max_rows=100000, verbose=True)`
- **Input:** A SAM or BAM file and the filter settings.
- **Output:** `runs` with `chrom`, `start`, `end`, and `depth`, plus `run_count`.
- **Example code:**

```python
from dnakit.annotation import alignment_depth

result = alignment_depth("aln.sam")
print([(row["start"], row["end"], row["depth"]) for row in result["runs"]])
```

- **Example output:**

```text
[(10, 12, 1), (12, 15, 2), (15, 16, 1)]
```

### 14) `dnakit.annotation.fastqc_summary()`

- **Purpose:** Scan a FASTQ file once and report FastQC-style basic statistics and module statuses.
- **Method:** Call polars-bio `fastqc`; GC percentage excludes N bases, and the per-sequence GC distribution is returned.
- **API:** `fastqc_summary(path, *, modules=None, verbose=True)`
- **Input:** A FASTQ file; `modules` optionally limits the modules computed.
- **Output:** `basic_statistics`, `module_status`, and `per_sequence_gc`.
- **Example code:**

```python
from dnakit.annotation import fastqc_summary

result = fastqc_summary("reads.fq")
print(result["basic_statistics"])
print(result["module_status"]["basic_stats"])
```

- **Example output:**

```text
{'gc_pct': 60.0, 'max_len': 6, 'min_len': 5, 'n_seq': 2, 'total_bases': 11}
PASS
```

### 15) `dnakit.annotation.read_contact_pairs()`

- **Purpose:** Read Hi-C read pairs in 4DN Pairs format.
- **Method:** Call polars-bio `read_pairs`; positions become 0-based and column names follow the file header.
- **API:** `read_contact_pairs(path, *, max_rows=100000, verbose=True)`
- **Input:** A Pairs file path.
- **Output:** `pairs` and `pair_count`.
- **Example code:**

```python
from dnakit.annotation import read_contact_pairs

result = read_contact_pairs("contacts.pairs")
print([(row["readID"], row["pos1"], row["pos2"]) for row in result["pairs"]])
```

- **Example output:**

```text
[('r1', 99, 499), ('r2', 199, 39)]
```

### 16) `dnakit.annotation.read_bed_intervals()`

- **Purpose:** Read a BED file with one record per interval, keeping the name and score columns.
- **Method:** Call polars-bio `read_bed` with 0-based half-open coordinates; rows are sorted by coordinate.
- **API:** `read_bed_intervals(path, *, max_rows=100000, verbose=True)`
- **Input:** A BED file path.
- **Output:** `intervals` and `interval_count`.
- **Example code:**

```python
from dnakit.annotation import read_bed_intervals

result = read_bed_intervals("a.bed")
print(result["interval_count"], [(row["name"], row["start"], row["end"]) for row in result["intervals"]][:2])
```

- **Example output:**

```text
4 [('a', 100, 200), ('b', 150, 300)]
```

### 17) `dnakit.annotation.read_cooler_matrix()`

- **Purpose:** Read a Cooler (.cool or .mcool) Hi-C contact matrix with one record per non-zero bin pair.
- **Method:** Call polars-bio `read_cool`; `.mcool` files need `resolution` to pick a resolution.
- **API:** `read_cooler_matrix(path, *, resolution=None, include_weights=False, max_rows=100000, verbose=True)`
- **Input:** A Cooler file; `resolution` and `include_weights` are optional.
- **Output:** `contacts` with `chrom1`, `start1`, `chrom2`, `start2`, and `count`, plus `contact_count` and `total_counts`.
- **Example code:**

```python
import pathlib

import cooler
import pandas as pd

from dnakit.annotation import read_cooler_matrix

bins = pd.DataFrame({"chrom": ["chr1", "chr1"], "start": [0, 100], "end": [100, 200]})
pixels = pd.DataFrame({"bin1_id": [0, 0, 1], "bin2_id": [0, 1, 1], "count": [3, 2, 5]})
pathlib.Path("contacts.cool").unlink(missing_ok=True)
cooler.create_cooler("contacts.cool", bins, pixels, ordered=False)
result = read_cooler_matrix("contacts.cool")
print(result["contact_count"], result["total_counts"])
print([(row["start1"], row["start2"], row["count"]) for row in result["contacts"]])
```

- **Example output:**

```text
3 10
[(0, 0, 3), (0, 100, 2), (100, 100, 5)]
```

### 18) `dnakit.annotation.read_stockholm_alignments()`

- **Purpose:** Read a Stockholm multiple-sequence alignment with one record per aligned sequence.
- **Method:** Call polars-bio `read_sto`; gap characters are kept in the aligned sequences.
- **API:** `read_stockholm_alignments(path, *, max_rows=100000, verbose=True)`
- **Input:** A Stockholm file (`.sto` or `.stk`).
- **Output:** `alignments` with `alignment_id`, `name`, and `sequence`, plus `sequence_count` and `alignment_columns`.
- **Example code:**

```python
from dnakit.annotation import read_stockholm_alignments

result = read_stockholm_alignments("align.sto")
print(result["sequence_count"], result["alignment_columns"])
print([(row["name"], row["sequence"]) for row in result["alignments"]])
```

- **Example output:**

```text
2 9
[('seq1', 'ACGT-ACGT'), ('seq2', 'ACGTTACGA')]
```

## gtars: genome sequences, region sets, partitions, LOLA, tokens, and VRS

### 19) `dnakit.annotation.region_gc_content()`

- **Purpose:** Compute the GC fraction of each BED region in a reference genome.
- **Method:** Read the FASTA with gtars `GenomeAssembly` and call `calc_gc_content`; N bases are not counted as GC.
- **API:** `region_gc_content(regions, genome_fasta, *, ignore_unknown_chromosomes=False, max_rows=100000, verbose=True)`
- **Input:** A BED region file and a FASTA genome; `ignore_unknown_chromosomes` decides whether regions on chromosomes missing from the FASTA are skipped.
- **Output:** `regions` with `chrom`, `start`, `end`, and `gc_fraction`, plus `mean_gc_fraction`.
- **Example code:**

```python
from dnakit.annotation import region_gc_content

result = region_gc_content("gc.bed", "genome.fa")
print(result["mean_gc_fraction"])
print([(row["chrom"], row["start"], row["end"], row["gc_fraction"]) for row in result["regions"]])
```

- **Example output:**

```text
0.5
[('chr1', 0, 10, 0.5), ('chr1', 10, 14, 0.5), ('chr2', 0, 4, 0.0), ('chr2', 8, 12, 1.0)]
```

### 20) `dnakit.annotation.region_dinucleotide_frequencies()`

- **Purpose:** Count, or give as percentages, the 16 dinucleotides inside each BED region.
- **Method:** Call gtars `calc_dinucl_freq`; dinucleotide names are uppercased and pooled totals are returned.
- **API:** `region_dinucleotide_frequencies(regions, genome_fasta, *, raw_counts=False, ignore_unknown_chromosomes=False, max_rows=100000, verbose=True)`
- **Input:** A BED region file and a FASTA genome; `raw_counts=True` returns counts, otherwise percentages.
- **Output:** `regions` with per-region `frequencies`, `pooled` totals, `dinucleotides`, and `value_unit`.
- **Example code:**

```python
from dnakit.annotation import region_dinucleotide_frequencies

result = region_dinucleotide_frequencies("gc.bed", "genome.fa", raw_counts=True)
print(result["value_unit"], result["regions"][0]["label"])
print(result["regions"][0]["frequencies"]["AC"], result["regions"][0]["frequencies"]["CG"])
```

- **Example output:**

```text
count chr1_0_10
3.0 2.0
```

### 21) `dnakit.annotation.partition_region_counts()`

- **Purpose:** Classify BED regions as core promoter, proximal promoter, exon, intron, or intergenic and count them.
- **Method:** Build the partitions with gtars `PartitionList.from_gtf` and call `calc_partitions`; with `chrom_sizes`, call `calc_expected_partitions` and return chi-square p-values.
- **API:** `partition_region_counts(regions, gene_annotation, *, core_promoter_bp=100, proximal_promoter_bp=2000, protein_coding_only=True, chrom_sizes=None, bp_proportion=False, verbose=True)`
- **Input:** A BED region file and a GTF annotation; promoter lengths and protein-coding filtering can be set.
- **Output:** `partitions` with observed counts (and expected counts and p-values when `chrom_sizes` is given), plus `total`.
- **Example code:**

```python
from dnakit.annotation import partition_region_counts

observed = partition_region_counts("regions.bed", "genes.gtf")
print(observed["total"], [(row["partition"], row["observed"]) for row in observed["partitions"]])
expected = partition_region_counts("regions.bed", "genes.gtf", chrom_sizes={"chr1": 10000})
print(expected["mode"], [round(row["expected"], 2) for row in expected["partitions"]])
```

- **Example output:**

```text
4 [('promoterCore', 2), ('promoterProx', 0), ('exon', 0), ('intron', 1), ('intergenic', 1)]
expected [0.08, 1.12, 0.56, 0.24, 2.0]
```

### 22) `dnakit.annotation.signal_matrix_summary()`

- **Purpose:** Take the signal of labeled regions that overlap the query regions and summarize each condition.
- **Method:** Read the signal table with gtars `SignalMatrix.from_tsv` and call `calc_summary_signal`; the first column holds `chrom_start_end` labels and the other columns hold one condition each.
- **API:** `signal_matrix_summary(regions, signal_table, *, verbose=True)`
- **Input:** A BED region file and a signal table (TSV).
- **Output:** `regions` with per-region `values`, `condition_names`, and `matrix_stats`.
- **Example code:**

```python
from dnakit.annotation import signal_matrix_summary

result = signal_matrix_summary("gc.bed", "signal.tsv")
print(result["condition_names"])
print([(row["label"], row["values"]) for row in result["regions"]])
```

- **Example output:**

```text
['condA', 'condB']
[('chr1_0_10', {'condA': 1.0, 'condB': 2.0}), ('chr1_10_14', {'condA': 1.0, 'condB': 2.0}), ('chr2_0_4', {'condA': 3.0, 'condB': 4.0}), ('chr2_8_12', {'condA': 3.0, 'condB': 4.0})]
```

### 23) `dnakit.annotation.consensus_regions()`

- **Purpose:** Count how many of the input BED files contain each distinct region.
- **Method:** Read each file as a gtars `RegionSet`, call `consensus`, and sort the result by coordinate.
- **API:** `consensus_regions(region_files, *, max_rows=100000, verbose=True)`
- **Input:** At least two BED file paths.
- **Output:** `regions` with `chrom`, `start`, `end`, and `count`, plus `file_count` and `region_count`.
- **Example code:**

```python
from dnakit.annotation import consensus_regions

result = consensus_regions(["set1.bed", "set2.bed"])
print([(row["chrom"], row["start"], row["end"], row["count"]) for row in result["regions"]])
```

- **Example output:**

```text
[('chr1', 100, 200, 2), ('chr1', 300, 400, 1), ('chr1', 500, 600, 1)]
```

### 24) `dnakit.annotation.lola_enrichment()`

- **Purpose:** Test whether the user regions overlap each database region set more than chance (LOLA-style enrichment).
- **Method:** Build the database with gtars `RegionDB.from_bed_files` and call `run_lola`; user and database regions must lie inside the given universe.
- **API:** `lola_enrichment(user_regions, universe, database_files, *, min_overlap=1, direction="enrichment", max_rows=100000, verbose=True)`
- **Input:** A user BED file, a universe BED3 file, and a list of database BED files whose names are unique.
- **Output:** `results` with `p_value_log10`, `odds_ratio`, `support`, and contingency counts b, c, d for each database set.
- **Example code:**

```python
from dnakit.annotation import lola_enrichment

result = lola_enrichment("user.bed", "universe.bed", ["db1.bed", "db2.bed"])
print([(row["db_name"], row["support"], row["b"], row["c"], row["d"]) for row in result["results"]])
```

- **Example output:**

```text
[('db1.bed', 1, 0, 1, 1), ('db2.bed', 1, 1, 1, 0)]
```

### 25) `dnakit.annotation.tokenize_fragments()`

- **Purpose:** Map fragment-file records onto universe regions (BED3) and give each barcode a list of region tokens.
- **Method:** Build a gtars `Tokenizer.from_bed` and call `tokenize_fragment_file`; a token id is the row number of the universe region.
- **API:** `tokenize_fragments(fragments, universe, *, max_rows=100000, verbose=True)`
- **Input:** A fragment file with chrom, start, end, and barcode columns, and a BED3 universe file.
- **Output:** `barcodes` with `barcode`, `token_ids`, and `token_count`, plus `universe_size`.
- **Example code:**

```python
from dnakit.annotation import tokenize_fragments

result = tokenize_fragments("frags.tsv", "universe.bed")
print(result["universe_size"], [(row["barcode"], row["token_ids"]) for row in result["barcodes"]])
```

- **Example output:**

```text
3 [('AAA', [0, 1]), ('BBB', [2])]
```

### 26) `dnakit.annotation.digest_fasta_sequences()`

- **Purpose:** Compute GA4GH sha512t24u and MD5 digests for every sequence in a FASTA file, plus a collection digest.
- **Method:** Call gtars `refget.digest_fasta`; digests follow the GA4GH refget specification.
- **API:** `digest_fasta_sequences(fasta, *, max_rows=100000, verbose=True)`
- **Input:** An uncompressed FASTA file path.
- **Output:** `sequences` with `name`, `length`, `sha512t24u`, and `md5`, plus `collection_digest` and `sequence_count`.
- **Example code:**

```python
from dnakit.annotation import digest_fasta_sequences

result = digest_fasta_sequences("genome.fa")
print(result["collection_digest"])
print([(row["name"], row["length"], row["sha512t24u"]) for row in result["sequences"]])
```

- **Example output:**

```text
rM-CtSsxdCWLtxc1sQpCDYBRFxBRXGP2
[('chr1', 14, 'bxzNj2Cq11mJJ1GuetSIcq2lO4sKlIds'), ('chr2', 12, 'yPlv_G-_BuLFIFYisu1dBipO_I3_bLNz')]
```

### 27) `dnakit.annotation.fasta_index_records()`

- **Purpose:** Return the index records (sequence name and length) of a FASTA file, matching samtools faidx names and lengths.
- **Method:** Call gtars `refget.compute_fai`; only uncompressed FASTA files are supported.
- **API:** `fasta_index_records(fasta, *, max_rows=100000, verbose=True)`
- **Input:** An uncompressed FASTA file path.
- **Output:** `records` with `name` and `length`, plus `sequence_count`.
- **Example code:**

```python
from dnakit.annotation import fasta_index_records

result = fasta_index_records("genome.fa")
print([(row["name"], row["length"]) for row in result["records"]])
```

- **Example output:**

```text
[('chr1', 14), ('chr2', 12)]
```

### 28) `dnakit.annotation.sequence_digest()`

- **Purpose:** Compute the GA4GH sha512t24u or MD5 digest of one DNA sequence.
- **Method:** The sequence is checked against the DNA alphabet and uppercased, then passed to gtars `sha512t24u_digest` or `md5_digest`.
- **API:** `sequence_digest(sequence, *, algorithm="sha512t24u", verbose=True)`
- **Input:** A DNA sequence string; `algorithm` is `sha512t24u` or `md5`.
- **Output:** `digest`, `algorithm`, and `length`.
- **Example code:**

```python
from dnakit.annotation import sequence_digest

print(sequence_digest("acgt")["digest"])
print(sequence_digest("ACGT", algorithm="md5")["digest"])
```

- **Example output:**

```text
aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2
f1f8f4bf413b16ad135722aa4591043e
```

### 29) `dnakit.annotation.normalize_variant_allele()`

- **Purpose:** Trim bases shared by both alleles and expand the variant over adjacent repeats to get a canonical form.
- **Method:** Check that `ref_allele` matches the reference at `start`, then call gtars `vrs.normalize_allele`; for example, deleting one T from TTTT is reported over the whole run as TTT.
- **API:** `normalize_variant_allele(sequence, start, ref_allele, alt_allele, *, verbose=True)`
- **Input:** The reference sequence, a 0-based start, the reference allele, and the alternate allele.
- **Output:** `start` and `end` (0-based, `end` exclusive) and `allele`.
- **Example code:**

```python
from dnakit.annotation import normalize_variant_allele

snv = normalize_variant_allele("ACGTACGTAC", 2, "G", "T")
print(snv["start"], snv["end"], snv["allele"])
deletion = normalize_variant_allele("ACGTTTTA", 2, "GT", "G")
print(deletion["start"], deletion["end"], repr(deletion["allele"]))
```

- **Example output:**

```text
2 3 T
3 7 'TTT'
```

### 30) `dnakit.annotation.vrs_allele_id()`

- **Purpose:** Compute the GA4GH VRS identifier (`ga4gh:VA.<digest>`) of one allele on a sequence.
- **Method:** Call gtars `vrs.vrs_id`; with `prefixed=False` the bare `vrs.vrs_digest` is returned instead.
- **API:** `vrs_allele_id(sequence_digest_value, start, end, alt, *, prefixed=True, verbose=True)`
- **Input:** The sequence's sha512t24u digest, a 0-based `start`, an exclusive `end`, and the alternate allele `alt`.
- **Output:** `identifier`, `prefixed`, and the digest and coordinates used.
- **Example code:**

```python
from dnakit.annotation import sequence_digest, vrs_allele_id

digest = sequence_digest("ACGTACGTACGGTT")["digest"]
print(vrs_allele_id(digest, 1, 2, "G")["identifier"])
print(vrs_allele_id(digest, 1, 2, "G", prefixed=False)["identifier"])
```

- **Example output:**

```text
ga4gh:VA.BwLTTHJnaEI7FfDFxp8JjUw08nmooAjb
BwLTTHJnaEI7FfDFxp8JjUw08nmooAjb
```

### 31) `dnakit.annotation.parse_hgvs_variant()`

- **Purpose:** Parse one HGVS expression such as `NC_000001.11:g.12345A>G` into its accession, reference type, and edits.
- **Method:** Call gtars `vrs.hgvs.parse_hgvs` and convert the result to a JSON-compatible dictionary.
- **API:** `parse_hgvs_variant(text, *, verbose=True)`
- **Input:** An HGVS string.
- **Output:** `accession`, `reference_type`, and the full `variant` dictionary.
- **Example code:**

```python
from dnakit.annotation import parse_hgvs_variant

result = parse_hgvs_variant("NC_000001.11:g.12345A>G")
print(result["accession"], result["reference_type"])
print(result["variant"]["pos_edit"]["edit"])
```

- **Example output:**

```text
NC_000001.11 g
{'kind': 'substitution', 'ref': 'A', 'alt': 'G'}
```

## genomicranges: GTF ranges, promoter windows, and binning

### 32) `dnakit.annotation.read_gtf_ranges()`

- **Purpose:** Read a GTF file with genomicranges, with one record per feature and every attribute column kept.
- **Method:** Call genomicranges `read_gtf`, convert its 1-based starts to 0-based, and restore the GTF end, which upstream stores minus one.
- **API:** `read_gtf_ranges(path, *, max_rows=100000, verbose=True)`
- **Input:** A GTF file path.
- **Output:** `features` with `chrom`, `start`, `end`, `strand`, and attribute columns, plus `feature_count`.
- **Example code:**

```python
from dnakit.annotation import read_gtf_ranges

result = read_gtf_ranges("genes.gtf")
print(result["feature_count"])
print([(row["feature"], row["gene_id"], row["start"], row["end"], row["strand"]) for row in result["features"]])
```

- **Example output:**

```text
5
[('exon', 'G1', 1000, 1200, '+'), ('gene', 'G1', 1000, 2000, '+'), ('exon', 'G1', 1800, 2000, '+'), ('exon', 'G2', 5000, 6000, '-'), ('gene', 'G2', 5000, 6000, '-')]
```

### 33) `dnakit.annotation.upstream_downstream_regions()`

- **Purpose:** Build promoter (TSS) or terminator (TES) windows around each BED interval.
- **Method:** Call genomicranges `compute_up_down` with 1-based inclusive ranges and integer strand codes (-1 for minus). Upstream returns widths, not ends, so the wrapper converts them to window ends and clips starts at 0.
- **API:** `upstream_downstream_regions(source, *, upstream=1000, downstream=0, site="TSS", max_rows=100000, verbose=True)`
- **Input:** A BED file; `upstream` and `downstream` are base counts; `site` is `TSS` or `TES`.
- **Output:** `windows` with `chrom`, `start`, `end`, `strand`, and the source coordinates, plus `window_count` and `site`.
- **Example code:**

```python
from dnakit.annotation import upstream_downstream_regions

tss = upstream_downstream_regions("regions.bed", upstream=100, downstream=50)
print([(row["source_start"], row["start"], row["end"]) for row in tss["windows"]])
tes = upstream_downstream_regions("regions.bed", upstream=100, downstream=50, site="TES")
print([(row["source_start"], row["start"], row["end"]) for row in tes["windows"]][:2])
```

- **Example output:**

```text
[(950, 850, 1000), (1500, 1400, 1550), (5950, 5850, 6000), (9000, 8900, 9050)]
[(950, 909, 1059), (1500, 1499, 1649)]
```

### 34) `dnakit.annotation.split_interval_bins()`

- **Purpose:** Tile the half-open range [start, end) into consecutive bins of `bin_size` bases; the last bin may be shorter.
- **Method:** Call genomicranges.utils `split_intervals` with the inclusive end (end-1), convert its offsets to half-open bins, and check that the bins cover the range exactly.
- **API:** `split_interval_bins(start, end, bin_size, *, verbose=True)`
- **Input:** Integers `start`, `end` (exclusive), and `bin_size`.
- **Output:** `bins` with `index`, `start`, and `end`, plus `bin_count`.
- **Example code:**

```python
from dnakit.annotation import split_interval_bins

result = split_interval_bins(0, 100, 30)
print([(row["start"], row["end"]) for row in result["bins"]])
```

- **Example output:**

```text
[(0, 30), (30, 60), (60, 90), (90, 100)]
```

### 35) `dnakit.annotation.slide_interval_windows()`

- **Purpose:** Slide a window of `width` bases by `step` across the half-open range [start, end).
- **Method:** Call genomicranges.utils `slide_intervals` and convert its inclusive results to half-open windows; a range shorter than `width` gives one window.
- **API:** `slide_interval_windows(start, end, width, step, *, verbose=True)`
- **Input:** Integers `start`, `end` (exclusive), `width`, and `step`.
- **Output:** `windows` with `index`, `start`, and `end`, plus `window_count`.
- **Example code:**

```python
from dnakit.annotation import slide_interval_windows

result = slide_interval_windows(0, 100, 30, 10)
print(result["window_count"], [(row["start"], row["end"]) for row in result["windows"]][:3])
```

- **Example output:**

```text
8 [(0, 30), (10, 40), (20, 50)]
```

## Coverage and skips

The table summarizes the status of every public function in each package inventory; the batch ledger records the reason for each one.

| Package | Inventory functions | Status counts |
| --- | --- | --- |
| polars_bio | 211 | duplicate-of 3, family-variant 17, internal 94, not-dna 10, other 68, plot 1, wrapped 18 |
| ncls | 5 | already-wrapped 1, internal 3, not-dna 1 |
| genomicranges | 22 | duplicate-of 3, family-variant 1, internal 8, needs-network 4, not-dna 2, wrapped 4 |
| gtars | 29 | family-variant 5, other 11, wrapped 13 |
