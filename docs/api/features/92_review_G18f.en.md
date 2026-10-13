# G18f review functions: polars-bio and gtars

This page documents the G18f wrappers. The eight polars-bio functions and six gtars functions are all exported from `dnakit.annotation`, for example `from dnakit.annotation import lola_check_universe`.

- Local functions: all 14 wrappers read or write local files only. They need no network, and the tests use no network either.
- Example input files are declared with the `Example file` label. Binary inputs (CRAM, BAM, bigBed, GTOK) are made in the example code, so each example runs as written.
- mygene and myvariant have no DNA-level function in this group. Their modules export nothing; see the unwrapped section below.

## Existing and unwrapped functions

- polars-bio: 68 `other` rows. 8 are wrapped (the functions on this page). 35 are family variants of those functions, such as the lazy `scan_*`, `sink_*`, and `write_bam` / `write_cram` / `write_sam` forms, and the 16 `register_*` forms of `query_genomic_tables`. 6 are session or global-state helpers (`get_option`, `set_option`, `set_loglevel`, `set_source_metadata`, `register_view`, `from_polars`). 15 need a fixture writer that is not installed: BGEN, PGEN, and VCF Zarr readers and describers (`needs-data:bgen-fixture`, `needs-data:pgen-fixture`, `needs-data:vcf-zarr-store`), plus `register_cram`, which takes no reference argument (`needs-data:reference-fasta`). 2 are generic table readers (`read_table`, `scan_table`; `not-dna`). `read_bigwig` and `scan_bigwig` duplicate the existing `read_pybigwig_intervals` (`duplicate-of`).
- gtars: 11 `other` rows. 6 are wrapped. `read_tokens_from_gtok_as_strings` is a family variant of `read_region_tokens`. `hgvs_to_vrs_id` needs a transcript store (`needs-data:cdot-transcript-store`). `median_abs_distance` is a generic median (`not-dna`). `force_unlock_store` and `store_lock_status` administer RefgetStore locks (`internal`).
- mygene: 10 `other` rows, all client classes or factories, cache and endpoint settings (`internal`), or metadata listings (`get_fields`, `metadata`; `not-dna`). No function returns gene records that are not already wrapped.
- myvariant: 10 `other` rows with the same breakdown as mygene. No function returns variant records that are not already wrapped.

## 1) `read_file_header`

- **Purpose:** Return the header of a VCF, BCF, SAM, or BAM file as polars-bio parses it, without reading records.
- **Method:** Opens the file with the matching lazy `scan_*` reader and reads `get_metadata` from the frame. The internal `_datafusion_table_name` field is removed.
- **API:** `read_file_header(path[required], file_format[optional], verbose[optional])`
- **Input:** A local VCF, BCF, SAM, or BAM file. `file_format` is `vcf` (default), `bcf`, `sam`, or `bam`.
- **Output:** A dictionary with `file_format`, `coordinate_system_zero_based`, `header` (format-specific fields such as `sample_names`, `contigs`, or `sort_order`), and `provenance`.

**Example file:** `sample.vcf`

```text
##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s1	s2
chr1	10	rs1	A	G	50	PASS	DP=12	GT	0/1	1/1
chr1	25	.	C	T	40	PASS	DP=8	GT	0/0	0/1
```

**Example code:**

```python
from dnakit.annotation import read_file_header

result = read_file_header("sample.vcf", file_format="vcf", verbose=False)
print(result["header"]["sample_names"])
print(result["header"]["contigs"])
```

**Example output:**

```text
['s1', 's2']
[{'id': 'chr1', 'length': 1000}]
```

## 2) `describe_genomic_file`

- **Purpose:** Describe the schema or header entries of a VCF, BCF, SAM, BAM, CRAM, cool, or Stockholm (STO) file.
- **Method:** Calls the polars-bio `describe_*` function for the format. For SAM, BAM, and CRAM, tag types come from a sample of records. CRAM needs `reference_path`, and `reference_path` is rejected for every other format.
- **API:** `describe_genomic_file(path[required], file_format[optional], reference_path[optional], verbose[optional])`
- **Input:** A local file. `file_format` is one of `vcf`, `bcf`, `sam`, `bam`, `cram`, `cool`, or `sto`.
- **Output:** A dictionary with `file_format`, `entry_count`, `entries` (one dictionary per entry; the keys depend on the format), and `provenance`.

**Example file:** `sample.vcf`

```text
##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s1	s2
chr1	10	rs1	A	G	50	PASS	DP=12	GT	0/1	1/1
chr1	25	.	C	T	40	PASS	DP=8	GT	0/0	0/1
```

**Example code:**

```python
from dnakit.annotation import describe_genomic_file

result = describe_genomic_file("sample.vcf", file_format="vcf", verbose=False)
for entry in result["entries"]:
    print(entry["name"], entry["field_type"], entry["data_type"])
```

**Example output:**

```text
DP INFO Integer
genotypes FORMAT Struct
```

## 3) `read_cram_records`

- **Purpose:** Read CRAM alignments, one record per read, decoded against a reference FASTA.
- **Method:** Calls polars-bio `read_cram` with the reference and reports coordinates as 0-based half-open intervals. The reference needs a `.fai` index, which `pysam.faidx` creates.
- **API:** `read_cram_records(path[required], reference_path[required], tag_fields[optional], max_rows[optional], verbose[optional])`
- **Input:** A CRAM file and its reference FASTA. The example writes the CRAM from the SAM file with pysam, since CRAM is binary.
- **Output:** A dictionary with `record_count`, `alignments` (name, chrom, start, end, flags, cigar, and the sequence fields), and `provenance`. The example reads SAM POS 5 as start 4 and 40 with `5M2I3M` as end 47.

**Example file:** `reads.sam`

```text
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:73
@RG	ID:rg1	SM:s1
r1	0	chr1	5	60	10M	*	0	0	ACGTACGTTA	IIIIIIIIII	RG:Z:rg1
r2	16	chr1	20	40	8M	*	0	0	GCCGATTA	IIIIIIII	RG:Z:rg1
r3	0	chr1	40	30	5M2I3M	*	0	0	ACGTAGCTAG	IIIIIIIIII	NM:i:2
```

**Example file:** `ref.fa`

```text
>chr1
ACGTACGTTAGCCGATTAGCATGCAATGCCGTAGCTAGCTAGGATCGATCGTAGCTAGCTAGCTACGATCGAT
```

**Example code:**

```python
import pathlib

import pysam
from dnakit.annotation import read_cram_records

pysam.faidx("ref.fa")
pathlib.Path("reads.cram").unlink(missing_ok=True)
pysam.view("-C", "-T", "ref.fa", "-o", "reads.cram", "reads.sam", catch_stdout=False)
result = read_cram_records("reads.cram", "ref.fa", verbose=False)
for row in result["alignments"]:
    print(row["name"], row["start"], row["end"], row["cigar"])
```

**Example output:**

```text
r1 4 14 10M
r2 19 27 8M
r3 39 47 5M2I3M
```

## 4) `read_bigbed_intervals`

- **Purpose:** Read every interval of a bigBed file together with its stored fields.
- **Method:** Calls polars-bio `read_bigbed` with `schema='rest'` and 0-based coordinates. The tab-separated text after `end` is split into `fields`.
- **API:** `read_bigbed_intervals(path[required], max_rows[optional], verbose[optional])`
- **Input:** A local bigBed file. The example writes one with pybigtools, because bigBed is binary.
- **Output:** A dictionary with `interval_count`, `intervals` (chrom, start, end, and `fields`), and `provenance`.

**Example code:**

```python
import pathlib

import pybigtools
from dnakit.annotation import read_bigbed_intervals

pathlib.Path("features.bb").unlink(missing_ok=True)
writer = pybigtools.open("features.bb", "w")
writer.write(
    {"chr1": 100, "chr2": 100},
    [
        ("chr1", 0, 10, "featA\t500\t+"),
        ("chr1", 20, 35, "featB\t300\t-"),
        ("chr2", 5, 15, "featC\t100\t+"),
    ],
)
result = read_bigbed_intervals("features.bb", verbose=False)
for row in result["intervals"]:
    print(row["chrom"], row["start"], row["end"], row["fields"][0])
```

**Example output:**

```text
chr1 0 10 featA
chr1 20 35 featB
chr2 5 15 featC
```

## 5) `read_alignment_sequences`

- **Purpose:** Read an A2M or A3M multiple-sequence alignment, one record per aligned sequence.
- **Method:** Calls polars-bio `read_a2m` or `read_a3m`. Gap characters and A3M lowercase insertions are kept as stored. `alignment_columns` is set only when every sequence has the same length.
- **API:** `read_alignment_sequences(path[required], file_format[optional], max_rows[optional], verbose[optional])`
- **Input:** A local A2M or A3M file with `>` headers. `file_format` is `a2m` (default) or `a3m`.
- **Output:** A dictionary with `file_format`, `sequence_count`, `alignment_columns`, `sequences` (name, description, sequence), and `provenance`.

**Example file:** `align.a2m`

```text
>seq1
ACGT-ACGT
>seq2
AC-TACGTT
```

**Example code:**

```python
from dnakit.annotation import read_alignment_sequences

result = read_alignment_sequences("align.a2m", file_format="a2m", verbose=False)
print(result["sequence_count"], result["alignment_columns"])
for row in result["sequences"]:
    print(row["name"], row["sequence"])
```

**Example output:**

```text
2 9
seq1 ACGT-ACGT
seq2 AC-TACGTT
```

## 6) `write_fastq_reads`

- **Purpose:** Write reads with names, DNA sequences, and quality strings to a new FASTQ file.
- **Method:** Checks each name (no whitespace), each sequence (DNA letters, uppercased), and each `quality_scores` string (Phred+33 characters, one per base) before calling polars-bio `write_fastq`. An existing file is never overwritten.
- **API:** `write_fastq_reads(reads[required], destination[required], verbose[optional])`
- **Input:** A list of dictionaries with `name`, `sequence`, and `quality_scores`, plus an optional `description`. The records returned by `read_fastq_reads` can be passed directly.
- **Output:** A dictionary with `records_written`, `destination`, and `provenance`.

**Example code:**

```python
import pathlib

from dnakit.annotation import write_fastq_reads

reads = [
    {"name": "r1", "description": "first", "sequence": "acgtn", "quality_scores": "IIII#"},
    {"name": "r2", "description": None, "sequence": "GG", "quality_scores": "II"},
]
pathlib.Path("out.fastq").unlink(missing_ok=True)
result = write_fastq_reads(reads, "out.fastq", verbose=False)
print(result["records_written"])
print(pathlib.Path("out.fastq").read_text(encoding="utf-8"))
```

**Example output:**

```text
2
@r1 first
ACGTN
+
IIII#
@r2
GG
+
II
```

## 7) `write_alignment_file`

- **Purpose:** Convert alignments between SAM, BAM, and CRAM and keep the record order.
- **Method:** Reads the source with polars-bio `read_sam`, `read_bam`, or `read_cram`, then writes with `write_sam`, `write_bam`, or `write_cram`. CRAM input or output needs `reference_path`. The writer sets the `@HD SO` value itself.
- **API:** `write_alignment_file(source[required], destination[required], source_format[optional], output_format[optional], reference_path[optional], verbose[optional])`
- **Input:** A local SAM, BAM, or CRAM file. `source_format` and `output_format` are `sam` (default source), `bam` (default output), or `cram`.
- **Output:** A dictionary with `source_format`, `output_format`, `records_written`, `destination`, and `provenance`.

**Example file:** `reads.sam`

```text
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:73
@RG	ID:rg1	SM:s1
r1	0	chr1	5	60	10M	*	0	0	ACGTACGTTA	IIIIIIIIII	RG:Z:rg1
r2	16	chr1	20	40	8M	*	0	0	GCCGATTA	IIIIIIII	RG:Z:rg1
r3	0	chr1	40	30	5M2I3M	*	0	0	ACGTAGCTAG	IIIIIIIIII	NM:i:2
```

**Example code:**

```python
import pathlib

import pysam
from dnakit.annotation import write_alignment_file

pathlib.Path("reads_out.bam").unlink(missing_ok=True)
result = write_alignment_file(
    "reads.sam",
    "reads_out.bam",
    source_format="sam",
    output_format="bam",
    verbose=False,
)
print(result["output_format"], result["records_written"])
with pysam.AlignmentFile("reads_out.bam", "rb", check_sq=False) as handle:
    print([record.reference_start + 1 for record in handle.fetch(until_eof=True)])
```

**Example output:**

```text
bam 3
[5, 20, 40]
```

## 8) `query_genomic_tables`

- **Purpose:** Run one read-only SELECT or WITH query over local genomic files.
- **Method:** Each file is registered under a lowercase table name in polars-bio's session for this call only, then the query runs. Tables are removed in a `finally` block, even when the query fails. Other SQL statements and table names already registered are rejected. Rows keep the query order when the SQL has `ORDER BY`; otherwise they are sorted.
- **API:** `query_genomic_tables(tables[required], query[required], max_rows[optional], verbose[optional])`; `tables` maps names to `(path, file_format)` pairs. Formats: `bed`, `gtf`, `gff`, `vcf`, `bcf`, `sam`, `bam`, `fasta`, `fastq`, `a2m`, `a3m`, `sto`, `pairs`, `cool`, `bigwig`, `bigbed`.
- **Input:** Local files and one SQL statement. The query uses DataFusion SQL; `end` is a reserved word, so it is quoted.
- **Output:** A dictionary with `table_names`, `columns`, `row_count`, `rows`, and `provenance`. BED and VCF starts are 0-based, so the VCF POS 10 appears as 9.

**Example file:** `iv.bed`

```text
chr1	0	10	featA	500	+
chr1	20	35	featB	300	-
chr2	5	15	featC	100	+
```

**Example file:** `sample.vcf`

```text
##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s1	s2
chr1	10	rs1	A	G	50	PASS	DP=12	GT	0/1	1/1
chr1	25	.	C	T	40	PASS	DP=8	GT	0/0	0/1
```

**Example code:**

```python
from dnakit.annotation import query_genomic_tables

result = query_genomic_tables(
    {"peaks": ("iv.bed", "bed"), "vars": ("sample.vcf", "vcf")},
    'SELECT chrom, start, "end", name FROM peaks ORDER BY start',
    verbose=False,
)
for row in result["rows"]:
    print(row["chrom"], row["start"], row["end"], row["name"])
```

**Example output:**

```text
chr1 0 10 featA
chr2 5 15 featC
chr1 20 35 featB
```

## 9) `lola_restricted_universe`

- **Purpose:** Build the restricted LOLA universe from several user region sets.
- **Method:** Calls gtars `lola.build_restricted_universe` on BED3 user sets. Overlapping regions are cut into non-overlapping pieces, so the result can serve as the universe for `lola_enrichment`.
- **API:** `lola_restricted_universe(user_regions[required], max_rows[optional], verbose[optional])`
- **Input:** A list of BED files, one per user set. Each row needs at least chrom, start, and end, with 0 <= start < end.
- **Output:** A dictionary with `user_set_count`, `region_count`, `regions` (chrom, start, end), and `provenance`. In the example, chr1 [0, 20) is cut at 5 and 10.

**Example file:** `setA.bed`

```text
chr1	0	10
chr1	5	20
```

**Example file:** `setB.bed`

```text
chr1	0	10
chr2	1	4
```

**Example code:**

```python
from dnakit.annotation import lola_restricted_universe

result = lola_restricted_universe(["setA.bed", "setB.bed"], verbose=False)
for region in result["regions"]:
    print(region["chrom"], region["start"], region["end"])
```

**Example output:**

```text
chr1 0 5
chr1 5 10
chr1 10 20
chr2 1 4
```

## 10) `lola_check_universe`

- **Purpose:** Check how well a universe BED file covers each user region set, with LOLA's warnings.
- **Method:** Calls gtars `lola.check_universe`. `coverage` is the fraction of user regions that overlap the universe. `many_to_many` counts user regions overlapping more than one universe region, which `lola_redefine_user_sets` removes.
- **API:** `lola_check_universe(user_regions[required], universe[required], max_rows[optional], verbose[optional])`
- **Input:** A list of user-set BED files and one universe BED file.
- **Output:** A dictionary with `user_set_count`, `user_sets` (user_set, total_regions, regions_in_universe, coverage, many_to_many), `warnings`, and `provenance`.

**Example file:** `setA.bed`

```text
chr1	0	10
chr1	5	20
```

**Example file:** `uni.bed`

```text
chr1	0	10
chr1	5	20
chr2	0	50
chr1	100	110
```

**Example code:**

```python
from dnakit.annotation import lola_check_universe

result = lola_check_universe(["setA.bed"], "uni.bed", verbose=False)
for item in result["user_sets"]:
    print(
        item["user_set"],
        item["total_regions"],
        item["regions_in_universe"],
        item["coverage"],
        item["many_to_many"],
    )
print(len(result["warnings"]))
```

**Example output:**

```text
0 2 2 1.0 2
1
```

## 11) `lola_redefine_user_sets`

- **Purpose:** Redefine each user region set as the universe regions it overlaps.
- **Method:** Calls gtars `lola.redefine_user_sets`. Each user set is replaced by the universe regions it touches, sorted by chrom, start, and end.
- **API:** `lola_redefine_user_sets(user_regions[required], universe[required], max_rows[optional], verbose[optional])`
- **Input:** A list of user-set BED files and one universe BED file.
- **Output:** A dictionary with `user_set_count`, `user_sets` (user_set and regions), and `provenance`.

**Example file:** `setA.bed`

```text
chr1	0	10
chr1	5	20
```

**Example file:** `setB.bed`

```text
chr1	0	10
chr2	1	4
```

**Example file:** `uni.bed`

```text
chr1	0	10
chr1	5	20
chr2	0	50
chr1	100	110
```

**Example code:**

```python
from dnakit.annotation import lola_redefine_user_sets

result = lola_redefine_user_sets(["setA.bed", "setB.bed"], "uni.bed", verbose=False)
for item in result["user_sets"]:
    spans = [(region["chrom"], region["start"], region["end"]) for region in item["regions"]]
    print(item["user_set"], spans)
```

**Example output:**

```text
0 [('chr1', 0, 10), ('chr1', 5, 20)]
1 [('chr1', 0, 10), ('chr1', 5, 20), ('chr2', 0, 50)]
```

## 12) `write_region_tokens`

- **Purpose:** Write region token ids to a new GTOK binary file.
- **Method:** Checks that every id is an integer from 0 to 2**32 - 1, then calls gtars `utils.write_tokens_to_gtok`. Ids are stored in the given order.
- **API:** `write_region_tokens(tokens[required], destination[required], verbose[optional])`
- **Input:** A non-empty list of non-negative integers, such as the token ids of a universe from `tokenize_fragments`.
- **Output:** A dictionary with `token_count`, `destination`, and `provenance`.

**Example code:**

```python
import pathlib

from dnakit.annotation import write_region_tokens

pathlib.Path("token_ids.gtok").unlink(missing_ok=True)
result = write_region_tokens([10, 0, 7], "token_ids.gtok", verbose=False)
print(result["token_count"], pathlib.Path(result["destination"]).name)
```

**Example output:**

```text
3 token_ids.gtok
```

## 13) `read_region_tokens`

- **Purpose:** Read region token ids from a GTOK binary file.
- **Method:** Calls gtars `utils.read_tokens_from_gtok`. Ids come back in the order they were stored. The string form gtars offers (`read_tokens_from_gtok_as_strings`) holds the same ids as text, so it is not wrapped separately.
- **API:** `read_region_tokens(path[required], max_rows[optional], verbose[optional])`
- **Input:** A local GTOK file, for example one written by `write_region_tokens`.
- **Output:** A dictionary with `token_count`, `tokens` (integers in stored order), and `provenance`.

**Example code:**

```python
import pathlib

from dnakit.annotation import read_region_tokens, write_region_tokens

pathlib.Path("tokens.gtok").unlink(missing_ok=True)
write_region_tokens([3, 1, 2], "tokens.gtok", verbose=False)
result = read_region_tokens("tokens.gtok", verbose=False)
print(result["token_count"], result["tokens"])
```

**Example output:**

```text
3 [3, 1, 2]
```

## 14) `vrs_location_digest`

- **Purpose:** Return the GA4GH VRS SequenceLocation digest of one interval on a sequence.
- **Method:** Calls gtars `vrs.location_digest` with the sequence digest (GA4GH sha512t24u), a 0-based start, and an exclusive end. The identifier is `ga4gh:SL.<digest>` unless `prefixed=False`.
- **API:** `vrs_location_digest(sequence_digest_value[required], start[required], end[required], prefixed[optional], verbose[optional])`
- **Input:** A 32-character sha512t24u sequence digest (for example from `sequence_digest`), plus integers with 0 <= start < end.
- **Output:** A dictionary with `identifier`, `digest`, `prefixed`, `sequence_digest`, `start`, `end`, and `provenance`.

**Example code:**

```python
from gtars.refget import sha512t24u_digest

from dnakit.annotation import vrs_location_digest

digest = str(sha512t24u_digest("ACGTACGTAC"))
result = vrs_location_digest(digest, 0, 4, verbose=False)
print(result["identifier"])
bare = vrs_location_digest(digest, 0, 4, prefixed=False, verbose=False)
print(bare["digest"] == result["digest"])
```

**Example output:**

```text
ga4gh:SL.xT-rAZC83CmO1DXxq6GGIm327TB1xI52
True
```
