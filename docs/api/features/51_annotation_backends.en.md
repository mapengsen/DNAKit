# Optional interval, annotation, and region-set backends

This group of functions in `dnakit.annotation` adapts optional packages: [bioframe](https://github.com/open2c/bioframe), [pybedtools](https://github.com/daler/pybedtools), [gtfparse](https://github.com/openvax/gtfparse), [gffutils](https://github.com/daler/gffutils), [pyBigWig](https://github.com/deeptools/pyBigWig), [pyranges1](https://github.com/pyranges/pyranges1), [polars-bio](https://github.com/biodatageeks/polars-bio), [gtars](https://github.com/databio/gtars), [ncls](https://github.com/pyranges/ncls), and [genomicranges](https://github.com/BiocPy/genomicranges). They compute merges, nearest intervals, coverage, intersections, complements, and sequence extraction for BED intervals; GTF and GFF parsing and database queries; bigWig value reads; Jaccard similarity, gaps, and tokenization for region sets; and interval index queries. Every function works only on explicit local files. None of them uses the network, aligns sequences, calls peaks, or downloads annotation databases.

The dependencies are optional and are imported only when the matching function is called. If a package is missing, the function raises an error that includes the install command. The pybedtools functions also need the `bedtools` executable installed separately.

## Example data

The examples below share these files. `a.bed` and `b.bed` are two small interval sets, `seqs.bed` and `ref.fa` are for sequence extraction, `genes.gff3` and `genes.gtf` are for annotation, and `universe.bed` is a BED3 region set for tokenization. The bigWig file and the gffutils databases are created by the example code at run time.

**Example file:** `a.bed`

```text
chr1	10	20	a	5	+
chr1	15	30	b	7	+
chr1	40	50	c	1	-
chr2	5	9	d	2	+
```

**Example file:** `b.bed`

```text
chr1	12	18	x
chr1	45	60	y
chr2	100	120	z
```

**Example file:** `seqs.bed`

```text
chr1	10	20	p1	0	+
chr1	0	5	p2	0	-
chr2	5	9	p3	0	+
```

**Example file:** `ref.fa`

```text
>chr1
ACGTACGTAC
GGGGTTTTCC
>chr2
ATATATATAT
```

**Example file:** `universe.bed`

```text
chr1	12	18
chr1	45	60
chr2	100	120
```

**Example file:** `genes.gff3`

```text
##gff-version 3
chr1	src	gene	100	500	.	+	.	ID=gene1;Name=GENE1
chr1	src	mRNA	100	500	.	+	.	ID=tx1;Parent=gene1;Name=GENE1-201
chr1	src	exon	100	200	.	+	.	ID=exon1;Parent=tx1
chr1	src	exon	300	500	.	+	.	ID=exon2;Parent=tx1
chr1	src	gene	800	900	.	-	.	ID=gene2;Name=GENE2
chr1	src	mRNA	800	900	.	-	.	ID=tx2;Parent=gene2
chr1	src	exon	800	900	.	-	.	ID=exon3;Parent=tx2
```

**Example file:** `genes.gtf`

```text
chr1	ens	gene	100	500	.	+	.	gene_id "ENSG1"; gene_name "GENE1"; gene_biotype "protein_coding";
chr1	ens	transcript	100	500	.	+	.	gene_id "ENSG1"; transcript_id "ENST1"; gene_name "GENE1"; gene_biotype "protein_coding"; transcript_biotype "protein_coding"; transcript_version "2";
chr1	ens	exon	100	200	.	+	.	gene_id "ENSG1"; transcript_id "ENST1"; exon_number "1"; gene_name "GENE1"; gene_biotype "protein_coding"; transcript_biotype "protein_coding";
chr1	ens	exon	300	500	.	+	.	gene_id "ENSG1"; transcript_id "ENST1"; exon_number "2"; gene_name "GENE1"; gene_biotype "protein_coding"; transcript_biotype "protein_coding";
chr1	ens	gene	800	900	.	-	.	gene_id "ENSG2"; gene_name "GENE2"; gene_biotype "lncRNA";
chr1	ens	transcript	800	900	.	-	.	gene_id "ENSG2"; transcript_id "ENST2"; gene_name "GENE2"; gene_biotype "lncRNA"; transcript_biotype "lncRNA";
chr1	ens	exon	800	900	.	-	.	gene_id "ENSG2"; transcript_id "ENST2"; exon_number "1"; gene_name "GENE2"; gene_biotype "lncRNA"; transcript_biotype "lncRNA";
```

## 1) `dnakit.annotation.merge_bioframe_intervals()`

- **Purpose:** Merge overlapping or adjacent BED intervals.
- **Method:** Call `bioframe.merge`. With `min_distance=0`, bookended intervals are merged too. With a larger value, intervals separated by no more than that gap are merged.
- **API:** `dnakit.annotation.merge_bioframe_intervals(source, *, min_distance=0, max_rows=100000)`.
- **Input:** A BED3 to BED6 file. `min_distance` is a non-negative integer. `max_rows` limits the number of input rows.
- **Output:** An `intervals` list. Each item has `chrom`, `start`, `end`, and `n_intervals`, plus `provenance`.
- **Example code:**

```python
from dnakit.annotation import merge_bioframe_intervals

result = merge_bioframe_intervals("a.bed")
print([(row["chrom"], row["start"], row["end"], row["n_intervals"]) for row in result["intervals"]])
```

- **Example output:**

```text
[('chr1', 10, 30, 2), ('chr1', 40, 50, 1), ('chr2', 5, 9, 1)]
```

## 2) `dnakit.annotation.closest_bioframe_intervals()`

- **Purpose:** For each source interval, find the k closest intervals in another file.
- **Method:** Call `bioframe.closest`. Distances use half-open coordinates, so a separated pair has distance `start_b - end_a`; overlapping pairs have distance 0. `ignore_overlaps=True` skips overlapping intervals.
- **API:** `dnakit.annotation.closest_bioframe_intervals(source, other, *, k=1, ignore_overlaps=False, max_rows=100000)`.
- **Input:** Two BED files. `k` is between 1 and 100. `ignore_overlaps` is a boolean.
- **Output:** A `records` list. Each item has the source columns, the closest interval's columns (suffixed with `_` for `other`), and `distance`.
- **Example code:**

```python
from dnakit.annotation import closest_bioframe_intervals

result = closest_bioframe_intervals("a.bed", "b.bed", k=1)
print([(row["name"], row["name_"], row["distance"]) for row in result["records"]])
```

- **Example output:**

```text
[('a', 'x', 0), ('b', 'x', 0), ('c', 'y', 0), ('d', 'z', 91)]
```

## 3) `dnakit.annotation.coverage_bioframe_intervals()`

- **Purpose:** Count how many bases of another file cover each source interval.
- **Method:** Call `bioframe.coverage`. The `coverage` column is measured in bases.
- **API:** `dnakit.annotation.coverage_bioframe_intervals(source, other, *, max_rows=100000)`.
- **Input:** Two BED files.
- **Output:** A `records` list. Each source interval has a `coverage` value.
- **Example code:**

```python
from dnakit.annotation import coverage_bioframe_intervals

result = coverage_bioframe_intervals("a.bed", "b.bed")
print([(row["name"], row["coverage"]) for row in result["records"]])
```

- **Example output:**

```text
[('a', 6), ('b', 3), ('c', 5), ('d', 0)]
```

## 4) `dnakit.annotation.merge_pybedtools_intervals()`

- **Purpose:** Merge BED intervals with BEDTools.
- **Method:** pybedtools sorts the file, then calls `bedtools merge -d`. Intervals separated by no more than `distance` are merged. The output keeps only the chrom, start, and end columns.
- **API:** `dnakit.annotation.merge_pybedtools_intervals(source, *, distance=0, max_rows=100000)`.
- **Input:** A BED file. `distance` is a non-negative integer. The `bedtools` executable is required.
- **Output:** An `intervals` list. Each item has `chrom`, `start`, and `end`.
- **Example code:**

```python
from dnakit.annotation import merge_pybedtools_intervals

result = merge_pybedtools_intervals("a.bed", distance=10)
print([(row["chrom"], row["start"], row["end"]) for row in result["intervals"]])
```

- **Example output:**

```text
[('chr1', 10, 50), ('chr2', 5, 9)]
```

## 5) `dnakit.annotation.closest_pybedtools_intervals()`

- **Purpose:** For each source interval, find the k closest intervals in another file and return the BEDTools distance.
- **Method:** Call `bedtools closest -d -t first`. Unlike `bioframe`, BEDTools adds 1 to the distance of separated intervals, giving `start_b - end_a + 1`.
- **API:** `dnakit.annotation.closest_pybedtools_intervals(source, other, *, k=1, max_rows=100000)`.
- **Input:** Two BED files. `k` is between 1 and 100. The `bedtools` executable is required.
- **Output:** A `records` list. Each item has `a` (the source fields), `b` (the closest fields, or `None` when missing), and `distance`.
- **Example code:**

```python
from dnakit.annotation import closest_pybedtools_intervals

result = closest_pybedtools_intervals("a.bed", "b.bed", k=1)
print([(row["a"][3], row["b"][3], row["distance"]) for row in result["records"]])
```

- **Example output:**

```text
[('a', 'x', 0), ('b', 'x', 0), ('c', 'y', 0), ('d', 'z', 92)]
```

## 6) `dnakit.annotation.coverage_pybedtools_intervals()`

- **Purpose:** Count the bases of each source interval covered by another file, and the covered fraction.
- **Method:** Call `bedtools coverage`. `fraction` is the number of covered bases divided by the interval length.
- **API:** `dnakit.annotation.coverage_pybedtools_intervals(source, other, *, max_rows=100000)`.
- **Input:** Two BED files. The `bedtools` executable is required.
- **Output:** A `records` list. Each item has `a`, `overlaps`, `bases_covered`, `length`, and `fraction`.
- **Example code:**

```python
from dnakit.annotation import coverage_pybedtools_intervals

result = coverage_pybedtools_intervals("a.bed", "b.bed")
print([(row["a"][3], row["bases_covered"], row["fraction"]) for row in result["records"]])
```

- **Example output:**

```text
[('a', 6, 0.6), ('b', 3, 0.2), ('c', 5, 0.5), ('d', 0, 0.0)]
```

## 7) `dnakit.annotation.extract_pybedtools_sequences()`

- **Purpose:** Extract FASTA sequences for BED intervals.
- **Method:** Call `bedtools getfasta -tab`. With `use_strand=True`, intervals on the minus strand are reverse-complemented. BEDTools skips intervals that run past the end of a sequence; this function raises an error instead of dropping them silently. If the FASTA file has no `.fai` index, BEDTools creates `<fasta>.fai` next to it.
- **API:** `dnakit.annotation.extract_pybedtools_sequences(source, fasta, *, use_strand=False, max_intervals=1000, max_sequence_length=100000)`.
- **Input:** A BED file and a FASTA file. Each interval must be no longer than `max_sequence_length`. The `bedtools` executable is required.
- **Output:** A `sequences` list. Each item has `name` (`chrom:start-end`), `sequence`, and `length`.
- **Example code:**

```python
from dnakit.annotation import extract_pybedtools_sequences

result = extract_pybedtools_sequences("seqs.bed", "ref.fa", use_strand=True)
print([(row["name"], row["sequence"]) for row in result["sequences"]])
```

- **Example output:**

```text
[('chr1:10-20(+)', 'GGGGTTTTCC'), ('chr1:0-5(-)', 'TACGT'), ('chr2:5-9(+)', 'TATA')]
```

## 8) `dnakit.annotation.read_gtfparse_table()`

- **Purpose:** Parse a GTF file into records, with one column for each attribute key.
- **Method:** Call `gtfparse.read_gtf`. `features` keeps only the listed feature types. With `normalize_gencode_biotypes=True`, GENCODE's `gene_type` and `transcript_type` are renamed to `gene_biotype` and `transcript_biotype`. Attributes that gtfparse leaves empty are reported as `None`.
- **API:** `dnakit.annotation.read_gtfparse_table(source, *, features=(), normalize_gencode_biotypes=False, max_rows=100000)`.
- **Input:** A GTF file (it may be `.gz`). `max_rows` limits the number of data lines in the file.
- **Output:** `columns`, `records`, `record_count`, and `provenance`. Coordinates are GTF 1-based inclusive values.
- **Example code:**

```python
from dnakit.annotation import read_gtfparse_table

result = read_gtfparse_table("genes.gtf", features=["exon"])
print(result["record_count"], [(row["gene_name"], row["start"], row["end"]) for row in result["records"]])
```

- **Example output:**

```text
3 [('GENE1', 100, 200), ('GENE1', 300, 500), ('GENE2', 800, 900)]
```

## 9) `dnakit.annotation.select_gtfparse_records()`

- **Purpose:** Select GTF records whose gene ID, gene name, or transcript ID appears in the given lists.
- **Method:** After reading the GTF, keep a record when its `gene_id`, `gene_name`, or `transcript_id` matches any listed value. At least one list must be given.
- **API:** `dnakit.annotation.select_gtfparse_records(source, *, gene_ids=(), gene_names=(), transcript_ids=(), features=(), max_rows=100000)`.
- **Input:** A GTF file, three lists of IDs or names, and an optional `features` filter.
- **Output:** `records`, `record_count`, `gene_ids`, and `transcript_ids`.
- **Example code:**

```python
from dnakit.annotation import select_gtfparse_records

result = select_gtfparse_records("genes.gtf", gene_names=["GENE2"])
print(result["gene_ids"], [(row["feature"], row["transcript_id"]) for row in result["records"]])
```

- **Example output:**

```text
['ENSG2'] [('gene', None), ('transcript', 'ENST2'), ('exon', 'ENST2')]
```

## 10) `dnakit.annotation.summarize_gtfparse_annotation()`

- **Purpose:** Count feature types, genes, transcripts, gene biotypes, and chromosomes in a GTF file.
- **Method:** Genes and transcripts are counted by distinct `gene_id` and `transcript_id`. Each gene contributes once to the biotype counts.
- **API:** `dnakit.annotation.summarize_gtfparse_annotation(source, *, max_rows=100000)`.
- **Input:** A GTF file.
- **Output:** `feature_counts`, `gene_count`, `transcript_count`, `gene_biotype_counts`, and `chromosome_counts`.
- **Example code:**

```python
from dnakit.annotation import summarize_gtfparse_annotation

result = summarize_gtfparse_annotation("genes.gtf")
print(result["feature_counts"], result["gene_biotype_counts"])
```

- **Example output:**

```text
{'gene': 2, 'transcript': 2, 'exon': 3} {'protein_coding': 1, 'lncRNA': 1}
```

## 11) `dnakit.annotation.create_gffutils_database()`

- **Purpose:** Create a gffutils SQLite database from a GFF3 or GTF file.
- **Method:** Call `gffutils.create_db`; gffutils detects the format. The database path must not exist yet, so an existing database is never overwritten. The input size is limited by `max_input_bytes`.
- **API:** `dnakit.annotation.create_gffutils_database(source, database, *, max_input_bytes=100000000)`.
- **Input:** A GFF3 or GTF file. `database` is a new database path whose parent directory must exist.
- **Output:** `database`, `input_format`, `featuretype_counts`, and `feature_count`. Coordinates are GFF 1-based inclusive values.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.annotation import create_gffutils_database

target = Path(tempfile.mkdtemp()) / "genes.db"
result = create_gffutils_database("genes.gff3", target)
print(result["input_format"], result["featuretype_counts"], result["feature_count"])
```

- **Example output:**

```text
gff3 {'exon': 3, 'gene': 2, 'mRNA': 2} 7
```

## 12) `dnakit.annotation.list_gffutils_children()`

- **Purpose:** List the child features of one feature in an existing gffutils database.
- **Method:** Call `FeatureDB.children`. `level=1` returns direct children and `level=2` returns grandchildren. A GTF-derived database lists each transcript as its own child; this function removes that self-link.
- **API:** `dnakit.annotation.list_gffutils_children(database, parent_id, *, featuretype=None, level=1, max_features=1000)`.
- **Input:** A gffutils database path, a parent feature ID, an optional `featuretype`, and `level` between 1 and 10.
- **Output:** A `records` list. Each item has `id`, `featuretype`, `seqid`, `start`, `end`, and `strand`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.annotation import create_gffutils_database, list_gffutils_children

database = Path(tempfile.mkdtemp()) / "genes.db"
create_gffutils_database("genes.gff3", database)
result = list_gffutils_children(database, "gene1", level=2)
print([(row["id"], row["featuretype"], row["start"], row["end"]) for row in result["records"]])
```

- **Example output:**

```text
[('exon1', 'exon', 100, 200), ('exon2', 'exon', 300, 500)]
```

## 13) `dnakit.annotation.count_gffutils_featuretypes()`

- **Purpose:** Count the features of each feature type in an existing gffutils database.
- **Method:** Read `featuretypes()` from the database, then call `count_features_of_type` for each type.
- **API:** `dnakit.annotation.count_gffutils_featuretypes(database)`.
- **Input:** A gffutils database path.
- **Output:** A `featuretype_counts` dictionary and the total `feature_count`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.annotation import count_gffutils_featuretypes, create_gffutils_database

database = Path(tempfile.mkdtemp()) / "genes.db"
create_gffutils_database("genes.gff3", database)
result = count_gffutils_featuretypes(database)
print(result["featuretype_counts"], result["feature_count"])
```

- **Example output:**

```text
{'exon': 3, 'gene': 2, 'mRNA': 2} 7
```

## 14) `dnakit.annotation.read_pybigwig_values()`

- **Purpose:** Read the bigWig value at every base of a region.
- **Method:** Call `pyBigWig`'s `values`. Bases without data are returned as `None`. The region length is limited by `max_bases`.
- **API:** `dnakit.annotation.read_pybigwig_values(source, chromosome, start, end, *, max_bases=100000)`.
- **Input:** A bigWig file, a chromosome name, and the 0-based half-open region `[start, end)`. `end` must not exceed the chromosome length.
- **Output:** `values` (one entry per base), `covered_bases`, and `length`.
- **Example code:**

```python
import pyBigWig

from dnakit.annotation import read_pybigwig_values

handle = pyBigWig.open("signal.bw", "w")
handle.addHeader([("chr1", 100), ("chr2", 50)])
handle.addEntries(["chr1", "chr1"], [10, 20], ends=[15, 25], values=[1.0, 3.0])
handle.addEntries(["chr2"], [0], ends=[10], values=[2.0])
handle.close()
result = read_pybigwig_values("signal.bw", "chr1", 8, 30)
print(result["covered_bases"], result["values"][:8])
```

- **Example output:**

```text
10 [None, None, 1.0, 1.0, 1.0, 1.0, 1.0, None]
```

## 15) `dnakit.annotation.read_pybigwig_intervals()`

- **Purpose:** Return the stored bigWig intervals that overlap a region.
- **Method:** Call `pyBigWig`'s `intervals`. The returned intervals are not clipped to the query region.
- **API:** `dnakit.annotation.read_pybigwig_intervals(source, chromosome, start, end, *, max_bases=1000000, max_intervals=10000)`.
- **Input:** A bigWig file, a chromosome name, and a 0-based half-open region. `max_bases` and `max_intervals` are explicit limits.
- **Output:** A `records` list. Each item has `start`, `end`, and `value`.
- **Example code:**

```python
import pyBigWig

from dnakit.annotation import read_pybigwig_intervals

handle = pyBigWig.open("signal.bw", "w")
handle.addHeader([("chr1", 100), ("chr2", 50)])
handle.addEntries(["chr1", "chr1"], [10, 20], ends=[15, 25], values=[1.0, 3.0])
handle.addEntries(["chr2"], [0], ends=[10], values=[2.0])
handle.close()
result = read_pybigwig_intervals("signal.bw", "chr1", 8, 30)
print([(row["start"], row["end"], row["value"]) for row in result["records"]])
```

- **Example output:**

```text
[(10, 15, 1.0), (20, 25, 3.0)]
```

## 16) `dnakit.annotation.summarize_pybigwig_bed_intervals()`

- **Purpose:** Compute one bigWig statistic for every interval in a BED file.
- **Method:** Call `pyBigWig`'s `stats` for each interval. `statistic` is `mean`, `min`, `max`, `coverage`, `std`, or `sum`. When an interval has no data, `value` is `None`.
- **API:** `dnakit.annotation.summarize_pybigwig_bed_intervals(source, bed, *, statistic="mean", max_intervals=10000)`.
- **Input:** A bigWig file and a BED file. Every interval must lie within the bigWig chromosome lengths.
- **Output:** A `records` list. Each item has `chrom`, `start`, `end`, `name`, and `value`.
- **Example code:**

```python
import pyBigWig

from dnakit.annotation import summarize_pybigwig_bed_intervals

handle = pyBigWig.open("signal.bw", "w")
handle.addHeader([("chr1", 100), ("chr2", 50)])
handle.addEntries(["chr1", "chr1"], [10, 20], ends=[15, 25], values=[1.0, 3.0])
handle.addEntries(["chr2"], [0], ends=[10], values=[2.0])
handle.close()
result = summarize_pybigwig_bed_intervals("signal.bw", "seqs.bed", statistic="sum")
print([(row["name"], row["value"]) for row in result["records"]])
```

- **Example output:**

```text
[('p1', 5.0), ('p2', None), ('p3', 8.0)]
```

## 17) `dnakit.annotation.intersect_pyranges1_intervals()`

- **Purpose:** Return the parts of source intervals that overlap intervals in another file.
- **Method:** Call `pyranges1`'s `intersect_overlaps`. The output intervals are clipped to the overlap and keep pyranges column names such as `Chromosome`, `Start`, and `End`.
- **API:** `dnakit.annotation.intersect_pyranges1_intervals(source, other, *, max_rows=100000)`.
- **Input:** Two BED files.
- **Output:** A `records` list.
- **Example code:**

```python
from dnakit.annotation import intersect_pyranges1_intervals

result = intersect_pyranges1_intervals("a.bed", "b.bed")
print([(row["Chromosome"], row["Start"], row["End"], row["Name"]) for row in result["records"]])
```

- **Example output:**

```text
[('chr1', 12, 18, 'a'), ('chr1', 15, 18, 'b'), ('chr1', 45, 50, 'c')]
```

## 18) `dnakit.annotation.nearest_pyranges1_intervals()`

- **Purpose:** For each source interval, find the k nearest intervals in another file.
- **Method:** Call `pyranges1`'s `nearest_ranges`. The `Distance` column follows the BEDTools convention, which adds 1 for separated intervals.
- **API:** `dnakit.annotation.nearest_pyranges1_intervals(source, other, *, k=1, max_rows=100000)`.
- **Input:** Two BED files. `k` is between 1 and 100.
- **Output:** A `records` list. Each item has the source columns, the `other` columns with a `_b` suffix, and `Distance`.
- **Example code:**

```python
from dnakit.annotation import nearest_pyranges1_intervals

result = nearest_pyranges1_intervals("a.bed", "b.bed", k=1)
print([row["Distance"] for row in result["records"]])
```

- **Example output:**

```text
[0, 0, 0, 92]
```

## 19) `dnakit.annotation.complement_pyranges1_intervals()`

- **Purpose:** Return the parts of each chromosome that no source interval covers.
- **Method:** Call `pyranges1`'s `complement_ranges`. Strands are ignored, and the leading gap from position 0 is included. Chromosome lengths come from `chrom_sizes`.
- **API:** `dnakit.annotation.complement_pyranges1_intervals(source, *, chrom_sizes, max_rows=100000)`.
- **Input:** A BED file. `chrom_sizes` is a `{chromosome: length}` mapping.
- **Output:** A `records` list. Each item has `Chromosome`, `Start`, and `End`.
- **Example code:**

```python
from dnakit.annotation import complement_pyranges1_intervals

result = complement_pyranges1_intervals("a.bed", chrom_sizes={"chr1": 100, "chr2": 50})
print([(row["Chromosome"], row["Start"], row["End"]) for row in result["records"]])
```

- **Example output:**

```text
[('chr1', 0, 10), ('chr1', 30, 40), ('chr1', 50, 100), ('chr2', 0, 5), ('chr2', 9, 50)]
```

## 20) `dnakit.annotation.merge_polars_bio_intervals()`

- **Purpose:** Merge BED intervals with polars-bio.
- **Method:** Read the BED with `use_zero_based=True` so the coordinates match the BED file. polars-bio merges only gaps that are strictly smaller than `min_dist`, so this function passes `min_distance + 1`. The result is the same rule as the other merge functions: gaps no larger than `min_distance` are merged.
- **API:** `dnakit.annotation.merge_polars_bio_intervals(source, *, min_distance=0, max_rows=100000)`.
- **Input:** A BED file. `min_distance` is a non-negative integer.
- **Output:** An `intervals` list. Each item has `chrom`, `start`, `end`, and `n_intervals`.
- **Example code:**

```python
from dnakit.annotation import merge_polars_bio_intervals

result = merge_polars_bio_intervals("a.bed", min_distance=10)
print([(row["chrom"], row["start"], row["end"], row["n_intervals"]) for row in result["intervals"]])
```

- **Example output:**

```text
[('chr1', 10, 50, 3), ('chr2', 5, 9, 1)]
```

## 21) `dnakit.annotation.nearest_polars_bio_intervals()`

- **Purpose:** For each source interval, find the k nearest intervals in another file.
- **Method:** Call `polars_bio.nearest`. With `overlap=True`, overlapping intervals have distance 0. Distances use 0-based half-open coordinates, so a separated pair has distance `start_b - end_a`.
- **API:** `dnakit.annotation.nearest_polars_bio_intervals(source, other, *, k=1, max_rows=100000)`.
- **Input:** Two BED files. `k` is between 1 and 100.
- **Output:** A `records` list. Each item has the `_1` and `_2` suffixed columns and `distance`.
- **Example code:**

```python
from dnakit.annotation import nearest_polars_bio_intervals

result = nearest_polars_bio_intervals("a.bed", "b.bed", k=1)
print([row["distance"] for row in result["records"]])
```

- **Example output:**

```text
[0, 0, 0, 91]
```

## 22) `dnakit.annotation.coverage_polars_bio_intervals()`

- **Purpose:** Count the bases of each source interval covered by another file.
- **Method:** Call `polars_bio.coverage`. Both inputs are read with 0-based coordinates.
- **API:** `dnakit.annotation.coverage_polars_bio_intervals(source, other, *, max_rows=100000)`.
- **Input:** Two BED files.
- **Output:** A `records` list. Each source interval has a `coverage` value.
- **Example code:**

```python
from dnakit.annotation import coverage_polars_bio_intervals

result = coverage_polars_bio_intervals("a.bed", "b.bed")
print([row["coverage"] for row in result["records"]])
```

- **Example output:**

```text
[6, 3, 5, 0]
```

## 23) `dnakit.annotation.jaccard_gtars_regions()`

- **Purpose:** Compute the base-pair Jaccard similarity of two BED region sets.
- **Method:** Call `gtars`'s `RegionSet.jaccard`. The value is the number of bases covered by both sets divided by the number covered by either set.
- **API:** `dnakit.annotation.jaccard_gtars_regions(source, other, *, max_rows=100000)`.
- **Input:** Two BED files.
- **Output:** `jaccard`, `source_regions`, and `other_regions`.
- **Example code:**

```python
from dnakit.annotation import jaccard_gtars_regions

result = jaccard_gtars_regions("a.bed", "b.bed")
print(result["jaccard"])
```

- **Example output:**

```text
0.171875
```

## 24) `dnakit.annotation.gaps_gtars_regions()`

- **Purpose:** Return the gaps in a region set for each chromosome.
- **Method:** Call `gtars`'s `RegionSet.gaps`. Gaps include the leading segment, the spaces between regions, and the trailing segment. A chromosome with no regions gives one gap spanning its whole length.
- **API:** `dnakit.annotation.gaps_gtars_regions(source, chrom_sizes, *, max_rows=100000)`.
- **Input:** A BED file. `chrom_sizes` is a `{chromosome: length}` mapping.
- **Output:** A `gaps` list. Each item has `chrom`, `start`, and `end` (0-based half-open).
- **Example code:**

```python
from dnakit.annotation import gaps_gtars_regions

result = gaps_gtars_regions("a.bed", {"chr1": 100, "chr2": 50})
print([(row["chrom"], row["start"], row["end"]) for row in result["gaps"]])
```

- **Example output:**

```text
[('chr1', 0, 10), ('chr1', 30, 40), ('chr1', 50, 100), ('chr2', 0, 5), ('chr2', 9, 50)]
```

## 25) `dnakit.annotation.tokenize_gtars_regions()`

- **Purpose:** Map BED regions onto a BED3 universe of regions and return one token for each overlap.
- **Method:** Read the universe with `gtars.tokenizers.Tokenizer.from_bed`, then call `tokenize`. Each query region that overlaps a universe region produces a token. Query regions with no overlap are dropped, and repeated tokens are kept.
- **API:** `dnakit.annotation.tokenize_gtars_regions(source, universe, *, max_rows=100000)`.
- **Input:** A query BED file. `universe` must be a BED3 file.
- **Output:** `tokens` (in the form `chr1:12-18`), `token_ids`, and `token_count`.
- **Example code:**

```python
from dnakit.annotation import tokenize_gtars_regions

result = tokenize_gtars_regions("a.bed", "universe.bed")
print(result["tokens"])
```

- **Example output:**

```text
['chr1:12-18', 'chr1:12-18', 'chr1:45-60']
```

## 26) `dnakit.annotation.query_ncls_overlaps()`

- **Purpose:** Return every overlapping pair between source intervals and query intervals.
- **Method:** Build one NCLS index per chromosome, using each source row number as its ID, then call `all_overlaps_both`. Coordinates are 0-based half-open, so intervals that only touch do not overlap.
- **API:** `dnakit.annotation.query_ncls_overlaps(source, queries, *, max_rows=100000, max_pairs=100000)`.
- **Input:** A source BED file and a query BED file. `max_pairs` limits the number of output pairs.
- **Output:** A `records` list. Each item has `chrom`, `query_index`, `target_index`, and the coordinates of both intervals. Indexes are zero-based row numbers within each file.
- **Example code:**

```python
from dnakit.annotation import query_ncls_overlaps

result = query_ncls_overlaps("a.bed", "b.bed")
print([(row["query_index"], row["target_index"]) for row in result["records"]])
```

- **Example output:**

```text
[(0, 0), (0, 1), (1, 2)]
```

## 27) `dnakit.annotation.count_ncls_overlaps()`

- **Purpose:** Count the source intervals that overlap each query interval.
- **Method:** Use the same index as `query_ncls_overlaps` and count matches per query row.
- **API:** `dnakit.annotation.count_ncls_overlaps(source, queries, *, max_rows=100000)`.
- **Input:** A source BED file and a query BED file.
- **Output:** A `records` list. Each query interval has an `overlap_count`, in input order.
- **Example code:**

```python
from dnakit.annotation import count_ncls_overlaps

result = count_ncls_overlaps("a.bed", "b.bed")
print([(row["name"], row["overlap_count"]) for row in result["records"]])
```

- **Example output:**

```text
[('x', 2), ('y', 1), ('z', 0)]
```

## 28) `dnakit.annotation.find_ncls_containing_intervals()`

- **Purpose:** Return the source intervals that fully contain each query interval.
- **Method:** Use the same index as `query_ncls_overlaps`, but call `all_containments_both`. A source interval matches when its start is not after the query start and its end is not before the query end.
- **API:** `dnakit.annotation.find_ncls_containing_intervals(source, queries, *, max_rows=100000, max_pairs=100000)`.
- **Input:** A source BED file and a query BED file.
- **Output:** A `records` list. Each item has `query_index` and `target_index`.
- **Example code:**

```python
from dnakit.annotation import find_ncls_containing_intervals

result = find_ncls_containing_intervals("a.bed", "b.bed")
print([(row["query_index"], row["target_index"]) for row in result["records"]])
```

- **Example output:**

```text
[(0, 0)]
```

## 29) `dnakit.annotation.reduce_genomicranges_intervals()`

- **Purpose:** Merge overlapping or adjacent intervals.
- **Method:** Convert the BED to `genomicranges.GenomicRanges` (1-based inclusive), call `reduce`, and convert back to 0-based half-open coordinates. With `min_gap_width=1`, overlapping and bookended intervals are merged; with `0`, only overlapping intervals are merged. `use_strand=False` ignores strand.
- **API:** `dnakit.annotation.reduce_genomicranges_intervals(source, *, min_gap_width=1, use_strand=False, max_rows=100000)`.
- **Input:** A BED file.
- **Output:** An `intervals` list. Each item has `chrom`, `start`, `end`, and `strand`. Metadata columns are not kept.
- **Example code:**

```python
from dnakit.annotation import reduce_genomicranges_intervals

result = reduce_genomicranges_intervals("a.bed")
print([(row["chrom"], row["start"], row["end"]) for row in result["intervals"]])
```

- **Example output:**

```text
[('chr1', 10, 30), ('chr1', 40, 50), ('chr2', 5, 9)]
```

## 30) `dnakit.annotation.intersect_genomicranges_intervals()`

- **Purpose:** Return the positions covered by both interval sets.
- **Method:** Call `GenomicRanges.intersect`. The result is the set intersection of the covered positions, not a pairwise intersection of intervals.
- **API:** `dnakit.annotation.intersect_genomicranges_intervals(source, other, *, use_strand=False, max_rows=100000)`.
- **Input:** Two BED files.
- **Output:** An `intervals` list. Each item has `chrom`, `start`, `end`, and `strand`.
- **Example code:**

```python
from dnakit.annotation import intersect_genomicranges_intervals

result = intersect_genomicranges_intervals("a.bed", "b.bed")
print([(row["chrom"], row["start"], row["end"]) for row in result["intervals"]])
```

- **Example output:**

```text
[('chr1', 12, 18), ('chr1', 45, 50)]
```

## 31) `dnakit.annotation.flank_genomicranges_intervals()`

- **Purpose:** Create flanking intervals of a fixed width beside each interval.
- **Method:** Call `GenomicRanges.flank`. `side="start"` gives the left flank, or the upstream flank when `use_strand=True`. `side="end"` gives the right flank, or the downstream flank. Flanks are not clipped to the chromosome length, and a flank that starts before position 0 raises an error.
- **API:** `dnakit.annotation.flank_genomicranges_intervals(source, width, *, side="start", use_strand=False, max_rows=100000)`.
- **Input:** A BED file. `width` is a positive integer. `side` is `"start"` or `"end"`.
- **Output:** An `intervals` list. Each item has `chrom`, `start`, `end`, and `strand`.
- **Example code:**

```python
from dnakit.annotation import flank_genomicranges_intervals

result = flank_genomicranges_intervals("a.bed", 5)
print([(row["chrom"], row["start"], row["end"]) for row in result["intervals"]])
```

- **Example output:**

```text
[('chr1', 5, 10), ('chr1', 10, 15), ('chr1', 35, 40), ('chr2', 0, 5)]
```

## Installing the optional dependencies

```bash
python -m pip install "dnakit[bio-bioframe]"        # bioframe
python -m pip install "dnakit[bio-pybedtools]"      # pybedtools; also needs the bedtools executable
python -m pip install "dnakit[bio-gtfparse]"        # gtfparse
python -m pip install "dnakit[bio-gffutils]"        # gffutils
python -m pip install "dnakit[bio-pybigwig]"        # pyBigWig
python -m pip install "dnakit[bio-pyranges1]"       # pyranges1
python -m pip install "dnakit[bio-polarsbio]"       # polars-bio
python -m pip install "dnakit[bio-gtars]"           # gtars
python -m pip install "dnakit[bio-ncls]"            # ncls and numpy
python -m pip install "dnakit[bio-genomicranges]"   # genomicranges and pandas
```

## How these functions are verified

The test file is `others/tests/unit/test_annotation_backends.py`, with 57 tests. They check the following:

- Each public function calls the real upstream package at least once. Results are compared with hand-calculated values or with direct upstream calls. Examples include bioframe and pybedtools base coverage, the gtars Jaccard value of 11/64, and ncls overlaps and containment on random intervals checked against brute-force enumeration.
- The BEDTools distance is one larger than the bioframe distance for separated intervals, and the polars-bio merge distance rule is checked against the other merge functions.
- gffutils GTF self-links, gtfparse empty attributes, and bigWig bases without data each have a test.
- Each module has a test showing that a missing optional package raises `BackendUnavailableError` with the install extra in the message.
- Non-boolean and out-of-range `k` or `max_rows`, reversed intervals, BED files without intervals, intervals outside a FASTA sequence, and a non-BED3 tokenization universe are all rejected.

No network services are wrapped in this batch, so there is no `DNAKIT_LIVE_TESTS` test.

## Current scope

This batch wraps 31 functions from 10 packages. All of them are DNA or gene interval operations on local files. The following operations are not wrapped:

- bioframe `complement`, `subtract`, and `cluster`, and pybedtools `intersect`: `intersect_pybedtools_bed` already covers pybedtools intersection, and the rest fall outside this batch's range.
- gffutils `FeatureDB` queries by region: `query_gffutils_database` already covers them.
- pyranges1 `join_overlaps` and `merge_overlaps`, and polars-bio `overlap`, `subtract`, and `cluster`: `overlap_polars_bio_bed` and `merge_pyranges_bed` already cover overlaps and merging, and the rest are outside this batch.
- pyBigWig write methods, and gtars `lola`, `refget`, and `genomic_distributions`: writing and these analysis modules are outside this read-only batch.
- gtfparse `write_gtf` and the gffutils GFF writer: writing output files is outside this read-only batch.
- ncls `k_overlaps_both` and `first_overlap_both`: they repeat the overlap queries above.
