# Feature 69: G3 interval, annotation, and bigWig wrappers

This page documents the 24 G3 DNAKit wrappers: bioframe (6), pybedtools (4), gffutils (4), gtfparse (3), pyBigWig (3), and pyranges1 (4). Each function is exported from `dnakit.annotation`, reads local files only, and never touches the network.

Import each function from the package, for example `from dnakit.annotation import overlap_bioframe_intervals`; the API line gives the full signature. The pybedtools functions need the `bedtools` executable (installed here at `/usr/bin/bedtools`). Optional packages are imported only when a function runs; a missing package raises an error with the install hint. Example outputs come from real runs.

## Example data

The examples below share these files. The pyBigWig example writes its bigWig to a temporary directory at run time, so it needs no extra file.

BED6 intervals for bioframe, pybedtools, and pyranges1 examples.

**Example file:** `bf_a.bed`

```text
chr1	10	20	a	5	+
chr1	15	30	b	7	+
chr1	40	50	c	1	-
chr2	5	9	d	2	+
```

BED4 intervals used as the second input of overlap-type examples.

**Example file:** `bf_b.bed`

```text
chr1	15	18	x
chr1	45	60	y
chr2	100	120	z
```

One BEDPE record (10 columns) for the tagging example.

**Example file:** `pbt_pairs.bedpe`

```text
chr1	10	20	chr1	100	110	p1	0	+	-
```

BED4 query intervals named q1 and q2.

**Example file:** `pbt_q.bed`

```text
chr1	12	18	q1
chr1	101	105	q2
```

Small GFF3 file: one gene with one mRNA and two exons, and one more gene.

**Example file:** `gff_genes.gff3`

```text
##gff-version 3
chr1	src	gene	100	500	.	+	.	ID=gene1;Name=GENE1
chr1	src	mRNA	100	500	.	+	.	ID=tx1;Parent=gene1;Name=GENE1-201
chr1	src	exon	100	200	.	+	.	ID=exon1;Parent=tx1
chr1	src	exon	300	500	.	+	.	ID=exon2;Parent=tx1
chr2	src	gene	50	90	.	-	.	ID=gene2
```

Small GTF file: gene, transcript, and one exon.

**Example file:** `gtp_genes.gtf`

```text
chr1	src	gene	100	500	.	+	.	gene_id "g1"; gene_name "GENE1";
chr1	src	transcript	100	500	.	+	.	gene_id "g1"; transcript_id "t1";
chr1	src	exon	100	200	.	+	.	gene_id "g1"; transcript_id "t1"; exon_number 1;
```

bedGraph rows (chrom, start, end, value) for the bigWig examples.

**Example file:** `bw_cov.bedgraph`

```text
chr1	0	10	1.5
chr1	10	20	3
chr2	5	9	2
```

## 1) `dnakit.annotation.overlap_bioframe_intervals()`

- **Purpose:** Find the pairs of overlapping intervals in two BED files.
- **Method:** Call `bioframe.overlap`. With `how="inner"` only matched rows are kept. `left`, `right`, and `outer` also keep unmatched rows, whose partner columns are empty.
- **API:** `dnakit.annotation.overlap_bioframe_intervals(source: 'str | Path', other: 'str | Path', *, how: 'str' = 'inner', max_rows: 'int' = 100000, verbose: bool = True)`
- **Input:** Two BED3 to BED6 files. `how` is inner, left, right, or outer. `max_rows` limits the rows read from each file.
- **Output:** A `records` list with one pair per item. Left columns have no suffix and right columns end with `_`. Also `record_count` and `provenance`.
- **Example code:**

```python
from dnakit.annotation import overlap_bioframe_intervals

result = overlap_bioframe_intervals("bf_a.bed", "bf_b.bed")
print(result["record_count"])
print([(row["chrom"], row["start"], row["end"]) for row in result["records"]])
```

- **Example output:**

```text
3
[('chr1', 10, 20), ('chr1', 15, 30), ('chr1', 40, 50)]
```

## 2) `dnakit.annotation.subtract_bioframe_intervals()`

- **Purpose:** Remove the bases covered by another BED file from each interval and return the remaining pieces.
- **Method:** Call `bioframe.subtract`. Removed bases are gone, and one interval can split into several pieces.
- **API:** `dnakit.annotation.subtract_bioframe_intervals(source: 'str | Path', other: 'str | Path', *, max_rows: 'int' = 100000, verbose: bool = True)`
- **Input:** Two BED3 to BED6 files. `max_rows` limits the rows read from each file.
- **Output:** An `intervals` list with the remaining pieces, plus `interval_count` and `provenance`.
- **Example code:**

```python
from dnakit.annotation import subtract_bioframe_intervals

result = subtract_bioframe_intervals("bf_a.bed", "bf_b.bed")
print([(row["chrom"], row["start"], row["end"]) for row in result["intervals"]])
```

- **Example output:**

```text
[('chr1', 10, 15), ('chr1', 18, 20), ('chr1', 18, 30), ('chr1', 40, 45), ('chr2', 5, 9)]
```

## 3) `dnakit.annotation.expand_bioframe_intervals()`

- **Purpose:** Grow each interval by `pad` bases on both sides.
- **Method:** Call `bioframe.expand`. The start moves left by `pad` and the end moves right by `pad`. Results are not clipped at 0.
- **API:** `dnakit.annotation.expand_bioframe_intervals(source: 'str | Path', *, pad: 'int', max_rows: 'int' = 100000, verbose: bool = True)`
- **Input:** One BED3 to BED6 file. `pad` is a non-negative integer. `max_rows` limits input rows.
- **Output:** An `intervals` list, `interval_count`, the `pad` value, and `provenance`.
- **Example code:**

```python
from dnakit.annotation import expand_bioframe_intervals

result = expand_bioframe_intervals("bf_a.bed", pad=5)
print([(row["chrom"], row["start"], row["end"]) for row in result["intervals"]])
```

- **Example output:**

```text
[('chr1', 5, 25), ('chr1', 10, 35), ('chr1', 35, 55), ('chr2', 0, 14)]
```

## 4) `dnakit.annotation.cluster_bioframe_intervals()`

- **Purpose:** Assign each interval to a cluster. Overlapping intervals, or intervals within `min_dist` of each other, share a cluster.
- **Method:** Call `bioframe.cluster`. Intervals on one chromosome are merged in order into clusters, and the cluster columns (`cluster`, `cluster_start`, `cluster_end`) are added.
- **API:** `dnakit.annotation.cluster_bioframe_intervals(source: 'str | Path', *, min_dist: 'int' = 0, max_rows: 'int' = 100000, verbose: bool = True)`
- **Input:** One BED3 to BED6 file. `min_dist` is a non-negative integer. `max_rows` limits input rows.
- **Output:** A `records` list with one row per input interval, including the cluster id, plus `record_count` and `provenance`.
- **Example code:**

```python
from dnakit.annotation import cluster_bioframe_intervals

result = cluster_bioframe_intervals("bf_a.bed")
print([(row["start"], row["end"], row["cluster"]) for row in result["records"]])
```

- **Example output:**

```text
[(10, 20, 0), (15, 30, 0), (40, 50, 1), (5, 9, 2)]
```

## 5) `dnakit.annotation.binnify_bioframe_chromosomes()`

- **Purpose:** Divide each chromosome of a length mapping into fixed-width bins for genome binning.
- **Method:** Call `bioframe.binnify`. Each chromosome is split from 0, and the last bin may be shorter than `binsize`.
- **API:** `dnakit.annotation.binnify_bioframe_chromosomes(chromsizes: 'Mapping[str, int]', *, binsize: 'int', max_rows: 'int' = 100000, verbose: bool = True)`
- **Input:** A mapping from chromosome name to length. `binsize` is a positive integer. `max_rows` limits the number of bins.
- **Output:** A `bins` list (chrom, start, end), `bin_count`, `binsize`, and `provenance`.
- **Example code:**

```python
from dnakit.annotation import binnify_bioframe_chromosomes

result = binnify_bioframe_chromosomes({"chr1": 25, "chr2": 10}, binsize=10)
print(result["bin_count"])
print([(row["chrom"], row["start"], row["end"]) for row in result["bins"]])
```

- **Example output:**

```text
4
[('chr1', 0, 10), ('chr1', 10, 20), ('chr1', 20, 25), ('chr2', 0, 10)]
```

## 6) `dnakit.annotation.parse_bioframe_region()`

- **Purpose:** Parse a UCSC-style region string such as `chr1:101-200`.
- **Method:** Call `bioframe.parse_region_string`. Thousands separators are removed. Numbers are returned as written, with no 1-based to 0-based shift.
- **API:** `dnakit.annotation.parse_bioframe_region(region: 'str', *, verbose: bool = True)`
- **Input:** A string of at most 200 characters.
- **Output:** `chrom`, `start`, and `end` fields, plus `provenance`. Strings that cannot be parsed raise `ConfigurationError`.
- **Example code:**

```python
from dnakit.annotation import parse_bioframe_region

result = parse_bioframe_region("chr1:1,101-2,000")
print(result["chrom"], result["start"], result["end"])
```

- **Example output:**

```text
chr1 1101 2000
```

## 7) `dnakit.annotation.five_prime_pybedtools_coordinates()`

- **Purpose:** Return a window around the 5' end of each feature.
- **Method:** Call `pybedtools.featurefuncs.five_prime` through `BedTool.each`. The 5' end is the start on + and the end on -.
- **API:** `dnakit.annotation.five_prime_pybedtools_coordinates(source: 'str | Path', *, upstream: 'int' = 500, downstream: 'int' = 500, max_rows: 'int' = 100000, verbose: bool = True)`
- **Input:** One BED3 to BED6 file. `upstream` and `downstream` are non-negative integers. `max_rows` limits input rows. Requires the BEDTools executable.
- **Output:** A `records` list (chrom, start, end, name, strand), `record_count`, and `provenance`.
- **Example code:**

```python
from dnakit.annotation import five_prime_pybedtools_coordinates

result = five_prime_pybedtools_coordinates("bf_a.bed", upstream=5, downstream=5)
print([(row["start"], row["end"]) for row in result["records"]])
```

- **Example output:**

```text
[(5, 15), (10, 20), (45, 55), (0, 10)]
```

## 8) `dnakit.annotation.midpoint_pybedtools_coordinates()`

- **Purpose:** Return the 1-bp midpoint of each BED feature.
- **Method:** Call `pybedtools.featurefuncs.midpoint`. The midpoint is `floor((start + end) / 2)`.
- **API:** `dnakit.annotation.midpoint_pybedtools_coordinates(source: 'str | Path', *, max_rows: 'int' = 100000, verbose: bool = True)`
- **Input:** One BED3 to BED6 file. `max_rows` limits input rows. Requires the BEDTools executable.
- **Output:** A `records` list where each item is a 1-bp interval, plus `record_count` and `provenance`.
- **Example code:**

```python
from dnakit.annotation import midpoint_pybedtools_coordinates

result = midpoint_pybedtools_coordinates("bf_a.bed")
print([(row["start"], row["end"]) for row in result["records"]])
```

- **Example output:**

```text
[(15, 16), (22, 23), (45, 46), (7, 8)]
```

## 9) `dnakit.annotation.tag_pybedtools_bedpe()`

- **Purpose:** Tag each end of a BEDPE record with the named query intervals it overlaps.
- **Method:** Call `pybedtools.contrib.long_range_interaction.tag_bedpe`. Each end is intersected with the BED4 queries, and rows with no overlap are dropped.
- **API:** `dnakit.annotation.tag_pybedtools_bedpe(bedpe: 'str | Path', queries: 'str | Path', *, max_rows: 'int' = 100000, verbose: bool = True)`
- **Input:** A 10-column BEDPE file and a BED4 query file whose column 4 is the name. Requires the BEDTools executable.
- **Output:** A `pairs` list. Each item has `name`, `left`, and `right` (each end's `query` and `overlap_bp`), plus `pair_count` and `provenance`.
- **Example code:**

```python
from dnakit.annotation import tag_pybedtools_bedpe

result = tag_pybedtools_bedpe("pbt_pairs.bedpe", "pbt_q.bed")
print(result["pairs"])
```

- **Example output:**

```text
[{'name': 'p1', 'left': [{'query': 'q1', 'overlap_bp': 6}], 'right': [{'query': 'q2', 'overlap_bp': 4}]}]
```

## 10) `dnakit.annotation.write_pybedtools_chromsizes()`

- **Purpose:** Write a two-column chromosome-size file from a length mapping.
- **Method:** Call `pybedtools.chromsizes_to_file`. Rows are sorted by chromosome name so the file is deterministic. The output must be a new file.
- **API:** `dnakit.annotation.write_pybedtools_chromsizes(chromsizes: 'dict[str, int]', output: 'str | Path', *, verbose: bool = True)`
- **Input:** A mapping from chromosome name to length. `output` is a new file path in an existing directory.
- **Output:** `path`, a `chromosomes` list (chrom, length), `chromosome_count`, and `provenance`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.annotation import write_pybedtools_chromsizes

out = Path(tempfile.mkdtemp()) / "sizes.txt"
result = write_pybedtools_chromsizes({"chr2": 50, "chr1": 100}, out)
print(result["chromosomes"])
```

- **Example output:**

```text
[{'chrom': 'chr1', 'length': 100}, {'chrom': 'chr2', 'length': 50}]
```

## 11) `dnakit.annotation.to_bed12_gffutils_transcript()`

- **Purpose:** Convert one transcript and its child features to a BED12 line.
- **Method:** Fetch the transcript from the database, then call `gffutils.convert.to_bed12`. The BED start is 0-based, and blocks come from the `child_type` features (exon by default).
- **API:** `dnakit.annotation.to_bed12_gffutils_transcript(database: 'str | Path', transcript_id: 'str', *, child_type: 'str' = 'exon', name_field: 'str' = 'ID', verbose: bool = True)`
- **Input:** A gffutils database path that already exists. `transcript_id` is the transcript ID. `child_type` and `name_field` are strings.
- **Output:** A `bed12` tab-separated string without a newline, plus `transcript_id` and `provenance`. An unknown ID raises `ConfigurationError`.
- **Example code:**

```python
import tempfile
from pathlib import Path

import gffutils

db = Path(tempfile.mkdtemp()) / "genes.db"
gffutils.create_db("gff_genes.gff3", str(db), force=True, verbose=False, merge_strategy="create_unique")
from dnakit.annotation import to_bed12_gffutils_transcript

print(to_bed12_gffutils_transcript(db, "tx1")["bed12"])
```

- **Example output:**

```text
chr1	99	500	tx1	.	+	100	500	0,0,0	2	101,201	0,200
```

## 12) `dnakit.annotation.tss_gffutils_intervals()`

- **Purpose:** List the 1-bp transcription start site (TSS) of each transcript in a database.
- **Method:** Call `gffutils.pybedtools_integration.tsses` with BED6 output. The start is 0-based.
- **API:** `dnakit.annotation.tss_gffutils_intervals(database: 'str | Path', *, verbose: bool = True)`
- **Input:** An existing gffutils database path. Requires pybedtools and the BEDTools executable.
- **Output:** A `records` list (chrom, start, end, name, strand), `record_count`, and `provenance`.
- **Example code:**

```python
import tempfile
from pathlib import Path

import gffutils

db = Path(tempfile.mkdtemp()) / "genes.db"
gffutils.create_db("gff_genes.gff3", str(db), force=True, verbose=False, merge_strategy="create_unique")
from dnakit.annotation import tss_gffutils_intervals

print(tss_gffutils_intervals(db)["records"])
```

- **Example output:**

```text
[{'chrom': 'chr1', 'start': 99, 'end': 100, 'name': 'tx1', 'strand': '+'}]
```

## 13) `dnakit.annotation.parse_gffutils_feature_line()`

- **Purpose:** Parse one nine-column GFF3 line into feature fields and attributes.
- **Method:** Call `gffutils.feature.feature_from_line` in strict mode. Coordinates are kept as stored in GFF (1-based, inclusive).
- **API:** `dnakit.annotation.parse_gffutils_feature_line(line: 'str', *, verbose: bool = True)`
- **Input:** One tab-separated nine-column GFF3 line of at most 65536 characters.
- **Output:** `seqid`, `source`, `featuretype`, `start`, `end`, `score`, `strand`, `frame`, and `attributes` (sorted keys, list values), plus `provenance`. A malformed line raises `ConfigurationError`.
- **Example code:**

```python
from dnakit.annotation import parse_gffutils_feature_line

line = "chr1\tsrc\texon\t300\t500\t.\t+\t.\tID=exon2;Parent=tx1"
result = parse_gffutils_feature_line(line)
print(result["featuretype"], result["start"], result["end"], result["attributes"])
```

- **Example output:**

```text
exon 300 500 {'ID': ['exon2'], 'Parent': ['tx1']}
```

## 14) `dnakit.annotation.inspect_gffutils_source()`

- **Purpose:** Count feature types, chromosomes, attribute keys, and features in a GFF or GTF file.
- **Method:** Call `gffutils.inspect.inspect` without logging. Count dictionaries are sorted by key.
- **API:** `dnakit.annotation.inspect_gffutils_source(path: 'str | Path', *, verbose: bool = True)`
- **Input:** A path to a GFF or GTF file.
- **Output:** The `featuretype`, `chrom`, and `attribute_keys` count dictionaries, `feature_count`, and `provenance`.
- **Example code:**

```python
from dnakit.annotation import inspect_gffutils_source

result = inspect_gffutils_source("gff_genes.gff3")
print(result["featuretype"], result["feature_count"])
```

- **Example output:**

```text
{'exon': 2, 'gene': 2, 'mRNA': 1} 5
```

## 15) `dnakit.annotation.parse_gtfparse_fixed_columns()`

- **Purpose:** Read the nine fixed GTF columns and keep the raw attribute text.
- **Method:** Call `gtfparse.parse_gtf` with `split_attributes=False`. Coordinates are 1-based inclusive.
- **API:** `dnakit.annotation.parse_gtfparse_fixed_columns(source: 'str | Path', *, max_rows: 'int' = 100000, verbose: bool = True)`
- **Input:** A GTF file. `max_rows` limits the number of rows.
- **Output:** A `records` list (seqname, source, feature, start, end, score, strand, frame, attribute), `record_count`, and `provenance`.
- **Example code:**

```python
from dnakit.annotation import parse_gtfparse_fixed_columns

result = parse_gtfparse_fixed_columns("gtp_genes.gtf")
print(result["record_count"], result["records"][2]["feature"], result["records"][2]["start"])
```

- **Example output:**

```text
3 exon 100
```

## 16) `dnakit.annotation.expand_gtfparse_attributes()`

- **Purpose:** Read a GTF file and expand its attribute column into one column per key.
- **Method:** Call `gtfparse.parse_gtf_and_expand_attributes`. Absent attribute keys are None.
- **API:** `dnakit.annotation.expand_gtfparse_attributes(source: 'str | Path', *, max_rows: 'int' = 100000, verbose: bool = True)`
- **Input:** A GTF file. `max_rows` limits the number of rows.
- **Output:** The sorted `columns`, `records`, `record_count`, and `provenance`.
- **Example code:**

```python
from dnakit.annotation import expand_gtfparse_attributes

result = expand_gtfparse_attributes("gtp_genes.gtf")
print(result["columns"])
print(result["records"][2]["transcript_id"], result["records"][2]["exon_number"])
```

- **Example output:**

```text
['end', 'exon_number', 'feature', 'frame', 'gene_id', 'gene_name', 'score', 'seqname', 'source', 'start', 'strand', 'transcript_id']
t1 1
```

## 17) `dnakit.annotation.write_gtfparse_table()`

- **Purpose:** Parse a GTF file and write it back out as a new GTF file.
- **Method:** Call `gtfparse.parse_gtf` and `gtfparse.write_gtf`. The `frame` column is written as parsed, so a `.` frame comes back as 0. The output must be a new file.
- **API:** `dnakit.annotation.write_gtfparse_table(source: 'str | Path', output: 'str | Path', *, verbose: bool = True)`
- **Input:** A source GTF file. `output` is a new file path in an existing directory.
- **Output:** `path`, `line_count`, and `provenance`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.annotation import write_gtfparse_table

out = Path(tempfile.mkdtemp()) / "round.gtf"
result = write_gtfparse_table("gtp_genes.gtf", out)
print(result["line_count"])
print(out.read_text(encoding="utf-8").splitlines()[2])
```

- **Example output:**

```text
3
chr1	src	exon	100	200	.	+	0	gene_id "g1"; transcript_id "t1"; exon_number 1;
```

## 18) `dnakit.annotation.read_bigwig_header()`

- **Purpose:** Read chromosome lengths and a header summary from a bigWig file.
- **Method:** Open the file read-only with `pyBigWig.open`, then read `chroms()` and `header()`. Files that are not bigWig are refused.
- **API:** `dnakit.annotation.read_bigwig_header(path: 'str | Path', *, verbose: bool = True)`
- **Input:** A bigWig file path.
- **Output:** A `chromosomes` list (chrom, length, sorted by name), `chromosome_count`, a `header` dictionary, and `provenance`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.annotation import read_bigwig_header, write_bigwig_intervals

out = Path(tempfile.mkdtemp()) / "cov.bw"
write_bigwig_intervals("bw_cov.bedgraph", out, chromsizes={"chr1": 100, "chr2": 50})
result = read_bigwig_header(out)
print(result["chromosomes"])
print(result["header"]["nBasesCovered"])
```

- **Example output:**

```text
[{'chrom': 'chr1', 'length': 100}, {'chrom': 'chr2', 'length': 50}]
24
```

## 19) `dnakit.annotation.write_bigwig_intervals()`

- **Purpose:** Write a new bigWig file from bedGraph rows (chrom, start, end, value).
- **Method:** Open the output with `pyBigWig.open(path, 'w')`, write the header with `addHeader`, and write the intervals with `addEntries`. Intervals must not overlap and must lie inside their chromosome.
- **API:** `dnakit.annotation.write_bigwig_intervals(bedgraph: 'str | Path', output: 'str | Path', *, chromsizes: 'Mapping[str, int]', verbose: bool = True)`
- **Input:** A bedGraph file. `chromsizes` is a length mapping that must cover every chromosome in the file. `output` is a new file path.
- **Output:** `path`, `interval_count`, `chromosomes`, and `provenance`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.annotation import write_bigwig_intervals

out = Path(tempfile.mkdtemp()) / "cov.bw"
result = write_bigwig_intervals("bw_cov.bedgraph", out, chromsizes={"chr1": 100, "chr2": 50})
print(result["interval_count"], result["chromosomes"])
```

- **Example output:**

```text
3 ['chr1', 'chr2']
```

## 20) `dnakit.annotation.summarize_bigwig_bins()`

- **Purpose:** Summarize a bigWig region into equal-width bins.
- **Method:** Call `stats(chrom, start, end, type=statistic, nBins=bins)` on a `pyBigWig` handle. Empty bins are returned as None.
- **API:** `dnakit.annotation.summarize_bigwig_bins(path: 'str | Path', *, chrom: 'str', start: 'int', end: 'int', bins: 'int', statistic: 'str' = 'mean', verbose: bool = True)`
- **Input:** A bigWig path. `chrom`, `start`, and `end` must lie inside the chromosome. `bins` is a positive integer. `statistic` is one of coverage, max, mean, min, std, or sum.
- **Output:** A `values` list (one item per bin), plus `chrom`, `start`, `end`, `bins`, `statistic`, and `provenance`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.annotation import summarize_bigwig_bins, write_bigwig_intervals

out = Path(tempfile.mkdtemp()) / "cov.bw"
write_bigwig_intervals("bw_cov.bedgraph", out, chromsizes={"chr1": 100, "chr2": 50})
result = summarize_bigwig_bins(out, chrom="chr1", start=0, end=20, bins=2)
print(result["values"])
```

- **Example output:**

```text
[1.5, 3.0]
```

## 21) `dnakit.annotation.tile_pyranges1_genome()`

- **Purpose:** Tile each chromosome of a length mapping into fixed-width intervals.
- **Method:** Call `pyranges1.tile_genome`. Each chromosome is tiled from 0. `full_last_tile` is passed to pyranges1 unchanged.
- **API:** `dnakit.annotation.tile_pyranges1_genome(chromsizes: 'Mapping[str, int]', *, tile_size: 'int', full_last_tile: 'bool' = False, max_rows: 'int' = 100000, verbose: bool = True)`
- **Input:** A length mapping, a positive integer `tile_size`, a boolean `full_last_tile`, and `max_rows` for the tile count.
- **Output:** A `tiles` list (chrom, start, end), `tile_count`, `tile_size`, and `provenance`.
- **Example code:**

```python
from dnakit.annotation import tile_pyranges1_genome

result = tile_pyranges1_genome({"chr1": 25, "chr2": 10}, tile_size=10)
print(result["tile_count"])
print([(row["chrom"], row["start"], row["end"]) for row in result["tiles"]])
```

- **Example output:**

```text
4
[('chr1', 0, 10), ('chr1', 10, 20), ('chr1', 20, 25), ('chr2', 0, 10)]
```

## 22) `dnakit.annotation.count_pyranges1_overlaps()`

- **Purpose:** Count the intervals of another BED file that overlap each source interval.
- **Method:** Call `pyranges1.count_overlaps`, with the source intervals as features and the other file as the counted set.
- **API:** `dnakit.annotation.count_pyranges1_overlaps(source: 'str | Path', other: 'str | Path', *, max_rows: 'int' = 100000, verbose: bool = True)`
- **Input:** Two BED files. `max_rows` limits the rows read from each file.
- **Output:** A `records` list (chrom, start, end, strand, and the `overlaps` count), plus `record_count` and `provenance`.
- **Example code:**

```python
from dnakit.annotation import count_pyranges1_overlaps

result = count_pyranges1_overlaps("bf_a.bed", "bf_b.bed")
print([(row["chrom"], row["start"], row["overlaps"]) for row in result["records"]])
```

- **Example output:**

```text
[('chr1', 10, 1), ('chr1', 15, 1), ('chr1', 40, 1), ('chr2', 5, 0)]
```

## 23) `dnakit.annotation.parse_pyranges1_string()`

- **Purpose:** Parse a headed multi-line table (Chromosome, Start, End, optional Strand) into intervals.
- **Method:** Call `pyranges1.from_string`. The first line is the header and must contain Chromosome, Start, and End. Strand is optional.
- **API:** `dnakit.annotation.parse_pyranges1_string(text: 'str', *, verbose: bool = True)`
- **Input:** A text block of at most 65536 characters.
- **Output:** A `records` list (chrom, start, end, and strand when present), `record_count`, and `provenance`. Text that cannot be parsed raises `ConfigurationError`.
- **Example code:**

```python
from dnakit.annotation import parse_pyranges1_string

text = "Chromosome Start End Strand\nchr1 10 20 +\nchr1 40 50 -\n"
result = parse_pyranges1_string(text)
print(result["records"])
```

- **Example output:**

```text
[{'chrom': 'chr1', 'start': 10, 'end': 20, 'strand': '+'}, {'chrom': 'chr1', 'start': 40, 'end': 50, 'strand': '-'}]
```

## 24) `dnakit.annotation.jaccard_pyranges1_intervals()`

- **Purpose:** Compute the base-pair Jaccard index of two BED files.
- **Method:** Call `pyranges1.stats.jaccard`. The numerator is the number of bases shared by both files. The denominator is the number of bases covered by either file.
- **API:** `dnakit.annotation.jaccard_pyranges1_intervals(source: 'str | Path', other: 'str | Path', *, strand_behavior: 'str' = 'auto', max_rows: 'int' = 100000, verbose: bool = True)`
- **Input:** Two BED files. `strand_behavior` is auto, same, opposite, or ignore. `max_rows` limits the rows read from each file.
- **Output:** `jaccard` (a float), `strand_behavior`, and `provenance`.
- **Example code:**

```python
from dnakit.annotation import jaccard_pyranges1_intervals

result = jaccard_pyranges1_intervals("bf_a.bed", "bf_b.bed")
print(round(result["jaccard"], 4))
```

- **Example output:**

```text
0.125
```
