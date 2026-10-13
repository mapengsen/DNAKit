# 93 G18g: interval and annotation functions from pybedtools, gffutils and gtfparse

This page describes the 22 functions wrapped for group G18g. They are the DNA-level public functions of pybedtools, gffutils and gtfparse that run offline on local inputs. Network downloads, plotting and UCSC binary tools are not wrapped. Every example output was produced by running the example. The wrappers live in the `pybedtools_review_functions`, `gffutils_review_functions` and `gtfparse_review_functions` modules of `dnakit.annotation`, and the examples import from those modules. The modules are not yet re-exported by the `__init__` of `dnakit.annotation`, so the API lines give the full module path.

## Coverage

| Package | Status | Count | other and cap before this pass |
| --- | --- | ---: | ---: |
| pybedtools | internal | 231 | 62 |
| pybedtools | wrapped | 13 |  |
| pybedtools | needs-binary | 7 |  |
| pybedtools | plot | 7 |  |
| pybedtools | needs-network | 5 |  |
| pybedtools | duplicate-of | 4 |  |
| pybedtools | family-variant | 4 |  |
| pybedtools | deprecated | 1 |  |
| pybedtools | not-dna | 1 |  |
| gffutils | internal | 103 | 35 |
| gffutils | wrapped | 13 |  |
| gffutils | family-variant | 10 |  |
| gffutils | duplicate-of | 5 |  |
| gffutils | already-wrapped | 1 |  |
| gffutils | cli | 1 |  |
| gtfparse | duplicate-of | 7 | 5 |
| gtfparse | wrapped | 5 |  |
| gtfparse | internal | 3 |  |
| gtfparse | family-variant | 2 |  |
| gtfparse | already-wrapped | 1 |  |

Status notes: `wrapped` marks the functions on this page; `duplicate-of` points to the wrapped name of the same upstream function; `family-variant` is chosen by an argument of a wrapped function; `internal` covers state, global configuration, temporary files and constructor-level helpers; `needs-binary` needs samtools, which is not installed here; `plot` is for plotting or display; `not-dna` is a network reachability check. Before this pass the packages had 102 `other` and `cap` rows (pybedtools 62, gffutils 35, gtfparse 5); all of them now have a new status.

## Wrapped functions and aliases

- `pybedtools.cbedtools.overlap` and `pybedtools.overlap` are the same upstream function; the alias row points to `overlap_pybedtools_bp`.
- `gffutils.DataIterator` and `gffutils.iterators.DataIterator` are the same function; the alias row points to `read_gffutils_records`.
- The nine gffutils `merge_criteria` functions and `overlap_any_threshold` are selected by the `criteria` argument of `merge_gffutils_features`, so they have status `family-variant`.
- The gtfparse module-level aliases point to the wrapped names of the same upstream functions.

## 1) `overlap_pybedtools_bp`

- **Purpose:** Signed overlap length of two 0-based half-open intervals, and whether they overlap.
- **Method:** Calls pybedtools.overlap. That upstream function takes both starts before both ends; the wrapper maps its arguments accordingly. A negative value is the gap between the intervals, and touching intervals give 0.
- **API:** `dnakit.annotation.overlap_pybedtools_bp(start1, end1, start2, end2, *, verbose=True)`
- **Input:** start1, end1, start2 and end2 are non-negative integers with start < end.
- **Output:** overlap_bp (signed integer) and overlaps (boolean).

- **Example code:**

```python
from dnakit.annotation import overlap_pybedtools_bp

for start1, end1, start2, end2 in [(0, 10, 5, 15), (0, 10, 10, 20), (0, 10, 20, 30)]:
    result = overlap_pybedtools_bp(start1, end1, start2, end2, verbose=False)
    print(result["overlap_bp"], result["overlaps"])
```

- **Example output:**

```text
5 True
0 False
-10 False
```

## 2) `normalize_pybedtools_scores_to_length`

- **Purpose:** Divide the numeric value in a BED column by the interval length and multiply by a scalar, for example to turn read counts into a per-kb rate.
- **Method:** Calls featurefuncs.normalized_to_length for each interval: value × scalar ÷ (end − start). The defaults index=4 and scalar=0.001 give a per-kb rate.
- **API:** `dnakit.annotation.normalize_pybedtools_scores_to_length(source, *, index=4, scalar=0.001, max_rows=100000, verbose=True)`
- **Input:** A BED3 to BED6 file; the column at index (0-based) must hold finite numbers.
- **Output:** records (chrom, start, end, name, score, strand), sorted by coordinate.

- **Example file:** `bg.bedgraph`

```text
chr1	0	100	10
chr1	100	200	20
```

- **Example code:**

```python
from dnakit.annotation import normalize_pybedtools_scores_to_length

result = normalize_pybedtools_scores_to_length("bg.bedgraph", index=3, scalar=1000, verbose=False)
for record in result["records"]:
    print(record["chrom"], record["start"], record["end"], record["name"])
```

- **Example output:**

```text
chr1 0 100 100.0
chr1 100 200 200.0
```

## 3) `scale_pybedtools_bedgraph_values`

- **Purpose:** Multiply the value column of a bedGraph by a scalar, for example to normalise or convert signal-track units.
- **Method:** Calls featurefuncs.bedgraph_scale for each interval, multiplying column 4 by scalar and writing the float back as text.
- **API:** `dnakit.annotation.scale_pybedtools_bedgraph_values(source, *, scalar, max_rows=100000, verbose=True)`
- **Input:** A bedGraph with exactly four columns; column 4 must hold finite numbers.
- **Output:** records; the scaled value is returned in name, sorted by coordinate.

- **Example file:** `bg.bedgraph`

```text
chr1	0	100	10
chr1	100	200	20
```

- **Example code:**

```python
from dnakit.annotation import scale_pybedtools_bedgraph_values

result = scale_pybedtools_bedgraph_values("bg.bedgraph", scalar=0.5, verbose=False)
for record in result["records"]:
    print(record["chrom"], record["start"], record["end"], record["name"])
```

- **Example output:**

```text
chr1 0 100 5.0
chr1 100 200 10.0
```

## 4) `center_pybedtools_coordinates`

- **Purpose:** Take a fixed-width window centred on each interval's midpoint, for example to standardise peak centres.
- **Method:** Calls featurefuncs.center. The window is 2 × (width // 2) bp; intervals shorter than width are unchanged; width=1 gives the base at the midpoint. Upstream clamps a window start to 1 rather than 0, and the wrapper keeps that behaviour.
- **API:** `dnakit.annotation.center_pybedtools_coordinates(source, *, width=100, max_rows=100000, verbose=True)`
- **Input:** A BED3 to BED6 file; width is an integer from 1 to 1e9.
- **Output:** records and width, sorted by coordinate.

- **Example file:** `peaks.bed`

```text
chr1	100	200	peakA	5	+
chr1	300	340	peakB	9	-
chr2	50	60	peakC	2	+
```

- **Example code:**

```python
from dnakit.annotation import center_pybedtools_coordinates

result = center_pybedtools_coordinates("peaks.bed", width=40, verbose=False)
for record in result["records"]:
    print(record["name"], record["start"], record["end"])
```

- **Example output:**

```text
peakA 130 170
peakB 300 340
peakC 50 60
```

## 5) `three_prime_pybedtools_coordinates`

- **Purpose:** Take a window around each interval's 3' end, for example around polyA sites or 3' UTR boundaries.
- **Method:** Calls featurefuncs.three_prime. The 3' end is stop on + and start on −. upstream counts bases on the feature side of that end, and downstream counts bases beyond it.
- **API:** `dnakit.annotation.three_prime_pybedtools_coordinates(source, *, upstream=500, downstream=500, max_rows=100000, verbose=True)`
- **Input:** A BED3 to BED6 file; upstream and downstream are non-negative integers.
- **Output:** records, upstream and downstream, sorted by coordinate.

- **Example file:** `peaks.bed`

```text
chr1	100	200	peakA	5	+
chr1	300	340	peakB	9	-
chr2	50	60	peakC	2	+
```

- **Example code:**

```python
from dnakit.annotation import three_prime_pybedtools_coordinates

result = three_prime_pybedtools_coordinates("peaks.bed", upstream=10, downstream=20, verbose=False)
for record in result["records"]:
    print(record["name"], record["start"], record["end"])
```

- **Example output:**

```text
peakA 190 220
peakB 280 310
peakC 50 80
```

## 6) `filter_pybedtools_intervals_by_length`

- **Purpose:** Keep BED intervals by length, with a lower bound, an upper bound, or both.
- **Method:** Chains BedTool.filter with featurefuncs.greater_than (length > bound) and featurefuncs.less_than (length < bound). Both thresholds are strict.
- **API:** `dnakit.annotation.filter_pybedtools_intervals_by_length(source, *, greater_than=None, less_than=None, max_rows=100000, verbose=True)`
- **Input:** A BED3 to BED6 file; give at least one of greater_than and less_than.
- **Output:** records, greater_than and less_than, sorted by coordinate.

- **Example file:** `peaks.bed`

```text
chr1	100	200	peakA	5	+
chr1	300	340	peakB	9	-
chr2	50	60	peakC	2	+
```

- **Example code:**

```python
from dnakit.annotation import filter_pybedtools_intervals_by_length

result = filter_pybedtools_intervals_by_length("peaks.bed", greater_than=30, less_than=100, verbose=False)
for record in result["records"]:
    print(record["name"], record["end"] - record["start"])
```

- **Example output:**

```text
peakB 40
```

## 7) `gff_to_bed_pybedtools`

- **Purpose:** Convert GFF features to BED6 so that they can be combined with other interval tools.
- **Method:** Calls featurefuncs.gff2bed on each line. name_field picks an attribute name or a column index for the name; by default the ID is used, and a missing attribute gives '.'.
- **API:** `dnakit.annotation.gff_to_bed_pybedtools(source, *, name_field=None, max_rows=100000, verbose=True)`
- **Input:** A nine-column GFF file; start and end are 1-based integers with start ≤ end.
- **Output:** records with BED6 fields; start is 0-based.

- **Example file:** `g.gff`

```text
chr1	src	gene	101	200	.	+	.	ID=g1;Name=GENE1
chr1	src	mRNA	101	180	5	+	.	ID=t1;Parent=g1
```

- **Example code:**

```python
from dnakit.annotation import gff_to_bed_pybedtools

result = gff_to_bed_pybedtools("g.gff", name_field="Name", verbose=False)
for record in result["records"]:
    print(record["chrom"], record["start"], record["end"], record["name"], record["score"], record["strand"])
```

- **Example output:**

```text
chr1 100 180 . 5 +
chr1 100 200 GENE1 . +
```

## 8) `bed_to_gff_pybedtools`

- **Purpose:** Write BED intervals as GFF records for downstream tools that need GFF.
- **Method:** Calls featurefuncs.bed2gff on each line: the start gains 1 to become 1-based, the name goes in the Name attribute, and source and feature are '.'.
- **API:** `dnakit.annotation.bed_to_gff_pybedtools(source, *, max_rows=100000, verbose=True)`
- **Input:** A BED3 to BED6 file.
- **Output:** records with the nine GFF fields, sorted by coordinate.

- **Example file:** `peaks.bed`

```text
chr1	100	200	peakA	5	+
chr1	300	340	peakB	9	-
chr2	50	60	peakC	2	+
```

- **Example code:**

```python
from dnakit.annotation import bed_to_gff_pybedtools

result = bed_to_gff_pybedtools("peaks.bed", verbose=False)
for record in result["records"]:
    print(record["seqid"], record["start"], record["end"], record["strand"], record["attributes"])
```

- **Example output:**

```text
chr1 101 200 + Name="peakA";
chr1 301 340 - Name="peakB";
chr2 51 60 + Name="peakC";
```

## 9) `cis_trans_pybedtools_interactions`

- **Purpose:** For chromatin interactions in BEDPE, find the named intervals hit by each end and their partners on the same end (cis) and on the other end (distal).
- **Method:** Calls tag_bedpe to intersect both BEDPE ends with named query sets, then cis_trans_interactions to summarise. At least two query sets are needed, because bedtools omits the name column when given one -b file.
- **API:** `dnakit.annotation.cis_trans_pybedtools_interactions(bedpe, queries, *, max_rows=100000, verbose=True)`
- **Input:** A ten-column BEDPE file; queries maps two to eight labels to BED4 files. Labels may use letters, digits, underscore, dot and hyphen.
- **Output:** records (target, cis, distal and label fields; '.' for a missing partner), sorted by field.

- **Example file:** `pairs.bedpe`

```text
chr1	1	10	chr1	50	90	pair1	5	+	-
chr1	2	15	chr1	200	210	pair2	1	+	+
```

- **Example file:** `tss.bed`

```text
chr1	5	6	gene1
chr1	60	61	gene2
chr1	88	89	gene3
```

- **Example file:** `pk.bed`

```text
chr1	3	4	peak1	50	.
```

- **Example code:**

```python
from dnakit.annotation import cis_trans_pybedtools_interactions

result = cis_trans_pybedtools_interactions(
    "pairs.bedpe", {"tss": "tss.bed", "pk": "pk.bed"}, verbose=False
)
print(result["record_count"])
for row in result["records"]:
    if row["label"] == "pair1":
        print(row["target_label"], row["target_name"], row["cis_name"], row["distal_name"])
```

- **Example output:**

```text
10
tss gene1 peak1 gene2
tss gene1 peak1 gene3
pk peak1 gene1 gene2
pk peak1 gene1 gene3
tss gene3 gene2 peak1
tss gene3 gene2 gene1
tss gene2 gene3 peak1
tss gene2 gene3 gene1
```

## 10) `random_jaccard_pybedtools_intervals`

- **Purpose:** Compare the Jaccard index of two BED files with one seeded random shuffle as a reference value.
- **Method:** Shuffles a once on the genome with bedtools shuffle (fixed -seed) and computes its Jaccard with b; also reports the unshuffled observed value. This is a single draw from the null model, not a p-value.
- **API:** `dnakit.annotation.random_jaccard_pybedtools_intervals(a, b, genome, *, seed=0, max_rows=100000, verbose=True)`
- **Input:** Two BED3 to BED6 files and a two-column chromosome-size file; seed is a non-negative integer.
- **Output:** observed and shuffled (intersection_bp, union_bp, jaccard, n_intersections).

- **Example file:** `a.bed`

```text
chr1	100	300	peakA
```

- **Example file:** `b.bed`

```text
chr1	150	350	b1
chr1	600	700	b2
```

- **Example file:** `genome.tsv`

```text
chr1	1000
chr2	500
```

- **Example code:**

```python
from dnakit.annotation import random_jaccard_pybedtools_intervals

result = random_jaccard_pybedtools_intervals("a.bed", "b.bed", "genome.tsv", seed=7, verbose=False)
print(result["observed"])
print(result["shuffled"])
```

- **Example output:**

```text
{'intersection_bp': 150, 'union_bp': 350, 'jaccard': 0.428571, 'n_intersections': 1}
{'intersection_bp': 12, 'union_bp': 488, 'jaccard': 0.0245902, 'n_intersections': 1}
```

## 11) `read_gffutils_records`

- **Purpose:** Read GFF3 or GTF features one by one, keeping coordinates, type and attributes.
- **Method:** Calls gffutils.DataIterator and keeps file order; attribute values become sorted lists; max_rows limits the number of records.
- **API:** `dnakit.annotation.read_gffutils_records(source, *, max_rows=100000, verbose=True)`
- **Input:** A GFF3 or GTF file.
- **Output:** records (seqid, source, featuretype, start, end, score, strand, frame, attributes).

- **Example file:** `small.gff3`

```text
##gff-version 3
chr1	ens	gene	101	400	.	+	.	ID=g1;Name=GENE1;gene_biotype=protein_coding
chr1	ens	mRNA	101	400	.	+	.	ID=t1;Parent=g1;Name=T1
chr1	ens	exon	101	200	.	+	.	ID=e1;Parent=t1
chr1	ens	exon	301	400	.	+	.	ID=e2;Parent=t1
chr1	ens	CDS	151	200	.	+	0	ID=c1;Parent=t1
chr1	ens	gene	501	700	.	-	.	ID=g2;Name=GENE2
chr1	ens	mRNA	501	700	.	-	.	ID=t2;Parent=g2
chr1	ens	exon	501	700	.	-	.	ID=e3;Parent=t2
chr2	ens	gene	11	40	.	+	.	ID=g3;Name=GENE3
chr2	ens	mRNA	11	40	.	+	.	ID=t3;Parent=g3
chr2	ens	exon	11	40	.	+	.	ID=e4;Parent=t3
```

- **Example code:**

```python
from dnakit.annotation import read_gffutils_records

result = read_gffutils_records("small.gff3", verbose=False)
print(result["record_count"])
first = result["records"][0]
print(first["featuretype"], first["start"], first["end"], first["attributes"]["Name"])
```

- **Example output:**

```text
11
gene 101 400 ['GENE1']
```

## 12) `bin_gffutils_interval`

- **Purpose:** Compute UCSC genomic bin numbers, for fast interval lookup of features.
- **Method:** Calls gffutils.bins.bins. one=True returns the smallest bin that holds the whole interval; one=False returns every bin that overlaps it.
- **API:** `dnakit.annotation.bin_gffutils_interval(start, stop, *, fmt='gff', one=True, verbose=True)`
- **Input:** Integer start and stop; with fmt='gff' start is 1-based and at least 1, and with fmt='bed' start is 0-based.
- **Output:** bins list in ascending order.

- **Example code:**

```python
from dnakit.annotation import bin_gffutils_interval

print(bin_gffutils_interval(100, 200, verbose=False)["bins"])
print(bin_gffutils_interval(1, 70000, one=False, verbose=False)["bins"])
```

- **Example output:**

```text
[4681]
[1, 9, 73, 585, 4681]
```

## 13) `asinterval_gffutils_feature`

- **Purpose:** Turn one database feature into the coordinates and name of a pybedtools interval.
- **Method:** Looks up the feature by ID (a GFF or GTF file is first loaded into an in-memory database) and calls gffutils.helpers.asinterval.
- **API:** `dnakit.annotation.asinterval_gffutils_feature(database, feature_id, *, verbose=True)`
- **Input:** A GFF3 or GTF file or a gffutils .db database, plus the feature ID.
- **Output:** chrom, start (0-based), end, name, score, strand, featuretype and fields.

- **Example file:** `small.gff3`

```text
##gff-version 3
chr1	ens	gene	101	400	.	+	.	ID=g1;Name=GENE1;gene_biotype=protein_coding
chr1	ens	mRNA	101	400	.	+	.	ID=t1;Parent=g1;Name=T1
chr1	ens	exon	101	200	.	+	.	ID=e1;Parent=t1
chr1	ens	exon	301	400	.	+	.	ID=e2;Parent=t1
chr1	ens	CDS	151	200	.	+	0	ID=c1;Parent=t1
chr1	ens	gene	501	700	.	-	.	ID=g2;Name=GENE2
chr1	ens	mRNA	501	700	.	-	.	ID=t2;Parent=g2
chr1	ens	exon	501	700	.	-	.	ID=e3;Parent=t2
chr2	ens	gene	11	40	.	+	.	ID=g3;Name=GENE3
chr2	ens	mRNA	11	40	.	+	.	ID=t3;Parent=g3
chr2	ens	exon	11	40	.	+	.	ID=e4;Parent=t3
```

- **Example code:**

```python
from dnakit.annotation import asinterval_gffutils_feature

result = asinterval_gffutils_feature("small.gff3", "g1", verbose=False)
print(result["chrom"], result["start"], result["end"], result["name"], result["strand"])
```

- **Example output:**

```text
chr1 100 400 g1 +
```

## 14) `features_gffutils_to_bedtool`

- **Purpose:** Convert database features, optionally of one type, into intervals for interval arithmetic.
- **Method:** Calls gffutils.pybedtools_integration.to_bedtool on features read in seqid and start order.
- **API:** `dnakit.annotation.features_gffutils_to_bedtool(database, *, featuretype=None, max_rows=100000, verbose=True)`
- **Input:** A GFF3 or GTF file or a .db database; featuretype is optional.
- **Output:** records (chrom, start, end, name, score, strand, featuretype), sorted by coordinate.

- **Example file:** `small.gff3`

```text
##gff-version 3
chr1	ens	gene	101	400	.	+	.	ID=g1;Name=GENE1;gene_biotype=protein_coding
chr1	ens	mRNA	101	400	.	+	.	ID=t1;Parent=g1;Name=T1
chr1	ens	exon	101	200	.	+	.	ID=e1;Parent=t1
chr1	ens	exon	301	400	.	+	.	ID=e2;Parent=t1
chr1	ens	CDS	151	200	.	+	0	ID=c1;Parent=t1
chr1	ens	gene	501	700	.	-	.	ID=g2;Name=GENE2
chr1	ens	mRNA	501	700	.	-	.	ID=t2;Parent=g2
chr1	ens	exon	501	700	.	-	.	ID=e3;Parent=t2
chr2	ens	gene	11	40	.	+	.	ID=g3;Name=GENE3
chr2	ens	mRNA	11	40	.	+	.	ID=t3;Parent=g3
chr2	ens	exon	11	40	.	+	.	ID=e4;Parent=t3
```

- **Example code:**

```python
from dnakit.annotation import features_gffutils_to_bedtool

result = features_gffutils_to_bedtool("small.gff3", featuretype="gene", verbose=False)
for record in result["records"]:
    print(record["chrom"], record["start"], record["end"], record["name"])
```

- **Example output:**

```text
chr1 100 400 g1
chr1 500 700 g2
chr2 10 40 g3
```

## 15) `infer_gffutils_attribute_dialect`

- **Purpose:** Infer whether attribute text is GFF3 or GTF, with its separators, quoting and key order.
- **Method:** Calls gffutils.helpers.infer_dialect.
- **API:** `dnakit.annotation.infer_gffutils_attribute_dialect(attributes, *, verbose=True)`
- **Input:** The attribute column text, non-empty and at most 65536 characters.
- **Output:** fmt (gff3 or gtf) and the dialect dictionary; sets become sorted lists.

- **Example code:**

```python
from dnakit.annotation import infer_gffutils_attribute_dialect

print(infer_gffutils_attribute_dialect("ID=g1;Name=x", verbose=False)["fmt"])
print(infer_gffutils_attribute_dialect('gene_id "A"; transcript_id "B";', verbose=False)["dialect"]["keyval separator"])
```

- **Example output:**

```text
gff3
```

## 16) `canonical_gffutils_transcript_sequences`

- **Purpose:** Pick one canonical transcript per gene and extract its exon sequence.
- **Method:** Calls gffutils.helpers.canonical_transcripts: the transcript with the longest CDS wins, or the longest transcript when there is no CDS. Sequences come from pyfaidx, reverse-complemented on the − strand. Upstream prints each choice; the wrapper hides that output. pyfaidx writes a .fai index next to the FASTA.
- **API:** `dnakit.annotation.canonical_gffutils_transcript_sequences(database, fasta, *, max_rows=100000, verbose=True)`
- **Input:** A GFF3 or GTF file (or .db) and a local FASTA covering every exon; every gene needs at least one transcript.
- **Output:** records (transcript_id, gene_ids, seqid, start, end, strand, length, sequence), sorted by transcript ID.

- **Example file:** `small.gff3`

```text
##gff-version 3
chr1	ens	gene	101	400	.	+	.	ID=g1;Name=GENE1;gene_biotype=protein_coding
chr1	ens	mRNA	101	400	.	+	.	ID=t1;Parent=g1;Name=T1
chr1	ens	exon	101	200	.	+	.	ID=e1;Parent=t1
chr1	ens	exon	301	400	.	+	.	ID=e2;Parent=t1
chr1	ens	CDS	151	200	.	+	0	ID=c1;Parent=t1
chr1	ens	gene	501	700	.	-	.	ID=g2;Name=GENE2
chr1	ens	mRNA	501	700	.	-	.	ID=t2;Parent=g2
chr1	ens	exon	501	700	.	-	.	ID=e3;Parent=t2
chr2	ens	gene	11	40	.	+	.	ID=g3;Name=GENE3
chr2	ens	mRNA	11	40	.	+	.	ID=t3;Parent=g3
chr2	ens	exon	11	40	.	+	.	ID=e4;Parent=t3
```

- **Example code:**

```python
from pathlib import Path
from dnakit.annotation import canonical_gffutils_transcript_sequences

pattern = ("ACGT" * 300)[:1000]
Path("ref.fa").write_text(f">chr1\n{pattern}\n>chr2\n{'G' * 100}\n", encoding="utf-8")
result = canonical_gffutils_transcript_sequences("small.gff3", "ref.fa", verbose=False)
for record in result["records"]:
    print(record["transcript_id"], record["gene_ids"], record["length"], record["sequence"][:12])
```

- **Example output:**

```text
t1 ['g1'] 50 GTACGTACGTAC
t2 ['g2'] 200 ACGTACGTACGT
t3 ['g3'] 30 GGGGGGGGGGGG
```

## 17) `gffutils_line_to_seqfeature`

- **Purpose:** Convert one GFF line to a Biopython SeqFeature for use with Biopython code.
- **Method:** Calls gffutils.biopython_integration.to_seqfeature and returns type, id, location and qualifiers.
- **API:** `dnakit.annotation.gffutils_line_to_seqfeature(line, *, verbose=True)`
- **Input:** One tab-separated GFF line with nine columns.
- **Output:** type, id, location (0-based start, exclusive end, strand) and qualifiers.

- **Example code:**

```python
from dnakit.annotation import gffutils_line_to_seqfeature

line = "chr1\tens\tgene\t101\t400\t.\t+\t.\tID=g1;Name=GENE1"
result = gffutils_line_to_seqfeature(line, verbose=False)
print(result["type"], result["location"], result["qualifiers"]["Name"])
```

- **Example output:**

```text
gene {'start': 100, 'end': 400, 'strand': 1} ['GENE1']
```

## 18) `seqfeature_to_gffutils_line`

- **Purpose:** Write an interval description as one GFF line; this is the reverse of the Biopython conversion.
- **Method:** Builds a Bio.SeqFeature and calls gffutils.biopython_integration.from_seqfeature to get a Feature, then writes its GFF line.
- **API:** `dnakit.annotation.seqfeature_to_gffutils_line(seqid, feature_type, start, end, *, strand='+', feature_id=None, source='.', score='.', frame='.', attributes=None, verbose=True)`
- **Input:** seqid, feature_type, a 0-based start, an inclusive end and a strand; optional feature_id, source, score, frame and an attribute mapping. Attribute names may not be GFF column names.
- **Output:** line (one GFF line), seqid, featuretype, start (1-based) and end.

- **Example code:**

```python
from dnakit.annotation import seqfeature_to_gffutils_line

result = seqfeature_to_gffutils_line(
    "chr1", "gene", 99, 400, feature_id="g9", attributes={"Name": "Z"}, verbose=False
)
print(result["line"])
```

- **Example output:**

```text
chr1	.	gene	100	400	.	+	.	ID=g9;Name=Z
```

## 19) `sanitize_gffutils_file`

- **Purpose:** Write a sanitised GFF copy: fix reversed coordinates, add a gid to each record, and drop features outside genes.
- **Method:** Copies the input to a new path, then calls gffutils.helpers.sanitize_gff_file with in_place, which rewrites only the copy. The output header carries a timestamp; the source file is not modified.
- **API:** `dnakit.annotation.sanitize_gffutils_file(source, output, *, verbose=True)`
- **Input:** A GFF3 file and an output path that does not exist yet (its directory must exist).
- **Output:** path, feature_count and gene_count.

- **Example file:** `small.gff3`

```text
##gff-version 3
chr1	ens	gene	101	400	.	+	.	ID=g1;Name=GENE1;gene_biotype=protein_coding
chr1	ens	mRNA	101	400	.	+	.	ID=t1;Parent=g1;Name=T1
chr1	ens	exon	101	200	.	+	.	ID=e1;Parent=t1
chr1	ens	exon	301	400	.	+	.	ID=e2;Parent=t1
chr1	ens	CDS	151	200	.	+	0	ID=c1;Parent=t1
chr1	ens	gene	501	700	.	-	.	ID=g2;Name=GENE2
chr1	ens	mRNA	501	700	.	-	.	ID=t2;Parent=g2
chr1	ens	exon	501	700	.	-	.	ID=e3;Parent=t2
chr2	ens	gene	11	40	.	+	.	ID=g3;Name=GENE3
chr2	ens	mRNA	11	40	.	+	.	ID=t3;Parent=g3
chr2	ens	exon	11	40	.	+	.	ID=e4;Parent=t3
```

- **Example code:**

```python
import os
from dnakit.annotation import sanitize_gffutils_file

if os.path.exists("sanitised.gff3"):
    os.remove("sanitised.gff3")
result = sanitize_gffutils_file("small.gff3", "sanitised.gff3", verbose=False)
print(result["feature_count"], result["gene_count"])
data = [line for line in open("sanitised.gff3", encoding="utf-8") if not line.startswith("#")]
print(data[0].rstrip("\n"))
```

- **Example output:**

```text
11 3
chr1	ens	gene	101	400	.	+	.	ID=g1;Name=GENE1;gene_biotype=protein_coding;gid=g1
```

## 20) `merge_gffutils_features`

- **Purpose:** Merge features of one type whose conditions hold, for example joining adjacent exons of one transcript.
- **Method:** Sorts by seqid and start, then calls FeatureDB.merge. criteria picks upstream merge_criteria; a feature joins a group only when every criterion holds, and threshold feeds the threshold criteria.
- **API:** `dnakit.annotation.merge_gffutils_features(source, *, featuretype, criteria=('seqid', 'overlap_end_inclusive', 'strand', 'feature_type'), threshold=0, max_rows=100000, verbose=True)`
- **Input:** A GFF3 or GTF file or .db, a featuretype, and optional criteria and threshold.
- **Output:** merged_count and merged (id, seqid, start, end, strand, member_ids).

- **Example file:** `small.gff3`

```text
##gff-version 3
chr1	ens	gene	101	400	.	+	.	ID=g1;Name=GENE1;gene_biotype=protein_coding
chr1	ens	mRNA	101	400	.	+	.	ID=t1;Parent=g1;Name=T1
chr1	ens	exon	101	200	.	+	.	ID=e1;Parent=t1
chr1	ens	exon	301	400	.	+	.	ID=e2;Parent=t1
chr1	ens	CDS	151	200	.	+	0	ID=c1;Parent=t1
chr1	ens	gene	501	700	.	-	.	ID=g2;Name=GENE2
chr1	ens	mRNA	501	700	.	-	.	ID=t2;Parent=g2
chr1	ens	exon	501	700	.	-	.	ID=e3;Parent=t2
chr2	ens	gene	11	40	.	+	.	ID=g3;Name=GENE3
chr2	ens	mRNA	11	40	.	+	.	ID=t3;Parent=g3
chr2	ens	exon	11	40	.	+	.	ID=e4;Parent=t3
```

- **Example code:**

```python
from dnakit.annotation import merge_gffutils_features

result = merge_gffutils_features(
    "small.gff3", featuretype="exon", criteria=["seqid", "strand", "feature_type"], verbose=False
)
print(result["merged_count"])
for record in result["merged"]:
    print(record["seqid"], record["start"], record["end"], record["member_ids"])
```

- **Example output:**

```text
3
chr1 101 400 ['e1', 'e2']
chr1 501 700 []
chr2 11 40 []
```

## 21) `expand_gtfparse_attribute_strings`

- **Purpose:** Expand the attribute text of a GTF ninth column into one column per key, for table-style analysis.
- **Method:** Calls gtfparse.attribute_parsing.expand_attribute_strings. Each input row gives one list entry; a repeated key in one row is joined with commas; an absent key is an empty string; a None row has no attributes.
- **API:** `dnakit.annotation.expand_gtfparse_attribute_strings(attribute_strings, *, quote_char="'", usecols=None, max_rows=100000, verbose=True)`
- **Input:** A list of attribute strings (or None), an optional quote character, and optional keys to keep.
- **Output:** row_count, column_count and columns (keys sorted).

- **Example code:**

```python
from dnakit.annotation import expand_gtfparse_attribute_strings

rows = ['gene_id "G1"; transcript_id "T1"; gene_name "A";', None, 'gene_id "G2"; tag "x"; tag "y";']
result = expand_gtfparse_attribute_strings(rows, verbose=False)
for key, values in result["columns"].items():
    print(key, values)
```

- **Example output:**

```text
gene_id ['G1', '', 'G2']
gene_name ['A', '', '']
tag ['', '', 'x,y']
transcript_id ['T1', '', '']
```

## 22) `create_missing_gtfparse_features`

- **Purpose:** Create missing transcript or gene rows when a GTF has only exon rows.
- **Method:** Calls gtfparse.parse_gtf_and_expand_attributes, then create_missing_features grouped by the key (for example transcript_id); start is the minimum and end the maximum of member rows. A feature type that already has rows gets no new rows.
- **API:** `dnakit.annotation.create_missing_gtfparse_features(source, *, feature='transcript', key='transcript_id', extra_columns=(), max_rows=100000, verbose=True)`
- **Input:** A GTF file, feature (transcript or gene), the key column and optional extra columns.
- **Output:** feature, key, created_count and records.

- **Example file:** `exons.gtf`

```text
chr1	HAVANA	exon	101	200	.	+	.	gene_id "G1"; transcript_id "T1"; gene_name "GENEA";
chr1	HAVANA	exon	301	400	.	+	.	gene_id "G1"; transcript_id "T1"; gene_name "GENEA";
chr1	HAVANA	exon	501	600	.	-	.	gene_id "G2"; transcript_id "T2"; gene_name "GENEB";
```

- **Example code:**

```python
from dnakit.annotation import create_missing_gtfparse_features

result = create_missing_gtfparse_features(
    "exons.gtf", feature="transcript", key="transcript_id", extra_columns=["gene_name"], verbose=False
)
for record in result["records"]:
    print(record["feature"], record["transcript_id"], record["start"], record["end"], record["strand"])
```

- **Example output:**

```text
transcript T1 101 400 +
transcript T2 501 600 -
```
