# G18a functions: bioframe intervals, regions, sequences, files and assemblies

This page documents the G18a wrappers in the `annotation` area. Each function lives in `dnakit.annotation` and can be imported directly, for example `from dnakit.annotation import validate_bioframe_bed`.
- All 36 wrappers in this group read only local input (BED, FASTA, SAM, tabix, bigBed, chrom.sizes and in-memory mappings) and need no network. UCSC downloads, remote centromere tables and UCSC mRNA fetches are outside this group.
- Example inputs are listed under the Example input files section below. The examples build `tab.bed.gz` and its `.tbi` with pysam; the bigBed example uses the `test.bigBed` file shipped with the pyBigWig test package.
- Coordinates are 0-based half-open unless stated. Region strings are read as written with no coordinate shift, matching the existing `parse_bioframe_region`.
- `digest_bioframe_fasta` corrects the cut positions reported by bioframe.digest; see that section.

## Existing and unwrapped functions

- This group handled the 153 `other` and `cap` ledger rows: 36 new wrappers; 62 aliases of the same upstream object marked `duplicate-of`; 44 internal helpers and constructors; 2 non-DNA colour formatter; 3 rows that need the pairix bindings; 6 family variants (`select_indices`, `select_labels`, `select_mask`, which use the same overlap test as `select_bioframe_region`).
- No `other` or `cap` rows remain. The ledger is `/root/scratch_dnakit/ledger/bioframe.json`.

## Example input files

- **Example file：** `a.bed`

```text
chr1	0	100	a	5	+
chr1	50	150	b	3	+
chr1	300	400	c	1	-
chr2	10	60	d	2	+
```

BED6, four intervals on two chromosomes

- **Example file：** `view.bed`

```text
chr1	0	200	v1
chr2	0	100	v2
```

BED4 view with names v1 and v2

- **Example file：** `labelled.bed`

```text
chr1	0	10	v1
chr2	0	5	v2
```

BED4 whose names match the view

- **Example file：** `inside.bed`

```text
chr1	10	20
chr2	0	50
```

BED3 lying inside the view

- **Example file：** `tiles.bed`

```text
chr1	0	100
chr1	100	200
```

BED3 with two abutting intervals

- **Example file：** `chr1_view.bed`

```text
chr1	0	200	v
```

BED4 single view region chr1:0-200

- **Example file：** `s.bed`

```text
chr10	5	9
chr2	1	3
chr1	7	8
chr2	0	1
```

BED3, unsorted

- **Example file：** `bad.bed`

```text
chr1	20	10	n
chr1	0	5	m
```

BED4, first row has start greater than end

- **Example file：** `b.bed`

```text
chr1	80	120	x
chr2	0	60	y
```

BED4 with two intervals

- **Example file：** `runs.bed`

```text
chr1	0	100	x
chr1	100	200	x
chr1	200	300	y
```

BED4; neighbours with the same name form one run

- **Example file：** `g.fa`

```text
>chr1
ACGTNNGGCCGATTACAGAATTCGC
>chr2
GGGGCCCCAAAATTTTGCGCGCAT
```

FASTA with two sequences; the .fai index is created on first use

- **Example file：** `small.bed`

```text
chr1	0	10	r1
chr1	10	20	r2
chr2	0	12	r3
```

BED4 with three short intervals

- **Example file：** `x.sam`

```text
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:1000
r1	0	chr1	11	60	8M	*	0	0	ACGTACGT	IIIIIIII	NM:i:0
r2	16	chr1	101	30	4M2D4M	*	0	0	ACGTACGT	IIIIIIII
```

SAM with two alignments; POS is 1-based in SAM

- **Example file：** `mrna.tsv`

```text
matches	misMatches	repMatches	nCount	qNumInsert	qBaseInsert	tNumInsert	tBaseInsert	strand	qName	qSize	qStart	qEnd	tName	tSize	tStart	tEnd	blockCount	blockSizes	qStarts	tStarts
100	0	0	0	0	0	0	0	+	NM1	500	0	100	chr1	1000	20	120	1	100,	0,	20,
```

Local mRNA table with a header; chr1 20-120

- **Example file：** `tab.bed`

```text
chr1	10	20	n
chr1	120	130	m
chr2	90	99	k
```

BED4; the example compresses it and builds a tabix index

- **Example file：** `sizes.txt`

```text
chr1	1000
chr2	500
chrUn_x	10
```

chrom.sizes file

## 1) `validate_bioframe_bed`

- **Purpose:** Check a BED3 to BED6 file against the bioframe bedframe rules.
- **Method:** Reads the file with `read_table` using the schema that matches its column count, then calls `is_bedframe(raise_errors=False)`. An invalid file gives `passed` false; nothing is raised and the file is not changed.
- **API:** `validate_bioframe_bed(source[required], max_rows[optional], verbose[optional])`
- **Input:** A local BED3 to BED6 file with at most `max_rows` rows.
- **Output:** `passed` (boolean), `record_count`, `provenance`.

**Example code:**

```python
from dnakit.annotation import validate_bioframe_bed

result = validate_bioframe_bed("a.bed", verbose=False)
print(result["passed"], result["record_count"])
```

**Example output:**

```text
True 4
```

## 2) `check_bioframe_cataloged`

- **Purpose:** Check that every interval name of a BED4+ file is a region name of a view.
- **Method:** Calls `is_cataloged(df_view_col="name")`. Only names are compared, not coordinates.
- **API:** `check_bioframe_cataloged(source[required], view[required], max_rows[optional], verbose[optional])`
- **Input:** A source BED4+ file and a view file (BED4+) whose names are the region names.
- **Output:** `passed`, `check` (always cataloged), `record_count`, `provenance`.

**Example code:**

```python
from dnakit.annotation import check_bioframe_cataloged

labelled = check_bioframe_cataloged("labelled.bed", "view.bed", verbose=False)
named = check_bioframe_cataloged("a.bed", "view.bed", verbose=False)
print(labelled["passed"], named["passed"])
```

**Example output:**

```text
True False
```

## 3) `check_bioframe_contained`

- **Purpose:** Check that every interval lies entirely inside some view region.
- **Method:** Calls `is_contained`. Intervals and views are compared as half-open ranges; an interval whose ends match the view ends is still contained.
- **API:** `check_bioframe_contained(source[required], view[required], max_rows[optional], verbose[optional])`
- **Input:** A source BED3 to BED6 file and a view file (BED4+).
- **Output:** `passed`, `check` (always contained), `record_count`, `provenance`.

**Example code:**

```python
from dnakit.annotation import check_bioframe_contained

inside = check_bioframe_contained("inside.bed", "view.bed", verbose=False)
outside = check_bioframe_contained("a.bed", "view.bed", verbose=False)
print(inside["passed"], outside["passed"])
```

**Example output:**

```text
True False
```

## 4) `check_bioframe_covering`

- **Purpose:** Check that the intervals of a source file cover every view region without gaps.
- **Method:** Calls `is_covering`. Passes only when the union of the intervals spans each view region with no gap.
- **API:** `check_bioframe_covering(source[required], view[required], max_rows[optional], verbose[optional])`
- **Input:** A source BED3 to BED6 file and a view file (BED4+).
- **Output:** `passed`, `check` (always covering), `record_count`, `provenance`.

**Example code:**

```python
from dnakit.annotation import check_bioframe_covering

result = check_bioframe_covering("tiles.bed", "chr1_view.bed", verbose=False)
print(result["passed"], result["record_count"])
```

**Example output:**

```text
True 2
```

## 5) `check_bioframe_overlapping`

- **Purpose:** Report whether any two intervals of a file share a base.
- **Method:** Calls `is_overlapping`. Intervals that only touch at one end (half-open) do not count as overlapping.
- **API:** `check_bioframe_overlapping(source[required], max_rows[optional], verbose[optional])`
- **Input:** A local BED3 to BED6 file.
- **Output:** `passed` (true means at least one overlap), `check`, `record_count`, `provenance`.

**Example code:**

```python
from dnakit.annotation import check_bioframe_overlapping

print(check_bioframe_overlapping("a.bed", verbose=False)["passed"])
```

**Example output:**

```text
True
```

## 6) `check_bioframe_sorted`

- **Purpose:** Report whether a BED file is already in bioframe sort order.
- **Method:** Calls `is_sorted`. `passed` is true when sorting would leave the file unchanged.
- **API:** `check_bioframe_sorted(source[required], max_rows[optional], verbose[optional])`
- **Input:** A local BED3 to BED6 file.
- **Output:** `passed`, `check` (always sorted), `record_count`, `provenance`.

**Example code:**

```python
from dnakit.annotation import check_bioframe_sorted

print(check_bioframe_sorted("s.bed", verbose=False)["passed"])
```

**Example output:**

```text
False
```

## 7) `check_bioframe_tiling`

- **Purpose:** Check that the intervals tile every view region end to end.
- **Method:** Calls `is_tiling(df_view_col=None)`. Intervals may neither overlap nor leave gaps.
- **API:** `check_bioframe_tiling(source[required], view[required], max_rows[optional], verbose[optional])`
- **Input:** A source BED3 to BED6 file and a view file (BED4+).
- **Output:** `passed`, `check` (always tiling), `record_count`, `provenance`.

**Example code:**

```python
from dnakit.annotation import check_bioframe_tiling

print(check_bioframe_tiling("tiles.bed", "chr1_view.bed", verbose=False)["passed"])
```

**Example output:**

```text
True
```

## 8) `check_bioframe_viewframe`

- **Purpose:** Check that a file is a valid view frame: unique names and valid intervals.
- **Method:** Calls `is_viewframe(raise_errors=False)`. A file without a name column (BED3) does not pass.
- **API:** `check_bioframe_viewframe(source[required], max_rows[optional], verbose[optional])`
- **Input:** A local file of BED3 or wider; a view should be BED4 to BED6.
- **Output:** `passed`, `check` (always viewframe), `record_count`, `provenance`.

**Example code:**

```python
from dnakit.annotation import check_bioframe_viewframe

print(check_bioframe_viewframe("view.bed", verbose=False)["passed"])
```

**Example output:**

```text
True
```

## 9) `resolve_bioframe_region`

- **Purpose:** Parse a region string into chrom, start and end; a bare chromosome name takes its end from chromsizes.
- **Method:** Calls `parse_region`. Numbers are read as written with commas removed, with no coordinate shift. A region without an end and without chromsizes is an error.
- **API:** `resolve_bioframe_region(region[required], chromsizes[optional], check_bounds[optional], verbose[optional])`
- **Input:** A region string such as `chr1:100-200` or `chr1`, and an optional chromsizes mapping.
- **Output:** `region`, `chrom`, `start`, `end`, `length`, `provenance`.

**Example code:**

```python
from dnakit.annotation import resolve_bioframe_region

span = resolve_bioframe_region("chr1:100-200", verbose=False)
whole = resolve_bioframe_region("chr1", chromsizes={"chr1": 1000}, verbose=False)
print(span["chrom"], span["start"], span["end"], whole["start"], whole["end"])
```

**Example output:**

```text
chr1 100 200 0 1000
```

## 10) `format_bioframe_region`

- **Purpose:** Write a chrom, start and end triple as a UCSC region string; the inverse of `resolve_bioframe_region`.
- **Method:** Calls `to_ucsc_string`, writing coordinates as given with no shift. Requires start less than end and a non-negative start.
- **API:** `format_bioframe_region(chrom[required], start[required], end[required], verbose[optional])`
- **Input:** A chromosome name and integer start and end.
- **Output:** `region`, `chrom`, `start`, `end`, `provenance`.

**Example code:**

```python
from dnakit.annotation import format_bioframe_region

print(format_bioframe_region("chr1", 0, 100, verbose=False)["region"])
```

**Example output:**

```text
chr1:0-100
```

## 11) `sanitize_bioframe_intervals`

- **Purpose:** Clean BED intervals: handle rows whose start is greater than end, and optionally drop null rows.
- **Method:** Calls `sanitize_bedframe`. `start_exceed_end_action` is None (raise on such rows), `flip` (swap the ends) or `drop` (remove the row).
- **API:** `sanitize_bioframe_intervals(source[required], start_exceed_end_action[optional], drop_null[optional], max_rows[optional], verbose[optional])`
- **Input:** A local BED3 to BED6 file.
- **Output:** `intervals`, `interval_count`, `start_exceed_end_action`, `drop_null`, `provenance`.

**Example code:**

```python
from dnakit.annotation import sanitize_bioframe_intervals

result = sanitize_bioframe_intervals("bad.bed", start_exceed_end_action="flip", verbose=False)
print([(row["start"], row["end"]) for row in result["intervals"]])
```

**Example output:**

```text
[(10, 20), (0, 5)]
```

## 12) `sort_bioframe_intervals`

- **Purpose:** Sort BED intervals by chromosome name as text, then by start.
- **Method:** Calls `sort_bedframe`. Chromosome names compare as text, so chr10 sorts before chr2.
- **API:** `sort_bioframe_intervals(source[required], max_rows[optional], verbose[optional])`
- **Input:** A local BED3 to BED6 file.
- **Output:** `intervals`, `interval_count`, `provenance`.

**Example code:**

```python
from dnakit.annotation import sort_bioframe_intervals

result = sort_bioframe_intervals("s.bed", verbose=False)
print([(row["chrom"], row["start"], row["end"]) for row in result["intervals"]])
```

**Example output:**

```text
[('chr1', 7, 8), ('chr10', 5, 9), ('chr2', 0, 1), ('chr2', 1, 3)]
```

## 13) `setdiff_bioframe_intervals`

- **Purpose:** Drop each interval that overlaps any interval of another file.
- **Method:** Calls `setdiff`. Whole intervals are removed; partly overlapping intervals are not clipped (use the existing `subtract_bioframe_intervals` to clip).
- **API:** `setdiff_bioframe_intervals(source[required], other[required], max_rows[optional], verbose[optional])`
- **Input:** Two local BED3 to BED6 files.
- **Output:** `intervals`, `interval_count`, `provenance`.

**Example code:**

```python
from dnakit.annotation import setdiff_bioframe_intervals

result = setdiff_bioframe_intervals("a.bed", "b.bed", verbose=False)
print([(row["chrom"], row["start"], row["end"]) for row in result["intervals"]])
```

**Example output:**

```text
[('chr1', 300, 400)]
```

## 14) `trim_bioframe_intervals`

- **Purpose:** Clip each interval to the view regions; intervals entirely outside the view are removed.
- **Method:** Calls `trim`. Rows that fall outside the view come back null from bioframe; the wrapper removes them and reports the count in `dropped_count`.
- **API:** `trim_bioframe_intervals(source[required], view[required], max_rows[optional], verbose[optional])`
- **Input:** A source BED3 to BED6 file and a view file (BED4+).
- **Output:** `intervals`, `interval_count`, `dropped_count`, `provenance`.

**Example code:**

```python
from dnakit.annotation import trim_bioframe_intervals

result = trim_bioframe_intervals("a.bed", "view.bed", verbose=False)
print(result["interval_count"], result["dropped_count"], [(row["start"], row["end"]) for row in result["intervals"]])
```

**Example output:**

```text
3 1 [(0, 100), (50, 150), (10, 60)]
```

## 15) `complement_bioframe_intervals`

- **Purpose:** Return the parts of each view region that no interval covers.
- **Method:** Calls `complement`. A view is required; without one bioframe uses a very large placeholder end, so the wrapper does not accept that input.
- **API:** `complement_bioframe_intervals(source[required], view[required], max_rows[optional], verbose[optional])`
- **Input:** A source BED3 to BED6 file and a view file (BED4+).
- **Output:** `gaps`, `gap_count`, `provenance`; each gap carries `view_region`, the view it belongs to.

**Example code:**

```python
from dnakit.annotation import complement_bioframe_intervals

result = complement_bioframe_intervals("b.bed", "view.bed", verbose=False)
print([(row["chrom"], row["start"], row["end"]) for row in result["gaps"]])
```

**Example output:**

```text
[('chr1', 0, 80), ('chr1', 120, 200), ('chr2', 60, 100)]
```

## 16) `count_bioframe_overlaps`

- **Purpose:** Count the intervals of another file that overlap each interval of the source.
- **Method:** Calls `count_overlaps`. Source row order is kept and a `count` column is added.
- **API:** `count_bioframe_overlaps(source[required], other[required], max_rows[optional], verbose[optional])`
- **Input:** Two local BED3 to BED6 files.
- **Output:** `intervals` (with `count`), `interval_count`, `provenance`.

**Example code:**

```python
from dnakit.annotation import count_bioframe_overlaps

result = count_bioframe_overlaps("a.bed", "b.bed", verbose=False)
print([(row["chrom"], row["start"], row["count"]) for row in result["intervals"]])
```

**Example output:**

```text
[('chr1', 0, 1), ('chr1', 50, 1), ('chr1', 300, 0), ('chr2', 10, 1)]
```

## 17) `select_bioframe_region`

- **Purpose:** Return the intervals that overlap a region string.
- **Method:** Calls `select`. A bare chromosome name selects every interval on that chromosome.
- **API:** `select_bioframe_region(source[required], region[required], max_rows[optional], verbose[optional])`
- **Input:** A local BED3 to BED6 file and a region string.
- **Output:** `region`, `intervals`, `interval_count`, `provenance`.

**Example code:**

```python
from dnakit.annotation import select_bioframe_region

result = select_bioframe_region("a.bed", "chr1:90-120", verbose=False)
print([row["name"] for row in result["intervals"]])
```

**Example output:**

```text
['a', 'b']
```

## 18) `assign_bioframe_views`

- **Purpose:** Label each interval with the name of the view region that contains it.
- **Method:** Calls `assign_view`. Intervals that overlap no view get a null `view_region`; with `drop_unassigned` true they are removed.
- **API:** `assign_bioframe_views(source[required], view[required], drop_unassigned[optional], max_rows[optional], verbose[optional])`
- **Input:** A source BED3 to BED6 file and a view file (BED4+).
- **Output:** `intervals` (with `view_region`), `interval_count`, `drop_unassigned`, `provenance`.

**Example code:**

```python
from dnakit.annotation import assign_bioframe_views

result = assign_bioframe_views("a.bed", "view.bed", verbose=False)
print([(row["name"], row["view_region"]) for row in result["intervals"]])
```

**Example output:**

```text
[('a', 'v1'), ('b', 'v1'), ('c', None), ('d', 'v2')]
```

## 19) `mark_bioframe_runs`

- **Purpose:** Number runs of neighbouring intervals that share one column value.
- **Method:** Calls `mark_runs`. Intervals must not overlap unless `allow_overlaps` is true. The column is name, score or strand and must exist in the file.
- **API:** `mark_bioframe_runs(source[required], column[required], allow_overlaps[optional], max_rows[optional], verbose[optional])`
- **Input:** A local BED file and a column name (name, score or strand).
- **Output:** `column`, `intervals` (with `run`), `interval_count`, `provenance`.

**Example code:**

```python
from dnakit.annotation import mark_bioframe_runs

result = mark_bioframe_runs("runs.bed", "name", verbose=False)
print([(row["name"], row["run"]) for row in result["intervals"]])
```

**Example output:**

```text
[('x', 0), ('x', 0), ('y', 1)]
```

## 20) `merge_bioframe_runs`

- **Purpose:** Merge neighbouring intervals that share one column value into single runs.
- **Method:** Calls `merge_runs`. The rules match `mark_bioframe_runs`.
- **API:** `merge_bioframe_runs(source[required], column[required], allow_overlaps[optional], max_rows[optional], verbose[optional])`
- **Input:** A local BED file and a column name (name, score or strand).
- **Output:** `column`, `intervals`, `interval_count`, `provenance`.

**Example code:**

```python
from dnakit.annotation import merge_bioframe_runs

result = merge_bioframe_runs("runs.bed", "name", verbose=False)
print([(row["chrom"], row["start"], row["end"]) for row in result["intervals"]])
```

**Example output:**

```text
[('chr1', 0, 200), ('chr1', 200, 300)]
```

## 21) `pair_bioframe_intervals_by_distance`

- **Purpose:** List unique interval pairs whose separation falls in a window.
- **Method:** Calls `pair_by_distance`. Separation is measured between midpoints by default and satisfies `min_sep <= distance < max_sep`. min_sep must be less than max_sep.
- **API:** `pair_bioframe_intervals_by_distance(source[required], min_sep[optional], max_sep[optional], min_intervening[optional], max_intervening[optional], relative_to[optional], max_rows[optional], verbose[optional])`
- **Input:** A local BED3 to BED6 file, with optional limits on the intervening distance.
- **Output:** `relative_to`, `pairs` (columns suffixed `_1` and `_2`), `pair_count`, `provenance`.

**Example code:**

```python
from dnakit.annotation import pair_bioframe_intervals_by_distance

result = pair_bioframe_intervals_by_distance("a.bed", min_sep=0, max_sep=280, verbose=False)
print(result["pair_count"], [(row["name_1"], row["name_2"]) for row in result["pairs"]])
```

**Example output:**

```text
2 [('a', 'b'), ('b', 'c')]
```

## 22) `sequence_gc_fraction`

- **Purpose:** Return the fraction of G and C in one DNA string.
- **Method:** Checks the characters with `sequence_text` (IUPAC codes allowed), then calls `seq_gc`. With `mapped_only`, N is excluded from the denominator; null when no base is counted.
- **API:** `sequence_gc_fraction(sequence[required], mapped_only[optional], verbose[optional])`
- **Input:** A DNA string, case-insensitive, at most ten million bases.
- **Output:** `length`, `gc_fraction`, `mapped_only`, `provenance`.

**Example code:**

```python
from dnakit.annotation import sequence_gc_fraction

result = sequence_gc_fraction("GGCCAATT", verbose=False)
print(result["gc_fraction"], result["length"])
```

**Example output:**

```text
0.5 8
```

## 23) `digest_bioframe_fasta`

- **Purpose:** Cut each FASTA sequence into restriction fragments.
- **Method:** Calls `bioframe.digest` (needs Biopython). bioframe adds 1 to Biopython's cut position, which moves each interior cut two bases to the right. The wrapper subtracts those two bases so EcoRI cuts after the G of G^AATTC. The test checks this against Biopython's `EcoRI.search` directly.
- **API:** `digest_bioframe_fasta(source[required], enzyme[required], max_rows[optional], verbose[optional])`
- **Input:** A local FASTA and a Biopython enzyme name such as EcoRI.
- **Output:** `enzyme`, `fragments`, `fragment_count`, `provenance`.

**Example code:**

```python
from dnakit.annotation import digest_bioframe_fasta

result = digest_bioframe_fasta("g.fa", "EcoRI", verbose=False)
print([(row["chrom"], row["start"], row["end"]) for row in result["fragments"]])
```

**Example output:**

```text
[('chr1', 0, 18), ('chr1', 18, 25), ('chr2', 0, 24)]
```

## 24) `frac_gc_bioframe_intervals`

- **Purpose:** Add the GC fraction of each BED interval of a FASTA sequence.
- **Method:** Calls `frac_gc` and renames its `GC` column to `frac_gc`. With `mapped_only`, N is excluded from the denominator. An interval whose chromosome is missing from the FASTA is an error.
- **API:** `frac_gc_bioframe_intervals(source[required], fasta[required], mapped_only[optional], max_rows[optional], verbose[optional])`
- **Input:** A BED3 to BED6 file and a FASTA file.
- **Output:** `mapped_only`, `intervals` (with `frac_gc`), `interval_count`, `provenance`.

**Example code:**

```python
from dnakit.annotation import frac_gc_bioframe_intervals

result = frac_gc_bioframe_intervals("small.bed", "g.fa", verbose=False)
print([round(row["frac_gc"], 3) for row in result["intervals"]])
```

**Example output:**

```text
[0.75, 0.3, 0.667]
```

## 25) `frac_mapped_bioframe_intervals`

- **Purpose:** Add the fraction of non-N bases of each BED interval of a FASTA sequence.
- **Method:** Calls `frac_mapped`. The denominator is the number of sequence bases inside the interval; an interval outside the sequence gets 0.
- **API:** `frac_mapped_bioframe_intervals(source[required], fasta[required], max_rows[optional], verbose[optional])`
- **Input:** A BED3 to BED6 file and a FASTA file.
- **Output:** `intervals` (with `frac_mapped`), `interval_count`, `provenance`.

**Example code:**

```python
from dnakit.annotation import frac_mapped_bioframe_intervals

result = frac_mapped_bioframe_intervals("small.bed", "g.fa", verbose=False)
print([round(row["frac_mapped"], 3) for row in result["intervals"]])
```

**Example output:**

```text
[0.8, 1.0, 1.0]
```

## 26) `frac_gene_coverage_bioframe`

- **Purpose:** Measure how many bases of each interval mRNA alignments cover, and how many alignments overlap it.
- **Method:** Calls `frac_gene_coverage`. The mRNA table must be a local tab-separated file with a header and the columns tName, tStart and tEnd. The upstream gives base counts (`coverage`) and overlap counts (`count`), not fractions.
- **API:** `frac_gene_coverage_bioframe(source[required], mrna[required], max_rows[optional], verbose[optional])`
- **Input:** A source BED3 to BED6 file and a local mRNA table (UCSC all_mrna layout with a header).
- **Output:** `intervals` (with `coverage` and `count`), `interval_count`, `provenance`.

**Example code:**

```python
from dnakit.annotation import frac_gene_coverage_bioframe

result = frac_gene_coverage_bioframe("a.bed", "mrna.tsv", verbose=False)
print([(row["name"], row["coverage"], row["count"]) for row in result["intervals"]])
```

**Example output:**

```text
[('a', 80, 1), ('b', 70, 1), ('c', 0, 0), ('d', 0, 0)]
```

## 27) `make_bioframe_chromarms`

- **Purpose:** Split each chromosome at its centromere midpoint into p and q arms.
- **Method:** Calls `make_chromarms`. Every chromosome needs a midpoint inside its length. Without one, bioframe keeps the whole chromosome as a single arm, so the wrapper rejects that input.
- **API:** `make_bioframe_chromarms(chromsizes[required], midpoints[required], max_rows[optional], verbose[optional])`
- **Input:** A chromosome-length mapping and a centromere midpoint mapping with the same chromosome set.
- **Output:** `arms` (each with `name`, such as chr1_p), `arm_count`, `provenance`.

**Example code:**

```python
from dnakit.annotation import make_bioframe_chromarms

result = make_bioframe_chromarms({"chr1": 1000, "chr2": 500}, {"chr1": 400, "chr2": 200}, verbose=False)
print([(row["name"], row["start"], row["end"]) for row in result["arms"]])
```

**Example output:**

```text
[('chr1_p', 0, 400), ('chr1_q', 400, 1000), ('chr2_p', 0, 200), ('chr2_q', 200, 500)]
```

## 28) `list_bioframe_assemblies`

- **Purpose:** List the genome assembly metadata that ships with bioframe.
- **Method:** Calls `assemblies_available`, which reads the YAML file bundled with the package and makes no network call. Results are sorted by provider and build.
- **API:** `list_bioframe_assemblies(verbose[optional])`
- **Input:** No input arguments.
- **Output:** `assemblies` (provider, provider_build, organism, default roles and units, url), `assembly_count`, `provenance`.

**Example code:**

```python
from dnakit.annotation import list_bioframe_assemblies

result = list_bioframe_assemblies(verbose=False)
hg38 = [row for row in result["assemblies"] if row["provider_build"] == "hg38"]
print(result["assembly_count"] > 10, hg38[0]["organism"], hg38[0]["provider"])
```

**Example output:**

```text
True homo sapiens ucsc
```

## 29) `describe_bioframe_assembly`

- **Purpose:** Return the sequence names and lengths of a bundled assembly, such as ucsc.hg38.
- **Method:** Calls `assembly_info`. Without roles and units the assembly defaults apply. An unknown or ambiguous name is an error.
- **API:** `describe_bioframe_assembly(name[required], roles[optional], units[optional], max_rows[optional], verbose[optional])`
- **Input:** An assembly name (provider.build works), and optional roles and units (None, `all`, or a list of strings).
- **Output:** `name`, `organism`, `provider`, `provider_build`, `release_year`, `url`, `sequences`, `sequence_count`, `provenance`.

**Example code:**

```python
from dnakit.annotation import describe_bioframe_assembly

result = describe_bioframe_assembly("ucsc.hg38", verbose=False)
first = result["sequences"][0]
print(result["sequence_count"], first["name"], first["length"])
```

**Example output:**

```text
25 chr1 248956422
```

## 30) `summarize_bioframe_fasta`

- **Purpose:** List the name and length of each FASTA sequence without reading the sequences into memory.
- **Method:** Calls `load_fasta`; lengths come from the pysam index. When the .fai index is missing, pysam writes it beside the file.
- **API:** `summarize_bioframe_fasta(source[required], max_rows[optional], verbose[optional])`
- **Input:** A local FASTA file.
- **Output:** `sequences` (name and length), `sequence_count`, `total_length`, `provenance`.

**Example code:**

```python
from dnakit.annotation import summarize_bioframe_fasta

result = summarize_bioframe_fasta("g.fa", verbose=False)
print(result["sequences"], result["total_length"])
```

**Example output:**

```text
[{'name': 'chr1', 'length': 25}, {'name': 'chr2', 'length': 24}] 49
```

## 31) `read_bioframe_chromsizes`

- **Purpose:** Read a local chrom.sizes or chromInfo file.
- **Method:** Calls `read_chromsizes` and keeps only the name and length. With `filter_chroms`, only chr1 to chr22, chrX, chrY and chrM are kept; `natsort` orders names naturally.
- **API:** `read_bioframe_chromsizes(source[required], filter_chroms[optional], natsort[optional], max_rows[optional], verbose[optional])`
- **Input:** A local chrom.sizes file. UCSC database names are not accepted because they need the network.
- **Output:** `chromsizes` (each with chrom and length), `chromosome_count`, `provenance`.

**Example code:**

```python
from dnakit.annotation import read_bioframe_chromsizes

result = read_bioframe_chromsizes("sizes.txt", verbose=False)
print(result["chromsizes"])
```

**Example output:**

```text
[{'chrom': 'chr1', 'length': 1000}, {'chrom': 'chr2', 'length': 500}]
```

## 32) `read_bioframe_tabix`

- **Purpose:** Read a bgzip-compressed BED file with a tabix index, whole or by region.
- **Method:** Calls `read_tabix` (pysam). A file without a .tbi or .csi index is an error; the message points to `pysam.tabix_index`. A region read needs chrom.
- **API:** `read_bioframe_tabix(source[required], chrom[optional], start[optional], end[optional], max_rows[optional], verbose[optional])`
- **Input:** A bgzip-compressed BED file and its tabix index. The example builds tab.bed.gz and the .tbi from tab.bed.
- **Output:** `intervals`, `interval_count`, `provenance`.

**Example code:**

```python
import pysam
from dnakit.annotation import read_bioframe_tabix

pysam.tabix_compress("tab.bed", "tab.bed.gz", force=True)
pysam.tabix_index("tab.bed.gz", preset="bed", force=True)
whole = read_bioframe_tabix("tab.bed.gz", verbose=False)
region = read_bioframe_tabix("tab.bed.gz", chrom="chr1", start=0, end=100, verbose=False)
print(whole["interval_count"], [row["name"] for row in region["intervals"]])
```

**Example output:**

```text
3 ['n']
```

## 33) `read_bioframe_alignments`

- **Purpose:** Read SAM or BAM alignment records.
- **Method:** Calls `read_alignments` (pysam). The returned POS is 0-based, one less than the SAM column. Region reads need an indexed BAM; a SAM file gives an error.
- **API:** `read_bioframe_alignments(source[required], chrom[optional], start[optional], end[optional], max_rows[optional], verbose[optional])`
- **Input:** A local SAM or BAM file, with optional chrom, start and end.
- **Output:** `alignments` (with QNAME, FLAG, RNAME, POS, MAPQ, CIGAR and more), `alignment_count`, `provenance`.

**Example code:**

```python
from dnakit.annotation import read_bioframe_alignments

result = read_bioframe_alignments("x.sam", verbose=False)
print(result["alignment_count"], [(row["QNAME"], row["POS"]) for row in result["alignments"]])
```

**Example output:**

```text
2 [('r1', 10), ('r2', 100)]
```

## 34) `read_bioframe_bigbed`

- **Purpose:** Read the intervals of one chromosome from a local bigBed file.
- **Method:** Calls `read_bigbed` with the pyBigWig engine. An unknown chromosome is an error. The example uses the `test.bigBed` file shipped with the pyBigWig test package, since a binary file cannot be embedded in the page.
- **API:** `read_bioframe_bigbed(source[required], chrom[required], start[optional], end[optional], max_rows[optional], verbose[optional])`
- **Input:** A local bigBed file, a chromosome name, and optional start and end.
- **Output:** `chrom`, `intervals` (with chrom, start, end and rest), `interval_count`, `provenance`.

**Example code:**

```python
import os
import pyBigWigTest
from dnakit.annotation import read_bioframe_bigbed

path = os.path.join(os.path.dirname(pyBigWigTest.__file__), "test.bigBed")
result = read_bioframe_bigbed(path, "chr1", verbose=False)
print(result["chrom"], result["interval_count"] > 0, result["intervals"][0]["start"])
```

**Example output:**

```text
chr1 True 10009333
```

## 35) `read_bioframe_table`

- **Purpose:** Read a tab-separated BED3 to BED12 table and check that its schema matches the column count.
- **Method:** Checks the column count of every row first (all rows must match), then calls `read_table`. Lines starting with # are skipped; track and browser lines are not.
- **API:** `read_bioframe_table(source[required], schema[optional], max_rows[optional], verbose[optional])`
- **Input:** A local table file. The schema may be omitted (inferred from the column count for bed3 to bed6, bed9 and bed12) or given explicitly.
- **Output:** `schema`, `rows`, `row_count`, `provenance`.

**Example code:**

```python
from dnakit.annotation import read_bioframe_table

result = read_bioframe_table("a.bed", verbose=False)
print(result["schema"], result["row_count"], result["rows"][0]["name"])
```

**Example output:**

```text
bed6 4 a
```

## 36) `write_bioframe_bed`

- **Purpose:** Re-write a BED file through `to_bed` to a new path after checking its fields.
- **Method:** Calls `to_bed`. The output path must not exist and its directory must exist. With `require_sorted` true, unsorted input is rejected. A failed write removes the partial output.
- **API:** `write_bioframe_bed(source[required], output[required], schema[optional], require_sorted[optional], max_rows[optional], verbose[optional])`
- **Input:** A source BED3 to BED6 file, a new output path, a schema (infer, or bed3 to bed12), and require_sorted.
- **Output:** `output`, `schema`, `record_count`, `provenance`. The output has no header and is tab-separated.

**Example code:**

```python
import pathlib
from dnakit.annotation import write_bioframe_bed

out = pathlib.Path("sorted_out.bed")
out.unlink(missing_ok=True)
result = write_bioframe_bed("a.bed", out, verbose=False)
print(result["record_count"], out.read_text(encoding="utf-8").splitlines()[0].split("\t")[:3])
```

**Example output:**

```text
4 ['chr1', '0', '100']
```
