# pybedtools: conversions, format checks and interval transforms

This group in `dnakit.annotation` wraps 15 operations of [pybedtools](https://daler.github.io/pybedtools/) that DNAKit did not have: BED12 to BED6, BAM to BED, BED or BEDPE to BAM, an interval table as records, tabix contig names, BAM/CRAM/BGZF detection, an overlap test of one interval against a file, slop, shift, fixed windows, concatenate-and-merge, truncation to chromosome bounds, removal of invalid records, introns from gene and exon features, and BAM tagging with interval files.

pybedtools drives the [BEDTools](https://bedtools.readthedocs.io/) programs, so BEDTools must be installed (tested with BEDTools 2.27.1 and pybedtools 0.12.1; install with `pip install "dnakit[bio-pybedtools]"` and `conda install -c bioconda bedtools`). Coordinates are 0-based half-open as in BED. Every function takes file paths and returns a JSON-friendly dict with the first `max_rows` records; `output_path` (a path that does not exist yet) receives the full result. The examples import from the module `dnakit.annotation.pybedtools_p33_functions`.

Other rows of the pybedtools work list are reached through functions that DNAKit already has: window search `bedtools_window`, shuffling `bedtools_shuffle`, multi-file intersection `bedtools_multiinter`, BEDPE pairs against BEDPE pairs `bedtools_pairtopair`, bedGraph union `bedtools_unionbedg`, multi-BAM coverage `bedtools_multicov`, and score normalisation `normalize_pybedtools_scores_to_length`.

## Example data

All examples use the files below, written into the working folder before the code runs.

**Example file:** `a.bed`

```text
chr1	100	200	a1	5	+
chr1	500	600	a2	6	-
chr2	10	90	a3	7	+
```

**Example file:** `b.bed`

```text
chr1	150	250	b1	1	+
chr1	900	950	b2	2	-
chr2	100	150	b3	3	+
```

**Example file:** `genome.sizes`

```text
chr1	2000
chr2	1000
```

**Example file:** `bed12.bed`

```text
chr1	100	400	g1	0	+	100	400	0	2	50,60	0,240
```

**Example file:** `pairs.bedpe`

```text
chr1	100	150	chr1	300	350	p1	5	+	-
chr2	10	60	chr2	200	250	p2	6	+	+
```

**Example file:** `genes.gff`

```text
chr1	t	gene	101	600	.	+	.	ID=g1
chr1	t	exon	101	200	.	+	.	Parent=g1
chr1	t	exon	301	400	.	+	.	Parent=g1
chr1	t	exon	501	600	.	+	.	Parent=g1
```

**Example file:** `wide.bed`

```text
chr1	1900	2500	w1
chr3	1	5	w2
chr2	5	50	w3
```

**Example file:** `bad.bed`

```text
chr1	10	100
chr1	100	10
chr1	5	6
```

**Example file:** `make_bam.py`

```python
import pysam

header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 2000}, {"SN": "chr2", "LN": 1000}]}
with pysam.AlignmentFile("reads.bam", "wb", header=header) as out:
    for i, (ref, pos) in enumerate([(0, 120), (0, 520), (0, 1500), (1, 20)]):
        read = pysam.AlignedSegment()
        read.query_name = f"read{i}"
        read.query_sequence = "A" * 40
        read.flag = 0
        read.reference_id = ref
        read.reference_start = pos
        read.mapping_quality = 60
        read.cigar = [(0, 40)]
        read.query_qualities = pysam.qualitystring_to_array("I" * 40)
        out.write(read)
```

## 1) `dnakit.annotation.convert_bed12_to_bed6_pybedtools()`

- **Purpose:** Split every BED12 record into one BED6 record per block.
- **Method:** Call `BedTool.bed6()`, which runs `bedtools bed12tobed6`. With `number_names=True` the 1-based block number is written into the score column (`-n`).
- **API:** `dnakit.annotation.convert_bed12_to_bed6_pybedtools(bed12_path: 'str | Path', *, number_names: 'bool' = False, output_path: 'str | Path | None' = None, max_rows: 'int' = 1000)`.
- **Input:** A BED12 file path, optional `number_names`, `output_path` and `max_rows`.
- **Output:** `n_records` and `records`, each a list of the six BED fields as strings.
- **Example code:**

```python
from dnakit.annotation.pybedtools_p33_functions import convert_bed12_to_bed6_pybedtools

result = convert_bed12_to_bed6_pybedtools("bed12.bed", verbose=False)
print(result["n_records"])
for record in result["records"]:
    print(record)
numbered = convert_bed12_to_bed6_pybedtools("bed12.bed", number_names=True, verbose=False)
print([r[4] for r in numbered["records"]])
```

- **Example output:**

```text
2
['chr1', '100', '150', 'g1', '0', '+']
['chr1', '340', '400', 'g1', '0', '+']
['1', '2']
```

## 2) `dnakit.annotation.bam_to_bed_pybedtools()`

- **Purpose:** Convert BAM alignments to BED, BED12 or split BED.
- **Method:** Call `BedTool.bam_to_bed()` (`bedtools bamtobed`). `bed12`, `split`, `split_deletions`, `edit_distance`, `tag`, `color` and `cigar` map to `-bed12`, `-split`, `-splitD`, `-ed`, `-tag`, `-color` and `-cigar`.
- **API:** `dnakit.annotation.bam_to_bed_pybedtools(bam_path: 'str | Path', *, bed12: 'bool' = False, split: 'bool' = False, split_deletions: 'bool' = False, edit_distance: 'bool' = False, tag: 'str | None' = None, color: 'str | None' = None, cigar: 'bool' = False, output_path: 'str | Path | None' = None, max_rows: 'int' = 1000)`.
- **Input:** A BAM file path and the options above.
- **Output:** `n_records` and `records` (chrom, start, end, read name, mapping quality, strand).
- **Example code:**

```python
from dnakit.annotation.pybedtools_p33_functions import bam_to_bed_pybedtools

exec(open("make_bam.py").read())
result = bam_to_bed_pybedtools("reads.bam", verbose=False)
print(result["n_records"])
for record in result["records"]:
    print(record)
```

- **Example output:**

```text
4
['chr1', '120', '160', 'read0', '60', '+']
['chr1', '520', '560', 'read1', '60', '+']
['chr1', '1500', '1540', 'read2', '60', '+']
['chr2', '20', '60', 'read3', '60', '+']
```

## 3) `dnakit.annotation.convert_intervals_to_bam_pybedtools()`

- **Purpose:** Write BED or BEDPE intervals as an unaligned BAM file.
- **Method:** `input_format="bed"` calls `BedTool.to_bam()` (`bedtools bedtobam`); `input_format="bedpe"` calls `BedTool.bedpe_to_bam()` (`bedpeToBam`). `mapq` sets the mapping quality, `bed12` turns BED12 blocks into CIGAR blocks (BED only) and `uncompressed` writes an uncompressed BAM.
- **API:** `dnakit.annotation.convert_intervals_to_bam_pybedtools(intervals_path: 'str | Path', genome_path: 'str | Path', output_path: 'str | Path', *, input_format: 'str' = 'bed', mapq: 'int' = 255, bed12: 'bool' = False, uncompressed: 'bool' = False)`.
- **Input:** An interval file, a chromosome-size file, a new output BAM path and the options above.
- **Output:** `output_path`, `n_alignments` (counted with pysam when installed) and `input_format`.
- **Example code:**

```python
from dnakit.annotation.pybedtools_p33_functions import convert_intervals_to_bam_pybedtools

result = convert_intervals_to_bam_pybedtools("a.bed", "genome.sizes", "a.bam", verbose=False)
print(result["input_format"], result["n_alignments"])
pe = convert_intervals_to_bam_pybedtools("pairs.bedpe", "genome.sizes", "pairs.bam", input_format="bedpe", verbose=False)
print(pe["input_format"], pe["n_alignments"])
```

- **Example output:**

```text
bed 3
bedpe 4
```

## 4) `dnakit.annotation.intervals_to_records_pybedtools()`

- **Purpose:** Parse BED, GFF, VCF or SAM records into named columns and typed values.
- **Method:** Call `BedTool.to_dataframe()`. Column names follow pybedtools (`chrom start end name score strand ...` for BED); `disable_auto_names=True` numbers the columns instead.
- **API:** `dnakit.annotation.intervals_to_records_pybedtools(path: 'str | Path', *, disable_auto_names: 'bool' = False, max_rows: 'int' = 1000)`.
- **Input:** An interval file path, optional `disable_auto_names` and `max_rows`.
- **Output:** `columns`, `dtypes`, `n_records` and `records` (dicts).
- **Example code:**

```python
from dnakit.annotation.pybedtools_p33_functions import intervals_to_records_pybedtools

result = intervals_to_records_pybedtools("a.bed", verbose=False)
print(result["columns"])
print(result["dtypes"]["start"], result["n_records"])
print(result["records"][0])
```

- **Example output:**

```text
['chrom', 'start', 'end', 'name', 'score', 'strand']
int64 3
{'chrom': 'chr1', 'start': 100, 'end': 200, 'name': 'a1', 'score': 5, 'strand': '+'}
```

## 5) `dnakit.annotation.list_tabix_contigs_pybedtools()`

- **Purpose:** List the contig names of a bgzip-compressed, tabix-indexed file.
- **Method:** Open the file with `pysam.TabixFile` (the `.tbi` or `.csi` index must be next to it), as `BedTool.tabix_contigs()` does.
- **API:** `dnakit.annotation.list_tabix_contigs_pybedtools(path: 'str | Path')`.
- **Input:** A bgzip file path with its index. Needs pysam.
- **Output:** `contigs` and `n_contigs`.
- **Example code:**

```python
from dnakit.annotation.pybedtools_p33_functions import list_tabix_contigs_pybedtools

import pysam

compressed = pysam.tabix_index("a.bed", preset="bed", force=True, keep_original=True)
result = list_tabix_contigs_pybedtools(compressed, verbose=False)
print(result["contigs"], result["n_contigs"])
```

- **Example output:**

```text
['chr1', 'chr2'] 2
```

## 6) `dnakit.annotation.detect_alignment_file_format_pybedtools()`

- **Purpose:** Detect BAM, CRAM, BGZF and gzip files from their magic bytes.
- **Method:** Call `pybedtools.helpers.isBAM`, `isCRAM`, `isBGZIP` and `isGZIP`. A BAM file is also BGZF-compressed, so `is_bgzip` is true for BAM; `format` names the most specific match.
- **API:** `dnakit.annotation.detect_alignment_file_format_pybedtools(path: 'str | Path')`.
- **Input:** A file path.
- **Output:** `is_bam`, `is_cram`, `is_bgzip`, `is_gzip` and `format` (`bam`, `cram`, `bgzf`, `gzip` or `other`).
- **Example code:**

```python
from dnakit.annotation.pybedtools_p33_functions import detect_alignment_file_format_pybedtools

exec(open("make_bam.py").read())
for path in ("reads.bam", "a.bed"):
    result = detect_alignment_file_format_pybedtools(path, verbose=False)
    print(path, result["format"], result["is_bam"], result["is_bgzip"])
```

- **Example output:**

```text
reads.bam bam True True
a.bed other False False
```

## 7) `dnakit.annotation.interval_has_hits_pybedtools()`

- **Purpose:** Test whether one query interval overlaps any interval of a file.
- **Method:** Build a pybedtools interval from `chrom`, `start`, `end` and `strand`, then call `BedTool.any_hits()`. `same_strand` requires equal strands and `min_overlap_fraction` the fraction of the query that must overlap.
- **API:** `dnakit.annotation.interval_has_hits_pybedtools(intervals_path: 'str | Path', chrom: 'str', start: 'int', end: 'int', *, strand: 'str' = '.', same_strand: 'bool' = False, min_overlap_fraction: 'float' = 0.0)`.
- **Input:** An interval file path, the query coordinates (BED style) and the options above.
- **Output:** `has_hits` and the `query`.
- **Example code:**

```python
from dnakit.annotation.pybedtools_p33_functions import interval_has_hits_pybedtools

print(interval_has_hits_pybedtools("a.bed", "chr1", 150, 160, verbose=False)["has_hits"])
print(interval_has_hits_pybedtools("a.bed", "chr1", 300, 400, verbose=False)["has_hits"])
print(interval_has_hits_pybedtools("a.bed", "chr1", 150, 160, strand="-", same_strand=True, verbose=False)["has_hits"])
```

- **Example output:**

```text
True
False
False
```

## 8) `dnakit.annotation.slop_pybedtools_intervals()`

- **Purpose:** Extend intervals by flanking bases, clipped at the chromosome ends.
- **Method:** Call `BedTool.slop()` (`bedtools slop`). `both` extends both sides; `left` and `right` (given together) extend each side, upstream and downstream when `strand_aware`. `as_fraction` reads the values as fractions of the interval length.
- **API:** `dnakit.annotation.slop_pybedtools_intervals(bed_path: 'str | Path', genome_path: 'str | Path', *, both: 'float | None' = None, left: 'float | None' = None, right: 'float | None' = None, strand_aware: 'bool' = False, as_fraction: 'bool' = False, header: 'bool' = False, output_path: 'str | Path | None' = None, max_rows: 'int' = 1000)`.
- **Input:** A BED file, a chromosome-size file and the extension options.
- **Output:** `n_records` and `records`.
- **Example code:**

```python
from dnakit.annotation.pybedtools_p33_functions import slop_pybedtools_intervals

result = slop_pybedtools_intervals("a.bed", "genome.sizes", both=20, verbose=False)
for record in result["records"]:
    print(record[:3])
stranded = slop_pybedtools_intervals("a.bed", "genome.sizes", left=10, right=0, strand_aware=True, verbose=False)
print(stranded["records"][1][:3])
```

- **Example output:**

```text
['chr1', '80', '220']
['chr1', '480', '620']
['chr2', '0', '110']
['chr1', '500', '610']
```

## 9) `dnakit.annotation.shift_pybedtools_intervals()`

- **Purpose:** Move intervals along the chromosome, clipped at the chromosome ends.
- **Method:** Call `BedTool.shift()` (`bedtools shift`). `shift` moves every interval; `plus_shift` and `minus_shift` move plus and minus strand intervals separately. `as_fraction` reads shifts as fractions of the interval length.
- **API:** `dnakit.annotation.shift_pybedtools_intervals(bed_path: 'str | Path', genome_path: 'str | Path', *, shift: 'float | None' = None, plus_shift: 'float | None' = None, minus_shift: 'float | None' = None, as_fraction: 'bool' = False, header: 'bool' = False, output_path: 'str | Path | None' = None, max_rows: 'int' = 1000)`.
- **Input:** A BED file, a chromosome-size file and the shift options.
- **Output:** `n_records` and `records`.
- **Example code:**

```python
from dnakit.annotation.pybedtools_p33_functions import shift_pybedtools_intervals

result = shift_pybedtools_intervals("a.bed", "genome.sizes", shift=30, verbose=False)
for record in result["records"]:
    print(record[:3])
split = shift_pybedtools_intervals("a.bed", "genome.sizes", plus_shift=5, minus_shift=-5, verbose=False)
print([r[1] for r in split["records"]])
```

- **Example output:**

```text
['chr1', '130', '230']
['chr1', '530', '630']
['chr2', '40', '120']
['105', '495', '15']
```

## 10) `dnakit.annotation.make_windows_pybedtools()`

- **Purpose:** Cut a genome or BED intervals into fixed-size or fixed-number windows.
- **Method:** Call `BedTool.window_maker()` (`bedtools makewindows`). `source_type` is `genome` (chromosome-size file) or `bed`. Give `window_size` (with optional `step_size` for sliding windows) or `number_of_windows`; `id_naming` is `winnum` or `srcwinnum`.
- **API:** `dnakit.annotation.make_windows_pybedtools(source_path: 'str | Path', *, source_type: 'str' = 'genome', window_size: 'int | None' = None, step_size: 'int | None' = None, number_of_windows: 'int | None' = None, reverse_numbering: 'bool' = False, id_naming: 'str | None' = None, output_path: 'str | Path | None' = None, max_rows: 'int' = 1000)`.
- **Input:** A chromosome-size file or BED file and the window options.
- **Output:** `n_records` and `records`.
- **Example code:**

```python
from dnakit.annotation.pybedtools_p33_functions import make_windows_pybedtools

result = make_windows_pybedtools("genome.sizes", window_size=700, verbose=False)
print(result["n_records"], result["records"][:2])
bed = make_windows_pybedtools("a.bed", source_type="bed", number_of_windows=2, id_naming="winnum", verbose=False)
print(bed["records"][:2])
```

- **Example output:**

```text
5 [['chr1', '0', '700'], ['chr1', '700', '1400']]
[['chr1', '100', '150', '1'], ['chr1', '150', '200', '2']]
```

## 11) `dnakit.annotation.cat_pybedtools_intervals()`

- **Purpose:** Concatenate interval files and, by default, merge the result.
- **Method:** Call `BedTool.cat()`. With `postmerge` the files are concatenated, sorted and merged (`merge_distance` is the maximum gap, `same_strand` merges per strand) and only the first three columns are kept; without it the files are concatenated as they are.
- **API:** `dnakit.annotation.cat_pybedtools_intervals(paths: 'Sequence[str | Path]', *, postmerge: 'bool' = True, force_truncate: 'bool' = False, merge_distance: 'int' = 0, same_strand: 'bool' = False, output_path: 'str | Path | None' = None, max_rows: 'int' = 1000)`.
- **Input:** Two or more interval files and the options above.
- **Output:** `n_records` and `records`.
- **Example code:**

```python
from dnakit.annotation.pybedtools_p33_functions import cat_pybedtools_intervals

merged = cat_pybedtools_intervals(["a.bed", "b.bed"], verbose=False)
print(merged["records"])
plain = cat_pybedtools_intervals(["a.bed", "b.bed"], postmerge=False, verbose=False)
print(plain["n_records"])
```

- **Example output:**

```text
[['chr1', '100', '250'], ['chr1', '500', '600'], ['chr1', '900', '950'], ['chr2', '10', '90'], ['chr2', '100', '150']]
6
```

## 12) `dnakit.annotation.truncate_pybedtools_to_chrom()`

- **Purpose:** Clip intervals to the chromosome bounds of a size file.
- **Method:** Call `BedTool.truncate_to_chrom()` with the sizes of the file as `(0, length)` bounds. Intervals beyond a chromosome end are cut and intervals on unknown chromosomes are dropped.
- **API:** `dnakit.annotation.truncate_pybedtools_to_chrom(bed_path: 'str | Path', genome_path: 'str | Path', *, output_path: 'str | Path | None' = None, max_rows: 'int' = 1000)`.
- **Input:** A BED file and a chromosome-size file.
- **Output:** `n_records` and `records`.
- **Example code:**

```python
from dnakit.annotation.pybedtools_p33_functions import truncate_pybedtools_to_chrom

result = truncate_pybedtools_to_chrom("wide.bed", "genome.sizes", verbose=False)
for record in result["records"]:
    print(record)
```

- **Example output:**

```text
['chr1', '1900', '2000', 'w1']
['chr2', '5', '50', 'w3']
```

## 13) `dnakit.annotation.remove_invalid_pybedtools_intervals()`

- **Purpose:** Drop records with invalid coordinates.
- **Method:** Call `BedTool.remove_invalid()`, which removes records whose start is greater than the end, negative coordinates and malformed lines. `n_removed` counts the dropped data lines.
- **API:** `dnakit.annotation.remove_invalid_pybedtools_intervals(bed_path: 'str | Path', *, output_path: 'str | Path | None' = None, max_rows: 'int' = 1000)`.
- **Input:** A BED file path.
- **Output:** `n_records`, `records` and `n_removed`.
- **Example code:**

```python
from dnakit.annotation.pybedtools_p33_functions import remove_invalid_pybedtools_intervals

result = remove_invalid_pybedtools_intervals("bad.bed", verbose=False)
print(result["n_records"], result["n_removed"])
print(result["records"])
```

- **Example output:**

```text
2 1
[['chr1', '10', '100'], ['chr1', '5', '6']]
```

## 14) `dnakit.annotation.introns_pybedtools()`

- **Purpose:** Derive intron and UTR intervals from GFF gene and exon features, or from BED12.
- **Method:** Call `BedTool.introns()`: every gene is paired with the same-strand exons inside it; gaps between exons become `intron`, the parts before the first and after the last exon become `utr5` or `utr3`. Overlapping genes can give records with start greater than end; those are dropped and counted in `n_invalid_dropped`.
- **API:** `dnakit.annotation.introns_pybedtools(annotation_path: 'str | Path', *, gene: 'str' = 'gene', exon: 'str' = 'exon', output_path: 'str | Path | None' = None, max_rows: 'int' = 1000)`.
- **Input:** A GFF file (or BED12) and the feature types `gene` and `exon`.
- **Output:** `n_records`, `records` (BED6 with the gene id and the kind as name) and `n_invalid_dropped`.
- **Example code:**

```python
from dnakit.annotation.pybedtools_p33_functions import introns_pybedtools

result = introns_pybedtools("genes.gff", verbose=False)
for record in result["records"]:
    print(record)
print(result["n_invalid_dropped"])
```

- **Example output:**

```text
['chr1', '200', '300', 'g1', 'intron', '+']
['chr1', '400', '500', 'g1', 'intron', '+']
0
```

## 15) `dnakit.annotation.tag_bam_pybedtools()`

- **Purpose:** Add a BAM tag naming the interval files (or interval names) that each alignment overlaps.
- **Method:** Call `BedTool.tag_bam()` (`bedtools tag`). `labels` gives one label per interval file; `report` selects labels, interval `names`, `scores` or `intervals` (chrom:start-end) as the tag value. The tagged alignments are written to a new BAM.
- **API:** `dnakit.annotation.tag_bam_pybedtools(bam_path: 'str | Path', interval_files: 'Sequence[str | Path]', output_path: 'str | Path', *, labels: 'Sequence[str] | None' = None, tag: 'str' = 'YB', same_strand: 'bool' = False, opposite_strand: 'bool' = False, min_overlap_fraction: 'float | None' = None, report: 'str' = 'labels')`.
- **Input:** A BAM file, a list of interval files, a new output BAM path and the tag options.
- **Output:** `output_path`, `n_alignments`, `tag` and `labels`.
- **Example code:**

```python
from dnakit.annotation.pybedtools_p33_functions import tag_bam_pybedtools

import pysam

exec(open("make_bam.py").read())
result = tag_bam_pybedtools("reads.bam", ["a.bed"], "tagged.bam", labels=["A"], tag="XA", verbose=False)
print(result["n_alignments"], result["tag"])
with pysam.AlignmentFile("tagged.bam", "rb") as handle:
    print([(read.query_name, read.get_tag("XA") if read.has_tag("XA") else None) for read in handle])
```

- **Example output:**

```text
4 XA
[('read0', 'A'), ('read1', 'A'), ('read2', None), ('read3', 'A')]
```
