# pybedtools: interval statistics, randomisation and reports

This group in `dnakit.annotation` wraps 14 further operations of [pybedtools](https://daler.github.io/pybedtools/) on top of [BEDTools](https://bedtools.readthedocs.io/): `map`, `annotate`, genome coverage, total covered bases, BEDPE pairs against BED, spacing, Fisher's exact test, nucleotide content, the randomised overlap distribution (count or base pairs), its significance summary and pairwise matrix, a feature-track plot, an IGV batch script and a UCSC link page.

The same requirements as in the previous document apply (BEDTools on `PATH`, `dnakit[bio-pybedtools]`; tested with BEDTools 2.27.1 and pybedtools 0.12.1). Inputs of `map`, `spacing` and `fisher` are sorted first because BEDTools needs sorted input. The randomisation functions take a seed: shuffle number `i` uses `bedtools shuffle -seed seed+i`, which reproduces pybedtools `debug=True` for `seed=0` and makes the result repeatable. The plot function draws on a matplotlib `Figure` (`dnakit[bio-plots]`) and returns the image path plus the drawn intervals. Nothing contacts the network; the UCSC link page and the IGV script only contain text.

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

**Example file:** `pairs.bedpe`

```text
chr1	100	150	chr1	300	350	p1	5	+	-
chr2	10	60	chr2	200	250	p2	6	+	+
```

**Example file:** `ref.fa`

```text
>chr1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
>chr2
AACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTT
```

## 1) `dnakit.annotation.map_pybedtools_values()`

- **Purpose:** Summarise a column of B over the B intervals that overlap each A interval.
- **Method:** Call `BedTool.map()` (`bedtools map`) after sorting A and B. `operation` is one of sum, min, max, absmin, absmax, mean, median, mode, antimode, stdev, sstdev, collapse, distinct, distinct_sort_num, distinct_sort_num_desc, count, count_distinct, first or last; `null` is written for A intervals without overlap.
- **API:** `dnakit.annotation.map_pybedtools_values(a_path: 'str | Path', b_path: 'str | Path', *, column: 'int' = 5, operation: 'str' = 'sum', null: 'str' = '.', delimiter: 'str' = ',', min_overlap_a: 'float | None' = None, min_overlap_b: 'float | None' = None, reciprocal: 'bool' = False, same_strand: 'bool' = False, opposite_strand: 'bool' = False, split: 'bool' = False, output_path: 'str | Path | None' = None, max_rows: 'int' = 1000)`.
- **Input:** Two interval files, the 1-based `column` of B, the `operation` and overlap options (`min_overlap_a`, `min_overlap_b`, `reciprocal`, `same_strand`, `opposite_strand`, `split`).
- **Output:** `n_records` and `records`: the A fields followed by the summary value.
- **Example code:**

```python
from dnakit.annotation.pybedtools_p33_functions import map_pybedtools_values

result = map_pybedtools_values("a.bed", "b.bed", column=5, operation="sum", verbose=False)
for record in result["records"]:
    print(record)
```

- **Example output:**

```text
['chr1', '100', '200', 'a1', '5', '+', '1']
['chr1', '500', '600', 'a2', '6', '-', '.']
['chr2', '10', '90', 'a3', '7', '+', '.']
```

## 2) `dnakit.annotation.annotate_pybedtools_intervals()`

- **Purpose:** Report how much of each A interval is covered by features of several files.
- **Method:** Call `BedTool.annotate()` (`bedtools annotate`). One column group per feature file: the fraction covered (default), overlap counts (`counts`) or both (`both`). `names` labels the files.
- **API:** `dnakit.annotation.annotate_pybedtools_intervals(a_path: 'str | Path', feature_files: 'Sequence[str | Path]', *, names: 'Sequence[str] | None' = None, counts: 'bool' = False, both: 'bool' = False, same_strand: 'bool' = False, opposite_strand: 'bool' = False, output_path: 'str | Path | None' = None, max_rows: 'int' = 1000)`.
- **Input:** An A file, a list of feature files and the strand and report options.
- **Output:** `n_records` and `records`: the A fields followed by one value per feature file.
- **Example code:**

```python
from dnakit.annotation.pybedtools_p33_functions import annotate_pybedtools_intervals

result = annotate_pybedtools_intervals("a.bed", ["b.bed"], counts=True, verbose=False)
for record in result["records"]:
    print(record)
```

- **Example output:**

```text
['chr1', '100', '200', 'a1', '5', '+', '1']
['chr1', '500', '600', 'a2', '6', '-', '0']
['chr2', '10', '90', 'a3', '7', '+', '0']
```

## 3) `dnakit.annotation.genome_coverage_pybedtools()`

- **Purpose:** Compute the coverage of a genome by BED intervals or BAM alignments.
- **Method:** Call `BedTool.genome_coverage()` (`bedtools genomecov`). `report` is `histogram` (default), `bedgraph`, `bedgraph_all`, `per_base` or `per_base_nonzero`. Interval input is sorted first; BAM input must be coordinate-sorted and needs no size file.
- **API:** `dnakit.annotation.genome_coverage_pybedtools(path: 'str | Path', genome_path: 'str | Path | None' = None, *, report: 'str' = 'histogram', scale: 'float' = 1.0, strand: 'str | None' = None, five_prime: 'bool' = False, three_prime: 'bool' = False, split: 'bool' = False, max_depth: 'int | None' = None, output_path: 'str | Path | None' = None, max_rows: 'int' = 1000)`.
- **Input:** A BED or BAM file, a chromosome-size file for BED input, and options (`scale`, `strand`, `five_prime`, `three_prime`, `split`, `max_depth`).
- **Output:** `n_records`, `records` and `report`.
- **Example code:**

```python
from dnakit.annotation.pybedtools_p33_functions import genome_coverage_pybedtools

graph = genome_coverage_pybedtools("a.bed", "genome.sizes", report="bedgraph", verbose=False)
print(graph["records"])
hist = genome_coverage_pybedtools("a.bed", "genome.sizes", verbose=False)
print(hist["records"][:3])
```

- **Example output:**

```text
[['chr1', '100', '200', '1'], ['chr1', '500', '600', '1'], ['chr2', '10', '90', '1']]
[['chr1', '0', '1800', '2000', '0.9'], ['chr1', '1', '200', '2000', '0.1'], ['chr2', '0', '920', '1000', '0.92']]
```

## 4) `dnakit.annotation.total_coverage_pybedtools()`

- **Purpose:** Count the bases covered by at least one interval.
- **Method:** Call `BedTool.total_coverage()`, which merges the intervals and adds up the merged lengths.
- **API:** `dnakit.annotation.total_coverage_pybedtools(path: 'str | Path')`.
- **Input:** An interval file path.
- **Output:** `total_bases` and `n_intervals`.
- **Example code:**

```python
from dnakit.annotation.pybedtools_p33_functions import total_coverage_pybedtools

result = total_coverage_pybedtools("a.bed", verbose=False)
print(result["total_bases"], result["n_intervals"])
```

- **Example output:**

```text
280 3
```

## 5) `dnakit.annotation.pair_to_bed_pybedtools()`

- **Purpose:** Keep the BEDPE pairs that overlap BED intervals.
- **Method:** Call `BedTool.pair_to_bed()` (`bedtools pairtobed`). `pair_type` is either, neither, xor, both, notboth, ispan, ospan, notispan or notospan; `min_overlap_fraction`, `same_strand` and `opposite_strand` refine the overlap.
- **API:** `dnakit.annotation.pair_to_bed_pybedtools(bedpe_path: 'str | Path', bed_path: 'str | Path', *, pair_type: 'str' = 'either', min_overlap_fraction: 'float | None' = None, same_strand: 'bool' = False, opposite_strand: 'bool' = False, output_path: 'str | Path | None' = None, max_rows: 'int' = 1000)`.
- **Input:** A BEDPE file, a BED file and the options above.
- **Output:** `n_records` and `records` (BEDPE rows).
- **Example code:**

```python
from dnakit.annotation.pybedtools_p33_functions import pair_to_bed_pybedtools

hit = pair_to_bed_pybedtools("pairs.bedpe", "a.bed", verbose=False)
print(hit["n_records"], hit["records"][0][6])
none = pair_to_bed_pybedtools("pairs.bedpe", "b.bed", pair_type="neither", verbose=False)
print(none["n_records"])
```

- **Example output:**

```text
2 p1
2
```

## 6) `dnakit.annotation.spacing_pybedtools_intervals()`

- **Purpose:** Report the gap to the previous interval for every interval.
- **Method:** Sort the input, then call `BedTool.spacing()` (`bedtools spacing`). The first interval of a chromosome gets `.` and overlapping or touching intervals get a negative or zero value as BEDTools writes them.
- **API:** `dnakit.annotation.spacing_pybedtools_intervals(bed_path: 'str | Path', *, output_path: 'str | Path | None' = None, max_rows: 'int' = 1000)`.
- **Input:** A BED file path.
- **Output:** `n_records` and `records`: the BED fields followed by the spacing.
- **Example code:**

```python
from dnakit.annotation.pybedtools_p33_functions import spacing_pybedtools_intervals

result = spacing_pybedtools_intervals("a.bed", verbose=False)
print([record[-1] for record in result["records"]])
```

- **Example output:**

```text
['.', '300', '.']
```

## 7) `dnakit.annotation.fisher_pybedtools_overlap()`

- **Purpose:** Test with Fisher's exact test whether two interval sets overlap more than expected by chance.
- **Method:** Sort A and B, then call `BedTool.fisher()` (`bedtools fisher`) and parse its report. The number of possible intervals is estimated from the genome size file.
- **API:** `dnakit.annotation.fisher_pybedtools_overlap(a_path: 'str | Path', b_path: 'str | Path', genome_path: 'str | Path', *, min_overlap_a: 'float | None' = None, reciprocal: 'bool' = False, same_strand: 'bool' = False, opposite_strand: 'bool' = False, split: 'bool' = False)`.
- **Input:** Two interval files, a chromosome-size file and overlap options.
- **Output:** `table` (2x2 counts), `left_tail`, `right_tail`, `two_tail`, `ratio` and the raw `report`.
- **Example code:**

```python
from dnakit.annotation.pybedtools_p33_functions import fisher_pybedtools_overlap

result = fisher_pybedtools_overlap("a.bed", "b.bed", "genome.sizes", verbose=False)
print(result["table"]["in -a"]["in -b"], result["table"]["not in -a"]["not in -b"])
print(result["two_tail"] < 0.001, result["right_tail"] < 0.001)
```

- **Example output:**

```text
1 13
False False
```

## 8) `dnakit.annotation.nucleotide_content_pybedtools()`

- **Purpose:** Profile the nucleotide composition of intervals in a FASTA file.
- **Method:** Call `BedTool.nucleotide_content()` (`bedtools nuc`). Each interval gets AT and GC fractions, base counts (A, C, G, T, N, other) and the sequence length; `pattern` adds a pattern count, `include_sequence` the sequence, `strand_aware` reverse-complements minus-strand intervals first.
- **API:** `dnakit.annotation.nucleotide_content_pybedtools(bed_path: 'str | Path', fasta_path: 'str | Path', *, strand_aware: 'bool' = False, pattern: 'str | None' = None, case_insensitive: 'bool' = False, include_sequence: 'bool' = False, output_path: 'str | Path | None' = None, max_rows: 'int' = 1000)`.
- **Input:** A BED file, a FASTA file and the options above.
- **Output:** `header` (the column names written by BEDTools), `n_records` and `records`.
- **Example code:**

```python
from dnakit.annotation.pybedtools_p33_functions import nucleotide_content_pybedtools

result = nucleotide_content_pybedtools("a.bed", "ref.fa", verbose=False)
print(result["header"][6:])
print(result["records"][0][6:])
```

- **Example output:**

```text
['7_pct_at', '8_pct_gc', '9_num_A', '10_num_C', '11_num_G', '12_num_T', '13_num_N', '14_num_oth', '15_seq_len']
['0.500000', '0.500000', '25', '25', '25', '25', '0', '0', '100']
```

## 9) `dnakit.annotation.random_overlap_distribution_pybedtools()`

- **Purpose:** Compute the distribution of overlaps between randomly shuffled A and fixed B intervals.
- **Method:** Shuffle A within the genome `iterations` times with `bedtools shuffle -seed seed+i` and intersect with B. `metric="count"` counts overlapping A intervals (`BedTool.randomintersection`); `metric="bp"` sums the overlapping base pairs (`BedTool.randomintersection_bp`). The overlap of the unshuffled A is returned for comparison.
- **API:** `dnakit.annotation.random_overlap_distribution_pybedtools(a_path: 'str | Path', b_path: 'str | Path', genome_path: 'str | Path', *, iterations: 'int' = 100, metric: 'str' = 'count', seed: 'int' = 0)`.
- **Input:** Two interval files, a chromosome-size file, `iterations` (up to 10000), `metric` and `seed`.
- **Output:** `actual`, `distribution` (one value per iteration), `metric`, `iterations` and `seed`.
- **Example code:**

```python
from dnakit.annotation.pybedtools_p33_functions import random_overlap_distribution_pybedtools

count = random_overlap_distribution_pybedtools("a.bed", "b.bed", "genome.sizes", iterations=5, verbose=False)
print(count["actual"], count["distribution"])
bp = random_overlap_distribution_pybedtools("a.bed", "b.bed", "genome.sizes", iterations=5, metric="bp", verbose=False)
print(bp["actual"], bp["distribution"])
```

- **Example output:**

```text
1 [0, 1, 2, 0, 1]
50 [0, 50, 25, 0, 84]
```

## 10) `dnakit.annotation.random_overlap_summary_pybedtools()`

- **Purpose:** Assess the significance of the overlap of A with B by shuffling A.
- **Method:** Compute the randomised overlap counts as above and summarise them as `BedTool.randomstats()` does: the actual count, the median of the randomised counts, the fractions of randomised counts above and below the actual one, the actual count divided by the median, the percentile of the actual count and the 2.5th and 97.5th percentiles.
- **API:** `dnakit.annotation.random_overlap_summary_pybedtools(a_path: 'str | Path', b_path: 'str | Path', genome_path: 'str | Path', *, iterations: 'int' = 100, seed: 'int' = 0, include_distribution: 'bool' = False)`.
- **Input:** Two interval files, a chromosome-size file, `iterations`, `seed` and `include_distribution`.
- **Output:** `actual`, `self`, `other`, `frac_randomized_above_actual`, `frac_randomized_below_actual`, `median_randomized`, `normalized`, `percentile`, `lower_2.5th`, `upper_97.5th` (and `distribution`).
- **Example code:**

```python
from dnakit.annotation.pybedtools_p33_functions import random_overlap_summary_pybedtools

result = random_overlap_summary_pybedtools("a.bed", "b.bed", "genome.sizes", iterations=20, verbose=False)
print(result["actual"], result["self"], result["other"])
print(result["median_randomized"], result["percentile"])
print(round(result["normalized"], 3))
```

- **Example output:**

```text
1 3 3
1.0 72.5
1.0
```

## 11) `dnakit.annotation.random_overlap_matrix_pybedtools()`

- **Purpose:** Compute randomisation overlap statistics for every ordered pair of several interval files.
- **Method:** For each pair (A, B) shuffle A `iterations` times and summarise as in `random_overlap_summary_pybedtools`, like `IntersectionMatrix.create_matrix()`. File names without extension label the matrix and must be unique.
- **API:** `dnakit.annotation.random_overlap_matrix_pybedtools(bed_paths: 'Sequence[str | Path]', genome_path: 'str | Path', *, iterations: 'int' = 100, seed: 'int' = 0)`.
- **Input:** 2 to 12 interval files, a chromosome-size file, `iterations` and `seed`.
- **Output:** `names` and `matrix`, a mapping `name_a -> name_b -> summary`.
- **Example code:**

```python
from dnakit.annotation.pybedtools_p33_functions import random_overlap_matrix_pybedtools

result = random_overlap_matrix_pybedtools(["a.bed", "b.bed"], "genome.sizes", iterations=5, verbose=False)
print(result["names"])
for first in result["names"]:
    print(first, [result["matrix"][first][second]["actual"] for second in result["names"]])
```

- **Example output:**

```text
['a', 'b']
a [3, 1]
b [1, 3]
```

## 12) `dnakit.annotation.plot_interval_tracks_pybedtools()`

- **Purpose:** Draw several BED files as stacked feature tracks over one region and save a PNG.
- **Method:** Read the intervals of each file that fall in `chrom:start-end` and draw them as bars on a matplotlib `Figure`, like `TrackCollection.plot()`. `label_intervals` writes the BED names above the bars.
- **API:** `dnakit.annotation.plot_interval_tracks_pybedtools(tracks: 'Mapping[str, str | Path]', chrom: 'str', start: 'int', end: 'int', output_path: 'str | Path', *, label_intervals: 'bool' = False)`.
- **Input:** A mapping of track label to BED file, the region, a new PNG path and `label_intervals`.
- **Output:** `image_path`, `region` and `tracks` (the intervals drawn per track).
- **Example code:**

```python
from dnakit.annotation.pybedtools_p33_functions import plot_interval_tracks_pybedtools

result = plot_interval_tracks_pybedtools({"a": "a.bed", "b": "b.bed"}, "chr1", 0, 1000, "tracks.png", verbose=False)
print(open(result["image_path"], "rb").read(4))
print({label: len(items) for label, items in result["tracks"].items()})
```

- **Example output:**

```text
b'\x89PNG'
{'a': 2, 'b': 2}
```

## 13) `dnakit.annotation.build_igv_script_pybedtools()`

- **Purpose:** Create an IGV batch script that takes a snapshot of every interval.
- **Method:** Call `BedTool.igv()` (`bedtools igv`) and return the script text. `slop` pads each locus, `sort` and `collapse` set the display, `name_in_file` names snapshots by BED name. IGV is not started.
- **API:** `dnakit.annotation.build_igv_script_pybedtools(bed_path: 'str | Path', *, snapshot_dir: 'str' = '/tmp/snapshots', image_format: 'str' = 'png', slop: 'int' = 0, sort: 'str | None' = None, collapse: 'bool' = False, name_in_file: 'bool' = False, output_path: 'str | Path | None' = None)`.
- **Input:** A BED file and the snapshot options (`snapshot_dir`, `image_format`, `slop`, ...).
- **Output:** `script`, `n_commands` and `output_path`.
- **Example code:**

```python
from dnakit.annotation.pybedtools_p33_functions import build_igv_script_pybedtools

result = build_igv_script_pybedtools("a.bed", snapshot_dir="/snaps", slop=5, verbose=False)
print(result["n_commands"])
print(result["script"].splitlines()[:3])
```

- **Example output:**

```text
3
['snapshotDirectory /snaps', 'goto chr1:95-205', 'snapshot chr1_100_200_slop5.png']
```

## 14) `dnakit.annotation.build_ucsc_links_page_pybedtools()`

- **Purpose:** Write an HTML table that links each interval to its UCSC Genome Browser position.
- **Method:** Run `bedtools links` (as `BedTool.links()`). The page is built offline; the links point at `base_url`/cgi-bin/hgTracks for `organism` and `database`.
- **API:** `dnakit.annotation.build_ucsc_links_page_pybedtools(bed_path: 'str | Path', *, base_url: 'str' = 'http://genome.ucsc.edu', organism: 'str' = 'human', database: 'str' = 'hg19', output_path: 'str | Path | None' = None)`.
- **Input:** A BED file, `base_url` (http or https), `organism` and `database`.
- **Output:** `html`, `n_links` and `output_path`.
- **Example code:**

```python
from dnakit.annotation.pybedtools_p33_functions import build_ucsc_links_page_pybedtools

result = build_ucsc_links_page_pybedtools("a.bed", database="hg38", verbose=False)
print(result["n_links"])
print([line.strip() for line in result["html"].splitlines() if "<a href=" in line][0])
```

- **Example output:**

```text
3
<a href=http://genome.ucsc.edu/cgi-bin/hgTracks?org=human&db=hg38&position=chr1:100-200>chr1:101-200</a>
```
