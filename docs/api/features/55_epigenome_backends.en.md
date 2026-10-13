# Optional epigenome backends: contact matrices, pair statistics, coverage and methylation

This group of functions in `dnakit.epigenome` is a set of adapters over optional packages: [cooler](https://github.com/open2c/cooler), [cooltools](https://github.com/open2c/cooltools), [pairtools](https://github.com/mirnylab/pairtools), [MACS3](https://github.com/macs3-project/MACS), [methylpy](https://github.com/yupenghe/methylpy) and [deeptools](https://github.com/deeptools/deepTools). They work only on explicit files or numbers. None of them aligns raw reads, queries a reference database, or runs a complete upstream pipeline (for example, `cooler cload` from BAM input, a bisulfite aligner, or MACS3's `callpeak` statistical model). The dependencies are optional and are imported only when the matching function is called. If a package is missing, the function raises an error that includes the install command.

## Example data

The examples below use these small files. `tiny.pairs` has 7 contacts that fall into 6 bin pairs, with 1-based positions. `typed.pairs` has a `pair_type` column. `dup.pairs` contains one group of PCR duplicates. `tiny.allc` has 6 cytosine sites. The files only demonstrate the interface and their numbers do not describe real data. Binary inputs (cooler, BAM and bigWig files) are created by the example code itself.

**Example file:** `tiny.pairs`

```text
## pairs format v1.0
#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2
r1	chr1	10	chr1	150	+	-
r2	chr1	20	chr1	160	+	-
r3	chr1	15	chr1	250	+	-
r4	chr1	510	chr1	520	+	-
r5	chr1	120	chr2	30	+	-
r6	chr2	40	chr2	90	+	-
r7	chr2	200	chr2	420	-	+
```

**Example file:** `typed.pairs`

```text
## pairs format v1.0
#sorted: chr1-pos1-chr2-pos2
#chromsize: chr1 1000
#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type
a	chr1	10	chr1	150	+	-	UU
b	chr1	20	chr1	300	+	+	UU
c	chr1	500	chr1	510	-	+	UR
d	chr1	600	chr2	40	+	-	UU
e	!	0	chr1	700	+	-	NU
```

**Example file:** `dup.pairs`

```text
## pairs format v1.0
#sorted: chr1-pos1-chr2-pos2
#chromsize: chr1 1000
#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type
a	chr1	10	chr1	150	+	-	UU
b	chr1	11	chr1	151	+	-	UU
c	chr1	14	chr1	152	+	-	UU
d	chr1	300	chr1	400	+	-	UU
e	chr1	300	chr1	400	-	+	UU
f	chr1	500	chr1	600	+	-	UU
g	!	0	chr1	700	+	-	NU
```

**Example file:** `tiny.allc`

```text
chrL	1	+	CGN	0	10	0
chrL	2	+	CGG	1	10	0
chrL	3	-	CHH	0	4	0
chr1	5	+	CGN	2	2	1
chr1	9	+	CHG	1	2	1
chr1	12	+	CGN	3	10	1
```

## 1) `dnakit.epigenome.build_contact_matrix()`

- **Purpose:** Bin a `.pairs` contact list into a new `.cool` matrix.
- **Method:** Read columns 2 to 5 (`chrom1`, `pos1`, `chrom2`, `pos2`); positions are 1-based. cooler's pairs pipeline bins them at `binsize` and sums counts for each bin pair, and the lower triangle is reflected into the upper triangle. The file is read in chunks. An existing `cool_path` is never overwritten.
- **API:** `dnakit.epigenome.build_contact_matrix(pairs_path, cool_path, *, chromsizes, binsize, max_bins=100000, chunksize=1000000)`.
- **Input:** A pairs file path; a new `cool_path`; `chromsizes` as `{chromosome: length}`; `binsize` in bp; `max_bins` limits the total number of bins.
- **Output:** A dictionary with `cool_path`, `binsize`, `chromosomes`, `bins`, `nonzero_pixels`, `total_contacts` and `provenance`.
- **Example code:**

```python
import tempfile
from pathlib import Path
from dnakit.epigenome import build_contact_matrix

work = Path(tempfile.mkdtemp())
result = build_contact_matrix("tiny.pairs", work / "tiny.cool", chromsizes={"chr1": 1000, "chr2": 600}, binsize=100)
print(result["nonzero_pixels"], result["total_contacts"], result["bins"])
```

- **Example output:**

```text
6 7 16
```

## 2) `dnakit.epigenome.fetch_contact_matrix()`

- **Purpose:** Read the dense contact matrix for one region.
- **Method:** Call cooler's `Cooler.matrix(balance=...).fetch(region)`. With `balance=True`, the stored weights are used; if there is no weight column the call fails. A region spanning more bins than `max_region_bins` is rejected, which bounds the output size.
- **API:** `dnakit.epigenome.fetch_contact_matrix(cool_path, region, *, balance=False, max_region_bins=200)`.
- **Input:** A `.cool` file path; a UCSC-style region such as `chr1` or `chr1:0-300`; whether to use balanced weights.
- **Output:** `matrix` (list of lists; missing values are `null`), `bins` (chromosome, start and end of each row), `binsize` and `balanced`.
- **Example code:**

```python
import tempfile
from pathlib import Path
from dnakit.epigenome import build_contact_matrix, fetch_contact_matrix

work = Path(tempfile.mkdtemp())
cool = build_contact_matrix("tiny.pairs", work / "tiny.cool", chromsizes={"chr1": 1000, "chr2": 600}, binsize=100)["cool_path"]
block = fetch_contact_matrix(cool, "chr1:0-300")
print(block["matrix"])
```

- **Example output:**

```text
[[0.0, 2.0, 1.0], [2.0, 0.0, 0.0], [1.0, 0.0, 0.0]]
```

## 3) `dnakit.epigenome.balance_contact_matrix()`

- **Purpose:** Compute iterative-correction weights and store them as the `weight` column of a `.cool` file.
- **Method:** Call cooler's `balance_cooler`. The weights are computed first and checked for finite values. If every bin is filtered out, the call fails and the file is not changed. An existing `weight` column is replaced.
- **API:** `dnakit.epigenome.balance_contact_matrix(cool_path, *, cis_only=False, ignore_diags=2, min_nnz=10, min_count=0, mad_max=5, tol=1e-5, max_iters=200)`.
- **Input:** A `.cool` file path and cooler's filter and convergence parameters. `min_nnz` and `min_count` decide which bins take part in balancing.
- **Output:** `finite_weights` (count of usable weights), `bins`, `scale`, `converged`, `variance`, `parameters` and `provenance`.
- **Example code:**

```python
import tempfile
from pathlib import Path
from dnakit.epigenome import build_contact_matrix, balance_contact_matrix

work = Path(tempfile.mkdtemp())
cool = build_contact_matrix("tiny.pairs", work / "tiny.cool", chromsizes={"chr1": 1000, "chr2": 600}, binsize=100)["cool_path"]
result = balance_contact_matrix(cool, min_nnz=1)
print(result["finite_weights"], result["converged"], result["scale"])
```

- **Example output:**

```text
6 True 1.0
```

## 4) `dnakit.epigenome.calculate_coverage()`

- **Purpose:** Compute, for each bin, the sum of cis contacts and the sum of all contacts.
- **Method:** Call cooltools' `coverage`. A diagonal pixel of the symmetric matrix counts once on each side. `cis` includes only partners on the same chromosome, and `total` includes all partners. With `balanced=True`, the `weight` column is used.
- **API:** `dnakit.epigenome.calculate_coverage(cool_path, *, balanced=False, ignore_diags=0, max_bins=100000)`.
- **Input:** A `.cool` file path; `balanced`; `ignore_diags`, the number of diagonals to skip; `max_bins`, which limits the number of bins.
- **Output:** A `bins` list where each item has `chrom`, `start`, `end`, `cis` and `total`, plus `balanced`, `ignore_diags` and `provenance`.
- **Example code:**

```python
import tempfile
from pathlib import Path
from dnakit.epigenome import build_contact_matrix, calculate_coverage

work = Path(tempfile.mkdtemp())
cool = build_contact_matrix("tiny.pairs", work / "tiny.cool", chromsizes={"chr1": 1000, "chr2": 600}, binsize=100)["cool_path"]
rows = calculate_coverage(cool, ignore_diags=0)["bins"]
print([(row["start"], row["cis"], row["total"]) for row in rows[:3]])
```

- **Example output:**

```text
[(0, 3.0, 3.0), (100, 2.0, 3.0), (200, 1.0, 1.0)]
```

## 5) `dnakit.epigenome.calculate_expected_contacts()`

- **Purpose:** Compute the mean contact frequency at each genomic distance within a chromosome, P(s).
- **Method:** Call cooltools' `expected_cis` with smoothing turned off. For each diagonal, `sum` is the total contact count and `mean` is `sum / valid_pixels`, where `valid_pixels` is the number of valid pixels on that diagonal. With `balanced=True`, balanced values are used.
- **API:** `dnakit.epigenome.calculate_expected_contacts(cool_path, *, balanced=False, ignore_diags=2, max_bins=100000)`.
- **Input:** A `.cool` file path; `balanced`; `ignore_diags`; `max_bins`.
- **Output:** A `rows` list where each item has `chrom`, `distance_bins`, `distance_bp`, `valid_pixels`, `sum` and `mean`, plus `provenance`.
- **Example code:**

```python
import tempfile
from pathlib import Path
from dnakit.epigenome import build_contact_matrix, calculate_expected_contacts

work = Path(tempfile.mkdtemp())
cool = build_contact_matrix("tiny.pairs", work / "tiny.cool", chromsizes={"chr1": 1000, "chr2": 600}, binsize=100)["cool_path"]
rows = calculate_expected_contacts(cool, ignore_diags=0)["rows"]
print(rows[1])
```

- **Example output:**

```text
{'chrom': 'chr1', 'distance_bins': 1, 'distance_bp': 100, 'valid_pixels': 9, 'sum': 2.0, 'mean': 0.2222222222222222}
```

## 6) `dnakit.epigenome.calculate_insulation_score()`

- **Purpose:** Compute the diamond insulation score and mark insulating boundaries.
- **Method:** Call cooltools' `insulation`. The input must be a balanced `.cool` file with a `weight` column. Boundaries are thresholded with the Li method (`threshold_li` from scikit-image). With `ignore_diags=None`, the number of diagonals recorded during balancing is used.
- **API:** `dnakit.epigenome.calculate_insulation_score(cool_path, window_bp, *, ignore_diags=None, max_bins=100000)`.
- **Input:** A balanced `.cool` file path; `window_bp`, the diamond window size in bp; `max_bins`, which limits the number of bins.
- **Output:** `window_bp`, `boundary_count` and a `bins` list. Each bin has `chrom`, `start`, `end`, `is_bad_bin`, `log2_insulation_score`, `valid_pixels`, `boundary_strength` and `is_boundary`.
- **Example code:**

```python
import tempfile
from pathlib import Path
from dnakit.epigenome import build_contact_matrix, balance_contact_matrix, calculate_insulation_score

work = Path(tempfile.mkdtemp())
lines = ["## pairs format v1.0", "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2"]
count = 0
for a in range(12):
    for b in range(a, 12):
        if b - a <= 2 and (a < 6) == (b < 6):
            for _ in range(3):
                count += 1
                lines.append(f"r{count}\tchr1\t{a * 100 + 10}\tchr1\t{b * 100 + 10}\t+\t-")
(work / "dense.pairs").write_text("\n".join(lines) + "\n", encoding="utf-8")
cool = build_contact_matrix(work / "dense.pairs", work / "dense.cool", chromsizes={"chr1": 1200}, binsize=100)["cool_path"]
balance_contact_matrix(cool, min_nnz=1)
result = calculate_insulation_score(cool, 300)
print(result["boundary_count"], [round(row["log2_insulation_score"], 3) for row in result["bins"][4:8]])
```

- **Example output:**

```text
0 [-0.585, -1.585, -1.585, -0.585]
```

## 7) `dnakit.epigenome.read_pairs_header()`

- **Purpose:** Read the column names, chromosome sizes and sort order from a `.pairs` header.
- **Method:** pairtools' `get_header` separates the header lines that start with `#`. `extract_column_names` and `extract_chromsizes` parse them. Without `#chromsize:` lines, `chromsizes` is an empty dictionary.
- **API:** `dnakit.epigenome.read_pairs_header(pairs_path)`.
- **Input:** A `.pairs` file path.
- **Output:** `header_lines`, `columns`, `sorted` (the `#sorted:` value, or `null`), `chromsizes` and `provenance`.
- **Example code:**

```python
from dnakit.epigenome import read_pairs_header

result = read_pairs_header("dup.pairs")
print(result["columns"], result["chromsizes"], result["sorted"])
```

- **Example output:**

```text
['readID', 'chrom1', 'pos1', 'chrom2', 'pos2', 'strand1', 'strand2', 'pair_type'] {'chr1': 1000} chr1-pos1-chr2-pos2
```

## 8) `dnakit.epigenome.summarize_pairs()`

- **Purpose:** Count mapped pairs, single-sided pairs, cis and trans contacts, and pairs of each `pair_type`.
- **Method:** Use pairtools' `PairCounter`. The file is read in blocks of `chunksize` rows and each block goes to `add_pairs_from_dataframe`. The file must have a `pair_type` column, and the call fails without one.
- **API:** `dnakit.epigenome.summarize_pairs(pairs_path, *, chunksize=100000)`.
- **Input:** A `.pairs` file path with a `pair_type` column; `chunksize`, the number of rows per block.
- **Output:** `total`, `total_mapped`, `total_single_sided_mapped`, `total_unmapped`, `cis`, `trans`, `pair_types` and `provenance`.
- **Example code:**

```python
from dnakit.epigenome import summarize_pairs

result = summarize_pairs("typed.pairs")
print(result["total"], result["cis"], result["trans"], result["pair_types"])
```

- **Example output:**

```text
5 3 1 {'UU': 3, 'UR': 1, 'NU': 1}
```

## 9) `dnakit.epigenome.deduplicate_pairs()`

- **Purpose:** Identify PCR duplicate pairs and count them. Optionally write the deduplicated pairs file.
- **Method:** Call pairtools' `streaming_dedup` with the scipy backend. Pairs whose two ends both lie within `max_mismatch` bp form a cluster, and one pair per cluster is kept. Pairs on opposite strands are not duplicates. Single-sided pairs are counted only; they are written only when `unmapped_path` is given. `duplicate_fraction` is `total_dups / total_mapped`.
- **API:** `dnakit.epigenome.deduplicate_pairs(pairs_path, *, output_path=None, unmapped_path=None, max_mismatch=3, method="max", chunksize=10000, carryover=100)`.
- **Input:** A `.pairs` file path with strand columns; optional `output_path` and `unmapped_path`, which must be new files; `max_mismatch`; `method`, either `max` or `sum`.
- **Output:** `total`, `total_mapped`, `total_single_sided_mapped`, `total_nodups`, `total_dups`, `duplicate_fraction`, `output`, `unmapped_output` and `provenance`.
- **Example code:**

```python
from dnakit.epigenome import deduplicate_pairs

result = deduplicate_pairs("dup.pairs", max_mismatch=3)
print(result["total"], result["total_dups"], result["total_nodups"], result["duplicate_fraction"])
```

- **Example output:**

```text
7 2 4 0.3333333333333333
```

## 10) `dnakit.epigenome.build_fragment_pileup()`

- **Purpose:** Accumulate fragments on one chromosome into piecewise-constant coverage.
- **Method:** Call MACS3's `Signal.Pileup.quick_pileup`. The contiguous segments start at position 0 and end at the largest fragment end. `baseline_value` is the minimum value of each segment, and `scale_factor` multiplies the coverage. Starts and ends are sorted separately, and every fragment end must be greater than its start.
- **API:** `dnakit.epigenome.build_fragment_pileup(starts, ends, *, scale_factor=1.0, baseline_value=0.0, max_fragments=1000000)`.
- **Input:** Equal-length lists of fragment starts and ends (0-based half-open); a scale factor; a baseline value; the maximum number of fragments.
- **Output:** `fragments`, `segments` (each with `start`, `end` and `value`) and `provenance`.
- **Example code:**

```python
from dnakit.epigenome import build_fragment_pileup

result = build_fragment_pileup([10, 50, 100], [40, 90, 140])
print([(row["start"], row["end"], row["value"]) for row in result["segments"]])
```

- **Example output:**

```text
[(0, 10, 0.0), (10, 40, 1.0), (40, 50, 0.0), (50, 90, 1.0), (90, 100, 0.0), (100, 140, 1.0)]
```

## 11) `dnakit.epigenome.build_extended_tag_pileup()`

- **Purpose:** Compute coverage from single-end tag positions extended by a fixed length on both sides.
- **Method:** Call MACS3's `Signal.Pileup.naive_quick_pileup`. Each tag covers `[position - extension, position + extension)`, and any part below 0 is clipped to 0. Positions are sorted before the calculation.
- **API:** `dnakit.epigenome.build_extended_tag_pileup(positions, extension, *, max_tags=1000000)`.
- **Input:** A list of tag positions (0-based); the extension length; the maximum number of tags.
- **Output:** `tags`, `extension`, `segments` (piecewise coverage) and `provenance`.
- **Example code:**

```python
from dnakit.epigenome import build_extended_tag_pileup

result = build_extended_tag_pileup([100, 200, 205], 50)
print([(row["start"], row["end"], row["value"]) for row in result["segments"]])
```

- **Example output:**

```text
[(0, 50, 0.0), (50, 155, 1.0), (155, 250, 2.0), (250, 255, 1.0)]
```

## 12) `dnakit.epigenome.call_fragment_peaks()`

- **Purpose:** Find summits above a threshold in fragment coverage.
- **Method:** Compute coverage with `quick_pileup`, then call MACS3's `naive_call_peaks`. Adjacent segments with coverage at least `min_value` are merged when the gap between them is at most `max_gap`. A merged region is kept when it is at least `min_length` long, and the summit is the position of the highest value. This is a naive rule. It is not MACS3's statistical peak model and reports no p-values.
- **API:** `dnakit.epigenome.call_fragment_peaks(starts, ends, *, min_value, max_gap=50, min_length=200, max_fragments=1000000)`.
- **Input:** Lists of fragment starts and ends; `min_value`; `max_gap`; `min_length`; the maximum number of fragments.
- **Output:** `fragments`, `min_value`, `max_gap`, `min_length`, `peaks` (each with `summit` and `height`) and `provenance`.
- **Example code:**

```python
from dnakit.epigenome import call_fragment_peaks

starts = [100] * 5 + [400] * 2 + [1000] * 4
ends = [300] * 5 + [500] * 2 + [1200] * 4
print(call_fragment_peaks(starts, ends, min_value=1.0)["peaks"])
```

- **Example output:**

```text
[{'summit': 200, 'height': 5.0}, {'summit': 1100, 'height': 4.0}]
```

## 13) `dnakit.epigenome.filter_allc_sites()`

- **Purpose:** Filter allc cytosine records by chromosome, cytosine context and coverage.
- **Method:** Call methylpy's `utilities.filter_allc_files`. A context such as `CGN` or `CHH` expands to the specific trinucleotide contexts. Coverage is the total read count in column 6 of the allc file.
- **API:** `dnakit.epigenome.filter_allc_sites(allc_path, output_path, *, chromosomes=None, mc_type=None, min_coverage=0, max_coverage=None)`.
- **Input:** An allc file path (plain text); a new output path; an optional chromosome list; `mc_type`; a coverage range.
- **Output:** `input_sites`, `kept_sites`, `output`, `chromosomes`, `mc_type`, `min_coverage`, `max_coverage` and `provenance`.
- **Example code:**

```python
import tempfile
from pathlib import Path
from dnakit.epigenome import filter_allc_sites

work = Path(tempfile.mkdtemp())
result = filter_allc_sites("tiny.allc", work / "filtered.allc", chromosomes=["chrL", "chr1"], mc_type="CGN", min_coverage=2)
print(result["input_sites"], result["kept_sites"])
```

- **Example output:**

```text
6 4
```

## 14) `dnakit.epigenome.calculate_non_conversion_rate()`

- **Purpose:** Estimate the bisulfite non-conversion rate from an unmethylated control region, for example lambda DNA.
- **Method:** Call methylpy's `calculate_non_conversion_rate`: the sum of methylated reads in the control region divided by the sum of all cytosine reads, across all contexts. methylpy writes its index beside the file it reads, so the function links the file into a temporary directory first, and no `.idx` file is created next to the input. Chromosome names are matched exactly as they appear in the allc file.
- **API:** `dnakit.epigenome.calculate_non_conversion_rate(allc_path, control_region)`.
- **Input:** An allc file path; a control region such as `chrom` or `chrom:start-end` (1-based, inclusive).
- **Output:** `control_region`, `non_conversion_rate` and `provenance`.
- **Example code:**

```python
from dnakit.epigenome import calculate_non_conversion_rate

print(calculate_non_conversion_rate("tiny.allc", "chrL")["non_conversion_rate"])
```

- **Example output:**

```text
0.041666666666666664
```

## 15) `dnakit.epigenome.calculate_site_binomial_pvalues()`

- **Purpose:** Run a one-sided binomial test at each cytosine to check whether methylation exceeds the non-conversion rate.
- **Method:** Call methylpy's `allc_run_binom_tests`. The p-value is `P(X >= methylated reads)`, where `X ~ Binomial(coverage, non_conversion)`. Sites with no methylated reads, or with coverage below `min_coverage`, get p = 2.0 and are not tested. The output file is sorted by p-value with the unix `sort` command.
- **API:** `dnakit.epigenome.calculate_site_binomial_pvalues(allc_path, output_path, *, non_conversion, min_coverage=1, sort_buffer="64M", max_report=20)`.
- **Input:** An allc file path; a new output path (letters, digits and `._/@+=,:-` only); `non_conversion` (0 to 1); `min_coverage`; `sort_buffer` such as `64M`; `max_report`, the number of lowest p-values to return.
- **Output:** `sites_written`, `tested_sites`, `mc_class_counts`, `lowest_p_values` (each with `chrom`, `pos`, `strand`, `mc_class`, `methylated`, `coverage` and `p_value`), `non_conversion_rate` and `provenance`.
- **Example code:**

```python
import tempfile
from pathlib import Path
from dnakit.epigenome import calculate_site_binomial_pvalues

work = Path(tempfile.mkdtemp())
result = calculate_site_binomial_pvalues("tiny.allc", work / "pvalues.tsv", non_conversion=0.01, min_coverage=2, max_report=2)
print(result["tested_sites"], result["lowest_p_values"][0])
```

- **Example output:**

```text
4 {'chrom': 'chr1', 'pos': 5, 'strand': '+', 'mc_class': 'CGN', 'methylated': 2, 'coverage': 2, 'p_value': 0.0001}
```

## 16) `dnakit.epigenome.summarize_bigwig_bins()`

- **Purpose:** Compute the mean signal of one or more bigWig files in fixed-width genome bins.
- **Method:** Call deeptools' `getScorePerBigWigBin.getScorePerBin`. Only chromosomes with the same name and length in every file are used. A bin's value is the mean over the bases in it that have data, and a bin without data is `null`. The last bin of a chromosome may be shorter.
- **API:** `dnakit.epigenome.summarize_bigwig_bins(bigwig_files, *, bin_size=1000, max_bins=100000, number_of_processors=1)`.
- **Input:** One or more bigWig paths; `bin_size` in bp; `max_bins`, which limits the total number of bins; `number_of_processors`.
- **Output:** `bin_size`, `files`, `chromosomes`, `bins` (each with `chrom`, `start`, `end` and one `values` entry per file) and `provenance`.
- **Example code:**

```python
import tempfile
from pathlib import Path
import pyBigWig
from dnakit.epigenome import summarize_bigwig_bins

work = Path(tempfile.mkdtemp())
bw = pyBigWig.open(str(work / "a.bw"), "w")
bw.addHeader([("chr1", 1000)])
bw.addEntries(["chr1", "chr1"], [0, 500], ends=[500, 1000], values=[2.0, 4.0])
bw.close()
result = summarize_bigwig_bins([work / "a.bw"], bin_size=200)
print([row["values"] for row in result["bins"]])
```

- **Example output:**

```text
[[2.0], [2.0], [3.0], [4.0], [4.0]]
```

## 17) `dnakit.epigenome.count_bam_reads_per_bin()`

- **Purpose:** Count the alignments that overlap each fixed-width genome bin.
- **Method:** Call deeptools' `countReadsPerBin.CountReadsPerBin` without read extension. Each alignment is counted on its own span, so one alignment can count in several bins it crosses. Alignments with mapping quality below `min_mapping_quality` are dropped. The input BAM files must have a `.bai` index.
- **API:** `dnakit.epigenome.count_bam_reads_per_bin(bam_files, *, bin_size=1000, max_bins=100000, min_mapping_quality=0, number_of_processors=1)`.
- **Input:** One or more indexed BAM paths; `bin_size`; `max_bins`; `min_mapping_quality` (0 to 255); `number_of_processors`.
- **Output:** `bin_size`, `files`, `chromosomes`, `bins` (each with `chrom`, `start`, `end` and one `counts` entry per file) and `provenance`.
- **Example code:**

```python
import tempfile
from pathlib import Path
import pysam
from dnakit.epigenome import count_bam_reads_per_bin

work = Path(tempfile.mkdtemp())
header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
bam = work / "reads.bam"
with pysam.AlignmentFile(str(bam), "wb", header=header) as out:
    for index, start in enumerate([100, 110, 120, 300, 310, 320, 400, 410]):
        read = pysam.AlignedSegment()
        read.query_name = f"q{index}"
        read.query_sequence = "A" * 50
        read.query_qualities = pysam.qualitystring_to_array("I" * 50)
        read.flag = 0
        read.reference_id = 0
        read.reference_start = start
        read.mapping_quality = 60
        read.cigar = [(0, 50)]
        out.write(read)
pysam.index(str(bam))
result = count_bam_reads_per_bin([bam], bin_size=100)
print([row["counts"][0] for row in result["bins"][:5]])
```

- **Example output:**

```text
[0, 3, 0, 3, 2]
```

## 18) `dnakit.epigenome.estimate_fragment_lengths()`

- **Purpose:** Estimate the fragment length and read length distributions by sampling an indexed BAM file.
- **Method:** Call deeptools' `get_read_and_fragment_length`. The fragment length is the template length of read1 in proper pairs. If the file has no proper pairs, all mapped alignments are used. The estimate comes from sampled genome regions, not from every alignment.
- **API:** `dnakit.epigenome.estimate_fragment_lengths(bam_file, *, number_of_processors=1)`.
- **Input:** An indexed BAM path; `number_of_processors`.
- **Output:** `bam`, `fragment_length` and `read_length` (each a dictionary with `sample_size`, `min`, `median`, `mean`, `std` and several percentiles, or `null` when the file has no data), plus `provenance`.
- **Example code:**

```python
import tempfile
from pathlib import Path
import pysam
from dnakit.epigenome import estimate_fragment_lengths

work = Path(tempfile.mkdtemp())
def make_read(name, flag, start, tlen=0):
    read = pysam.AlignedSegment()
    read.query_name = name
    read.query_sequence = "A" * 50
    read.query_qualities = pysam.qualitystring_to_array("I" * 50)
    read.flag = flag
    read.reference_id = 0
    read.reference_start = start
    read.mapping_quality = 60
    read.cigar = [(0, 50)]
    read.next_reference_id = 0
    read.template_length = tlen
    return read

header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
reads = []
for index, start in enumerate([100, 110, 120, 400, 410]):
    reads.append(make_read(f"q{index}", 99, start, 250))
    reads.append(make_read(f"q{index}", 147, start + 200, -250))
reads.sort(key=lambda read: read.reference_start)
bam = work / "pairs.bam"
with pysam.AlignmentFile(str(bam), "wb", header=header) as out:
    for read in reads:
        out.write(read)
pysam.index(str(bam))
result = estimate_fragment_lengths(bam)
print(result["fragment_length"]["sample_size"], result["fragment_length"]["median"], result["read_length"]["median"])
```

- **Example output:**

```text
5 250.0 50.0
```

## Installing the optional dependencies

```bash
python -m pip install "dnakit[bio-cooler]"          # cooler: build_contact_matrix / fetch_contact_matrix / balance_contact_matrix
python -m pip install "dnakit[bio-cooltools]"       # cooltools: calculate_* on .cool files
python -m pip install "dnakit[bio-pairtools]"       # pairtools: read_pairs_header / summarize_pairs / deduplicate_pairs
python -m pip install "dnakit[bio-macs3]"           # MACS3: build_* pileups and call_fragment_peaks
python -m pip install "dnakit[bio-methylpy]"        # methylpy: allc filtering and binomial tests
python -m pip install "dnakit[bio-deeptools]"       # deeptools; also needs pyBigWig and pysam (bio-pybigwig, bio-reads)
```

Project pages: [cooler](https://github.com/open2c/cooler), [cooltools](https://github.com/open2c/cooltools), [pairtools](https://github.com/mirnylab/pairtools), [MACS3](https://github.com/macs3-project/MACS), [methylpy](https://github.com/yupenghe/methylpy), [deeptools](https://github.com/deeptools/deepTools), [pyBigWig](https://github.com/deeptools/pyBigWig), [pysam](https://github.com/pysam-developers/pysam).

methylpy 1.4.7 imports `pkg_resources` at load time, which setuptools 80 and earlier provide. If that module is missing from the environment, install `setuptools<81`.

## How these functions are verified

The test file is `others/tests/unit/test_epigenome_backends.py`, with 34 tests. They check the following:

- Contact matrices: `build_contact_matrix()` counts match a hand count, including pixels (0,1)=2, (0,2)=1, (1,10)=1, (5,5)=1, (10,10)=1 and (11,14)=1. `fetch_contact_matrix()` equals a direct `Cooler.matrix().fetch()`, and the weights written by `balance_contact_matrix()` equal a direct `balance_cooler` call.
- cooltools: `calculate_coverage()` cis and total sums match a hand sum and a direct `cooltools.coverage` call. Diagonal means from `calculate_expected_contacts()` match a hand calculation and a direct `expected_cis` call. `calculate_insulation_score()` matches a direct `insulation` call and finds the synthetic boundary near 2 kb.
- pairtools: `summarize_pairs()` counts, the duplicate chain a-b-c, the non-duplicate pairs on opposite strands, and the output for single-sided pairs are checked by hand.
- MACS3: segment coverage, extended-tag coverage and summits match hand results.
- methylpy: the non-conversion rate matches the hand ratio 1/24. Binomial p-values match the hand formulas (0.01^2 and 1 - 0.99^2) and scipy. The allc filter output is checked by hand.
- deeptools: bigWig bin means, the handling of missing bases, BAM bin counts (including the MAPQ filter) and fragment lengths match hand results.
- Each module names its `bio-*` install hint when the package is missing. Missing inputs, existing outputs, bin counts above `max_bins` and unsafe shell output paths are all rejected.
- The tests do not use the network. There is no `DNAKIT_LIVE_TESTS` test because this batch has no network service.

## Current scope

This batch wraps three operations from each of the six packages, for 18 functions in total. Not wrapped:

- cooler: `cooler cload` from BAM and `.hic` input, `.mcool` multi-resolution files, `zoomify`, `merge` and `coarsen`.
- cooltools: compartment eigenvectors, dot calling, saddle, pileup/snipping and directionality.
- pairtools: `parse` (SAM to pairs), `sort`, `select`, `split`, `restrict` and `flip`, plus the YAML statistics output and filter expressions. Only uncompressed plain-text pairs were tested as input.
- MACS3: the `callpeak` statistical model. Only the naive peak rule is wrapped. BAM-based `callpeak` is already available through the CLI in `dnakit.ngs.call_chip_peaks()`, so it is not wrapped again.
- methylpy: the alignment pipeline, BAM-based `call_methylated_sites`, DMR calling (DMRfind) and regional methylation levels (`get_methylation_levels_DMRfind`).
- deeptools: `bamCoverage` writing bigWig files, `computeMatrix` with heatmaps, correlation and PCA, and GC bias correction.
- Relation to existing adapters: `dnakit.integrations.summarize_pybigwig_region()` gives single-region pyBigWig statistics, which differ from the deeptools window matrix here. `dnakit.epigenome.extract_methylation_methyldackel()` calls the MethylDackel CLI and is unrelated to the methylpy functions.
