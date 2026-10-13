# G10: file and matrix functions from cooler, cooltools and pairtools

This group wraps 42 public functions from three upstream packages that work on DNA or genomic coordinates and run offline on local inputs: `cooler_functions.py` has 11 functions (chromosome sizes, genomic intervals, FASTA lengths, and contact-matrix file operations: reading, coarsening, merging and multi-resolution output); `cooltools_functions.py` has 20 functions (expected contacts, compartments, saddle plots, virtual 4C, pileup, viewframes, and dense-matrix correction and decomposition); `pairtools_functions.py` has 11 functions (restriction fragments, scaling curves, genomic regions, read mismatches, SAM alignment parsing, .pairs headers and haplotype phase).

The functions live in the three modules of `dnakit.epigenome` (for example `dnakit.epigenome`). Each one imports its upstream package only when called; a missing optional dependency raises `BackendUnavailableError`. Inputs are checked first, and upstream `ValueError`, `KeyError` and `OSError` are reported as `ConfigurationError`. Results are plain dictionaries, NaN is returned as `None`, and grouped results are sorted explicitly. `fetch_chrom_sizes` needs network access; the default tests use a mocked client, and its example only shows the argument check.

The examples build contact matrices with the existing `build_contact_matrix()` from `contacts.pairs` below and balance them with `balance_contact_matrix()`; each example writes its outputs to `g10_out/` in the current directory.

## Example files

The files below are written to the example directory before the examples run.

**Example file:** `chrom.sizes`

```text
chr1	1000
chr2	800
```

**Example file:** `contacts.pairs`

```text
## pairs format v1.0
#chromsize: chr1 1000
#chromsize: chr2 800
#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type
r0	chr2	649	chr1	86	+	-	UU
r1	chr1	182	chr1	447	+	-	UU
r2	chr1	801	chr1	843	+	-	UU
r3	chr2	266	chr2	319	+	-	UU
r4	chr1	433	chr1	446	+	-	UU
r5	chr1	691	chr1	699	+	-	UU
r6	chr2	587	chr2	642	+	-	UU
r7	chr1	887	chr1	919	+	-	UU
r8	chr1	517	chr1	672	+	-	UU
r9	chr2	139	chr2	384	+	-	UU
r10	chr1	738	chr1	847	+	-	UU
r11	chr1	320	chr2	519	+	-	UU
r12	chr2	695	chr2	696	+	-	UU
r13	chr1	293	chr1	433	+	-	UU
r14	chr1	943	chr1	999	+	-	UU
r15	chr2	239	chr2	318	+	-	UU
r16	chr1	663	chr1	734	+	-	UU
r17	chr1	585	chr1	656	+	-	UU
r18	chr2	380	chr2	639	+	-	UU
r19	chr1	31	chr1	201	+	-	UU
r20	chr1	254	chr1	429	+	-	UU
r21	chr2	73	chr2	441	+	-	UU
r22	chr1	928	chr2	166	+	-	UU
r23	chr1	248	chr1	396	+	-	UU
r24	chr2	239	chr2	495	+	-	UU
r25	chr1	654	chr1	817	+	-	UU
r26	chr1	219	chr1	358	+	-	UU
r27	chr2	4	chr2	217	+	-	UU
r28	chr1	683	chr1	730	+	-	UU
r29	chr1	957	chr1	999	+	-	UU
r30	chr2	607	chr2	631	+	-	UU
r31	chr1	593	chr1	679	+	-	UU
r32	chr1	849	chr1	912	+	-	UU
r33	chr2	309	chr1	147	+	-	UU
r34	chr1	799	chr1	975	+	-	UU
r35	chr1	292	chr1	366	+	-	UU
r36	chr2	419	chr2	558	+	-	UU
r37	chr1	562	chr1	999	+	-	UU
r38	chr1	931	chr1	954	+	-	UU
r39	chr2	158	chr2	375	+	-	UU
```

**Example file:** `expected.tsv`

```text
region1	region2	dist	n_valid	count.avg	balanced.avg
chr1	chr1	0	10	1.2	0.9
chr1	chr1	1	9	0.8	0.7
chr1	chr1	2	8	0.5	0.6
chr2	chr2	0	6	1.0	0.8
```

**Example file:** `features.bed`

```text
chr1	300	400
chr1	500	600
chr2	200	300
```

**Example file:** `genome.fa`

```text
>chr1
ACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCA
>chr2
GGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATT
```

**Example file:** `rearrange_view.tsv`

```text
chrom	start	end	name	new_chrom	strand
chr1	0	1000	chr1	chrX	+
chr2	0	800	chr2	chrX	-
```

**Example file:** `regions.bed`

```text
chr1	0	150	x
chr1	200	300	y
chr2	10	20	z
```

**Example file:** `sample.sam`

```text
@HD	VN:1.6	SO:queryname
@SQ	SN:chr1	LN:1000
@SQ	SN:chr2	LN:800
r1	65	chr1	101	60	10M	=	101	0	ACGTACGTAC	IIIIIIIIII	NM:i:0	AS:i:10	XS:i:0
r1	129	chr2	301	60	10M	=	301	0	ACGTACGTAC	IIIIIIIIII	NM:i:0	AS:i:10	XS:i:0
```

**Example file:** `track.bed`

```text
chr1	0	100	-1.000000
chr1	100	200	-0.777778
chr1	200	300	-0.555556
chr1	300	400	-0.333333
chr1	400	500	-0.111111
chr1	500	600	0.111111
chr1	600	700	0.333333
chr1	700	800	0.555556
chr1	800	900	0.777778
chr1	900	1000	1.000000
chr2	0	100	1.000000
chr2	100	200	0.714286
chr2	200	300	0.428571
chr2	300	400	0.142857
chr2	400	500	-0.142857
chr2	500	600	-0.428571
chr2	600	700	-0.714286
chr2	700	800	-1.000000
```

**Example file:** `view.bed`

```text
chr1	0	1000	chr1
chr2	0	800	chr2
```


## cooler functions (11)

## 1) `dnakit.epigenome.read_chrom_sizes()`

- **Purpose:** Read chromosome lengths from a UCSC-style chrom.sizes or chromInfo file.
- **Method:** Calls cooler.util.read_chromsizes. By default only chr1-22, chrX, chrY and chrM names are kept; `all_names=True` keeps every name, and `name_patterns` takes custom regular expressions. Rows keep the file order.
- **API:** `dnakit.epigenome.read_chrom_sizes(path, *, all_names=False, name_patterns=None, verbose=True)`.
- **Input:** `path`: a two-column (name, length) file or a UCSC chromInfo file; `all_names`; `name_patterns` (a list of regular expressions).
- **Output:** A dictionary with `chromosome_count` and `chromosomes` (each row has `name` and `length`), plus `provenance`.
- **Example code:**

```python
from dnakit.epigenome import read_chrom_sizes

result = read_chrom_sizes("chrom.sizes", all_names=True, verbose=False)
for row in result["chromosomes"]:
    print(row["name"], row["length"])
```

- **Example output:**

```text
chr1 1000
chr2 800
```

## 2) `dnakit.epigenome.fetch_chrom_sizes()`

- **Purpose:** Download the chromosome-length table of a UCSC genome assembly (needs network).
- **Method:** First checks that `database` is a UCSC name of 2 to 64 letters, digits or underscores, then requests `chromInfo.txt.gz` through cooler.util.fetch_chromsizes. Network errors become `BackendUnavailableError`. Tests use a mocked HTTP client, so the default test run makes no network calls.
- **API:** `dnakit.epigenome.fetch_chrom_sizes(database, *, all_names=False, name_patterns=None, verbose=True)`.
- **Input:** `database`, for example `hg38`; `all_names` and `name_patterns` as in `read_chrom_sizes`.
- **Output:** The same dictionary as `read_chrom_sizes`, with `source` giving the download address.
- **Example code:**

```python
from dnakit.epigenome import fetch_chrom_sizes
from dnakit.exceptions import ConfigurationError

try:
    fetch_chrom_sizes("hg 38", verbose=False)
except ConfigurationError as exc:
    print("rejected:", str(exc).split("] ", 1)[-1])
```

- **Example output:**

```text
rejected: database must be a UCSC assembly name such as hg38.
```

## 3) `dnakit.epigenome.parse_genomic_region()`

- **Purpose:** Parse a region string such as chr1:100-200 into a half-open interval (0-based start, 1-based end).
- **Method:** Calls cooler.util.parse_region. A chromosome name alone gives the whole chromosome when `chromsizes` is given. Ranges outside the chromosome are rejected.
- **API:** `dnakit.epigenome.parse_genomic_region(region, *, chromsizes=None, verbose=True)`.
- **Input:** `region`: a string; `chromsizes`: an optional mapping of chromosome lengths.
- **Output:** A dictionary with `chrom`, `start`, `end` and `length`.
- **Example code:**

```python
from dnakit.epigenome import parse_genomic_region

region = parse_genomic_region("chr1:100-200", verbose=False)
print(region["chrom"], region["start"], region["end"], region["length"])
whole = parse_genomic_region("chr2", chromsizes={"chr1": 1000, "chr2": 800}, verbose=False)
print(whole["start"], whole["end"])
```

- **Example output:**

```text
chr1 100 200 100
0 800
```

## 4) `dnakit.epigenome.query_bed_intervals()`

- **Purpose:** Return the intervals of one region from a sorted, non-overlapping BED file.
- **Method:** Checks start and end coordinates, sorting and overlap within each chromosome, then calls cooler.util.bedslice. Intervals that overlap the region are returned, sorted by start and end.
- **API:** `dnakit.epigenome.query_bed_intervals(bed_path, region, *, chromsizes, verbose=True)`.
- **Input:** `bed_path`: a BED file with at least chrom, start and end; `region`; `chromsizes`.
- **Output:** A dictionary with `interval_count` and `intervals` (each row has `chrom`, `start` and `end`).
- **Example code:**

```python
from dnakit.epigenome import query_bed_intervals

result = query_bed_intervals("regions.bed", "chr1:100-250", chromsizes={"chr1": 1000, "chr2": 800}, verbose=False)
for row in result["intervals"]:
    print(row["chrom"], row["start"], row["end"])
```

- **Example output:**

```text
chr1 0 150
chr1 200 300
```

## 5) `dnakit.epigenome.rename_cooler_chromosomes()`

- **Purpose:** Copy a cooler file and rename its chromosomes in the copy; the input file is unchanged.
- **Method:** Checks that the new names are unique and that the old names exist, copies the file, then calls cooler.rename_chroms on the copy.
- **API:** `dnakit.epigenome.rename_cooler_chromosomes(cool_path, output_path, rename_map, *, verbose=True)`.
- **Input:** `cool_path`: the input .cool file; `output_path`: a new file that must not exist; `rename_map`: old names mapped to new names.
- **Output:** A dictionary with `chromosomes_before`, `chromosomes_after` and `renamed`.
- **Example code:**

```python
import shutil
from pathlib import Path
from dnakit.epigenome.backends import build_contact_matrix

out = Path("g10_out")
shutil.rmtree(out, ignore_errors=True)
out.mkdir()
build_contact_matrix("contacts.pairs", out / "in.cool", chromsizes={"chr1": 1000, "chr2": 800}, binsize=100, verbose=False)
from dnakit.epigenome import rename_cooler_chromosomes

result = rename_cooler_chromosomes(out / "in.cool", out / "renamed.cool", {"chr1": "1"}, verbose=False)
print(result["chromosomes_before"], result["chromosomes_after"])
```

- **Example output:**

```text
['chr1', 'chr2'] ['1', 'chr2']
```

## 6) `dnakit.epigenome.coarsen_contact_matrix()`

- **Purpose:** Lower the resolution of a contact matrix by an integer factor (sum of k-by-k pixel blocks).
- **Method:** Calls cooler.coarsen_cooler one chromosome block at a time. The output bin size is `factor` times the input size, and the total number of contacts is unchanged.
- **API:** `dnakit.epigenome.coarsen_contact_matrix(cool_path, output_path, *, factor, chunksize=10000000, verbose=True)`.
- **Input:** `cool_path`; `output_path` (a new file); `factor` (an integer from 2 to 1000); `chunksize`.
- **Output:** A dictionary with `input_binsize`, `output_binsize`, `output_bins` and contact counts.
- **Example code:**

```python
import shutil
from pathlib import Path
from dnakit.epigenome.backends import build_contact_matrix

out = Path("g10_out")
shutil.rmtree(out, ignore_errors=True)
out.mkdir()
build_contact_matrix("contacts.pairs", out / "in.cool", chromsizes={"chr1": 1000, "chr2": 800}, binsize=100, verbose=False)
from dnakit.epigenome import coarsen_contact_matrix

result = coarsen_contact_matrix(out / "in.cool", out / "coarse.cool", factor=2, verbose=False)
print(result["input_binsize"], result["output_binsize"], result["output_total_contacts"])
```

- **Example output:**

```text
100 200 40
```

## 7) `dnakit.epigenome.merge_contact_matrices()`

- **Purpose:** Add two or more coolers that have identical bins into one file.
- **Method:** Checks that there are 2 to 64 distinct inputs, then calls cooler.merge_coolers. The output contact total is the sum of the inputs.
- **API:** `dnakit.epigenome.merge_contact_matrices(cool_paths, output_path, *, mergebuf=10000000, verbose=True)`.
- **Input:** `cool_paths`: a list of file paths; `output_path` (a new file); `mergebuf`.
- **Output:** A dictionary with `input_total_contacts`, `output_total_contacts` and `output_nonzero_pixels`.
- **Example code:**

```python
import shutil
from pathlib import Path
from dnakit.epigenome.backends import build_contact_matrix

out = Path("g10_out")
shutil.rmtree(out, ignore_errors=True)
out.mkdir()
build_contact_matrix("contacts.pairs", out / "in.cool", chromsizes={"chr1": 1000, "chr2": 800}, binsize=100, verbose=False)
from dnakit.epigenome.backends import build_contact_matrix
from dnakit.epigenome import merge_contact_matrices

build_contact_matrix("contacts.pairs", out / "second.cool", chromsizes={"chr1": 1000, "chr2": 800}, binsize=100, verbose=False)
result = merge_contact_matrices([out / "in.cool", out / "second.cool"], out / "merged.cool", verbose=False)
print(result["input_total_contacts"], result["output_total_contacts"])
```

- **Example output:**

```text
80 80
```

## 8) `dnakit.epigenome.zoomify_contact_matrix()`

- **Purpose:** Build a multi-resolution (.mcool) file from base coolers.
- **Method:** Calls cooler.zoomify_cooler, coarsening recursively; each resolution is stored under `/resolutions/<bin size>`. The result gives bin and contact counts for every level.
- **API:** `dnakit.epigenome.zoomify_contact_matrix(cool_paths, output_path, *, resolutions, chunksize=10000000, verbose=True)`.
- **Input:** `cool_paths`: one or more base files; `output_path` (a new file); `resolutions` (ascending, distinct bin sizes).
- **Output:** A dictionary whose `resolutions` list has `resolution`, `binsize`, `bins` and `total_contacts` for each level.
- **Example code:**

```python
import shutil
from pathlib import Path
from dnakit.epigenome.backends import build_contact_matrix

out = Path("g10_out")
shutil.rmtree(out, ignore_errors=True)
out.mkdir()
build_contact_matrix("contacts.pairs", out / "in.cool", chromsizes={"chr1": 1000, "chr2": 800}, binsize=100, verbose=False)
from dnakit.epigenome import zoomify_contact_matrix

result = zoomify_contact_matrix(out / "in.cool", out / "multi.mcool", resolutions=[100, 200, 400], verbose=False)
for level in result["resolutions"]:
    print(level["resolution"], level["bins"], level["total_contacts"])
```

- **Example output:**

```text
100 18 40
200 9 40
400 5 40
```

## 9) `dnakit.epigenome.load_fasta_lengths()`

- **Purpose:** Return the lengths of named FASTA records without loading the bases into memory.
- **Method:** Calls cooler.util.load_fasta, which opens a lazy pyfaidx index. Rows follow the order of `names`. A missing name is an error.
- **API:** `dnakit.epigenome.load_fasta_lengths(fasta_path, names, *, verbose=True)`.
- **Input:** `fasta_path`: a FASTA file; `names`: a list of sequence names.
- **Output:** A dictionary with `sequence_count` and `sequences` (each row has `name` and `length`).
- **Example code:**

```python
from dnakit.epigenome import load_fasta_lengths

result = load_fasta_lengths("genome.fa", ["chr2", "chr1"], verbose=False)
for row in result["sequences"]:
    print(row["name"], row["length"])
```

- **Example output:**

```text
chr2 800
chr1 1000
```

## 10) `dnakit.epigenome.read_cooler_pixels()`

- **Purpose:** Read the first nonzero pixels (upper-triangle contacts) of a cooler, optionally joined to bin coordinates.
- **Method:** Calls cooler.Cooler.pixels (api.pixels) and slices in stored order (bin1_id, bin2_id). With `join=True`, the chromosome, start and end of both bins are added.
- **API:** `dnakit.epigenome.read_cooler_pixels(cool_path, *, limit=100, join=False, verbose=True)`.
- **Input:** `cool_path`; `limit` (1 to 10000); `join` (a boolean).
- **Output:** A dictionary with `nonzero_pixels`, `returned`, `columns` and `pixels`.
- **Example code:**

```python
import shutil
from pathlib import Path
from dnakit.epigenome.backends import build_contact_matrix

out = Path("g10_out")
shutil.rmtree(out, ignore_errors=True)
out.mkdir()
build_contact_matrix("contacts.pairs", out / "in.cool", chromsizes={"chr1": 1000, "chr2": 800}, binsize=100, verbose=False)
from dnakit.epigenome import read_cooler_pixels

result = read_cooler_pixels(out / "in.cool", limit=3, join=True, verbose=False)
print(result["nonzero_pixels"], result["returned"])
for row in result["pixels"]:
    print(row["chrom1"], row["start1"], row["chrom2"], row["start2"], row["count"])
```

- **Example output:**

```text
28 3
chr1 0 chr1 200 1
chr1 0 chr2 600 1
chr1 100 chr1 400 1
```

## 11) `dnakit.epigenome.read_cooler_chroms()`

- **Purpose:** List the chromosomes stored in a cooler file, in matrix order.
- **Method:** Calls cooler.Cooler.chroms (api.chroms) to read the chromosome table.
- **API:** `dnakit.epigenome.read_cooler_chroms(cool_path, *, verbose=True)`.
- **Input:** `cool_path`: a .cool file.
- **Output:** A dictionary with `chromosome_count` and `chromosomes` (each row has `name` and `length`).
- **Example code:**

```python
import shutil
from pathlib import Path
from dnakit.epigenome.backends import build_contact_matrix

out = Path("g10_out")
shutil.rmtree(out, ignore_errors=True)
out.mkdir()
build_contact_matrix("contacts.pairs", out / "in.cool", chromsizes={"chr1": 1000, "chr2": 800}, binsize=100, verbose=False)
from dnakit.epigenome import read_cooler_chroms

result = read_cooler_chroms(out / "in.cool", verbose=False)
for row in result["chromosomes"]:
    print(row["name"], row["length"])
```

- **Example output:**

```text
chr1 1000
chr2 800
```

## cooltools functions (20)

## 12) `dnakit.epigenome.calculate_expected_trans()`

- **Purpose:** Average contact frequency for every pair of chromosomes (inter-chromosomal expected).
- **Method:** Calls cooltools.expected_trans with the balancing column (default `weight`). Rows are sorted by region1 and region2.
- **API:** `dnakit.epigenome.calculate_expected_trans(cool_path, *, weight='weight', verbose=True)`.
- **Input:** `cool_path`: a balanced .cool file; `weight`: the name of the weight column.
- **Output:** A dictionary with `pair_count` and `pairs` (each row has `region1`, `region2`, `count.sum`, `balanced.sum` and more).
- **Example code:**

```python
import shutil
from pathlib import Path
from dnakit.epigenome.backends import balance_contact_matrix, build_contact_matrix

out = Path("g10_out")
shutil.rmtree(out, ignore_errors=True)
out.mkdir()
build_contact_matrix("contacts.pairs", out / "bal.cool", chromsizes={"chr1": 1000, "chr2": 800}, binsize=100, verbose=False)
balance_contact_matrix(out / "bal.cool", min_nnz=1, verbose=False)
from dnakit.epigenome import calculate_expected_trans

result = calculate_expected_trans(out / "bal.cool", verbose=False)
print(result["pair_count"])
for row in result["pairs"]:
    print(row["region1"], row["region2"], row["count.sum"], round(row["balanced.sum"], 4))
```

- **Example output:**

```text
1
chr1 chr2 4.0 1.8912
```

## 13) `dnakit.epigenome.calculate_compartments()`

- **Purpose:** Cis A/B compartment eigenvalues and eigenvectors for each chromosome.
- **Method:** Calls cooltools.eigs_cis with whole chromosomes as the view. The ARPACK start vector and the C random seed are fixed so that repeated calls agree. The largest-magnitude entry of each eigenvector is positive.
- **API:** `dnakit.epigenome.calculate_compartments(cool_path, *, n_eigs=3, weight='weight', verbose=True)`.
- **Input:** `cool_path`; `n_eigs` (1 to 10); `weight`.
- **Output:** A dictionary with `eigenvalues` (one row per chromosome) and `eigenvectors` (one row per bin, with columns E1, E2 ...).
- **Example code:**

```python
import shutil
from pathlib import Path
from dnakit.epigenome.backends import balance_contact_matrix, build_contact_matrix

out = Path("g10_out")
shutil.rmtree(out, ignore_errors=True)
out.mkdir()
build_contact_matrix("contacts.pairs", out / "bal.cool", chromsizes={"chr1": 1000, "chr2": 800}, binsize=100, verbose=False)
balance_contact_matrix(out / "bal.cool", min_nnz=1, verbose=False)
from dnakit.epigenome import calculate_compartments

result = calculate_compartments(out / "bal.cool", n_eigs=2, verbose=False)
for row in result["eigenvalues"]:
    print(row["name"], round(row["eigval1"], 4), round(row["eigval2"], 4))
print(len(result["eigenvectors"]))
```

- **Example output:**

```text
chr1 8.0892 -6.6156
chr2 4.3574 -4.3272
18
```

## 14) `dnakit.epigenome.calculate_trans_compartments()`

- **Purpose:** Genome-wide trans compartment eigenvectors from inter-chromosomal contacts.
- **Method:** Calls cooltools.eigs_trans over all chromosomes, with the same fixed random inputs as above. Eigenvector signs follow the same convention.
- **API:** `dnakit.epigenome.calculate_trans_compartments(cool_path, *, n_eigs=3, weight='weight', verbose=True)`.
- **Input:** `cool_path`; `n_eigs`; `weight`.
- **Output:** A dictionary with `eigenvalues` and `eigenvectors`.
- **Example code:**

```python
import shutil
from pathlib import Path
from dnakit.epigenome.backends import balance_contact_matrix, build_contact_matrix

out = Path("g10_out")
shutil.rmtree(out, ignore_errors=True)
out.mkdir()
build_contact_matrix("contacts.pairs", out / "bal.cool", chromsizes={"chr1": 1000, "chr2": 800}, binsize=100, verbose=False)
balance_contact_matrix(out / "bal.cool", min_nnz=1, verbose=False)
from dnakit.epigenome import calculate_trans_compartments

result = calculate_trans_compartments(out / "bal.cool", n_eigs=1, verbose=False)
print(round(result["eigenvalues"][0]["eigval1"], 4), len(result["eigenvectors"]))
```

- **Example output:**

```text
8.0001 18
```

## 15) `dnakit.epigenome.calculate_saddle()`

- **Purpose:** Saddle plot: mean balanced contacts between bins of a genomic track, grouped by quantile.
- **Method:** Computes the expected values with cooltools.expected_cis or expected_trans, then calls cooltools.saddle with the quantile range `qrange`. The output matrix has two extra rows and columns, one ring for values below or above the quantile range.
- **API:** `dnakit.epigenome.calculate_saddle(cool_path, track_path, *, contact_type='cis', n_bins=5, qrange=(0.02, 0.98), min_diag=3, weight='weight', verbose=True)`.
- **Input:** `cool_path`; `track_path`: a track file with chrom, start, end and value columns; `n_bins`; `qrange` (low and high quantiles); `min_diag`.
- **Output:** A dictionary with `matrix_size`, `interaction_sum` and `interaction_count`.
- **Example code:**

```python
import shutil
from pathlib import Path
from dnakit.epigenome.backends import balance_contact_matrix, build_contact_matrix

out = Path("g10_out")
shutil.rmtree(out, ignore_errors=True)
out.mkdir()
build_contact_matrix("contacts.pairs", out / "bal.cool", chromsizes={"chr1": 1000, "chr2": 800}, binsize=100, verbose=False)
balance_contact_matrix(out / "bal.cool", min_nnz=1, verbose=False)
from dnakit.epigenome import calculate_saddle

result = calculate_saddle(out / "bal.cool", "track.bed", n_bins=3, qrange=(0.0, 1.0), min_diag=1, verbose=False)
print(result["matrix_size"])
print([round(value, 4) if value is not None else None for value in result["interaction_sum"][1]])
```

- **Example output:**

```text
5
[0.0, 6.902, 26.0064, 0.0, 0.0]
```

## 16) `dnakit.epigenome.sample_contacts()`

- **Purpose:** Draw a random subset of contacts and write it to a new cooler file.
- **Method:** With `count`, cooltools.sample performs exact sampling; with `fraction`, binomial sampling is used. Give exactly one. The seed is set only during the call, and numpy's global random state is restored afterwards.
- **API:** `dnakit.epigenome.sample_contacts(cool_path, output_path, *, count=None, fraction=None, seed=0, verbose=True)`.
- **Input:** `cool_path`; `output_path` (a new file); `count` or `fraction`; `seed`.
- **Output:** A dictionary with `method`, `seed`, `total_contacts_before` and `total_contacts_after`.
- **Example code:**

```python
import shutil
from pathlib import Path
from dnakit.epigenome.backends import balance_contact_matrix, build_contact_matrix

out = Path("g10_out")
shutil.rmtree(out, ignore_errors=True)
out.mkdir()
build_contact_matrix("contacts.pairs", out / "bal.cool", chromsizes={"chr1": 1000, "chr2": 800}, binsize=100, verbose=False)
balance_contact_matrix(out / "bal.cool", min_nnz=1, verbose=False)
from dnakit.epigenome import sample_contacts

result = sample_contacts(out / "bal.cool", out / "sampled.cool", count=20, seed=1, verbose=False)
print(result["method"], result["total_contacts_before"], result["total_contacts_after"])
```

- **Example output:**

```text
exact 40 20
```

## 17) `dnakit.epigenome.calculate_virtual_4c()`

- **Purpose:** Contact profile of one viewpoint against every bin of the genome (virtual 4C).
- **Method:** Parses a viewpoint such as chr1:0-100 and calls cooltools.virtual4c, which is marked experimental upstream. Missing values are returned as None.
- **API:** `dnakit.epigenome.calculate_virtual_4c(cool_path, viewpoint, *, weight='weight', verbose=True)`.
- **Input:** `cool_path`; `viewpoint`: a region string; `weight`.
- **Output:** A dictionary with `viewpoint`, `bin_count` and `profile` (one row per bin).
- **Example code:**

```python
import shutil
from pathlib import Path
from dnakit.epigenome.backends import balance_contact_matrix, build_contact_matrix

out = Path("g10_out")
shutil.rmtree(out, ignore_errors=True)
out.mkdir()
build_contact_matrix("contacts.pairs", out / "bal.cool", chromsizes={"chr1": 1000, "chr2": 800}, binsize=100, verbose=False)
balance_contact_matrix(out / "bal.cool", min_nnz=1, verbose=False)
from dnakit.epigenome import calculate_virtual_4c

result = calculate_virtual_4c(out / "bal.cool", "chr1:0-100", verbose=False)
print(result["bin_count"], result["viewpoint"]["chrom"], result["viewpoint"]["end"])
print([row["balanced"] is None or round(row["balanced"], 4) for row in result["profile"][:3]])
```

- **Example output:**

```text
18 chr1 100
[True, True, 0.5738]
```

## 18) `dnakit.epigenome.calculate_directionality_index()`

- **Purpose:** Directionality index (DI) for each bin, used to locate TAD boundaries.
- **Method:** Calls cooltools.directionality with a window of `window_bp` bases, which must be at least the bin size. Missing values are returned as None.
- **API:** `dnakit.epigenome.calculate_directionality_index(cool_path, *, window_bp=100000, weight='weight', verbose=True)`.
- **Input:** `cool_path`; `window_bp`; `weight`.
- **Output:** A dictionary with `bin_count` and `bins` (each row has `chrom`, `start`, `end` and the DI columns).
- **Example code:**

```python
import shutil
from pathlib import Path
from dnakit.epigenome.backends import balance_contact_matrix, build_contact_matrix

out = Path("g10_out")
shutil.rmtree(out, ignore_errors=True)
out.mkdir()
build_contact_matrix("contacts.pairs", out / "bal.cool", chromsizes={"chr1": 1000, "chr2": 800}, binsize=100, verbose=False)
balance_contact_matrix(out / "bal.cool", min_nnz=1, verbose=False)
from dnakit.epigenome import calculate_directionality_index

result = calculate_directionality_index(out / "bal.cool", window_bp=300, verbose=False)
print(result["bin_count"])
print([row["directionality_index_300"] is None or round(row["directionality_index_300"], 4) for row in result["bins"][:4]])
```

- **Example output:**

```text
18
[True, True, 0.0069, True]
```

## 19) `dnakit.epigenome.pileup_features()`

- **Purpose:** Average balanced contact snippet around each BED feature (pileup).
- **Method:** Calls cooltools.pileup to cut a snippet around each feature, then takes the mean of the finite values (NaN ignored). No expected-value normalisation is applied.
- **API:** `dnakit.epigenome.pileup_features(cool_path, features_path, *, flank_bp=100000, weight='weight', verbose=True)`.
- **Input:** `cool_path`; `features_path`: a BED file with at least chrom, start and end; `flank_bp`; `weight`.
- **Output:** A dictionary with `feature_count`, `features_with_data`, `snippet_size` and `mean_snippet`.
- **Example code:**

```python
import shutil
from pathlib import Path
from dnakit.epigenome.backends import balance_contact_matrix, build_contact_matrix

out = Path("g10_out")
shutil.rmtree(out, ignore_errors=True)
out.mkdir()
build_contact_matrix("contacts.pairs", out / "bal.cool", chromsizes={"chr1": 1000, "chr2": 800}, binsize=100, verbose=False)
balance_contact_matrix(out / "bal.cool", min_nnz=1, verbose=False)
from dnakit.epigenome import pileup_features

result = pileup_features(out / "bal.cool", "features.bed", flank_bp=200, verbose=False)
print(result["feature_count"], result["features_with_data"], result["snippet_size"])
```

- **Example output:**

```text
3 3 5
```

## 20) `dnakit.epigenome.calculate_logbin_expected()`

- **Purpose:** Expected contact frequency against distance, grouped into log-spaced bins with local slopes.
- **Method:** Computes the expected values with cooltools.expected_cis, then calls cooltools.api.expected.logbin_expected. Returns the binned table, the slope table and the bin edges.
- **API:** `dnakit.epigenome.calculate_logbin_expected(cool_path, *, bins_per_order_magnitude=10, min_nvalid=200, min_count=50, weight='weight', verbose=True)`.
- **Input:** `cool_path`; `bins_per_order_magnitude` (bins per decade); `min_nvalid`; `min_count`; `weight`.
- **Output:** A dictionary with `binned_rows`, `binned`, `slopes` and `bin_edges`.
- **Example code:**

```python
import shutil
from pathlib import Path
from dnakit.epigenome.backends import balance_contact_matrix, build_contact_matrix

out = Path("g10_out")
shutil.rmtree(out, ignore_errors=True)
out.mkdir()
build_contact_matrix("contacts.pairs", out / "bal.cool", chromsizes={"chr1": 1000, "chr2": 800}, binsize=100, verbose=False)
balance_contact_matrix(out / "bal.cool", min_nnz=1, verbose=False)
from dnakit.epigenome import calculate_logbin_expected

result = calculate_logbin_expected(out / "bal.cool", bins_per_order_magnitude=4, min_nvalid=1, min_count=1, verbose=False)
print(result["binned_rows"], result["bin_edges"])
```

- **Example output:**

```text
4 [1, 2, 3, 6]
```

## 21) `dnakit.epigenome.rearrange_contact_matrix()`

- **Purpose:** Reorder the bins and pixels of a cooler to follow a new viewframe order and write a new file.
- **Method:** Reads a tab-separated table with chrom, start, end, name and the new chromosome and strand columns, then calls cooltools.api.rearrange.rearrange_cooler. Regions on the `-` strand are reversed.
- **API:** `dnakit.epigenome.rearrange_contact_matrix(cool_path, view_path, output_path, *, new_chrom_column='new_chrom', orientation_column='strand', verbose=True)`.
- **Input:** `cool_path`; `view_path`: the viewframe table; `output_path` (a new file); `new_chrom_column`; `orientation_column`.
- **Output:** A dictionary with `chromosomes_before`, `chromosomes_after`, `output_bins` and `output_total_contacts`.
- **Example code:**

```python
import shutil
from pathlib import Path
from dnakit.epigenome.backends import balance_contact_matrix, build_contact_matrix

out = Path("g10_out")
shutil.rmtree(out, ignore_errors=True)
out.mkdir()
build_contact_matrix("contacts.pairs", out / "bal.cool", chromsizes={"chr1": 1000, "chr2": 800}, binsize=100, verbose=False)
balance_contact_matrix(out / "bal.cool", min_nnz=1, verbose=False)
from dnakit.epigenome import rearrange_contact_matrix

result = rearrange_contact_matrix(out / "bal.cool", "rearrange_view.tsv", out / "rearranged.cool", verbose=False)
print(result["chromosomes_before"], result["chromosomes_after"], result["output_total_contacts"])
```

- **Example output:**

```text
['chr1', 'chr2'] ['chrX'] 40
```

## 22) `dnakit.epigenome.align_windows_to_bins()`

- **Purpose:** Turn genomic loci into windows aligned to a fixed bin grid.
- **Method:** Takes `flank_bp` on each side of every locus and moves both window ends to bin edges, using cooltools.api.snipping.make_bin_aligned_windows. `lo` and `hi` are bin indices counted from `region_start_bp`.
- **API:** `dnakit.epigenome.align_windows_to_bins(binsize, chroms, centers_bp, *, flank_bp=0, region_start_bp=0, verbose=True)`.
- **Input:** `binsize`; `chroms`; `centers_bp` (the same length as chroms); `flank_bp`; `region_start_bp`.
- **Output:** A dictionary with `window_count` and `windows` (each row has `chrom`, `start`, `end`, `lo` and `hi`).
- **Example code:**

```python
from dnakit.epigenome import align_windows_to_bins

result = align_windows_to_bins(100, ["chr1"], [450], flank_bp=200, verbose=False)
window = result["windows"][0]
print(window["chrom"], window["start"], window["end"], window["lo"], window["hi"])
```

- **Example output:**

```text
chr1 200 700 2 7
```

## 23) `dnakit.epigenome.assign_features_to_view()`

- **Purpose:** Name the viewframe region that contains each BED feature.
- **Method:** Checks interval overlap with cooltools.lib.common.assign_view_auto. Features that overlap no region get `region` None and are counted in `unassigned_count`. Feature order is kept.
- **API:** `dnakit.epigenome.assign_features_to_view(features_path, view_path, *, verbose=True)`.
- **Input:** `features_path`: a BED file; `view_path`: a BED file whose fourth column is the region name.
- **Output:** A dictionary with `feature_count`, `unassigned_count` and `features` (each row has `chrom`, `start`, `end` and `region`).
- **Example code:**

```python
from dnakit.epigenome import assign_features_to_view

result = assign_features_to_view("features.bed", "view.bed", verbose=False)
for row in result["features"]:
    print(row["chrom"], row["start"], row["region"])
```

- **Example output:**

```text
chr1 300 chr1
chr1 500 chr1
chr2 200 chr2
```

## 24) `dnakit.epigenome.make_chromosome_view()`

- **Purpose:** Build a whole-chromosome viewframe (one region per chromosome) from the cooler chromosome sizes.
- **Method:** Calls cooltools.lib.common.make_cooler_view. With `ucsc_names=True` the region names use UCSC style.
- **API:** `dnakit.epigenome.make_chromosome_view(cool_path, *, ucsc_names=False, verbose=True)`.
- **Input:** `cool_path`; `ucsc_names`.
- **Output:** A dictionary with `region_count` and `regions` (each row has `chrom`, `start`, `end` and `name`).
- **Example code:**

```python
import shutil
from pathlib import Path
from dnakit.epigenome.backends import balance_contact_matrix, build_contact_matrix

out = Path("g10_out")
shutil.rmtree(out, ignore_errors=True)
out.mkdir()
build_contact_matrix("contacts.pairs", out / "bal.cool", chromsizes={"chr1": 1000, "chr2": 800}, binsize=100, verbose=False)
balance_contact_matrix(out / "bal.cool", min_nnz=1, verbose=False)
from dnakit.epigenome import make_chromosome_view

result = make_chromosome_view(out / "bal.cool", verbose=False)
for row in result["regions"]:
    print(row["name"], row["start"], row["end"])
```

- **Example output:**

```text
chr1 0 1000
chr2 0 800
```

## 25) `dnakit.epigenome.read_viewframe_bed()`

- **Purpose:** Read and check a viewframe BED file (unique names, no overlaps).
- **Method:** Calls cooltools.lib.io.read_viewframe_from_file. With `cool_path`, the regions are also checked against the chromosome lengths.
- **API:** `dnakit.epigenome.read_viewframe_bed(view_path, *, cool_path=None, verbose=True)`.
- **Input:** `view_path`: a BED file; `cool_path` (optional).
- **Output:** A dictionary with `region_count`, `regions` and `verified_against`.
- **Example code:**

```python
import shutil
from pathlib import Path
from dnakit.epigenome.backends import balance_contact_matrix, build_contact_matrix

out = Path("g10_out")
shutil.rmtree(out, ignore_errors=True)
out.mkdir()
build_contact_matrix("contacts.pairs", out / "bal.cool", chromsizes={"chr1": 1000, "chr2": 800}, binsize=100, verbose=False)
balance_contact_matrix(out / "bal.cool", min_nnz=1, verbose=False)
from dnakit.epigenome import read_viewframe_bed

result = read_viewframe_bed("view.bed", cool_path=out / "bal.cool", verbose=False)
print(result["region_count"], [row["name"] for row in result["regions"]])
```

- **Example output:**

```text
2 ['chr1', 'chr2']
```

## 26) `dnakit.epigenome.read_expected_table()`

- **Purpose:** Read an expected-contacts table in cooltools v1.0 format and check its column types.
- **Method:** Calls cooltools.lib.io.read_expected_from_file. A cis table has region1, region2, dist, n_valid and value columns. At most `limit` rows are returned.
- **API:** `dnakit.epigenome.read_expected_table(expected_path, *, contact_type='cis', limit=100, verbose=True)`.
- **Input:** `expected_path`: a tab-separated file; `contact_type` (`cis` or `trans`); `limit`.
- **Output:** A dictionary with `row_count`, `columns` and `rows`.
- **Example code:**

```python
from dnakit.epigenome import read_expected_table

result = read_expected_table("expected.tsv", limit=2, verbose=False)
print(result["row_count"], result["columns"])
print(result["rows"][1]["dist"], result["rows"][1]["balanced.avg"])
```

- **Example output:**

```text
4 ['region1', 'region2', 'dist', 'n_valid', 'count.avg', 'balanced.avg']
1 0.7
```

## 27) `dnakit.epigenome.observed_over_expected_matrix()`

- **Purpose:** Divide a dense symmetric contact matrix by its distance-dependent mean, giving an observed-over-expected (O/E) matrix.
- **Method:** Calls cooltools.lib.numutils.observed_over_expected. Diagonals are grouped into distance bins that grow exponentially, and each diagonal is divided by the mean of its bin. The input must be symmetric.
- **API:** `dnakit.epigenome.observed_over_expected_matrix(matrix, *, dist_bin_edge_ratio=1.03, verbose=True)`.
- **Input:** `matrix`: a square matrix (a list of lists); `dist_bin_edge_ratio` (greater than 1 and at most 2).
- **Output:** A dictionary with `observed_over_expected`, `distance_bin_edges`, `sum_pixels` and `n_pixels`.
- **Example code:**

```python
from dnakit.epigenome import observed_over_expected_matrix

matrix = [[1.0, 2.0, 3.0], [2.0, 1.0, 2.0], [3.0, 2.0, 1.0]]
result = observed_over_expected_matrix(matrix, verbose=False)
print([[round(value, 4) for value in row] for row in result["observed_over_expected"]])
```

- **Example output:**

```text
[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0], [1.0, 1.0, 1.0]]
```

## 28) `dnakit.epigenome.iterative_correction_dense()`

- **Purpose:** Iterative correction (matrix balancing) of a dense contact matrix, symmetric or asymmetric.
- **Method:** A symmetric matrix uses cooltools.lib.numutils.iterative_correction_symmetric, which accepts `ignore_diags`. An asymmetric matrix uses `iterative_correction_asymmetric`, which requires `ignore_diags` to be 0.
- **API:** `dnakit.epigenome.iterative_correction_dense(matrix, *, symmetric=True, max_iter=1000, ignore_diags=0, tol=1e-05, verbose=True)`.
- **Input:** `matrix`; `symmetric`; `max_iter`; `ignore_diags`; `tol`.
- **Output:** A dictionary with `converged`, `iterations`, `corrected` and `bias`.
- **Example code:**

```python
from dnakit.epigenome import iterative_correction_dense

matrix = [[1.0, 2.0, 3.0], [2.0, 4.0, 5.0], [3.0, 5.0, 6.0]]
result = iterative_correction_dense(matrix, max_iter=500, verbose=False)
print(result["converged"], [round(value, 4) for value in result["bias"]])
```

- **Example output:**

```text
True [0.5618, 1.0536, 1.3846]
```

## 29) `dnakit.epigenome.dense_cis_eigenvectors()`

- **Purpose:** Eigen-decomposition of one dense cis contact matrix, giving compartment eigenvalues and eigenvectors.
- **Method:** Calls cooltools.api.eigdecomp.cis_eig and drops the `ignore_diags` diagonals next to the main diagonal. The ARPACK start vector is fixed, and each eigenvector's largest-magnitude entry is positive.
- **API:** `dnakit.epigenome.dense_cis_eigenvectors(matrix, *, n_eigs=3, ignore_diags=2, clip_percentile=0.0, verbose=True)`.
- **Input:** `matrix`: a symmetric square matrix; `n_eigs`; `ignore_diags`; `clip_percentile`.
- **Output:** A dictionary with `eigenvalues` and `eigenvectors` (one row per eigenvector).
- **Example code:**

```python
from dnakit.epigenome import dense_cis_eigenvectors

matrix = [[4.0, 3.0, 1.0, 0.5], [3.0, 4.0, 1.0, 0.5], [1.0, 1.0, 4.0, 3.0], [0.5, 0.5, 3.0, 4.0]]
result = dense_cis_eigenvectors(matrix, n_eigs=1, ignore_diags=0, verbose=False)
print([round(value, 4) for value in result["eigenvalues"]])
```

- **Example output:**

```text
[0.8095]
```

## 30) `dnakit.epigenome.expected_from_dense_matrix()`

- **Purpose:** Count and sum the contacts on each diagonal of a dense matrix, the raw summary behind a cis expected curve.
- **Method:** Calls cooltools.api.expected.diagsum_from_array. Rows are sorted by distance.
- **API:** `dnakit.epigenome.expected_from_dense_matrix(matrix, *, ignore_diags=2, region_name=None, verbose=True)`.
- **Input:** `matrix`: a square matrix; `ignore_diags`; `region_name` (optional).
- **Output:** A dictionary with `diagonals` (each row has `dist`, `n_valid` and `balanced.sum`).
- **Example code:**

```python
from dnakit.epigenome import expected_from_dense_matrix

matrix = [[1.0, 2.0, 3.0], [2.0, 4.0, 5.0], [3.0, 5.0, 6.0]]
result = expected_from_dense_matrix(matrix, ignore_diags=0, verbose=False)
for row in result["diagonals"]:
    print(row["dist"], row["n_valid"], row["balanced.sum"])
```

- **Example output:**

```text
0 3 11.0
1 2 7.0
2 1 3.0
```

## 31) `dnakit.epigenome.digitize_track()`

- **Purpose:** Quantile-bin a genomic signal track into integer classes for saddle analysis.
- **Method:** Calls cooltools.api.saddle.digitize. Classes 1 to n cover the values inside the quantile range; 0 is below the range, n+1 is above it, and -1 marks missing values.
- **API:** `dnakit.epigenome.digitize_track(track_path, *, n_bins=5, qrange=(0.02, 0.98), verbose=True)`.
- **Input:** `track_path`: a file with chrom, start, end and value columns; `n_bins`; `qrange` (low and high quantiles).
- **Output:** A dictionary with `intervals` (each row has `chrom`, `start`, `end` and `class`) and `bin_edges`.
- **Example code:**

```python
from dnakit.epigenome import digitize_track

result = digitize_track("track.bed", n_bins=3, qrange=(0.0, 1.0), verbose=False)
print([row["class"] for row in result["intervals"][:8]])
```

- **Example output:**

```text
[1, 1, 1, 2, 2, 2, 2, 3]
```

## pairtools functions (11)

## 32) `dnakit.epigenome.find_restriction_fragment()`

- **Purpose:** Find the restriction fragment that contains a position.
- **Method:** Binary search over sorted cut sites, using pairtools.lib.restrict.find_rfrag. A fragment is the half-open interval between two consecutive sites. Positions outside the site range are rejected.
- **API:** `dnakit.epigenome.find_restriction_fragment(sites, chrom, position, *, verbose=True)`.
- **Input:** `sites`: a mapping from chromosome name to ascending site positions; `chrom`; `position`.
- **Output:** A dictionary with `fragment_index`, `start`, `end` and `length`.
- **Example code:**

```python
from dnakit.epigenome import find_restriction_fragment

sites = {"chr1": [0, 300, 700, 1000], "chr2": [0, 800]}
result = find_restriction_fragment(sites, "chr1", 450, verbose=False)
print(result["fragment_index"], result["start"], result["end"], result["length"])
```

- **Example output:**

```text
1 300 700 400
```

## 33) `dnakit.epigenome.compute_contact_scaling()`

- **Purpose:** Contact frequency against genomic distance (scaling curve) for a .pairs file.
- **Method:** Calls pairtools.lib.scaling.compute_scaling and adds up the results over regions and strand pairs for each distance bin. The contact frequency is the number of contacts divided by the number of possible base-pair positions.
- **API:** `dnakit.epigenome.compute_contact_scaling(pairs_path, *, chromsizes, min_distance=1, max_distance=1000000000, bins_per_decade=8, ignore_trans=False, verbose=True)`.
- **Input:** `pairs_path`; `chromsizes`; `min_distance`; `max_distance`; `bins_per_decade`; `ignore_trans`.
- **Output:** A dictionary with `cis_pairs`, `trans_pairs` and `bins` (each row has `min_dist`, `max_dist`, `n_pairs` and `contact_frequency`).
- **Example code:**

```python
from dnakit.epigenome import compute_contact_scaling

result = compute_contact_scaling("contacts.pairs", chromsizes={"chr1": 1000, "chr2": 800}, max_distance=1000, bins_per_decade=2, verbose=False)
print(result["cis_pairs"], result["trans_pairs"])
print(result["bins"][1]["min_dist"], result["bins"][1]["max_dist"], result["bins"][1]["n_pairs"])
```

- **Example output:**

```text
36 4
1 3 1
```

## 34) `dnakit.epigenome.normalize_scaling_curve()`

- **Purpose:** Divide a contact-frequency curve by its mean inside a distance window, so that the window mean is 1.
- **Method:** Calls pairtools.lib.scaling.norm_scaling_factor and norm_scaling. Distances must be unique and ascending.
- **API:** `dnakit.epigenome.normalize_scaling_curve(distances, frequencies, *, norm_window=(10000, 1000000000), verbose=True)`.
- **Input:** `distances`: a list of integers; `frequencies`: non-negative values of the same length; `norm_window`: the low and high distances.
- **Output:** A dictionary with `factor` and `points` (each row has `distance` and `frequency`).
- **Example code:**

```python
from dnakit.epigenome import normalize_scaling_curve

result = normalize_scaling_curve([1, 10, 100], [0.5, 0.25, 0.1], norm_window=(1, 100), verbose=False)
print(round(result["factor"], 4))
print([round(point["frequency"], 4) for point in result["points"]])
```

- **Example output:**

```text
0.2833
[1.7647, 0.8824, 0.3529]
```

## 35) `dnakit.epigenome.assign_genomic_regions()`

- **Purpose:** Decide which region contains each position, for example a stretch of a chromosome.
- **Method:** Calls pairtools.lib.scaling.assign_regs. The region index counts regions of the same chromosome after sorting by start. Positions in no region get None.
- **API:** `dnakit.epigenome.assign_genomic_regions(positions, regions, *, verbose=True)`.
- **Input:** `positions`: a list of mappings with `chrom` and `pos`; `regions`: a list of mappings with `chrom`, `start` and `end`.
- **Output:** A dictionary with `assigned_count` and `positions` (each row has `region_index`, `region_start` and `region_end`).
- **Example code:**

```python
from dnakit.epigenome import assign_genomic_regions

positions = [{"chrom": "chr1", "pos": 50}, {"chrom": "chr1", "pos": 1500}, {"chrom": "chr2", "pos": 10}, {"chrom": "chr3", "pos": 5}]
regions = [{"chrom": "chr1", "start": 0, "end": 1000}, {"chrom": "chr1", "start": 1000, "end": 2000}, {"chrom": "chr2", "start": 0, "end": 800}]
result = assign_genomic_regions(positions, regions, verbose=False)
print([row["region_index"] for row in result["positions"]], result["assigned_count"])
```

- **Example output:**

```text
[0, 1, 0, None] 3
```

## 36) `dnakit.epigenome.find_read_mismatches()`

- **Purpose:** List the read bases that differ from the reference, given the aligned pairs.
- **Method:** Calls pairtools.lib.parse_pysam.get_mismatches_c, which compares each read base with the aligned reference base (case-insensitive). Rows are sorted by read position and then reference position.
- **API:** `dnakit.epigenome.find_read_mismatches(sequence, qualities, aligned_pairs, *, verbose=True)`.
- **Input:** `sequence`: the read sequence (only A, C, G, T or N); `qualities`: PHRED values, one per base; `aligned_pairs`: (read position, reference position, reference base) triples.
- **Output:** A dictionary with `mismatch_count` and `mismatches` (each row has `ref_base`, `read_base`, `quality`, `ref_pos` and `read_pos`).
- **Example code:**

```python
from dnakit.epigenome import find_read_mismatches

pairs = [(0, 100, "A"), (1, 101, "C"), (2, 102, "G"), (3, 103, "T")]
result = find_read_mismatches("ACTT", [30, 31, 32, 33], pairs, verbose=False)
print(result["mismatch_count"], result["mismatches"][0])
```

- **Example output:**

```text
1 {'ref_base': 'G', 'read_base': 'T', 'quality': 32, 'ref_pos': 102, 'read_pos': 2}
```

## 37) `dnakit.epigenome.parse_sam_alignments()`

- **Purpose:** Parse the first records of a SAM or BAM file into strand, 5' end and 3' end coordinates.
- **Method:** Reads the file with pairtools' AlignmentFilePairtoolized and calls pairtools.lib.parse.parse_pysam_entry on each record. The input must have @SQ header lines.
- **API:** `dnakit.epigenome.parse_sam_alignments(sam_path, *, min_mapq=30, max_records=100, verbose=True)`.
- **Input:** `sam_path`; `min_mapq`; `max_records` (1 to 100000).
- **Output:** A dictionary with `record_count` and `alignments` (each row has `read_id`, `chrom`, `pos5`, `pos3`, `strand` and `mapq`).
- **Example code:**

```python
from dnakit.epigenome import parse_sam_alignments

result = parse_sam_alignments("sample.sam", min_mapq=30, max_records=2, verbose=False)
for row in result["alignments"]:
    print(row["read_id"], row["chrom"], row["pos5"], row["pos3"], row["strand"])
```

- **Example output:**

```text
r1 chr1 101 110 +
r1 chr2 301 310 +
```

## 38) `dnakit.epigenome.read_pairs_table()`

- **Purpose:** Read a .pairs file: header, chromosome sizes, row count and the first rows.
- **Method:** Calls pairtools.lib.pairsio.read_pairs in a single process. `row_count` counts every data line, while `rows` holds only the first `limit` lines.
- **API:** `dnakit.epigenome.read_pairs_table(pairs_path, *, limit=100, verbose=True)`.
- **Input:** `pairs_path`: a .pairs file; `limit` (1 to 1000).
- **Output:** A dictionary with `row_count`, `columns`, `chromsizes`, `header` and `rows`.
- **Example code:**

```python
from dnakit.epigenome import read_pairs_table

result = read_pairs_table("contacts.pairs", limit=2, verbose=False)
print(result["row_count"], result["chromsizes"])
print(result["rows"][0]["readID"], result["rows"][0]["chrom1"], result["rows"][0]["pos1"])
```

- **Example output:**

```text
40 {'chr1': 1000, 'chr2': 800}
r0 chr2 649
```

## 39) `dnakit.epigenome.make_pairs_header()`

- **Purpose:** Build the standard .pairs header lines for a genome.
- **Method:** Calls pairtools.lib.headerops.make_standard_pairsheader, which writes the version, shape, assembly name, chromosome lengths and column names.
- **API:** `dnakit.epigenome.make_pairs_header(*, chromsizes, assembly='unknown', columns=None, shape='upper triangle', verbose=True)`.
- **Input:** `chromsizes` (required); `assembly`; `columns` (optional); `shape` (`upper triangle` or `whole matrix`).
- **Output:** A dictionary with `chromosome_count` and `lines`.
- **Example code:**

```python
from dnakit.epigenome import make_pairs_header

result = make_pairs_header(chromsizes={"chr1": 1000, "chr2": 800}, assembly="hg38", verbose=False)
for line in result["lines"][:4]:
    print(line)
```

- **Example output:**

```text
## pairs format v1.0.0
#shape: upper triangle
#genome_assembly: hg38
#chromsize: chr1 1000
```

## 40) `dnakit.epigenome.subset_pairs_header()`

- **Purpose:** Keep only the chromsize and chromosomes header lines for the chosen chromosomes.
- **Method:** Calls pairtools.lib.headerops.subset_chroms_in_pairsheader. All other header lines are kept unchanged.
- **API:** `dnakit.epigenome.subset_pairs_header(header_lines, chroms, *, verbose=True)`.
- **Input:** `header_lines`: a list of header lines; `chroms`: the chromosome names to keep.
- **Output:** A dictionary with `input_line_count`, `output_line_count` and `lines`.
- **Example code:**

```python
from dnakit.epigenome import make_pairs_header, subset_pairs_header

header = make_pairs_header(chromsizes={"chr1": 1000, "chr2": 800}, verbose=False)["lines"]
result = subset_pairs_header(header, ["chr2"], verbose=False)
print(result["input_line_count"], result["output_line_count"])
print([line for line in result["lines"] if line.startswith("#chromsize")])
```

- **Example output:**

```text
6 5
['#chromsize: chr2 800']
```

## 41) `dnakit.epigenome.read_chrom_order()`

- **Purpose:** Number each chromosome in a chromsizes file in file order, starting at 1.
- **Method:** Calls pairtools.lib.headerops.get_chrom_order. With `sam_chroms`, only those chromosomes are numbered; the others are added at the end in sorted order.
- **API:** `dnakit.epigenome.read_chrom_order(chroms_path, *, sam_chroms=None, verbose=True)`.
- **Input:** `chroms_path`: a two-column (name, length) file; `sam_chroms` (optional).
- **Output:** A dictionary with `chromosome_count` and `chromosomes` (each row has `chrom` and `order`).
- **Example code:**

```python
from dnakit.epigenome import read_chrom_order

result = read_chrom_order("chrom.sizes", sam_chroms=["chr2", "chr1", "chrC"], verbose=False)
print([(row["chrom"], row["order"]) for row in result["chromosomes"]])
```

- **Example output:**

```text
[('chr1', 1), ('chr2', 2), ('chrC', 3)]
```

## 42) `dnakit.epigenome.assign_haplotype_phase()`

- **Purpose:** Decide which haplotype a read belongs to from its XA or XB alternative-hit tag.
- **Method:** Calls pairtools.lib.phase.phase_side_XB (`variant='XB'`) or phase_side_XA (`variant='XA'`, which needs `mismatches`). The suffix of the chromosome name, such as `_1` or `_2`, sets the haplotype number.
- **API:** `dnakit.epigenome.assign_haplotype_phase(chrom, *, tag_value, primary_score, secondary_score, variant='XB', mismatches=None, phase_suffixes=('_1', '_2'), verbose=True)`.
- **Input:** `chrom`; `tag_value`: the tag text; `primary_score` and `secondary_score`: the primary and secondary alignment scores; `variant`; `mismatches`; `phase_suffixes`.
- **Output:** A dictionary with `phase` (`0`, `1`, `.` or `!`), `chromosome_base` and `scores`.
- **Example code:**

```python
from dnakit.epigenome import assign_haplotype_phase

result = assign_haplotype_phase("chr1_1", tag_value="", primary_score=50, secondary_score=40, variant="XB", verbose=False)
print(result["phase"], result["chromosome_base"], result["scores"])
```

- **Example output:**

```text
0 chr1 [50, 40, -1]
```
