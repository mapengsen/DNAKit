# G18l: review functions from cooler, cooltools and pairtools

This group reviews the functions that were `other` or `cap` in the cooler, cooltools and pairtools ledgers. It wraps 24 functions that work on DNA or genomic coordinates and run offline on local inputs: 4 from cooler, 16 from cooltools and 4 from pairtools. The remaining functions are recorded in the ledger with a reason, such as `cli`, `internal`, `not-dna`, `needs-network`, `family-variant` or `duplicate-of`.
The functions live in the three review modules of `dnakit.epigenome` (`cooler_review_functions.py`, `cooltools_review_functions.py` and `pairtools_review_functions.py`). Each one imports its upstream package only when called; a missing optional dependency raises `BackendUnavailableError`. Inputs are checked first, and upstream `ValueError`, `KeyError` and `OSError` are reported as `ConfigurationError`. Every public function takes a keyword-only `verbose=True`, returns a plain dictionary and includes `provenance`.
The eigen-decomposition and decoy functions (`dense_trans_eigenvectors`, `replace_cis_with_decoys`) save and restore the numpy global state and the C library `random()` state during the call (glibc and BSD provide `initstate` and `setstate`), and they fix the ARPACK start vector. The patch is restored in `finally` but is not thread-safe. On platforms without `initstate` and `setstate` the C random state is not isolated. `call_chromatin_dots` needs an explicit kernel width at fine resolutions, so the example sets it.
The examples share one synthetic cooler (`synthetic.cool`, two chromosomes of 2000 bp and 1500 bp at 100 bp resolution, created by the setup example). Outputs are written to the current directory. The example output is produced by running the code.
## Example files

The files below are written to the example directory before the examples run.

**Example file:** `track.bedgraph`

```text
chr1	0	100	1.0
chr1	100	200	2.0
chr1	200	300	3.0
chr2	0	100	0.5
chr2	100	200	0.6
```

**Example file:** `features.bed`

```text
chr1	300	500
chr2	100	200
```

**Example file:** `features.bedpe`

```text
chr1	100	150	chr1	700	750
chr1	500	550	chr2	50	80
```

**Example file:** `view.bed`

```text
chr1	0	400	A
chr1	600	1000	B
chr2	0	800	X
```

**Example file:** `cells.pairs`

```text
## pairs format v1.0
#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type
r1	chr1	100	chr1	300	+	-	UU
r2	chr1	100	chr1	300	+	-	UU
r3	chr2	10	chr2	20	+	-	UU
r4	chr1	500	chr2	100	+	-	UU
```

**Example file:** `cov.pairs`

```text
## pairs format v1.0
#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type
r1	chr1	100	chr1	300	+	-	UU
r2	chr1	100	chr1	300	+	-	UU
r3	chr1	100	chr1	300	+	-	UU
r4	chr2	10	chr2	20	+	-	UU
r5	chr1	500	chr2	100	+	-	UU
```

**Example file:** `reads.sam`

```text
@HD	VN:1.6	SO:queryname
@SQ	SN:chr1	LN:1000
@SQ	SN:chr2	LN:800
r1	65	chr1	101	60	10M	=	101	0	ACGTACGTAC	IIIIIIIIII	NM:i:0	AS:i:10	XS:i:0
r1	129	chr2	301	60	10M	=	301	0	ACGTACGTAC	IIIIIIIIII	NM:i:0	AS:i:10	XS:i:0
r2	65	chr2	51	60	10M	=	51	0	ACGTACGTAC	IIIIIIIIII	NM:i:0	AS:i:10	XS:i:0
r2	129	chr2	201	60	10M	=	201	0	ACGTACGTAC	IIIIIIIIII	NM:i:0	AS:i:10	XS:i:0
```


## Setup

This example creates `synthetic.cool` and removes outputs left by earlier runs.

**Example code:**

```python
import os

import cooler
import numpy as np
import pandas as pd

for name in ["synthetic.cool", "cells.scool", "all.pairs", "strict.pairs", "removed.pairs",
             "reads.pairs"]:
    if os.path.exists(name):
        os.remove(name)
rng = np.random.default_rng(0)
bins = cooler.binnify(pd.Series({"chr1": 2000, "chr2": 1500}), 100)
rows = []
for i in range(len(bins)):
    for j in range(i, len(bins)):
        same = bins.chrom[i] == bins.chrom[j]
        count = rng.poisson(6.0 / (1 + abs(i - j)) if same else 0.5)
        if i == j:
            count += 3
        if count > 0:
            rows.append((i, j, int(count)))
rows.append((2, 9, 40))
pixels = pd.DataFrame(rows, columns=["bin1_id", "bin2_id", "count"])
pixels = pixels.sort_values(["bin1_id", "bin2_id"])
cooler.create_cooler("synthetic.cool", bins, pixels, dtypes={"count": "int"}, ensure_sorted=True)
cooler.balance_cooler(cooler.Cooler("synthetic.cool"), store=True)
clr = cooler.Cooler("synthetic.cool")
print("bins", len(clr.bins()[:]), "chromosomes", list(clr.chromnames))
```

**Example output:**

```text
bins 35 chromosomes ['chr1', 'chr2']
```

## cooler functions (4)

## 1) `dnakit.epigenome.create_single_cell_cooler()`

- **Purpose:** Write a single-cell (.scool) cooler from one .pairs file per cell.
- **Method:** Each cell's .pairs file is binned at `binsize`, upper-triangle pixels are summed per bin pair, and cooler.create_scool writes the file. The output is a new file and is never overwritten.
- **API:** `dnakit.epigenome.create_single_cell_cooler(output_path, pairs_by_cell, chromsizes, binsize, *, verbose=True)`.
- **Input:** `output_path`: a new .scool file; `pairs_by_cell`: cell name to .pairs file; `chromsizes`: chromosome lengths; `binsize`: bin size.
- **Output:** A dictionary with `bin_count`, `cell_count` and `cells` (each with `pair_count` and `nonzero_pixels`), plus `provenance`.
- **Example code:**

```python
from dnakit.epigenome import create_single_cell_cooler

result = create_single_cell_cooler(
    "cells.scool",
    {"cellA": "cells.pairs", "cellB": "cells.pairs"},
    {"chr1": 1000, "chr2": 800},
    100,
    verbose=False,
)
print(result["operation"], result["bin_count"], result["cell_count"])
for cell in result["cells"]:
    print(cell["name"], cell["pair_count"], cell["nonzero_pixels"])
```

- **Example output:**

```text
create_single_cell_cooler 18 2
cellA 4 3
cellB 4 3
```

## 2) `dnakit.epigenome.list_cooler_groups()`

- **Purpose:** List the groups that hold a cooler data collection in an HDF5 file.
- **Method:** Calls cooler.api.list_coolers; the result is sorted by group path. A plain .cool file gives `/`, and a multi-resolution file gives one group per resolution.
- **API:** `dnakit.epigenome.list_cooler_groups(path, *, verbose=True)`.
- **Input:** `path`: an HDF5 file (.cool, .mcool or .scool).
- **Output:** A dictionary with `group_count` and `groups`, plus `provenance`.
- **Example code:**

```python
from dnakit.epigenome import list_cooler_groups

print(list_cooler_groups("synthetic.cool", verbose=False)["groups"])
print(list_cooler_groups("cells.scool", verbose=False)["groups"])
```

- **Example output:**

```text
['/']
['/cells/cellA', '/cells/cellB']
```

## 3) `dnakit.epigenome.list_scool_cells()`

- **Purpose:** List the cells of a single-cell (.scool) cooler file.
- **Method:** Checks the file type with cooler.fileops.is_scool_file and rejects other files; then lists the cell groups sorted by path and gives each cell name.
- **API:** `dnakit.epigenome.list_scool_cells(path, *, verbose=True)`.
- **Input:** `path`: a single-cell cooler file.
- **Output:** A dictionary with `cell_count` and `cells` (each with `name` and `group`), plus `provenance`.
- **Example code:**

```python
from dnakit.epigenome import list_scool_cells

result = list_scool_cells("cells.scool", verbose=False)
print(result["cell_count"], [cell["name"] for cell in result["cells"]])
```

- **Example output:**

```text
2 ['cellA', 'cellB']
```

## 4) `dnakit.epigenome.infer_cooler_binsize()`

- **Purpose:** Infer the bin size of a cooler from its bin table.
- **Method:** Calls cooler.util.get_binsize: every bin except the last bin of each chromosome must have the same length. When the bins are not uniform, `binsize` is None.
- **API:** `dnakit.epigenome.infer_cooler_binsize(cool_path, *, verbose=True)`.
- **Input:** `cool_path`: a cooler file.
- **Output:** A dictionary with `binsize`, `uniform` and `bin_count`, plus `provenance`.
- **Example code:**

```python
from dnakit.epigenome import infer_cooler_binsize

result = infer_cooler_binsize("synthetic.cool", verbose=False)
print(result["binsize"], result["uniform"], result["bin_count"])
```

- **Example output:**

```text
100 True 35
```

## cooltools functions (16)

## 5) `dnakit.epigenome.call_chromatin_dots()`

- **Purpose:** Call chromatin dots with HiCCUPS-style kernels of a chosen width.
- **Method:** The cis expected values come from cooltools.expected_cis (balanced.avg), four kernels are built with numutils.get_kernel, and cooltools.dots runs. Default kernels are not recommended at fine resolutions, so the kernel width is given explicitly.
- **API:** `dnakit.epigenome.call_chromatin_dots(cool_path, *, weight='weight', kernel_half_width=3, donut_inner_half_width=1, max_loci_separation=10000000, n_lambda_bins=40, verbose=True)`.
- **Input:** `cool_path`: a balanced cooler; `kernel_half_width` and `donut_inner_half_width`: outer and inner half-widths of the kernels.
- **Output:** A dictionary with `dot_count` and `dots` (coordinates, count and per-kernel la_exp values and q-values), plus `provenance`.
- **Example code:**

```python
from dnakit.epigenome import call_chromatin_dots

result = call_chromatin_dots("synthetic.cool", verbose=False)
print(result["kernel_half_width"], result["donut_inner_half_width"], result["dot_count"])
```

- **Example output:**

```text
3 1 0
```

## 6) `dnakit.epigenome.align_track_to_cooler_bins()`

- **Purpose:** Align a bedGraph-like track to the bins of a cooler.
- **Method:** Calls cooltools.lib.common.align_track_with_cooler: it checks the bin size, merges the track with the bin table, and sets bad bins (non-finite weight) to None unless masking is turned off.
- **API:** `dnakit.epigenome.align_track_to_cooler_bins(track_path, cool_path, *, weight='weight', mask_bad_bins=True, verbose=True)`.
- **Input:** `track_path`: a TSV with chrom, start, end and value; `cool_path`: a balanced cooler.
- **Output:** A dictionary with `bin_count`, `valued_bin_count` and `bins`, plus `provenance`.
- **Example code:**

```python
from dnakit.epigenome import align_track_to_cooler_bins

result = align_track_to_cooler_bins("track.bedgraph", "synthetic.cool", verbose=False)
print(result["bin_count"], result["valued_bin_count"])
for row in result["bins"]:
    if row["value"] is not None:
        print(row["chrom"], row["start"], row["end"], row["value"])
```

- **Example output:**

```text
35 5
chr1 0 100 1.0
chr1 100 200 2.0
chr1 200 300 3.0
chr2 0 100 0.5
chr2 100 200 0.6
```

## 7) `dnakit.epigenome.make_track_view()`

- **Purpose:** Build a viewframe with one region per chromosome covered by a track.
- **Method:** Calls cooltools.lib.common.view_from_track: each chromosome takes its smallest start and largest end, and bioframe supplies the region names.
- **API:** `dnakit.epigenome.make_track_view(track_path, *, verbose=True)`.
- **Input:** `track_path`: a TSV with chrom, start and end.
- **Output:** A dictionary with `region_count` and `regions` (chrom, start, end and name), plus `provenance`.
- **Example code:**

```python
from dnakit.epigenome import make_track_view

result = make_track_view("track.bedgraph", verbose=False)
for row in result["regions"]:
    print(row["chrom"], row["start"], row["end"], row["name"])
```

- **Example output:**

```text
chr1 0 300 chr1
chr2 0 200 chr2
```

## 8) `dnakit.epigenome.mask_track_bad_bins()`

- **Purpose:** Set track values to missing for bins that the cooler marks as bad.
- **Method:** Calls cooltools.lib.common.mask_cooler_bad_bins: only rows whose chrom, start and end match a bin are kept, and bins with a non-finite weight get NaN.
- **API:** `dnakit.epigenome.mask_track_bad_bins(track_path, cool_path, *, weight='weight', verbose=True)`.
- **Input:** `track_path`: a TSV; `cool_path`: a cooler file; `weight`: the weight column name.
- **Output:** A dictionary with `row_count`, `masked_count` and `rows`, plus `provenance`.
- **Example code:**

```python
from dnakit.epigenome import mask_track_bad_bins

result = mask_track_bad_bins("track.bedgraph", "synthetic.cool", verbose=False)
print(result["row_count"], result["masked_count"])
```

- **Example output:**

```text
5 0
```

## 9) `dnakit.epigenome.assign_paired_features_to_view()`

- **Purpose:** Name the view region of each side of paired features.
- **Method:** Reads a six-column BEDPE file and a BED view, then calls cooltools.lib.common.assign_view_paired. A side that lies in no region is None.
- **API:** `dnakit.epigenome.assign_paired_features_to_view(features_path, view_path, *, drop_unassigned=False, verbose=True)`.
- **Input:** `features_path`: six-column BEDPE; `view_path`: BED view (a name column is used when present); `drop_unassigned`: whether to drop unnamed features.
- **Output:** A dictionary with `feature_count`, `assigned_count` and `features` (with `region1` and `region2`), plus `provenance`.
- **Example code:**

```python
from dnakit.epigenome import assign_paired_features_to_view

result = assign_paired_features_to_view("features.bedpe", "view.bed", verbose=False)
for row in result["features"]:
    print(row["chrom1"], row["start1"], row["chrom2"], row["start2"],
          row["region1"], row["region2"])
```

- **Example output:**

```text
chr1 100 chr1 700 A B
chr1 500 chr2 50 None X
```

## 10) `dnakit.epigenome.assign_bins_to_regions()`

- **Purpose:** Name the region that contains each bin index.
- **Method:** Sorts the regions by bin_start, rejects overlapping regions, then calls cooltools.lib.common.assign_regions_to_bins. A bin outside every region is None.
- **API:** `dnakit.epigenome.assign_bins_to_regions(bin_ids, regions, *, verbose=True)`.
- **Input:** `bin_ids`: a list of bin indices; `regions`: a list of mappings with name, bin_start and bin_end.
- **Output:** A dictionary with `bin_count`, `assigned_count` and `assignments` (bin_id and region), plus `provenance`.
- **Example code:**

```python
from dnakit.epigenome import assign_bins_to_regions

result = assign_bins_to_regions(
    [0, 4, 5, 11, 12],
    [
        {"name": "a", "bin_start": 0, "bin_end": 5},
        {"name": "b", "bin_start": 5, "bin_end": 12},
    ],
    verbose=False,
)
print([(row["bin_id"], row["region"]) for row in result["assignments"]])
```

- **Example output:**

```text
[(0, 'a'), (4, 'a'), (5, 'b'), (11, 'b'), (12, None)]
```

## 11) `dnakit.epigenome.dense_trans_eigenvectors()`

- **Purpose:** Compute trans compartment eigenvectors of a dense genome-wide matrix.
- **Method:** Calls cooltools.api.eigdecomp.trans_eig: cis blocks are removed, extreme trans values are clipped, fake cis values are filled in, and the matrix is decomposed. Each eigenvector's largest absolute entry is positive. The numpy global state and the C random() state are saved and restored around the call, and ARPACK gets a fixed start vector.
- **API:** `dnakit.epigenome.dense_trans_eigenvectors(matrix, partition, *, n_eigs=3, perc_top=99.95, perc_bottom=1.0, verbose=True)`.
- **Input:** `matrix`: a symmetric square matrix; `partition`: chromosome offsets from 0 to the matrix size; `n_eigs`: number of eigenvectors.
- **Output:** A dictionary with `eigenvalues` and `eigenvectors`, plus `provenance`.
- **Example code:**

```python
import numpy as np

from dnakit.epigenome import dense_trans_eigenvectors

rng = np.random.default_rng(5)
base = rng.random((12, 12))
matrix = (base + base.T).tolist()
result = dense_trans_eigenvectors(matrix, [0, 6, 12], n_eigs=2, verbose=False)
print(result["size"], result["partition"])
print([round(value, 4) for value in result["eigenvalues"]])
```

- **Example output:**

```text
12 [0, 6, 12]
[2.5948, -2.4452]
```

## 12) `dnakit.epigenome.pairwise_diagonal_sums()`

- **Purpose:** Sum the pixels of each diagonal within every chromosome block.
- **Method:** Calls cooltools.api.expected.diagsum_pairwise and sorts by region pair and distance. The two closest diagonals are ignored by default.
- **API:** `dnakit.epigenome.pairwise_diagonal_sums(cool_path, *, weight='weight', ignore_diags=2, verbose=True)`.
- **Input:** `cool_path`: a balanced cooler; `weight` and `ignore_diags`.
- **Output:** A dictionary with `row_count` and `rows` (region1, region2, dist, n_total, n_valid and count.sum), plus `provenance`.
- **Example code:**

```python
from dnakit.epigenome import pairwise_diagonal_sums

result = pairwise_diagonal_sums("synthetic.cool", verbose=False)
chr1_rows = [row for row in result["rows"] if row["region1"] == "chr1"]
print(result["row_count"], len(chr1_rows))
```

- **Example output:**

```text
35 20
```

## 13) `dnakit.epigenome.pairwise_block_sums()`

- **Purpose:** Sum the pixels of the rectangular block for every pair of chromosomes.
- **Method:** Calls cooltools.api.expected.blocksum_pairwise; the rows are sorted by region pair.
- **API:** `dnakit.epigenome.pairwise_block_sums(cool_path, *, weight='weight', verbose=True)`.
- **Input:** `cool_path`: a balanced cooler; `weight`.
- **Output:** A dictionary with `row_count` and `rows` (region1, region2, n_valid and count.sum), plus `provenance`.
- **Example code:**

```python
from dnakit.epigenome import pairwise_block_sums

result = pairwise_block_sums("synthetic.cool", verbose=False)
row = result["rows"][0]
print(result["row_count"], row["region1"], row["region2"], row["n_valid"], row["count.sum"])
```

- **Example output:**

```text
1 chr1 chr2 280 134
```

## 14) `dnakit.epigenome.smooth_contact_decay()`

- **Purpose:** Smooth the balanced contact-versus-distance curve in log space.
- **Method:** Builds the cis contact-versus-distance table with cooltools.expected_cis, then calls genomewide_smooth_cvd (mode="genomewide", all regions pooled) or per_region_smooth_cvd (mode="per_region", one region at a time). The suffix is always .smoothed.
- **API:** `dnakit.epigenome.smooth_contact_decay(cool_path, *, mode='genomewide', weight='weight', sigma_log10=0.1, window_sigma=5, points_per_sigma=10, verbose=True)`.
- **Input:** `cool_path`: a balanced cooler; `mode`; `sigma_log10`, `window_sigma` and `points_per_sigma`: Gaussian kernel settings.
- **Output:** A dictionary with `row_count` and `rows` (regions, distance and smoothed columns), plus `provenance`.
- **Example code:**

```python
from dnakit.epigenome import smooth_contact_decay

result = smooth_contact_decay("synthetic.cool", verbose=False)
row = result["rows"][2]
print(result["mode"], result["row_count"], row["dist"], round(row["balanced.avg.smoothed"], 6))
```

- **Example output:**

```text
genomewide 35 2 0.07744
```

## 15) `dnakit.epigenome.saddle_compartment_strength()`

- **Purpose:** Compute compartment strength from saddle sums and counts.
- **Method:** Calls cooltools.api.saddle.saddle_strength: for k bins at each corner it takes the ratio of the same-compartment mean (AA+BB) to the different-compartment mean (AB+BA).
- **API:** `dnakit.epigenome.saddle_compartment_strength(saddle_sums, saddle_counts, *, verbose=True)`.
- **Input:** `saddle_sums` and `saddle_counts`: square matrices of the same shape (None allowed).
- **Output:** A dictionary with `size` and `strengths` (each with k and strength), plus `provenance`.
- **Example code:**

```python
from dnakit.epigenome import saddle_compartment_strength

sums = [[1, 2, 3, 4], [2, 5, 6, 7], [3, 6, 8, 9], [4, 7, 9, 10]]
counts = [[2, 2, 2, 2] for _ in range(4)]
result = saddle_compartment_strength(sums, counts, verbose=False)
print([(row["k"], row["strength"]) for row in result["strengths"]])
```

- **Example output:**

```text
[(1, 1.375), (2, 1.15), (3, 1.03)]
```

## 16) `dnakit.epigenome.expand_features_to_bins()`

- **Purpose:** Expand BED or BEDPE features by a flank and snap them to bin edges.
- **Method:** Calls cooltools.api.snipping.expand_align_features: the centre is floored to the resolution, flank//resolution bins are added on each side, and the edges are multiplied back by the resolution.
- **API:** `dnakit.epigenome.expand_features_to_bins(features_path, *, flank, resolution, format='bed', verbose=True)`.
- **Input:** `features_path`: BED or BEDPE; `flank`: flank length in bp; `resolution`: resolution in bp; `format`: `bed` or `bedpe`.
- **Output:** A dictionary with `feature_count` and `features` (original and expanded coordinates), plus `provenance`.
- **Example code:**

```python
from dnakit.epigenome import expand_features_to_bins

result = expand_features_to_bins("features.bed", flank=100, resolution=50, verbose=False)
for row in result["features"]:
    print(row["chrom"], row["start"], row["end"])
```

- **Example output:**

```text
chr1 300 550
chr2 50 300
```

## 17) `dnakit.epigenome.adaptive_coarsegrain_matrix()`

- **Purpose:** Adaptively coarse-grain a contact matrix: sparse pixels are pooled with their neighbours until their counts are high enough.
- **Method:** Calls cooltools.lib.numutils.adaptive_coarsegrain: 2×2 pooling repeats until each pixel has at least `cutoff` raw counts or `max_levels` is reached. The output has the input shape.
- **API:** `dnakit.epigenome.adaptive_coarsegrain_matrix(matrix, counts, *, cutoff=5.0, max_levels=8, min_shape=8, verbose=True)`.
- **Input:** `matrix`: a contact matrix; `counts`: raw counts of the same shape; `cutoff`, `max_levels` and `min_shape`.
- **Output:** A dictionary with `size` and `matrix`, plus `provenance`.
- **Example code:**

```python
import numpy as np

from dnakit.epigenome import adaptive_coarsegrain_matrix

rng = np.random.default_rng(3)
matrix = rng.random((9, 9))
counts = np.full((9, 9), 2.0)
result = adaptive_coarsegrain_matrix(
    matrix.tolist(), counts.tolist(), cutoff=3, max_levels=2, min_shape=4, verbose=False
)
print(result["size"], len(result["matrix"]), len(result["matrix"][0]))
```

- **Example output:**

```text
9 9 9
```

## 18) `dnakit.epigenome.replace_cis_with_decoys()`

- **Purpose:** Replace the masked cis pixels of a contact matrix with decoy contacts.
- **Method:** Calls cooltools.lib.numutils.fake_cis: pixels with mask 1 take values from the C random() stream, and pixels with mask 0 stay unchanged. The C stream is switched to a private state seeded with `seed` for the call and restored afterwards.
- **API:** `dnakit.epigenome.replace_cis_with_decoys(matrix, mask, *, seed=0, verbose=True)`.
- **Input:** `matrix`: a square matrix; `mask`: a 0/1 matrix of the same shape; `seed`: random seed.
- **Output:** A dictionary with `replaced_pixels`, `seed` and `matrix`, plus `provenance`.
- **Example code:**

```python
import numpy as np

from dnakit.epigenome import replace_cis_with_decoys

rng = np.random.default_rng(4)
base = rng.random((12, 12))
matrix = (base + base.T).tolist()
mask = [[1 if i < 6 and j < 6 else 0 for j in range(12)] for i in range(12)]
result = replace_cis_with_decoys(matrix, mask, seed=0, verbose=False)
print(result["replaced_pixels"], result["seed"])
print(result["matrix"][6][6] == matrix[6][6], result["matrix"][0][0] == matrix[0][0])
```

- **Example output:**

```text
36 0
True False
```

## 19) `dnakit.epigenome.expected_full_table()`

- **Purpose:** Compute expected contacts for every pair of chromosome regions (distance decay in cis, block means in trans).
- **Method:** Calls expected_full from the cooltools sandbox (expected_full_fast when fast=True). The rows are sorted by region pair and distance.
- **API:** `dnakit.epigenome.expected_full_table(cool_path, *, weight='weight', fast=False, verbose=True)`.
- **Input:** `cool_path`: a balanced cooler; `weight`; `fast`: whether to use the faster variant.
- **Output:** A dictionary with `row_count` and `rows` (region1, region2, dist, n_valid, count.sum and expected), plus `provenance`.
- **Example code:**

```python
from dnakit.epigenome import expected_full_table

result = expected_full_table("synthetic.cool", verbose=False)
print(result["row_count"], result["rows"][0]["region1"], result["rows"][0]["region2"])
```

- **Example output:**

```text
36 chr1 chr1
```

## 20) `dnakit.epigenome.observed_over_expected_pixels()`

- **Purpose:** Compute observed over expected for every pixel of a cooler.
- **Method:** Computes the expected table with expected_full, labels each bin with its view region, and calls the sandbox obs_over_exp: oe = balanced count / expected.
- **API:** `dnakit.epigenome.observed_over_expected_pixels(cool_path, *, weight='weight', verbose=True)`.
- **Input:** `cool_path`: a balanced cooler; `weight`.
- **Output:** A dictionary with `pixel_count` and `pixels` (bin indices, balanced count, expected and oe), plus `provenance`.
- **Example code:**

```python
from dnakit.epigenome import observed_over_expected_pixels

result = observed_over_expected_pixels("synthetic.cool", verbose=False)
valid = [row for row in result["pixels"] if row["oe"] is not None]
print(result["pixel_count"], len(valid))
print(round(valid[0]["oe"], 4))
```

- **Example output:**

```text
337 272
0.3723
```

## pairtools functions (4)

## 21) `dnakit.epigenome.sam_to_pairs()`

- **Purpose:** Classify the SAM alignments of each read into Hi-C pairs and write a .pairs file.
- **Method:** Calls pairtools.lib.parse.streaming_classify (the regular parser). Only unique alignments are kept, and a standard .pairs header is written without the SAM text. The output is a new file.
- **API:** `dnakit.epigenome.sam_to_pairs(sam_path, output_path, *, min_mapq=40, max_molecule_size=500, max_inter_align_gap=20, walks_policy='mask', verbose=True)`.
- **Input:** `sam_path`: a SAM or BAM file with an @SQ header; `output_path`: a new .pairs file; `min_mapq`, `max_molecule_size`, `max_inter_align_gap` and `walks_policy`.
- **Output:** A dictionary with `pair_count`, `chromosome_count` and `columns`, plus `provenance`.
- **Example code:**

```python
from dnakit.epigenome import sam_to_pairs

result = sam_to_pairs("reads.sam", "reads.pairs", verbose=False)
print(result["pair_count"], result["chromosome_count"])
print(open("reads.pairs").read())
```

- **Example output:**

```text
2 2
## pairs format v1.0.0
#shape: upper triangle
#genome_assembly: unknown
#chromsize: chr1 1000
#chromsize: chr2 800
#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type
r1	chr1	101	chr2	301	+	+	UU
r2	chr2	51	chr2	201	+	+	UU
```

## 22) `dnakit.epigenome.filter_pairs_by_coverage()`

- **Purpose:** Remove pairs whose fragment ends are covered too often by nearby pairs.
- **Method:** Calls pairtools.lib.filterbycov.streaming_filterbycov: each end counts its partners within max_dist, combined by method (sum or max). Pairs above max_cov are removed and can be written to removed_path.
- **API:** `dnakit.epigenome.filter_pairs_by_coverage(pairs_path, output_path, *, method='max', max_dist=500, max_cov=8, removed_path=None, verbose=True)`.
- **Input:** `pairs_path`: a .pairs file; `output_path`: the kept pairs; `method`, `max_dist` (default 500) and `max_cov` (default 8).
- **Output:** A dictionary with `kept_count` and `removed_count`, plus `provenance`.
- **Example code:**

```python
from dnakit.epigenome import filter_pairs_by_coverage

everything = filter_pairs_by_coverage("cov.pairs", "all.pairs", max_cov=100, verbose=False)
print(everything["kept_count"], everything["removed_count"])
strict = filter_pairs_by_coverage(
    "cov.pairs",
    "strict.pairs",
    method="max",
    max_dist=500,
    max_cov=3,
    removed_path="removed.pairs",
    verbose=False,
)
print(strict["kept_count"], strict["removed_count"])
```

- **Example output:**

```text
5 None
1 4
```

## 23) `dnakit.epigenome.mark_sam_duplicate()`

- **Purpose:** Set the duplicate flag (bit 1024) of one SAM record and change its Yt:Z tag to DD.
- **Method:** Calls pairtools.lib.dedup.mark_sam_as_dup. Tabs are replaced by the pairsam field separator \x19 before the call and restored afterwards.
- **API:** `dnakit.epigenome.mark_sam_duplicate(sam_line, *, verbose=True)`.
- **Input:** `sam_line`: one SAM record with at least 11 columns.
- **Output:** A dictionary with `sam`, `flag_before`, `flag_after` and `duplicate_flag_set`, plus `provenance`.
- **Example code:**

```python
from dnakit.epigenome import mark_sam_duplicate

line = "r1\t99\tchr1\t101\t60\t10M\t=\t201\t110\tACGTACGTAC\tIIIIIIIIII\tYt:Z:UU"
result = mark_sam_duplicate(line, verbose=False)
print(result["flag_before"], result["flag_after"], result["duplicate_flag_set"])
print(result["sam"])
```

- **Example output:**

```text
99 1123 True
r1	1123	chr1	101	60	10M	=	201	110	ACGTACGTAC	IIIIIIIIII	Yt:Z:DD
```

## 24) `dnakit.epigenome.chromsizes_from_sam_header()`

- **Purpose:** Read the sequence names and lengths from the @SQ lines of a SAM or BAM header.
- **Method:** Opens the file with pairtools' pysam reader and calls headerops.get_chromsizes_from_pysam_header; the header order is kept.
- **API:** `dnakit.epigenome.chromsizes_from_sam_header(sam_path, *, verbose=True)`.
- **Input:** `sam_path`: a SAM or BAM file with an @SQ header.
- **Output:** A dictionary with `chromosome_count` and `chromosomes` (name and length), plus `provenance`.
- **Example code:**

```python
from dnakit.epigenome import chromsizes_from_sam_header

result = chromsizes_from_sam_header("reads.sam", verbose=False)
for row in result["chromosomes"]:
    print(row["name"], row["length"])
```

- **Example output:**

```text
chr1 1000
chr2 800
```
