# Optional regulatory-genomics backends: motif scans, fragment counts and model attributions

This group of functions in `dnakit.regulation` is a set of adapters over optional packages: [gimmemotifs](https://github.com/vanheeringen-lab/gimmemotifs), [snapATAC2](https://github.com/scverse/snapatac2) and [tangermeme](https://github.com/jmschrei/tangermeme). gimmemotifs reads PFM motifs, scans given sequences for sites and scores k-mers. snapATAC2 computes single-cell QC metrics, bin count matrices and peak count matrices from a five-column fragment file. tangermeme calls a PyTorch model that the user supplies, for predictions, motif marginalization and DeepLIFT/SHAP attributions. Every function works only on an explicit file, sequence set or model. None of them downloads a genome, background sequences or motif database, trains a model, or computes the significance of motif enrichment.

The dependencies are optional and are imported only when the matching function is called. If a package is missing, the function raises an error that includes the install command.

## Example data

The examples below use these three small files. They only demonstrate the interface, and their counts and scores do not describe real data.

**Example file:** `motifs.pfm`

```text
>GTTT_demo
0 0 10 0
0 0 0 10
0 0 0 10
0 0 0 10
>ACGT_demo
10 0 0 0
0 10 0 0
0 0 10 0
0 0 0 10
```

**Example file:** `fragments.tsv`

```text
chr1	100	300	cellA	1
chr1	120	330	cellA	2
chr1	600	800	cellA	1
chrM	10	200	cellA	1
chr1	650	900	cellB	1
chr2	50	250	cellB	1
chr2	1200	1400	cellB	3
```

**Example file:** `peaks.bed`

```text
chr1	50	400
chr1	590	950
chr2	1100	1500
```

## 1) `dnakit.regulation.read_motif_file()`

- **Purpose:** Read a PFM, TRANSFAC or JASPAR motif file and return each motif's ID, length, consensus, information content, score range and PFM counts.
- **Method:** Parse the file with gimmemotifs `read_motifs()`. The returned `pfm` keeps the counts stored in the file, so it can be passed directly to `scan_motif_hits()`. The file size and the number of motifs have explicit limits.
- **API:** `dnakit.regulation.read_motif_file(path, *, fmt="pfm", max_bytes=10000000, max_motifs=1000)`.
- **Input:** A text motif file path. `fmt` is `"pfm"`, `"transfac"` or `"jaspar"`.
- **Output:** A dictionary with `motif_count` and `motifs`. Each motif has `id`, `length`, `consensus`, `information_content`, `min_score`, `max_score` and `pfm`.
- **Example code:**

```python
from dnakit.regulation import read_motif_file

result = read_motif_file("motifs.pfm")
for motif in result["motifs"]:
    print(motif["id"], motif["consensus"], motif["length"])
```

- **Example output:**

```text
GTTT_demo GTTT 4
ACGT_demo ACGT 4
```

## 2) `dnakit.regulation.scan_motif_hits()`

- **Purpose:** Find sites of PFM motifs in given DNA sequences, reporting hits on both the forward and reverse strands.
- **Method:** Compute log-odds scores with gimmemotifs `Motif.scan_all()`: for each window, the sum of `ln(ppm / 0.25 + 0.01)` over motif positions. The threshold is a fraction of each motif's own score range: `cutoff = min + (max - min) * threshold`. Reverse-complement scanning is on by default. `max_hits_per_sequence` limits the hits of each motif in each sequence; going over the limit raises an error instead of truncating.
- **API:** `dnakit.regulation.scan_motif_hits(sequences, motifs, *, threshold=0.9, scan_rc=True, max_hits_per_sequence=100, max_sequences=1000, max_sequence_length=100000)`.
- **Input:** `DNA(...)` or `{sequence name: sequence}`. Sequences may contain only A, C, G and T. `motifs` is `{motif ID: PFM count matrix}`, with one row per position and columns A, C, G, T.
- **Output:** `hit_count` and `hits`. Each hit has `sequence`, `motif`, `start` (0-based forward coordinate), `end` (exclusive), `strand` and `score`.
- **Example code:**

```python
from dnakit.regulation import read_motif_file, scan_motif_hits

motifs = {motif["id"]: motif["pfm"] for motif in read_motif_file("motifs.pfm")["motifs"]}
result = scan_motif_hits(
    {"promoter": "GGAAACGGTTTAAC", "control": "CCCCCCCC"},
    motifs,
    threshold=0.9,
)
for hit in result["hits"]:
    print(hit["sequence"], hit["motif"], hit["start"], hit["strand"], round(hit["score"], 3))
```

- **Example output:**

```text
promoter GTTT_demo 2 - 5.554
promoter GTTT_demo 7 + 5.554
```

## 3) `dnakit.regulation.score_motif_kmers()`

- **Purpose:** Compute the log-odds score of equal-length k-mers against one PFM.
- **Method:** Call gimmemotifs `Motif.score_kmer()`. The function also returns the motif's maximum and minimum scores, which help convert scores into thresholds.
- **API:** `dnakit.regulation.score_motif_kmers(motif, kmers, *, max_kmers=10000)`.
- **Input:** One PFM count matrix (rows are positions, columns are A, C, G, T) and a list of k-mer strings. Each k-mer must have the motif length and contain only A, C, G and T.
- **Output:** `consensus`, `max_score`, `min_score` and `scores`. Each score has `kmer` and `score`.
- **Example code:**

```python
from dnakit.regulation import score_motif_kmers

result = score_motif_kmers(
    [[0, 0, 10, 0], [0, 0, 0, 10], [0, 0, 0, 10], [0, 0, 0, 10]],
    ["GTTT", "GTTA"],
)
print(result["consensus"], round(result["max_score"], 3))
for item in result["scores"]:
    print(item["kmer"], round(item["score"], 3))
```

- **Example output:**

```text
GTTT 5.554
GTTT 5.554
GTTA -0.4
```

## 4) `dnakit.regulation.read_fragment_qc()`

- **Purpose:** Compute per-barcode fragment count, duplicate fraction and mitochondrial fraction.
- **Method:** Read a plain-text, headerless, five-column paired-end fragment file (chromosome, start, end, barcode, read count). Each line is validated, then the file is imported with snapATAC2 `pp.import_fragments()`. The per-cell `n_fragment`, `frac_dup` and `frac_mito` values come from its `obs` table. `frac_dup` is duplicate reads divided by all reads; `frac_mito` is fragments on chrM/M divided by all fragments.
- **API:** `dnakit.regulation.read_fragment_qc(fragment_file, chrom_sizes, *, min_num_fragments=1, max_fragments=5000000, max_cells=100000)`.
- **Input:** A fragment file. `chrom_sizes` is `{chromosome: length}`, and every chromosome in the file must appear in it. `min_num_fragments` is passed to snapATAC2 and removes barcodes with fewer fragments.
- **Output:** `cell_count` and `cells`. Each cell has `barcode`, `n_fragment`, `frac_dup` and `frac_mito`.
- **Example code:**

```python
from dnakit.regulation import read_fragment_qc

result = read_fragment_qc(
    "fragments.tsv",
    {"chr1": 2000, "chr2": 2000, "chrM": 1000},
    min_num_fragments=1,
)
for cell in result["cells"]:
    print(cell["barcode"], cell["n_fragment"], round(cell["frac_dup"], 3), cell["frac_mito"])
```

- **Example output:**

```text
cellA 3 0.2 0.25
cellB 3 0.4 0.0
```

## 5) `dnakit.regulation.count_tile_matrix()`

- **Purpose:** Build a cell by fixed-width bin count matrix.
- **Method:** Call snapATAC2 `pp.add_tile_matrix()`. Bins are named `chrom:start-end` (0-based, half-open). chrM, chrY, M and Y are excluded by default. `counting_strategy` is `"fragment"` (each fragment counted once), `"insertion"` or `"paired-insertion"`.
- **API:** `dnakit.regulation.count_tile_matrix(fragment_file, chrom_sizes, *, bin_size=500, counting_strategy="fragment", exclude_chroms=("chrM", "chrY", "M", "Y"), min_num_fragments=1, max_fragments=5000000, max_cells=100000, max_entries=1000000)`.
- **Input:** The same fragment file and `chrom_sizes` as `read_fragment_qc()`. `bin_size` is a positive integer in base pairs.
- **Output:** `shape`, `barcodes`, `features`, `nonzero_count`, `total_count` and `entries`. `entries` are `[cell index, bin index, count]` triplets, sorted by row and column, with indices starting at 0.
- **Example code:**

```python
from dnakit.regulation import count_tile_matrix

result = count_tile_matrix(
    "fragments.tsv",
    {"chr1": 2000, "chr2": 2000, "chrM": 1000},
    bin_size=500,
    min_num_fragments=1,
)
print(result["shape"], result["total_count"])
print(result["features"][:4])
print(result["entries"])
```

- **Example output:**

```text
[2, 8] 6
['chr1:0-500', 'chr1:500-1000', 'chr1:1000-1500', 'chr1:1500-2000']
[[0, 0, 2], [0, 1, 1], [1, 1, 1], [1, 4, 1], [1, 6, 1]]
```

## 6) `dnakit.regulation.count_peak_fragments()`

- **Purpose:** Build a cell by peak count matrix.
- **Method:** Call snapATAC2 `pp.make_peak_matrix()`. Peaks come from the first three columns of a BED file (chromosome, start, end); extra columns are ignored. Peak names are `chrom:start-end`.
- **API:** `dnakit.regulation.count_peak_fragments(fragment_file, chrom_sizes, peak_file, *, counting_strategy="fragment", min_num_fragments=1, max_fragments=5000000, max_cells=100000, max_peaks=100000, max_entries=1000000)`.
- **Input:** The fragment file, `chrom_sizes` and a BED peak file. Chromosomes in the BED file must also appear in `chrom_sizes`.
- **Output:** The same fields as `count_tile_matrix()`, with peak names in `features`.
- **Example code:**

```python
from dnakit.regulation import count_peak_fragments

result = count_peak_fragments(
    "fragments.tsv",
    {"chr1": 2000, "chr2": 2000, "chrM": 1000},
    "peaks.bed",
    min_num_fragments=1,
)
print(result["features"])
print(result["entries"])
```

- **Example output:**

```text
['chr1:50-400', 'chr1:590-950', 'chr2:1100-1500']
[[0, 0, 2], [0, 1, 1], [1, 1, 1], [1, 2, 1]]
```

## 7) `dnakit.regulation.predict_sequences()`

- **Purpose:** Run batched predictions with tangermeme `predict()` on a set of equal-length DNA sequences.
- **Method:** Convert the sequences to one-hot arrays (in A, C, G, T order) and call the model in batches of `batch_size`. The model must be a `torch.nn.Module` in eval mode; otherwise the function raises an error, because dropout and batch norm would change the outputs. Only a single output tensor is supported.
- **API:** `dnakit.regulation.predict_sequences(model, sequences, *, batch_size=32, max_sequences=256, max_sequence_length=10000)`.
- **Input:** A PyTorch model and sequences in the same format as `scan_motif_hits()`. All sequences must have the same length and contain only A, C, G and T.
- **Output:** `sequence_count`, `shape` and `predictions` (a nested list whose first dimension is the sequence).
- **Example code:**

```python
import torch
from dnakit.regulation import predict_sequences

torch.manual_seed(0)
model = torch.nn.Sequential(torch.nn.Flatten(), torch.nn.Linear(4 * 6, 1)).eval()
result = predict_sequences(model, {"a": "ACGTAC", "b": "TTTTGG"})
print(result["shape"])
print([round(row[0], 4) for row in result["predictions"]])
```

- **Example output:**

```text
[2, 1]
[0.094, 0.1866]
```

## 8) `dnakit.regulation.marginalize_motif_effect()`

- **Purpose:** Marginalization analysis: compare model outputs before and after the same motif is written into each background sequence.
- **Method:** Call tangermeme `marginalize()`. The motif overwrites each sequence; when `start` is None it is placed at the centre, at position `length // 2 - motif length // 2`. `delta` is the output after writing minus the output before.
- **API:** `dnakit.regulation.marginalize_motif_effect(model, sequences, motif, *, start=None, max_sequences=256, max_sequence_length=10000)`.
- **Input:** A model, sequences and one ACGT motif string. `start` is the 0-based position of the motif in each sequence and must fit inside the sequence.
- **Output:** `motif`, `start`, `shape`, `y_before`, `y_after` and `delta`.
- **Example code:**

```python
import torch
from dnakit.regulation import marginalize_motif_effect


class CountA(torch.nn.Module):
    """Outputs the number of A bases in each sequence (one-hot channel 0)."""

    def forward(self, inputs):
        return inputs[:, 0, :].sum(dim=1, keepdim=True)


model = CountA().eval()
result = marginalize_motif_effect(model, {"c": "CCCCCCCCCC", "ac": "ACCCCCCCCC"}, "AAAA")
print(result["start"], result["y_before"], result["y_after"], result["delta"])
```

- **Example output:**

```text
3 [[0.0], [1.0]] [[4.0], [5.0]] [[4.0], [4.0]]
```

## 9) `dnakit.regulation.deep_lift_shap_attribution()`

- **Purpose:** Compute per-position, per-base DeepLIFT/SHAP attributions for one model output.
- **Method:** Call tangermeme `deep_lift_shap()`. The references are dinucleotide-shuffled copies of each input sequence, drawn with a fixed `random_state` so results are reproducible. `max_positions` limits the number of sequences times their length.
- **API:** `dnakit.regulation.deep_lift_shap_attribution(model, sequences, *, target=0, n_shuffles=20, random_state=0, max_sequences=256, max_sequence_length=10000, max_positions=200000)`.
- **Input:** A model in eval mode with a single output tensor, sequences, and the output column `target` to explain.
- **Output:** `shape` and `attributions`, indexed by sequence, position and base (A, C, G, T).
- **Example code:**

```python
import torch
from dnakit.regulation import deep_lift_shap_attribution


class CountA(torch.nn.Module):
    """Outputs the number of A bases in each sequence (one-hot channel 0)."""

    def forward(self, inputs):
        return inputs[:, 0, :].sum(dim=1, keepdim=True)


model = CountA().eval()
result = deep_lift_shap_attribution(model, {"s": "ACGTTGCAAC"}, n_shuffles=4, random_state=3)
print(result["shape"])
print([round(value, 2) for value in result["attributions"][0][0]])
```

- **Example output:**

```text
[1, 10, 4]
[0.0, -0.0, -0.0, -0.0]
```

## Installing the optional dependencies

```bash
python -m pip install "dnakit[bio-gimmemotifs]"     # gimmemotifs: read_motif_file / scan_motif_hits / score_motif_kmers
python -m pip install "dnakit[bio-snapatac2]"       # snapATAC2: read_fragment_qc / count_tile_matrix / count_peak_fragments
python -m pip install "dnakit[bio-tangermeme]"      # tangermeme and PyTorch: predict_sequences / marginalize_motif_effect / deep_lift_shap_attribution
```

Project pages: [gimmemotifs](https://github.com/vanheeringen-lab/gimmemotifs), [snapATAC2](https://github.com/scverse/snapatac2), [tangermeme](https://github.com/jmschrei/tangermeme).

## How these functions are verified

The test file is `others/tests/unit/test_regulation_backends.py`, with 20 tests. They check:

- gimmemotifs: PFM, TRANSFAC and JASPAR files give the same counts; reverse-complement hit coordinates match a hand calculation; log-odds scores match a hand sum of `ln(ppm / 0.25 + 0.01)` and a direct `Motif.score_kmer()` call; going over `max_hits_per_sequence` raises an error.
- snapatac2: fragment QC (fragments per cell, duplicate fraction, mitochondrial fraction) matches hand values and a direct `pp.import_fragments()` call; the bin count matrix matches hand counts and a direct `add_tile_matrix()` call; the peak count matrix matches hand counts.
- tangermeme: predictions match a direct model call; the before and after outputs of motif marginalization match hand counts; DeepLIFT attributions match a direct `deep_lift_shap()` call and repeat exactly for the same seed.
- Invalid files, separators, chromosome names, coordinates and row limits are rejected. A missing package raises an install hint. A `BaseException` raised upstream (for example a pyo3 panic) becomes a `BackendExecutionError`; the test uses a stand-in exception to check this.

## Current scope

This batch wraps DNA-level operations from the three packages. Wrapped: gimmemotifs reading of PFM, TRANSFAC and JASPAR files, sequence scanning with `scan_all()`, and k-mer scoring; snapATAC2 fragment QC, bin count matrices and peak count matrices; tangermeme batched prediction, motif marginalization and DeepLIFT/SHAP attributions.

Not wrapped, with reasons:

- gimmemotifs `Scanner` and FPR-based thresholds: they need a genome or background sequence file. This batch only scans offline.
- gimmemotifs motif discovery (MEME, HOMER and similar): these rely on external programs. The motif writers are not wrapped either.
- snapATAC2 single-end fragments (`is_paired=False`): they need an extra strand column. This batch accepts only the five-column paired-end format. gzip input, gene activity matrices, TSS enrichment (which needs gene annotation), peak calling, dimensionality reduction and clustering are also not wrapped.
- tangermeme saturation mutagenesis, variant effects, ablation, sequence design and readers such as `read_meme`: these can be added in later batches as needed.
