# gimmemotifs, tangermeme and snapATAC2 regulatory functions

This page covers 33 wrapped functions in `dnakit.regulation`: 12 from gimmemotifs (PFM construction and comparison, motif scanning and enrichment statistics, clustering, peak and interval handling, promoter windows), 13 from tangermeme (model-based mutagenesis, variant effects, ablation, motif spacing and design, seqlets, k-mers, locus extraction, background loci, information content and motif placement) and 8 from snapATAC2 (TSS enrichment, FRiP, fragment-length and per-chromosome counts, gene activity, region-to-gene links, motif enrichment and peak merging). Every function reads local files, explicit sequences or a PyTorch model that the caller supplies. None downloads a genome or motif database, and none trains a model. The dependencies are optional: [gimmemotifs](https://github.com/vanheeringen-lab/gimmemotifs), [snapATAC2](https://github.com/scverse/snapatac2) and [tangermeme](https://github.com/jmschrei/tangermeme) (tangermeme also needs PyTorch).

If a package is missing, the function raises an error with the `bio-gimmemotifs`, `bio-tangermeme` or `bio-snapatac2` install hint. Upstream packages are imported only when the function is called.

The examples import from the module files (`dnakit.regulation.<module>`); package-level exports are wired separately.

## Example data

The examples below share these small files. They only demonstrate the interface, and their counts and scores do not describe real data. The example code writes its own genome FASTA.

**Example file:** `fragments.tsv`

```text
chr1	100	300	cellA	1
chr1	120	330	cellA	2
chr1	600	800	cellA	1
chr2	50	250	cellA	1
chr1	650	900	cellB	1
chr2	50	250	cellB	1
chr2	1200	1400	cellB	3
chr1	2000	2200	cellB	1
```

**Example file:** `genes.gtf`

```text
chr1	test	gene	90	400	.	+	.	gene_id "G1"; gene_name "GENEA";
chr1	test	transcript	90	400	.	+	.	gene_id "G1"; transcript_id "T1"; gene_name "GENEA";
chr1	test	gene	1000	1300	.	-	.	gene_id "G2"; gene_name "GENEB";
chr1	test	transcript	1000	1300	.	-	.	gene_id "G2"; transcript_id "T2"; gene_name "GENEB";
chr2	test	gene	400	1600	.	+	.	gene_id "G3"; gene_name "GENEC";
chr2	test	transcript	400	1600	.	+	.	gene_id "G3"; transcript_id "T3"; gene_name "GENEC";
```

**Example file:** `peaks.bed`

```text
chr1	100	400
chr1	1000	1300
chr2	500	900
```

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

**Example file:** `motifs.meme`

```text
MEME version 4

ALPHABET= ACGT

strands: + -

Background letter frequencies
A 0.25 C 0.25 G 0.25 T 0.25

MOTIF GTAT_demo
letter-probability matrix: alength= 4 w= 4 nsites= 20 E= 0
0 0 1 0
0 0 0 1
1 0 0 0
0 0 0 1
```

**Example file:** `foreground.fa`

```text
>f1
ACGTTGATAAGCCGTAGGATAAGTTCAG
>f2
TTGATAAGGGCCATTAGAAGTCCTAGGA
```

**Example file:** `background.fa`

```text
>b1
GGGCCCAAATTTGGGCCCAAATTTGGGC
>b2
ATATATCGCGCGTATATAGCGCGTATAT
```

**Example file:** `peaks.narrowPeak`

```text
chr1	100	400	p1	10	+	200	5.0	3.0	150
chr1	1000	1300	p2	8	+	120	4.0	2.0	100
chr2	500	900	p3	9	+	300	4.5	2.5	200
```

**Example file:** `summits.bed`

```text
chr1	250	251	summit1	9.5
chr2	600	601	summit2	7.0
```

**Example file:** `genes6.bed`

```text
chr1	500	900	G1	0	+
chr1	2000	2500	G2	0	-
chr2	400	800	G3	0	+
```

**Example file:** `loci.bed`

```text
chr1	100	120
chr2	500	520
```

## 1) `dnakit.regulation.motif_from_consensus()`

- **Purpose:** Build a PFM count matrix from an IUPAC consensus string.
- **Method:** Calls gimmemotifs `motif_from_consensus()`; each position totals `n_sites`, and ambiguous codes split the count across their bases.
- **API:** `dnakit.regulation.motif_from_consensus(consensus, *, n_sites=1200, verbose=True)`.
- **Input:** An IUPAC nucleotide consensus string of at most 100 symbols.
- **Output:** `consensus`, `length` and `pfm` (one row per position, columns A, C, G, T).
- **Example code:**

```python
from dnakit.regulation import motif_from_consensus

result = motif_from_consensus("GATAAG", n_sites=100)
print(result["consensus"], result["length"])
print(result["pfm"][0], result["pfm"][2])
```

- **Example output:**

```text
GATAAG 6
[0, 0, 100, 0] [0, 0, 0, 100]
```

## 2) `dnakit.regulation.motif_from_alignment()`

- **Purpose:** Build a PFM and probability matrix from aligned, equal-length sequences.
- **Method:** Calls gimmemotifs `motif_from_align()`; `pfm` holds the base counts of the alignment and `ppm` the probabilities gimmemotifs derives (with a small pseudocount).
- **API:** `dnakit.regulation.motif_from_alignment(sequences, *, verbose=True)`.
- **Input:** A list of two or more equal-length DNA strings.
- **Output:** `consensus`, `length`, `pfm`, `ppm` and `sequence_count`.
- **Example code:**

```python
from dnakit.regulation import motif_from_alignment

result = motif_from_alignment(["ACGT", "ACGA", "ACGT"])
print(result["consensus"], result["length"])
print(result["pfm"][3])
```

- **Example output:**

```text
ACGw 4
[1, 0, 0, 2]
```

## 3) `dnakit.regulation.compare_motif_positions()`

- **Purpose:** Compare two equal-length PFMs position by position.
- **Method:** Both matrices are normalised per row and given gimmemotifs' pseudocount, then each position is scored with `akl` (average Kullback-Leibler similarity), `chisq` or `ssd` (sum of squared distances).
- **API:** `dnakit.regulation.compare_motif_positions(pfm_a, pfm_b, *, metric="akl", verbose=True)`.
- **Input:** Two equal-length PFMs (four non-negative counts per row) and the metric name.
- **Output:** `position_scores`, `mean_score`, `metric` and `length`.
- **Example code:**

```python
from dnakit.regulation import compare_motif_positions

pfm_a = [[10, 0, 0, 0], [0, 10, 0, 0], [0, 0, 10, 0]]
pfm_b = [[10, 0, 0, 0], [0, 0, 10, 0], [0, 0, 0, 10]]
akl = compare_motif_positions(pfm_a, pfm_b)
ssd = compare_motif_positions(pfm_a, pfm_b, metric="ssd")
print([round(v, 3) for v in akl["position_scores"]])
print([round(v, 3) for v in ssd["position_scores"]])
```

- **Example output:**

```text
[10.0, -3.815, -3.815]
[2.0, 0.0, 0.0]
```

## 4) `dnakit.regulation.motif_score_correlation()`

- **Purpose:** Best Pearson correlation between the scores of two PFMs.
- **Method:** Calls gimmemotifs `seqcor()`, which scans both motifs over a 7-mer de Bruijn sequence and compares their score profiles.
- **API:** `dnakit.regulation.motif_score_correlation(pfm_a, pfm_b, *, verbose=True)`.
- **Input:** Two PFMs.
- **Output:** `score` (correlation), `position` (offset of the second motif) and `strand` (+1 or -1).
- **Example code:**

```python
from dnakit.regulation import motif_score_correlation

pfm_a = [[0, 0, 10, 0], [0, 0, 0, 10], [10, 0, 0, 0], [0, 0, 0, 10]]
pfm_b = [[10, 0, 0, 0], [0, 10, 0, 0], [0, 0, 10, 0], [0, 0, 0, 10]]
result = motif_score_correlation(pfm_a, pfm_b)
print(round(result["score"], 4), result["position"], result["strand"])
```

- **Example output:**

```text
0.5012 -2 -1
```

## 5) `dnakit.regulation.scan_best_motif_hits()`

- **Purpose:** Find the best hit of each motif in each FASTA sequence.
- **Method:** Calls gimmemotifs `scan_to_best_match()`, comparing log-odds scores (the sum of `ln(ppm / 0.25 + 0.01)`) on both strands. `start` is a forward-strand coordinate, so minus-strand hits also use the forward start.
- **API:** `dnakit.regulation.scan_best_motif_hits(sequences, motifs, *, max_sequences=1000, verbose=True)`.
- **Input:** A FASTA file and a PFM file.
- **Output:** `hits`: one item per sequence and motif with `sequence`, `motif`, `score`, `start`, `end` and `strand`, in FASTA order then motif ID.
- **Example code:**

```python
from dnakit.regulation import scan_best_motif_hits

result = scan_best_motif_hits("foreground.fa", "motifs.pfm")
print(result["sequence_count"], result["motif_count"], result["hit_count"])
for hit in result["hits"]:
    print(hit["sequence"], hit["motif"], hit["start"], hit["end"], hit["strand"])
```

- **Example output:**

```text
2 2 4
f1 ACGT_demo 0 4 +
f1 GTTT_demo 2 6 +
f2 ACGT_demo 17 21 +
f2 GTTT_demo 2 6 +
```

## 6) `dnakit.regulation.motif_enrichment_stats()`

- **Purpose:** Enrichment statistics for each motif in foreground and background FASTA files.
- **Method:** Calls gimmemotifs `calc_stats()` for ROC AUC, PR AUC, recall and enrichment at a false-positive rate. No z-score or GC correction is applied, so no genome is needed.
- **API:** `dnakit.regulation.motif_enrichment_stats(foreground, background, motifs, *, verbose=True)`.
- **Input:** Foreground FASTA, background FASTA and a PFM file.
- **Output:** `motifs`: one item per motif with `motif` and `stats` (a mapping from statistic name to value), sorted by motif ID.
- **Example code:**

```python
from dnakit.regulation import motif_enrichment_stats

result = motif_enrichment_stats("foreground.fa", "background.fa", "motifs.pfm")
for motif in result["motifs"]:
    print(motif["motif"], round(motif["stats"]["roc_auc"], 4))
```

- **Example output:**

```text
ACGT_demo_ACGT 0.875
GTTT_demo_GTTT 0.25
```

## 7) `dnakit.regulation.cluster_motif_set()`

- **Purpose:** Hierarchically cluster PFM motifs by similarity.
- **Method:** Calls gimmemotifs `cluster_motifs()`, scoring every pair with `metric` and `match`, and merges pairs above `threshold`. Each cluster gives its members and the ID of its averaged representative.
- **API:** `dnakit.regulation.cluster_motif_set(motifs, *, metric="wic", match="total", threshold=0.95, verbose=True)`.
- **Input:** A PFM file with at least two motifs; `metric` is wic, pcc, ed, distance or chisq; `match` is total or subtotal.
- **Output:** `clusters`: each item has `representative` and `members`, sorted by the first member.
- **Example code:**

```python
from dnakit.regulation import cluster_motif_set

result = cluster_motif_set("motifs.pfm", threshold=0.95)
for cluster in result["clusters"]:
    print(cluster["representative"], cluster["members"])
```

- **Example output:**

```text
ACGT_demo ['ACGT_demo']
GTTT_demo ['GTTT_demo']
```

## 8) `dnakit.regulation.read_peak_summits()`

- **Purpose:** Read peaks from a MACS2 summits.bed or narrowPeak file.
- **Method:** Calls gimmemotifs `read_peak_file_to_df()`. Summit rows are 1 bp; narrowPeak `start` is the peak start plus the summit offset.
- **API:** `dnakit.regulation.read_peak_summits(path, *, max_peaks=100000, verbose=True)`.
- **Input:** A summits.bed or narrowPeak file; other BED files are rejected.
- **Output:** `peaks` (columns such as `chrom`, `start`, `end`, `name`, `value`, `log_value`) and `columns`.
- **Example code:**

```python
from dnakit.regulation import read_peak_summits

result = read_peak_summits("summits.bed")
for peak in result["peaks"]:
    print(peak["chrom"], peak["start"], peak["end"], peak["name"], round(peak["log_value"], 4))
```

- **Example output:**

```text
chr1 250 251 summit1 2.3514
chr2 600 601 summit2 2.0794
```

## 9) `dnakit.regulation.summit_intervals_from_narrowpeak()`

- **Purpose:** Fixed-width intervals centred on narrowPeak summits.
- **Method:** Calls gimmemotifs `narrowpeak_to_bed()`; each interval starts at `peak start + summit - size // 2` and is `size` bases long.
- **API:** `dnakit.regulation.summit_intervals_from_narrowpeak(path, *, size, verbose=True)`.
- **Input:** A narrowPeak file (column 10 is the summit offset) and the width `size`.
- **Output:** `intervals`: each item has `chrom`, `start`, `end` and `columns`.
- **Example code:**

```python
from dnakit.regulation import summit_intervals_from_narrowpeak

result = summit_intervals_from_narrowpeak("peaks.narrowPeak", size=50)
for item in result["intervals"]:
    print(item["chrom"], item["start"], item["end"])
```

- **Example output:**

```text
chr1 225 275
chr1 1075 1125
chr2 675 725
```

## 10) `dnakit.regulation.equal_size_intervals()`

- **Purpose:** Resize BED intervals to a fixed width around their midpoints.
- **Method:** Calls gimmemotifs `write_equalsize_bedfile()`. With `keep_columns=False` only the first six BED columns are kept, which matches the output of the upstream `peakfile2bedfile()`.
- **API:** `dnakit.regulation.equal_size_intervals(path, *, size, keep_columns=True, verbose=True)`.
- **Input:** A BED file with at least three columns and the width `size`.
- **Output:** `intervals`: each item has `chrom`, `start`, `end` and `columns` (the remaining columns).
- **Example code:**

```python
from dnakit.regulation import equal_size_intervals

result = equal_size_intervals("peaks.narrowPeak", size=100, keep_columns=False)
for item in result["intervals"]:
    print(item["chrom"], item["start"], item["end"], item["columns"][0])
```

- **Example output:**

```text
chr1 200 300 p1
chr1 1100 1200 p2
chr2 650 750 p3
```

## 11) `dnakit.regulation.random_genomic_regions()`

- **Purpose:** Random fixed-width intervals from a genome FASTA.
- **Method:** Calls gimmemotifs `create_random_genomic_bedfile()`. A fixed `seed` makes the result reproducible; intervals are sorted by chromosome and start.
- **API:** `dnakit.regulation.random_genomic_regions(genome, size, n, *, seed=None, verbose=True)`.
- **Input:** A genome FASTA, the width `size`, the count `n` and an optional `seed`.
- **Output:** `intervals`: each item has `chrom`, `start` and `end`.
- **Example code:**

```python
open("genome.fa", "w").write(">chr1\n" + "ACGTTGATAAGCCGTAGGAT" * 100 + "\n" + ">chr2\n" + "TTGATAAGGGCCATTAGAAG" * 100 + "\n")
from dnakit.regulation import random_genomic_regions

result = random_genomic_regions("genome.fa", 100, 3, seed=1)
for item in result["intervals"]:
    print(item["chrom"], item["start"], item["end"])
```

- **Example output:**

```text
chr1 484 584
chr2 854 954
chr2 941 1041
```

## 12) `dnakit.regulation.promoter_regions()`

- **Purpose:** Sample promoter windows from a six-column gene BED file.
- **Method:** Calls gimmemotifs `create_promoter_bedfile()`: plus-strand genes give the `size` bases before the start, minus-strand genes the bases after the end. The global random state is saved and restored around the call.
- **API:** `dnakit.regulation.promoter_regions(genes, size, n, *, seed=None, verbose=True)`.
- **Input:** A six-column gene BED file (chromosome, start, end, name, score, strand), the window width and the count `n`.
- **Output:** `intervals`: each item has `chrom`, `start`, `end` and `strand`.
- **Example code:**

```python
from dnakit.regulation import promoter_regions

result = promoter_regions("genes6.bed", 100, 3, seed=1)
for item in result["intervals"]:
    print(item["chrom"], item["start"], item["end"], item["strand"])
```

- **Example output:**

```text
chr1 400 500 +
chr1 2500 2600 -
chr2 300 400 +
```

## 13) `dnakit.regulation.saturation_mutagenesis()`

- **Purpose:** Attribute every base in a window by testing all single-base substitutions.
- **Method:** Calls tangermeme `saturation_mutagenesis()`. Each position's mutant-minus-reference output is centred across the four bases; with `hypothetical=False` only the observed base is kept.
- **API:** `dnakit.regulation.saturation_mutagenesis(model, sequences, *, start=0, end=None, target=None, hypothetical=False, batch_size=32, verbose=True)`.
- **Input:** A PyTorch model in eval mode and equal-length DNA sequences (a mapping from name to sequence).
- **Output:** `attributions`, indexed by sequence, base (A, C, G, T) and position over the window `[start, end)`.
- **Example code:**

```python
import torch

class Linear(torch.nn.Module):
    def __init__(self):
        super().__init__()
        self.register_buffer("weights", torch.sin(torch.arange(80, dtype=torch.float32)).reshape(4, 20))

    def forward(self, x, *args):
        return (x * self.weights).sum(dim=(1, 2)).unsqueeze(-1)

model = Linear().eval()
seqs = {"s1": "ACGTTGATAAGCCGTAGGAT", "s2": "TTGATAAGGGCCATTAGAAG"}

from dnakit.regulation import saturation_mutagenesis

result = saturation_mutagenesis(model, seqs, start=0, end=5)
print(result["shape"], result["sequences"])
print([round(v, 4) for v in result["attributions"][0][0]])
```

- **Example output:**

```text
[2, 4, 5] ['s1', 's2']
[-0.3383, 0.0, 0.0, 0.0, -0.0]
```

## 14) `dnakit.regulation.substitution_effect()`

- **Purpose:** Model output change after single-base substitutions.
- **Method:** Calls tangermeme `substitution_effect()`; `delta` is the output after minus the output before. Several substitutions in one sequence are applied together.
- **API:** `dnakit.regulation.substitution_effect(model, sequences, substitutions, *, verbose=True)`.
- **Input:** A model, a sequence mapping and a list of `{"sequence", "position", "base"}` items.
- **Output:** `y_before`, `y_after` and `delta` (one row per sequence).
- **Example code:**

```python
import torch

class Linear(torch.nn.Module):
    def __init__(self):
        super().__init__()
        self.register_buffer("weights", torch.sin(torch.arange(80, dtype=torch.float32)).reshape(4, 20))

    def forward(self, x, *args):
        return (x * self.weights).sum(dim=(1, 2)).unsqueeze(-1)

model = Linear().eval()
seqs = {"s1": "ACGTTGATAAGCCGTAGGAT", "s2": "TTGATAAGGGCCATTAGAAG"}

from dnakit.regulation import substitution_effect

result = substitution_effect(model, seqs, [{"sequence": "s1", "position": 3, "base": "G"}])
print([round(row[0], 4) for row in result["delta"]])
```

- **Example output:**

```text
[-0.9991, 0.0]
```

## 15) `dnakit.regulation.insertion_effect()`

- **Purpose:** Model output change after inserting bases.
- **Method:** Calls tangermeme `insertion_effect()`. The inserted sequence is trimmed back to the model length from the right, or from the left when `left=True`.
- **API:** `dnakit.regulation.insertion_effect(model, sequences, insertions, *, left=False, verbose=True)`.
- **Input:** A model, a sequence mapping and `{"sequence", "position", "base"}` insertions (position is the gap before that original base).
- **Output:** `y_before`, `y_after`, `delta` and `left`.
- **Example code:**

```python
import torch

class Linear(torch.nn.Module):
    def __init__(self):
        super().__init__()
        self.register_buffer("weights", torch.sin(torch.arange(80, dtype=torch.float32)).reshape(4, 20))

    def forward(self, x, *args):
        return (x * self.weights).sum(dim=(1, 2)).unsqueeze(-1)

model = Linear().eval()
seqs = {"s1": "ACGTTGATAAGCCGTAGGAT", "s2": "TTGATAAGGGCCATTAGAAG"}

from dnakit.regulation import insertion_effect

result = insertion_effect(model, seqs, [{"sequence": "s2", "position": 4, "base": "C"}])
print([round(row[0], 4) for row in result["y_after"]])
```

- **Example output:**

```text
[0.7906, -4.9901]
```

## 16) `dnakit.regulation.deletion_effect()`

- **Purpose:** Model output change after deleting bases.
- **Method:** Calls tangermeme `deletion_effect()`. Sequences must be longer than the model input, because each deletion is trimmed back to the model length; `left` chooses the left or right window.
- **API:** `dnakit.regulation.deletion_effect(model, sequences, deletions, *, left=False, verbose=True)`.
- **Input:** Equal-length sequences longer than the model input, and deletions `{"sequence", "position"}`.
- **Output:** `y_before`, `y_after` and `delta`.
- **Example code:**

```python
import torch

class Linear(torch.nn.Module):
    def __init__(self):
        super().__init__()
        self.register_buffer("weights", torch.sin(torch.arange(80, dtype=torch.float32)).reshape(4, 20))

    def forward(self, x, *args):
        return (x * self.weights).sum(dim=(1, 2)).unsqueeze(-1)

model = Linear().eval()
seqs = {"s1": "ACGTTGATAAGCCGTAGGAT", "s2": "TTGATAAGGGCCATTAGAAG"}

from dnakit.regulation import deletion_effect

long_seqs = {"s1": seqs["s1"] + "A", "s2": seqs["s2"] + "A"}
result = deletion_effect(model, long_seqs, [{"sequence": "s1", "position": 5}])
print([round(row[0], 4) for row in result["y_after"]])
```

- **Example output:**

```text
[-1.1958, -8.5005]
```

## 17) `dnakit.regulation.ablate_region()`

- **Purpose:** Shuffle a region repeatedly and measure the model output.
- **Method:** Calls tangermeme `ablate()`. `shuffle` is `dinucleotide` (keeps dinucleotide counts) or `shuffle` (random permutation). A fixed `random_state` makes the shuffles reproducible.
- **API:** `dnakit.regulation.ablate_region(model, sequences, *, start, end, n=20, shuffle="dinucleotide", random_state=0, verbose=True)`.
- **Input:** A model, equal-length sequences, the region `[start, end)` and the shuffle count `n`.
- **Output:** `y_before` (one row per sequence) and `y_after` (sequence, shuffle, output).
- **Example code:**

```python
import torch

class Linear(torch.nn.Module):
    def __init__(self):
        super().__init__()
        self.register_buffer("weights", torch.sin(torch.arange(80, dtype=torch.float32)).reshape(4, 20))

    def forward(self, x, *args):
        return (x * self.weights).sum(dim=(1, 2)).unsqueeze(-1)

model = Linear().eval()
seqs = {"s1": "ACGTTGATAAGCCGTAGGAT", "s2": "TTGATAAGGGCCATTAGAAG"}

from dnakit.regulation import ablate_region

result = ablate_region(model, seqs, start=2, end=14, n=2, random_state=5)
print(len(result["y_after"]), len(result["y_after"][0]))
print([round(v[0], 4) for v in result["y_after"][0]])
```

- **Example output:**

```text
2 2
[1.9214, -1.6687]
```

## 18) `dnakit.regulation.motif_pair_spacing()`

- **Purpose:** Model output for motifs placed at several spacings.
- **Method:** Calls tangermeme `space()`. Each row of `spacings` gives the gaps between adjacent motifs; `start` is the position of the first motif (centred by default).
- **API:** `dnakit.regulation.motif_pair_spacing(model, sequences, motifs, spacings, *, start=None, verbose=True)`.
- **Input:** Motif strings (ACGT), spacing rows with one non-negative gap per adjacent pair, and a model.
- **Output:** `y_before` and `y_after` (spacing row, sequence, output).
- **Example code:**

```python
import torch

class Linear(torch.nn.Module):
    def __init__(self):
        super().__init__()
        self.register_buffer("weights", torch.sin(torch.arange(80, dtype=torch.float32)).reshape(4, 20))

    def forward(self, x, *args):
        return (x * self.weights).sum(dim=(1, 2)).unsqueeze(-1)

model = Linear().eval()
seqs = {"s1": "ACGTTGATAAGCCGTAGGAT", "s2": "TTGATAAGGGCCATTAGAAG"}

from dnakit.regulation import motif_pair_spacing

result = motif_pair_spacing(model, seqs, ["GG", "TT"], [[0], [2]], start=2)
print(result["spacing_count"], len(result["y_after"][0]))
print([round(v[0], 4) for v in result["y_after"][1]])
```

- **Example output:**

```text
2 2
[-7.6877, -10.1997]
```

## 19) `dnakit.regulation.design_greedy_substitution()`

- **Purpose:** Greedily insert motifs into one sequence to move the model output towards a target.
- **Method:** Calls tangermeme `greedy_substitution()`. Each round tries every position and candidate motif (with reverse complements if requested), keeps the lowest-error substitution, and stops at `max_iter` or when the error is below `tol`.
- **API:** `dnakit.regulation.design_greedy_substitution(model, sequence, target, *, motifs=None, reverse_complement=True, max_iter=10, tol=0.001, verbose=True)`.
- **Input:** One DNA sequence, a target output list and an optional motif list.
- **Output:** `designed_sequence`, `input_sequence`, `target` and `length`.
- **Example code:**

```python
import torch

class Linear(torch.nn.Module):
    def __init__(self):
        super().__init__()
        self.register_buffer("weights", torch.sin(torch.arange(80, dtype=torch.float32)).reshape(4, 20))

    def forward(self, x, *args):
        return (x * self.weights).sum(dim=(1, 2)).unsqueeze(-1)

model = Linear().eval()
seqs = {"s1": "ACGTTGATAAGCCGTAGGAT", "s2": "TTGATAAGGGCCATTAGAAG"}

from dnakit.regulation import design_greedy_substitution

result = design_greedy_substitution(model, seqs["s1"], [5.0], motifs=["GGG"], max_iter=2)
print(result["designed_sequence"], result["length"])
```

- **Example output:**

```text
ACGTTGATAAGGGGTAGGGG 20
```

## 20) `dnakit.regulation.call_seqlets()`

- **Purpose:** Call seqlets (contiguous high-attribution regions) from attribution rows.
- **Method:** Calls tangermeme `recursive_seqlets()`, which keeps regions with empirical p-values below `threshold`. Results are sorted by example, start and end.
- **API:** `dnakit.regulation.call_seqlets(attributions, *, threshold=0.01, min_seqlet_len=4, max_seqlet_len=25, additional_flanks=0, verbose=True)`.
- **Input:** A list with one row per example and one value per position.
- **Output:** `seqlets`: each item has `example_index`, `start`, `end`, `attribution` and `p_value`.
- **Example code:**

```python
import numpy as np
from dnakit.regulation import call_seqlets

rng = np.random.RandomState(0)
values = rng.randn(20, 300) * 0.5
values[0, 100:112] = 8.0
result = call_seqlets(values.tolist(), threshold=0.01)
print(result["seqlet_count"])
for item in result["seqlets"][:2]:
    print(item["example_index"], item["start"], item["end"])
```

- **Example output:**

```text
10
0 97 115
4 151 155
```

## 21) `dnakit.regulation.kmer_counts()`

- **Purpose:** Count contiguous k-mers per sequence, or sum per-position scores over them.
- **Method:** Calls tangermeme `kmers()`. With `scores`, each k-mer occurrence adds the sum of the scores it spans; without them each occurrence counts 1. Only k-mers that occur are returned.
- **API:** `dnakit.regulation.kmer_counts(sequences, k, *, scores=None, verbose=True)`.
- **Input:** A mapping from name to sequence, `k` (1 to 8) and optional per-position scores.
- **Output:** `entries`: each item has `sequence`, `kmer` and `value`, sorted by sequence and k-mer.
- **Example code:**

```python
from dnakit.regulation import kmer_counts

counts = kmer_counts({"s": "ACGT"}, 2)
print([(item["kmer"], item["value"]) for item in counts["entries"]])
weighted = kmer_counts({"s": "ACGT"}, 2, scores={"s": [1.0, 2.0, 3.0, 4.0]})
print([item["value"] for item in weighted["entries"]])
```

- **Example output:**

```text
[('AC', 1.0), ('CG', 1.0), ('GT', 1.0)]
[3.0, 5.0, 7.0]
```

## 22) `dnakit.regulation.extract_loci_sequences()`

- **Purpose:** Cut fixed-width sequences around BED loci from a genome FASTA.
- **Method:** Calls tangermeme `extract_loci()`. Each window is centred on the locus midpoint and is `in_window` bases wide.
- **API:** `dnakit.regulation.extract_loci_sequences(loci, genome, *, in_window=200, out_window=1000, max_jitter=0, verbose=True)`.
- **Input:** A BED file of loci (at least three columns) and a genome FASTA.
- **Output:** `loci`: each item has `chrom`, `start`, `end` and `sequence`.
- **Example code:**

```python
open("genome.fa", "w").write(">chr1\n" + "ACGTTGATAAGCCGTAGGAT" * 100 + "\n" + ">chr2\n" + "TTGATAAGGGCCATTAGAAG" * 100 + "\n")

from dnakit.regulation import extract_loci_sequences

result = extract_loci_sequences("loci.bed", "genome.fa", in_window=10, out_window=8)
for item in result["loci"]:
    print(item["chrom"], item["start"], item["end"], item["sequence"])
```

- **Example output:**

```text
chr1 100 120 GATAAGCCGT
chr2 500 520 AAGGGCCATT
```

## 23) `dnakit.regulation.matched_background_loci()`

- **Purpose:** Sample background loci whose GC content matches the input loci.
- **Method:** Calls tangermeme `extract_matching_loci()`: candidate windows are binned by GC content, windows with too many N are dropped, and samples follow the input distribution. A fixed `random_state` makes the sample reproducible.
- **API:** `dnakit.regulation.matched_background_loci(loci, genome, *, in_window=200, out_window=1000, max_n_perc=0.1, gc_bin_width=0.02, random_state=0, verbose=True)`.
- **Input:** A BED file of loci, a genome FASTA and the window width.
- **Output:** `loci` (`chrom`, `start`, `end`), `background_count` and `input_locus_count`.
- **Example code:**

```python
open("genome.fa", "w").write(">chr1\n" + "ACGTTGATAAGCCGTAGGAT" * 100 + "\n" + ">chr2\n" + "TTGATAAGGGCCATTAGAAG" * 100 + "\n")

from dnakit.regulation import matched_background_loci

result = matched_background_loci("loci.bed", "genome.fa", in_window=40, out_window=20)
print(result["input_locus_count"], result["background_count"])
print(all(item["end"] - item["start"] == 40 for item in result["loci"]))
```

- **Example output:**

```text
2 2
True
```

## 24) `dnakit.regulation.pwm_information_content()`

- **Purpose:** Information content (bits) of each PWM position.
- **Method:** Calls tangermeme `information_content()`, which gives `2 - entropy` for each position.
- **API:** `dnakit.regulation.pwm_information_content(pwm, *, verbose=True)`.
- **Input:** A list of rows with four probabilities (A, C, G, T).
- **Output:** `information_content` (one value per position) and `total`.
- **Example code:**

```python
from dnakit.regulation import pwm_information_content

result = pwm_information_content([[0.25, 0.25, 0.25, 0.25], [1.0, 0.0, 0.0, 0.0], [0.5, 0.5, 0.0, 0.0]])
print(result["information_content"], result["total"])
```

- **Example output:**

```text
[0.0, 2.0, 1.0] 3.0
```

## 25) `dnakit.regulation.place_motif()`

- **Purpose:** Write a motif into each sequence by substitution or insertion.
- **Method:** Calls tangermeme `substitute()` or `insert()`. `mode="insert"` trims the result back to the original length from the right; `start=None` centres the motif.
- **API:** `dnakit.regulation.place_motif(sequences, motif, *, start=None, mode="substitute", verbose=True)`.
- **Input:** A mapping from name to sequence, an ACGT motif and the `mode`.
- **Output:** `sequences`: each item has `name` and `sequence`.
- **Example code:**

```python
from dnakit.regulation import place_motif

seqs = {"s1": "ACGTTGATAAGCCGTAGGAT"}
print(place_motif(seqs, "GGG", start=5)["sequences"][0]["sequence"])
print(place_motif(seqs, "GGG", start=5, mode="insert")["sequences"][0]["sequence"])
```

- **Example output:**

```text
ACGTTGGGAAGCCGTAGGAT
ACGTTGGGGATAAGCCGTAG
```

## 26) `dnakit.regulation.tss_enrichment()`

- **Purpose:** Transcription start site enrichment (TSSe) per cell.
- **Method:** Imports the fragments into a temporary snapATAC2 file, calls `metrics.tsse()` with TSSs from the GTF transcripts, and deletes the file.
- **API:** `dnakit.regulation.tss_enrichment(fragment_file, chrom_sizes, gene_gtf, *, exclude_chroms=("chrM", "M"), min_num_fragments=1, max_fragments=5000000, max_cells=100000, verbose=True)`.
- **Input:** A five-column fragment file, chromosome sizes and a GTF gene annotation.
- **Output:** `cells` (`barcode` and `tsse` per cell), `library_tsse` and `frac_overlap_TSS`.
- **Example code:**

```python
from dnakit.regulation import tss_enrichment

result = tss_enrichment("fragments.tsv", {"chr1": 3000, "chr2": 2000}, "genes.gtf")
for cell in result["cells"]:
    print(cell["barcode"], round(cell["tsse"], 4))
print(round(result["library_tsse"], 4))
```

- **Example output:**

```text
cellA 0.0
cellB 0.8658
0.8264
```

## 27) `dnakit.regulation.fragment_in_peaks_fraction()`

- **Purpose:** Fraction of each cell's fragments that fall in peaks (FRiP).
- **Method:** Calls snapATAC2 `metrics.frip()`. With `normalized=True` the count is divided by the cell's total fragments; `count_as_insertion=True` counts Tn5 insertion sites.
- **API:** `dnakit.regulation.fragment_in_peaks_fraction(fragment_file, chrom_sizes, peak_file, *, normalized=True, count_as_insertion=False, min_num_fragments=1, max_fragments=5000000, max_cells=100000, max_peaks=100000, verbose=True)`.
- **Input:** A five-column fragment file, chromosome sizes and a peak BED file (first three columns are coordinates).
- **Output:** `cells`: each item has `barcode` and `frip`.
- **Example code:**

```python
from dnakit.regulation import fragment_in_peaks_fraction

result = fragment_in_peaks_fraction("fragments.tsv", {"chr1": 3000, "chr2": 2000}, "peaks.bed")
for cell in result["cells"]:
    print(cell["barcode"], round(cell["frip"], 4))
```

- **Example output:**

```text
cellA 0.5
cellB 0.0
```

## 28) `dnakit.regulation.fragment_size_distribution()`

- **Purpose:** Dataset-wide fragment length distribution.
- **Method:** Calls snapATAC2 `metrics.frag_size_distr()`; `counts[i]` is the number of fragments of length i bp. Longer fragments are not recorded.
- **API:** `dnakit.regulation.fragment_size_distribution(fragment_file, chrom_sizes, *, max_recorded_size=1000, min_num_fragments=1, max_fragments=5000000, max_cells=100000, verbose=True)`.
- **Input:** A five-column fragment file and chromosome sizes.
- **Output:** `counts`, `fragment_count` and `max_recorded_size`.
- **Example code:**

```python
from dnakit.regulation import fragment_size_distribution

result = fragment_size_distribution("fragments.tsv", {"chr1": 3000, "chr2": 2000}, max_recorded_size=300)
print(result["fragment_count"], result["counts"][200], result["counts"][210])
```

- **Example output:**

```text
8 6 1
```

## 29) `dnakit.regulation.chromosome_fragment_counts()`

- **Purpose:** Per-cell fragment summary for each chromosome.
- **Method:** Calls snapATAC2 `metrics.summary_by_chrom()`; `mode` is `count` (fragments), `sum` (sum of read counts) or `mean`.
- **API:** `dnakit.regulation.chromosome_fragment_counts(fragment_file, chrom_sizes, *, mode="count", min_num_fragments=1, max_fragments=5000000, max_cells=100000, verbose=True)`.
- **Input:** A five-column fragment file, chromosome sizes and the summary mode.
- **Output:** `chromosomes` (sorted by name) and `cells` (`counts` per cell, in the chromosome order).
- **Example code:**

```python
from dnakit.regulation import chromosome_fragment_counts

result = chromosome_fragment_counts("fragments.tsv", {"chr1": 3000, "chr2": 2000})
print(result["chromosomes"])
for cell in result["cells"]:
    print(cell["barcode"], cell["counts"])
```

- **Example output:**

```text
['chr1', 'chr2']
cellA [3.0, 1.0]
cellB [2.0, 2.0]
```

## 30) `dnakit.regulation.gene_activity_counts()`

- **Purpose:** Cell-by-gene fragment counts over gene windows (gene activity matrix).
- **Method:** Calls snapATAC2 `pp.make_gene_matrix()`. Each gene's window covers `upstream` bases before the start, the gene body (optional) and `downstream` bases after the end, with genes taken from the GTF.
- **API:** `dnakit.regulation.gene_activity_counts(fragment_file, chrom_sizes, gene_gtf, *, upstream=2000, downstream=0, include_gene_body=True, min_num_fragments=1, max_fragments=5000000, max_cells=100000, max_entries=1000000, verbose=True)`.
- **Input:** A five-column fragment file, chromosome sizes, a GTF annotation and the window settings.
- **Output:** `genes`, `entries` (`[cell index, gene index, count]`), `nonzero_count` and `total_count`.
- **Example code:**

```python
from dnakit.regulation import gene_activity_counts

result = gene_activity_counts("fragments.tsv", {"chr1": 3000, "chr2": 2000}, "genes.gtf", upstream=500)
print(result["genes"], result["total_count"])
print(result["entries"][:3])
```

- **Example output:**

```text
['GENEA', 'GENEB', 'GENEC'] 5
[[0, 0, 2], [0, 2, 1], [1, 2, 2]]
```

## 31) `dnakit.regulation.regulatory_regions_to_genes()`

- **Purpose:** Link regions to genes within a distance window.
- **Method:** Calls snapATAC2 `tl.init_network_from_annotation()`. A region is linked when it lies within `upstream` bases before a gene's start or `downstream` bases after its end; distances come from snapATAC2.
- **API:** `dnakit.regulation.regulatory_regions_to_genes(regions, gene_gtf, *, upstream=250000, downstream=250000, id_type="gene_name", coding_gene_only=True, verbose=True)`.
- **Input:** Regions as `chr:start-end` strings and a GTF annotation.
- **Output:** `edges`: each item has `region`, `gene` and `distance`, sorted by region then gene.
- **Example code:**

```python
from dnakit.regulation import regulatory_regions_to_genes

result = regulatory_regions_to_genes(["chr1:100-300", "chr1:1000-1200", "chr2:500-700"], "genes.gtf", upstream=250, downstream=250)
for edge in result["edges"]:
    print(edge["region"], edge["gene"], edge["distance"])
```

- **Example output:**

```text
chr1:100-300 GENEA 11
chr1:1000-1200 GENEB 99
chr2:500-700 GENEC 101
```

## 32) `dnakit.regulation.region_motif_enrichment()`

- **Purpose:** Test motif enrichment in region groups against background.
- **Method:** Calls snapATAC2 `tl.motif_enrichment()` with sequences from the genome FASTA; `method` is `binomial`, `hypergeometric` or None (snapATAC2 default).
- **API:** `dnakit.regulation.region_motif_enrichment(regions, motifs, genome, *, background=None, method=None, verbose=True)`.
- **Input:** Region groups (a mapping from group name to regions), a MEME motif file, a genome FASTA and optional background regions.
- **Output:** `groups`: each group has `group` and `rows` (log2 fold change, p-value and adjusted p-value per motif).
- **Example code:**

```python
open("genome.fa", "w").write(">chr1\n" + "ACGTTGATAAGCCGTAGGAT" * 100 + "\n" + ">chr2\n" + "TTGATAAGGGCCATTAGAAG" * 100 + "\n")

from dnakit.regulation import region_motif_enrichment

result = region_motif_enrichment({"g1": ["chr1:100-300", "chr1:1000-1200"], "g2": ["chr2:500-700"]}, "motifs.meme", "genome.fa")
for group in result["groups"]:
    print(group["group"], group["rows"][0]["id"], round(group["rows"][0]["p-value"], 4))
```

- **Example output:**

```text
g1 GTAT_demo 0.3333
g2 GTAT_demo 0.3333
```

## 33) `dnakit.regulation.merge_peak_sets()`

- **Purpose:** Merge group-specific peaks into one non-overlapping set.
- **Method:** Calls snapATAC2 `tl.merge_peaks()`. Each peak is centred on its summit and extended `half_width` bases on both sides; overlaps are resolved by snapATAC2.
- **API:** `dnakit.regulation.merge_peak_sets(peaks_by_group, chrom_sizes, *, half_width=250, verbose=True)`.
- **Input:** MACS2-style peak rows per group (`chrom`, `start`, `end`, `name`, `score`, `strand`, `signal_value`, `p_value`, `q_value`, `summit`) and chromosome sizes.
- **Output:** `peaks`: each item has `chrom`, `start`, `end` and `groups` (the groups that contain the peak).
- **Example code:**

```python
from dnakit.regulation import merge_peak_sets

peak = {"chrom": "chr1", "start": 100, "end": 400, "name": "p", "score": 10, "strand": ".", "signal_value": 5.0, "p_value": 3.0, "q_value": 2.0, "summit": 200}
result = merge_peak_sets({"A": [peak], "B": [dict(peak, chrom="chr2", start=60, end=300, summit=100)]}, {"chr1": 3000, "chr2": 2000}, half_width=50)
for item in result["peaks"]:
    print(item["chrom"], item["start"], item["end"], item["groups"])
```

- **Example output:**

```text
chr1 250 351 ['A']
chr2 110 211 ['B']
```

## How these functions are verified

The test file is `others/tests/unit/test_full_G12.py`, with 37 tests and at least one test per wrapped function. The main checks are:

- gimmemotifs: PFM counts from a consensus and from an alignment match hand counts; scan coordinates on the forward and reverse strands match sites placed by hand; `calc_stats()`, `scan_to_best_match()`, `cluster_motifs()` and `seqcor()` match direct upstream calls; promoter window coordinates match hand calculations, and the global random state is restored after the call.
- tangermeme: the substitution output change equals the difference of two weights; the saturation attribution equals the observed base's weight minus the mean over the four bases; insertion, deletion, ablation, spacing and design results match direct upstream calls; k-mer counts and weighted sums match hand values; locus sequences match slices of the genome; information content and motif placement match hand values.
- snapATAC2: TSS enrichment, FRiP, fragment lengths and chromosome counts match hand values or direct upstream calls; the gene activity matrix matches a direct `make_gene_matrix()` call; region-to-gene distances match a direct call; motif enrichment p-values match a direct call; merged peak coordinates match hand values.
- When any of the three packages is missing, the error message includes the matching install hint (simulated with monkeypatch).

## Current scope

The wrapped operations are listed in the 33 sections above. The remaining functions are recorded in the module ledger, grouped as follows:

- Already covered by existing wrappers: motif reading in PFM, TRANSFAC and JASPAR formats (`read_motif_file`), snapATAC2 `import_fragments`, `add_tile_matrix` and `make_peak_matrix`, and tangermeme `predict`, `marginalize` and `deep_lift_shap`.
- Same as existing dnakit functions: `rc`, `as_fasta`, one-hot encoding, reverse complement, VCF reading, dinucleotide shuffling and GC content.
- Not wrapped because they need an external program (de novo motif tools, gffread), network access (JASPAR, snapATAC2 datasets and genome downloads) or extra data files (GC indexes, ortholog databases, genome-based backgrounds). Gimmemotifs functions that need background data are marked with a `needs-data:` reason in the ledger.
- Plotting and report functions are not wrapped.
- DNA-related functions not selected in this pass (for example tangermeme `gapped_kmers`, `tfmodisco_seqlets` and `annotate_seqlets`, and snapATAC2 `add_tf_binding`, `macs3`, `export_fragments` and `export_coverage`) are recorded as `other`, with the reason in the ledger.
