# G18i: gimmemotifs, tangermeme and snapATAC2 regulatory review functions

This page covers the 36 wrapped functions that review task G18i adds to `dnakit.regulation`: 9 from gimmemotifs (summit merging, peak coverage, motif scanning to a BED file, region scoring, motif activity, ROC and KS statistics, PFM alignment, FASTA record counts), 19 from tangermeme (sequence editing and shuffling, random sequences, PWM consensus, coordinate conversion, FASTA output, gapped k-mers, annotation counts, model-based ablation and marginalisation, sequence design, PISA attribution and seqlets) and 8 from snapATAC2 (fragment export, coverage tracks, BAM to fragments, contact and value import, 10x QC metrics, MACS3 peaks and TF binding edges).

Every function reads local files, explicit sequences, matrices or a PyTorch model that the caller supplies. None downloads a genome or database. The dependencies are optional: when a package is missing, the function raises an error with the `bio-gimmemotifs`, `bio-tangermeme` or `bio-snapatac2` install hint. Upstream packages are imported only when the function is called.

The examples import from the module files (`dnakit.regulation.<module>`); package-level exports are wired separately.


## 1) `dnakit.regulation.combine_summit_files()`

**Example file:** `gm_summits1.bed`

```text
chr1	100	101	p1	9.0
chr2	500	501	p2	7.0
```

**Example file:** `gm_summits2.bed`

```text
chr1	110	111	p3	8.0
chr1	2000	2001	p4	6.0
```

- **Purpose:** Merge several MACS2-style summit files and keep the strongest summit in each window.
- **Method:** Calls gimmemotifs `combine_peaks()` with a local genome FASTA as the chromosome source; `window` sets the merge width and `scale_value` rescales values per file before taking the maximum.
- **API:** `dnakit.regulation.combine_summit_files(summit_files, genome_fasta, *, window=200, scale_value=False, verbose=True)`.
- **Input:** A list of 1 to 50 summit file paths, a FASTA genome path, the window width (1 to 100000) and the scaling flag.
- **Output:** `summit_count` and `summits` (sorted `chrom`, `start` and `end` records).
- **Example code:**

```python
from pathlib import Path
Path("gm_genome.fa").write_text(
    ">chr1\n" + "ACGTTGATAAGCCGTAGGATAAGTTCAGTTGATAAGGGCCATTAGAAGTCCTAGGA" * 70 + "\n"
    + ">chr2\n" + "GGGCCCAAATTTGGGCCCAAATTTGATAAGGGCCCAAATTTGGGCCCAAATTTGGGC" * 60 + "\n"
)
from dnakit.regulation import combine_summit_files

result = combine_summit_files(["gm_summits1.bed", "gm_summits2.bed"], "gm_genome.fa", window=200)
print(result["summit_count"], result["window"])
for item in result["summits"]:
    print(item["chrom"], item["start"], item["end"])
```

- **Example output:**

```text
3 200
chr1 100 101
chr1 2000 2001
chr2 500 501
```

## 2) `dnakit.regulation.peak_coverage_table()`

**Example file:** `gm_summits1.bed`

```text
chr1	100	101	p1	9.0
chr2	500	501	p2	7.0
```

- **Purpose:** Read the mean bigWig coverage in a window around each peak summit.
- **Method:** Calls gimmemotifs `coverage_table()`; each peak is widened to a `window`-base interval centred on the summit, `log_transform` applies log(1 + x), and a positive `top` keeps the highest-variance peaks.
- **API:** `dnakit.regulation.peak_coverage_table(peak_file, coverage_files, *, window=200, log_transform=True, top=-1, verbose=True)`.
- **Input:** A peak file path, 1 to 50 local bigWig paths, the window width, the log-transform flag and `top` (-1 keeps all peaks).
- **Output:** `tracks` (track names), `peaks` (interval labels) and `values` (one row per peak, one column per track).
- **Example code:**

```python
import pyBigWig
bw = pyBigWig.open("gm_cov.bw", "w")
bw.addHeader([("chr1", 4000), ("chr2", 3000)])
bw.addEntries(["chr1"], [0], ends=[4000], values=[2.0])
bw.close()
from dnakit.regulation import peak_coverage_table

result = peak_coverage_table("gm_summits1.bed", ["gm_cov.bw"], window=200)
print(result["tracks"], result["peaks"])
print([round(v, 4) for v in result["values"][0]], [round(v, 4) for v in result["values"][1]])
```

- **Example output:**

```text
['gm_cov'] ['chr1:0-200', 'chr2:400-600']
[1.0986] [0.0]
```

## 3) `dnakit.regulation.scan_motifs_to_file()`

**Example file:** `gm_motifs.pfm`

```text
>GATA
0 0 10 0
0 0 0 10
10 0 0 0
0 10 0 0
```

- **Purpose:** Write motif hits above a threshold to a BED file and return the sorted hit list.
- **Method:** Calls gimmemotifs `scan_to_file()` with log-odds scores, without z-normalisation or GC correction; `cutoff` is a fraction from 0 to 1 of the threshold, and 0 keeps every hit. The output file must not exist.
- **API:** `dnakit.regulation.scan_motifs_to_file(sequences_fasta, motifs_pfm, output_bed, *, cutoff=0.0, verbose=True)`.
- **Input:** A FASTA sequence file, a PFM file, an output BED path that does not exist yet, and the threshold fraction.
- **Output:** `hit_count` and `hits` (each hit has `chrom`, `start`, `end`, `motif`, `score` and `strand`).
- **Example code:**

```python
from pathlib import Path
Path("gm_genome.fa").write_text(
    ">chr1\n" + "ACGTTGATAAGCCGTAGGATAAGTTCAGTTGATAAGGGCCATTAGAAGTCCTAGGA" * 70 + "\n"
    + ">chr2\n" + "GGGCCCAAATTTGGGCCCAAATTTGATAAGGGCCCAAATTTGGGCCCAAATTTGGGC" * 60 + "\n"
)
from pathlib import Path
from dnakit.regulation import scan_motifs_to_file

Path("gm_hits.bed").unlink(missing_ok=True)
result = scan_motifs_to_file("gm_genome.fa", "gm_motifs.pfm", "gm_hits.bed", cutoff=0.0)
print(result["hit_count"] > 0, result["output"])
hit = result["hits"][0]
print(hit["chrom"], hit["start"], hit["end"], hit["motif"], round(hit["score"], 3), hit["strand"])
```

- **Example output:**

```text
True gm_hits.bed
chr1 13 17 GATA -0.4 +
```

## 4) `dnakit.regulation.scan_region_scores()`

**Example file:** `gm_regions.txt`

```text
chr1:100-300
chr2:400-600
```

- **Purpose:** Score every motif in every region and return a region-by-motif table of best scores.
- **Method:** Calls gimmemotifs `scan_regionfile_to_table()` with sequences from a FASTA genome; no z-normalisation or GC correction is applied, so no GC background is needed. Each chromosome must be longer than 1,280 bases.
- **API:** `dnakit.regulation.scan_region_scores(regions_file, genome_fasta, motifs_pfm, *, verbose=True)`.
- **Input:** A file with one `chrom:start-end` region per line, a FASTA genome and a PFM file.
- **Output:** `regions`, `motifs` and `scores` (one row per region, columns matching `motifs`).
- **Example code:**

```python
from pathlib import Path
Path("gm_genome.fa").write_text(
    ">chr1\n" + "ACGTTGATAAGCCGTAGGATAAGTTCAGTTGATAAGGGCCATTAGAAGTCCTAGGA" * 70 + "\n"
    + ">chr2\n" + "GGGCCCAAATTTGGGCCCAAATTTGATAAGGGCCCAAATTTGGGCCCAAATTTGGGC" * 60 + "\n"
)
from dnakit.regulation import scan_region_scores

result = scan_region_scores("gm_regions.txt", "gm_genome.fa", "gm_motifs.pfm")
print(result["motifs"], result["regions"])
print([round(v, 4) for v in result["scores"][0]])
```

- **Example output:**

```text
['GATA'] ['chr1:100-300', 'chr2:400-600']
[-0.4004]
```

## 5) `dnakit.regulation.motif_activity_scores()`

**Example file:** `gm_labels.tsv`

```text
label	group
r0	1
r1	1
r2	0
r3	0
r4	1
r5	0
```

**Example file:** `gm_counts.tsv`

```text
region	GATA	TTGA
r0	3	0
r1	2	1
r2	0	0
r3	0	1
r4	4	0
r5	0	2
```

- **Purpose:** Compare motif counts in foreground and background regions with a hypergeometric test to get a per-motif activity.
- **Method:** Calls gimmemotifs `moap()` with the `hypergeom` method on a precomputed count table (`scoring="count"`); no genome scan is run.
- **API:** `dnakit.regulation.motif_activity_scores(labels_file, motif_counts_file, *, method="hypergeom", verbose=True)`.
- **Input:** A labels table with a header, indexed by region (1 foreground, 0 background), and a non-negative integer count table with the same index.
- **Output:** `motifs` (one row per motif), `classes` (label classes) and `values` (activity per motif and class).
- **Example code:**

```python
from dnakit.regulation import motif_activity_scores

result = motif_activity_scores("gm_labels.tsv", "gm_counts.tsv")
print(result["classes"], result["motifs"])
for motif, row in zip(result["motifs"], result["values"]):
    print(motif, [round(v, 3) for v in row])
```

- **Example output:**

```text
['0', '1'] ['GATA', 'TTGA']
GATA [-0.0, 0.699]
TTGA [0.0, -0.0]
```

## 6) `dnakit.regulation.motif_score_metric()`

- **Purpose:** Compute one ROC, PR or enrichment statistic for foreground and background motif scores.
- **Method:** Calls the gimmemotifs `rocmetrics` function named by `metric`: `roc_auc`, `pr_auc`, `max_fmeasure`, `mncp`, `max_enrichment`, `recall_at_fdr`, `matches_at_fpr`, `phyper_at_fpr`, `fraction_fpr`, `score_at_fpr`, `enr_at_fpr`, `roc_auc_xlim` or `roc_values`. Non-finite results are returned as None.
- **API:** `dnakit.regulation.motif_score_metric(foreground, background, *, metric="roc_auc", fpr=0.01, fdr=0.1, minbg=2, verbose=True)`.
- **Input:** Foreground and background score lists (finite numbers), the metric name, and the `fpr`, `fdr` and `minbg` parameters.
- **Output:** `value` (a scalar statistic), or `curve` (`fpr` and `tpr` for `roc_values`).
- **Example code:**

```python
from dnakit.regulation import motif_score_metric

roc = motif_score_metric([3, 4, 5, 6], [1, 2, 3, 4], metric="roc_auc")
print(roc["metric"], round(roc["value"], 4))
curve = motif_score_metric([3, 4, 5, 6], [1, 2, 3, 4], metric="roc_values")
print(len(curve["curve"]["fpr"]), len(curve["curve"]["tpr"]))
```

- **Example output:**

```text
roc_auc 0.875
5 5
```

## 7) `dnakit.regulation.motif_position_ks()`

- **Purpose:** Kolmogorov-Smirnov test of motif hit positions, to check whether hits sit toward one end of a sequence.
- **Method:** Calls gimmemotifs `ks_pvalue()` or `ks_significance()`; `output` "pvalue" returns the p-value and "significance" returns -log10(p).
- **API:** `dnakit.regulation.motif_position_ks(foreground_positions, background_positions=None, *, output="pvalue", verbose=True)`.
- **Input:** Foreground hit positions, optional background positions, and the output type.
- **Output:** `output`, `value` (None when not finite) and the sample counts of both sets.
- **Example code:**

```python
from dnakit.regulation import motif_position_ks

fg = [10, 20, 30, 40, 50]
bg = [100, 110, 120, 130, 140]
pvalue = motif_position_ks(fg, bg, output="pvalue")
print(pvalue["output"], round(pvalue["value"], 4))
significance = motif_position_ks(fg, bg, output="significance")
print(round(significance["value"], 4))
```

- **Example output:**

```text
pvalue 0.9616
0.017
```

## 8) `dnakit.regulation.equal_length_motif_pair()`

- **Purpose:** Align two PFMs with a position shift and pad them to equal length.
- **Method:** Counts are first normalised to probabilities per row, then gimmemotifs `make_equal_length()` is called. `pos` shifts `pfm_b` against `pfm_a`; missing positions get 0.25 per base. `truncate` removes the padded ends of one or both matrices.
- **API:** `dnakit.regulation.equal_length_motif_pair(pfm_a, pfm_b, *, pos=0, truncate=None, verbose=True)`.
- **Input:** Two PFMs (four A, C, G, T counts per row), an integer shift (-100 to 100) and the truncation option.
- **Output:** `length` and the aligned `pfm_a` and `pfm_b` (probabilities).
- **Example code:**

```python
from dnakit.regulation import equal_length_motif_pair

result = equal_length_motif_pair([[10, 0, 0, 0], [0, 10, 0, 0]], [[0, 0, 10, 0]], pos=0)
print(result["length"])
print(result["pfm_b"])
```

- **Example output:**

```text
2
[[0.0, 0.0, 1.0, 0.0], [0.25, 0.25, 0.25, 0.25]]
```

## 9) `dnakit.regulation.fasta_record_count()`

**Example file:** `gm_records.fa`

```text
>a
ACGT
>b
GGCC
>c
TTAA
```

- **Purpose:** Count the sequence records in a FASTA file.
- **Method:** Calls gimmemotifs `number_of_seqs_in_file()`, which counts the records that start with `>`.
- **API:** `dnakit.regulation.fasta_record_count(path, *, verbose=True)`.
- **Input:** 个
- **Output:** r
- **Example code:**

```python
from dnakit.regulation import fasta_record_count

result = fasta_record_count("gm_records.fa")
print(result["record_count"], result["file"])
```

- **Example output:**

```text
3 gm_records.fa
```

## 10) `dnakit.regulation.shuffle_region()`

- **Purpose:** Shuffle the bases inside one region of a sequence and keep everything outside it.
- **Method:** Calls tangermeme `ersatz.shuffle()` on a one-hot matrix over `[start, end)`, repeated `n` times with a fixed `random_state`.
- **API:** `dnakit.regulation.shuffle_region(sequence, start, end, *, n=1, random_state=None, verbose=True)`.
- **Input:** A DNA string, a 0-based half-open interval `[start, end)`, the copy count `n` (1 to 1000) and an optional random seed.
- **Output:** `sequences` (`n` shuffled sequences) with `start` and `end`.
- **Example code:**

```python
from dnakit.regulation import shuffle_region

result = shuffle_region("ACGTACGTAA", 2, 6, n=2, random_state=0)
for sequence in result["sequences"]:
    print(sequence[:2], sorted(sequence[2:6]), sequence[6:])
```

- **Example output:**

```text
AC ['A', 'C', 'G', 'T'] GTAA
AC ['A', 'C', 'G', 'T'] GTAA
```

## 11) `dnakit.regulation.randomize_region()`

- **Purpose:** Replace the bases of a region with bases drawn from given A, C, G and T probabilities.
- **Method:** Calls tangermeme `ersatz.randomize()`; `probs` is normalised and applied at every position of the region, uniform by default.
- **API:** `dnakit.regulation.randomize_region(sequence, start, end, *, probs=None, n=1, random_state=None, verbose=True)`.
- **Input:** A DNA string, an interval, four non-negative probabilities (A, C, G, T), the copy count and an optional seed.
- **Output:** `sequences`, the normalised `probs`, and the interval.
- **Example code:**

```python
from dnakit.regulation import randomize_region

result = randomize_region("ACGTACGTAA", 2, 6, probs=[0, 0, 1, 0], random_state=0)
print(result["sequences"][0])
print(result["probs"])
```

- **Example output:**

```text
ACGGGGGTAA
[0.0, 0.0, 1.0, 0.0]
```

## 12) `dnakit.regulation.delete_region()`

- **Purpose:** Delete one stretch of bases from a sequence and return the shorter sequence.
- **Method:** Calls tangermeme `ersatz.delete()`, which removes the one-hot columns of `[start, end)`.
- **API:** `dnakit.regulation.delete_region(sequence, start, end, *, verbose=True)`.
- **Input:** A DNA string and a half-open interval `[start, end)`.
- **Output:** `sequence`, `length` and `deleted`.
- **Example code:**

```python
from dnakit.regulation import delete_region

result = delete_region("ACGTACGTAA", 2, 6)
print(result["sequence"], result["length"], result["deleted"])
```

- **Example output:**

```text
ACGTAA 6 4
```

## 13) `dnakit.regulation.random_one_hot_sequences()`

- **Purpose:** Draw random DNA sequences from given base probabilities, for tests and benchmarks.
- **Method:** Calls tangermeme `utils.random_one_hot()`; the four bases are equally likely when `probs` is omitted.
- **API:** `dnakit.regulation.random_one_hot_sequences(count, length, *, probs=None, random_state=None, verbose=True)`.
- **Input:** The number of sequences (1 to 1000), the length (1 to 100000), optional four probabilities and an optional seed.
- **Output:** `sequences`, with `count`, `length` and the normalised `probs`.
- **Example code:**

```python
from dnakit.regulation import random_one_hot_sequences

print(random_one_hot_sequences(2, 5, probs=[1, 0, 0, 0])["sequences"])
print(random_one_hot_sequences(2, 6, random_state=7)["sequences"])
```

- **Example output:**

```text
['AAAAA', 'AAAAA']
['ATCGTG', 'GACCGT']
```

## 14) `dnakit.regulation.place_motif_set()`

- **Purpose:** Write several motifs into one sequence at given spacings.
- **Method:** Calls tangermeme `ersatz.multisubstitute()`; the first motif starts at `start` and each later motif starts `spacing[i]` bases after the previous one ends, substituted in place.
- **API:** `dnakit.regulation.place_motif_set(sequence, motifs, spacing, *, start, verbose=True)`.
- **Input:** A DNA string, 2 to 10 motifs, a spacing list with one fewer value than the motifs, and the start position.
- **Output:** `sequence`, `start` and `spacing`.
- **Example code:**

```python
from dnakit.regulation import place_motif_set

result = place_motif_set("AAAAAAAAAA", ["GG", "TT"], [2], start=1)
print(result["sequence"])
```

- **Example output:**

```text
AGGAATTAAA
```

## 15) `dnakit.regulation.local_dinucleotide_shuffle()`

- **Purpose:** Dinucleotide-shuffle inside local bins, keeping dinucleotide counts per bin.
- **Method:** Calls tangermeme `ersatz.local_dinucleotide_shuffle()`; sequences are cut into `bin_size` bins and shuffled inside each bin; equal-length sequences are processed together.
- **API:** `dnakit.regulation.local_dinucleotide_shuffle(sequences, *, n=1, bin_size=1024, min_bin_size=512, random_state=None, verbose=True)`.
- **Input:** A list of equal-length DNA strings, the shuffles per sequence, the bin size and a seed.
- **Output:** `sequences` (`n` shuffles per input sequence), with `bin_size` and `min_bin_size`.
- **Example code:**

```python
from dnakit.regulation import local_dinucleotide_shuffle

from collections import Counter
seq = "ACGTTGATAAGCCGTAGGAT" * 75
result = local_dinucleotide_shuffle([seq], n=2, bin_size=512, min_bin_size=256, random_state=0)
print(len(result["sequences"]), len(result["sequences"][0]), len(result["sequences"][0][0]))
print(Counter(result["sequences"][0][0]) == Counter(seq))
```

- **Example output:**

```text
1 2 1500
True
```

## 16) `dnakit.regulation.pwm_consensus_sequence()`

- **Purpose:** Convert a position weight matrix to its consensus DNA sequence.
- **Method:** Calls tangermeme `utils.pwm_consensus()` and takes the base with the largest value in each row.
- **API:** `dnakit.regulation.pwm_consensus_sequence(pwm, *, verbose=True)`.
- **Input:** A matrix with four non-negative numbers (A, C, G, T) per row; no row may be all zero.
- **Output:** `consensus` and `length`.
- **Example code:**

```python
from dnakit.regulation import pwm_consensus_sequence

print(pwm_consensus_sequence([[0.7, 0.1, 0.1, 0.1], [0.1, 0.1, 0.7, 0.1], [0.1, 0.1, 0.1, 0.7]])["consensus"])
```

- **Example output:**

```text
AGT
```

## 17) `dnakit.regulation.example_to_sequence_coords()`

- **Purpose:** Convert coordinates inside each example to genome coordinates.
- **Method:** Calls tangermeme `utils.example_to_fasta_coords()`; the genome position is the locus start plus the example offset, optionally centred with `window`.
- **API:** `dnakit.regulation.example_to_sequence_coords(examples, loci, *, window=None, one_indexed=False, verbose=True)`.
- **Input:** Example rows (`example_idx`, `start`, `end`) and locus rows (`chrom`, `start`, `end`).
- **Output:** `intervals`, each with `chrom`, `start` and `end`.
- **Example code:**

```python
from dnakit.regulation import example_to_sequence_coords

result = example_to_sequence_coords([{"example_idx": 0, "start": 5, "end": 10}], [{"chrom": "chr1", "start": 100, "end": 200}])
print(result["intervals"])
```

- **Example output:**

```text
[{'chrom': 'chr1', 'start': 105, 'end': 110}]
```

## 18) `dnakit.regulation.one_hot_to_fasta_file()`

- **Purpose:** Write DNA sequences to a FASTA file.
- **Method:** Encodes the sequences as one-hot and calls tangermeme `io.one_hot_to_fasta()` to write them; the output file must not exist.
- **API:** `dnakit.regulation.one_hot_to_fasta_file(sequences, output, *, headers=None, verbose=True)`.
- **Input:** A list of equal-length DNA strings, an output path and optional record names.
- **Output:** `output`, `record_count` and `headers`.
- **Example code:**

```python
from dnakit.regulation import one_hot_to_fasta_file

from pathlib import Path
Path("tm_out.fa").unlink(missing_ok=True)
result = one_hot_to_fasta_file(["ACGT", "TTGA"], "tm_out.fa", headers=["a", "b"])
print(result["record_count"], result["headers"])
print(open("tm_out.fa").read().strip())
```

- **Example output:**

```text
2 ['a', 'b']
> a
ACGT
> b
TTGA
```

## 19) `dnakit.regulation.gapped_kmer_features()`

- **Purpose:** Count gapped k-mers (k-mers with a few skipped positions) in each sequence.
- **Method:** Calls tangermeme `kmers.gapped_kmers()` and returns the non-zero counts as `[sequence index, feature index, count]`. Sequences must be longer than the largest window.
- **API:** `dnakit.regulation.gapped_kmer_features(sequences, *, min_k=2, max_k=2, max_gap=1, max_len=3, max_gkmers=5, verbose=True)`.
- **Input:** A list of equal-length DNA strings, with the k, gap and window limits.
- **Output:** `counts` (a list of triples) and `feature_count`.
- **Example code:**

```python
from dnakit.regulation import gapped_kmer_features

result = gapped_kmer_features(["ACGTTGATAAGCCGTAGGAT", "TTGATAAGGGCCATTAGAAG"])
print(result["feature_count"], result["counts"][:4])
```

- **Example output:**

```text
125 [[0, 8, 2], [0, 9, 2], [0, 17, 2], [0, 21, 2]]
```

## 20) `dnakit.regulation.ablate_annotations()`

- **Purpose:** Shuffle each annotated span on its own and record the model output before and after.
- **Method:** Calls tangermeme `ablate_annotations()`: each span is shuffled `n` times inside its sequence with its base composition kept, `random_state` fixes the draws, and the model runs on CPU.
- **API:** `dnakit.regulation.ablate_annotations(model, sequences, annotations, *, n=1, random_state=None, verbose=True)`.
- **Input:** A torch model in eval mode, a list of equal-length DNA strings, annotation rows with `example_idx`, `start` and `end` (0-based, end exclusive), `n` (1 to 1000) and an optional `random_state`.
- **Output:** `before` (one row per annotation), `after` (one row per annotation and shuffle), `annotation_count` and `n`.
- **Example code:**

```python
import torch

class CountC(torch.nn.Module):
    def forward(self, x):
        return x[:, 1, :].sum(dim=-1, keepdim=True)
from dnakit.regulation import ablate_annotations

result = ablate_annotations(CountC().eval(), ["ACGTACGTAA", "TTCCGGAATT"], [{"example_idx": 0, "start": 0, "end": 4}, {"example_idx": 1, "start": 2, "end": 6}], n=2, random_state=0)
print(result["before"])
print(result["after"])
```

- **Example output:**

```text
[[2.0], [2.0]]
[[[2.0], [2.0]], [[2.0], [2.0]]]
```

## 21) `dnakit.regulation.annotate_seqlets()`

- **Purpose:** Match each seqlet to its most similar motifs with TOMTOM and report p-values.
- **Method:** Calls tangermeme `annotate_seqlets()`, which uses the TOMTOM code in memelite (pure Python and numba, no MEME binary); each PWM row is scaled to sum to one, and the `n_nearest` most significant hits are kept.
- **API:** `dnakit.regulation.annotate_seqlets(sequences, seqlets, motifs, *, n_nearest=1, verbose=True)`.
- **Input:** A list of equal-length DNA strings; seqlet rows (`example_idx`, `start`, `end`); a dict from motif name to a PWM with one row of A, C, G, T per position; `n_nearest` (1 to the number of motifs). Use at least 25 motifs.
- **Output:** `matches`: one list per seqlet, each entry with `motif` (None when there is no significant hit) and `p_value`; also `seqlet_count` and `motif_count`.
- **Example code:**

```python
from dnakit.regulation import annotate_seqlets

import numpy

database = {"GATA": [[0, 0, 1, 0], [1, 0, 0, 0], [0, 0, 0, 1], [1, 0, 0, 0]]}
rng = numpy.random.RandomState(7)
for index in range(30):
    database[f"rand{index}"] = rng.dirichlet(numpy.ones(4), size=6).tolist()
result = annotate_seqlets(["TTCCGGAATTGATACCGGTTAACCGG"], [{"example_idx": 0, "start": 10, "end": 14}], database, n_nearest=2)
best = result["matches"][0][0]
print(result["motif_count"], best["motif"], best["p_value"] < 0.001)
```

- **Example output:**

```text
31 GATA True
```

## 22) `dnakit.regulation.count_seqlet_annotations()`

- **Purpose:** Count how often each annotation occurs in each example.
- **Method:** Calls tangermeme `count_annotations()` on `[example_idx, annotation_idx]` rows; counts are int64 so they cannot overflow the uint8 default; `dim` None returns the matrix, while 0 or 1 sums along that axis.
- **API:** `dnakit.regulation.count_seqlet_annotations(pairs, *, dim=None, verbose=True)`.
- **Input:** A list of 1 to 10000 non-negative integer `[example_idx, annotation_idx]` rows and an optional `dim` (None, 0 or 1).
- **Output:** `counts` (a matrix or a flat list) and `dim`.
- **Example code:**

```python
from dnakit.regulation import count_seqlet_annotations

pairs = [[0, 0], [0, 1], [0, 1], [1, 0]]
print(count_seqlet_annotations(pairs)["counts"])
print(count_seqlet_annotations(pairs, dim=1)["counts"])
```

- **Example output:**

```text
[[1, 2], [1, 0]]
[3, 1]
```

## 23) `dnakit.regulation.pairwise_annotation_counts()`

- **Purpose:** Count annotation pairs within an example, and optionally bin them by the gap between the two annotations.
- **Method:** With `max_distance` None it calls tangermeme `pairwise_annotations()`, counting each distinct pair once per example in a symmetric matrix; with an integer `max_distance` it calls `pairwise_annotations_spacing()` and returns an annotation by annotation by distance array, where distance is the number of bases between the end of the left span and the start of the right one.
- **API:** `dnakit.regulation.pairwise_annotation_counts(annotations, *, max_distance=None, symmetric=True, verbose=True)`.
- **Input:** With `max_distance` None, `[example_idx, annotation_idx]` rows; otherwise mapping rows with `example_idx`, `annotation`, `start` and `end`, where spans in one example must not overlap. `symmetric` is True or False.
- **Output:** `mode` (`pairs` or `spacing`), `max_distance`, `counts` (nested lists) and `row_count`.
- **Example code:**

```python
from dnakit.regulation import pairwise_annotation_counts

print(pairwise_annotation_counts([[0, 0], [0, 1], [0, 2], [0, 2]])["counts"])
rows = [{"example_idx": 0, "annotation": 0, "start": 0, "end": 4}, {"example_idx": 0, "annotation": 1, "start": 6, "end": 10}]
spaced = pairwise_annotation_counts(rows, max_distance=5)["counts"]
print(spaced[0][1][2], spaced[1][0][2])
```

- **Example output:**

```text
[[0, 1, 1], [1, 0, 1], [1, 1, 0]]
1 1
```

## 24) `dnakit.regulation.screen_sequence_design()`

- **Purpose:** Screen random DNA sequences for those whose model output is closest to a target.
- **Method:** Calls tangermeme `screen()`: each batch draws `batch_size` random one-hot sequences and scores them, keeping the `n_best` with the smallest squared error; it stops when the worst kept candidate is within `tol` or after `max_iter` batches, so it always ends. `squared_error` is recomputed from the model output of the returned sequences.
- **API:** `dnakit.regulation.screen_sequence_design(model, length, target, *, n_best=1, tol=1e-3, max_iter=100, batch_size=32, random_state=None, verbose=True)`.
- **Input:** A model in eval mode, a sequence length (1 to 1000), a target vector with one number per output, `n_best` (1 to 16), `tol`, `max_iter` (1 to 10000), `batch_size` and an optional `random_state`. With `n_best` above 1, tangermeme fails if two candidates have exactly equal loss, so keep `n_best` at 1 for models with discrete outputs.
- **Output:** `sequences` (best first), `predictions`, `squared_error`, `length` and `n_best`.
- **Example code:**

```python
import torch

class CountA(torch.nn.Module):
    def forward(self, x):
        return x[:, 0, :].sum(dim=-1, keepdim=True)

from dnakit.regulation import screen_sequence_design

result = screen_sequence_design(CountA().eval(), 6, [0.0], random_state=0)
print(result["sequences"], result["squared_error"])
```

- **Example output:**

```text
['TTCTCG'] [0.0]
```

## 25) `dnakit.regulation.greedy_marginalize_design()`

- **Purpose:** Build a construct by greedily placing motifs so that the model output moves toward a target.
- **Method:** Calls tangermeme `greedy_marginalize()`: the first motif starts at `len // 2 - 1`, later motifs are searched within `max_spacing` bases of the placed ones, and a motif is kept only if the mean squared error drops by more than 0.001, up to `max_iter` motifs. With `reverse_complement` the reverse complements are also tried. Uncovered positions are written as `N` and the ends are trimmed.
- **API:** `dnakit.regulation.greedy_marginalize_design(model, sequence, target, motifs, *, max_spacing=12, reverse_complement=True, max_iter=10, batch_size=32, verbose=True)`.
- **Input:** A model in eval mode, a background sequence (at least `2 * (max_spacing + longest motif)` bases), a target vector, 1 to 20 motifs, `max_spacing` (0 to 100), `reverse_complement`, `max_iter` (1 to 100) and `batch_size`.
- **Output:** `sequence` (the construct, which may contain `N`) and `length`.
- **Example code:**

```python
import torch

class CountT(torch.nn.Module):
    def forward(self, x):
        return x[:, 3, :].sum(dim=-1, keepdim=True)

from dnakit.regulation import greedy_marginalize_design

result = greedy_marginalize_design(CountT().eval(), "A" * 60, [4.0], ["TT"], reverse_complement=False, max_spacing=3)
print(result["sequence"], result["length"])
```

- **Example output:**

```text
TTNNNTT 7
```

## 26) `dnakit.regulation.marginalize_annotations()`

- **Purpose:** Insert each annotated span into the centre of background sequences and record the model output before and after.
- **Method:** Calls tangermeme `marginalize_annotations()`: the span replaces the bases of each background sequence at `len // 2 - span // 2`. Results are indexed by annotation and then by background sequence.
- **API:** `dnakit.regulation.marginalize_annotations(model, sequences, background, annotations, *, verbose=True)`.
- **Input:** A model in eval mode, a list of equal-length source sequences, a list of equal-length background sequences, and annotation rows with `example_idx`, `start` and `end`. No span may be longer than the background.
- **Output:** `insert_start`, `before` and `after`; `before` and `after` have the shape annotation by background sequence by output.
- **Example code:**

```python
import torch

class CountC(torch.nn.Module):
    def forward(self, x):
        return x[:, 1, :].sum(dim=-1, keepdim=True)
from dnakit.regulation import marginalize_annotations

result = marginalize_annotations(CountC().eval(), ["TTCCCCTTTT"], ["AAAAAAAAAA"], [{"example_idx": 0, "start": 2, "end": 6}])
print(result["insert_start"], result["before"], result["after"])
```

- **Example output:**

```text
[3] [[[0.0]]] [[[4.0]]]
```

## 27) `dnakit.regulation.pisa_attribution()`

- **Purpose:** Attribute each model output to the bases of each sequence with PISA.
- **Method:** Calls tangermeme `pisa()` on CPU: each sequence is compared with `n_shuffles` dinucleotide shuffles and gradient-times-input is averaged. With `hypothetical` False the result is multiplied by the observed one-hot bases, so only the base present at each position is non-zero. NumPy's global random state is restored after the call.
- **API:** `dnakit.regulation.pisa_attribution(model, sequences, *, n_shuffles=20, hypothetical=False, batch_size=32, random_state=None, verbose=True)`.
- **Input:** A model in eval mode, a list of equal-length DNA strings (at most 100000 bases in total), `n_shuffles` (1 to 200), `hypothetical`, `batch_size` and an optional `random_state`.
- **Output:** `attributions` (nested lists shaped sequences by outputs by 4 by length) and `shape`.
- **Example code:**

```python
import torch
from dnakit.regulation import pisa_attribution

layer = torch.nn.Linear(4 * 30, 2)
with torch.no_grad():
    layer.weight.copy_(torch.linspace(-1.0, 1.0, 2 * 4 * 30).reshape(2, 4 * 30))
    layer.bias.zero_()
model = torch.nn.Sequential(torch.nn.Flatten(), layer).eval()
result = pisa_attribution(model, ["ACGTTGATAAGCCGTAGGATAAGTTCAGTT"], n_shuffles=3, random_state=0)
print(result["shape"])
```

- **Example output:**

```text
[1, 2, 4, 30]
```

## 28) `dnakit.regulation.tfmodisco_seqlets()`

- **Purpose:** Extract seqlets from attribution scores with the TF-MoDISco thresholding rule.
- **Method:** Calls tangermeme `tfmodisco_seqlets()`: windows of `window_size` are summed, a threshold is chosen from a null model to reach about `target_fdr`, and each seqlet keeps `flank` bases on both sides. `attribution` sums the seqlet after trimming, and the result is deterministic.
- **API:** `dnakit.regulation.tfmodisco_seqlets(attributions, *, window_size=21, flank=10, target_fdr=0.2, min_passing_frac=0.03, max_passing_frac=0.2, weak_threshold_for_counting_sign=0.8, verbose=True)`.
- **Input:** A two-dimensional attribution matrix (one row per example and one column per position, rows at least `window_size` long) and the window, flank and threshold settings.
- **Output:** `seqlets` (`example_idx`, `start`, `end` (exclusive) and `attribution`) and `seqlet_count`.
- **Example code:**

```python
import numpy
from dnakit.regulation import tfmodisco_seqlets

rng = numpy.random.RandomState(0)
attr = rng.normal(0, 0.1, size=(3, 60))
attr[0, 20:41] += 1.0
attr[1, 25:46] += 1.0
attr[2, 10:31] += 1.0
result = tfmodisco_seqlets(attr.tolist())
for row in result["seqlets"]:
    print(row["example_idx"], row["start"], row["end"])
```

- **Example output:**

```text
0 10 51
1 15 56
2 0 41
```

## 29) `dnakit.regulation.export_fragment_files()`

**Example file:** `sa_fragments.tsv`

```text
chr1	100	300	cellA	1
chr1	400	600	cellA	1
chr1	700	900	cellA	1
chr2	50	250	cellB	1
chr1	1000	1200	cellB	1
chr2	1200	1400	cellC	1
chr1	2000	2200	cellC	1
```

- **Purpose:** Split a fragment file into one five-column fragment file per cell group.
- **Method:** Imports the fragments with snapATAC2 (every barcode kept), then calls `export_fragments()`; file names are `{prefix}{group}{suffix}`. Each output line is one stored fragment with read count 1. Group names must start with a letter or digit, and existing output files are refused.
- **API:** `dnakit.regulation.export_fragment_files(fragment_file, chrom_sizes, groups, out_dir, *, prefix='', suffix='.bed', verbose=True)`.
- **Input:** A five-column fragment file, a chromosome-length mapping, a barcode-to-group mapping that names every barcode, an existing output directory, and an optional prefix and suffix.
- **Output:** `group_count` and `groups` (per group: `group`, `file` and `fragment_count`, sorted by group name).
- **Example code:**

```python
import shutil
from pathlib import Path

from dnakit.regulation import export_fragment_files

shutil.rmtree("sa_export_out", ignore_errors=True)
Path("sa_export_out").mkdir()
result = export_fragment_files(
    "sa_fragments.tsv",
    {"chr1": 3000, "chr2": 2000},
    {"cellA": "G1", "cellB": "G2", "cellC": "G1"},
    "sa_export_out",
    prefix="sa_",
)
print(result["group_count"])
for group in result["groups"]:
    print(group["group"], group["file"], group["fragment_count"])
```

- **Example output:**

```text
2
G1 sa_G1.bed 5
G2 sa_G2.bed 2
```

## 30) `dnakit.regulation.export_coverage_tracks()`

**Example file:** `sa_fragments.tsv`

```text
chr1	100	300	cellA	1
chr1	400	600	cellA	1
chr1	700	900	cellA	1
chr2	50	250	cellB	1
chr1	1000	1200	cellB	1
chr2	1200	1400	cellC	1
chr1	2000	2200	cellC	1
```

- **Purpose:** Write one genome-wide coverage track per cell group, as bedGraph or bigWig.
- **Method:** Imports the fragments, then calls `export_coverage()` with fragment counting and no normalization. `bin_size` sets the bin width and an empty `max_fragment_length` applies no length filter. Files are named `{prefix}{group}.bedgraph` or `.bw`.
- **API:** `dnakit.regulation.export_coverage_tracks(fragment_file, chrom_sizes, groups, out_dir, *, bin_size=10, normalization=None, output_format='bedgraph', max_fragment_length=None, prefix='', verbose=True)`.
- **Input:** A fragment file, a chromosome-length mapping, a barcode-to-group mapping, an existing output directory, a bin width (1 to 10000000), a normalization (None, RPKM, CPM or BPM), an output format and a fragment length limit.
- **Output:** `track_count`, `output_format`, `bin_size`, `normalization` and `tracks` (per group: `group` and `file`).
- **Example code:**

```python
import shutil
from pathlib import Path

from dnakit.regulation import export_coverage_tracks

shutil.rmtree("sa_coverage_out", ignore_errors=True)
Path("sa_coverage_out").mkdir()
result = export_coverage_tracks(
    "sa_fragments.tsv",
    {"chr1": 3000, "chr2": 2000},
    {"cellA": "G1", "cellB": "G2", "cellC": "G1"},
    "sa_coverage_out",
    bin_size=100,
)
print(result["track_count"], [track["file"] for track in result["tracks"]])
print(Path("sa_coverage_out/G2.bedgraph").read_text(encoding="utf-8").splitlines()[0])
```

- **Example output:**

```text
2 ['G1.bedgraph', 'G2.bedgraph']
chr1	1000	1200	1
```

## 31) `dnakit.regulation.fragments_from_bam()`

**Example file:** `sa_reads.sam`

```text
@HD	VN:1.6	SO:unsorted
@SQ	SN:chr1	LN:3000
@SQ	SN:chr2	LN:2000
r1	99	chr1	101	60	10M	=	391	300	AAAAAAAAAA	IIIIIIIIII	CB:Z:AAAA
r1	147	chr1	391	60	10M	=	101	-300	AAAAAAAAAA	IIIIIIIIII	CB:Z:AAAA
r2	99	chr1	601	60	10M	=	791	200	AAAAAAAAAA	IIIIIIIIII	CB:Z:AAAA
r2	147	chr1	791	60	10M	=	601	-200	AAAAAAAAAA	IIIIIIIIII	CB:Z:AAAA
r3	99	chr2	51	60	10M	=	241	200	AAAAAAAAAA	IIIIIIIIII	CB:Z:CCCC
r3	147	chr2	241	60	10M	=	51	-200	AAAAAAAAAA	IIIIIIIIII	CB:Z:CCCC
```

- **Purpose:** Convert a barcode-tagged BAM file into a five-column fragment file and return the QC fractions.
- **Method:** Calls `snapATAC2.pp.make_fragment_file()` with the CB tag as barcode, a minimum MAPQ of 30 and Tn5 shifts (+4 at the start, -5 at the end). The example converts the SAM text to BAM with `pysam.view` first.
- **API:** `dnakit.regulation.fragments_from_bam(bam_file, output_fragment_file, *, barcode_tag='CB', is_paired=True, min_mapq=30, mitochondrial_chroms=('chrM', 'M'), verbose=True)`.
- **Input:** A paired-end BAM file, a fragment output path that does not yet exist, a two-character SAM tag, a MAPQ floor (0 to 255) and a list of mitochondrial chromosome names.
- **Output:** `fragment_count`, `barcode_count`, `qc` (snapATAC2 fractions) and `output_file`.
- **Example code:**

```python
from pathlib import Path

import pysam

from dnakit.regulation import fragments_from_bam

Path("sa_reads.bam").unlink(missing_ok=True)
Path("sa_bam_fragments.tsv").unlink(missing_ok=True)
Path("sa_reads.bam").write_bytes(pysam.view("-b", "sa_reads.sam"))
result = fragments_from_bam("sa_reads.bam", "sa_bam_fragments.tsv")
print(result["fragment_count"], result["barcode_count"])
print(Path("sa_bam_fragments.tsv").read_text(encoding="utf-8"), end="")
```

- **Example output:**

```text
3 2
chr1	104	395	AAAA	1
chr1	604	795	AAAA	1
chr2	54	245	CCCC	1
```

## 32) `dnakit.regulation.import_contact_pairs()`

**Example file:** `sa_contacts.tsv`

```text
cellA	chr1	10	chr1	40	1
cellA	chr1	120	chr1	950	2
cellB	chr2	300	chr2	700	1
cellB	chr1	50	chr2	800	1
```

- **Purpose:** Bin single-cell chromatin contact pairs into an hdf5 AnnData file.
- **Method:** Each row is `barcode, chrom1, pos1, chrom2, pos2, count`, with 0-based positions that must lie inside the chromosome lengths. Calls `snapATAC2.pp.import_contacts()`, with cells as rows and bin pairs as columns.
- **API:** `dnakit.regulation.import_contact_pairs(contact_file, chrom_sizes, output_h5ad, *, bin_size=500000, sorted_by_barcode=False, verbose=True)`.
- **Input:** A contact pair file, a chromosome-length mapping, an output path ending in `.h5ad` that does not yet exist, a bin width (1 to 10000000) and the sorting flag.
- **Output:** `cell_count`, `feature_count` (number of bin pairs), `contact_pair_count`, `contact_count` (sum of counts), `chromosomes` and `output_file`.
- **Example code:**

```python
from pathlib import Path

from dnakit.regulation import import_contact_pairs

Path("sa_contacts.h5ad").unlink(missing_ok=True)
result = import_contact_pairs(
    "sa_contacts.tsv",
    {"chr1": 1000, "chr2": 1000},
    "sa_contacts.h5ad",
    bin_size=500,
)
print(result["cell_count"], result["feature_count"], result["contact_count"])
```

- **Example output:**

```text
2 16 5
```

## 33) `dnakit.regulation.import_base_pair_values()`

**Example file:** `sa_cellA.tsv`

```text
chrom	pos	methyl	unmethyl
chr1	10	3	7
chr1	20	0	4
chr2	5	2	0
```

**Example file:** `sa_cellB.tsv`

```text
chrom	pos	methyl	unmethyl
chr1	10	1	1
```

- **Purpose:** Import one per-cell base-level count table into an hdf5 AnnData file, for example WGBS methylation counts.
- **Method:** The input directory may hold only `.tsv` files, one per cell, and a cell is named by its file stem. Each header must be `chrom, pos, methyl, unmethyl`, with each position once. Calls `snapATAC2.pp.import_values()` and checks each cell's value count against its file.
- **API:** `dnakit.regulation.import_base_pair_values(input_dir, chrom_sizes, output_h5ad, *, verbose=True)`.
- **Input:** A directory holding one `.tsv` file per cell and nothing else, a chromosome-length mapping, and an output path ending in `.h5ad` that does not yet exist.
- **Output:** `cell_count`, `total_value_count`, `reference_length`, `cells` sorted by cell (each with `cell` and `value_count`) and `output_file`.
- **Example code:**

```python
import shutil
from pathlib import Path

from dnakit.regulation import import_base_pair_values

shutil.rmtree("sa_values_dir", ignore_errors=True)
Path("sa_values_dir").mkdir()
Path("sa_values.h5ad").unlink(missing_ok=True)
for name in ("sa_cellA.tsv", "sa_cellB.tsv"):
    shutil.copy(name, "sa_values_dir")
result = import_base_pair_values(
    "sa_values_dir",
    {"chr1": 100, "chr2": 50},
    "sa_values.h5ad",
)
print(result["cell_count"], result["total_value_count"])
for cell in result["cells"]:
    print(cell["cell"], cell["value_count"])
```

- **Example output:**

```text
2 4
sa_cellA 3
sa_cellB 1
```

## 34) `dnakit.regulation.tenx_atac_metrics()`

**Example file:** `sa_make_10x_bam.py`

```text
import random

import pysam


def write_bam(path):
    """30 cells with 300 read pairs inside a peak, and 120 background barcodes with 3 pairs."""

    header = {
        "HD": {"VN": "1.6", "SO": "unsorted"},
        "SQ": [{"SN": "chr1", "LN": 3000}, {"SN": "chr2", "LN": 2000}],
    }
    rng = random.Random(7)
    reads = []
    for index in range(30):
        for k in range(300):
            reads.append((f"C{index:03d}_{k}", f"C{index:03d}", 0, rng.randint(900, 1400), rng.randint(160, 240)))
    for index in range(120):
        for k in range(3):
            tid = rng.choice([0, 1])
            reads.append((f"B{index:03d}_{k}", f"B{index:03d}", tid, rng.randint(0, 1700), rng.randint(160, 240)))
    rng.shuffle(reads)
    with pysam.AlignmentFile(path, "wb", header=header) as out:
        for name, barcode, tid, start, tlen in reads:
            for read1 in (True, False):
                record = pysam.AlignedSegment()
                record.query_name = name
                record.query_sequence = "A" * 50
                record.flag = 0x1 | 0x2 | (0x40 if read1 else 0x80) | (0 if read1 else 0x10)
                record.reference_id = tid
                record.reference_start = start if read1 else start + tlen - 50
                record.mapping_quality = 60
                record.cigarstring = "50M"
                record.next_reference_id = tid
                record.next_reference_start = start + tlen - 50 if read1 else start
                record.template_length = tlen if read1 else -tlen
                record.query_qualities = pysam.qualitystring_to_array("I" * 50)
                record.set_tag("CB", barcode, value_type="Z")
                out.write(record)
```

**Example file:** `sa_genes.gtf`

```text
chr1	t	gene	1000	1500	.	+	.	gene_id "G1"; gene_name "GATA1";
chr1	t	transcript	1000	1500	.	+	.	gene_id "G1"; transcript_id "T1"; gene_name "GATA1";
chr1	t	exon	1000	1500	.	+	.	gene_id "G1"; transcript_id "T1"; gene_name "GATA1";
```

**Example file:** `sa_peaks.bed`

```text
chr1	900	1400	p1	10	.
```

- **Purpose:** Run the 10x ATAC quality-control recipe on a barcoded BAM and return sequencing, cell, complexity, mapping and targeting metrics.
- **Method:** Calls `snapATAC2.pp.recipe_10x_metrics()`, which converts the BAM to fragments, imports them, computes TSS enrichment from the local GTF and derives peak metrics from the BED peaks. The example BAM comes from `sa_make_10x_bam.py`. MACS3 peak calling is not run; peaks must be supplied as a BED file.
- **API:** `dnakit.regulation.tenx_atac_metrics(bam_file, output_fragment_file, output_h5ad, chrom_sizes, gene_gtf, peaks_bed, *, barcode_tag='CB', verbose=True)`.
- **Input:** A BAM with a CB tag, two output paths that do not yet exist (a fragment file and an `.h5ad`), a chromosome-length mapping, a local GTF gene annotation and a BED peak file.
- **Output:** `metrics` (a nested dictionary with the snapATAC2 recipe keys), `fragment_file` and `h5ad_file`.
- **Example code:**

```python
from pathlib import Path

from sa_make_10x_bam import write_bam

from dnakit.regulation import tenx_atac_metrics

for name in ("sa_10x.bam", "sa_10x_fragments.tsv", "sa_10x.h5ad"):
    Path(name).unlink(missing_ok=True)
write_bam("sa_10x.bam")
result = tenx_atac_metrics(
    "sa_10x.bam",
    "sa_10x_fragments.tsv",
    "sa_10x.h5ad",
    {"chr1": 3000, "chr2": 2000},
    "sa_genes.gtf",
    "sa_peaks.bed",
)
metrics = result["metrics"]
print(metrics["Sequencing"]["Sequenced_read_pairs"], metrics["Targeting"]["Number_of_peaks"], metrics["Cells"]["Number_of_cells"])
```

- **Example output:**

```text
9360.0 1 31
```

## 35) `dnakit.regulation.call_peaks_macs3()`

- **Purpose:** Call open chromatin peaks with MACS3 separately for each cell group.
- **Method:** Imports the fragments into a temporary snapATAC2 file, calls `snapATAC2.tl.macs3()` for each group and sorts each group's peaks by chromosome and position. The example generates its fragments in code and prints the first peak of each group.
- **API:** `dnakit.regulation.call_peaks_macs3(fragment_file, chrom_sizes, groups, *, qvalue=0.05, shift=-100, extsize=200, verbose=True)`.
- **Input:** A fragment file, a chromosome-length mapping, a group mapping that names every barcode, and the q-value threshold (above 0 and at most 1), shift and extension length.
- **Output:** `group_count` and `groups` (per group: `group`, `peak_count` and `peaks`, with MACS3 columns such as `chrom`, `start`, `end` and `q_value`).
- **Example code:**

```python
import random

from dnakit.regulation import call_peaks_macs3

rng = random.Random(3)
rows = []
for _ in range(60):
    start = 1000 + rng.randint(0, 200)
    rows.append(f"chr1\t{start}\t{start + 150}\tcellA\t1")
for _ in range(60):
    start = 500 + rng.randint(0, 100)
    rows.append(f"chr2\t{start}\t{start + 150}\tcellB\t1")
for _ in range(20):
    start = rng.randint(0, 1800)
    rows.append(f"chr1\t{start}\t{start + 150}\tcellC\t1")
with open("sa_clustered.tsv", "w", encoding="utf-8") as handle:
    handle.write("\n".join(rows) + "\n")
result = call_peaks_macs3(
    "sa_clustered.tsv",
    {"chr1": 3000, "chr2": 2000},
    {"cellA": "G1", "cellB": "G2", "cellC": "G1"},
)
for group in result["groups"]:
    peak = group["peaks"][0]
    print(group["group"], group["peak_count"], peak["chrom"], peak["start"], peak["end"])
```

- **Example output:**

```text
G1 1 chr1 934 1420
G2 1 chr2 415 834
```

## 36) `dnakit.regulation.add_tf_binding_edges()`

**Example file:** `sa_network.gtf`

```text
chr1	t	gene	900	1500	.	+	.	gene_id "G1"; gene_name "GATA1";
chr1	t	transcript	900	1500	.	+	.	gene_id "G1"; transcript_id "T1"; gene_name "GATA1";
chr1	t	exon	900	1500	.	+	.	gene_id "G1"; transcript_id "T1"; gene_name "GATA1";
chr2	t	gene	400	800	.	+	.	gene_id "G2"; gene_name "OTHER";
chr2	t	transcript	400	800	.	+	.	gene_id "G2"; transcript_id "T2"; gene_name "OTHER";
chr2	t	exon	400	800	.	+	.	gene_id "G2"; transcript_id "T2"; gene_name "OTHER";
```

**Example file:** `sa_motifs.meme`

```text
MEME version 4

ALPHABET= ACGT

strands: + -

Background letter frequencies
A 0.25 C 0.25 G 0.25 T 0.25

MOTIF GATA1
letter-probability matrix: alength= 4 w= 10 nsites= 20 E= 0
 0 0 1 0
 1 0 0 0
 0 0 0 1
 1 0 0 0
 1 0 0 0
 0 0 1 0
 0 0 1 0
 1 0 0 0
 0 0 0 1
 0 1 0 0
```

- **Purpose:** Build a region-to-gene network from a GTF and add motif-to-region binding edges by scanning a local genome.
- **Method:** Builds the network with `snapATAC2.tl.init_network_from_annotation()` using gene names and 250 kb flanks, then calls `snapATAC2.tl.add_tf_binding()` to scan the MEME motifs in the region sequences at a p-value threshold. Only motif-to-region edges are returned. The example writes its genome in code, with the motif at chr1:1100-1110.
- **API:** `dnakit.regulation.add_tf_binding_edges(regions, gtf_file, motif_file, genome_fasta, *, pvalue=1e-5, upstream=250000, downstream=250000, coding_gene_only=True, verbose=True)`.
- **Input:** A list of `chrom:start-end` regions, a local GTF, a MEME motif file, a local FASTA genome, and the p-value threshold (above 0 and at most 1) with the flank distances.
- **Output:** `motif_count`, `region_count`, `gene_count`, `binding_edge_count` and `binding_edges` sorted by motif and region.
- **Example code:**

```python
from pathlib import Path

from dnakit.regulation import add_tf_binding_edges

sequence = "A" * 1100 + "GATAAGGATC" + "A" * 1890
Path("sa_genome.fa").write_text(">chr1\n" + sequence + "\n>chr2\n" + "C" * 2000 + "\n", encoding="utf-8")
result = add_tf_binding_edges(
    ["chr1:1000-1200", "chr1:2000-2200", "chr2:500-700"],
    "sa_network.gtf",
    "sa_motifs.meme",
    "sa_genome.fa",
)
print(result["motif_count"], result["region_count"], result["gene_count"])
for edge in result["binding_edges"]:
    print(edge["motif"], edge["region"])
```

- **Example output:**

```text
1 3 2
GATA1 chr1:1000-1200
```

## How these functions are verified

The test file is `others/tests/unit/test_review_G18i.py`. Each wrapped function has at least one test that matches a direct upstream call or a hand calculation. Missing-package install hints are tested with monkeypatch. Random functions take an explicit seed and do not change global random state.


## Current scope

- gimmemotifs: 17 functions need external data (genomepy genomes and annotations, GC indexes, ortholog databases, OrthoFinder output, ROC report tables or genome-name backgrounds); 2 need de novo output from MEME or HOMER; `run_maelstrom` is a command-line pipeline; the rest are output or input helpers, aliases and plot data. `scan_to_file` with an `fpr` threshold needs a background file, so the wrapper exposes only `cutoff`.
- tangermeme: `chunk`, `unchunk`, `extract_signal` and `entropy` are array helpers and are not wrapped on their own; `ersatz.random_one_hot`, `design.screen.random_one_hot`, `space.multisubstitute` and `ablate.shuffle` are duplicates of wrapped functions.
- snapATAC2: `link_tf_to_gene`, `prune_network` and `pagerank` run on an in-memory regulatory network with no sequence or interval step, so they are recorded as `not-dna`.

