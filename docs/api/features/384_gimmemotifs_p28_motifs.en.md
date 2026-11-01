# GimmeMotifs motif comparison, matrix transformation and scans

These functions in `dnakit.regulation` wrap the motif-level routines of [GimmeMotifs](https://github.com/vanheeringen-lab/gimmemotifs) 0.18.4: motif comparison, the transformation and export of motif matrices, and motif scans over FASTA sequences. A motif is given as a PFM, a list with one row per position and the counts of A, C, G and T in each row (probability rows are accepted too, as in GimmeMotifs). Every function calls the GimmeMotifs routine of the same name and returns a JSON-serialisable dictionary with `status`, the result and a `provenance` entry that names the routine and the GimmeMotifs version. Nothing is downloaded.

GimmeMotifs is an optional dependency (`pip install "dnakit[bio-gimmemotifs]"`) and is imported only when a function is called. The examples import from the module `dnakit.regulation.gimmemotifs_p28_motif_functions`; after the area package is wired the same names are available from `dnakit.regulation`.

These functions cover the GimmeMotifs routines `MotifComparer.compare_motifs`, `MotifComparer.pvalue`, `MotifComparer.generate_score_dist`, `comparison.akl` and `ssd` (already available as `compare_motif_positions`), `Motif.pcc`, `max_pcc`, `ic`, `matrix_ic`, `pwm_min_score`, `pwm_max_score`, `rc`, `shuffle`, `pfm_to_ppm`, `to_motevo`, `to_meme`, `to_ppm`, `to_consensusv2`, `average_motifs`, `sample`, `pwm_scan_score`, `pwm_scan_to_gff`, `consensus_scan` and `utils.ppm_pseudocount`.

## Example data

The scan examples read one three-record FASTA file; the other examples pass the motif matrices `A` (4 positions, consensus ACGT) and `B` (5 positions) directly.

**Example file:** `seqs.fa`

```text
>s1
TTACGTTTACGTAA
>s2
GGGGACGTCCCC
>s3
CCCCCCCCCC
```

**Example file:** `set.pfm`

```text
>a
10 0 0 0
0 10 0 0
0 0 10 0
0 0 0 10
>b
0 10 0 0
0 0 10 0
0 0 0 10
10 0 0 0
>c
5 5 0 0
0 5 5 0
0 0 5 5
5 0 0 5
>d
10 0 0 0
10 0 0 0
0 0 10 0
0 0 10 0
0 10 0 0
```

## 1) `dnakit.regulation.gimmemotifs_compare_motif_pair()`

- **Purpose:** Compare two motifs over all alignments and strands and return the best similarity score, or its p-value.
- **Method:** Build two `Motif` objects and call `MotifComparer.compare_motifs`. `metric` is `wic` (weighted information content, the default), `pcc`, `ed`, `distance`, `chisq`, `akl`, `ssd` or `seqcor`. `match` is `total` (the whole match, including unaligned positions), `subtotal` or `partial` (overlapping positions only); `combine` is `mean` or `sum`. With `pvalue=True` the score is converted by `MotifComparer.pvalue` to the upper tail of a normal distribution fitted per pair of motif lengths. GimmeMotifs ships that distribution only for `wic`, `total` and `mean`; for any other setting pass `score_distribution`, a file written by `gimmemotifs_motif_score_distribution`.
- **API:** `dnakit.regulation.gimmemotifs_compare_motif_pair(pfm_a, pfm_b, *, metric="wic", match="total", combine="mean", pvalue=False, score_distribution=None)`.
- **Input:** Two PFMs. The options above; `score_distribution` only together with `pvalue=True`.
- **Output:** `value` (score or p-value, see `value_kind`), `position` (offset of the second motif against the first), `strand` (1 or -1), the motif lengths and `provenance`.
- **Example code:**

```python
from dnakit.regulation.gimmemotifs_p28_motif_functions import gimmemotifs_compare_motif_pair

A = [[10, 0, 0, 0], [0, 10, 0, 0], [0, 0, 8, 2], [1, 1, 1, 7]]
B = [[0, 10, 0, 0], [0, 0, 10, 0], [1, 1, 1, 7], [5, 5, 0, 0], [0, 0, 0, 10]]

for kwargs in ({}, {"metric": "pcc"}, {"metric": "akl", "match": "partial"}, {"pvalue": True}):
    result = gimmemotifs_compare_motif_pair(A, B, **kwargs)
    print(kwargs, result["value_kind"], round(result["value"], 4), result["position"], result["strand"])
```

- **Example output:**

```text
{} score -1.6806 1 1
{'metric': 'pcc'} score 0.4948 1 1
{'metric': 'akl', 'match': 'partial'} score 9.4234 -2 -1
{'pvalue': True} pvalue 0.1521 1 1
```

## 2) `dnakit.regulation.gimmemotifs_motif_score_distribution()`

- **Purpose:** Compute the score distribution of all motif pairs by length, the table behind motif-similarity p-values.
- **Method:** Read the motifs of a PFM file and run `MotifComparer.generate_score_dist` serially: every motif is compared with every motif, and for each pair of motif lengths the mean and standard deviation of the scores are written as `length1<TAB>length2<TAB>mean<TAB>sd`. Upstream writes the file into the GimmeMotifs installation; here it goes to `output_file`, which must not exist, so the installed tables stay untouched. `metric` is `wic`, `pcc`, `ed` or `distance`; `match` is `total` or `subtotal`.
- **API:** `dnakit.regulation.gimmemotifs_motif_score_distribution(motifs_pfm, output_file, *, metric="wic", match="total", combine="mean")`.
- **Input:** A PFM file with 1 to 1,000 motifs and a new output path.
- **Output:** The `distribution` rows (`length_a`, `length_b`, `mean`, `sd`), `n_motifs` and the written file in `outputs`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.regulation.gimmemotifs_p28_motif_functions import gimmemotifs_compare_motif_pair, gimmemotifs_motif_score_distribution

A = [[10, 0, 0, 0], [0, 10, 0, 0], [0, 0, 8, 2], [1, 1, 1, 7]]
B = [[0, 10, 0, 0], [0, 0, 10, 0], [1, 1, 1, 7], [5, 5, 0, 0], [0, 0, 0, 10]]

out = Path(tempfile.mkdtemp()) / "pcc_dist.txt"
dist = gimmemotifs_motif_score_distribution("set.pfm", out, metric="pcc")
print(dist["n_motifs"], [(r["length_a"], r["length_b"]) for r in dist["distribution"]])
scored = gimmemotifs_compare_motif_pair(A, B, metric="pcc", pvalue=True, score_distribution=out)
print(scored["value_kind"], 0 <= scored["value"] <= 1)
```

- **Example output:**

```text
4 [(4, 4), (4, 5), (5, 4), (5, 5)]
pvalue True
```

## 3) `dnakit.regulation.gimmemotifs_motif_pcc()`

- **Purpose:** Pearson correlation of two motifs' position probabilities, at one alignment or at the best of all alignments.
- **Method:** With `position` (the offset of the second motif against the first; -2 puts it two positions to the left) the function calls `Motif.pcc`: the two probability matrices are aligned, the unaligned ends are cut off, and the Pearson correlations of the four-value rows are summed. Without `position` it calls `Motif.max_pcc`, which tries every offset, also for the reverse complement of the second motif when `revcomp` is true, and keeps the highest sum.
- **API:** `dnakit.regulation.gimmemotifs_motif_pcc(pfm_a, pfm_b, *, position=None, revcomp=True)`.
- **Input:** Two PFMs, an optional integer offset and `revcomp`.
- **Output:** `mode` (`pcc` or `max_pcc`), `score`, `position` and `strand` (1 or -1).
- **Example code:**

```python
from dnakit.regulation.gimmemotifs_p28_motif_functions import gimmemotifs_motif_pcc

A = [[10, 0, 0, 0], [0, 10, 0, 0], [0, 0, 8, 2], [1, 1, 1, 7]]
B = [[0, 10, 0, 0], [0, 0, 10, 0], [1, 1, 1, 7], [5, 5, 0, 0], [0, 0, 0, 10]]

best = gimmemotifs_motif_pcc(A, B)
print(best["mode"], round(best["score"], 4), best["position"], best["strand"])
fixed = gimmemotifs_motif_pcc(A, B, position=-1)
print(fixed["mode"], round(fixed["score"], 4), fixed["position"])
```

- **Example output:**

```text
max_pcc 2.9685 1 1
pcc -0.4291 -1
```

## 4) `dnakit.regulation.gimmemotifs_motif_ic_similarity()`

- **Purpose:** Information-content similarity of two motifs, together with the information content of each motif.
- **Method:** `mode="matrix"` calls `Motif.matrix_ic`: both strands and every offset are scored with the information content shared by the two probability rows relative to the uniform background, and the best score, its offset and its strand are returned. `mode="aligned"` calls `Motif.ic` for the single alignment at `position`; positions that overlap nothing are compared with the background and weighted by `bg_factor`. `Motif.information_content` (bits, 2 per fully conserved position) is reported for both motifs. GimmeMotifs `matrix_ic` raises a broadcasting error for some pairs in which the second motif is longer than the first (an upstream defect); the function reports that error with a hint to swap the motifs or use `aligned`.
- **API:** `dnakit.regulation.gimmemotifs_motif_ic_similarity(pfm_a, pfm_b, *, mode="matrix", position=0, bg_factor=1.0)`.
- **Input:** Two PFMs, `mode`, an integer `position` and a non-negative `bg_factor`.
- **Output:** `score`, `position`, `strand`, `information_content_a` and `information_content_b`.
- **Example code:**

```python
from dnakit.regulation.gimmemotifs_p28_motif_functions import gimmemotifs_motif_ic_similarity

A = [[10, 0, 0, 0], [0, 10, 0, 0], [0, 0, 8, 2], [1, 1, 1, 7]]
B = [[0, 10, 0, 0], [0, 0, 10, 0], [1, 1, 1, 7], [5, 5, 0, 0], [0, 0, 0, 10]]

best = gimmemotifs_motif_ic_similarity(A, B)
print(best["mode"], round(best["score"], 4), best["position"], best["strand"])
print(round(best["information_content_a"], 4), round(best["information_content_b"], 4))
one = gimmemotifs_motif_ic_similarity(A, B, mode="aligned", position=1)
print(one["mode"], round(one["score"], 4))
```

- **Example output:**

```text
matrix 3.5532 1 1
5.909 7.6267
aligned -14.7849
```

## 5) `dnakit.regulation.gimmemotifs_motif_score_range()`

- **Purpose:** Minimum and maximum PWM (log-odds) score a motif can reach.
- **Method:** Build the `Motif` and read `min_score` and `max_score`, the values behind `Motif.pwm_min_score` and `pwm_max_score` (deprecated upstream in favour of the properties). The log-odds of a cell is `log(p / 0.25 + 0.01)`; the scores are the sums of the smallest and the largest value of each position.
- **API:** `dnakit.regulation.gimmemotifs_motif_score_range(pfm)`.
- **Input:** One PFM.
- **Output:** `min_score`, `max_score` and `length`.
- **Example code:**

```python
from dnakit.regulation.gimmemotifs_p28_motif_functions import gimmemotifs_motif_score_range

A = [[10, 0, 0, 0], [0, 10, 0, 0], [0, 0, 8, 2], [1, 1, 1, 7]]
B = [[0, 10, 0, 0], [0, 0, 10, 0], [1, 1, 1, 7], [5, 5, 0, 0], [0, 0, 0, 10]]

result = gimmemotifs_motif_score_range(A)
print(result["length"], round(result["min_score"], 4), round(result["max_score"], 4))
```

- **Example output:**

```text
4 -14.5894 4.9759
```

## 6) `dnakit.regulation.gimmemotifs_motif_consensus()`

- **Purpose:** IUPAC consensus sequence of a motif.
- **Method:** Call `Motif.to_consensusv2`, which returns the GimmeMotifs consensus: a base in upper case when it dominates a position, an IUPAC code in upper case for two or three bases that together dominate, a lower-case code for weaker ambiguity and `n` for an uninformative position.
- **API:** `dnakit.regulation.gimmemotifs_motif_consensus(pfm)`.
- **Input:** One PFM.
- **Output:** `consensus` and `length`.
- **Example code:**

```python
from dnakit.regulation.gimmemotifs_p28_motif_functions import gimmemotifs_motif_consensus

A = [[10, 0, 0, 0], [0, 10, 0, 0], [0, 0, 8, 2], [1, 1, 1, 7]]
B = [[0, 10, 0, 0], [0, 0, 10, 0], [1, 1, 1, 7], [5, 5, 0, 0], [0, 0, 0, 10]]

print(gimmemotifs_motif_consensus(A)["consensus"])
print(gimmemotifs_motif_consensus(B)["consensus"])
print(gimmemotifs_motif_consensus([[3, 1, 1, 1], [2, 2, 2, 2]])["consensus"])
```

- **Example output:**

```text
ACGT
CGTmT
nn
```

## 7) `dnakit.regulation.gimmemotifs_motif_reverse_complement()`

- **Purpose:** Reverse-complement a motif matrix.
- **Method:** Call `Motif.rc`: the order of the positions is reversed and the A, C, G, T columns are reversed, so that A pairs with T and C with G. The motif is rebuilt from the counts.
- **API:** `dnakit.regulation.gimmemotifs_motif_reverse_complement(pfm)`.
- **Input:** One PFM.
- **Output:** `pfm` (counts), `ppm` (probabilities) and `consensus` of the reverse complement.
- **Example code:**

```python
from dnakit.regulation.gimmemotifs_p28_motif_functions import gimmemotifs_motif_reverse_complement

A = [[10, 0, 0, 0], [0, 10, 0, 0], [0, 0, 8, 2], [1, 1, 1, 7]]
B = [[0, 10, 0, 0], [0, 0, 10, 0], [1, 1, 1, 7], [5, 5, 0, 0], [0, 0, 0, 10]]

result = gimmemotifs_motif_reverse_complement(A)
print(result["pfm"])
print(result["consensus"])
```

- **Example output:**

```text
[[7.0, 1.0, 1.0, 1.0], [2.0, 8.0, 0.0, 0.0], [0.0, 0.0, 10.0, 0.0], [0.0, 0.0, 0.0, 10.0]]
ACGT
```

## 8) `dnakit.regulation.gimmemotifs_motif_shuffle()`

- **Purpose:** Shuffle the positions of a motif, the usual control motif of the same composition.
- **Method:** Call `Motif.shuffle`, which permutes the rows of the matrix with numpy and keeps the order of the A, C, G, T columns. Upstream draws from numpy's global random state; `seed` seeds it before the call (and the previous state is restored afterwards), so the permutation is reproducible.
- **API:** `dnakit.regulation.gimmemotifs_motif_shuffle(pfm, *, seed=None)`.
- **Input:** One PFM and an optional integer `seed`.
- **Output:** `pfm`, `ppm` and `seed`.
- **Example code:**

```python
from dnakit.regulation.gimmemotifs_p28_motif_functions import gimmemotifs_motif_shuffle

A = [[10, 0, 0, 0], [0, 10, 0, 0], [0, 0, 8, 2], [1, 1, 1, 7]]
B = [[0, 10, 0, 0], [0, 0, 10, 0], [1, 1, 1, 7], [5, 5, 0, 0], [0, 0, 0, 10]]

first = gimmemotifs_motif_shuffle(A, seed=1)
again = gimmemotifs_motif_shuffle(A, seed=1)
print(first["pfm"])
print(first["pfm"] == again["pfm"])
```

- **Example output:**

```text
[[1.0, 1.0, 1.0, 7.0], [0.0, 0.0, 8.0, 2.0], [10.0, 0.0, 0.0, 0.0], [0.0, 10.0, 0.0, 0.0]]
True
```

## 9) `dnakit.regulation.gimmemotifs_pfm_to_ppm()`

- **Purpose:** Convert motif counts to frequencies with a pseudocount.
- **Method:** Call `Motif.pfm_to_ppm`: every cell becomes `(count + pseudo) / (row sum + 4 * pseudo)`, so each row sums to 1 and no value is zero.
- **API:** `dnakit.regulation.gimmemotifs_pfm_to_ppm(pfm, *, pseudo=0.001)`.
- **Input:** One PFM of counts and a non-negative `pseudo`.
- **Output:** `ppm`, one row of four probabilities per position.
- **Example code:**

```python
from dnakit.regulation.gimmemotifs_p28_motif_functions import gimmemotifs_pfm_to_ppm

A = [[10, 0, 0, 0], [0, 10, 0, 0], [0, 0, 8, 2], [1, 1, 1, 7]]
B = [[0, 10, 0, 0], [0, 0, 10, 0], [1, 1, 1, 7], [5, 5, 0, 0], [0, 0, 0, 10]]

ppm = gimmemotifs_pfm_to_ppm(A)["ppm"]
print([round(v, 5) for v in ppm[0]])
print([round(v, 5) for v in ppm[2]])
print(round(sum(ppm[3]), 6))
```

- **Example output:**

```text
[0.9997, 0.0001, 0.0001, 0.0001]
[0.0001, 0.0001, 0.79978, 0.20002]
1.0
```

## 10) `dnakit.regulation.gimmemotifs_ppm_pseudocount()`

- **Purpose:** Add a pseudocount to a probability matrix and renormalise each row.
- **Method:** Call `gimmemotifs.utils.ppm_pseudocount`: `(p + pseudo) / sum(p + pseudo)` for each row. The default `1e-6` is the value GimmeMotifs adds before comparing motifs so that no probability is exactly zero.
- **API:** `dnakit.regulation.gimmemotifs_ppm_pseudocount(ppm, *, pseudo=1e-06)`.
- **Input:** A probability matrix (rows of four values) and a non-negative `pseudo`.
- **Output:** `ppm` with the pseudocount applied.
- **Example code:**

```python
from dnakit.regulation.gimmemotifs_p28_motif_functions import gimmemotifs_ppm_pseudocount

result = gimmemotifs_ppm_pseudocount([[1, 0, 0, 0], [0.25, 0.25, 0.25, 0.25]])
print([f"{v:.9f}" for v in result["ppm"][0]])
print(result["ppm"][1])
```

- **Example output:**

```text
['0.999997000', '0.000001000', '0.000001000', '0.000001000']
[0.25, 0.25, 0.25, 0.25]
```

## 11) `dnakit.regulation.gimmemotifs_format_motif()`

- **Purpose:** Write a motif as text in a GimmeMotifs export format.
- **Method:** `output_format` selects the GimmeMotifs method: `ppm` (`Motif.to_ppm`, a `>name` line and tab-separated probabilities with `precision` decimals), `pfm` (`Motif.to_pfm`, counts), `meme` (`Motif.to_meme`, a MEME `letter-probability matrix` block), `motevo` (`Motif.to_motevo`, the MotEvo count format) or `transfac` (`Motif.to_transfac`). The text is returned and not written to a file.
- **API:** `dnakit.regulation.gimmemotifs_format_motif(pfm, *, output_format="ppm", motif_id="motif", precision=4)`.
- **Input:** One PFM, the format, a single-line `motif_id` and `precision` from 0 to 12.
- **Output:** `text`, `output_format` and `motif_id`.
- **Example code:**

```python
from dnakit.regulation.gimmemotifs_p28_motif_functions import gimmemotifs_format_motif

A = [[10, 0, 0, 0], [0, 10, 0, 0], [0, 0, 8, 2], [1, 1, 1, 7]]
B = [[0, 10, 0, 0], [0, 0, 10, 0], [1, 1, 1, 7], [5, 5, 0, 0], [0, 0, 0, 10]]

for fmt in ("ppm", "motevo", "meme"):
    print(gimmemotifs_format_motif(A, output_format=fmt, motif_id="mA", precision=2)["text"])
```

- **Example output:**

```text
>mA
1.00	0.00	0.00	0.00
0.00	1.00	0.00	0.00
0.00	0.00	0.80	0.20
0.10	0.10	0.10	0.70
//
NA mA
P0	A	C	G	T
0	10	0	0	0
1	0	10	0	0
2	0	0	8	2
3	1	1	1	7
//
MOTIF mA
BL   MOTIF mA width=0 seqs=0
letter-probability matrix: alength= 4 w= 4 nsites= 10 E= 0
0.9997	0.0001	0.0001	0.0001
0.0001	0.9997	0.0001	0.0001
0.0001	0.0001	0.7998	0.2
0.1	0.1001	0.1001	0.6998
```

## 12) `dnakit.regulation.gimmemotifs_average_motifs()`

- **Purpose:** Average two motifs at a given alignment into a new motif.
- **Method:** Call `Motif.average_motifs`. `position` is the offset of the second motif against the first (-2 places it two positions to the left) and `orientation` is 1 for the second motif as given or -1 for its reverse complement. Positions covered by one motif only are completed with the other motif's frequencies (scaled to the count of the first), or with the uniform background when `include_bg` is true; the counts of the aligned matrices are then added.
- **API:** `dnakit.regulation.gimmemotifs_average_motifs(pfm_a, pfm_b, *, position=0, orientation=1, include_bg=False)`.
- **Input:** Two PFMs, an integer `position`, `orientation` of 1 or -1 and `include_bg`.
- **Output:** `pfm`, `ppm` and `consensus` of the averaged motif.
- **Example code:**

```python
from dnakit.regulation.gimmemotifs_p28_motif_functions import gimmemotifs_average_motifs

A = [[10, 0, 0, 0], [0, 10, 0, 0], [0, 0, 8, 2], [1, 1, 1, 7]]
B = [[0, 10, 0, 0], [0, 0, 10, 0], [1, 1, 1, 7], [5, 5, 0, 0], [0, 0, 0, 10]]

result = gimmemotifs_average_motifs(A, B, position=-1, orientation=-1)
print(result["pfm"])
print(result["consensus"])
```

- **Example output:**

```text
[[20.0, 0.0, 0.0, 0.0], [10.0, 0.0, 5.0, 5.0], [7.0, 11.0, 1.0, 1.0], [0.0, 10.0, 8.0, 2.0], [1.0, 1.0, 11.0, 7.0]]
Awmsk
```

## 13) `dnakit.regulation.gimmemotifs_sample_motif_sequences()`

- **Purpose:** Sample random sequences that follow a motif's probability matrix.
- **Method:** Call `Motif.sample` with a numpy `default_rng(seed)`: each position of each sequence is drawn independently from the probabilities of that position. A position with probability 1 for one base always gives that base.
- **API:** `dnakit.regulation.gimmemotifs_sample_motif_sequences(pfm, n_sequences, *, seed=None)`.
- **Input:** One PFM, the number of sequences (1 to 1,000,000) and an optional integer `seed`.
- **Output:** `sequences`, `n_sequences` and `seed`.
- **Example code:**

```python
from dnakit.regulation.gimmemotifs_p28_motif_functions import gimmemotifs_sample_motif_sequences

A = [[10, 0, 0, 0], [0, 10, 0, 0], [0, 0, 8, 2], [1, 1, 1, 7]]
B = [[0, 10, 0, 0], [0, 0, 10, 0], [1, 1, 1, 7], [5, 5, 0, 0], [0, 0, 0, 10]]

result = gimmemotifs_sample_motif_sequences(A, 6, seed=4)
print(result["sequences"])
print(all(s.startswith("AC") for s in result["sequences"]))
```

- **Example output:**

```text
['ACGC', 'ACGT', 'ACTT', 'ACGT', 'ACTT', 'ACTC']
True
```

## 14) `dnakit.regulation.gimmemotifs_pwm_scan()`

- **Purpose:** Scan FASTA sequences with one motif and report scores, positions, full hits or a GFF file.
- **Method:** The cutoff is a fraction (0 to 1) of the way from the motif's minimum to its maximum PWM score. `output_format="scores"` calls `Motif.pwm_scan_score` (the score of each reported match; defaults `cutoff=0`, `nreport=1`, so the best match of every sequence), `"positions"` calls `Motif.scan` (match positions), `"hits"` calls `Motif.scan_all` (position, score and strand of each match) and `"gff"` calls `Motif.pwm_scan_to_gff`, which also writes the hits to a new GFF file (`gff_path`) with the motif name `motif_id`. The defaults for the last three are `cutoff=0.9`, `nreport=50`. `scan_rc` includes the reverse strand. Sequence names are the complete FASTA header lines.
- **API:** `dnakit.regulation.gimmemotifs_pwm_scan(sequences_fasta, pfm, *, output_format="scores", cutoff=None, nreport=None, scan_rc=True, gff_path=None, motif_id="motif")`.
- **Input:** A FASTA file of at most 100,000 sequences, one PFM, the options above.
- **Output:** `matches` per sequence, the parameters used and, for `gff`, the file in `outputs`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.regulation.gimmemotifs_p28_motif_functions import gimmemotifs_pwm_scan

A = [[10, 0, 0, 0], [0, 10, 0, 0], [0, 0, 8, 2], [1, 1, 1, 7]]
B = [[0, 10, 0, 0], [0, 0, 10, 0], [1, 1, 1, 7], [5, 5, 0, 0], [0, 0, 0, 10]]

best = gimmemotifs_pwm_scan("seqs.fa", A)
print({name: [round(v, 3) for v in values] for name, values in best["matches"].items()})
hits = gimmemotifs_pwm_scan("seqs.fa", A, output_format="hits", cutoff=0.95)
print(hits["matches"]["s1"])
gff = Path(tempfile.mkdtemp()) / "hits.gff"
gimmemotifs_pwm_scan("seqs.fa", A, output_format="gff", cutoff=0.95, gff_path=gff, motif_id="mA")
print(gff.read_text().splitlines()[0])
```

- **Example output:**

```text
{'s1': [4.976], 's2': [4.976], 's3': [-8.634]}
[{'position': 2, 'score': 4.975905319526848, 'strand': 1}, {'position': 8, 'score': 4.975905319526848, 'strand': 1}, {'position': 2, 'score': 4.975905319526848, 'strand': -1}, {'position': 8, 'score': 4.975905319526848, 'strand': -1}]
s1	pfmscan	misc_feature	2	6	4.976	+	.	motif_name "mA" ; motif_instance "ACGT"
```

## 15) `dnakit.regulation.gimmemotifs_consensus_scan()`

- **Purpose:** Scan FASTA sequences for a motif's IUPAC consensus as an exact pattern.
- **Method:** Call `Motif.consensus_scan`: every symbol of the consensus (`Motif.to_consensusv2`) is expanded to a character class and the resulting regular expression is searched in each sequence on the forward strand. The middle coordinate of each match is reported (0-based, so it can end in .5).
- **API:** `dnakit.regulation.gimmemotifs_consensus_scan(sequences_fasta, pfm)`.
- **Input:** A FASTA file of at most 100,000 sequences and one PFM.
- **Output:** `matches` (match middles per sequence) and the `consensus` that was searched.
- **Example code:**

```python
from dnakit.regulation.gimmemotifs_p28_motif_functions import gimmemotifs_consensus_scan

A = [[10, 0, 0, 0], [0, 10, 0, 0], [0, 0, 8, 2], [1, 1, 1, 7]]
B = [[0, 10, 0, 0], [0, 0, 10, 0], [1, 1, 1, 7], [5, 5, 0, 0], [0, 0, 0, 10]]

result = gimmemotifs_consensus_scan("seqs.fa", A)
print(result["consensus"], result["matches"])
```

- **Example output:**

```text
ACGT {'s1': [4.0, 10.0], 's2': [6.0], 's3': []}
```
