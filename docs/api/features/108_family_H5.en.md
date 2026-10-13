# 108 Family-variant review (H5): microbiome, regulation, epigenome, variants, predictions

This page records how review group H5 checked the `family-variant` ledger rows: 122 rows in 17 packages. Each row was tested for being a true alias, re-export or branch reachable through a parameter of an existing wrapper. After the check, 40 rows became `duplicate-of` (same object or same code), 51 stay `family-variant` (the parameter exists and the result equals the direct upstream call), and 31 became `wrapped` through 26 new functions. The new functions compute locally and use no network; inputs are sketches, number lists or local files. Optional extras: `bio-sourmash`, `bio-gimmemotifs`, `bio-tangermeme`, `bio-macs3`, `bio-cooltools`, `bio-pairtools`, `bio-deeptools`, `bio-methylpy`, `bio-varcode` and `bio-reads` (pysam). The modules are `dnakit.microbiome`, `dnakit.regulation`, `dnakit.regulation`, `dnakit.epigenome`, `dnakit.epigenome`, `dnakit.epigenome`, `dnakit.epigenome`, `dnakit.epigenome`, `dnakit.variants` and `dnakit.variants`; the package exports are wired by the lead, so the examples import from the modules.

## Review result

| Package | Before | duplicate-of | family-variant (verified) | wrapped |
| --- | --- | --- | --- | --- |
| sourmash | 36 | 19 | 5 | 12 |
| gimmemotifs | 16 | 0 | 13 | 3 |
| MACS3 | 14 | 7 | 2 | 5 |
| grelu | 10 | 2 | 8 | 0 |
| cooltools | 8 | 1 | 5 | 2 |
| alphagenome | 7 | 0 | 7 | 0 |
| pairtools | 5 | 2 | 1 | 2 |
| tangermeme | 5 | 2 | 1 | 2 |
| deeptools | 5 | 4 | 0 | 1 |
| varcode | 4 | 1 | 1 | 2 |
| bed_reader | 3 | 1 | 2 | 0 |
| methylpy | 2 | 0 | 1 | 1 |
| hgvs | 2 | 0 | 2 | 0 |
| vcfpy | 2 | 0 | 2 | 0 |
| pysam | 1 | 0 | 0 | 1 |
| dendropy | 1 | 1 | 0 | 0 |
| enformer_pytorch | 1 | 0 | 1 | 0 |
| **Total** | 122 | 40 | 51 | 31 |

## 1) `dnakit.microbiome.write_sourmash_signature_archive()`

- **Purpose:** Write MinHash sketches to a sourmash zip collection or SQLite index.
- **Method:** The output name selects the format (.zip and .sig.zip give zip, .sqldb gives SQLite); sourmash.save_load.SaveSignaturesToLocation does the writing, and abundances are kept when a sketch has them. The existing write_sourmash_signatures writes only .sig and .sig.gz.
- **API:** `dnakit.microbiome.write_sourmash_signature_archive(sketches, output, *, verbose=True)`.
- **Input:** sketches is a list of sketch mappings (ksize, scaled, hashes, optional id and abundances); output is a path that does not exist yet.
- **Output:** format (zip or sqlite), signatures (count), output and provenance; the file reads back with read_sourmash_signatures.
- **Example code:**

```python
a = {"id": "a", "ksize": 21, "scaled": 1, "hashes": [1, 2, 3, 4, 5, 6]}
b = {"id": "b", "ksize": 21, "scaled": 1, "hashes": [4, 5, 6, 7, 8, 9]}
c = {"id": "c", "ksize": 21, "scaled": 1, "hashes": [20, 21, 22]}
import tempfile
from pathlib import Path

from dnakit.microbiome import read_sourmash_signatures
from dnakit.microbiome import write_sourmash_signature_archive

with tempfile.TemporaryDirectory() as folder:
    for name in ("sketches.sig.zip", "sketches.sqldb"):
        result = write_sourmash_signature_archive([a, b], Path(folder) / name, verbose=False)
        back = read_sourmash_signatures(result["output"], verbose=False)
        print(result["format"], result["signatures"], [s["id"] for s in back["sketches"]])
        print(back["sketches"][1]["hashes"])
```

- **Example output:**

```text
zip 2 ['a', 'b']
[4, 5, 6, 7, 8, 9]
sqlite 2 ['a', 'b']
[4, 5, 6, 7, 8, 9]
```

## 2) `dnakit.microbiome.parse_sourmash_signature_json()`

- **Purpose:** Parse sourmash signature JSON text into sketches without reading a file.
- **Method:** With require_single=False it calls sourmash.signature.load_signatures_from_json; with True it calls load_one_signature_from_json, which needs exactly one matching signature. ksize selects a k-mer size. The output is sorted by id, ksize and md5sum.
- **API:** `dnakit.microbiome.parse_sourmash_signature_json(json_text, *, ksize=None, require_single=False, max_signatures=1000, verbose=True)`.
- **Input:** json_text is a JSON string; ksize is optional; require_single is a boolean; max_signatures is a limit.
- **Output:** signatures (count), sketches (id, ksize, scaled, hash_count, hashes, md5sum, and abundances when tracked) and provenance.
- **Example code:**

```python
a = {"id": "a", "ksize": 21, "scaled": 1, "hashes": [1, 2, 3, 4, 5, 6]}
b = {"id": "b", "ksize": 21, "scaled": 1, "hashes": [4, 5, 6, 7, 8, 9]}
c = {"id": "c", "ksize": 21, "scaled": 1, "hashes": [20, 21, 22]}
from dnakit.microbiome import write_sourmash_signatures
from dnakit.microbiome import parse_sourmash_signature_json
import tempfile
from pathlib import Path

k31 = {"id": "k31", "ksize": 31, "scaled": 1, "hashes": [1, 2, 3]}
with tempfile.TemporaryDirectory() as folder:
    path = Path(folder) / "all.sig"
    write_sourmash_signatures([a, b, k31], path, verbose=False)
    text = path.read_text()
every = parse_sourmash_signature_json(text, verbose=False)
print(every["signatures"], [s["id"] for s in every["sketches"]])
one = parse_sourmash_signature_json(text, ksize=31, require_single=True, verbose=False)
print(one["sketches"][0]["id"], one["sketches"][0]["ksize"], one["sketches"][0]["hashes"])
```

- **Example output:**

```text
3 ['a', 'b', 'k31']
k31 31 [1, 2, 3]
```

## 3) `dnakit.microbiome.select_sourmash_signatures()`

- **Purpose:** Select DNA signatures from sourmash files by k-mer size, name and md5.
- **Method:** expect='many' calls sourmash_args.load_many_signatures; 'one' calls load_one_signature for each source and needs exactly one; 'query' calls load_query_signature and needs a single source with one signature left after selection. name_pattern is a case-insensitive regular expression searched in the name, file name and md5.
- **API:** `dnakit.microbiome.select_sourmash_signatures(sources, *, ksize=None, name_pattern=None, md5_prefix=None, expect='many', max_signatures=1000, verbose=True)`.
- **Input:** sources is a list of file paths; ksize, name_pattern and md5_prefix are optional filters; expect is many, one or query.
- **Output:** signatures (count), sketches (with source and md5sum, sorted by source, id and md5sum) and provenance.
- **Example code:**

```python
a = {"id": "a", "ksize": 21, "scaled": 1, "hashes": [1, 2, 3, 4, 5, 6]}
b = {"id": "b", "ksize": 21, "scaled": 1, "hashes": [4, 5, 6, 7, 8, 9]}
c = {"id": "c", "ksize": 21, "scaled": 1, "hashes": [20, 21, 22]}
import tempfile
from pathlib import Path

from dnakit.microbiome import write_sourmash_signatures
from dnakit.microbiome import select_sourmash_signatures

k31 = {"id": "k31", "ksize": 31, "scaled": 1, "hashes": [1, 2, 3]}
with tempfile.TemporaryDirectory() as folder:
    path = Path(folder) / "all.sig"
    write_sourmash_signatures([a, b, k31], path, verbose=False)
    def ids(**options):
        result = select_sourmash_signatures([path], verbose=False, **options)
        return [item["id"] for item in result["sketches"]]
    print(ids())
    print(ids(ksize=21))
    print(ids(name_pattern="^B"))
    print(ids(expect="one", ksize=31))
    print(ids(expect="query", ksize=31))
```

- **Example output:**

```text
['a', 'b', 'k31']
['a', 'b']
['b']
['k31']
['k31']
```

## 4) `dnakit.microbiome.compare_minhash_abundance_matrix()`

- **Purpose:** Compare sketches pairwise, with optional abundance weighting and several processes.
- **Method:** Calls sourmash.compare.compare_all_pairs: n_jobs=1 uses compare_serial, otherwise compare_parallel. Sketches with abundances give the angular (cosine) similarity unless ignore_abundance is true; otherwise the Jaccard index is used. return_ani=True returns the estimated ANI.
- **API:** `dnakit.microbiome.compare_minhash_abundance_matrix(sketches, *, ignore_abundance=False, return_ani=False, n_jobs=1, verbose=True)`.
- **Input:** sketches (at least 2, same ksize and scaled; abundances on all or none); ignore_abundance, return_ani, and n_jobs (1 to 16).
- **Output:** measure (jaccard, angular_similarity or ani), abundance_weighted, ids, matrix and provenance.
- **Example code:**

```python
from dnakit.microbiome import compare_minhash_abundance_matrix

hashes = [1, 2, 3, 4]
a = {"id": "a", "ksize": 21, "scaled": 1, "hashes": hashes, "abundances": [1, 1, 1, 1]}
b = {"id": "b", "ksize": 21, "scaled": 1, "hashes": hashes, "abundances": [1, 1, 2, 2]}
result = compare_minhash_abundance_matrix([a, b], verbose=False)
print(result["measure"], round(result["matrix"][0][1], 4))
flat = compare_minhash_abundance_matrix([a, b], ignore_abundance=True, verbose=False)
print(flat["measure"], flat["matrix"][0][1])
parallel = compare_minhash_abundance_matrix([a, b], n_jobs=2, verbose=False)
print(parallel["matrix"] == result["matrix"])
```

- **Example output:**

```text
angular_similarity 0.7952
jaccard 1.0
True
```

## 5) `dnakit.microbiome.search_minhash_database_ranked()`

- **Purpose:** Rank database sketches by similarity to a query sketch, as the sourmash search command does.
- **Method:** For mode jaccard, containment and max_containment it calls sourmash.search.search_databases_with_flat_query; for mode='abundance' it calls search_databases_with_abund_query, which needs abundances on the query and every database sketch and scores by angular similarity. Matches below threshold are dropped and the rest are sorted by score, then id.
- **API:** `dnakit.microbiome.search_minhash_database_ranked(query, database, *, mode='jaccard', threshold=0.08, verbose=True)`.
- **Input:** query and database (sketch mappings); mode; threshold (0 to 1, default 0.08).
- **Output:** query_id, mode, match_count, matches (id, md5, score) and provenance.
- **Example code:**

```python
a = {"id": "a", "ksize": 21, "scaled": 1, "hashes": list(range(1, 2001))}
b = {"id": "b", "ksize": 21, "scaled": 1, "hashes": list(range(1001, 3001))}
c = {"id": "c", "ksize": 21, "scaled": 1, "hashes": list(range(5001, 5501))}
from dnakit.microbiome import search_minhash_database_ranked

for mode in ("jaccard", "containment", "max_containment"):
    result = search_minhash_database_ranked(a, [a, b, c], mode=mode, verbose=False)
    print(mode, [(m["id"], round(m["score"], 4)) for m in result["matches"]])
```

- **Example output:**

```text
jaccard [('a', 1.0), ('b', 0.3333)]
containment [('a', 1.0), ('b', 0.5)]
max_containment [('a', 1.0), ('b', 0.5)]
```

## 6) `dnakit.microbiome.search_sourmash_sbt_index()`

- **Purpose:** Search a saved sourmash sequence bloom tree (SBT) index with a query sketch.
- **Method:** Loads the index with sourmash.sbtmh.load_sbt_index and calls sourmash.sbtmh.search_sbt_index. The score is the Jaccard index, so the matches equal those of search_minhash_database_ranked(mode='jaccard') on the same sketches.
- **API:** `dnakit.microbiome.search_sourmash_sbt_index(index, query, *, threshold=0.08, verbose=True)`.
- **Input:** index is a .sbt.zip or .sbt.json file (as written by sourmash index); query is a sketch with the ksize and scaled of the index; threshold.
- **Output:** query_id, match_count, matches (id, md5, jaccard, sorted by jaccard then id) and provenance.
- **Example code:**

```python
a = {"id": "a", "ksize": 21, "scaled": 1, "hashes": list(range(1, 2001))}
b = {"id": "b", "ksize": 21, "scaled": 1, "hashes": list(range(1001, 3001))}
c = {"id": "c", "ksize": 21, "scaled": 1, "hashes": list(range(5001, 5501))}
import tempfile
from pathlib import Path

import sourmash
from sourmash.sbtmh import create_sbt_index

from dnakit.microbiome import search_sourmash_sbt_index

tree = create_sbt_index()
for sketch in (a, b, c):
    minhash = sourmash.MinHash(n=0, ksize=21, scaled=1)
    minhash.add_many(sketch["hashes"])
    tree.insert(sourmash.SourmashSignature(minhash, name=sketch["id"]))
with tempfile.TemporaryDirectory() as folder:
    path = Path(folder) / "index.sbt.zip"
    tree.save(str(path))
    result = search_sourmash_sbt_index(path, a, threshold=0.1, verbose=False)
print([(m["id"], round(m["jaccard"], 4)) for m in result["matches"]])
```

- **Example output:**

```text
[('a', 1.0), ('b', 0.3333)]
```

## 7) `dnakit.microbiome.downsample_minhash_num_sketch()`

- **Purpose:** Downsample a num-based MinHash sketch to fewer hashes.
- **Method:** Calls sourmash.minhash.flatten_and_downsample_num, which keeps the smallest num hashes. The new num cannot exceed the sketch's own num. The existing downsample_minhash_sketch handles only scaled sketches.
- **API:** `dnakit.microbiome.downsample_minhash_num_sketch(sketch, *, num, verbose=True)`.
- **Input:** sketch has ksize, num and hashes (optional id); num is the new hash count.
- **Output:** input_num, num, hash_count, hashes and provenance.
- **Example code:**

```python
from dnakit.microbiome import downsample_minhash_num_sketch

sketch = {"id": "n", "ksize": 21, "num": 6, "hashes": [34, 3, 21, 5, 13, 8]}
result = downsample_minhash_num_sketch(sketch, num=3, verbose=False)
print(result["input_num"], result["num"], result["hashes"])
```

- **Example output:**

```text
6 3 [3, 5, 8]
```

## 8) `dnakit.regulation.normalize_peak_widths()`

- **Purpose:** Resize peaks to a fixed window around their midpoint, honouring the strand.
- **Method:** Calls gimmemotifs.preprocessing.peakfile2bedfile. The midpoint is (start + end) // 2; the + strand (the default) becomes midpoint ± window // 2; for other strands the upstream writes the two coordinates the other way round, so start is larger than end and the row carries start_after_end=True. Rows whose start would be below 0 are dropped, not clamped. Missing name, score and strand become an empty string, 0 and +.
- **API:** `dnakit.regulation.normalize_peak_widths(path, *, window=200, verbose=True)`.
- **Input:** path is a tab-separated BED-like file (at least 3 columns); window is the width.
- **Output:** intervals (chrom, start, end, name, score, strand, start_after_end), input_rows, dropped_rows and provenance.

**Example file:** `peaks.bed`

```text
chr1	100	300	p1	5	+
chr1	1000	1100	p2	7	-
chr2	10	50
chr3	500	700	p4
```

- **Example code:**

```python
from dnakit.regulation import normalize_peak_widths

result = normalize_peak_widths("peaks.bed", window=100, verbose=False)
for item in result["intervals"]:
    print(item["chrom"], item["start"], item["end"], item["name"], item["strand"], item["start_after_end"])
print(result["input_rows"], result["dropped_rows"])
```

- **Example output:**

```text
chr1 150 250 p1 + False
chr1 1100 1000 p2 - True
chr3 550 650 p4 + False
4 1
```

## 9) `dnakit.regulation.count_motif_matches_at_fpr()`

- **Purpose:** Count foreground and background sequences scoring above the score cutoff at a false-positive rate.
- **Method:** Calls gimmemotifs.rocmetrics.matches_at_fpr: the cutoff is the (100 - fpr * 100) percentile of the background scores, and the result is the number of foreground and of background scores at or above it. motif_score_metric returns one number only and cannot return these two counts.
- **API:** `dnakit.regulation.count_motif_matches_at_fpr(foreground, background, *, fpr=0.01, verbose=True)`.
- **Input:** foreground and background are lists of numbers; fpr is in (0, 1].
- **Output:** foreground_matches, background_matches, their fractions, the two totals, fpr and provenance.
- **Example code:**

```python
from dnakit.regulation import count_motif_matches_at_fpr

background = [float(value) for value in range(100)]
foreground = [95.0, 99.0, 50.0, 80.0, 98.5]
result = count_motif_matches_at_fpr(foreground, background, fpr=0.1, verbose=False)
print(result["foreground_matches"], result["background_matches"], result["foreground_fraction"])
```

- **Example output:**

```text
3 10 0.6
```

## 10) `dnakit.regulation.calculate_partial_motif_roc_auc()`

- **Purpose:** Compute the area under the motif-score ROC curve up to a false-positive-rate limit.
- **Method:** Calls gimmemotifs.rocmetrics.roc_auc_xlim with the foreground and background score lists (the upstream parameters are named x_bla and y_bla but hold the scores). The area is the integral of the true-positive rate over false-positive rates 0 to fpr_limit, so it is at most fpr_limit. motif_score_metric passes ROC curve points to this function, which makes gimmemotifs fail.
- **API:** `dnakit.regulation.calculate_partial_motif_roc_auc(foreground, background, *, fpr_limit=0.1, verbose=True)`.
- **Input:** foreground and background are lists of finite numbers; fpr_limit is in (0, 1].
- **Output:** area, area_fraction_of_limit, the two counts, fpr_limit and provenance.
- **Example code:**

```python
from dnakit.regulation import calculate_partial_motif_roc_auc

foreground = [3.0, 4.0, 5.0]
background = [1.0, 2.0, 3.5, 4.5]
result = calculate_partial_motif_roc_auc(foreground, background, fpr_limit=0.5, verbose=False)
print(result["area"], result["area_fraction_of_limit"])
```

- **Example output:**

```text
0.25 0.5
```

## 11) `dnakit.regulation.design_beam_substitution()`

- **Purpose:** Beam search over motif substitutions that move the model output towards a target.
- **Method:** Calls tangermeme.design.beam_substitution. Each round tiles every motif (and its reverse complement when reverse_complement=True) at every position of every sequence in the beam, scores the full sequences by mean squared error to target, and keeps the beam_size best distinct ones. It stops when the best loss improves by no more than tol or after max_iter rounds. beam_size=1 equals design_greedy_substitution.
- **API:** `dnakit.regulation.design_beam_substitution(model, sequence, target, *, motifs=None, reverse_complement=True, beam_size=4, n_best=1, max_iter=10, tol=0.001, verbose=True)`.
- **Input:** model is a PyTorch model in eval mode; sequence is one DNA string; target has one number per output; motifs is optional; beam_size, n_best (at most beam_size), max_iter, tol.
- **Output:** designed_sequences (lowest loss first), designed_sequence (the best), input_sequence, beam_size and provenance.
- **Example code:**

```python
import torch

from dnakit.regulation import design_greedy_substitution
from dnakit.regulation import design_beam_substitution

length = 20
weights = torch.tensor(
    [[((base * 7 + pos * 3) % 11) / 4.0 - 1.0 for pos in range(length)] for base in range(4)]
)


class Linear(torch.nn.Module):
    def __init__(self):
        super().__init__()
        self.register_buffer("weights", weights)

    def forward(self, x, *args):
        return (x * self.weights).sum(dim=(1, 2)).unsqueeze(-1)


model = Linear().eval()
sequence = "ACGTTGATAAGCCGTAGGAT"
beam = design_beam_substitution(
    model, sequence, [5.0], motifs=["GGG", "AC"], beam_size=4, n_best=2, max_iter=3, verbose=False
)
print(beam["designed_sequences"])
single = design_beam_substitution(model, sequence, [5.0], motifs=["GGG"], beam_size=1, max_iter=2, verbose=False)
greedy = design_greedy_substitution(model, sequence, [5.0], motifs=["GGG"], max_iter=2, verbose=False)
print(single["designed_sequence"] == greedy["designed_sequence"])
```

- **Example output:**

```text
['ACCCCGATAAGCCGTACCCT', 'ACCCCGATACCCCGTACCCT']
True
```

## 12) `dnakit.regulation.insert_motif()`

- **Purpose:** Insert a motif into each sequence and keep the longer result.
- **Method:** Calls tangermeme.ersatz.insert. The motif is inserted before base start (0-based; the centre of each sequence when start is None), so every output is len(motif) bases longer than its input and no base is lost. place_motif(mode='insert') trims the result back to the input length from the right, which drops the last bases.
- **API:** `dnakit.regulation.insert_motif(sequences, motif, *, start=None, verbose=True)`.
- **Input:** sequences is a DNA string or a mapping of names to strings (same length); motif is DNA; start is optional.
- **Output:** sequences (name, sequence), input_length, output_length, start and provenance.
- **Example code:**

```python
from dnakit.regulation import place_motif
from dnakit.regulation import insert_motif

sequence = "ACGTTGATAAGCCGTAGGAT"
longer = insert_motif({"s": sequence}, "GG", start=3, verbose=False)
print(longer["sequences"][0]["sequence"], longer["input_length"], longer["output_length"])
trimmed = place_motif({"s": sequence}, "GG", start=3, mode="insert", verbose=False)
print(trimmed["sequences"][0]["sequence"])
```

- **Example output:**

```text
ACGGGTTGATAAGCCGTAGGAT 20 22
ACGGGTTGATAAGCCGTAGG
```

## 13) `dnakit.epigenome.build_counted_fragment_pileup()`

- **Purpose:** Pile up paired-end fragments that each stand for several identical reads.
- **Method:** Calls MACS3.Signal.PileupV2.pileup_from_LRC_as_list: fragment i covers [left, right) and adds counts[i] to each base. left_sorted=True promises the fragments are sorted by left position (it is checked) so MACS3 skips its sort. The layout matches build_paired_fragment_pileup: positions are segment ends and values the pileup of the segment ending there.
- **API:** `dnakit.epigenome.build_counted_fragment_pileup(left_positions, right_positions, counts, *, scale_factor=1.0, baseline_value=0.0, left_sorted=False, verbose=True)`.
- **Input:** left_positions, right_positions and counts (0 to 65535) are integer lists of equal length; scale_factor, baseline_value, left_sorted.
- **Output:** fragments, total_count, positions, values and provenance.
- **Example code:**

```python
from dnakit.epigenome import build_counted_fragment_pileup

left = [10, 20, 40, 100]
right = [50, 60, 90, 150]
counts = [2, 1, 3, 1]
result = build_counted_fragment_pileup(left, right, counts, verbose=False)
print(result["positions"])
print(result["values"])
```

- **Example output:**

```text
[10, 20, 40, 50, 60, 90, 100, 150]
[0.0, 2.0, 3.0, 6.0, 4.0, 3.0, 0.0, 1.0]
```

## 14) `dnakit.epigenome.build_fragment_end_pileup()`

- **Purpose:** Pile up windows centred on both ends of count-weighted fragments.
- **Method:** Calls MACS3.Signal.PileupV2.pileup_from_LRC_centers_as_list. Every fragment contributes two windows of width window with weight counts[i]: one starting at left - window // 2 and one at right - window // 2. The fragment body is not covered, so the pileup shows where fragments start and end.
- **API:** `dnakit.epigenome.build_fragment_end_pileup(left_positions, right_positions, counts, *, window, scale_factor=1.0, baseline_value=0.0, verbose=True)`.
- **Input:** left_positions, right_positions and counts are integer lists of equal length; window is the width; scale_factor, baseline_value.
- **Output:** fragments, window, total_count, positions, values and provenance.
- **Example code:**

```python
from dnakit.epigenome import build_fragment_end_pileup

result = build_fragment_end_pileup([100], [200], [2], window=20, verbose=False)
print(result["positions"])
print(result["values"])
```

- **Example output:**

```text
[90, 110, 190, 210]
[0.0, 2.0, 0.0, 2.0]
```

## 15) `dnakit.epigenome.build_hmmratac_fragment_pileup()`

- **Purpose:** Pile up paired-end fragments weighted by their length, as HMMRATAC does.
- **Method:** Calls MACS3.Signal.PileupV2.pileup_from_LR_hmmratac. length_weights maps a fragment length (right - left, written as a decimal string) to the weight added to every base of that fragment; every length present needs a weight, and the weights can come from hmmratac_fragment_weights.
- **API:** `dnakit.epigenome.build_hmmratac_fragment_pileup(left_positions, right_positions, length_weights, *, verbose=True)`.
- **Input:** left_positions and right_positions are integer lists of equal length; length_weights maps length strings to weights.
- **Output:** fragments, positions, values and provenance.
- **Example code:**

```python
from dnakit.epigenome import build_hmmratac_fragment_pileup

left = [10, 20, 40, 100]
right = [50, 60, 90, 150]
result = build_hmmratac_fragment_pileup(left, right, {"40": 0.5, "50": 2.0}, verbose=False)
print(result["positions"])
print(result["values"])
```

- **Example output:**

```text
[10, 20, 40, 50, 60, 90, 100, 150]
[0.0, 0.5, 1.0, 3.0, 2.5, 2.0, 0.0, 2.0]
```

## 16) `dnakit.epigenome.calculate_allele_specific_model_likelihood()`

- **Purpose:** Log-likelihood and BIC of the heterozygous model with an allele-specific ratio.
- **Method:** Calls MACS3.Signal.VariantStat.CalModel_Heter_AS. The treatment reads (_t) are fitted with a free allele ratio capped at max_allowed_ratio; the control reads (_c) are fitted without allele preference. Each argument lists the base qualities of the reads supporting the top allele (top1) or the second allele (top2). calculate_allele_model_likelihood covers only the models without a free ratio.
- **API:** `dnakit.epigenome.calculate_allele_specific_model_likelihood(top1_bq_t, top1_bq_c, top2_bq_t, top2_bq_c, *, max_allowed_ratio=0.99, verbose=True)`.
- **Input:** top1_bq_t, top1_bq_c, top2_bq_t and top2_bq_c are integer lists of equal length; max_allowed_ratio is in [0.5, 1].
- **Output:** reads, max_allowed_ratio, log_likelihood, bic and provenance.
- **Example code:**

```python
from dnakit.epigenome import calculate_allele_specific_model_likelihood

reads = ([30, 32, 35, 28, 30, 33], [2, 3, 1, 2, 3, 2], [0] * 6, [0] * 6)
result = calculate_allele_specific_model_likelihood(*reads, verbose=False)
print(round(result["log_likelihood"], 4), round(result["bic"], 4))
```

- **Example output:**

```text
-11.2048 29.8643
```

## 17) `dnakit.epigenome.calculate_allele_specific_genotype_quality()`

- **Purpose:** Phred-like score for the allele-specific model against a second model.
- **Method:** Calls MACS3.Signal.VariantStat.calculate_GQ_heterASsig. With L1 = 1 and L2 = min(1, exp(lnL2 - lnL1)) (at least 1e-110) the score is int(-10 * log10(L2 / (L1 + L2))), or 255 when the ratio is not above 1e-110. It is 3 when the two log-likelihoods are equal and grows as model 1 fits better. The three-model calculate_genotype_quality uses a different formula.
- **API:** `dnakit.epigenome.calculate_allele_specific_genotype_quality(log_likelihood_1, log_likelihood_2, *, verbose=True)`.
- **Input:** log_likelihood_1 and log_likelihood_2 are natural-log likelihoods.
- **Output:** score (an integer) and provenance.
- **Example code:**

```python
from dnakit.epigenome import calculate_allele_specific_genotype_quality

for first, second in ((-10.0, -15.0), (-10.0, -10.0), (-15.0, -10.0), (-10.0, -300.0)):
    score = calculate_allele_specific_genotype_quality(first, second, verbose=False)["score"]
    print(first, second, score)
```

- **Example output:**

```text
-10.0 -15.0 21
-10.0 -10.0 3
-15.0 -10.0 3
-10.0 -300.0 255
```

## 18) `dnakit.epigenome.pairwise_region_sums()`

- **Purpose:** Diagonal sums within each region pair plus one total row for each pair of chromosomes.
- **Method:** Calls cooltools.sandbox.obs_over_exp_cooler.sum_pairwise on the chromosome view of the cooler. Rows with region1 equal to region2 and dist of 0 or more equal the rows of pairwise_diagonal_sums; every pair of different regions also has one row with dist = -1 holding the block total (count.sum) and the number of valid pixels (n_valid); n_total is empty there.
- **API:** `dnakit.epigenome.pairwise_region_sums(cool_path, *, weight='weight', ignore_diags=2, verbose=True)`.
- **Input:** cool_path is a balanced cooler; weight is the weight column; ignore_diags is the number of ignored diagonals.
- **Output:** rows (region1, region2, dist, n_total, n_valid, count.sum), row_count, block_row_count and provenance.
- **Example code:**

```python
import tempfile
from pathlib import Path

import cooler
import numpy as np
import pandas as pd

from dnakit.epigenome import pairwise_region_sums

rng = np.random.default_rng(0)
bins = cooler.util.binnify(pd.Series({"chr1": 1500, "chr2": 2000}, name="length"), 100)
rows = []
for i in range(len(bins)):
    for j in range(i, len(bins)):
        same = bins.chrom[i] == bins.chrom[j]
        count = rng.poisson(6.0 / (1 + abs(i - j)) if same else 0.5) + (3 if i == j else 0)
        if count > 0:
            rows.append((i, j, int(count)))
pixels = pd.DataFrame(rows, columns=["bin1_id", "bin2_id", "count"])
with tempfile.TemporaryDirectory() as folder:
    path = Path(folder) / "toy.cool"
    cooler.create_cooler(str(path), bins, pixels, dtypes={"count": "int"}, ensure_sorted=True)
    cooler.balance_cooler(cooler.Cooler(str(path)), store=True)
    result = pairwise_region_sums(path, verbose=False)
block = [row for row in result["rows"] if row["dist"] == -1][0]
print(result["row_count"], result["block_row_count"])
print(block["region1"], block["region2"], block["count.sum"], block["n_valid"])
```

- **Example output:**

```text
36 1
chr1 chr2 144.0 300.0
```

## 19) `dnakit.epigenome.compute_region_contact_scaling()`

- **Purpose:** Raw contact counts and possible-pair areas per log-spaced distance bin for one region pair.
- **Method:** Calls cooltools.sandbox.pairs_scaling_functions.compute_scaling. Contact i joins positions1[i] and positions2[i] on one chromosome. Only contacts with positions1 inside region1 and positions2 inside region2 (both [start, end); region2 defaults to region1) are counted, in the distance bin of positions2 - positions1, for distances from min_distance up to but not including max_distance. The bin area is the integral of the number of possible position pairs, and frequency = observed / area. Unlike compute_contact_scaling the input is positions, not a pairs file, and no chromosome sizes are needed.
- **API:** `dnakit.epigenome.compute_region_contact_scaling(positions1, positions2, *, region1, region2=None, min_distance=10, max_distance=10000000, n_bins=50, verbose=True)`.
- **Input:** positions1 and positions2 are integer lists of equal length; region1 (and optional region2) is [start, end]; min_distance, max_distance, n_bins.
- **Output:** bins (min_dist, max_dist, observed, area, frequency), contacts_counted and provenance.
- **Example code:**

```python
from dnakit.epigenome import compute_region_contact_scaling

result = compute_region_contact_scaling(
    [100, 150, 200, 1000, 5000, 200],
    [130, 400, 1000, 1500, 9000, 260],
    region1=[0, 10000],
    min_distance=10,
    max_distance=10000,
    n_bins=5,
    verbose=False,
)
for item in result["bins"]:
    print(item["min_dist"], item["max_dist"], item["observed"], round(item["area"]))
print(result["contacts_counted"])
```

- **Example output:**

```text
10 56 1 458482
56 316 2 2551640
316 1778 2 13089286
1778 10000 1 33800642
6
```

## 20) `dnakit.epigenome.deduplicate_pairs_online()`

- **Purpose:** Mark PCR duplicate pairs with pairtools' Cython online algorithm.
- **Method:** Calls pairtools.lib.dedup.streaming_dedup_cython. A pair is a duplicate when both sides lie within max_mismatch bp of a pair that was kept earlier (method='max' uses the larger and 'sum' the summed difference of the two positions). Duplicates are not chained: deduplicate_pairs joins pairs that are close to a duplicate into one cluster, so on a run of pairs spaced 2 bp apart it keeps fewer pairs. With mark_duplicates=True the pair type of duplicates written to duplicates_path becomes DD.
- **API:** `dnakit.epigenome.deduplicate_pairs_online(pairs_path, *, output_path=None, duplicates_path=None, unmapped_path=None, max_mismatch=3, method='max', mark_duplicates=False, verbose=True)`.
- **Input:** pairs_path is a sorted .pairs file; output_path, duplicates_path and unmapped_path are optional paths that do not exist yet; max_mismatch, method, mark_duplicates.
- **Output:** total, total_mapped, total_nodups, total_dups, duplicate_fraction, unique_read_ids, duplicate_read_ids (each at most 10000) and provenance.

**Example file:** `dup.pairs`

```text
## pairs format v1.0.0
#shape: upper triangle
#chromsize: chr1 10000
#chromsize: chr2 10000
#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type
r1	chr1	100	chr2	500	+	-	UU
r5	chr1	100	chr2	502	+	-	UU
r2	chr1	102	chr2	500	+	-	UU
r3	chr1	104	chr2	500	+	-	UU
r4	chr1	106	chr2	500	+	-	UU
r6	chr1	3000	chr1	4000	+	-	UU
r8	chr1	3000	chr1	4000	-	+	UU
r7	chr1	3001	chr1	4002	+	-	UU
```

- **Example code:**

```python
from dnakit.epigenome import deduplicate_pairs
from dnakit.epigenome import deduplicate_pairs_online

online = deduplicate_pairs_online("dup.pairs", max_mismatch=2, verbose=False)
print(online["total_nodups"], online["total_dups"])
print(online["unique_read_ids"])
print(sorted(online["duplicate_read_ids"]))
chained = deduplicate_pairs("dup.pairs", max_mismatch=2, verbose=False)
print(chained["total_nodups"], chained["total_dups"])
```

- **Example output:**

```text
4 4
['r1', 'r3', 'r6', 'r8']
['r2', 'r4', 'r5', 'r7']
3 5
```

## 21) `dnakit.epigenome.normalize_scaling_curve_unity()`

- **Purpose:** Scale a contact-frequency curve so that its area under the curve is 1.
- **Method:** Calls pairtools.lib.scaling.unity_norm_scaling. Every frequency is divided by sum(bin_length * frequency) over the finite bins whose geometric midpoint sqrt(low * high) lies strictly inside norm_range; norm_range=None uses all finite bins. normalize_scaling_curve instead divides by the mean frequency in a window.
- **API:** `dnakit.epigenome.normalize_scaling_curve_unity(bin_edges, frequencies, *, norm_range=(10000, 1000000000), verbose=True)`.
- **Input:** bin_edges are n + 1 ascending integer edges; frequencies are n non-negative values; norm_range is [low, high] or None.
- **Output:** bins (min_dist, max_dist, midpoint, frequency), norm_range and provenance.
- **Example code:**

```python
from dnakit.epigenome import normalize_scaling_curve_unity

edges = [10, 100, 1000, 10000]
frequencies = [0.5, 0.1, 0.01]
result = normalize_scaling_curve_unity(edges, frequencies, norm_range=None, verbose=False)
print([round(item["frequency"], 6) for item in result["bins"]])
area = sum((high - low) * f for low, high, f in zip(edges[:-1], edges[1:], frequencies))
print(area)
```

- **Example output:**

```text
[0.002222, 0.000444, 4.4e-05]
225.0
```

## 22) `dnakit.epigenome.chop_regions_from_middle()`

- **Purpose:** Cut a fixed number of bases on each side of the centre of a set of exons.
- **Method:** Calls deeptools.heatmapper.chopRegionsFromMiddle. The centre is the point after total_exon_length // 2 exonic bases; an exon that spans it is split. The last left bases before the centre and the first right bases after it are returned as lists of [start, end] intervals; when a side has fewer bases than requested, pad_left or pad_right holds the missing count. Unlike chop_regions_for_profile the intervals are cut from the centre and there is no body.
- **API:** `dnakit.epigenome.chop_regions_from_middle(exons, *, left=0, right=0, verbose=True)`.
- **Input:** exons is a sorted, non-overlapping list of [start, end]; left and right are base counts.
- **Output:** left_bins, right_bins, pad_left, pad_right and provenance.
- **Example code:**

```python
from dnakit.epigenome import chop_regions_from_middle

exons = [[0, 200], [300, 400], [800, 900]]
result = chop_regions_from_middle(exons, left=150, right=80, verbose=False)
print(result["left_bins"], result["right_bins"], result["pad_left"], result["pad_right"])
short = chop_regions_from_middle([[0, 50]], left=100, right=10, verbose=False)
print(short["left_bins"], short["right_bins"], short["pad_left"], short["pad_right"])
```

- **Example output:**

```text
[[50, 200]] [[300, 380]] 0 0
[[0, 25]] [[25, 35]] 75 0
```

## 23) `dnakit.epigenome.convert_paired_fastq_to_bisulfite()`

- **Purpose:** Convert the reads of one mate file of a paired-end library for bisulfite alignment.
- **Method:** Calls methylpy.call_mc_pe.convert_reads_pe. For read 1 (is_read2=False) every C becomes T, as in convert_fastq_to_bisulfite; for read 2 (is_read2=True) every G becomes A. In both cases the positions of the original C (read 1) or G (read 2) bases are encoded in the read name after a !, and spaces in the name also become !. The output file must not exist yet.
- **API:** `dnakit.epigenome.convert_paired_fastq_to_bisulfite(fastq_path, output_path, *, is_read2=False, buffer_line_number=100000, verbose=True)`.
- **Input:** fastq_path is a FASTQ file; output_path is a new file; is_read2 is a boolean; buffer_line_number.
- **Output:** reads, converted_base, is_read2, input, output and provenance.

**Example file:** `reads.fq`

```text
@r1 x
ACCGTTCA
+
IIIIIIII
@r2
GGATCC
+
IIIIII
```

- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import convert_paired_fastq_to_bisulfite

with tempfile.TemporaryDirectory() as folder:
    for mate in (False, True):
        output = Path(folder) / f"mate_{int(mate) + 1}.fq"
        result = convert_paired_fastq_to_bisulfite("reads.fq", output, is_read2=mate, verbose=False)
        lines = output.read_text().splitlines()
        print(result["converted_base"], result["reads"], lines[1], lines[5])
```

- **Example output:**

```text
C to T 2 ATTGTTTA GGATTT
G to A 2 ACCATTCA AAATCC
```

## 24) `dnakit.variants.load_maf_table()`

- **Purpose:** Read a TCGA-style MAF file as a table with all of its columns.
- **Method:** Calls varcode.maf.load_maf_dataframe. The file needs the 19 standard MAF columns in the standard order (names are matched ignoring case and renamed to the standard spelling) and may have more. load_maf_variants keeps only contig, start, ref and alt; this function also returns gene symbol, variant class, sample barcodes and every other column. Empty cells are null. max_rows limits the rows read.
- **API:** `dnakit.variants.load_maf_table(path, *, max_rows=None, verbose=True)`.
- **Input:** path is a local MAF file; max_rows is optional.
- **Output:** row_count, columns, rows (one mapping per row) and provenance.

**Example file:** `sample.maf`

```text
Hugo_Symbol	Entrez_Gene_Id	Center	NCBI_Build	Chromosome	Start_Position	End_Position	Strand	Variant_Classification	Variant_Type	Reference_Allele	Tumor_Seq_Allele1	Tumor_Seq_Allele2	dbSNP_RS	dbSNP_Val_Status	Tumor_Sample_Barcode	Matched_Norm_Sample_Barcode	Match_Norm_Seq_Allele1	Match_Norm_Seq_Allele2
TP53			GRCh38	17	7676154	7676154	+	Missense_Mutation	SNP	C	C	T			TCGA-AA-0001-01A			
```

- **Example code:**

```python
from dnakit.variants import load_maf_table

result = load_maf_table("sample.maf", verbose=False)
row = result["rows"][0]
print(result["row_count"], len(result["columns"]))
print(row["Hugo_Symbol"], row["Variant_Classification"], row["Start_Position"], row["Entrez_Gene_Id"])
```

- **Example output:**

```text
1 19
TP53 Missense_Mutation 7676154 None
```

## 25) `dnakit.variants.read_vcf_table()`

- **Purpose:** Read a VCF as a raw text table with varcode.
- **Method:** Calls varcode.vcf.read_vcf_into_dataframe. Every column is the text of the VCF field: ALT stays comma separated and INFO (only with include_info=True) is not parsed. When sample_names is given (it needs include_info=True) the FORMAT column and those sample columns are added. varcode only labels the first sample columns in file order, so the names must be the first samples of the VCF header. Each sample cell is the whole colon-separated string, so genotype quality and depth are kept, unlike read_vcf.
- **API:** `dnakit.variants.read_vcf_table(path, *, include_info=False, sample_names=None, verbose=True)`.
- **Input:** path is a VCF (optionally .gz); include_info is a boolean; sample_names is optional.
- **Output:** row_count, columns, rows and provenance.

**Example file:** `calls.vcf`

```text
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
1	100	.	A	G,T	50	PASS	DP=10;AF=0.1,0.2	GT:GQ	0/1:35	1/2:40
1	200	rs1	CT	C	60	PASS	DP=20	GT:GQ	0/0:99	0/1:12
```

- **Example code:**

```python
from dnakit.variants import read_vcf_table

result = read_vcf_table("calls.vcf", include_info=True, sample_names=["S1", "S2"], verbose=False)
print(result["columns"])
first = result["rows"][0]
print(first["ALT"], first["INFO"], first["FORMAT"], first["S1"], first["S2"])
```

- **Example output:**

```text
['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'S1', 'S2']
G,T DP=10;AF=0.1,0.2 GT:GQ 0/1:35 1/2:40
```

## 26) `dnakit.variants.export_vcf_to_phasing_format()`

- **Purpose:** Export VCF genotypes to the IMPUTE2 or SHAPEIT text formats with bcftools convert.
- **Method:** fmt='gensample' writes <prefix>.gen.gz (three genotype probabilities per sample and site) and <prefix>.samples; 'hapsample' writes <prefix>.hap.gz (two haplotype alleles per sample, so phased genotypes are expected) and <prefix>.samples; 'haplegendsample' writes .hap.gz, .legend.gz and .samples. Sites are named chrom:pos_ref_alt. The bcftools build bundled with pysam is called with catch_stdout=False, because otherwise it writes no files. The output files must not exist yet.
- **API:** `dnakit.variants.export_vcf_to_phasing_format(source, output_prefix, *, fmt='gensample', max_rows=100000, verbose=True)`.
- **Input:** source is a VCF path; output_prefix is the output prefix (its directory must exist); fmt; max_rows limits the returned rows.
- **Output:** files, sample_ids, site_count, rows (lines of the main data file), rows_truncated and provenance.

**Example file:** `ph.vcf`

```text
##fileformat=VCFv4.2
##contig=<ID=1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
1	100	rs1	A	G	.	PASS	.	GT	0|1	1|1
1	200	rs2	C	T	.	PASS	.	GT	0|0	1|0
```

- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants import export_vcf_to_phasing_format

with tempfile.TemporaryDirectory() as folder:
    for fmt in ("gensample", "hapsample", "haplegendsample"):
        result = export_vcf_to_phasing_format("ph.vcf", Path(folder) / fmt, fmt=fmt, verbose=False)
        print(fmt, [Path(name).name for name in result["files"]], result["sample_ids"])
        print(result["rows"][0])
```

- **Example output:**

```text
gensample ['gensample.gen.gz', 'gensample.samples'] ['S1', 'S2']
1:100_A_G 1:100_A_G 100 A G 0 1 0 0 0 1
hapsample ['hapsample.hap.gz', 'hapsample.samples'] ['S1', 'S2']
1 1:100_A_G 100 A G 0 1 1 1
haplegendsample ['haplegendsample.hap.gz', 'haplegendsample.legend.gz', 'haplegendsample.samples'] ['S1', 'S2']
0 1 1 1
```

## Main wrong claims

- gimmemotifs: `motif_score_metric(metric='matches_at_fpr')` raises TypeError (`float()` of a list), because upstream returns two counts, and `metric='roc_auc_xlim'` passes ROC curve points to a function that expects the score lists and fails with BackendExecutionError. Both rows became `wrapped`. `peakfile2bedfile` is not `equal_size_intervals(keep_columns=False)`: it is strand aware (start larger than end on non-plus strands), drops rows with start below 0 where the other clamps to 0, and always writes six columns.
- sourmash: `nodegraph.extract_nodegraph_info` has the same code as khmer's and both round the table size to hundreds, so the ledger remark about no rounding was wrong; `nodegraph.calc_expected_collisions` gives the same hash sizes, occupied bins and false-positive rate as `expected_kmer_table_collisions` on one FASTA file (now `duplicate-of`). `compare_all_pairs` brings abundance-weighted angular similarity and ANI output, and `search_databases_with_flat_query` and `abund_query` rank by score; no existing wrapper covers them, so they are wrapped.
- pairtools: `streaming_dedup_cython` cannot be reached through `deduplicate_pairs` (which fixes the scipy backend), and its result differs: on 8 pairs with max_mismatch 2 the online algorithm keeps 4 and the clustering one keeps 3. `unity_norm_scaling` divides by the area under the curve, which `norm_window` (window mean) does not do.
- cooltools: `sandbox.sum_pairwise` is not a different set of diagonals; it returns the diagonal rows plus one total row for every pair of different chromosomes. `compute_scaling` takes position lists and returns raw counts and areas; it is not a raw-count variant of `compute_contact_scaling`.
- tangermeme: `place_motif(mode='insert')` trims the inserted sequences back to the input length and loses the last bases; the beam width and n_best of `beam_substitution` do not exist in `design_greedy_substitution`.
- MACS3: `pileup_from_PN` equals `build_shifted_read_pileup(five_shift=0, three_shift=extsize)` apart from a leading zero point at a negative position; `se_all_in_one_pileup` gives the same output as `pileup_from_PN_shifted`. `CalModel_Heter_AS` and `calculate_GQ_heterASsig` cannot be reached through existing parameters.
- varcode: the tables of `load_maf_dataframe` and `read_vcf_into_dataframe` keep every column and the full FORMAT strings, while the existing wrappers keep contig, start, ref and alt or parsed genotypes. Upstream `sample_names` only labels the first sample columns in file order, so a later sample name returns wrong data; the wrapper checks the header.
- pysam: the ledger said `bcftools convert` GEN and HAP export wrote no files. The default `catch_stdout=True` swallows the output; with `catch_stdout=False` all three formats are written.
- methylpy: `convert_reads_pe` converts read 2 from G to A, which `convert_fastq_to_bisulfite` cannot do.
- deeptools: `chopRegionsFromMiddle` is not a variant of `chopRegions`; it cuts flanks on both sides of the exon centre and has no body.

## How these functions are verified

- Test file `others/tests/unit/test_family_H5.py`: 55 tests, all offline. Each new function has at least one hand calculation or an item-by-item comparison with the direct upstream call. For the 51 kept `family-variant` rows each package has a test showing that the stated parameter reaches the upstream function with an equal result, and alias rows are checked with `is`.
- All 26 new functions pass the agent registry's parameter analysis (`analyze_callable`); mapping parameters have string keys.
