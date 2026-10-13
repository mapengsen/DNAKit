# Optional alignment backends: minimap2, edlib, parasail, and pyhmmer

This group of functions in `dnakit.alignment` adapts optional packages: [mappy](https://github.com/lh3/minimap2) (the Python binding of minimap2), [edlib](https://github.com/Martinsos/edlib), [parasail](https://github.com/jeffdaily/parasail), and [pyhmmer](https://github.com/althonos/pyhmmer). They compute DNA mapping positions, edit distances, scored alignments, and profile HMM searches. They do not download reference databases or call variant callers or assemblers.

The dependencies are optional and are imported only when the matching function is called. If a package is missing, the function raises an error that includes the install command.

## Example data

The `reference.fa` file below has only two DNA sequences. It is for demonstrating the interface. Alignment results apply only to these short sequences and do not describe a real genome.

**Example file:** `reference.fa`

```text
>chr1
CAAGTGTCCGGCCATGTCCCTCTCAGATTCGAACACAGAACAGAGGATGGACACGGGCGGAGCCAGACGAGTCCTTGAAGCGACGCTGAAATTCCGTCATAGATCATGTGAGACTTTCTG
>chr2
GTCGTGTTCAACAGCCGTTGTTCTCTCCAATGGGAGAGGGGCACTTTTTTATCGAGGATACCGCTGTGCAATTAGTGCGTAGCCTAAATCTGGAATACAT
```

## 1) `dnakit.alignment.map_query_sequences_paf()`

- **Purpose:** Map a set of named DNA query sequences to a FASTA or `.mmi` reference with minimap2 and return PAF-style records.
- **Method:** Call `mappy.Aligner.map()` for each query, then convert each hit to PAF fields. `qstart`, `qend`, `tstart`, and `tend` are zero-based half-open coordinates. `nmatch` is the number of matching bases, and `alen` is the alignment block length. `edit_distance` comes from the NM tag, and `cigar` is minimap2's CIGAR. When `strand` is `-`, the query was aligned as its reverse complement.
- **API:** `dnakit.alignment.map_query_sequences_paf(reference, queries, *, preset="map-ont", max_queries=10000, max_hits_per_query=100)`.
- **Input:** `reference` is a FASTA file or a `.mmi` index written by `build_minimap2_index()`. `queries` is `{name: sequence}` or a DNA object and may contain only A, C, G, T, and N. `preset` is `sr`, `map-ont`, `map-pb`, `map-hifi`, `asm5`, `asm10`, or `asm20`. An error is raised if there are more than `max_queries` queries or if one query has more than `max_hits_per_query` hits.
- **Output:** A `records` list. Each record contains `qname`, `qlen`, `qstart`, `qend`, `strand`, `tname`, `tlen`, `tstart`, `tend`, `nmatch`, `alen`, `mapq`, `primary`, `edit_distance`, and `cigar`. The result also has `query_count`, `hit_count`, and `provenance`.
- **Example code:**

```python
from dnakit.alignment import map_query_sequences_paf

queries = {
    "read1": "TCTCAGATTCGAACACAGAACAGAGGATGGACACGGGCGGAGCCAGACGAGTCCTTGAAG",
    "read2": "ATGACGGAATTTCAGCGTCGCTTCAAGGACTCGTCTGGCTCCGCCCGTGTCCATCCTCTG",
}
result = map_query_sequences_paf("reference.fa", queries)
for record in result["records"]:
    print(record["qname"], record["tname"], record["tstart"], record["tend"], record["strand"], record["cigar"])
```

- **Example output:**

```text
read1 chr1 20 80 + 60M
read2 chr1 40 100 - 60M
```

## 2) `dnakit.alignment.build_minimap2_index()`

- **Purpose:** Write a minimap2 `.mmi` index for a FASTA reference and list its contig names.
- **Method:** Call `mappy.Aligner(reference, preset=..., fn_idx_out=...)` to write the index. If the output path already exists, the function raises an error unless `overwrite=True`.
- **API:** `dnakit.alignment.build_minimap2_index(reference, output_path, *, preset="map-ont", overwrite=False)`.
- **Input:** A FASTA reference and an output path whose parent directory already exists. `preset` accepts the same values as in `map_query_sequences_paf()`.
- **Output:** `output` (the index path), `contig_count`, `contigs` (the contig names), and `index_bytes`.
- **Example code:**

```python
from dnakit.alignment import build_minimap2_index

result = build_minimap2_index("reference.fa", "reference.mmi", overwrite=True)
print(result["contig_count"], result["contigs"])
```

- **Example output:**

```text
2 ['chr1', 'chr2']
```

## 3) `dnakit.alignment.fetch_minimap2_reference_region()`

- **Purpose:** Extract a specified interval from one contig of a reference.
- **Method:** Slice the sequence store that mappy loads. Coordinates are zero-based and half-open. mappy silently truncates intervals that run past the end of a contig, so this function checks the returned length and raises an error on a mismatch.
- **API:** `dnakit.alignment.fetch_minimap2_reference_region(reference, contig, start, end, *, max_reference_bytes=2000000000)`.
- **Input:** A FASTA or `.mmi` reference, a contig name, and integer coordinates with `0 <= start < end`. The function raises an error if the reference file is larger than `max_reference_bytes`.
- **Output:** `sequence`, `length`, `contig`, `start`, and `end`.
- **Example code:**

```python
from dnakit.alignment import fetch_minimap2_reference_region

result = fetch_minimap2_reference_region("reference.fa", "chr2", 0, 12)
print(result["sequence"], result["length"])
```

- **Example output:**

```text
GTCGTGTTCAAC 12
```

## 4) `dnakit.alignment.align_dna_edlib()`

- **Purpose:** Align two DNA sequences with Edlib using global (`NW`) or infix (`HW`) edit-distance alignment, and return an extended CIGAR with base counts.
- **Method:** Call `edlib.align(task="path")`. The CIGAR uses `=`, `X`, `I`, and `D`: `I` is a base present in the query and absent from the target, and `D` is a base present in the target and absent from the query. With `iupac=True`, an IUPAC ambiguity code equals each base it represents. By default, it equals only itself. In `HW` mode, the query is aligned to one substring of the target, and `target_start` and `target_end` give that substring.
- **API:** `dnakit.alignment.align_dna_edlib(query, target, *, mode="NW", k=-1, iupac=False, max_length=100000)`.
- **Input:** Two DNA sequences, which may contain IUPAC symbols. `mode` is `NW` or `HW`. `k` is the largest allowed edit distance, and `-1` means no limit.
- **Output:** `found`, `edit_distance`, `cigar`, `query_start`, `query_end`, `target_start`, `target_end`, `matches`, `mismatches`, `insertions`, `deletions`, and `identity`. If the edit distance exceeds `k`, `found` is `False`, and the distance and CIGAR are `None`.
- **Example code:**

```python
from dnakit.alignment import align_dna_edlib

for iupac in (False, True):
    result = align_dna_edlib("ACRT", "ACGT", iupac=iupac)
    print(iupac, result["edit_distance"], result["cigar"])
```

- **Example output:**

```text
False 1 2=1X1=
True 0 4=
```

## 5) `dnakit.alignment.locate_approximate_pattern()`

- **Purpose:** Find the approximate matches with the lowest edit distance in a text. This suits short primers or motifs.
- **Method:** Call `edlib.align(mode="HW", task="locations", k=max_edits)`. The result lists the match spans that reach the lowest edit distance. It does not list every position within `max_edits`. Spans are zero-based and half-open.
- **API:** `dnakit.alignment.locate_approximate_pattern(pattern, text, *, max_edits=0, iupac=False, max_length=100000)`.
- **Input:** A DNA pattern and a DNA text. `max_edits` is an integer from 0 to 1000.
- **Output:** `found`, `edit_distance`, `locations` (each item has `start` and `end`), and `location_count`.
- **Example code:**

```python
from dnakit.alignment import locate_approximate_pattern

result = locate_approximate_pattern("ACGT", "TTACGTTTACGTT", max_edits=0)
print(result["edit_distance"], result["locations"])
```

- **Example output:**

```text
0 [{'start': 2, 'end': 6}, {'start': 8, 'end': 12}]
```

## 6) `dnakit.alignment.calculate_edit_distance_matrix()`

- **Purpose:** Compute the global edit-distance matrix for every pair of named DNA sequences.
- **Method:** For each pair, call `edlib.align(mode="NW", task="distance")`. The matrix is symmetric, and its diagonal is 0.
- **API:** `dnakit.alignment.calculate_edit_distance_matrix(sequences, *, iupac=False, max_sequences=100, max_length=5000)`.
- **Input:** `{name: sequence}` or a DNA object. The number of sequences and the length of each sequence have explicit limits.
- **Output:** `names`, `matrix` (in the order of `names`), and `pair_count`.
- **Example code:**

```python
from dnakit.alignment import calculate_edit_distance_matrix

result = calculate_edit_distance_matrix({"a": "ACGT", "b": "ACGT", "c": "ACGA", "d": "TTTT"})
print(result["names"], result["matrix"])
```

- **Example output:**

```text
['a', 'b', 'c', 'd'] [[0, 0, 1, 3], [0, 0, 1, 3], [1, 1, 0, 4], [3, 3, 4, 0]]
```

## 7) `dnakit.alignment.align_dna_parasail()`

- **Purpose:** Align two DNA sequences with parasail in global (`global`), local (`local`), or semi-global (`semiglobal`) mode, using a substitution matrix and affine gap penalties.
- **Method:** With `matrix="custom"`, `match` and `mismatch` build an ACGTN substitution matrix. With `matrix="dnafull"`, parasail's built-in EDNAFULL matrix is used, and `match` and `mismatch` are ignored. A gap of length k costs `gap_open + (k - 1) * gap_extend`. In local and semi-global modes, the spans cover only the part where both sequences have residues, and they exclude free end gaps. The `cigar` covers the whole traceback, so in semi-global mode it may include free gaps at either end.
- **API:** `dnakit.alignment.align_dna_parasail(query, target, *, mode="global", matrix="custom", match=5, mismatch=-4, gap_open=10, gap_extend=1, max_cells=25000000)`.
- **Input:** Only A, C, G, T, and N are allowed. `mode` is `global`, `local`, or `semiglobal`. `max_cells` limits `len(query) * len(target)`.
- **Output:** `score`, `found`, `cigar`, `query_aligned`, `target_aligned`, `query_start`, `query_end`, `target_start`, `target_end`, `matches`, `mismatches`, `insertions`, `deletions`, and `identity`.
- **Example code:**

```python
from dnakit.alignment import align_dna_parasail

result = align_dna_parasail(
    "TTTGATTACATTT", "GGGGATTACAGGG", mode="local", match=2, mismatch=-1, gap_open=5, gap_extend=1
)
print(result["score"], result["cigar"], result["query_start"], result["query_end"])
print(result["target_start"], result["target_end"], result["identity"])
```

- **Example output:**

```text
14 7= 3 10
3 10 1.0
```

## 8) `dnakit.alignment.score_dna_parasail()`

- **Purpose:** Compute only the optimal score and the end cell of the optimal alignment, without a traceback.
- **Method:** Call parasail's scan variants (`*_scan_sat`). These do not store a traceback matrix, so they use less memory than `align_dna_parasail()`. If the local score is 0, the end coordinates are `None`.
- **API:** `dnakit.alignment.score_dna_parasail(query, target, *, mode="global", matrix="custom", match=5, mismatch=-4, gap_open=10, gap_extend=1, max_cells=25000000)`.
- **Input:** The same as `align_dna_parasail()`.
- **Output:** `score`, `query_end`, and `target_end`, all zero-based half-open coordinates.
- **Example code:**

```python
from dnakit.alignment import score_dna_parasail

options = {"match": 2, "mismatch": -1, "gap_open": 5, "gap_extend": 1}
result = score_dna_parasail("GATTACA", "TTTGATTACATTT", mode="semiglobal", **options)
print(result["score"], result["query_end"], result["target_end"])
```

- **Example output:**

```text
14 7 10
```

## 9) `dnakit.alignment.score_dna_parasail_matrix()`

- **Purpose:** Compute the parasail score for every pair of named DNA sequences.
- **Method:** Compute the upper triangle and mirror it to the lower triangle. Scores follow the rules of `score_dna_parasail()`.
- **API:** `dnakit.alignment.score_dna_parasail_matrix(sequences, *, mode="global", matrix="custom", match=5, mismatch=-4, gap_open=10, gap_extend=1, max_sequences=50, max_length=2000, max_cells=25000000)`.
- **Input:** `{name: sequence}` or a DNA object. The number of sequences, the length of each sequence, and the DP cells for each pair have explicit limits.
- **Output:** `names`, `scores` (in the order of `names`), and `pair_count`.
- **Example code:**

```python
from dnakit.alignment import score_dna_parasail_matrix

sequences = {"x": "ACGTACGT", "y": "ACGTTACGT", "z": "TTTTGGGG"}
result = score_dna_parasail_matrix(
    sequences, mode="local", match=2, mismatch=-1, gap_open=5, gap_extend=1
)
print(result["names"], result["scores"])
```

- **Example output:**

```text
['x', 'y', 'z'] [[16, 11, 2], [11, 18, 4], [2, 4, 16]]
```

## 10) `dnakit.alignment.build_dna_hmm()`

- **Purpose:** Build a HMMER3 profile HMM from an aligned DNA file, and optionally write it to disk.
- **Method:** Read the multiple alignment with pyhmmer in the DNA alphabet, then build the HMM with a `Builder` and `Background` that use default settings. The HMM name defaults to the file name without its extension. Pass `name` to set it explicitly.
- **API:** `dnakit.alignment.build_dna_hmm(alignment, *, name=None, output_path=None, overwrite=False, max_sequences=10000, max_columns=100000)`.
- **Input:** A DNA alignment file in which every sequence has the same length. Gaps are written as `-`. `output_path` is an optional HMMER3 output path.
- **Output:** `name`, `length` (the model length, which equals the number of alignment columns), `nseq`, `effective_nseq`, `consensus`, and `output`.
- **Example file:** `msa.fa`

```text
>seqA
ACGTACGTTAGCAAGGCTTA
>seqB
ACGTACGATAGCAAGGCTTA
>seqC
ACGTTCGTTAGCAAGG-TTA
>seqD
ACGTACGTTAGCTAGGCTTA
```

- **Example code:**

```python
from dnakit.alignment import build_dna_hmm

result = build_dna_hmm("msa.fa", output_path="msa.hmm", overwrite=True)
print(result["name"], result["length"], result["nseq"], result["consensus"])
```

- **Example output:**

```text
msa 20 4 AcgtacgttAgcaAggcttA
```

## 11) `dnakit.alignment.search_dna_hmm()`

- **Purpose:** Search a DNA sequence file with a DNA profile HMM, checking both the forward and reverse strands.
- **Method:** Call pyhmmer's `nhmmer` using its long-target pipeline. `strand` is `both` (the default), `forward` (top strand only), or `reverse` (bottom strand only). `evalue` is the reporting threshold. Output coordinates are 1-based and inclusive, with `start <= end`. The direction is given by `strand`.
- **API:** `dnakit.alignment.search_dna_hmm(profile, sequences, *, evalue=10.0, strand="both", max_hits=1000, max_sequences=100000, max_residues=100000000)`.
- **Input:** A HMMER3 file that contains exactly one DNA profile, such as a file written by `build_dna_hmm()`, and a FASTA or FASTQ sequence file. Protein profiles are rejected.
- **Output:** A `hits` list. Each reported domain has `target`, `start`, `end`, `strand`, `hmm_start`, `hmm_end`, `score`, `evalue`, and `sequence_evalue`. The result also has `hit_count` and `sequences_searched`.
- **Example file:** `targets.fa`

```text
>t_fwd
GGGGACGTACGTTAGCAAGGCTTAGGGG
>t_rev
CCCCTAAGCCTTGCTAACGTACGTCCCC
>t_none
TTTTTTTTTTTTTTTTTTTTTTTTTTTT
>t_mut
CCCCACGTACGATAGCAAGGCTTACCCC
```

- **Example code:**

```python
from dnakit.alignment import search_dna_hmm

result = search_dna_hmm("msa.hmm", "targets.fa")
for hit in result["hits"]:
    print(hit["target"], hit["start"], hit["end"], hit["strand"])
```

- **Example output:**

```text
t_fwd 5 24 +
t_rev 5 24 -
t_mut 5 24 +
```

## 12) `dnakit.alignment.read_hmm_profile()`

- **Purpose:** Read each profile's name, alphabet, length, consensus, and sequence counts from a HMMER3 file, without searching.
- **Method:** Read the profiles one at a time with `pyhmmer.plan7.HMMFile`, and limit how many profiles are read.
- **API:** `dnakit.alignment.read_hmm_profile(path, *, max_profiles=1000)`.
- **Input:** The path to a HMMER3 file. The file may contain several profiles.
- **Output:** A `profiles` list. Each item has `name`, `accession`, `description`, `alphabet`, `length`, `consensus`, `nseq`, and `effective_nseq`. The result also has `profile_count`.
- **Example code:**

```python
from dnakit.alignment import read_hmm_profile

profile = read_hmm_profile("msa.hmm")["profiles"][0]
print(profile["name"], profile["alphabet"], profile["length"], profile["nseq"])
```

- **Example output:**

```text
msa DNA 20 4
```

## Installing the optional dependencies

```bash
python -m pip install "dnakit[bio-mappy]"           # mappy (minimap2): map_query_sequences_paf and related functions
python -m pip install "dnakit[bio-edlib]"           # edlib: align_dna_edlib and related functions
python -m pip install "dnakit[bio-parasail]"        # parasail: align_dna_parasail and related functions
python -m pip install "dnakit[bio-pyhmmer]"         # pyhmmer: build_dna_hmm and related functions
```

Project pages: [minimap2](https://github.com/lh3/minimap2), [edlib](https://github.com/Martinsos/edlib), [parasail](https://github.com/jeffdaily/parasail), [pyhmmer](https://github.com/althonos/pyhmmer).

## How these functions are verified

The test file is `others/tests/unit/test_alignment_backends.py`, with 30 tests. They check the following:

- mappy: the coordinates and MAPQ of an exact forward hit match a direct mappy call. The strand and coordinates of a reverse-complement hit match hand calculation. The NM tag and matching-base count for one substitution match hand calculation. A written index can be reloaded, and an interval past the end of a contig is rejected.
- edlib: the CIGAR and edit distance of `4=1D4=` match a direct edlib call. The IUPAC option and the `HW` target span match hand calculation. Locate spans match direct edlib `locations`, and the distance matrix matches pairwise calls.
- parasail: scores of 20, 11, 10, 4, and 14 match hand calculation. Local and semi-global spans match hand results. Scores match direct parasail calls. The score matrix matches pairwise calls and is symmetric.
- pyhmmer: the model length and consensus match a direct pyhmmer build, and the sequence count matches the input. The hit coordinates for the three sequences (forward, reverse, and one substitution) are 5 to 24, which matches hand calculation. The `strand` filter works, and protein profiles are rejected.
- Each module has a test that a missing optional package raises an error with the install hint. Other tests cover unknown arguments, values above explicit limits, and booleans passed as integers.

## Current scope

This batch covers in-memory minimap2 mapping, minimap2 index writing and reading, Edlib NW and HW alignment and distance matrices, parasail global, local, and semi-global alignment and scoring, and pyhmmer DNA profile construction, nhmmer search, and metadata reading. The following are not included yet:

- minimap2 splice alignment, the `cs` tag, and arbitrary command-line options. Mapping results report only the CIGAR and NM.
- Edlib's `SHW` (prefix) mode, and batch alignment of several queries against one target.
- IUPAC scoring in parasail, protein substitution matrices such as BLOSUM, and the striped and profile interfaces.
- pyhmmer protein searches (`hmmsearch` and `hmmscan`), `jackhmmer`, `phmmer`, and HMM `cutoffs` settings.
