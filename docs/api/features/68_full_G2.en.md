# Alignment extensions: parasail, edlib and pyhmmer

This group extends the parasail, edlib and pyhmmer adapters already in `dnakit.alignment`. parasail's 5479 inventory entries (including the `bindings_v1` and `bindings_v2` mirrors) fall into 13 families of algorithm by end-gap mode. Each family has one representative, and SIMD method and lane width are parameters that select the variant. This page wraps ten semi-global end-gap families, banded global scores and SSW local alignment, plus one edlib NICE text view and one pyhmmer DNA profile alignment. Every function works on sequences or local files passed in explicitly, with no network access.

The examples import from the module that defines each function.

## Example data

The file below is only used by the `pyhmmer_align_to_hmm()` example to build a DNA profile. The other examples use inline strings.

**Example file:** `aligned.fa`

```text
>a
ACGTACGTAC
>b
ACGTTCGTAC
>c
ACGAACGTAC
```

## Relation to existing functions and upstream entries

| Upstream package | Inventory entries | Wrapped here | Handling of the other entries |
| --- | --- | --- | --- |
| mappy | 3 | none | `revcomp` and `fastx_read` duplicate `reverse_complement_dna` and `read_fastx`; `verbose` is a log setting |
| edlib | 2 | `edlib_nice_alignment` | `align` is already wrapped by `align_dna_edlib` and related functions |
| parasail | 5479 | 12 representatives (10 semi-global families, `nw_banded`, `ssw`) | the standard `nw`, `sw` and `sg` families are covered by `align_dna_parasail` and `score_dna_parasail`; the rest are family members, binding mirrors or configuration |
| pyhmmer | 18 | `pyhmmer_align_to_hmm` | `nhmmer` is covered by `search_dna_hmm`; `hmmsearch` duplicates it; `phmmer`, `jackhmmer` and `hmmscan` are protein searches; `hmmpress` only writes an index |

## 1) `dnakit.alignment.parasail_sg_qb_align()`

- **Purpose：** Semi-global alignment in which the query begin is free: an unaligned query prefix costs nothing.
- **Method：** Calls parasail's `sg_qb_trace` family; the representative is `sg_qb_trace_scan_sat`. `method` (`scan`, `striped`, `diag`) and `width` (8, 16, 32, 64 or `sat`) select the SIMD member and raise when a lane saturates. The traceback includes the free end-gap columns, so `cigar` covers the whole path; coordinates and identity count only the aligned core.
- **API：** `dnakit.alignment.parasail_sg_qb_align(query, target, *, match=5, mismatch=-4, gap_open=10, gap_extend=1, method="scan", width="sat", max_cells=25000000, verbose=True)`.
- **Input：** `query` and `target`: strings of ACGTN or DNA objects; `match` and `mismatch`: base scores; `gap_open` and `gap_extend`: affine gaps, a gap of length k costs `gap_open + (k - 1) * gap_extend`; `method` and `width`: SIMD member; `max_cells`: cap on dynamic-programming cells (default 25 million).
- **Output：** A dict with `score`, `found`, `free_end_gaps`, `query_start`, `query_end`, `target_start`, `target_end` (zero-based half-open), `query_aligned`, `target_aligned`, `cigar`, `matches`, `mismatches`, `insertions`, `deletions`, `alignment_columns`, `identity` and `provenance`.
- **Example code：**

```python
from dnakit.alignment import parasail_sg_qb_align

result = parasail_sg_qb_align("TTTTACGTACGGA", "ACGTACGGA")
print(result["score"], result["free_end_gaps"], result["cigar"])
print(result["query_start"], result["query_end"], result["target_start"], result["target_end"])
```

- **Example output：**

```text
45 ['query_begin'] 4I9=
4 13 0 9
```

## 2) `dnakit.alignment.parasail_sg_qe_align()`

- **Purpose：** Semi-global alignment in which the query end is free: an unaligned query suffix costs nothing.
- **Method：** Calls parasail's `sg_qe_trace` family; the representative is `sg_qe_trace_scan_sat`. `method` (`scan`, `striped`, `diag`) and `width` (8, 16, 32, 64 or `sat`) select the SIMD member and raise when a lane saturates. The traceback includes the free end-gap columns, so `cigar` covers the whole path; coordinates and identity count only the aligned core.
- **API：** `dnakit.alignment.parasail_sg_qe_align(query, target, *, match=5, mismatch=-4, gap_open=10, gap_extend=1, method="scan", width="sat", max_cells=25000000, verbose=True)`.
- **Input：** `query` and `target`: strings of ACGTN or DNA objects; `match` and `mismatch`: base scores; `gap_open` and `gap_extend`: affine gaps, a gap of length k costs `gap_open + (k - 1) * gap_extend`; `method` and `width`: SIMD member; `max_cells`: cap on dynamic-programming cells (default 25 million).
- **Output：** A dict with `score`, `found`, `free_end_gaps`, `query_start`, `query_end`, `target_start`, `target_end` (zero-based half-open), `query_aligned`, `target_aligned`, `cigar`, `matches`, `mismatches`, `insertions`, `deletions`, `alignment_columns`, `identity` and `provenance`.
- **Example code：**

```python
from dnakit.alignment import parasail_sg_qe_align

result = parasail_sg_qe_align("ACGTACGGATTTT", "ACGTACGGA")
print(result["score"], result["free_end_gaps"], result["cigar"])
print(result["query_start"], result["query_end"], result["target_start"], result["target_end"])
```

- **Example output：**

```text
45 ['query_end'] 9=4I
0 9 0 9
```

## 3) `dnakit.alignment.parasail_sg_qx_align()`

- **Purpose：** Semi-global alignment in which both query ends are free, so the query may overhang the target.
- **Method：** Calls parasail's `sg_qx_trace` family; the representative is `sg_qx_trace_scan_sat`. `method` (`scan`, `striped`, `diag`) and `width` (8, 16, 32, 64 or `sat`) select the SIMD member and raise when a lane saturates. The traceback includes the free end-gap columns, so `cigar` covers the whole path; coordinates and identity count only the aligned core.
- **API：** `dnakit.alignment.parasail_sg_qx_align(query, target, *, match=5, mismatch=-4, gap_open=10, gap_extend=1, method="scan", width="sat", max_cells=25000000, verbose=True)`.
- **Input：** `query` and `target`: strings of ACGTN or DNA objects; `match` and `mismatch`: base scores; `gap_open` and `gap_extend`: affine gaps, a gap of length k costs `gap_open + (k - 1) * gap_extend`; `method` and `width`: SIMD member; `max_cells`: cap on dynamic-programming cells (default 25 million).
- **Output：** A dict with `score`, `found`, `free_end_gaps`, `query_start`, `query_end`, `target_start`, `target_end` (zero-based half-open), `query_aligned`, `target_aligned`, `cigar`, `matches`, `mismatches`, `insertions`, `deletions`, `alignment_columns`, `identity` and `provenance`.
- **Example code：**

```python
from dnakit.alignment import parasail_sg_qx_align

result = parasail_sg_qx_align("TTTTACGTACGGATTTT", "ACGTACGGA")
print(result["score"], result["free_end_gaps"], result["cigar"])
print(result["query_start"], result["query_end"], result["target_start"], result["target_end"])
```

- **Example output：**

```text
45 ['query_begin', 'query_end'] 4I9=4I
4 13 0 9
```

## 4) `dnakit.alignment.parasail_sg_db_align()`

- **Purpose：** Semi-global alignment in which the database (target) begin is free.
- **Method：** Calls parasail's `sg_db_trace` family; the representative is `sg_db_trace_scan_sat`. `method` (`scan`, `striped`, `diag`) and `width` (8, 16, 32, 64 or `sat`) select the SIMD member and raise when a lane saturates. The traceback includes the free end-gap columns, so `cigar` covers the whole path; coordinates and identity count only the aligned core.
- **API：** `dnakit.alignment.parasail_sg_db_align(query, target, *, match=5, mismatch=-4, gap_open=10, gap_extend=1, method="scan", width="sat", max_cells=25000000, verbose=True)`.
- **Input：** `query` and `target`: strings of ACGTN or DNA objects; `match` and `mismatch`: base scores; `gap_open` and `gap_extend`: affine gaps, a gap of length k costs `gap_open + (k - 1) * gap_extend`; `method` and `width`: SIMD member; `max_cells`: cap on dynamic-programming cells (default 25 million).
- **Output：** A dict with `score`, `found`, `free_end_gaps`, `query_start`, `query_end`, `target_start`, `target_end` (zero-based half-open), `query_aligned`, `target_aligned`, `cigar`, `matches`, `mismatches`, `insertions`, `deletions`, `alignment_columns`, `identity` and `provenance`.
- **Example code：**

```python
from dnakit.alignment import parasail_sg_db_align

result = parasail_sg_db_align("ACGTACGGA", "TTTTACGTACGGA")
print(result["score"], result["free_end_gaps"], result["cigar"])
print(result["query_start"], result["query_end"], result["target_start"], result["target_end"])
```

- **Example output：**

```text
45 ['database_begin'] 4D9=
0 9 4 13
```

## 5) `dnakit.alignment.parasail_sg_de_align()`

- **Purpose：** Semi-global alignment in which the database (target) end is free.
- **Method：** Calls parasail's `sg_de_trace` family; the representative is `sg_de_trace_scan_sat`. `method` (`scan`, `striped`, `diag`) and `width` (8, 16, 32, 64 or `sat`) select the SIMD member and raise when a lane saturates. The traceback includes the free end-gap columns, so `cigar` covers the whole path; coordinates and identity count only the aligned core.
- **API：** `dnakit.alignment.parasail_sg_de_align(query, target, *, match=5, mismatch=-4, gap_open=10, gap_extend=1, method="scan", width="sat", max_cells=25000000, verbose=True)`.
- **Input：** `query` and `target`: strings of ACGTN or DNA objects; `match` and `mismatch`: base scores; `gap_open` and `gap_extend`: affine gaps, a gap of length k costs `gap_open + (k - 1) * gap_extend`; `method` and `width`: SIMD member; `max_cells`: cap on dynamic-programming cells (default 25 million).
- **Output：** A dict with `score`, `found`, `free_end_gaps`, `query_start`, `query_end`, `target_start`, `target_end` (zero-based half-open), `query_aligned`, `target_aligned`, `cigar`, `matches`, `mismatches`, `insertions`, `deletions`, `alignment_columns`, `identity` and `provenance`.
- **Example code：**

```python
from dnakit.alignment import parasail_sg_de_align

result = parasail_sg_de_align("ACGTACGGA", "ACGTACGGATTTT")
print(result["score"], result["free_end_gaps"], result["cigar"])
print(result["query_start"], result["query_end"], result["target_start"], result["target_end"])
```

- **Example output：**

```text
45 ['database_end'] 9=4D
0 9 0 9
```

## 6) `dnakit.alignment.parasail_sg_dx_align()`

- **Purpose：** Semi-global alignment in which both target ends are free: the query fits inside the target (fitting).
- **Method：** Calls parasail's `sg_dx_trace` family; the representative is `sg_dx_trace_scan_sat`. `method` (`scan`, `striped`, `diag`) and `width` (8, 16, 32, 64 or `sat`) select the SIMD member and raise when a lane saturates. The traceback includes the free end-gap columns, so `cigar` covers the whole path; coordinates and identity count only the aligned core.
- **API：** `dnakit.alignment.parasail_sg_dx_align(query, target, *, match=5, mismatch=-4, gap_open=10, gap_extend=1, method="scan", width="sat", max_cells=25000000, verbose=True)`.
- **Input：** `query` and `target`: strings of ACGTN or DNA objects; `match` and `mismatch`: base scores; `gap_open` and `gap_extend`: affine gaps, a gap of length k costs `gap_open + (k - 1) * gap_extend`; `method` and `width`: SIMD member; `max_cells`: cap on dynamic-programming cells (default 25 million).
- **Output：** A dict with `score`, `found`, `free_end_gaps`, `query_start`, `query_end`, `target_start`, `target_end` (zero-based half-open), `query_aligned`, `target_aligned`, `cigar`, `matches`, `mismatches`, `insertions`, `deletions`, `alignment_columns`, `identity` and `provenance`.
- **Example code：**

```python
from dnakit.alignment import parasail_sg_dx_align

result = parasail_sg_dx_align("ACGTACGGA", "TTTTACGTACGGATTTT")
print(result["score"], result["free_end_gaps"], result["cigar"])
print(result["query_start"], result["query_end"], result["target_start"], result["target_end"])
```

- **Example output：**

```text
45 ['database_begin', 'database_end'] 4D9=4D
0 9 4 13
```

## 7) `dnakit.alignment.parasail_sg_qb_db_align()`

- **Purpose：** Semi-global alignment in which the query begin or the database begin may be skipped, but not both at once; both sequence ends are aligned.
- **Method：** Calls parasail's `sg_qb_db_trace` family; the representative is `sg_qb_db_trace_scan_sat`. `method` (`scan`, `striped`, `diag`) and `width` (8, 16, 32, 64 or `sat`) select the SIMD member and raise when a lane saturates. The traceback includes the free end-gap columns, so `cigar` covers the whole path; coordinates and identity count only the aligned core.
- **API：** `dnakit.alignment.parasail_sg_qb_db_align(query, target, *, match=5, mismatch=-4, gap_open=10, gap_extend=1, method="scan", width="sat", max_cells=25000000, verbose=True)`.
- **Input：** `query` and `target`: strings of ACGTN or DNA objects; `match` and `mismatch`: base scores; `gap_open` and `gap_extend`: affine gaps, a gap of length k costs `gap_open + (k - 1) * gap_extend`; `method` and `width`: SIMD member; `max_cells`: cap on dynamic-programming cells (default 25 million).
- **Output：** A dict with `score`, `found`, `free_end_gaps`, `query_start`, `query_end`, `target_start`, `target_end` (zero-based half-open), `query_aligned`, `target_aligned`, `cigar`, `matches`, `mismatches`, `insertions`, `deletions`, `alignment_columns`, `identity` and `provenance`.
- **Example code：**

```python
from dnakit.alignment import parasail_sg_qb_db_align

result = parasail_sg_qb_db_align("TTTTACGTACGGA", "ACGTACGGA")
print(result["score"], result["free_end_gaps"], result["cigar"])
print(result["query_start"], result["query_end"], result["target_start"], result["target_end"])
```

- **Example output：**

```text
45 ['query_begin', 'database_begin'] 4I9=
4 13 0 9
```

## 8) `dnakit.alignment.parasail_sg_qb_de_align()`

- **Purpose：** Semi-global alignment with a free query begin and database end: a query suffix overlapping a target prefix.
- **Method：** Calls parasail's `sg_qb_de_trace` family; the representative is `sg_qb_de_trace_scan_sat`. `method` (`scan`, `striped`, `diag`) and `width` (8, 16, 32, 64 or `sat`) select the SIMD member and raise when a lane saturates. The traceback includes the free end-gap columns, so `cigar` covers the whole path; coordinates and identity count only the aligned core.
- **API：** `dnakit.alignment.parasail_sg_qb_de_align(query, target, *, match=5, mismatch=-4, gap_open=10, gap_extend=1, method="scan", width="sat", max_cells=25000000, verbose=True)`.
- **Input：** `query` and `target`: strings of ACGTN or DNA objects; `match` and `mismatch`: base scores; `gap_open` and `gap_extend`: affine gaps, a gap of length k costs `gap_open + (k - 1) * gap_extend`; `method` and `width`: SIMD member; `max_cells`: cap on dynamic-programming cells (default 25 million).
- **Output：** A dict with `score`, `found`, `free_end_gaps`, `query_start`, `query_end`, `target_start`, `target_end` (zero-based half-open), `query_aligned`, `target_aligned`, `cigar`, `matches`, `mismatches`, `insertions`, `deletions`, `alignment_columns`, `identity` and `provenance`.
- **Example code：**

```python
from dnakit.alignment import parasail_sg_qb_de_align

result = parasail_sg_qb_de_align("TTTTACGTACGGA", "ACGTACGGATTTT")
print(result["score"], result["free_end_gaps"], result["cigar"])
print(result["query_start"], result["query_end"], result["target_start"], result["target_end"])
```

- **Example output：**

```text
45 ['query_begin', 'database_end'] 4I9=4D
4 13 0 9
```

## 9) `dnakit.alignment.parasail_sg_qe_db_align()`

- **Purpose：** Semi-global alignment with a free query end and database begin: a query prefix overlapping a target suffix.
- **Method：** Calls parasail's `sg_qe_db_trace` family; the representative is `sg_qe_db_trace_scan_sat`. `method` (`scan`, `striped`, `diag`) and `width` (8, 16, 32, 64 or `sat`) select the SIMD member and raise when a lane saturates. The traceback includes the free end-gap columns, so `cigar` covers the whole path; coordinates and identity count only the aligned core.
- **API：** `dnakit.alignment.parasail_sg_qe_db_align(query, target, *, match=5, mismatch=-4, gap_open=10, gap_extend=1, method="scan", width="sat", max_cells=25000000, verbose=True)`.
- **Input：** `query` and `target`: strings of ACGTN or DNA objects; `match` and `mismatch`: base scores; `gap_open` and `gap_extend`: affine gaps, a gap of length k costs `gap_open + (k - 1) * gap_extend`; `method` and `width`: SIMD member; `max_cells`: cap on dynamic-programming cells (default 25 million).
- **Output：** A dict with `score`, `found`, `free_end_gaps`, `query_start`, `query_end`, `target_start`, `target_end` (zero-based half-open), `query_aligned`, `target_aligned`, `cigar`, `matches`, `mismatches`, `insertions`, `deletions`, `alignment_columns`, `identity` and `provenance`.
- **Example code：**

```python
from dnakit.alignment import parasail_sg_qe_db_align

result = parasail_sg_qe_db_align("ACGTACGGATTTT", "TTTTACGTACGGA")
print(result["score"], result["free_end_gaps"], result["cigar"])
print(result["query_start"], result["query_end"], result["target_start"], result["target_end"])
```

- **Example output：**

```text
45 ['query_end', 'database_begin'] 4D9=4I
0 9 4 13
```

## 10) `dnakit.alignment.parasail_sg_qe_de_align()`

- **Purpose：** Semi-global alignment in which the query end or the database end may be skipped, but not both at once; both sequence starts are aligned.
- **Method：** Calls parasail's `sg_qe_de_trace` family; the representative is `sg_qe_de_trace_scan_sat`. `method` (`scan`, `striped`, `diag`) and `width` (8, 16, 32, 64 or `sat`) select the SIMD member and raise when a lane saturates. The traceback includes the free end-gap columns, so `cigar` covers the whole path; coordinates and identity count only the aligned core.
- **API：** `dnakit.alignment.parasail_sg_qe_de_align(query, target, *, match=5, mismatch=-4, gap_open=10, gap_extend=1, method="scan", width="sat", max_cells=25000000, verbose=True)`.
- **Input：** `query` and `target`: strings of ACGTN or DNA objects; `match` and `mismatch`: base scores; `gap_open` and `gap_extend`: affine gaps, a gap of length k costs `gap_open + (k - 1) * gap_extend`; `method` and `width`: SIMD member; `max_cells`: cap on dynamic-programming cells (default 25 million).
- **Output：** A dict with `score`, `found`, `free_end_gaps`, `query_start`, `query_end`, `target_start`, `target_end` (zero-based half-open), `query_aligned`, `target_aligned`, `cigar`, `matches`, `mismatches`, `insertions`, `deletions`, `alignment_columns`, `identity` and `provenance`.
- **Example code：**

```python
from dnakit.alignment import parasail_sg_qe_de_align

result = parasail_sg_qe_de_align("ACGTACGGATTTT", "ACGTACGGA")
print(result["score"], result["free_end_gaps"], result["cigar"])
print(result["query_start"], result["query_end"], result["target_start"], result["target_end"])
```

- **Example output：**

```text
45 ['query_end', 'database_end'] 9=4I
0 9 0 9
```

## 11) `dnakit.alignment.parasail_nw_banded_score()`

- **Purpose：** Global score restricted to a band of half-width `band` around the main diagonal. A path that strays too far scores lower than the full global alignment.
- **Method：** Calls parasail's `nw_banded` and returns only the score (no traceback). When the lengths differ, parasail widens the band so that the final cell stays reachable.
- **API：** `dnakit.alignment.parasail_nw_banded_score(query, target, *, band=10, match=5, mismatch=-4, gap_open=10, gap_extend=1, max_cells=25000000, verbose=True)`.
- **Input：** `query` and `target`: ACGTN sequences; `band`: half-width (1 to 100000); scoring as above; `max_cells`: cap on cells inside the band (length of the query times 2 × band + 1).
- **Output：** A dict with `score`, `band`, `query_length`, `target_length` and `provenance`.
- **Example code：**

```python
from dnakit.alignment import parasail_nw_banded_score

query = "ACGTTGCAAGTCCATGGATCCAGTTACGGTAACGTTAGCCA"
target = query[:12] + "TTTTTT" + query[12:30] + query[36:]
print(parasail_nw_banded_score(query, target, band=1)["score"])
print(parasail_nw_banded_score(query, target, band=6)["score"])
```

- **Example output：**

```text
81
145
```

## 12) `dnakit.alignment.parasail_ssw_align()`

- **Purpose：** Best local DNA alignment from Striped Smith-Waterman (SSW), with score, coordinates, CIGAR and aligned text; suited to locating short reads.
- **Method：** Calls parasail's `ssw`. SSW keeps 16-bit scores, so parasail returns no result for scores near 32767, and this function raises. The SSW operation codes are decoded into `=`, `X`, `I` and `D`, which build the aligned text.
- **API：** `dnakit.alignment.parasail_ssw_align(query, target, *, match=5, mismatch=-4, gap_open=10, gap_extend=1, max_cells=25000000, verbose=True)`.
- **Input：** `query` and `target`: ACGTN sequences; scoring and gap parameters as above; `max_cells`: cap on dynamic-programming cells.
- **Output：** A dict with `score`, `found`, `query_start`, `query_end`, `target_start`, `target_end` (zero-based half-open), `query_aligned`, `target_aligned`, `cigar`, the count fields, `identity` and `provenance`. When no hit is found, `found` is False.
- **Example code：**

```python
from dnakit.alignment import parasail_ssw_align

result = parasail_ssw_align("ACGTTGCAAGTCCATG", "ACGTTGCAAGGGTCCATG")
print(result["score"], result["cigar"])
print(result["query_aligned"])
print(result["target_aligned"])
```

- **Example output：**

```text
69 9=2D7=
ACGTTGCAA--GTCCATG
ACGTTGCAAGGGTCCATG
```

## 13) `dnakit.alignment.edlib_nice_alignment()`

- **Purpose：** Align two DNA sequences with Edlib and return the three NICE text rows: query, match line and target.
- **Method：** Calls `edlib.align(task="path")` and then `edlib.getNiceAlignment()`. The edit distance and extended CIGAR match `align_dna_edlib()`. `mode` is `NW` (global) or `HW` (query global, target local).
- **API：** `dnakit.alignment.edlib_nice_alignment(query, target, *, mode="NW", k=-1, iupac=False, gap_symbol="-", max_length=100000, verbose=True)`.
- **Input：** `query` and `target`: DNA sequences; `mode`: `NW` or `HW`; `k`: edit-distance cap (-1 means no cap); `iupac`: make IUPAC ambiguity codes equal the bases they cover; `gap_symbol`: one printable non-alphanumeric, non-space character; `max_length`: length limit per sequence.
- **Output：** A dict with `found`, `edit_distance`, `query_aligned`, `matched_aligned`, `target_aligned`, `cigar`, `target_start`, `target_end`, `identity` and `provenance`. When the distance exceeds `k`, `found` is False.
- **Example code：**

```python
from dnakit.alignment import edlib_nice_alignment

result = edlib_nice_alignment("ACGTTGCA", "GGACGTAGCAGG")
print(result["edit_distance"], result["cigar"])
print(result["query_aligned"])
print(result["matched_aligned"])
print(result["target_aligned"])
```

- **Example output：**

```text
5 2D4=1X3=2D
--ACGTTGCA--
--||||.|||--
GGACGTAGCAGG
```

## 14) `dnakit.alignment.pyhmmer_align_to_hmm()`

- **Purpose：** Align named DNA sequences to one DNA profile HMM and return one aligned row per sequence, in model columns.
- **Method：** Reads the single DNA profile in a HMMER3 file, digitises the sequences with its alphabet and calls `pyhmmer.hmmalign`. Upper-case letters are match columns, lower-case letters are insert residues outside the model, and `.` is a gap; `trim=True` keeps only model columns. Protein profiles raise an error.
- **API：** `dnakit.alignment.pyhmmer_align_to_hmm(profile, sequences, *, trim=False, all_consensus_cols=True, max_sequences=10000, max_length=100000, verbose=True)`.
- **Input：** `profile`: path to a HMMER3 file with one DNA profile; `sequences`: a DNA object or a name-to-sequence mapping; `trim`: keep only model columns; `all_consensus_cols`: keep all consensus columns; `max_sequences` and `max_length`: limits.
- **Output：** A dict with `profile`, `profile_length`, `consensus`, `sequence_count`, `alignment_columns`, `rows` (each with `name` and `aligned`) and `provenance`.
- **Example code：**

```python
from dnakit.alignment.backends import build_dna_hmm
from dnakit.alignment import pyhmmer_align_to_hmm

build_dna_hmm("aligned.fa", name="probe", output_path="probe.hmm", overwrite=True)
result = pyhmmer_align_to_hmm("probe.hmm", {"s1": "GGACGTACGTACGG", "s2": "ACGTTCGTAC"}, trim=True)
print(result["profile"], result["profile_length"], result["consensus"])
for row in result["rows"]:
    print(row["name"], row["aligned"])
```

- **Example output：**

```text
probe 10 acgtacgtac
s1 ACGTACGTAC
s2 ACGTTCGTAC
```

## Installing the optional dependencies

`parasail`, `edlib` and `pyhmmer` are optional. Each is imported only when its function runs, and a missing package raises an error with the install command: `dnakit[bio-parasail]`, `dnakit[bio-edlib]`, `dnakit[bio-pyhmmer]`.

## Verification

`others/tests/unit/test_full_G2.py` checks each wrapper against a direct upstream call or a hand calculation. No function here needs the network.

## Scope

- The semi-global families expose the `trace` output only. The `table` and `rowcol` dynamic-programming matrices are debugging output and are not wrapped here.
- The `mappy.Aligner` mapping and indexing API is not in the inventory; `map_query_sequences_paf()` and related functions in `backends.py` already cover it.
- In the semi-global families, the two begins (or the two ends) of a pair cannot both be skipped in one alignment: the dynamic program starts on the first row or column and ends on the last row or column. `qb_db` and `qe_de` therefore skip one side only.
- SIMD variants are chosen with `method` and `width`; the profile-input forms (`*_profile_*`) are not exposed.
