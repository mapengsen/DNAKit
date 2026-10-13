# DnaChisel, DnaCauldron and pydna sequence functions

This group wraps module-level sequence functions from DnaChisel, DnaCauldron and pydna as `dnakit.molecular_lab` functions. Each function imports its package only when it runs; a missing package raises `BackendUnavailableError` with the `dnakit[bio-dnachisel]`, `dnakit[bio-dnacauldron]` or `dnakit[bio-assembly]` install hint.

The functions only process sequences or files passed explicitly. `fetch_genbank_record()` needs network access, so its example uses a stand-in in place of the download. The example code imports from the module files; after the functions are exported from `dnakit.molecular_lab`, they can also be imported from the package. Functions that already have wrappers (such as `reverse_translate_protein()`, `design_pcr_primers()`, `find_recombination_sites()` and `optimize_dnachisel_record()`) are not repeated here.

## 1) `dnakit.molecular_lab.generate_constrained_dna_sequence()`

- **Purpose:** Generate a random DNA sequence that satisfies DnaChisel constraints.
- **Method:** Calls DnaChisel's `random_compatible_dna_sequence`, retrying until the constraints hold, then checks every constraint with `DnaOptimizationProblem`. Constraint names are limited to `AvoidPattern`, `AvoidStopCodons`, `AvoidHairpins` and `EnforceGCContent`. A given `seed` makes the result reproducible.
- **API:** `dnakit.molecular_lab.generate_constrained_dna_sequence(sequence_length, *, constraints=(), seed=None, max_random_iters=5000, verbose=True)`
- **Input:** An integer length from 1 to 20000; `constraints` is a list of mappings, each with a `name` and that constraint's arguments, for example `{"name": "AvoidPattern", "pattern": "GGTCTC"}`.
- **Output:** `sequence`, `sequence_length`, `seed`, `constraint_names`, `all_constraints_pass` and `provenance`.
- **Example code:**

```python
from dnakit.molecular_lab import generate_constrained_dna_sequence

result = generate_constrained_dna_sequence(
    30,
    constraints=[{"name": "AvoidPattern", "pattern": "GGTCTC"}],
    seed=1,
)
print(result["sequence"])
print(result["all_constraints_pass"])
```

- **Example output:**

```text
TGAAGTGTGAATAGTACTCACTCAGACATC
True
```

## 2) `dnakit.molecular_lab.find_sequence_difference_segments()`

- **Purpose:** Find the contiguous segments where two equal-length DNA sequences differ.
- **Method:** Calls DnaChisel's `sequences_differences_segments` for a position-by-position comparison and merges adjacent differing positions into one segment. No alignment is made, and insertions and deletions are not handled.
- **API:** `dnakit.molecular_lab.find_sequence_difference_segments(sequence_a, sequence_b, *, verbose=True)`
- **Input:** Two equal-length sequences of A, C, G and T; whitespace is removed and case is ignored.
- **Output:** `segments` (a list of `[start, end)` pairs), `segment_count`, `differing_positions`, `sequence_length` and `provenance`.
- **Example code:**

```python
from dnakit.molecular_lab import find_sequence_difference_segments

result = find_sequence_difference_segments("ACGTTACGTA", "ACGAAACGTC")
print(result["segments"])
print(result["differing_positions"])
```

- **Example output:**

```text
[[3, 5], [9, 10]]
3
```

## 3) `dnakit.molecular_lab.list_common_restriction_enzymes()`

- **Purpose:** List restriction enzymes with a given recognition-site length and enough suppliers from DnaChisel's bundled enzyme table.
- **Method:** Calls DnaChisel's `list_common_enzymes`, filters by `site_length` and `min_suppliers`, and sorts the names. No network access is used.
- **API:** `dnakit.molecular_lab.list_common_restriction_enzymes(*, site_length=6, min_suppliers=1, verbose=True)`
- **Input:** `site_length` is an integer from 4 to 12; `min_suppliers` is an integer from 0 to 20.
- **Output:** `enzyme_count`, the sorted `enzymes` list, plus `site_length`, `min_suppliers` and `provenance`.
- **Example code:**

```python
from dnakit.molecular_lab import list_common_restriction_enzymes

result = list_common_restriction_enzymes(site_length=6, min_suppliers=5)
print(result["enzyme_count"])
print(result["enzymes"][:5])
```

- **Example output:**

```text
37
['AatII', 'AccI', 'ApaI', 'BalI', 'BamHI']
```

## 4) `dnakit.molecular_lab.group_nearby_positions()`

- **Purpose:** Group a list of positions by distance, for example to cluster nearby variants or mutation sites.
- **Method:** Calls DnaChisel's `group_nearby_indices`. Neighbouring positions join a group while their difference is below `max_gap`; `max_group_spread` limits the distance from a group's first to its last position. Either rule can be omitted (`None` means no limit). Output is sorted by position.
- **API:** `dnakit.molecular_lab.group_nearby_positions(positions, *, max_gap=None, max_group_spread=None, verbose=True)`
- **Input:** A list of unique non-negative integers, at most 100000; `max_gap` and `max_group_spread` are non-negative integers or `None`.
- **Output:** `groups` (each group sorted), `group_count`, `position_count`, `max_gap`, `max_group_spread` and `provenance`.
- **Example code:**

```python
from dnakit.molecular_lab import group_nearby_positions

result = group_nearby_positions([30, 1, 11, 2, 10, 3], max_gap=2)
print(result["groups"])
```

- **Example output:**

```text
[[1, 2, 3], [10, 11], [30]]
```

## 5) `dnakit.molecular_lab.autoselect_golden_gate_enzyme()`

- **Purpose:** Decide which Type IIS enzyme a set of Golden Gate parts was designed for.
- **Method:** Each sequence becomes a linear record that is passed to DnaCauldron's `autoselect_enzyme`, which picks the candidate whose cut pattern is closest to two sites per part. If no candidate fits, DnaCauldron still returns its best candidate, so callers should check the cut sites themselves.
- **API:** `dnakit.molecular_lab.autoselect_golden_gate_enzyme(parts, *, enzymes=("BsmBI", "BsaI", "BbsI", "AarI", "SapI"), verbose=True)`
- **Input:** 1 to 200 DNA strings, each flanked by enzyme sites; `enzymes` is a subset of candidate names without repeats.
- **Output:** `enzyme`, `candidate_enzymes`, `part_count` and `provenance`.
- **Example code:**

```python
parts = [
    "GGTCTCAATGCTTTAAACCCGGGTTTTGCTTGGAGACC",
    "GGTCTCAGCTTAAAACCCCGGGGAAAACGTAGGAGACC",
    "GGTCTCACGTACCCCAAAATTTTGGGGATGCGGAGACC",
]

from dnakit.molecular_lab import autoselect_golden_gate_enzyme

print(autoselect_golden_gate_enzyme(parts)["enzyme"])
```

- **Example output:**

```text
BsaI
```

## 6) `dnakit.molecular_lab.list_type2s_digestion_overhangs()`

- **Purpose:** List the sticky-end overhangs that Type IIS digestion creates.
- **Method:** Calls DnaCauldron's `list_digestion_overhangs`. The result is unique and sorted; DnaCauldron's `"None"` placeholder for missing ends is removed and counted separately. With `parts_only` true, only overhangs on part fragments are kept.
- **API:** `dnakit.molecular_lab.list_type2s_digestion_overhangs(sequences, *, enzyme="auto", parts_only=True, verbose=True)`
- **Input:** 1 to 200 DNA strings; `enzyme` is `"auto"` or one of `BsmBI`, `BsaI`, `BbsI`, `AarI`, `SapI`.
- **Output:** `overhangs`, `overhang_count`, `none_placeholder_count`, `enzyme`, `parts_only`, `sequence_count` and `provenance`.
- **Example code:**

```python
parts = [
    "GGTCTCAATGCTTTAAACCCGGGTTTTGCTTGGAGACC",
    "GGTCTCAGCTTAAAACCCCGGGGAAAACGTAGGAGACC",
    "GGTCTCACGTACCCCAAAATTTTGGGGATGCGGAGACC",
]

from dnakit.molecular_lab import list_type2s_digestion_overhangs

result = list_type2s_digestion_overhangs(parts, enzyme="BsaI")
print(result["overhangs"])
print(result["overhang_count"])
```

- **Example output:**

```text
['ATGC', 'CGTA', 'GCTT']
3
```

## 7) `dnakit.molecular_lab.crop_dna_record_with_features()`

- **Purpose:** Crop a DNA segment and keep the features that overlap the cropped window.
- **Method:** DnaCauldron's `crop_record_with_saddling_features` decides the cropped sequence and which features are kept. DnaCauldron 2.0.12 leaves partly included features in absolute coordinates, so this function clips and shifts those coordinates to the window. Features outside the window are dropped.
- **API:** `dnakit.molecular_lab.crop_dna_record_with_features(sequence, *, start, end, features=(), verbose=True)`
- **Input:** A DNA string; `start` and `end` form a 0-based half-open window; `features` is a list of mappings with `type`, `start`, `end` and an optional `label`.
- **Output:** `sequence`, `sequence_length`, `features` (coordinates relative to the new sequence), `feature_count`, `source_length` and `provenance`.
- **Example code:**

```python
from dnakit.molecular_lab import crop_dna_record_with_features

result = crop_dna_record_with_features(
    "ACGT" * 50,
    start=10,
    end=120,
    features=[
        {"type": "CDS", "start": 0, "end": 40, "label": "a"},
        {"type": "misc_feature", "start": 60, "end": 100},
        {"type": "misc_feature", "start": 150, "end": 180},
    ],
)
print(result["sequence_length"])
print(result["features"])
```

- **Example output:**

```text
110
[{'type': 'CDS', 'start': 0, 'end': 30, 'label': 'a'}, {'type': 'misc_feature', 'start': 50, 'end': 90, 'label': ''}]
```

## 8) `dnakit.molecular_lab.simulate_pcr_amplicon()`

- **Purpose:** Simulate the PCR product of a linear template with a primer pair.
- **Method:** Calls pydna's `amplify.pcr`. Each primer must anneal over at least `min_annealing` bases. The product runs from the 5' end of the forward primer to the 5' end of the reverse primer. A template without a product raises an error.
- **API:** `dnakit.molecular_lab.simulate_pcr_amplicon(template, forward_primer, reverse_primer, *, min_annealing=14, verbose=True)`
- **Input:** Template DNA (at most 20000 bp); two primer sequences (at most 200 nt each); `min_annealing` is an integer from 8 to 60.
- **Output:** `product_sequence`, `product_length`, `template_length`, `forward_primer`, `reverse_primer` and `provenance`.
- **Example code:**

```python
from dnakit.molecular_lab import simulate_pcr_amplicon

template = "AAAAA" + "ACGTTGCAAGGCTTAGCAT" + "C" * 16 + "GCAAGGCTTAGCATTGGACC" + "TTTTT"
result = simulate_pcr_amplicon(template, "ACGTTGCAAGGCTTAGCAT", "GGTCCAATGCTAAGCCTTGC")
print(result["product_length"])
print(result["product_sequence"])
```

- **Example output:**

```text
55
ACGTTGCAAGGCTTAGCATCCCCCCCCCCCCCCCCGCAAGGCTTAGCATTGGACC
```

## 9) `dnakit.molecular_lab.assemble_pydna_fragments()`

- **Purpose:** Enumerate the products that a set of DNA fragments can assemble into.
- **Method:** Calls pydna's `assembly2` algorithms. Overlap methods (`gibson`, `in_fusion`, `fusion_pcr`, `in_vivo`) search for end overlaps of at least `min_overlap`; enzyme methods (`golden_gate`, `restriction_ligation`) digest with the given `Bio.Restriction` enzymes and ligate. Products are sorted by length and sequence, and at most `max_products` are listed.
- **API:** `dnakit.molecular_lab.assemble_pydna_fragments(fragments, *, method="gibson", min_overlap=25, enzymes=(), circular_only=False, max_products=50, verbose=True)`
- **Input:** 2 to 20 DNA fragments, at most 200000 bp in total; enzyme methods need `enzymes`, for example `["BsaI"]`; overlap methods must leave `enzymes` empty.
- **Output:** `product_count`, `products` (each with sequence, length and circularity), `method`, `min_overlap`, `enzymes` and `provenance`.
- **Example code:**

```python
parts = [
    "GGTCTCAATGCTTTAAACCCGGGTTTTGCTTGGAGACC",
    "GGTCTCAGCTTAAAACCCCGGGGAAAACGTAGGAGACC",
    "GGTCTCACGTACCCCAAAATTTTGGGGATGCGGAGACC",
]

from dnakit.molecular_lab import assemble_pydna_fragments

result = assemble_pydna_fragments(
    parts, method="golden_gate", enzymes=["BsaI"], circular_only=True
)
print(result["product_count"], result["products"][0]["length"], result["products"][0]["circular"])
print(result["products"][0]["sequence"])
```

- **Example output:**

```text
1 60 True
ATGCTTTAAACCCGGGTTTTGCTTAAAACCCCGGGGAAAACGTACCCCAAAATTTTGGGG
```

## 10) `dnakit.molecular_lab.find_common_substrings()`

- **Purpose:** Find every common substring of at least `min_length` bases between two sequences.
- **Method:** Calls pydna's `common_sub_strings`. Each match is `(x_start, y_start, length)`; matches are sorted by position and at most 1000 are listed.
- **API:** `dnakit.molecular_lab.find_common_substrings(sequence_x, sequence_y, *, min_length=25, verbose=True)`
- **Input:** Two DNA sequences (at most 20000 bp each); `min_length` is an integer from 1 to 20000.
- **Output:** `matches`, `match_count`, `min_length` and `provenance`.
- **Example code:**

```python
from dnakit.molecular_lab import find_common_substrings

result = find_common_substrings("ACGTTGCAAGGCTTAGCAT", "TTGCAAGGCTTAGCATCC", min_length=5)
print(result["matches"])
```

- **Example output:**

```text
[{'x_start': 3, 'y_start': 0, 'length': 16}]
```

## 11) `dnakit.molecular_lab.find_terminal_overlaps()`

- **Purpose:** Find where the end of `sequence_x` matches the start of `sequence_y`, as used for Gibson-style joins.
- **Method:** Calls pydna's `common_sub_strings.terminal_overlap`. A returned match should sit at the 3' end of `sequence_x` and the 5' start of `sequence_y`.
- **API:** `dnakit.molecular_lab.find_terminal_overlaps(sequence_x, sequence_y, *, min_length=15, verbose=True)`
- **Input:** Two DNA sequences (at most 20000 bp each); `min_length` is an integer from 1 to 20000.
- **Output:** `matches`, `match_count`, `min_length` and `provenance`.
- **Example code:**

```python
from dnakit.molecular_lab import find_terminal_overlaps

result = find_terminal_overlaps("ACGTTGCAAGGCTTAGCAT", "TTGCAAGGCTTAGCATCC", min_length=5)
print(result["match_count"])
print(result["matches"])
```

- **Example output:**

```text
1
[{'x_start': 3, 'y_start': 0, 'length': 16}]
```

## 12) `dnakit.molecular_lab.primer_pair_hybridization()`

- **Purpose:** Compute the duplex products two primers can form with each other, for example to check for primer dimers.
- **Method:** Calls pydna's `oligonucleotide_hybridization` and `oligonucleotide_hybridization_overhangs`. Two primers pair when their shared complementary region is at least `minimal_annealing` bases; pairs with mismatches raise an error.
- **API:** `dnakit.molecular_lab.primer_pair_hybridization(forward_primer, reverse_primer, *, minimal_annealing=10, verbose=True)`
- **Input:** Two primer sequences written 5' to 3' (at most 200 nt each); `minimal_annealing` is an integer from 4 to 60.
- **Output:** `overhangs` (possible overlap lengths as returned by pydna, with pydna's sign convention), `product_count`, `products` and `provenance`.
- **Example code:**

```python
from dnakit.molecular_lab import primer_pair_hybridization

result = primer_pair_hybridization("GGGGGACGTTGCAAGGCTTAGCATG", "CATGCTAAGCCTTGCAACGTCCCCC")
print(result["overhangs"])
print(result["products"])
```

- **Example output:**

```text
[0]
[{'sequence': 'GGGGGACGTTGCAAGGCTTAGCATG', 'length': 25}]
```

## 13) `dnakit.molecular_lab.primer_melting_temperature()`

- **Purpose:** Compute a primer's melting temperature with pydna's default nearest-neighbour parameters.
- **Method:** Calls pydna's `tm.tm_default` in strict mode, which rejects characters other than A, C, G and T. The value is rounded to three decimals.
- **API:** `dnakit.molecular_lab.primer_melting_temperature(sequence, *, verbose=True)`
- **Input:** A DNA string of 6 to 200 bases.
- **Output:** `tm_celsius`, `length`, `sequence` and `provenance`.
- **Example code:**

```python
from dnakit.molecular_lab import primer_melting_temperature

result = primer_melting_temperature("ACGTTGCAAGGCTTAGCATG")
print(result["tm_celsius"])
```

- **Example output:**

```text
62.197
```

## 14) `dnakit.molecular_lab.parse_primer_sequences()`

- **Purpose:** Read primer sequences from FASTA-formatted text or a file.
- **Method:** Calls pydna's `parsers.parse_primers`. `source` is either an existing file path or the text itself; names come from FASTA headers and the input order is kept.
- **API:** `dnakit.molecular_lab.parse_primer_sequences(source, *, verbose=True)`
- **Input:** A file path or text; at most 1000 primers, and the file or text must be 1 MB or smaller.
- **Output:** `primers` (each with name, sequence and length), `primer_count`, `source` and `provenance`.
- **Example file:** `primers.fa`

```text
>p1
ACGTACGTACGT
>p2
GGGGCCCCAAAA
```
- **Example code:**

```python
from dnakit.molecular_lab import parse_primer_sequences

result = parse_primer_sequences("primers.fa")
print(result["primer_count"])
print(result["primers"])
```

- **Example output:**

```text
2
[{'name': 'p1', 'sequence': 'ACGTACGTACGT', 'length': 12}, {'name': 'p2', 'sequence': 'GGGGCCCCAAAA', 'length': 12}]
```

## 15) `dnakit.molecular_lab.fetch_genbank_record()`

- **Purpose:** Download one GenBank nucleotide record by accession, for gene- or sequence-level lookups.
- **Method:** Calls pydna's `genbank.genbank`, which downloads through NCBI Entrez. The accession must look like `CS570233.1` or `NC_000913.3`. `email` is the contact address NCBI asks for; it is passed to the client and not stored in the result. This function needs network access. The example replaces the download with a stand-in, and the default tests use monkeypatch.
- **API:** `dnakit.molecular_lab.fetch_genbank_record(accession, *, email, verbose=True)`
- **Input:** An accession string; `email` is a contact address containing `@`.
- **Output:** `accession`, `name`, `description`, `topology`, `length`, `sequence` and `provenance`.
- **Example code:**

```python
from pydna.dseqrecord import Dseqrecord
import pydna.genbank
from dnakit.molecular_lab import fetch_genbank_record


def offline_stand_in(accession, *args, email=None, **kwargs):
    return Dseqrecord("ACGTACGTAA", name="DEMO")


# Replaces the NCBI download for this offline example only.
pydna.genbank.genbank = offline_stand_in
result = fetch_genbank_record("NC_000913.3", email="reader@example.org")
print(result["accession"], result["length"], result["topology"])
print(result["sequence"])
```

- **Example output:**

```text
NC_000913.3 10 linear
ACGTACGTAA
```
