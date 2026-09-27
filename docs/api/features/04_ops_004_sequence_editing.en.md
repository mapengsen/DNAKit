# Sequence Editing

Insert, delete or replace DNA fragments, cut duplex DNA, trim sequence ends or mask selected regions, handling coordinates and annotations according to each operation.

## 1) Insert, delete and substitute sequences

- **Purpose:** Insert, delete or replace bases at the specified position, return the edited sequence, operation records and coordinate changes, and process the affected features synchronously.
- **Method:** Split at the requested position or interval, modify the fragment, and update affected coordinates.
- **API:** `dnakit.ops.insert(dna[required], position[required], fragment[required])`, `dnakit.ops.delete(dna[required], start[required], end[required])`, `dnakit.ops.substitute(dna[required], start[required], end[required], fragment[required])`.
- **Input:** A single-record `DNA` and a position/interval; insertion or replacement also requires new DNA fragment. All three functions return a new `DNA`.
- **Output:** All three functions return new edited `DNA`; read the sequence from `symbols`. The input is unchanged.
- **Example code:**

```python
from dnakit import DNA
from dnakit.ops import delete, insert, substitute

dna = DNA("AACCGG", id="seq-1")
inserted = insert(dna, 2, "TT")
deleted = delete(dna, 2, 4)
replaced = substitute(dna, 2, 4, "TN")
print(inserted.symbols)
print(deleted.symbols)
print(replaced.symbols)
```

- **Example output:**

```text
AATTCCGG
AAGG
AATNGG
```

<!-- BEGIN GENERATED API: digest -->

<span id="digest"></span>

## 2) Simulate restriction digestion and resulting fragments {#dnakit.databases.dna_virtual_digest}

- **Purpose:** Simulate restriction digestion and resulting fragments.
- **Method:** Match restriction sites and split linear or circular sequences at their cleavage offsets.
- **API：** `dnakit.databases.dna_virtual_digest(*, operation=None, sequence, enzymes=None, circular=None, verbose=True)`
- **Input:** Provide `sequence`.
- **Output:** Returns a mapping with `fragments`, `n_fragments`, `enzymes_used`, `cut_sites` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.dna_virtual_digest(
    operation='virtual_digest',
    sequence='ATGGAGGAGCCGCAGTCAGATCCTAGCGTTGAATTCGGATCCAAGCTTACGTACGT',
    enzymes=['EcoRI', 'BamHI', 'HindIII'],
    circular=False,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "fragments": [
    {
      "sequence": "<string>",
      "length": "<integer>",
      "start": "<integer>",
      "end": "<integer>"
    }
  ],
  "n_fragments": "<integer>",
  "enzymes_used": [
    "<string>"
  ],
  "cut_sites": [
    {
      "enzyme": "<string>",
      "position": "<integer>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/04_ops_004_sequence_editing.md#dnakit.databases.dna_virtual_digest-parameters)

<!-- END GENERATED API: digest -->

## 3) Cut duplex DNA at strand-specific positions {#bio-52-cut}

- **Purpose:** Cut duplex DNA and inspect blunt or sticky fragment ends.
- **Method:** Store both 5′→3′ strands and their offset, then create fragments at the specified cut positions.
- **API:** `dnakit.bioinformatics.cloning.Duplex.from_sequence(sequence, circular=False)`; `Duplex.cut(watson_position, crick_position)`
- **Input:** A linear complementary duplex or A/C/G/T sequence; both cut positions use left-to-right Watson coordinates and must lie inside the strands. Nicks, internal gaps and mismatches are unsupported.
- **Output:** Returns two cut `Duplex` fragments and their end information.
- **Example code:**

```python
from dnakit.bioinformatics.cloning import Duplex

molecule = Duplex.from_sequence("AACCGGTT")
left, right = molecule.cut(3, 5)
print(left.ends)
```

- **Example output:**

```text
{'left': {'polarity': 'blunt', 'sequence': ''}, 'right': {'polarity': '5prime', 'sequence': 'CG'}}
```

[Dependencies and implementation notes](../../faq/bioinformatics.md#optimization).

## 4) Trimming: trim sequence ends {#trimming}

Trimming removes specified lengths from the left and right ends of a DNA sequence, returning a new, shortened sequence.

- **Purpose:** Remove specified lengths from both ends and update annotations; the caller supplies the lengths to trim.
- **Method:** Remove left and right bases from the ends, retain the middle segment and transform annotations according to policy.
- **API:** `dnakit.ops.trim(dna[required], left[optional], right[optional], feature_policy[optional], letter_annotation_policy[optional])`.
- **Input:** `left` is the length to remove from the left end, `right` is the length to remove from the right end; both are non-negative integers.
- **Output:** Returns a new trimmed `DNA` without modifying the input in place.
- **Example code:**

```python
from dnakit import DNA
from dnakit.ops import trim

dna = DNA("AACCGGTT", id="seq-1")
trimmed = trim(dna, left=1, right=2)
print(trimmed.symbols)
```

- **Example output:**

```text
ACCGG
```

[FAQ: trim, delete and subsequence](../../faq/tutorial_basics.md#trimming-comparison)

## 5) Masking: mask sequence regions {#masking}

Masking replaces bases in one or more selected DNA sequence regions with a specified character (`N` by default), preserving the original sequence length and coordinates.

- **Purpose:** Mask low-quality, low-complexity, or specified regions without changing sequence length and coordinates to prevent these positions from affecting search, statistics, or modeling.
- **Method:** Replace bases inside each half-open interval with symbol while preserving length and applying annotation policies.
- **API:** `dnakit.ops.mask(dna[required], intervals[required], symbol[optional], feature_policy[optional], letter_annotation_policy[optional])`.
- **Input:** `intervals` is one or more 0-based half-open intervals `(start, end)`; `symbol` defaults to `N`.
- **Output:** Returns a new masked `DNA` of the same length without modifying the input in place.
- **Example code:**

```python
from dnakit import DNA
from dnakit.ops import mask

dna = DNA("AACCGGTT", id="seq-1")
masked = mask(dna, [(2, 4), (6, 8)])
print(masked.symbols)
```

- **Example output:**

```text
AANNGGNN
```
