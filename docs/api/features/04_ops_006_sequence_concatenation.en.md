# OPS-006 Sequence splicing

Splice multiple DNA sequences or perform de-overlapping splicing based on precise overlapping regions.

## 1. Ordinary sequence splicing

- **Function:** Connect multiple DNA fragments in the input order, insert linker or structured gap at the junction, and return the coordinates of each fragment in the new sequence, suitable for constructing combined sequences.
- **API:** `dnakit.ops.concat(sequences[required], linker[optional], gap[optional])`.
- **Input:** At least two `DNASequence` or normalized strings; `linker` and `gap` cannot be used at the same time.
- **Sample code:**

```python
from dnakit import DNASequence
from dnakit.ops import concat

joined = concat(
    [DNASequence("AA"), DNASequence("CC")],
    linker="T",
)
print(joined.symbols)
```

- **Example results:**

```text
AATCC
```

## 2. Remove overlapping splicing

- **Function:** Complete the connection based on the overlapping relationship between the ends of two pieces of DNA, so that the matched overlapping region is retained only once, and the joint position and final sequence are reported.
- **API:** `dnakit.ops.concat_overlap(sequences[required], min_overlap[optional], max_overlap[optional])`.
- **Input:** Must be two linear `DNASequence`s with no explicit gap or a normalized string.
- **Rule:** Automatically select the longest "left fragment suffix = right fragment prefix" exact overlap; `min_overlap` defaults to 1, `max_overlap` can limit the search length. IUPAC characters are compared literally and no compatibility inferences are made.
- **Sample code:**

```python
from dnakit import DNASequence
from dnakit.ops import concat_overlap

joined = concat_overlap(
    [DNASequence("AAACCC"), DNASequence("CCCGG")],
)
print(joined.symbols)
```

- **Example results:**

```text
AAACCCGG
```
