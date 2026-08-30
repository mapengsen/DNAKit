# OPS-007 Masking

Masking uses specified characters to replace the sequence of one or more target intervals in the DNA sequence with any character N, retaining the original sequence length and coordinates.

- **Function:** Mask low-quality, low-complexity, or specified regions without changing sequence length and coordinates to prevent these positions from affecting search, statistics, or modeling.
- **API:** `dnakit.ops.mask(dna[required], intervals[required], symbol[optional], feature_policy[optional], letter_annotation_policy[optional])`.
- **Input:** `intervals` is one or more 0-based half-open intervals `(start, end)`; `symbol` defaults to `N`.
- **Returns:** A new `DNA` with the same length, without modifying the input in place.

<span id="_1"></span>**Example**

```python
from dnakit import DNA
from dnakit.ops import mask

dna = DNA("AACCGGTT", id="seq-1")
masked = mask(dna, [(2, 4), (6, 8)])
print(masked.symbols)
```

```text
AANNGGNN
```
