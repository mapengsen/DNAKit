# OPS-007 Trimming

Trimming removes specified lengths from the left and right ends of a DNA sequence, returning a new, shortened sequence.

- **Function:** Identify and delete adapters, primers or low-quality bases from both ends of the sequence, and return retained fragments and trimming coordinates for preprocessing of reads and amplification products.
- **API:** `dnakit.ops.trim(dna[required], left[optional], right[optional], feature_policy[optional], letter_annotation_policy[optional])`.
- **Input:** `left` is the length to remove from the left end, `right` is the length to remove from the right end; both are non-negative integers.
- **Returns:** New `DNA`, which does not modify the input in place.

<span id="_1"></span>**Differences from arbitrary interval operations**

| Requirements | APIs that should be used |
|------------------------| ----------------------------------|
| Delete bases from the left and right ends of the sequence | `trim(dna, left=..., right=...)` |
| Delete any range in the sequence | `delete(dna, start, end)` |
| Only keep any interval in the sequence | `subsequence(dna, start, end)` |

`trim()` refers specifically to trimming at both ends; deletion of any internal positions is sequence editing and is handled by `delete()`.

<span id="_2"></span>**Example**

```python
from dnakit import DNA
from dnakit.ops import trim

dna = DNA("AACCGGTT", id="seq-1")
trimmed = trim(dna, left=1, right=2)
print(trimmed.symbols)
```

```text
ACCGG
```
