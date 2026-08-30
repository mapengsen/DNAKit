# OPS-003 subsequence extraction

Extracts a specified fragment from a DNA sequence.

- **Function:** Intercept DNA fragments according to the start and end coordinates and retain the source range. It also supports interception of circular sequences across origins for extracting genes, functional regions or specified windows.
- **API:** `dnakit.ops.subsequence(sequence[required], start[required], end[required], allow_gaps[optional])`, `dnakit.ops.circular_subsequence(sequence[required], start[required], end[required], allow_gaps[optional])`.
- **Input:** Required `DNASequence`, `start`, `end`; Optional `allow_gaps=True` Keep or truncate known gaps.
- **Sample code:**

```python
from dnakit import DNASequence
from dnakit.ops import subsequence

seq = DNASequence("AACCGG")
selected = subsequence(seq, 1, 5)
print(selected.symbols)
```

- **Example results:**

```text
ACCG
```
