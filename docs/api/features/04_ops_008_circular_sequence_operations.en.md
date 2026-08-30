# OPS-008 Circular sequence operation

Rotate the origin of a circular DNA, or extract a sequence fragment that spans the circular origin.

- **Function:** Respecify the storage starting point without changing the content of the circular DNA, or intercept the fragments spanning the origin, and correctly convert the old and new coordinates.
- **Normal API:** `dnakit.ops.rotate(dna[required], offset[required])`, `dnakit.ops.canonical_origin(dna[required])`, both return new `DNA` and synchronize comments.
- **Advanced Serial Value API:** `dnakit.ops.circular_subsequence(sequence[required], start[required], end[required], allow_gaps[optional])`.
- **Input:** Required Single record `DNA` explicitly declared as circular; rotate requires an offset. When only extracting sequence values, pass `dna.sequence` to `circular_subsequence()`.
- **Sample code:**

```python
from dnakit import DNA
from dnakit.ops import canonical_origin, circular_subsequence, rotate

dna = DNA("GATTACA", topology="circular")
rotated = rotate(dna, 2)
canonical = canonical_origin(dna)
wrapped = circular_subsequence(dna.sequence, 5, 2)
print(rotated.symbols)
print(canonical.symbols)
print(wrapped.symbols)
```

- **Example results:**

```text
TTACAGA
ACAGATT
CAGA
```
