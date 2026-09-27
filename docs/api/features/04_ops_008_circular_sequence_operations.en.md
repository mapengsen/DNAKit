# Circular sequence operation

Rotate the origin of a circular DNA, or extract a sequence fragment that spans the circular origin.

## 1) Rotate circular sequences and extract origin-spanning regions

- **Purpose:** Respecify the storage starting point without changing the content of the circular DNA, or intercept the fragments spanning the origin, and correctly convert the old and new coordinates.
- **Method:** Rotate the stored origin by the offset; canonical_origin selects a canonical origin and transforms annotations.
- **API:** `dnakit.ops.rotate(dna[required], offset[required])`, `dnakit.ops.canonical_origin(dna[required])`, both return new `DNA` and synchronize comments.
- **Input:** Required Single record `DNA` explicitly declared as circular; rotate requires an offset. When only extracting sequence values, pass `dna.sequence` to `circular_subsequence()`.
- **Output:** Returns the rotated, canonicalized, or origin-spanning extracted `DNA` sequence.
- **Example code:**

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

- **Example output:**

```text
TTACAGA
ACAGATT
CAGA
```

[FAQ: details](../../faq/tutorial_basics.md#compatibility-04_ops_008_circular_sequence_operations)
