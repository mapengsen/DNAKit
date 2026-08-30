# OPS-004 Sequence Editing

Insert, delete or replace DNA fragments at specified positions to obtain edited sequences and record coordinate changes.

- **Function:** Insert, delete or replace bases at the specified position, return the edited sequence, operation records and coordinate changes, and process the affected features synchronously.
- **API:** `dnakit.ops.insert(dna[required], position[required], fragment[required])`, `dnakit.ops.delete(dna[required], start[required], end[required])`, `dnakit.ops.substitute(dna[required], start[required], end[required], fragment[required])`.
- **Input:** Regular user leaflet record `DNA` and position/interval; insertion or replacement also requires new DNA fragment. All three functions return a new `DNA`.
- **Sample code:**

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

- **Example results:**

```text
AATTCCGG
AAGG
AATNGG
```
