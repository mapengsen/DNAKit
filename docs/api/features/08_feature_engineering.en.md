# DNAFingerprint

DNAKit currently offers two methods of DNA fingerprint calculation. Each bit has only `0` or `1`, which can be directly used for Tanimoto/Jaccard similarity, retrieval, clustering and deduplication.

<a id="1-hashed-k-mer"></a>

## 1) Hashed k-mer bit fingerprint

- **Purpose:** Extracts k-mers occurring in a sequence and stably hashes them into fixed-length 0/1 bit vectors for use in Tanimoto/Jaccard similarity, clustering and retrieval, while avoiding excessively large explicit `4^k` feature spaces.
- **Method:** Hash k-mers into n_bits positions and set observed positions to 1; collisions can map distinct k-mers to one bit.
- **API:** `dnakit.fingerprints.hashed_kmer_fingerprint(value[required], k[required], n_bits[optional], canonical[optional], seed[optional], representation[optional], ambiguity_policy[optional], overlapping[optional], cross_gaps[optional])`
- **Input:** A DNA sequence and k; optional n_bits, seed and canonical settings.
- **Output:** Returns a `BitFingerprintResult`; the example shows `dimension`, `values`.
- **Common parameters:** Default `n_bits=2048`, `canonical=True`, `seed=0`; `representation="sparse"` only saves bits with value `1`.

- **Example code:**

```python
from dnakit import DNASequence
from dnakit.fingerprints import hashed_kmer_fingerprint

result = hashed_kmer_fingerprint(
    DNASequence("ACGTAC"),
    k=3,
    n_bits=16,
    seed=7,
    representation="sparse",
)

print(result.dimension)
print(dict(result.values))
```

- **Example output:**

```text
16
{'bit:1': 1, 'bit:10': 1}
```

## 2) Panel existence fingerprint

- **Purpose:** Convert the user-defined naming pattern panel into an interpretable 0/1 bit vector. Each bit indicates whether the corresponding motif or recognition sequence exists, making it easy to compare and filter sequences according to known functional patterns.
- **Method:** Scan named patterns in fixed panel order; set each bit to 1 when at least one match exists, otherwise 0.
- **API:** `dnakit.fingerprints.panel_fingerprint(value[required], panel[required], mode[optional], overlapping[optional], representation[optional], max_panel_size[optional], max_matches_per_pattern[optional])`
- **Input:** A DNA sequence and a named pattern panel with fixed order.
- **Output:** Returns a `BitFingerprintResult`; the example shows `feature_names`, `dense_values`.
- **Mode:** Default is `mode="iupac"`, `mode="exact"` can also be used; by default, both positive and negative chains are scanned.

- **Example code:**

```python
from dnakit import DNASequence
from dnakit.fingerprints import panel_fingerprint

panel = {
    "start": "ATG",
    "EcoRI": "GAATTC",
    "TATA_box": "TATAWAWR",
}
result = panel_fingerprint(DNASequence("ATGGAATTC"), panel)

print(result.feature_names)
print(result.dense_values())
```

- **Example output:**

```text
('panel:EcoRI', 'panel:TATA_box', 'panel:start')
(1, 0, 1)
```

- **Explanation:** The corresponding modes of `EcoRI` and `start` exist, so the corresponding bits are `1`; `TATA_box` does not exist, so it is `0`.

[FAQ: details](../../faq/tutorial_basics.md#fingerprint-terms)
