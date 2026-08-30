# DNAFingerprint

DNAKit currently offers two methods of DNA fingerprint calculation. Each bit has only `0` or `1`, which can be directly used for Tanimoto/Jaccard similarity, retrieval, clustering and deduplication.

## 1) Hashed k-mer bit fingerprint

- **What it does:** Extracts k-mers occurring in a sequence and stably hashes them into fixed-length 0/1 bit vectors for use in Tanimoto/Jaccard similarity, clustering and retrieval, while avoiding excessively large explicit `4^k` feature spaces.
- **Features:** Functionally closest to the Morgan/ECFP bit fingerprint in the molecule; no matter how big `k` is, the output length is always `n_bits`.
- **API:** `dnakit.fingerprints.hashed_kmer_fingerprint(value[required], k[required], n_bits[optional], canonical[optional], seed[optional], representation[optional], ambiguity_policy[optional], overlapping[optional], cross_gaps[optional])`
- **Common parameters:** Default `n_bits=2048`, `canonical=True`, `seed=0`; `representation="sparse"` only saves bits with value `1`.

**Example:**

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

**Output:**

```text
16
{'bit:1': 1, 'bit:10': 1}
```

## 2) Panel existence fingerprint

- **Function:** Convert the user-defined naming pattern panel into an interpretable 0/1 bit vector. Each bit indicates whether the corresponding motif or recognition sequence exists, making it easy to compare and filter sequences according to known functional patterns.
- **Features:** Similar to MACCS Keys in molecules, the bits have clear meanings and are suitable for detecting known motifs, promoter patterns or enzyme digestion recognition sequences.
- **API:** `dnakit.fingerprints.panel_fingerprint(value[required], panel[required], mode[optional], overlapping[optional], representation[optional], max_panel_size[optional], max_matches_per_pattern[optional])`
- **Mode:** Default is `mode="iupac"`, `mode="exact"` can also be used; by default, both positive and negative chains are scanned.

**Example:**

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

Output:

```text
('panel:EcoRI', 'panel:TATA_box', 'panel:start')
(1, 0, 1)
```

- **Explanation:** The corresponding modes of `EcoRI` and `start` exist, so the corresponding bits are `1`; `TATA_box` does not exist, so it is `0`.
