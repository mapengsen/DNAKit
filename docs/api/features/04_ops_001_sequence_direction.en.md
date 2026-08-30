# OPS-001 sequence completion

Reverse, complement, or reverse complement the DNA sequence to obtain a new sequence after the direction is reversed.

## 1) OPS-001.1 reverse order

- **Function:** Only reverse the character order of the DNA sequence without replacing the bases, used to check or construct the reverse position sequence.
- **API:** `dnakit.ops.reverse(sequence[required])`.
- **Input:** One required `DNASequence`; no additional required parameters.
- **Sample code:**

```python
from dnakit import DNASequence
from dnakit.ops import reverse

seq = DNASequence("AATCGC")
print(reverse(seq).symbols)
```

- **Example results:**

```text
CGCTAA
```

## 2) OPS-001.2 complementary

- **Function:** Replace each DNA base with a complementary base according to base pairing rules, but maintain the original position order, and use it to construct the corresponding complementary chain.
- **API:** `dnakit.ops.complement(sequence[required])`.
- **Input:** One required `DNASequence`; no additional required parameters.
- **Sample code:**

```python
from dnakit import DNASequence
from dnakit.ops import complement

seq = DNASequence("AATCGC")
print(complement(seq).symbols)
```

- **Example results:**

```text
TTAGCG
```

## 3) OPS-001.3 Reverse complementation

- **Function:** Simultaneously reverse the sequence order and replace it with complementary bases to obtain a 5′→3′ representation of the other strand, which can be used for reverse strand search, direction unification and primer analysis.
- **API:** `dnakit.ops.reverse_complement(dna[required], feature_policy[optional])`.
- **Input:** Normal user leaflet record `DNA`; will synchronize feature, strand and base-by-base annotations and return new `DNA`.
- **Sample code:**

```python
from dnakit import DNA
from dnakit.ops import reverse_complement

dna = DNA("AATCGC", id="seq-1")
print(reverse_complement(dna).symbols)
```

- **Example results:**

```text
GCGATT
```
