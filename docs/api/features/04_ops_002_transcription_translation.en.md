# OPS-002 Transcription and Translation

Transcribe DNA sequences into RNA or translate into protein sequences according to the specified reading frame and code table.

## 1) OPS-002.1 Transcription

- **Function:** Convert thymine `T` in DNA to uracil `U` according to the specified chain direction to obtain an uppercase RNA sequence for checking transcript products.
- **API:** `dnakit.ops.transcribe(sequence[required], strand[optional])`.
- **Input:** One `DNASequence` is required; `strand` optional `"forward"` or `"reverse"`, forward chain is used by default.
- **Sample code:**

```python
from dnakit import DNASequence
from dnakit.ops import transcribe

seq = DNASequence("ATGGGCTAA")
print(transcribe(seq))
```

- **Example results:**

```text
AUGGGCUAA
```

## 2) OPS-002.2 Translation

- **Function:** Translate DNA or RNA into protein sequences according to the specified chain direction, reading frame and genetic code table, used to check coding regions and amino acid sequences.
- **API:** `dnakit.ops.translate(sequence[required], frame[optional], table[optional], strand[optional], stop_policy[optional], unknown_policy[optional], incomplete_policy[optional])`.
- **Input:** Required `DNASequence` or normalized uppercase DNA/RNA string; optional `strand`, `frame`, `table`, stop codon policy, ambiguous codon policy, and incomplete codon policy.
- **Sample code:**

```python
from dnakit import DNASequence
from dnakit.ops import translate

seq = DNASequence("ATGGGCTAA")
print(translate(seq))
```

- **Example results:**

```text
MG*
```
