# Transcription and Translation

Transcribe DNA sequences into RNA or translate into protein sequences according to the specified reading frame and code table.

## 1) Transcription {#1-ops-0021-transcription}

- **Purpose:** Convert thymine `T` in DNA to uracil `U` according to the specified chain direction to obtain an uppercase RNA sequence for checking transcript products.
- **Method:** Replace T with U; reverse-strand transcription first takes the reverse complement.
- **API:** `dnakit.ops.transcribe(sequence[required], strand[optional])`.
- **Input:** One `DNASequence` is required; `strand` optional `"forward"` or `"reverse"`, forward chain is used by default.
- **Output:** Returns the transcribed RNA string.
- **Example code:**

```python
from dnakit import DNASequence
from dnakit.ops import transcribe

seq = DNASequence("ATGGGCTAA")
print(transcribe(seq))
```

- **Example output:**

```text
AUGGGCUAA
```

## 2) Translation {#2-ops-0022-translation}

- **Purpose:** Translate DNA or RNA into protein sequences according to the specified chain direction, reading frame and genetic code table, used to check coding regions and amino acid sequences.
- **Method:** Read codons in the selected frame and code table, applying stop and ambiguity policies.
- **API:** `dnakit.ops.translate(sequence[required], frame[optional], table[optional], strand[optional], stop_policy[optional], unknown_policy[optional], incomplete_policy[optional])`.
- **Input:** Required `DNASequence` or normalized uppercase DNA/RNA string; optional `strand`, `frame`, `table`, stop codon policy, ambiguous codon policy, and incomplete codon policy.
- **Output:** Returns the translated amino-acid sequence as a string.
- **Example code:**

```python
from dnakit import DNASequence
from dnakit.ops import translate

seq = DNASequence("ATGGGCTAA")
print(translate(seq))
```

- **Example output:**

```text
MG*
```

<!-- BEGIN GENERATED API: translation -->

<span id="translation"></span>

## 3) Translate forward reading frames and find the longest ORF {#dnakit.databases.dna_translate_reading_frames}

- **Purpose:** Translate forward reading frames and find the longest ORF.
- **Method:** Translate the three forward frames with the standard genetic code and compare uninterrupted protein segments.
- **API：** `dnakit.databases.dna_translate_reading_frames(*, operation=None, sequence, frame=None, genetic_code=None, verbose=True)`
- **Input:** Provide `sequence`.
- **Output:** Returns a mapping with `sequence_length`, `frames`, `best_frame`, `best_orf` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.dna_translate_reading_frames(
    operation='translate_reading_frames',
    sequence='ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "sequence_length": "<integer>",
  "frames": "<object>",
  "best_frame": "<integer>",
  "best_orf": "<string>",
  "best_orf_length_aa": "<integer>"
}
```

[Full parameters and usage notes](../../faq/api/04_ops_002_transcription_translation.md#dnakit.databases.dna_translate_reading_frames-parameters)

<!-- END GENERATED API: translation -->
