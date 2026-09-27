# sequence completion

Reverse, complement, or reverse complement the DNA sequence to obtain a new sequence after the direction is reversed.

## 1) reverse order {#1-ops-0011-reverse-order}

- **Purpose:** Only reverse the character order of the DNA sequence without replacing the bases, used to check or construct the reverse position sequence.
- **Method:** Read bases from end to start without changing their symbols.
- **API:** `dnakit.ops.reverse(sequence[required])`.
- **Input:** One required `DNASequence`; no additional required parameters.
- **Output:** Returns the processed `DNASequence`; read its symbols from `symbols`.
- **Example code:**

```python
from dnakit import DNASequence
from dnakit.ops import reverse

seq = DNASequence("AATCGC")
print(reverse(seq).symbols)
```

- **Example output:**

```text
CGCTAA
```

## 2) complementary {#2-ops-0012-complementary}

- **Purpose:** Replace each DNA base with a complementary base according to base pairing rules, but maintain the original position order, and use it to construct the corresponding complementary chain.
- **Method:** Replace each base with its DNA complement without reordering positions.
- **API:** `dnakit.ops.complement(sequence[required])`.
- **Input:** One required `DNASequence`; no additional required parameters.
- **Output:** Returns the processed `DNASequence`; read its symbols from `symbols`.
- **Example code:**

```python
from dnakit import DNASequence
from dnakit.ops import complement

seq = DNASequence("AATCGC")
print(complement(seq).symbols)
```

- **Example output:**

```text
TTAGCG
```

## 3) Reverse complementation {#3-ops-0013-reverse-complementation}

- **Purpose:** Simultaneously reverse the sequence order and replace it with complementary bases to obtain a 5′→3′ representation of the other strand, which can be used for reverse strand search, direction unification and primer analysis.
- **Method:** Reverse and complement bases while remapping strand, features and per-base annotations.
- **API:** `dnakit.ops.reverse_complement(value, output_path=None, *, feature_policy="preserve", read_config=None, write_config=None, sequence_column=None, show_progress=True)`.
- **Input:** DNA text, a file path, or one `DNA`, `DNASequence`, or `DNARecord`. Text returns a string; files return `WriteResult`. Input and output paths must differ; existing outputs are protected. [FAQ: details](../../faq/tutorial_integrations.md#reverse-complement).
- **Output:** Returns a reverse-complemented `DNASequence` or `DNA`; file input can be written to the selected path.
- **Example code:**

```python
from dnakit import DNA
from dnakit.ops import reverse_complement

dna = DNA("AATCGC", id="seq-1")
print(reverse_complement(dna).symbols)
```

- **Example output:**

```text
GCGATT
```

<!-- BEGIN GENERATED API: sequence-direction -->

<span id="sequence-direction"></span>

## 4) Compute a DNA reverse complement {#dnakit.databases.sequence_reverse_complement}

- **Purpose:** Compute a DNA reverse complement.
- **Method:** Complement A↔T, C↔G and N↔N, then reverse the sequence.
- **API：** `dnakit.databases.sequence_reverse_complement(*, operation=None, sequence, verbose=True)`
- **Input:** Provide `sequence`.
- **Output:** Returns a mapping with `data` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.sequence_reverse_complement(
    operation='reverse_complement',
    sequence='ATGCGATCG',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "data": {
    "original": "<string>",
    "reverse_complement": "<string>",
    "length": "<integer>"
  }
}
```

[Full parameters and usage notes](../../faq/api/04_ops_001_sequence_direction.md#dnakit.databases.sequence_reverse_complement-parameters)

<!-- END GENERATED API: sequence-direction -->

## 5) Reverse-complement duplex DNA {#bio-52-reverse}

- **Purpose:** Reverse the orientation of duplex DNA while preserving its duplex and end structure.
- **Method:** Swap the Watson and Crick strands and adjust their relative offset.
- **API:** `dnakit.bioinformatics.cloning.Duplex.from_sequence(sequence, circular=False)`; `Duplex.reverse_complement()`
- **Input:** A fully complementary duplex or A/C/G/T sequence; nicks, internal gaps and mismatches are unsupported.
- **Output:** Returns the reverse-complemented `Duplex` with both strands and end structures.
- **Example code:**

```python
from dnakit.bioinformatics.cloning import Duplex

molecule = Duplex.from_sequence("AACCGGTT")
print(molecule.reverse_complement().watson)
```

- **Example output:**

```text
AACCGGTT
```

The example sequence equals its reverse complement.

[Dependencies and implementation notes](../../faq/bioinformatics.md#optimization).
