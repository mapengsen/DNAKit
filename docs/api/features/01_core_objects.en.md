# 1. Core data object

Ordinary users only need to remember one core object: `dnakit.DNA(...)`. This entry is used to input one sequence or multiple sequences, ID, topology, metadata and feature, and the return type is always `DNA`.

```python
import dnakit

# One sequence
dna = dnakit.DNA("ACGT")

# One sequence with attached information
annotated = dnakit.DNA(
    "ACGT",
    id="seq-1",
    topology="circular",
    metadata={"species": "synthetic"},
    features=[{"type": "motif", "start": 1, "end": 3}],
)

# Multiple sequences, still returning DNA
dataset = dnakit.DNA(["ACGT", "TTAA"])
detailed = dnakit.DNA(
    [
        {"sequence": "ACGT", "id": "seq-1"},
        {"sequence": "TTAA", "id": "seq-2", "topology": "circular"},
    ]
)

print(type(dna).__name__, type(dataset).__name__)
print(dataset.ids, dataset[0].symbols)
```

```text
DNA DNA
('sequence_1', 'sequence_2') ACGT
```

There are only three usage rules: a string represents a sequence; a string list represents multiple sequences; when multiple sequences require different additional information, use a dictionary list containing `sequence`. `data[0]` and `data[1:3]` still return `DNA`. The following `DNASequence`, `DNARecord` and `DNASet` are internally clearly layered and compatible with old code objects, and do not need to be constructed separately for normal use.

## 1) CORE-001 DNA sequence object

- **Function:** Save DNA sequence values, alphabets, linear or circular states, single-stranded or double-stranded states, and Gap information as basic sequence objects for descriptors, searches, alignments, and other calculations.
- **Normal API:** `dnakit.DNA(data[required], alphabet[optional], topology[optional], strandedness[optional])`.
- **Advanced Compatibility:** `dnakit.DNASequence(parts[required], alphabet[optional], topology[optional], strandedness[optional])`.
- **Input:** `DNA` can directly receive raw strings and automatically normalize them; advanced `DNASequence` only receives normalized sequences.
- **Sample code:**

```python
import dnakit

seq = dnakit.DNA("ACGT", topology="linear", strandedness="double")
print(seq.symbols)
print(seq.symbol_length)
print(seq.record_count)
```

- **Example results:**

```text
ACGT
4
1
```

## 2) CORE-002 DNA recording object

- **Function:** Attach ID, description, functional area and sample information to a DNA sequence so that the calculation results can be traced back to the specific record and its source.
- **Normal API:** `dnakit.DNA(data[required], id[optional], description[optional], features[optional], metadata[optional], letter_annotations[optional])`.
- **Advanced Compatibility:** `dnakit.DNARecord(sequence[required], id[required], ...)`.
- **Input:** Normal entry only requires the sequence; the rest of the information is optional parameters in the same call. Automatically generated `sequence_1` when no ID is provided.
- **Sample code:**

```python
import dnakit

record = dnakit.DNA(
    "ACGT",
    id="seq-1",
    description="Example sequence",
    metadata={"species": "human"},
)
print(record.id, record.symbols, record.metadata["species"])
```

- **Example results:**

```text
seq-1 ACGT human
```

## 3) CORE-003 DNA dataset object

- **Function:** Use the same `DNA` object to manage one or more records in a fixed order, and supports subscript and slice selection, serving as a unified data input entrance for ordinary users.
- **Normal API:** `dnakit.DNA(data[required], name[optional], source[optional], version[optional], collection_metadata[optional])`; use `dnakit.read(..., mode="dna")` for file reading.
- **Advanced Compatibility:** `dnakit.DNASet(...)`, `DNASet.from_records(...)`, `DNASet.from_sequences(...)` and `read_set(...)`.
- **Input:** Use a string list for simple multiple sequences; use a dictionary list when each record requires different information.
- **Sample code:**

```python
import dnakit

dataset = dnakit.DNA(
    ["AC", "GT"],
    name="demo",
)
print(dataset.ids)
print(dataset[1].symbols)
```

- **Example results:**

```text
('sequence_1', 'sequence_2')
GT
```

## 4) CORE-004 characteristic object

- **What it does:** A feature is an additional annotation for a sequence in `DNA` that tags ORFs, motifs, restriction enzyme sites, and repeats, and preserves regional meaning when editing, exporting, and visualizing; it is not another DNA data object.
- **Normal API:** Pass the dictionary directly in `dnakit.DNA(..., features=[...])`; use `type`, `start`, `end` for simple intervals, and the remaining fields are optional.
- **Advanced Compatibility Objects:** `dnakit.DNAFeature(...)`, `dnakit.Interval(...)`; use `CompoundLocation`, `UnresolvedLocation` for compound or unresolved positions.
- **Input:** Common feature dictionary must fill in `type`, `start`, `end`; `location` can also be used instead of `start/end`.
- **Sample code:**

```python
import dnakit

dna = dnakit.DNA(
    "ACGT",
    id="seq-1",
    features=[
        {
            "type": "motif",
            "start": 1,
            "end": 3,
            "id": "m1",
            "strand": "forward",
            "label": "Example site",
        }
    ],
)
feature = dna.features[0]
print(feature.type, feature.location, feature.strand.value)
```

- **Example results:**

```text
motif Interval(start=1, end=3) forward
```

## 5) CORE-005 Gap object

- **Function:** Explicitly save gaps of known or unknown length in the sequence to prevent missing areas from being mistaken as continuous bases during coordinate calculation or sequence splicing.
- **API:** `dnakit.DNA(parts[required], ...)`, `dnakit.Gap(length[required], kind[optional], crossable[optional], evidence[optional], metadata[optional])`, `dnakit.GapKind`; `DNASequence.from_fragments(...)` is the advanced compatible entry.
- **Input:** The required length is a positive integer, or `None` represents an unknown length; optional `kind`, `crossable`, `evidence`, `metadata`.
- **Sample code:**

```python
import dnakit

gap = dnakit.Gap(500, kind="scaffold", crossable=False, evidence=("paired-ends",))
seq = dnakit.DNA(["AC", gap, "GT"])
print(seq.symbol_length)
print(seq.coordinate_span)
```

- **Example results:**

```text
4
504
```

## 6) CORE-006 sequence type declaration

- **Function:** Declares the characters allowed in the sequence, linear or circular morphology, and single-chain or double-chain type, for input verification and subsequent algorithms to select the correct processing rules.
- **API:** `dnakit.DNA(data[required], alphabet[optional], topology[optional], strandedness[optional])`; `DNASequence`, `DNAAlphabet`, `Topology` and `Strandedness` are reserved as advanced types.
- **Input:** Required DNA content; optional `alphabet="strict"|"iupac"`, `topology="linear"|"circular"`, `strandedness="single"|"double"`.
- **Sample code:**

```python
import dnakit

seq = dnakit.DNA(
    "ACGN",
    alphabet="iupac",
    topology="circular",
    strandedness="double",
)
print(seq.alphabet.value, seq.topology.value, seq.strandedness.value)
```

- **Example results:**

```text
iupac circular double
```

## 7) CORE-007 coordinate system

- **Function:** Convert the coordinates of different file formats into the 0-based half-open interval used by DNAKit to avoid one-bit deviation during interval interception and format conversion.
- **API:** `dnakit.core.ExternalInterval(start[required], end[required], system[required], strand[optional])`, `dnakit.Interval(start[required], end[required])`, `dnakit.core.CompoundLocation(parts[required])`, `dnakit.core.import_location(external[required], sequence_length[optional])`, `dnakit.core.export_location(location[required], target_system[required], sequence_length[optional])`, `dnakit.core.reverse_strand_location(location[required], sequence_length[required])`.
- **Input:** Required starting point, end point, source or target coordinate system; optional `strand`, `sequence_length`. The sequence length must be given when crossing the ring origin.
- **Sample code:**

```python
from dnakit.core import ExternalInterval, export_location, import_location

external = ExternalInterval(2, 8, system="1-based-closed", strand="forward")
internal = import_location(external)
(converted,) = export_location(internal, target_system="0-based-half-open")
print(internal)
print(converted.start, converted.end)
```

- **Example results:**

```text
Interval(start=1, end=8)
1 8
```
