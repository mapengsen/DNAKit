# 1. Core data object

Use `dnakit.DNA(...)` for one or more DNA sequences.

## 1) DNA sequence object {#1-core-001-dna-sequence-object}

- **Purpose:** Save DNA sequence values, alphabets, linear or circular states, single-stranded or double-stranded states, and Gap information as basic sequence objects for descriptors, searches, alignments, and other calculations.
- **Method:** Normalize the string and store sequence, topology and strandedness in one DNA object.
- **API:** `dnakit.DNA(data[required], alphabet[optional], topology[optional], strandedness[optional])`.
- **Input:** `DNA` can directly receive raw strings and automatically normalize them; advanced `DNASequence` only receives normalized sequences.
- **Output:** `DNA(...)` returns a `DNASequence` with sequence symbols and length.
- **Example code:**

```python
import dnakit

seq = dnakit.DNA("ACGT", topology="linear", strandedness="double")
print(seq.symbols)
print(seq.symbol_length)
print(seq.record_count)
```

- **Example output:**

```text
ACGT
4
1
```

## 2) DNA recording object {#2-core-002-dna-recording-object}

- **Purpose:** Attach ID, description, functional area and sample information to a DNA sequence so that the calculation results can be traced back to the specific record and its source.
- **Method:** Attach identifiers, descriptions, interval annotations and metadata to one sequence.
- **API:** `dnakit.DNA(data[required], id[optional], description[optional], features[optional], metadata[optional], letter_annotations[optional])`.
- **Input:** Normal entry only requires the sequence; the rest of the information is optional parameters in the same call. Automatically generated `sequence_1` when no ID is provided.
- **Output:** `DNA(...)` returns a `DNARecord` with an ID, sequence, features, and metadata.
- **Example code:**

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

- **Example output:**

```text
seq-1 ACGT human
```

## 3) DNA dataset object {#3-core-003-dna-dataset-object}

- **Purpose:** Use the same `DNA` object to manage one or more records in a fixed order, and supports subscript and slice selection, serving as a unified data input entrance for ordinary users.
- **Method:** Build records in input order; integer indexing and slicing both return DNA.
- **API:** `dnakit.DNA(data[required], name[optional], source[optional], version[optional], collection_metadata[optional])`; use `dnakit.read(..., mode="dna")` for file reading.
- **Input:** Use a string list for simple multiple sequences; use a dictionary list when each record requires different information.
- **Output:** `DNA(...)` returns a `DNASet` with ordered records and their IDs.
- **Example code:**

```python
import dnakit

dataset = dnakit.DNA(
    ["AC", "GT"],
    name="demo",
)
print(dataset.ids)
print(dataset[1].symbols)
```

- **Example output:**

```text
('sequence_1', 'sequence_2')
GT
```

## 4) characteristic object {#4-core-004-characteristic-object}

- **Purpose:** A feature is an additional annotation for a sequence in `DNA` that tags ORFs, motifs, restriction enzyme sites, and repeats, and preserves regional meaning when editing, exporting, and visualizing; it is not another DNA data object.
- **Method:** Convert feature dictionaries into typed annotations with 0-based half-open coordinates.
- **API:** Pass the dictionary directly in `dnakit.DNA(..., features=[...])`; use `type`, `start`, `end` for simple intervals, and the remaining fields are optional.
- **Input:** Common feature dictionary must fill in `type`, `start`, `end`; `location` can also be used instead of `start/end`.
- **Output:** Returns a `DNARecord` containing `DNAFeature` entries with type, interval, and strand.
- **Example code:**

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

- **Example output:**

```text
motif Interval(start=1, end=3) forward
```

## 5) Gap object {#5-core-005-gap-object}

- **Purpose:** Explicitly save gaps of known or unknown length in the sequence to prevent missing areas from being mistaken as continuous bases during coordinate calculation or sequence splicing.
- **Method:** Represent known or unknown lengths with explicit Gap objects without converting them into bases.
- **API:** `dnakit.DNA(parts[required], ...)`, `dnakit.Gap(length[required], kind[optional], crossable[optional], evidence[optional], metadata[optional])`, `dnakit.GapKind`; `DNASequence.from_fragments(...)` is the advanced compatible entry.
- **Input:** The required length is a positive integer, or `None` represents an unknown length; optional `kind`, `crossable`, `evidence`, `metadata`.
- **Output:** `Gap` stores a gap length; a `DNASequence` containing it exposes symbol length and coordinate span.
- **Example code:**

```python
import dnakit

gap = dnakit.Gap(500, kind="scaffold", crossable=False, evidence=("paired-ends",))
seq = dnakit.DNA(["AC", gap, "GT"])
print(seq.symbol_length)
print(seq.coordinate_span)
```

- **Example output:**

```text
4
504
```

## 6) sequence type declaration {#6-core-006-sequence-type-declaration}

- **Purpose:** Declares the characters allowed in the sequence, linear or circular morphology, and single-chain or double-chain type, for input verification and subsequent algorithms to select the correct processing rules.
- **Method:** Use the declared alphabet, topology and strandedness to validate the object state.
- **API:** `dnakit.DNA(data[required], alphabet[optional], topology[optional], strandedness[optional])`; `DNASequence`, `DNAAlphabet`, `Topology` and `Strandedness` are reserved as advanced types.
- **Input:** Required DNA content; optional `alphabet="strict"|"iupac"`, `topology="linear"|"circular"`, `strandedness="single"|"double"`.
- **Output:** Returns a `DNASequence` with `alphabet`, `topology`, and `strandedness` declarations.
- **Example code:**

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

- **Example output:**

```text
iupac circular double
```

## 7) coordinate system {#7-core-007-coordinate-system}

- **Purpose:** Convert the coordinates of different file formats into the 0-based half-open interval used by DNAKit to avoid one-bit deviation during interval interception and format conversion.
- **Method:** Convert endpoints from the source coordinate system into internal start-inclusive, end-exclusive intervals.
- **API:** `dnakit.core.ExternalInterval(start[required], end[required], system[required], strand[optional])`, `dnakit.Interval(start[required], end[required])`, `dnakit.core.CompoundLocation(parts[required])`, `dnakit.core.import_location(external[required], sequence_length[optional])`, `dnakit.core.export_location(location[required], target_system[required], sequence_length[optional])`, `dnakit.core.reverse_strand_location(location[required], sequence_length[required])`.
- **Input:** Required starting point, end point, source or target coordinate system; optional `strand`, `sequence_length`. The sequence length must be given when crossing the ring origin.
- **Output:** Returns an interval object with converted `start` and `end` coordinates.
- **Example code:**

```python
from dnakit.core import ExternalInterval, export_location, import_location

external = ExternalInterval(2, 8, system="1-based-closed", strand="forward")
internal = import_location(external)
(converted,) = export_location(internal, target_system="0-based-half-open")
print(internal)
print(converted.start, converted.end)
```

- **Example output:**

```text
Interval(start=1, end=8)
1 8
```

[FAQ: details](../../faq/tutorial_basics.md#core-objects)

[FAQ: details](../../faq/tutorial_basics.md#compatibility-01_core_objects)
