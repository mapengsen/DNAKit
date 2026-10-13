# Sequence analysis

This category contains **52** documented candidate capabilities; 52 have a locally callable formal entry (0 partial due to binary-format dependencies; 0 bounded external backends; 0 dependency-adapter only). The documented functions are exported through the formal DNAKit namespace; functions with a compatible JSON call plan are also exposed through the Agent/MCP registry.

## 1) apply_reverse_complement_augmentation

- **Purpose:** Document a reusable apply reverse complement augmentation capability for Sequence analysis workflows. Candidate CAP-0583 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the apply reverse complement augmentation operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/sequence/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequence.apply_reverse_complement_augmentation(sequence[required], include_original[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `sequence`; optional arguments: `include_original`. Sequence arguments accept DNA, RNA, or protein text according to the operation.

- **Output:** The local `apply_reverse_complement_augmentation` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequence import apply_reverse_complement_augmentation

result = apply_reverse_complement_augmentation(
    sequence="ATGGCC"
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0583",
  "operation": "apply_reverse_complement_augmentation",
  "result_type": "list",
  "result_preview": "['ATGGCC', 'GGCCAT']",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 2) back_transcribe_rna_to_dna

- **Purpose:** Document a reusable back transcribe rna to dna capability for Sequence analysis workflows. Candidate CAP-0584 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the back transcribe rna to dna operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/sequence/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequence.back_transcribe_rna_to_dna(sequence[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `sequence`; optional arguments: no. Sequence arguments accept DNA, RNA, or protein text according to the operation.

- **Output:** The local `back_transcribe_rna_to_dna` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequence import back_transcribe_rna_to_dna

result = back_transcribe_rna_to_dna(
    sequence="AUGGCC"
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0584",
  "operation": "back_transcribe_rna_to_dna",
  "result_type": "str",
  "result_preview": "'ATGGCC'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 3) build_minimizer_index

- **Purpose:** Document a reusable build minimizer index capability for Sequence analysis workflows. Candidate CAP-0585 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the build minimizer index operation, and return structured results with provenance and diagnostics. The dependency-free implementation is available at `src/dnakit/sequence/__init__.py` and is exported through the formal DNAKit namespace. It operates on explicit DNA strings without a reference database, trained model, or external executable. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequence.build_minimizer_index(sequence[required], k[optional], window[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `sequence`; optional arguments: `k`, `window`. Sequence arguments accept DNA, RNA, or protein text according to the operation.

- **Output:** The local `build_minimizer_index` function returns its deterministic Python result from explicit DNA strings. It does not query a reference database, load a trained model, or invoke an external executable.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequence import build_minimizer_index

result = build_minimizer_index(
    sequence="ATGGCC",
    k=3,
    window=3
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0585",
  "operation": "build_minimizer_index",
  "result_type": "dict",
  "result_preview": "{'ATG': [0], 'GCC': [3]}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 4) calculate_sequence_metrics

- **Purpose:** Document a reusable calculate sequence metrics capability for Sequence analysis workflows. Candidate CAP-0586 comes from 27 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate sequence metrics operation, and return structured results with provenance and diagnostics. The optional `metric` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/sequence/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequence.calculate_sequence_metrics(sequence[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `sequence`; optional arguments: no. Sequence arguments accept DNA, RNA, or protein text according to the operation.

- **Output:** The local `calculate_sequence_metrics` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequence import calculate_sequence_metrics

result = calculate_sequence_metrics(
    sequence="ATGGCC"
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0586",
  "operation": "calculate_sequence_metrics",
  "result_type": "dict",
  "result_preview": "{'length': 6, 'gc_fraction': 0.6666666666666666, 'gc_percent': 66.66666666666666, 'n_count': 0, 'ambiguous_count': 0, 'entropy_bits': 1.9182958340544896, 'composition': {'A': 1, 'C': 2, 'G': 2, 'T': 1}}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 5) cluster_sequences_by_identity

- **Purpose:** Document a reusable cluster sequences by identity capability for Sequence analysis workflows. Candidate CAP-0587 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the cluster sequences by identity operation, and return structured results with provenance and diagnostics. The dependency-free implementation is available at `src/dnakit/sequence/__init__.py` and is exported through the formal DNAKit namespace. It operates on explicit DNA strings without a reference database, trained model, or external executable. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequence.cluster_sequences_by_identity(sequences[required], threshold[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `sequences`; optional arguments: `threshold`. Sequence arguments accept DNA, RNA, or protein text according to the operation.

- **Output:** The local `cluster_sequences_by_identity` function returns its deterministic Python result from explicit DNA strings. It does not query a reference database, load a trained model, or invoke an external executable.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequence import cluster_sequences_by_identity

result = cluster_sequences_by_identity(
    sequences=["ACGT", "ACGA"],
    threshold=0.5
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0587",
  "operation": "cluster_sequences_by_identity",
  "result_type": "list",
  "result_preview": "[[0, 1]]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 6) compare_motif_alleles

- **Purpose:** Document a reusable compare motif alleles capability for Sequence analysis workflows. Candidate CAP-0588 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the compare motif alleles operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/sequence/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequence.compare_motif_alleles(left[required], right[required], motif[required], start[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `left`, `right`, `motif`; optional arguments: `start`. Sequence arguments accept DNA, RNA, or protein text according to the operation.

- **Output:** The local `compare_motif_alleles` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequence import compare_motif_alleles

result = compare_motif_alleles(
    left="ACGT",
    right="AGGT",
    motif="GCC"
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0588",
  "operation": "compare_motif_alleles",
  "result_type": "dict",
  "result_preview": "{'start': 0, 'end': 3, 'motif': 'GCC', 'left': 'ACG', 'right': 'AGG', 'left_matches': False, 'right_matches': False, 'changed': True}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 7) compare_sequence_sketches

- **Purpose:** Document a reusable compare sequence sketches capability for Sequence analysis workflows. Candidate CAP-0589 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the compare sequence sketches operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/sequence/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequence.compare_sequence_sketches(left[required], right[required], k[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `left`, `right`; optional arguments: `k`. Sequence arguments accept DNA, RNA, or protein text according to the operation.

- **Output:** The local `compare_sequence_sketches` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequence import compare_sequence_sketches

result = compare_sequence_sketches(
    left="ACGT",
    right="AGGT",
    k=3
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0589",
  "operation": "compare_sequence_sketches",
  "result_type": "dict",
  "result_preview": "{'k': 3, 'left_size': 2, 'right_size': 2, 'intersection_size': 0, 'jaccard': 0.0}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 8) complement_dna

- **Purpose:** Document a reusable complement dna capability for Sequence analysis workflows. Candidate CAP-0590 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the complement dna operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/sequence/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `no_exact_public_match`. Existing-name evidence: `complement`.

- **API:** `dnakit.sequence.complement_dna(sequence[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `sequence`; optional arguments: no. Sequence arguments accept DNA, RNA, or protein text according to the operation.

- **Output:** The local `complement_dna` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequence import complement_dna

result = complement_dna(
    sequence="ATGGCC"
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0590",
  "operation": "complement_dna",
  "result_type": "str",
  "result_preview": "'TACCGG'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 9) count_canonical_kmers

- **Purpose:** Document a reusable count canonical kmers capability for Sequence analysis workflows. Candidate CAP-0591 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the count canonical kmers operation, and return structured results with provenance and diagnostics. The dependency-free implementation is available at `src/dnakit/sequence/__init__.py` and is exported through the formal DNAKit namespace. It operates on explicit DNA strings without a reference database, trained model, or external executable. Review basis: `no_exact_public_match`. Existing-name evidence: `canonical_kmer`.

- **API:** `dnakit.sequence.count_canonical_kmers(sequence[required], k[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `sequence`; optional arguments: `k`. Sequence arguments accept DNA, RNA, or protein text according to the operation.

- **Output:** The local `count_canonical_kmers` function returns its deterministic Python result from explicit DNA strings. It does not query a reference database, load a trained model, or invoke an external executable.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequence import count_canonical_kmers

result = count_canonical_kmers(
    sequence="ATGGCC",
    k=3
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0591",
  "operation": "count_canonical_kmers",
  "result_type": "dict",
  "result_preview": "{'ATG': 1, 'CCA': 1, 'GCC': 2}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 10) count_kmers

- **Purpose:** Document a reusable count kmers capability for Sequence analysis workflows. Candidate CAP-0592 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the count kmers operation, and return structured results with provenance and diagnostics. The dependency-free implementation is available at `src/dnakit/sequence/__init__.py` and is exported through the formal DNAKit namespace. It operates on explicit DNA strings without a reference database, trained model, or external executable. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequence.count_kmers(sequence[required], k[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `sequence`; optional arguments: `k`. Sequence arguments accept DNA, RNA, or protein text according to the operation.

- **Output:** The local `count_kmers` function returns its deterministic Python result from explicit DNA strings. It does not query a reference database, load a trained model, or invoke an external executable.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequence import count_kmers

result = count_kmers(
    sequence="ATGGCC",
    k=3
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0592",
  "operation": "count_kmers",
  "result_type": "dict",
  "result_preview": "{'ATG': 1, 'GCC': 1, 'GGC': 1, 'TGG': 1}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 11) count_kmers_batch

- **Purpose:** Document a reusable count kmers batch capability for Sequence analysis workflows. Candidate CAP-0593 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the count kmers batch operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/sequence/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequence.count_kmers_batch(sequences[required], k[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `sequences`; optional arguments: `k`. Sequence arguments accept DNA, RNA, or protein text according to the operation.

- **Output:** The local `count_kmers_batch` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequence import count_kmers_batch

result = count_kmers_batch(
    sequences=["ACGT", "ACGA"],
    k=3
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0593",
  "operation": "count_kmers_batch",
  "result_type": "list",
  "result_preview": "[{'ACG': 1, 'CGT': 1}, {'ACG': 1, 'CGA': 1}]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 12) deduplicate_sequences

- **Purpose:** Document a reusable deduplicate sequences capability for Sequence analysis workflows. Candidate CAP-0594 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the deduplicate sequences operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/sequence/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequence.deduplicate_sequences(sequences[required], canonical[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `sequences`; optional arguments: `canonical`. Sequence arguments accept DNA, RNA, or protein text according to the operation.

- **Output:** The local `deduplicate_sequences` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequence import deduplicate_sequences

result = deduplicate_sequences(
    sequences=["ACGT", "ACGA"]
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0594",
  "operation": "deduplicate_sequences",
  "result_type": "list",
  "result_preview": "['ACGT', 'ACGA']",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 13) digest_sequence_in_silico

- **Purpose:** Document a reusable digest sequence in silico capability for Sequence analysis workflows. Candidate CAP-0595 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the digest sequence in silico operation, and return structured results with provenance and diagnostics. The dependency-free implementation is available at `src/dnakit/sequence/__init__.py` and is exported through the formal DNAKit namespace. It operates on explicit DNA strings without a reference database, trained model, or external executable. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequence.digest_sequence_in_silico(sequence[required], enzymes[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `sequence`, `enzymes`; optional arguments: no. Sequence arguments accept DNA, RNA, or protein text according to the operation.

- **Output:** The local `digest_sequence_in_silico` function returns its deterministic Python result from explicit DNA strings. It does not query a reference database, load a trained model, or invoke an external executable.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequence import digest_sequence_in_silico

result = digest_sequence_in_silico(
    sequence="ATGGCC",
    enzymes={"EcoRI": "GAATTC"}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0595",
  "operation": "digest_sequence_in_silico",
  "result_type": "dict",
  "result_preview": "{'sites': {'EcoRI': []}, 'cut_positions': [0, 6], 'fragment_lengths': [6]}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 14) extract_cds_from_transcript

- **Purpose:** Document a reusable extract cds from transcript capability for Sequence analysis workflows. Candidate CAP-0596 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the extract cds from transcript operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/sequence/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequence.extract_cds_from_transcript(transcript[required], intervals[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `transcript`, `intervals`; optional arguments: no. Sequence arguments accept DNA, RNA, or protein text according to the operation.

- **Output:** The local `extract_cds_from_transcript` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequence import extract_cds_from_transcript

result = extract_cds_from_transcript(
    transcript="ATGGCC",
    intervals=[(0, 6)]
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0596",
  "operation": "extract_cds_from_transcript",
  "result_type": "str",
  "result_preview": "'ATGGCC'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 15) generate_kmer_table

- **Purpose:** Document a reusable generate kmer table capability for Sequence analysis workflows. Candidate CAP-0598 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the generate kmer table operation, and return structured results with provenance and diagnostics. The dependency-free implementation is available at `src/dnakit/sequence/__init__.py` and is exported through the formal DNAKit namespace. It operates on explicit DNA strings without a reference database, trained model, or external executable. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequence.generate_kmer_table(sequences[required], k[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `sequences`; optional arguments: `k`. Sequence arguments accept DNA, RNA, or protein text according to the operation.

- **Output:** The local `generate_kmer_table` function returns its deterministic Python result from explicit DNA strings. It does not query a reference database, load a trained model, or invoke an external executable.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequence import generate_kmer_table

result = generate_kmer_table(
    sequences=["ACGT", "ACGA"],
    k=3
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0598",
  "operation": "generate_kmer_table",
  "result_type": "list",
  "result_preview": "[{'sequence_index': 0, 'kmer': 'ACG', 'count': 1}, {'sequence_index': 0, 'kmer': 'CGT', 'count': 1}, {'sequence_index': 1, 'kmer': 'ACG', 'count': 1}, {'sequence_index': 1, 'kmer': 'CGA', 'count': 1}]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 16) generate_markov_sequence

- **Purpose:** Document a reusable generate markov sequence capability for Sequence analysis workflows. Candidate CAP-0599 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the generate markov sequence operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/sequence/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequence.generate_markov_sequence(sequence[required], length[required], order[optional], seed[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `sequence`, `length`; optional arguments: `order`, `seed`. Sequence arguments accept DNA, RNA, or protein text according to the operation.

- **Output:** The local `generate_markov_sequence` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequence import generate_markov_sequence

result = generate_markov_sequence(
    sequence="ATGGCC",
    length=12,
    seed=7
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0599",
  "operation": "generate_markov_sequence",
  "result_type": "str",
  "result_preview": "'ATGCCCCCCCCC'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 17) generate_random_dna

- **Purpose:** Document a reusable generate random dna capability for Sequence analysis workflows. Candidate CAP-0600 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the generate random dna operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/sequence/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequence.generate_random_dna(length[required], seed[optional], alphabet[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `length`; optional arguments: `seed`, `alphabet`. Sequence arguments accept DNA, RNA, or protein text according to the operation.

- **Output:** The local `generate_random_dna` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequence import generate_random_dna

result = generate_random_dna(
    length=12,
    seed=7
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0600",
  "operation": "generate_random_dna",
  "result_type": "str",
  "result_preview": "'GGGGGGGGGGGG'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 18) mask_low_complexity

- **Purpose:** Document a reusable mask low complexity capability for Sequence analysis workflows. Candidate CAP-0601 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the mask low complexity operation, and return structured results with provenance and diagnostics. The dependency-free implementation is available at `src/dnakit/sequence/__init__.py` and is exported through the formal DNAKit namespace. It operates on explicit DNA strings without a reference database, trained model, or external executable. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequence.mask_low_complexity(sequence[required], window[optional], entropy_threshold[optional], mask[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `sequence`; optional arguments: `window`, `entropy_threshold`, `mask`. Sequence arguments accept DNA, RNA, or protein text according to the operation.

- **Output:** The local `mask_low_complexity` function returns its deterministic Python result from explicit DNA strings. It does not query a reference database, load a trained model, or invoke an external executable.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequence import mask_low_complexity

result = mask_low_complexity(
    sequence="ATGGCC",
    window=3
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0601",
  "operation": "mask_low_complexity",
  "result_type": "str",
  "result_preview": "'ANNNNN'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 19) merge_overlapping_windows

- **Purpose:** Document a reusable merge overlapping windows capability for Sequence analysis workflows. Candidate CAP-0602 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the merge overlapping windows operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/sequence/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequence.merge_overlapping_windows(windows[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `windows`; optional arguments: no. Sequence arguments accept DNA, RNA, or protein text according to the operation.

- **Output:** The local `merge_overlapping_windows` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequence import merge_overlapping_windows

result = merge_overlapping_windows(
    windows=[(0, 3), (2, 5)]
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0602",
  "operation": "merge_overlapping_windows",
  "result_type": "list",
  "result_preview": "[[0, 5]]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 20) mutate_sequence_at_positions

- **Purpose:** Document a reusable mutate sequence at positions capability for Sequence analysis workflows. Candidate CAP-0603 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the mutate sequence at positions operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/sequence/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequence.mutate_sequence_at_positions(sequence[required], mutations[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `sequence`, `mutations`; optional arguments: no. Sequence arguments accept DNA, RNA, or protein text according to the operation.

- **Output:** The local `mutate_sequence_at_positions` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequence import mutate_sequence_at_positions

result = mutate_sequence_at_positions(
    sequence="ATGGCC",
    mutations={1: "T"}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0603",
  "operation": "mutate_sequence_at_positions",
  "result_type": "str",
  "result_preview": "'ATGGCC'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 21) normalize_dna_sequence

- **Purpose:** Document a reusable normalize dna sequence capability for Sequence analysis workflows. Candidate CAP-0604 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the normalize dna sequence operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/sequence/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequence.normalize_dna_sequence(sequence[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `sequence`; optional arguments: no. Sequence arguments accept DNA, RNA, or protein text according to the operation.

- **Output:** The local `normalize_dna_sequence` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequence import normalize_dna_sequence

result = normalize_dna_sequence(
    sequence="ATGGCC"
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0604",
  "operation": "normalize_dna_sequence",
  "result_type": "str",
  "result_preview": "'ATGGCC'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 22) normalize_protein_sequence

- **Purpose:** Document a reusable normalize protein sequence capability for Sequence analysis workflows. Candidate CAP-0605 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the normalize protein sequence operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/sequence/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequence.normalize_protein_sequence(sequence[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `sequence`; optional arguments: no. Sequence arguments accept DNA, RNA, or protein text according to the operation.

- **Output:** The local `normalize_protein_sequence` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequence import normalize_protein_sequence

result = normalize_protein_sequence(
    sequence="ATGGCC"
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0605",
  "operation": "normalize_protein_sequence",
  "result_type": "str",
  "result_preview": "'ATGGCC'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 23) normalize_rna_sequence

- **Purpose:** Document a reusable normalize rna sequence capability for Sequence analysis workflows. Candidate CAP-0606 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the normalize rna sequence operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/sequence/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequence.normalize_rna_sequence(sequence[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `sequence`; optional arguments: no. Sequence arguments accept DNA, RNA, or protein text according to the operation.

- **Output:** The local `normalize_rna_sequence` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequence import normalize_rna_sequence

result = normalize_rna_sequence(
    sequence="AUGGCC"
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0606",
  "operation": "normalize_rna_sequence",
  "result_type": "str",
  "result_preview": "'AUGGCC'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 24) one_hot_encode_dna

- **Purpose:** Document a reusable one hot encode dna capability for Sequence analysis workflows. Candidate CAP-0607 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the one hot encode dna operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/sequence/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `no_exact_public_match`. Existing-name evidence: `one_hot_encode`.

- **API:** `dnakit.sequence.one_hot_encode_dna(sequence[required], include_ambiguous[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `sequence`; optional arguments: `include_ambiguous`. Sequence arguments accept DNA, RNA, or protein text according to the operation.

- **Output:** The local `one_hot_encode_dna` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequence import one_hot_encode_dna

result = one_hot_encode_dna(
    sequence="ATGGCC"
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0607",
  "operation": "one_hot_encode_dna",
  "result_type": "list",
  "result_preview": "[[1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 1, 0, 0]]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 25) one_hot_encode_rna

- **Purpose:** Document a reusable one hot encode rna capability for Sequence analysis workflows. Candidate CAP-0608 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the one hot encode rna operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/sequence/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `no_exact_public_match`. Existing-name evidence: `one_hot_encode`.

- **API:** `dnakit.sequence.one_hot_encode_rna(sequence[required], include_ambiguous[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `sequence`; optional arguments: `include_ambiguous`. Sequence arguments accept DNA, RNA, or protein text according to the operation.

- **Output:** The local `one_hot_encode_rna` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequence import one_hot_encode_rna

result = one_hot_encode_rna(
    sequence="AUGGCC"
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0608",
  "operation": "one_hot_encode_rna",
  "result_type": "list",
  "result_preview": "[[1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 1, 0, 0]]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 26) pad_sequence

- **Purpose:** Document a reusable pad sequence capability for Sequence analysis workflows. Candidate CAP-0609 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the pad sequence operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/sequence/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `no_exact_public_match`. Existing-name evidence: `ena_sequence|plot_sequence|sequence`.

- **API:** `dnakit.sequence.pad_sequence(sequence[required], length[required], symbol[optional], side[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `sequence`, `length`; optional arguments: `symbol`, `side`. Sequence arguments accept DNA, RNA, or protein text according to the operation.

- **Output:** The local `pad_sequence` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequence import pad_sequence

result = pad_sequence(
    sequence="ATGGCC",
    length=12
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0609",
  "operation": "pad_sequence",
  "result_type": "str",
  "result_preview": "'ATGGCCNNNNNN'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 27) predict_coding_regions

- **Purpose:** Document a reusable predict coding regions capability for Sequence analysis workflows. Candidate CAP-0610 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the predict coding regions operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/sequence/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequence.predict_coding_regions(sequence[required], min_aa[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `sequence`; optional arguments: `min_aa`. Sequence arguments accept DNA, RNA, or protein text according to the operation.

- **Output:** The local `predict_coding_regions` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequence import predict_coding_regions

result = predict_coding_regions(
    sequence="ATGGCC"
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0610",
  "operation": "predict_coding_regions",
  "result_type": "list",
  "result_preview": "[]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 28) remove_sequence_gaps

- **Purpose:** Document a reusable remove sequence gaps capability for Sequence analysis workflows. Candidate CAP-0611 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the remove sequence gaps operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/sequence/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequence.remove_sequence_gaps(sequence[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `sequence`; optional arguments: no. Sequence arguments accept DNA, RNA, or protein text according to the operation.

- **Output:** The local `remove_sequence_gaps` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequence import remove_sequence_gaps

result = remove_sequence_gaps(
    sequence="ATGGCC"
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0611",
  "operation": "remove_sequence_gaps",
  "result_type": "str",
  "result_preview": "'ATGGCC'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 29) remove_sequence_whitespace

- **Purpose:** Document a reusable remove sequence whitespace capability for Sequence analysis workflows. Candidate CAP-0612 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the remove sequence whitespace operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/sequence/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequence.remove_sequence_whitespace(sequence[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `sequence`; optional arguments: no. Sequence arguments accept DNA, RNA, or protein text according to the operation.

- **Output:** The local `remove_sequence_whitespace` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequence import remove_sequence_whitespace

result = remove_sequence_whitespace(
    sequence="ATGGCC"
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0612",
  "operation": "remove_sequence_whitespace",
  "result_type": "str",
  "result_preview": "'ATGGCC'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 30) replace_ambiguous_bases

- **Purpose:** Document a reusable replace ambiguous bases capability for Sequence analysis workflows. Candidate CAP-0613 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the replace ambiguous bases operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/sequence/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequence.replace_ambiguous_bases(sequence[required], replacement[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `sequence`; optional arguments: `replacement`. Sequence arguments accept DNA, RNA, or protein text according to the operation.

- **Output:** The local `replace_ambiguous_bases` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequence import replace_ambiguous_bases

result = replace_ambiguous_bases(
    sequence="ATGGCC"
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0613",
  "operation": "replace_ambiguous_bases",
  "result_type": "str",
  "result_preview": "'ATGGCC'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 31) reverse_complement_dna

- **Purpose:** Document a reusable reverse complement dna capability for Sequence analysis workflows. Candidate CAP-0614 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the reverse complement dna operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/sequence/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `no_exact_public_match`. Existing-name evidence: `reverse_complement|reverse_complement_record|reverse_complement_search`.

- **API:** `dnakit.sequence.reverse_complement_dna(sequence[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `sequence`; optional arguments: no. Sequence arguments accept DNA, RNA, or protein text according to the operation.

- **Output:** The local `reverse_complement_dna` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequence import reverse_complement_dna

result = reverse_complement_dna(
    sequence="ATGGCC"
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0614",
  "operation": "reverse_complement_dna",
  "result_type": "str",
  "result_preview": "'GGCCAT'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 32) reverse_sequence

- **Purpose:** Document a reusable reverse sequence capability for Sequence analysis workflows. Candidate CAP-0615 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the reverse sequence operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/sequence/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequence.reverse_sequence(sequence[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `sequence`; optional arguments: no. Sequence arguments accept DNA, RNA, or protein text according to the operation.

- **Output:** The local `reverse_sequence` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequence import reverse_sequence

result = reverse_sequence(
    sequence="ATGGCC"
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0615",
  "operation": "reverse_sequence",
  "result_type": "str",
  "result_preview": "'CCGGTA'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 33) scan_pwm_batch

- **Purpose:** Document a reusable scan pwm batch capability for Sequence analysis workflows. Candidate CAP-0617 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the scan pwm batch operation, and return structured results with provenance and diagnostics. The dependency-free implementation is available at `src/dnakit/sequence/__init__.py` and is exported through the formal DNAKit namespace. It operates on explicit DNA strings without a reference database, trained model, or external executable. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequence.scan_pwm_batch(sequences[required], pwm[required], threshold[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `sequences`, `pwm`; optional arguments: `threshold`. Sequence arguments accept DNA, RNA, or protein text according to the operation.

- **Output:** The local `scan_pwm_batch` function returns its deterministic Python result from explicit DNA strings. It does not query a reference database, load a trained model, or invoke an external executable.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequence import scan_pwm_batch

result = scan_pwm_batch(
    sequences=["ACGT", "ACGA"],
    pwm=[{"A": 1, "C": 0, "G": 0, "T": 0}],
    threshold=0.5
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0617",
  "operation": "scan_pwm_batch",
  "result_type": "list",
  "result_preview": "[[{'start': 0, 'end': 1, 'sequence': 'A', 'score': 1.0}], [{'start': 0, 'end': 1, 'sequence': 'A', 'score': 1.0}, {'start': 3, 'end': 4, 'sequence': 'A', 'score': 1.0}]]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 34) score_motif_hit

- **Purpose:** Document a reusable score motif hit capability for Sequence analysis workflows. Candidate CAP-0618 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the score motif hit operation, and return structured results with provenance and diagnostics. The dependency-free implementation is available at `src/dnakit/sequence/__init__.py` and is exported through the formal DNAKit namespace. It operates on explicit DNA strings without a reference database, trained model, or external executable. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequence.score_motif_hit(sequence[required], motif[required], position[optional], match_score[optional], mismatch_score[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `sequence`, `motif`; optional arguments: `position`, `match_score`, `mismatch_score`. Sequence arguments accept DNA, RNA, or protein text according to the operation.

- **Output:** The local `score_motif_hit` function returns its deterministic Python result from explicit DNA strings. It does not query a reference database, load a trained model, or invoke an external executable.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequence import score_motif_hit

result = score_motif_hit(
    sequence="ATGGCC",
    motif="GCC"
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0618",
  "operation": "score_motif_hit",
  "result_type": "dict",
  "result_preview": "{'position': 0, 'motif': 'GCC', 'observed': 'ATG', 'score': 0.0, 'matches': False}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 35) select_representative_sequences

- **Purpose:** Document a reusable select representative sequences capability for Sequence analysis workflows. Candidate CAP-0619 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the select representative sequences operation, and return structured results with provenance and diagnostics. The dependency-free implementation is available at `src/dnakit/sequence/__init__.py` and is exported through the formal DNAKit namespace. It operates on explicit DNA strings without a reference database, trained model, or external executable. Review basis: `no_exact_public_match`. Existing-name evidence: `select_representatives`.

- **API:** `dnakit.sequence.select_representative_sequences(sequences[required], threshold[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `sequences`; optional arguments: `threshold`. Sequence arguments accept DNA, RNA, or protein text according to the operation.

- **Output:** The local `select_representative_sequences` function returns its deterministic Python result from explicit DNA strings. It does not query a reference database, load a trained model, or invoke an external executable.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequence import select_representative_sequences

result = select_representative_sequences(
    sequences=["ACGT", "ACGA"],
    threshold=0.5
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0619",
  "operation": "select_representative_sequences",
  "result_type": "list",
  "result_preview": "[0]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 36) simulate_indels

- **Purpose:** Document a reusable simulate indels capability for Sequence analysis workflows. Candidate CAP-0620 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the simulate indels operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/sequence/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequence.simulate_indels(sequence[required], count[optional], seed[optional], max_length[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `sequence`; optional arguments: `count`, `seed`, `max_length`. Sequence arguments accept DNA, RNA, or protein text according to the operation.

- **Output:** The local `simulate_indels` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequence import simulate_indels

result = simulate_indels(
    sequence="ATGGCC",
    seed=7
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0620",
  "operation": "simulate_indels",
  "result_type": "str",
  "result_preview": "'ATTGGCC'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 37) simulate_point_mutations

- **Purpose:** Document a reusable simulate point mutations capability for Sequence analysis workflows. Candidate CAP-0621 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the simulate point mutations operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/sequence/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequence.simulate_point_mutations(sequence[required], rate[optional], seed[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `sequence`; optional arguments: `rate`, `seed`. Sequence arguments accept DNA, RNA, or protein text according to the operation.

- **Output:** The local `simulate_point_mutations` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequence import simulate_point_mutations

result = simulate_point_mutations(
    sequence="ATGGCC",
    rate=0.5,
    seed=7
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0621",
  "operation": "simulate_point_mutations",
  "result_type": "str",
  "result_preview": "'CATTTA'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 38) simulate_recombination

- **Purpose:** Document a reusable simulate recombination capability for Sequence analysis workflows. Candidate CAP-0622 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the simulate recombination operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/sequence/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequence.simulate_recombination(left[required], right[required], crossover[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `left`, `right`; optional arguments: `crossover`. Sequence arguments accept DNA, RNA, or protein text according to the operation.

- **Output:** The local `simulate_recombination` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequence import simulate_recombination

result = simulate_recombination(
    left="ACGT",
    right="AGGT"
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0622",
  "operation": "simulate_recombination",
  "result_type": "str",
  "result_preview": "'ACGT'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 39) simulate_sequence_evolution

- **Purpose:** Document a reusable simulate sequence evolution capability for Sequence analysis workflows. Candidate CAP-0623 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the simulate sequence evolution operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/sequence/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequence.simulate_sequence_evolution(sequence[required], generations[optional], mutation_rate[optional], seed[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `sequence`; optional arguments: `generations`, `mutation_rate`, `seed`. Sequence arguments accept DNA, RNA, or protein text according to the operation.

- **Output:** The local `simulate_sequence_evolution` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequence import simulate_sequence_evolution

result = simulate_sequence_evolution(
    sequence="ATGGCC",
    seed=7
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0623",
  "operation": "simulate_sequence_evolution",
  "result_type": "list",
  "result_preview": "['ATGGCC', 'ATGGCC']",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 40) split_sequence_windows

- **Purpose:** Document a reusable split sequence windows capability for Sequence analysis workflows. Candidate CAP-0624 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the split sequence windows operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/sequence/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequence.split_sequence_windows(sequence[required], window[required], step[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `sequence`, `window`; optional arguments: `step`. Sequence arguments accept DNA, RNA, or protein text according to the operation.

- **Output:** The local `split_sequence_windows` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequence import split_sequence_windows

result = split_sequence_windows(
    sequence="ATGGCC",
    window=3,
    step=3
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0624",
  "operation": "split_sequence_windows",
  "result_type": "list",
  "result_preview": "[{'start': 0, 'end': 3, 'sequence': 'ATG'}, {'start': 3, 'end': 6, 'sequence': 'GCC'}]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 41) tokenize_dna_bpe

- **Purpose:** Document a reusable tokenize dna bpe capability for Sequence analysis workflows. Candidate CAP-0625 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the tokenize dna bpe operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/sequence/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequence.tokenize_dna_bpe(sequence[required], merges[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `sequence`; optional arguments: `merges`. Sequence arguments accept DNA, RNA, or protein text according to the operation.

- **Output:** The local `tokenize_dna_bpe` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequence import tokenize_dna_bpe

result = tokenize_dna_bpe(
    sequence="ATGGCC"
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0625",
  "operation": "tokenize_dna_bpe",
  "result_type": "list",
  "result_preview": "['A', 'T', 'G', 'G', 'C', 'C']",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 42) tokenize_dna_kmer

- **Purpose:** Document a reusable tokenize dna kmer capability for Sequence analysis workflows. Candidate CAP-0626 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the tokenize dna kmer operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/sequence/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequence.tokenize_dna_kmer(sequence[required], k[optional], stride[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `sequence`; optional arguments: `k`, `stride`. Sequence arguments accept DNA, RNA, or protein text according to the operation.

- **Output:** The local `tokenize_dna_kmer` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequence import tokenize_dna_kmer

result = tokenize_dna_kmer(
    sequence="ATGGCC",
    k=3
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0626",
  "operation": "tokenize_dna_kmer",
  "result_type": "list",
  "result_preview": "['ATG', 'TGG', 'GGC', 'GCC']",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 43) transcribe_dna_to_rna

- **Purpose:** Document a reusable transcribe dna to rna capability for Sequence analysis workflows. Candidate CAP-0627 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the transcribe dna to rna operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/sequence/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequence.transcribe_dna_to_rna(sequence[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `sequence`; optional arguments: no. Sequence arguments accept DNA, RNA, or protein text according to the operation.

- **Output:** The local `transcribe_dna_to_rna` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequence import transcribe_dna_to_rna

result = transcribe_dna_to_rna(
    sequence="AUGGCC"
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0627",
  "operation": "transcribe_dna_to_rna",
  "result_type": "str",
  "result_preview": "'AUGGCC'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 44) translate_codon

- **Purpose:** Document a reusable translate codon capability for Sequence analysis workflows. Candidate CAP-0628 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the translate codon operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/sequence/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequence.translate_codon(codon[required], unknown[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `codon`; optional arguments: `unknown`. Sequence arguments accept DNA, RNA, or protein text according to the operation.

- **Output:** The local `translate_codon` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequence import translate_codon

result = translate_codon(
    codon="ATG"
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0628",
  "operation": "translate_codon",
  "result_type": "str",
  "result_preview": "'M'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 45) translate_dna_to_protein

- **Purpose:** Document a reusable translate dna to protein capability for Sequence analysis workflows. Candidate CAP-0629 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the translate dna to protein operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/sequence/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequence.translate_dna_to_protein(sequence[required], frame[optional], stop_at_stop[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `sequence`; optional arguments: `frame`, `stop_at_stop`. Sequence arguments accept DNA, RNA, or protein text according to the operation.

- **Output:** The local `translate_dna_to_protein` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequence import translate_dna_to_protein

result = translate_dna_to_protein(
    sequence="ATGGCC",
    frame=0
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0629",
  "operation": "translate_dna_to_protein",
  "result_type": "str",
  "result_preview": "'MA'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 46) translate_orf

- **Purpose:** Document a reusable translate orf capability for Sequence analysis workflows. Candidate CAP-0630 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the translate orf operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/sequence/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `no_exact_public_match`. Existing-name evidence: `translate`.

- **API:** `dnakit.sequence.translate_orf(sequence[required], start_codon[optional], stop_codons[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `sequence`; optional arguments: `start_codon`, `stop_codons`. Sequence arguments accept DNA, RNA, or protein text according to the operation.

- **Output:** The local `translate_orf` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequence import translate_orf

result = translate_orf(
    sequence="ATGGCC"
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0630",
  "operation": "translate_orf",
  "result_type": "list",
  "result_preview": "[]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 47) translate_rna_to_protein

- **Purpose:** Document a reusable translate rna to protein capability for Sequence analysis workflows. Candidate CAP-0631 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the translate rna to protein operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/sequence/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequence.translate_rna_to_protein(sequence[required], frame[optional], stop_at_stop[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `sequence`; optional arguments: `frame`, `stop_at_stop`. Sequence arguments accept DNA, RNA, or protein text according to the operation.

- **Output:** The local `translate_rna_to_protein` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequence import translate_rna_to_protein

result = translate_rna_to_protein(
    sequence="AUGGCC",
    frame=0
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0631",
  "operation": "translate_rna_to_protein",
  "result_type": "str",
  "result_preview": "'MA'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 48) translate_six_frames

- **Purpose:** Document a reusable translate six frames capability for Sequence analysis workflows. Candidate CAP-0632 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the translate six frames operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/sequence/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequence.translate_six_frames(sequence[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `sequence`; optional arguments: no. Sequence arguments accept DNA, RNA, or protein text according to the operation.

- **Output:** The local `translate_six_frames` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequence import translate_six_frames

result = translate_six_frames(
    sequence="ATGGCC"
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0632",
  "operation": "translate_six_frames",
  "result_type": "dict",
  "result_preview": "{'+1': 'MA', '+2': 'W', '+3': 'G', '-1': 'GH', '-2': 'A', '-3': 'P'}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 49) trim_sequence

- **Purpose:** Document a reusable trim sequence capability for Sequence analysis workflows. Candidate CAP-0633 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the trim sequence operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/sequence/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequence.trim_sequence(sequence[required], left[optional], right[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `sequence`; optional arguments: `left`, `right`. Sequence arguments accept DNA, RNA, or protein text according to the operation.

- **Output:** The local `trim_sequence` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequence import trim_sequence

result = trim_sequence(
    sequence="ATGGCC"
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0633",
  "operation": "trim_sequence",
  "result_type": "str",
  "result_preview": "'ATGGCC'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 50) validate_dna_alphabet

- **Purpose:** Document a reusable validate dna alphabet capability for Sequence analysis workflows. Candidate CAP-0634 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the validate dna alphabet operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/sequence/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequence.validate_dna_alphabet(sequence[required], allow_ambiguous[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `sequence`; optional arguments: `allow_ambiguous`. Sequence arguments accept DNA, RNA, or protein text according to the operation.

- **Output:** The local `validate_dna_alphabet` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequence import validate_dna_alphabet

result = validate_dna_alphabet(
    sequence="ATGGCC"
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0634",
  "operation": "validate_dna_alphabet",
  "result_type": "dict",
  "result_preview": "{'valid': True, 'length': 6, 'invalid_symbols': [], 'alphabet': 'iupac'}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 51) validate_protein_alphabet

- **Purpose:** Document a reusable validate protein alphabet capability for Sequence analysis workflows. Candidate CAP-0635 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the validate protein alphabet operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/sequence/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequence.validate_protein_alphabet(sequence[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `sequence`; optional arguments: no. Sequence arguments accept DNA, RNA, or protein text according to the operation.

- **Output:** The local `validate_protein_alphabet` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequence import validate_protein_alphabet

result = validate_protein_alphabet(
    sequence="ATGGCC"
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0635",
  "operation": "validate_protein_alphabet",
  "result_type": "dict",
  "result_preview": "{'valid': True, 'length': 6, 'invalid_symbols': [], 'alphabet': 'protein'}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 52) validate_rna_alphabet

- **Purpose:** Document a reusable validate rna alphabet capability for Sequence analysis workflows. Candidate CAP-0636 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the validate rna alphabet operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/sequence/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequence.validate_rna_alphabet(sequence[required], allow_ambiguous[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `sequence`; optional arguments: `allow_ambiguous`. Sequence arguments accept DNA, RNA, or protein text according to the operation.

- **Output:** The local `validate_rna_alphabet` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequence import validate_rna_alphabet

result = validate_rna_alphabet(
    sequence="AUGGCC"
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0636",
  "operation": "validate_rna_alphabet",
  "result_type": "dict",
  "result_preview": "{'valid': True, 'length': 6, 'invalid_symbols': [], 'alphabet': 'iupac'}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```
