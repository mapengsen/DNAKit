# pydna primer screening and sequence annotation

These DNAKit interfaces call their documented upstream implementations. Optional dependencies are imported only when needed. Coordinate conventions, defaults and parameter choices are specified per function.

Validation scope: Recombinase site scanning needs pyahocorasick. No package was installed into the validation environment.

## `dnakit.molecular_lab.pydna_p34_functions.pydna_dna_mass_from_amount`

**Purpose:** Mass in grams of an amount of a double-stranded molecule (pydna ``Dseqrecord.m``).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def pydna_dna_mass_from_amount(dna: str | dict[str, Any], *, amount_mol: float=5e-14, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Mass in grams of an amount of a double-stranded molecule (pydna ``Dseqrecord.m``).

The mass is the molecular weight of the molecule (``Dseq.mw``, g/mol) times ``amount_mol``
(mol). pydna's default amount is 5e-14 mol (0.05 pmol).
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.pydna_p34_functions.pydna_find_peptide`

**Purpose:** Find a peptide in the three forward frames of a DNA sequence (pydna ``find_aminoacids``).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def pydna_find_peptide(dna: str, peptide: str, *, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Find a peptide in the three forward frames of a DNA sequence (pydna ``find_aminoacids``).

Returns the DNA slice (0-based start, end) of the first match, or None. ``find_aa`` is the
alias. The sequence is searched in its three forward frames (the reverse strand is not).
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.pydna_p34_functions.pydna_map_sanger_traces`

**Purpose:** Map ABI Sanger trace files onto a sequence (pydna ``Dseqrecord.map_trace_files``).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def pydna_map_sanger_traces(sequence: str | dict[str, Any], trace_glob: str, *, limit: int=25, target_start: int | None=None, target_end: int | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Map ABI Sanger trace files onto a sequence (pydna ``Dseqrecord.map_trace_files``).

``trace_glob`` is a file pattern such as ``/data/reads/*.ab1``. A trace maps when it shares a
stretch of at least ``limit`` bases with the sequence (either strand). With ``target_start`` and
``target_end`` only traces that cover that region count as matching. The result lists the matching
files with the mapped coordinates.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.pydna_p34_functions.pydna_sync_circular_origin`

**Purpose:** Rotate (and if needed reverse-complement) a circular sequence to start where a reference starts.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def pydna_sync_circular_origin(sequence: str, reference: str, *, limit: int=25, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Rotate (and if needed reverse-complement) a circular sequence to start where a reference starts.

This is pydna ``Dseqrecord.synced``. The two sequences must share a stretch of at least
``limit`` bases (``limit`` is lowered for short sequences, as pydna does).
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.pydna_p34_functions.pydna_orfs_to_features`

**Purpose:** Annotate the open reading frames of both strands as CDS features with their translation.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def pydna_orfs_to_features(sequence: str | dict[str, Any], *, min_size: int=300, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Annotate the open reading frames of both strands as CDS features with their translation.

This is pydna ``Dseqrecord.orfs_to_features``. ``min_size`` is the minimum ORF length in bases
(pydna default 300). Each feature has the start, end and strand, the translation and the SEGUID
checksums of the DNA and the protein.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.pydna_p34_functions.pydna_clean_genbank_text`

**Purpose:** Repair malformed GenBank text so that it can be parsed (pydna ``genbankfixer.gbtext_clean``).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def pydna_clean_genbank_text(genbank_text: str, *, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Repair malformed GenBank text so that it can be parsed (pydna ``genbankfixer.gbtext_clean``).

Returns the repaired GenBank text and the parsed JSON representation (name, size, topology,
sequence, features) that pydna builds on the way.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.pydna_p34_functions.pydna_primer_pairs_by_size`

**Purpose:** Find primer pairs from a primer list whose PCR product size lies in a range.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def pydna_primer_pairs_by_size(template: str | dict[str, Any], primers: list[Any], *, short: int=500, long: int=2000, limit: int=16, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Find primer pairs from a primer list whose PCR product size lies in a range.

This is pydna ``primer_screen.primer_pairs``. Product size includes the primers; only primers
that bind in one place each are used. ``limit`` is the length of the 3' end used to find the
binding sites. A circular template also yields products across the origin.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.pydna_p34_functions.pydna_flanking_primer_pairs`

**Purpose:** Find primer pairs that flank a target region of the template.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def pydna_flanking_primer_pairs(template: str | dict[str, Any], primers: list[Any], *, target_start: int, target_end: int, limit: int=16, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Find primer pairs that flank a target region of the template.

This is pydna ``primer_screen.flanking_primer_pairs``: the forward primer binds at or before
``target_start`` and the reverse primer at or after ``target_end`` (0-based template
coordinates).
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.pydna_p34_functions.pydna_diagnostic_primer_sets`

**Purpose:** Select primer pairs or triplets whose PCR products have different sizes on each template.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def pydna_diagnostic_primer_sets(templates: list[Any], primers: list[Any], *, primers_per_set: int=2, min_relative_size_difference: float=0.2, short: int=500, long: int=1500, limit: int=16, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Select primer pairs or triplets whose PCR products have different sizes on each template.

``primers_per_set`` is 2 (pydna ``primer_screen.diff_primer_pairs``, a primer pair that gives a
distinct product size from every template) or 3 (``primer_screen.diff_primer_triplets``, a shared
primer plus one primer specific to each template). Such sets verify clones or genetic
modifications by product size. Needs at least two templates. Two product sizes count as distinct
when they differ by at least ``min_relative_size_difference`` (default 0.2, pydna's rule of 20 %
of the larger product) of the larger one.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.pydna_p34_functions.pydna_find_recombinase_sites`

**Purpose:** Find the two recognition sites of a site-specific recombinase on both strands.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def pydna_find_recombinase_sites(sequence: str | dict[str, Any], site1: str, site2: str, *, site1_name: str='site1', site2_name: str='site2', verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Find the two recognition sites of a site-specific recombinase on both strands.

This is pydna ``recombinase.Recombinase.find`` and ``Recombinase.annotate``
(and the same methods of ``RecombinaseCollection``). The homology core of each site is written in
lowercase, as in ``ATGCCCTAAaaTT`` / ``AAaaTTTTTTTCCCT``, and the two cores must match. The result
lists the 0-based locations per site and the annotated GenBank features.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.pydna_p34_functions.pydna_terminal_codon_usage`

**Purpose:** Look up the usage value of the first and last codon of a coding sequence.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def pydna_terminal_codon_usage(cds: str, *, organism: str='sce', codon: str='both', verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Look up the usage value of the first and last codon of a coding sequence.

``codon`` is ``start`` (pydna ``Seq.startcodon``), ``stop`` (``Seq.stopcodon``) or ``both``. The
value is the entry of pydna's start-codon or stop-codon usage table for ``organism`` (only
``sce``, Saccharomyces cerevisiae, is available) and is None when the codon is not in the table.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.pydna_p34_functions.pydna_cds_expression_report`

**Purpose:** Summarise a coding sequence for expression: CAI, GC, start and stop codon, N-end, rare codons.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def pydna_cds_expression_report(cds: str, *, organism: str='sce', verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Summarise a coding sequence for expression: CAI, GC, start and stop codon, N-end, rare codons.

This is the report of pydna ``Seq.express``: length in codons, codon adaptation index (geometric
mean of the organism's codon weights, as ``pydna_codon_adaptation_index``), GC fraction, start and
stop codon values, the N-end rule value of the second amino acid, the count of each rare codon
and the fraction of rare codons. ``organism`` is ``sce``.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.pydna_p34_functions.pydna_is_open_reading_frame`

**Purpose:** Check whether a sequence is an open reading frame on its 5'->3' strand (pydna ``isorf``).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def pydna_is_open_reading_frame(sequence: str, *, table: int=1, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Check whether a sequence is an open reading frame on its 5'->3' strand (pydna ``isorf``).

True when it starts with a start codon, ends with a stop codon, is a whole number of codons and
has no internal stop codon, under the NCBI translation ``table`` (default 1, the standard code).
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.pydna_p34_functions.pydna_seguid_stamp`

**Purpose:** Compute the SEGUID checksum stamp of a record and write it into the record's COMMENT.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def pydna_seguid_stamp(sequence: str | dict[str, Any], *, timestamp: str | None=None, tool: str='pydna', comment: str='', verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Compute the SEGUID checksum stamp of a record and write it into the record's COMMENT.

This is pydna ``SeqRecord.stamp``. The stamp is ``cdseguid=...`` for a circular molecule and
``ldseguid=...`` for a linear one; the comment reads ``<tool> <stamp> <timestamp> <comment>``.
``timestamp`` defaults to the current time (ISO format); pass a fixed text for reproducible
output.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.pydna_p34_functions.pydna_iupac_site_to_regex`

**Purpose:** Convert a degenerate (IUPAC) DNA site into a regular expression (pydna ``compute_regex_site``).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def pydna_iupac_site_to_regex(site: str, *, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Convert a degenerate (IUPAC) DNA site into a regular expression (pydna ``compute_regex_site``).

Ambiguity codes become character classes (``R`` becomes ``[AG]``) and the pattern is case
insensitive. The result also reports whether the pattern is the same for the reverse
complement.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## Executed example

**Example code:**

```python
from dnakit.molecular_lab.pydna_p34_functions import pydna_is_open_reading_frame
r = pydna_is_open_reading_frame('ATGAAATAA', verbose=False)
print(r['is_orf'])
```

**Example output:**

```text
True
```
