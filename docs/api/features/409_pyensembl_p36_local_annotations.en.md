# PyEnsembl local annotation queries

These DNAKit interfaces call their documented upstream implementations. Optional dependencies are imported only when needed. Coordinate conventions, defaults and parameter choices are specified per function.

## `dnakit.annotation.pyensembl_p36_functions.load_hgnc_gene_aliases`

**Purpose:** Load the HGNC gene-name alias table (symbols, aliases, previous symbols) mapped to Ensembl gene IDs.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def load_hgnc_gene_aliases(hgnc_tsv: str | Path | None=None, *, download: bool=False, source_url: str | None=None, cache_directory: str | Path | None=None, names: Sequence[str] | None=None, max_names: int=100000, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Load the HGNC gene-name alias table (symbols, aliases, previous symbols) mapped to Ensembl gene IDs.

With ``hgnc_tsv`` the file is read offline (``GeneNameAliases.from_hgnc``; plain or gzip complete-set
TSV). With ``download=True`` the current HGNC complete set (or ``source_url``) is fetched and cached
(``GeneNameAliases.download_hgnc``; network). ``names`` restricts the result to those symbols; otherwise
at most ``max_names`` entries are returned and ``total_names`` gives the full size.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.annotation.pyensembl_p36_functions.fetch_genomic_sequence`

**Purpose:** Read genomic DNA for a 1-based inclusive interval from a local genome FASTA (``Genome.sequence``).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def fetch_genomic_sequence(genome_fasta: str | Path, contig: str, start: int, end: int, *, strand: str='+', mask: str='upper', verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Read genomic DNA for a 1-based inclusive interval from a local genome FASTA (``Genome.sequence``).

``strand="-"`` returns the reverse complement. ``mask="upper"`` upper-cases the result, ``mask="raw"``
keeps soft masking (pyensembl's default is ``upper``). The contig may be named as in the FASTA.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.annotation.pyensembl_p36_functions.fetch_transcript_sequence`

**Purpose:** Look up the cDNA sequences of transcript IDs in a local transcript FASTA (``Genome.transcript_sequence``).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def fetch_transcript_sequence(transcript_fasta: str | Path, transcript_ids: Sequence[str], *, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Look up the cDNA sequences of transcript IDs in a local transcript FASTA (``Genome.transcript_sequence``).

IDs absent from the FASTA map to ``null`` in ``sequences`` and are listed in ``missing``. A versioned ID
such as ``ENST00000269305.8`` is accepted when the FASTA records that version.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.annotation.pyensembl_p36_functions.query_annotation_at_locus`

**Purpose:** List the transcripts, exons or protein IDs that overlap a locus (``Genome.transcripts_at_locus``,

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def query_annotation_at_locus(gtf: str | Path, contig: str, start: int, end: int | None=None, *, feature: str='transcript', strand: str | None=None, output: str='records', cache_dir: str | Path | None=None, max_results: int=1000, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
List the transcripts, exons or protein IDs that overlap a locus (``Genome.transcripts_at_locus``,
``exons_at_locus``, ``protein_ids_at_locus`` and their ``*_ids_at_locus`` / ``*_names_at_locus`` forms).

``feature`` is ``transcript``, ``exon`` or ``protein``; ``output`` is ``records`` (objects as
dictionaries), ``ids``, or ``names`` (transcripts only). ``end`` defaults to ``start``. For genes use the
existing ``query_genes_at_locus``.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.annotation.pyensembl_p36_functions.merge_gene_intervals`

**Purpose:** Union of all gene loci on a contig as sorted non-overlapping intervals (``Genome.merged_gene_intervals``).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def merge_gene_intervals(gtf: str | Path, contig: str, *, strand: str | None=None, cache_dir: str | Path | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Union of all gene loci on a contig as sorted non-overlapping intervals (``Genome.merged_gene_intervals``).

Adjacent intervals (``end + 1 == next start``) are merged; ``strand`` restricts to one strand.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.annotation.pyensembl_p36_functions.find_nearest_transcript`

**Purpose:** Find the transcript whose locus is nearest to a position or interval, even if none overlaps

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def find_nearest_transcript(gtf: str | Path, contig: str, position: int, *, end: int | None=None, strand: str | None=None, cache_dir: str | Path | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Find the transcript whose locus is nearest to a position or interval, even if none overlaps
(``Genome.nearest_transcript``).

``distance`` is 0 when the transcript overlaps the position (``position`` to ``end``).
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.annotation.pyensembl_p36_functions.get_annotation_by_id`

**Purpose:** Get the annotation of one gene, transcript or exon by ID (``Genome.gene_by_id``, ``transcript_by_id``,

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def get_annotation_by_id(gtf: str | Path, feature: str, identifier: str, *, cache_dir: str | Path | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Get the annotation of one gene, transcript or exon by ID (``Genome.gene_by_id``, ``transcript_by_id``,
``exon_by_id``; the name and locus look-ups ``gene_name_of_gene_id``, ``locus_of_*_id`` and
``transcript_name_of_transcript_id`` are included in the record).

``feature`` is ``gene``, ``transcript`` or ``exon``.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.annotation.pyensembl_p36_functions.list_annotation_contigs`

**Purpose:** List the contig names that carry genes in an annotation (``Genome.contigs``).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def list_annotation_contigs(gtf: str | Path, *, cache_dir: str | Path | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
List the contig names that carry genes in an annotation (``Genome.contigs``).
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.annotation.pyensembl_p36_functions.list_annotation_features`

**Purpose:** List genes, transcripts, exons or protein IDs, filtered by contig, strand or biotype

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def list_annotation_features(gtf: str | Path, feature: str='gene', *, contig: str | None=None, strand: str | None=None, biotype: str | None=None, output: str='ids', cache_dir: str | Path | None=None, max_results: int=1000, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
List genes, transcripts, exons or protein IDs, filtered by contig, strand or biotype
(``Genome.gene_ids`` / ``gene_names`` / ``genes``, ``transcript_ids`` / ``transcript_names`` /
``transcripts``, ``exon_ids`` / ``exons``, ``protein_ids``).

``feature`` is ``gene``, ``transcript``, ``exon`` or ``protein``; ``output`` is ``ids``, ``names``
(genes and transcripts) or ``records`` (genes, transcripts, exons). ``biotype`` applies to genes and
transcripts only (pyensembl has no exon or protein biotype).
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.annotation.pyensembl_p36_functions.related_annotation_features`

**Purpose:** Follow the relations between genes, transcripts, exons and proteins in an annotation.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def related_annotation_features(gtf: str | Path, from_kind: str, identifier: str, to_kind: str, *, cache_dir: str | Path | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Follow the relations between genes, transcripts, exons and proteins in an annotation.

``from_kind`` is ``gene_id``, ``gene_name``, ``transcript_id``, ``transcript_name``, ``exon_id`` or
``protein_id``; ``to_kind`` is ``gene``, ``transcript`` or ``exon``. Supported pairs: protein to gene
(``gene_by_protein_id``) and to transcript (``transcript_by_protein_id``); transcript to gene
(``gene_name_of_transcript_id`` / ``_name``) and to exons (``exon_ids_of_transcript_id`` / ``_name``);
transcript name to transcripts (``transcripts_by_name``); exon to gene (``gene_name_of_exon_id``) and to
transcripts (``transcript_ids_of_exon_id``); gene to transcripts (``transcript_ids_of_gene_id`` /
``_gene_name``) and to exons (``exon_ids_of_gene_id`` / ``_gene_name``).
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.annotation.pyensembl_p36_functions.map_position_to_transcript_offset`

**Purpose:** Convert a genomic position to the 0-based offset in the spliced transcript sequence

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def map_position_to_transcript_offset(gtf: str | Path, transcript_id: str, position: int, *, cache_dir: str | Path | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Convert a genomic position to the 0-based offset in the spliced transcript sequence
(``Transcript.spliced_offset``).

The offset counts exon bases only, from the transcript's 5' end (strand-aware). A position that lies
in an intron or outside every exon is reported with ``offset: null`` and the reason.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.annotation.pyensembl_p36_functions.genomic_locus_offset`

**Purpose:** Offset of a position (or a position range) from the stranded start of a locus

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def genomic_locus_offset(contig: str, start: int, end: int, strand: str, *, position: int | None=None, range_start: int | None=None, range_end: int | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Offset of a position (or a position range) from the stranded start of a locus
(``Locus.offset``, ``Locus.offset_range``).

On the plus strand the offset of position 13 in 10..20 is 3, on the minus strand 7. Give ``position``
or both ``range_start`` and ``range_end``; positions must lie inside the locus.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.annotation.pyensembl_p36_functions.compare_genomic_loci`

**Purpose:** Distance and overlap length of two loci (``Locus.distance_to_locus`` / ``distance_to_interval`` and

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def compare_genomic_loci(first: Mapping[str, Any], second: Mapping[str, Any], *, ignore_strand: bool=False, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Distance and overlap length of two loci (``Locus.distance_to_locus`` / ``distance_to_interval`` and
``Locus.overlap_length``).

Each locus is ``{"contig", "start", "end", "strand"}`` (1-based inclusive). ``distance`` is 0 for
overlapping loci and ``null`` for loci on different contigs (or strands, unless ``ignore_strand``);
``overlap_length`` is the number of shared bases (0 when they cannot overlap).
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## Executed example

**Example code:**

```python
from dnakit.annotation.pyensembl_p36_functions import genomic_locus_offset
r = genomic_locus_offset('chr1', 10, 20, '-', position=13, verbose=False)
print(r['offset'])
```

**Example output:**

```text
7
```
