# DNA Cauldron records and assembly graphs

These DNAKit interfaces call their documented upstream implementations. Optional dependencies are imported only when needed. Coordinate conventions, defaults and parameter choices are specified per function.

## `dnakit.molecular_lab.dnacauldron_p35_functions.build_dnacauldron_record`

**Purpose:** Build an annotated DNA record from a raw DNA string with DNA Cauldron.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def build_dnacauldron_record(sequence: str, *, record_id: str='record', topology: str='linear', features: list[dict[str, Any]] | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Build an annotated DNA record from a raw DNA string with DNA Cauldron.

``topology`` is ``linear`` or ``circular`` (DNA Cauldron's ``set_record_topology``). Each item of
``features`` is ``{"type", "start", "end", "strand", "label"}`` and is added with DNA Cauldron's
``annotate_record`` (type default ``misc_feature``, strand default 1, span default the whole record).
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.dnacauldron_p35_functions.annotate_dnacauldron_record`

**Purpose:** Add feature annotations to a DNA record and/or set its topology (circular or linear).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def annotate_dnacauldron_record(record: Mapping[str, Any], *, features: list[dict[str, Any]] | None=None, topology: str | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Add feature annotations to a DNA record and/or set its topology (circular or linear).

``record`` is a record dictionary (see the module description). Each item of ``features`` is
``{"type", "start", "end", "strand", "label"}`` and is added with DNA Cauldron's ``annotate_record``;
``topology`` (``linear`` or ``circular``) is applied with ``set_record_topology``. The record you pass is
not changed; the annotated record is returned.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.dnacauldron_p35_functions.write_dnacauldron_records`

**Purpose:** Write DNA records to GenBank or FASTA files, one file per record, with DNA Cauldron's ``write_record``.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def write_dnacauldron_records(records: list[dict[str, Any]], output_dir: str, *, output_format: str='genbank', verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Write DNA records to GenBank or FASTA files, one file per record, with DNA Cauldron's ``write_record``.

``output_format`` is ``genbank`` (files ``<id>.gb``, as ``BackboneChoice.write_final_records`` writes) or
``fasta`` (files ``<id>.fa``). Files are created in ``output_dir`` and are never overwritten. As in DNA
Cauldron, GenBank locus names are cut to 20 characters; the file name keeps the full id.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.dnacauldron_p35_functions.get_dnacauldron_repository_records`

**Purpose:** Retrieve sequence records by ID from a DNA Cauldron sequence repository built from local files.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def get_dnacauldron_repository_records(record_ids: list[str], *, files: list[str] | None=None, folder: str | None=None, collection: str='parts', use_file_names_as_ids: bool=True, topology: str='default_to_linear', verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Retrieve sequence records by ID from a DNA Cauldron sequence repository built from local files.

The repository is filled from ``files`` (FASTA, GenBank, SnapGene or zip) and/or a ``folder`` through
``SequenceRepository.import_records`` (``collection`` names the collection; ``topology`` is ``linear``,
``circular``, ``default_to_linear`` or ``default_to_circular``). ``record_ids`` lists the records to
return (``get_record`` for one ID, ``get_records`` for several); an ID that is missing raises
``ConfigurationError`` with DNA Cauldron's name suggestions.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.dnacauldron_p35_functions.digest_dnacauldron_sticky_ends`

**Purpose:** Digest DNA with one or several restriction enzymes into sticky-end fragments (DNA Cauldron).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def digest_dnacauldron_sticky_ends(sequence: str, enzymes: list[str], *, circular: bool=False, method: str='sequence', verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Digest DNA with one or several restriction enzymes into sticky-end fragments (DNA Cauldron).

``method`` is ``sequence`` (``StickyEndSeq.list_from_sequence_digestion``) or ``record``
(``StickyEndFragment.list_from_record_digestion``, the route used by assembly mixes). Each fragment
lists its middle ``sequence`` and the left and right protruding ends with their strand. A sequence with
no site returns one fragment without ends. ``enzymes`` are Biopython restriction enzyme names.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.dnacauldron_p35_functions.transform_dnacauldron_sticky_sequence`

**Purpose:** Reverse-complement, flatten or slice a sticky-end sequence.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def transform_dnacauldron_sticky_sequence(fragment: dict[str, Any], *, operation: str='reverse_complement', start: int | None=None, end: int | None=None, discard_sticky_ends: bool=False, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Reverse-complement, flatten or slice a sticky-end sequence.

``operation`` is ``reverse_complement`` (the ends are swapped and reverse-complemented), ``to_standard``
(the plain DNA string left end + middle + right end; with ``discard_sticky_ends`` only the middle) or
``slice`` (``start`` and ``end`` are zero-based coordinates on the middle sequence; the overhangs are
dropped, as ``StickyEndSeq.slice_seq`` does).
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.dnacauldron_p35_functions.check_dnacauldron_sticky_end_compatibility`

**Purpose:** Check whether sticky-end fragments clip together in the given order.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def check_dnacauldron_sticky_end_compatibility(fragments: list[dict[str, Any]], *, circular: bool=False, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Check whether sticky-end fragments clip together in the given order.

Fragment ``i`` clips with fragment ``i+1`` when its right end has the same bases as the next left end on
the opposite strand (``StickyEnd.will_clip_directly_with``; empty ends never clip). With
``circular=True`` the last fragment is also checked against the first.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.dnacauldron_p35_functions.sticky_fragment_to_dnacauldron_record`

**Purpose:** Convert a sticky-end fragment or one of its overhangs to an annotated record.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def sticky_fragment_to_dnacauldron_record(fragment: dict[str, Any], *, part: str='fragment', record_id: str='fragment', verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Convert a sticky-end fragment or one of its overhangs to an annotated record.

``part`` is ``fragment`` (left overhang + middle + right overhang, the overhangs annotated with their
strand, as ``StickyEndFragment.as_biopython_record``), ``left_end`` or ``right_end``
(``StickyEnd.as_biopython_record``, one record labelled ``(+) strand`` or ``(-) strand``).
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.dnacauldron_p35_functions.check_dnacauldron_homology`

**Purpose:** Validate a homology between two sequences by length, melting temperature and edit distance.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def check_dnacauldron_homology(sequence: str, other_sequence: str | None=None, *, min_size: int=15, max_size: int | None=80, min_tm: float=0, max_tm: float | None=None, max_distance: int=0, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Validate a homology between two sequences by length, melting temperature and edit distance.

This is ``HomologyChecker.check_homology`` with DNA Cauldron's defaults (15 to 80 bases, no Tm bounds,
edit distance 0). ``other_sequence`` is optional: without it only the length and Tm are checked.
The melting temperature is DNA Cauldron's rule 2 degrees per A/T and 4 per G/C. With
``max_distance`` above 0 small edits are accepted (Levenshtein distance).
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.dnacauldron_p35_functions.compute_dnacauldron_slots`

**Purpose:** Compute the slot of each fragment of a restriction-ligation mix (standardized overhang pairs).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def compute_dnacauldron_slots(parts: dict[str, str], *, enzymes: list[str] | None=None, topology: str='auto', filter_site_fragments: bool=True, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Compute the slot of each fragment of a restriction-ligation mix (standardized overhang pairs).

A fragment's slot is the alphabetically smaller of its (left, right) overhangs and the reverse
complement of the reversed pair. Parts that share a slot are interchangeable in an assembly. ``parts``
maps names to DNA strings (circular by default); ``enzymes`` default to ``["BsmBI"]``.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.dnacauldron_p35_functions.compute_dnacauldron_slots_graph`

**Purpose:** Compute the slots graph of a restriction-ligation mix as nodes and edges.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def compute_dnacauldron_slots_graph(parts: dict[str, str], *, enzymes: list[str] | None=None, topology: str='auto', filter_site_fragments: bool=True, with_overhangs: bool=True, directed: bool=True, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Compute the slots graph of a restriction-ligation mix as nodes and edges.

A slot node is written ``LEFT-RIGHT`` (its two overhangs). With ``with_overhangs`` the overhangs are
nodes too, and with ``directed`` the edges carry DNA Cauldron's arrow directions; without overhangs two
slots are linked when they share an overhang.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.dnacauldron_p35_functions.compute_dnacauldron_connections_graph`

**Purpose:** Compute which fragments of an assembly mix clip together and in which order.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def compute_dnacauldron_connections_graph(parts: dict[str, str], *, mix_type: str='restriction_ligation', enzymes: list[str] | None=None, topology: str='auto', filter_site_fragments: bool=True, bridging_oligos: dict[str, str] | None=None, min_size: int=15, max_size: int | None=80, min_tm: float=0, max_tm: float | None=None, max_distance: int=0, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Compute which fragments of an assembly mix clip together and in which order.

``mix_type`` is ``restriction_ligation`` (digestion by ``enzymes``, default ``["BsmBI"]``),
``homology`` (end homologies checked with ``min_size``/``max_size``/``min_tm``/``max_tm``/
``max_distance``, Gibson style) or ``ligase_cycling`` (needs ``bridging_oligos``). The result lists the
fragments (index, label, overhangs) and the directed edges ``[i, j]``: fragment ``i`` clips before ``j``.
Fragments rejected by the site filter are left out of the edges.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.dnacauldron_p35_functions.plot_dnacauldron_connections_graph`

**Purpose:** Plot the fragment connection graph of an assembly mix to a new PNG file.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def plot_dnacauldron_connections_graph(parts: dict[str, str], output_path: str, *, mix_type: str='restriction_ligation', enzymes: list[str] | None=None, topology: str='auto', filter_site_fragments: bool=True, bridging_oligos: dict[str, str] | None=None, min_size: int=15, max_size: int | None=80, min_tm: float=0, max_tm: float | None=None, max_distance: int=0, figure_size: float=12.0, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Plot the fragment connection graph of an assembly mix to a new PNG file.

The mix is built as in ``compute_dnacauldron_connections_graph``. DNA Cauldron draws the
de-duplicated graph (``plot_connections_graph``): one node per fragment, labelled with the overhangs and
the part name, one arrow per allowed junction. ``figure_size`` is the figure edge in inches.
The result gives the PNG path and the drawn nodes and edges.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.dnacauldron_p35_functions.plot_dnacauldron_slots_graph`

**Purpose:** Plot the slots graph of a restriction-ligation mix to a new PNG file.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def plot_dnacauldron_slots_graph(parts: dict[str, str], output_path: str, *, enzymes: list[str] | None=None, topology: str='auto', filter_site_fragments: bool=True, with_overhangs: bool=False, directed: bool=True, highlighted_parts: list[str] | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Plot the slots graph of a restriction-ligation mix to a new PNG file.

Slots are numbered nodes with a legend of the parts that fill each slot; ``highlighted_parts`` marks
parts in blue; ``with_overhangs`` adds the overhangs as nodes (``StickyEndAssemblyMix.plot_slots_graph``).
The result gives the PNG path and the slots graph data.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.dnacauldron_p35_functions.enumerate_dnacauldron_circular_assemblies`

**Purpose:** Enumerate or sample the circular constructs an assembly mix can form.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def enumerate_dnacauldron_circular_assemblies(parts: dict[str, str], *, mix_type: str='restriction_ligation', enzymes: list[str] | None=None, topology: str='auto', filter_site_fragments: bool=True, bridging_oligos: dict[str, str] | None=None, min_size: int=15, max_size: int | None=80, min_tm: float=0, max_tm: float | None=None, max_distance: int=0, output: str='records', sampling: str='exhaustive', max_results: int=20, annotate_homologies: bool=False, random_seed: int | None=None, staling_cutoff: int=100, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Enumerate or sample the circular constructs an assembly mix can form.

``output`` is ``fragments_sets`` (the ordered fragment labels of each construct,
``compute_circular_fragments_sets``) or ``records`` (the assembled circular sequences,
``compute_circular_assemblies``). ``sampling`` is ``exhaustive`` (networkx cycle order, at most
``max_results``) or ``random`` (``compute_random_circular_fragments_sets``; ``random_seed`` makes it
repeatable; ``staling_cutoff`` stops when only repeats are found, reported as ``exhausted``). Each set
is standardized, so rotations and reverse complements of one construct are reported once.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.dnacauldron_p35_functions.enumerate_dnacauldron_linear_assemblies`

**Purpose:** Enumerate the linear constructs of an assembly mix (shortest path between fragments).

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def enumerate_dnacauldron_linear_assemblies(parts: dict[str, str], *, mix_type: str='restriction_ligation', enzymes: list[str] | None=None, topology: str='auto', filter_site_fragments: bool=True, bridging_oligos: dict[str, str] | None=None, min_size: int=15, max_size: int | None=80, min_tm: float=0, max_tm: float | None=None, max_distance: int=0, min_parts: int=2, max_results: int=20, annotate_homologies: bool=False, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Enumerate the linear constructs of an assembly mix (shortest path between fragments).

This is ``compute_linear_assemblies``: for every pair of connected fragments the shortest path through
the connection graph is assembled, assemblies with fewer than ``min_parts`` fragments are dropped and
equivalent chains are reported once. For sticky-end constructs ``left_end`` and ``right_end`` give the
free overhangs and ``sequence`` includes them. DNA Cauldron calls this feature "a bit undertested".
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.dnacauldron_p35_functions.transform_dnacauldron_fragment_chain`

**Purpose:** Reverse-complement or standardize the circular fragment chains of an assembly mix.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def transform_dnacauldron_fragment_chain(parts: dict[str, str], *, operation: str='standardized', mix_type: str='restriction_ligation', enzymes: list[str] | None=None, topology: str='auto', filter_site_fragments: bool=True, bridging_oligos: dict[str, str] | None=None, min_size: int=15, max_size: int | None=80, min_tm: float=0, max_tm: float | None=None, max_distance: int=0, max_results: int=20, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Reverse-complement or standardize the circular fragment chains of an assembly mix.

The raw cycles of the connection graph are wrapped in DNA Cauldron ``FragmentChain`` objects. With
``operation="reverse_complement"`` each chain is replaced by the inverted list of reversed fragments;
with ``"standardized"`` it is given its canonical rotation and orientation (backbone first when one part
is marked as backbone, otherwise at most half reversed nucleotides, rotated to start at the largest
fragment). The result lists, per cycle, the fragment labels before and after.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## Executed example

**Example code:**

```python
from dnakit.molecular_lab.dnacauldron_p35_functions import build_dnacauldron_record
r = build_dnacauldron_record('ACGT', record_id='part', verbose=False)
print(r['record']['id'], r['record']['sequence'])
```

**Example output:**

```text
part ACGT
```
