# pydna molecule ends and cloning history

These DNAKit interfaces call their documented upstream implementations. Optional dependencies are imported only when needed. Coordinate conventions, defaults and parameter choices are specified per function.

## `dnakit.molecular_lab.pydna_p34_functions.pydna_fill_in_ends`

**Purpose:** Fill in recessed 3' ends (5' overhangs) with a polymerase that has no exonuclease activity.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def pydna_fill_in_ends(dna: str | dict[str, Any], *, nucleotides: str | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Fill in recessed 3' ends (5' overhangs) with a polymerase that has no exonuclease activity.

This is the exo-minus Klenow fill-in of pydna ``Dseq.fill_in``. ``nucleotides`` limits the
dNTPs that are present (any combination of A, C, G, T); None means all four. ``dna`` is a string
or a ``{watson, crick, ovhg, circular}`` object. The result holds the new molecule.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.pydna_p34_functions.pydna_mung_bean_blunt`

**Purpose:** Blunt a molecule by removing single-stranded overhangs, as mung bean nuclease does.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def pydna_mung_bean_blunt(dna: str | dict[str, Any], *, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Blunt a molecule by removing single-stranded overhangs, as mung bean nuclease does.

Covers pydna ``Dseq.nibble_to_blunt`` and its alias ``Dseq.mung``: 5' and 3' overhangs are
both degraded, so the result is blunt whatever kind of overhang the input had.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.pydna_p34_functions.pydna_t4_polymerase_polish`

**Purpose:** Polish ends with T4 DNA polymerase: fill in 5' overhangs and chew back 3' overhangs.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def pydna_t4_polymerase_polish(dna: str | dict[str, Any], *, nucleotides: str | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Polish ends with T4 DNA polymerase: fill in 5' overhangs and chew back 3' overhangs.

This is pydna ``Dseq.T4``. ``nucleotides`` restricts the dNTPs present (A, C, G, T in any
combination); None means all four. With fewer dNTPs the exonuclease activity stops at the first
base that could be refilled, so partial sticky ends can be produced.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.pydna_p34_functions.pydna_exonuclease_five_prime`

**Purpose:** Resect one end of a molecule with a 5'->3' exonuclease, removing ``bases`` nucleotides.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def pydna_exonuclease_five_prime(dna: str | dict[str, Any], *, end: str='left', bases: int=1, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Resect one end of a molecule with a 5'->3' exonuclease, removing ``bases`` nucleotides.

``end`` is ``left`` (pydna ``Dseq.nibble_five_prime_left``) or ``right``
(``Dseq.nibble_five_prime_right``). The 5' strand of that end is shortened, which exposes a 3'
overhang.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.pydna_p34_functions.pydna_exonuclease_three_prime`

**Purpose:** Resect one end of a molecule with a 3'->5' exonuclease, removing ``bases`` nucleotides.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def pydna_exonuclease_three_prime(dna: str | dict[str, Any], *, end: str='left', bases: int=1, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Resect one end of a molecule with a 3'->5' exonuclease, removing ``bases`` nucleotides.

``end`` is ``left`` (pydna ``Dseq.nibble_three_prime_left``) or ``right``
(``Dseq.nibble_three_prime_right``). The 3' strand of that end is shortened, which exposes a 5'
overhang.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.pydna_p34_functions.pydna_terminal_transferase_tail`

**Purpose:** Add a tail of nucleotides to the 3' ends, as terminal transferase (TdT, Taq tailing) does.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def pydna_terminal_transferase_tail(dna: str | dict[str, Any], *, nucleotides: str='a', verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Add a tail of nucleotides to the 3' ends, as terminal transferase (TdT, Taq tailing) does.

This is pydna ``Dseq.terminal_transferase`` (and ``Dseqrecord.terminal_transferase``).
``nucleotides`` is the string of bases that is added to each 3' end; the default ``a`` gives the
single A overhang used in T/A cloning.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.pydna_p34_functions.pydna_user_enzyme_treatment`

**Purpose:** Simulate USER enzyme treatment: each uracil is excised and leaves a gap in its strand.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def pydna_user_enzyme_treatment(dna: str | dict[str, Any], *, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Simulate USER enzyme treatment: each uracil is excised and leaves a gap in its strand.

This is pydna ``Dseq.user``. Write the uracil as ``U`` in the strand; short stretches
between the gap and the end are released, which creates a 3' overhang.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.pydna_p34_functions.pydna_cast_to_double_strand`

**Purpose:** Convert single-stranded ends into double-stranded sequence by adding the complement.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def pydna_cast_to_double_strand(dna: str | dict[str, Any], *, end: str='both', verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Convert single-stranded ends into double-stranded sequence by adding the complement.

``end`` is ``both`` (pydna ``Dseq.cast_to_ds``), ``left`` (``Dseq.cast_to_ds_left``) or
``right`` (``Dseq.cast_to_ds_right``). The result is blunt on the converted ends.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.pydna_p34_functions.pydna_melt_short_strands`

**Purpose:** Release short single-stranded stretches from a duplex and report where they are.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def pydna_melt_short_strands(dna: str | dict[str, Any], *, length: int, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Release short single-stranded stretches from a duplex and report where they are.

Single-stranded pieces up to ``length`` nucleotides that are separated by nicks or gaps are shed
without changing the length of the remaining molecule. This covers pydna
``Dseq.get_ss_meltsites`` (the ``sites``), ``Dseq.melt`` (``melt_fragments``: the shed strands and the duplex pieces left when short double-
stranded stretches also melt) and ``Dseq.melt_ss_dna`` (``remaining`` and ``shed``).
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.pydna_p34_functions.pydna_shed_strand_segments`

**Purpose:** Separate chosen segments of the top (watson) and bottom (crick) strand from a molecule.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def pydna_shed_strand_segments(dna: str | dict[str, Any], *, watson_cutpairs: list[list[int]] | None=None, crick_cutpairs: list[list[int]] | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Separate chosen segments of the top (watson) and bottom (crick) strand from a molecule.

Each cut pair is ``[start, end]`` in the coordinates pydna ``Dseq.shed_ss_dna`` uses. The result
has the new molecule and the single-stranded pieces that were shed.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.pydna_p34_functions.pydna_restriction_cutters`

**Purpose:** List the restriction enzymes that cut a molecule a given number of times.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def pydna_restriction_cutters(dna: str | dict[str, Any], *, cutter_class: str='once', cuts: int=1, enzymes: list[str] | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
List the restriction enzymes that cut a molecule a given number of times.

``cutter_class`` selects the group: ``none`` (pydna ``Dseq.no_cutters``), ``once``
(``unique_cutters`` / ``once_cutters``), ``twice`` (``twice_cutters``), ``n`` (``n_cutters``
with ``cuts`` sites) or ``any`` (``cutters``, at least one site). ``enzymes`` is a list of
Biopython enzyme names; the default is pydna's batch of commercially available enzymes
(``CommOnly``). The result also holds the number of cuts of every selected enzyme.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.pydna_p34_functions.pydna_restriction_ligation_assembly`

**Purpose:** Digest fragments with restriction enzymes, ligate the compatible ends, and flag partial digests.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def pydna_restriction_ligation_assembly(fragments: list[Any], enzymes: list[str], *, allow_blunt: bool=True, circular_only: bool=False, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Digest fragments with restriction enzymes, ligate the compatible ends, and flag partial digests.

Runs pydna ``assembly2.restriction_ligation_assembly`` and, for every product, pydna
``assembly2.partially_digested_seqs_in_restriction_ligation_assembly``: the input fragments
that took part in a product without being cut between two sites of the enzymes. Fragments are
strings or ``{sequence, circular, name}`` objects; ``enzymes`` are Biopython names. With
``circular_only`` only circular products are returned.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.pydna_p34_functions.pydna_pcr_alignment_figure`

**Purpose:** Run an in-silico PCR and draw the primer-template alignment as text.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def pydna_pcr_alignment_figure(template: str | dict[str, Any], forward_primer: str, reverse_primer: str, *, limit: int=13, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Run an in-silico PCR and draw the primer-template alignment as text.

Uses pydna ``amplify.pcr`` and ``Amplicon.figure``. The figure shows the forward primer on the
top strand and the reverse primer on the bottom strand, with tails marked. ``limit`` is the
number of 3' bases that must anneal (pydna default 13). Exactly one product must form.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.pydna_p34_functions.pydna_assembly_figure`

**Purpose:** Assemble fragments by their shared end sequences and draw the junctions.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def pydna_assembly_figure(fragments: list[Any], *, limit: int=25, figure: str='compact', product_index: int=0, circular: bool=False, output_path: str | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Assemble fragments by their shared end sequences and draw the junctions.

Uses pydna ``assembly.Assembly`` and its contigs. ``figure`` is ``compact`` (``Contig.figure``,
fragment names with the overlap lengths), ``detailed`` (``Contig.detailed_figure``, the aligned
overlap sequences) or ``matplotlib`` (``Contig.figure_mpl``, written to ``output_path`` as
a PNG). ``limit`` is the minimum overlap length. ``circular`` chooses circular instead of linear
products; ``product_index`` selects one of several products (longest first).
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.pydna_p34_functions.pydna_build_cloning_strategy`

**Purpose:** Run a cloning step with pydna and export the product and its history as OpenCloning JSON.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def pydna_build_cloning_strategy(fragments: list[Any], *, operation: str='none', primers: list[Any] | None=None, enzymes: list[str] | None=None, limit: int=25, product_index: int=0, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Run a cloning step with pydna and export the product and its history as OpenCloning JSON.

``operation`` is ``none`` (export the fragments as they are), ``pcr`` (``fragments[0]`` is the
template; ``primers`` are the forward and reverse primer), ``gibson`` (Gibson assembly of the
fragments with minimum overlap ``limit``; for ``pcr`` the primer anneal limit is ``limit`` capped at 14) or ``restriction_ligation`` (``enzymes`` required). The
export is pydna ``CloningStrategy.from_dseqrecords`` followed by ``model_dump_json``; the JSON can
be read back with ``pydna_cloning_strategy_to_records``.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.pydna_p34_functions.pydna_cloning_strategy_to_records`

**Purpose:** Read an OpenCloning cloning-strategy JSON text and return the sequences it contains.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def pydna_cloning_strategy_to_records(strategy_json: str, *, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Read an OpenCloning cloning-strategy JSON text and return the sequences it contains.

Uses pydna ``opencloning_models.CloningStrategy.to_dseqrecords`` (which also reads the
``TextFileSequence`` entries through ``read_dseqrecord_from_text_file_sequence``). Each record
reports its name, sequence, topology and the cloning history as text.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.pydna_p34_functions.pydna_cloning_history`

**Purpose:** Show, validate or normalise the cloning history of the sequences in an OpenCloning strategy.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def pydna_cloning_history(strategy_json: str, *, action: str='display', verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Show, validate or normalise the cloning history of the sequences in an OpenCloning strategy.

``action`` is ``display`` (the history tree as text, pydna ``Dseqrecord.history``),
``validate`` (replay every recorded operation and compare the result with the stored sequence,
``Dseqrecord.validate_history``) or ``normalize`` (replay the operations and return the rebuilt
sequences, ``Dseqrecord.normalize_history``). Records without history are reported as such.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## Executed example

**Example code:**

```python
from dnakit.molecular_lab.pydna_p34_functions import pydna_fill_in_ends
r = pydna_fill_in_ends('ACGT', verbose=False)
print(r['result']['watson'], r['result']['crick'])
```

**Example output:**

```text
ACGT ACGT
```
