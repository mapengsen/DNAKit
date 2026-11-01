# DnaChisel rules and DendroPy population models

These DNAKit interfaces call their documented upstream implementations. Optional dependencies are imported only when needed. Coordinate conventions, defaults and parameter choices are specified per function.

Validation scope: BLAST/Bowtie/primer3-dependent rules need their corresponding programs or packages and caller databases. No external program or remote resource is invoked by the examples.

## `dnakit.molecular_lab.dnachisel_p38_functions.evaluate_dnachisel_design_rule`

**Purpose:** Score a DNA sequence against one DnaChisel design rule, and optionally optimise it for the rule.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def evaluate_dnachisel_design_rule(sequence: str, rule: str, *, parameters: Mapping[str, Any] | None=None, circular: bool=False, role: str='objective', optimize: bool=False, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Score a DNA sequence against one DnaChisel design rule, and optionally optimise it for the rule.

``rule`` selects the specification: ``AllowPrimer`` (primer-friendly region: Tm window, no repeats, homologies
or heterodimers), ``AvoidBlastMatches`` (BLAST+ ``blast_db`` or ``sequences``), ``AvoidHeterodimerization``
(``other_primers_sequences``, ``tmax``), ``AvoidMatches`` (``match_length`` and a Bowtie ``bowtie_index``),
``EnforceChoice`` (``choices``), ``EnforcePatternOccurence`` (``pattern``, ``occurences``),
``EnforceRegionsCompatibility`` (``locations`` and ``condition`` ``different`` or ``min_hamming_distance`` with
``min_distance``), ``EnforceSequence`` (degenerate ``sequence`` at ``location``), ``EnforceTerminalGCContent``
(``window_size``, ``mini``, ``maxi``), ``SequenceLengthBounds`` (``min_length``, ``max_length``),
``UniquifyAllKmers`` (``k``), ``AvoidRareCodons`` (``min_frequency``, ``species``), ``MatchTargetCodonUsage``
(``species``) and ``HarmonizeRCA`` (``original_species``, ``species``). ``parameters`` holds the DnaChisel
keyword arguments of that rule; ``location`` is ``[start, end]`` or ``[start, end, strand]``.

``role`` is ``objective`` (the sequence is only scored; default) or ``constraint`` (DnaChisel first resolves
constraints that restrict nucleotides, which can change the sequence). With ``optimize=True`` the rule is used
as an objective (or constraint) and the sequence is optimised with DnaChisel's random-mutation optimiser; the
result holds the evaluation before and after and the new sequence. Optimisation is stochastic in DnaChisel.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.dnachisel_p38_functions.build_pssm_pattern_from_sequences`

**Purpose:** Build a PSSM pattern from same-length DNA sequences with pseudocounts, using

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def build_pssm_pattern_from_sequences(sequences: Sequence[str], *, name: str='unnamed', pseudocounts: Any='jaspar', threshold: float | None=None, relative_threshold: float | None=None, scan_sequence: str | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Build a PSSM pattern from same-length DNA sequences with pseudocounts, using
``dnachisel.MotifPssmPattern.from_sequences`` (and ``apply_pseudocounts``).

``pseudocounts`` is ``jaspar`` (automatic, Biopython), a mapping ``{"A": .., "C": .., "G": .., "T": ..}`` or
null for none. A match threshold is given absolute (``threshold``) or relative (``relative_threshold``, 0 to
1 between the minimum and maximum possible score). With ``scan_sequence`` the forward-strand matches are
returned. Returns the log-odds matrix by base.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.molecular_lab.dnachisel_p38_functions.read_pssm_patterns_from_motif_file`

**Purpose:** Read motifs from a JASPAR, MEME or TRANSFAC file and return them as PSSM patterns, using

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def read_pssm_patterns_from_motif_file(motifs_path: str, file_format: str, *, pseudocounts: Any='jaspar', threshold: float | None=None, relative_threshold: float | None=None, scan_sequence: str | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Read motifs from a JASPAR, MEME or TRANSFAC file and return them as PSSM patterns, using
``dnachisel.MotifPssmPattern.list_from_file`` (Bio.motifs parsers).

``file_format`` is ``jaspar``, ``meme`` (MEME XML, as read by Biopython) or ``TRANSFAC``. Pseudocounts and thresholds are as in
``build_pssm_pattern_from_sequences``; with ``scan_sequence`` the forward-strand matches of every motif are
listed.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.phylogeny.dendropy_p38_functions.read_multi_phylip_alignments`

**Purpose:** Read a multi-block PHYLIP file (several alignments in one file) with dendropy's ``multiphylip`` reader.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def read_multi_phylip_alignments(path: str, *, data_type: str='dna', strict: bool=False, interleaved: bool=False, limit: int=100, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Read a multi-block PHYLIP file (several alignments in one file) with dendropy's ``multiphylip`` reader.

Every block starts with its own ``ntax nchar`` header line. ``data_type`` is ``dna`` or ``rna``; ``strict``
selects strict (10-character names) rather than relaxed PHYLIP; ``interleaved`` reads interleaved blocks.
Returns each alignment's taxa, length and sequences (at most ``limit`` blocks).
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.phylogeny.dendropy_p38_functions.calculate_hky85_rate_matrix`

**Purpose:** HKY85 instantaneous nucleotide rate matrix Q and the branch-length correction factor, using

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def calculate_hky85_rate_matrix(*, base_frequencies: Sequence[float]=(0.25, 0.25, 0.25, 0.25), kappa: float=1.0, rate: float=1.0, branch_length: float | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
HKY85 instantaneous nucleotide rate matrix Q and the branch-length correction factor, using
``dendropy.model.discrete.Hky85.qmatrix`` and ``corrected_substitution_rate``.

``base_frequencies`` are ordered A, C, G, T; ``kappa`` is the transition/transversion rate ratio; ``rate`` is
the relative rate. ``corrected_substitution_rate`` is the factor that converts a branch length to expected
substitutions per site (``rate * 0.5 / f``). With ``branch_length`` the matrix of transition probabilities
over that length (``pij``) is returned too.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.phylogeny.dendropy_p38_functions.sample_stationary_nucleotide_sequence`

**Purpose:** Sample a root nucleotide sequence from the stationary base distribution, using

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def sample_stationary_nucleotide_sequence(length: int, *, base_frequencies: Sequence[float]=(0.25, 0.25, 0.25, 0.25), seed: int | None=None, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Sample a root nucleotide sequence from the stationary base distribution, using
``dendropy.model.discrete.NucleotideCharacterEvolutionModel.stationary_sample``.

``base_frequencies`` are ordered A, C, G, T; ``seed`` makes the draw reproducible.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.phylogeny.dendropy_p38_functions.simulate_fragmented_population_sequences`

**Purpose:** Simulate DNA sequences from fragmented (structured) populations, using

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def simulate_fragmented_population_sequences(species_name: str, *, div_time_gens: int, num_desc_pops: int=2, mutrate_per_site_per_generation: float=1e-07, desc_pop_size: int=10000, samples_per_pop: int=10, seq_len: int=2000, seed: int | None=None, output_path: str | None=None, output_format: str='fasta', verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Simulate DNA sequences from fragmented (structured) populations, using
``dendropy.simulate.popgensim.FragmentedPopulations.generate_sequences``.

A population tree with ``num_desc_pops`` descendant populations of ``desc_pop_size`` diverging
``div_time_gens`` generations ago is built, a gene tree is drawn with ``samples_per_pop`` genes per
population, and ``seq_len`` sites are evolved on the mutation tree (dendropy's built-in HKY85, equal base
frequencies). The sequences are returned and, with ``output_path``, written in ``output_format`` (fasta,
nexus, nexml or phylip). ``seed`` makes the run reproducible.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.phylogeny.dendropy_p38_functions.calculate_folded_site_frequency_spectrum`

**Purpose:** Folded (minor-allele) site frequency spectrum of a DNA alignment, using

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def calculate_folded_site_frequency_spectrum(sequences: Mapping[str, str], *, pad_to_unfolded_length: bool=False, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Folded (minor-allele) site frequency spectrum of a DNA alignment, using
``DiscreteCharacterMatrix.folded_site_frequency_spectrum``.

``sequences`` maps IDs to aligned DNA strings. Entry i counts the sites whose minor allele occurs in i
sequences; the vector has ceil(N/2) entries, or N+1 with ``pad_to_unfolded_length``, whose first entry then
counts monomorphic sites (as dendropy does).
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.phylogeny.dendropy_p38_functions.calculate_population_pair_statistics`

**Purpose:** Two-population pairwise summary statistics, using

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def calculate_population_pair_statistics(population1: Mapping[str, str], population2: Mapping[str, str], *, ignore_uncertain: bool=True, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Two-population pairwise summary statistics, using
``dendropy.calculate.popgenstat.PopulationPairSummaryStatistics`` (Hickerson et al. 2006, Wakeley 1996).

Each population maps IDs to aligned DNA strings of the same length. Returns average pairwise differences
(overall, between, within, net), segregating sites, Watterson's theta, Wakeley's psi and Tajima's D. Needs
at least two sequences per population.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.phylogeny.dendropy_p38_functions.parse_genbank_xml_records`

**Purpose:** Parse GenBank INSDSeq XML records offline, using dendropy's GenBank interop classes.

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def parse_genbank_xml_records(xml_path: str, *, molecule: str='dna', include_features: bool=True, fasta_labels: str='defline', limit: int=100, verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Parse GenBank INSDSeq XML records offline, using dendropy's GenBank interop classes.

``xml_path`` is an NCBI efetch ``rettype=gbc, retmode=xml`` file (INSDSet/INSDSeq XML, the format dendropy's
GenBank interop reads; plain ``gb`` GBSeq XML is not read by dendropy). For each record the function returns
the metadata (accession, version, organism, definition, taxonomy, dates, other sequence IDs), the sequence,
its FASTA text (``GenBankAccessionRecord.as_fasta``; ``fasta_labels`` is ``defline`` for the NCBI-style
``gi|..|gb|..|`` header or ``accession_organism`` for a composed label) and, with ``include_features``, the
feature table (key, location, qualifiers). ``molecule`` is ``dna`` or ``rna``.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.phylogeny.dendropy_p38_functions.genbank_xml_to_alignment`

**Purpose:** Build a nucleotide character matrix from GenBank INSDSeq XML records and write it, using

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def genbank_xml_to_alignment(xml_path: str, output_path: str, *, output_format: str='nexml', molecule: str='dna', label_components: Sequence[str]=('accession', 'organism'), label_separator: str='_', verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Build a nucleotide character matrix from GenBank INSDSeq XML records and write it, using
``GenBankResourceStore.generate_char_matrix``.

Taxon labels are composed from ``label_components`` (record attributes such as ``accession``, ``organism``,
``gi``) joined by ``label_separator``. The records need not be aligned: write ``fasta`` for unaligned sequences.
(dendropy 5.1's GenBank annotation helpers ``as_annotation`` and ``as_reference_annotation`` raise ``TypeError``
and are not used; metadata and feature tables are available from ``parse_genbank_xml_records``.)
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.phylogeny.dendropy_p38_functions.concatenate_alignments`

**Purpose:** Concatenate alignments into a partitioned supermatrix, using

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def concatenate_alignments(alignment_paths: Sequence[str], output_path: str, *, input_format: str='fasta', output_format: str='nexus', data_type: str='dna', verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Concatenate alignments into a partitioned supermatrix, using
``dendropy.datamodel.charmatrixmodel.CharacterMatrix.concatenate`` (and ``concatenate_from_paths``).

Every alignment must contain the same taxa (labels are matched by name). The parts are recorded as
character subsets named after the input files; the ``nexus`` output keeps them as ``charset`` blocks. The
result reports each part's column range (0-based, half-open).
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## `dnakit.phylogeny.dendropy_p38_functions.extract_alignment_columns`

**Purpose:** Extract alignment columns by index or range, using ``CharacterMatrix.export_character_indices`` (and

**Method:** Validate inputs through the explicit wrapper and call the corresponding backend.

**API:**

```text
def extract_alignment_columns(alignment_path: str, output_path: str, *, columns: Sequence[int] | None=None, ranges: Sequence[Sequence[int]] | None=None, input_format: str='fasta', output_format: str='fasta', data_type: str='dna', verbose: bool=True) -> dict[str, Any]:
```

**Input:** Required arguments, optional parameters and native defaults appear in the signature. Parameter and result details follow.

```text
Extract alignment columns by index or range, using ``CharacterMatrix.export_character_indices`` (and
``export_character_subset`` for named subsets).

``columns`` lists 0-based column indices; ``ranges`` lists ``[start, end)`` intervals (0-based, half-open);
both can be given and are merged in alignment order. When the alignment is NEXUS with ``charset`` blocks,
``input_format="nexus"`` keeps them and ``ranges`` can be replaced by naming nothing: the function
extracts only the requested columns.
```

**Output:** A JSON-serialisable dictionary; file-producing functions return paths. Missing optional dependencies raise BackendUnavailableError with an installation hint.

## Executed example

**Example code:**

```python
from dnakit.molecular_lab.dnachisel_p38_functions import evaluate_dnachisel_design_rule
from dnakit.phylogeny.dendropy_p38_functions import calculate_hky85_rate_matrix
r = evaluate_dnachisel_design_rule('ACGT', 'SequenceLengthBounds', parameters={'min_length': 2, 'max_length': 6}, verbose=False)
print(r['evaluation']['score'])
print(calculate_hky85_rate_matrix(verbose=False)['qmatrix'][0])
```

**Example output:**

```text
0.0
[-1.0, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333]
```
