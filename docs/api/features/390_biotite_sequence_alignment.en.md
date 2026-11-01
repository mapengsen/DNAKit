# Biotite sequence and nucleic-acid operations

Biotite >=2,<3. Examples use local fixture files from the source checkout; no network is used.

## align_banded_pairwise

**Purpose：** Pairwise alignment restricted to a diagonal band (``align_banded``).

**Method：** Pairwise alignment restricted to a diagonal band (``align_banded``).

``band`` is ``[lower, upper]`` diagonals (``j - i``). ``local=True`` gives a local
alignment. The default scores are biotite's: match 1, mismatch -1, linear gap -1.

**API：** `dnakit.alignment.biotite_p30_functions.align_banded_pairwise(sequence_1: 'str', sequence_2: 'str', *, band: 'Sequence[int]', seq_type: 'str' = 'dna', matrix: 'str | None' = None, match: 'int' = 1, mismatch: 'int' = -1, gap_penalty: 'int' = -1, gap_extension: 'int | None' = None, local: 'bool' = False, max_number: 'int' = 1000, verbose: bool = True) -> 'dict[str, Any]'`

**Input：** Use the explicitly typed parameters in the signature; optional settings retain the documented upstream defaults.

**Output：** JSON-compatible result with status, computed values or image path, and provenance.

**Example code:**

```python
from dnakit.alignment.biotite_p30_functions import align_banded_pairwise
result = align_banded_pairwise(sequence_1='ACGTACGTTAGC', sequence_2='TTACGTACGTAAGC', band=[-3, 3], local=True, verbose=False)
print(result["status"])
```

**Example output:**

```text
completed
```

## align_local_gapped_from_seed

**Purpose：** Gapped local alignment grown from a seed with the X-drop rule (``align_local_gapped``).

**Method：** Gapped local alignment grown from a seed with the X-drop rule (``align_local_gapped``).

**API：** `dnakit.alignment.biotite_p30_functions.align_local_gapped_from_seed(sequence_1: 'str', sequence_2: 'str', *, seed: 'Sequence[int]', threshold: 'int', seq_type: 'str' = 'dna', matrix: 'str | None' = None, match: 'int' = 1, mismatch: 'int' = -1, gap_penalty: 'int' = -1, gap_extension: 'int | None' = None, max_number: 'int' = 1, direction: 'str' = 'both', verbose: bool = True) -> 'dict[str, Any]'`

**Input：** Use the explicitly typed parameters in the signature; optional settings retain the documented upstream defaults.

**Output：** JSON-compatible result with status, computed values or image path, and provenance.

**Example code:**

```python
from dnakit.alignment.biotite_p30_functions import align_local_gapped_from_seed
result = align_local_gapped_from_seed(sequence_1='ACGTACGTTAGC', sequence_2='TTACGTACGTAAGC', seed=[4, 6], threshold=5, verbose=False)
print(result["status"])
```

**Example output:**

```text
completed
```

## align_multiple_progressive

**Purpose：** Progressive multiple sequence alignment with a guide tree (``align_multiple``).

**Method：** Progressive multiple sequence alignment with a guide tree (``align_multiple``).

Returns the gapped sequences (input order), the guide-tree leaf order for display, the
guide-tree edges and the pairwise distance matrix used to build the tree.

**API：** `dnakit.alignment.biotite_p30_functions.align_multiple_progressive(sequences: 'Sequence[str]', *, seq_type: 'str' = 'dna', matrix: 'str | None' = None, match: 'int' = 1, mismatch: 'int' = -1, gap_penalty: 'int' = -1, gap_extension: 'int | None' = None, terminal_penalty: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`

**Input：** Use the explicitly typed parameters in the signature; optional settings retain the documented upstream defaults.

**Output：** JSON-compatible result with status, computed values or image path, and provenance.

**Example code:**

```python
from dnakit.alignment.biotite_p30_functions import align_multiple_progressive
result = align_multiple_progressive(sequences=['ACGTACGTTAGCATCG', 'TTACGTACGTAAGC', 'ACGTACGTTAGCATCG'], verbose=False)
print(result["status"])
```

**Example output:**

```text
completed
```

## align_ungapped_pairwise

**Purpose：** Align two sequences position by position without gaps (``align_ungapped``).

**Method：** Align two sequences position by position without gaps (``align_ungapped``).

The sequences must have the same length; the score is the sum of substitution scores.

**API：** `dnakit.alignment.biotite_p30_functions.align_ungapped_pairwise(sequence_1: 'str', sequence_2: 'str', *, seq_type: 'str' = 'dna', matrix: 'str | None' = None, match: 'int' = 1, mismatch: 'int' = -1, verbose: bool = True) -> 'dict[str, Any]'`

**Input：** Use the explicitly typed parameters in the signature; optional settings retain the documented upstream defaults.

**Output：** JSON-compatible result with status, computed values or image path, and provenance.

**Example code:**

```python
from dnakit.alignment.biotite_p30_functions import align_ungapped_pairwise
result = align_ungapped_pairwise(sequence_1='ACGT', sequence_2='ACGA', verbose=False)
print(result["status"])
```

**Example output:**

```text
completed
```

## alignment_evalue_parameters

**Purpose：** Estimate the Gumbel parameters lambda and K of a scoring scheme by sampling.

**Method：** Estimate the Gumbel parameters lambda and K of a scoring scheme by sampling.

Random sequences are drawn with the given symbol ``frequencies`` (uniform by default) and
aligned locally (``EValueEstimator.from_samples``). ``matrix`` names a Biotite matrix;
otherwise the ``match``/``mismatch`` scores are used for DNA. ``seed`` makes it reproducible.

**API：** `dnakit.alignment.biotite_p30_functions.alignment_evalue_parameters(*, seq_type: 'str' = 'dna', matrix: 'str | None' = None, match: 'int' = 1, mismatch: 'int' = -1, gap_penalty: 'int' = -5, gap_extension: 'int | None' = None, frequencies: 'Sequence[float] | None' = None, sample_length: 'int' = 1000, sample_size: 'int' = 1000, seed: 'int' = 0, verbose: bool = True) -> 'dict[str, Any]'`

**Input：** Use the explicitly typed parameters in the signature; optional settings retain the documented upstream defaults.

**Output：** JSON-compatible result with status, computed values or image path, and provenance.

**Example code:**

```python
from dnakit.alignment.biotite_p30_functions import alignment_evalue_parameters
result = alignment_evalue_parameters(sample_length=100, sample_size=100, seed=0, verbose=False)
print(result["status"])
```

**Example output:**

```text
completed
```

## alignment_from_strings

**Purpose：** Build an alignment (trace of residue indices) from gapped strings.

**Method：** Build an alignment (trace of residue indices) from gapped strings.

Covers ``Alignment.from_strings`` and ``Alignment.trace_from_strings``. In the trace, ``-1``
marks a gap; the ungapped sequences are returned too.

**API：** `dnakit.alignment.biotite_p30_functions.alignment_from_strings(aligned: 'Sequence[str]', *, seq_type: 'str' = 'dna', gap_character: 'str' = '-', verbose: bool = True) -> 'dict[str, Any]'`

**Input：** Use the explicitly typed parameters in the signature; optional settings retain the documented upstream defaults.

**Output：** JSON-compatible result with status, computed values or image path, and provenance.

**Example code:**

```python
from dnakit.alignment.biotite_p30_functions import alignment_from_strings
result = alignment_from_strings(aligned=['ACGTACGT', 'AC-TACGT', 'ACGTTCGT'], verbose=False)
print(result["status"])
```

**Example output:**

```text
completed
```

## alignment_log_evalue

**Purpose：** Natural logarithm of the E-value of an alignment score (``EValueEstimator.log_evalue``).

**Method：** Natural logarithm of the E-value of an alignment score (``EValueEstimator.log_evalue``).

Give ``lam`` and ``k`` (from ``alignment_evalue_parameters``) or leave both None to estimate
them by sampling for the given scoring scheme. Also returns the E-value itself.

**API：** `dnakit.alignment.biotite_p30_functions.alignment_log_evalue(score: 'int', seq1_length: 'int', seq2_length: 'int', *, lam: 'float | None' = None, k: 'float | None' = None, seq_type: 'str' = 'dna', matrix: 'str | None' = None, match: 'int' = 1, mismatch: 'int' = -1, gap_penalty: 'int' = -5, gap_extension: 'int | None' = None, seed: 'int' = 0, verbose: bool = True) -> 'dict[str, Any]'`

**Input：** Use the explicitly typed parameters in the signature; optional settings retain the documented upstream defaults.

**Output：** JSON-compatible result with status, computed values or image path, and provenance.

**Example code:**

```python
from dnakit.alignment.biotite_p30_functions import alignment_log_evalue
result = alignment_log_evalue(score=30, seq1_length=100, seq2_length=100, lam=1.0, k=0.1, verbose=False)
print(result["status"])
```

**Example output:**

```text
completed
```

## alignment_pairwise_identity

**Purpose：** Matrix of pairwise sequence identities of an alignment (``get_pairwise_sequence_identity``).

**Method：** Matrix of pairwise sequence identities of an alignment (``get_pairwise_sequence_identity``).

``mode``: ``all`` divides by all columns, ``not_terminal`` by the columns without terminal
gaps of the pair (biotite default), ``shortest`` by the length of the shorter sequence.

**API：** `dnakit.alignment.biotite_p30_functions.alignment_pairwise_identity(aligned: 'Sequence[str]', *, seq_type: 'str' = 'dna', mode: 'str' = 'not_terminal', verbose: bool = True) -> 'dict[str, Any]'`

**Input：** Use the explicitly typed parameters in the signature; optional settings retain the documented upstream defaults.

**Output：** JSON-compatible result with status, computed values or image path, and provenance.

**Example code:**

```python
from dnakit.alignment.biotite_p30_functions import alignment_pairwise_identity
result = alignment_pairwise_identity(aligned=['ACGTACGT', 'AC-TACGT', 'ACGTTCGT'], mode='all', verbose=False)
print(result["status"])
```

**Example output:**

```text
completed
```

## alignment_remove_gaps

**Purpose：** Drop every alignment column that contains a gap in any sequence (``remove_gaps``).

**Method：** Drop every alignment column that contains a gap in any sequence (``remove_gaps``).

**API：** `dnakit.alignment.biotite_p30_functions.alignment_remove_gaps(aligned: 'Sequence[str]', *, seq_type: 'str' = 'dna', verbose: bool = True) -> 'dict[str, Any]'`

**Input：** Use the explicitly typed parameters in the signature; optional settings retain the documented upstream defaults.

**Output：** JSON-compatible result with status, computed values or image path, and provenance.

**Example code:**

```python
from dnakit.alignment.biotite_p30_functions import alignment_remove_gaps
result = alignment_remove_gaps(aligned=['ACGTACGT', 'AC-TACGT', 'ACGTTCGT'], verbose=False)
print(result["status"])
```

**Example output:**

```text
completed
```

## alignment_remove_terminal_gaps

**Purpose：** Trim the leading and trailing columns that contain a gap (``remove_terminal_gaps``).

**Method：** Trim the leading and trailing columns that contain a gap (``remove_terminal_gaps``).

Also reports the number of trimmed columns at each end (``find_terminal_gaps``).

**API：** `dnakit.alignment.biotite_p30_functions.alignment_remove_terminal_gaps(aligned: 'Sequence[str]', *, seq_type: 'str' = 'dna', verbose: bool = True) -> 'dict[str, Any]'`

**Input：** Use the explicitly typed parameters in the signature; optional settings retain the documented upstream defaults.

**Output：** JSON-compatible result with status, computed values or image path, and provenance.

**Example code:**

```python
from dnakit.alignment.biotite_p30_functions import alignment_remove_terminal_gaps
result = alignment_remove_terminal_gaps(aligned=['--ACGT', 'TTACG-'], verbose=False)
print(result["status"])
```

**Example output:**

```text
completed
```

## alignment_score

**Purpose：** Score an existing pairwise alignment of gapped strings (``biotite...align.score``).

**Method：** Score an existing pairwise alignment of gapped strings (``biotite...align.score``).

**API：** `dnakit.alignment.biotite_p30_functions.alignment_score(aligned: 'Sequence[str]', *, seq_type: 'str' = 'dna', matrix: 'str | None' = None, match: 'int' = 1, mismatch: 'int' = -1, gap_penalty: 'int' = -1, gap_extension: 'int | None' = None, terminal_penalty: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`

**Input：** Use the explicitly typed parameters in the signature; optional settings retain the documented upstream defaults.

**Output：** JSON-compatible result with status, computed values or image path, and provenance.

**Example code:**

```python
from dnakit.alignment.biotite_p30_functions import alignment_score
result = alignment_score(aligned=['ACGTACGT', 'AC-TACGT'], verbose=False)
print(result["status"])
```

**Example output:**

```text
completed
```

## alignment_sequence_identity

**Purpose：** Overall sequence identity of an alignment, a single number (``get_sequence_identity``).

**Method：** Overall sequence identity of an alignment, a single number (``get_sequence_identity``).

``mode`` is as in ``alignment_pairwise_identity``.

**API：** `dnakit.alignment.biotite_p30_functions.alignment_sequence_identity(aligned: 'Sequence[str]', *, seq_type: 'str' = 'dna', mode: 'str' = 'not_terminal', verbose: bool = True) -> 'dict[str, Any]'`

**Input：** Use the explicitly typed parameters in the signature; optional settings retain the documented upstream defaults.

**Output：** JSON-compatible result with status, computed values or image path, and provenance.

**Example code:**

```python
from dnakit.alignment.biotite_p30_functions import alignment_sequence_identity
result = alignment_sequence_identity(aligned=['ACGTACGT', 'AC-TACGT', 'ACGTTCGT'], verbose=False)
print(result["status"])
```

**Example output:**

```text
completed
```

## codon_table_is_start_codon

**Purpose：** Tell whether each codon is a start codon of an NCBI genetic code table.

**Method：** Tell whether each codon is a start codon of an NCBI genetic code table.

Wraps ``CodonTable.is_start_codon``. ``codons`` are DNA triplets (A, C, G, T); ``table`` is
the NCBI translation table number (1 standard, 11 bacterial, ...).

**API：** `dnakit.alignment.biotite_p30_functions.codon_table_is_start_codon(codons: 'Sequence[str]', *, table: 'int' = 1, verbose: bool = True) -> 'dict[str, Any]'`

**Input：** Use the explicitly typed parameters in the signature; optional settings retain the documented upstream defaults.

**Output：** JSON-compatible result with status, computed values or image path, and provenance.

**Example code:**

```python
from dnakit.alignment.biotite_p30_functions import codon_table_is_start_codon
result = codon_table_is_start_codon(codons=['ATG', 'GTG', 'AAA'], table=11, verbose=False)
print(result["status"])
```

**Example output:**

```text
completed
```

## kmer_index

**Purpose：** Index the positions of every k-mer of a set of sequences.

**Method：** Index the positions of every k-mer of a set of sequences.

``bucketed=False`` uses ``KmerTable.from_sequences`` (exact k-mer to position lists);
``bucketed=True`` uses ``BucketKmerTable.from_sequences`` (hash buckets, ``n_buckets``).
``spacing`` is a spaced-seed pattern such as ``"1101"``. Each position is
``[sequence_index, offset]``. At most ``max_kmers`` distinct k-mers are listed.

**API：** `dnakit.alignment.biotite_p30_functions.kmer_index(sequences: 'Sequence[str]', *, k: 'int' = 4, seq_type: 'str' = 'dna', spacing: 'str | None' = None, bucketed: 'bool' = False, n_buckets: 'int | None' = None, max_kmers: 'int' = 1000, verbose: bool = True) -> 'dict[str, Any]'`

**Input：** Use the explicitly typed parameters in the signature; optional settings retain the documented upstream defaults.

**Output：** JSON-compatible result with status, computed values or image path, and provenance.

**Example code:**

```python
from dnakit.alignment.biotite_p30_functions import kmer_index
result = kmer_index(sequences=['ACGTACGTTAGC', 'TTACGTA'], k=4, verbose=False)
print(result["status"])
```

**Example output:**

```text
completed
```

## kmer_match

**Purpose：** Find k-mer matches (alignment seeds) between a query and reference sequences.

**Method：** Find k-mer matches (alignment seeds) between a query and reference sequences.

The references are indexed once (``KmerTable`` or, with ``bucketed=True``,
``BucketKmerTable``) and the query is matched against them. With ``similarity_threshold``
similar k-mers whose substitution score reaches the threshold also match
(``ScoreThresholdRule``; ``matrix`` names the substitution matrix, default NUC for DNA and
BLOSUM62 for protein). Each match is ``[query_position, reference_index,
reference_position]``.

**API：** `dnakit.alignment.biotite_p30_functions.kmer_match(reference_sequences: 'Sequence[str]', query: 'str', *, k: 'int' = 4, seq_type: 'str' = 'dna', spacing: 'str | None' = None, similarity_threshold: 'int | None' = None, matrix: 'str | None' = None, bucketed: 'bool' = False, max_matches: 'int' = 10000, verbose: bool = True) -> 'dict[str, Any]'`

**Input：** Use the explicitly typed parameters in the signature; optional settings retain the documented upstream defaults.

**Output：** JSON-compatible result with status, computed values or image path, and provenance.

**Example code:**

```python
from dnakit.alignment.biotite_p30_functions import kmer_match
result = kmer_match(reference_sequences=['ACGTACGTTAGC', 'TTACGTA'], query='CGTACG', k=4, verbose=False)
print(result["status"])
```

**Example output:**

```text
completed
```

## mincode_selection

**Purpose：** Select a subset of the k-mers of a sequence with the mincode method.

**Method：** Select a subset of the k-mers of a sequence with the mincode method.

Keeps the k-mers whose code is below a threshold set by ``compression`` (the desired factor
by which the k-mer count is reduced; ``MincodeSelector``).

**API：** `dnakit.alignment.biotite_p30_functions.mincode_selection(sequence: 'str', *, k: 'int' = 5, compression: 'float' = 2.0, seq_type: 'str' = 'dna', verbose: bool = True) -> 'dict[str, Any]'`

**Input：** Use the explicitly typed parameters in the signature; optional settings retain the documented upstream defaults.

**Output：** JSON-compatible result with status, computed values or image path, and provenance.

**Example code:**

```python
from dnakit.alignment.biotite_p30_functions import mincode_selection
result = mincode_selection(sequence='ACGTACGTTAGCATCG', k=4, compression=2.0, verbose=False)
print(result["status"])
```

**Example output:**

```text
completed
```

## plot_sequence_alignment

**Purpose：** Draw an alignment to an image file and return the file path with the plotted data.

**Method：** Draw an alignment to an image file and return the file path with the plotted data.

``coloring`` selects the symbol plotter: ``type`` colours by symbol type
(``plot_alignment_type_based``), ``similarity`` by column similarity
(``plot_alignment_similarity_based``), ``values`` shades each symbol by a per-position value
(``plot_alignment_array``; ``values`` is a columns-by-sequences matrix of floats).

**API：** `dnakit.alignment.biotite_p30_functions.plot_sequence_alignment(aligned: 'Sequence[str]', output_path: 'str', *, seq_type: 'str' = 'dna', coloring: 'str' = 'type', values: 'Sequence[Sequence[float]] | None' = None, labels: 'Sequence[str] | None' = None, symbols_per_line: 'int' = 50, show_numbers: 'bool' = False, dpi: 'int' = 150, verbose: bool = True) -> 'dict[str, Any]'`

**Input：** Use the explicitly typed parameters in the signature; optional settings retain the documented upstream defaults.

**Output：** JSON-compatible result with status, computed values or image path, and provenance.

**Example code:**

```python
from dnakit.alignment.biotite_p30_functions import plot_sequence_alignment
result = plot_sequence_alignment(aligned=['ACGTACGT', 'AC-TACGT', 'ACGTTCGT'], output_path='plot_sequence_alignment.png', coloring='type', verbose=False)
print(result["status"])
```

**Example output:**

```text
completed
```

## seed_extension_ungapped

**Purpose：** Extend a seed without gaps in both sequences with the X-drop rule (``SeedExtension``).

**Method：** Extend a seed without gaps in both sequences with the X-drop rule (``SeedExtension``).

``seed`` is ``[position_in_sequence_1, position_in_sequence_2]``; extension stops when the
score falls ``threshold`` below the best score seen. ``direction`` is ``both``,
``upstream`` or ``downstream``.

**API：** `dnakit.alignment.biotite_p30_functions.seed_extension_ungapped(sequence_1: 'str', sequence_2: 'str', *, seed: 'Sequence[int]', threshold: 'int', seq_type: 'str' = 'dna', matrix: 'str | None' = None, match: 'int' = 1, mismatch: 'int' = -1, direction: 'str' = 'both', verbose: bool = True) -> 'dict[str, Any]'`

**Input：** Use the explicitly typed parameters in the signature; optional settings retain the documented upstream defaults.

**Output：** JSON-compatible result with status, computed values or image path, and provenance.

**Example code:**

```python
from dnakit.alignment.biotite_p30_functions import seed_extension_ungapped
result = seed_extension_ungapped(sequence_1='ACGTACGTTAGC', sequence_2='TTACGTACGTAAGC', seed=[4, 6], threshold=5, verbose=False)
print(result["status"])
```

**Example output:**

```text
completed
```

## sequence_profile_from_alignment

**Purpose：** Count the symbols of every alignment column (a position-frequency matrix).

**Method：** Count the symbols of every alignment column (a position-frequency matrix).

``aligned`` are equal-length gapped strings. Returns the symbols (including the gap column)
and one count row per alignment column (``SequenceProfile.from_alignment``).

**API：** `dnakit.alignment.biotite_p30_functions.sequence_profile_from_alignment(aligned: 'Sequence[str]', *, seq_type: 'str' = 'dna', verbose: bool = True) -> 'dict[str, Any]'`

**Input：** Use the explicitly typed parameters in the signature; optional settings retain the documented upstream defaults.

**Output：** JSON-compatible result with status, computed values or image path, and provenance.

**Example code:**

```python
from dnakit.alignment.biotite_p30_functions import sequence_profile_from_alignment
result = sequence_profile_from_alignment(aligned=['ACGTACGT', 'AC-TACGT', 'ACGTTCGT'], verbose=False)
print(result["status"])
```

**Example output:**

```text
completed
```

## sequence_profile_probability

**Purpose：** Probability of ``sequence`` under the profile of an alignment.

**Method：** Probability of ``sequence`` under the profile of an alignment.

The profile is built from ``aligned`` (``SequenceProfile.from_alignment``); the sequence must
have one symbol per alignment column. ``pseudocount`` is added to every symbol count.

**API：** `dnakit.alignment.biotite_p30_functions.sequence_profile_probability(aligned: 'Sequence[str]', sequence: 'str', *, seq_type: 'str' = 'dna', pseudocount: 'int' = 0, verbose: bool = True) -> 'dict[str, Any]'`

**Input：** Use the explicitly typed parameters in the signature; optional settings retain the documented upstream defaults.

**Output：** JSON-compatible result with status, computed values or image path, and provenance.

**Example code:**

```python
from dnakit.alignment.biotite_p30_functions import sequence_profile_probability
result = sequence_profile_probability(aligned=['ACGTACGT', 'AC-TACGT', 'ACGTTCGT'], sequence='ACGTACGT', verbose=False)
print(result["status"])
```

**Example output:**

```text
completed
```

## sequence_profile_score

**Purpose：** Log-odds score of ``sequence`` against the weight matrix of an alignment profile.

**Method：** Log-odds score of ``sequence`` against the weight matrix of an alignment profile.

``background_frequencies`` is one float or one frequency per alphabet symbol; None means
uniform. ``pseudocount`` is added to every symbol count (``SequenceProfile.sequence_score``).

**API：** `dnakit.alignment.biotite_p30_functions.sequence_profile_score(aligned: 'Sequence[str]', sequence: 'str', *, seq_type: 'str' = 'dna', background_frequencies: 'float | Sequence[float] | None' = None, pseudocount: 'int' = 0, verbose: bool = True) -> 'dict[str, Any]'`

**Input：** Use the explicitly typed parameters in the signature; optional settings retain the documented upstream defaults.

**Output：** JSON-compatible result with status, computed values or image path, and provenance.

**Example code:**

```python
from dnakit.alignment.biotite_p30_functions import sequence_profile_score
result = sequence_profile_score(aligned=['ACGTACGT', 'AC-TACGT', 'ACGTTCGT'], sequence='ACGTACGT', pseudocount=1, verbose=False)
print(result["status"])
```

**Example output:**

```text
completed
```

## similar_kmers

**Purpose：** List all k-mers whose substitution score against ``kmer`` reaches ``threshold``.

**Method：** List all k-mers whose substitution score against ``kmer`` reaches ``threshold``.

Uses ``ScoreThresholdRule.similar_kmers``. ``matrix`` is a Biotite matrix name (default NUC
for DNA, BLOSUM62 for protein).

**API：** `dnakit.alignment.biotite_p30_functions.similar_kmers(kmer: 'str', *, threshold: 'int', seq_type: 'str' = 'dna', matrix: 'str | None' = None, max_results: 'int' = 10000, verbose: bool = True) -> 'dict[str, Any]'`

**Input：** Use the explicitly typed parameters in the signature; optional settings retain the documented upstream defaults.

**Output：** JSON-compatible result with status, computed values or image path, and provenance.

**Example code:**

```python
from dnakit.alignment.biotite_p30_functions import similar_kmers
result = similar_kmers(kmer='ACGT', threshold=3, verbose=False)
print(result["status"])
```

**Example output:**

```text
completed
```

## syncmer_selection

**Purpose：** Select the syncmers (k-mers whose smallest s-mer is at one of the offsets) of a sequence.

**Method：** Select the syncmers (k-mers whose smallest s-mer is at one of the offsets) of a sequence.

Uses ``CachedSyncmerSelector``. ``s`` is the s-mer length (s < k); ``offset`` the allowed
positions of the minimal s-mer inside the k-mer. Returns positions and k-mer strings.

**API：** `dnakit.alignment.biotite_p30_functions.syncmer_selection(sequence: 'str', *, k: 'int' = 5, s: 'int' = 3, seq_type: 'str' = 'dna', offset: 'Sequence[int]' = (0,), verbose: bool = True) -> 'dict[str, Any]'`

**Input：** Use the explicitly typed parameters in the signature; optional settings retain the documented upstream defaults.

**Output：** JSON-compatible result with status, computed values or image path, and provenance.

**Example code:**

```python
from dnakit.alignment.biotite_p30_functions import syncmer_selection
result = syncmer_selection(sequence='ACGTACGTTAGCATCG', k=4, s=2, verbose=False)
print(result["status"])
```

**Example output:**

```text
completed
```

## write_a3m_alignment

**Purpose：** Write an alignment of one query against targets as an A3M file (``set_a3m_alignments``).

**Method：** Write an alignment of one query against targets as an A3M file (``set_a3m_alignments``).

``aligned[0]`` is the query and every other string a target aligned to it. A3M lower-cases
target residues that are inserted relative to the query and drops the gaps of the query.

**API：** `dnakit.alignment.biotite_p30_functions.write_a3m_alignment(aligned: 'Sequence[str]', output_path: 'str', *, seq_type: 'str' = 'protein', query_label: 'str' = 'query', target_labels: 'Sequence[str] | None' = None, verbose: bool = True) -> 'dict[str, Any]'`

**Input：** Use the explicitly typed parameters in the signature; optional settings retain the documented upstream defaults.

**Output：** JSON-compatible result with status, computed values or image path, and provenance.

**Example code:**

```python
from dnakit.alignment.biotite_p30_functions import write_a3m_alignment
result = write_a3m_alignment(aligned=['ACGTACGT', 'AC-TACGT', 'ACGTTCGT'], output_path='write_a3m_alignment.a3m', seq_type='dna', verbose=False)
print(result["status"])
```

**Example output:**

```text
completed
```

