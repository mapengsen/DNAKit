# Seqeron Sequence MCP (`seqeron-sequence`)

DNAKit's MCP tool index keeps **30** DNA/gene tools of this server (35 listed, read on 2026-10-16). They run once you install the server; DNAKit does not bundle them. See [calling a tool that was found](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found).

| Item | Value |
| --- | --- |
| Source | [AlexanderV/Seqeron @ `7ee2db6`](https://github.com/AlexanderV/Seqeron/tree/7ee2db67549fbc97b82b73e4a38b5cf73170737b) |
| License | MIT |
| Command | `dotnet run --project src/Seqeron/Mcp/Seqeron.Mcp.Sequence` |
| Transports | stdio |
| Installation | Source-based: clone the pinned commit, install the .NET 10 SDK (net10.0 target), set workdir to the checkout root and executable_path to the dotnet executable; the first dotnet run builds the project. No environment variables are read. This entry covers only the Seqeron.Mcp.Sequence server; the repository ships ten more servers (Parsers, Alignment, Analysis, Annotation, MolTools, Population, ...) that are not cataloged. |
| Effects | Pure local computation on sequences passed as tool arguments; the Sequence server's code makes no network calls and reads or writes no files. |

| Tool | Category | Summary |
| --- | --- | --- |
| [`complement_base`](#complement_base) | `sequence` | Returns the Watson-Crick complement of a single DNA or RNA nucleotide base. |
| [`complexity_compression_ratio`](#complexity_compression_ratio) | `sequence` | Computes normalized Lempel-Ziv (LZ76) complexity of a sequence; lower values indicate more repetitive sequence. |
| [`complexity_dust_score`](#complexity_dust_score) | `sequence` | Calculates the DUST low-complexity score of a DNA sequence from triplet counts (as used by dustmasker/sdust). |
| [`complexity_kmer_entropy`](#complexity_kmer_entropy) | `sequence` | Calculates Shannon entropy of overlapping k-mers in a DNA sequence with optional bias correction and normalization. |
| [`complexity_linguistic`](#complexity_linguistic) | `sequence` | Calculates DNA linguistic complexity as the ratio of observed to possible subwords up to a maximum word length. |
| [`complexity_mask_low`](#complexity_mask_low) | `sequence` | Masks low-complexity regions in a DNA sequence with SDUST or dustmasker-compatible DUST, with hard or soft masking. |
| [`complexity_shannon`](#complexity_shannon) | `sequence` | Calculates Shannon entropy in bits per base of a DNA sequence (maximum 2 bits). |
| [`dna_reverse_complement`](#dna_reverse_complement) | `sequence` | Returns the reverse complement of a DNA sequence. |
| [`dna_validate`](#dna_validate) | `sequence` | Validates that a sequence contains only A, C, G and T and reports the result. |
| [`gc_content`](#gc_content) | `sequence` | Calculates the GC content percentage of a DNA or RNA sequence. |
| [`is_valid_dna`](#is_valid_dna) | `sequence` | Quickly checks whether a sequence contains only A, T, G and C characters. |
| [`is_valid_rna`](#is_valid_rna) | `sequence` | Quickly checks whether a sequence contains only A, U, G and C characters. |
| [`iupac_code`](#iupac_code) | `sequence` | Returns the IUPAC ambiguity code representing a given set of nucleotide bases. |
| [`iupac_match`](#iupac_match) | `sequence` | Checks whether two IUPAC nucleotide codes can represent the same base. |
| [`iupac_matches`](#iupac_matches) | `sequence` | Checks whether a specific nucleotide matches an IUPAC ambiguity code. |
| [`kmer_analyze`](#kmer_analyze) | `sequence` | Computes Jellyfish-style k-mer statistics (total, distinct, singletons, counts, entropy) with optional canonical and count filters. |
| [`kmer_count`](#kmer_count) | `sequence` | Counts k-mer frequencies in a sequence with optional ACGT-only and canonical (reverse-complement) modes. |
| [`kmer_distance`](#kmer_distance) | `sequence` | Calculates an alignment-free k-mer distance between two sequences using Euclidean, d2star, Jensen-Shannon or other metrics. |
| [`kmer_entropy`](#kmer_entropy) | `sequence` | Calculates Shannon entropy of k-mer frequencies in a sequence. |
| [`linguistic_complexity`](#linguistic_complexity) | `sequence` | Calculates linguistic complexity (0-1) of a sequence from k-mer diversity up to a maximum k. |
| [`melting_temperature`](#melting_temperature) | `sequence` | Calculates DNA melting temperature (Tm) using the Wallace rule or a GC-content formula. |
| [`molecular_weight_nucleotide`](#molecular_weight_nucleotide) | `sequence` | Calculates the molecular weight in daltons of a DNA or RNA sequence. |
| [`nucleotide_composition`](#nucleotide_composition) | `sequence` | Calculates nucleotide counts (A, T, G, C, U) and GC content of a DNA or RNA sequence. |
| [`rna_from_dna`](#rna_from_dna) | `sequence` | Transcribes a DNA sequence to RNA by replacing T with U. |
| [`rna_validate`](#rna_validate) | `sequence` | Validates that a sequence contains only A, C, G and U and reports the result. |
| [`shannon_entropy`](#shannon_entropy) | `sequence` | Calculates the Shannon entropy of a sequence's symbol distribution. |
| [`summarize_sequence`](#summarize_sequence) | `sequence` | Generates summary statistics for a DNA or RNA sequence: composition, GC content, entropy, complexity and melting temperature. |
| [`thermodynamics`](#thermodynamics) | `sequence` | Calculates nearest-neighbor thermodynamic parameters (dH, dS, dG, Tm) of a DNA duplex at given Na+ and primer concentrations. |
| [`translate_dna`](#translate_dna) | `sequence` | Translates a DNA sequence to protein with the standard genetic code in a chosen reading frame, optionally to the first stop. |
| [`translate_rna`](#translate_rna) | `sequence` | Translates an RNA sequence to protein with the standard genetic code in a chosen reading frame, optionally to the first stop. |

## `complement_base` {#complement_base}

Returns the Watson-Crick complement of a single DNA or RNA nucleotide base.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `nucleotide`
- **Keywords:** complement, nucleotide, base pairing, watson-crick

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `nucleotide` | string | yes | The nucleotide base (A, T, G, C, or U) |

??? note "Upstream description"

    Get the Watson-Crick complement of a single nucleotide base (A↔T, C↔G for DNA; A↔U for RNA).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'seqeron-sequence',
    'complement_base',
    {'nucleotide': '...'},
)
```

## `complexity_compression_ratio` {#complexity_compression_ratio}

Computes normalized Lempel-Ziv (LZ76) complexity of a sequence; lower values indicate more repetitive sequence.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `sequence`
- **Keywords:** lempel-ziv, lz76, sequence complexity, repeats

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sequence` | string | yes | The sequence to analyze |

??? note "Upstream description"

    Normalized Lempel-Ziv (LZ76) complexity c/(n/log_b n) (Lempel &amp; Ziv 1976; Zhang et al. 2009): c = number of components of the LZ76 exhaustive history, n = length, b = number of distinct symbols (at least 2). Not a compressor's compression ratio (the tool name is historical). About 1 for random sequences (finite sequences may exceed 1); lower values indicate more repetitive/less complex sequences.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'seqeron-sequence',
    'complexity_compression_ratio',
    {'sequence': '...'},
)
```

## `complexity_dust_score` {#complexity_dust_score}

Calculates the DUST low-complexity score of a DNA sequence from triplet counts (as used by dustmasker/sdust).

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `sequence`
- **Keywords:** dust, low complexity, dustmasker, triplets

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sequence` | string | yes | The DNA sequence to analyze |
| `wordSize` | integer | no | Word size; must be 3 (DUST is defined for triplets only). Kept for compatibility. Default `3` |

??? note "Upstream description"

    Calculate the DUST low-complexity score (Morgulis et al. 2006; the score thresholded by NCBI dustmasker and lh3/sdust): sum over triplets of c(c-1)/2 divided by (number of triplets - 1). Higher scores indicate lower complexity.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'seqeron-sequence',
    'complexity_dust_score',
    {'sequence': '...'},
)
```

## `complexity_kmer_entropy` {#complexity_kmer_entropy}

Calculates Shannon entropy of overlapping k-mers in a DNA sequence with optional bias correction and normalization.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `sequence`
- **Keywords:** k-mer, shannon entropy, block entropy, miller-madow, sequence complexity

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sequence` | string | yes | The DNA sequence to analyze |
| `k` | integer | no | K-mer size (default: 2 for dinucleotides) Default `2` |
| `correction` | string | no | Bias correction: 'none' (plug-in, default), 'millerMadow' (+ (D-1)/(2N) nats, D = observed k-mers), 'grassberger' (Grassberger 2003: ln N - (1/N) sum n_i G(n_i)). Default `none` |
| `normalize` | boolean | no | Divide the (corrected) entropy by log2 N, N = number of k-mers (0 when N &lt;= 1). Default false. Default `False` |

??? note "Upstream description"

    Calculate the Shannon entropy (bits) of the overlapping k-mer distribution (block entropy) for DNA complexity analysis; optional finite-sample bias correction (Miller-Madow 1955 or Grassberger 2003) and normalisation by log2 of the number of k-mers (BBTools EntropyTracker 0-1 scale).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'seqeron-sequence',
    'complexity_kmer_entropy',
    {'sequence': '...'},
)
```

## `complexity_linguistic` {#complexity_linguistic}

Calculates DNA linguistic complexity as the ratio of observed to possible subwords up to a maximum word length.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `sequence`
- **Keywords:** linguistic complexity, subwords, sequence complexity

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sequence` | string | yes | The DNA sequence to analyze |
| `maxWordLength` | integer | no | Maximum word length to consider (default: 10) Default `10` |
| `alphabetSize` | integer \| null | no | Optional fixed alphabet size a (Troyanskaya et al. 2002; 4 for DNA as in Rosalind LING); must be &gt;= the number of distinct symbols. Default: inferred from the sequence. Default `None` |

??? note "Upstream description"

    Calculate DNA linguistic complexity as ratio of observed to possible subwords. LC = 1.0 for maximum complexity.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'seqeron-sequence',
    'complexity_linguistic',
    {'sequence': '...'},
)
```

## `complexity_mask_low` {#complexity_mask_low}

Masks low-complexity regions in a DNA sequence with SDUST or dustmasker-compatible DUST, with hard or soft masking.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `sequence`
- **Keywords:** dust, sdust, dustmasker, low complexity masking, soft mask

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sequence` | string | yes | The DNA sequence to mask (A/C/G/T plus IUPAC codes such as N; case-insensitive) |
| `windowSize` | integer | no | SDUST window length in bases (default: 64, &gt;= 3) Default `64` |
| `threshold` | number | no | DUST threshold; intervals scoring strictly above it are masked (default: 2.0) Default `2.0` |
| `maskChar` | string | no | Character to use for masking (default: 'N'; ignored when softMask is true) Default `N` |
| `linker` | integer | no | dustmasker linker: merge masked intervals separated by fewer than linker unmasked bases (1-32, default: 1 = sdust/dustmasker default) Default `1` |
| `softMask` | boolean | no | Soft-mask: lower-case masked bases, upper-case the rest (dustmasker -outfmt fasta). Default: false Default `False` |
| `engine` | string | no | DUST engine: 'sdust' (default: lh3/sdust port, non-ACGT symbols split the scan) or 'dustmasker' (exact NCBI dustmasker 2.12.0 parity: IUPAC codes scanned as bases, only N runs longer than the window plus leading/trailing N runs cut the scan and are masked; window 8-64, 10*threshold an integer 2-64). Default `sdust` |

??? note "Upstream description"

    Mask low-complexity regions in a DNA sequence with the symmetric DUST algorithm (SDUST; Morgulis et al. 2006, identical to lh3/sdust). N and other IUPAC codes are accepted; each maximal A/C/G/T run is scanned independently (sdust's documented contract: "N effectively breaks input into pieces of independent sequences"; sdust's code itself carries its scoring window across N). Optional dustmasker linker merge and soft (lower-case) masking. engine='dustmasker' reproduces NCBI dustmasker 2.12.0 output exactly (symdust core, IUPAC scanned as bases, long/terminal N runs masked).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'seqeron-sequence',
    'complexity_mask_low',
    {'sequence': '...'},
)
```

## `complexity_shannon` {#complexity_shannon}

Calculates Shannon entropy in bits per base of a DNA sequence (maximum 2 bits).

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `sequence`
- **Keywords:** shannon entropy, bits per base, sequence complexity

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sequence` | string | yes | The DNA sequence to analyze |

??? note "Upstream description"

    Calculate DNA Shannon entropy (bits per base). Maximum entropy for DNA is 2 bits.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'seqeron-sequence',
    'complexity_shannon',
    {'sequence': '...'},
)
```

## `dna_reverse_complement` {#dna_reverse_complement}

Returns the reverse complement of a DNA sequence.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `sequence`
- **Keywords:** reverse complement, strand, complement

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sequence` | string | yes | The DNA sequence |

??? note "Upstream description"

    Get the reverse complement of a DNA sequence. A↔T, C↔G, then reversed.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'seqeron-sequence',
    'dna_reverse_complement',
    {'sequence': '...'},
)
```

## `dna_validate` {#dna_validate}

Validates that a sequence contains only A, C, G and T and reports the result.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `sequence`
- **Keywords:** validation, nucleotide alphabet, acgt

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sequence` | string | yes | The DNA sequence to validate |

??? note "Upstream description"

    Validate a DNA sequence. Returns whether the sequence contains only valid nucleotides (A, C, G, T).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'seqeron-sequence',
    'dna_validate',
    {'sequence': '...'},
)
```

## `gc_content` {#gc_content}

Calculates the GC content percentage of a DNA or RNA sequence.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `sequence`
- **Keywords:** gc percent, nucleotide composition, base content

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sequence` | string | yes | The DNA or RNA sequence |

??? note "Upstream description"

    Calculate the GC content (percentage of G and C nucleotides) of a DNA/RNA sequence.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'seqeron-sequence',
    'gc_content',
    {'sequence': '...'},
)
```

## `is_valid_dna` {#is_valid_dna}

Quickly checks whether a sequence contains only A, T, G and C characters.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `sequence`
- **Keywords:** validation, dna alphabet, acgt

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sequence` | string | yes | The sequence to validate |

??? note "Upstream description"

    Quick check if a sequence contains only valid DNA characters (A, T, G, C). Faster than dna_validate but returns less information.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'seqeron-sequence',
    'is_valid_dna',
    {'sequence': '...'},
)
```

## `is_valid_rna` {#is_valid_rna}

Quickly checks whether a sequence contains only A, U, G and C characters.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `sequence`
- **Keywords:** validation, rna alphabet, acgu

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sequence` | string | yes | The sequence to validate |

??? note "Upstream description"

    Quick check if a sequence contains only valid RNA characters (A, U, G, C). Faster than rna_validate but returns less information.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'seqeron-sequence',
    'is_valid_rna',
    {'sequence': '...'},
)
```

## `iupac_code` {#iupac_code}

Returns the IUPAC ambiguity code representing a given set of nucleotide bases.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `bases`
- **Keywords:** iupac, ambiguity code, degenerate base, purine

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `bases` | string | yes | Nucleotide bases to encode (e.g., 'AG' for purine R) |

??? note "Upstream description"

    Get the IUPAC ambiguity code that represents a set of nucleotide bases.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'seqeron-sequence',
    'iupac_code',
    {'bases': '...'},
)
```

## `iupac_match` {#iupac_match}

Checks whether two IUPAC nucleotide codes can represent the same base.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `code1`, `code2`
- **Keywords:** iupac, ambiguity code, degenerate base

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `code1` | string | yes | First IUPAC code |
| `code2` | string | yes | Second IUPAC code |

??? note "Upstream description"

    Check if two IUPAC codes can represent the same nucleotide base.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'seqeron-sequence',
    'iupac_match',
    {'code1': '...', 'code2': '...'},
)
```

## `iupac_matches` {#iupac_matches}

Checks whether a specific nucleotide matches an IUPAC ambiguity code.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `nucleotide`, `iupacCode`
- **Keywords:** iupac, ambiguity code, nucleotide match

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `nucleotide` | string | yes | The nucleotide to check (A, C, G, T) |
| `iupacCode` | string | yes | The IUPAC code to match against |

??? note "Upstream description"

    Check if a specific nucleotide matches an IUPAC ambiguity code.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'seqeron-sequence',
    'iupac_matches',
    {'nucleotide': '...', 'iupacCode': '...'},
)
```

## `kmer_analyze` {#kmer_analyze}

Computes Jellyfish-style k-mer statistics (total, distinct, singletons, counts, entropy) with optional canonical and count filters.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `sequence`
- **Keywords:** k-mer, jellyfish, canonical k-mers, distinct k-mers, singletons

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sequence` | string | yes | The sequence to analyze |
| `k` | integer | no | K-mer length (default: 3) Default `3` |
| `lowerCount` | integer | no | Ignore k-mers with count below this value (Jellyfish stats -L; default 0). Default `0` |
| `upperCount` | integer | no | Ignore k-mers with count above this value (Jellyfish stats -U; default unbounded). Default `2147483647` |
| `canonical` | boolean | no | Canonical k-mers min(k-mer, reverse complement) (jellyfish count -C); implies acgtOnly. Default false. Default `False` |
| `acgtOnly` | boolean | no | Skip windows containing a non-ACGT symbol (Jellyfish convention). Default false. Default `False` |

??? note "Upstream description"

    Comprehensive k-mer analysis (Jellyfish stats fields): total, distinct (uniqueKmers/distinctKmers), singleton (count-1) k-mers, min/max/mean count, and Shannon entropy; optional Jellyfish -L/-U count filters and canonical (count -C) / acgtOnly modes (same as analyze_kmers).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'seqeron-sequence',
    'kmer_analyze',
    {'sequence': '...'},
)
```

## `kmer_count` {#kmer_count}

Counts k-mer frequencies in a sequence with optional ACGT-only and canonical (reverse-complement) modes.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `sequence`
- **Keywords:** k-mer, counting, canonical, jellyfish, frequency

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sequence` | string | yes | The sequence to analyze |
| `k` | integer | no | K-mer length (default: 3) Default `3` |
| `canonical` | boolean | no | Canonical counting (jellyfish count -C): key = lexicographically smaller of the k-mer and its reverse complement; implies acgtOnly. Default false. Default `False` |
| `acgtOnly` | boolean | no | Skip every window containing a symbol other than A/C/G/T (case-insensitive), as Jellyfish does. Default false (all symbols counted literally). Default `False` |

??? note "Upstream description"

    Count k-mer (substring of length k) frequencies in a sequence. Returns a dictionary of k-mers and their counts. Optional Jellyfish modes (same as the Analysis server's count_kmers): acgtOnly skips windows containing a non-ACGT symbol; canonical keys each k-mer by min(k-mer, reverse complement) (jellyfish count -C, implies acgtOnly).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'seqeron-sequence',
    'kmer_count',
    {'sequence': '...'},
)
```

## `kmer_distance` {#kmer_distance}

Calculates an alignment-free k-mer distance between two sequences using Euclidean, d2star, Jensen-Shannon or other metrics.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `sequence1`, `sequence2`
- **Keywords:** alignment-free, k-mer, d2star, sequence similarity, jensen-shannon

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sequence1` | string | yes | First sequence |
| `sequence2` | string | yes | Second sequence |
| `k` | integer | no | K-mer length (default: 3) Default `3` |
| `metric` | string | no | Metric: euclidean (default), squared_euclidean_counts, manhattan, chebyshev, canberra, cosine, d2, d2star, d2shepherd (alias d2s), jensen_shannon (alias js), euclidean_counts, ev (alias evolutionary). Default `euclidean` |
| `markovOrder` | integer | no | Background Markov order r (0 &lt;= r &lt; k) for d2star/d2shepherd, or -1 to choose each sequence's order by BIC; default 0. Must be 0 for the other metrics. Default `0` |
| `bothStrands` | boolean | no | d2star/d2shepherd: CAFE -R both-strand mode (count of w + count of its reverse complement, background averaged over both); ev: spaced's reverse-complement mode (seq1 on both strands). Default false. Default `False` |

??? note "Upstream description"

    Calculate k-mer based distance between two sequences using Euclidean distance of k-mer frequencies. Lower values indicate more similar sequences. Optional metric (same as the Analysis server's kmer_distance): euclidean (default), squared_euclidean_counts, manhattan, chebyshev, canberra, cosine, d2 (a similarity), d2star / d2shepherd (background-adjusted, k &lt;= 12, Markov order markovOrder, bothStrands = CAFE -R), jensen_shannon, euclidean_counts, ev (spaced -d EV evolutionary distance with the contiguous pattern 1^k; bothStrands = spaced's reverse-complement mode).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'seqeron-sequence',
    'kmer_distance',
    {'sequence1': '...', 'sequence2': '...'},
)
```

## `kmer_entropy` {#kmer_entropy}

Calculates Shannon entropy of k-mer frequencies in a sequence.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `sequence`
- **Keywords:** k-mer, shannon entropy, sequence complexity

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sequence` | string | yes | The sequence to analyze |
| `k` | integer | no | K-mer length (default: 2) Default `2` |

??? note "Upstream description"

    Calculate Shannon entropy based on k-mer frequencies. Higher values indicate more complexity.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'seqeron-sequence',
    'kmer_entropy',
    {'sequence': '...'},
)
```

## `linguistic_complexity` {#linguistic_complexity}

Calculates linguistic complexity (0-1) of a sequence from k-mer diversity up to a maximum k.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `sequence`
- **Keywords:** linguistic complexity, k-mer diversity, sequence complexity

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sequence` | string | yes | The sequence to analyze |
| `maxK` | integer | no | Maximum k-mer length to consider (default: 6) Default `6` |
| `alphabetSize` | integer \| null | no | Optional fixed alphabet size a for the denominator min(a^i, N-i+1) (4 for DNA as in Rosalind LING, 20 for protein); default: inferred from the sequence. Default `None` |

??? note "Upstream description"

    Calculate linguistic complexity of a sequence based on k-mer diversity. Values range from 0 to 1.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'seqeron-sequence',
    'linguistic_complexity',
    {'sequence': '...'},
)
```

## `melting_temperature` {#melting_temperature}

Calculates DNA melting temperature (Tm) using the Wallace rule or a GC-content formula.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `sequence`
- **Keywords:** tm, melting temperature, wallace rule, oligo, primer

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sequence` | string | yes | The DNA sequence |
| `useWallaceRule` | boolean | no | Use Wallace rule for short oligos (default: true) Default `True` |

??? note "Upstream description"

    Calculate the melting temperature (Tm) of a DNA sequence using Wallace rule or GC formula.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'seqeron-sequence',
    'melting_temperature',
    {'sequence': '...'},
)
```

## `molecular_weight_nucleotide` {#molecular_weight_nucleotide}

Calculates the molecular weight in daltons of a DNA or RNA sequence.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `sequence`
- **Keywords:** molecular weight, daltons, oligonucleotide

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sequence` | string | yes | The DNA or RNA sequence |
| `isDna` | boolean | no | True for DNA, false for RNA (default: true) Default `True` |

??? note "Upstream description"

    Calculate the molecular weight of a DNA or RNA sequence in Daltons (Da).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'seqeron-sequence',
    'molecular_weight_nucleotide',
    {'sequence': '...'},
)
```

## `nucleotide_composition` {#nucleotide_composition}

Calculates nucleotide counts (A, T, G, C, U) and GC content of a DNA or RNA sequence.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `sequence`
- **Keywords:** base composition, nucleotide counts, gc content

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sequence` | string | yes | The DNA or RNA sequence to analyze |

??? note "Upstream description"

    Calculate nucleotide composition (A, T, G, C, U counts) and GC content of a DNA/RNA sequence.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'seqeron-sequence',
    'nucleotide_composition',
    {'sequence': '...'},
)
```

## `rna_from_dna` {#rna_from_dna}

Transcribes a DNA sequence to RNA by replacing T with U.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `sequence`
- **Keywords:** transcription, t to u, rna

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sequence` | string | yes | The DNA sequence to transcribe |

??? note "Upstream description"

    Transcribe DNA to RNA by replacing T (thymine) with U (uracil).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'seqeron-sequence',
    'rna_from_dna',
    {'sequence': '...'},
)
```

## `rna_validate` {#rna_validate}

Validates that a sequence contains only A, C, G and U and reports the result.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `sequence`
- **Keywords:** validation, rna alphabet, acgu

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sequence` | string | yes | The RNA sequence to validate |

??? note "Upstream description"

    Validate an RNA sequence. Returns whether the sequence contains only valid nucleotides (A, C, G, U).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'seqeron-sequence',
    'rna_validate',
    {'sequence': '...'},
)
```

## `shannon_entropy` {#shannon_entropy}

Calculates the Shannon entropy of a sequence's symbol distribution.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `sequence`
- **Keywords:** shannon entropy, sequence complexity, randomness

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sequence` | string | yes | The sequence to analyze |

??? note "Upstream description"

    Calculate Shannon entropy of a sequence. Higher values indicate more complexity/randomness.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'seqeron-sequence',
    'shannon_entropy',
    {'sequence': '...'},
)
```

## `summarize_sequence` {#summarize_sequence}

Generates summary statistics for a DNA or RNA sequence: composition, GC content, entropy, complexity and melting temperature.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `sequence`
- **Keywords:** sequence statistics, gc content, entropy, tm, composition

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sequence` | string | yes | The DNA or RNA sequence to analyze |

??? note "Upstream description"

    Generate comprehensive summary statistics for a DNA/RNA sequence including composition, GC content, entropy, complexity, and Tm.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'seqeron-sequence',
    'summarize_sequence',
    {'sequence': '...'},
)
```

## `thermodynamics` {#thermodynamics}

Calculates nearest-neighbor thermodynamic parameters (dH, dS, dG, Tm) of a DNA duplex at given Na+ and primer concentrations.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `sequence`
- **Keywords:** nearest-neighbor, delta g, duplex, tm, primer, oligo

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sequence` | string | yes | The DNA sequence |
| `naConcentration` | number | no | Na+ concentration in M (default: 0.05 = 50mM) Default `0.05` |
| `primerConcentration` | number | no | Primer concentration in M (default: 0.00000025 = 250nM) Default `2.5e-07` |

??? note "Upstream description"

    Calculate thermodynamic properties (ΔH, ΔS, ΔG, Tm) of a DNA duplex using the nearest-neighbor method.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'seqeron-sequence',
    'thermodynamics',
    {'sequence': '...'},
)
```

## `translate_dna` {#translate_dna}

Translates a DNA sequence to protein with the standard genetic code in a chosen reading frame, optionally to the first stop.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `sequence`
- **Keywords:** translation, genetic code, reading frame, codon, orf

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sequence` | string | yes | The DNA sequence to translate |
| `frame` | integer | no | Reading frame (0, 1, or 2, default: 0) Default `0` |
| `toFirstStop` | boolean | no | Stop at first stop codon (default: false) Default `False` |

??? note "Upstream description"

    Translate a DNA sequence to protein using the standard genetic code.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'seqeron-sequence',
    'translate_dna',
    {'sequence': '...'},
)
```

## `translate_rna` {#translate_rna}

Translates an RNA sequence to protein with the standard genetic code in a chosen reading frame, optionally to the first stop.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `sequence`
- **Keywords:** translation, genetic code, reading frame, codon, mrna

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sequence` | string | yes | The RNA sequence to translate |
| `frame` | integer | no | Reading frame (0, 1, or 2, default: 0) Default `0` |
| `toFirstStop` | boolean | no | Stop at first stop codon (default: false) Default `False` |

??? note "Upstream description"

    Translate an RNA sequence to protein using the standard genetic code.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'seqeron-sequence',
    'translate_rna',
    {'sequence': '...'},
)
```
