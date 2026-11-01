# Seqeron Sequence MCP (`seqeron-sequence`)

DNAKit 的 MCP 工具索引收录了这个服务器的 **30** 个 DNA/基因相关工具（服务器共列出 35 个，2026-10-16 读取）。工具在用户自己安装服务器后运行，DNAKit 本身不包含它们。调用方式见 [调用索引中的工具](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found)。下面的工具说明是上游英文原文。

| 项目 | 内容 |
| --- | --- |
| 源码 | [AlexanderV/Seqeron @ `7ee2db6`](https://github.com/AlexanderV/Seqeron/tree/7ee2db67549fbc97b82b73e4a38b5cf73170737b) |
| 许可证 | MIT |
| 启动命令 | `dotnet run --project src/Seqeron/Mcp/Seqeron.Mcp.Sequence` |
| 传输方式 | stdio |
| 安装 | Source-based: clone the pinned commit, install the .NET 10 SDK (net10.0 target), set workdir to the checkout root and executable_path to the dotnet executable; the first dotnet run builds the project. No environment variables are read. This entry covers only the Seqeron.Mcp.Sequence server; the repository ships ten more servers (Parsers, Alignment, Analysis, Annotation, MolTools, Population, ...) that are not cataloged. |
| 副作用 | Pure local computation on sequences passed as tool arguments; the Sequence server's code makes no network calls and reads or writes no files. |

| 工具 | 类别 | 简介 |
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

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`nucleotide`
- **关键词：**complement, nucleotide, base pairing, watson-crick

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `nucleotide` | string | 是 | The nucleotide base (A, T, G, C, or U) |

??? note "上游说明（英文原文）"

    Get the Watson-Crick complement of a single nucleotide base (A↔T, C↔G for DNA; A↔U for RNA).

安装服务器后通过 DNAKit 调用：

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

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`sequence`
- **关键词：**lempel-ziv, lz76, sequence complexity, repeats

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `sequence` | string | 是 | The sequence to analyze |

??? note "上游说明（英文原文）"

    Normalized Lempel-Ziv (LZ76) complexity c/(n/log_b n) (Lempel &amp; Ziv 1976; Zhang et al. 2009): c = number of components of the LZ76 exhaustive history, n = length, b = number of distinct symbols (at least 2). Not a compressor's compression ratio (the tool name is historical). About 1 for random sequences (finite sequences may exceed 1); lower values indicate more repetitive/less complex sequences.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`sequence`
- **关键词：**dust, low complexity, dustmasker, triplets

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `sequence` | string | 是 | The DNA sequence to analyze |
| `wordSize` | integer | 否 | Word size; must be 3 (DUST is defined for triplets only). Kept for compatibility. 默认 `3` |

??? note "上游说明（英文原文）"

    Calculate the DUST low-complexity score (Morgulis et al. 2006; the score thresholded by NCBI dustmasker and lh3/sdust): sum over triplets of c(c-1)/2 divided by (number of triplets - 1). Higher scores indicate lower complexity.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`sequence`
- **关键词：**k-mer, shannon entropy, block entropy, miller-madow, sequence complexity

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `sequence` | string | 是 | The DNA sequence to analyze |
| `k` | integer | 否 | K-mer size (default: 2 for dinucleotides) 默认 `2` |
| `correction` | string | 否 | Bias correction: 'none' (plug-in, default), 'millerMadow' (+ (D-1)/(2N) nats, D = observed k-mers), 'grassberger' (Grassberger 2003: ln N - (1/N) sum n_i G(n_i)). 默认 `none` |
| `normalize` | boolean | 否 | Divide the (corrected) entropy by log2 N, N = number of k-mers (0 when N &lt;= 1). Default false. 默认 `False` |

??? note "上游说明（英文原文）"

    Calculate the Shannon entropy (bits) of the overlapping k-mer distribution (block entropy) for DNA complexity analysis; optional finite-sample bias correction (Miller-Madow 1955 or Grassberger 2003) and normalisation by log2 of the number of k-mers (BBTools EntropyTracker 0-1 scale).

安装服务器后通过 DNAKit 调用：

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

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`sequence`
- **关键词：**linguistic complexity, subwords, sequence complexity

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `sequence` | string | 是 | The DNA sequence to analyze |
| `maxWordLength` | integer | 否 | Maximum word length to consider (default: 10) 默认 `10` |
| `alphabetSize` | integer \| null | 否 | Optional fixed alphabet size a (Troyanskaya et al. 2002; 4 for DNA as in Rosalind LING); must be &gt;= the number of distinct symbols. Default: inferred from the sequence. 默认 `None` |

??? note "上游说明（英文原文）"

    Calculate DNA linguistic complexity as ratio of observed to possible subwords. LC = 1.0 for maximum complexity.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`sequence`
- **关键词：**dust, sdust, dustmasker, low complexity masking, soft mask

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `sequence` | string | 是 | The DNA sequence to mask (A/C/G/T plus IUPAC codes such as N; case-insensitive) |
| `windowSize` | integer | 否 | SDUST window length in bases (default: 64, &gt;= 3) 默认 `64` |
| `threshold` | number | 否 | DUST threshold; intervals scoring strictly above it are masked (default: 2.0) 默认 `2.0` |
| `maskChar` | string | 否 | Character to use for masking (default: 'N'; ignored when softMask is true) 默认 `N` |
| `linker` | integer | 否 | dustmasker linker: merge masked intervals separated by fewer than linker unmasked bases (1-32, default: 1 = sdust/dustmasker default) 默认 `1` |
| `softMask` | boolean | 否 | Soft-mask: lower-case masked bases, upper-case the rest (dustmasker -outfmt fasta). Default: false 默认 `False` |
| `engine` | string | 否 | DUST engine: 'sdust' (default: lh3/sdust port, non-ACGT symbols split the scan) or 'dustmasker' (exact NCBI dustmasker 2.12.0 parity: IUPAC codes scanned as bases, only N runs longer than the window plus leading/trailing N runs cut the scan and are masked; window 8-64, 10*threshold an integer 2-64). 默认 `sdust` |

??? note "上游说明（英文原文）"

    Mask low-complexity regions in a DNA sequence with the symmetric DUST algorithm (SDUST; Morgulis et al. 2006, identical to lh3/sdust). N and other IUPAC codes are accepted; each maximal A/C/G/T run is scanned independently (sdust's documented contract: "N effectively breaks input into pieces of independent sequences"; sdust's code itself carries its scoring window across N). Optional dustmasker linker merge and soft (lower-case) masking. engine='dustmasker' reproduces NCBI dustmasker 2.12.0 output exactly (symdust core, IUPAC scanned as bases, long/terminal N runs masked).

安装服务器后通过 DNAKit 调用：

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

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`sequence`
- **关键词：**shannon entropy, bits per base, sequence complexity

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `sequence` | string | 是 | The DNA sequence to analyze |

??? note "上游说明（英文原文）"

    Calculate DNA Shannon entropy (bits per base). Maximum entropy for DNA is 2 bits.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`sequence`
- **关键词：**reverse complement, strand, complement

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `sequence` | string | 是 | The DNA sequence |

??? note "上游说明（英文原文）"

    Get the reverse complement of a DNA sequence. A↔T, C↔G, then reversed.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`sequence`
- **关键词：**validation, nucleotide alphabet, acgt

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `sequence` | string | 是 | The DNA sequence to validate |

??? note "上游说明（英文原文）"

    Validate a DNA sequence. Returns whether the sequence contains only valid nucleotides (A, C, G, T).

安装服务器后通过 DNAKit 调用：

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

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`sequence`
- **关键词：**gc percent, nucleotide composition, base content

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `sequence` | string | 是 | The DNA or RNA sequence |

??? note "上游说明（英文原文）"

    Calculate the GC content (percentage of G and C nucleotides) of a DNA/RNA sequence.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`sequence`
- **关键词：**validation, dna alphabet, acgt

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `sequence` | string | 是 | The sequence to validate |

??? note "上游说明（英文原文）"

    Quick check if a sequence contains only valid DNA characters (A, T, G, C). Faster than dna_validate but returns less information.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`sequence`
- **关键词：**validation, rna alphabet, acgu

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `sequence` | string | 是 | The sequence to validate |

??? note "上游说明（英文原文）"

    Quick check if a sequence contains only valid RNA characters (A, U, G, C). Faster than rna_validate but returns less information.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`bases`
- **关键词：**iupac, ambiguity code, degenerate base, purine

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `bases` | string | 是 | Nucleotide bases to encode (e.g., 'AG' for purine R) |

??? note "上游说明（英文原文）"

    Get the IUPAC ambiguity code that represents a set of nucleotide bases.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`code1`, `code2`
- **关键词：**iupac, ambiguity code, degenerate base

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `code1` | string | 是 | First IUPAC code |
| `code2` | string | 是 | Second IUPAC code |

??? note "上游说明（英文原文）"

    Check if two IUPAC codes can represent the same nucleotide base.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`nucleotide`, `iupacCode`
- **关键词：**iupac, ambiguity code, nucleotide match

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `nucleotide` | string | 是 | The nucleotide to check (A, C, G, T) |
| `iupacCode` | string | 是 | The IUPAC code to match against |

??? note "上游说明（英文原文）"

    Check if a specific nucleotide matches an IUPAC ambiguity code.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`sequence`
- **关键词：**k-mer, jellyfish, canonical k-mers, distinct k-mers, singletons

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `sequence` | string | 是 | The sequence to analyze |
| `k` | integer | 否 | K-mer length (default: 3) 默认 `3` |
| `lowerCount` | integer | 否 | Ignore k-mers with count below this value (Jellyfish stats -L; default 0). 默认 `0` |
| `upperCount` | integer | 否 | Ignore k-mers with count above this value (Jellyfish stats -U; default unbounded). 默认 `2147483647` |
| `canonical` | boolean | 否 | Canonical k-mers min(k-mer, reverse complement) (jellyfish count -C); implies acgtOnly. Default false. 默认 `False` |
| `acgtOnly` | boolean | 否 | Skip windows containing a non-ACGT symbol (Jellyfish convention). Default false. 默认 `False` |

??? note "上游说明（英文原文）"

    Comprehensive k-mer analysis (Jellyfish stats fields): total, distinct (uniqueKmers/distinctKmers), singleton (count-1) k-mers, min/max/mean count, and Shannon entropy; optional Jellyfish -L/-U count filters and canonical (count -C) / acgtOnly modes (same as analyze_kmers).

安装服务器后通过 DNAKit 调用：

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

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`sequence`
- **关键词：**k-mer, counting, canonical, jellyfish, frequency

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `sequence` | string | 是 | The sequence to analyze |
| `k` | integer | 否 | K-mer length (default: 3) 默认 `3` |
| `canonical` | boolean | 否 | Canonical counting (jellyfish count -C): key = lexicographically smaller of the k-mer and its reverse complement; implies acgtOnly. Default false. 默认 `False` |
| `acgtOnly` | boolean | 否 | Skip every window containing a symbol other than A/C/G/T (case-insensitive), as Jellyfish does. Default false (all symbols counted literally). 默认 `False` |

??? note "上游说明（英文原文）"

    Count k-mer (substring of length k) frequencies in a sequence. Returns a dictionary of k-mers and their counts. Optional Jellyfish modes (same as the Analysis server's count_kmers): acgtOnly skips windows containing a non-ACGT symbol; canonical keys each k-mer by min(k-mer, reverse complement) (jellyfish count -C, implies acgtOnly).

安装服务器后通过 DNAKit 调用：

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

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`sequence1`, `sequence2`
- **关键词：**alignment-free, k-mer, d2star, sequence similarity, jensen-shannon

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `sequence1` | string | 是 | First sequence |
| `sequence2` | string | 是 | Second sequence |
| `k` | integer | 否 | K-mer length (default: 3) 默认 `3` |
| `metric` | string | 否 | Metric: euclidean (default), squared_euclidean_counts, manhattan, chebyshev, canberra, cosine, d2, d2star, d2shepherd (alias d2s), jensen_shannon (alias js), euclidean_counts, ev (alias evolutionary). 默认 `euclidean` |
| `markovOrder` | integer | 否 | Background Markov order r (0 &lt;= r &lt; k) for d2star/d2shepherd, or -1 to choose each sequence's order by BIC; default 0. Must be 0 for the other metrics. 默认 `0` |
| `bothStrands` | boolean | 否 | d2star/d2shepherd: CAFE -R both-strand mode (count of w + count of its reverse complement, background averaged over both); ev: spaced's reverse-complement mode (seq1 on both strands). Default false. 默认 `False` |

??? note "上游说明（英文原文）"

    Calculate k-mer based distance between two sequences using Euclidean distance of k-mer frequencies. Lower values indicate more similar sequences. Optional metric (same as the Analysis server's kmer_distance): euclidean (default), squared_euclidean_counts, manhattan, chebyshev, canberra, cosine, d2 (a similarity), d2star / d2shepherd (background-adjusted, k &lt;= 12, Markov order markovOrder, bothStrands = CAFE -R), jensen_shannon, euclidean_counts, ev (spaced -d EV evolutionary distance with the contiguous pattern 1^k; bothStrands = spaced's reverse-complement mode).

安装服务器后通过 DNAKit 调用：

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

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`sequence`
- **关键词：**k-mer, shannon entropy, sequence complexity

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `sequence` | string | 是 | The sequence to analyze |
| `k` | integer | 否 | K-mer length (default: 2) 默认 `2` |

??? note "上游说明（英文原文）"

    Calculate Shannon entropy based on k-mer frequencies. Higher values indicate more complexity.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`sequence`
- **关键词：**linguistic complexity, k-mer diversity, sequence complexity

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `sequence` | string | 是 | The sequence to analyze |
| `maxK` | integer | 否 | Maximum k-mer length to consider (default: 6) 默认 `6` |
| `alphabetSize` | integer \| null | 否 | Optional fixed alphabet size a for the denominator min(a^i, N-i+1) (4 for DNA as in Rosalind LING, 20 for protein); default: inferred from the sequence. 默认 `None` |

??? note "上游说明（英文原文）"

    Calculate linguistic complexity of a sequence based on k-mer diversity. Values range from 0 to 1.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`sequence`
- **关键词：**tm, melting temperature, wallace rule, oligo, primer

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `sequence` | string | 是 | The DNA sequence |
| `useWallaceRule` | boolean | 否 | Use Wallace rule for short oligos (default: true) 默认 `True` |

??? note "上游说明（英文原文）"

    Calculate the melting temperature (Tm) of a DNA sequence using Wallace rule or GC formula.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`sequence`
- **关键词：**molecular weight, daltons, oligonucleotide

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `sequence` | string | 是 | The DNA or RNA sequence |
| `isDna` | boolean | 否 | True for DNA, false for RNA (default: true) 默认 `True` |

??? note "上游说明（英文原文）"

    Calculate the molecular weight of a DNA or RNA sequence in Daltons (Da).

安装服务器后通过 DNAKit 调用：

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

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`sequence`
- **关键词：**base composition, nucleotide counts, gc content

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `sequence` | string | 是 | The DNA or RNA sequence to analyze |

??? note "上游说明（英文原文）"

    Calculate nucleotide composition (A, T, G, C, U counts) and GC content of a DNA/RNA sequence.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`sequence`
- **关键词：**transcription, t to u, rna

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `sequence` | string | 是 | The DNA sequence to transcribe |

??? note "上游说明（英文原文）"

    Transcribe DNA to RNA by replacing T (thymine) with U (uracil).

安装服务器后通过 DNAKit 调用：

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

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`sequence`
- **关键词：**validation, rna alphabet, acgu

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `sequence` | string | 是 | The RNA sequence to validate |

??? note "上游说明（英文原文）"

    Validate an RNA sequence. Returns whether the sequence contains only valid nucleotides (A, C, G, U).

安装服务器后通过 DNAKit 调用：

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

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`sequence`
- **关键词：**shannon entropy, sequence complexity, randomness

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `sequence` | string | 是 | The sequence to analyze |

??? note "上游说明（英文原文）"

    Calculate Shannon entropy of a sequence. Higher values indicate more complexity/randomness.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`sequence`
- **关键词：**sequence statistics, gc content, entropy, tm, composition

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `sequence` | string | 是 | The DNA or RNA sequence to analyze |

??? note "上游说明（英文原文）"

    Generate comprehensive summary statistics for a DNA/RNA sequence including composition, GC content, entropy, complexity, and Tm.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`sequence`
- **关键词：**nearest-neighbor, delta g, duplex, tm, primer, oligo

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `sequence` | string | 是 | The DNA sequence |
| `naConcentration` | number | 否 | Na+ concentration in M (default: 0.05 = 50mM) 默认 `0.05` |
| `primerConcentration` | number | 否 | Primer concentration in M (default: 0.00000025 = 250nM) 默认 `2.5e-07` |

??? note "上游说明（英文原文）"

    Calculate thermodynamic properties (ΔH, ΔS, ΔG, Tm) of a DNA duplex using the nearest-neighbor method.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`sequence`
- **关键词：**translation, genetic code, reading frame, codon, orf

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `sequence` | string | 是 | The DNA sequence to translate |
| `frame` | integer | 否 | Reading frame (0, 1, or 2, default: 0) 默认 `0` |
| `toFirstStop` | boolean | 否 | Stop at first stop codon (default: false) 默认 `False` |

??? note "上游说明（英文原文）"

    Translate a DNA sequence to protein using the standard genetic code.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`sequence`
- **关键词：**translation, genetic code, reading frame, codon, mrna

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `sequence` | string | 是 | The RNA sequence to translate |
| `frame` | integer | 否 | Reading frame (0, 1, or 2, default: 0) 默认 `0` |
| `toFirstStop` | boolean | 否 | Stop at first stop codon (default: false) 默认 `False` |

??? note "上游说明（英文原文）"

    Translate an RNA sequence to protein using the standard genetic code.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'seqeron-sequence',
    'translate_rna',
    {'sequence': '...'},
)
```
