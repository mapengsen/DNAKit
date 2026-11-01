# MCP Sequence Simulation (`sequence-simulation`)

DNAKit 的 MCP 工具索引收录了这个服务器的 **5** 个 DNA/基因相关工具（服务器共列出 6 个，2026-10-16 读取）。工具在用户自己安装服务器后运行，DNAKit 本身不包含它们。调用方式见 [调用索引中的工具](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found)。下面的工具说明是上游英文原文。

| 项目 | 内容 |
| --- | --- |
| 源码 | [coreymhudson/mcp-sequence-simulation @ `dea87b2`](https://github.com/coreymhudson/mcp-sequence-simulation/tree/dea87b2a45594ff969eb5f05826642f99d42ac29) |
| 许可证 | BSD-2-Clause |
| 启动命令 | `node dist/server.js` |
| 传输方式 | stdio |
| 安装 | Install separately from the pinned source checkout with npm install &amp;&amp; npm run build (Node.js with npm; package.json declares no engines field), then set workdir to the checkout and keep dist/server.js; no environment variables are used. |
| 副作用 | Pure local computation that returns simulated sequences, trees and FASTQ text in the tool response; the source contains no network calls, file writes or environment access. |

| 工具 | 类别 | 简介 |
| --- | --- | --- |
| [`evolve_sequence`](#evolve_sequence) | `sequence` | Simulates evolution of a starting sequence over generations with mutation rate, population size and selection on a fitness function such as GC content. |
| [`generate_dna_sequence`](#generate_dna_sequence) | `sequence` | Generates random DNA sequences of given length and GC content using random, Markov or codon-biased models, as FASTA or plain text. |
| [`mutate_sequence`](#mutate_sequence) | `sequence` | Applies random substitutions, insertions and deletions to a DNA or protein sequence with a transition/transversion bias for DNA. |
| [`simulate_fastq_file`](#simulate_fastq_file) | `sequence` | Simulates single- or paired-end FASTQ reads from a reference DNA sequence with coverage, insert size, error rate and quality models (NEAT-style). |
| [`simulate_phylogeny`](#simulate_phylogeny) | `sequence` | Simulates a phylogenetic tree and evolves a root sequence along it using JC69, K80, HKY85 or GTR models, outputting FASTA, NEXUS or PHYLIP. |

## `evolve_sequence` {#evolve_sequence}

Simulates evolution of a starting sequence over generations with mutation rate, population size and selection on a fitness function such as GC content.

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`sequence`, `generations`, `populationSize`, `mutationRate`
- **关键词：**evolution simulation, mutation rate, selection, gc content, lineage, population

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `sequence` | string | 是 | Starting sequence |
| `generations` | number | 是 | Number of generations to simulate |
| `populationSize` | number | 是 | Population size for each generation |
| `mutationRate` | number | 是 | Mutation rate per generation (0-1) |
| `selectionPressure` | number | 否 | Selection pressure (0-1), 0 = no selection, 1 = strong selection |
| `fitnessFunction` | string | 否 | Fitness function: 'gc-content', 'length', 'hydrophobic', or 'custom' 可选值：`gc-content`, `length`, `hydrophobic`, `custom` |
| `targetValue` | number | 否 | Target value for fitness function (e.g., target GC content) |
| `seed` | number | 否 | Random seed for reproducible results (optional) |
| `trackLineages` | boolean | 否 | Track individual lineages through generations 默认 `False` |
| `outputFormat` | string | 否 | Output format: 'summary', 'detailed', or 'fasta' 可选值：`summary`, `detailed`, `fasta` |

??? note "上游说明（英文原文）"

    Simulate evolution of sequences over multiple generations

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'sequence-simulation',
    'evolve_sequence',
    {'sequence': '...',
     'generations': 0.0,
     'populationSize': 0.0,
     'mutationRate': 0.0},
)
```

## `generate_dna_sequence` {#generate_dna_sequence}

Generates random DNA sequences of given length and GC content using random, Markov or codon-biased models, as FASTA or plain text.

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`length`
- **关键词：**random dna, synthetic sequence, gc content, markov, codon bias, fasta

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `length` | number | 是 | Length of the DNA sequence to generate |
| `gcContent` | number | 否 | GC content (0-1), default is 0.5 |
| `count` | number | 否 | Number of sequences to generate, default is 1 |
| `seed` | number | 否 | Random seed for reproducible results (optional) |
| `model` | string | 否 | Generation model: 'random', 'markov', or 'codon-biased' 可选值：`random`, `markov`, `codon-biased` |
| `outputFormat` | string | 否 | Output format: 'fasta' or 'plain' 可选值：`fasta`, `plain` |

??? note "上游说明（英文原文）"

    Generate random DNA sequences with specified parameters

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'sequence-simulation',
    'generate_dna_sequence',
    {'length': 0.0},
)
```

## `mutate_sequence` {#mutate_sequence}

Applies random substitutions, insertions and deletions to a DNA or protein sequence with a transition/transversion bias for DNA.

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`sequence`, `sequenceType`
- **关键词：**mutation, substitution, indel, transition transversion, in silico mutagenesis

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `sequence` | string | 是 | Input sequence to mutate |
| `sequenceType` | string | 是 | Type of sequence: 'dna' or 'protein' 可选值：`dna`, `protein` |
| `substitutionRate` | number | 否 | Substitution mutation rate (0-1), default 0.01 |
| `insertionRate` | number | 否 | Insertion mutation rate (0-1), default 0.001 |
| `deletionRate` | number | 否 | Deletion mutation rate (0-1), default 0.001 |
| `transitionBias` | number | 否 | Transition vs transversion bias for DNA (default 2.0) |
| `iterations` | number | 否 | Number of mutation iterations, default 1 |
| `seed` | number | 否 | Random seed for reproducible results (optional) |
| `outputFormat` | string | 否 | Output format: 'fasta' or 'plain' 可选值：`fasta`, `plain` |

??? note "上游说明（英文原文）"

    Apply mutations to DNA or protein sequences

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'sequence-simulation',
    'mutate_sequence',
    {'sequence': '...', 'sequenceType': 'dna'},
)
```

## `simulate_fastq_file` {#simulate_fastq_file}

Simulates single- or paired-end FASTQ reads from a reference DNA sequence with coverage, insert size, error rate and quality models (NEAT-style).

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`referenceSequence`, `readLength`, `coverage`
- **关键词：**fastq, read simulation, ngs, coverage, paired-end, quality scores, neat

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `referenceSequence` | string | 是 | Reference DNA sequence to generate reads from |
| `readLength` | number | 是 | Length of each sequencing read |
| `coverage` | number | 是 | Target sequencing coverage depth |
| `readType` | string | 否 | Type of reads to generate 可选值：`single-end`, `paired-end` |
| `insertSize` | number | 否 | Mean insert size for paired-end reads (ignored for single-end) |
| `insertSizeStd` | number | 否 | Standard deviation of insert size for paired-end reads |
| `errorRate` | number | 否 | Base calling error rate (0-1) |
| `qualityModel` | string | 否 | Quality score model to use 可选值：`illumina`, `454`, `ion-torrent`, `pacbio` |
| `mutationRate` | number | 否 | Rate of true mutations to introduce (0-1) |
| `seed` | number | 否 | Random seed for reproducible results (optional) |
| `outputFormat` | string | 否 | Output format for the reads 可选值：`fastq`, `json` |

??? note "上游说明（英文原文）"

    Simulate FASTQ sequencing reads with realistic quality scores and error models based on NEAT methodology. Implementation inspired by: Stephens et al. (2016) 'Simulating Next-Generation Sequencing Datasets from Empirical Mutation and Sequencing Models.' PLOS ONE 11(11): e0167047. https://doi.org/10.1371/journal.pone.0167047

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'sequence-simulation',
    'simulate_fastq_file',
    {'referenceSequence': '...', 'readLength': 0.0, 'coverage': 0.0},
)
```

## `simulate_phylogeny` {#simulate_phylogeny}

Simulates a phylogenetic tree and evolves a root sequence along it using JC69, K80, HKY85 or GTR models, outputting FASTA, NEXUS or PHYLIP.

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`rootSequence`
- **关键词：**phylogenetic tree, newick, substitution model, gtr, hky85, molecular clock

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `rootSequence` | string | 是 | Root sequence for phylogenetic simulation |
| `treeStructure` | string | 否 | Newick format tree or 'random' for random tree |
| `numTaxa` | number | 否 | Number of taxa for random tree generation |
| `mutationRate` | number | 否 | Mutation rate per branch length unit |
| `branchLengthVariation` | number | 否 | Variation in branch lengths (0-1) |
| `molecularClock` | boolean | 否 | Use molecular clock (uniform rates) 默认 `True` |
| `substitutionModel` | string | 否 | Substitution model: 'JC69', 'K80', 'HKY85', or 'GTR' 可选值：`JC69`, `K80`, `HKY85`, `GTR` |
| `seed` | number | 否 | Random seed for reproducible results (optional) |
| `outputFormat` | string | 否 | Output format: 'fasta', 'nexus', or 'phylip' 可选值：`fasta`, `nexus`, `phylip` |

??? note "上游说明（英文原文）"

    Simulate phylogenetic tree and sequence evolution

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'sequence-simulation',
    'simulate_phylogeny',
    {'rootSequence': '...'},
)
```
