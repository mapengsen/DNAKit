# MCP Sequence Simulation (`sequence-simulation`)

DNAKit's MCP tool index keeps **5** DNA/gene tools of this server (6 listed, read on 2026-10-16). They run once you install the server; DNAKit does not bundle them. See [calling a tool that was found](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found).

| Item | Value |
| --- | --- |
| Source | [coreymhudson/mcp-sequence-simulation @ `dea87b2`](https://github.com/coreymhudson/mcp-sequence-simulation/tree/dea87b2a45594ff969eb5f05826642f99d42ac29) |
| License | BSD-2-Clause |
| Command | `node dist/server.js` |
| Transports | stdio |
| Installation | Install separately from the pinned source checkout with npm install &amp;&amp; npm run build (Node.js with npm; package.json declares no engines field), then set workdir to the checkout and keep dist/server.js; no environment variables are used. |
| Effects | Pure local computation that returns simulated sequences, trees and FASTQ text in the tool response; the source contains no network calls, file writes or environment access. |

| Tool | Category | Summary |
| --- | --- | --- |
| [`evolve_sequence`](#evolve_sequence) | `sequence` | Simulates evolution of a starting sequence over generations with mutation rate, population size and selection on a fitness function such as GC content. |
| [`generate_dna_sequence`](#generate_dna_sequence) | `sequence` | Generates random DNA sequences of given length and GC content using random, Markov or codon-biased models, as FASTA or plain text. |
| [`mutate_sequence`](#mutate_sequence) | `sequence` | Applies random substitutions, insertions and deletions to a DNA or protein sequence with a transition/transversion bias for DNA. |
| [`simulate_fastq_file`](#simulate_fastq_file) | `sequence` | Simulates single- or paired-end FASTQ reads from a reference DNA sequence with coverage, insert size, error rate and quality models (NEAT-style). |
| [`simulate_phylogeny`](#simulate_phylogeny) | `sequence` | Simulates a phylogenetic tree and evolves a root sequence along it using JC69, K80, HKY85 or GTR models, outputting FASTA, NEXUS or PHYLIP. |

## `evolve_sequence` {#evolve_sequence}

Simulates evolution of a starting sequence over generations with mutation rate, population size and selection on a fitness function such as GC content.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `sequence`, `generations`, `populationSize`, `mutationRate`
- **Keywords:** evolution simulation, mutation rate, selection, gc content, lineage, population

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sequence` | string | yes | Starting sequence |
| `generations` | number | yes | Number of generations to simulate |
| `populationSize` | number | yes | Population size for each generation |
| `mutationRate` | number | yes | Mutation rate per generation (0-1) |
| `selectionPressure` | number | no | Selection pressure (0-1), 0 = no selection, 1 = strong selection |
| `fitnessFunction` | string | no | Fitness function: 'gc-content', 'length', 'hydrophobic', or 'custom' One of `gc-content`, `length`, `hydrophobic`, `custom` |
| `targetValue` | number | no | Target value for fitness function (e.g., target GC content) |
| `seed` | number | no | Random seed for reproducible results (optional) |
| `trackLineages` | boolean | no | Track individual lineages through generations Default `False` |
| `outputFormat` | string | no | Output format: 'summary', 'detailed', or 'fasta' One of `summary`, `detailed`, `fasta` |

??? note "Upstream description"

    Simulate evolution of sequences over multiple generations

Call it through DNAKit once the server is installed:

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

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `length`
- **Keywords:** random dna, synthetic sequence, gc content, markov, codon bias, fasta

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `length` | number | yes | Length of the DNA sequence to generate |
| `gcContent` | number | no | GC content (0-1), default is 0.5 |
| `count` | number | no | Number of sequences to generate, default is 1 |
| `seed` | number | no | Random seed for reproducible results (optional) |
| `model` | string | no | Generation model: 'random', 'markov', or 'codon-biased' One of `random`, `markov`, `codon-biased` |
| `outputFormat` | string | no | Output format: 'fasta' or 'plain' One of `fasta`, `plain` |

??? note "Upstream description"

    Generate random DNA sequences with specified parameters

Call it through DNAKit once the server is installed:

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

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `sequence`, `sequenceType`
- **Keywords:** mutation, substitution, indel, transition transversion, in silico mutagenesis

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sequence` | string | yes | Input sequence to mutate |
| `sequenceType` | string | yes | Type of sequence: 'dna' or 'protein' One of `dna`, `protein` |
| `substitutionRate` | number | no | Substitution mutation rate (0-1), default 0.01 |
| `insertionRate` | number | no | Insertion mutation rate (0-1), default 0.001 |
| `deletionRate` | number | no | Deletion mutation rate (0-1), default 0.001 |
| `transitionBias` | number | no | Transition vs transversion bias for DNA (default 2.0) |
| `iterations` | number | no | Number of mutation iterations, default 1 |
| `seed` | number | no | Random seed for reproducible results (optional) |
| `outputFormat` | string | no | Output format: 'fasta' or 'plain' One of `fasta`, `plain` |

??? note "Upstream description"

    Apply mutations to DNA or protein sequences

Call it through DNAKit once the server is installed:

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

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `referenceSequence`, `readLength`, `coverage`
- **Keywords:** fastq, read simulation, ngs, coverage, paired-end, quality scores, neat

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `referenceSequence` | string | yes | Reference DNA sequence to generate reads from |
| `readLength` | number | yes | Length of each sequencing read |
| `coverage` | number | yes | Target sequencing coverage depth |
| `readType` | string | no | Type of reads to generate One of `single-end`, `paired-end` |
| `insertSize` | number | no | Mean insert size for paired-end reads (ignored for single-end) |
| `insertSizeStd` | number | no | Standard deviation of insert size for paired-end reads |
| `errorRate` | number | no | Base calling error rate (0-1) |
| `qualityModel` | string | no | Quality score model to use One of `illumina`, `454`, `ion-torrent`, `pacbio` |
| `mutationRate` | number | no | Rate of true mutations to introduce (0-1) |
| `seed` | number | no | Random seed for reproducible results (optional) |
| `outputFormat` | string | no | Output format for the reads One of `fastq`, `json` |

??? note "Upstream description"

    Simulate FASTQ sequencing reads with realistic quality scores and error models based on NEAT methodology. Implementation inspired by: Stephens et al. (2016) 'Simulating Next-Generation Sequencing Datasets from Empirical Mutation and Sequencing Models.' PLOS ONE 11(11): e0167047. https://doi.org/10.1371/journal.pone.0167047

Call it through DNAKit once the server is installed:

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

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `rootSequence`
- **Keywords:** phylogenetic tree, newick, substitution model, gtr, hky85, molecular clock

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `rootSequence` | string | yes | Root sequence for phylogenetic simulation |
| `treeStructure` | string | no | Newick format tree or 'random' for random tree |
| `numTaxa` | number | no | Number of taxa for random tree generation |
| `mutationRate` | number | no | Mutation rate per branch length unit |
| `branchLengthVariation` | number | no | Variation in branch lengths (0-1) |
| `molecularClock` | boolean | no | Use molecular clock (uniform rates) Default `True` |
| `substitutionModel` | string | no | Substitution model: 'JC69', 'K80', 'HKY85', or 'GTR' One of `JC69`, `K80`, `HKY85`, `GTR` |
| `seed` | number | no | Random seed for reproducible results (optional) |
| `outputFormat` | string | no | Output format: 'fasta', 'nexus', or 'phylip' One of `fasta`, `nexus`, `phylip` |

??? note "Upstream description"

    Simulate phylogenetic tree and sequence evolution

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'sequence-simulation',
    'simulate_phylogeny',
    {'rootSequence': '...'},
)
```
