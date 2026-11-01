# ToolUniverse: molecular_tools (`tooluniverse`)

DNAKit's MCP tool index keeps **845** DNA/gene tools of this server (2927 listed, read on 2026-10-16). They run once you install the server; DNAKit does not bundle them. See [calling a tool that was found](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found).

This page lists the 23 tools in category `molecular_tools`; other categories: [`sequence`](tooluniverse-sequence.md) (44), [`genome_annotation`](tooluniverse-genome_annotation.md) (253), [`variant`](tooluniverse-variant.md) (155), [`association`](tooluniverse-association.md) (138), [`expression`](tooluniverse-expression.md) (56), [`regulation_epigenome`](tooluniverse-regulation_epigenome.md) (54), [`sequencing_data`](tooluniverse-sequencing_data.md) (91), [`prediction_model`](tooluniverse-prediction_model.md) (17), [`breeding_population`](tooluniverse-breeding_population.md) (14).

| Item | Value |
| --- | --- |
| Source | [mims-harvard/ToolUniverse @ `9959ebe`](https://github.com/mims-harvard/ToolUniverse/tree/9959ebee7638773314362f1746980a15f5b17e5b) |
| License | Apache-2.0 |
| Command | `tooluniverse` |
| Transports | stdio |
| Installation | Install the pinned ToolUniverse distribution and required provider dependencies separately. Use the installed tooluniverse console script; unlike the README uvx --refresh example, this bridge does not refresh or install anything. |
| Effects | Starts the compact ToolUniverse MCP catalog. Selected tools may contact services, load models, execute programs or write files according to their provider definitions. |

| Tool | Category | Summary |
| --- | --- | --- |
| [`Addgene_get_plasmid`](#Addgene_get_plasmid) | `molecular_tools` | Retrieves an Addgene plasmid record with cloning, backbone, insert and resistance details. |
| [`Addgene_search_depositors`](#Addgene_search_depositors) | `molecular_tools` | Searches Addgene depositors by name or institution with plasmid counts and examples. |
| [`Addgene_search_plasmids`](#Addgene_search_plasmids) | `molecular_tools` | Searches the Addgene plasmid catalog by name, gene, organism or vector type. |
| [`CodonUsage_get_optimal_codons`](#CodonUsage_get_optimal_codons) | `molecular_tools` | Returns the most frequently used codon per amino acid for an organism (by NCBI taxid) from the Codon Usage Database. |
| [`CodonUsage_get_table`](#CodonUsage_get_table) | `molecular_tools` | Retrieves the full 64-codon usage table (counts, per-thousand frequency, fraction) for an organism by NCBI taxid. |
| [`DNA_codon_optimize`](#DNA_codon_optimize) | `molecular_tools` | Codon-optimizes an amino acid sequence for human, E. coli, mouse or yeast and reports CAI and GC. |
| [`DNA_crispr_guide_design`](#DNA_crispr_guide_design) | `molecular_tools` | Finds and ranks SpCas9 NGG-PAM guide RNAs on both strands of a DNA sequence by rule-based heuristics. |
| [`DNA_find_restriction_sites`](#DNA_find_restriction_sites) | `molecular_tools` | Finds restriction enzyme recognition sites in a DNA sequence using NEB and Biopython libraries. |
| [`DNA_gibson_design`](#DNA_gibson_design) | `molecular_tools` | Computes Gibson Assembly overlaps between ordered DNA fragments. |
| [`DNA_golden_gate_assemble`](#DNA_golden_gate_assemble) | `molecular_tools` | Computes the product of a Golden Gate reaction from input parts and a Type IIS enzyme. |
| [`DNA_golden_gate_design`](#DNA_golden_gate_design) | `molecular_tools` | Designs Golden Gate parts with unique 4-bp overhangs and BsaI/BbsI sites. |
| [`DNA_primer_design`](#DNA_primer_design) | `molecular_tools` | Designs forward and reverse PCR primers for a target region using SantaLucia thermodynamics. |
| [`DNA_virtual_digest`](#DNA_virtual_digest) | `molecular_tools` | Performs a virtual restriction digest of a linear or circular DNA sequence. |
| [`IDT_analyze_oligo`](#IDT_analyze_oligo) | `molecular_tools` | Calculates oligo Tm, GC content, molecular weight and extinction coefficient via IDT OligoAnalyzer. |
| [`IDT_check_self_dimer`](#IDT_check_self_dimer) | `molecular_tools` | Checks an oligonucleotide for self-dimer structures with delta-G and alignments via IDT OligoAnalyzer. |
| [`NEB_Tm_calculate`](#NEB_Tm_calculate) | `molecular_tools` | Calculates primer Tm and annealing temperature with polymerase-specific salt corrections via the NEB Tm Calculator. |
| [`NEB_Tm_list_polymerases`](#NEB_Tm_list_polymerases) | `molecular_tools` | Lists NEB DNA polymerase products available for Tm calculation. |
| [`REBASE_get_enzyme`](#REBASE_get_enzyme) | `molecular_tools` | Returns REBASE restriction enzyme recognition site, cut position, ends, methylation sensitivity and suppliers. |
| [`REBASE_list_isoschizomers`](#REBASE_list_isoschizomers) | `molecular_tools` | Lists restriction enzymes sharing a recognition specificity with the named enzyme, commercial first. |
| [`REBASE_search_by_site`](#REBASE_search_by_site) | `molecular_tools` | Finds restriction enzymes recognizing a site or cutting a given DNA sequence, with cut counts. |
| [`SynBioHub_get_collections`](#SynBioHub_get_collections) | `molecular_tools` | Lists public SynBioHub collections of genetic parts and designs such as the iGEM Registry. |
| [`SynBioHub_get_part`](#SynBioHub_get_part) | `molecular_tools` | Gets SBOL details of a SynBioHub genetic part including description, DNA sequence, length and functional roles. |
| [`SynBioHub_search_parts`](#SynBioHub_search_parts) | `molecular_tools` | Searches the SynBioHub repository for genetic parts and designs by keyword. |

## `Addgene_get_plasmid` {#Addgene_get_plasmid}

Retrieves an Addgene plasmid record with cloning, backbone, insert and resistance details.

- **Category:** `molecular_tools` (Plasmids, primers, CRISPR guides and cloning)
- **Required arguments:** `plasmid_id`
- **Keywords:** plasmid, vector backbone, insert, cloning, promoter
- **Needs API keys:** `ADDGENE_API_KEY`
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | Operation type One of `get_plasmid` Default `get_plasmid` |
| `plasmid_id` | integer | yes | Addgene plasmid catalog ID (e.g., 39296 for pSpCas9(BB)-2A-Puro, 48138 for pX330) |

??? note "Upstream description"

    Get detailed information about a specific Addgene plasmid by its catalog ID. Returns full record: name, description, cloning details (backbone, promoter, vector types), insert information (gene, species, tags, mutations), bacterial resistance, growth conditions (strain, temperature), article references (DOI, PubMed), and depositor comments. Requires ADDGENE_API_KEY environment variable.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Addgene_get_plasmid', 'arguments': {'plasmid_id': 0}},
)
```

## `Addgene_search_depositors` {#Addgene_search_depositors}

Searches Addgene depositors by name or institution with plasmid counts and examples.

- **Category:** `molecular_tools` (Plasmids, primers, CRISPR guides and cloning)
- **Required arguments:** none
- **Keywords:** plasmid depositor, principal investigator, addgene lab
- **Needs API keys:** `ADDGENE_API_KEY`
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | Operation type One of `search_depositors` Default `search_depositors` |
| `name` | string \| null | no | PI or depositor name to search (e.g., 'Feng Zhang', 'Jennifer Doudna') |
| `institution` | string \| null | no | Institution name to search (e.g., 'MIT', 'Broad Institute') |

??? note "Upstream description"

    Search for Addgene depositors (PIs / principal investigators) by name or institution. Returns unique depositor names with their plasmid counts and example plasmids. At least one of name or institution is required. Requires ADDGENE_API_KEY environment variable. Use for: finding plasmids from a specific lab, browsing depositor contributions.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Addgene_search_depositors', 'arguments': {}},
)
```

## `Addgene_search_plasmids` {#Addgene_search_plasmids}

Searches the Addgene plasmid catalog by name, gene, organism or vector type.

- **Category:** `molecular_tools` (Plasmids, primers, CRISPR guides and cloning)
- **Required arguments:** `query`
- **Keywords:** plasmid search, vector, expression system, insert gene
- **Needs API keys:** `ADDGENE_API_KEY`
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | Operation type One of `search_plasmids` Default `search_plasmids` |
| `query` | string | yes | Plasmid name or keyword to search (e.g., 'pSpCas9', 'GFP', 'lentiviral'). Searches the name field. |
| `organism` | string \| null | no | Filter by organism/species (e.g., 'Human', 'Mouse', 'E. coli'). Maps to the species field. |
| `vector_type` | string \| null | no | Filter by vector type (e.g., 'AAV', 'Lentiviral', 'Retroviral', 'Mammalian Expression'). Maps to vector_types field. |
| `limit` | integer | no | Maximum number of results to return (default 10, max 100) Default `10` |

??? note "Upstream description"

    Search Addgene's plasmid catalog by name, gene, species, or vector type. Addgene is a nonprofit global plasmid repository with 100,000+ plasmids. Returns plasmid name, ID, purpose, depositor, genes, species, vector types, and expression systems. Requires ADDGENE_API_KEY environment variable. Use for: finding CRISPR plasmids, AAV vectors, expression constructs, shRNA libraries.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Addgene_search_plasmids', 'arguments': {'query': '...'}},
)
```

## `CodonUsage_get_optimal_codons` {#CodonUsage_get_optimal_codons}

Returns the most frequently used codon per amino acid for an organism (by NCBI taxid) from the Codon Usage Database.

- **Category:** `molecular_tools` (Plasmids, primers, CRISPR guides and cloning)
- **Required arguments:** `taxid`
- **Keywords:** codon optimization, codon bias, taxid, expression host, preferred codons
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `taxid` | integer \| string | yes | NCBI taxonomy identifier, e.g. 9606 (human) or 83333 (E. coli K-12). |
| `include_stop_codons` | boolean \| null | no | Include the stop codon group (reported as 'End'). Default false. |

??? note "Upstream description"

    Return the most-used codon for each amino acid in an organism, which is the reference table needed to codon-optimize a sequence for expression in that host. Example: taxid=83333 gives the preferred codons for E. coli K-12, which differ from human at ten amino acids. Use with DNA_codon_optimize when the target host is outside its four built-in species. Optimizing every position to the single best codon maximizes adaptation but can deplete tRNA pools, so the returned fraction supports proportional sampling instead.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'CodonUsage_get_optimal_codons',
     'arguments': {'taxid': 0}},
)
```

## `CodonUsage_get_table` {#CodonUsage_get_table}

Retrieves the full 64-codon usage table (counts, per-thousand frequency, fraction) for an organism by NCBI taxid.

- **Category:** `molecular_tools` (Plasmids, primers, CRISPR guides and cloning)
- **Required arguments:** `taxid`
- **Keywords:** codon usage, codon frequency, taxid, codon bias, synonymous codons
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `taxid` | integer \| string | yes | NCBI taxonomy identifier, e.g. 9606 (human), 83333 (E. coli K-12), 4932 (yeast), 7227 (fruit fly), 3702 (Arabidopsis). Tables are often held at strain level, so 562 for E. coli has no table but 83333 does. |
| `amino_acid` | string \| null | no | Restrict to one amino acid using its three-letter code, e.g. 'Leu', 'Ala', or 'End' for stop codons. |

??? note "Upstream description"

    Retrieve the codon usage table for an organism from the Codon Usage Database, indexed by NCBI taxonomy identifier. Returns all 64 codons with their raw count, frequency per thousand codons, and the fraction of their amino acid's usage they account for. DNA_codon_optimize carries hardcoded tables for only four species; this supplies the reference data for any sequenced organism. Example: taxid=9606 for human, or taxid=9606 with amino_acid='Leu' for just the leucine codons.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'CodonUsage_get_table', 'arguments': {'taxid': 0}},
)
```

## `DNA_codon_optimize` {#DNA_codon_optimize}

Codon-optimizes an amino acid sequence for human, E. coli, mouse or yeast and reports CAI and GC.

- **Category:** `molecular_tools` (Plasmids, primers, CRISPR guides and cloning)
- **Required arguments:** `sequence`
- **Keywords:** codon optimization, cai, back-translation, expression host
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | Operation (optional; defaults to 'codon_optimize' for this tool). One of `codon_optimize` Default `codon_optimize` |
| `sequence` | string | yes | Amino acid sequence in single-letter code (e.g., 'MEPVDDLPL'). Stop codon (*) is optional at the end. |
| `species` | string \| null | no | Target expression organism. Options: 'human' (default), 'ecoli', 'mouse', 'yeast'. One of `human`, `ecoli`, `mouse`, `yeast`, `None` Default `human` |

??? note "Upstream description"

    Codon-optimize an amino acid sequence for expression in a target species (human, ecoli, mouse, or yeast) using species-specific codon frequency tables. Returns the optimized DNA sequence, GC content, Codon Adaptation Index (CAI), and sequence length. Use for: synthetic gene design, heterologous protein expression, improving codon usage for a target host organism.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'DNA_codon_optimize', 'arguments': {'sequence': '...'}},
)
```

## `DNA_crispr_guide_design` {#DNA_crispr_guide_design}

Finds and ranks SpCas9 NGG-PAM guide RNAs on both strands of a DNA sequence by rule-based heuristics.

- **Category:** `molecular_tools` (Plasmids, primers, CRISPR guides and cloning)
- **Required arguments:** `sequence`
- **Keywords:** crispr, sgrna, spcas9, pam, guide design
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | Operation (optional; defaults to 'crispr_guide_design' for this tool). One of `crispr_guide_design` Default `crispr_guide_design` |
| `sequence` | string | yes | DNA target sequence (A, T, G, C only; no N), at least 23 bp. |
| `gc_min` | number \| null | no | Minimum acceptable protospacer GC%, 0-100 (default 30). Default `30.0` |
| `gc_max` | number \| null | no | Maximum acceptable protospacer GC%, 0-100 (default 70). Default `70.0` |
| `limit` | integer \| null | no | Max candidate guides to return, 1-100 (default 10). Default `10` |

??? note "Upstream description"

    Find candidate SpCas9 (NGG PAM) guide RNAs in a target DNA sequence, on both strands. Ranks 20 bp protospacers by a transparent rule-based heuristic (GC content, PAM-proximal G, Pol III terminator and homopolymer avoidance) -- not a trained on-target efficiency model, and no off-target search against any genome is performed. Use for a first-pass candidate list; verify with CRISPOR or the Broad GPP sgRNA designer before ordering.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'DNA_crispr_guide_design',
     'arguments': {'sequence': '...'}},
)
```

## `DNA_find_restriction_sites` {#DNA_find_restriction_sites}

Finds restriction enzyme recognition sites in a DNA sequence using NEB and Biopython libraries.

- **Category:** `molecular_tools` (Plasmids, primers, CRISPR guides and cloning)
- **Required arguments:** `sequence`
- **Keywords:** restriction enzyme, recognition site, neb, type iis
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | Operation (optional; defaults to 'find_restriction_sites' for this tool). One of `find_restriction_sites` Default `find_restriction_sites` |
| `sequence` | string | yes | DNA sequence (A, T, G, C, N only). Case insensitive. Spaces and newlines are ignored. |
| `enzymes` | array \| null | no | Optional list of specific enzyme names to check (e.g., \['EcoRI', 'BamHI'\]). If null, checks all 24 NEB enzymes. |

??? note "Upstream description"

    Find restriction enzyme recognition sites in a DNA sequence using the NEB enzyme library (25 curated common enzymes plus Biopython's ~600-enzyme library, including the Type IIS enzymes used for Golden Gate/MoClo: BsaI, BbsI, Esp3I/BsmBI, SapI). Returns positions, recognition sequences, and cut counts for each enzyme that has sites. Use for: cloning strategy design, vector linearization planning, molecular biology workflows.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'DNA_find_restriction_sites',
     'arguments': {'sequence': '...'}},
)
```

## `DNA_gibson_design` {#DNA_gibson_design}

Computes Gibson Assembly overlaps between ordered DNA fragments.

- **Category:** `molecular_tools` (Plasmids, primers, CRISPR guides and cloning)
- **Required arguments:** `fragments`
- **Keywords:** gibson assembly, overlap, cloning, dna fragments
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | Operation (optional; defaults to 'gibson_design' for this tool). One of `gibson_design` Default `gibson_design` |
| `fragments` | array of string | yes | List of DNA fragment sequences (at least 2). Each must be longer than overlap_length. |
| `overlap_length` | integer \| null | no | Overlap length in bp for Gibson Assembly (default: 20 bp). Must be at least 1. Default `20` |

??? note "Upstream description"

    Design Gibson Assembly overlaps for seamless DNA fragment assembly. For each fragment, computes the 5' overhang (from the previous fragment's 3' end) and 3' overhang (from the next fragment's 5' end) based on the specified overlap length (default 20 bp). Returns the modified fragment sequences ready for Gibson Assembly cloning. Use for: seamless cloning, multi-fragment assembly, vector backbone insertion.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'DNA_gibson_design', 'arguments': {'fragments': []}},
)
```

## `DNA_golden_gate_assemble` {#DNA_golden_gate_assemble}

Computes the product of a Golden Gate reaction from input parts and a Type IIS enzyme.

- **Category:** `molecular_tools` (Plasmids, primers, CRISPR guides and cloning)
- **Required arguments:** `fragments`
- **Keywords:** golden gate, type iis, bsai, assembly product, plasmid
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | Operation (optional; defaults to 'golden_gate_assemble' for this tool). One of `golden_gate_assemble` Default `golden_gate_assemble` |
| `fragments` | array of string | yes | The DNA sequences combined in the reaction (plasmids or linear parts), A/T/G/C only. At least 2. |
| `enzyme` | string \| null | no | Type IIS enzyme used in the reaction, e.g. 'BsaI' (default), 'BbsI', 'Esp3I'/'BsmBI', 'SapI'. Default `BsaI` |
| `circular` | boolean \| null | no | Whether the inputs are circular plasmids (default true). Set false for linear parts. Default `True` |
| `labels` | array \| null | no | Optional names for the inputs (e.g. plasmid names), echoed in the assembly order. |

??? note "Upstream description"

    Work out the product of a Golden Gate reaction that has already been performed: given the plasmids/parts that were combined and the Type IIS enzyme used (BsaI, BbsI, Esp3I/BsmBI, SapI), digest each input, discard fragments that keep a recognition site (they are re-cut in the reaction), and ligate the rest by matching 4-bp overhangs into the circular product. Returns the assembled sequence, its length, and the fragment order with the overhang at each junction. Use to answer 'I combined these plasmids with &lt;enzyme&gt; — what does the resulting plasmid contain/express?'. This is the reverse of DNA_golden_gate_design, which designs parts rather than analysing a finished reaction.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'DNA_golden_gate_assemble',
     'arguments': {'fragments': []}},
)
```

## `DNA_golden_gate_design` {#DNA_golden_gate_design}

Designs Golden Gate parts with unique 4-bp overhangs and BsaI/BbsI sites.

- **Category:** `molecular_tools` (Plasmids, primers, CRISPR guides and cloning)
- **Required arguments:** `parts`
- **Keywords:** golden gate, overhang, bsai, bbsi, cloning
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | Operation (optional; defaults to 'golden_gate_design' for this tool). One of `golden_gate_design` Default `golden_gate_design` |
| `parts` | array of string | yes | List of DNA part sequences to assemble (at least 2). A, T, G, C only. |
| `enzyme` | string \| null | no | Type IIS restriction enzyme to use. Options: 'BsaI' (default) or 'BbsI'. One of `BsaI`, `BbsI`, `None` Default `BsaI` |

??? note "Upstream description"

    Design Golden Gate Assembly parts with Type IIS restriction enzyme sites (BsaI or BbsI). Assigns unique non-palindromic 4-bp overhangs to each junction and builds the full synthetic part sequences including enzyme recognition sites. Use for: scarless multi-part assembly, modular cloning (MoClo), Golden Gate cloning workflows.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'DNA_golden_gate_design', 'arguments': {'parts': []}},
)
```

## `DNA_primer_design` {#DNA_primer_design}

Designs forward and reverse PCR primers for a target region using SantaLucia thermodynamics.

- **Category:** `molecular_tools` (Plasmids, primers, CRISPR guides and cloning)
- **Required arguments:** `sequence`
- **Keywords:** pcr primer, tm, nearest-neighbor, amplicon
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | Operation (optional; defaults to 'primer_design' for this tool). One of `primer_design` Default `primer_design` |
| `sequence` | string | yes | DNA template sequence (A, T, G, C, N only). Must be at least 200 bp for good primer design. |
| `target_start` | integer \| null | no | 0-based start position of the target region to amplify. Defaults to 0. |
| `target_end` | integer \| null | no | 0-based end position of the target region to amplify. Defaults to sequence length. |
| `tm_target` | number \| null | no | Target melting temperature in Celsius (default: 60.0°C). Default `60.0` |
| `product_size_min` | integer \| null | no | Minimum acceptable PCR product size in bp (default: 100). Default `100` |
| `product_size_max` | integer \| null | no | Maximum acceptable PCR product size in bp (default: 1000). Default `1000` |

??? note "Upstream description"

    Design forward and reverse PCR primers for a target DNA region using the SantaLucia 1998 nearest-neighbor thermodynamic model. Selects primers with Tm closest to the target (default 60°C), GC content 40-60%, length 18-25 bp, and no 3'-end repeats. Use for: PCR amplification, genotyping, diagnostic assay design, amplicon sequencing.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'DNA_primer_design', 'arguments': {'sequence': '...'}},
)
```

## `DNA_virtual_digest` {#DNA_virtual_digest}

Performs a virtual restriction digest of a linear or circular DNA sequence.

- **Category:** `molecular_tools` (Plasmids, primers, CRISPR guides and cloning)
- **Required arguments:** `sequence`
- **Keywords:** restriction digest, fragments, neb, circular plasmid
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | Operation (optional; defaults to 'virtual_digest' for this tool). One of `virtual_digest` Default `virtual_digest` |
| `sequence` | string | yes | DNA sequence (A, T, G, C, N only). Case insensitive. |
| `enzymes` | array \| null | no | List of enzyme names to use for digestion (e.g., \['EcoRI', 'BamHI'\]). If null, uses all 24 NEB enzymes. |
| `circular` | boolean \| null | no | Treat sequence as circular DNA (e.g., plasmid). Default: false (linear). Default `False` |

??? note "Upstream description"

    Perform a virtual restriction enzyme digest on a DNA sequence using the NEB enzyme library (25 curated enzymes plus Biopython's ~600-enzyme library, including Type IIS enzymes such as BsaI, BbsI, Esp3I/BsmBI and SapI, which cut outside their recognition site). Returns all resulting fragments with their sequences, lengths, and positions, plus a list of cut sites. Supports both linear and circular DNA. Use for: cloning planning, gel electrophoresis prediction, verifying construct integrity.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'DNA_virtual_digest', 'arguments': {'sequence': '...'}},
)
```

## `IDT_analyze_oligo` {#IDT_analyze_oligo}

Calculates oligo Tm, GC content, molecular weight and extinction coefficient via IDT OligoAnalyzer.

- **Category:** `molecular_tools` (Plasmids, primers, CRISPR guides and cloning)
- **Required arguments:** `sequence`
- **Keywords:** oligo, primer, melting temperature, tm, gc content
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sequence` | string | yes | Oligonucleotide sequence (5' to 3'). Must contain only A, T, C, G (DNA) or A, U, C, G (RNA) bases. Length: 5-200 bases. Example: 'AGCGGATAACAATTTCACACAGGA' |
| `na_concentration_mm` | number \| null | no | Sodium (Na+) concentration in millimolar (mM). Affects Tm calculation. Typical range: 0-1000 mM. Default: 50 mM. Standard PCR buffer is ~50 mM. |
| `mg_concentration_mm` | number \| null | no | Magnesium (Mg2+) concentration in millimolar (mM). Stabilizes duplexes, raising Tm. Typical range: 0-20 mM. Default: 0 mM. Standard PCR uses 1.5-3.0 mM. |
| `dntps_concentration_mm` | number \| null | no | dNTPs concentration in millimolar (mM). dNTPs chelate Mg2+, reducing effective Mg2+. Typical range: 0-10 mM. Default: 0 mM. Standard PCR uses 0.2-0.8 mM. |
| `oligo_concentration_um` | number \| null | no | Oligo concentration in micromolar (uM). Affects Tm slightly. Default: 0.25 uM. |
| `oligo_type` | string \| null | no | Type of oligonucleotide: 'DNA' or 'RNA'. Default: 'DNA'. |

??? note "Upstream description"

    Analyze an oligonucleotide using the IDT OligoAnalyzer API. Calculates melting temperature (Tm), GC content, molecular weight, molar extinction coefficient, and nmol/OD260 and ug/OD260 yields. Tm is calculated using nearest-neighbor thermodynamic parameters (SantaLucia 1998) under user-specified salt conditions (Na+, Mg2+, dNTPs). Supports DNA and RNA oligos up to 200 bases. Use this BEFORE ordering primers/probes to verify Tm, check GC balance, and estimate synthesis yields. Complements NEB_Tm_calculate (which provides polymerase-specific Tm with salt corrections) and DNA_calculate_gc_content (local, no Tm). For self-dimer/homodimer assessment, use IDT_check_self_dimer.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'IDT_analyze_oligo', 'arguments': {'sequence': '...'}},
)
```

## `IDT_check_self_dimer` {#IDT_check_self_dimer}

Checks an oligonucleotide for self-dimer structures with delta-G and alignments via IDT OligoAnalyzer.

- **Category:** `molecular_tools` (Plasmids, primers, CRISPR guides and cloning)
- **Required arguments:** `sequence`
- **Keywords:** primer dimer, self-complementarity, delta g, oligo
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sequence` | string | yes | Oligonucleotide sequence (5' to 3'). Must contain only A, T, C, G, U bases. Length: 5-200 bases. Example: 'AGCGGATAACAATTTCACACAGGA' |
| `na_concentration_mm` | number \| null | no | Sodium (Na+) concentration in millimolar (mM). Affects dimer stability. Default: 50 mM. |
| `mg_concentration_mm` | number \| null | no | Magnesium (Mg2+) concentration in millimolar (mM). Default: 0 mM. |
| `temperature_celsius` | number \| null | no | Temperature in Celsius for dimer stability calculation. Use your expected annealing temperature. Default: 25. |

??? note "Upstream description"

    Check an oligonucleotide for self-complementarity and homodimer (self-dimer) formation using the IDT OligoAnalyzer API. Returns predicted self-dimer structures with their thermodynamic stability (delta-G in kcal/mol), number of base pairs, and alignment diagrams. A more negative delta-G indicates a more stable (problematic) dimer. Risk assessment: delta-G &gt; -3.0 kcal/mol is low risk, -3.0 to -6.0 is moderate risk, &lt; -6.0 is high risk. Use this tool AFTER designing primers (DNA_primer_design or manual) and BEFORE ordering, to verify primers won't form stable homodimers that reduce PCR efficiency. For basic oligo properties (Tm, GC%, MW), use IDT_analyze_oligo instead.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'IDT_check_self_dimer', 'arguments': {'sequence': '...'}},
)
```

## `NEB_Tm_calculate` {#NEB_Tm_calculate}

Calculates primer Tm and annealing temperature with polymerase-specific salt corrections via the NEB Tm Calculator.

- **Category:** `molecular_tools` (Plasmids, primers, CRISPR guides and cloning)
- **Required arguments:** `primer_sequence`
- **Keywords:** primer tm, annealing temperature, pcr, nearest-neighbor
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `primer_sequence` | string | yes | Forward primer DNA sequence (5' to 3'). Must contain only A, T, C, G bases and be at least 10 nt long. Example: 'AGCGGATAACAATTTCACACAGGA' |
| `primer_sequence_2` | string \| null | no | Optional reverse primer DNA sequence (5' to 3'). When provided, the API also calculates annealing temperature (Ta). Example: 'CCTCACTAAAGGGAACAAAAGCTG' |
| `polymerase` | string \| null | no | NEB polymerase product code. Determines buffer composition for Tm calculation. Common codes: 'q5-0' (Q5 HF, default), 'phusion-0' (Phusion HF Buffer), 'taq-0' (Taq Standard Buffer), 'onetaq-0' (OneTaq Standard Buffer), 'custom' (specify monovalent_salt_mm). Use NEB_Tm_list_polymerases for full list. Default: 'q5-0'. |
| `primer_concentration` | integer \| null | no | Primer concentration in nanomolar (nM). Affects Tm calculation. Common values: 200-500 nM. Default: 250. |
| `monovalent_salt_mm` | integer \| null | no | Monovalent salt concentration in mM. Only used when polymerase is 'custom'. Typical range: 20-100 mM. Default: 50. |

??? note "Upstream description"

    Calculate primer melting temperature (Tm) and annealing temperature (Ta) using the NEB Tm Calculator API. Uses thermodynamic nearest-neighbor calculations (SantaLucia 1998) with polymerase-specific salt corrections for accurate PCR temperature predictions. Supports all NEB DNA polymerases including Q5, Phusion, Taq, OneTaq, and custom buffer conditions. When two primers are provided, also calculates the recommended annealing temperature (Ta). Complements the local DNA_primer_design tool by providing manufacturer-validated Tm values, though results can diverge by 10+°C between this tool, DNA_primer_design, and IDT_analyze_oligo for the same primer since each uses different nearest-neighbor parameter sets and salt-correction models -- treat cross-tool Tm comparisons as approximate, not interchangeable. Common polymerase codes: 'q5-0' (Q5 High-Fidelity), 'phusion-0' (Phusion HF Buffer), 'taq-0' (Taq Standard Buffer), 'onetaq-0' (OneTaq Standard Buffer). Use NEB_Tm_list_polymerases to see all available product codes.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'NEB_Tm_calculate',
     'arguments': {'primer_sequence': '...'}},
)
```

## `NEB_Tm_list_polymerases` {#NEB_Tm_list_polymerases}

Lists NEB DNA polymerase products available for Tm calculation.

- **Category:** `molecular_tools` (Plasmids, primers, CRISPR guides and cloning)
- **Required arguments:** none
- **Keywords:** polymerase, q5, phusion, taq
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `filter` | string \| null | no | Optional keyword to filter products by name, product code, or catalog number. Case-insensitive. Examples: 'Q5', 'Phusion', 'hot start', 'master mix', 'M0530'. |

??? note "Upstream description"

    List all NEB DNA polymerase products available for Tm calculation. Returns product codes (prodcode), product names, and NEB catalog numbers for all 54 supported polymerases/kits. Includes Q5, Phusion, Taq, OneTaq, LongAmp, Vent, Deep Vent, and specialty polymerases. Use the returned prodcode values with NEB_Tm_calculate to get polymerase-specific melting temperatures. Can optionally filter by keyword (e.g., 'Q5', 'Phusion', 'hot start').

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'NEB_Tm_list_polymerases', 'arguments': {}},
)
```

## `REBASE_get_enzyme` {#REBASE_get_enzyme}

Returns REBASE restriction enzyme recognition site, cut position, ends, methylation sensitivity and suppliers.

- **Category:** `molecular_tools` (Plasmids, primers, CRISPR guides and cloning)
- **Required arguments:** `name`
- **Keywords:** recognition site, cut site, sticky ends, neb
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `name` | string | yes | Enzyme name, e.g. 'EcoRI', 'BamHI', 'SmaI'. Matched case-insensitively. |

??? note "Upstream description"

    Look up a restriction enzyme in REBASE, the reference catalogue of ~6,000 enzymes from New England Biolabs. Returns the recognition sequence with its cut point, whether the cut leaves blunt or sticky ends, the source organism, methylation sensitivity, commercial availability, and isoschizomers. Example: name='EcoRI' returns G^AATTC. The DNA_find_restriction_sites and DNA_virtual_digest calculators work from 25 curated enzymes; use this to look up anything outside that set.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'REBASE_get_enzyme', 'arguments': {'name': '...'}},
)
```

## `REBASE_list_isoschizomers` {#REBASE_list_isoschizomers}

Lists restriction enzymes sharing a recognition specificity with the named enzyme, commercial first.

- **Category:** `molecular_tools` (Plasmids, primers, CRISPR guides and cloning)
- **Required arguments:** `name`
- **Keywords:** isoschizomer, recognition site, commercial enzyme
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `name` | string | yes | Enzyme name to find isoschizomers of, e.g. 'EcoRI'. |

??? note "Upstream description"

    List every restriction enzyme sharing a prototype with the named enzyme, i.e. all enzymes with the same recognition specificity. Commercially available enzymes are listed first, so this answers 'which enzyme can I actually buy that cuts like this one'. Example: name='EcoRI'.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'REBASE_list_isoschizomers', 'arguments': {'name': '...'}},
)
```

## `REBASE_search_by_site` {#REBASE_search_by_site}

Finds restriction enzymes recognizing a site or cutting a given DNA sequence, with cut counts.

- **Category:** `molecular_tools` (Plasmids, primers, CRISPR guides and cloning)
- **Required arguments:** none
- **Keywords:** recognition site, digest, iupac, cloning
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `site` | string \| null | no | Exact recognition sequence to match, e.g. 'GAATTC'. Mutually exclusive with sequence. |
| `sequence` | string \| null | no | DNA sequence to find cutters for, e.g. a plasmid region. Mutually exclusive with site. |
| `only_commercial` | boolean \| null | no | If true, return only commercially available enzymes. |
| `limit` | integer \| null | no | Maximum enzymes to return (default 25, max 200). |

??? note "Upstream description"

    Find restriction enzymes by what they recognize. Give site to find every enzyme reading an exact recognition sequence (site='GAATTC'), or give sequence to find every enzyme that would cut a DNA sequence, with the number of cut sites each produces. IUPAC ambiguity codes in enzyme sites are expanded, so degenerate recognition sequences match correctly. Set only_commercial=true to restrict to enzymes you can buy.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'REBASE_search_by_site', 'arguments': {}},
)
```

## `SynBioHub_get_collections` {#SynBioHub_get_collections}

Lists public SynBioHub collections of genetic parts and designs such as the iGEM Registry.

- **Category:** `molecular_tools` (Plasmids, primers, CRISPR guides and cloning)
- **Required arguments:** none
- **Keywords:** synthetic biology, biobricks, igem, genetic parts, sbol
- **Read-only:** yes (declared by the server)

??? note "Upstream description"

    List all public collections available on SynBioHub. Collections are curated sets of genetic parts and designs, including the iGEM Parts Registry, CIDAR Lab parts, Free Genes feature libraries, and organism-specific collections. Returns collection names, descriptions, URIs, and member counts. Use this to discover available part libraries before searching for specific parts. Note: SynBioHub anonymous REST access was disabled upstream. Every /public/... read returns 401 'Login required'. The tool returns a structured error; users with a SynBioHub account can supply the auth token via SYNBIOHUB_API_TOKEN.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'SynBioHub_get_collections', 'arguments': {}},
)
```

## `SynBioHub_get_part` {#SynBioHub_get_part}

Gets SBOL details of a SynBioHub genetic part including description, DNA sequence, length and functional roles.

- **Category:** `molecular_tools` (Plasmids, primers, CRISPR guides and cloning)
- **Required arguments:** none
- **Keywords:** synthetic biology, sbol, dna sequence, promoter, biobrick
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `display_id` | string \| null | no | BioBrick/part display ID from the iGEM collection. Examples: 'BBa_E0040' (GFP), 'BBa_R0010' (lac promoter), 'BBa_B0034' (RBS). If not from iGEM, use part_uri instead. |
| `part_uri` | string \| null | no | Full SynBioHub URI for the part (from search results). Example: 'https://synbiohub.org/public/igem/BBa_E0040/1'. Use this for parts not in iGEM collection. |

??? note "Upstream description"

    Get detailed SBOL (Synthetic Biology Open Language) information for a specific genetic part from SynBioHub. Returns the part's title, description, DNA sequence, sequence length, functional roles, creation date, and provenance. Can look up by full URI or by iGEM BioBrick display ID (e.g., BBa_E0040 for GFP). Example: get details for BBa_E0040 (GFP coding sequence, 720bp) including its full DNA sequence. Note: SynBioHub anonymous REST access was disabled upstream. Every /public/... read returns 401 'Login required'. The tool returns a structured error; users with a SynBioHub account can supply the auth token via SYNBIOHUB_API_TOKEN.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'SynBioHub_get_part', 'arguments': {}},
)
```

## `SynBioHub_search_parts` {#SynBioHub_search_parts}

Searches the SynBioHub repository for genetic parts and designs by keyword.

- **Category:** `molecular_tools` (Plasmids, primers, CRISPR guides and cloning)
- **Required arguments:** `query`
- **Keywords:** synthetic biology, igem registry, promoter, terminator, biobricks
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string | yes | Search query for finding genetic parts. Can include gene names, part functions, BioBrick IDs, or biological descriptions. Examples: 'GFP', 'lac promoter', 'BBa_E0040', 'terminator', 'riboswitch', 'CRISPR', 'T7 promoter'. |
| `offset` | integer \| null | no | Offset for pagination (default: 0). |
| `limit` | integer \| null | no | Maximum number of results to return (default: 10, max: 50). |

??? note "Upstream description"

    Search the SynBioHub synthetic biology parts repository for genetic parts and designs by keyword. SynBioHub hosts the iGEM Registry of Standard Biological Parts (20,000+ BioBricks) and other public SBOL collections. Parts include promoters, coding sequences (CDS/genes), terminators, ribosome binding sites (RBS), reporters (GFP, RFP, LacZ), regulatory elements, and composite genetic devices. Returns part IDs, names, descriptions, SBOL types, and functional roles. Example: search 'GFP' finds green fluorescent protein parts; search 'lac promoter' finds lac operon regulatory parts. Note: SynBioHub anonymous REST access was disabled upstream. Every /public/... read returns 401 'Login required'. The tool returns a structured error; users with a SynBioHub account can supply the auth token via SYNBIOHUB_API_TOKEN.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'SynBioHub_search_parts', 'arguments': {'query': '...'}},
)
```
