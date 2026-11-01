# ToolUniverse: molecular_tools (`tooluniverse`)

DNAKit 的 MCP 工具索引收录了这个服务器的 **845** 个 DNA/基因相关工具（服务器共列出 2927 个，2026-10-16 读取）。工具在用户自己安装服务器后运行，DNAKit 本身不包含它们。调用方式见 [调用索引中的工具](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found)。下面的工具说明是上游英文原文。

本页列出类别 `molecular_tools` 的 23 个工具，其他类别：[`sequence`](tooluniverse-sequence.md) (44), [`genome_annotation`](tooluniverse-genome_annotation.md) (253), [`variant`](tooluniverse-variant.md) (155), [`association`](tooluniverse-association.md) (138), [`expression`](tooluniverse-expression.md) (56), [`regulation_epigenome`](tooluniverse-regulation_epigenome.md) (54), [`sequencing_data`](tooluniverse-sequencing_data.md) (91), [`prediction_model`](tooluniverse-prediction_model.md) (17), [`breeding_population`](tooluniverse-breeding_population.md) (14)。

| 项目 | 内容 |
| --- | --- |
| 源码 | [mims-harvard/ToolUniverse @ `9959ebe`](https://github.com/mims-harvard/ToolUniverse/tree/9959ebee7638773314362f1746980a15f5b17e5b) |
| 许可证 | Apache-2.0 |
| 启动命令 | `tooluniverse` |
| 传输方式 | stdio |
| 安装 | Install the pinned ToolUniverse distribution and required provider dependencies separately. Use the installed tooluniverse console script; unlike the README uvx --refresh example, this bridge does not refresh or install anything. |
| 副作用 | Starts the compact ToolUniverse MCP catalog. Selected tools may contact services, load models, execute programs or write files according to their provider definitions. |

| 工具 | 类别 | 简介 |
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

- **类别：**`molecular_tools`（Plasmids, primers, CRISPR guides and cloning）
- **必填参数：**`plasmid_id`
- **关键词：**plasmid, vector backbone, insert, cloning, promoter
- **需要 API 密钥：**`ADDGENE_API_KEY`
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `operation` | string | 否 | Operation type 可选值：`get_plasmid` 默认 `get_plasmid` |
| `plasmid_id` | integer | 是 | Addgene plasmid catalog ID (e.g., 39296 for pSpCas9(BB)-2A-Puro, 48138 for pX330) |

??? note "上游说明（英文原文）"

    Get detailed information about a specific Addgene plasmid by its catalog ID. Returns full record: name, description, cloning details (backbone, promoter, vector types), insert information (gene, species, tags, mutations), bacterial resistance, growth conditions (strain, temperature), article references (DOI, PubMed), and depositor comments. Requires ADDGENE_API_KEY environment variable.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`molecular_tools`（Plasmids, primers, CRISPR guides and cloning）
- **必填参数：**无
- **关键词：**plasmid depositor, principal investigator, addgene lab
- **需要 API 密钥：**`ADDGENE_API_KEY`
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `operation` | string | 否 | Operation type 可选值：`search_depositors` 默认 `search_depositors` |
| `name` | string \| null | 否 | PI or depositor name to search (e.g., 'Feng Zhang', 'Jennifer Doudna') |
| `institution` | string \| null | 否 | Institution name to search (e.g., 'MIT', 'Broad Institute') |

??? note "上游说明（英文原文）"

    Search for Addgene depositors (PIs / principal investigators) by name or institution. Returns unique depositor names with their plasmid counts and example plasmids. At least one of name or institution is required. Requires ADDGENE_API_KEY environment variable. Use for: finding plasmids from a specific lab, browsing depositor contributions.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`molecular_tools`（Plasmids, primers, CRISPR guides and cloning）
- **必填参数：**`query`
- **关键词：**plasmid search, vector, expression system, insert gene
- **需要 API 密钥：**`ADDGENE_API_KEY`
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `operation` | string | 否 | Operation type 可选值：`search_plasmids` 默认 `search_plasmids` |
| `query` | string | 是 | Plasmid name or keyword to search (e.g., 'pSpCas9', 'GFP', 'lentiviral'). Searches the name field. |
| `organism` | string \| null | 否 | Filter by organism/species (e.g., 'Human', 'Mouse', 'E. coli'). Maps to the species field. |
| `vector_type` | string \| null | 否 | Filter by vector type (e.g., 'AAV', 'Lentiviral', 'Retroviral', 'Mammalian Expression'). Maps to vector_types field. |
| `limit` | integer | 否 | Maximum number of results to return (default 10, max 100) 默认 `10` |

??? note "上游说明（英文原文）"

    Search Addgene's plasmid catalog by name, gene, species, or vector type. Addgene is a nonprofit global plasmid repository with 100,000+ plasmids. Returns plasmid name, ID, purpose, depositor, genes, species, vector types, and expression systems. Requires ADDGENE_API_KEY environment variable. Use for: finding CRISPR plasmids, AAV vectors, expression constructs, shRNA libraries.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`molecular_tools`（Plasmids, primers, CRISPR guides and cloning）
- **必填参数：**`taxid`
- **关键词：**codon optimization, codon bias, taxid, expression host, preferred codons
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `taxid` | integer \| string | 是 | NCBI taxonomy identifier, e.g. 9606 (human) or 83333 (E. coli K-12). |
| `include_stop_codons` | boolean \| null | 否 | Include the stop codon group (reported as 'End'). Default false. |

??? note "上游说明（英文原文）"

    Return the most-used codon for each amino acid in an organism, which is the reference table needed to codon-optimize a sequence for expression in that host. Example: taxid=83333 gives the preferred codons for E. coli K-12, which differ from human at ten amino acids. Use with DNA_codon_optimize when the target host is outside its four built-in species. Optimizing every position to the single best codon maximizes adaptation but can deplete tRNA pools, so the returned fraction supports proportional sampling instead.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`molecular_tools`（Plasmids, primers, CRISPR guides and cloning）
- **必填参数：**`taxid`
- **关键词：**codon usage, codon frequency, taxid, codon bias, synonymous codons
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `taxid` | integer \| string | 是 | NCBI taxonomy identifier, e.g. 9606 (human), 83333 (E. coli K-12), 4932 (yeast), 7227 (fruit fly), 3702 (Arabidopsis). Tables are often held at strain level, so 562 for E. coli has no table but 83333 does. |
| `amino_acid` | string \| null | 否 | Restrict to one amino acid using its three-letter code, e.g. 'Leu', 'Ala', or 'End' for stop codons. |

??? note "上游说明（英文原文）"

    Retrieve the codon usage table for an organism from the Codon Usage Database, indexed by NCBI taxonomy identifier. Returns all 64 codons with their raw count, frequency per thousand codons, and the fraction of their amino acid's usage they account for. DNA_codon_optimize carries hardcoded tables for only four species; this supplies the reference data for any sequenced organism. Example: taxid=9606 for human, or taxid=9606 with amino_acid='Leu' for just the leucine codons.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`molecular_tools`（Plasmids, primers, CRISPR guides and cloning）
- **必填参数：**`sequence`
- **关键词：**codon optimization, cai, back-translation, expression host
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `operation` | string | 否 | Operation (optional; defaults to 'codon_optimize' for this tool). 可选值：`codon_optimize` 默认 `codon_optimize` |
| `sequence` | string | 是 | Amino acid sequence in single-letter code (e.g., 'MEPVDDLPL'). Stop codon (*) is optional at the end. |
| `species` | string \| null | 否 | Target expression organism. Options: 'human' (default), 'ecoli', 'mouse', 'yeast'. 可选值：`human`, `ecoli`, `mouse`, `yeast`, `None` 默认 `human` |

??? note "上游说明（英文原文）"

    Codon-optimize an amino acid sequence for expression in a target species (human, ecoli, mouse, or yeast) using species-specific codon frequency tables. Returns the optimized DNA sequence, GC content, Codon Adaptation Index (CAI), and sequence length. Use for: synthetic gene design, heterologous protein expression, improving codon usage for a target host organism.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`molecular_tools`（Plasmids, primers, CRISPR guides and cloning）
- **必填参数：**`sequence`
- **关键词：**crispr, sgrna, spcas9, pam, guide design
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `operation` | string | 否 | Operation (optional; defaults to 'crispr_guide_design' for this tool). 可选值：`crispr_guide_design` 默认 `crispr_guide_design` |
| `sequence` | string | 是 | DNA target sequence (A, T, G, C only; no N), at least 23 bp. |
| `gc_min` | number \| null | 否 | Minimum acceptable protospacer GC%, 0-100 (default 30). 默认 `30.0` |
| `gc_max` | number \| null | 否 | Maximum acceptable protospacer GC%, 0-100 (default 70). 默认 `70.0` |
| `limit` | integer \| null | 否 | Max candidate guides to return, 1-100 (default 10). 默认 `10` |

??? note "上游说明（英文原文）"

    Find candidate SpCas9 (NGG PAM) guide RNAs in a target DNA sequence, on both strands. Ranks 20 bp protospacers by a transparent rule-based heuristic (GC content, PAM-proximal G, Pol III terminator and homopolymer avoidance) -- not a trained on-target efficiency model, and no off-target search against any genome is performed. Use for a first-pass candidate list; verify with CRISPOR or the Broad GPP sgRNA designer before ordering.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`molecular_tools`（Plasmids, primers, CRISPR guides and cloning）
- **必填参数：**`sequence`
- **关键词：**restriction enzyme, recognition site, neb, type iis
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `operation` | string | 否 | Operation (optional; defaults to 'find_restriction_sites' for this tool). 可选值：`find_restriction_sites` 默认 `find_restriction_sites` |
| `sequence` | string | 是 | DNA sequence (A, T, G, C, N only). Case insensitive. Spaces and newlines are ignored. |
| `enzymes` | array \| null | 否 | Optional list of specific enzyme names to check (e.g., \['EcoRI', 'BamHI'\]). If null, checks all 24 NEB enzymes. |

??? note "上游说明（英文原文）"

    Find restriction enzyme recognition sites in a DNA sequence using the NEB enzyme library (25 curated common enzymes plus Biopython's ~600-enzyme library, including the Type IIS enzymes used for Golden Gate/MoClo: BsaI, BbsI, Esp3I/BsmBI, SapI). Returns positions, recognition sequences, and cut counts for each enzyme that has sites. Use for: cloning strategy design, vector linearization planning, molecular biology workflows.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`molecular_tools`（Plasmids, primers, CRISPR guides and cloning）
- **必填参数：**`fragments`
- **关键词：**gibson assembly, overlap, cloning, dna fragments
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `operation` | string | 否 | Operation (optional; defaults to 'gibson_design' for this tool). 可选值：`gibson_design` 默认 `gibson_design` |
| `fragments` | array of string | 是 | List of DNA fragment sequences (at least 2). Each must be longer than overlap_length. |
| `overlap_length` | integer \| null | 否 | Overlap length in bp for Gibson Assembly (default: 20 bp). Must be at least 1. 默认 `20` |

??? note "上游说明（英文原文）"

    Design Gibson Assembly overlaps for seamless DNA fragment assembly. For each fragment, computes the 5' overhang (from the previous fragment's 3' end) and 3' overhang (from the next fragment's 5' end) based on the specified overlap length (default 20 bp). Returns the modified fragment sequences ready for Gibson Assembly cloning. Use for: seamless cloning, multi-fragment assembly, vector backbone insertion.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`molecular_tools`（Plasmids, primers, CRISPR guides and cloning）
- **必填参数：**`fragments`
- **关键词：**golden gate, type iis, bsai, assembly product, plasmid
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `operation` | string | 否 | Operation (optional; defaults to 'golden_gate_assemble' for this tool). 可选值：`golden_gate_assemble` 默认 `golden_gate_assemble` |
| `fragments` | array of string | 是 | The DNA sequences combined in the reaction (plasmids or linear parts), A/T/G/C only. At least 2. |
| `enzyme` | string \| null | 否 | Type IIS enzyme used in the reaction, e.g. 'BsaI' (default), 'BbsI', 'Esp3I'/'BsmBI', 'SapI'. 默认 `BsaI` |
| `circular` | boolean \| null | 否 | Whether the inputs are circular plasmids (default true). Set false for linear parts. 默认 `True` |
| `labels` | array \| null | 否 | Optional names for the inputs (e.g. plasmid names), echoed in the assembly order. |

??? note "上游说明（英文原文）"

    Work out the product of a Golden Gate reaction that has already been performed: given the plasmids/parts that were combined and the Type IIS enzyme used (BsaI, BbsI, Esp3I/BsmBI, SapI), digest each input, discard fragments that keep a recognition site (they are re-cut in the reaction), and ligate the rest by matching 4-bp overhangs into the circular product. Returns the assembled sequence, its length, and the fragment order with the overhang at each junction. Use to answer 'I combined these plasmids with &lt;enzyme&gt; — what does the resulting plasmid contain/express?'. This is the reverse of DNA_golden_gate_design, which designs parts rather than analysing a finished reaction.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`molecular_tools`（Plasmids, primers, CRISPR guides and cloning）
- **必填参数：**`parts`
- **关键词：**golden gate, overhang, bsai, bbsi, cloning
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `operation` | string | 否 | Operation (optional; defaults to 'golden_gate_design' for this tool). 可选值：`golden_gate_design` 默认 `golden_gate_design` |
| `parts` | array of string | 是 | List of DNA part sequences to assemble (at least 2). A, T, G, C only. |
| `enzyme` | string \| null | 否 | Type IIS restriction enzyme to use. Options: 'BsaI' (default) or 'BbsI'. 可选值：`BsaI`, `BbsI`, `None` 默认 `BsaI` |

??? note "上游说明（英文原文）"

    Design Golden Gate Assembly parts with Type IIS restriction enzyme sites (BsaI or BbsI). Assigns unique non-palindromic 4-bp overhangs to each junction and builds the full synthetic part sequences including enzyme recognition sites. Use for: scarless multi-part assembly, modular cloning (MoClo), Golden Gate cloning workflows.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`molecular_tools`（Plasmids, primers, CRISPR guides and cloning）
- **必填参数：**`sequence`
- **关键词：**pcr primer, tm, nearest-neighbor, amplicon
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `operation` | string | 否 | Operation (optional; defaults to 'primer_design' for this tool). 可选值：`primer_design` 默认 `primer_design` |
| `sequence` | string | 是 | DNA template sequence (A, T, G, C, N only). Must be at least 200 bp for good primer design. |
| `target_start` | integer \| null | 否 | 0-based start position of the target region to amplify. Defaults to 0. |
| `target_end` | integer \| null | 否 | 0-based end position of the target region to amplify. Defaults to sequence length. |
| `tm_target` | number \| null | 否 | Target melting temperature in Celsius (default: 60.0°C). 默认 `60.0` |
| `product_size_min` | integer \| null | 否 | Minimum acceptable PCR product size in bp (default: 100). 默认 `100` |
| `product_size_max` | integer \| null | 否 | Maximum acceptable PCR product size in bp (default: 1000). 默认 `1000` |

??? note "上游说明（英文原文）"

    Design forward and reverse PCR primers for a target DNA region using the SantaLucia 1998 nearest-neighbor thermodynamic model. Selects primers with Tm closest to the target (default 60°C), GC content 40-60%, length 18-25 bp, and no 3'-end repeats. Use for: PCR amplification, genotyping, diagnostic assay design, amplicon sequencing.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`molecular_tools`（Plasmids, primers, CRISPR guides and cloning）
- **必填参数：**`sequence`
- **关键词：**restriction digest, fragments, neb, circular plasmid
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `operation` | string | 否 | Operation (optional; defaults to 'virtual_digest' for this tool). 可选值：`virtual_digest` 默认 `virtual_digest` |
| `sequence` | string | 是 | DNA sequence (A, T, G, C, N only). Case insensitive. |
| `enzymes` | array \| null | 否 | List of enzyme names to use for digestion (e.g., \['EcoRI', 'BamHI'\]). If null, uses all 24 NEB enzymes. |
| `circular` | boolean \| null | 否 | Treat sequence as circular DNA (e.g., plasmid). Default: false (linear). 默认 `False` |

??? note "上游说明（英文原文）"

    Perform a virtual restriction enzyme digest on a DNA sequence using the NEB enzyme library (25 curated enzymes plus Biopython's ~600-enzyme library, including Type IIS enzymes such as BsaI, BbsI, Esp3I/BsmBI and SapI, which cut outside their recognition site). Returns all resulting fragments with their sequences, lengths, and positions, plus a list of cut sites. Supports both linear and circular DNA. Use for: cloning planning, gel electrophoresis prediction, verifying construct integrity.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`molecular_tools`（Plasmids, primers, CRISPR guides and cloning）
- **必填参数：**`sequence`
- **关键词：**oligo, primer, melting temperature, tm, gc content
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `sequence` | string | 是 | Oligonucleotide sequence (5' to 3'). Must contain only A, T, C, G (DNA) or A, U, C, G (RNA) bases. Length: 5-200 bases. Example: 'AGCGGATAACAATTTCACACAGGA' |
| `na_concentration_mm` | number \| null | 否 | Sodium (Na+) concentration in millimolar (mM). Affects Tm calculation. Typical range: 0-1000 mM. Default: 50 mM. Standard PCR buffer is ~50 mM. |
| `mg_concentration_mm` | number \| null | 否 | Magnesium (Mg2+) concentration in millimolar (mM). Stabilizes duplexes, raising Tm. Typical range: 0-20 mM. Default: 0 mM. Standard PCR uses 1.5-3.0 mM. |
| `dntps_concentration_mm` | number \| null | 否 | dNTPs concentration in millimolar (mM). dNTPs chelate Mg2+, reducing effective Mg2+. Typical range: 0-10 mM. Default: 0 mM. Standard PCR uses 0.2-0.8 mM. |
| `oligo_concentration_um` | number \| null | 否 | Oligo concentration in micromolar (uM). Affects Tm slightly. Default: 0.25 uM. |
| `oligo_type` | string \| null | 否 | Type of oligonucleotide: 'DNA' or 'RNA'. Default: 'DNA'. |

??? note "上游说明（英文原文）"

    Analyze an oligonucleotide using the IDT OligoAnalyzer API. Calculates melting temperature (Tm), GC content, molecular weight, molar extinction coefficient, and nmol/OD260 and ug/OD260 yields. Tm is calculated using nearest-neighbor thermodynamic parameters (SantaLucia 1998) under user-specified salt conditions (Na+, Mg2+, dNTPs). Supports DNA and RNA oligos up to 200 bases. Use this BEFORE ordering primers/probes to verify Tm, check GC balance, and estimate synthesis yields. Complements NEB_Tm_calculate (which provides polymerase-specific Tm with salt corrections) and DNA_calculate_gc_content (local, no Tm). For self-dimer/homodimer assessment, use IDT_check_self_dimer.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`molecular_tools`（Plasmids, primers, CRISPR guides and cloning）
- **必填参数：**`sequence`
- **关键词：**primer dimer, self-complementarity, delta g, oligo
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `sequence` | string | 是 | Oligonucleotide sequence (5' to 3'). Must contain only A, T, C, G, U bases. Length: 5-200 bases. Example: 'AGCGGATAACAATTTCACACAGGA' |
| `na_concentration_mm` | number \| null | 否 | Sodium (Na+) concentration in millimolar (mM). Affects dimer stability. Default: 50 mM. |
| `mg_concentration_mm` | number \| null | 否 | Magnesium (Mg2+) concentration in millimolar (mM). Default: 0 mM. |
| `temperature_celsius` | number \| null | 否 | Temperature in Celsius for dimer stability calculation. Use your expected annealing temperature. Default: 25. |

??? note "上游说明（英文原文）"

    Check an oligonucleotide for self-complementarity and homodimer (self-dimer) formation using the IDT OligoAnalyzer API. Returns predicted self-dimer structures with their thermodynamic stability (delta-G in kcal/mol), number of base pairs, and alignment diagrams. A more negative delta-G indicates a more stable (problematic) dimer. Risk assessment: delta-G &gt; -3.0 kcal/mol is low risk, -3.0 to -6.0 is moderate risk, &lt; -6.0 is high risk. Use this tool AFTER designing primers (DNA_primer_design or manual) and BEFORE ordering, to verify primers won't form stable homodimers that reduce PCR efficiency. For basic oligo properties (Tm, GC%, MW), use IDT_analyze_oligo instead.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`molecular_tools`（Plasmids, primers, CRISPR guides and cloning）
- **必填参数：**`primer_sequence`
- **关键词：**primer tm, annealing temperature, pcr, nearest-neighbor
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `primer_sequence` | string | 是 | Forward primer DNA sequence (5' to 3'). Must contain only A, T, C, G bases and be at least 10 nt long. Example: 'AGCGGATAACAATTTCACACAGGA' |
| `primer_sequence_2` | string \| null | 否 | Optional reverse primer DNA sequence (5' to 3'). When provided, the API also calculates annealing temperature (Ta). Example: 'CCTCACTAAAGGGAACAAAAGCTG' |
| `polymerase` | string \| null | 否 | NEB polymerase product code. Determines buffer composition for Tm calculation. Common codes: 'q5-0' (Q5 HF, default), 'phusion-0' (Phusion HF Buffer), 'taq-0' (Taq Standard Buffer), 'onetaq-0' (OneTaq Standard Buffer), 'custom' (specify monovalent_salt_mm). Use NEB_Tm_list_polymerases for full list. Default: 'q5-0'. |
| `primer_concentration` | integer \| null | 否 | Primer concentration in nanomolar (nM). Affects Tm calculation. Common values: 200-500 nM. Default: 250. |
| `monovalent_salt_mm` | integer \| null | 否 | Monovalent salt concentration in mM. Only used when polymerase is 'custom'. Typical range: 20-100 mM. Default: 50. |

??? note "上游说明（英文原文）"

    Calculate primer melting temperature (Tm) and annealing temperature (Ta) using the NEB Tm Calculator API. Uses thermodynamic nearest-neighbor calculations (SantaLucia 1998) with polymerase-specific salt corrections for accurate PCR temperature predictions. Supports all NEB DNA polymerases including Q5, Phusion, Taq, OneTaq, and custom buffer conditions. When two primers are provided, also calculates the recommended annealing temperature (Ta). Complements the local DNA_primer_design tool by providing manufacturer-validated Tm values, though results can diverge by 10+°C between this tool, DNA_primer_design, and IDT_analyze_oligo for the same primer since each uses different nearest-neighbor parameter sets and salt-correction models -- treat cross-tool Tm comparisons as approximate, not interchangeable. Common polymerase codes: 'q5-0' (Q5 High-Fidelity), 'phusion-0' (Phusion HF Buffer), 'taq-0' (Taq Standard Buffer), 'onetaq-0' (OneTaq Standard Buffer). Use NEB_Tm_list_polymerases to see all available product codes.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`molecular_tools`（Plasmids, primers, CRISPR guides and cloning）
- **必填参数：**无
- **关键词：**polymerase, q5, phusion, taq
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `filter` | string \| null | 否 | Optional keyword to filter products by name, product code, or catalog number. Case-insensitive. Examples: 'Q5', 'Phusion', 'hot start', 'master mix', 'M0530'. |

??? note "上游说明（英文原文）"

    List all NEB DNA polymerase products available for Tm calculation. Returns product codes (prodcode), product names, and NEB catalog numbers for all 54 supported polymerases/kits. Includes Q5, Phusion, Taq, OneTaq, LongAmp, Vent, Deep Vent, and specialty polymerases. Use the returned prodcode values with NEB_Tm_calculate to get polymerase-specific melting temperatures. Can optionally filter by keyword (e.g., 'Q5', 'Phusion', 'hot start').

安装服务器后通过 DNAKit 调用：

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

- **类别：**`molecular_tools`（Plasmids, primers, CRISPR guides and cloning）
- **必填参数：**`name`
- **关键词：**recognition site, cut site, sticky ends, neb
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `name` | string | 是 | Enzyme name, e.g. 'EcoRI', 'BamHI', 'SmaI'. Matched case-insensitively. |

??? note "上游说明（英文原文）"

    Look up a restriction enzyme in REBASE, the reference catalogue of ~6,000 enzymes from New England Biolabs. Returns the recognition sequence with its cut point, whether the cut leaves blunt or sticky ends, the source organism, methylation sensitivity, commercial availability, and isoschizomers. Example: name='EcoRI' returns G^AATTC. The DNA_find_restriction_sites and DNA_virtual_digest calculators work from 25 curated enzymes; use this to look up anything outside that set.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`molecular_tools`（Plasmids, primers, CRISPR guides and cloning）
- **必填参数：**`name`
- **关键词：**isoschizomer, recognition site, commercial enzyme
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `name` | string | 是 | Enzyme name to find isoschizomers of, e.g. 'EcoRI'. |

??? note "上游说明（英文原文）"

    List every restriction enzyme sharing a prototype with the named enzyme, i.e. all enzymes with the same recognition specificity. Commercially available enzymes are listed first, so this answers 'which enzyme can I actually buy that cuts like this one'. Example: name='EcoRI'.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`molecular_tools`（Plasmids, primers, CRISPR guides and cloning）
- **必填参数：**无
- **关键词：**recognition site, digest, iupac, cloning
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `site` | string \| null | 否 | Exact recognition sequence to match, e.g. 'GAATTC'. Mutually exclusive with sequence. |
| `sequence` | string \| null | 否 | DNA sequence to find cutters for, e.g. a plasmid region. Mutually exclusive with site. |
| `only_commercial` | boolean \| null | 否 | If true, return only commercially available enzymes. |
| `limit` | integer \| null | 否 | Maximum enzymes to return (default 25, max 200). |

??? note "上游说明（英文原文）"

    Find restriction enzymes by what they recognize. Give site to find every enzyme reading an exact recognition sequence (site='GAATTC'), or give sequence to find every enzyme that would cut a DNA sequence, with the number of cut sites each produces. IUPAC ambiguity codes in enzyme sites are expanded, so degenerate recognition sequences match correctly. Set only_commercial=true to restrict to enzymes you can buy.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`molecular_tools`（Plasmids, primers, CRISPR guides and cloning）
- **必填参数：**无
- **关键词：**synthetic biology, biobricks, igem, genetic parts, sbol
- **只读：**是（服务器声明）

??? note "上游说明（英文原文）"

    List all public collections available on SynBioHub. Collections are curated sets of genetic parts and designs, including the iGEM Parts Registry, CIDAR Lab parts, Free Genes feature libraries, and organism-specific collections. Returns collection names, descriptions, URIs, and member counts. Use this to discover available part libraries before searching for specific parts. Note: SynBioHub anonymous REST access was disabled upstream. Every /public/... read returns 401 'Login required'. The tool returns a structured error; users with a SynBioHub account can supply the auth token via SYNBIOHUB_API_TOKEN.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`molecular_tools`（Plasmids, primers, CRISPR guides and cloning）
- **必填参数：**无
- **关键词：**synthetic biology, sbol, dna sequence, promoter, biobrick
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `display_id` | string \| null | 否 | BioBrick/part display ID from the iGEM collection. Examples: 'BBa_E0040' (GFP), 'BBa_R0010' (lac promoter), 'BBa_B0034' (RBS). If not from iGEM, use part_uri instead. |
| `part_uri` | string \| null | 否 | Full SynBioHub URI for the part (from search results). Example: 'https://synbiohub.org/public/igem/BBa_E0040/1'. Use this for parts not in iGEM collection. |

??? note "上游说明（英文原文）"

    Get detailed SBOL (Synthetic Biology Open Language) information for a specific genetic part from SynBioHub. Returns the part's title, description, DNA sequence, sequence length, functional roles, creation date, and provenance. Can look up by full URI or by iGEM BioBrick display ID (e.g., BBa_E0040 for GFP). Example: get details for BBa_E0040 (GFP coding sequence, 720bp) including its full DNA sequence. Note: SynBioHub anonymous REST access was disabled upstream. Every /public/... read returns 401 'Login required'. The tool returns a structured error; users with a SynBioHub account can supply the auth token via SYNBIOHUB_API_TOKEN.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`molecular_tools`（Plasmids, primers, CRISPR guides and cloning）
- **必填参数：**`query`
- **关键词：**synthetic biology, igem registry, promoter, terminator, biobricks
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `query` | string | 是 | Search query for finding genetic parts. Can include gene names, part functions, BioBrick IDs, or biological descriptions. Examples: 'GFP', 'lac promoter', 'BBa_E0040', 'terminator', 'riboswitch', 'CRISPR', 'T7 promoter'. |
| `offset` | integer \| null | 否 | Offset for pagination (default: 0). |
| `limit` | integer \| null | 否 | Maximum number of results to return (default: 10, max: 50). |

??? note "上游说明（英文原文）"

    Search the SynBioHub synthetic biology parts repository for genetic parts and designs by keyword. SynBioHub hosts the iGEM Registry of Standard Biological Parts (20,000+ BioBricks) and other public SBOL collections. Parts include promoters, coding sequences (CDS/genes), terminators, ribosome binding sites (RBS), reporters (GFP, RFP, LacZ), regulatory elements, and composite genetic devices. Returns part IDs, names, descriptions, SBOL types, and functional roles. Example: search 'GFP' finds green fluorescent protein parts; search 'lac promoter' finds lac operon regulatory parts. Note: SynBioHub anonymous REST access was disabled upstream. Every /public/... read returns 401 'Login required'. The tool returns a structured error; users with a SynBioHub account can supply the auth token via SYNBIOHUB_API_TOKEN.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'SynBioHub_search_parts', 'arguments': {'query': '...'}},
)
```
