# ToolUniverse: sequence (`tooluniverse`)

DNAKit's MCP tool index keeps **845** DNA/gene tools of this server (2927 listed, read on 2026-10-16). They run once you install the server; DNAKit does not bundle them. See [calling a tool that was found](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found).

This page lists the 44 tools in category `sequence`; other categories: [`genome_annotation`](tooluniverse-genome_annotation.md) (253), [`variant`](tooluniverse-variant.md) (155), [`association`](tooluniverse-association.md) (138), [`expression`](tooluniverse-expression.md) (56), [`regulation_epigenome`](tooluniverse-regulation_epigenome.md) (54), [`sequencing_data`](tooluniverse-sequencing_data.md) (91), [`prediction_model`](tooluniverse-prediction_model.md) (17), [`molecular_tools`](tooluniverse-molecular_tools.md) (23), [`breeding_population`](tooluniverse-breeding_population.md) (14).

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
| [`BLAST_nucleotide_search`](#BLAST_nucleotide_search) | `sequence` | Runs remote NCBI blastn for a nucleotide sequence and returns hits. |
| [`BOLDSystems_get_record`](#BOLDSystems_get_record) | `sequence` | Retrieves a BOLD DNA barcode record by process ID. |
| [`BOLDSystems_search_by_bin`](#BOLDSystems_search_by_bin) | `sequence` | Lists BOLD DNA barcode records sharing a Barcode Index Number. |
| [`BOLDSystems_search_by_taxon`](#BOLDSystems_search_by_taxon) | `sequence` | Searches BOLD DNA barcode records by taxon and country. |
| [`dbfetch_fetch_batch`](#dbfetch_fetch_batch) | `sequence` | Fetches multiple EBI Dbfetch entries by IDs from a chosen database in a chosen format such as FASTA. |
| [`dbfetch_fetch_entry`](#dbfetch_fetch_entry) | `sequence` | Fetches a single EBI Dbfetch entry by database and ID in FASTA, XML or other formats. |
| [`DNA_calculate_gc_content`](#DNA_calculate_gc_content) | `sequence` | Calculates GC content and nucleotide composition of a DNA sequence. |
| [`DNA_find_orfs`](#DNA_find_orfs) | `sequence` | Finds ORFs in a DNA sequence on chosen strands, sorted by length with positions and frames. |
| [`DNA_reverse_complement`](#DNA_reverse_complement) | `sequence` | Returns the reverse complement of a DNA sequence. |
| [`DNA_translate_reading_frames`](#DNA_translate_reading_frames) | `sequence` | Translates a DNA sequence to protein in the three forward reading frames using a chosen genetic code and reports the longest ORF per frame. |
| [`DNA_translate_sequence`](#DNA_translate_sequence) | `sequence` | Translates a DNA coding sequence to protein using a chosen genetic code. |
| [`EBI_build_phylogenetic_tree`](#EBI_build_phylogenetic_tree) | `sequence` | Builds a neighbour-joining or UPGMA Newick tree from an aligned FASTA via EMBL-EBI simple_phylogeny. |
| [`EBI_msa_align`](#EBI_msa_align) | `sequence` | Runs de novo multiple sequence alignment (Clustal Omega, MUSCLE, MAFFT etc.) via EMBL-EBI. |
| [`EBI_pairwise_align`](#EBI_pairwise_align) | `sequence` | Aligns two sequences with EMBOSS needle, stretcher or water and reports identity and score. |
| [`EBI_translate_sequence`](#EBI_translate_sequence) | `sequence` | Translates nucleotide to protein, back-translates, or runs six-frame translation with ORFs via EMBOSS. |
| [`ena_get_entry`](#ena_get_entry) | `sequence` | Gets ENA entry information for an EMBL/GenBank accession. |
| [`ena_get_entry_history`](#ena_get_entry_history) | `sequence` | Gets version history for an ENA EMBL/GenBank entry. |
| [`ena_get_entry_summary`](#ena_get_entry_summary) | `sequence` | Gets a detailed ENA entry summary for an EMBL/GenBank accession. |
| [`ena_get_sequence_embl`](#ena_get_sequence_embl) | `sequence` | Gets an EMBL-formatted nucleotide record from ENA by accession. |
| [`ena_get_sequence_fasta`](#ena_get_sequence_fasta) | `sequence` | Gets a nucleotide sequence in FASTA format from ENA by accession. |
| [`ensembl_get_alignment`](#ensembl_get_alignment) | `sequence` | Gets EPO/PECAN multi-species genomic alignments for a region from Ensembl. |
| [`ensembl_get_sequence`](#ensembl_get_sequence) | `sequence` | Gets sequence for an Ensembl gene, transcript or region in genomic, CDS, cDNA or peptide form. |
| [`EnsemblSeq_get_id_sequence`](#EnsemblSeq_get_id_sequence) | `sequence` | Gets protein, cDNA/CDS or genomic sequence for an Ensembl ID. |
| [`EnsemblSeq_get_region_sequence`](#EnsemblSeq_get_region_sequence) | `sequence` | Gets the reference DNA sequence for a chromosomal region from Ensembl. |
| [`IMGT_get_germline_gene_fasta`](#IMGT_get_germline_gene_fasta) | `sequence` | Retrieves IMGT/GENE-DB germline IG/TR V/D/J-region reference sequences in FASTA for a gene type and species. |
| [`IMGT_get_sequence`](#IMGT_get_sequence) | `sequence` | Gets an immunoglobulin/T cell receptor sequence from IMGT/LIGM-DB by accession in FASTA or EMBL format. |
| [`NCBI_fetch_accessions`](#NCBI_fetch_accessions) | `sequence` | Converts NCBI Nucleotide UIDs to GenBank/EMBL/RefSeq accession numbers. |
| [`NCBI_get_sequence`](#NCBI_get_sequence) | `sequence` | Retrieves nucleotide sequences from NCBI by accession in FASTA, GenBank or EMBL format. |
| [`NCBI_search_nucleotide`](#NCBI_search_nucleotide) | `sequence` | Searches NCBI Nucleotide (GenBank/EMBL/RefSeq) by organism, gene, strain, keywords or sequence type. |
| [`Pathoplexus_count_sequences`](#Pathoplexus_count_sequences) | `sequence` | Counts and aggregates open pathogen genome sequences in Pathoplexus by country, lineage or other fields. |
| [`Pathoplexus_get_sequences_fasta`](#Pathoplexus_get_sequences_fasta) | `sequence` | Downloads unaligned or aligned nucleotide or amino-acid FASTA sequences for a Pathoplexus organism. |
| [`phykit_batch_analysis`](#phykit_batch_analysis) | `sequence` | Runs PhyKIT functions (treeness, saturation, parsimony-informative sites) over a directory of alignments or trees. |
| [`Rfam_get_alignment`](#Rfam_get_alignment) | `sequence` | Returns the Rfam seed multiple sequence alignment in Stockholm, Pfam or FASTA format. |
| [`Rfam_get_covariance_model`](#Rfam_get_covariance_model) | `sequence` | Returns the Infernal covariance model for an Rfam family. |
| [`Rfam_get_tree_data`](#Rfam_get_tree_data) | `sequence` | Returns the NHX phylogenetic tree for an Rfam family seed alignment. |
| [`Rfam_search_sequence`](#Rfam_search_sequence) | `sequence` | Searches an RNA sequence against all Rfam families with Infernal cmscan. |
| [`RNAcentral_get_sequence`](#RNAcentral_get_sequence) | `sequence` | Returns the FASTA nucleotide sequence of an RNAcentral ncRNA. |
| [`Sequence_count_residues`](#Sequence_count_residues) | `sequence` | Counts occurrences and positions of a given nucleotide or amino acid residue in a sequence string or a UniProt-fetched protein sequence. |
| [`Sequence_dn_ds`](#Sequence_dn_ds) | `sequence` | Computes dN/dS (Ka/Ks) between two codon-aligned coding sequences with Nei-Gojobori and Jukes-Cantor correction. |
| [`Sequence_gc_content`](#Sequence_gc_content) | `sequence` | Calculates GC count, fraction and percentage and full nucleotide composition of a DNA or RNA sequence. |
| [`Sequence_reverse_complement`](#Sequence_reverse_complement) | `sequence` | Returns the reverse complement of a DNA sequence (A/T/C/G/N). |
| [`Sequence_stats`](#Sequence_stats) | `sequence` | Computes length, composition and type-specific metrics (GC% for nucleotides, molecular weight for proteins) for a sequence, auto-detecting its type. |
| [`SGD_get_sequence`](#SGD_get_sequence) | `sequence` | Gets genomic DNA, coding DNA and protein sequences plus genomic coordinates for a yeast gene in reference strain S288C from SGD. |
| [`UCSC_get_sequence`](#UCSC_get_sequence) | `sequence` | Gets DNA sequence for a genomic region of a UCSC genome assembly given a locus or chrom/start/end. |

## `BLAST_nucleotide_search` {#BLAST_nucleotide_search}

Runs remote NCBI blastn for a nucleotide sequence and returns hits.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `sequence`
- **Keywords:** blastn, nucleotide alignment, homology search, ncbi
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sequence` | string | yes | DNA sequence to search (minimum 10 nucleotides) |
| `database` | string | no | Database to search. Options: 'nt' (comprehensive but slow), 'refseq_select_rna' (curated, faster), 'refseq_representative_genomes' (representative genomes) Default `nt` |
| `expect` | number | no | E-value threshold for reporting hits Default `10.0` |
| `hitlist_size` | integer | no | Maximum number of hits to return Default `50` |

??? note "Upstream description"

    Search nucleotide sequences using NCBI BLAST blastn against nucleotide databases. Requires Biopython (pip install biopython). IMPORTANT: This tool uses NCBI's remote BLAST service which is inherently SLOW - expect 5-30 minutes per query depending on server load. This is normal NCBI behavior, not a bug. For faster results, use smaller databases like 'refseq_select_rna' instead of 'nt'.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'BLAST_nucleotide_search',
     'arguments': {'sequence': '...'}},
)
```

## `BOLDSystems_get_record` {#BOLDSystems_get_record}

Retrieves a BOLD DNA barcode record by process ID.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `process_id`
- **Keywords:** dna barcode, coi, bold, taxonomy
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `process_id` | string | yes | BOLD process ID, e.g. 'ABRMM002-06'. |

??? note "Upstream description"

    Retrieve a single BOLD (Barcode of Life Data System) DNA barcode record by its BOLD process ID. Returns full taxonomy, BIN cluster, collection metadata, and COI-5P INSDC accession. Example: process_id='ABRMM002-06' returns a Panthera leo record collected at the Metro Toronto Zoo.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'BOLDSystems_get_record',
     'arguments': {'process_id': '...'}},
)
```

## `BOLDSystems_search_by_bin` {#BOLDSystems_search_by_bin}

Lists BOLD DNA barcode records sharing a Barcode Index Number.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `bin_id`
- **Keywords:** barcode index number, bin, species cluster, coi
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `bin_id` | string | yes | BOLD BIN identifier, e.g. 'BOLD:AAD6819' (the 'BOLD:' prefix is added automatically if omitted). |
| `limit` | integer \| null | no | Max records to return, 1-200. Default 30. |

??? note "Upstream description"

    Search BOLD (Barcode of Life Data System) for all DNA barcode records sharing a BIN (Barcode Index Number) -- BOLD's algorithmic, sequence-similarity-based species-cluster identifier, which often resolves identity where taxonomy alone cannot. Example: bin_id='BOLD:AAD6819' returns 51 specimen records of Panthera leo sharing that barcode cluster. This is a capability unique to BOLD; ToolUniverse's iDigBio tools have no genetic barcode or BIN data.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'BOLDSystems_search_by_bin',
     'arguments': {'bin_id': '...'}},
)
```

## `BOLDSystems_search_by_taxon` {#BOLDSystems_search_by_taxon}

Searches BOLD DNA barcode records by taxon and country.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `taxon_name`
- **Keywords:** dna barcode, taxon, biodiversity, coi
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `taxon_name` | string | yes | Taxon name matching rank, e.g. 'Panthera' (genus), 'Felidae' (family), 'Carnivora' (order), or 'Panthera leo' (species). |
| `rank` | string \| null | no | Taxonomic rank of taxon_name. Default 'genus'. One of `genus`, `family`, `order`, `species`, `None` |
| `country` | string \| null | no | Optional country/ocean name to narrow results, e.g. 'Canada'. |
| `limit` | integer \| null | no | Max records to return, 1-200. Default 30. |

??? note "Upstream description"

    Search BOLD (Barcode of Life Data System) DNA barcode records by taxon. rank='genus'/'family'/'order' filters server-side; rank='species' (e.g. 'Panthera leo') queries at genus level and filters client-side, since BOLD does not support server-side species filtering. Optionally narrow by country. Example: taxon_name='Panthera', rank='genus' returns all big-cat barcode records with BIN cluster IDs, COI-5P accession numbers, and collection metadata.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'BOLDSystems_search_by_taxon',
     'arguments': {'taxon_name': '...'}},
)
```

## `dbfetch_fetch_batch` {#dbfetch_fetch_batch}

Fetches multiple EBI Dbfetch entries by IDs from a chosen database in a chosen format such as FASTA.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `db`, `ids`
- **Keywords:** dbfetch, fasta, embl, ena, batch retrieval
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `db` | string | yes | Database name |
| `ids` | string \| array | yes | Comma-separated IDs or list of IDs (e.g., 'P05067,P04637' or \['P05067', 'P04637'\]) |
| `format` | string | no | Output format (default: 'fasta') Default `fasta` |

??? note "Upstream description"

    Fetch multiple database entries by IDs in batch. Supports comma-separated IDs or list of IDs. Returns entries in specified format.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'dbfetch_fetch_batch',
     'arguments': {'db': '...', 'ids': '...'}},
)
```

## `dbfetch_fetch_entry` {#dbfetch_fetch_entry}

Fetches a single EBI Dbfetch entry by database and ID in FASTA, XML or other formats.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `db`, `id`
- **Keywords:** dbfetch, fasta, embl, ena, entry retrieval
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `db` | string | yes | Database name (e.g., 'uniprotkb', 'pdb', 'embl', 'ensembl') |
| `id` | string | yes | Entry ID (e.g., UniProt accession 'P05067', PDB ID '1A2B') |
| `format` | string | no | Output format (e.g., 'fasta', 'xml', 'json', 'txt'). Default: 'fasta' One of `fasta`, `xml`, `json`, `txt`, `gff`, `embl`, `genbank` Default `fasta` |

??? note "Upstream description"

    Fetch a single database entry by ID from various databases (UniProt, PDB, etc.) in specified format. Supports FASTA, XML, JSON, and other formats.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'dbfetch_fetch_entry',
     'arguments': {'db': '...', 'id': '...'}},
)
```

## `DNA_calculate_gc_content` {#DNA_calculate_gc_content}

Calculates GC content and nucleotide composition of a DNA sequence.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `sequence`
- **Keywords:** gc content, nucleotide composition, base counts
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | Operation (optional; defaults to 'calculate_gc_content' for this tool). One of `calculate_gc_content` Default `calculate_gc_content` |
| `sequence` | string | yes | DNA sequence (A, T, G, C, N only). Case insensitive. |

??? note "Upstream description"

    Calculate GC content percentage and full nucleotide composition of a DNA sequence. Returns GC%, AT%, individual base counts (A, T, G, C, N), and interpretation (low/normal/high GC). Useful for primer design quality checks, cloning efficiency prediction, and sequence characterization. High GC (&gt;60%) sequences may require special PCR conditions.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'DNA_calculate_gc_content',
     'arguments': {'sequence': '...'}},
)
```

## `DNA_find_orfs` {#DNA_find_orfs}

Finds ORFs in a DNA sequence on chosen strands, sorted by length with positions and frames.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `sequence`
- **Keywords:** orf, open reading frame, start codon, stop codon
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | Operation (optional; defaults to 'find_orfs' for this tool). One of `find_orfs` Default `find_orfs` |
| `sequence` | string | yes | DNA sequence (A, T, G, C, N only) |
| `min_length` | integer | no | Minimum ORF length in nucleotides (default: 100 nt = ~33 amino acids) Default `100` |
| `strand` | string \| null | no | Which strand to search: 'forward', 'reverse', or 'both' (default) One of `forward`, `reverse`, `both`, `None` Default `both` |

??? note "Upstream description"

    Find open reading frames (ORFs) in a DNA sequence by scanning for ATG start codons followed by stop codons (TAA, TAG, TGA). Returns ORFs sorted by length with positions, strand, and reading frame. Supports forward strand, reverse strand, or both. Use for: gene prediction, identifying protein-coding regions, analyzing synthetic sequences.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'DNA_find_orfs', 'arguments': {'sequence': '...'}},
)
```

## `DNA_reverse_complement` {#DNA_reverse_complement}

Returns the reverse complement of a DNA sequence.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `sequence`
- **Keywords:** reverse complement, antisense, dna strand
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | Operation (optional; defaults to 'reverse_complement' for this tool). One of `reverse_complement` Default `reverse_complement` |
| `sequence` | string | yes | DNA sequence (A, T, G, C, N only). Case insensitive. |

??? note "Upstream description"

    Generate the reverse complement of a DNA sequence. Complements each base (A↔T, G↔C) then reverses the sequence. Essential for: designing reverse primers, understanding antisense strand, converting sequences between conventions. N bases are preserved as N in the reverse complement.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'DNA_reverse_complement',
     'arguments': {'sequence': '...'}},
)
```

## `DNA_translate_reading_frames` {#DNA_translate_reading_frames}

Translates a DNA sequence to protein in the three forward reading frames using a chosen genetic code and reports the longest ORF per frame.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `sequence`
- **Keywords:** orf finding, open reading frame, genetic code, codon, translation, protein from dna
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | Operation type One of `translate_reading_frames` Default `translate_reading_frames` |
| `sequence` | string | yes | DNA sequence (A, T, G, C only). Case insensitive. Spaces and newlines are stripped. |
| `frame` | string \| null | no | Reading frame to translate: '1', '2', '3', or 'all' (default). Frame 1 starts at position 0, frame 2 at position 1, frame 3 at position 2. One of `1`, `2`, `3`, `all`, `None` Default `all` |
| `genetic_code` | integer \| null | no | NCBI genetic code number. Only 1 (the standard code) is implemented; any other value is rejected with an error rather than translated with the standard table. Alternative code tables (2 = vertebrate mitochondrial, 3 = yeast mitochondrial, 11 = bacterial, etc.) are not supported by this tool. One of `1`, `None` Default `1` |

??? note "Upstream description"

    Translate a DNA sequence to protein in all 3 forward reading frames and identify the longest ORF. Uses the standard genetic code (NCBI Code 1). Returns protein sequences for each frame, longest ORF per frame, and the best frame with longest uninterrupted protein. Use for: reading frame analysis, finding the correct ORF in an unknown sequence, comparing all three frames.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'DNA_translate_reading_frames',
     'arguments': {'sequence': '...'}},
)
```

## `DNA_translate_sequence` {#DNA_translate_sequence}

Translates a DNA coding sequence to protein using a chosen genetic code.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `sequence`
- **Keywords:** translation, genetic code, codon, protein sequence
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | Operation (optional; defaults to 'translate_sequence' for this tool). One of `translate_sequence` Default `translate_sequence` |
| `sequence` | string | yes | DNA coding sequence starting with ATG (A, T, G, C only). Should be in-frame. |
| `codon_table` | string \| null | no | Genetic code to use. Currently only 'standard' (NCBI Code 1) is supported. One of `standard`, `None` Default `standard` |

??? note "Upstream description"

    Translate a DNA coding sequence to protein using the standard genetic code (NCBI Code 1). Returns protein sequence (single-letter amino acids), stop codon positions, and protein length. Trims at first stop codon by default. Use for: checking synthetic gene constructs, verifying reading frames, predicting protein sequence from genomic DNA.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'DNA_translate_sequence',
     'arguments': {'sequence': '...'}},
)
```

## `EBI_build_phylogenetic_tree` {#EBI_build_phylogenetic_tree}

Builds a neighbour-joining or UPGMA Newick tree from an aligned FASTA via EMBL-EBI simple_phylogeny.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `aligned_sequences`
- **Keywords:** phylogeny, newick, neighbour-joining, upgma, msa
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `aligned_sequences` | string | yes | Aligned sequences in FASTA format (&gt;=2 records, equal length with gaps as '-'). Typically the aligned_fasta from EBI_msa_align. |
| `clustering` | string | no | Tree-building method: 'Neighbour-joining' (default) or 'UPGMA'. Default `Neighbour-joining` |
| `distance_correction` | boolean | no | Apply Kimura distance correction for multiple substitutions (recommended for divergent sequences). Default `False` |
| `exclude_gaps` | boolean | no | Exclude alignment columns containing any gap before computing distances. Default `False` |

??? note "Upstream description"

    Build a phylogenetic tree from a multiple sequence alignment via EMBL-EBI simple_phylogeny. Takes an ALIGNED FASTA (for example the aligned_fasta output of EBI_msa_align) and returns a Newick tree by neighbour-joining (default) or UPGMA, optionally with Kimura distance correction and gap-column exclusion. Job-based: submits and polls to completion. Use after EBI_msa_align to infer evolutionary relationships among sequences; for raw unaligned sequences, run EBI_msa_align first.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'EBI_build_phylogenetic_tree',
     'arguments': {'aligned_sequences': '...'}},
)
```

## `EBI_msa_align` {#EBI_msa_align}

Runs de novo multiple sequence alignment (Clustal Omega, MUSCLE, MAFFT etc.) via EMBL-EBI.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `sequences`
- **Keywords:** msa, clustal omega, mafft, muscle, alignment
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sequences` | string | yes | Two or more sequences in FASTA format (each starting with '&gt;'). Example: '&gt;seq1\\nMVLSPADK...\\n&gt;seq2\\nMVHLTPEE...' |
| `method` | string | no | Alignment program: 'clustalo' (Clustal Omega, default, scales to many sequences), 'muscle', 'mafft', 'kalign', or 'tcoffee'. Default `clustalo` |
| `sequence_type` | string | no | Sequence molecule type: 'protein' (default), 'dna', or 'rna'. Ignored by muscle (auto-detected). Default `protein` |

??? note "Upstream description"

    Multiple sequence alignment (MSA) of user-provided sequences via EMBL-EBI Job Dispatcher. Unlike pre-computed alignment lookups (Rfam, Ensembl Compara), this aligns YOUR OWN protein/DNA/RNA sequences de novo. Submit &gt;=2 sequences in FASTA, pick an aligner (clustalo \[default\], muscle, mafft, kalign, tcoffee), and get back the alignment in FASTA and Clustal formats plus a guide tree. Job-based: submits and polls to completion (typically seconds to ~2 min). Use for aligning homologs, building an MSA before phylogenetics, conservation analysis, or designing primers/probes against conserved regions.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'EBI_msa_align', 'arguments': {'sequences': '...'}},
)
```

## `EBI_pairwise_align` {#EBI_pairwise_align}

Aligns two sequences with EMBOSS needle, stretcher or water and reports identity and score.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `sequence_a`, `sequence_b`
- **Keywords:** pairwise alignment, needleman-wunsch, smith-waterman, emboss
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sequence_a` | string | yes | First sequence, FASTA or raw. |
| `sequence_b` | string | yes | Second sequence, FASTA or raw. |
| `algorithm` | string \| null | no | 'needle' (global, default), 'stretcher' (global, long sequences), or 'matcher' (local). |
| `sequence_type` | string \| null | no | 'protein' (default) or 'dna'. |
| `matrix` | string \| null | no | Scoring matrix, e.g. 'EBLOSUM62' for protein or 'EDNAFULL' for DNA. |
| `gap_open` | number \| null | no | Gap opening penalty. |
| `gap_extend` | number \| null | no | Gap extension penalty. |

??? note "Upstream description"

    Align exactly two sequences via EMBL-EBI EMBOSS and return the alignment with percent identity, similarity, gaps, and score. Choose algorithm='needle' for Needleman-Wunsch global alignment (default), 'stretcher' for memory-efficient global alignment of long sequences, or 'matcher' for Smith-Waterman local alignment. Use EBI_msa_align instead when aligning three or more sequences. Jobs are submitted and polled, so a call takes several seconds.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'EBI_pairwise_align',
     'arguments': {'sequence_a': '...', 'sequence_b': '...'}},
)
```

## `EBI_translate_sequence` {#EBI_translate_sequence}

Translates nucleotide to protein, back-translates, or runs six-frame translation with ORFs via EMBOSS.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `sequence`
- **Keywords:** transeq, backtranseq, six-frame, orf, translation
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sequence` | string | yes | Input sequence, FASTA or raw. |
| `mode` | string \| null | no | 'dna_to_protein' (default), 'protein_to_dna', or 'six_frame'. |
| `frame` | string \| null | no | Reading frame for dna_to_protein: '1','2','3','-1','-2','-3', or '6' for all. |
| `codon_table` | string \| null | no | Genetic code, e.g. '0' standard, '2' vertebrate mitochondrial, '11' bacterial. |
| `min_orf_size` | integer \| null | no | Minimum ORF size in nucleotides for six_frame mode. |

??? note "Upstream description"

    Translate nucleotide to protein, back-translate protein to nucleotide, or produce a six-frame translation with ORF detection, via EMBL-EBI EMBOSS. Set mode='dna_to_protein' (default, transeq), 'protein_to_dna' (backtranseq), or 'six_frame' (sixpack). Supports alternative genetic codes through codon_table. Complements the local DNA_* helpers by handling alternative codon tables and six-frame ORF calling server-side.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'EBI_translate_sequence',
     'arguments': {'sequence': '...'}},
)
```

## `ena_get_entry` {#ena_get_entry}

Gets ENA entry information for an EMBL/GenBank accession.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `accession`
- **Keywords:** ena, embl, genbank, accession
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `accession` | string | yes | EMBL/GenBank accession number. NOT RefSeq (NC_*, NM_*, NP_*). Examples: 'U00096', 'AJ312385' |
| `expanded` | boolean | no | Return expanded metadata (default: false) Default `False` |

??? note "Upstream description"

    Get entry information from ENA by accession number. Supports EMBL/GenBank accessions only (e.g., U00096, AJ312385, M15390). Does NOT support RefSeq accessions (NC_*, NM_*, NP_* prefixes). Currently extracts metadata from FASTA header. For comprehensive metadata including project, sample, taxon, publications, and sequence length, use ena_get_entry_summary instead. Returns sequence header information with basic metadata.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ena_get_entry', 'arguments': {'accession': '...'}},
)
```

## `ena_get_entry_history` {#ena_get_entry_history}

Gets version history for an ENA EMBL/GenBank entry.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `accession`
- **Keywords:** ena, versions, sequence history
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `accession` | string | yes | EMBL/GenBank accession number. NOT RefSeq (NC_*, NM_*, NP_*). Examples: 'U00096', 'AJ312385' |

??? note "Upstream description"

    Get version history for an ENA entry by accession number. Supports EMBL/GenBank accessions only (e.g., U00096, AJ312385, M15390). Does NOT support RefSeq accessions (NC_*, NM_*, NP_* prefixes). Returns all versions of the entry with version numbers, dates, and status information. Each version includes sequence version number, first public date, last updated date, and current status.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ena_get_entry_history',
     'arguments': {'accession': '...'}},
)
```

## `ena_get_entry_summary` {#ena_get_entry_summary}

Gets a detailed ENA entry summary for an EMBL/GenBank accession.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `accession`
- **Keywords:** ena, entry summary, embl, genbank
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `accession` | string | yes | EMBL/GenBank accession number. NOT RefSeq (NC_*, NM_*, NP_*). Examples: 'U00096', 'AJ312385', 'M15390' |

??? note "Upstream description"

    Get comprehensive summary information for an ENA entry by accession number. Supports EMBL/GenBank accessions only (e.g., U00096, AJ312385, M15390). Does NOT support RefSeq accessions (NC_*, NM_*, NP_* prefixes). Returns detailed metadata including description, version, project ID, sample ID, molecule type, taxon ID, data type, publication references, sequence length, and status information. This provides more complete metadata than ena_get_entry.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ena_get_entry_summary',
     'arguments': {'accession': '...'}},
)
```

## `ena_get_sequence_embl` {#ena_get_sequence_embl}

Gets an EMBL-formatted nucleotide record from ENA by accession.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `accession`
- **Keywords:** ena, embl format, nucleotide record
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `accession` | string | yes | EMBL/GenBank accession number. NOT RefSeq (NC_*, NM_*, NP_*). Examples: 'U00096', 'AJ312385', 'M15390' |

??? note "Upstream description"

    Get nucleotide sequence in EMBL format from ENA by accession number. Supports EMBL/GenBank accessions only (e.g., U00096, AJ312385, M15390). Does NOT support RefSeq accessions (NC_*, NM_*, NP_* prefixes). Returns EMBL-formatted sequence record.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ena_get_sequence_embl',
     'arguments': {'accession': '...'}},
)
```

## `ena_get_sequence_fasta` {#ena_get_sequence_fasta}

Gets a nucleotide sequence in FASTA format from ENA by accession.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `accession`
- **Keywords:** ena, fasta, nucleotide sequence, download
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `accession` | string | yes | EMBL/GenBank accession number. Supported formats: U\#\#\#\#\#, M\#\#\#\#\#, AJ\#\#\#\#\#\#, AL\#\#\#\#\#\#, CP\#\#\#\#\#\#, etc. NOT RefSeq (NC_*, NM_*, NP_*). Examples: 'U00096', 'AJ312385', 'M15390' |
| `download` | boolean | no | Force download instead of inline display (default: false) Default `False` |

??? note "Upstream description"

    Get nucleotide sequence in FASTA format from ENA by accession number. Supports EMBL/GenBank accessions only (e.g., U00096, AJ312385, M15390). Does NOT support RefSeq accessions (NC_*, NM_*, NP_* prefixes). Returns sequence with header information.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ena_get_sequence_fasta',
     'arguments': {'accession': '...'}},
)
```

## `ensembl_get_alignment` {#ensembl_get_alignment}

Gets EPO/PECAN multi-species genomic alignments for a region from Ensembl.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `species`, `region`
- **Keywords:** genomic alignment, epo, conservation, comparative genomics
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `species` | string | yes | Source species (e.g., 'human', 'homo_sapiens') |
| `region` | string | yes | Genomic region in format 'chromosome:start-end' (e.g., '2:106000000-106100000') |
| `species_set_group` | string | no | Predefined species group for multiple alignment. Options: 'mammals', 'primates', 'fish', 'birds', 'sauropsids'. Required for EPO method. Default `mammals` |
| `method` | string | no | Alignment method: EPO (high quality multiple), EPO_EXTENDED (extended coverage), PECAN (low coverage) One of `EPO`, `EPO_EXTENDED`, `PECAN` Default `EPO` |
| `compact` | boolean | no | Compact output format Default `True` |
| `display_species_set` | string | no | Comma-separated list to filter which species to display in results (optional, e.g., 'mouse,rat') |

??? note "Upstream description"

    Get genomic alignments between species for a genomic region. Returns multiple alignments showing conserved sequences across species groups. Supports EPO and PECAN alignment methods for predefined species groups (mammals, primates, fish, birds, sauropsids). Useful for studying sequence conservation.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ensembl_get_alignment',
     'arguments': {'species': '...', 'region': '...'}},
)
```

## `ensembl_get_sequence` {#ensembl_get_sequence}

Gets sequence for an Ensembl gene, transcript or region in genomic, CDS, cDNA or peptide form.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `id`
- **Keywords:** ensembl, genomic sequence, cdna, cds
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `id` | string | yes | Ensembl gene/transcript ID (e.g., 'ENSG00000139618' or 'ENST00000380152') |
| `type` | string | no | Sequence type: 'genomic' (genomic DNA), 'cds' (coding sequence), 'cdna' (cDNA), 'protein' (protein/peptide sequence -- 'peptide' is also accepted and mapped to 'protein', the name Ensembl's own REST API actually expects) One of `genomic`, `cds`, `cdna`, `protein`, `peptide` Default `genomic` |
| `species` | string | no | Species name (optional, defaults to human) |
| `multiple_sequences` | boolean | no | Return multiple sequences if available (e.g., all transcripts) Default `True` |

??? note "Upstream description"

    Get DNA or protein sequence for a gene, transcript, or genomic region. Returns sequence in JSON format with metadata. Supports multiple sequence types: genomic, CDS, cDNA, and peptide sequences.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ensembl_get_sequence', 'arguments': {'id': '...'}},
)
```

## `EnsemblSeq_get_id_sequence` {#EnsemblSeq_get_id_sequence}

Gets protein, cDNA/CDS or genomic sequence for an Ensembl ID.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `ensembl_id`
- **Keywords:** ensembl, sequence, cdna, cds, protein
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `ensembl_id` | string | yes | Ensembl stable ID. Examples: 'ENSP00000269305' (TP53 protein, 393 aa), 'ENST00000269305' (TP53-201 transcript), 'ENSG00000141510' (TP53 gene). Protein IDs start with ENSP, transcript with ENST, gene with ENSG. |
| `type` | string | no | Sequence type to return. Options: 'protein' (amino acid), 'cdna' (spliced cDNA), 'cds' (coding sequence), 'genomic' (unspliced DNA). Default: 'protein' for ENSP IDs, auto-detected for others. Default `protein` |

??? note "Upstream description"

    Get the amino acid or nucleotide sequence for an Ensembl ID (protein, transcript, or gene). Retrieves protein sequences from ENSP IDs, cDNA/CDS from ENST IDs, or genomic DNA from ENSG IDs. Example: ENSP00000269305 returns the 393 amino acid protein sequence of TP53 canonical isoform. ENST00000269305 with type='cdna' returns the 2512 bp cDNA sequence of TP53-201.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'EnsemblSeq_get_id_sequence',
     'arguments': {'ensembl_id': '...'}},
)
```

## `EnsemblSeq_get_region_sequence` {#EnsemblSeq_get_region_sequence}

Gets the reference DNA sequence for a chromosomal region from Ensembl.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `region`
- **Keywords:** reference genome, dna sequence, region
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `region` | string | yes | Genomic region in format 'chr:start-end' or 'chr:start..end:strand'. Examples: '17:7668421-7668520' (TP53 locus, 100bp), '7:140424943-140425043' (BRAF locus, 100bp), 'X:15560138-15602945' (MECP2 gene). Strand: 1 for forward, -1 for reverse. |
| `species` | string | no | Ensembl species name. Default: 'homo_sapiens'. Examples: 'mus_musculus', 'danio_rerio', 'drosophila_melanogaster'. Default `homo_sapiens` |

??? note "Upstream description"

    Get the DNA nucleotide sequence for a specific genomic region from Ensembl. Retrieves the reference genome sequence for given chromosome coordinates. Useful for examining regulatory regions, intergenic sequences, or any arbitrary genomic interval. Example: region '17:7668421-7668520' in homo_sapiens returns 100 bp of the TP53 locus. Supports any Ensembl species and coordinate strand.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'EnsemblSeq_get_region_sequence',
     'arguments': {'region': '...'}},
)
```

## `IMGT_get_germline_gene_fasta` {#IMGT_get_germline_gene_fasta}

Retrieves IMGT/GENE-DB germline IG/TR V/D/J-region reference sequences in FASTA for a gene type and species.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** none
- **Keywords:** imgt, germline, v(d)j, immunoglobulin, tcr, fasta
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | Default `get_germline_fasta` |
| `gene_type` | string | no | IG/TR gene type/locus: IGHV, IGHD, IGHJ, IGKV, IGLV, TRAV, TRBV, TRBD, TRBJ, etc. (alias: gene, query) |
| `gene` | string | no | Alias for gene_type. |
| `species` | string | no | Species name (default: Homo sapiens). Examples: 'Mus musculus', 'Homo sapiens'. Default `Homo sapiens` |
| `label` | string | no | IMGT GENElect label controlling the region set (default '7.2' = all V/D/J-REGION reference alleles). Default `7.2` |

??? note "Upstream description"

    Retrieve IMGT germline immunoglobulin (IG) / T-cell-receptor (TR) gene reference sequences in FASTA from IMGT/GENE-DB for a given gene type and species. Returns labeled FASTA records (V/D/J-REGION) e.g. '&gt;AB019441\|IGHV(II)-1-1*01\|Homo sapiens\|...'. Use for antibody/TCR repertoire reference sets and germline alignment. gene_type is an IG/TR locus (IGHV, IGHD, IGHJ, IGKV, IGLV, TRAV, TRBV, TRBD, TRBJ, ...). No authentication required.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'IMGT_get_germline_gene_fasta', 'arguments': {}},
)
```

## `IMGT_get_sequence` {#IMGT_get_sequence}

Gets an immunoglobulin/T cell receptor sequence from IMGT/LIGM-DB by accession in FASTA or EMBL format.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `accession`
- **Keywords:** imgt, immunoglobulin, embl, fasta, dbfetch
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | Default `get_sequence` |
| `accession` | string | yes | IMGT/LIGM-DB accession or EMBL/GenBank accession |
| `format` | string | no | Output format: fasta or embl (default: fasta) One of `fasta`, `embl` Default `fasta` |

??? note "Upstream description"

    Get immunoglobulin/T cell receptor sequence from IMGT by accession number. IMGT is the reference for IG/TR/MHC gene sequences. Returns sequence in FASTA or EMBL format. Uses EBI DBFetch. No authentication required.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'IMGT_get_sequence', 'arguments': {'accession': '...'}},
)
```

## `NCBI_fetch_accessions` {#NCBI_fetch_accessions}

Converts NCBI Nucleotide UIDs to GenBank/EMBL/RefSeq accession numbers.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `uids`
- **Keywords:** genbank, uid, accession, refseq
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | any | no | Operation type (fixed: fetch_accession) Default `fetch_accession` |
| `uids` | any | yes | GenBank UID(s) from NCBI_search_nucleotide. Can be single UID or array of UIDs. |

??? note "Upstream description"

    Convert GenBank UIDs to accession numbers (U00096, NC_045512, etc.). Takes UIDs from NCBI_search_nucleotide and returns actual accession numbers that can be used with database tools. Returns both EMBL/GenBank accessions (U*, M*, AJ*) that work with ENA Browser tools, and RefSeq accessions (NC_*, NM_*) for NCBI tools. Use this as second step after NCBI_search_nucleotide to get usable accession numbers.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'NCBI_fetch_accessions', 'arguments': {'uids': '...'}},
)
```

## `NCBI_get_sequence` {#NCBI_get_sequence}

Retrieves nucleotide sequences from NCBI by accession in FASTA, GenBank or EMBL format.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `accession`
- **Keywords:** fasta, genbank, nucleotide sequence, refseq
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | any | no | Operation type (fixed: fetch_sequence) Default `fetch_sequence` |
| `accession` | string | yes | NCBI accession number (e.g., 'U00096', 'NC_045512', 'NM_000546'). Works with GenBank, EMBL, and RefSeq accessions. |
| `format` | string | no | Sequence format: 'fasta' for FASTA, 'gb' for GenBank, 'embl' for EMBL, 'gp' for GenPept protein. Default: fasta One of `fasta`, `gb`, `embl`, `gp` Default `fasta` |

??? note "Upstream description"

    Retrieve DNA/RNA sequence data from NCBI by accession number. Returns sequences in specified format (fasta, gb, embl). Works with both EMBL/GenBank accessions (U00096, M15390) and RefSeq accessions (NC_045512, NM_000546). Use this to get actual sequence data after finding accessions with NCBI_search_nucleotide. Returns full sequence with annotations depending on format.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'NCBI_get_sequence', 'arguments': {'accession': '...'}},
)
```

## `NCBI_search_nucleotide` {#NCBI_search_nucleotide}

Searches NCBI Nucleotide (GenBank/EMBL/RefSeq) by organism, gene, strain, keywords or sequence type.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** none
- **Keywords:** genbank search, complete genome, mrna, nucleotide
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | any | no | Operation type (fixed: search) Default `search` |
| `organism` | string | no | Organism name to search (e.g., 'Escherichia coli', 'Homo sapiens'). Searches in \[Organism\] field. |
| `gene` | string | no | Gene name to search (e.g., 'BRCA1', 'TP53'). Searches in \[Gene\] field. |
| `strain` | string | no | Strain name to filter (e.g., 'K-12', 'MG1655'). Searches in \[Strain\] field. |
| `keywords` | string | no | Keywords to search in titles (e.g., 'complete genome', 'mitochondrial'). Searches in \[Title\] field. |
| `seq_type` | string | no | Sequence type filter: 'complete_genome' for complete genomes, 'mrna' for mRNA sequences, 'refseq' for RefSeq curated sequences One of `complete_genome`, `mrna`, `refseq` |
| `query` | string | no | Free-form search query if not using structured filters. Use NCBI query syntax (e.g., '(Escherichia coli\[Organism\]) AND (K-12\[Strain\])') |
| `limit` | integer | no | Maximum number of results to return (default: 20, max: 100) Default `20` |
| `sort` | string | no | Sort order for results (default: relevance) One of `relevance`, `pub_date`, `title` Default `relevance` |

??? note "Upstream description"

    Search NCBI Nucleotide database (GenBank/EMBL/RefSeq) by organism, gene name, or keywords. Returns GenBank UIDs which can be used to fetch accession numbers. Supports filtering by organism, gene, strain, sequence type (complete genome, mRNA, RefSeq). Use this to discover DNA/RNA sequences for specific organisms or genes. Returns both EMBL/GenBank accessions (U00096, M15390) and RefSeq accessions (NC_*, NM_*). Follow up with NCBI_fetch_accessions to get actual accession numbers.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'NCBI_search_nucleotide', 'arguments': {}},
)
```

## `Pathoplexus_count_sequences` {#Pathoplexus_count_sequences}

Counts and aggregates open pathogen genome sequences in Pathoplexus by country, lineage or other fields.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `organism`
- **Keywords:** pathogen genome, genomic surveillance, lineage, lapis
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `organism` | string | yes | Pathoplexus organism slug: 'west-nile', 'ebola-zaire', 'ebola-sudan', 'cchf', 'mpox'. |
| `country` | string \| null | no | Filter by country (geoLocCountry), e.g. 'USA'. |
| `lineage` | string \| null | no | Filter by lineage. |
| `group_by` | string \| null | no | Metadata field to group counts by, e.g. 'geoLocCountry', 'lineage', 'sampleCollectionDate'. |

??? note "Upstream description"

    Count and aggregate open pathogen genome sequences in Pathoplexus (LAPIS) for genomic-epidemiology surveillance of emerging pathogens. Returns sequence counts for an organism, optionally filtered by country and lineage and grouped by a metadata field (e.g. geoLocCountry, lineage, sampleCollectionDate) to see geographic/temporal distribution. Known organisms: 'west-nile', 'ebola-zaire', 'ebola-sudan', 'cchf' (Crimean-Congo hemorrhagic fever), 'mpox'. No API key required.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Pathoplexus_count_sequences',
     'arguments': {'organism': '...'}},
)
```

## `Pathoplexus_get_sequences_fasta` {#Pathoplexus_get_sequences_fasta}

Downloads unaligned or aligned nucleotide or amino-acid FASTA sequences for a Pathoplexus organism.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `organism`
- **Keywords:** fasta, pathogen genome, alignment, nucleotide
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `organism` | string | yes | Pathoplexus organism slug: 'west-nile', 'ebola-zaire', 'ebola-sudan', 'cchf', 'mpox'. |
| `country` | string \| null | no | Filter by country (geoLocCountry), e.g. 'USA'. |
| `lineage` | string \| null | no | Filter by lineage. |
| `sequence_type` | string \| null | no | 'nucleotide' (default) or 'aminoAcid'. |
| `aligned` | boolean \| null | no | Return aligned sequences instead of unaligned (default false). |
| `limit` | integer \| null | no | Number of sequences to download (default 1, max 100). |

??? note "Upstream description"

    Download the actual genome sequences (FASTA) from Pathoplexus/LAPIS for a surveilled pathogen. Returns unaligned (default) or aligned nucleotide or amino-acid sequences as FASTA text, optionally filtered by country and lineage. This is the only way in ToolUniverse to retrieve the nucleotide/protein sequence of Pathoplexus genomes — the other tools expose only counts and mutations. Known organisms: 'west-nile', 'ebola-zaire', 'ebola-sudan', 'cchf', 'mpox'. No API key required.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Pathoplexus_get_sequences_fasta',
     'arguments': {'organism': '...'}},
)
```

## `phykit_batch_analysis` {#phykit_batch_analysis}

Runs PhyKIT functions (treeness, saturation, parsimony-informative sites) over a directory of alignments or trees.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** none
- **Keywords:** phylogenomics, alignment, parsimony informative, treeness
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | single: one file, batch: all files in directory, gap_percentage: compute gap % across alignments One of `single`, `batch`, `gap_percentage` Default `single` |
| `function` | string | no | PhyKIT function to run. 'treeness_over_rcv' (treeness/RCV) is distinct from 'treeness' and needs both an alignment and a tree. One of `treeness`, `saturation`, `dvmc`, `long_branch_score`, `total_tree_length`, `parsimony_informative`, `treeness_over_rcv` |
| `file` | string \| null | no | Path to single tree/alignment file (for operation=single) |
| `directory` | string \| null | no | Directory containing tree/alignment files (for operation=batch or gap_percentage) |
| `extension` | string | no | File extension to match (default: .treefile for trees, .fa for alignments) |
| `tree_directory` | string \| null | no | Directory with paired tree files (for saturation function) |
| `tree_extension` | string | no | Tree file extension (default: .treefile) |
| `per_tree_stat` | string | no | How to summarize per-tree values for long_branch_score (default: mean) One of `mean`, `median` |

??? note "Upstream description"

    Run PhyKIT phylogenetics functions (treeness, saturation, dvmc, long_branch_score, total_tree_length, parsimony_informative) on a directory of tree/alignment files. Returns summary statistics (mean, median, min, max). For long_branch_score, aggregates per-taxon values within each tree before cross-tree summary. Use for batch phylogenetic analyses comparing gene trees across organism groups.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'phykit_batch_analysis', 'arguments': {}},
)
```

## `Rfam_get_alignment` {#Rfam_get_alignment}

Returns the Rfam seed multiple sequence alignment in Stockholm, Pfam or FASTA format.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `family_id`
- **Keywords:** stockholm, msa, ncrna, seed alignment
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | Operation type One of `get_alignment` Default `get_alignment` |
| `family_id` | string | yes | Required: Rfam accession or family name |
| `format` | string | no | Alignment format: stockholm (structured), pfam (single-line), fasta (gapped), fastau (ungapped) One of `stockholm`, `pfam`, `fasta`, `fastau` Default `stockholm` |
| `gzip` | boolean | no | Compress output with gzip Default `False` |

??? note "Upstream description"

    Get RNA family seed alignment in multiple formats. Returns curated multiple sequence alignment used to build covariance model. Formats: stockholm (default, with structure annotation), pfam (stockholm single-line), fasta (gapped), fastau (ungapped). Optional gzip compression for large alignments. Use for: sequence analysis, structure prediction, phylogenetic studies, homology detection. Alignments show conserved positions and base pairing.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Rfam_get_alignment', 'arguments': {'family_id': '...'}},
)
```

## `Rfam_get_covariance_model` {#Rfam_get_covariance_model}

Returns the Infernal covariance model for an Rfam family.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `family_id`
- **Keywords:** infernal, cmsearch, covariance model, ncrna
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | Operation type One of `get_covariance_model` Default `get_covariance_model` |
| `family_id` | string | yes | Required: Rfam accession or family name |

??? note "Upstream description"

    Get Infernal covariance model (CM) for RNA family. Returns CM file in Infernal format, which captures both sequence and structure conservation. Can be used with Infernal tools (cmsearch, cmalign, cmbuild) for RNA sequence analysis. Use for: homology searches, structural alignment, custom RNA annotation, sequence classification. CMs are more sensitive than profile HMMs for structured RNAs.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Rfam_get_covariance_model',
     'arguments': {'family_id': '...'}},
)
```

## `Rfam_get_tree_data` {#Rfam_get_tree_data}

Returns the NHX phylogenetic tree for an Rfam family seed alignment.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `family_id`
- **Keywords:** nhx, phylogeny, seed alignment
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | Operation type One of `get_tree_data` Default `get_tree_data` |
| `family_id` | string | yes | Required: Rfam accession or family name |

??? note "Upstream description"

    Get phylogenetic tree for RNA family in NHX format (New Hampshire eXtended). Returns tree structure based on seed alignment showing evolutionary relationships. Tree includes species names and sequence information. Can be visualized with tree viewers (FigTree, iTOL, ETE). Use for: phylogenetic analysis, understanding RNA family evolution, species distribution analysis, taxonomic studies. NHX format includes extended annotations beyond standard Newick.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Rfam_get_tree_data', 'arguments': {'family_id': '...'}},
)
```

## `Rfam_search_sequence` {#Rfam_search_sequence}

Searches an RNA sequence against all Rfam families with Infernal cmscan.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `sequence`
- **Keywords:** cmscan, infernal, rna sequence, homology
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | Operation type One of `search_sequence` Default `search_sequence` |
| `sequence` | string | yes | Required: RNA sequence to search (can be multi-line, no FASTA header needed) |
| `max_wait_seconds` | integer | no | Maximum time to wait for results (default: 120s). Job continues running if timeout reached. Default `120` |

??? note "Upstream description"

    Search RNA sequence against all Rfam families using Infernal cmscan. Submits sequence to batch search server, polls for completion, and returns hits with family matches, positions, scores, and E-values. Typically completes in 4-5 seconds for short sequences. Results available for 1 week. Use for: annotating novel RNA sequences, identifying RNA family membership, finding conserved elements, structural RNA discovery. Returns family IDs, alignment details, and significance scores.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Rfam_search_sequence', 'arguments': {'sequence': '...'}},
)
```

## `RNAcentral_get_sequence` {#RNAcentral_get_sequence}

Returns the FASTA nucleotide sequence of an RNAcentral ncRNA.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `urs_id`
- **Keywords:** fasta, rna sequence, urs
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | Fixed to 'sequence' for this tool. One of `sequence` Default `sequence` |
| `urs_id` | string | yes | RNAcentral URS identifier, e.g. 'URS00003B7674'. A trailing '_taxid' suffix is accepted and stripped. |

??? note "Upstream description"

    Get the nucleotide sequence (FASTA) of a non-coding RNA from RNAcentral. Provide an RNAcentral URS identifier (e.g. 'URS00003B7674'). Returns the FASTA header, the full RNA sequence string, and its length.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'RNAcentral_get_sequence', 'arguments': {'urs_id': '...'}},
)
```

## `Sequence_count_residues` {#Sequence_count_residues}

Counts occurrences and positions of a given nucleotide or amino acid residue in a sequence string or a UniProt-fetched protein sequence.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `residue`
- **Keywords:** nucleotide count, residue positions, composition, amino acid, cysteine
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | Operation type One of `count_residues` Default `count_residues` |
| `sequence` | string | no | Sequence string (DNA/RNA/protein). Optional if uniprot_id is provided. |
| `uniprot_id` | string | no | UniProt accession to fetch sequence from (e.g., 'P00533' for EGFR). Optional if sequence is provided. |
| `residue` | string | yes | Single character residue to count (e.g., 'C' for cysteine) |

??? note "Upstream description"

    Count occurrences of a specific amino acid or nucleotide residue in a sequence. Accepts either a raw sequence string or a UniProt accession (fetches protein sequence from UniProt REST API). Reports count, fraction, percentage, and 1-based positions. Useful for cysteine counting (disulfide bond estimation), charged residue analysis, or nucleotide composition.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Sequence_count_residues',
     'arguments': {'residue': '...'}},
)
```

## `Sequence_dn_ds` {#Sequence_dn_ds}

Computes dN/dS (Ka/Ks) between two codon-aligned coding sequences with Nei-Gojobori and Jukes-Cantor correction.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** none
- **Keywords:** dn/ds, ka/ks, selection, nei-gojobori, coding sequence
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `seq1` | string \| null | no | First coding sequence (in-frame, codon-aligned to seq2). DNA/RNA letters. |
| `seq2` | string \| null | no | Second coding sequence (same length/frame as seq1). |
| `fasta1_path` | string \| null | no | Alternative to seq1: path to a single-record FASTA |
| `fasta2_path` | string \| null | no | Alternative to seq2: path to a single-record FASTA |

??? note "Upstream description"

    Compute dN/dS (Ka/Ks) between two coding sequences via the Nei-Gojobori (1986) estimator with Jukes-Cantor correction — the standard test for selection: dN/dS &gt; 1 = positive/diversifying, &lt;&lt; 1 = purifying, ~ 1 = near-neutral. Give two in-frame, codon-aligned coding sequences inline ('seq1'/'seq2') or as single-record FASTA files. Returns dN, dS, the ratio, per-site counts, and an interpretation. Pure Python, no API key. Works for any organisms/genes.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Sequence_dn_ds', 'arguments': {}},
)
```

## `Sequence_gc_content` {#Sequence_gc_content}

Calculates GC count, fraction and percentage and full nucleotide composition of a DNA or RNA sequence.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `sequence`
- **Keywords:** gc percent, nucleotide composition, base composition, primer design, dna
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | Operation type One of `gc_content` Default `gc_content` |
| `sequence` | string | yes | DNA or RNA sequence string |

??? note "Upstream description"

    Calculate GC content of a DNA or RNA sequence. Reports GC count, GC fraction, GC percentage, and full nucleotide composition. Ambiguous bases (N) are excluded from the denominator. Useful for primer design, genome characterization, and codon usage analysis.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Sequence_gc_content', 'arguments': {'sequence': '...'}},
)
```

## `Sequence_reverse_complement` {#Sequence_reverse_complement}

Returns the reverse complement of a DNA sequence (A/T/C/G/N).

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `sequence`
- **Keywords:** complementary strand, strand orientation, primer design, probe design, dna
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | Operation type One of `reverse_complement` Default `reverse_complement` |
| `sequence` | string | yes | DNA sequence (A/T/C/G/N only) |

??? note "Upstream description"

    Return the reverse complement of a DNA sequence. Supports A, T, C, G, and N bases. Useful for primer design, probe design, and strand orientation analysis.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Sequence_reverse_complement',
     'arguments': {'sequence': '...'}},
)
```

## `Sequence_stats` {#Sequence_stats}

Computes length, composition and type-specific metrics (GC% for nucleotides, molecular weight for proteins) for a sequence, auto-detecting its type.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** none
- **Keywords:** sequence length, gc content, composition, molecular weight, dna rna protein
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | Operation type One of `stats` Default `stats` |
| `sequence` | string | no | Sequence string. Optional if uniprot_id is provided. |
| `uniprot_id` | string | no | UniProt accession to fetch protein sequence from. Optional if sequence is provided. |

??? note "Upstream description"

    Compute basic statistics for a DNA, RNA, or protein sequence. Auto-detects sequence type, reports length, full composition, and type-specific metrics (GC% for nucleotides, estimated monoisotopic molecular weight for proteins -- for the conventional average mass reported by ExPASy ProtParam/UniProt, use ProtParam_calculate instead). Accepts sequence string or UniProt accession.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Sequence_stats', 'arguments': {}},
)
```

## `SGD_get_sequence` {#SGD_get_sequence}

Gets genomic DNA, coding DNA and protein sequences plus genomic coordinates for a yeast gene in reference strain S288C from SGD.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `sgd_id`
- **Keywords:** yeast, genomic dna, coding sequence, s288c, genomic coordinates
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sgd_id` | string | yes | SGD locus identifier. Examples: 'S000000364' (CDC28), 'S000005564' (ACT1). |

??? note "Upstream description"

    Get genomic, coding, and protein sequences plus genomic coordinates for a yeast gene from SGD (reference strain S288C). Returns the genomic DNA start/end/strand and full residue strings for genomic_dna, coding_dna, and protein. Example: S000000364 (CDC28) returns genomic_dna start=560078 end=560974 with residues 'ATGAGCGGTGAATTAGCAAAT...'.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'SGD_get_sequence', 'arguments': {'sgd_id': '...'}},
)
```

## `UCSC_get_sequence` {#UCSC_get_sequence}

Gets DNA sequence for a genomic region of a UCSC genome assembly given a locus or chrom/start/end.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `genome`
- **Keywords:** genomic dna, reference genome, hg38, locus, chromosome region
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `genome` | string | yes | UCSC genome assembly identifier. Examples: 'hg38' (human), 'mm39' (mouse), 'dm6' (fly), 'sacCer3' (yeast). |
| `chrom` | string | no | Chromosome name. Examples: 'chr1', 'chr17', 'chrX', 'chrM'. |
| `start` | integer | no | Start position. Interpreted per 'coordinate_system' (default 0-based, inclusive). Example: 7668421. |
| `end` | integer | no | End position. Exclusive when coordinate_system is 0-based, inclusive when 1-based. Span must be &lt;= 100000 bp. Example: 7668521. |
| `region` | string | no | Genomic region written the usual way, 1-based INCLUSIVE on both ends, e.g. 'chr14:89000000-89000100' (101 bp). This is the form shown by the UCSC browser, IGV and samtools, so you can paste a locus verbatim. PREFERRED over chrom/start/end: it removes the 0-based/1-based conversion that otherwise silently drops the first base. |
| `coordinate_system` | string | no | Which convention chrom/start/end are given in. '0-based' (default) is half-open, matching the UCSC API: length = end - start. '1-based' is inclusive, matching a written locus: length = end - start + 1. Ignored when 'region' is supplied. One of `0-based`, `1-based` Default `0-based` |

??? note "Upstream description"

    Get DNA sequence from the UCSC Genome Browser for a genomic region. Pass 'region' as a written locus (e.g. 'chr14:89000000-89000100', 1-based inclusive, 101 bp) to paste coordinates verbatim, or give chrom/start/end with 'coordinate_system'. Returns the nucleotide sequence plus the resolved region in both conventions. Useful for primer design, sequence analysis and regulatory element inspection.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'UCSC_get_sequence', 'arguments': {'genome': '...'}},
)
```
