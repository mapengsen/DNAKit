# Optional sequence functions: DNA and gene-level wrappers over Biopython and scikit-bio

This page describes 26 functions in `dnakit.sequence` and `dnakit.sequence`: 20 from Biopython and 6 from scikit-bio. They call DNA- and gene-level features of the upstream packages on explicit strings, mappings or local files. None of them accesses the network or runs an external program. In this version, import them from those two modules, for example `from dnakit.sequence import read_gaf_annotations`.

The dependencies are optional and are imported only when the matching function is called. If a package is missing, the function raises an error that includes the install command. Other upstream functions, such as protein-structure, ecology, remote-retrieval and duplicate implementations, are not wrapped; their disposition is recorded in the ledger.
## Example data

The examples below use these small files. They only demonstrate the interface and their values are not real data. Inputs for the other examples are written in the code.

**Example file:** `primersearch.txt`

```text
Primer name P1
Amplimer 1
	Sequence: ATGCGATTTGCAAGCTTAGC
	Amplimer length: 42 bp
	Forward primer: 0 mismatches
	Reverse primer: 0 mismatches
```

**Example file:** `contig.ace`

```text
AS 1 2

CO Contig1 12 2 4 U
ACGTACGTACGT

BQ
30 30 30 30 30 30 30 30 30 30 30 30

AF read1 U 1
AF read2 U 2

RD read1 12 0 0
ACGTACGTACGT

QA 1 12 1 12
DS CHROMAT_FILE read1 PHD_FILE read1.phd.1 TIME Thu Jan 1 00:00:00 2024

RD read2 12 0 0
ACGTACGTACGT

QA 1 12 1 12
DS CHROMAT_FILE read2 PHD_FILE read2.phd.1 TIME Thu Jan 1 00:00:00 2024

```

**Example file:** `read.phd`

```text
BEGIN_SEQUENCE read1
BEGIN_COMMENT
CHROMAT_FILE: read1.ab1
PHRED_VERSION: 0.980528.c
CALL_METHOD: phred
END_COMMENT
BEGIN_DNA
A 10 1
C 20 4
G 30 2
T 40 3
END_DNA
END_SEQUENCE
```

**Example file:** `pop.genepop`

```text
Example GenePop file
Locus1
Locus2
Pop
ind1 , 0101 0203
ind2 , 0102 0000
Pop
ind3 , 0201 0303
ind4 , 0000 0303
```

**Example file:** `annotations.gaf`

```text
!gaf-version: 2.2
UniProtKB	P12345	GENEA		GO:0003700	PMID:1234567	IDA		F	TF A	GENEA|TFA	protein	taxon:9606	20240101	UniProt
UniProtKB	P67890	GENEB	involved_in	GO:0006355	PMID:7654321	IEA		P	Protein B	GENEB	protein	taxon:9606	20240102	UniProt
```

**Example file:** `gene.txt`

```text
ENTRY       hsa:7157                    CDS       T01001
SYMBOL      TP53
NAME        (RefSeq) tumor protein p53
ORTHOLOGY   K04442  tumor protein p53
ORGANISM    hsa  Homo sapiens (human)
POSITION    17:7661779..7687538
DBLINKS     NCBI-GeneID: 7157
///
```

**Example file:** `alignments.psl`

```text
psLayout version 3

match	mis-	rep.	N's	Q gap	Q gap	T gap	T gap	strand	Q name	Q size	Q start	Q end	T name	T size	T start	T end	block count	blockSizes	qStarts	tStarts
     	match	match	   	count	bases	count	bases	      									
---------------------------------------------------------------------------------------------------------------
20	0	0	0	0	0	0	0	+	read1	20	0	20	chr1	1000	100	120	1	20,	0,	100,
```

**Example file:** `primer3.txt`

```text
# PRIMER3 test
#                      Start  Len   Tm     GC%   Sequence
     PRODUCT SIZE: 60
     FORWARD PRIMER    1  20  60.00  50.00  ACGTACGTACGTACGTACGT
     REVERSE PRIMER   60  20  60.00  50.00  TGCATGCATGCATGCATGCA
```

**Example file:** `reads.fq`

```text
@r1
ACGT
+
IIII
```

## 1) `dnakit.sequence.gc_content_by_codon_position()`

- **Purpose:** Calculate total G+C content and the G+C percentage at codon positions 1, 2 and 3.
- **Method:** Call Biopython `Bio.SeqUtils.GC123`. Whitespace is removed and letters are uppercased. IUPAC ambiguity symbols count toward the length but not as G or C.
- **API:** `dnakit.sequence.gc_content_by_codon_position(sequence[required], max_length[optional], verbose[optional])`.
- **Input:** One DNA string no longer than `max_length` (default 1,000,000).
- **Output:** `gc_percent` (with `total`, `first`, `second` and `third`, as percentages), `sequence_length` and `provenance`.
- **Example code:**

```python
from dnakit.sequence import gc_content_by_codon_position

result = gc_content_by_codon_position("ATGCGTACGTTAGC", verbose=False)
print(result["gc_percent"])
```

- **Example output:**

```text
{'total': 50.0, 'first': 40.0, 'second': 60.0, 'third': 50.0}
```

## 2) `dnakit.sequence.gc_skew_windows()`

- **Purpose:** Compute the GC skew (G-C)/(G+C) for consecutive windows, for example to look for replication origins or strand bias.
- **Method:** Call `Bio.SeqUtils.GC_skew`. The last window can be shorter. A window with no G or C has skew 0.0. Window coordinates are 1-based and inclusive.
- **API:** `dnakit.sequence.gc_skew_windows(sequence[required], window[optional], max_length[optional], verbose[optional])`.
- **Input:** One DNA string. `window` is a positive integer no longer than the sequence.
- **Output:** `window`, `window_count` and `windows`. Each window has `start`, `end` and `skew`.
- **Example code:**

```python
from dnakit.sequence import gc_skew_windows

result = gc_skew_windows("GGGGCCCCAAAATTTTGGGG", window=5, verbose=False)
for item in result["windows"]:
    print(item["start"], item["end"], item["skew"])
```

- **Example output:**

```text
1 5 0.6
6 10 -1.0
11 15 0.0
16 20 1.0
```

## 3) `dnakit.sequence.local_composition_complexity()`

- **Purpose:** Compute the local composition complexity (LCC) of every sliding window. Lower values mean the window is dominated by fewer bases.
- **Method:** Call `Bio.SeqUtils.lcc.lcc_mult` with window length `window`. Only A, C, G and T are accepted. `values[i]` describes the window that starts at position i + 1.
- **API:** `dnakit.sequence.local_composition_complexity(sequence[required], window[optional], max_length[optional], verbose[optional])`.
- **Input:** One ACGT string. `window` must not exceed the sequence length.
- **Output:** `value_count`, `values`, `minimum`, `maximum`, `mean` and `provenance`.
- **Example code:**

```python
from dnakit.sequence import local_composition_complexity

result = local_composition_complexity("ACGTACGTAAAAAAAACGT", window=5, verbose=False)
print(result["value_count"], [round(value, 4) for value in result["values"][:4]])
```

- **Example output:**

```text
15 [0.961, 0.961, 0.961, 0.961]
```

## 4) `dnakit.sequence.ambiguous_codon_amino_acids()`

- **Purpose:** List every amino acid that one IUPAC-ambiguous codon can encode. A codon that can only encode a stop is reported separately.
- **Method:** Call `Bio.Data.CodonTable.list_possible_proteins` with an NCBI genetic code (table 1 by default). A codon that mixes amino acids and stops raises an error.
- **API:** `dnakit.sequence.ambiguous_codon_amino_acids(codon[required], genetic_code[optional], verbose[optional])`.
- **Input:** A three-symbol IUPAC DNA string. `genetic_code` is an integer from 1 to 33.
- **Output:** `amino_acids` (sorted), `amino_acid_count`, `stop_only` and `provenance`.
- **Example code:**

```python
from dnakit.sequence import ambiguous_codon_amino_acids

for codon in ["GCN", "TTY", "TAR"]:
    result = ambiguous_codon_amino_acids(codon, verbose=False)
    print(codon, result["amino_acids"], result["stop_only"])
```

- **Example output:**

```text
GCN ['A'] False
TTY ['F'] False
TAR [] True
```

## 5) `dnakit.sequence.expand_ambiguous_codons()`

- **Purpose:** Derive the IUPAC-ambiguous codons that are fully covered by a set of concrete codons, such as TAR from TAA and TAG.
- **Method:** Call `Bio.Data.CodonTable.list_ambiguous_codons`. A codon is added only when every concrete codon it covers is present, so TRR is not added.
- **API:** `dnakit.sequence.expand_ambiguous_codons(codons[required], max_codons[optional], verbose[optional])`.
- **Input:** A sequence of three-base ACGT codons. Duplicates are removed.
- **Output:** `input_codons`, `added_codons`, `expanded_codons` (all sorted) and `provenance`.
- **Example code:**

```python
from dnakit.sequence import expand_ambiguous_codons

result = expand_ambiguous_codons(["TGA", "TAA", "TAG"], verbose=False)
print(result["added_codons"])
print(result["expanded_codons"])
```

- **Example output:**

```text
['TAR', 'TRA']
['TAA', 'TAG', 'TAR', 'TGA', 'TRA']
```

## 6) `dnakit.sequence.sequence_checksum()`

- **Purpose:** Compute a SEGUID, CRC32, CRC64 or GCG checksum of one sequence, for recording and matching sequence identity.
- **Method:** Call the matching function in `Bio.SeqUtils.CheckSum`. Whitespace is removed and letters are uppercased first. SEGUID and CRC64 return strings; CRC32 and GCG return integers.
- **API:** `dnakit.sequence.sequence_checksum(sequence[required], method[optional], max_length[optional], verbose[optional])`.
- **Input:** One IUPAC DNA string. `method` is `seguid`, `crc32`, `crc64` or `gcg`.
- **Output:** `method`, `sequence_length`, `value` and `provenance`.
- **Example code:**

```python
from dnakit.sequence import sequence_checksum

for method in ["seguid", "crc32", "crc64", "gcg"]:
    print(method, sequence_checksum("ATGCGT", method=method, verbose=False)["value"])
```

- **Example output:**

```text
seguid Yrhkea/OcXcnrtHsYwY41PbWdEk
crc32 527699539
crc64 CRC-71A87EB861ADB000
gcg 1573
```

## 7) `dnakit.sequence.mcdonald_kreitman_test()`

- **Purpose:** Run the McDonald-Kreitman neutrality test, which compares within-species polymorphism with between-species fixed differences.
- **Method:** Call `Bio.Align.analysis.mktest`. Synonymous and nonsynonymous changes are counted per codon, and the G test gives the p-value. Every cell of the 2x2 table must be non-zero, otherwise a configuration error is raised.
- **API:** `dnakit.sequence.mcdonald_kreitman_test(sequences[required], species[required], genetic_code[optional], max_length[optional], verbose[optional])`.
- **Input:** `sequences` maps IDs to coding sequences, and `species` maps the same IDs to species labels. Sequences must be equal length, a multiple of 3, ACGT only, from at least two species.
- **Output:** `p_value`, `sequence_count`, `codon_count`, `species`, `genetic_code` and `provenance`.
- **Example code:**

```python
from dnakit.sequence import mcdonald_kreitman_test

sequences = {"a1": "GCTGAAAAATTT", "a2": "GCCGAAAAATTA",
             "b1": "GCTGAGCAATTT", "b2": "GCTGAGCAATTT"}
species = {"a1": "A", "a2": "A", "b1": "B", "b2": "B"}
result = mcdonald_kreitman_test(sequences, species, verbose=False)
print(result["species"], result["codon_count"], round(result["p_value"], 4))
```

- **Example output:**

```text
['A', 'B'] 4 1.0
```

## 8) `dnakit.sequence.pairwise_dn_ds_matrix()`

- **Purpose:** Compute pairwise dN (nonsynonymous rate) and dS (synonymous rate) matrices for a set of coding sequences.
- **Method:** Call `Bio.Align.analysis.calculate_dn_ds_matrix` with NG86 or LWL85. Biopython reports saturation as -1; this function reports it as `null` and counts it in `saturated_pair_count`.
- **API:** `dnakit.sequence.pairwise_dn_ds_matrix(sequences[required], method[optional], genetic_code[optional], max_length[optional], verbose[optional])`.
- **Input:** A mapping from IDs to coding sequences. Sequences must be equal length, a multiple of 3, ACGT only, and at least two.
- **Output:** `ids`, `dn`, `ds` (full symmetric matrices with zero diagonals), `saturated_pair_count` and `provenance`.
- **Example code:**

```python
from dnakit.sequence import pairwise_dn_ds_matrix

sequences = {
    "s1": "ATGGCTGCCAAATTTGGGCTGAAAGATCGTAAGCCTTTCGCA",
    "s2": "ATGGCAGCCAAATTCGGGCTGAAAGATCGCAAGCCTTTCGCA",
    "s3": "ATGGCTGCCAAATTTGGGCTGAAAGATCGTAAGCCTTTCGCT",
}
result = pairwise_dn_ds_matrix(sequences, verbose=False)
print(result["ids"])
print([round(value, 4) for value in result["ds"][0]])
print(result["saturated_pair_count"])
```

- **Example output:**

```text
['s1', 's2', 's3']
[0.0, 0.4197, 0.1156]
0
```

## 9) `dnakit.sequence.build_codon_alignment()`

- **Purpose:** Lay out DNA codons on the column structure of a protein alignment to make a codon alignment.
- **Method:** Call `Bio.codonalign.build`. Gaps in the protein become three-base codon gaps. Each DNA codon must match its amino acid.
- **API:** `dnakit.sequence.build_codon_alignment(protein_alignment[required], nucleotide_sequences[required], genetic_code[optional], max_length[optional], verbose[optional])`.
- **Input:** `protein_alignment` is an equal-length protein alignment (gaps allowed). `nucleotide_sequences` are ungapped ACGT sequences. The IDs must match.
- **Output:** `ids`, `alignment_length`, `sequences` (ID to codon-aligned string) and `provenance`.
- **Example code:**

```python
from dnakit.sequence import build_codon_alignment

protein = {"s1": "MAA-KF", "s2": "MAGKF-", "s3": "MAAKY-"}
nucleotides = {"s1": "ATGGCTGCCAAATTT", "s2": "ATGGCAGGCAAATTT", "s3": "ATGGCTGCCAAATAT"}
result = build_codon_alignment(protein, nucleotides, verbose=False)
print(result["alignment_length"])
for identifier, row in result["sequences"].items():
    print(identifier, row)
```

- **Example output:**

```text
18
s1 ATGGCTGCC---AAATTT
s2 ATGGCAGGCAAATTT---
s3 ATGGCTGCCAAATAT---
```

## 10) `dnakit.sequence.parse_primersearch_output()`

- **Purpose:** Read a saved EMBOSS primersearch report and return the amplicon sequence, length and details for each primer.
- **Method:** Call `Bio.Emboss.PrimerSearch.read`. Only a saved report is parsed; EMBOSS is not run. Primer names are sorted.
- **API:** `dnakit.sequence.parse_primersearch_output(path[required], verbose[optional])`.
- **Input:** The path to a primersearch output file.
- **Output:** `primer_count`, `amplicon_count`, `primers` (primer name to amplicon list) and `provenance`.
- **Example code:**

```python
from dnakit.sequence import parse_primersearch_output

result = parse_primersearch_output("primersearch.txt", verbose=False)
amplicon = result["primers"]["P1"][0]
print(amplicon["sequence"], amplicon["length"])
print(amplicon["details"])
```

- **Example output:**

```text
ATGCGATTTGCAAGCTTAGC 42
['Forward primer: 0 mismatches', 'Reverse primer: 0 mismatches']
```

## 11) `dnakit.sequence.read_ace_assembly()`

- **Purpose:** Read an ACE assembly file into contigs, the reads in each contig and the position of each read on its contig.
- **Method:** Call `Bio.Sequencing.Ace.read`. Quality and alignment clipping positions keep the 1-based coordinates written in the file.
- **API:** `dnakit.sequence.read_ace_assembly(path[required], verbose[optional])`.
- **Input:** The path to an ACE file.
- **Output:** `contig_count`, `read_count`, `contigs` (each with `name`, `length`, `sequence`, `reads` and `placements`) and `provenance`.
- **Example code:**

```python
from dnakit.sequence import read_ace_assembly

result = read_ace_assembly("contig.ace", verbose=False)
contig = result["contigs"][0]
print(contig["name"], contig["length"], result["read_count"])
print(contig["placements"])
```

- **Example output:**

```text
Contig1 12 2
[{'read': 'read1', 'strand': 'U', 'padded_start': 1}, {'read': 'read2', 'strand': 'U', 'padded_start': 2}]
```

## 12) `dnakit.sequence.read_phd_sequence()`

- **Purpose:** Read one PHRED base-call (PHD) record and return its bases, per-base qualities and peak positions.
- **Method:** Call `Bio.Sequencing.Phd.read`. A line without a peak location gives `null`. Empty comment values are dropped.
- **API:** `dnakit.sequence.read_phd_sequence(path[required], verbose[optional])`.
- **Input:** The path to a PHD file that holds one record.
- **Output:** `file_name`, `sequence`, `length`, `qualities`, `peak_locations`, `comments` and `provenance`.
- **Example code:**

```python
from dnakit.sequence import read_phd_sequence

result = read_phd_sequence("read.phd", verbose=False)
print(result["file_name"], result["sequence"], result["qualities"], result["peak_locations"])
```

- **Example output:**

```text
read1 ACGT [10, 20, 30, 40] [1, 4, 2, 3]
```

## 13) `dnakit.sequence.read_genepop_genotypes()`

- **Purpose:** Read a GenePop genotype file into loci, populations and two-allele genotypes for population-genetic summaries.
- **Method:** Call `Bio.PopGen.GenePop.read`. Missing genotypes `0000` become `[null, null]`. Other allele codes keep their file values.
- **API:** `dnakit.sequence.read_genepop_genotypes(path[required], verbose[optional])`.
- **Input:** The path to a GenePop file with two alleles per locus.
- **Output:** `loci`, `marker_length`, `population_count`, `individual_count`, `populations` (with individual names and genotypes) and `provenance`.
- **Example code:**

```python
from dnakit.sequence import read_genepop_genotypes

result = read_genepop_genotypes("pop.genepop", verbose=False)
print(result["loci"], result["population_count"], result["individual_count"])
print(result["populations"][0]["individuals"][1])
```

- **Example output:**

```text
['Locus1', 'Locus2'] 2 4
{'name': 'ind2', 'genotypes': [[1, 2], [None, None]]}
```

## 14) `dnakit.sequence.read_gaf_annotations()`

- **Purpose:** Read a GAF 1.0 or 2.x gene-annotation file, return every record and count GO terms and aspects.
- **Method:** Call `Bio.UniProt.GOA.gafiterator`. Empty qualifiers, references and synonyms are dropped. If the record count exceeds `max_records`, an error is raised rather than the result being truncated.
- **API:** `dnakit.sequence.read_gaf_annotations(path[required], max_records[optional], verbose[optional])`.
- **Input:** The path to a GAF file. `max_records` defaults to 100,000.
- **Output:** `record_count`, `go_term_counts`, `aspect_counts`, `records` and `provenance`.
- **Example code:**

```python
from dnakit.sequence import read_gaf_annotations

result = read_gaf_annotations("annotations.gaf", verbose=False)
print(result["record_count"], result["go_term_counts"], result["aspect_counts"])
```

- **Example output:**

```text
2 {'GO:0003700': 1, 'GO:0006355': 1} {'F': 1, 'P': 1}
```

## 15) `dnakit.sequence.parse_kegg_gene_entries()`

- **Purpose:** Parse a saved KEGG GENE flat file into names, orthology, organism, chromosome position and database links.
- **Method:** Call `Bio.KEGG.Gene.parse`. Only a local file is read; no KEGG REST request is made.
- **API:** `dnakit.sequence.parse_kegg_gene_entries(path[required], verbose[optional])`.
- **Input:** The path to a KEGG GENE flat file.
- **Output:** `entry_count`, `entries` (each with `entry`, `names`, `orthology`, `organism`, `position`, `dblinks`) and `provenance`.
- **Example code:**

```python
from dnakit.sequence import parse_kegg_gene_entries

result = parse_kegg_gene_entries("gene.txt", verbose=False)
entry = result["entries"][0]
print(entry["entry"], entry["organism"]["code"], entry["orthology"])
```

- **Example output:**

```text
hsa:7157 hsa [{'id': 'K04442', 'name': 'tumor protein p53'}]
```

## 16) `dnakit.sequence.load_substitution_matrix()`

- **Purpose:** Load a substitution matrix bundled with Biopython, including DNA alignment scoring matrices such as NUC.4.4.
- **Method:** Call `Bio.Align.substitution_matrices.load`. The name must appear in Biopython's list. Rows and columns follow the alphabet order.
- **API:** `dnakit.sequence.load_substitution_matrix(name[required], verbose[optional])`.
- **Input:** A matrix name such as `NUC.4.4` or `BLOSUM62`.
- **Output:** `name`, `alphabet`, `size`, `symmetric`, `matrix` and `provenance`.
- **Example code:**

```python
from dnakit.sequence import load_substitution_matrix

result = load_substitution_matrix("NUC.4.4", verbose=False)
print(result["alphabet"], result["symmetric"])
print(result["matrix"][0][:4])
```

- **Example output:**

```text
ATGCSWRYKMBVHDN True
[5.0, -4.0, -4.0, -4.0]
```

## 17) `dnakit.sequence.convert_quality_score()`

- **Purpose:** Convert one quality score between the Solexa and PHRED scales, for older FASTQ data.
- **Method:** Call the formulas in `Bio.SeqIO.QualityIO`. PHRED = 10·log10(10^(Q/10) + 1); the Solexa conversion uses the inverse formula. A PHRED input cannot be negative.
- **API:** `dnakit.sequence.convert_quality_score(value[required], direction[optional], verbose[optional])`.
- **Input:** One finite number. `direction` is `solexa_to_phred` (default) or `phred_to_solexa`.
- **Output:** `direction`, `input`, `output` and `provenance`.
- **Example code:**

```python
from dnakit.sequence import convert_quality_score

print(round(convert_quality_score(10.0, verbose=False)["output"], 4))
print(round(convert_quality_score(10.0, direction="phred_to_solexa", verbose=False)["output"], 4))
```

- **Example output:**

```text
10.4139
9.5424
```

## 18) `dnakit.sequence.read_psl_alignments()`

- **Purpose:** Read PSL pairwise genomic alignments, such as BLAT or minimap2 output, into block coordinates on the target and query.
- **Method:** Call `Bio.Align.parse(format='psl')`. `target_coordinates` and `query_coordinates` are the block boundaries on each sequence; consecutive boundaries with equal steps are ungapped.
- **API:** `dnakit.sequence.read_psl_alignments(path[required], max_alignments[optional], verbose[optional])`.
- **Input:** The path to a PSL file. `max_alignments` defaults to 100,000.
- **Output:** `alignment_count`, `alignments` (each with `target_id`, `query_id`, start and end coordinates and block coordinates) and `provenance`.
- **Example code:**

```python
from dnakit.sequence import read_psl_alignments

result = read_psl_alignments("alignments.psl", verbose=False)
print(result["alignment_count"])
print(result["alignments"][0])
```

- **Example output:**

```text
1
{'target_id': 'chr1', 'query_id': 'read1', 'target_start': 100, 'target_end': 120, 'query_start': 0, 'query_end': 20, 'target_coordinates': [100, 120], 'query_coordinates': [0, 20]}
```

## 19) `dnakit.sequence.bootstrap_consensus_tree()`

- **Purpose:** Resample an aligned DNA set by bootstrap and return the majority-rule consensus tree in Newick format.
- **Method:** Call `Bio.Phylo.Consensus.bootstrap_consensus` with NJ or UPGMA and the identity, blastn or trans distance model. The random seed is fixed during the call and the previous state is restored afterwards.
- **API:** `dnakit.sequence.bootstrap_consensus_tree(sequences[required], replicates[optional], method[optional], distance_model[optional], seed[optional], max_length[optional], verbose[optional])`.
- **Input:** A mapping from IDs to equal-length aligned ACGT sequences, from 3 to 200 sequences. `replicates` is from 1 to 1,000.
- **Output:** `ids`, `method`, `distance_model`, `replicates`, `seed`, `newick` and `provenance`.
- **Example code:**

```python
from dnakit.sequence import bootstrap_consensus_tree

rows = {
    "s1": "ACGTACGTAAGGCCTTACGT",
    "s2": "ACGTACCTAAGGCATTACGA",
    "s3": "ACGTTCGTAAGGCCTAACGT",
    "s4": "ACGTTCGAAAGGCCTAACGA",
}
result = bootstrap_consensus_tree(rows, replicates=5, seed=1, verbose=False)
print(result["newick"])
```

- **Example output:**

```text
(s3:0.025,s4:0.075,(s1:0.025,s2:0.125)100.00:0.08):0;
```

## 20) `dnakit.sequence.parse_primer3_output()`

- **Purpose:** Read one Primer3 design report into primer pairs, product sizes and optional internal oligos.
- **Method:** Call `Bio.Emboss.Primer3.read`. Only a saved report is parsed; Primer3 is not run. The report must contain exactly one target sequence.
- **API:** `dnakit.sequence.parse_primer3_output(path[required], verbose[optional])`.
- **Input:** The path to a Primer3 output file.
- **Output:** `pair_count`, `primer_pairs` (each with product size, forward and reverse primers and optional `internal`), `comments` and `provenance`.
- **Example code:**

```python
from dnakit.sequence import parse_primer3_output

result = parse_primer3_output("primer3.txt", verbose=False)
pair = result["primer_pairs"][0]
print(pair["product_size"], pair["forward"]["sequence"], pair["reverse"]["start"])
```

- **Example output:**

```text
60 ACGTACGTACGTACGTACGT 60
```

## 21) `dnakit.sequence.nucleotide_distance()`

- **Purpose:** Compute an evolutionary distance between two aligned, gap-free DNA sequences: p-distance, JC69, K2P, F81, F84, TN93, LogDet or paralinear.
- **Method:** Call the matching function in `skbio.sequence.distance`. If scikit-bio returns NaN, for example for a saturated distance, the result has `distance: null` and `status: undefined`.
- **API:** `dnakit.sequence.nucleotide_distance(seq1[required], seq2[required], model[optional], verbose[optional])`.
- **Input:** Two equal-length IUPAC DNA strings without gap symbols. `model` names the distance.
- **Output:** `model`, `sequence_length`, `distance`, `status` and `provenance`.
- **Example code:**

```python
from dnakit.sequence import nucleotide_distance

left = "ACGTACGTAAGGCCTTACGT"
right = "ACGTACCTAAGGCATTACGA"
for model in ["pdist", "jc69", "tn93"]:
    result = nucleotide_distance(left, right, model=model, verbose=False)
    print(model, round(result["distance"], 6))
```

- **Example output:**

```text
pdist 0.15
jc69 0.167358
tn93 0.170588
```

## 22) `dnakit.sequence.nucleotide_transition_matrix()`

- **Purpose:** Turn an evolutionary distance into a 4 x 4 nucleotide substitution probability matrix, for simulation or likelihood work.
- **Method:** Call the JC69, K2P, F81, HKY85 or TN93 functions in `skbio.sequence.transition` at float64 precision. The relative `kappa` parameters must lie in (0, 1], which is scikit-bio's convention.
- **API:** `dnakit.sequence.nucleotide_transition_matrix(distance[required], model[optional], kappa[optional], kappa_purine[optional], kappa_pyrimidine[optional], base_frequencies[optional], verbose[optional])`.
- **Input:** A non-negative distance. `base_frequencies` gives positive frequencies for A, C, G and T that sum to 1 (equal by default).
- **Output:** `model`, `distance`, `alphabet` (`ACGT`), `base_frequencies`, `matrix` and `provenance`.
- **Example code:**

```python
from dnakit.sequence import nucleotide_transition_matrix

result = nucleotide_transition_matrix(0.2, verbose=False)
print(result["alphabet"], [round(value, 4) for value in result["matrix"][0]])
```

- **Example output:**

```text
ACGT [0.8244, 0.0585, 0.0585, 0.0585]
```

## 23) `dnakit.sequence.alignment_distance_matrix()`

- **Purpose:** Compute a pairwise distance matrix between the rows of a DNA multiple alignment.
- **Method:** Call `skbio.alignment.align_dists`. Gap symbols `-` are kept and scikit-bio's `shared_by_all` default applies.
- **API:** `dnakit.sequence.alignment_distance_matrix(aligned_sequences[required], metric[optional], verbose[optional])`.
- **Input:** A mapping from IDs to equal-length aligned sequences (gaps allowed), from 2 to 500 rows. `metric` is hamming, pdist, jc69, k2p, f81, f84, tn93, logdet or paralin.
- **Output:** `metric`, `ids`, `matrix` (a symmetric square matrix) and `provenance`.
- **Example code:**

```python
from dnakit.sequence import alignment_distance_matrix

rows = {"x": "ACGTACGTAA", "y": "ACGTTCGTAC", "z": "ACGAACGTAA"}
result = alignment_distance_matrix(rows, metric="hamming", verbose=False)
print(result["ids"])
print(result["matrix"])
```

- **Example output:**

```text
['x', 'y', 'z']
[[0.0, 0.2, 0.1], [0.2, 0.0, 0.3], [0.1, 0.3, 0.0]]
```

## 24) `dnakit.sequence.multiple_nucleotide_alignment()`

- **Purpose:** Progressively align several DNA sequences with scikit-bio, without MAFFT, MUSCLE or Clustal binaries.
- **Method:** Call `skbio.alignment.multi_align_nucl` with configurable match, mismatch and gap scores. Output rows follow the sorted input IDs. Removing the gaps gives back the input sequences.
- **API:** `dnakit.sequence.multiple_nucleotide_alignment(sequences[required], match_score[optional], mismatch_score[optional], gap_cost[optional], free_ends[optional], verbose[optional])`.
- **Input:** A mapping from IDs to unaligned IUPAC DNA sequences, from 2 to 100 sequences.
- **Output:** `ids`, `alignment` (ID to gapped sequence), `alignment_length`, the scoring parameters and `provenance`.
- **Example code:**

```python
from dnakit.sequence import multiple_nucleotide_alignment

rows = {"s1": "ACGTACGTAAGGCC", "s2": "ACGTACCTAAGGCA", "s3": "ACGTTCGTAAGGCCTA"}
result = multiple_nucleotide_alignment(rows, verbose=False)
print(result["alignment_length"])
for identifier, row in result["alignment"].items():
    print(identifier, row)
```

- **Example output:**

```text
16
s1 ACGTACGTAAGGCC--
s2 ACGTACCTAAGGCA--
s3 ACGTTCGTAAGGCCTA
```

## 25) `dnakit.sequence.alignment_score()`

- **Purpose:** Score a DNA alignment with match, mismatch and gap penalties, for comparing alternative alignments.
- **Method:** Call `skbio.alignment.align_score` on a mapping from IDs to equal-length sequences. `free_ends` controls whether terminal gaps are penalized.
- **API:** `dnakit.sequence.alignment_score(aligned_sequences[required], match_score[optional], mismatch_score[optional], gap_cost[optional], free_ends[optional], verbose[optional])`.
- **Input:** A mapping from IDs to equal-length aligned sequences (gaps allowed), at least two.
- **Output:** `ids`, `alignment_length`, `score`, the scoring parameters and `provenance`.
- **Example code:**

```python
from dnakit.sequence import alignment_score

print(alignment_score({"x": "ACGTA", "y": "ACGTA"}, verbose=False)["score"])
print(alignment_score({"x": "AC-TA", "y": "ACGTA"}, verbose=False)["score"])
```

- **Example output:**

```text
5.0
2.0
```

## 26) `dnakit.sequence.sniff_sequence_format()`

- **Purpose:** Detect the format of a local sequence file and the reader options scikit-bio suggests, before reading it.
- **Method:** Call `skbio.io.sniff`. Only the start of the file is read and the file is not converted. An unrecognized file raises a configuration error.
- **API:** `dnakit.sequence.sniff_sequence_format(path[required], verbose[optional])`.
- **Input:** The path to a local sequence file.
- **Output:** `format`, `reader_options`, `source` and `provenance`.
- **Example code:**

```python
from dnakit.sequence import sniff_sequence_format

result = sniff_sequence_format("reads.fq", verbose=False)
print(result["format"], result["reader_options"])
```

- **Example output:**

```text
fastq {}
```

## Installing the optional dependencies

```bash
python -m pip install "dnakit[bio-formats]"   # Biopython
python -m pip install "dnakit[bio-skibio]"    # scikit-bio
```

## Verification

The unit tests are in `others/tests/unit/test_full_G1a.py`. Each function has at least one test that uses local temporary files and calls the real upstream package.
