# H3 family-variant review: Bio, skbio, screed, cutadapt, umi_tools, khmer

This page records the row-by-row review of 108 `family-variant` ledger rows and the 12 wrappers added as a result (10 from Biopython, 1 from scikit-bio, 1 from khmer). The new functions live in `dnakit.sequence`, `dnakit.sequence` and `dnakit.microbiome`; they are exported from the area packages and the examples import from there.

- Offline: every function reads call arguments or local files only, uses no network and writes no file.
- Rule: an alias, re-export or the same object is `duplicate-of`; `family-variant` is kept only when a parameter of an existing wrapper really gives the variant behaviour (the parameter was checked with `inspect.signature` and the result by calling it); the rest became `wrapped`, `deprecated`, `broken-upstream` or `internal`.

## Review result

| Package | family-variant before | duplicate-of | family-variant (verified) | wrapped | Other |
| --- | ---: | ---: | ---: | ---: | --- |
| Bio | 56 | 9 | 24 | 17 | deprecated 3, broken-upstream 1, internal 2 |
| skbio | 18 | 5 | 12 | 1 | none |
| screed | 5 | 3 | 2 | 0 | none |
| cutadapt | 13 | 6 | 2 | 0 | internal 5 |
| umi_tools | 8 | 2 | 4 | 0 | internal 2 |
| khmer | 8 | 0 | 5 | 3 | none |

Main findings (claims that were wrong):

- `Bio.Blast.parse` and `read` read BLAST XML2, while the SearchIO parser behind `parse_search_results` returns an empty result for an XML2 file without an error; `Bio.Blast.write(fmt='XML2')` fails in Biopython 1.88 because `XML2Writer` has no `_write_query_id`.
- `Bio.motifs.pfm.write` is not the JASPAR pfm layout; its text equals `write_motifs(file_format='clusterbuster')`.
- `Bio.UniProt.GOA.gafbyproteiniterator` drops the last protein group (P1,P1,P2,P3 gives only the first two groups) and is marked `broken-upstream`; `as_tab`, `as_fastq_illumina` and `as_fastq_solexa` are deprecated by Biopython.
- `khmer.load_counttable`: `query_saved_countgraph` silently returns 0 on a counttable file (another hash layout), so it is not a variant; the same holds for nodetable and nodegraph, and `inspect_khmer_graph` reads only the header.
- `skbio.sequence.distance.jc69_correct` takes p-distances; `nucleotide_distance(model='jc69')` cannot correct a p-distance measured elsewhere.
- Multi-record readers (several trees, alignment blocks, PHD or Primer3 records, GenBank header fields) cannot be reached by a parameter of the single-record wrappers, so new functions were added.
- `cutadapt.align.naive/py/slow_edit_environment` are test references and `umi_tools.network.recursive_search` keeps state in a function attribute; both are `internal`.

The existing `duplicate-of` claims for `Bio.Phylo.parse` and `Bio.Blast.parse` are outside this group's rows, but the findings above show they are not fully equivalent (several trees, XML2).

## 1) `read_alignment_file`

- **Purpose:** Read pairwise alignments from a BED, Clustal, FASTA, MAF, PHYLIP, PSL or SAM file.
- **Method:** Calls Bio.Align.parse; with `exactly_one=True` it calls Bio.Align.read and raises `ConfigurationError` unless the file holds exactly one alignment. Each alignment reports sequence names, aligned columns, the score (when the format stores one), block coordinates (one row per sequence) and the gapped rows (when the file carries the sequences; PSL, SAM and BED do not, so this is null). The existing `read_psl_alignments` reads PSL only.
- **API:** `read_alignment_file(path[required], file_format[optional], exactly_one[optional], max_alignments[optional], verbose[optional])`
- **Input:** A file path, a format name, whether exactly one alignment is required and an alignment count limit.
- **Output:** A dictionary with `alignment_count`, `alignments` (`sequence_ids`, `aligned_columns`, `score`, `coordinates`, `aligned_rows`) and `provenance`.

**Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import read_alignment_file
from dnakit.sequence import write_pairwise_alignment_file

path = Path(tempfile.mkdtemp()) / "pair.maf"
write_pairwise_alignment_file("ACGTACGGTT", "ACGTACGTT", path, file_format="maf", verbose=False)
result = read_alignment_file(path, file_format="maf", exactly_one=True, verbose=False)
item = result["alignments"][0]
print(result["alignment_count"], item["sequence_ids"], item["score"])
print(item["aligned_rows"])
print(item["coordinates"])
```

**Example output:**

```text
1 ['query', 'target'] 16.0
['ACGTACGGTT', 'ACGTAC-GTT']
[[0, 6, 7, 10], [0, 6, 6, 9]]
```

## 2) `read_alignment_blocks`

- **Purpose:** Read every alignment block of a MAF or FASTA-M10 file.
- **Method:** Calls Bio.AlignIO.parse. A MAF file has one block per aligned region and a FASTA-M10 file one block per hit; the existing `read_alignment` uses AlignIO.read and raises 'More than one record' on several blocks. Each block lists the row names, gapped sequences and row annotations (the MAF `start`, `size`, `strand` and `srcSize`) and the block annotations.
- **API:** `read_alignment_blocks(path[required], file_format[optional], max_blocks[optional], verbose[optional])`
- **Input:** A file path, a format name (`maf` or `fasta-m10`) and a block count limit.
- **Output:** A dictionary with `block_count`, `blocks` (`rows`, `alignment_length`, `annotations`) and `provenance`.

**Example file:** `blocks.maf`

```text
##maf version=1 scoring=none
a score=23262.0
s hg18.chr7    27578828 38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s panTro1.chr6 28741140 38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s baboon         116834 38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG

a score=5062.0
s hg18.chr7    27699739 6 + 158545518 TAAAGA
s panTro1.chr6 28862317 6 + 161576975 TAAAGA
s baboon         241163 6 +   4622798 TAAAGA

```

**Example code:**

```python
from dnakit.sequence import read_alignment_blocks

result = read_alignment_blocks("blocks.maf", verbose=False)
print(result["block_count"])
for block in result["blocks"]:
    first = block["rows"][0]
    print(block["row_count"], block["alignment_length"], first["id"], first["annotations"]["start"])
```

**Example output:**

```text
2
3 42 hg18.chr7 27578828
3 6 hg18.chr7 27699739
```

## 3) `read_substitution_matrix_file`

- **Purpose:** Read a custom nucleotide substitution matrix file.
- **Method:** Calls Bio.Align.substitution_matrices.read. The file has optional `#` comment lines, a header row of symbols and one row per symbol. The alphabet may contain only A, C, G, T, U and IUPAC symbols, otherwise `ConfigurationError` is raised. Bundled matrices are still loaded by name with `load_substitution_matrix`.
- **API:** `read_substitution_matrix_file(path[required], verbose[optional])`
- **Input:** A matrix file path.
- **Output:** A dictionary with `alphabet`, `size`, `symmetric`, `matrix`, `comments` and `provenance`.

**Example file:** `nuc.mat`

```text
# demo nucleotide matrix
   A  C  G  T
A  2 -3 -3 -3
C -3  2 -3 -3
G -3 -3  2 -3
T -3 -3 -3  2
```

**Example code:**

```python
from dnakit.sequence import read_substitution_matrix_file

result = read_substitution_matrix_file("nuc.mat", verbose=False)
print(result["alphabet"], result["symmetric"], result["comments"])
print(result["matrix"][0])
```

**Example output:**

```text
ACGT True ['demo nucleotide matrix']
[2.0, -3.0, -3.0, -3.0]
```

## 4) `read_blast_records`

- **Purpose:** Read BLAST XML or XML2 output, one record per query.
- **Method:** Calls Bio.Blast.parse; with `exactly_one=True` it calls Bio.Blast.read. Each record gives the query, program, version and database, and for every hit its high-scoring pairs (score, bit score, E-value, identity count and the target and query coordinates, 0-based and half-open). The SearchIO parser behind `parse_search_results` returns an empty result, without an error, for the XML2 layout.
- **API:** `read_blast_records(path[required], exactly_one[optional], max_records[optional], max_hits[optional], verbose[optional])`
- **Input:** A file path, whether exactly one query is required, and record and hit count limits.
- **Output:** A dictionary with `record_count`, `records` (`hits`, `hsps`) and `provenance`.

**Example file:** `blast_small.xml`

```text
<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
<BlastOutput>
  <BlastOutput_program>blastn</BlastOutput_program>
  <BlastOutput_version>BLASTN 2.13.0+</BlastOutput_version>
  <BlastOutput_db>demo_db</BlastOutput_db>
  <BlastOutput_query-ID>Query_1</BlastOutput_query-ID>
  <BlastOutput_query-def>query_one demo query</BlastOutput_query-def>
  <BlastOutput_query-len>20</BlastOutput_query-len>
    <BlastOutput_iterations>
    <Iteration>
      <Iteration_iter-num>1</Iteration_iter-num>
      <Iteration_query-ID>Query_1</Iteration_query-ID>
      <Iteration_query-def>query_one demo query</Iteration_query-def>
      <Iteration_query-len>20</Iteration_query-len>
      <Iteration_hits>
        <Hit>
          <Hit_num>1</Hit_num>
          <Hit_id>gnl|BL_ORD_ID|0</Hit_id>
          <Hit_def>subject_one demo subject</Hit_def>
          <Hit_accession>0</Hit_accession>
          <Hit_len>60</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>37.4</Hsp_bit-score>
              <Hsp_score>20</Hsp_score>
              <Hsp_evalue>1.2e-05</Hsp_evalue>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>20</Hsp_query-to>
              <Hsp_hit-from>11</Hsp_hit-from>
              <Hsp_hit-to>30</Hsp_hit-to>
              <Hsp_query-frame>1</Hsp_query-frame>
              <Hsp_hit-frame>1</Hsp_hit-frame>
              <Hsp_identity>20</Hsp_identity>
              <Hsp_positive>20</Hsp_positive>
              <Hsp_gaps>0</Hsp_gaps>
              <Hsp_align-len>20</Hsp_align-len>
              <Hsp_qseq>ACGTACGGTTAACCGGATAC</Hsp_qseq>
              <Hsp_hseq>ACGTACGGTTAACCGGATAC</Hsp_hseq>
              <Hsp_midline>||||||||||||||||||||</Hsp_midline>
            </Hsp>
          </Hit_hsps>
        </Hit>
      </Iteration_hits>
          </Iteration>
  </BlastOutput_iterations>
</BlastOutput>
```

**Example code:**

```python
from dnakit.sequence import read_blast_records

result = read_blast_records("blast_small.xml", verbose=False)
record = result["records"][0]
hit = record["hits"][0]
hsp = hit["hsps"][0]
print(record["query_id"], record["query_length"], record["program"], record["database"])
print(hit["target_id"], hit["target_length"])
print(hsp["evalue"], hsp["bit_score"], hsp["target_coordinates"], hsp["query_coordinates"])
```

**Example output:**

```text
Query_1 20 blastn demo_db
gnl|BL_ORD_ID|0 60
1.2e-05 37.4 [10, 30] [0, 20]
```

## 5) `parse_primer3_records`

- **Purpose:** Read every record of an EMBOSS eprimer3 or Primer3 report (one per target sequence).
- **Method:** Calls Bio.Emboss.Primer3.parse. The existing `parse_primer3_output` uses Primer3.read and raises `ConfigurationError` on a multi-target report. Each record gives its header comments and the primer pairs (forward, reverse and an optional internal oligo). An empty file, or one without records, raises `ConfigurationError` (the upstream parse loops forever on an empty stream, so it is checked first).
- **API:** `parse_primer3_records(path[required], max_records[optional], verbose[optional])`
- **Input:** A file path and a record count limit.
- **Output:** A dictionary with `record_count`, `pair_count`, `records` and `provenance`.

**Example file:** `primer3_multi.txt`

```text
# EPRIMER3 RELEASE 5.0.0
#                      Start  Len   Tm     GC%   Sequence

   1 PRODUCT SIZE: 60
     FORWARD PRIMER    1  20  60.00  50.00  ACGTACGTACGTACGTACGT
     REVERSE PRIMER   60  20  60.00  50.00  TGCATGCATGCATGCATGCA

# EPRIMER3 RELEASE 5.0.0
#                      Start  Len   Tm     GC%   Sequence

   1 PRODUCT SIZE: 80
     FORWARD PRIMER    5  18  59.50  55.00  GGCCTTAAGGCCTTAAGG
     REVERSE PRIMER   80  18  59.80  50.00  CCGGAATTCCGGAATTCC
```

**Example code:**

```python
from dnakit.sequence import parse_primer3_records

result = parse_primer3_records("primer3_multi.txt", verbose=False)
print(result["record_count"], result["pair_count"])
for record in result["records"]:
    pair = record["primer_pairs"][0]
    print(pair["product_size"], pair["forward"]["sequence"], pair["reverse"]["start"])
```

**Example output:**

```text
2 2
60 ACGTACGTACGTACGTACGT 60
80 GGCCTTAAGGCCTTAAGG 80
```

## 6) `read_genbank_record_fields`

- **Purpose:** Read the header fields and feature table of GenBank records.
- **Method:** Calls Bio.GenBank.parse; with `exactly_one=True` it calls Bio.GenBank.read. The low-level parser keeps fields that `read_sequence_records` drops: the LOCUS line (name, length, molecule type, division, date), accessions, version, keywords, source and organism, taxonomy, the number of references, the comment and, for every feature, the location text and qualifiers.
- **API:** `read_genbank_record_fields(path[required], exactly_one[optional], max_records[optional], max_bases[optional], verbose[optional])`
- **Input:** A file path, whether exactly one record is required, and record and base count limits.
- **Output:** A dictionary with `record_count`, `total_bases`, `records` and `provenance`.

**Example file:** `two.gb`

```text
LOCUS       SEQ1                      18 bp    DNA     linear   BCT 01-JAN-2020
DEFINITION  Demo sequence one.
ACCESSION   SEQ1
VERSION     SEQ1.1
KEYWORDS    demo; test.
SOURCE      Escherichia coli
  ORGANISM  Escherichia coli
            Bacteria; Proteobacteria; Gammaproteobacteria.
FEATURES             Location/Qualifiers
     source          1..18
                     /organism="Escherichia coli"
     gene            1..12
                     /gene="abcA"
ORIGIN
        1 atgaaacgta ttgcgtaa
//
LOCUS       SEQ2                      12 bp    DNA     linear   BCT 02-JAN-2020
DEFINITION  Demo sequence two.
ACCESSION   SEQ2
VERSION     SEQ2.1
KEYWORDS    .
SOURCE      Bacillus subtilis
  ORGANISM  Bacillus subtilis
            Bacteria; Firmicutes; Bacilli.
FEATURES             Location/Qualifiers
     CDS             1..12
                     /gene="xyzB"
ORIGIN
        1 atggcgtacg tt
//
```

**Example code:**

```python
from dnakit.sequence import read_genbank_record_fields

result = read_genbank_record_fields("two.gb", verbose=False)
print(result["record_count"], result["total_bases"])
for record in result["records"]:
    print(record["locus"], record["data_file_division"], record["date"], record["organism"])
    print(record["taxonomy"], [f["type"] for f in record["features"]])
```

**Example output:**

```text
2 30
SEQ1 BCT 01-JAN-2020 Escherichia coli
['Bacteria', 'Proteobacteria', 'Gammaproteobacteria'] ['source', 'gene']
SEQ2 BCT 02-JAN-2020 Bacillus subtilis
['Bacteria', 'Firmicutes', 'Bacilli'] ['CDS']
```

## 7) `read_phylogeny_collection`

- **Purpose:** Read every tree of a Newick, NEXUS, NeXML, phyloXML or CDAO file.
- **Method:** Calls Bio.Phylo.parse, which reaches the parsers of NexusIO, NeXMLIO, PhyloXMLIO and CDAOIO. `read_tree` reads one tree and raises on several, so bootstrap or posterior tree sets need this function. Each tree reports its name, rooted flag, sorted tip names, total branch length and a Newick string. A file without trees raises `ConfigurationError`.
- **API:** `read_phylogeny_collection(path[required], file_format[optional], max_trees[optional], verbose[optional])`
- **Input:** A file path, a format name and a tree count limit.
- **Output:** A dictionary with `tree_count`, `trees` and `provenance`.

**Example file:** `trees.nexus`

```text
#NEXUS
BEGIN TAXA;
DIMENSIONS NTAX=4;
TAXLABELS A B C D;
END;
BEGIN TREES;
TREE t1 = ((A:0.1,B:0.2):0.05,(C:0.3,D:0.4):0.06);
TREE t2 = ((A:0.1,C:0.2):0.05,(B:0.3,D:0.4):0.06);
END;
```

**Example code:**

```python
from dnakit.sequence import read_phylogeny_collection

result = read_phylogeny_collection("trees.nexus", file_format="nexus", verbose=False)
print(result["tree_count"])
for tree in result["trees"]:
    print(tree["name"], tree["tips"], round(tree["total_branch_length"], 2))
    print(tree["newick"])
```

**Example output:**

```text
2
t1 ['A', 'B', 'C', 'D'] 1.11
((A:0.1,B:0.2):0.05,(C:0.3,D:0.4):0.06):0;
t2 ['A', 'B', 'C', 'D'] 1.11
((A:0.1,C:0.2):0.05,(B:0.3,D:0.4):0.06):0;
```

## 8) `bootstrap_replicate_trees`

- **Purpose:** Build the bootstrap replicate trees of an aligned set of DNA sequences.
- **Method:** Calls Bio.Phylo.Consensus.bootstrap_trees: alignment columns are resampled with replacement and one tree (`nj` or `upgma`) is built per replicate. `bootstrap_consensus_tree` returns only the majority-rule consensus; this function returns the replicate trees, which `bootstrap_branch_support` can count. The global random state is seeded with `seed` for the call and restored afterwards.
- **API:** `bootstrap_replicate_trees(sequences[required], replicates[optional], method[optional], distance_model[optional], seed[optional], max_length[optional], verbose[optional])`
- **Input:** An equal-length DNA alignment (3 to 200 rows), the replicate count, method, distance model and seed.
- **Output:** A dictionary with `trees` (Newick list), `ids`, `replicates`, `seed` and `provenance`.

**Example code:**

```python
import random

from dnakit.sequence import bootstrap_replicate_trees

alignment = {
    "a": "ACGTACGTAC",
    "b": "ACGTACGAAC",
    "c": "ACGAACGTTC",
    "d": "TCGTACGTAC",
}
random.seed(5)
expected = random.random()
random.seed(5)
first = bootstrap_replicate_trees(alignment, replicates=4, seed=7, verbose=False)
second = bootstrap_replicate_trees(alignment, replicates=4, seed=7, verbose=False)
print(len(first["trees"]), first["trees"] == second["trees"], first["ids"])
print(all(name in tree for tree in first["trees"] for name in "abcd"))
print(random.random() == expected)
```

**Example output:**

```text
4 True ['a', 'b', 'c', 'd']
True
True
```

## 9) `read_phd_records`

- **Purpose:** Read every PHRED base-call record of a PHD file.
- **Method:** Calls Bio.Sequencing.Phd.parse. `read_phd_sequence` uses Phd.read and accepts exactly one record. Each record gives the file name, the base calls, the per-base PHRED qualities, the peak locations (null when the file has none) and the comment fields.
- **API:** `read_phd_records(path[required], max_records[optional], verbose[optional])`
- **Input:** A file path and a record count limit.
- **Output:** A dictionary with `record_count`, `records` and `provenance`.

**Example file:** `reads.phd`

```text
BEGIN_SEQUENCE read1
BEGIN_COMMENT
CHROMAT_FILE: read1.ab1
PHRED_VERSION: 0.980528.c
END_COMMENT
BEGIN_DNA
A 10 1
C 20 4
END_DNA
END_SEQUENCE

BEGIN_SEQUENCE read2
BEGIN_COMMENT
CHROMAT_FILE: read2.ab1
PHRED_VERSION: 0.980528.c
END_COMMENT
BEGIN_DNA
G 30 2
T 40 3
A 50 5
END_DNA
END_SEQUENCE
```

**Example code:**

```python
from dnakit.sequence import read_phd_records

result = read_phd_records("reads.phd", verbose=False)
print(result["record_count"])
for record in result["records"]:
    print(record["file_name"], record["sequence"], record["qualities"], record["peak_locations"])
```

**Example output:**

```text
2
read1 AC [10, 20] [1, 4]
read2 GTA [30, 40, 50] [2, 3, 5]
```

## 10) `read_goa_gene_product_file`

- **Purpose:** Read a GOA gene-product association (GPA) or gene-product information (GPI) file.
- **Method:** `file_format='gpa'` calls Bio.UniProt.GOA.gpa_iterator (versions 1.0 and 1.1: gene product, GO term, evidence and reference per line); `'gpi'` calls gpi_iterator (versions 1.0, 1.1 and 1.2: symbol, name, synonyms, type, taxon and cross references). The version comes from the required `!gpa-version` or `!gpi-version` header line. `read_gaf_annotations` reads the GAF layout only.
- **API:** `read_goa_gene_product_file(path[required], file_format[optional], max_records[optional], verbose[optional])`
- **Input:** A file path, a format name (`gpa` or `gpi`) and a record count limit.
- **Output:** A dictionary with `version`, `record_count`, `records` (GOA field names; list fields are lists) and `provenance`.

**Example file:** `prod.gpa`

```text
!gpa-version: 1.1
UniProtKB	P12345	enables	GO:0003700	PMID:1234567	ECO:0000314			20240101	UniProt		
UniProtKB	P67890	involved_in	GO:0006355	PMID:7654321|PMID:1	ECO:0000501	GO:1		20240102	UniProt		
```

**Example code:**

```python
from dnakit.sequence import read_goa_gene_product_file

result = read_goa_gene_product_file("prod.gpa", verbose=False)
print(result["version"], result["record_count"])
for record in result["records"]:
    print(record["DB_Object_ID"], record["GO_ID"], record["Qualifier"], record["DB:Reference"], record["With"])
```

**Example output:**

```text
1.1 2
P12345 GO:0003700 ['enables'] ['PMID:1234567'] []
P67890 GO:0006355 ['involved_in'] ['PMID:7654321', 'PMID:1'] ['GO:1']
```

## 11) `correct_jc69_distance`

- **Purpose:** Apply the Jukes-Cantor (JC69) correction to observed p-distances.
- **Method:** Calls skbio.sequence.distance.jc69_correct: D = -a ln(1 - p/a) with a = (n-1)/n for an alphabet of n symbols (4 for DNA). Each p must lie in [0, 1]. A p at or above a (saturation) has no finite distance; it is null in `corrected` and the status is `partial`. `nucleotide_distance(model='jc69')` starts from two sequences and cannot correct a p-distance measured elsewhere.
- **API:** `correct_jc69_distance(p_distances[required], alphabet_size[optional], verbose[optional])`
- **Input:** A list of p-distances and the alphabet size.
- **Output:** A dictionary with `corrected`, `saturated_count`, `saturation_limit`, `status` and `provenance`.

**Example code:**

```python
from dnakit.sequence import correct_jc69_distance

result = correct_jc69_distance([0.0, 0.1, 0.5, 0.8], verbose=False)
print(result["status"], result["saturation_limit"], result["saturated_count"])
print([None if v is None else round(v, 6) for v in result["corrected"]])
```

**Example output:**

```text
partial 0.75 1
[0.0, 0.107326, 0.823959, None]
```

## 12) `query_saved_khmer_table`

- **Purpose:** Look up k-mers in a saved khmer counttable, nodetable or nodegraph.
- **Method:** `table_type` must match the class that wrote the file, because khmer stores no marker for it: `counttable` calls khmer.load_counttable and returns counts; `nodetable` and `nodegraph` call load_nodetable and load_nodegraph and return 1 when the k-mer was added and 0 otherwise. `query_saved_countgraph` reads a counttable file with another hash layout and silently returns 0. A k-mer and its reverse complement share one entry; the result is sorted by k-mer.
- **API:** `query_saved_khmer_table(source[required], kmers[required], table_type[optional], verbose[optional])`
- **Input:** A file path, a list of k-mers (length equal to the table ksize) and the table type.
- **Output:** A dictionary with `values`, `value_kind` (`count` or `presence`), `ksize`, `hashsizes` and `provenance`.

**Example code:**

```python
import tempfile
from pathlib import Path

import khmer

from dnakit.microbiome import query_saved_khmer_table
from dnakit.microbiome import query_saved_countgraph

table = khmer.Counttable(5, 1000, 4)
table.consume("ACGTACGGATCGATCGGCTAGCTAGG")
table.consume("ACGTA")
path = Path(tempfile.mkdtemp()) / "demo.ct"
table.save(str(path))
result = query_saved_khmer_table(path, ["ACGTA", "TACGT", "GGGGG"], verbose=False)
print(result["value_kind"], result["ksize"], result["values"])
print(query_saved_countgraph(path, ["ACGTA"], verbose=False)["counts"])
```

**Example output:**

```text
count 5 {'ACGTA': 2, 'GGGGG': 0, 'TACGT': 2}
{'ACGTA': 0}
```
