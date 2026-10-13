# G18o functions: Biopython, scikit-bio, cutadapt and umi_tools review wrappers

This page documents the 26 wrappers that review group G18o adds in the `sequence` and `sequencing_qc` areas. The functions are exported from their area packages, for example `from dnakit.sequence import write_motifs`.

- Offline: every function reads only call arguments or local files (SAM, BAM, FASTA, QUAL, BLAST tables and Newick strings) and never uses the network.
- Functions that write files create a new file only; an existing output path raises `ConfigurationError`, so the examples can be run repeatedly.
- Inputs needed by an example are given with the `Example file` label, and the example code reads them from the run directory; examples that need intermediate files use a temporary directory.
- pyslow5 1.5.0 can crash the interpreter on a corrupt SLOW5/BLOW5 file. G18o adds no pyslow5 function; the existing `read_slow5_signal` and related backends check the SLOW5/BLOW5 magic bytes before opening a file, and that check must stay.
- parasail's `ssw_init` segfaults with score_size 16. G18o does not wrap it; a future wrapper must guard score_size first.

## Ledger status (88 other and cap rows)

| Package | Other/cap before | Wrapped | Other statuses |
| --- | ---: | ---: | --- |
| Bio | 39 | 13 rows (14 functions, including `write_gaf_file`) | family-variant 11, duplicate-of 1, not-dna 6, internal 6, deprecated 2 |
| skbio | 10 | 5 | internal 4, not-dna 1 |
| parasail | 13 | 0 | internal 7, not-dna 6 (CPU probes and a timer) |
| mappy | 1 | 0 | internal 1 (`mappy.verbose`) |
| pyhmmer | 1 | 0 | internal 1 (`hmmpress`) |
| cutadapt | 11 | 1 | cli 4, internal 5, duplicate-of 1 |
| umi_tools | 11 | 6 | family-variant 1, cli 1, internal 3 |
| pyfastx | 1 | 0 | internal 1 (`version`) |
| pyslow5 | 1 | 0 | already-wrapped 1 (existing `read_slow5_signal` and related) |

Main reasons for leaving rows unwrapped: `not-dna` covers Entrez, MEDLINE, SCAD and SFF metadata readers and parasail's CPU probes; `internal` covers helpers, global setters (which change process state) and report formatters; `cli` covers command-line reports and the RSEM script steps; `deprecated` covers the retired UniGene records.

## 1) `write_pairwise_alignment_file`

- **Purpose:** Align two DNA sequences with Biopython and write the best alignment to a new file.
- **Method:** Scores are fixed at match 2, mismatch -1, gap open -2 and gap extend -1, with global alignment by default; `mode="local"` gives a local alignment. Formats are fasta, clustal, maf, phylip, psl, sam or bed. The output file must not exist.
- **API:** `write_pairwise_alignment_file(query[required], target[required], output_path[required], file_format[optional], mode[optional], query_id[optional], target_id[optional], verbose[optional])`
- **Input:** Two DNA strings (IUPAC symbols allowed), a new file path, a format name, the alignment mode and sequence names.
- **Output:** A dictionary with `score`, `aligned_query`, `aligned_target`, `alignment_count`, `output_path` and `provenance`.

**Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import write_pairwise_alignment_file

out = Path(tempfile.mkdtemp()) / "pair.fa"
result = write_pairwise_alignment_file("ACGTACGT", "ACGAACGT", out, verbose=False)
print(result["score"], result["aligned_query"], result["aligned_target"])
print(out.read_text().strip())
```

**Example output:**

```text
13.0 ACGTACGT ACGAACGT
>query
ACGTACGT
>target
ACGAACGT
```

## 2) `convert_search_results`

- **Purpose:** Convert a local BLAST, BLAT or HMMER search result to another format.
- **Method:** Calls Bio.SearchIO.convert. The input and output formats must be in an allowed list. When upstream cannot write the requested output from the input (for example blast-tab cannot give blast-xml), the call raises `ConfigurationError`.
- **API:** `convert_search_results(source[required], source_format[required], output_path[required], output_format[required], verbose[optional])`
- **Input:** A search-result file (for example a BLAST tabular file without a header), the input format name and a new output path.
- **Output:** A dictionary with `query_count`, `hit_count`, `hsp_count`, `fragment_count`, `output_path` and `provenance`.

**Example file:** `hits.tsv`

```text
q1	t1	100.00	8	0	0	1	8	1	8	1e-05	16.0
q2	t2	87.50	8	1	0	1	8	1	8	2e-03	12.5
```

**Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import convert_search_results

out = Path(tempfile.mkdtemp()) / "converted.tsv"
result = convert_search_results("hits.tsv", "blast-tab", out, "blast-tab", verbose=False)
print(result["query_count"], result["hit_count"], result["hsp_count"])
print(out.read_text().strip())
```

**Example output:**

```text
2 2 2
q1	t1	100.00	8	0	0	1	8	1	8	1e-05	16.0
q2	t2	87.50	8	1	0	1	8	1	8	0.002	12.5
```

## 3) `codon_back_table`

- **Purpose:** Return one codon for each amino acid of an NCBI genetic code, with stop written as `*`.
- **Method:** Calls Bio.Data.CodonTable.make_back_table. The table keeps only the last codon of each amino acid in sorted order, so the result is deterministic but lists one synonymous codon, not all of them.
- **API:** `codon_back_table(genetic_code[optional], stop_codon[optional], verbose[optional])`
- **Input:** An NCBI genetic code number from 1 to 33 (an unambiguous DNA table) and a stop codon of that table.
- **Output:** A dictionary with `table_name`, `codon_count`, `back_table` (amino acid to codon) and `provenance`.

**Example code:**

```python
from dnakit.sequence import codon_back_table

result = codon_back_table(1, verbose=False)
back = result["back_table"]
print(result["table_name"], result["codon_count"])
print(back["M"], back["L"], back["*"])
```

**Example output:**

```text
Standard 61
ATG TTG TAA
```

## 4) `adam_consensus_tree`

- **Purpose:** Compute the Adam consensus of several Newick trees.
- **Method:** Each tree is parsed into a Bio.Phylo tree, then Bio.Phylo.Consensus.adam_consensus is called. The output is the upstream result: consensus branch lengths are set to 0 by upstream, and conflicting topologies follow upstream's algorithm.
- **API:** `adam_consensus_tree(newick_trees[required], verbose[optional])`
- **Input:** A sequence of Newick strings, from 1 to 1000 trees, with the same tip names.
- **Output:** A dictionary with `tree_count`, `newick` and `provenance`.

**Example code:**

```python
from dnakit.sequence import adam_consensus_tree

newick = "((A:1,B:1):1,(C:1,D:1):1);"
result = adam_consensus_tree([newick, newick], verbose=False)
print(result["tree_count"], result["newick"])
```

**Example output:**

```text
2 ((C:1,D:1):0,(A:1,B:1):0):0;
```

## 5) `bootstrap_branch_support`

- **Purpose:** Annotate a target tree with how many bootstrap trees support each internal branch.
- **Method:** Calls Bio.Phylo.Consensus.get_support. Support is the percentage of replicate trees that contain the branch, written on the internal node. For example, 2 of 3 replicates give 66.67.
- **API:** `bootstrap_branch_support(target_newick[required], replicate_newicks[required], verbose[optional])`
- **Input:** The target Newick string and 1 to 1000 replicate trees with the same tip names.
- **Output:** A dictionary with `replicate_count`, `scale` (always percent), `newick` and `provenance`.

**Example code:**

```python
from dnakit.sequence import bootstrap_branch_support

target = "((A:0.1,B:0.2):0.3,(C:0.1,D:0.2):0.3);"
replicates = [target, "((A:0.1,C:0.2):0.3,(B:0.1,D:0.2):0.3);", target]
result = bootstrap_branch_support(target, replicates, verbose=False)
print(result["replicate_count"], result["scale"])
print(result["newick"])
```

**Example output:**

```text
3 percent
((A:0.1,B:0.2)66.67:0.3,(C:0.1,D:0.2)66.67:0.3)100.00:0;
```

## 6) `tree_to_graph_edges`

- **Purpose:** Return the nodes and parent-to-child edges of a Newick tree as plain lists, without networkx or igraph.
- **Method:** Nodes are numbered in pre-order from 0 at the root, and each parent-to-child pair is an edge. The edge set equals the one from Bio.Phylo.to_networkx and Bio.Phylo.to_igraph.
- **API:** `tree_to_graph_edges(newick[required], verbose[optional])`
- **Input:** One Newick string.
- **Output:** A dictionary with `node_count`, `edge_count`, `tip_count`, `nodes` (id, name, branch length, tip flag), `edges` (source, target, branch length) and `provenance`.

**Example code:**

```python
from dnakit.sequence import tree_to_graph_edges

result = tree_to_graph_edges("((A:0.1,B:0.2):0.3,C:0.5);", verbose=False)
print(result["node_count"], result["edge_count"], result["tip_count"])
for edge in result["edges"]:
    print(edge["source"], edge["target"], edge["branch_length"])
```

**Example output:**

```text
5 4 3
0 1 0.3
0 4 0.5
1 2 0.1
1 3 0.2
```

## 7) `complement_rna_sequence`

- **Purpose:** Return the complement of an RNA string, keeping U.
- **Method:** Calls Bio.Seq.complement_rna (A pairs with U, C with G). IUPAC symbols are allowed; a T in the input is rejected.
- **API:** `complement_rna_sequence(sequence[required], verbose[optional])`
- **Input:** One RNA string (ACGU and IUPAC symbols).
- **Output:** A dictionary with `input`, `sequence`, `length` and `provenance`.

**Example code:**

```python
from dnakit.sequence import complement_rna_sequence

print(complement_rna_sequence("ACGU", verbose=False)["sequence"])
```

**Example output:**

```text
UGCA
```

## 8) `reverse_complement_rna_sequence`

- **Purpose:** Return the reverse complement of an RNA string, keeping U.
- **Method:** Calls Bio.Seq.reverse_complement_rna: the string is reversed and then complemented.
- **API:** `reverse_complement_rna_sequence(sequence[required], verbose[optional])`
- **Input:** One RNA string (ACGU and IUPAC symbols).
- **Output:** A dictionary with `input`, `sequence`, `length` and `provenance`.

**Example code:**

```python
from dnakit.sequence import reverse_complement_rna_sequence

print(reverse_complement_rna_sequence("AAGC", verbose=False)["sequence"])
```

**Example output:**

```text
GCUU
```

## 9) `read_fasta_qual_pairs`

- **Purpose:** Read a matched FASTA and QUAL file pair as sequences with PHRED qualities.
- **Method:** Calls Bio.SeqIO.QualityIO.PairedFastaQualIterator, which pairs the two files by position. Qualities are integers. More than `max_records` records raise an error.
- **API:** `read_fasta_qual_pairs(fasta_path[required], qual_path[required], max_records[optional], verbose[optional])`
- **Input:** A local FASTA file and a local QUAL file.
- **Output:** A dictionary with `record_count`, `records` (each with id, description, sequence and qualities) and `provenance`.

**Example file:** `reads.fa`

```text
>r1 desc
ACGT
>r2
GGCC
```

**Example file:** `reads.qual`

```text
>r1 desc
30 20 10 5
>r2
40 40 40 40
```

**Example code:**

```python
from dnakit.sequence import read_fasta_qual_pairs

result = read_fasta_qual_pairs("reads.fa", "reads.qual", verbose=False)
print(result["record_count"])
for record in result["records"]:
    print(record["id"], record["sequence"], record["qualities"])
```

**Example output:**

```text
2
r1 ACGT [30, 20, 10, 5]
r2 GGCC [40, 40, 40, 40]
```

## 10) `write_qual_file`

- **Purpose:** Write the PHRED qualities of DNA records to a new QUAL file.
- **Method:** Each record needs one quality per base, with values from 0 to 93. The file is written with Biopython's `qual` format writer.
- **API:** `write_qual_file(records[required], output_path[required], verbose[optional])`
- **Input:** A list of records (each with `id`, `sequence` and `qualities`) and a new output path.
- **Output:** A dictionary with `output_path`, `record_count` and `provenance`.

**Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import write_qual_file

out = Path(tempfile.mkdtemp()) / "reads.qual"
record = {"id": "r1", "sequence": "ACGT", "qualities": [30, 20, 10, 5]}
result = write_qual_file([record], out, verbose=False)
print(result["record_count"])
print(out.read_text().strip())
```

**Example output:**

```text
1
>r1
30 20 10 5
```

## 11) `melting_temperature_gc`

- **Purpose:** Estimate the melting temperature of an unambiguous DNA primer from its GC content.
- **Method:** Calls Bio.SeqUtils.MeltingTemp.Tm_GC with check and strict both true, so ambiguous bases such as N raise an error. `valueset` selects one of coefficient sets 1 to 7; salt concentrations are in mM; `saltcorr` is a salt-correction method from 0 to 7. This is not the nearest-neighbour model; see the existing `melting_temperature`.
- **API:** `melting_temperature_gc(sequence[required], valueset[optional], na[optional], k[optional], tris[optional], mg[optional], dntps[optional], saltcorr[optional], mismatch[optional], verbose[optional])`
- **Input:** One unambiguous DNA string and the salt and correction options.
- **Output:** A dictionary with `length`, `tm_celsius`, `valueset`, `saltcorr`, `salt_mM` and `provenance`.

**Example code:**

```python
from dnakit.sequence import melting_temperature_gc

result = melting_temperature_gc("ACGTACGTACGTACGTAAGC", verbose=False)
print(result["length"], round(result["tm_celsius"], 2))
```

**Example output:**

```text
20 50.4
```

## 12) `write_gaf_file`

- **Purpose:** Write gene ontology (GO) annotations to a new GAF 2.0 file.
- **Method:** Each record is a mapping with GAF 2.0 field names; a missing field is written empty, and a list value is joined with `|`. Text containing a tab, a newline or `|` is rejected. Biopython's `writerec` writes each line. This is the representative for `writerec` and `writebyproteinrec`.
- **API:** `write_gaf_file(records[required], output_path[required], verbose[optional])`
- **Input:** A list of GAF records and a new output path.
- **Output:** A dictionary with `output_path`, `record_count`, `gaf_version` (always 2.0) and `provenance`.

**Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import write_gaf_file

record = {
    "DB": "UniProtKB",
    "DB_Object_ID": "P12345",
    "DB_Object_Symbol": "GENE1",
    "GO_ID": "GO:0005634",
    "Evidence": "IDA",
    "Aspect": "C",
    "Taxon_ID": ["taxon:9606"],
}
out = Path(tempfile.mkdtemp()) / "annotations.gaf"
result = write_gaf_file([record], out, verbose=False)
fields = out.read_text().rstrip("\n").split("\t")
print(result["record_count"], result["gaf_version"], len(fields))
print(fields[0], fields[1], fields[2], fields[4], fields[12])
```

**Example output:**

```text
1 2.0 17
UniProtKB P12345 GENE1 GO:0005634 taxon:9606
```

## 13) `read_mast_report`

- **Purpose:** Read a local MEME Suite MAST XML report and return its motifs, search settings and per-sequence hit diagrams.
- **Method:** Calls Bio.motifs.mast.read. A diagram lists the hits in order: numbers are gap lengths, `[+1]` is a forward hit of motif 1, and `[-1]` is a reverse-complement hit. A malformed report raises `ConfigurationError`.
- **API:** `read_mast_report(path[required], verbose[optional])`
- **Input:** One local MAST XML file.
- **Output:** A dictionary with `version`, `database`, `alphabet`, `strand_handling`, `motifs`, `sequences` (names and diagrams) and `provenance`.

**Example file:** `mast.xml`

```text
<?xml version="1.0" encoding="UTF-8"?>
<mast version="5.5.0">
<sequence_dbs><sequence_db source="seqs.fa"/></sequence_dbs>
<alphabet name="DNA"/>
<settings strand_handling="combine"/>
<motifs>
<motif id="MOTIF1" alt="m1" length="4"/>
<motif id="MOTIF2" alt="m2" length="3"/>
</motifs>
<sequences>
<sequence name="s1" length="20">
<seg start="1"><hit idx="0" pos="3" rc="n"/></seg>
</sequence>
<sequence name="s2" length="20"/>
</sequences>
</mast>
```

**Example code:**

```python
from dnakit.sequence import read_mast_report

result = read_mast_report("mast.xml", verbose=False)
print(result["version"], result["alphabet"], result["strand_handling"])
print([motif["name"] for motif in result["motifs"]])
for item in result["sequences"]:
    print(item["name"], item["diagram"])
```

**Example output:**

```text
5.5.0 DNA combine
['1', '2']
s1 2-[+1]-14
s2 20
```

## 14) `write_motifs`

- **Purpose:** Format DNA motifs (several equal-length sites per motif) as JASPAR, PFM, TRANSFAC or ClusterBuster text.
- **Method:** Each motif is built with Bio.motifs.create and formatted with Bio.motifs.write. Names are written in sorted order; `pfm` holds exactly one motif. The result has the text, each motif's length and consensus.
- **API:** `write_motifs(motifs[required], file_format[required], verbose[optional])`
- **Input:** A mapping from motif name to its list of sites, for example `{"M1": ["ACGT", "ACGA"]}`, and a format name.
- **Output:** A dictionary with `file_format`, `motif_count`, `motifs` (name, length, consensus), `text` and `provenance`.

**Example code:**

```python
from dnakit.sequence import write_motifs

result = write_motifs({"M1": ["ACGT", "ACGA", "ACGT"]}, "jaspar", verbose=False)
print(result["motif_count"], result["motifs"][0]["consensus"])
print(result["text"].strip())
```

**Example output:**

```text
1 ACGT
>M1 M1
A [  3.00   0.00   0.00   1.00]
C [  0.00   3.00   0.00   0.00]
G [  0.00   0.00   3.00   0.00]
T [  0.00   0.00   0.00   2.00]
```

## 15) `build_bme_tree`

- **Purpose:** Build a balanced minimum evolution tree from the p-distances of aligned DNA.
- **Method:** The p-distance is the number of differing sites divided by the number of sites where both sequences have an unambiguous base (gaps and ambiguity codes are ignored). Then skbio.tree.bme is called. A pair with no comparable site raises an error.
- **API:** `build_bme_tree(sequences[required], neg_as_zero[optional], verbose[optional])`
- **Input:** A mapping from taxon name to an equal-length aligned DNA string (gaps `-` allowed), 3 to 500 taxa.
- **Output:** A dictionary with `taxa` (sorted by name), `distance` (always p-distance), `newick` and `provenance`.

**Example code:**

```python
from dnakit.sequence import build_bme_tree

sequences = {"A": "ACGTACGTAC", "B": "ACGTACGAAC", "C": "TCGTACGTAC", "D": "TCGTTCGTAG"}
result = build_bme_tree(sequences, verbose=False)
print(result["taxa"], result["distance"])
print(result["newick"])
```

**Example output:**

```text
['A', 'B', 'C', 'D'] p-distance
(B:0.1,(C:0.0,D:0.2):0.1,A:0.0);
```

## 16) `build_gme_tree`

- **Purpose:** Build a greedy minimum evolution tree from the p-distances of aligned DNA.
- **Method:** Uses the same p-distances as `build_bme_tree`, then calls skbio.tree.gme.
- **API:** `build_gme_tree(sequences[required], neg_as_zero[optional], verbose[optional])`
- **Input:** As for `build_bme_tree`.
- **Output:** A dictionary with `taxa`, `distance`, `newick` and `provenance`.

**Example code:**

```python
from dnakit.sequence import build_gme_tree

sequences = {"A": "ACGTACGTAC", "B": "ACGTACGAAC", "C": "TCGTACGTAC", "D": "TCGTTCGTAG"}
result = build_gme_tree(sequences, verbose=False)
print(result["taxa"], result["distance"])
print(result["newick"])
```

**Example output:**

```text
['A', 'B', 'C', 'D'] p-distance
(B:0.1,(C:0.0,D:0.2):0.1,A:0.0);
```

## 17) `improve_tree_nni`

- **Purpose:** Apply nearest-neighbour interchanges (NNI) to a Newick tree so that it fits the p-distances of aligned DNA better.
- **Method:** Tip names must match the sequence names exactly, and the tree needs at least four tips because a three-tip tree has no internal edge to swap. Calls skbio.tree.nni; an upstream ValueError becomes `ConfigurationError`.
- **API:** `improve_tree_nni(newick[required], sequences[required], balanced[optional], neg_as_zero[optional], verbose[optional])`
- **Input:** A Newick string and the aligned DNA mapping.
- **Output:** A dictionary with `tip_count`, `balanced`, `neg_as_zero`, `distance`, `newick` and `provenance`.

**Example code:**

```python
from dnakit.sequence import improve_tree_nni

sequences = {"A": "ACGTACGTAC", "B": "ACGTACGAAC", "C": "TCGTACGTAC", "D": "TCGTTCGTAG"}
result = improve_tree_nni("((A:1,C:1):1,(B:1,D:1):1);", sequences, verbose=False)
print(result["tip_count"], result["balanced"])
print(result["newick"])
```

**Example output:**

```text
4 True
(B:0.1,(C:0.0,D:0.2):0.1,A:0.0);
```

## 18) `tree_path_distances`

- **Purpose:** Compare several trees by their tip-to-tip path lengths.
- **Method:** Calls skbio.tree.path_dists, computed exactly with no sampling. `metric` is euclidean, cityblock or sqeuclidean; `use_length` chooses branch lengths, otherwise each edge counts as 1.
- **API:** `tree_path_distances(newicks[required], ids[optional], metric[optional], use_length[optional], shared_by_all[optional], verbose[optional])`
- **Input:** 2 to 200 Newick strings and an optional list of labels (unique, no spaces).
- **Output:** A dictionary with `ids`, `matrix`, `metric`, `use_length`, `shared_by_all` and `provenance`.

**Example code:**

```python
from dnakit.sequence import tree_path_distances

newicks = ["((A:1,B:1):1,(C:1,D:1):1);", "((A:1,C:1):1,(B:1,D:1):1);"]
result = tree_path_distances(newicks, ids=["t1", "t2"], verbose=False)
print(result["ids"])
print(result["matrix"])
```

**Example output:**

```text
['t1', 't2']
[[0.0, 4.0], [4.0, 0.0]]
```

## 19) `weighted_robinson_foulds_distances`

- **Purpose:** Compare several trees by weighted Robinson-Foulds distances.
- **Method:** Calls skbio.tree.wrf_dists. Each split carries its branch length; `metric` combines the length differences of shared splits with the lengths of unshared ones. `include_tips` counts terminal branches, and `rooted` compares rooted trees.
- **API:** `weighted_robinson_foulds_distances(newicks[required], ids[optional], metric[optional], rooted[optional], include_tips[optional], shared_by_all[optional], verbose[optional])`
- **Input:** 2 to 200 Newick strings and an optional list of labels.
- **Output:** A dictionary with `ids`, `matrix`, `metric`, `rooted`, `include_tips`, `shared_by_all` and `provenance`.

**Example code:**

```python
from dnakit.sequence import weighted_robinson_foulds_distances

newicks = ["((A:1,B:1):1,(C:1,D:1):1);", "((A:1,C:1):1,(B:1,D:1):1);"]
result = weighted_robinson_foulds_distances(newicks, ids=["t1", "t2"], verbose=False)
print(result["metric"], result["matrix"])
```

**Example output:**

```text
cityblock [[0.0, 4.0], [4.0, 0.0]]
```

## 20) `adapter_untrimmed_interval`

- **Purpose:** Return the interval of a read that remains after removing one or more adapters.
- **Method:** Adapters are applied in the given order. Each one is matched only to the part left after the previous adapter, as cutadapt does for chained adapters. The untrimmed interval is `cutadapt.adapters.remainder` of those matches; with no match the whole read remains.
- **API:** `adapter_untrimmed_interval(sequence[required], adapters[required], position[optional], max_errors[optional], min_overlap[optional], verbose[optional])`
- **Input:** One DNA string, 1 to 10 adapter sequences and the adapter position (`3prime`, `5prime` or `anywhere`).
- **Output:** A dictionary with `start` and `stop` (0-based, half-open), `length`, `untrimmed_sequence`, `adapter_matches` (adapter, start, stop and errors for each match) and `provenance`.

**Example code:**

```python
from dnakit.sequencing_qc import adapter_untrimmed_interval

result = adapter_untrimmed_interval("ACGTACGTTTAGATCGGAA", ["AGATCGG"], verbose=False)
print(result["start"], result["stop"], result["untrimmed_sequence"])
match = result["adapter_matches"][0]
print(match["start"], match["stop"], match["errors"])
```

**Example output:**

```text
0 10 ACGTACGTTT
10 17 0
```

## 21) `barcode_from_bam_tag`

- **Purpose:** Extract a UMI, and optionally a cell barcode, from the string tags of one alignment.
- **Method:** An in-memory pysam record is built and passed to `umi_tools.sam_methods.get_barcode_tag`. The UMI comes from `umi_tag` (default RX), cut at `umi_tag_split` and stripped of `umi_tag_delim`; the cell barcode is handled the same way, with `cell_tag_split` defaulting to `-` to drop a 10x GEM suffix. A missing tag raises an error.
- **API:** `barcode_from_bam_tag(tags[required], umi_tag[optional], umi_tag_split[optional], umi_tag_delim[optional], cell_barcode[optional], cell_tag[optional], cell_tag_split[optional], cell_tag_delim[optional], verbose[optional])`
- **Input:** A mapping from two-letter SAM tag names to string values.
- **Output:** A dictionary with `umi`, `cell_barcode` (None unless requested), `umi_tag`, `cell_tag` and `provenance`.

**Example code:**

```python
from dnakit.sequencing_qc import barcode_from_bam_tag

umi = barcode_from_bam_tag({"RX": "ACGTAC-1"}, umi_tag_split="-", verbose=False)
cell = barcode_from_bam_tag(
    {"RX": "ACGTAC-1", "CB": "AAACCT-1"}, cell_barcode=True, cell_tag="CB", verbose=False
)
print(umi["umi"], umi["cell_barcode"])
print(cell["umi"], cell["cell_barcode"])
```

**Example output:**

```text
ACGTAC None
ACGTAC-1 AAACCT
```

## 22) `metacontig_contig_map`

- **Purpose:** Map each gene (metacontig) to the references of a local SAM or BAM file that belong to it.
- **Method:** Calls `umi_tools.sam_methods.getMetaContig2contig`. The gene-transcript file has two tab-separated columns, and only references listed in the alignment header are kept. As in umi_tools, lines starting with `#` are skipped and the first empty line ends the file.
- **API:** `metacontig_contig_map(bam_path[required], gene_transcript_path[required], verbose[optional])`
- **Input:** A local SAM or BAM file and a two-column gene-transcript file.
- **Output:** A dictionary with `metacontig_count`, `metacontigs` (gene to sorted list of references) and `provenance`.

**Example file:** `reads.sam`

```text
@HD	VN:1.6	SO:coordinate
@SQ	SN:tx1	LN:100
@SQ	SN:tx2	LN:100
@SQ	SN:tx3	LN:100
r1	0	tx1	1	60	5M	*	0	0	ACGTA	IIIII	RX:Z:ACGTAC-1
r2	0	tx2	3	60	4M	*	0	0	ACGT	IIII	RX:Z:GGGGTT-1
r3	16	tx2	10	60	3S5M	*	0	0	ACGTACGT	IIIIIIII
```

**Example file:** `gene_to_transcript.tsv`

```text
GENE1	tx1
GENE1	tx2
GENE2	tx3
```

**Example code:**

```python
from dnakit.sequencing_qc import metacontig_contig_map

result = metacontig_contig_map("reads.sam", "gene_to_transcript.tsv", verbose=False)
print(result["metacontig_count"])
for gene, transcripts in result["metacontigs"].items():
    print(gene, transcripts)
```

**Example output:**

```text
2
GENE1 ['tx1', 'tx2']
GENE2 ['tx3']
```

## 23) `tag_reads_by_metacontig`

- **Purpose:** Count the alignments of an indexed BAM file in gene (metacontig) order.
- **Method:** Calls `umi_tools.sam_methods.metafetcher`, which tags each record in memory with `metatag`; no file is written. Counts are per gene, so a reference shared by two genes is counted for both. The BAM needs a `.bai` or `.csi` index.
- **API:** `tag_reads_by_metacontig(bam_path[required], gene_transcript_path[required], metatag[optional], verbose[optional])`
- **Input:** A local indexed BAM file and a gene-transcript file.
- **Output:** A dictionary with `metacontig_count`, `read_count`, `reads_per_metacontig` (records per gene) and `provenance`.

**Example code:**

```python
import tempfile
from pathlib import Path

import pysam

from dnakit.sequencing_qc import tag_reads_by_metacontig

bam = Path(tempfile.mkdtemp()) / "reads.bam"
with pysam.AlignmentFile("reads.sam", "r") as source:
    with pysam.AlignmentFile(str(bam), "wb", template=source) as target:
        for record in source:
            target.write(record)
pysam.index(str(bam))
result = tag_reads_by_metacontig(bam, "gene_to_transcript.tsv", verbose=False)
print(result["metacontig_count"], result["read_count"], result["reads_per_metacontig"])
```

**Example output:**

```text
2 3 {'GENE1': 3, 'GENE2': 0}
```

## 24) `soft_clip_adjusted_position`

- **Purpose:** Compute the UMI deduplication position of an alignment from its CIGAR, and whether it is spliced.
- **Method:** Calls `umi_tools.sam_methods.get_read_position`. For a forward read the position is the start minus any leading soft clip; for a reverse read it is the end plus any trailing soft clip. A CIGAR with an `N`, or a 5' soft clip longer than `soft_clip_threshold` (default 4, the umi_tools default), counts as spliced.
- **API:** `soft_clip_adjusted_position(cigar[required], reference_start[optional], is_reverse[optional], soft_clip_threshold[optional], verbose[optional])`
- **Input:** A CIGAR string, the 0-based reference start and the strand.
- **Output:** A dictionary with `start`, `position`, `is_spliced`, `splice_offset` (offset of the first junction from the read start, or None) and `provenance`.

**Example code:**

```python
from dnakit.sequencing_qc import soft_clip_adjusted_position

fwd = soft_clip_adjusted_position("3S5M2N2M", reference_start=10, verbose=False)
rev = soft_clip_adjusted_position("5M2S", reference_start=10, is_reverse=True, verbose=False)
print(fwd["position"], fwd["splice_offset"])
print(rev["start"], rev["position"], rev["is_spliced"])
```

**Example output:**

```text
7 8
10 17 False
```

## 25) `parse_prefixed_umi_name`

- **Purpose:** Read the UMI, and optionally a cell barcode, from a read name written by umi_tools extract.
- **Method:** Calls `umi_tools.sam_methods.get_barcode_umis`. The name is split at `:`; a `UMI_<sequence>` field gives the UMI and a `CELL_<sequence>` field gives the cell barcode when `cell_barcode=True`. A name without a UMI field raises an error.
- **API:** `parse_prefixed_umi_name(read_name[required], cell_barcode[optional], verbose[optional])`
- **Input:** One read name string.
- **Output:** A dictionary with `umi`, `cell_barcode`, `cell_barcode_requested` and `provenance`.

**Example code:**

```python
from dnakit.sequencing_qc import parse_prefixed_umi_name

result = parse_prefixed_umi_name("HWI:1:UMI_ACGT:CELL_TTAA", cell_barcode=True, verbose=False)
print(result["umi"], result["cell_barcode"])
```

**Example output:**

```text
ACGT TTAA
```

## 26) `resolve_above_knee_barcodes`

- **Purpose:** Discard, or correct, whitelist barcodes that lie within a few edits of a more common barcode.
- **Method:** Calls `umi_tools.whitelist_methods.errorDetectAboveThreshold`. Barcodes are checked from the lowest count upward, and one with a near miss (within `errors` edits) among the more common barcodes is discarded. With `resolution_method="correct"`, a single substitution near miss is recorded in `true_to_false_map`. The input mappings are not modified.
- **API:** `resolve_above_knee_barcodes(cell_barcode_counts[required], cell_whitelist[required], true_to_false_map[required], errors[optional], resolution_method[optional], verbose[optional])`
- **Input:** A barcode count mapping, a whitelist list and a true-to-false barcode mapping.
- **Output:** A dictionary with `kept`, `discarded`, `true_to_false_map`, `resolution_method`, `errors` and `provenance`.

**Example code:**

```python
from dnakit.sequencing_qc import resolve_above_knee_barcodes

result = resolve_above_knee_barcodes(
    {"AAAA": 10, "AAAT": 5}, ["AAAA", "AAAT"], {}, verbose=False
)
print(result["kept"], result["discarded"])
```

**Example output:**

```text
['AAAA'] ['AAAT']
```
