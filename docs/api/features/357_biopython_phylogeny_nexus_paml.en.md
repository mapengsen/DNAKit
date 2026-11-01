# Biopython phylogeny, NEXUS, tree queries, dN/dS and PAML

This group in `dnakit.phylogeny` wraps [Biopython](https://biopython.org) (`Bio.Nexus`, `Bio.Phylo`, `Bio.Align.analysis`, `Bio.codonalign`, `Bio.Phylo.PAML`; tested with Biopython 1.88) and the PAML programs. Alignments are mappings of taxon name to aligned sequence; trees are Newick strings. Biopython is imported only when a function runs; install it with `python -m pip install "dnakit[bio-formats]"`. Functions with randomness take a `seed` and restore the global random state afterwards. `run_paml` needs PAML (`codeml`, `baseml`, `yn00`) on `PATH`, for example `conda install -c bioconda paml`. The examples import from the module `dnakit.phylogeny.biopython_p19_phylo_functions`.

## Example data

`codons.phy` is a 4-sequence, 600-base in-frame coding alignment in PHYLIP format (simulated: one random sequence and three mutated copies). `pair_yn00.out` is the real PAML 4.10.10 `yn00` results file for the pair `a`, `d` of that alignment.

**Example file:** `codons.phy`

```text
 4 600
a  GGATCACAGTCTACACTGCTCACTCCAACCCCGGCCCCTGAGTCCGAGGAGAGGGTGCTTCAGAGTATGTATACCACTGGGGCTGATACGGCGGAGGGCACGTCAATACGGTTCAATGCCCTACTGCATGCTCTTGTGGTTCATCTGCATGGAGAGGGTGGGCATGGGTGGGGGTGCTGGCCCGTGATCTGGACCTCCCATCCACAGCTCATTGTACCGAGTGTAGAGAGGGGCTTGTCCTTCCAGATAGCGTTTCTGTTTCGGTGTAGGTGCGCTTCGACTATGCTACTGCGGTTAACGGGGATGGCAAGTACATTTTTTCGTAGATGTGCCTTGCTAACGAAAGTATTAAACACGTCCCTCACAATAGAATCAGCTTTGGACGCGCGACGGCCGTTCCAGAAAATCTTTGAATACTCAATCCTGCGGGTTCGGGCTCCTAAAACCCATGCTTTGTGTTACCCAGTTCGAGCGCATAGGGAATTCAGGTCCACACATGGCTGGATCCCCATGATATTCAAGAACTATACATTAAGTGCTACCTCCAGAACACATGTTTCAGTCACGGCTTGCCATCATCGATCACGGAATGTAGCATCA
b  GGATCACAGTCTACACTGCTCACTCCAAACCCGGCCCCTGAGTCCGAGGAGAGGGTGCTTCAGAGAATGTATACCACTGGGGCTGATACGGCGGAGGCCACGTCAATACGGTTCAATGCCCTACTGCATGCTCTTGTGGTCCATCTGCATGGAGAGGGTGGGCATGGATGGGGGTGCTGGCCCGTGATCTGGACCTCCCATCCACAGCTCATTGTACCGAGTGTAGAGAGGGGTTTATCCTTCCAGATAGCGTTTCTGTTTCGGTGTAGGTGCGCTTCGACTATGCTACTGCGGTTAACGGGGATGGCAAGTAAATTTTTTCGTGGATCTGCCTTGCTAACGAAAGTATTAAACACGTCCCTCACCGTAGAATCAGCTTTGGACGCGCGACGGCCGTTCCAGAAAATCTTTGAATACTCAATCCTGCGGGTTCGGGCACCTAAATCCCATGCTTTGTGTTACCCAGTTCGAGCGCATAGGGAATTCAGGTCCACACATGGCTGGATCCCGATGATATTCAAGAACTATACATTAAGTGCTACCGCCAGAACACATGTTTCAGTCACGGCTTGCCATCATCGATCACGGAATGTAGCATCA
c  GGATCACAGTTTACACTGCTGACTCCAACCCCGACCCCTGAGTCCGAGGAGAGGGTGCTTCAGCGTATGTATACCACTGGGGCTGATACAACGGAGAGCACGTCAAAACGGTTCAATGCCCTACTGCATGCTCTTGTGGTTCATCTGCATGGAGAGGGTGGGCAAGGGTGGGGGTGCTGGCCCGTGATCTGGACCTCCCATCCACAGCTCATAGTACCGATTGTAGAGAGGGGCTTGTCCTTCCCGAGAGCGTTTCTGTTTCGGTGTAGGTGCGCTTCGACAACGCTAATCCGGTTAACGGGGATGGCAACTACTTTTTATCGTAGATGTGCCTTGCTAACGAAAGTATTAAACACGTCCGTCACAATAGAATCAGCTTTGGACGCGCGATGGCCGTTCCAGGAAACCTTTGAATACTCAATCCTGCGGCTTCGGGCTCCTAAATCCCATGCCTTGTGTTTCCCAGTTCGAGCGCATAGGGTATTCAGGTCCACACATGGCTGTATCCCCATGTTATTCAAGAAATATACATTAAGTGCTACCTCCAGAACACATGTTTCAGTCACGGCTTGCCATCATCGATCACGGAATGTAGCATCA
d  GGATCACAGTCTACACGGCTCACTCCAACCCCGGCCCCTGAGTCCGAGGAGTGGGTGCTTCAAACTATGTATACCACTCGGGCTGATACGTCGGAGGGGACGTTAATACGGTTCCATGCCATAATGCATGCTCTTGTGGTTTATCTGCAAGGAGAGGGTGGGCATGGGTGGGGGTGCTGGCCTGTGATCTGGACCTCCAAGCCACAGCTCACTGTACAGAGTGTAGAGAGGGGCTTGTCCTTCCAGATAGCATTGCTGTTTCGGTGTAGGTGCGCTTCGACTATGCTTCTGCGGTTACCGGGAATGGCACGTGCATTTTTTCGTAGATGCTTCCTGCTAACGAAATTATTTAACACGTCGCTCACAATAGAATCAGCTTTGGAGTCGCGACGGCCGTTCAAGGCTATCTTTGAATACTCAATCCTGTGGGCCCGGCCTCCTATACCCCATGCTTTGTGTTACCCAGTTCCAGCGCATAGCGAATTCAGGTACACCCATGGCGCTATCACCATGATATTCAAGAACTATACAATAAGTGCTATCTCCAGAACACATGTTTCAGCCACAGCTTGCCATCATCGATCACGGAATGTAGCATCA
```

**Example file:** `pair_yn00.out`

```text
YN00         seq.phy

ns =   2	ls = 200

Codon position x base (3x4) table for each sequence.

a                   
position  1:    T:0.23000    C:0.25500    A:0.26000    G:0.25500
position  2:    T:0.28000    C:0.31000    A:0.19500    G:0.21500
position  3:    T:0.25500    C:0.20000    A:0.22500    G:0.32000

d                   
position  1:    T:0.25000    C:0.25500    A:0.25500    G:0.24000
position  2:    T:0.28000    C:0.32000    A:0.19500    G:0.20500
position  3:    T:0.25500    C:0.20000    A:0.23000    G:0.31500

Average
position  1:    T:0.24000    C:0.25500    A:0.25750    G:0.24750
position  2:    T:0.28000    C:0.31500    A:0.19500    G:0.21000
position  3:    T:0.25500    C:0.20000    A:0.22750    G:0.31750

Codon usage for each species
----------------------------------------------------------------------
Phe TTT   5   5 | Ser TCT   1   1 | Tyr TAT   2   3 | Cys TGT   3   2
    TTC   5   6 |     TCC   6   4 |     TAC   2   3 |     TGC   3   4
Leu TTA   3   3 |     TCA   7   6 | *** TAA   0   0 | *** TGA   0   0
    TTG   4   4 |     TCG   1   4 |     TAG   0   0 | Trp TGG   4   5
----------------------------------------------------------------------
Leu CTT   2   3 | Pro CCT   2   4 | His CAT  11   9 | Arg CGT   1   2
    CTC   3   3 |     CCC   2   1 |     CAC   0   0 |     CGC   0   0
    CTA   3   1 |     CCA   3   4 | Gln CAA   0   2 |     CGA   3   2
    CTG   6   5 |     CCG   3   3 |     CAG   5   4 |     CGG   7   8
----------------------------------------------------------------------
Ile ATT   1   0 | Thr ACT   3   5 | Asn AAT   2   1 | Ser AGT   4   2
    ATC   4   5 |     ACC   5   5 |     AAC   2   2 |     AGC   0   1
    ATA   4   7 |     ACA   6   5 | Lys AAA   3   1 | Arg AGA   2   2
Met ATG   4   5 |     ACG   6   4 |     AAG   1   3 |     AGG   5   3
----------------------------------------------------------------------
Val GTT   4   3 | Ala GCT   8   9 | Asp GAT   1   1 | Gly GGT   1   1
    GTC   1   0 |     GCC   3   4 |     GAC   1   0 |     GGC   3   2
    GTA   4   3 |     GCA   2   4 | Glu GAA   3   3 |     GGA   2   3
    GTG   3   3 |     GCG   4   1 |     GAG   6   7 |     GGG   5   4
----------------------------------------------------------------------


Sums
------------------------------------------------------
Phe TTT  10 | Ser TCT   2 | Tyr TAT   5 | Cys TGT   5
    TTC  11 |     TCC  10 |     TAC   5 |     TGC   7
Leu TTA   6 |     TCA  13 | *** TAA   0 | *** TGA   0
    TTG   8 |     TCG   5 |     TAG   0 | Trp TGG   9
------------------------------------------------------
Leu CTT   5 | Pro CCT   6 | His CAT  20 | Arg CGT   3
    CTC   6 |     CCC   3 |     CAC   0 |     CGC   0
    CTA   4 |     CCA   7 | Gln CAA   2 |     CGA   5
    CTG  11 |     CCG   6 |     CAG   9 |     CGG  15
------------------------------------------------------
Ile ATT   1 | Thr ACT   8 | Asn AAT   3 | Ser AGT   6
    ATC   9 |     ACC  10 |     AAC   4 |     AGC   1
    ATA  11 |     ACA  11 | Lys AAA   4 | Arg AGA   4
Met ATG   9 |     ACG  10 |     AAG   4 |     AGG   8
------------------------------------------------------
Val GTT   7 | Ala GCT  17 | Asp GAT   2 | Gly GGT   2
    GTC   1 |     GCC   7 |     GAC   1 |     GGC   5
    GTA   7 |     GCA   6 | Glu GAA   6 |     GGA   5
    GTG   6 |     GCG   5 |     GAG  13 |     GGG   9
------------------------------------------------------



(A) Nei-Gojobori (1986) method



Nei & Gojobori 1986. dN/dS (dN, dS)
(Note: This matrix is not used in later ML. analysis.
Use runmode = -2 for ML pairwise comparison.)

a                   
d                    1.1566 (0.1033 0.0894)


(B) Yang & Nielsen (2000) method

Yang Z, Nielsen R (2000) Estimating synonymous and nonsynonymous substitution rates under realistic evolutionary models. Mol. Biol. Evol. 17:32-43

(equal weighting of pathways)

seq. seq.     S       N        t   kappa   omega     dN +- SE    dS +- SE

   2    1   163.2   436.8   0.2998  1.1043  1.2262 0.1052 +- 0.0164  0.0858 +- 0.0242


(C) LWL85, LPB93 & LWLm methods

Li W.-H., C.-I. Wu, Luo (1985) A new method for estimating synonymous and nonsynonymous rates of nucleotide substitutions considering the relative likelihood of nucleotide and codon changes. Mol. Biol. Evol. 2: 150-174.
Li W-H (1993) Unbiased estimation of the rates of synonymous and nonsynonymous substitution. J. Mol. Evol. 36:96-99
Pamilo P, Bianchi NO (1993) Evolution of the Zfx and Zfy genes - rates and interdependence between the genes. Mol. Biol. Evol. 10:271-281
Yang Z (2006) Computational Molecular Evolution. Oxford University Press, Oxford. Eqs. 2.12 & 2.13

2 (d) vs. 1 (a)

L(i):      378.5     113.0     108.5  sum=    600.0
Ns(i):    9.5000    3.5000    5.0000  sum=  18.0000
Nv(i):   23.5000    9.5000    5.0000  sum=  38.0000
A(i):     0.0264    0.0329    0.0502
B(i):     0.0663    0.0920    0.0483
LWL85:  dS =  0.0986 dN =  0.1002 w = 1.0165 S =  146.2 N =  453.8
LWL85m: dS =  0.0868 dN =  0.1048 w = 1.2080 S =  166.1 N =  433.9 (rho = 0.510)
LPB93:  dS =  0.0897 dN =  0.0986 w = 1.0991
```

## 1) `dnakit.phylogeny.biopython_p19_phylo_functions.nexus_terminal_gaps_to_missing()`

- **Purpose:** Replace the leading and trailing gaps of every aligned sequence by the missing-data symbol.
- **Method:** Writes the alignment as a NEXUS data block, reads it with `Bio.Nexus.Nexus.Nexus` and calls `terminal_gap_to_missing`. Internal gaps stay gaps; with `skip_n=False` terminal runs of N are treated like gaps too.
- **API:** `dnakit.phylogeny.biopython_p19_phylo_functions.nexus_terminal_gaps_to_missing(alignment, *, missing='?', skip_n=True, datatype='dna')`.
- **Input:** A mapping of taxon name to aligned sequence (gap `-`, missing `?`; equal lengths), the missing symbol and `datatype` (`dna`, `rna` or `nucleotide`).
- **Output:** `taxa` and the converted `alignment`.
- **Example code:**

```python
from dnakit.phylogeny.biopython_p19_phylo_functions import nexus_terminal_gaps_to_missing

result = nexus_terminal_gaps_to_missing({"A": "--ACGT--A", "B": "?ACGTTTAA", "C": "-ACG--TTA"})
for name, row in result["alignment"].items():
    print(name, row)
```

- **Example output:**

```text
A ??ACGT--A
B ?ACGTTTAA
C ?ACG--TTA
```

## 2) `dnakit.phylogeny.biopython_p19_phylo_functions.nexus_combine_matrices()`

- **Purpose:** Concatenate several character matrices into one supermatrix.
- **Method:** Reads each alignment as a NEXUS matrix and calls `Bio.Nexus.Nexus.combine`. Taxa that are absent from a data set receive `?` for its columns; one character set per data set records its columns.
- **API:** `dnakit.phylogeny.biopython_p19_phylo_functions.nexus_combine_matrices(matrices, *, datatype='dna')`.
- **Input:** A mapping of data-set name to alignment (taxon name to aligned sequence) and `datatype`.
- **Output:** `taxa`, `character_count`, the combined `alignment`, `charsets` and the `partition` (0-based column indices).
- **Example code:**

```python
from dnakit.phylogeny.biopython_p19_phylo_functions import nexus_combine_matrices

result = nexus_combine_matrices({
    "g1": {"A": "--ACGT--A", "B": "?ACGTTTAA", "C": "-ACG--TTA"},
    "g2": {"A": "ACGTAA", "D": "TTTTGG"},
})
print(result["taxa"], result["character_count"])
for name, row in result["alignment"].items():
    print(name, row)
print(result["charsets"])
```

- **Example output:**

```text
['A', 'B', 'C', 'D'] 15
A --ACGT--AACGTAA
B ?ACGTTTAA??????
C -ACG--TTA??????
D ?????????TTTTGG
{'g1': [0, 1, 2, 3, 4, 5, 6, 7, 8], 'g2': [9, 10, 11, 12, 13, 14]}
```

## 3) `dnakit.phylogeny.biopython_p19_phylo_functions.nexus_write_partition_files()`

- **Purpose:** Write one NEXUS data file for each character partition of an alignment.
- **Method:** Calls `Nexus.write_nexus_data_partitions` with the partition (for example codon positions). Each file `<file_stem>_<partition>.data` holds only the columns of its partition; `exclude` removes columns and `delete` removes taxa. The output directory must not exist.
- **API:** `dnakit.phylogeny.biopython_p19_phylo_functions.nexus_write_partition_files(alignment, charpartition, output_dir, *, file_stem='alignment', exclude=None, delete=None, interleave=False, datatype='dna')`.
- **Input:** An alignment mapping, `charpartition` (partition name to 0-based column indices), a new `output_dir`, `file_stem`, `exclude`, `delete`, `interleave` and `datatype`.
- **Output:** `files` (partition name to path) and `partition_sizes`.
- **Example code:**

```python
import os, tempfile
from dnakit.phylogeny.biopython_p19_phylo_functions import nexus_write_partition_files

target = os.path.join(tempfile.mkdtemp(), "parts")
result = nexus_write_partition_files(
    {"A": "--ACGT--A", "B": "?ACGTTTAA", "C": "-ACG--TTA"},
    {"first": [0, 1, 2], "rest": [3, 4, 5, 6, 7, 8]},
    target,
)
print(result["partition_sizes"])
print(sorted(os.path.basename(path) for path in result["files"].values()))
print(open(result["files"]["first"]).read())
```

- **Example output:**

```text
{'first': 3, 'rest': 6}
['alignment_first.data', 'alignment_rest.data']
#NEXUS
[
Partition: first
]
begin data;
dimensions ntax=3 nchar=3;
format datatype=dna missing=? gap=-;
matrix
A --A
B ?AC
C -AC
;
end;
```

## 4) `dnakit.phylogeny.biopython_p19_phylo_functions.nexus_weighted_stepmatrix()`

- **Purpose:** Compute a weighted-parsimony step matrix from the character states observed in an alignment.
- **Method:** Follows Wheeler (1990): counts how often two states occur together in a site, normalises the counts and takes the negative logarithm (`Bio.Nexus.Nexus.Nexus.weighted_stepmatrix`). A pair that never co-occurs gets `inf`.
- **API:** `dnakit.phylogeny.biopython_p19_phylo_functions.nexus_weighted_stepmatrix(alignment, *, name='your_name_here', exclude=None, delete=None, datatype='dna')`.
- **Input:** An alignment mapping with at least two taxa, the matrix `name`, `exclude` (0-based columns), `delete` (taxa) and `datatype`.
- **Output:** `symbols`, the NEXUS `stepmatrix` text and `weights` (`null` stands for `inf`).
- **Example code:**

```python
from dnakit.phylogeny.biopython_p19_phylo_functions import nexus_weighted_stepmatrix

result = nexus_weighted_stepmatrix({"A": "--ACGT--A", "B": "?ACGTTTAA", "C": "-ACG--TTA"}, name="steps")
print(result["symbols"])
print("\n".join(line.rstrip() for line in result["stepmatrix"].splitlines()))
```

- **Example output:**

```text
['A', 'C', 'G', 'T', '-']
usertype steps stepmatrix=5
        A        C        G        T        -
[A]     .       1.39     inf.     1.39     inf.
[C]    1.39      .       1.39     inf.     inf.
[G]    inf.     1.39      .       1.39     inf.
[T]    1.39     inf.     1.39      .       inf.
[-]    inf.     inf.     inf.     inf.      .
;
```

## 5) `dnakit.phylogeny.biopython_p19_phylo_functions.nexus_bootstrap_alignment()`

- **Purpose:** Resample the columns of an alignment with replacement.
- **Method:** Calls `Bio.Nexus.Nexus.Nexus.bootstrap` after removing excluded columns and deleted taxa. With `seed` the draw is reproducible; the global random state is restored afterwards.
- **API:** `dnakit.phylogeny.biopython_p19_phylo_functions.nexus_bootstrap_alignment(alignment, *, seed=None, exclude=None, delete=None, datatype='dna')`.
- **Input:** An alignment mapping, `seed`, `exclude`, `delete` and `datatype`.
- **Output:** The bootstrapped `alignment` (same length as the input after exclusions).
- **Example code:**

```python
from dnakit.phylogeny.biopython_p19_phylo_functions import nexus_bootstrap_alignment

result = nexus_bootstrap_alignment({"A": "ACGTACGTAC", "B": "ACGTACGTTC", "C": "ACGAACGTTC"}, seed=7)
for name, row in result["alignment"].items():
    print(name, row)
```

- **Example output:**

```text
A CGGACACCCA
B CGGACTCCCA
C CGGACTCCCA
```

## 6) `dnakit.phylogeny.biopython_p19_phylo_functions.biopython_alignment_distance_matrix()`

- **Purpose:** Compute the pairwise distance matrix of a DNA alignment with Biopython's `DistanceCalculator`.
- **Method:** `identity` is the fraction of differing aligned positions; `blastn`, `megablast`, `trans` and `hoxd70` derive the distance from a nucleotide substitution matrix (`Bio.Phylo.TreeConstruction.DistanceCalculator.get_distance`).
- **API:** `dnakit.phylogeny.biopython_p19_phylo_functions.biopython_alignment_distance_matrix(alignment, *, model='identity')`.
- **Input:** An alignment mapping with at least two sequences and `model`.
- **Output:** `ids` and the full symmetric `matrix`.
- **Example code:**

```python
from dnakit.phylogeny.biopython_p19_phylo_functions import biopython_alignment_distance_matrix

aln = {"a": "ACGTACGTAC", "b": "ACGTACGTTC", "c": "ACGAACGTTC", "d": "TCGAACGATC"}
for model in ("identity", "blastn"):
    result = biopython_alignment_distance_matrix(aln, model=model)
    print(model, result["ids"], [round(value, 3) for value in result["matrix"][3]])
```

- **Example output:**

```text
identity ['a', 'b', 'c', 'd'] [0.4, 0.3, 0.2, 0.0]
blastn ['a', 'b', 'c', 'd'] [0.72, 0.54, 0.36, 0.0]
```

## 7) `dnakit.phylogeny.biopython_p19_phylo_functions.parsimony_tree_search()`

- **Purpose:** Search for a maximum-parsimony tree by nearest-neighbour interchange.
- **Method:** Runs Biopython's `ParsimonyTreeConstructor` with `NNITreeSearcher` and `ParsimonyScorer` (Fitch, unit costs). Without a starting tree it starts from a neighbour-joining tree on identity distances.
- **API:** `dnakit.phylogeny.biopython_p19_phylo_functions.parsimony_tree_search(alignment, *, starting_tree=None)`.
- **Input:** An alignment mapping with at least three sequences and an optional Newick `starting_tree` with branch lengths and the same taxa.
- **Output:** The `newick` tree and its `parsimony_score`.
- **Example code:**

```python
from dnakit.phylogeny.biopython_p19_phylo_functions import parsimony_tree_search

result = parsimony_tree_search({"a": "ACGTACGTAC", "b": "ACGTACGTTC", "c": "ACGAACGTTC", "d": "TCGAACGATC"})
print(result["newick"])
print(result["parsimony_score"], result["starting_tree"])
```

- **Example output:**

```text
(d:0.1625,((a:0.075,b:0.05)Inner1:0.025,c:0.05)Inner2:0.0875)Inner3:0;
4.0 neighbour joining (identity distance)
```

## 8) `dnakit.phylogeny.biopython_p19_phylo_functions.strict_consensus_tree()`

- **Purpose:** Build the strict consensus of several trees: the clades present in every tree.
- **Method:** Matches clades by tip name. Biopython's `strict_consensus` compares tip positions and is only correct when all trees list the tips in the same order; this function gives the same tree in that case and the correct tree otherwise.
- **API:** `dnakit.phylogeny.biopython_p19_phylo_functions.strict_consensus_tree(newick_trees)`.
- **Input:** A list of at least two Newick strings with the same taxa.
- **Output:** `newick`, `tree_count` and `shared_clade_count`.
- **Example code:**

```python
from dnakit.phylogeny.biopython_p19_phylo_functions import strict_consensus_tree

result = strict_consensus_tree(["((A,B),(C,D),E);", "((A,B),(C,E),D);"])
print(result["newick"], result["shared_clade_count"])
```

- **Example output:**

```text
((A:0,B:0):0,C:0,D:0,E:0):0; 1
```

## 9) `dnakit.phylogeny.biopython_p19_phylo_functions.collapse_tree_genera()`

- **Purpose:** Collapse every subtree whose tips belong to one genus.
- **Method:** The genus is the first word of the taxon name (`Genus_species`). Uses `Bio.Nexus.Trees.Tree.collapse_genera`; a collapsed subtree becomes one tip named `<Genus> <collapsed>`.
- **API:** `dnakit.phylogeny.biopython_p19_phylo_functions.collapse_tree_genera(newick, *, space_equals_underscore=True)`.
- **Input:** A Newick string and `space_equals_underscore`.
- **Output:** The collapsed `newick` and the tip counts before and after.
- **Example code:**

```python
from dnakit.phylogeny.biopython_p19_phylo_functions import collapse_tree_genera

result = collapse_tree_genera("((Homo_sapiens:1,Homo_erectus:1):1,(Pan_troglodytes:1,Gorilla_gorilla:1):1);")
print(result["newick"])
print(result["terminal_count_before"], result["terminal_count_after"])
```

- **Example output:**

```text
(Homo <collapsed>:1.00000,(Pan_troglodytes:1.00000,Gorilla_gorilla:1.00000)0.00000:1.00000);
4 3
```

## 10) `dnakit.phylogeny.biopython_p19_phylo_functions.split_tree_clade()`

- **Purpose:** Split a named clade into new descendant clades.
- **Method:** Calls `Bio.Phylo.BaseTree.TreeMixin.split`. The new children are named `<target>0`, `<target>1`, ... and get `branch_length`; a tip becomes an internal node.
- **API:** `dnakit.phylogeny.biopython_p19_phylo_functions.split_tree_clade(newick, target, *, count=2, branch_length=1.0)`.
- **Input:** A Newick string, the unique `target` clade name, `count` and `branch_length`.
- **Output:** `new_clades` and the edited `newick`.
- **Example code:**

```python
from dnakit.phylogeny.biopython_p19_phylo_functions import split_tree_clade

result = split_tree_clade("((A:1,B:1):1,C:2);", "C", count=3, branch_length=0.5)
print(result["new_clades"])
print(result["newick"])
```

- **Example output:**

```text
['C0', 'C1', 'C2']
((A:1,B:1):1,(C0:0.5,C1:0.5,C2:0.5)C:2):0;
```

## 11) `dnakit.phylogeny.biopython_p19_phylo_functions.random_tree()`

- **Purpose:** Generate a random bifurcating tree on the given taxa.
- **Method:** Calls `Bio.Phylo.BaseTree.Tree.randomized`: it repeatedly splits a random tip and shuffles the names. This is not a coalescent or birth-death model. With `seed` the tree is reproducible.
- **API:** `dnakit.phylogeny.biopython_p19_phylo_functions.random_tree(taxa, *, branch_length=1.0, branch_stdev=None, seed=None)`.
- **Input:** A list of 2 to 10000 distinct taxon names, `branch_length`, optional `branch_stdev` and `seed`.
- **Output:** The `newick` tree.
- **Example code:**

```python
from dnakit.phylogeny.biopython_p19_phylo_functions import random_tree

print(random_tree(["A", "B", "C", "D", "E"], seed=3)["newick"])
print(random_tree(["A", "B", "C", "D", "E"], seed=3)["newick"])
```

- **Example output:**

```text
(((C:1,(A:1,E:1)n001:1)n00:1,D:1)n0:1,B:1):0;
(((C:1,(A:1,E:1)n001:1)n00:1,D:1)n0:1,B:1):0;
```

## 12) `dnakit.phylogeny.biopython_p19_phylo_functions.tree_is_bifurcating()`

- **Purpose:** Test whether a tree is strictly bifurcating.
- **Method:** Calls `Bio.Phylo.BaseTree.TreeMixin.is_bifurcating`: every internal node has two children; a root with three children (unrooted form) is accepted.
- **API:** `dnakit.phylogeny.biopython_p19_phylo_functions.tree_is_bifurcating(newick)`.
- **Input:** A Newick string.
- **Output:** `is_bifurcating`, `terminal_count` and the set of child counts of the internal nodes.
- **Example code:**

```python
from dnakit.phylogeny.biopython_p19_phylo_functions import tree_is_bifurcating

for text in ("((A,B),(C,D));", "(A,B,C,D);", "((A,B,C),D);"):
    result = tree_is_bifurcating(text)
    print(text, result["is_bifurcating"], result["internal_child_counts"])
```

- **Example output:**

```text
((A,B),(C,D)); True [2]
(A,B,C,D); False [4]
((A,B,C),D); False [2, 3]
```

## 13) `dnakit.phylogeny.biopython_p19_phylo_functions.tree_is_parent_of()`

- **Purpose:** Test whether one named clade contains another.
- **Method:** Calls `Bio.Phylo.BaseTree.TreeMixin.is_parent_of` on the clade `ancestor` with the clade `descendant`. As in Biopython a clade counts as its own descendant.
- **API:** `dnakit.phylogeny.biopython_p19_phylo_functions.tree_is_parent_of(newick, ancestor, descendant)`.
- **Input:** A Newick string and two unique clade names.
- **Output:** `is_parent_of` (true or false).
- **Example code:**

```python
from dnakit.phylogeny.biopython_p19_phylo_functions import tree_is_parent_of

tree = "((A,B)AB,(C,D)CD)root;"
print(tree_is_parent_of(tree, "AB", "A")["is_parent_of"])
print(tree_is_parent_of(tree, "AB", "C")["is_parent_of"])
```

- **Example output:**

```text
True
False
```

## 14) `dnakit.phylogeny.biopython_p19_phylo_functions.find_tree_clades()`

- **Purpose:** Find clades by name, tip status or branch length.
- **Method:** Uses `find_clades` (all matches) or `find_any` (`first_only=True`) of `Bio.Phylo`. `name` is a regular expression that must match the whole name; `terminal` selects tips or internal nodes; `branch_length` is compared by equality within 1e-9.
- **API:** `dnakit.phylogeny.biopython_p19_phylo_functions.find_tree_clades(newick, *, name=None, terminal=None, branch_length=None, order='preorder', first_only=False, max_results=10000)`.
- **Input:** A Newick string, `name`, `terminal`, `branch_length`, `order` (`preorder`, `postorder` or `level`), `first_only` and `max_results`.
- **Output:** `match_count` and a `clades` list with name, branch length, confidence, tip status and tip count.
- **Example code:**

```python
from dnakit.phylogeny.biopython_p19_phylo_functions import find_tree_clades

tree = "((A_x:1,A_y:2)0.9:1,(C:1,D:1)0.8:2,E:3);"
print([c["name"] for c in find_tree_clades(tree, name="A.*")["clades"]])
print(find_tree_clades(tree, terminal=False)["match_count"])
print(find_tree_clades(tree, branch_length=3)["clades"])
```

- **Example output:**

```text
['A_x', 'A_y']
3
[{'name': 'E', 'branch_length': 3.0, 'confidence': None, 'terminal': True, 'terminal_count': 1}]
```

## 15) `dnakit.phylogeny.biopython_p19_phylo_functions.trees_are_identical()`

- **Purpose:** Test whether two trees have the same topology and taxa.
- **Method:** Calls `Bio.Nexus.Trees.Tree.is_identical`, which compares the nested sets of children; branch lengths and the order of children are ignored.
- **API:** `dnakit.phylogeny.biopython_p19_phylo_functions.trees_are_identical(newick_a, newick_b)`.
- **Input:** Two Newick strings.
- **Output:** `identical` (true or false).
- **Example code:**

```python
from dnakit.phylogeny.biopython_p19_phylo_functions import trees_are_identical

print(trees_are_identical("((A,B),(C,D));", "((D,C),(B,A));")["identical"])
print(trees_are_identical("((A,B),(C,D));", "((A,C),(B,D));")["identical"])
```

- **Example output:**

```text
True
False
```

## 16) `dnakit.phylogeny.biopython_p19_phylo_functions.trees_are_compatible()`

- **Purpose:** Check two trees for conflicting supported clades.
- **Method:** Calls `Bio.Nexus.Trees.Tree.is_compatible`: internal branches with support at least `threshold` are compared and every pair of clades that overlap without nesting is a conflict. Support must be written `)support:length`, for example `(A,B)90:1`. With `strict=True` the taxa must be equal.
- **API:** `dnakit.phylogeny.biopython_p19_phylo_functions.trees_are_compatible(newick_a, newick_b, *, threshold, strict=True)`.
- **Input:** Two Newick strings, `threshold` and `strict`.
- **Output:** `compatible`, `conflict_count` and the `conflicts` (clades, supports and shared taxa).
- **Example code:**

```python
from dnakit.phylogeny.biopython_p19_phylo_functions import trees_are_compatible

same = trees_are_compatible("((A,B)90:1,(C,D)90:1,E);", "((A,B)95:1,(C,D)80:1,E);", threshold=50)
print(same["compatible"], same["conflict_count"])
clash = trees_are_compatible("((A,B)90:1,(C,D)90:1,E);", "((A,C)90:1,(B,D)90:1,E);", threshold=50)
print(clash["compatible"], clash["conflict_count"])
print(clash["conflicts"][0])
```

- **Example output:**

```text
True 0
False 4
{'clade_a': ['A', 'B'], 'support_a': 90.0, 'clade_b': ['A', 'C'], 'support_b': 90.0, 'shared': ['A']}
```

## 17) `dnakit.phylogeny.biopython_p19_phylo_functions.convert_tree_support()`

- **Purpose:** Convert clade support values.
- **Method:** `absolute_to_relative` divides clade counts (as written by PHYLIP `consense`) by the number of replicates (`convert_absolute_support`). `branchlength_to_support` moves values that were read as branch lengths (for example PAUP bootstrap trees) into the support field (`branchlength2support`). Both are `Bio.Nexus.Trees.Tree` methods.
- **API:** `dnakit.phylogeny.biopython_p19_phylo_functions.convert_tree_support(newick, *, mode='absolute_to_relative', replicates=None)`.
- **Input:** A Newick string, `mode` and `replicates` (required for `absolute_to_relative`).
- **Output:** `newick` with `support:length` labels and the `nodes` with taxon, support and branch length.
- **Example code:**

```python
from dnakit.phylogeny.biopython_p19_phylo_functions import convert_tree_support

result = convert_tree_support("((A:1,B:2)90:1,(C:1,D:1)80:2,E:3);", replicates=100)
print(result["newick"])
moved = convert_tree_support("((A:0.9,B:2)0.9:1,(C:1,D:1)0.8:2,E:3);", mode="branchlength_to_support")
print(moved["newick"])
```

- **Example output:**

```text
((A:1.00000,B:2.00000)0.90:1.00000,(C:1.00000,D:1.00000)0.80:2.00000,E:3.00000);
((A:0.00000,B:0.00000)1.00:0.00000,(C:0.00000,D:0.00000)2.00:0.00000,E:0.00000);
```

## 18) `dnakit.phylogeny.biopython_p19_phylo_functions.codon_dn_ds_estimates()`

- **Purpose:** Estimate pairwise dN and dS for a codon alignment.
- **Method:** Calls `Bio.Align.analysis.calculate_dn_ds_matrix`. `NG86` is Nei-Gojobori, `LWL85` Li-Wu-Luo, `YN00` Yang-Nielsen and `ML` the maximum-likelihood method of Goldman and Yang. A saturated pair (Biopython returns -1) is reported as `null`. The default `YN00` differs from the `NG86` default of `pairwise_dn_ds_matrix`.
- **API:** `dnakit.phylogeny.biopython_p19_phylo_functions.codon_dn_ds_estimates(codon_alignment, *, method='YN00', genetic_code=1)`.
- **Input:** An in-frame, gap-free coding alignment (A, C, G, T only) as a mapping, `method` and `genetic_code` (NCBI table number).
- **Output:** `ids` and the full symmetric `dn` and `ds` matrices.
- **Example code:**

```python
rows = dict(line.split() for line in open("codons.phy").read().splitlines()[1:])
from dnakit.phylogeny.biopython_p19_phylo_functions import codon_dn_ds_estimates

for method in ("NG86", "LWL85", "YN00"):
    result = codon_dn_ds_estimates(rows, method=method)
    print(method, "dN(a,d) =", round(result["dn"][0][3], 4), "dS(a,d) =", round(result["ds"][0][3], 4))
```

- **Example output:**

```text
NG86 dN(a,d) = 0.1017 dS(a,d) = 0.0936
LWL85 dN(a,d) = 0.0981 dS(a,d) = 0.1017
YN00 dN(a,d) = 0.0992 dS(a,d) = 0.1019
```

## 19) `dnakit.phylogeny.biopython_p19_phylo_functions.codon_dn_ds_trees()`

- **Purpose:** Build a dN tree and a dS tree from a codon alignment.
- **Method:** Calls `Bio.codonalign.codonalignment.CodonAlignment.get_dn_ds_tree`: the pairwise dN and dS matrices (`NG86`, `LWL85`, `YN00` or `ML`) are turned into trees with UPGMA or neighbour joining.
- **API:** `dnakit.phylogeny.biopython_p19_phylo_functions.codon_dn_ds_trees(codon_alignment, *, dn_ds_method='NG86', tree_method='UPGMA', genetic_code=1)`.
- **Input:** A codon alignment mapping as above, `dn_ds_method`, `tree_method` (`UPGMA` or `NJ`) and `genetic_code`.
- **Output:** `dn_tree_newick` and `ds_tree_newick`.
- **Example code:**

```python
rows = dict(line.split() for line in open("codons.phy").read().splitlines()[1:])
from dnakit.phylogeny.biopython_p19_phylo_functions import codon_dn_ds_trees

result = codon_dn_ds_trees(rows, dn_ds_method="NG86", tree_method="UPGMA")
print(result["dn_tree_newick"])
print(result["ds_tree_newick"])
```

- **Example output:**

```text
(d:0.066962759,(c:0.033571855,(b:0.010170115,a:0.010170115)Inner1:0.02340174)Inner2:0.033390904)Inner3:0;
(d:0.066604305,(b:0.03646134,(c:0.02201535,a:0.02201535)Inner1:0.01444599)Inner2:0.030142964)Inner3:0;
```

## 20) `dnakit.phylogeny.biopython_p19_phylo_functions.run_paml()`

- **Purpose:** Run the PAML programs `codeml`, `baseml` or `yn00` on a sequence file.
- **Method:** Copies the alignment (and tree) into a new working directory, writes the control file from Biopython's option table (`Bio.Phylo.PAML.<program>`), runs the program with the control file as its only argument and parses the results file. A bare Newick tree gets the `<taxa> 1` header line that PAML expects. Needs PAML on `PATH` (`conda install -c bioconda paml`) or `executable`.
- **API:** `dnakit.phylogeny.biopython_p19_phylo_functions.run_paml(program, alignment_path, working_dir, *, tree_path=None, options=None, executable=None, timeout_seconds=3600.0)`.
- **Input:** `program`, the PHYLIP/PAML `alignment_path`, a new `working_dir`, an optional `tree_path` (not for `yn00`), `options` by PAML name (unknown names are refused), `executable` and `timeout_seconds`.
- **Output:** `status`, the program `version`, `return_code`, `results_file`, the working-directory `files` and the parsed `results`.
- **Example code:**

```python
import os, tempfile
from dnakit.phylogeny.biopython_p19_phylo_functions import run_paml

work = os.path.join(tempfile.mkdtemp(), "yn00")
result = run_paml("yn00", "codons.phy", work, options={"icode": 0, "verbose": 0})
print(result["status"], result["version"], result["return_code"])
print(sorted(result["results"]))
print(result["results"]["a"]["b"]["YN00"])
```

- **Example output:**

```text
completed 4.10.10 0
['a', 'b', 'c', 'd']
{'S': 165.2, 'N': 434.8, 't': 0.0817, 'kappa': 1.1311, 'omega': 0.4805, 'dN': 0.021, 'dN SE': 0.007, 'dS': 0.0437, 'dS SE': 0.0167}
```

## 21) `dnakit.phylogeny.biopython_p19_phylo_functions.read_paml_results()`

- **Purpose:** Parse a PAML results file of `codeml`, `baseml` or `yn00`.
- **Method:** Calls `Bio.Phylo.PAML.<program>.read` (`codeml.read`, `baseml.read`, `yn00.read`) and converts the result to plain JSON types. `yn00` gives the Nei-Gojobori, Yang-Nielsen, LWL85 and LPB93 tables per pair, `codeml` the site-model results and `baseml` the model parameters and tree lengths.
- **API:** `dnakit.phylogeny.biopython_p19_phylo_functions.read_paml_results(program, results_file)`.
- **Input:** `program` and the path of a results file.
- **Output:** The parsed `results`.
- **Example code:**

```python
from dnakit.phylogeny.biopython_p19_phylo_functions import read_paml_results

result = read_paml_results("yn00", "pair_yn00.out")
print(result["results"]["a"]["d"]["YN00"])
print(result["results"]["a"]["d"]["NG86"])
```

- **Example output:**

```text
{'S': 163.2, 'N': 436.8, 't': 0.2998, 'kappa': 1.1043, 'omega': 1.2262, 'dN': 0.1052, 'dN SE': 0.0164, 'dS': 0.0858, 'dS SE': 0.0242}
{'omega': 1.1566, 'dN': 0.1033, 'dS': 0.0894}
```
