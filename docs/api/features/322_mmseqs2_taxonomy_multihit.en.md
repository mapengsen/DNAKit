# MMseqs2 taxonomy assignment and multi-hit search

These functions in `dnakit.alignment` wrap the command-line program [MMseqs2](https://github.com/soedinglab/MMseqs2): each runs one MMseqs2 module with a fixed argument list and no shell. The user installs MMseqs2 (conda or a release binary) and puts it on PATH or passes `mmseqs_path`; the tested build is commit 76de5d8ae114750adc59288516d947b6a72a4b05. Databases are path prefixes; options left unset are not passed, so MMseqs2's defaults apply; new outputs must not exist yet. The examples import by module path.

## Example data

**Example file:** `proteins.fa`

```text
>fam0_m0
MMFPCDVENWCTHCDQQDIDVQCWEIWCWWPCICVFLQFVEWLVGEWWHNEVDWCYHSVQ
>fam0_m1
MMFPCDVLNWCTHCDQQDDDVQCWEIWCWWPLICVFLQRVEWLVGEWWHNEVDWCYHSVF
>fam0_m2
MMFPCQVENWCTHCDQQDIDVQCWEIWCWWPCIHVFLQFVEWLVEEWWHNEVDWCYHSKQ
>fam1_m0
MGNIVVTMIYHIPIHTSNAAKSKHYNRNNDIEISHMHSYYASNDEPHSGQMDPRPDGGFA
>fam1_m1
MGNIVVTMIYEIPIHTSNAAKSKHYNRNNDIEISRMHSYYASTDEPHHGQMDPRPDGGFA
>fam1_m2
MGNIVVTMIYHIPIHTSNAAKSKHYNSNNDIEISHNVSYYASNDEPHKGQMDPRPDGGFC
>fam2_m0
MGQDKADKDYIDKERAMVQKYFCTIEGKCGHLLTHLRTGKNAKCAATVHTSIREQSVPTL
>fam2_m1
MGQDKADKFKIPKERAMVQKYFCTIEGKPKHLLKHLRTGKTAKCAATVHTSIREQSVPTL
>fam2_m2
MGWDKADCDYIDKERAMDQKYFSTKEGKCGHLLTHFRTGKNAKCHATVHTSIREQRTPTL
>fam3_m0
MHHDWDFTKNFYTKENISSPAGASRPLFQNPMEMAMMPEHALKNDPPWDNQKCKECLFIK
>fam3_m1
MHHDWPFTKNQYTKEVISSPAGASRPLFQNPMEMAMMPMHALKNQPPWDNQKCKECLIIK
>fam3_m2
MHHDWDFTKNTYTKENISSPAGASFPLFQNPMEMTMMPEVALKNDPPWDNQSCKICLFMK
```

**Example file:** `contigs.fa`

```text
>contig0
TTTCCTCATGCAATTCAAAACCATGTCCGTAATGTAGGCGATGATGTTTCCTTGTGATGTTGAAAATTGGTGTACTCATTGTGATCAACAAGATATTGATGTTCAATGTTGGGAAATTTGGTGTTGGTGGCCTTGTATTTGTGTTTTTCTGCAATTTGTTGAATGGCTGGTTGGTGAATGGTGGCATAATGAAGTTGATTGGTGTTATCATTCTGTTCAATAAAAATAGTAAACCATTTTACGGAGGATACCAAATTCCTCCTTTAAGCAAAACCACCATCAGGACGAGGATCCATTTGACCAGAATGAGGTTCATCATTAGAAGCATAATAAGAATGCATATGAGAAATTTCAATATCATTATTACGATTATAATGTTTAGATTTAGCAGCATTAGAAGTATGAATAGGAATATGATAAATCATAGTAACAACAATATTACCCATTATTCAGGACCTAACCTGAGGTAAACCAGGTCTCTCCGCC
>contig1
CCCTTATAAAAGCTGTTGCACCTAGCCAAGTTCAACGGCAATGGGTCAAGATAAAGCTGATAAAGATTATATTGATAAAGAACGTGCTATGGTTCAAAAATATTTTTGTACTATTGAAGGTAAATGTGGTCATCTGCTGACTCATCTGCGTACTGGTAAAAATGCTAAATGTGCTGCTACTGTTCATACTTCTATTCGTGAACAATCTGTTCCTACTCTGTAAGCTGCAATGGAAATAGGCAATGACGGATATATATTAAAAATTATTTAATAAACAGACATTCTTTACATTTTTGATTATCCCAAGGAGGATCATTTTTCAGAGCATGTTCAGGCATCATAGCCATTTCCATAGGATTTTGAAACAGAGGACGAGAAGCACCAGCAGGAGAAGAAATATTTTCTTTAGTATAAAAATTTTTAGTAAAATCCCAATCATGATGCATGTGTTTTAAGATACATTGAGGCCCGTTCGTGCTCCTCGCC
>contig2
CTGAAGCATTGCTTTGTGAAGAGGGACTTCAGCCAATAGAATGATGTTTCCTTGTGATGTTCTGAATTGGTGTACTCATTGTGATCAACAAGATGATGATGTTCAATGTTGGGAAATTTGGTGTTGGTGGCCTCTGATTTGTGTTTTTCTGCAACGTGTTGAATGGCTGGTTGGTGAATGGTGGCATAATGAAGTTGATTGGTGTTATCATTCTGTTTTTTAACCTGCATACCGGCTCATTCTTCATGTGCAACCTAGGGAGATTAAGCAAAACCACCATCAGGACGAGGATCCATTTGACCATGATGAGGTTCATCAGTAGAAGCATAATAAGAATGCATACGAGAAATTTCAATATCATTATTACGATTATAATGTTTAGATTTAGCAGCATTAGAAGTATGAATAGGAATTTCATAAATCATAGTAACAACAATATTACCCATATGTGTACATACGCTCTTACTGCGGTCGCGTCTAATAATA
>contig3
TACATTTGCTTCGTTGACTAGCAACCCAGGGCTATAGCTAATGGGTCAAGATAAAGCTGATAAATTTAAAATTCCTAAAGAACGTGCTATGGTTCAAAAATATTTTTGTACTATTGAAGGTAAACCTAAACATCTGCTGAAACATCTGCGTACTGGTAAAACTGCTAAATGTGCTGCTACTGTTCATACTTCTATTCGTGAACAATCTGTTCCTACTCTGTAATTCCCCCCGCGGCCCACCCAGTATTCCTAACGGAGCATAATTATTTAATAATCAGACATTCTTTACATTTTTGATTATCCCAAGGAGGTTGATTTTTCAGAGCATGCATAGGCATCATAGCCATTTCCATAGGATTTTGAAACAGAGGACGAGAAGCACCAGCAGGAGAAGAAATAACTTCTTTAGTATATTGATTTTTAGTAAAAGGCCAATCATGATGCATATCCCACCCGAACTAAGTTTGTCGAACCTTGGTCCAAGAT
```

**Example file:** `taxmap.tsv`

```text
fam0_m0	561
fam0_m1	561
fam0_m2	561
fam1_m0	561
fam1_m1	561
fam1_m2	561
fam2_m0	9605
fam2_m1	9605
fam2_m2	9605
fam3_m0	9605
fam3_m1	9605
fam3_m2	9605
```

**Example file:** `nodes.dmp`

```text
1	|	1	|	no rank	|		|
2	|	1	|	superkingdom	|		|
2759	|	1	|	superkingdom	|		|
561	|	2	|	genus	|		|
9605	|	2759	|	genus	|		|
```

**Example file:** `names.dmp`

```text
1	|	root	|		|	scientific name	|
2	|	Bacteria	|		|	scientific name	|
2759	|	Eukaryota	|		|	scientific name	|
561	|	Escherichia	|		|	scientific name	|
9605	|	Homo	|		|	scientific name	|
```

**Example file:** `merged.dmp`

```text
9999	|	561	|
```

**Example file:** `set_a.fa`

```text
>contig0
TTTCCTCATGCAATTCAAAACCATGTCCGTAATGTAGGCGATGATGTTTCCTTGTGATGTTGAAAATTGGTGTACTCATTGTGATCAACAAGATATTGATGTTCAATGTTGGGAAATTTGGTGTTGGTGGCCTTGTATTTGTGTTTTTCTGCAATTTGTTGAATGGCTGGTTGGTGAATGGTGGCATAATGAAGTTGATTGGTGTTATCATTCTGTTCAATAAAAATAGTAAACCATTTTACGGAGGATACCAAATTCCTCCTTTAAGCAAAACCACCATCAGGACGAGGATCCATTTGACCAGAATGAGGTTCATCATTAGAAGCATAATAAGAATGCATATGAGAAATTTCAATATCATTATTACGATTATAATGTTTAGATTTAGCAGCATTAGAAGTATGAATAGGAATATGATAAATCATAGTAACAACAATATTACCCATTATTCAGGACCTAACCTGAGGTAAACCAGGTCTCTCCGCC
>contig1
CCCTTATAAAAGCTGTTGCACCTAGCCAAGTTCAACGGCAATGGGTCAAGATAAAGCTGATAAAGATTATATTGATAAAGAACGTGCTATGGTTCAAAAATATTTTTGTACTATTGAAGGTAAATGTGGTCATCTGCTGACTCATCTGCGTACTGGTAAAAATGCTAAATGTGCTGCTACTGTTCATACTTCTATTCGTGAACAATCTGTTCCTACTCTGTAAGCTGCAATGGAAATAGGCAATGACGGATATATATTAAAAATTATTTAATAAACAGACATTCTTTACATTTTTGATTATCCCAAGGAGGATCATTTTTCAGAGCATGTTCAGGCATCATAGCCATTTCCATAGGATTTTGAAACAGAGGACGAGAAGCACCAGCAGGAGAAGAAATATTTTCTTTAGTATAAAAATTTTTAGTAAAATCCCAATCATGATGCATGTGTTTTAAGATACATTGAGGCCCGTTCGTGCTCCTCGCC
```

**Example file:** `set_b.fa`

```text
>contig2
CTGAAGCATTGCTTTGTGAAGAGGGACTTCAGCCAATAGAATGATGTTTCCTTGTGATGTTCTGAATTGGTGTACTCATTGTGATCAACAAGATGATGATGTTCAATGTTGGGAAATTTGGTGTTGGTGGCCTCTGATTTGTGTTTTTCTGCAACGTGTTGAATGGCTGGTTGGTGAATGGTGGCATAATGAAGTTGATTGGTGTTATCATTCTGTTTTTTAACCTGCATACCGGCTCATTCTTCATGTGCAACCTAGGGAGATTAAGCAAAACCACCATCAGGACGAGGATCCATTTGACCATGATGAGGTTCATCAGTAGAAGCATAATAAGAATGCATACGAGAAATTTCAATATCATTATTACGATTATAATGTTTAGATTTAGCAGCATTAGAAGTATGAATAGGAATTTCATAAATCATAGTAACAACAATATTACCCATATGTGTACATACGCTCTTACTGCGGTCGCGTCTAATAATA
>contig3
TACATTTGCTTCGTTGACTAGCAACCCAGGGCTATAGCTAATGGGTCAAGATAAAGCTGATAAATTTAAAATTCCTAAAGAACGTGCTATGGTTCAAAAATATTTTTGTACTATTGAAGGTAAACCTAAACATCTGCTGAAACATCTGCGTACTGGTAAAACTGCTAAATGTGCTGCTACTGTTCATACTTCTATTCGTGAACAATCTGTTCCTACTCTGTAATTCCCCCCGCGGCCCACCCAGTATTCCTAACGGAGCATAATTATTTAATAATCAGACATTCTTTACATTTTTGATTATCCCAAGGAGGTTGATTTTTCAGAGCATGCATAGGCATCATAGCCATTTCCATAGGATTTTGAAACAGAGGACGAGAAGCACCAGCAGGAGAAGAAATAACTTCTTTAGTATATTGATTTTTAGTAAAAGGCCAATCATGATGCATATCCCACCCGAACTAAGTTTGTCGAACCTTGGTCCAAGAT
```

## 1) `dnakit.alignment.mmseqs_createtaxdb()`

- **Purpose:** Annotate a sequence database with taxonomy labels (``mmseqs createtaxdb``).
- **Method:** Runs the MMseqs2 module with a fixed argument list and no shell; options you do not set are not passed, so MMseqs2's defaults apply; output prefixes and files must not exist yet.
- **API:** `dnakit.alignment.mmseqs_createtaxdb(sequence_db, *, ncbi_tax_dump=None, tax_mapping_file=None, tax_mapping_mode=None, tax_db_mode=None, threads=1, mmseqs_path=None, verbose=True)`.
- **Input:** `sequence_db`: prefix of the sequence database (modified in place by adding the taxonomy files)..
- **Output:** A dict with `status`, `version`, `command`, `outputs` (new files or database prefixes), `entry_counts` or `line_counts`, `log_tail` and `provenance`; some functions add `records`, `rows` or `clusters`.
- **Example code:**

```python
import shutil
from pathlib import Path

import dnakit.alignment.mmseqs2_p07_align_cli_functions as al
import dnakit.alignment.mmseqs2_p07_taxonomy_cli_functions as tx
import dnakit.alignment.mmseqs2_p07_workflow_cli_functions as wf

shutil.rmtree("work", ignore_errors=True)
Path("work").mkdir()
wf.mmseqs_createdb("proteins.fa", "work/ref", verbose=False)
wf.mmseqs_createdb("contigs.fa", "work/ctg", verbose=False)
r = tx.mmseqs_createtaxdb("work/ref", ncbi_tax_dump=".", tax_mapping_file="taxmap.tsv", verbose=False)
al.mmseqs_extractorfs("work/ctg", "work/orfs", translate=True, min_length=30, verbose=False)
al.mmseqs_prefilter("work/orfs", "work/ref", "work/orfpf", verbose=False)
al.mmseqs_rescorediagonal("work/orfs", "work/ref", "work/orfpf", "work/orfaln", rescore_mode=2, add_backtrace=True, verbose=False)
print([Path(p).name for p in r["taxonomy_files"]])
```

- **Example output:**

```text
['ref_mapping', 'ref_taxonomy']
```

## 2) `dnakit.alignment.mmseqs_taxonomy()`

- **Purpose:** Taxonomic classification of sequences (``mmseqs easy-taxonomy``, ``mmseqs taxonomy``).
- **Method:** Runs the MMseqs2 module with a fixed argument list and no shell; options you do not set are not passed, so MMseqs2's defaults apply; output prefixes and files must not exist yet.
- **API:** `dnakit.alignment.mmseqs_taxonomy(query, target_db, output, *, input_kind='fasta', lca_mode=None, majority=None, vote_mode=None, ..., threads=1, mmseqs_path=None, verbose=True)`.
- **Input:** `query`: see the docstring; `target_db`: see the docstring; `output`: see the docstring.
- **Output:** A dict with `status`, `version`, `command`, `outputs` (new files or database prefixes), `entry_counts` or `line_counts`, `log_tail` and `provenance`; some functions add `records`, `rows` or `clusters`.
- **Example code:**

```python
import dnakit.alignment.mmseqs2_p07_taxonomy_cli_functions as tx

r = tx.mmseqs_taxonomy("contigs.fa", "work/ref", "work/tax", verbose=False)
print([(row["query"], row["name"]) for row in r["records"]])
d = tx.mmseqs_taxonomy("work/orfs", "work/ref", "work/taxdb", input_kind="database", tax_output_mode=2, verbose=False)
print(d["entry_counts"])
```

- **Example output:**

```text
[('contig0', 'Escherichia'), ('contig1', 'Homo'), ('contig2', 'Escherichia'), ('contig3', 'unclassified')]
{'taxonomy_db': 30}
```

## 3) `dnakit.alignment.mmseqs_taxonomyreport()`

- **Purpose:** Create a Kraken-style or Krona taxonomy report from a taxonomy result (`mmseqs taxonomyreport`).
- **Method:** Runs the MMseqs2 module with a fixed argument list and no shell; options you do not set are not passed, so MMseqs2's defaults apply; output prefixes and files must not exist yet.
- **API:** `dnakit.alignment.mmseqs_taxonomyreport(sequence_db, result_db, output, *, report_mode=None, threads=1, mmseqs_path=None, verbose=True)`.
- **Input:** `sequence_db`: taxonomy-annotated sequence database.; `result_db`: taxonomy result (or any result/sequence) database.; `output`: new report file..
- **Output:** A dict with `status`, `version`, `command`, `outputs` (new files or database prefixes), `entry_counts` or `line_counts`, `log_tail` and `provenance`; some functions add `records`, `rows` or `clusters`.
- **Example code:**

```python
import dnakit.alignment.mmseqs2_p07_taxonomy_cli_functions as tx

r = tx.mmseqs_taxonomyreport("work/ref", "work/taxdb", "work/report.txt", verbose=False)
for row in r["rows"]:
    print(row[3:])
```

- **Example output:**

```text
['no rank', '0', 'unclassified']
['no rank', '1', 'root']
['superkingdom', '2', '  Bacteria']
['genus', '561', '    Escherichia']
['superkingdom', '2759', '  Eukaryota']
['genus', '9605', '    Homo']
```

## 4) `dnakit.alignment.mmseqs_addtaxonomy()`

- **Purpose:** Add taxonomic labels to a result database (`mmseqs addtaxonomy`).
- **Method:** Runs the MMseqs2 module with a fixed argument list and no shell; options you do not set are not passed, so MMseqs2's defaults apply; output prefixes and files must not exist yet.
- **API:** `dnakit.alignment.mmseqs_addtaxonomy(target_db, result_db, annotated_db, *, tax_lineage=None, lca_ranks=None, pick_id_from=None, threads=1, mmseqs_path=None, verbose=True)`.
- **Input:** `target_db`: taxonomy-annotated target sequence database.; `result_db`: result (alignment) database.; `annotated_db`: new result database with taxon columns..
- **Output:** A dict with `status`, `version`, `command`, `outputs` (new files or database prefixes), `entry_counts` or `line_counts`, `log_tail` and `provenance`; some functions add `records`, `rows` or `clusters`.
- **Example code:**

```python
import dnakit.alignment.mmseqs2_p07_taxonomy_cli_functions as tx

r = tx.mmseqs_addtaxonomy("work/ref", "work/orfaln", "work/annotated", tax_lineage=1, verbose=False)
print(r["entry_counts"])
```

- **Example output:**

```text
{'annotated_db': 8}
```

## 5) `dnakit.alignment.mmseqs_filtertaxdb()`

- **Purpose:** Keep the entries of a taxonomy result whose taxon matches a taxon expression (`mmseqs filtertaxdb`).
- **Method:** Runs the MMseqs2 module with a fixed argument list and no shell; options you do not set are not passed, so MMseqs2's defaults apply; output prefixes and files must not exist yet.
- **API:** `dnakit.alignment.mmseqs_filtertaxdb(target_db, taxonomy_db, filtered_db, *, taxon_list=None, threads=1, mmseqs_path=None, verbose=True)`.
- **Input:** `target_db`: taxonomy-annotated target database.; `taxonomy_db`: taxonomy result database.; `filtered_db`: new taxonomy result database..
- **Output:** A dict with `status`, `version`, `command`, `outputs` (new files or database prefixes), `entry_counts` or `line_counts`, `log_tail` and `provenance`; some functions add `records`, `rows` or `clusters`.
- **Example code:**

```python
import dnakit.alignment.mmseqs2_p07_taxonomy_cli_functions as tx

r = tx.mmseqs_filtertaxdb("work/ref", "work/taxdb", "work/only_561", taxon_list="561", verbose=False)
print(r["entry_counts"])
```

- **Example output:**

```text
{'filtered_db': 30}
```

## 6) `dnakit.alignment.mmseqs_filtertaxseqdb()`

- **Purpose:** Keep the sequences of a taxonomy-annotated database that belong to a taxon expression (`mmseqs filtertaxseqdb`).
- **Method:** Runs the MMseqs2 module with a fixed argument list and no shell; options you do not set are not passed, so MMseqs2's defaults apply; output prefixes and files must not exist yet.
- **API:** `dnakit.alignment.mmseqs_filtertaxseqdb(sequence_db, filtered_db, *, taxon_list=None, subdb_mode=None, threads=1, mmseqs_path=None, verbose=True)`.
- **Input:** `sequence_db`: taxonomy-annotated sequence database.; `filtered_db`: new sequence database..
- **Output:** A dict with `status`, `version`, `command`, `outputs` (new files or database prefixes), `entry_counts` or `line_counts`, `log_tail` and `provenance`; some functions add `records`, `rows` or `clusters`.
- **Example code:**

```python
import dnakit.alignment.mmseqs2_p07_taxonomy_cli_functions as tx

r = tx.mmseqs_filtertaxseqdb("work/ref", "work/only_homo", taxon_list="9605", verbose=False)
print(r["entry_counts"])
```

- **Example output:**

```text
{'filtered_db': 6}
```

## 7) `dnakit.alignment.mmseqs_majoritylca()`

- **Purpose:** Lowest common ancestor of the hits of every query by majority voting (`mmseqs majoritylca`).
- **Method:** Runs the MMseqs2 module with a fixed argument list and no shell; options you do not set are not passed, so MMseqs2's defaults apply; output prefixes and files must not exist yet.
- **API:** `dnakit.alignment.mmseqs_majoritylca(target_db, result_db, taxonomy_db, *, majority=None, vote_mode=None, lca_ranks=None, blacklist=None, ..., threads=1, mmseqs_path=None, verbose=True)`.
- **Input:** `target_db`: taxonomy-annotated target sequence database.; `result_db`: alignment or prefilter result database.; `taxonomy_db`: new taxonomy result database..
- **Output:** A dict with `status`, `version`, `command`, `outputs` (new files or database prefixes), `entry_counts` or `line_counts`, `log_tail` and `provenance`; some functions add `records`, `rows` or `clusters`.
- **Example code:**

```python
import dnakit.alignment.mmseqs2_p07_taxonomy_cli_functions as tx

r = tx.mmseqs_majoritylca("work/ref", "work/orfaln", "work/mlca", majority=0.6, verbose=False)
print(r["entry_counts"])
```

- **Example output:**

```text
{'taxonomy_db': 30}
```

## 8) `dnakit.alignment.mmseqs_lcaalign()`

- **Purpose:** Gapped alignment tailored to the lowest-common-ancestor computation (`mmseqs lcaalign`).
- **Method:** Runs the MMseqs2 module with a fixed argument list and no shell; options you do not set are not passed, so MMseqs2's defaults apply; output prefixes and files must not exist yet.
- **API:** `dnakit.alignment.mmseqs_lcaalign(query_db, target_db, result_db, alignment_db, *, e_value=None, min_seq_id=None, coverage=None, cov_mode=None, ..., threads=1, mmseqs_path=None, verbose=True)`.
- **Input:** `query_db`: query sequence database.; `target_db`: target sequence database.; `result_db`: prefilter result database.; `alignment_db`: new alignment database..
- **Output:** A dict with `status`, `version`, `command`, `outputs` (new files or database prefixes), `entry_counts` or `line_counts`, `log_tail` and `provenance`; some functions add `records`, `rows` or `clusters`.
- **Example code:**

```python
import dnakit.alignment.mmseqs2_p07_taxonomy_cli_functions as tx

r = tx.mmseqs_lcaalign("work/orfs", "work/ref", "work/orfpf", "work/lcaaln", add_backtrace=True, verbose=False)
print(r["entry_counts"])
```

- **Example output:**

```text
{'alignment_db': 30}
```

## 9) `dnakit.alignment.mmseqs_aggregate_taxonomy()`

- **Purpose:** Aggregate taxon labels of a set of sequences into one (``mmseqs aggregatetax``, ``aggregatetaxweights``).
- **Method:** Runs the MMseqs2 module with a fixed argument list and no shell; options you do not set are not passed, so MMseqs2's defaults apply; output prefixes and files must not exist yet.
- **API:** `dnakit.alignment.mmseqs_aggregate_taxonomy(sequence_db, set_to_sequence_db, taxonomy_per_sequence_db, output, *, alignment_per_sequence_db=None, swap_headers=False, majority=None, vote_mode=None, ..., threads=1, mmseqs_path=None, verbose=True)`.
- **Input:** `sequence_db`: taxonomy-annotated target sequence database.; `set_to_sequence_db`: database that maps every set key to its member keys. For ORFs of contigs it is the transposed header database of the ORF database; pass ``<orf_db>_h`` with ``swap_headers=True`` and the facade transposes it (``mmseqs swapdb``) first.; `taxonomy_per_sequence_db`: taxonomy result of the sequences.; `output`: prefix of the new taxonomy result database per set..
- **Output:** A dict with `status`, `version`, `command`, `outputs` (new files or database prefixes), `entry_counts` or `line_counts`, `log_tail` and `provenance`; some functions add `records`, `rows` or `clusters`.
- **Example code:**

```python
import dnakit.alignment.mmseqs2_p07_taxonomy_cli_functions as tx
import dnakit.alignment.mmseqs2_p07_workflow_cli_functions as wf

r = tx.mmseqs_aggregate_taxonomy("work/ref", "work/orfs_h", "work/taxdb", "work/per_contig", swap_headers=True, verbose=False)
w = tx.mmseqs_aggregate_taxonomy("work/ref", "work/orfs_h", "work/taxdb", "work/per_contig_w", alignment_per_sequence_db="work/taxdb_aln", swap_headers=True, majority=0.6, verbose=False)
t = wf.mmseqs_createtsv("work/ctg", "work/per_contig", "work/per_contig.tsv", target_column=0, verbose=False)
print(r["entry_counts"], w["entry_counts"])
for row in t["rows"]:
    print(row[0], row[2], row[4])
```

- **Example output:**

```text
{'taxonomy_db': 4} {'taxonomy_db': 4}
contig0 561 Escherichia
contig1 9605 Homo
contig2 561 Escherichia
contig3 0 unclassified
```

## 10) `dnakit.alignment.mmseqs_multihitdb()`

- **Purpose:** Build a grouped sequence database (sets of ORFs per contig) for multi-hit search (`mmseqs multihitdb`).
- **Method:** Runs the MMseqs2 module with a fixed argument list and no shell; options you do not set are not passed, so MMseqs2's defaults apply; output prefixes and files must not exist yet.
- **API:** `dnakit.alignment.mmseqs_multihitdb(fasta_files, set_db, *, dbtype=None, shuffle=None, createdb_mode=None, min_length=None, ..., threads=1, mmseqs_path=None, verbose=True)`.
- **Input:** `fasta_files`: nucleotide FASTA file(s); every file becomes one set (for example a genome), its ORFs are the members.; `set_db`: prefix of the new set database..
- **Output:** A dict with `status`, `version`, `command`, `outputs` (new files or database prefixes), `entry_counts` or `line_counts`, `log_tail` and `provenance`; some functions add `records`, `rows` or `clusters`.
- **Example code:**

```python
import dnakit.alignment.mmseqs2_p07_taxonomy_cli_functions as tx

a = tx.mmseqs_multihitdb(["set_a.fa", "set_b.fa"], "work/sets1", min_length=30, verbose=False)
tx.mmseqs_multihitdb(["set_b.fa", "set_a.fa"], "work/sets2", min_length=30, verbose=False)
print(a["entry_counts"])
```

- **Example output:**

```text
{'set_db': 30}
```

## 11) `dnakit.alignment.mmseqs_multihitsearch()`

- **Purpose:** Search a grouped query set database against a grouped target set database (`mmseqs multihitsearch`).
- **Method:** Runs the MMseqs2 module with a fixed argument list and no shell; options you do not set are not passed, so MMseqs2's defaults apply; output prefixes and files must not exist yet.
- **API:** `dnakit.alignment.mmseqs_multihitsearch(query_set_db, target_set_db, result_db, *, e_value=None, min_seq_id=None, coverage=None, cov_mode=None, ..., threads=1, mmseqs_path=None, verbose=True)`.
- **Input:** `query_set_db`: query set database from multihitdb.; `target_set_db`: target set database from multihitdb.; `result_db`: new result database..
- **Output:** A dict with `status`, `version`, `command`, `outputs` (new files or database prefixes), `entry_counts` or `line_counts`, `log_tail` and `provenance`; some functions add `records`, `rows` or `clusters`.
- **Example code:**

```python
import dnakit.alignment.mmseqs2_p07_taxonomy_cli_functions as tx

r = tx.mmseqs_multihitsearch("work/sets1", "work/sets2", "work/hits", verbose=False)
print(r["entry_counts"])
```

- **Example output:**

```text
{'result_db': 2}
```

## 12) `dnakit.alignment.mmseqs_besthitperset()`

- **Purpose:** For every sequence set keep the best hit of the sequence-level result and update its p-value (`mmseqs besthitperset`).
- **Method:** Runs the MMseqs2 module with a fixed argument list and no shell; options you do not set are not passed, so MMseqs2's defaults apply; output prefixes and files must not exist yet.
- **API:** `dnakit.alignment.mmseqs_besthitperset(query_set_db, target_set_db, result_db, best_hit_db, *, simple_best_hit=None, threads=1, mmseqs_path=None, verbose=True)`.
- **Input:** `query_set_db`: query set database from mmseqs_multihitdb.; `target_set_db`: target set database from mmseqs_multihitdb.; `result_db`: sequence-level search or alignment result of the two set databases.; `best_hit_db`: new result database with one best hit per set pair..
- **Output:** A dict with `status`, `version`, `command`, `outputs` (new files or database prefixes), `entry_counts` or `line_counts`, `log_tail` and `provenance`; some functions add `records`, `rows` or `clusters`.
- **Example code:**

```python
import dnakit.alignment.mmseqs2_p07_align_cli_functions as al
import dnakit.alignment.mmseqs2_p07_taxonomy_cli_functions as tx

al.mmseqs_prefilter("work/sets1", "work/sets2", "work/spf", verbose=False)
al.mmseqs_rescorediagonal("work/sets1", "work/sets2", "work/spf", "work/saln", rescore_mode=2, verbose=False)
r = tx.mmseqs_besthitperset("work/sets1", "work/sets2", "work/saln", "work/best", verbose=False)
print(r["entry_counts"])
```

- **Example output:**

```text
{'best_hit_db': 30}
```

## 13) `dnakit.alignment.mmseqs_combinepvalperset()`

- **Purpose:** Combine the hit p-values of every pair of sequence sets (`mmseqs combinepvalperset`).
- **Method:** Runs the MMseqs2 module with a fixed argument list and no shell; options you do not set are not passed, so MMseqs2's defaults apply; output prefixes and files must not exist yet.
- **API:** `dnakit.alignment.mmseqs_combinepvalperset(query_set_db, target_set_db, result_db, pvalue_db, *, alpha=None, aggregation_mode=None, threads=1, mmseqs_path=None, verbose=True)`.
- **Input:** `query_set_db`: query set database from mmseqs_multihitdb.; `target_set_db`: target set database from mmseqs_multihitdb.; `result_db`: set-level result, for example the output of mmseqs_multihitsearch.; `pvalue_db`: new p-value database..
- **Output:** A dict with `status`, `version`, `command`, `outputs` (new files or database prefixes), `entry_counts` or `line_counts`, `log_tail` and `provenance`; some functions add `records`, `rows` or `clusters`.
- **Example code:**

```python
import dnakit.alignment.mmseqs2_p07_taxonomy_cli_functions as tx

r = tx.mmseqs_combinepvalperset("work/sets1", "work/sets2", "work/hits", "work/pvalues", verbose=False)
print(r["entry_counts"])
```

- **Example output:**

```text
{'pvalue_db': 2}
```

## 14) `dnakit.alignment.mmseqs_mergeresultsbyset()`

- **Purpose:** Merge the results of the ORFs of every set (contig) back onto the set (`mmseqs mergeresultsbyset`).
- **Method:** Runs the MMseqs2 module with a fixed argument list and no shell; options you do not set are not passed, so MMseqs2's defaults apply; output prefixes and files must not exist yet.
- **API:** `dnakit.alignment.mmseqs_mergeresultsbyset(set_to_member_db, result_db, merged_db, *, , threads=1, mmseqs_path=None, verbose=True)`.
- **Input:** `set_to_member_db`: the <set_db>_set_to_member database that mmseqs_multihitdb wrote.; `result_db`: sequence-level result to merge, for example from mmseqs_besthitperset.; `merged_db`: new result database with one entry per set..
- **Output:** A dict with `status`, `version`, `command`, `outputs` (new files or database prefixes), `entry_counts` or `line_counts`, `log_tail` and `provenance`; some functions add `records`, `rows` or `clusters`.
- **Example code:**

```python
import dnakit.alignment.mmseqs2_p07_taxonomy_cli_functions as tx

r = tx.mmseqs_mergeresultsbyset("work/sets1_set_to_member", "work/best", "work/merged", verbose=False)
print(r["entry_counts"])
```

- **Example output:**

```text
{'merged_db': 2}
```
