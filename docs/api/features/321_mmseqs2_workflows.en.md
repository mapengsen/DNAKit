# MMseqs2 database, search, clustering and export workflows

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

**Example file:** `features.gff3`

```text
##gff-version 3
contig0	toy	CDS	41	223	.	+	0	ID=cds0
contig0	toy	CDS	264	446	.	-	0	ID=cds1
contig1	toy	CDS	41	223	.	+	0	ID=cds2
contig1	toy	CDS	264	446	.	-	0	ID=cds3
contig2	toy	CDS	41	223	.	+	0	ID=cds4
contig2	toy	CDS	264	446	.	-	0	ID=cds5
contig3	toy	CDS	41	223	.	+	0	ID=cds6
contig3	toy	CDS	264	446	.	-	0	ID=cds7
```

## 1) `dnakit.alignment.mmseqs_createdb()`

- **Purpose:** Convert FASTA/FASTQ files into an MMseqs2 sequence database (`mmseqs createdb`).
- **Method:** Runs the MMseqs2 module with a fixed argument list and no shell; options you do not set are not passed, so MMseqs2's defaults apply; output prefixes and files must not exist yet.
- **API:** `dnakit.alignment.mmseqs_createdb(fasta_files, database, *, dbtype=None, shuffle=None, createdb_mode=None, write_lookup=None, ..., threads=1, mmseqs_path=None, verbose=True)`.
- **Input:** `fasta_files`: FASTA or FASTQ file(s), plain, gzip or bzip2.; `database`: prefix of the new sequence database..
- **Output:** A dict with `status`, `version`, `command`, `outputs` (new files or database prefixes), `entry_counts` or `line_counts`, `log_tail` and `provenance`; some functions add `records`, `rows` or `clusters`.
- **Example code:**

```python
import shutil
from pathlib import Path

import dnakit.alignment.mmseqs2_p07_align_cli_functions as al
import dnakit.alignment.mmseqs2_p07_workflow_cli_functions as wf

shutil.rmtree("work", ignore_errors=True)
Path("work").mkdir()
result = wf.mmseqs_createdb("proteins.fa", "work/prot", verbose=False)
wf.mmseqs_createdb("contigs.fa", "work/nt", dbtype=2, verbose=False)
al.mmseqs_prefilter("work/prot", "work/prot", "work/pf", verbose=False)
al.mmseqs_rescorediagonal("work/prot", "work/prot", "work/pf", "work/aln", rescore_mode=2, add_backtrace=True, verbose=False)
print(result["status"], result["entry_counts"])
```

- **Example output:**

```text
completed {'database': 12}
```

## 2) `dnakit.alignment.mmseqs_map()`

- **Purpose:** Map nearly identical sequences of a query database onto a target database (`mmseqs map`).
- **Method:** Runs the MMseqs2 module with a fixed argument list and no shell; options you do not set are not passed, so MMseqs2's defaults apply; output prefixes and files must not exist yet.
- **API:** `dnakit.alignment.mmseqs_map(query_db, target_db, alignment_db, *, e_value=None, min_seq_id=None, coverage=None, cov_mode=None, ..., threads=1, mmseqs_path=None, verbose=True)`.
- **Input:** `query_db`: query sequence database.; `target_db`: target sequence database.; `alignment_db`: new alignment database..
- **Output:** A dict with `status`, `version`, `command`, `outputs` (new files or database prefixes), `entry_counts` or `line_counts`, `log_tail` and `provenance`; some functions add `records`, `rows` or `clusters`.
- **Example code:**

```python
import dnakit.alignment.mmseqs2_p07_workflow_cli_functions as wf

r = wf.mmseqs_map("work/prot", "work/prot", "work/map", verbose=False)
print(r["entry_counts"])
```

- **Example output:**

```text
{'alignment_db': 12}
```

## 3) `dnakit.alignment.mmseqs_convertalis()`

- **Purpose:** Export an alignment database to a BLAST-tab, SAM or custom-column table (`mmseqs convertalis`).
- **Method:** Runs the MMseqs2 module with a fixed argument list and no shell; options you do not set are not passed, so MMseqs2's defaults apply; output prefixes and files must not exist yet.
- **API:** `dnakit.alignment.mmseqs_convertalis(query_db, target_db, alignment_db, output, *, format_mode=None, format_output=None, search_type=None, translation_table=None, ..., threads=1, mmseqs_path=None, verbose=True)`.
- **Input:** `query_db`: query sequence database.; `target_db`: target sequence database.; `alignment_db`: alignment database.; `output`: new output table..
- **Output:** A dict with `status`, `version`, `command`, `outputs` (new files or database prefixes), `entry_counts` or `line_counts`, `log_tail` and `provenance`; some functions add `records`, `rows` or `clusters`.
- **Example code:**

```python
import dnakit.alignment.mmseqs2_p07_workflow_cli_functions as wf

r = wf.mmseqs_convertalis("work/prot", "work/prot", "work/aln", "work/aln.tsv", format_output="query,target,fident,bits", verbose=False)
print(r["line_counts"])
print(r["records"][1])
```

- **Example output:**

```text
{'output': 36}
{'query': 'fam0_m0', 'target': 'fam0_m2', 'fident': '0.933', 'bits': '173'}
```

## 4) `dnakit.alignment.mmseqs_clust()`

- **Purpose:** Cluster an alignment or prefilter result graph with set-cover, connected-component or greedy-by-length clustering (`mmseqs clust`).
- **Method:** Runs the MMseqs2 module with a fixed argument list and no shell; options you do not set are not passed, so MMseqs2's defaults apply; output prefixes and files must not exist yet.
- **API:** `dnakit.alignment.mmseqs_clust(sequence_db, result_db, cluster_db, *, cluster_mode=None, max_iterations=None, similarity_type=None, set_mode=None, threads=1, mmseqs_path=None, verbose=True)`.
- **Input:** `sequence_db`: sequence database.; `result_db`: alignment or prefilter result database.; `cluster_db`: new cluster database..
- **Output:** A dict with `status`, `version`, `command`, `outputs` (new files or database prefixes), `entry_counts` or `line_counts`, `log_tail` and `provenance`; some functions add `records`, `rows` or `clusters`.
- **Example code:**

```python
import dnakit.alignment.mmseqs2_p07_workflow_cli_functions as wf

r = wf.mmseqs_clust("work/prot", "work/aln", "work/clu", verbose=False)
print(r["entry_counts"])
```

- **Example output:**

```text
{'cluster_db': 4}
```

## 5) `dnakit.alignment.mmseqs_createtsv()`

- **Purpose:** Export a result database (clusters, alignments, taxonomy) to a tab-separated file (`mmseqs createtsv`).
- **Method:** Runs the MMseqs2 module with a fixed argument list and no shell; options you do not set are not passed, so MMseqs2's defaults apply; output prefixes and files must not exist yet.
- **API:** `dnakit.alignment.mmseqs_createtsv(query_db, result_db, output, *, target_db=None, first_seq_as_repr=None, target_column=None, full_header=None, ..., threads=1, mmseqs_path=None, verbose=True)`.
- **Input:** `query_db`: query sequence database.; `result_db`: result database.; `output`: new TSV file..
- **Output:** A dict with `status`, `version`, `command`, `outputs` (new files or database prefixes), `entry_counts` or `line_counts`, `log_tail` and `provenance`; some functions add `records`, `rows` or `clusters`.
- **Example code:**

```python
import dnakit.alignment.mmseqs2_p07_workflow_cli_functions as wf

r = wf.mmseqs_createtsv("work/prot", "work/clu", "work/clu.tsv", target_db="work/prot", verbose=False)
print(r["rows"][:3])
```

- **Example output:**

```text
[['fam0_m2', 'fam0_m2'], ['fam0_m2', 'fam0_m0'], ['fam0_m2', 'fam0_m1']]
```

## 6) `dnakit.alignment.mmseqs_convert2fasta()`

- **Purpose:** Write a sequence database as a FASTA file (`mmseqs convert2fasta`).
- **Method:** Runs the MMseqs2 module with a fixed argument list and no shell; options you do not set are not passed, so MMseqs2's defaults apply; output prefixes and files must not exist yet.
- **API:** `dnakit.alignment.mmseqs_convert2fasta(sequence_db, output, *, use_header_file=None, threads=1, mmseqs_path=None, verbose=True)`.
- **Input:** `sequence_db`: sequence database.; `output`: new FASTA file..
- **Output:** A dict with `status`, `version`, `command`, `outputs` (new files or database prefixes), `entry_counts` or `line_counts`, `log_tail` and `provenance`; some functions add `records`, `rows` or `clusters`.
- **Example code:**

```python
import dnakit.alignment.mmseqs2_p07_workflow_cli_functions as wf

r = wf.mmseqs_convert2fasta("work/prot", "work/prot.fa", verbose=False)
print(r["sequence_count"])
```

- **Example output:**

```text
12
```

## 7) `dnakit.alignment.mmseqs_result2repseq()`

- **Purpose:** Extract the representative sequence of every cluster of a result database (`mmseqs result2repseq`).
- **Method:** Runs the MMseqs2 module with a fixed argument list and no shell; options you do not set are not passed, so MMseqs2's defaults apply; output prefixes and files must not exist yet.
- **API:** `dnakit.alignment.mmseqs_result2repseq(sequence_db, result_db, representative_db, *, , threads=1, mmseqs_path=None, verbose=True)`.
- **Input:** `sequence_db`: sequence database.; `result_db`: cluster or result database.; `representative_db`: new sequence database of the representatives..
- **Output:** A dict with `status`, `version`, `command`, `outputs` (new files or database prefixes), `entry_counts` or `line_counts`, `log_tail` and `provenance`; some functions add `records`, `rows` or `clusters`.
- **Example code:**

```python
import dnakit.alignment.mmseqs2_p07_workflow_cli_functions as wf

r = wf.mmseqs_result2repseq("work/prot", "work/clu", "work/rep", verbose=False)
print(r["entry_counts"])
```

- **Example output:**

```text
{'representative_db': 4}
```

## 8) `dnakit.alignment.mmseqs_pickconsensusrep()`

- **Purpose:** Choose a new representative sequence per cluster based on the cluster consensus (`mmseqs pickconsensusrep`).
- **Method:** Runs the MMseqs2 module with a fixed argument list and no shell; options you do not set are not passed, so MMseqs2's defaults apply; output prefixes and files must not exist yet.
- **API:** `dnakit.alignment.mmseqs_pickconsensusrep(sequence_db, cluster_db, new_cluster_db, *, , threads=1, mmseqs_path=None, verbose=True)`.
- **Input:** `sequence_db`: sequence database.; `cluster_db`: cluster database.; `new_cluster_db`: new cluster database..
- **Output:** A dict with `status`, `version`, `command`, `outputs` (new files or database prefixes), `entry_counts` or `line_counts`, `log_tail` and `provenance`; some functions add `records`, `rows` or `clusters`.
- **Example code:**

```python
import dnakit.alignment.mmseqs2_p07_workflow_cli_functions as wf

r = wf.mmseqs_pickconsensusrep("work/prot", "work/clu", "work/clu_consensus", verbose=False)
print(r["entry_counts"])
```

- **Example output:**

```text
{'new_cluster_db': 4}
```

## 9) `dnakit.alignment.mmseqs_clusthash()`

- **Purpose:** Cluster equal-length sequences by hashing them (`mmseqs clusthash`); the result is an alignment database of identical or near-identical pairs.
- **Method:** Runs the MMseqs2 module with a fixed argument list and no shell; options you do not set are not passed, so MMseqs2's defaults apply; output prefixes and files must not exist yet.
- **API:** `dnakit.alignment.mmseqs_clusthash(sequence_db, alignment_db, *, min_seq_id=None, max_seq_len=None, threads=1, mmseqs_path=None, verbose=True)`.
- **Input:** `sequence_db`: sequence database.; `alignment_db`: new alignment database..
- **Output:** A dict with `status`, `version`, `command`, `outputs` (new files or database prefixes), `entry_counts` or `line_counts`, `log_tail` and `provenance`; some functions add `records`, `rows` or `clusters`.
- **Example code:**

```python
import dnakit.alignment.mmseqs2_p07_workflow_cli_functions as wf

r = wf.mmseqs_clusthash("work/prot", "work/hash", verbose=False)
print(r["entry_counts"])
```

- **Example output:**

```text
{'alignment_db': 12}
```

## 10) `dnakit.alignment.mmseqs_clusterupdate()`

- **Purpose:** Update an earlier clustering with new sequences (`mmseqs clusterupdate`); MMseqs2 deletes the old cluster database.
- **Method:** Runs the MMseqs2 module with a fixed argument list and no shell; options you do not set are not passed, so MMseqs2's defaults apply; output prefixes and files must not exist yet.
- **API:** `dnakit.alignment.mmseqs_clusterupdate(old_sequence_db, new_sequence_db, old_cluster_db, mapped_sequence_db, new_cluster_db, *, e_value=None, min_seq_id=None, coverage=None, cov_mode=None, ..., threads=1, mmseqs_path=None, verbose=True)`.
- **Input:** `old_sequence_db`: sequence database that was clustered before.; `new_sequence_db`: sequence database with the old and the new sequences.; `old_cluster_db`: previous cluster database.; `mapped_sequence_db`: new sequence database whose keys are mapped to the old ones.; `new_cluster_db`: new cluster database..
- **Output:** A dict with `status`, `version`, `command`, `outputs` (new files or database prefixes), `entry_counts` or `line_counts`, `log_tail` and `provenance`; some functions add `records`, `rows` or `clusters`.
- **Example code:**

```python
import dnakit.alignment.mmseqs2_p07_workflow_cli_functions as wf

wf.mmseqs_clust("work/prot", "work/aln", "work/old_clu", verbose=False)
r = wf.mmseqs_clusterupdate("work/prot", "work/prot", "work/old_clu", "work/mapped", "work/new_clu", verbose=False)
print(r["entry_counts"], r["old_cluster_db_removed"])
```

- **Example output:**

```text
{'mapped_sequence_db': 12, 'new_cluster_db': 4} True
```

## 11) `dnakit.alignment.mmseqs_summarizealis()`

- **Purpose:** Summarise the alignments of every query into one row with unique coverage, coverage and mean sequence identity (`mmseqs summarizealis`).
- **Method:** Runs the MMseqs2 module with a fixed argument list and no shell; options you do not set are not passed, so MMseqs2's defaults apply; output prefixes and files must not exist yet.
- **API:** `dnakit.alignment.mmseqs_summarizealis(alignment_db, summary_db, *, , threads=1, mmseqs_path=None, verbose=True)`.
- **Input:** `alignment_db`: alignment database.; `summary_db`: new summary database..
- **Output:** A dict with `status`, `version`, `command`, `outputs` (new files or database prefixes), `entry_counts` or `line_counts`, `log_tail` and `provenance`; some functions add `records`, `rows` or `clusters`.
- **Example code:**

```python
import dnakit.alignment.mmseqs2_p07_workflow_cli_functions as wf

r = wf.mmseqs_summarizealis("work/aln", "work/summary", verbose=False)
print(r["entry_counts"])
```

- **Example output:**

```text
{'summary_db': 12}
```

## 12) `dnakit.alignment.mmseqs_filterresult()`

- **Purpose:** Filter the pairwise alignment results of every query by redundancy, identity to the query, score and coverage (`mmseqs filterresult`).
- **Method:** Runs the MMseqs2 module with a fixed argument list and no shell; options you do not set are not passed, so MMseqs2's defaults apply; output prefixes and files must not exist yet.
- **API:** `dnakit.alignment.mmseqs_filterresult(query_db, target_db, result_db, filtered_db, *, max_seq_id=None, qid=None, qsc=None, cov=None, ..., threads=1, mmseqs_path=None, verbose=True)`.
- **Input:** `query_db`: query sequence database.; `target_db`: target sequence database.; `result_db`: alignment result database.; `filtered_db`: new filtered result database..
- **Output:** A dict with `status`, `version`, `command`, `outputs` (new files or database prefixes), `entry_counts` or `line_counts`, `log_tail` and `provenance`; some functions add `records`, `rows` or `clusters`.
- **Example code:**

```python
import dnakit.alignment.mmseqs2_p07_workflow_cli_functions as wf

r = wf.mmseqs_filterresult("work/prot", "work/prot", "work/aln", "work/filtered", qid="0.9", verbose=False)
print(r["entry_counts"])
```

- **Example output:**

```text
{'filtered_db': 12}
```

## 13) `dnakit.alignment.mmseqs_extractalignedregion()`

- **Purpose:** Cut the aligned region out of the query or the target sequences (`mmseqs extractalignedregion`).
- **Method:** Runs the MMseqs2 module with a fixed argument list and no shell; options you do not set are not passed, so MMseqs2's defaults apply; output prefixes and files must not exist yet.
- **API:** `dnakit.alignment.mmseqs_extractalignedregion(query_db, target_db, alignment_db, region_db, *, extract_mode=None, threads=1, mmseqs_path=None, verbose=True)`.
- **Input:** `query_db`: query sequence database.; `target_db`: target sequence database.; `alignment_db`: alignment database.; `region_db`: new sequence database of the aligned regions..
- **Output:** A dict with `status`, `version`, `command`, `outputs` (new files or database prefixes), `entry_counts` or `line_counts`, `log_tail` and `provenance`; some functions add `records`, `rows` or `clusters`.
- **Example code:**

```python
import dnakit.alignment.mmseqs2_p07_workflow_cli_functions as wf

r = wf.mmseqs_extractalignedregion("work/prot", "work/prot", "work/aln", "work/region", extract_mode=1, verbose=False)
print(r["entry_counts"])
```

- **Example output:**

```text
{'region_db': 36}
```

## 14) `dnakit.alignment.mmseqs_result2msa()`

- **Purpose:** Multiple sequence alignment per query (``mmseqs result2msa``, ``mmseqs result2dnamsa``).
- **Method:** Runs the MMseqs2 module with a fixed argument list and no shell; options you do not set are not passed, so MMseqs2's defaults apply; output prefixes and files must not exist yet.
- **API:** `dnakit.alignment.mmseqs_result2msa(query_db, target_db, result_db, msa_db, *, dna=False, msa_format_mode=None, filter_msa=None, filter_min_enable=None, ..., threads=1, mmseqs_path=None, verbose=True)`.
- **Input:** `query_db`: see the docstring; `target_db`: see the docstring; `result_db`: see the docstring; `msa_db`: see the docstring.
- **Output:** A dict with `status`, `version`, `command`, `outputs` (new files or database prefixes), `entry_counts` or `line_counts`, `log_tail` and `provenance`; some functions add `records`, `rows` or `clusters`.
- **Example code:**

```python
import dnakit.alignment.mmseqs2_p07_align_cli_functions as al
import dnakit.alignment.mmseqs2_p07_workflow_cli_functions as wf

r = wf.mmseqs_result2msa("work/prot", "work/prot", "work/aln", "work/msa", verbose=False)
al.mmseqs_prefilter("work/nt", "work/nt", "work/ntpf", verbose=False)
al.mmseqs_rescorediagonal("work/nt", "work/nt", "work/ntpf", "work/ntaln", rescore_mode=2, add_backtrace=True, verbose=False)
d = wf.mmseqs_result2msa("work/nt", "work/nt", "work/ntaln", "work/dnamsa", dna=True, verbose=False)
print(r["entry_counts"], d["entry_counts"])
```

- **Example output:**

```text
{'msa_db': 12} {'msa_db': 4}
```

## 15) `dnakit.alignment.mmseqs_gff2db()`

- **Purpose:** Extract the regions of GFF3 features from a sequence database into a new database (`mmseqs gff2db`).
- **Method:** Runs the MMseqs2 module with a fixed argument list and no shell; options you do not set are not passed, so MMseqs2's defaults apply; output prefixes and files must not exist yet.
- **API:** `dnakit.alignment.mmseqs_gff2db(gff_files, sequence_db, feature_db, *, gff_type=None, id_offset=None, threads=1, mmseqs_path=None, verbose=True)`.
- **Input:** `gff_files`: GFF3 file(s).; `sequence_db`: nucleotide sequence database of the genome.; `feature_db`: new sequence database of the feature sequences..
- **Output:** A dict with `status`, `version`, `command`, `outputs` (new files or database prefixes), `entry_counts` or `line_counts`, `log_tail` and `provenance`; some functions add `records`, `rows` or `clusters`.
- **Example code:**

```python
import dnakit.alignment.mmseqs2_p07_workflow_cli_functions as wf

r = wf.mmseqs_gff2db(["features.gff3"], "work/nt", "work/features", gff_type="CDS", verbose=False)
print(r["entry_counts"])
```

- **Example output:**

```text
{'feature_db': 8}
```

## 16) `dnakit.alignment.mmseqs_maskbygff()`

- **Purpose:** Mask the regions of selected GFF3 features in a sequence database (``mmseqs maskbygff``).
- **Method:** Runs the MMseqs2 module with a fixed argument list and no shell; options you do not set are not passed, so MMseqs2's defaults apply; output prefixes and files must not exist yet.
- **API:** `dnakit.alignment.mmseqs_maskbygff(gff_file, sequence_db, masked_db, *, gff_type=None, id_offset=None, names_as_ids=True, threads=1, mmseqs_path=None, verbose=True)`.
- **Input:** `gff_file`: GFF3 file. sequence_db: sequence database of the genome. masked_db: new database.; `sequence_db`: see the docstring; `masked_db`: see the docstring.
- **Output:** A dict with `status`, `version`, `command`, `outputs` (new files or database prefixes), `entry_counts` or `line_counts`, `log_tail` and `provenance`; some functions add `records`, `rows` or `clusters`.
- **Example code:**

```python
import dnakit.alignment.mmseqs2_p07_workflow_cli_functions as wf

r = wf.mmseqs_maskbygff("features.gff3", "work/nt", "work/masked", verbose=False)
print(r["entry_counts"])
```

- **Example output:**

```text
{'masked_db': 4}
```

## 17) `dnakit.alignment.mmseqs_createlinindex()`

- **Purpose:** Create the bottom-m-hash k-mer index of a sequence database (``mmseqs createlinindex``).
- **Method:** Runs the MMseqs2 module with a fixed argument list and no shell; options you do not set are not passed, so MMseqs2's defaults apply; output prefixes and files must not exist yet.
- **API:** `dnakit.alignment.mmseqs_createlinindex(sequence_db, *, kmer_size=None, kmer_per_seq=None, min_seq_id=None, search_type=None, ..., threads=1, mmseqs_path=None, verbose=True)`.
- **Input:** `sequence_db`: prefix of the sequence database to index..
- **Output:** A dict with `status`, `version`, `command`, `outputs` (new files or database prefixes), `entry_counts` or `line_counts`, `log_tail` and `provenance`; some functions add `records`, `rows` or `clusters`.
- **Example code:**

```python
import dnakit.alignment.mmseqs2_p07_workflow_cli_functions as wf

r = wf.mmseqs_createlinindex("work/nt", search_type=3, verbose=False)
print(r["entry_counts"])
```

- **Example output:**

```text
{'kmer_index_db': 18}
```

## 18) `dnakit.alignment.mmseqs_linsearch()`

- **Purpose:** Linear-time, less sensitive sequence search (``mmseqs easy-linsearch``, ``mmseqs linsearch``).
- **Method:** Runs the MMseqs2 module with a fixed argument list and no shell; options you do not set are not passed, so MMseqs2's defaults apply; output prefixes and files must not exist yet.
- **API:** `dnakit.alignment.mmseqs_linsearch(query, target, output, *, input_kind='fasta', search_type=None, e_value=None, min_seq_id=None, ..., threads=1, mmseqs_path=None, verbose=True)`.
- **Input:** `query`: see the docstring; `target`: see the docstring; `output`: see the docstring.
- **Output:** A dict with `status`, `version`, `command`, `outputs` (new files or database prefixes), `entry_counts` or `line_counts`, `log_tail` and `provenance`; some functions add `records`, `rows` or `clusters`.
- **Example code:**

```python
import dnakit.alignment.mmseqs2_p07_workflow_cli_functions as wf

r = wf.mmseqs_linsearch("contigs.fa", "contigs.fa", "work/lin.tsv", search_type=3, verbose=False)
print(sorted({row["query"] for row in r["records"]}))
d = wf.mmseqs_linsearch("work/nt", "work/nt", "work/lin_db", input_kind="database", search_type=3, verbose=False)
print(d["entry_counts"])
```

- **Example output:**

```text
['contig0', 'contig1', 'contig2', 'contig3']
{'alignment_db': 4}
```

## 19) `dnakit.alignment.mmseqs_linclust()`

- **Purpose:** Linear-time, less sensitive clustering (``mmseqs easy-linclust``, ``mmseqs linclust``).
- **Method:** Runs the MMseqs2 module with a fixed argument list and no shell; options you do not set are not passed, so MMseqs2's defaults apply; output prefixes and files must not exist yet.
- **API:** `dnakit.alignment.mmseqs_linclust(sequences, output, *, input_kind='fasta', min_seq_id=None, coverage=None, cov_mode=None, ..., threads=1, mmseqs_path=None, verbose=True)`.
- **Input:** `sequences`: see the docstring; `output`: see the docstring.
- **Output:** A dict with `status`, `version`, `command`, `outputs` (new files or database prefixes), `entry_counts` or `line_counts`, `log_tail` and `provenance`; some functions add `records`, `rows` or `clusters`.
- **Example code:**

```python
import dnakit.alignment.mmseqs2_p07_workflow_cli_functions as wf

r = wf.mmseqs_linclust("proteins.fa", "work/lc", verbose=False)
print(r["cluster_count"], r["sequence_count"])
```

- **Example output:**

```text
7 12
```

## 20) `dnakit.alignment.mmseqs_rbh()`

- **Purpose:** Reciprocal best-hit search of two sequence sets (``mmseqs easy-rbh``, ``mmseqs rbh``).
- **Method:** Runs the MMseqs2 module with a fixed argument list and no shell; options you do not set are not passed, so MMseqs2's defaults apply; output prefixes and files must not exist yet.
- **API:** `dnakit.alignment.mmseqs_rbh(query, target, output, *, input_kind='fasta', sensitivity=None, e_value=None, min_seq_id=None, ..., threads=1, mmseqs_path=None, verbose=True)`.
- **Input:** `query`: see the docstring; `target`: see the docstring; `output`: see the docstring.
- **Output:** A dict with `status`, `version`, `command`, `outputs` (new files or database prefixes), `entry_counts` or `line_counts`, `log_tail` and `provenance`; some functions add `records`, `rows` or `clusters`.
- **Example code:**

```python
import dnakit.alignment.mmseqs2_p07_workflow_cli_functions as wf

r = wf.mmseqs_rbh("proteins.fa", "proteins.fa", "work/rbh.tsv", verbose=False)
print(r["pair_count"], r["records"][0]["query"])
```

- **Example output:**

```text
12 fam0_m0
```
