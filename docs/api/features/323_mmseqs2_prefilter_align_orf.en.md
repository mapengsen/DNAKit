# MMseqs2 prefilters, alignment modules and ORF utilities

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

**Example file:** `cds.fa`

```text
>fam0_m0
ATGATGTTTCCTTGTGATGTTGAAAATTGGTGTACTCATTGTGATCAACAAGATATTGATGTTCAATGTTGGGAAATTTGGTGTTGGTGGCCTTGTATTTGTGTTTTTCTGCAATTTGTTGAATGGCTGGTTGGTGAATGGTGGCATAATGAAGTTGATTGGTGTTATCATTCTGTTCAA
>fam0_m1
ATGATGTTTCCTTGTGATGTTCTGAATTGGTGTACTCATTGTGATCAACAAGATGATGATGTTCAATGTTGGGAAATTTGGTGTTGGTGGCCTCTGATTTGTGTTTTTCTGCAACGTGTTGAATGGCTGGTTGGTGAATGGTGGCATAATGAAGTTGATTGGTGTTATCATTCTGTTTTT
>fam0_m2
ATGATGTTTCCTTGTCAAGTTGAAAATTGGTGTACTCATTGTGATCAACAAGATATTGATGTTCAATGTTGGGAAATTTGGTGTTGGTGGCCTTGTATTCATGTTTTTCTGCAATTTGTTGAATGGCTGGTTGAAGAATGGTGGCATAATGAAGTTGATTGGTGTTATCATTCTAAACAA
>fam1_m0
ATGGGTAATATTGTTGTTACTATGATTTATCATATTCCTATTCATACTTCTAATGCTGCTAAATCTAAACATTATAATCGTAATAATGATATTGAAATTTCTCATATGCATTCTTATTATGCTTCTAATGATGAACCTCATTCTGGTCAAATGGATCCTCGTCCTGATGGTGGTTTTGCT
>fam1_m1
ATGGGTAATATTGTTGTTACTATGATTTATGAAATTCCTATTCATACTTCTAATGCTGCTAAATCTAAACATTATAATCGTAATAATGATATTGAAATTTCTCGTATGCATTCTTATTATGCTTCTACTGATGAACCTCATCATGGTCAAATGGATCCTCGTCCTGATGGTGGTTTTGCT
>fam1_m2
ATGGGTAATATTGTTGTTACTATGATTTATCATATTCCTATTCATACTTCTAATGCTGCTAAATCTAAACATTATAATTCTAATAATGATATTGAAATTTCTCATAATGTTTCTTATTATGCTTCTAATGATGAACCTCATAAAGGTCAAATGGATCCTCGTCCTGATGGTGGTTTTTGT
>fam2_m0
ATGGGTCAAGATAAAGCTGATAAAGATTATATTGATAAAGAACGTGCTATGGTTCAAAAATATTTTTGTACTATTGAAGGTAAATGTGGTCATCTGCTGACTCATCTGCGTACTGGTAAAAATGCTAAATGTGCTGCTACTGTTCATACTTCTATTCGTGAACAATCTGTTCCTACTCTG
>fam2_m1
ATGGGTCAAGATAAAGCTGATAAATTTAAAATTCCTAAAGAACGTGCTATGGTTCAAAAATATTTTTGTACTATTGAAGGTAAACCTAAACATCTGCTGAAACATCTGCGTACTGGTAAAACTGCTAAATGTGCTGCTACTGTTCATACTTCTATTCGTGAACAATCTGTTCCTACTCTG
>fam2_m2
ATGGGTTGGGATAAAGCTGATTGTGATTATATTGATAAAGAACGTGCTATGGATCAAAAATATTTTTCTACTAAAGAAGGTAAATGTGGTCATCTGCTGACTCATTTTCGTACTGGTAAAAATGCTAAATGTCATGCTACTGTTCATACTTCTATTCGTGAACAACGTACTCCTACTCTG
>fam3_m0
ATGCATCATGATTGGGATTTTACTAAAAATTTTTATACTAAAGAAAATATTTCTTCTCCTGCTGGTGCTTCTCGTCCTCTGTTTCAAAATCCTATGGAAATGGCTATGATGCCTGAACATGCTCTGAAAAATGATCCTCCTTGGGATAATCAAAAATGTAAAGAATGTCTGTTTATTAAA
>fam3_m1
ATGCATCATGATTGGCCTTTTACTAAAAATCAATATACTAAAGAAGTTATTTCTTCTCCTGCTGGTGCTTCTCGTCCTCTGTTTCAAAATCCTATGGAAATGGCTATGATGCCTATGCATGCTCTGAAAAATCAACCTCCTTGGGATAATCAAAAATGTAAAGAATGTCTGATTATTAAA
>fam3_m2
ATGCATCATGATTGGGATTTTACTAAAAATACTTATACTAAAGAAAATATTTCTTCTCCTGCTGGTGCTTCTTTTCCTCTGTTTCAAAATCCTATGGAAATGACTATGATGCCTGAAGTTGCTCTGAAAAATGATCCTCCTTGGGATAATCAATCTTGTAAAATTTGTCTGTTTATGAAA
```

## 1) `dnakit.alignment.mmseqs_prefilter()`

- **Purpose:** k-mer prefilter with double consecutive diagonal hits (`mmseqs prefilter`).
- **Method:** Runs the MMseqs2 module with a fixed argument list and no shell; options you do not set are not passed, so MMseqs2's defaults apply; output prefixes and files must not exist yet.
- **API:** `dnakit.alignment.mmseqs_prefilter(query_db, target_db, prefilter_db, *, sensitivity=None, kmer_size=None, max_seqs=None, min_ungapped_score=None, ..., threads=1, mmseqs_path=None, verbose=True)`.
- **Input:** `query_db`: query sequence database.; `target_db`: target sequence database.; `prefilter_db`: new prefilter result database..
- **Output:** A dict with `status`, `version`, `command`, `outputs` (new files or database prefixes), `entry_counts` or `line_counts`, `log_tail` and `provenance`; some functions add `records`, `rows` or `clusters`.
- **Example code:**

```python
import shutil
from pathlib import Path

import dnakit.alignment.mmseqs2_p07_align_cli_functions as al
import dnakit.alignment.mmseqs2_p07_workflow_cli_functions as wf

shutil.rmtree("work", ignore_errors=True)
Path("work").mkdir()
wf.mmseqs_createdb("proteins.fa", "work/prot", verbose=False)
wf.mmseqs_createdb("contigs.fa", "work/nt", dbtype=2, verbose=False)
wf.mmseqs_createdb("proteins.fa", "work/prot0", shuffle=False, verbose=False)
wf.mmseqs_createdb("cds.fa", "work/cds", shuffle=False, verbose=False)
r = al.mmseqs_prefilter("work/prot", "work/prot", "work/pf", verbose=False)
print(r["entry_counts"])
```

- **Example output:**

```text
{'prefilter_db': 12}
```

## 2) `dnakit.alignment.mmseqs_ungappedprefilter()`

- **Purpose:** Prefilter that scores the optimal ungapped diagonal of every query-target pair (`mmseqs ungappedprefilter`).
- **Method:** Runs the MMseqs2 module with a fixed argument list and no shell; options you do not set are not passed, so MMseqs2's defaults apply; output prefixes and files must not exist yet.
- **API:** `dnakit.alignment.mmseqs_ungappedprefilter(query_db, target_db, prefilter_db, *, max_seqs=None, min_ungapped_score=None, comp_bias_corr=None, coverage=None, ..., threads=1, mmseqs_path=None, verbose=True)`.
- **Input:** `query_db`: query sequence database.; `target_db`: target sequence database.; `prefilter_db`: new prefilter result database..
- **Output:** A dict with `status`, `version`, `command`, `outputs` (new files or database prefixes), `entry_counts` or `line_counts`, `log_tail` and `provenance`; some functions add `records`, `rows` or `clusters`.
- **Example code:**

```python
import dnakit.alignment.mmseqs2_p07_align_cli_functions as al

r = al.mmseqs_ungappedprefilter("work/prot", "work/prot", "work/upf", verbose=False)
print(r["entry_counts"])
```

- **Example output:**

```text
{'prefilter_db': 12}
```

## 3) `dnakit.alignment.mmseqs_gappedprefilter()`

- **Purpose:** Prefilter that scores the optimal gapped Smith-Waterman alignment of every query-target pair, slow (`mmseqs gappedprefilter`).
- **Method:** Runs the MMseqs2 module with a fixed argument list and no shell; options you do not set are not passed, so MMseqs2's defaults apply; output prefixes and files must not exist yet.
- **API:** `dnakit.alignment.mmseqs_gappedprefilter(query_db, target_db, prefilter_db, *, max_seqs=None, min_ungapped_score=None, comp_bias_corr=None, coverage=None, ..., threads=1, mmseqs_path=None, verbose=True)`.
- **Input:** `query_db`: query sequence database.; `target_db`: target sequence database.; `prefilter_db`: new prefilter result database..
- **Output:** A dict with `status`, `version`, `command`, `outputs` (new files or database prefixes), `entry_counts` or `line_counts`, `log_tail` and `provenance`; some functions add `records`, `rows` or `clusters`.
- **Example code:**

```python
import dnakit.alignment.mmseqs2_p07_align_cli_functions as al

r = al.mmseqs_gappedprefilter("work/prot", "work/prot", "work/gpf", verbose=False)
print(r["entry_counts"])
```

- **Example output:**

```text
{'prefilter_db': 12}
```

## 4) `dnakit.alignment.mmseqs_kmermatcher()`

- **Purpose:** Find bottom-m-hashed k-mer matches within one sequence database (`mmseqs kmermatcher`); the building block of linear-time clustering.
- **Method:** Runs the MMseqs2 module with a fixed argument list and no shell; options you do not set are not passed, so MMseqs2's defaults apply; output prefixes and files must not exist yet.
- **API:** `dnakit.alignment.mmseqs_kmermatcher(sequence_db, prefilter_db, *, kmer_size=None, min_seq_id=None, coverage=None, cov_mode=None, ..., threads=1, mmseqs_path=None, verbose=True)`.
- **Input:** `sequence_db`: sequence database.; `prefilter_db`: new prefilter result database..
- **Output:** A dict with `status`, `version`, `command`, `outputs` (new files or database prefixes), `entry_counts` or `line_counts`, `log_tail` and `provenance`; some functions add `records`, `rows` or `clusters`.
- **Example code:**

```python
import dnakit.alignment.mmseqs2_p07_align_cli_functions as al

r = al.mmseqs_kmermatcher("work/prot", "work/kmm", verbose=False)
print(r["entry_counts"])
```

- **Example output:**

```text
{'prefilter_db': 12}
```

## 5) `dnakit.alignment.mmseqs_kmersearch()`

- **Purpose:** Find bottom-m-hashed k-mer matches between a query database and a k-mer index of a target database (`mmseqs kmersearch`).
- **Method:** Runs the MMseqs2 module with a fixed argument list and no shell; options you do not set are not passed, so MMseqs2's defaults apply; output prefixes and files must not exist yet.
- **API:** `dnakit.alignment.mmseqs_kmersearch(query_db, kmer_index_db, prefilter_db, *, coverage=None, cov_mode=None, kmer_per_seq=None, pick_n_sim_kmer=None, ..., threads=1, mmseqs_path=None, verbose=True)`.
- **Input:** `query_db`: query sequence database.; `kmer_index_db`: k-mer index of the target from mmseqs_createlinindex (the .linidx database).; `prefilter_db`: new prefilter result database..
- **Output:** A dict with `status`, `version`, `command`, `outputs` (new files or database prefixes), `entry_counts` or `line_counts`, `log_tail` and `provenance`; some functions add `records`, `rows` or `clusters`.
- **Example code:**

```python
import dnakit.alignment.mmseqs2_p07_align_cli_functions as al
import dnakit.alignment.mmseqs2_p07_workflow_cli_functions as wf

wf.mmseqs_createlinindex("work/nt", search_type=3, verbose=False)
r = al.mmseqs_kmersearch("work/nt", "work/nt.linidx", "work/kms", verbose=False)
print(r["entry_counts"])
```

- **Example output:**

```text
{'prefilter_db': 4}
```

## 6) `dnakit.alignment.mmseqs_rescorediagonal()`

- **Purpose:** Compute the sequence identity of every prefilter diagonal by Hamming distance or local or global alignment (`mmseqs rescorediagonal`).
- **Method:** Runs the MMseqs2 module with a fixed argument list and no shell; options you do not set are not passed, so MMseqs2's defaults apply; output prefixes and files must not exist yet.
- **API:** `dnakit.alignment.mmseqs_rescorediagonal(query_db, target_db, prefilter_db, result_db, *, e_value=None, min_seq_id=None, coverage=None, cov_mode=None, ..., threads=1, mmseqs_path=None, verbose=True)`.
- **Input:** `query_db`: query sequence database.; `target_db`: target sequence database.; `prefilter_db`: prefilter result database.; `result_db`: new result database..
- **Output:** A dict with `status`, `version`, `command`, `outputs` (new files or database prefixes), `entry_counts` or `line_counts`, `log_tail` and `provenance`; some functions add `records`, `rows` or `clusters`.
- **Example code:**

```python
import dnakit.alignment.mmseqs2_p07_align_cli_functions as al

r = al.mmseqs_rescorediagonal("work/prot", "work/prot", "work/pf", "work/aln", rescore_mode=2, add_backtrace=True, verbose=False)
print(r["entry_counts"])
```

- **Example output:**

```text
{'result_db': 12}
```

## 7) `dnakit.alignment.mmseqs_alignall()`

- **Purpose:** Gapped local all-versus-all alignment of the members inside every result entry (`mmseqs alignall`).
- **Method:** Runs the MMseqs2 module with a fixed argument list and no shell; options you do not set are not passed, so MMseqs2's defaults apply; output prefixes and files must not exist yet.
- **API:** `dnakit.alignment.mmseqs_alignall(sequence_db, result_db, alignment_db, *, e_value=None, min_seq_id=None, coverage=None, cov_mode=None, ..., threads=1, mmseqs_path=None, verbose=True)`.
- **Input:** `sequence_db`: sequence database.; `result_db`: result (for example cluster) database.; `alignment_db`: new alignment database..
- **Output:** A dict with `status`, `version`, `command`, `outputs` (new files or database prefixes), `entry_counts` or `line_counts`, `log_tail` and `provenance`; some functions add `records`, `rows` or `clusters`.
- **Example code:**

```python
import dnakit.alignment.mmseqs2_p07_align_cli_functions as al
import dnakit.alignment.mmseqs2_p07_workflow_cli_functions as wf

wf.mmseqs_clust("work/prot", "work/aln", "work/clu", verbose=False)
r = al.mmseqs_alignall("work/prot", "work/clu", "work/all", verbose=False)
print(r["entry_counts"])
```

- **Example output:**

```text
{'alignment_db': 4}
```

## 8) `dnakit.alignment.mmseqs_transitivealign()`

- **Purpose:** Transfer alignments through a shared sequence by transitivity (`mmseqs transitivealign`).
- **Method:** Runs the MMseqs2 module with a fixed argument list and no shell; options you do not set are not passed, so MMseqs2's defaults apply; output prefixes and files must not exist yet.
- **API:** `dnakit.alignment.mmseqs_transitivealign(sequence_db, alignment_db, transitive_db, *, e_value=None, min_seq_id=None, coverage=None, cov_mode=None, ..., threads=1, mmseqs_path=None, verbose=True)`.
- **Input:** `sequence_db`: sequence database.; `alignment_db`: alignment database to extend.; `transitive_db`: new alignment database..
- **Output:** A dict with `status`, `version`, `command`, `outputs` (new files or database prefixes), `entry_counts` or `line_counts`, `log_tail` and `provenance`; some functions add `records`, `rows` or `clusters`.
- **Example code:**

```python
import dnakit.alignment.mmseqs2_p07_align_cli_functions as al

r = al.mmseqs_transitivealign("work/prot", "work/aln", "work/trans", verbose=False)
print(r["entry_counts"])
```

- **Example output:**

```text
{'transitive_db': 12}
```

## 9) `dnakit.alignment.mmseqs_fwbw()`

- **Purpose:** Forward-backward (maximum accuracy) alignment of the pairs of an alignment database (`mmseqs fwbw`).
- **Method:** Runs the MMseqs2 module with a fixed argument list and no shell; options you do not set are not passed, so MMseqs2's defaults apply; output prefixes and files must not exist yet.
- **API:** `dnakit.alignment.mmseqs_fwbw(query_db, target_db, alignment_db, result_db, *, e_value=None, min_seq_id=None, coverage=None, cov_mode=None, ..., threads=1, mmseqs_path=None, verbose=True)`.
- **Input:** `query_db`: query sequence database.; `target_db`: target sequence database.; `alignment_db`: alignment database to realign.; `result_db`: new alignment database..
- **Output:** A dict with `status`, `version`, `command`, `outputs` (new files or database prefixes), `entry_counts` or `line_counts`, `log_tail` and `provenance`; some functions add `records`, `rows` or `clusters`.
- **Example code:**

```python
import dnakit.alignment.mmseqs2_p07_align_cli_functions as al

r = al.mmseqs_fwbw("work/prot", "work/prot", "work/aln", "work/fwbw", verbose=False)
print(r["entry_counts"])
```

- **Example output:**

```text
{'result_db': 12}
```

## 10) `dnakit.alignment.mmseqs_alignbykmer()`

- **Purpose:** Heuristic gapped local alignment seeded by k-mer matches (`mmseqs alignbykmer`).
- **Method:** Runs the MMseqs2 module with a fixed argument list and no shell; options you do not set are not passed, so MMseqs2's defaults apply; output prefixes and files must not exist yet.
- **API:** `dnakit.alignment.mmseqs_alignbykmer(query_db, target_db, result_db, alignment_db, *, kmer_size=None, spaced_kmer_mode=None, e_value=None, min_seq_id=None, ..., threads=1, mmseqs_path=None, verbose=True)`.
- **Input:** `query_db`: query sequence database.; `target_db`: target sequence database.; `result_db`: prefilter result database.; `alignment_db`: new alignment database..
- **Output:** A dict with `status`, `version`, `command`, `outputs` (new files or database prefixes), `entry_counts` or `line_counts`, `log_tail` and `provenance`; some functions add `records`, `rows` or `clusters`.
- **Example code:**

```python
import dnakit.alignment.mmseqs2_p07_align_cli_functions as al

r = al.mmseqs_alignbykmer("work/prot", "work/prot", "work/pf", "work/abk", verbose=False)
print(r["entry_counts"])
```

- **Example output:**

```text
{'alignment_db': 12}
```

## 11) `dnakit.alignment.mmseqs_extractorfs()`

- **Purpose:** Six-frame extraction of open reading frames from nucleotide sequences (`mmseqs extractorfs`).
- **Method:** Runs the MMseqs2 module with a fixed argument list and no shell; options you do not set are not passed, so MMseqs2's defaults apply; output prefixes and files must not exist yet.
- **API:** `dnakit.alignment.mmseqs_extractorfs(sequence_db, orf_db, *, min_length=None, max_length=None, max_gaps=None, contig_start_mode=None, ..., threads=1, mmseqs_path=None, verbose=True)`.
- **Input:** `sequence_db`: nucleotide sequence database.; `orf_db`: new ORF database (DNA, or protein with translate=True)..
- **Output:** A dict with `status`, `version`, `command`, `outputs` (new files or database prefixes), `entry_counts` or `line_counts`, `log_tail` and `provenance`; some functions add `records`, `rows` or `clusters`.
- **Example code:**

```python
import dnakit.alignment.mmseqs2_p07_align_cli_functions as al

r = al.mmseqs_extractorfs("work/nt", "work/orfs", translate=True, min_length=30, verbose=False)
print(r["entry_counts"])
```

- **Example output:**

```text
{'orf_db': 30}
```

## 12) `dnakit.alignment.mmseqs_extractframes()`

- **Purpose:** Extract the reading frames of nucleotide sequences, translated or not (`mmseqs extractframes`).
- **Method:** Runs the MMseqs2 module with a fixed argument list and no shell; options you do not set are not passed, so MMseqs2's defaults apply; output prefixes and files must not exist yet.
- **API:** `dnakit.alignment.mmseqs_extractframes(sequence_db, frame_db, *, forward_frames=None, reverse_frames=None, translation_table=None, translate=None, ..., threads=1, mmseqs_path=None, verbose=True)`.
- **Input:** `sequence_db`: nucleotide sequence database.; `frame_db`: new database of the frames..
- **Output:** A dict with `status`, `version`, `command`, `outputs` (new files or database prefixes), `entry_counts` or `line_counts`, `log_tail` and `provenance`; some functions add `records`, `rows` or `clusters`.
- **Example code:**

```python
import dnakit.alignment.mmseqs2_p07_align_cli_functions as al

r = al.mmseqs_extractframes("work/nt", "work/frames", forward_frames="1,2", reverse_frames="1", translate=True, verbose=False)
print(r["entry_counts"])
```

- **Example output:**

```text
{'frame_db': 12}
```

## 13) `dnakit.alignment.mmseqs_orftocontig()`

- **Purpose:** Write the ORF locations on their contigs in alignment format (`mmseqs orftocontig`).
- **Method:** Runs the MMseqs2 module with a fixed argument list and no shell; options you do not set are not passed, so MMseqs2's defaults apply; output prefixes and files must not exist yet.
- **API:** `dnakit.alignment.mmseqs_orftocontig(contig_db, orf_db, alignment_db, *, , threads=1, mmseqs_path=None, verbose=True)`.
- **Input:** `contig_db`: nucleotide sequence database of the contigs.; `orf_db`: database of the extracted ORFs (extractorfs output); its header database <orf_db>_h is read.; `alignment_db`: new alignment database..
- **Output:** A dict with `status`, `version`, `command`, `outputs` (new files or database prefixes), `entry_counts` or `line_counts`, `log_tail` and `provenance`; some functions add `records`, `rows` or `clusters`.
- **Example code:**

```python
import dnakit.alignment.mmseqs2_p07_align_cli_functions as al

r = al.mmseqs_orftocontig("work/nt", "work/orfs", "work/orf2contig", verbose=False)
print(r["entry_counts"])
```

- **Example output:**

```text
{'alignment_db': 30}
```

## 14) `dnakit.alignment.mmseqs_recoverlongestorf()`

- **Purpose:** Recover the longest ORF per query contig of a result for taxonomic assignment (`mmseqs recoverlongestorf`).
- **Method:** Runs the MMseqs2 module with a fixed argument list and no shell; options you do not set are not passed, so MMseqs2's defaults apply; output prefixes and files must not exist yet.
- **API:** `dnakit.alignment.mmseqs_recoverlongestorf(orf_db, result_db, output, *, , threads=1, mmseqs_path=None, verbose=True)`.
- **Input:** `orf_db`: ORF database.; `result_db`: result database of the ORFs.; `output`: new TSV file..
- **Output:** A dict with `status`, `version`, `command`, `outputs` (new files or database prefixes), `entry_counts` or `line_counts`, `log_tail` and `provenance`; some functions add `records`, `rows` or `clusters`.
- **Example code:**

```python
import dnakit.alignment.mmseqs2_p07_align_cli_functions as al

al.mmseqs_prefilter("work/orfs", "work/prot", "work/orfpf", verbose=False)
al.mmseqs_rescorediagonal("work/orfs", "work/prot", "work/orfpf", "work/orfaln", rescore_mode=2, verbose=False)
r = al.mmseqs_recoverlongestorf("work/orfs", "work/orfaln", "work/longest.tsv", verbose=False)
print(r["line_counts"])
```

- **Example output:**

```text
{'output': 0}
```

## 15) `dnakit.alignment.mmseqs_splitsequence()`

- **Purpose:** Split long sequences into overlapping pieces (`mmseqs splitsequence`).
- **Method:** Runs the MMseqs2 module with a fixed argument list and no shell; options you do not set are not passed, so MMseqs2's defaults apply; output prefixes and files must not exist yet.
- **API:** `dnakit.alignment.mmseqs_splitsequence(sequence_db, split_db, *, sequence_overlap=None, max_seq_len=None, sequence_split_mode=None, headers_split_mode=None, ..., threads=1, mmseqs_path=None, verbose=True)`.
- **Input:** `sequence_db`: sequence database.; `split_db`: new sequence database of the pieces..
- **Output:** A dict with `status`, `version`, `command`, `outputs` (new files or database prefixes), `entry_counts` or `line_counts`, `log_tail` and `provenance`; some functions add `records`, `rows` or `clusters`.
- **Example code:**

```python
import dnakit.alignment.mmseqs2_p07_align_cli_functions as al

r = al.mmseqs_splitsequence("work/nt", "work/split", max_seq_len=300, sequence_overlap=30, verbose=False)
print(r["entry_counts"])
```

- **Example output:**

```text
{'split_db': 8}
```

## 16) `dnakit.alignment.mmseqs_masksequence()`

- **Purpose:** Soft-mask low-complexity regions of a sequence database with tantan (`mmseqs masksequence`).
- **Method:** Runs the MMseqs2 module with a fixed argument list and no shell; options you do not set are not passed, so MMseqs2's defaults apply; output prefixes and files must not exist yet.
- **API:** `dnakit.alignment.mmseqs_masksequence(sequence_db, masked_db, *, mask=None, mask_prob=None, mask_lower_case=None, mask_n_repeat=None, threads=1, mmseqs_path=None, verbose=True)`.
- **Input:** `sequence_db`: sequence database.; `masked_db`: new sequence database with the masked residues in lower case..
- **Output:** A dict with `status`, `version`, `command`, `outputs` (new files or database prefixes), `entry_counts` or `line_counts`, `log_tail` and `provenance`; some functions add `records`, `rows` or `clusters`.
- **Example code:**

```python
import dnakit.alignment.mmseqs2_p07_align_cli_functions as al

r = al.mmseqs_masksequence("work/prot", "work/masked", verbose=False)
print(r["entry_counts"])
```

- **Example output:**

```text
{'masked_db': 12}
```

## 17) `dnakit.alignment.mmseqs_offsetalignment()`

- **Purpose:** Shift the alignments of ORFs or split sequences by their offset on the original sequence (`mmseqs offsetalignment`).
- **Method:** Runs the MMseqs2 module with a fixed argument list and no shell; options you do not set are not passed, so MMseqs2's defaults apply; output prefixes and files must not exist yet.
- **API:** `dnakit.alignment.mmseqs_offsetalignment(query_db, query_orf_db, target_db, target_orf_db, alignment_db, offset_db, *, search_type=None, chain_alignments=None, merge_query=None, threads=1, mmseqs_path=None, verbose=True)`.
- **Input:** `query_db`: original query sequence database.; `query_orf_db`: query ORF or piece database.; `target_db`: original target sequence database.; `target_orf_db`: target ORF or piece database.; `alignment_db`: alignment database of the ORFs.; `offset_db`: new alignment database in original coordinates..
- **Output:** A dict with `status`, `version`, `command`, `outputs` (new files or database prefixes), `entry_counts` or `line_counts`, `log_tail` and `provenance`; some functions add `records`, `rows` or `clusters`.
- **Example code:**

```python
import dnakit.alignment.mmseqs2_p07_align_cli_functions as al

r = al.mmseqs_offsetalignment("work/nt", "work/orfs", "work/prot", "work/prot", "work/orfaln", "work/offset", verbose=False)
print(r["entry_counts"])
```

- **Example output:**

```text
{'offset_db': 4}
```

## 18) `dnakit.alignment.mmseqs_proteinaln2nucl()`

- **Purpose:** Map protein alignments back to nucleotide alignments (`mmseqs proteinaln2nucl`).
- **Method:** Runs the MMseqs2 module with a fixed argument list and no shell; options you do not set are not passed, so MMseqs2's defaults apply; output prefixes and files must not exist yet.
- **API:** `dnakit.alignment.mmseqs_proteinaln2nucl(nucleotide_query_db, nucleotide_target_db, protein_query_db, protein_target_db, alignment_db, nucleotide_alignment_db, *, gap_open=None, gap_extend=None, threads=1, mmseqs_path=None, verbose=True)`.
- **Input:** `nucleotide_query_db`: nucleotide query database.; `nucleotide_target_db`: nucleotide target database.; `protein_query_db`: protein query database (translated).; `protein_target_db`: protein target database (translated).; `alignment_db`: protein alignment database.; `nucleotide_alignment_db`: new nucleotide alignment database..
- **Output:** A dict with `status`, `version`, `command`, `outputs` (new files or database prefixes), `entry_counts` or `line_counts`, `log_tail` and `provenance`; some functions add `records`, `rows` or `clusters`.
- **Example code:**

```python
import dnakit.alignment.mmseqs2_p07_align_cli_functions as al

al.mmseqs_prefilter("work/prot0", "work/prot0", "work/pf0", verbose=False)
al.mmseqs_rescorediagonal("work/prot0", "work/prot0", "work/pf0", "work/aln0", rescore_mode=2, add_backtrace=True, verbose=False)
r = al.mmseqs_proteinaln2nucl("work/cds", "work/cds", "work/prot0", "work/prot0", "work/aln0", "work/nucl_aln", verbose=False)
print(r["entry_counts"])
```

- **Example output:**

```text
{'nucleotide_alignment_db': 12}
```
