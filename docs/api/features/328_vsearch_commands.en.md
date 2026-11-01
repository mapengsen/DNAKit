# VSEARCH: search, clustering, chimeras, dereplication, read filtering and merging

This group of functions in `dnakit.fasta_fastq` is a set of typed adapters over [VSEARCH](https://github.com/torognes/vsearch) (tested with 2.32.0). Each function runs one VSEARCH command with a fixed, typed list of options and without a shell; the options of the commands that differ only by a parameter (the chimera detection methods, the dereplication variants, the read filters, the UDB reports) are one function with that parameter. They cover global alignment and search, UNOISE3 denoising, chimera detection, dereplication and rereplication, sorting, FASTQ statistics, read filtering, pair merging, joining and synchronising, record and subsequence extraction, orientation, restriction cutting, scrambling, SFF conversion, SINTAX classification and UDB indexes. They do not expose every VSEARCH option, and they do not cluster (`dnakit.fasta_fastq.cluster_vsearch()` does).

VSEARCH is dual-licensed (BSD-2-Clause or GPL-3.0). DNAKit runs it as a separate process and copies none of its code.

The functions are in the module `dnakit.fasta_fastq.vsearch_p09_cli_functions`. VSEARCH is not a Python package: install it with `conda install -c bioconda vsearch`, `apt install vsearch` or a release binary, and put it on PATH or pass `vsearch_path`. If it is missing, the functions raise `BackendUnavailableError` with an install hint, before any output directory is created. Every function writes into a new `output_directory` (it must not exist yet) under documented file names, validates its inputs (`ConfigurationError`) and raises `BackendExecutionError` with the log tail when VSEARCH fails. Optional parameters default to `None` and are then not passed, so the VSEARCH defaults apply. VSEARCH rejects an option that its command does not accept.

## Example data

All examples read these small files from the working directory. `toy.fa` has two different 90-base sequences with abundance annotations, a variant of the first, a chimera of the two and an exact duplicate of the first. `site.fa` has one record with an EcoRI site (`GAATTC`). `tax.fa` is a two-record taxonomy database. `reads.fq` has four reads whose last bases have low quality, and `r1.fq`/`r2.fq` are three read pairs with a 30-base overlap (`r1_orphan.fq` has one extra forward read). Every example puts its output into a new temporary directory, and the examples assume that `vsearch` is on PATH.

**Example file:** `toy.fa`

```text
>a;size=50
TTTCCTCATGCAATTCAAAACCATGTCCGTAATGTAGGCGAAATAGTAAACCATTTTACGGAGGATACCAAATTCCTCCTTATTCAGGAC
>b;size=30
CTAACCTGAGGTAAACCAGGTCTCTCCGCCCCCTTATAAAAGCTGTTGCACCTAGCCAAGTTCAACGGCAGCTGCAATGGAAATAGGCAA
>a2;size=5
TTTCCTCATGCAATTCAAAACCATGTCCGTAATGTAGGCGAAATAGTAAACCATTTTACGGAGGATACCAAATTCCTCCTTATTCAGGAT
>chim;size=3
TTTCCTCATGCAATTCAAAACCATGTCCGTAATGTAGGCGAAATATTGCACCTAGCCAAGTTCAACGGCAGCTGCAATGGAAATAGGCAA
>dup1;size=2
TTTCCTCATGCAATTCAAAACCATGTCCGTAATGTAGGCGAAATAGTAAACCATTTTACGGAGGATACCAAATTCCTCCTTATTCAGGAC
```

**Example file:** `site.fa`

```text
>withsite
TTTCCTCATGCAATTCAAAACCATGTCCGTAATGTAGGCGGAATTCCTAACCTGAGGTAAACCAGGTCTCTCCGCCCCCTTATAAA
>nosite
CTAACCTGAGGTAAACCAGGTCTCTCCGCCCCCTTATAAAAGCTGTTGCA
```

**Example file:** `tax.fa`

```text
>t1;tax=d:Bacteria,p:Firmicutes,c:Bacilli
TTTCCTCATGCAATTCAAAACCATGTCCGTAATGTAGGCGAAATAGTAAACCATTTTACGGAGGATACCAAATTCCTCCTTATTCAGGAC
>t2;tax=d:Bacteria,p:Proteobacteria,c:Gamma
CTAACCTGAGGTAAACCAGGTCTCTCCGCCCCCTTATAAAAGCTGTTGCACCTAGCCAAGTTCAACGGCAGCTGCAATGGAAATAGGCAA
```

**Example file:** `reads.fq`

```text
@read0
TGACGGATATATATTAAAAAGTGTTTTAAGATACATTGAG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII((((((((((
@read1
GCCCGTTCGTGCTCCTCGCCCTGAAGCATTGCTTTGTGAAGAGGGACTTC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII((((((((((((((((((((
@read2
AGCCAATAGACCTGCATACCGGCTCATTCTTCATGTGCAACCTAGGGAGAATGTGTACAT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII((((((((((((((((((((((((((((((
@read3
ACGCTCTTACTGCGGTCGCGTCTAATAATATACATTTGCTTCGTTGACTAGCAACCCAGGGCTATAGCTA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII((((((((((((((((((((((((((((((((((((((((
```

**Example file:** `r1.fq`

```text
@pair0/1
TTCCCCCCGCGGCCCACCCAGTATTCCTAACGGAGCATAAATCCCACCCG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@pair1/1
GGTCCAAGATCGGGACTCGGTCTCCAGGTAAGACGGGCTCATTCATAAAC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@pair2/1
TCTATTTGTGGGTGGGAACACTTAGTAGACTTGCAATCCAATTACAGCAG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
```

**Example file:** `r2.fq`

```text
@pair0/2
AAGGTTCGACAAACTTAGTTCGGGTGGGATTTATGCTCCGTTAGGAATAC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@pair1/2
AGATTATACCCCTTAGTAACGTTTATGAATGAGCCCGTCTTACCTGGAGA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@pair2/2
GCGCCCCTAGGCGCACAAGACTGCTGTAATTGGATTGCAAGTCTACTAAG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
```

**Example file:** `r1_orphan.fq`

```text
@pair0/1
TTCCCCCCGCGGCCCACCCAGTATTCCTAACGGAGCATAAATCCCACCCG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@pair1/1
GGTCCAAGATCGGGACTCGGTCTCCAGGTAAGACGGGCTCATTCATAAAC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@pair2/1
TCTATTTGTGGGTGGGAACACTTAGTAGACTTGCAATCCAATTACAGCAG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@lonely/1
CCCAAAGGTAAACGAACCGTTGCGGTCAATCTTGTCGCGGCTGATGAATT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
```

## 1) `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_allpairs_global()`

- **Purpose:** Align every pair of sequences of a FASTA file globally and report the identity of each pair.
- **Method:** Run `vsearch --allpairs_global` with `--alnout` and `--blast6out`. `accept_all=True` adds `--acceptall`, so that all pairs are reported; with `accept_all=False` only pairs that reach `min_identity` (`--id`) are. `identity_definition` is `--iddef`, `match_score`, `mismatch_score`, `gap_open` and `gap_extend` are the scoring options. The BLAST-like table is parsed into `records`.
- **API:** `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_allpairs_global(input_fasta: 'str | Path', output_directory: 'str | Path', *, accept_all: 'bool' = True, min_identity: 'float | None' = None, identity_definition: 'int | None' = None, match_score: 'int | None' = None, mismatch_score: 'int | None' = None, gap_open: 'str | None' = None, gap_extend: 'str | None' = None, threads: 'int | None' = None, vsearch_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** A FASTA file, a new output directory, and the optional parameters above.
- **Output:** A dictionary with `status`, `version`, `command`, `outputs` (`alignments` = `alignments.aln`, `hits` = `hits.tsv`), `records` (query, target, identity, alignment length, mismatches, gap opens, positions, e-value, bit score), `record_count` and `log_tail`.
- **Example code:**

```python
import tempfile
from dnakit.fasta_fastq.vsearch_p09_cli_functions import vsearch_allpairs_global

result = vsearch_allpairs_global("toy.fa", tempfile.mkdtemp() + "/run", verbose=False)
best = max(result["records"], key=lambda row: row["identity"])
print(result["status"], result["record_count"])
print(best["query"], best["target"], best["identity"])
```

- **Example output:**

```text
completed 10
a;size=50 dup1;size=2 100.0
```

## 2) `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_search_global()`

- **Purpose:** Search query sequences against a FASTA or UDB database by global alignment.
- **Method:** Run `vsearch --usearch_global` with `--db`, the required `--id` (`min_identity`) and a BLAST-like table. Search settings map to `--strand`, `--iddef`, `--maxaccepts`, `--maxrejects`, `--minwordmatches`, `--wordlength`, `--query_cov`, `--target_cov`, `--top_hits_only`, `--output_no_hits`, `--qmask`, `--dbmask` and `--sizein`. The flags `write_alignments`, `write_uc`, `write_otu_table`, `write_sam` and `write_matched_queries` add further output files. The single-sequence and batch search calls of the VSEARCH library API are this operation.
- **API:** `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_search_global(query_fasta: 'str | Path', database: 'str | Path', output_directory: 'str | Path', *, min_identity: 'float', strand: "Literal['plus', 'both'] | None" = None, identity_definition: 'int | None' = None, max_accepts: 'int | None' = None, max_rejects: 'int | None' = None, min_word_matches: 'int | None' = None, word_length: 'int | None' = None, query_coverage: 'float | None' = None, target_coverage: 'float | None' = None, top_hits_only: 'bool' = False, output_no_hits: 'bool' = False, query_mask: "Literal['none', 'dust', 'soft'] | None" = None, database_mask: "Literal['none', 'dust', 'soft'] | None" = None, abundance_in: 'bool' = False, write_alignments: 'bool' = False, write_uc: 'bool' = False, write_otu_table: 'bool' = False, write_sam: 'bool' = False, write_matched_queries: 'bool' = False, threads: 'int | None' = None, vsearch_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** Query FASTA, database (FASTA or UDB), a new output directory, and `min_identity` between 0 and 1.
- **Output:** `status`, `outputs` (`hits` and the optional files), `records` (BLAST-like rows), `record_count`, plus the common fields.
- **Example code:**

```python
import tempfile
from dnakit.fasta_fastq.vsearch_p09_cli_functions import vsearch_search_global

result = vsearch_search_global("toy.fa", "toy.fa", tempfile.mkdtemp() + "/run", min_identity=0.9, abundance_in=True, write_otu_table=True, verbose=False)
for row in sorted(result["records"], key=lambda row: row["query"]):
    print(row["query"], row["target"], row["identity"])
print(sorted(result["outputs"]))
```

- **Example output:**

```text
a2;size=5 a2;size=5 100.0
a;size=50 a;size=50 100.0
b;size=30 b;size=30 100.0
chim;size=3 chim;size=3 100.0
dup1;size=2 a;size=50 100.0
['hits', 'otu_table']
```

## 3) `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_cluster_unoise()`

- **Purpose:** Denoise amplicon reads with the UNOISE3 algorithm and keep one centroid per real biological sequence.
- **Method:** Run `vsearch --cluster_unoise` with `--centroids` and `--uc`. `min_abundance` is `--minsize` (VSEARCH default 8), `unoise_alpha` is `--unoise_alpha` (2.0), `identity` is `--id`; `abundance_in`/`abundance_out` handle the `;size=` annotation and `relabel` renames the centroids. `write_otu_table` adds `--otutabout`.
- **API:** `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_cluster_unoise(input_fasta: 'str | Path', output_directory: 'str | Path', *, min_abundance: 'int | None' = None, unoise_alpha: 'float | None' = None, identity: 'float | None' = None, abundance_in: 'bool' = False, abundance_out: 'bool' = False, relabel: 'str | None' = None, write_otu_table: 'bool' = False, threads: 'int | None' = None, vsearch_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** A FASTA file whose labels carry `;size=` abundances and a new output directory.
- **Output:** `outputs` (`centroids` = `centroids.fasta`, `clusters` = `clusters.uc`, optional `otu_table`), `centroid_count` and the common fields.
- **Example code:**

```python
import tempfile
from dnakit.fasta_fastq.vsearch_p09_cli_functions import vsearch_cluster_unoise

result = vsearch_cluster_unoise("toy.fa", tempfile.mkdtemp() + "/run", min_abundance=2, abundance_in=True, abundance_out=True, verbose=False)
print(result["status"], result["centroid_count"])
print([line for line in open(result["outputs"]["centroids"]) if line.startswith(">")])
```

- **Example output:**

```text
completed 3
['>a;size=57\n', '>b;size=30\n', '>chim;size=3\n']
```

## 4) `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_detect_chimeras()`

- **Purpose:** Detect chimeric sequences with UCHIME, UCHIME2/3, UCHIME against a reference, or the new de novo algorithm.
- **Method:** `method` selects the VSEARCH command: `uchime_denovo`, `uchime2_denovo`, `uchime3_denovo`, `uchime_ref` (needs `database`) or `chimeras_denovo` (long exact sequences; options `chimeras_diff_pct`, `chimeras_length_min`, `chimeras_parents_max`, `chimeras_parts`). `abundance_skew` is `--abskew`; `min_score`, `min_divergence`, `min_diffs`, `no_vote_pseudocount`, `no_vote_weight` are `--minh`, `--mindiv`, `--mindiffs`, `--dn`, `--xn`. Parameters that belong to another method are refused before the program runs. The 18-column report is parsed into `records`.
- **API:** `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_detect_chimeras(input_fasta: 'str | Path', output_directory: 'str | Path', *, method: "Literal['uchime_denovo', 'uchime2_denovo', 'uchime3_denovo', 'uchime_ref', 'chimeras_denovo']" = 'uchime_denovo', database: 'str | Path | None' = None, abundance_skew: 'float | None' = None, abundance_in: 'bool' = False, abundance_out: 'bool' = False, min_score: 'float | None' = None, min_divergence: 'float | None' = None, min_diffs: 'int | None' = None, no_vote_pseudocount: 'float | None' = None, no_vote_weight: 'float | None' = None, exclude_same_label: 'bool' = False, exclude_same_sequence: 'bool' = False, chimeras_diff_pct: 'float | None' = None, chimeras_length_min: 'int | None' = None, chimeras_parents_max: 'int | None' = None, chimeras_parts: 'int | None' = None, relabel: 'str | None' = None, threads: 'int | None' = None, vsearch_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** A FASTA file with `;size=` abundances (`abundance_in=True`), a new output directory and, for `uchime_ref`, a reference FASTA.
- **Output:** `outputs` (`chimeras`, `nonchimeras`, `report`, `alignments`, and `borderline` for the UCHIME methods), `records` (score, query, parents, identities, votes, divergence, `chimeric` = Y/N/?), `chimera_count`, `nonchimera_count` and `method`.
- **Example code:**

```python
import tempfile
from dnakit.fasta_fastq.vsearch_p09_cli_functions import vsearch_detect_chimeras

result = vsearch_detect_chimeras("toy.fa", tempfile.mkdtemp() + "/run", method="uchime_denovo", abundance_in=True, verbose=False)
print(result["chimera_count"], result["nonchimera_count"])
print([row["query"] for row in result["records"] if row["chimeric"] == "Y"])
```

- **Example output:**

```text
1 4
['chim;size=3']
```

## 5) `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_cut_restriction_sites()`

- **Purpose:** Cut sequences at a restriction-enzyme site.
- **Method:** Run `vsearch --cut` with `--cut_pattern`, a site written with `^` and `_` marking the cut on the forward and the reverse strand (for example `G^AATT_C`). Fragments are written for the forward and the reverse strand; sequences without the site are written to the discarded files.
- **API:** `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_cut_restriction_sites(input_fasta: 'str | Path', output_directory: 'str | Path', *, cut_pattern: 'str', vsearch_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** A FASTA file, a new output directory and `cut_pattern` (nucleotide letters including both `^` and `_`).
- **Output:** `outputs` (`fragments_forward`, `fragments_reverse`, `discarded_forward`, `discarded_reverse`), `forward_fragment_count`, `reverse_fragment_count`.
- **Example code:**

```python
import tempfile
from dnakit.fasta_fastq.vsearch_p09_cli_functions import vsearch_cut_restriction_sites

result = vsearch_cut_restriction_sites("site.fa", tempfile.mkdtemp() + "/run", cut_pattern="G^AATT_C", verbose=False)
print(result["forward_fragment_count"], result["reverse_fragment_count"])
print([line.strip() for line in open(result["outputs"]["fragments_forward"]) if line.startswith(">")])
```

- **Example output:**

```text
2 2
['>withsite', '>withsite']
```

## 6) `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_dereplicate_fulllength()`

- **Purpose:** Collapse identical sequences and count them (full-length dereplication).
- **Method:** `method='fulllength'` runs `--derep_fulllength`, `'smallmem'` runs `--derep_smallmem` (less memory, no UC file) and `'fastx_uniques'` runs `--fastx_uniques` for FASTA or FASTQ input (`fastq_output=True` keeps the qualities). `strand`, `abundance_in`, `abundance_out`, `min_unique_size`, `max_unique_size`, `top_n`, `relabel` and `relabel_mode` map to `--strand`, `--sizein`, `--sizeout`, `--minuniquesize`, `--maxuniquesize`, `--topn`, `--relabel` and the digest relabel switches.
- **API:** `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_dereplicate_fulllength(input_fasta: 'str | Path', output_directory: 'str | Path', *, method: "Literal['fulllength', 'smallmem', 'fastx_uniques']" = 'fulllength', fastq_output: 'bool' = False, strand: "Literal['plus', 'both'] | None" = None, abundance_in: 'bool' = False, abundance_out: 'bool' = False, min_unique_size: 'int | None' = None, max_unique_size: 'int | None' = None, top_n: 'int | None' = None, relabel: 'str | None' = None, relabel_mode: "Literal['keep', 'md5', 'sha1', 'self'] | None" = None, vsearch_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** A FASTA (or FASTQ for `fastx_uniques`) file and a new output directory.
- **Output:** `outputs` (`unique` and, except for `smallmem`, `clusters`) and `unique_count`.
- **Example code:**

```python
import tempfile
from dnakit.fasta_fastq.vsearch_p09_cli_functions import vsearch_dereplicate_fulllength

result = vsearch_dereplicate_fulllength("toy.fa", tempfile.mkdtemp() + "/run", abundance_in=True, abundance_out=True, verbose=False)
print(result["unique_count"])
print([line.strip() for line in open(result["outputs"]["unique"]) if line.startswith(">")])
```

- **Example output:**

```text
4
['>a;size=52', '>b;size=30', '>a2;size=5', '>chim;size=3']
```

## 7) `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_dereplicate_by_id()`

- **Purpose:** Dereplicate on both the identifier and the sequence.
- **Method:** Run `vsearch --derep_id`. Only records with the same label and the same sequence are collapsed, so equal sequences with different labels stay apart. The options are the dereplication options of the previous function.
- **API:** `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_dereplicate_by_id(input_fasta: 'str | Path', output_directory: 'str | Path', *, strand: "Literal['plus', 'both'] | None" = None, abundance_in: 'bool' = False, abundance_out: 'bool' = False, min_unique_size: 'int | None' = None, max_unique_size: 'int | None' = None, top_n: 'int | None' = None, relabel: 'str | None' = None, relabel_mode: "Literal['keep', 'md5', 'sha1', 'self'] | None" = None, vsearch_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** A FASTA file and a new output directory.
- **Output:** `outputs` (`unique`, `clusters`) and `unique_count`.
- **Example code:**

```python
import tempfile
from dnakit.fasta_fastq.vsearch_p09_cli_functions import vsearch_dereplicate_by_id

result = vsearch_dereplicate_by_id("toy.fa", tempfile.mkdtemp() + "/run", verbose=False)
print(result["unique_count"])
```

- **Example output:**

```text
5
```

## 8) `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_dereplicate_prefix()`

- **Purpose:** Merge sequences that are prefixes of longer sequences.
- **Method:** Run `vsearch --derep_prefix`. A sequence identical to the start of a longer one is counted with the longer sequence. `threads` is passed as `--threads`.
- **API:** `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_dereplicate_prefix(input_fasta: 'str | Path', output_directory: 'str | Path', *, strand: "Literal['plus', 'both'] | None" = None, abundance_in: 'bool' = False, abundance_out: 'bool' = False, min_unique_size: 'int | None' = None, max_unique_size: 'int | None' = None, top_n: 'int | None' = None, relabel: 'str | None' = None, relabel_mode: "Literal['keep', 'md5', 'sha1', 'self'] | None" = None, threads: 'int | None' = None, vsearch_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** A FASTA file and a new output directory.
- **Output:** `outputs` (`unique`, `clusters`) and `unique_count`.
- **Example code:**

```python
import tempfile
from dnakit.fasta_fastq.vsearch_p09_cli_functions import vsearch_dereplicate_prefix

prefix = tempfile.mkdtemp() + "/prefixes.fa"
lines = open("toy.fa").read().splitlines()
open(prefix, "w").write(">long\n" + lines[1] + "\n>short\n" + lines[1][:60] + "\n>other\n" + lines[3] + "\n")
result = vsearch_dereplicate_prefix(prefix, tempfile.mkdtemp() + "/run", abundance_out=True, verbose=False)
print(result["unique_count"])
print([line.strip() for line in open(result["outputs"]["unique"]) if line.startswith(">")])
```

- **Example output:**

```text
2
['>long;size=2', '>other;size=1']
```

## 9) `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_rereplicate()`

- **Purpose:** Expand every record into as many copies as its `;size=` abundance says.
- **Method:** Run `vsearch --rereplicate`. `strip_abundance=True` adds `--xsize` to remove the size annotation from the labels.
- **API:** `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_rereplicate(input_fasta: 'str | Path', output_directory: 'str | Path', *, strip_abundance: 'bool' = False, vsearch_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** A FASTA file with `;size=` annotations and a new output directory.
- **Output:** `outputs` (`rereplicated`) and `record_count`.
- **Example code:**

```python
import tempfile
from dnakit.fasta_fastq.vsearch_p09_cli_functions import vsearch_rereplicate

result = vsearch_rereplicate("toy.fa", tempfile.mkdtemp() + "/run", verbose=False)
print(result["record_count"])
```

- **Example output:**

```text
90
```

## 10) `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_sort_by_size()`

- **Purpose:** Sort FASTA records by decreasing abundance.
- **Method:** Run `vsearch --sortbysize`. `min_size` and `max_size` drop records outside the abundance range, `top_n` keeps the first n records, `abundance_out` writes `;size=` and `relabel`/`relabel_mode` rename the records.
- **API:** `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_sort_by_size(input_fasta: 'str | Path', output_directory: 'str | Path', *, min_size: 'int | None' = None, max_size: 'int | None' = None, top_n: 'int | None' = None, abundance_out: 'bool' = False, relabel: 'str | None' = None, relabel_mode: "Literal['keep', 'md5', 'sha1', 'self'] | None" = None, vsearch_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** A FASTA file with `;size=` annotations and a new output directory.
- **Output:** `outputs` (`sorted`) and `record_count`.
- **Example code:**

```python
import tempfile
from dnakit.fasta_fastq.vsearch_p09_cli_functions import vsearch_sort_by_size

result = vsearch_sort_by_size("toy.fa", tempfile.mkdtemp() + "/run", min_size=3, verbose=False)
print(result["record_count"])
print([line.strip() for line in open(result["outputs"]["sorted"]) if line.startswith(">")])
```

- **Example output:**

```text
4
['>a;size=50', '>b;size=30', '>a2;size=5', '>chim;size=3']
```

## 11) `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_fastq_eestats()`

- **Purpose:** Per-position quality-score and expected-error statistics of a FASTQ file.
- **Method:** Run `vsearch --fastq_eestats` and parse the table. `fastq_ascii`, `fastq_qmin` and `fastq_qmax` set the quality encoding.
- **API:** `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_fastq_eestats(input_fastq: 'str | Path', output_directory: 'str | Path', *, fastq_ascii: 'int | None' = None, fastq_qmin: 'int | None' = None, fastq_qmax: 'int | None' = None, vsearch_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** A FASTQ file and a new output directory.
- **Output:** `outputs` (`statistics` = `eestats.tsv`), `records` (one row per position with `pos`, `recs`, and min/low/median/mean/high/max of Q, Pe and EE) and `position_count`.
- **Example code:**

```python
import tempfile
from dnakit.fasta_fastq.vsearch_p09_cli_functions import vsearch_fastq_eestats

result = vsearch_fastq_eestats("reads.fq", tempfile.mkdtemp() + "/run", verbose=False)
print(result["position_count"])
print(result["records"][0]["pos"], result["records"][0]["recs"], result["records"][-1]["recs"])
```

- **Example output:**

```text
70
1 4 1
```

## 12) `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_fastq_eestats2()`

- **Purpose:** How many reads survive each length cutoff and maximum expected error.
- **Method:** Run `vsearch --fastq_eestats2`. `ee_cutoffs` are the MaxEE values (default 0.5, 1, 2) and `length_cutoffs` is `min,max,increment` (default `50,*,50`).
- **API:** `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_fastq_eestats2(input_fastq: 'str | Path', output_directory: 'str | Path', *, ee_cutoffs: 'Sequence[float] | None' = None, length_cutoffs: 'str | None' = None, fastq_ascii: 'int | None' = None, fastq_qmin: 'int | None' = None, fastq_qmax: 'int | None' = None, vsearch_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** A FASTQ file, a new output directory and the optional cutoffs.
- **Output:** `outputs` (`statistics`), `records` (per length cutoff a `retained` list with `max_ee`, `reads`, `percent`), and the summary `reads`, `max_length`, `mean_length`.
- **Example code:**

```python
import tempfile
from dnakit.fasta_fastq.vsearch_p09_cli_functions import vsearch_fastq_eestats2

result = vsearch_fastq_eestats2("reads.fq", tempfile.mkdtemp() + "/run", ee_cutoffs=[1, 3], length_cutoffs="30,70,20", verbose=False)
print(result["reads"], result["max_length"])
for row in result["records"]:
    print(row["length"], [item["reads"] for item in row["retained"]])
```

- **Example output:**

```text
4 70
30 [4, 4]
50 [0, 0]
70 [0, 0]
```

## 13) `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_fastq_stats()`

- **Purpose:** Detailed statistics of a FASTQ file: read lengths, quality scores and truncation.
- **Method:** Run `vsearch --fastq_stats`; the report goes to `stats.txt` through `--log` and the read-length and quality-score tables are parsed from it.
- **API:** `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_fastq_stats(input_fastq: 'str | Path', output_directory: 'str | Path', *, fastq_ascii: 'int | None' = None, fastq_qmin: 'int | None' = None, fastq_qmax: 'int | None' = None, vsearch_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** A FASTQ file and a new output directory.
- **Output:** `outputs` (`report`), `read_length_distribution`, `quality_distribution`, `read_count`, `mean_length` and `bases`.
- **Example code:**

```python
import tempfile
from dnakit.fasta_fastq.vsearch_p09_cli_functions import vsearch_fastq_stats

result = vsearch_fastq_stats("reads.fq", tempfile.mkdtemp() + "/run", verbose=False)
print(result["read_count"], result["mean_length"])
print(result["read_length_distribution"][0])
```

- **Example output:**

```text
4 55.0
{'length': 70, 'reads': 1, 'percent': 25.0, 'cumulative_percent': 25.0}
```

## 14) `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_filter_reads()`

- **Purpose:** Trim and filter FASTA or FASTQ records by quality, expected error, length, N content and abundance.
- **Method:** Run `vsearch --fastx_filter`. One function covers the criteria `min_quality` (`--fastq_minqual`), `truncate_quality` (`--fastq_truncqual`), `truncate_ee` and `truncate_ee_rate` (`--fastq_truncee`, `--fastq_truncee_rate`), `max_ee`, `max_ee_rate`, `min_length`, `max_length`, `max_ns`, `strip_left`, `strip_right`, `truncate_length`, `truncate_length_keep`, and the abundance filters `min_size` and `max_size`. `reverse_path` filters read pairs together, `write_discarded` keeps the rejected records.
- **API:** `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_filter_reads(input_path: 'str | Path', output_directory: 'str | Path', *, output_format: "Literal['fasta', 'fastq'] | None" = None, reverse_path: 'str | Path | None' = None, min_quality: 'int | None' = None, truncate_quality: 'int | None' = None, truncate_ee: 'float | None' = None, truncate_ee_rate: 'float | None' = None, max_ee: 'float | None' = None, max_ee_rate: 'float | None' = None, min_length: 'int | None' = None, max_length: 'int | None' = None, max_ns: 'int | None' = None, strip_left: 'int | None' = None, strip_right: 'int | None' = None, truncate_length: 'int | None' = None, truncate_length_keep: 'int | None' = None, min_size: 'int | None' = None, max_size: 'int | None' = None, write_discarded: 'bool' = False, include_expected_errors: 'bool' = False, relabel: 'str | None' = None, relabel_mode: "Literal['keep', 'md5', 'sha1', 'self'] | None" = None, abundance_out: 'bool' = False, fastq_ascii: 'int | None' = None, fastq_qmin: 'int | None' = None, fastq_qmax: 'int | None' = None, vsearch_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** A FASTQ or FASTA file, a new output directory and the criteria. The output format follows the input unless `output_format` is given.
- **Output:** `outputs` (`filtered`, optionally `filtered_reverse`, `discarded`, `discarded_reverse`), `passed_count` and, with `write_discarded`, `discarded_count`.
- **Example code:**

```python
import tempfile
from dnakit.fasta_fastq.vsearch_p09_cli_functions import vsearch_filter_reads

result = vsearch_filter_reads("reads.fq", tempfile.mkdtemp() + "/run", truncate_quality=20, min_length=30, write_discarded=True, verbose=False)
print(result["passed_count"], result["discarded_count"])
lengths = [len(line.strip()) for line in open(result["outputs"]["filtered"]).readlines()[1::4]]
print(lengths)
```

- **Example output:**

```text
4 0
[30, 30, 30, 30]
```

## 15) `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_sync_pairs()`

- **Purpose:** Put the reads of two paired-end files in the same order and set aside the reads without a mate.
- **Method:** Run `vsearch --fastx_syncpairs` with the four output files for synced and orphan reads. `read_separators` are the characters before the mate number (default `/`).
- **API:** `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_sync_pairs(forward_path: 'str | Path', reverse_path: 'str | Path', output_directory: 'str | Path', *, output_format: "Literal['fasta', 'fastq'] | None" = None, read_separators: 'str | None' = None, vsearch_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** Forward and reverse FASTQ (or FASTA) files and a new output directory.
- **Output:** `outputs` (`synced_forward`, `synced_reverse`, `orphans_forward`, `orphans_reverse`), `synced_pairs`, `forward_orphans`, `reverse_orphans`.
- **Example code:**

```python
import tempfile
from dnakit.fasta_fastq.vsearch_p09_cli_functions import vsearch_sync_pairs

result = vsearch_sync_pairs("r1_orphan.fq", "r2.fq", tempfile.mkdtemp() + "/run", verbose=False)
print(result["synced_pairs"], result["forward_orphans"], result["reverse_orphans"])
```

- **Example output:**

```text
3 1 0
```

## 16) `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_merge_pairs()`

- **Purpose:** Merge overlapping paired-end reads into one sequence.
- **Method:** Run `vsearch --fastq_mergepairs`. The parameters are `allow_stagger`, `max_diffs`, `max_diff_pct`, `max_ee`, `max_ns`, `min_length`, `min_overlap`, `min_merge_length`, `max_merge_length`, `truncate_quality`, `label_suffix`, `include_expected_errors`, the quality-range options and `write_error_table` (`--eetabbedout`).
- **API:** `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_merge_pairs(forward_fastq: 'str | Path', reverse_fastq: 'str | Path', output_directory: 'str | Path', *, allow_stagger: 'bool' = False, max_diffs: 'int | None' = None, max_diff_pct: 'float | None' = None, max_ee: 'float | None' = None, max_ns: 'int | None' = None, min_length: 'int | None' = None, min_overlap: 'int | None' = None, min_merge_length: 'int | None' = None, max_merge_length: 'int | None' = None, truncate_quality: 'int | None' = None, label_suffix: 'str | None' = None, include_expected_errors: 'bool' = False, write_error_table: 'bool' = False, fastq_ascii: 'int | None' = None, fastq_qmin: 'int | None' = None, fastq_qmax: 'int | None' = None, fastq_qmin_out: 'int | None' = None, fastq_qmax_out: 'int | None' = None, threads: 'int | None' = None, vsearch_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** Forward and reverse FASTQ files and a new output directory.
- **Output:** `outputs` (`merged`, `notmerged_forward`, `notmerged_reverse`, optional `errors`), `merged_count`, `not_merged_count`.
- **Example code:**

```python
import tempfile
from dnakit.fasta_fastq.vsearch_p09_cli_functions import vsearch_merge_pairs

result = vsearch_merge_pairs("r1.fq", "r2.fq", tempfile.mkdtemp() + "/run", verbose=False)
print(result["merged_count"], result["not_merged_count"])
print([len(line.strip()) for line in open(result["outputs"]["merged"]).readlines()[1::4]])
```

- **Example output:**

```text
3 0
[70, 70, 70]
```

## 17) `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_join_pairs()`

- **Purpose:** Join paired-end reads into one sequence with a gap between them.
- **Method:** Run `vsearch --fastq_join`. The reverse read is reverse-complemented and appended after `pad_sequence` (`--join_padgap`, default `NNNNNNNN`) with the qualities `pad_quality` (`--join_padgapq`, default `IIIIIIII`).
- **API:** `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_join_pairs(forward_fastq: 'str | Path', reverse_fastq: 'str | Path', output_directory: 'str | Path', *, pad_sequence: 'str | None' = None, pad_quality: 'str | None' = None, write_fasta: 'bool' = False, vsearch_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** Forward and reverse FASTQ files and a new output directory.
- **Output:** `outputs` (`joined`, optional `joined_fasta`) and `joined_count`.
- **Example code:**

```python
import tempfile
from dnakit.fasta_fastq.vsearch_p09_cli_functions import vsearch_join_pairs

result = vsearch_join_pairs("r1.fq", "r2.fq", tempfile.mkdtemp() + "/run", verbose=False)
sequence = open(result["outputs"]["joined"]).readlines()[1].strip()
print(result["joined_count"], len(sequence), sequence[50:58])
```

- **Example output:**

```text
3 108 NNNNNNNN
```

## 18) `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_extract_records()`

- **Purpose:** Extract records from a FASTA or FASTQ file by header label.
- **Method:** Run `vsearch --fastx_getseqs`; one `label` is the single-record `--fastx_getseq` case. Selectors are `label`, `labels_file`, `label_word`, `label_words_file` (with `label_field`); `substring_match` allows a match anywhere in the header; `write_unmatched` also writes the other records.
- **API:** `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_extract_records(input_path: 'str | Path', output_directory: 'str | Path', *, label: 'str | None' = None, labels_file: 'str | Path | None' = None, label_word: 'str | None' = None, label_words_file: 'str | Path | None' = None, label_field: 'str | None' = None, substring_match: 'bool' = False, output_format: "Literal['fasta', 'fastq'] | None" = None, write_unmatched: 'bool' = False, abundance_in: 'bool' = False, abundance_out: 'bool' = False, relabel: 'str | None' = None, relabel_mode: "Literal['keep', 'md5', 'sha1', 'self'] | None" = None, vsearch_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** A FASTA or FASTQ file, a new output directory and at least one selector.
- **Output:** `outputs` (`extracted`, optional `unmatched`) and `extracted_count`.
- **Example code:**

```python
import tempfile
from dnakit.fasta_fastq.vsearch_p09_cli_functions import vsearch_extract_records

result = vsearch_extract_records("toy.fa", tempfile.mkdtemp() + "/run", label="b;size=30", write_unmatched=True, verbose=False)
print(result["extracted_count"])
print([line.strip() for line in open(result["outputs"]["unmatched"]) if line.startswith(">")])
```

- **Example output:**

```text
1
['>a;size=50', '>a2;size=5', '>chim;size=3', '>dup1;size=2']
```

## 19) `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_extract_subsequence()`

- **Purpose:** Extract a subsequence of one record by label and positions.
- **Method:** Run `vsearch --fastx_getsubseq` with `--label`, `--subseq_start` and `--subseq_end` (1-based, inclusive). `label_suffix` is appended to the label.
- **API:** `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_extract_subsequence(input_path: 'str | Path', output_directory: 'str | Path', *, label: 'str', start: 'int', end: 'int', substring_match: 'bool' = False, output_format: "Literal['fasta', 'fastq'] | None" = None, label_suffix: 'str | None' = None, vsearch_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** A FASTA or FASTQ file, a new output directory, `label`, `start` and `end`.
- **Output:** `outputs` (`subsequence`) and `extracted_count`.
- **Example code:**

```python
import tempfile
from dnakit.fasta_fastq.vsearch_p09_cli_functions import vsearch_extract_subsequence

result = vsearch_extract_subsequence("toy.fa", tempfile.mkdtemp() + "/run", label="a;size=50", start=1, end=20, verbose=False)
print(open(result["outputs"]["subsequence"]).read().split()[1])
```

- **Example output:**

```text
TTTCCTCATGCAATTCAAAA
```

## 20) `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_orient_sequences()`

- **Purpose:** Orient sequences to the strand of a reference database.
- **Method:** Run `vsearch --orient` with the reference as `--db`; sequences that match on the opposite strand are reverse-complemented. `word_length`, `query_mask` and `database_mask` map to `--wordlength`, `--qmask` and `--dbmask`.
- **API:** `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_orient_sequences(input_path: 'str | Path', reference_fasta: 'str | Path', output_directory: 'str | Path', *, output_format: "Literal['fasta', 'fastq'] | None" = None, word_length: 'int | None' = None, query_mask: "Literal['none', 'dust', 'soft'] | None" = None, database_mask: "Literal['none', 'dust', 'soft'] | None" = None, threads: 'int | None' = None, vsearch_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** A FASTA or FASTQ file, a reference FASTA and a new output directory.
- **Output:** `outputs` (`oriented`, `undetermined`, `orientation`), `records` (query, strand `+`/`-`/`?`, forward and reverse matches), `record_count`.
- **Example code:**

```python
import tempfile
from dnakit.fasta_fastq.vsearch_p09_cli_functions import vsearch_orient_sequences

comp = str.maketrans("ACGT", "TGCA")
forward = open("toy.fa").read().splitlines()[1]
flipped = tempfile.mkdtemp() + "/flipped.fa"
open(flipped, "w").write(">flipped\n" + forward.translate(comp)[::-1] + "\n>straight\n" + forward + "\n")
result = vsearch_orient_sequences(flipped, "toy.fa", tempfile.mkdtemp() + "/run", verbose=False)
print([(row["query"], row["strand"]) for row in result["records"]])
```

- **Example output:**

```text
[('flipped', '-'), ('straight', '+')]
```

## 21) `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_scramble_sequences()`

- **Purpose:** Randomise the nucleotide order within each sequence.
- **Method:** Run `vsearch --scramble`. `kmer_size` (`--scramble_kmer`) preserves the counts of words up to that size and `random_seed` (`--randseed`) makes the result reproducible.
- **API:** `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_scramble_sequences(input_path: 'str | Path', output_directory: 'str | Path', *, output_format: "Literal['fasta', 'fastq'] | None" = None, kmer_size: 'int | None' = None, random_seed: 'int | None' = None, label_suffix: 'str | None' = None, vsearch_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** A FASTA or FASTQ file and a new output directory.
- **Output:** `outputs` (`scrambled`) and `record_count`.
- **Example code:**

```python
import tempfile
from dnakit.fasta_fastq.vsearch_p09_cli_functions import vsearch_scramble_sequences

result = vsearch_scramble_sequences("toy.fa", tempfile.mkdtemp() + "/run", random_seed=3, verbose=False)
original = open("toy.fa").read().splitlines()[1]
scrambled = "".join(open(result["outputs"]["scrambled"]).read().split(">")[1].split("\n")[1:])
print(result["record_count"], sorted(original) == sorted(scrambled), original != scrambled)
```

- **Example output:**

```text
5 True True
```

## 22) `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_convert_sff()`

- **Purpose:** Convert an SFF flowgram file (454, Ion Torrent) to FASTQ.
- **Method:** Run `vsearch --sff_convert`. `clip_ends=True` applies the clip points stored in the file (`--sff_clip`); `fastq_ascii_out`, `fastq_qmin_out` and `fastq_qmax_out` set the output quality encoding.
- **API:** `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_convert_sff(input_sff: 'str | Path', output_directory: 'str | Path', *, clip_ends: 'bool' = False, fastq_ascii_out: 'int | None' = None, fastq_qmin_out: 'int | None' = None, fastq_qmax_out: 'int | None' = None, vsearch_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** An SFF file and a new output directory.
- **Output:** `outputs` (`converted`) and `record_count`.
- **Example code:**

```python
import tempfile
from dnakit.fasta_fastq.vsearch_p09_cli_functions import vsearch_convert_sff

import struct

flows, key = "TACG" * 3, "TCAG"
head_len = 31 + len(flows) + len(key)
head_len += -head_len % 8
blob = struct.pack(">IIQIIHHHB", 0x2E736666, 1, 0, 0, 1, head_len, len(key), len(flows), 1) + flows.encode() + key.encode()
blob += bytes(head_len - len(blob))
bases = "ACGTACGTAC"
read_head = struct.pack(">HHIHHHH", 24, 2, len(bases), 1, len(bases), 0, 0) + b"r1"
read_head += bytes(24 - len(read_head))
data = struct.pack(">12H", *[100] * 12) + bytes([1] * 10) + bases.encode() + bytes([30] * 10)
data += bytes(-len(data) % 8)
sff = tempfile.mkdtemp() + "/toy.sff"
open(sff, "wb").write(blob + read_head + data)
result = vsearch_convert_sff(sff, tempfile.mkdtemp() + "/run", verbose=False)
print(open(result["outputs"]["converted"]).read().split())
```

- **Example output:**

```text
['@r1', 'ACGTACGTAC', '+', '??????????']
```

## 23) `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_sintax()`

- **Purpose:** Classify sequences taxonomically with the SINTAX algorithm.
- **Method:** Run `vsearch --sintax` against a FASTA or UDB database whose labels carry `;tax=d:...,p:...` annotations. `cutoff` is `--sintax_cutoff`, `random_seed` is `--randseed` and `random_ties` is `--sintax_random`.
- **API:** `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_sintax(input_path: 'str | Path', database: 'str | Path', output_directory: 'str | Path', *, cutoff: 'float | None' = None, random_seed: 'int | None' = None, random_ties: 'bool' = False, threads: 'int | None' = None, vsearch_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** Query FASTA or FASTQ, the taxonomy database and a new output directory.
- **Output:** `outputs` (`classification`), `records` (query, per-rank hits with confidence, strand, taxonomy above the cutoff) and `record_count`.
- **Example code:**

```python
import tempfile
from dnakit.fasta_fastq.vsearch_p09_cli_functions import vsearch_sintax

result = vsearch_sintax("toy.fa", "tax.fa", tempfile.mkdtemp() + "/run", cutoff=0.5, random_seed=1, verbose=False)
for row in result["records"][:2]:
    print(row["query"], row["taxonomy"])
```

- **Example output:**

```text
a;size=50 d:Bacteria,p:Firmicutes,c:Bacilli
b;size=30 d:Bacteria,p:Proteobacteria,c:Gamma
```

## 24) `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_make_udb()`

- **Purpose:** Build a UDB word-index database from a FASTA file.
- **Method:** Run `vsearch --makeudb_usearch`. `word_length` (`--wordlength`, default 8) and `database_mask` (`--dbmask`, default dust) are stored in the index; `hardmask` masks with N.
- **API:** `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_make_udb(input_fasta: 'str | Path', output_directory: 'str | Path', *, word_length: 'int | None' = None, database_mask: "Literal['none', 'dust', 'soft'] | None" = None, hardmask: 'bool' = False, vsearch_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** A FASTA file and a new output directory.
- **Output:** `outputs` (`database` = `database.udb`) and `database_bytes`.
- **Example code:**

```python
import tempfile
from dnakit.fasta_fastq.vsearch_p09_cli_functions import vsearch_make_udb

result = vsearch_make_udb("toy.fa", tempfile.mkdtemp() + "/run", word_length=6, verbose=False)
print(result["status"], result["database_bytes"] > 0)
```

- **Example output:**

```text
completed True
```

## 25) `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_udb_to_fasta()`

- **Purpose:** Export the sequences of a UDB file as FASTA.
- **Method:** Run `vsearch --udb2fasta`.
- **API:** `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_udb_to_fasta(input_udb: 'str | Path', output_directory: 'str | Path', *, vsearch_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** A UDB file and a new output directory.
- **Output:** `outputs` (`sequences`) and `record_count`.
- **Example code:**

```python
import tempfile
from dnakit.fasta_fastq.vsearch_p09_cli_functions import vsearch_make_udb, vsearch_udb_to_fasta

built = vsearch_make_udb("toy.fa", tempfile.mkdtemp() + "/udb", verbose=False)
result = vsearch_udb_to_fasta(built["outputs"]["database"], tempfile.mkdtemp() + "/run", verbose=False)
print(result["record_count"])
```

- **Example output:**

```text
5
```

## 26) `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_udb_report()`

- **Purpose:** Describe a UDB file: summary information or k-mer word statistics.
- **Method:** `report='info'` runs `vsearch --udbinfo` (number of sequences, alphabet, word width, index settings, parsed into `fields`); `report='stats'` runs `vsearch --udbstats` (word-size distribution and the most frequent words). The report is read through `--log`.
- **API:** `dnakit.fasta_fastq.vsearch_p09_cli_functions.vsearch_udb_report(input_udb: 'str | Path', output_directory: 'str | Path', *, report: "Literal['info', 'stats']" = 'info', vsearch_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Input:** A UDB file, a new output directory and `report`.
- **Output:** `outputs` (`report`), `report_text` and, for `info`, `fields`.
- **Example code:**

```python
import tempfile
from dnakit.fasta_fastq.vsearch_p09_cli_functions import vsearch_make_udb, vsearch_udb_report

built = vsearch_make_udb("toy.fa", tempfile.mkdtemp() + "/udb", word_length=6, verbose=False)
info = vsearch_udb_report(built["outputs"]["database"], tempfile.mkdtemp() + "/info", verbose=False)
print(info["fields"]["Seqs"], info["fields"]["Word width"])
stats = vsearch_udb_report(built["outputs"]["database"], tempfile.mkdtemp() + "/stats", report="stats", verbose=False)
print(stats["report_text"].splitlines()[0])
```

- **Example output:**

```text
5 6
450 nt in 5 seqs, min 90, max 90, avg 90
```
