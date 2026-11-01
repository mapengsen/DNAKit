# sourmash: lineage databases, taxonomy reports and signature utilities

This group of functions in `dnakit.microbiome` wraps [sourmash](https://github.com/sourmash-bio/sourmash), the FracMinHash toolkit for comparing and classifying sequence data. The functions work on local files that you supply: LCA lineage databases, signature files, SBT and SQLite indexes, `sourmash compare` matrices, and `sourmash gather` and lineage CSV files. They load and describe lineage databases, classify and summarise signatures on the lineage tree, compare lineage tables, write the report of an LCA index build, turn gather results into Krona, CSV, text, Kraken-style, CAMI and LIN-group reports, and handle signatures and SBT indexes (sketch from a list, grep, inflate, combine, migrate, convert). Nothing is downloaded, and no function builds a database from sequence reads.

sourmash is published under the BSD-3-Clause licence (checked in the installed package metadata). The dependency is optional and is imported only when a function is called. If it is missing, the function raises `BackendUnavailableError` with the install command. The examples import from `dnakit.microbiome`.

## Example data

The examples use a tiny synthetic data set instead of real genomes, so every number can be checked by hand: three genomes at k=31 and scaled=1000 whose hashes are small integers (hashes 6 to 10 are shared by the two *Escherichia* genomes). The examples that need a comparison matrix, an SBT index or a SQLite index create it in a temporary folder with the NumPy and sourmash APIs, and every output file is written to a temporary folder.

LCA database with three genomes: GCF_A (Escherichia coli, hashes 1 to 10), GCF_B (Escherichia fergusonii, hashes 6 to 15) and GCF_C (Bacillus subtilis, hashes 21 to 30).

**Example file:** `db.lca.json`

```text
{"version": "2.1", "type": "sourmash_lca", "license": "CC0", "ksize": 31, "scaled": 1000, "moltype": "DNA", "lid_to_lineage": {"0": [["superkingdom", "Bacteria"], ["phylum", "Pseudomonadota"], ["class", "Gammaproteobacteria"], ["order", "Enterobacterales"], ["family", "Enterobacteriaceae"], ["genus", "Escherichia"], ["species", "Escherichia coli"]], "1": [["superkingdom", "Bacteria"], ["phylum", "Pseudomonadota"], ["class", "Gammaproteobacteria"], ["order", "Enterobacterales"], ["family", "Enterobacteriaceae"], ["genus", "Escherichia"], ["species", "Escherichia fergusonii"]], "2": [["superkingdom", "Bacteria"], ["phylum", "Bacillota"], ["class", "Bacilli"], ["order", "Bacillales"], ["family", "Bacillaceae"], ["genus", "Bacillus"], ["species", "Bacillus subtilis"]]}, "hashval_to_idx": {"1": [0], "2": [0], "3": [0], "4": [0], "5": [0], "6": [0, 1], "7": [0, 1], "8": [0, 1], "9": [0, 1], "10": [0, 1], "11": [1], "12": [1], "13": [1], "14": [1], "15": [1], "21": [2], "22": [2], "23": [2], "24": [2], "25": [2], "26": [2], "27": [2], "28": [2], "29": [2], "30": [2]}, "ident_to_name": {"GCF_A": "GCF_A Escherichia coli", "GCF_B": "GCF_B Escherichia fergusonii", "GCF_C": "GCF_C Bacillus subtilis"}, "ident_to_idx": {"GCF_A": 0, "GCF_B": 1, "GCF_C": 2}, "idx_to_lid": {"0": 0, "1": 1, "2": 2}}
```

Query signature `query1` with the hashes 1 to 8, 22 and 23.

**Example file:** `query.sig`

```text
[{"class":"sourmash_signature","email":"","hash_function":"0.murmur64","filename":null,"name":"query1","license":"CC0","signatures":[{"num":0,"ksize":31,"seed":42,"max_hash":18446744073709552,"mins":[1,2,3,4,5,6,7,8,22,23],"md5sum":"80370ff3293372ea9e60f59626212d6a","molecule":"DNA"}],"version":0.4}]
```

Query signature `query2` with the hashes 11 to 15.

**Example file:** `query2.sig`

```text
[{"class":"sourmash_signature","email":"","hash_function":"0.murmur64","filename":null,"name":"query2","license":"CC0","signatures":[{"num":0,"ksize":31,"seed":42,"max_hash":18446744073709552,"mins":[11,12,13,14,15],"md5sum":"556d99975b863d98753e9796ed9cd6a5","molecule":"DNA"}],"version":0.4}]
```

The three genome signatures of the database as one signature file.

**Example file:** `genomes.sig`

```text
[{"class":"sourmash_signature","email":"","hash_function":"0.murmur64","filename":null,"name":"GCF_A Escherichia coli","license":"CC0","signatures":[{"num":0,"ksize":31,"seed":42,"max_hash":18446744073709552,"mins":[1,2,3,4,5,6,7,8,9,10],"md5sum":"41276f1cde0088f537f94eb1456e7dc7","molecule":"DNA"}],"version":0.4},{"class":"sourmash_signature","email":"","hash_function":"0.murmur64","filename":null,"name":"GCF_B Escherichia fergusonii","license":"CC0","signatures":[{"num":0,"ksize":31,"seed":42,"max_hash":18446744073709552,"mins":[6,7,8,9,10,11,12,13,14,15],"md5sum":"28cd92c421ea26a1c81c441fcbd134b4","molecule":"DNA"}],"version":0.4},{"class":"sourmash_signature","email":"","hash_function":"0.murmur64","filename":null,"name":"GCF_C Bacillus subtilis","license":"CC0","signatures":[{"num":0,"ksize":31,"seed":42,"max_hash":18446744073709552,"mins":[21,22,23,24,25,26,27,28,29,30],"md5sum":"44c0573131a0d93f2277857df1d1b155","molecule":"DNA"}],"version":0.4}]
```

A signature `reads_abund` that tracks abundance.

**Example file:** `abund.sig`

```text
[{"class":"sourmash_signature","email":"","hash_function":"0.murmur64","filename":null,"name":"reads_abund","license":"CC0","signatures":[{"num":0,"ksize":31,"seed":42,"max_hash":18446744073709552,"mins":[1,2,3,21,30],"md5sum":"e3fb529d70ad3c826c56d757251058fd","abundances":[5,3,1,2,7],"molecule":"DNA"}],"version":0.4}]
```

A gather CSV for `query1` with two matches.

**Example file:** `gather.csv`

```text
query_name,name,f_unique_weighted,f_unique_to_query,unique_intersect_bp,remaining_bp,query_md5,query_filename,query_bp,ksize,scaled,query_n_hashes,sum_weighted_found,total_weighted_hashes
query1,GCF_A Escherichia coli,0.8,0.8,8000,2000,80370ff3,query.sig,10000,31,1000,10,8,10
query1,GCF_C Bacillus subtilis,0.2,0.2,2000,0,80370ff3,query.sig,10000,31,1000,10,10,10
```

Lineage table for the three genomes.

**Example file:** `tax.csv`

```text
ident,superkingdom,phylum,class,order,family,genus,species
GCF_A,Bacteria,Pseudomonadota,Gammaproteobacteria,Enterobacterales,Enterobacteriaceae,Escherichia,Escherichia coli
GCF_B,Bacteria,Pseudomonadota,Gammaproteobacteria,Enterobacterales,Enterobacteriaceae,Escherichia,Escherichia fergusonii
GCF_C,Bacteria,Bacillota,Bacilli,Bacillales,Bacillaceae,Bacillus,Bacillus subtilis
```

The same lineages with an NCBI `taxpath` column.

**Example file:** `taxid.csv`

```text
ident,superkingdom,phylum,class,order,family,genus,species,taxpath
GCF_A,Bacteria,Pseudomonadota,Gammaproteobacteria,Enterobacterales,Enterobacteriaceae,Escherichia,Escherichia coli,2|1224|1236|91347|543|561|562
GCF_B,Bacteria,Pseudomonadota,Gammaproteobacteria,Enterobacterales,Enterobacteriaceae,Escherichia,Escherichia fergusonii,2|1224|1236|91347|543|561|564
GCF_C,Bacteria,Bacillota,Bacilli,Bacillales,Bacillaceae,Bacillus,Bacillus subtilis,2|1239|91061|1385|186817|1386|1423
```

LIN lineages for the three genomes.

**Example file:** `lin.csv`

```text
ident,lin
GCF_A,0;0;0
GCF_B,0;0;1
GCF_C,1;0;0
```

Two lingroups.

**Example file:** `lingroup.csv`

```text
name,lin
group_ecoli_like,0
group_bacillus_like,1
```

A classify-style table for the comparison.

**Example file:** `classify.csv`

```text
ID,status,superkingdom,phylum,class,order,family,genus,species
q1,found,Bacteria,Pseudomonadota,Gammaproteobacteria,Enterobacterales,Enterobacteriaceae,Escherichia,Escherichia coli
q2,found,Bacteria,Pseudomonadota,Gammaproteobacteria,Enterobacterales,Enterobacteriaceae,Escherichia,
q3,found,Bacteria,Bacillota,Bacilli,Bacillales,Bacillaceae,Bacillus,Bacillus subtilis
q4,found,Bacteria,Bacillota,Bacilli,Bacillales,Bacillaceae,Bacillus,Bacillus subtilis
```

A custom lineage table for the comparison.

**Example file:** `custom.csv`

```text
ID,superkingdom,phylum,class,order,family,genus,species
q1,Bacteria,Pseudomonadota,Gammaproteobacteria,Enterobacterales,Enterobacteriaceae,Escherichia,Escherichia coli
q2,Bacteria,Pseudomonadota,Gammaproteobacteria,Enterobacterales,Enterobacteriaceae,Escherichia,Escherichia fergusonii
q3,Bacteria,Pseudomonadota,Gammaproteobacteria,Enterobacterales,Enterobacteriaceae,Escherichia,Escherichia coli
q5,Archaea,Euryarchaeota,Methanobacteria,Methanobacteriales,Methanobacteriaceae,Methanobacterium,Methanobacterium formicicum
```

The list for `sketch_fromfile_signatures()`.

**Example file:** `genomes.csv`

```text
name,genome_filename,protein_filename
genome1,genome1.fa,
```

A 93-base genome.

**Example file:** `genome1.fa`

```text
>genome1
ATGGCGTACGTTAGCCTAGGATCCGATTACAGGCTTAAGCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCGGGTACCGAGCTCGAATTC
```

## 1) `dnakit.microbiome.check_and_load_gather_csvs_summary()`

- **Purpose:** Load several `sourmash gather` CSV files against lineage tables and report one summary row per query, failing early on empty files or missing identifiers.
- **Method:** Loads the lineage tables with sourmash's `MultiLineageDB.load()` and passes the gather files to `sourmash.tax.tax_utils.check_and_load_gather_csvs()`. That upstream function rejects empty gather files and, with `fail_on_missing_taxonomy=True`, match identifiers that have no lineage. DNAKit then reduces each loaded query to a summary row. `force=True` is passed to both the lineage loader and the gather loader: an empty or unreadable gather file is then skipped instead of raising `ValueError`, but a query that appears in two gather files still raises it.
- **API:** `check_and_load_gather_csvs_summary(gather_csvs: 'Sequence[str]', taxonomy_files: 'Sequence[str]', *, fail_on_missing_taxonomy: 'bool' = False, force: 'bool' = False, keep_full_identifiers: 'bool' = False, keep_identifier_versions: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package sourmash (extra `dnakit[bio-sourmash]`).
- **Input:** `gather_csvs` and `taxonomy_files` are non-empty lists of existing file paths (a bare string is refused). Lineage tables are CSV or sqlite files in sourmash lineage format (an `ident` column plus one column per rank). The gather CSV needs the sourmash 4.4+ columns `query_name`, `name`, `f_unique_weighted`, `f_unique_to_query`, `unique_intersect_bp`, `remaining_bp`, `query_md5`, `query_filename`, `query_bp`, `ksize` and `scaled`. `keep_full_identifiers` and `keep_identifier_versions` control how the match name is turned into a lineage identifier.
- **Output:** A dictionary with `status`, `backend`, the resolved `gather_csvs` paths, and `queries`: one row per query with `query_name`, `query_md5`, `query_filename`, `ksize`, `scaled`, `match_count`, `missed_identifiers` and `skipped_identifiers`.
- **Example code:**

```python
from dnakit.microbiome import check_and_load_gather_csvs_summary

result = check_and_load_gather_csvs_summary(["gather.csv"], ["tax.csv"], verbose=False)
print(result["status"], len(result["gather_csvs"]))
for query in result["queries"]:
    print(query["query_name"], query["match_count"], query["missed_identifiers"], query["scaled"])
```

- **Example output:**

```text
completed 1
query1 2 0 1000
```

## 2) `dnakit.microbiome.classify_signature_lineage()`

- **Purpose:** Assign every signature in one file to the lowest-common-ancestor lineage of its k-mers, like `sourmash lca classify`.
- **Method:** Loads the LCA databases with `sourmash.lca.lca_utils.load_databases()`, loads the signatures of the file at the databases' k-mer size, downsamples a finer signature to the database `scaled` value, and calls `sourmash.lca.command_classify.classify_signature()` for each one. Lineages supported by fewer than `threshold` k-mers are ignored. With `majority=True` the most frequent lineage is used instead of the lowest common ancestor. A query that is coarser than the database raises `ConfigurationError`.
- **API:** `classify_signature_lineage(signature: 'str', lca_databases: 'Sequence[str]', *, threshold: 'int' = 5, majority: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package sourmash (extra `dnakit[bio-sourmash]`).
- **Input:** `signature` is an existing signature file (it may hold several signatures). `lca_databases` is a non-empty list of existing LCA database files that share one k-mer size. `threshold` is a non-negative integer (sourmash default 5). `majority` is a boolean.
- **Output:** A dictionary with `status`, `backend`, `ksize`, `scaled`, `threshold`, `majority` and `classifications`: one row per signature with `name`, `status` (`found`, `disagree` or `nomatch`) and `lineage`, a list of `{"rank", "name"}` mappings from the most general rank to the most specific.
- **Example code:**

```python
from dnakit.microbiome import classify_signature_lineage

for threshold in (5, 2):
    result = classify_signature_lineage("query.sig", ["db.lca.json"], threshold=threshold, verbose=False)
    row = result["classifications"][0]
    print(threshold, row["name"], row["status"], row["lineage"][-1]["rank"], row["lineage"][-1]["name"])
```

- **Example output:**

```text
5 query1 found species Escherichia coli
2 query1 disagree superkingdom Bacteria
```

## 3) `dnakit.microbiome.gather_lineage_assignments()`

- **Purpose:** List, for every k-mer hash of one signature, the lineages that the LCA databases assign it to.
- **Method:** Loads the databases with `load_databases()`, requires exactly one signature at the databases' k-mer size in the file, downsamples it to the database `scaled` value, and calls `sourmash.lca.lca_utils.gather_assignments()` on its hashes. Hash values are turned into strings and the lineages of each hash are sorted, so the result is stable.
- **API:** `gather_lineage_assignments(signature: 'str', lca_databases: 'Sequence[str]', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package sourmash (extra `dnakit[bio-sourmash]`).
- **Input:** `signature` is an existing signature file that holds exactly one signature with the databases' k-mer size (otherwise `ConfigurationError`). `lca_databases` is a non-empty list of existing LCA database files.
- **Output:** A dictionary with `status`, `backend`, `signature`, `ksize`, `scaled`, `assigned_hash_count` and `assignments`: one row per hash found in a database, with `hashval` (a string) and `lineages` (a list of lineages). Hashes that no database contains are not listed.
- **Example code:**

```python
from dnakit.microbiome import gather_lineage_assignments

result = gather_lineage_assignments("query.sig", ["db.lca.json"], verbose=False)
print(result["ksize"], result["scaled"], result["assigned_hash_count"])
for row in result["assignments"]:
    species = [lineage[-1]["name"] for lineage in row["lineages"]]
    print(row["hashval"], species)
```

- **Example output:**

```text
31 1000 10
1 ['Escherichia coli']
2 ['Escherichia coli']
3 ['Escherichia coli']
4 ['Escherichia coli']
5 ['Escherichia coli']
6 ['Escherichia coli', 'Escherichia fergusonii']
7 ['Escherichia coli', 'Escherichia fergusonii']
8 ['Escherichia coli', 'Escherichia fergusonii']
22 ['Bacillus subtilis']
23 ['Bacillus subtilis']
```

## 4) `dnakit.microbiome.lca_rank_counts()`

- **Purpose:** Count, for every lineage, the k-mer hashes of the databases whose lowest common ancestor is that lineage.
- **Method:** Loads the databases and calls `sourmash.lca.command_rankinfo.make_lca_counts()`, which walks all hashes of the databases and records the lowest common ancestor of each. DNAKit serialises the lineages and sorts the rows by count (largest first) and then by lineage.
- **API:** `lca_rank_counts(lca_databases: 'Sequence[str]', *, min_num: 'int' = 0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package sourmash (extra `dnakit[bio-sourmash]`).
- **Input:** `lca_databases` is a non-empty list of existing LCA database files. `min_num` is a non-negative integer: lineages with fewer hashes are dropped (default 0). The databases are used at their own `scaled` value.
- **Output:** A dictionary with `status`, `backend`, `ksize`, `scaled`, `min_num` and `counts`: rows with `lineage` (list of `{"rank", "name"}`) and `count`.
- **Example code:**

```python
from dnakit.microbiome import lca_rank_counts

result = lca_rank_counts(["db.lca.json"], verbose=False)
for row in result["counts"]:
    print(row["lineage"][-1]["rank"], row["lineage"][-1]["name"], row["count"])
```

- **Example output:**

```text
species Bacillus subtilis 10
genus Escherichia 5
species Escherichia coli 5
species Escherichia fergusonii 5
```

## 5) `dnakit.microbiome.load_comparison_matrix()`

- **Purpose:** Read the similarity matrix and the labels written by `sourmash compare -o BASE`.
- **Method:** Checks that `basefile` and `<basefile>.labels.txt` exist and calls `sourmash.fig.load_matrix_and_labels()`. The matrix is converted to floats and returned as nested lists.
- **API:** `load_comparison_matrix(basefile: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages sourmash (extra `dnakit[bio-sourmash]`) and numpy (extra `dnakit[bio-analysis]`).
- **Input:** `basefile` is the path given to `sourmash compare -o`: a NumPy `.npy` array stored under that exact name (no extension is added), with the label file next to it as `<basefile>.labels.txt`.
- **Output:** A dictionary with `status`, `backend`, `basefile`, `size` (the matrix side length), `labels` and `matrix` (a list of rows).
- **Example code:**

```python
import tempfile
from pathlib import Path

import numpy

from dnakit.microbiome import load_comparison_matrix

work = Path(tempfile.mkdtemp())
base = work / "cmp"
with open(base, "wb") as handle:
    numpy.save(handle, numpy.array([[1.0, 0.5, 0.0], [0.5, 1.0, 0.25], [0.0, 0.25, 1.0]]))
(work / "cmp.labels.txt").write_text("genome_A\ngenome_B\ngenome_C\n")
result = load_comparison_matrix(str(base), verbose=False)
print(result["size"], result["labels"])
print(result["matrix"][1])
```

- **Example output:**

```text
3 ['genome_A', 'genome_B', 'genome_C']
[0.5, 1.0, 0.25]
```

## 6) `dnakit.microbiome.load_gather_results_summary()`

- **Purpose:** Read one `sourmash gather` CSV, attach lineages from lineage tables and summarise each query.
- **Method:** Loads the lineage tables with `MultiLineageDB.load()` and calls `sourmash.tax.tax_utils.load_gather_results()`. Each query is reduced to a summary row. With `fail_on_missing_taxonomy=True`, a match without a lineage is an error; otherwise it is counted in `missed_identifiers`.
- **API:** `load_gather_results_summary(gather_csv: 'str', taxonomy_files: 'Sequence[str]', *, fail_on_missing_taxonomy: 'bool' = False, keep_full_identifiers: 'bool' = False, keep_identifier_versions: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package sourmash (extra `dnakit[bio-sourmash]`).
- **Input:** `gather_csv` is an existing gather CSV with the sourmash 4.4+ columns (see function 1). `taxonomy_files` is a non-empty list of CSV or sqlite lineage tables. `keep_full_identifiers` and `keep_identifier_versions` control how a match name is reduced to an identifier.
- **Output:** A dictionary with `status`, `backend`, `gather_csv`, `columns` (the header of the gather CSV) and `queries` with the same row fields as function 1.
- **Example code:**

```python
from dnakit.microbiome import load_gather_results_summary

result = load_gather_results_summary("gather.csv", ["tax.csv"], verbose=False)
print(len(result["columns"]), result["columns"][:3])
for query in result["queries"]:
    print(query["query_name"], query["query_md5"], query["match_count"], query["ksize"])
```

- **Example output:**

```text
14 ['query_name', 'name', 'f_unique_weighted']
query1 80370ff3 2 31
```

## 7) `dnakit.microbiome.load_lca_database()`

- **Purpose:** Load one LCA database and report its k-mer size, scaled value and size.
- **Method:** Calls `sourmash.lca.lca_db.load_single_database()` on the file and counts the signatures and the hashes it holds.
- **API:** `load_lca_database(path: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package sourmash (extra `dnakit[bio-sourmash]`).
- **Input:** `path` is an existing LCA database file (JSON, or a sqlite LCA database) written by `sourmash lca index`.
- **Output:** A dictionary with `status`, `backend`, `path`, `ksize`, `scaled`, `signature_count` and `hash_count`.
- **Example code:**

```python
from dnakit.microbiome import load_lca_database

result = load_lca_database("db.lca.json", verbose=False)
print(result["ksize"], result["scaled"], result["signature_count"], result["hash_count"])
```

- **Example output:**

```text
31 1000 3 25
```

## 8) `dnakit.microbiome.load_lca_databases()`

- **Purpose:** Load several LCA databases at once and describe each one.
- **Method:** Calls `sourmash.lca.lca_utils.load_databases()`, which requires that all databases share one k-mer size, and reports the shared k-mer size and `scaled` value together with a summary of each database. With `scaled`, the databases are loaded at that coarser resolution.
- **API:** `load_lca_databases(lca_databases: 'Sequence[str]', *, scaled: 'int | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package sourmash (extra `dnakit[bio-sourmash]`).
- **Input:** `lca_databases` is a non-empty list of existing LCA database files. `scaled` is an optional positive integer. A value larger than the database's own downsamples the databases; a smaller or equal value leaves them at their own `scaled`.
- **Output:** A dictionary with `status`, `backend`, `ksize`, `scaled` and `databases`: one row per file with `path`, `ksize`, `scaled`, `signature_count` and `hash_count`.
- **Example code:**

```python
from dnakit.microbiome import load_lca_databases

result = load_lca_databases(["db.lca.json"], verbose=False)
print(result["ksize"], result["scaled"], result["databases"][0]["hash_count"])
coarse = load_lca_databases(["db.lca.json"], scaled=2000, verbose=False)
print(coarse["scaled"], coarse["databases"][0]["scaled"])
```

- **Example output:**

```text
31 1000 25
2000 2000
```

## 9) `dnakit.microbiome.load_sbt_index_summary()`

- **Purpose:** Open a sourmash SBT index (`.sbt.zip`) and list its signatures.
- **Method:** Calls `sourmash.sbtmh.load_sbt_index()` (without the version warning) and iterates over the signatures of the index to collect names, k-mer sizes and `scaled` values.
- **API:** `load_sbt_index_summary(path: 'str', *, cache_size: 'int | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package sourmash (extra `dnakit[bio-sourmash]`).
- **Input:** `path` is an existing SBT index file. `cache_size` is an optional cache size passed to sourmash.
- **Output:** A dictionary with `status`, `backend`, `path`, `index_type` (always `SBT`), `signature_count`, `signature_names`, `ksizes` and `scaled_values` (both sorted lists).
- **Example code:**

```python
import tempfile
from pathlib import Path

import sourmash
from sourmash.sbtmh import create_sbt_index

from dnakit.microbiome import load_sbt_index_summary

path = Path(tempfile.mkdtemp()) / "genomes.sbt.zip"
tree = create_sbt_index(1000)
for signature in sourmash.load_file_as_signatures("genomes.sig"):
    tree.insert(signature)
tree.save(str(path))
result = load_sbt_index_summary(str(path), verbose=False)
print(result["index_type"], result["signature_count"], result["ksizes"], result["scaled_values"])
print(sorted(result["signature_names"]))
```

- **Example output:**

```text
SBT 3 [31] [1000]
['GCF_A Escherichia coli', 'GCF_B Escherichia fergusonii', 'GCF_C Bacillus subtilis']
```

## 10) `dnakit.microbiome.load_sqlite_index_summary()`

- **Purpose:** Open a sourmash SQLite signature collection and list its signatures.
- **Method:** Calls `sourmash.index.sqlite_index.load_sqlite_index()` and iterates over the signatures to collect names, k-mer sizes and `scaled` values. `index_type` is the class name of the object that sourmash returns.
- **API:** `load_sqlite_index_summary(path: 'str', *, request_manifest: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package sourmash (extra `dnakit[bio-sourmash]`).
- **Input:** `path` is an existing SQLite signature collection. `request_manifest` is passed to sourmash.
- **Output:** A dictionary with `status`, `backend`, `path`, `index_type`, `signature_count`, `signature_names`, `ksizes` and `scaled_values`.
- **Example code:**

```python
import tempfile
from pathlib import Path

import sourmash
from sourmash.index.sqlite_index import SqliteIndex

from dnakit.microbiome import load_sqlite_index_summary

path = Path(tempfile.mkdtemp()) / "genomes.sqldb"
index = SqliteIndex.create(str(path))
for signature in sourmash.load_file_as_signatures("genomes.sig"):
    index.insert(signature)
index.commit()
index.close()
result = load_sqlite_index_summary(str(path), verbose=False)
print(result["index_type"], result["signature_count"], result["ksizes"], result["scaled_values"])
print(sorted(result["signature_names"]))
```

- **Example output:**

```text
SqliteIndex 3 [31] [1000]
['GCF_A Escherichia coli', 'GCF_B Escherichia fergusonii', 'GCF_C Bacillus subtilis']
```

## 11) `dnakit.microbiome.plot_comparison_matrix_composite()`

- **Purpose:** Draw a dendrogram with the reordered heatmap for a `sourmash compare` matrix and save it as an image.
- **Method:** Loads the matrix with `sourmash.fig.load_matrix_and_labels()`, draws it with `sourmash.fig.plot_composite_matrix()` (single-linkage clustering, YlGnBu colours), saves the figure at 150 dpi and closes it. The leaf order of the dendrogram and the reordered matrix are returned. A matrix outside [0, 1] is an error unless `force=True`, which rescales it.
- **API:** `plot_comparison_matrix_composite(basefile: 'str', output_image: 'str', *, show_labels: 'bool' = True, vmin: 'float' = 0.0, vmax: 'float' = 1.0, force: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages sourmash (extra `dnakit[bio-sourmash]`), matplotlib (extra `dnakit[bio-plots]`) and scipy (extra `dnakit[bio-analysis]`).
- **Input:** `basefile` as in function 5. `output_image` is a new file name ending in `.png`, `.pdf` or `.svg` in an existing folder. `show_labels` toggles the labels. `vmin` must be smaller than `vmax` (defaults 0.0 and 1.0).
- **Output:** A dictionary with `status`, `backend`, `image`, `labels` (in dendrogram order), `matrix` (reordered rows), `vmin` and `vmax`.
- **Example code:**

```python
import tempfile
from pathlib import Path

import numpy

from dnakit.microbiome import plot_comparison_matrix_composite

work = Path(tempfile.mkdtemp())
base = work / "cmp"
with open(base, "wb") as handle:
    numpy.save(handle, numpy.array([[1.0, 0.1, 0.8], [0.1, 1.0, 0.2], [0.8, 0.2, 1.0]]))
(work / "cmp.labels.txt").write_text("genome_A\ngenome_B\ngenome_C\n")
image = work / "cmp.png"
result = plot_comparison_matrix_composite(str(base), str(image), verbose=False)
print(result["labels"])
print([[round(value, 2) for value in row] for row in result["matrix"]])
print(image.read_bytes()[:4] == b"\x89PNG", image.stat().st_size > 1000)
```

- **Example output:**

```text
['genome_B', 'genome_A', 'genome_C']
[[1.0, 0.1, 0.2], [0.1, 1.0, 0.8], [0.2, 0.8, 1.0]]
True True
```

## 12) `dnakit.microbiome.summarize_lineage_counts()`

- **Purpose:** Count the k-mers of one signature on the lineage tree: each lineage and all its ancestors receive the k-mers assigned below them.
- **Method:** Loads the databases, requires exactly one signature at their k-mer size, downsamples it to the database `scaled` value, and calls `sourmash.lca.command_summarize.summarize()`. Lineages with fewer than `threshold` k-mers are not counted. The rows are sorted by count (largest first), then by lineage.
- **API:** `summarize_lineage_counts(signature: 'str', lca_databases: 'Sequence[str]', *, threshold: 'int' = 5, ignore_abundance: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package sourmash (extra `dnakit[bio-sourmash]`).
- **Input:** `signature` is an existing file with exactly one signature at the databases' k-mer size. `lca_databases` is a non-empty list of LCA database files. `threshold` is a non-negative integer (sourmash default 5). With `ignore_abundance=True` every hash is counted once even if the signature tracks abundance.
- **Output:** A dictionary with `status`, `backend`, `signature`, `ksize`, `scaled`, `threshold`, `ignore_abundance` and `counts`: rows with `lineage` and `count`.
- **Example code:**

```python
from dnakit.microbiome import summarize_lineage_counts

result = summarize_lineage_counts("query.sig", ["db.lca.json"], threshold=2, verbose=False)
for row in result["counts"]:
    last = row["lineage"][-1]
    print(last["rank"], last["name"], row["count"])
```

- **Example output:**

```text
superkingdom Bacteria 10
phylum Pseudomonadota 8
class Gammaproteobacteria 8
order Enterobacterales 8
family Enterobacteriaceae 8
genus Escherichia 8
species Escherichia coli 5
phylum Bacillota 2
class Bacilli 2
order Bacillales 2
family Bacillaceae 2
genus Bacillus 2
species Bacillus subtilis 2
```

## 13) `dnakit.microbiome.classify_lca_signatures()`

- **Purpose:** Classify every signature of several files against LCA databases and write the result as a CSV file.
- **Method:** Runs the computation of `sourmash lca classify`: loads the databases (optionally at a coarser `scaled`), reads each query file at the databases' k-mer size, downsamples signatures when needed, and calls `classify_signature()` for each. The CSV is written row by row with the columns `ID`, `status` and one column per taxonomic rank.
- **API:** `classify_lca_signatures(query_files: 'Sequence[str]', lca_databases: 'Sequence[str]', output_csv: 'str', *, threshold: 'int' = 5, majority: 'bool' = False, scaled: 'int | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package sourmash (extra `dnakit[bio-sourmash]`).
- **Input:** `query_files` and `lca_databases` are non-empty lists of existing files. `output_csv` is a new file in an existing folder. `threshold` is a non-negative integer (default 5). `majority` uses the most frequent lineage. `scaled` is an optional positive integer for loading the databases at a coarser resolution.
- **Output:** A dictionary with `status`, `backend`, `output_csv`, `ksize`, `scaled`, `threshold`, `majority`, `classification_count` and `classifications`: the same rows as the CSV, each with `id`, `status` and `lineage`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.microbiome import classify_lca_signatures

output = Path(tempfile.mkdtemp()) / "classified.csv"
result = classify_lca_signatures(
    ["query.sig", "query2.sig"], ["db.lca.json"], str(output), verbose=False
)
print(result["classification_count"], result["ksize"], result["scaled"])
print(output.read_text().splitlines()[0])
for row in result["classifications"]:
    print(row["id"], row["status"], row["lineage"][-1]["name"])
```

- **Example output:**

```text
2 31 1000
ID,status,superkingdom,phylum,class,order,family,genus,species,strain
query1 found Escherichia coli
query2 found Escherichia fergusonii
```

## 14) `dnakit.microbiome.compare_lineage_assignment_csvs()`

- **Purpose:** Compare the lineage that two tables assign to each identifier, like `sourmash lca compare-csv`.
- **Method:** Reads the first table with `sourmash.lca.command_index.load_taxonomy_assignments()` starting at column 3 (a classify-style CSV) and the second with the given delimiter and start column. For every identifier present in both tables with different lineages, it builds a lineage tree and calls `find_lca()`: the verdict is `compatible` when one lineage is an ancestor of the other and `incompatible` otherwise. Incompatible lineages are counted at the rank of the common ancestor (or at the first rank when there is none).
- **API:** `compare_lineage_assignment_csvs(classify_csv: 'str', custom_csv: 'str', *, start_column: 'int' = 2, tabs: 'bool' = False, no_headers: 'bool' = False, force: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package sourmash (extra `dnakit[bio-sourmash]`).
- **Input:** `classify_csv` and `custom_csv` are existing files. `start_column` is the 1-based column where the lineage starts in the second table (at least 2, default 2). `tabs` selects a tab-separated second table. `no_headers` says the second table has no header row. `force` is passed to sourmash and keeps going over problem rows.
- **Output:** A dictionary with `status`, `backend`, the two paths, `classify_rows`, `custom_rows`, `classify_lineage_count`, `custom_lineage_count`, `missing_in_custom`, `missing_in_classify`, `common_count`, `different_count`, `compatible_count`, `incompatible_count`, `incompatible_by_rank` (one entry per rank) and `differences` (rows with `identifier`, `verdict` and `lca`).
- **Example code:**

```python
from dnakit.microbiome import compare_lineage_assignment_csvs

result = compare_lineage_assignment_csvs("classify.csv", "custom.csv", verbose=False)
print(result["classify_rows"], result["custom_rows"], result["common_count"])
print(result["missing_in_custom"], result["missing_in_classify"])
print(result["compatible_count"], result["incompatible_count"])
for row in result["differences"]:
    print(row["identifier"], row["verdict"], [pair["name"] for pair in row["lca"]][-1])
print([item for item in result["incompatible_by_rank"] if item["count"]])
```

- **Example output:**

```text
4 4 3
['q4'] ['q5']
1 1
q2 compatible Escherichia fergusonii
q3 incompatible Bacteria
[{'rank': 'superkingdom', 'count': 1}]
```

## 15) `dnakit.microbiome.lca_rank_info()`

- **Purpose:** Count the database k-mer hashes whose lowest common ancestor falls on each taxonomic rank, like `sourmash lca rankinfo`.
- **Method:** Loads the databases (optionally at a coarser `scaled`), collects the per-lineage counts with `make_lca_counts()`, adds them up by the rank of the last entry of each lineage, and reports a count and a percentage for every rank in taxonomic order. The percentages are not rounded.
- **API:** `lca_rank_info(lca_databases: 'Sequence[str]', *, min_num: 'int' = 0, scaled: 'int | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package sourmash (extra `dnakit[bio-sourmash]`).
- **Input:** `lca_databases` is a non-empty list of existing LCA database files. `min_num` is a non-negative integer (default 0). `scaled` is an optional positive integer.
- **Output:** A dictionary with `status`, `backend`, `ksize`, `scaled`, `min_num`, `total`, `note` (`None`, or a message when no hash has a lineage) and `ranks`: one row per rank with `rank`, `count` and `percent`. `ranks` is empty when the total is 0.
- **Example code:**

```python
from dnakit.microbiome import lca_rank_info

result = lca_rank_info(["db.lca.json"], verbose=False)
print(result["total"], result["note"])
for row in result["ranks"]:
    if row["count"]:
        print(row["rank"], row["count"], round(row["percent"], 1))
```

- **Example output:**

```text
25 None
genus 5 20.0
species 20 80.0
```

## 16) `dnakit.microbiome.summarize_lca_signatures()`

- **Purpose:** Summarise the k-mers of every signature of several files on the lineage tree, optionally as a CSV file, like `sourmash lca summarize`.
- **Method:** Loads the databases, reads every signature of the query files one by one with `load_singletons_and_count()`, assigns the hashes with `summarize()` and gives each lineage and all its ancestors the number of hashes. When `output_csv` is given, the rows are written with sourmash's own `output_csv()` (columns `count`, one per rank, `filename`, `sig_name`, `sig_md5`, `total_counts`).
- **API:** `summarize_lca_signatures(query_files: 'Sequence[str]', lca_databases: 'Sequence[str]', *, output_csv: 'str | None' = None, threshold: 'int' = 5, ignore_abundance: 'bool' = False, scaled: 'int | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package sourmash (extra `dnakit[bio-sourmash]`).
- **Input:** `query_files` and `lca_databases` are non-empty lists of existing files. `output_csv` is an optional new file in an existing folder. `threshold` is a non-negative integer (default 5). `ignore_abundance=True` counts each hash once. `scaled` is an optional positive integer for the databases.
- **Output:** A dictionary with `status`, `backend`, `output_csv` (a path or `None`), `ksize`, `scaled`, `threshold`, `ignore_abundance`, `signature_count` and `signatures`: one row per signature with `filename`, `name`, `md5`, `total` (the sum of abundances, or the number of hashes with `ignore_abundance`) and `counts` (rows with `lineage` and `count`).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.microbiome import summarize_lca_signatures

output = Path(tempfile.mkdtemp()) / "summary.csv"
result = summarize_lca_signatures(
    ["query.sig", "query2.sig"], ["db.lca.json"], output_csv=str(output), threshold=2, verbose=False
)
print(result["signature_count"], output.read_text().splitlines()[0])
for signature in result["signatures"]:
    top = max(signature["counts"], key=lambda row: (len(row["lineage"]), row["count"]))
    print(signature["name"], signature["total"], len(signature["counts"]), top["lineage"][-1]["name"], top["count"])
```

- **Example output:**

```text
2 count,superkingdom,phylum,class,order,family,genus,species,strain,filename,sig_name,sig_md5,total_counts
query1 10.0 13 Escherichia coli 5
query2 5.0 7 Escherichia fergusonii 5
```

## 17) `dnakit.microbiome.write_lca_index_report()`

- **Purpose:** Write the anomaly report of an LCA index build as a text file.
- **Method:** Converts the unused lineages to sourmash `LineagePair` tuples and calls `sourmash.lca.command_index.generate_report()`, the function behind the `--report` file of `sourmash lca index`. It does not build an index and reads no signature.
- **API:** `write_lca_index_report(duplicate_signatures: 'Sequence[str]', identifiers_without_lineage: 'Sequence[str]', identifiers_without_signature: 'Sequence[str]', unused_identifiers: 'Sequence[str]', unused_lineages: 'Sequence[Sequence[Mapping[str, str]]]', output_report: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package sourmash (extra `dnakit[bio-sourmash]`).
- **Input:** The first four arguments are lists of identifiers (duplicate signatures, identifiers without a lineage, identifiers without a signature, unused identifiers). `unused_lineages` is a list of lineages, each a list of `{"rank", "name"}` mappings from the most general rank to the most specific. `output_report` is a new file in an existing folder.
- **Output:** A dictionary with `status`, `backend`, `report` (the path) and the counts `duplicate_count`, `no_lineage_count`, `no_signature_count`, `unused_identifier_count` and `unused_lineage_count`. The report file lists each group with its number and entries; unused lineages are written as `;`-separated names.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.microbiome import write_lca_index_report

lineage = [{"rank": "superkingdom", "name": "Bacteria"}, {"rank": "phylum", "name": "Bacillota"}]
report = Path(tempfile.mkdtemp()) / "report.txt"
result = write_lca_index_report(
    ["GCF_A"], ["GCF_X"], ["GCF_Y"], ["GCF_Z"], [lineage], str(report), verbose=False
)
print({key: value for key, value in result.items() if key.endswith("_count")})
print(report.read_text().rstrip())
```

- **Example output:**

```text
{'duplicate_count': 1, 'no_lineage_count': 1, 'no_signature_count': 1, 'unused_identifier_count': 1, 'unused_lineage_count': 1}
Duplicate signatures: 1
GCF_A
----
Unused identifiers: 1
GCF_Z
----
No lineage provided for these identifiers: 1
GCF_X
----
No signatures found for these identifiers: 1
GCF_Y
----
Unused lineages: 1
Bacteria;Bacillota;;;;;;
```

## 18) `dnakit.microbiome.sketch_fromfile_signatures()`

- **Purpose:** Sketch every genome listed in CSV files and save the signatures into one new collection.
- **Method:** Runs `sourmash sketch fromfile CSV... --output-signatures OUT --license LICENSE` with one `-p` per parameter string, plus `--already-done` and `--check-sequence` when requested. Relative genome paths in the CSV are resolved from the current working directory because the command runs there. The format of the collection follows the suffix of `output_signatures` (for example `.zip`).
- **API:** `sketch_fromfile_signatures(csv_paths: 'Sequence[str]', output_signatures: 'str', *, param_strings: 'Sequence[str]' = ('dna,k=31,scaled=1000',), already_done: 'Sequence[str]' = (), license: 'str' = 'CC0', check_sequence: 'bool' = False, sourmash_path: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package sourmash (extra `dnakit[bio-sourmash]`). The function runs the `sourmash` command, which is installed with the package; it is looked up on PATH, then next to the Python interpreter, or taken from `sourmash_path`. If it is not found, `BackendUnavailableError` is raised.
- **Input:** `csv_paths` is a non-empty list of CSV files with the columns `name`, `genome_filename` and `protein_filename` (a protein path may be empty). `output_signatures` is a new file in an existing folder. `param_strings` is a non-empty list of sourmash parameter strings such as `dna,k=31,scaled=1000` (the default). `already_done` lists existing collections whose signatures are not recalculated. `license` is the signature licence (default `CC0`). `check_sequence` makes sourmash check the sequence.
- **Output:** A dictionary with `status` (`completed` when the command exits with code 0 and every expected output exists, otherwise `failed`), `backend`, `command` (the arguments, with the program shown as `sourmash`), `return_code`, `outputs` (the expected output paths that exist) and `log_tail` (the last 4000 characters of the command output).
- **Example code:**

```python
import tempfile
from pathlib import Path

import sourmash

from dnakit.microbiome import sketch_fromfile_signatures

output = Path(tempfile.mkdtemp()) / "sketched.zip"
result = sketch_fromfile_signatures(
    ["genomes.csv"], str(output), param_strings=["dna,k=21,scaled=1"], verbose=False
)
print(result["status"], result["return_code"], len(result["outputs"]))
for signature in sourmash.load_file_as_signatures(str(output)):
    print(signature.name, signature.minhash.ksize, signature.minhash.scaled, len(signature.minhash))
```

- **Example output:**

```text
completed 0 1
genome1 21 1 73
```

## 19) `dnakit.microbiome.grep_signatures()`

- **Purpose:** Keep the signatures whose name matches a regular expression and save them to a new file.
- **Method:** Runs `sourmash sig grep PATTERN SIGNATURE... -o OUT`, adding `-i`, `-v`, `-k`, the molecule-type flag and `--csv` for the options that are set.
- **API:** `grep_signatures(pattern: 'str', signatures: 'Sequence[str]', output: 'str', *, ignore_case: 'bool' = False, invert_match: 'bool' = False, ksize: 'int | None' = None, moltype: 'str | None' = None, csv_manifest: 'str | None' = None, sourmash_path: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package sourmash (extra `dnakit[bio-sourmash]`). The function runs the `sourmash` command, which is installed with the package; it is looked up on PATH, then next to the Python interpreter, or taken from `sourmash_path`. If it is not found, `BackendUnavailableError` is raised.
- **Input:** `pattern` is a non-empty regular expression. `signatures` is a non-empty list of existing signature files. `output` is a new file in an existing folder. `ignore_case` matches both cases. `invert_match` keeps the signatures that do not match. `ksize` is a positive integer. `moltype` is one of `DNA`, `RNA`, `protein`, `dayhoff` or `hp`. `csv_manifest` is a new file that receives a CSV manifest of the kept signatures.
- **Output:** A dictionary with `status` (`completed` when the command exits with code 0 and every expected output exists, otherwise `failed`), `backend`, `command` (the arguments, with the program shown as `sourmash`), `return_code`, `outputs` (the expected output paths that exist) and `log_tail` (the last 4000 characters of the command output).
- **Example code:**

```python
import tempfile
from pathlib import Path

import sourmash

from dnakit.microbiome import grep_signatures

work = Path(tempfile.mkdtemp())
result = grep_signatures("escherichia", ["genomes.sig"], str(work / "kept.sig"), ignore_case=True, verbose=False)
print(result["status"], [Path(p).name for p in result["outputs"]])
print([s.name for s in sourmash.load_file_as_signatures(str(work / "kept.sig"))])
result = grep_signatures("escherichia", ["genomes.sig"], str(work / "rest.sig"), ignore_case=True, invert_match=True, verbose=False)
print([s.name for s in sourmash.load_file_as_signatures(str(work / "rest.sig"))])
```

- **Example output:**

```text
completed ['kept.sig']
['GCF_A Escherichia coli', 'GCF_B Escherichia fergusonii']
['GCF_C Bacillus subtilis']
```

## 20) `dnakit.microbiome.inflate_signatures()`

- **Purpose:** Give other signatures the abundances of a reference signature that tracks abundance.
- **Method:** Runs `sourmash sig inflate REFERENCE OTHER... -o OUT` with `-k KSIZE`, or with `--dna` when no `ksize` is given. Upstream, each other signature is inflated with `MinHash.inflate()` against the reference, so only the hashes present in both signatures keep an abundance and the output signatures track abundance.
- **API:** `inflate_signatures(signature_from: 'str', other_signatures: 'Sequence[str]', output: 'str', *, ksize: 'int | None' = None, sourmash_path: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package sourmash (extra `dnakit[bio-sourmash]`). The function runs the `sourmash` command, which is installed with the package; it is looked up on PATH, then next to the Python interpreter, or taken from `sourmash_path`. If it is not found, `BackendUnavailableError` is raised.
- **Input:** `signature_from` is an existing signature file with abundances. `other_signatures` is a non-empty list of existing signature files. `output` is a new file in an existing folder. `ksize` is a positive integer.
- **Output:** A dictionary with `status` (`completed` when the command exits with code 0 and every expected output exists, otherwise `failed`), `backend`, `command` (the arguments, with the program shown as `sourmash`), `return_code`, `outputs` (the expected output paths that exist) and `log_tail` (the last 4000 characters of the command output).
- **Example code:**

```python
import tempfile
from pathlib import Path

import sourmash

from dnakit.microbiome import inflate_signatures

output = Path(tempfile.mkdtemp()) / "inflated.sig"
result = inflate_signatures("abund.sig", ["genomes.sig"], str(output), ksize=31, verbose=False)
print(result["status"], result["return_code"])
for signature in sourmash.load_file_as_signatures(str(output)):
    print(signature.name, signature.minhash.track_abundance, sorted(signature.minhash.hashes.items()))
```

- **Example output:**

```text
completed 0
GCF_A Escherichia coli True [(1, 5), (2, 3), (3, 1)]
GCF_B Escherichia fergusonii True []
GCF_C Bacillus subtilis True [(21, 2), (30, 7)]
```

## 21) `dnakit.microbiome.combine_sbt_indexes()`

- **Purpose:** Merge several SBT indexes into one new SBT index.
- **Method:** Runs `sourmash sbt_combine OUT INPUT... -x BF_SIZE`.
- **API:** `combine_sbt_indexes(output_sbt: 'str', input_sbts: 'Sequence[str]', *, bf_size: 'float' = 100000.0, sourmash_path: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package sourmash (extra `dnakit[bio-sourmash]`). The function runs the `sourmash` command, which is installed with the package; it is looked up on PATH, then next to the Python interpreter, or taken from `sourmash_path`. If it is not found, `BackendUnavailableError` is raised.
- **Input:** `output_sbt` is a new file in an existing folder, normally ending in `.sbt.zip`. `input_sbts` is a non-empty list of existing SBT index files. `bf_size` is the positive Bloom filter size used when the index is saved (default 100000).
- **Output:** A dictionary with `status` (`completed` when the command exits with code 0 and every expected output exists, otherwise `failed`), `backend`, `command` (the arguments, with the program shown as `sourmash`), `return_code`, `outputs` (the expected output paths that exist) and `log_tail` (the last 4000 characters of the command output).
- **Example code:**

```python
import tempfile
from pathlib import Path

import sourmash
from sourmash.sbtmh import create_sbt_index

from dnakit.microbiome import combine_sbt_indexes

work = Path(tempfile.mkdtemp())
signatures = list(sourmash.load_file_as_signatures("genomes.sig"))
inputs = []
for number, chunk in enumerate((signatures[:1], signatures[1:])):
    tree = create_sbt_index(1000)
    for signature in chunk:
        tree.insert(signature)
    inputs.append(str(work / f"part{number}.sbt.zip"))
    tree.save(inputs[-1])
result = combine_sbt_indexes(str(work / "all.sbt.zip"), inputs, verbose=False)
print(result["status"], result["return_code"], [Path(p).name for p in result["outputs"]])
merged = sourmash.load_file_as_index(str(work / "all.sbt.zip"))
print(sorted(signature.name for signature in merged.signatures()))
```

- **Example output:**

```text
completed 0 ['all.sbt.zip']
['GCF_A Escherichia coli', 'GCF_B Escherichia fergusonii', 'GCF_C Bacillus subtilis']
```

## 22) `dnakit.microbiome.migrate_sbt_index()`

- **Purpose:** Rewrite an SBT index in the latest format and keep the original file unchanged.
- **Method:** Creates `output_dir` if needed, copies the index into it under the same file name (the storage tag of a `.sbt.zip` index depends on its file name, so a renamed copy would lose its signatures), and runs `sourmash migrate` on the copy. If the index already has the latest format, the copy is simply rewritten.
- **API:** `migrate_sbt_index(sbt_path: 'str', output_dir: 'str', *, sourmash_path: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package sourmash (extra `dnakit[bio-sourmash]`). The function runs the `sourmash` command, which is installed with the package; it is looked up on PATH, then next to the Python interpreter, or taken from `sourmash_path`. If it is not found, `BackendUnavailableError` is raised.
- **Input:** `sbt_path` is an existing self-contained `.sbt.zip` index. `output_dir` is the folder for the copy; the file name inside it must not exist yet.
- **Output:** A dictionary with `status` (`completed` when the command exits with code 0 and every expected output exists, otherwise `failed`), `backend`, `command` (the arguments, with the program shown as `sourmash`), `return_code`, `outputs` (the expected output paths that exist) and `log_tail` (the last 4000 characters of the command output).
- **Example code:**

```python
import tempfile
from pathlib import Path

import sourmash
from sourmash.sbtmh import create_sbt_index

from dnakit.microbiome import migrate_sbt_index

work = Path(tempfile.mkdtemp())
tree = create_sbt_index(1000)
for signature in sourmash.load_file_as_signatures("genomes.sig"):
    tree.insert(signature)
tree.save(str(work / "genomes.sbt.zip"))
result = migrate_sbt_index(str(work / "genomes.sbt.zip"), str(work / "migrated"), verbose=False)
print(result["status"], result["return_code"], [Path(p).name for p in result["outputs"]])
print(Path(result["outputs"][0]).parent.name)
migrated = sourmash.load_file_as_index(result["outputs"][0])
print(len(list(migrated.signatures())))
```

- **Example output:**

```text
completed 0 ['genomes.sbt.zip']
migrated
3
```

## 23) `dnakit.microbiome.convert_sbt_storage()`

- **Purpose:** Convert the storage of an SBT index to another backend, working on a copy.
- **Method:** Copies the index to `output_sbt` and runs `sourmash storage convert OUT -b BACKEND(PATH)` on the copy. The target path is resolved to an absolute path and must not exist yet.
- **API:** `convert_sbt_storage(sbt_path: 'str', output_sbt: 'str', backend: 'str', *, sourmash_path: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package sourmash (extra `dnakit[bio-sourmash]`). The function runs the `sourmash` command, which is installed with the package; it is looked up on PATH, then next to the Python interpreter, or taken from `sourmash_path`. If it is not found, `BackendUnavailableError` is raised.
- **Input:** `sbt_path` is an existing SBT index file. `output_sbt` is a new file in an existing folder. `backend` is `zip(PATH)` or `fs(PATH)`, where PATH is a new file or folder that receives the signature data; any other form raises `ConfigurationError`.
- **Output:** A dictionary with `status` (`completed` when the command exits with code 0 and every expected output exists, otherwise `failed`), `backend`, `command` (the arguments, with the program shown as `sourmash`), `return_code`, `outputs` (the expected output paths that exist) and `log_tail` (the last 4000 characters of the command output). Both the converted index and the new storage target are expected outputs.
- **Example code:**

```python
import tempfile
from pathlib import Path

import sourmash
from sourmash.sbtmh import create_sbt_index

from dnakit.microbiome import convert_sbt_storage

work = Path(tempfile.mkdtemp())
tree = create_sbt_index(1000)
for signature in sourmash.load_file_as_signatures("genomes.sig"):
    tree.insert(signature)
tree.save(str(work / "genomes.sbt.zip"))
result = convert_sbt_storage(
    str(work / "genomes.sbt.zip"), str(work / "converted.sbt.zip"), f"zip({work / 'signatures.zip'})", verbose=False
)
print(result["status"], result["return_code"], [Path(p).name for p in result["outputs"]])
```

- **Example output:**

```text
completed 0 ['converted.sbt.zip', 'signatures.zip']
```

## 24) `dnakit.microbiome.krona_taxonomy_profile()`

- **Purpose:** Write a Krona table with the fraction of each query assigned to each lineage at one rank.
- **Method:** Runs `sourmash tax metagenome -g GATHER... -t TAXONOMY... -F krona -o BASE --rank RANK`. The output file is `BASE.krona.tsv`; sourmash aggregates the lineages at the rank.
- **API:** `krona_taxonomy_profile(gather_csvs: 'Sequence[str]', taxonomy_csvs: 'Sequence[str]', output_base: 'str', *, rank: 'str', sourmash_path: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package sourmash (extra `dnakit[bio-sourmash]`). The function runs the `sourmash` command, which is installed with the package; it is looked up on PATH, then next to the Python interpreter, or taken from `sourmash_path`. If it is not found, `BackendUnavailableError` is raised.
- **Input:** `gather_csvs` and `taxonomy_csvs` are non-empty lists of existing files: gather CSVs from `sourmash gather` and lineage CSVs (an `ident` column plus one column per rank). `output_base` is a path prefix in an existing folder; the output file is `output_base` plus the suffix given above and must not exist yet. `rank` is required and must be one of `superkingdom`, `phylum`, `class`, `order`, `family`, `genus`, `species` or `strain`.
- **Output:** A dictionary with `status` (`completed` when the command exits with code 0 and every expected output exists, otherwise `failed`), `backend`, `command` (the arguments, with the program shown as `sourmash`), `return_code`, `outputs` (the expected output paths that exist) and `log_tail` (the last 4000 characters of the command output).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.microbiome import krona_taxonomy_profile

base = Path(tempfile.mkdtemp()) / "profile"
result = krona_taxonomy_profile(["gather.csv"], ["tax.csv"], str(base), rank="genus", verbose=False)
print(result["status"], [Path(p).name for p in result["outputs"]])
print(Path(result["outputs"][0]).read_text().rstrip())
```

- **Example output:**

```text
completed ['profile.krona.tsv']
fraction	superkingdom	phylum	class	order	family	genus
0.8	Bacteria	Pseudomonadota	Gammaproteobacteria	Enterobacterales	Enterobacteriaceae	Escherichia
0.2	Bacteria	Bacillota	Bacilli	Bacillales	Bacillaceae	Bacillus
```

## 25) `dnakit.microbiome.csv_taxonomy_summary()`

- **Purpose:** Write the per-rank summary of gather results as a CSV file.
- **Method:** Runs `sourmash tax metagenome -g GATHER... -t TAXONOMY... -F csv_summary -o BASE`. The output file is `BASE.summarized.csv`, with the columns `query_name`, `rank`, `fraction`, `lineage`, the query identifiers and the matching statistics.
- **API:** `csv_taxonomy_summary(gather_csvs: 'Sequence[str]', taxonomy_csvs: 'Sequence[str]', output_base: 'str', *, sourmash_path: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package sourmash (extra `dnakit[bio-sourmash]`). The function runs the `sourmash` command, which is installed with the package; it is looked up on PATH, then next to the Python interpreter, or taken from `sourmash_path`. If it is not found, `BackendUnavailableError` is raised.
- **Input:** `gather_csvs` and `taxonomy_csvs` are non-empty lists of existing files: gather CSVs from `sourmash gather` and lineage CSVs (an `ident` column plus one column per rank). `output_base` is a path prefix in an existing folder; the output file is `output_base` plus the suffix given above and must not exist yet.
- **Output:** A dictionary with `status` (`completed` when the command exits with code 0 and every expected output exists, otherwise `failed`), `backend`, `command` (the arguments, with the program shown as `sourmash`), `return_code`, `outputs` (the expected output paths that exist) and `log_tail` (the last 4000 characters of the command output).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.microbiome import csv_taxonomy_summary

base = Path(tempfile.mkdtemp()) / "profile"
result = csv_taxonomy_summary(["gather.csv"], ["tax.csv"], str(base), verbose=False)
print(result["status"], [Path(p).name for p in result["outputs"]])
lines = Path(result["outputs"][0]).read_text().splitlines()
print(lines[0])
for line in lines[1:]:
    fields = line.split(",")
    if fields[1] == "genus":
        print(fields[0], fields[1], fields[2], fields[3].split(";")[-1])
```

- **Example output:**

```text
completed ['profile.summarized.csv']
query_name,rank,fraction,lineage,query_md5,query_filename,f_weighted_at_rank,bp_match_at_rank,query_ani_at_rank,total_weighted_hashes
query1 genus 0.8 Escherichia
query1 genus 0.2 Bacillus
```

## 26) `dnakit.microbiome.human_taxonomy_summary()`

- **Purpose:** Write a human-readable table of the lineages found in each query.
- **Method:** Runs `sourmash tax metagenome -g GATHER... -t TAXONOMY... -F human -o BASE`, with `--rank RANK` when a rank is given. The output file is `BASE.human.txt`. Without a rank, sourmash reports the lineages it finds.
- **API:** `human_taxonomy_summary(gather_csvs: 'Sequence[str]', taxonomy_csvs: 'Sequence[str]', output_base: 'str', *, rank: 'str | None' = None, sourmash_path: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package sourmash (extra `dnakit[bio-sourmash]`). The function runs the `sourmash` command, which is installed with the package; it is looked up on PATH, then next to the Python interpreter, or taken from `sourmash_path`. If it is not found, `BackendUnavailableError` is raised.
- **Input:** `gather_csvs` and `taxonomy_csvs` are non-empty lists of existing files: gather CSVs from `sourmash gather` and lineage CSVs (an `ident` column plus one column per rank). `output_base` is a path prefix in an existing folder; the output file is `output_base` plus the suffix given above and must not exist yet. `rank` is optional and, when given, one of the eight ranks listed for function 24.
- **Output:** A dictionary with `status` (`completed` when the command exits with code 0 and every expected output exists, otherwise `failed`), `backend`, `command` (the arguments, with the program shown as `sourmash`), `return_code`, `outputs` (the expected output paths that exist) and `log_tail` (the last 4000 characters of the command output).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.microbiome import human_taxonomy_summary

base = Path(tempfile.mkdtemp()) / "profile"
result = human_taxonomy_summary(["gather.csv"], ["tax.csv"], str(base), rank="genus", verbose=False)
print(result["status"], [Path(p).name for p in result["outputs"]])
print(Path(result["outputs"][0]).read_text().rstrip())
```

- **Example output:**

```text
completed ['profile.human.txt']
sample name    proportion   cANI   lineage
-----------    ----------   ----   -------
query1            80.0%     99.3%  Bacteria;Pseudomonadota;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia
query1            20.0%     94.9%  Bacteria;Bacillota;Bacilli;Bacillales;Bacillaceae;Bacillus
```

## 27) `dnakit.microbiome.lineage_taxonomy_summary()`

- **Purpose:** Write the fraction of each query assigned to each lineage at one rank as a table with one column per query.
- **Method:** Runs `sourmash tax metagenome -g GATHER... -t TAXONOMY... -F lineage_summary -o BASE --rank RANK`. The output file is `BASE.lineage_summary.tsv`.
- **API:** `lineage_taxonomy_summary(gather_csvs: 'Sequence[str]', taxonomy_csvs: 'Sequence[str]', output_base: 'str', *, rank: 'str', sourmash_path: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package sourmash (extra `dnakit[bio-sourmash]`). The function runs the `sourmash` command, which is installed with the package; it is looked up on PATH, then next to the Python interpreter, or taken from `sourmash_path`. If it is not found, `BackendUnavailableError` is raised.
- **Input:** `gather_csvs` and `taxonomy_csvs` are non-empty lists of existing files: gather CSVs from `sourmash gather` and lineage CSVs (an `ident` column plus one column per rank). `output_base` is a path prefix in an existing folder; the output file is `output_base` plus the suffix given above and must not exist yet. `rank` is required and is one of the eight ranks listed for function 24.
- **Output:** A dictionary with `status` (`completed` when the command exits with code 0 and every expected output exists, otherwise `failed`), `backend`, `command` (the arguments, with the program shown as `sourmash`), `return_code`, `outputs` (the expected output paths that exist) and `log_tail` (the last 4000 characters of the command output).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.microbiome import lineage_taxonomy_summary

base = Path(tempfile.mkdtemp()) / "profile"
result = lineage_taxonomy_summary(["gather.csv"], ["tax.csv"], str(base), rank="genus", verbose=False)
print(result["status"], [Path(p).name for p in result["outputs"]])
print(Path(result["outputs"][0]).read_text().rstrip())
```

- **Example output:**

```text
completed ['profile.lineage_summary.tsv']
lineage	query1
Bacteria;Bacillota;Bacilli;Bacillales;Bacillaceae;Bacillus	0.2
Bacteria;Pseudomonadota;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia	0.8
```

## 28) `dnakit.microbiome.bioboxes_taxonomy_profile()`

- **Purpose:** Write a CAMI bioboxes taxonomic profile of the gather results.
- **Method:** Checks that the header of every lineage CSV has a `taxpath` column, then runs `sourmash tax metagenome -g GATHER... -t TAXONOMY... -F bioboxes -o BASE`. The output file is `BASE.bioboxes.profile`. The `taxpath` column holds the NCBI taxids of the lineage, separated by `|`, one per rank. In sourmash 4.9.4 a reported lineage without a taxid makes the command fail with a `TypeError`, and the function then returns status `failed` with the log tail.
- **API:** `bioboxes_taxonomy_profile(gather_csvs: 'Sequence[str]', taxonomy_csvs: 'Sequence[str]', output_base: 'str', *, sourmash_path: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package sourmash (extra `dnakit[bio-sourmash]`). The function runs the `sourmash` command, which is installed with the package; it is looked up on PATH, then next to the Python interpreter, or taken from `sourmash_path`. If it is not found, `BackendUnavailableError` is raised.
- **Input:** `gather_csvs` and `taxonomy_csvs` are non-empty lists of existing files: gather CSVs from `sourmash gather` and lineage CSVs (an `ident` column plus one column per rank). `output_base` is a path prefix in an existing folder; the output file is `output_base` plus the suffix given above and must not exist yet. Every lineage CSV needs a `taxpath` column, otherwise `ConfigurationError` is raised.
- **Output:** A dictionary with `status` (`completed` when the command exits with code 0 and every expected output exists, otherwise `failed`), `backend`, `command` (the arguments, with the program shown as `sourmash`), `return_code`, `outputs` (the expected output paths that exist) and `log_tail` (the last 4000 characters of the command output).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.microbiome import bioboxes_taxonomy_profile

base = Path(tempfile.mkdtemp()) / "profile"
result = bioboxes_taxonomy_profile(["gather.csv"], ["taxid.csv"], str(base), verbose=False)
print(result["status"], [Path(p).name for p in result["outputs"]])
lines = Path(result["outputs"][0]).read_text().splitlines()
print(lines[1:5])
for line in lines[6:]:
    fields = line.split("\t")
    if fields[1] == "species":
        print(fields[0], fields[1], fields[2], fields[4])
```

- **Example output:**

```text
completed ['profile.bioboxes.profile']
['@SampleID:query1', '@Version:0.10.0', '@Ranks:superkingdom|phylum|class|order|family|genus|species|strain', '@__program__:sourmash']
562 species 2|1224|1236|91347|543|561|562 80.00
1423 species 2|1239|91061|1385|186817|1386|1423 20.00
```

## 29) `dnakit.microbiome.kreport_taxonomy_report()`

- **Purpose:** Write a Kraken-style report of the gather results.
- **Method:** Runs `sourmash tax metagenome -g GATHER... -t TAXONOMY... -F kreport -o BASE`. The output file is `BASE.kreport.txt`, with the columns percentage, base pairs contained, base pairs assigned, rank code, taxid and name. The taxid column is empty, and an unclassified part of the query is reported as `U`.
- **API:** `kreport_taxonomy_report(gather_csvs: 'Sequence[str]', taxonomy_csvs: 'Sequence[str]', output_base: 'str', *, sourmash_path: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package sourmash (extra `dnakit[bio-sourmash]`). The function runs the `sourmash` command, which is installed with the package; it is looked up on PATH, then next to the Python interpreter, or taken from `sourmash_path`. If it is not found, `BackendUnavailableError` is raised.
- **Input:** `gather_csvs` and `taxonomy_csvs` are non-empty lists of existing files: gather CSVs from `sourmash gather` and lineage CSVs (an `ident` column plus one column per rank). `output_base` is a path prefix in an existing folder; the output file is `output_base` plus the suffix given above and must not exist yet. The gather CSV must carry a non-zero `total_weighted_hashes` column (written by sourmash 4.5 and later); otherwise sourmash stops with a `ValueError`, the status is `failed`, and an empty output file is left behind.
- **Output:** A dictionary with `status` (`completed` when the command exits with code 0 and every expected output exists, otherwise `failed`), `backend`, `command` (the arguments, with the program shown as `sourmash`), `return_code`, `outputs` (the expected output paths that exist) and `log_tail` (the last 4000 characters of the command output).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.microbiome import kreport_taxonomy_report

base = Path(tempfile.mkdtemp()) / "profile"
result = kreport_taxonomy_report(["gather.csv"], ["tax.csv"], str(base), verbose=False)
print(result["status"], [Path(p).name for p in result["outputs"]])
for line in Path(result["outputs"][0]).read_text().splitlines():
    if line.split("\t")[3] in ("D", "S"):
        print(line)
```

- **Example output:**

```text
completed ['profile.kreport.txt']
100.00	10000	0	D		Bacteria
80.00	8000	8000	S		Escherichia coli
20.00	2000	2000	S		Bacillus subtilis
```

## 30) `dnakit.microbiome.lingroup_taxonomy_report()`

- **Purpose:** Summarise the gather results by LIN lingroup.
- **Method:** Checks that the header of every lineage CSV has a `lin` column, then runs `sourmash tax metagenome -g GATHER... -t TAXONOMY... -F lingroup -o BASE --lins --lingroup LINGROUP_CSV`. The output file is `BASE.lingroup.tsv`, with the columns `name`, `lin`, `percent_containment` and `num_bp_contained`.
- **API:** `lingroup_taxonomy_report(gather_csvs: 'Sequence[str]', lin_taxonomy_csvs: 'Sequence[str]', lingroup_csv: 'str', output_base: 'str', *, sourmash_path: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package sourmash (extra `dnakit[bio-sourmash]`). The function runs the `sourmash` command, which is installed with the package; it is looked up on PATH, then next to the Python interpreter, or taken from `sourmash_path`. If it is not found, `BackendUnavailableError` is raised.
- **Input:** `gather_csvs` is a non-empty list of gather CSVs and `lin_taxonomy_csvs` a non-empty list of lineage CSVs that carry an `ident` and a `lin` column (otherwise `ConfigurationError`). `lingroup_csv` is an existing CSV with the columns `name` and `lin`, giving the LIN prefix of each group. `output_base` is a path prefix in an existing folder; `BASE.lingroup.tsv` must not exist yet. As for function 29, the gather CSV needs a non-zero `total_weighted_hashes` column.
- **Output:** A dictionary with `status` (`completed` when the command exits with code 0 and every expected output exists, otherwise `failed`), `backend`, `command` (the arguments, with the program shown as `sourmash`), `return_code`, `outputs` (the expected output paths that exist) and `log_tail` (the last 4000 characters of the command output).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.microbiome import lingroup_taxonomy_report

base = Path(tempfile.mkdtemp()) / "profile"
result = lingroup_taxonomy_report(["gather.csv"], ["lin.csv"], "lingroup.csv", str(base), verbose=False)
print(result["status"], [Path(p).name for p in result["outputs"]])
header, *rows = Path(result["outputs"][0]).read_text().splitlines()
print(header)
for row in sorted(rows):
    print(row)
```

- **Example output:**

```text
completed ['profile.lingroup.tsv']
name	lin	percent_containment	num_bp_contained
group_bacillus_like	1	20.00	2000
group_ecoli_like	0	80.00	8000
```
