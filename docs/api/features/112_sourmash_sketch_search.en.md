# sourmash: sketch, compare, search and index DNA signatures

This group of functions in `dnakit.microbiome` is a set of typed facades over the command line of [sourmash](https://github.com/sourmash-bio/sourmash), which sketches DNA into small MinHash signatures and compares, searches, indexes and classifies them. The group covers sketching a FASTA file (DNA and six-frame translated protein), pairwise comparison and plots, search, gather, prefetch and categorize against a database, set operations on signatures (merge, intersect, subtract, downsample, flatten, filter), bookkeeping (describe, manifest, split, extract, rename, picklist checks), Mash import and export, and lowest-common-ancestor (LCA) taxonomy databases. Each function builds a fixed `sourmash` command line without a shell, runs it in the folder of its output, and returns a JSON-serialisable dictionary; output files must be new paths. None of them downloads data or contacts a server.

sourmash is BSD-3-Clause licensed (checked from the package metadata), so there is no copyleft caveat. Import the functions from `dnakit.microbiome`, for example `from dnakit.microbiome import sketch_dna`. The dependency is optional and the program is looked up only when a function is called: if `sourmash` is not on PATH and no `sourmash_path` is given, the function raises `BackendUnavailableError` with the install command. All examples pass `sourmash_path` found with `shutil.which`, fall back to the `sourmash` script next to the Python interpreter, work in a temporary folder and pass `verbose=False`.

## Example data

The examples use two single-record FASTA files of 60 and 65 bases. The first 35 bases are identical, so at k-mer size 21 the two sketches share 15 of their 40 and 45 k-mers (made with `scaled=1`, which keeps every hash). The examples sketch these files themselves into a temporary folder.

**Example file:** `a.fa`

```text
>seqA
GCTAAAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCG
```

**Example file:** `b.fa`

```text
>seqB
GCTAAAGACAATTACATAACATACACGTCAGCACGCTTAAGGGTTAAGTAAGTGTGATGCATACG
```

## 1) `dnakit.microbiome.check_signatures_against_picklist()`

- **Purpose:** Check which rows of a picklist are present in signature files and list the missing ones.
- **Method:** Run `sourmash sig check <signatures> --picklist <picklist> -o <missing csv> --no-require-manifest`, adding `-m <manifest>` when `output_manifest` is given. Exit codes 0 and 1 both count as completed. The function additionally reads the missing-rows CSV and sets `all_picklist_rows_found` only when the exit code is 0 and the CSV has no data row.
- **API:** `dnakit.microbiome.check_signatures_against_picklist(signatures: 'Sequence[str]', picklist: 'str', output_missing_csv: 'str', *, output_manifest: 'str | None' = None, sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `sourmash` on PATH (or its location in `sourmash_path`). It is provided by the optional package sourmash 4.9 or later (extra `dnakit[bio-sourmash]`); without it, `BackendUnavailableError` is raised.
- **Input:** `signatures` is a non-empty list of existing signature files. `picklist` is passed to sourmash unchanged, in the form `file.csv:column:coltype` (for example `names.csv:name:name`); it is not checked by DNAKit, and a relative path is resolved against the folder of `output_missing_csv`, so an absolute path is safest. `output_missing_csv` (and `output_manifest`, if given) must not exist yet.
- **Output:** The standard result dictionary (here `status` is `completed` for exit code 0 or 1, and `outputs` lists the missing-rows CSV and the optional manifest) plus `all_picklist_rows_found`, a boolean.
- **Example code:**

```python
import shutil, sys, tempfile
from pathlib import Path

from dnakit.microbiome import sketch_dna, check_signatures_against_picklist

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())
sketch_dna("a.fa", work / "a.sig", scaled=1, name="sigA first", sourmash_path=tool, verbose=False)
sketch_dna("b.fa", work / "b.sig", scaled=1, name="sigB second", sourmash_path=tool, verbose=False)

(work / "wanted.csv").write_text("name\nsigA first\nnot a signature\n")
result = check_signatures_against_picklist(
    [work / "a.sig", work / "b.sig"], f"{work / 'wanted.csv'}:name:name", work / "missing.csv",
    sourmash_path=tool, verbose=False,
)
print(result["status"], result["return_code"], result["all_picklist_rows_found"])
print((work / "missing.csv").read_text().splitlines())
```

- **Example output:**

```text
completed 0 False
['name', 'not a signature']
```

## 2) `dnakit.microbiome.classify_with_lca()`

- **Purpose:** Assign a taxonomic lineage to each query signature using lowest-common-ancestor (LCA) databases.
- **Method:** Run `sourmash lca classify --query <queries> --db <databases> --threshold <n> -o <csv>`, adding `--majority` when `majority` is true. Each query row reports the lineage supported by the hashes found in the databases.
- **API:** `dnakit.microbiome.classify_with_lca(queries: 'Sequence[str]', databases: 'Sequence[str]', output_csv: 'str', *, threshold: 'int' = 5, majority: 'bool' = False, sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `sourmash` on PATH (or its location in `sourmash_path`). It is provided by the optional package sourmash 4.9 or later (extra `dnakit[bio-sourmash]`); without it, `BackendUnavailableError` is raised.
- **Input:** `queries` and `databases` are non-empty lists of existing files; the databases are LCA JSON files such as those written by `index_lca_database()`. `threshold` is the minimum number of hashes that must support a lineage (default 5; it is not checked by DNAKit). `majority` uses majority-vote classification instead of the lowest common ancestor. `output_csv` must not exist yet.
- **Output:** The standard result dictionary: `status` (`completed` when sourmash exits with 0 and every expected output file exists, otherwise `failed`), `backend` (`sourmash`), `command`, `return_code`, `outputs` (the existing output files) and `log_tail` (the last 4000 characters of sourmash's output, which contain terminal control codes). The CSV has the columns `ID`, `status`, and one column per taxonomic rank.
- **Example code:**

```python
import shutil, sys, tempfile
from pathlib import Path
import csv

from dnakit.microbiome import sketch_dna, classify_with_lca, index_lca_database

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())
sketch_dna("a.fa", work / "a.sig", scaled=1, name="sigA first", sourmash_path=tool, verbose=False)
sketch_dna("b.fa", work / "b.sig", scaled=1, name="sigB second", sourmash_path=tool, verbose=False)

(work / "lineages.csv").write_text(
    "accession,superkingdom,phylum,class,order,family,genus,species\n"
    "sigA,Bacteria,Firmicutes,Bacilli,Order1,Family1,Genusone,Genusone alpha\n"
    "sigB,Bacteria,Firmicutes,Bacilli,Order1,Family1,Genusone,Genusone beta\n"
)
index_lca_database([work / "a.sig", work / "b.sig"], work / "lineages.csv", work / "lca.lca.json", scaled=1, sourmash_path=tool, verbose=False)
result = classify_with_lca([work / "a.sig"], [work / "lca.lca.json"], work / "lca.csv", threshold=1, sourmash_path=tool, verbose=False)
rows = list(csv.DictReader(open(work / "lca.csv")))
print(result["status"], len(rows), rows[0]["status"], rows[0]["genus"], rows[0]["species"])
```

- **Example output:**

```text
completed 1 found Genusone Genusone alpha
```

## 3) `dnakit.microbiome.collect_signature_manifests()`

- **Purpose:** Collect the manifests of many signature files into one standalone manifest database.
- **Method:** Run `sourmash sig collect <signatures> -o <manifest> --no-require-manifest`. The output is a SQLite manifest (`.sqlmf`) that points at the original files; the signatures themselves are not copied.
- **API:** `dnakit.microbiome.collect_signature_manifests(signatures: 'Sequence[str]', output_manifest: 'str', *, sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `sourmash` on PATH (or its location in `sourmash_path`). It is provided by the optional package sourmash 4.9 or later (extra `dnakit[bio-sourmash]`); without it, `BackendUnavailableError` is raised.
- **Input:** `signatures` is a non-empty list of existing signature files. `output_manifest` must not exist yet; use the `.sqlmf` extension so that sourmash writes SQLite.
- **Output:** The standard result dictionary: `status` (`completed` when sourmash exits with 0 and every expected output file exists, otherwise `failed`), `backend` (`sourmash`), `command`, `return_code`, `outputs` (the existing output files) and `log_tail` (the last 4000 characters of sourmash's output, which contain terminal control codes).
- **Example code:**

```python
import shutil, sys, tempfile
from pathlib import Path
import sqlite3

from dnakit.microbiome import sketch_dna, collect_signature_manifests

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())
sketch_dna("a.fa", work / "a.sig", scaled=1, name="sigA first", sourmash_path=tool, verbose=False)
sketch_dna("b.fa", work / "b.sig", scaled=1, name="sigB second", sourmash_path=tool, verbose=False)

result = collect_signature_manifests([work / "a.sig", work / "b.sig"], work / "collection.sqlmf", sourmash_path=tool, verbose=False)
connection = sqlite3.connect(work / "collection.sqlmf")
print(result["status"], connection.execute("select count(*) from sourmash_sketches").fetchone()[0])
connection.close()
```

- **Example output:**

```text
completed 2
```

## 4) `dnakit.microbiome.index_lca_database()`

- **Purpose:** Build a lowest-common-ancestor (LCA) taxonomy database from signatures and a lineage CSV.
- **Method:** Run `sourmash lca index <lineages csv> <database> <signatures> -k <ksize> --scaled <scaled> --dna`, adding `--split-identifiers` when `split_identifiers` is true. With the split option, only the first word of each signature name (up to the first space or dot) is matched to the lineage identifier. sourmash appends `.lca.json` to a database name that does not already end in `.lca.json` or `.lca.json.gz`, but DNAKit looks for exactly the path you gave, so a name such as `db.json` is reported as `failed` (return code 0, no `outputs`) although `db.json.lca.json` was written. Give a name that ends in `.lca.json`.
- **API:** `dnakit.microbiome.index_lca_database(signatures: 'Sequence[str]', lineages_csv: 'str', output_database: 'str', *, ksize: 'int' = 21, scaled: 'int' = 1000, split_identifiers: 'bool' = True, sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `sourmash` on PATH (or its location in `sourmash_path`). It is provided by the optional package sourmash 4.9 or later (extra `dnakit[bio-sourmash]`); without it, `BackendUnavailableError` is raised.
- **Input:** `signatures` is a non-empty list of existing signature files, `lineages_csv` an existing CSV whose first column is an identifier column (`accession`, `ident` or `identifiers`) followed by rank columns (`superkingdom,phylum,class,order,family,genus,species`, optionally `strain`). `ksize` and `scaled` must match the signatures (defaults 21 and 1000). `output_database` must not exist yet and should end in `.lca.json`.
- **Output:** The standard result dictionary: `status` (`completed` when sourmash exits with 0 and every expected output file exists, otherwise `failed`), `backend` (`sourmash`), `command`, `return_code`, `outputs` (the existing output files) and `log_tail` (the last 4000 characters of sourmash's output, which contain terminal control codes).
- **Example code:**

```python
import shutil, sys, tempfile
from pathlib import Path
import json

from dnakit.microbiome import sketch_dna, index_lca_database

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())
sketch_dna("a.fa", work / "a.sig", scaled=1, name="sigA first", sourmash_path=tool, verbose=False)
sketch_dna("b.fa", work / "b.sig", scaled=1, name="sigB second", sourmash_path=tool, verbose=False)

(work / "lineages.csv").write_text(
    "accession,superkingdom,phylum,class,order,family,genus,species\n"
    "sigA,Bacteria,Firmicutes,Bacilli,Order1,Family1,Genusone,Genusone alpha\n"
    "sigB,Bacteria,Firmicutes,Bacilli,Order1,Family1,Genusone,Genusone beta\n"
)
result = index_lca_database([work / "a.sig", work / "b.sig"], work / "lineages.csv", work / "lca.lca.json", scaled=1, sourmash_path=tool, verbose=False)
data = json.loads((work / "lca.lca.json").read_text())
print(result["status"], data["ksize"], data["scaled"], len(data["lid_to_lineage"]), len(data["hashval_to_idx"]))
```

- **Example output:**

```text
completed 21 1 2 70
```

## 5) `dnakit.microbiome.ingest_mash_csv()`

- **Purpose:** Convert Mash CSV hash dumps into a sourmash signature file.
- **Method:** Run `sourmash sig ingest --csv <files> -o <signature>`. Each input line becomes one `num` sketch.
- **API:** `dnakit.microbiome.ingest_mash_csv(csv_files: 'Sequence[str]', output_signature: 'str', *, sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `sourmash` on PATH (or its location in `sourmash_path`). It is provided by the optional package sourmash 4.9 or later (extra `dnakit[bio-sourmash]`); without it, `BackendUnavailableError` is raised.
- **Input:** `csv_files` is a non-empty list of existing files. In sourmash 4.9 each line must have exactly five comma-separated fields: `murmur64`, `42`, the k-mer size, a name, and all hashes as integers separated by single spaces. `output_signature` must not exist yet.
- **Output:** The standard result dictionary: `status` (`completed` when sourmash exits with 0 and every expected output file exists, otherwise `failed`), `backend` (`sourmash`), `command`, `return_code`, `outputs` (the existing output files) and `log_tail` (the last 4000 characters of sourmash's output, which contain terminal control codes).
- **Example code:**

```python
import shutil, sys, tempfile
from pathlib import Path
from sourmash import load_file_as_signatures

from dnakit.microbiome import sketch_dna, ingest_mash_csv

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())

(work / "mash.csv").write_text("murmur64,42,21,mash sample,101 202 303 404 505\n")
result = ingest_mash_csv([work / "mash.csv"], work / "mash.sig", sourmash_path=tool, verbose=False)
(sig,) = load_file_as_signatures(str(work / "mash.sig"))
print(result["status"], sig.minhash.ksize, sig.minhash.num, sorted(sig.minhash.hashes))
```

- **Example output:**

```text
completed 21 5 [101, 202, 303, 404, 505]
```

## 6) `dnakit.microbiome.categorize_signatures()`

- **Purpose:** Assign each query signature to its best match in a database.
- **Method:** Run `sourmash categorize <database> <queries> -k <ksize> --dna --threshold <t> --csv <csv>`. The CSV has no header line; the columns are the query file, the query name, the best match name and the similarity.
- **API:** `dnakit.microbiome.categorize_signatures(database: 'str', queries: 'Sequence[str]', output_csv: 'str', *, ksize: 'int' = 21, threshold: 'float' = 0.0, sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `sourmash` on PATH (or its location in `sourmash_path`). It is provided by the optional package sourmash 4.9 or later (extra `dnakit[bio-sourmash]`); without it, `BackendUnavailableError` is raised.
- **Input:** `database` is an existing signature collection, `queries` a non-empty list of existing signature files. `threshold` is a similarity between 0 and 1 (default 0.0), `ksize` defaults to 21. `output_csv` must not exist yet.
- **Output:** The standard result dictionary: `status` (`completed` when sourmash exits with 0 and every expected output file exists, otherwise `failed`), `backend` (`sourmash`), `command`, `return_code`, `outputs` (the existing output files) and `log_tail` (the last 4000 characters of sourmash's output, which contain terminal control codes).
- **Example code:**

```python
import shutil, sys, tempfile
from pathlib import Path
import csv

from dnakit.microbiome import sketch_dna, categorize_signatures

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())
sketch_dna("a.fa", work / "a.sig", scaled=1, name="sigA first", sourmash_path=tool, verbose=False)
sketch_dna("b.fa", work / "b.sig", scaled=1, name="sigB second", sourmash_path=tool, verbose=False)

result = categorize_signatures(work / "b.sig", [work / "a.sig"], work / "categories.csv", threshold=0.1, sourmash_path=tool, verbose=False)
rows = list(csv.reader(open(work / "categories.csv")))
print(result["status"], len(rows), rows[0][1], rows[0][2], round(float(rows[0][3]), 3))
```

- **Example output:**

```text
completed 1 sigA first sigB second 0.214
```

## 7) `dnakit.microbiome.export_signature_to_mash()`

- **Purpose:** Export a signature to the Mash JSON dump format.
- **Method:** Run `sourmash sig export <signature> -k <ksize> -o <json>`. The JSON holds `kmer`, `sketchSize`, `hashType`, `hashBits`, `hashSeed` and a `sketches` list with the hash values.
- **API:** `dnakit.microbiome.export_signature_to_mash(signature: 'str', output_json: 'str', *, ksize: 'int' = 21, sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `sourmash` on PATH (or its location in `sourmash_path`). It is provided by the optional package sourmash 4.9 or later (extra `dnakit[bio-sourmash]`); without it, `BackendUnavailableError` is raised.
- **Input:** `signature` is an existing signature file containing a sketch of the requested `ksize` (default 21). `output_json` must not exist yet.
- **Output:** The standard result dictionary: `status` (`completed` when sourmash exits with 0 and every expected output file exists, otherwise `failed`), `backend` (`sourmash`), `command`, `return_code`, `outputs` (the existing output files) and `log_tail` (the last 4000 characters of sourmash's output, which contain terminal control codes).
- **Example code:**

```python
import shutil, sys, tempfile
from pathlib import Path
import json

from dnakit.microbiome import sketch_dna, export_signature_to_mash

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())
sketch_dna("a.fa", work / "a.sig", scaled=1, name="sigA first", sourmash_path=tool, verbose=False)
sketch_dna("b.fa", work / "b.sig", scaled=1, name="sigB second", sourmash_path=tool, verbose=False)

result = export_signature_to_mash(work / "a.sig", work / "a_mash.json", sourmash_path=tool, verbose=False)
data = json.loads((work / "a_mash.json").read_text())
print(result["status"], data["kmer"], data["sketchSize"], len(data["sketches"][0]["hashes"]))
```

- **Example output:**

```text
completed 21 40 40
```

## 8) `dnakit.microbiome.extract_signature_kmers()`

- **Purpose:** Write the sequences that contain the k-mers hashed in a signature.
- **Method:** Run `sourmash sig kmers --signatures <signature> --sequences <files> -k <ksize> --save-sequences <output>`. The output holds the sequences that contain k-mers whose hashes are in the signature; a signature made with a large `scaled` contains only a small fraction of all hashes, so the example uses `scaled=1`.
- **API:** `dnakit.microbiome.extract_signature_kmers(signature: 'str', sequences: 'Sequence[str]', output_sequences: 'str', *, ksize: 'int' = 21, sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `sourmash` on PATH (or its location in `sourmash_path`). It is provided by the optional package sourmash 4.9 or later (extra `dnakit[bio-sourmash]`); without it, `BackendUnavailableError` is raised.
- **Input:** `signature` is an existing signature file, `sequences` a non-empty list of existing FASTA or FASTQ files, `ksize` defaults to 21. `output_sequences` must not exist yet.
- **Output:** The standard result dictionary: `status` (`completed` when sourmash exits with 0 and every expected output file exists, otherwise `failed`), `backend` (`sourmash`), `command`, `return_code`, `outputs` (the existing output files) and `log_tail` (the last 4000 characters of sourmash's output, which contain terminal control codes).
- **Example code:**

```python
import shutil, sys, tempfile
from pathlib import Path

from dnakit.microbiome import sketch_dna, extract_signature_kmers

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())
sketch_dna("a.fa", work / "a.sig", scaled=1, name="sigA first", sourmash_path=tool, verbose=False)
sketch_dna("b.fa", work / "b.sig", scaled=1, name="sigB second", sourmash_path=tool, verbose=False)

result = extract_signature_kmers(work / "b.sig", ["a.fa"], work / "matches.fa", sourmash_path=tool, verbose=False)
print(result["status"], (work / "matches.fa").read_text().splitlines()[0])
```

- **Example output:**

```text
completed >seqA
```

## 9) `dnakit.microbiome.multigather_signatures()`

- **Purpose:** Run a gather search for every query signature against the databases.
- **Method:** Run `sourmash multigather --query <queries> --db <databases> -k <ksize> --dna --threshold-bp <bp> --output-dir <dir>`. DNAKit creates the output folder if it is missing, lists its contents before and after, and reports the new files as `outputs`; the status depends only on the exit code.
- **API:** `dnakit.microbiome.multigather_signatures(queries: 'Sequence[str]', databases: 'Sequence[str]', output_dir: 'str', *, ksize: 'int' = 21, threshold_bp: 'int' = 0, sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `sourmash` on PATH (or its location in `sourmash_path`). It is provided by the optional package sourmash 4.9 or later (extra `dnakit[bio-sourmash]`); without it, `BackendUnavailableError` is raised.
- **Input:** `queries` and `databases` are non-empty lists of existing signature files. `output_dir` is created if needed and may already contain files. `ksize` defaults to 21 and `threshold_bp` (minimum overlap in base pairs) to 0.
- **Output:** The standard result dictionary, with `status` set to `completed` when the exit code is 0, and `outputs` listing the files that appeared in `output_dir` (for each query, named after its source file: a `.csv`, a `.matches.sig` and an `.unassigned.sig` file).
- **Example code:**

```python
import shutil, sys, tempfile
from pathlib import Path

from dnakit.microbiome import sketch_dna, multigather_signatures

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())
sketch_dna("a.fa", work / "a.sig", scaled=1, name="sigA first", sourmash_path=tool, verbose=False)
sketch_dna("b.fa", work / "b.sig", scaled=1, name="sigB second", sourmash_path=tool, verbose=False)

result = multigather_signatures([work / "a.sig"], [work / "b.sig"], work / "gathers", sourmash_path=tool, verbose=False)
print(result["status"], [Path(name).suffix for name in result["outputs"]])
```

- **Example output:**

```text
completed ['.csv', '.sig', '.sig']
```

## 10) `dnakit.microbiome.overlap_signatures()`

- **Purpose:** Report the Jaccard similarity, containment and shared and unique hash counts of two signatures.
- **Method:** Run `sourmash sig overlap <a> <b>` and parse every output line of the form `name: number` into a dictionary of floats. Each file must contain exactly one signature, otherwise sourmash fails and `status` is `failed`. Note that the `size` and `sum hashes` lines are printed once per signature, so the dictionary keeps the values of the second signature.
- **API:** `dnakit.microbiome.overlap_signatures(signature_a: 'str', signature_b: 'str', *, sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `sourmash` on PATH (or its location in `sourmash_path`). It is provided by the optional package sourmash 4.9 or later (extra `dnakit[bio-sourmash]`); without it, `BackendUnavailableError` is raised.
- **Input:** `signature_a` and `signature_b` are existing signature files with one signature each and the same k-mer size and scaled value.
- **Output:** `status`, `backend`, `return_code`, `values` (floats keyed by sourmash's labels such as `jaccard similarity`, `number of hashes in common`, `only in first`, `only in second` and `total (union)`) and `report` (the last 4000 characters of the raw output).
- **Example code:**

```python
import shutil, sys, tempfile
from pathlib import Path

from dnakit.microbiome import sketch_dna, overlap_signatures

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())
sketch_dna("a.fa", work / "a.sig", scaled=1, name="sigA first", sourmash_path=tool, verbose=False)
sketch_dna("b.fa", work / "b.sig", scaled=1, name="sigB second", sourmash_path=tool, verbose=False)

result = overlap_signatures(work / "a.sig", work / "b.sig", sourmash_path=tool, verbose=False)
print(result["status"])
values = result["values"]
for key in ("jaccard similarity", "number of hashes in common", "only in first", "only in second", "total (union)"):
    print(key, values[key])
```

- **Example output:**

```text
completed
jaccard similarity 0.21429
number of hashes in common 15.0
only in first 25.0
only in second 30.0
total (union) 70.0
```

## 11) `dnakit.microbiome.plot_comparison_matrix()`

- **Purpose:** Draw the clustered similarity heat map of a sourmash compare matrix as a PNG.
- **Method:** Run `sourmash plot <matrix> --output-dir <dir>`. sourmash reads the `.labels.txt` file that `compare_signatures()` writes next to the matrix, clusters the samples and draws the heat map with matplotlib and scipy (both installed with sourmash). DNAKit creates the folder if needed and reports the new files.
- **API:** `dnakit.microbiome.plot_comparison_matrix(comparison_matrix: 'str', output_dir: 'str', *, sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `sourmash` on PATH (or its location in `sourmash_path`). It is provided by the optional package sourmash 4.9 or later (extra `dnakit[bio-sourmash]`); without it, `BackendUnavailableError` is raised and, for plotting, matplotlib and scipy, which are installed with sourmash.
- **Input:** `comparison_matrix` is an existing matrix file written by `compare_signatures()` (its `.labels.txt` file must be beside it). `output_dir` is created if needed.
- **Output:** The standard result dictionary, with `outputs` listing all new files, `status` `completed` only when the exit code is 0 and a PNG was written, and `image`, the path of the first new PNG (or `None`).
- **Example code:**

```python
import shutil, sys, tempfile
from pathlib import Path

from dnakit.microbiome import sketch_dna, compare_signatures, plot_comparison_matrix

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())
sketch_dna("a.fa", work / "a.sig", scaled=1, name="sigA first", sourmash_path=tool, verbose=False)
sketch_dna("b.fa", work / "b.sig", scaled=1, name="sigB second", sourmash_path=tool, verbose=False)

compare_signatures([work / "a.sig", work / "b.sig"], work / "cmp", sourmash_path=tool, verbose=False)
result = plot_comparison_matrix(work / "cmp", work / "plots", sourmash_path=tool, verbose=False)
print(result["status"], Path(result["image"]).suffix, Path(result["image"]).stat().st_size > 1000)
```

- **Example output:**

```text
completed .png True
```

## 12) `dnakit.microbiome.prefetch_signatures()`

- **Purpose:** Select the database signatures that overlap a query signature.
- **Method:** Run `sourmash prefetch <query> <databases> -k <ksize> --dna --threshold-bp <bp> -o <csv>`. Prefetch is a fast first pass that lists every database signature sharing hashes with the query, without the greedy decomposition of `gather`.
- **API:** `dnakit.microbiome.prefetch_signatures(query: 'str', databases: 'Sequence[str]', output_csv: 'str', *, ksize: 'int' = 21, threshold_bp: 'int' = 0, sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `sourmash` on PATH (or its location in `sourmash_path`). It is provided by the optional package sourmash 4.9 or later (extra `dnakit[bio-sourmash]`); without it, `BackendUnavailableError` is raised.
- **Input:** `query` is an existing signature file, `databases` a non-empty list of existing signature collections. `ksize` defaults to 21 and `threshold_bp` to 0. `output_csv` must not exist yet.
- **Output:** The standard result dictionary: `status` (`completed` when sourmash exits with 0 and every expected output file exists, otherwise `failed`), `backend` (`sourmash`), `command`, `return_code`, `outputs` (the existing output files) and `log_tail` (the last 4000 characters of sourmash's output, which contain terminal control codes). The CSV has one row per overlapping database signature.
- **Example code:**

```python
import shutil, sys, tempfile
from pathlib import Path
import csv

from dnakit.microbiome import sketch_dna, prefetch_signatures

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())
sketch_dna("a.fa", work / "a.sig", scaled=1, name="sigA first", sourmash_path=tool, verbose=False)
sketch_dna("b.fa", work / "b.sig", scaled=1, name="sigB second", sourmash_path=tool, verbose=False)

result = prefetch_signatures(work / "a.sig", [work / "b.sig"], work / "prefetch.csv", sourmash_path=tool, verbose=False)
rows = list(csv.DictReader(open(work / "prefetch.csv")))
print(result["status"], len(rows), rows[0]["match_name"], rows[0]["intersect_bp"])
```

- **Example output:**

```text
completed 1 sigB second 15
```

## 13) `dnakit.microbiome.sketch_translated_dna()`

- **Purpose:** Translate DNA in all six frames and sketch the protein k-mers.
- **Method:** Run `sourmash sketch translate <fasta> -p k=<ksize>,scaled=<scaled> -o <signature>`. The signature holds a protein sketch; for `k=5` sourmash reports ksize 5 for the protein sketch.
- **API:** `dnakit.microbiome.sketch_translated_dna(fasta: 'str', output_signature: 'str', *, ksize: 'int' = 10, scaled: 'int' = 200, sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `sourmash` on PATH (or its location in `sourmash_path`). It is provided by the optional package sourmash 4.9 or later (extra `dnakit[bio-sourmash]`); without it, `BackendUnavailableError` is raised.
- **Input:** `fasta` is an existing FASTA file, `ksize` (default 10) and `scaled` (default 200) are positive integers; otherwise `ConfigurationError` is raised. `output_signature` must not exist yet.
- **Output:** The standard result dictionary: `status` (`completed` when sourmash exits with 0 and every expected output file exists, otherwise `failed`), `backend` (`sourmash`), `command`, `return_code`, `outputs` (the existing output files) and `log_tail` (the last 4000 characters of sourmash's output, which contain terminal control codes).
- **Example code:**

```python
import shutil, sys, tempfile
from pathlib import Path
from sourmash import load_file_as_signatures

from dnakit.microbiome import sketch_dna, sketch_translated_dna

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())

result = sketch_translated_dna("a.fa", work / "protein.sig", ksize=5, scaled=1, sourmash_path=tool, verbose=False)
(sig,) = load_file_as_signatures(str(work / "protein.sig"))
print(result["status"], sig.minhash.moltype, sig.minhash.ksize, len(sig.minhash))
```

- **Example output:**

```text
completed protein 5 92
```

## 14) `dnakit.microbiome.compare_signatures()`

- **Purpose:** Compare signatures pairwise and write the similarity (or containment) matrix.
- **Method:** Run `sourmash compare <signatures> -k <ksize> -o <matrix>`, adding `--containment` when `containment` is true. sourmash writes the matrix as a NumPy file (readable with `numpy.load`) and the sample order next to it as `<matrix>.labels.txt`.
- **API:** `dnakit.microbiome.compare_signatures(signatures: 'Sequence[str]', output_matrix: 'str', *, ksize: 'int' = 21, containment: 'bool' = False, sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `sourmash` on PATH (or its location in `sourmash_path`). It is provided by the optional package sourmash 4.9 or later (extra `dnakit[bio-sourmash]`); without it, `BackendUnavailableError` is raised.
- **Input:** `signatures` is a list of at least two existing signature files that contain a sketch of `ksize` (default 21); fewer than two raises `ConfigurationError`. `output_matrix` must not exist yet.
- **Output:** The standard result dictionary: `status` (`completed` when sourmash exits with 0 and every expected output file exists, otherwise `failed`), `backend` (`sourmash`), `command`, `return_code`, `outputs` (the existing output files) and `log_tail` (the last 4000 characters of sourmash's output, which contain terminal control codes). `outputs` contains the matrix file; the labels file is not listed.
- **Example code:**

```python
import shutil, sys, tempfile
from pathlib import Path
import numpy as np

from dnakit.microbiome import sketch_dna, compare_signatures

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())
sketch_dna("a.fa", work / "a.sig", scaled=1, name="sigA first", sourmash_path=tool, verbose=False)
sketch_dna("b.fa", work / "b.sig", scaled=1, name="sigB second", sourmash_path=tool, verbose=False)

result = compare_signatures([work / "a.sig", work / "b.sig"], work / "cmp", sourmash_path=tool, verbose=False)
matrix = np.load(work / "cmp")
print(result["status"], matrix.shape, matrix.round(3).tolist())
```

- **Example output:**

```text
completed (2, 2) [[1.0, 0.214], [0.214, 1.0]]
```

## 15) `dnakit.microbiome.concatenate_signatures()`

- **Purpose:** Concatenate signatures into one file without merging them.
- **Method:** Run `sourmash sig cat <signatures> -o <signature>`. Every signature stays a separate entry in the output.
- **API:** `dnakit.microbiome.concatenate_signatures(signatures: 'Sequence[str]', output_signature: 'str', *, sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `sourmash` on PATH (or its location in `sourmash_path`). It is provided by the optional package sourmash 4.9 or later (extra `dnakit[bio-sourmash]`); without it, `BackendUnavailableError` is raised.
- **Input:** `signatures` is a non-empty list of existing signature files. `output_signature` must not exist yet.
- **Output:** The standard result dictionary: `status` (`completed` when sourmash exits with 0 and every expected output file exists, otherwise `failed`), `backend` (`sourmash`), `command`, `return_code`, `outputs` (the existing output files) and `log_tail` (the last 4000 characters of sourmash's output, which contain terminal control codes).
- **Example code:**

```python
import shutil, sys, tempfile
from pathlib import Path
from sourmash import load_file_as_signatures

from dnakit.microbiome import sketch_dna, concatenate_signatures

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())
sketch_dna("a.fa", work / "a.sig", scaled=1, name="sigA first", sourmash_path=tool, verbose=False)
sketch_dna("b.fa", work / "b.sig", scaled=1, name="sigB second", sourmash_path=tool, verbose=False)

result = concatenate_signatures([work / "a.sig", work / "b.sig"], work / "both.sig", sourmash_path=tool, verbose=False)
sigs = load_file_as_signatures(str(work / "both.sig"))
print(result["status"], [sig.name for sig in sigs])
```

- **Example output:**

```text
completed ['sigA first', 'sigB second']
```

## 16) `dnakit.microbiome.extract_signatures()`

- **Purpose:** Extract the signatures with one k-mer size from a collection.
- **Method:** Run `sourmash sig extract <collection> -k <ksize> -o <signature>` with `--dna` (default) or `--protein` when `dna` is false.
- **API:** `dnakit.microbiome.extract_signatures(collection: 'Sequence[str]', output_signature: 'str', *, ksize: 'int', dna: 'bool' = True, sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `sourmash` on PATH (or its location in `sourmash_path`). It is provided by the optional package sourmash 4.9 or later (extra `dnakit[bio-sourmash]`); without it, `BackendUnavailableError` is raised.
- **Input:** `collection` is a non-empty list of existing signature files, `ksize` is required (keyword only), `dna` chooses nucleotide or protein sketches. `output_signature` must not exist yet.
- **Output:** The standard result dictionary: `status` (`completed` when sourmash exits with 0 and every expected output file exists, otherwise `failed`), `backend` (`sourmash`), `command`, `return_code`, `outputs` (the existing output files) and `log_tail` (the last 4000 characters of sourmash's output, which contain terminal control codes).
- **Example code:**

```python
import shutil, sys, tempfile
from pathlib import Path
from sourmash import load_file_as_signatures

from dnakit.microbiome import sketch_dna, extract_signatures

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())
sketch_dna("a.fa", work / "a.sig", scaled=1, name="sigA first", sourmash_path=tool, verbose=False)
sketch_dna("b.fa", work / "b.sig", scaled=1, name="sigB second", sourmash_path=tool, verbose=False)

sketch_dna("b.fa", work / "multi.sig", ksizes=[21, 31], scaled=1, sourmash_path=tool, verbose=False)
result = extract_signatures([work / "multi.sig"], work / "k31.sig", ksize=31, sourmash_path=tool, verbose=False)
(sig,) = load_file_as_signatures(str(work / "k31.sig"))
print(result["status"], sig.minhash.ksize, len(sig.minhash))
```

- **Example output:**

```text
completed 31 35
```

## 17) `dnakit.microbiome.filter_signature_abundance()`

- **Purpose:** Keep the hashes whose abundance lies in an inclusive range.
- **Method:** Run `sourmash sig filter <signature> -o <output>` with `-m <min_abundance>` and/or `-M <max_abundance>`. It only makes sense for signatures that track abundance; the facade `sketch_dna()` makes signatures without abundance, so the example builds an abundance-tracking signature with the sourmash Python API.
- **API:** `dnakit.microbiome.filter_signature_abundance(signature: 'str', output_signature: 'str', *, min_abundance: 'int | None' = None, max_abundance: 'int | None' = None, sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `sourmash` on PATH (or its location in `sourmash_path`). It is provided by the optional package sourmash 4.9 or later (extra `dnakit[bio-sourmash]`); without it, `BackendUnavailableError` is raised.
- **Input:** `signature` is an existing signature file; at least one of `min_abundance` and `max_abundance` must be given, otherwise `ConfigurationError` is raised. `output_signature` must not exist yet.
- **Output:** The standard result dictionary: `status` (`completed` when sourmash exits with 0 and every expected output file exists, otherwise `failed`), `backend` (`sourmash`), `command`, `return_code`, `outputs` (the existing output files) and `log_tail` (the last 4000 characters of sourmash's output, which contain terminal control codes).
- **Example code:**

```python
import shutil, sys, tempfile
from pathlib import Path
from sourmash import MinHash, SourmashSignature, load_file_as_signatures
from sourmash.sourmash_args import SaveSignaturesToLocation

from dnakit.microbiome import sketch_dna, filter_signature_abundance

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())

minhash = MinHash(n=0, ksize=21, scaled=1, track_abundance=True)
minhash.add_sequence("GCTAAAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCG")
minhash.add_sequence("GCTAAAGACAATTACATAACATACACGTCAGCACGA")
with SaveSignaturesToLocation(str(work / "abund.sig")) as saver:
    saver.add(SourmashSignature(minhash, name="abund"))
result = filter_signature_abundance(work / "abund.sig", work / "kept.sig", min_abundance=2, sourmash_path=tool, verbose=False)
(sig,) = load_file_as_signatures(str(work / "kept.sig"))
print(result["status"], len(minhash), len(sig.minhash), sorted(set(sig.minhash.hashes.values())))
```

- **Example output:**

```text
completed 40 16 [2]
```

## 18) `dnakit.microbiome.manifest_signature()`

- **Purpose:** Write a CSV manifest describing every signature in a file.
- **Method:** Run `sourmash sig manifest <signature> -o <csv>`. The CSV starts with a `# SOURMASH-MANIFEST-VERSION` comment line, then a header with columns such as `md5`, `ksize`, `moltype`, `num`, `scaled`, `n_hashes`, `with_abundance` and `name`.
- **API:** `dnakit.microbiome.manifest_signature(signature: 'str', output_csv: 'str', *, sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `sourmash` on PATH (or its location in `sourmash_path`). It is provided by the optional package sourmash 4.9 or later (extra `dnakit[bio-sourmash]`); without it, `BackendUnavailableError` is raised.
- **Input:** `signature` is an existing signature file. `output_csv` must not exist yet.
- **Output:** The standard result dictionary: `status` (`completed` when sourmash exits with 0 and every expected output file exists, otherwise `failed`), `backend` (`sourmash`), `command`, `return_code`, `outputs` (the existing output files) and `log_tail` (the last 4000 characters of sourmash's output, which contain terminal control codes).
- **Example code:**

```python
import shutil, sys, tempfile
from pathlib import Path
import csv

from dnakit.microbiome import sketch_dna, manifest_signature

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())
sketch_dna("a.fa", work / "a.sig", scaled=1, name="sigA first", sourmash_path=tool, verbose=False)
sketch_dna("b.fa", work / "b.sig", scaled=1, name="sigB second", sourmash_path=tool, verbose=False)

result = manifest_signature(work / "a.sig", work / "manifest.csv", sourmash_path=tool, verbose=False)
lines = (work / "manifest.csv").read_text().splitlines()
rows = list(csv.DictReader(lines[1:]))
print(result["status"], len(rows), rows[0]["name"], rows[0]["ksize"], rows[0]["n_hashes"])
```

- **Example output:**

```text
completed 1 sigA first 21 40
```

## 19) `dnakit.microbiome.split_signatures()`

- **Purpose:** Split signature files into one file per signature with standardised names.
- **Method:** Run `sourmash sig split <signatures> --outdir <dir>`. DNAKit creates the folder if needed (it must be new or empty), and `outputs` is the sorted list of `*.sig` files in it; `status` is `completed` only when the exit code is 0 and at least one `.sig` file exists.
- **API:** `dnakit.microbiome.split_signatures(signatures: 'Sequence[str]', output_directory: 'str', *, sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `sourmash` on PATH (or its location in `sourmash_path`). It is provided by the optional package sourmash 4.9 or later (extra `dnakit[bio-sourmash]`); without it, `BackendUnavailableError` is raised.
- **Input:** `signatures` is a non-empty list of existing signature files. `output_directory` must not exist or must be empty, otherwise `ConfigurationError` is raised.
- **Output:** The standard result dictionary: `status` (`completed` when sourmash exits with 0 and every expected output file exists, otherwise `failed`), `backend` (`sourmash`), `command`, `return_code`, `outputs` (the existing output files) and `log_tail` (the last 4000 characters of sourmash's output, which contain terminal control codes). Here `outputs` lists every `.sig` file written to the folder.
- **Example code:**

```python
import shutil, sys, tempfile
from pathlib import Path

from dnakit.microbiome import sketch_dna, split_signatures

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())
sketch_dna("a.fa", work / "a.sig", scaled=1, name="sigA first", sourmash_path=tool, verbose=False)
sketch_dna("b.fa", work / "b.sig", scaled=1, name="sigB second", sourmash_path=tool, verbose=False)

result = split_signatures([work / "a.sig", work / "b.sig"], work / "parts", sourmash_path=tool, verbose=False)
print(result["status"], len(result["outputs"]), all(name.endswith(".sig") for name in result["outputs"]))
```

- **Example output:**

```text
completed 2 True
```

## 20) `dnakit.microbiome.downsample_signature()`

- **Purpose:** Reduce the resolution of a signature to a larger scaled value.
- **Method:** Run `sourmash sig downsample <signature> --scaled <n> -o <output>`. Only hashes below `2^64 / n` remain, so a signature can be downsampled but never made finer.
- **API:** `dnakit.microbiome.downsample_signature(signature: 'str', output_signature: 'str', *, scaled: 'int', sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `sourmash` on PATH (or its location in `sourmash_path`). It is provided by the optional package sourmash 4.9 or later (extra `dnakit[bio-sourmash]`); without it, `BackendUnavailableError` is raised.
- **Input:** `signature` is an existing signature file, `scaled` (required, keyword only) a positive integer larger than the signature's own value. `output_signature` must not exist yet.
- **Output:** The standard result dictionary: `status` (`completed` when sourmash exits with 0 and every expected output file exists, otherwise `failed`), `backend` (`sourmash`), `command`, `return_code`, `outputs` (the existing output files) and `log_tail` (the last 4000 characters of sourmash's output, which contain terminal control codes).
- **Example code:**

```python
import shutil, sys, tempfile
from pathlib import Path
from sourmash import load_file_as_signatures

from dnakit.microbiome import sketch_dna, downsample_signature

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())
sketch_dna("a.fa", work / "a.sig", scaled=1, name="sigA first", sourmash_path=tool, verbose=False)
sketch_dna("b.fa", work / "b.sig", scaled=1, name="sigB second", sourmash_path=tool, verbose=False)

result = downsample_signature(work / "a.sig", work / "a_scaled4.sig", scaled=4, sourmash_path=tool, verbose=False)
(sig,) = load_file_as_signatures(str(work / "a_scaled4.sig"))
print(result["status"], sig.minhash.scaled, len(sig.minhash))
```

- **Example output:**

```text
completed 4 10
```

## 21) `dnakit.microbiome.flatten_signature()`

- **Purpose:** Remove abundance information from a signature.
- **Method:** Run `sourmash sig flatten <signature> -o <output>`. The hashes stay; their counts are dropped.
- **API:** `dnakit.microbiome.flatten_signature(signature: 'str', output_signature: 'str', *, sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `sourmash` on PATH (or its location in `sourmash_path`). It is provided by the optional package sourmash 4.9 or later (extra `dnakit[bio-sourmash]`); without it, `BackendUnavailableError` is raised.
- **Input:** `signature` is an existing signature file. `output_signature` must not exist yet.
- **Output:** The standard result dictionary: `status` (`completed` when sourmash exits with 0 and every expected output file exists, otherwise `failed`), `backend` (`sourmash`), `command`, `return_code`, `outputs` (the existing output files) and `log_tail` (the last 4000 characters of sourmash's output, which contain terminal control codes).
- **Example code:**

```python
import shutil, sys, tempfile
from pathlib import Path
from sourmash import MinHash, SourmashSignature, load_file_as_signatures
from sourmash.sourmash_args import SaveSignaturesToLocation

from dnakit.microbiome import sketch_dna, flatten_signature

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())

minhash = MinHash(n=0, ksize=21, scaled=1, track_abundance=True)
minhash.add_sequence("GCTAAAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCG")
minhash.add_sequence("GCTAAAGACAATTACATAACATACACGTCAGCACGA")
with SaveSignaturesToLocation(str(work / "abund.sig")) as saver:
    saver.add(SourmashSignature(minhash, name="abund"))
result = flatten_signature(work / "abund.sig", work / "flat.sig", sourmash_path=tool, verbose=False)
(sig,) = load_file_as_signatures(str(work / "flat.sig"))
print(result["status"], minhash.track_abundance, sig.minhash.track_abundance, len(sig.minhash))
```

- **Example output:**

```text
completed True False 40
```

## 22) `dnakit.microbiome.rename_signature()`

- **Purpose:** Give a signature a new name.
- **Method:** Run `sourmash sig rename <signature> <new name> -o <output>`.
- **API:** `dnakit.microbiome.rename_signature(signature: 'str', new_name: 'str', output_signature: 'str', *, sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `sourmash` on PATH (or its location in `sourmash_path`). It is provided by the optional package sourmash 4.9 or later (extra `dnakit[bio-sourmash]`); without it, `BackendUnavailableError` is raised.
- **Input:** `signature` is an existing signature file, `new_name` a non-empty string (blank raises `ConfigurationError`). `output_signature` must not exist yet.
- **Output:** The standard result dictionary: `status` (`completed` when sourmash exits with 0 and every expected output file exists, otherwise `failed`), `backend` (`sourmash`), `command`, `return_code`, `outputs` (the existing output files) and `log_tail` (the last 4000 characters of sourmash's output, which contain terminal control codes).
- **Example code:**

```python
import shutil, sys, tempfile
from pathlib import Path
from sourmash import load_file_as_signatures

from dnakit.microbiome import sketch_dna, rename_signature

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())
sketch_dna("a.fa", work / "a.sig", scaled=1, name="sigA first", sourmash_path=tool, verbose=False)
sketch_dna("b.fa", work / "b.sig", scaled=1, name="sigB second", sourmash_path=tool, verbose=False)

result = rename_signature(work / "a.sig", "renamed sample", work / "renamed.sig", sourmash_path=tool, verbose=False)
(sig,) = load_file_as_signatures(str(work / "renamed.sig"))
print(result["status"], sig.name)
```

- **Example output:**

```text
completed renamed sample
```

## 23) `dnakit.microbiome.signature_file_info()`

- **Purpose:** Report file-level information about a signature file.
- **Method:** Run `sourmash sig fileinfo <signature>` and return its text: the file type, whether it has a manifest, the number of signatures, the total hashes and a summary per k-mer size and molecule type.
- **API:** `dnakit.microbiome.signature_file_info(signature: 'str', *, sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `sourmash` on PATH (or its location in `sourmash_path`). It is provided by the optional package sourmash 4.9 or later (extra `dnakit[bio-sourmash]`); without it, `BackendUnavailableError` is raised.
- **Input:** `signature` is an existing signature file.
- **Output:** `status` (`completed` when the exit code is 0), `backend`, `return_code` and `description`, the last 4000 characters of sourmash's text output (which may contain terminal control codes).
- **Example code:**

```python
import shutil, sys, tempfile
from pathlib import Path

from dnakit.microbiome import sketch_dna, signature_file_info

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())
sketch_dna("a.fa", work / "a.sig", scaled=1, name="sigA first", sourmash_path=tool, verbose=False)
sketch_dna("b.fa", work / "b.sig", scaled=1, name="sigB second", sourmash_path=tool, verbose=False)

result = signature_file_info(work / "a.sig", sourmash_path=tool, verbose=False)
print(result["status"])
print([line.strip() for line in result["description"].splitlines() if "total hashes" in line or "sketches with" in line])
```

- **Example output:**

```text
completed
['total hashes: 40', '1 sketches with DNA, k=21, scaled=1                40 total hashes']
```

## 24) `dnakit.microbiome.intersect_signatures()`

- **Purpose:** Keep the k-mers shared by all signatures.
- **Method:** Run `sourmash sig intersect <signatures> -o <output>`.
- **API:** `dnakit.microbiome.intersect_signatures(signatures: 'Sequence[str]', output_signature: 'str', *, sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `sourmash` on PATH (or its location in `sourmash_path`). It is provided by the optional package sourmash 4.9 or later (extra `dnakit[bio-sourmash]`); without it, `BackendUnavailableError` is raised.
- **Input:** `signatures` is a list of at least two existing signature files with the same k-mer size and scaled value; fewer than two raises `ConfigurationError`. `output_signature` must not exist yet.
- **Output:** The standard result dictionary: `status` (`completed` when sourmash exits with 0 and every expected output file exists, otherwise `failed`), `backend` (`sourmash`), `command`, `return_code`, `outputs` (the existing output files) and `log_tail` (the last 4000 characters of sourmash's output, which contain terminal control codes).
- **Example code:**

```python
import shutil, sys, tempfile
from pathlib import Path
from sourmash import load_file_as_signatures

from dnakit.microbiome import sketch_dna, intersect_signatures

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())
sketch_dna("a.fa", work / "a.sig", scaled=1, name="sigA first", sourmash_path=tool, verbose=False)
sketch_dna("b.fa", work / "b.sig", scaled=1, name="sigB second", sourmash_path=tool, verbose=False)

result = intersect_signatures([work / "a.sig", work / "b.sig"], work / "shared.sig", sourmash_path=tool, verbose=False)
(sig,) = load_file_as_signatures(str(work / "shared.sig"))
print(result["status"], len(sig.minhash))
```

- **Example output:**

```text
completed 15
```

## 25) `dnakit.microbiome.merge_signatures()`

- **Purpose:** Merge several signatures of the same k-mer size into one.
- **Method:** Run `sourmash sig merge <signatures> -o <output>`, adding `--flatten` when `flatten` is true (removes all abundances so that a mix of abundance and non-abundance signatures can be merged). The result holds the union of the hashes.
- **API:** `dnakit.microbiome.merge_signatures(signatures: 'Sequence[str]', output_signature: 'str', *, flatten: 'bool' = False, sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `sourmash` on PATH (or its location in `sourmash_path`). It is provided by the optional package sourmash 4.9 or later (extra `dnakit[bio-sourmash]`); without it, `BackendUnavailableError` is raised.
- **Input:** `signatures` is a list of at least two existing signature files with the same k-mer size; fewer than two raises `ConfigurationError`. `output_signature` must not exist yet.
- **Output:** The standard result dictionary: `status` (`completed` when sourmash exits with 0 and every expected output file exists, otherwise `failed`), `backend` (`sourmash`), `command`, `return_code`, `outputs` (the existing output files) and `log_tail` (the last 4000 characters of sourmash's output, which contain terminal control codes).
- **Example code:**

```python
import shutil, sys, tempfile
from pathlib import Path
from sourmash import load_file_as_signatures

from dnakit.microbiome import sketch_dna, merge_signatures

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())
sketch_dna("a.fa", work / "a.sig", scaled=1, name="sigA first", sourmash_path=tool, verbose=False)
sketch_dna("b.fa", work / "b.sig", scaled=1, name="sigB second", sourmash_path=tool, verbose=False)

result = merge_signatures([work / "a.sig", work / "b.sig"], work / "merged.sig", sourmash_path=tool, verbose=False)
(sig,) = load_file_as_signatures(str(work / "merged.sig"))
print(result["status"], len(sig.minhash))
```

- **Example output:**

```text
completed 70
```

## 26) `dnakit.microbiome.subtract_signatures()`

- **Purpose:** Remove the k-mers of other signatures from one signature.
- **Method:** Run `sourmash sig subtract <signature> <others> -o <output>`.
- **API:** `dnakit.microbiome.subtract_signatures(signature: 'str', others: 'Sequence[str]', output_signature: 'str', *, sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `sourmash` on PATH (or its location in `sourmash_path`). It is provided by the optional package sourmash 4.9 or later (extra `dnakit[bio-sourmash]`); without it, `BackendUnavailableError` is raised.
- **Input:** `signature` is an existing signature file, `others` a non-empty list of existing signature files (empty raises `ConfigurationError`). `output_signature` must not exist yet.
- **Output:** The standard result dictionary: `status` (`completed` when sourmash exits with 0 and every expected output file exists, otherwise `failed`), `backend` (`sourmash`), `command`, `return_code`, `outputs` (the existing output files) and `log_tail` (the last 4000 characters of sourmash's output, which contain terminal control codes).
- **Example code:**

```python
import shutil, sys, tempfile
from pathlib import Path
from sourmash import load_file_as_signatures

from dnakit.microbiome import sketch_dna, subtract_signatures

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())
sketch_dna("a.fa", work / "a.sig", scaled=1, name="sigA first", sourmash_path=tool, verbose=False)
sketch_dna("b.fa", work / "b.sig", scaled=1, name="sigB second", sourmash_path=tool, verbose=False)

result = subtract_signatures(work / "a.sig", [work / "b.sig"], work / "only_a.sig", sourmash_path=tool, verbose=False)
(sig,) = load_file_as_signatures(str(work / "only_a.sig"))
print(result["status"], len(sig.minhash))
```

- **Example output:**

```text
completed 25
```

## 27) `dnakit.microbiome.describe_signature()`

- **Purpose:** Describe the contents of a signature file: name, k-mer size, scaled value and size.
- **Method:** Run `sourmash sig describe <signature>` and return its text. For each signature it shows the name, source file, md5, a line with `k=`, molecule type, `num`, `scaled`, seed and abundance flag, and the number of hashes.
- **API:** `dnakit.microbiome.describe_signature(signature: 'str', *, sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `sourmash` on PATH (or its location in `sourmash_path`). It is provided by the optional package sourmash 4.9 or later (extra `dnakit[bio-sourmash]`); without it, `BackendUnavailableError` is raised.
- **Input:** `signature` is an existing signature file.
- **Output:** `status` (`completed` when the exit code is 0), `backend`, `return_code` and `description`, the last 4000 characters of sourmash's text output (which may contain terminal control codes).
- **Example code:**

```python
import shutil, sys, tempfile
from pathlib import Path

from dnakit.microbiome import sketch_dna, describe_signature

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())
sketch_dna("a.fa", work / "a.sig", scaled=1, name="sigA first", sourmash_path=tool, verbose=False)
sketch_dna("b.fa", work / "b.sig", scaled=1, name="sigB second", sourmash_path=tool, verbose=False)

result = describe_signature(work / "a.sig", sourmash_path=tool, verbose=False)
print(result["status"])
print([line for line in result["description"].splitlines() if line.startswith(("signature:", "k=", "size:"))])
```

- **Example output:**

```text
completed
['signature: sigA first', 'k=21 molecule=DNA num=0 scaled=1 seed=42 track_abundance=0', 'size: 40']
```

## 28) `dnakit.microbiome.gather_signatures()`

- **Purpose:** Decompose a query signature into the best database matches.
- **Method:** Run `sourmash gather <query> <database> -k <ksize> -o <csv> --threshold-bp <bp>`. Gather repeatedly picks the database signature with the largest remaining overlap and removes its hashes from the query, so the rows are non-overlapping matches. `threshold_bp` is the minimum overlap reported; the default 50000 is sourmash's default; with a tiny test signature no match reaches it, sourmash then writes no CSV, and `status` is `failed`. The example therefore sets 0.
- **API:** `dnakit.microbiome.gather_signatures(query: 'str', database: 'str', output_csv: 'str', *, ksize: 'int' = 21, threshold_bp: 'int' = 50000, sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `sourmash` on PATH (or its location in `sourmash_path`). It is provided by the optional package sourmash 4.9 or later (extra `dnakit[bio-sourmash]`); without it, `BackendUnavailableError` is raised.
- **Input:** `query` is an existing signature file, `database` an existing signature collection, `ksize` defaults to 21, `threshold_bp` must not be negative (`ConfigurationError` otherwise). `output_csv` must not exist yet.
- **Output:** The standard result dictionary: `status` (`completed` when sourmash exits with 0 and every expected output file exists, otherwise `failed`), `backend` (`sourmash`), `command`, `return_code`, `outputs` (the existing output files) and `log_tail` (the last 4000 characters of sourmash's output, which contain terminal control codes). The CSV has one row per gather match, with columns such as `intersect_bp`, `f_orig_query`, `f_match` and `name`.
- **Example code:**

```python
import shutil, sys, tempfile
from pathlib import Path
import csv

from dnakit.microbiome import sketch_dna, gather_signatures

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())
sketch_dna("a.fa", work / "a.sig", scaled=1, name="sigA first", sourmash_path=tool, verbose=False)
sketch_dna("b.fa", work / "b.sig", scaled=1, name="sigB second", sourmash_path=tool, verbose=False)

result = gather_signatures(work / "a.sig", work / "b.sig", work / "gather.csv", threshold_bp=0, sourmash_path=tool, verbose=False)
rows = list(csv.DictReader(open(work / "gather.csv")))
print(result["status"], len(rows), rows[0]["name"], rows[0]["intersect_bp"], round(float(rows[0]["f_orig_query"]), 3))
```

- **Example output:**

```text
completed 1 sigB second 15 0.375
```

## 29) `dnakit.microbiome.index_signatures()`

- **Purpose:** Build a searchable signature database (SBT) from signature files.
- **Method:** Run `sourmash index <database> <signatures> -k <ksize>`. sourmash appends `.sbt.zip` to the name, and DNAKit therefore checks and reports `<output_database>.sbt.zip`.
- **API:** `dnakit.microbiome.index_signatures(signatures: 'Sequence[str]', output_database: 'str', *, ksize: 'int' = 21, sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `sourmash` on PATH (or its location in `sourmash_path`). It is provided by the optional package sourmash 4.9 or later (extra `dnakit[bio-sourmash]`); without it, `BackendUnavailableError` is raised.
- **Input:** `signatures` is a non-empty list of existing signature files, `ksize` defaults to 21. Neither `output_database` nor `<output_database>.sbt.zip` may exist.
- **Output:** The standard result dictionary: `status` (`completed` when sourmash exits with 0 and every expected output file exists, otherwise `failed`), `backend` (`sourmash`), `command`, `return_code`, `outputs` (the existing output files) and `log_tail` (the last 4000 characters of sourmash's output, which contain terminal control codes). `database` is the path of the `.sbt.zip` file.
- **Example code:**

```python
import shutil, sys, tempfile
from pathlib import Path

from dnakit.microbiome import sketch_dna, index_signatures

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())
sketch_dna("a.fa", work / "a.sig", scaled=1, name="sigA first", sourmash_path=tool, verbose=False)
sketch_dna("b.fa", work / "b.sig", scaled=1, name="sigB second", sourmash_path=tool, verbose=False)

result = index_signatures([work / "a.sig", work / "b.sig"], work / "db", sourmash_path=tool, verbose=False)
print(result["status"], Path(result["database"]).name, len(result["outputs"]))
```

- **Example output:**

```text
completed db.sbt.zip 1
```

## 30) `dnakit.microbiome.search_signatures()`

- **Purpose:** Search a signature database for matches to a query signature.
- **Method:** Run `sourmash search <query> <database> -k <ksize> --threshold <t> -o <csv>`, adding `--containment` when `containment` is true. Matches are ranked by Jaccard similarity, or by containment of the query when requested.
- **API:** `dnakit.microbiome.search_signatures(query: 'str', database: 'str', output_csv: 'str', *, ksize: 'int' = 21, threshold: 'float' = 0.08, containment: 'bool' = False, sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `sourmash` on PATH (or its location in `sourmash_path`). It is provided by the optional package sourmash 4.9 or later (extra `dnakit[bio-sourmash]`); without it, `BackendUnavailableError` is raised.
- **Input:** `query` is an existing signature file, `database` an existing signature collection, `ksize` defaults to 21, `threshold` is between 0 and 1 (default 0.08; otherwise `ConfigurationError`). `output_csv` must not exist yet.
- **Output:** The standard result dictionary: `status` (`completed` when sourmash exits with 0 and every expected output file exists, otherwise `failed`), `backend` (`sourmash`), `command`, `return_code`, `outputs` (the existing output files) and `log_tail` (the last 4000 characters of sourmash's output, which contain terminal control codes). The CSV has one row per match, with columns such as `similarity`, `md5`, `name` and `filename`.
- **Example code:**

```python
import shutil, sys, tempfile
from pathlib import Path
import csv

from dnakit.microbiome import sketch_dna, search_signatures

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())
sketch_dna("a.fa", work / "a.sig", scaled=1, name="sigA first", sourmash_path=tool, verbose=False)
sketch_dna("b.fa", work / "b.sig", scaled=1, name="sigB second", sourmash_path=tool, verbose=False)

result = search_signatures(work / "a.sig", work / "b.sig", work / "search.csv", threshold=0.1, sourmash_path=tool, verbose=False)
rows = list(csv.DictReader(open(work / "search.csv")))
print(result["status"], len(rows), rows[0]["name"], round(float(rows[0]["similarity"]), 3))
```

- **Example output:**

```text
completed 1 sigB second 0.214
```

## 31) `dnakit.microbiome.sketch_dna()`

- **Purpose:** Compute a DNA MinHash signature of a FASTA file with scaled sketching.
- **Method:** Run `sourmash sketch dna <fasta> -o <signature> -p k=<k1>,k=<k2>,...,scaled=<n>`, adding `--name <name>` when `name` is given. All records of the file go into one signature per k-mer size, and abundance is not tracked.
- **API:** `dnakit.microbiome.sketch_dna(fasta: 'str', output_signature: 'str', *, ksizes: 'Sequence[int]' = (21,), scaled: 'int' = 1000, name: 'str | None' = None, sourmash_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `sourmash` on PATH (or its location in `sourmash_path`). It is provided by the optional package sourmash 4.9 or later (extra `dnakit[bio-sourmash]`); without it, `BackendUnavailableError` is raised.
- **Input:** `fasta` is an existing FASTA file; `ksizes` is a non-empty sequence of positive integers (default `(21,)`), `scaled` a positive integer (default 1000), otherwise `ConfigurationError`. sourmash warns for `scaled` below 100 but continues, which is why the examples use 1. `output_signature` must not exist yet.
- **Output:** The standard result dictionary: `status` (`completed` when sourmash exits with 0 and every expected output file exists, otherwise `failed`), `backend` (`sourmash`), `command`, `return_code`, `outputs` (the existing output files) and `log_tail` (the last 4000 characters of sourmash's output, which contain terminal control codes).
- **Example code:**

```python
import shutil, sys, tempfile
from pathlib import Path
from sourmash import load_file_as_signatures

from dnakit.microbiome import sketch_dna

tool = shutil.which("sourmash") or str(Path(sys.executable).with_name("sourmash"))
work = Path(tempfile.mkdtemp())

result = sketch_dna("a.fa", work / "a.sig", scaled=1, name="sigA first", sourmash_path=tool, verbose=False)
(sig,) = load_file_as_signatures(str(work / "a.sig"))
print(result["status"], result["return_code"], len(sig.minhash), sig.name)
```

- **Example output:**

```text
completed 0 40 sigA first
```
