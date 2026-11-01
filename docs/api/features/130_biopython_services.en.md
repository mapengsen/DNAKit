# Biopython: NCBI, KEGG and TogoWS services, trees and plots

This group of functions in `dnakit.sequence` is a set of adapters over [Biopython](https://biopython.org/) (package `Bio`, [source](https://github.com/biopython/biopython)). Nineteen of them query remote services and return the reply as a JSON-friendly dictionary: BLAST searches through NCBI QBLAST, the NCBI Entrez E-utilities (`einfo`, `esearch`, `efetch`, `esummary`, `elink`, `epost`), the KEGG REST service and the TogoWS service. The other three work locally: they draw a Newick tree as ASCII art or as a PNG with `Bio.Phylo`, and compute and plot the GC skew of a DNA sequence with `Bio.SeqUtils`. The module only wraps these upstream functions; it does not mirror any database and it does not interpret the biology of a reply.

Each network function checks its arguments and then calls one upstream function when it is invoked, so nothing is sent to a remote service before that. Sizes are bounded (at most 200 NCBI identifiers per request, 10 KEGG entries for `kegg_get_entries()`, 100 KEGG identifiers per conversion or link request and 1000 TogoWS identifiers per iterator call). An NCBI e-mail address is optional: it is set on the Biopython module only while the call runs and is not written to the result. Output files must be new paths; an existing file is refused. Errors raised by the remote service (unknown identifiers, HTTP errors, NCBI error replies) propagate unchanged. Check the terms of use of NCBI, KEGG and TogoWS before you send requests.

Biopython is released under the Biopython License Agreement (some files are dual licensed with the BSD 3-Clause licence), checked in the installed package metadata and its `LICENSE.rst`. It is a permissive licence, not the GPL. The dependencies are optional and are imported only when a function is called; a missing package raises an error with the install command. The examples import from `dnakit.sequence`.

## Example data

No input file is needed. The three offline examples pass their inputs directly and write their files to a temporary directory. The examples of the network functions call the function with an argument that is rejected by its validation, so that they run without internet access and print the real error; each of those sections also shows an unexecuted call as it would be made with network access.

## 1) `dnakit.sequence.blast_qblast_search()`

- **Purpose:** Run a BLAST search on the NCBI servers and summarise the hits.
- **Method:** Validate the arguments, then call `Bio.Blast.qblast()` with `format_type="XML"` (the module-level function of the newer `Bio.Blast` API), read the XML completely, optionally save it to `output_xml`, and parse it with the legacy `Bio.Blast.NCBIXML.read()` into a plain summary.
- **API:** `blast_qblast_search(program: 'str', database: 'str', sequence: 'str', *, hitlist_size: 'int' = 50, expect: 'float' = 10.0, entrez_query: 'str' = '(none)', output_xml: 'str | None' = None, email: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package Biopython (extra `dnakit[bio-formats]`) and internet access to NCBI BLAST (QBLAST). No key is required and an e-mail address is optional. A search is queued on the NCBI server and can take minutes. The example below passes an argument that is rejected before any request is sent, so it runs offline.
- **Input:** `program` is one of `blastn`, `blastp`, `blastx`, `tblastn`, `tblastx`. `database` is a short name such as `nt`. `sequence` is a non-empty string of at most 100000 characters (it is not checked for DNA letters, although the docstring says DNA query). `hitlist_size` is an integer from 1 to 500, `expect` is a positive number, `entrez_query` is at most 512 characters, `output_xml` is an optional new file path and `email` an optional address.
- **Output:** Every network function returns a dictionary with `status` (`completed`), `service`, `upstream` (the Biopython function that was called), `parameters` (the validated request parameters), `retrieved_utc`, `biopython_version`, `output_path` (a path or `None`) and `data`. `parameters` holds the query length instead of the query. `data` has `query_id`, `query`, `query_letters`, `database`, `alignment_count` and `alignments`; each alignment has `hit_id`, `hit_def`, `length` and `hsps` (`score`, `bits`, `expect`, `identities`, `align_length`, `query_start`, `query_end`, `subject_start`, `subject_end`). When `output_xml` is given the raw XML is also written there.
- **Example code:**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.sequence import blast_qblast_search

try:
    blast_qblast_search("blastz", "nt", "ACGT", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error)
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] program must be one of blastn, blastp, blastx, tblastn, tblastx.
```

Usage with the dependency installed:

```python
from dnakit.sequence import blast_qblast_search

result = blast_qblast_search(
    "blastn", "nt", "ATGGCGGCGCTGAGCGGTGGCGGTGGTGGCGGCGCGGAGCCGGGCCAGGCTCTGTTCAACGGGGACATGGAGCCCGAGGCCGGC",
    hitlist_size=5, email="you@example.org",
)
print(result["data"]["alignment_count"])
```

## 2) `dnakit.sequence.gc_skew_window_plot()`

- **Purpose:** Compute the GC skew of a DNA sequence in fixed windows and save a line plot of the skew and its running sum as a new PNG.
- **Method:** Check that the sequence holds only IUPAC DNA letters, compute the per-window values with `Bio.SeqUtils.GC_skew()` as (G - C) / (G + C), accumulate them with `itertools.accumulate`, and draw both lines with the matplotlib object API (`Figure` and `FigureCanvasAgg`, no pyplot state). It does not call `Bio.SeqUtils.xGC_skew()`, which needs a Tk display. A final partial window is dropped by `GC_skew()`.
- **API:** `gc_skew_window_plot(sequence: 'str', image_path: 'str', *, window: 'int' = 1000, dpi: 'int' = 150, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages Biopython (extra `dnakit[bio-formats]`) and matplotlib (extra `dnakit[bio-plots]`). It works offline and needs no display.
- **Input:** `sequence` is a non-empty string of IUPAC DNA letters, at most 10000000 characters. `image_path` must be a new path ending in `.png` in an existing folder. `window` is an integer from 1 to 10000000 (default 1000) and `dpi` an integer from 50 to 600.
- **Output:** A dictionary with `status` (`completed`), `image_path`, `format` (`png`) and `data`. `data` has `sequence_length`, `window`, `window_count`, `window_starts`, `gc_skew` and `accumulated_gc_skew`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import gc_skew_window_plot

image = Path(tempfile.mkdtemp()) / "skew.png"
result = gc_skew_window_plot("GGGGGGCCAAGGGCCCCCCTT", image, window=7, verbose=False)
data = result["data"]
print(result["status"], result["format"], image.read_bytes()[:4])
print(data["sequence_length"], data["window_count"], data["window_starts"])
print([round(value, 3) for value in data["gc_skew"]])
print([round(value, 3) for value in data["accumulated_gc_skew"]])
```

- **Example output:**

```text
completed png b'\x89PNG'
21 3 [0, 7, 14]
[0.714, 0.2, -1.0]
[0.714, 0.914, -0.086]
```

## 3) `dnakit.sequence.kegg_convert_identifiers()`

- **Purpose:** Convert KEGG identifiers to or from identifiers of an outside database.
- **Method:** Validate the arguments and call `Bio.KEGG.REST.kegg_conv()`. The tab-separated reply is split into rows.
- **API:** `kegg_convert_identifiers(target_db: 'str', source: 'str | Sequence[str]', *, option: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package Biopython (extra `dnakit[bio-formats]`) and internet access to the KEGG REST service (rest.kegg.jp). No key is required. The example below passes an argument that is rejected before any request is sent, so it runs offline; to run the function for real you need the network.
- **Input:** `target_db` is a short name such as `ncbi-geneid` or `hsa`. `source` is a database name or one identifier string, or a list of 1 to 100 identifiers (letters, digits and `_ . : | + -` only). `option` is `None`, `turtle` or `n-triple`.
- **Output:** Every network function returns a dictionary with `status` (`completed`), `service`, `upstream` (the Biopython function that was called), `parameters` (the validated request parameters), `retrieved_utc`, `biopython_version`, `output_path` (a path or `None`) and `data`. `data` has `rows` (a list of lists of strings) and `row_count`.
- **Example code:**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.sequence import kegg_convert_identifiers

try:
    kegg_convert_identifiers("ncbi-geneid", "hsa:10458", option="xml", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error)
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] option must be one of None, turtle, n-triple.
```

Usage with the dependency installed:

```python
from dnakit.sequence import kegg_convert_identifiers

result = kegg_convert_identifiers("ncbi-geneid", ["hsa:10458", "hsa:10459"])
print(result["data"]["rows"])
```

## 4) `dnakit.sequence.kegg_database_info()`

- **Purpose:** Describe a KEGG database, with its release and statistics.
- **Method:** Validate the name and call `Bio.KEGG.REST.kegg_info()`; the reply text is returned unchanged.
- **API:** `kegg_database_info(database: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package Biopython (extra `dnakit[bio-formats]`) and internet access to the KEGG REST service (rest.kegg.jp). No key is required. The example below passes an argument that is rejected before any request is sent, so it runs offline; to run the function for real you need the network.
- **Input:** `database` is a KEGG database or organism identifier such as `pathway` or `hsa`: letters, digits and `_ . : | + -` only, at most 256 characters.
- **Output:** Every network function returns a dictionary with `status` (`completed`), `service`, `upstream` (the Biopython function that was called), `parameters` (the validated request parameters), `retrieved_utc`, `biopython_version`, `output_path` (a path or `None`) and `data`. `data` has `text` and `line_count`.
- **Example code:**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.sequence import kegg_database_info

try:
    kegg_database_info("bad name", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error)
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] database must be an identifier without spaces, commas or slashes.
```

Usage with the dependency installed:

```python
from dnakit.sequence import kegg_database_info

result = kegg_database_info("pathway")
print(result["data"]["text"].splitlines()[0])
```

## 5) `dnakit.sequence.kegg_find_entries()`

- **Purpose:** Search a KEGG database by keyword, molecular formula, exact mass or molecular weight.
- **Method:** Validate the arguments and call `Bio.KEGG.REST.kegg_find()`. The tab-separated reply is split into rows.
- **API:** `kegg_find_entries(database: 'str', query: 'str', *, option: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package Biopython (extra `dnakit[bio-formats]`) and internet access to the KEGG REST service (rest.kegg.jp). No key is required. The example below passes an argument that is rejected before any request is sent, so it runs offline; to run the function for real you need the network.
- **Input:** `database` is a short name such as `compound` or `genes`. `query` is a non-empty string of at most 512 characters. `option` is `None`, `formula`, `exact_mass` or `mol_weight`.
- **Output:** Every network function returns a dictionary with `status` (`completed`), `service`, `upstream` (the Biopython function that was called), `parameters` (the validated request parameters), `retrieved_utc`, `biopython_version`, `output_path` (a path or `None`) and `data`. `data` has `rows` and `row_count`.
- **Example code:**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.sequence import kegg_find_entries

try:
    kegg_find_entries("compound", "C7H10O5", option="weight", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error)
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] option must be one of None, formula, exact_mass, mol_weight.
```

Usage with the dependency installed:

```python
from dnakit.sequence import kegg_find_entries

result = kegg_find_entries("compound", "C7H10O5", option="formula")
print(result["data"]["row_count"])
```

## 6) `dnakit.sequence.kegg_get_entries()`

- **Purpose:** Retrieve KEGG entries as text, in the format chosen by `option`.
- **Method:** Validate the arguments and call `Bio.KEGG.REST.kegg_get()`. The reply text is returned unchanged.
- **API:** `kegg_get_entries(entries: 'str | Sequence[str]', *, option: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package Biopython (extra `dnakit[bio-formats]`) and internet access to the KEGG REST service (rest.kegg.jp). No key is required. The example below passes an argument that is rejected before any request is sent, so it runs offline; to run the function for real you need the network.
- **Input:** `entries` is one identifier string or a list of 1 to 10 identifiers. `option` is `None`, `aaseq`, `ntseq`, `mol`, `kcf`, `kgml` or `json`; `image` is refused.
- **Output:** Every network function returns a dictionary with `status` (`completed`), `service`, `upstream` (the Biopython function that was called), `parameters` (the validated request parameters), `retrieved_utc`, `biopython_version`, `output_path` (a path or `None`) and `data`. `data` has `text` and `line_count`.
- **Example code:**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.sequence import kegg_get_entries

try:
    kegg_get_entries("hsa:10458", option="image", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error)
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] option must be one of None, aaseq, ntseq, mol, kcf, kgml, json; image is not supported.
```

Usage with the dependency installed:

```python
from dnakit.sequence import kegg_get_entries

result = kegg_get_entries(["hsa:10458"], option="aaseq")
print(result["data"]["text"].splitlines()[0])
```

## 7) `dnakit.sequence.kegg_link_entries()`

- **Purpose:** Find the KEGG entries that are linked to a database or to given entries.
- **Method:** Validate the arguments and call `Bio.KEGG.REST.kegg_link()`. The tab-separated reply is split into rows. The accepted options are the same as for `kegg_convert_identifiers()`.
- **API:** `kegg_link_entries(target_db: 'str', source: 'str | Sequence[str]', *, option: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package Biopython (extra `dnakit[bio-formats]`) and internet access to the KEGG REST service (rest.kegg.jp). No key is required. The example below passes an argument that is rejected before any request is sent, so it runs offline; to run the function for real you need the network.
- **Input:** `target_db` is a short name such as `pathway`. `source` is a database name or one identifier string, or a list of 1 to 100 identifiers. `option` is `None`, `turtle` or `n-triple`.
- **Output:** Every network function returns a dictionary with `status` (`completed`), `service`, `upstream` (the Biopython function that was called), `parameters` (the validated request parameters), `retrieved_utc`, `biopython_version`, `output_path` (a path or `None`) and `data`. `data` has `rows` and `row_count`.
- **Example code:**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.sequence import kegg_link_entries

try:
    kegg_link_entries("pathway", "hsa:10458", option="rdf", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error)
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] option must be one of None, turtle, n-triple.
```

Usage with the dependency installed:

```python
from dnakit.sequence import kegg_link_entries

result = kegg_link_entries("pathway", ["hsa:10458"])
print(result["data"]["rows"])
```

## 8) `dnakit.sequence.kegg_list_entries()`

- **Purpose:** List the entries of a KEGG database, or the pathways of one organism.
- **Method:** Validate the arguments and call `Bio.KEGG.REST.kegg_list()` with the database and the optional organism. The tab-separated reply is split into rows.
- **API:** `kegg_list_entries(database: 'str | Sequence[str]', *, organism: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package Biopython (extra `dnakit[bio-formats]`) and internet access to the KEGG REST service (rest.kegg.jp). No key is required. The example below passes an argument that is rejected before any request is sent, so it runs offline; to run the function for real you need the network.
- **Input:** `database` is a database name or one entry string, or a list of 1 to 100 identifiers. `organism` is `None` or a short name such as `hsa`; it is passed to KEGG as the second argument of the list request.
- **Output:** Every network function returns a dictionary with `status` (`completed`), `service`, `upstream` (the Biopython function that was called), `parameters` (the validated request parameters), `retrieved_utc`, `biopython_version`, `output_path` (a path or `None`) and `data`. `data` has `rows` and `row_count`.
- **Example code:**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.sequence import kegg_list_entries

try:
    kegg_list_entries("pathway", organism="hsa!", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error)
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] organism must be a short name such as 'pubmed', 'nt' or 'hsa'.
```

Usage with the dependency installed:

```python
from dnakit.sequence import kegg_list_entries

result = kegg_list_entries("pathway", organism="hsa")
print(result["data"]["row_count"])
```

## 9) `dnakit.sequence.ncbi_efetch()`

- **Purpose:** Fetch NCBI Entrez records as text, optionally into a new file.
- **Method:** Validate the arguments and call `Bio.Entrez.efetch()` with the identifiers joined by commas. The reply is read as text and, when `output_path` is given, written to that new file.
- **API:** `ncbi_efetch(db: 'str', ids: 'Sequence[str]', *, rettype: 'str' = 'fasta', retmode: 'str' = 'text', output_path: 'str | None' = None, email: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package Biopython (extra `dnakit[bio-formats]`) and internet access to the NCBI Entrez E-utilities. No key is required and an e-mail address is optional. The example below passes an argument that is rejected before any request is sent, so it runs offline; to run the function for real you need the network.
- **Input:** `db` is a short name such as `nuccore` or `pubmed`. `ids` is a list of 1 to 200 identifier strings. `rettype` (default `fasta`) and `retmode` (default `text`) are short names. `output_path` is `None` or a path that does not exist yet. `email` is `None` or an e-mail address.
- **Output:** Every network function returns a dictionary with `status` (`completed`), `service`, `upstream` (the Biopython function that was called), `parameters` (the validated request parameters), `retrieved_utc`, `biopython_version`, `output_path` (a path or `None`) and `data`. `data` has `text` and `characters`; `output_path` is the written file or `None`.
- **Example code:**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.sequence import ncbi_efetch

try:
    ncbi_efetch("nuccore", [], verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error)
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] ids must hold 1 to 200 identifiers.
```

Usage with the dependency installed:

```python
from dnakit.sequence import ncbi_efetch

result = ncbi_efetch("nuccore", ["NM_000546.6"], rettype="fasta", output_path="tp53.fa")
print(result["data"]["characters"])
```

## 10) `dnakit.sequence.ncbi_einfo()`

- **Purpose:** List the NCBI Entrez databases, or describe one database.
- **Method:** Validate the arguments and call `Bio.Entrez.einfo()` with or without `db`; the XML reply is parsed with `Bio.Entrez.read()` and converted to plain Python types.
- **API:** `ncbi_einfo(*, db: 'str | None' = None, email: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package Biopython (extra `dnakit[bio-formats]`) and internet access to the NCBI Entrez E-utilities. No key is required and an e-mail address is optional. The example below passes an argument that is rejected before any request is sent, so it runs offline; to run the function for real you need the network.
- **Input:** `db` is `None` (list all databases) or a short name such as `pubmed`. `email` is `None` or an e-mail address.
- **Output:** Every network function returns a dictionary with `status` (`completed`), `service`, `upstream` (the Biopython function that was called), `parameters` (the validated request parameters), `retrieved_utc`, `biopython_version`, `output_path` (a path or `None`) and `data`. `data` is the parsed reply: with no `db` it has `DbList`; with a `db` it has `DbInfo`.
- **Example code:**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.sequence import ncbi_einfo

try:
    ncbi_einfo(email="not-an-email", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error)
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] email must be an e-mail address when it is given.
```

Usage with the dependency installed:

```python
from dnakit.sequence import ncbi_einfo

result = ncbi_einfo()
print(len(result["data"]["DbList"]))
```

## 11) `dnakit.sequence.ncbi_elink()`

- **Purpose:** Find the NCBI Entrez records in one database that are linked to identifiers of another database.
- **Method:** Validate the arguments and call `Bio.Entrez.elink()`; the XML reply is parsed with `Bio.Entrez.read()` and converted to plain Python types.
- **API:** `ncbi_elink(dbfrom: 'str', db: 'str', ids: 'Sequence[str]', *, linkname: 'str | None' = None, email: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package Biopython (extra `dnakit[bio-formats]`) and internet access to the NCBI Entrez E-utilities. No key is required and an e-mail address is optional. The example below passes an argument that is rejected before any request is sent, so it runs offline; to run the function for real you need the network.
- **Input:** `dbfrom` and `db` are short names. `ids` is a list of 1 to 200 identifiers, each without spaces, commas or slashes. `linkname` is `None` or a short name. `email` is `None` or an e-mail address.
- **Output:** Every network function returns a dictionary with `status` (`completed`), `service`, `upstream` (the Biopython function that was called), `parameters` (the validated request parameters), `retrieved_utc`, `biopython_version`, `output_path` (a path or `None`) and `data`. `data` is the parsed reply, a list with one entry per request that has `LinkSetDb` when links were found.
- **Example code:**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.sequence import ncbi_elink

try:
    ncbi_elink("gene", "nuccore", ["1,2"], verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error)
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] ids must be an identifier without spaces, commas or slashes.
```

Usage with the dependency installed:

```python
from dnakit.sequence import ncbi_elink

result = ncbi_elink("gene", "nuccore", ["7157"], linkname="gene_nuccore_refseqrna")
print(len(result["data"]))
```

## 12) `dnakit.sequence.ncbi_epost()`

- **Purpose:** Post identifiers to the NCBI history server and return the WebEnv key.
- **Method:** Validate the arguments and call `Bio.Entrez.epost()`; the XML reply is parsed with `Bio.Entrez.read()` and converted to plain Python types.
- **API:** `ncbi_epost(db: 'str', ids: 'Sequence[str]', *, email: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package Biopython (extra `dnakit[bio-formats]`) and internet access to the NCBI Entrez E-utilities. No key is required and an e-mail address is optional. The example below passes an argument that is rejected before any request is sent, so it runs offline; to run the function for real you need the network.
- **Input:** `db` is a short name. `ids` is a list of 1 to 200 identifier strings. `email` is `None` or an e-mail address.
- **Output:** Every network function returns a dictionary with `status` (`completed`), `service`, `upstream` (the Biopython function that was called), `parameters` (the validated request parameters), `retrieved_utc`, `biopython_version`, `output_path` (a path or `None`) and `data`. `data` is the parsed reply with `WebEnv` and `QueryKey`.
- **Example code:**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.sequence import ncbi_epost

try:
    ncbi_epost("nuccore", [str(number) for number in range(201)], verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error)
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] ids must hold 1 to 200 identifiers.
```

Usage with the dependency installed:

```python
from dnakit.sequence import ncbi_epost

result = ncbi_epost("nuccore", ["NM_000546.6", "NM_007294.4"])
print(sorted(result["data"]))
```

## 13) `dnakit.sequence.ncbi_esearch()`

- **Purpose:** Search an NCBI Entrez database and return the total count and the matching UIDs.
- **Method:** Validate the arguments and call `Bio.Entrez.esearch()`; the XML reply is parsed with `Bio.Entrez.read()`, converted to plain Python types, and extended with `count` and `ids`.
- **API:** `ncbi_esearch(db: 'str', term: 'str', *, retmax: 'int' = 20, email: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package Biopython (extra `dnakit[bio-formats]`) and internet access to the NCBI Entrez E-utilities. No key is required and an e-mail address is optional. The example below passes an argument that is rejected before any request is sent, so it runs offline; to run the function for real you need the network.
- **Input:** `db` is a short name such as `pubmed`. `term` is a non-empty string of at most 4096 characters. `retmax` is an integer from 0 to 10000 (default 20). `email` is `None` or an e-mail address.
- **Output:** Every network function returns a dictionary with `status` (`completed`), `service`, `upstream` (the Biopython function that was called), `parameters` (the validated request parameters), `retrieved_utc`, `biopython_version`, `output_path` (a path or `None`) and `data`. `data` is the parsed reply (`Count`, `RetMax`, `IdList`, ...) plus `count` (an integer) and `ids` (a list of strings).
- **Example code:**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.sequence import ncbi_esearch

try:
    ncbi_esearch("pubmed", "BRCA1", retmax=10001, verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error)
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] retmax must be an integer from 0 to 10000.
```

Usage with the dependency installed:

```python
from dnakit.sequence import ncbi_esearch

result = ncbi_esearch("pubmed", "BRCA1[Title] AND 2020[PDAT]", retmax=5)
print(result["data"]["count"], result["data"]["ids"])
```

## 14) `dnakit.sequence.ncbi_esummary()`

- **Purpose:** Retrieve NCBI Entrez document summaries for a list of identifiers.
- **Method:** Validate the arguments and call `Bio.Entrez.esummary()`; the XML reply is parsed with `Bio.Entrez.read()` and converted to plain Python types.
- **API:** `ncbi_esummary(db: 'str', ids: 'Sequence[str]', *, email: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package Biopython (extra `dnakit[bio-formats]`) and internet access to the NCBI Entrez E-utilities. No key is required and an e-mail address is optional. The example below passes an argument that is rejected before any request is sent, so it runs offline; to run the function for real you need the network.
- **Input:** `db` is a short name. `ids` is a list of 1 to 200 identifier strings. `email` is `None` or an e-mail address.
- **Output:** Every network function returns a dictionary with `status` (`completed`), `service`, `upstream` (the Biopython function that was called), `parameters` (the validated request parameters), `retrieved_utc`, `biopython_version`, `output_path` (a path or `None`) and `data`. `data` is the parsed reply: its shape depends on the database (a list of summary records for most databases).
- **Example code:**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.sequence import ncbi_esummary

try:
    ncbi_esummary("pub med", ["1"], verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error)
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] db must be a short name such as 'pubmed', 'nt' or 'hsa'.
```

Usage with the dependency installed:

```python
from dnakit.sequence import ncbi_esummary

result = ncbi_esummary("pubmed", ["31452104"])
print(result["data"][0]["Title"])
```

## 15) `dnakit.sequence.ncbiwww_qblast_search()`

- **Purpose:** Run a BLAST search on the NCBI servers through the older `NCBIWWW` module and summarise the hits.
- **Method:** Same as `blast_qblast_search()` except that the upstream call is `Bio.Blast.NCBIWWW.qblast()` with `format_type="XML"`. The XML is read completely, optionally saved to `output_xml`, and parsed with `Bio.Blast.NCBIXML.read()` into a plain summary.
- **API:** `ncbiwww_qblast_search(program: 'str', database: 'str', sequence: 'str', *, hitlist_size: 'int' = 50, expect: 'float' = 10.0, entrez_query: 'str' = '(none)', output_xml: 'str | None' = None, email: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package Biopython (extra `dnakit[bio-formats]`) and internet access to NCBI BLAST (QBLAST). No key is required and an e-mail address is optional. A search is queued on the NCBI server and can take minutes. The example below passes an argument that is rejected before any request is sent, so it runs offline.
- **Input:** The same arguments and limits as `blast_qblast_search()`. `program` is one of `blastn`, `blastp`, `blastx`, `tblastn`, `tblastx`. `database` is a short name such as `nt`. `sequence` is a non-empty string of at most 100000 characters (it is not checked for DNA letters, although the docstring says DNA query). `hitlist_size` is an integer from 1 to 500, `expect` is a positive number, `entrez_query` is at most 512 characters, `output_xml` is an optional new file path and `email` an optional address.
- **Output:** Every network function returns a dictionary with `status` (`completed`), `service`, `upstream` (the Biopython function that was called), `parameters` (the validated request parameters), `retrieved_utc`, `biopython_version`, `output_path` (a path or `None`) and `data`. `parameters` holds the query length instead of the query. `data` has `query_id`, `query`, `query_letters`, `database`, `alignment_count` and `alignments`; each alignment has `hit_id`, `hit_def`, `length` and `hsps` (`score`, `bits`, `expect`, `identities`, `align_length`, `query_start`, `query_end`, `subject_start`, `subject_end`). When `output_xml` is given the raw XML is also written there.
- **Example code:**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.sequence import ncbiwww_qblast_search

try:
    ncbiwww_qblast_search("blastn", "nt", "   ", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error)
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] sequence must be a non-empty string.
```

Usage with the dependency installed:

```python
from dnakit.sequence import ncbiwww_qblast_search

result = ncbiwww_qblast_search(
    "blastn", "nt", "ATGGCGGCGCTGAGCGGTGGCGGTGGTGGCGGCGCGGAGCCGGGCCAGGCTCTGTTCAACGGGGACATGGAGCCCGAGGCCGGC",
    hitlist_size=5, expect=1e-5, output_xml="blast_result.xml",
)
print(result["data"]["alignment_count"])
```

## 16) `dnakit.sequence.phylo_draw_ascii_text()`

- **Purpose:** Draw a Newick tree as ASCII art and save the drawing as a new text file.
- **Method:** Parse the tree with `Bio.Phylo.read()` (format `newick`), draw it with `Bio.Phylo.draw_ascii()` into a string buffer of the requested width, write the text to `output_path` and return it.
- **API:** `phylo_draw_ascii_text(newick: 'str', output_path: 'str', *, column_width: 'int' = 80, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package Biopython (extra `dnakit[bio-formats]`). It works offline and does not need matplotlib.
- **Input:** `newick` is a non-empty Newick string of at most 5000000 characters. `output_path` must be a new file path in an existing folder. `column_width` is an integer from 20 to 1000 (default 80).
- **Output:** A dictionary with `status` (`completed`), `output_path` and `data`. `data` has `text` (the drawing), `line_count`, `tip_names` and `column_width`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import phylo_draw_ascii_text

target = Path(tempfile.mkdtemp()) / "tree.txt"
result = phylo_draw_ascii_text("((A:1,B:2)0.9:1,C:3);", target, column_width=40, verbose=False)
data = result["data"]
for line in data["text"].splitlines():
    if line.strip():
        print("|" + line.rstrip())
print(data["tip_names"], data["column_width"])
print(target.read_text() == data["text"])
```

- **Example output:**

```text
|              ___________ A
|  ___________|
|_|           |_______________________ B
| |
| |___________________________________ C
['A', 'B', 'C'] 40
True
```

## 17) `dnakit.sequence.phylo_draw_newick_png()`

- **Purpose:** Draw a Newick tree as a new PNG image and return the tree data.
- **Method:** Parse the tree with `Bio.Phylo.read()` (format `newick`), draw it with `Bio.Phylo.draw()` on a matplotlib `Figure` with the Agg canvas (6 x 4 inches, no pyplot state, nothing is shown) and save the PNG. The tree is also written back with `Bio.Phylo.write()` to give a normalised Newick string.
- **API:** `phylo_draw_newick_png(newick: 'str', image_path: 'str', *, show_confidence: 'bool' = True, dpi: 'int' = 150, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages Biopython (extra `dnakit[bio-formats]`) and matplotlib (extra `dnakit[bio-plots]`). It works offline and needs no display.
- **Input:** `newick` is a non-empty Newick string of at most 5000000 characters. `image_path` must be a new path ending in `.png` in an existing folder. `show_confidence` (default `True`) labels internal nodes that have a confidence value. `dpi` is an integer from 50 to 600 (default 150).
- **Output:** A dictionary with `status` (`completed`), `image_path`, `format` (`png`) and `data`. `data` has `newick` (normalised), `tip_names`, `terminal_count`, `clade_count` and `total_branch_length`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import phylo_draw_newick_png

image = Path(tempfile.mkdtemp()) / "tree.png"
result = phylo_draw_newick_png("((A:1,B:2)0.9:1,C:3);", image, verbose=False)
data = result["data"]
print(result["status"], result["format"], image.read_bytes()[:4])
print(data["newick"])
print(data["tip_names"], data["terminal_count"], data["clade_count"], data["total_branch_length"])
```

- **Example output:**

```text
completed png b'\x89PNG'
((A:1,B:2)0.90:1,C:3):0;
['A', 'B', 'C'] 3 5 7.0
```

## 18) `dnakit.sequence.togows_convert_data()`

- **Purpose:** Convert a record between file formats with the TogoWS service.
- **Method:** Validate the arguments and call `Bio.TogoWS.convert()`. The reply is read as text.
- **API:** `togows_convert_data(data: 'str', in_format: 'str', out_format: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package Biopython (extra `dnakit[bio-formats]`) and internet access to TogoWS (togows.dbcls.jp). No key is required. The example below passes an argument that is rejected before any request is sent, so it runs offline; to run the function for real you need the network.
- **Input:** `data` is a non-empty string of at most 10000000 characters. `in_format` and `out_format` are short format names that TogoWS supports (for example `genbank` and `fasta`).
- **Output:** Every network function returns a dictionary with `status` (`completed`), `service`, `upstream` (the Biopython function that was called), `parameters` (the validated request parameters), `retrieved_utc`, `biopython_version`, `output_path` (a path or `None`) and `data`. `parameters` holds the two formats and the input length. `data` has `text` and `characters`.
- **Example code:**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.sequence import togows_convert_data

try:
    togows_convert_data(">x\nACGT", "1bad", "fasta", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error)
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] in_format must be a short name such as 'pubmed', 'nt' or 'hsa'.
```

Usage with the dependency installed:

```python
from dnakit.sequence import togows_convert_data

genbank_text = open("record.gb").read()
result = togows_convert_data(genbank_text, "genbank", "fasta")
print(result["data"]["text"].splitlines()[0])
```

## 19) `dnakit.sequence.togows_count_matches()`

- **Purpose:** Count the matches of a query in a TogoWS database.
- **Method:** Validate the arguments and call `Bio.TogoWS.search_count()`; the count is returned as an integer.
- **API:** `togows_count_matches(database: 'str', query: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package Biopython (extra `dnakit[bio-formats]`) and internet access to TogoWS (togows.dbcls.jp). No key is required. The example below passes an argument that is rejected before any request is sent, so it runs offline; to run the function for real you need the network.
- **Input:** `database` is a short name such as `ncbi-nucleotide` or `pubmed`. `query` is a non-empty string of at most 1024 characters.
- **Output:** Every network function returns a dictionary with `status` (`completed`), `service`, `upstream` (the Biopython function that was called), `parameters` (the validated request parameters), `retrieved_utc`, `biopython_version`, `output_path` (a path or `None`) and `data`. `data` has `count`.
- **Example code:**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.sequence import togows_count_matches

try:
    togows_count_matches("ncbi-nucleotide", "  ", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error)
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] query must be a non-empty string.
```

Usage with the dependency installed:

```python
from dnakit.sequence import togows_count_matches

result = togows_count_matches("pubmed", "lung and cancer")
print(result["data"]["count"])
```

## 20) `dnakit.sequence.togows_fetch_entry()`

- **Purpose:** Fetch one TogoWS record, optionally in one format or reduced to one field.
- **Method:** Validate the arguments and call `Bio.TogoWS.entry()`. The reply is read as text.
- **API:** `togows_fetch_entry(database: 'str', identifier: 'str', *, output_format: 'str | None' = None, field: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package Biopython (extra `dnakit[bio-formats]`) and internet access to TogoWS (togows.dbcls.jp). No key is required. The example below passes an argument that is rejected before any request is sent, so it runs offline; to run the function for real you need the network.
- **Input:** `database` is a short name. `identifier` is one identifier of letters, digits and `_ . : | + -` (no slashes or commas). `output_format` and `field` are `None` or short names.
- **Output:** Every network function returns a dictionary with `status` (`completed`), `service`, `upstream` (the Biopython function that was called), `parameters` (the validated request parameters), `retrieved_utc`, `biopython_version`, `output_path` (a path or `None`) and `data`. `data` has `text` and `characters`.
- **Example code:**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.sequence import togows_fetch_entry

try:
    togows_fetch_entry("ncbi-nucleotide", "a/b", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error)
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] identifier must be an identifier without spaces, commas or slashes.
```

Usage with the dependency installed:

```python
from dnakit.sequence import togows_fetch_entry

result = togows_fetch_entry("ncbi-nucleotide", "X52280", output_format="fasta")
print(result["data"]["text"].splitlines()[0])
```

## 21) `dnakit.sequence.togows_iterate_ids()`

- **Purpose:** Collect up to a bounded number of identifiers for a TogoWS query.
- **Method:** Validate the arguments, call `Bio.TogoWS.search_iter()` and take at most `limit` identifiers from the iterator with `itertools.islice`.
- **API:** `togows_iterate_ids(database: 'str', query: 'str', *, limit: 'int' = 100, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package Biopython (extra `dnakit[bio-formats]`) and internet access to TogoWS (togows.dbcls.jp). No key is required. The example below passes an argument that is rejected before any request is sent, so it runs offline; to run the function for real you need the network.
- **Input:** `database` is a short name. `query` is a non-empty string of at most 1024 characters. `limit` is an integer from 1 to 1000 (default 100).
- **Output:** Every network function returns a dictionary with `status` (`completed`), `service`, `upstream` (the Biopython function that was called), `parameters` (the validated request parameters), `retrieved_utc`, `biopython_version`, `output_path` (a path or `None`) and `data`. `data` has `ids` and `returned_count`.
- **Example code:**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.sequence import togows_iterate_ids

try:
    togows_iterate_ids("ncbi-nucleotide", "x", limit=1001, verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error)
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] limit must be an integer from 1 to 1000.
```

Usage with the dependency installed:

```python
from dnakit.sequence import togows_iterate_ids

result = togows_iterate_ids("pubmed", "lung and cancer", limit=25)
print(result["data"]["returned_count"])
```

## 22) `dnakit.sequence.togows_search_ids()`

- **Purpose:** Search a TogoWS database and return one page of identifiers, or the raw reply in a given format.
- **Method:** Validate the arguments and call `Bio.TogoWS.search()`. Without `output_format` the reply is split on whitespace into identifiers; with it the raw text is returned.
- **API:** `togows_search_ids(database: 'str', query: 'str', *, offset: 'int | None' = None, limit: 'int | None' = None, output_format: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package Biopython (extra `dnakit[bio-formats]`) and internet access to TogoWS (togows.dbcls.jp). No key is required. The example below passes an argument that is rejected before any request is sent, so it runs offline; to run the function for real you need the network.
- **Input:** `database` is a short name. `query` is a non-empty string of at most 1024 characters. `offset` is `None` or an integer from 1 to 1000000000 (TogoWS counts from 1). `limit` is `None` or an integer from 1 to 1000000. `output_format` is `None` or a short name such as `json`.
- **Output:** Every network function returns a dictionary with `status` (`completed`), `service`, `upstream` (the Biopython function that was called), `parameters` (the validated request parameters), `retrieved_utc`, `biopython_version`, `output_path` (a path or `None`) and `data`. `data` has `ids` and `count` or, with `output_format`, `text`.
- **Example code:**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.sequence import togows_search_ids

try:
    togows_search_ids("ncbi-nucleotide", "x", offset=0, verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error)
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] offset must be an integer from 1 to 1000000000.
```

Usage with the dependency installed:

```python
from dnakit.sequence import togows_search_ids

result = togows_search_ids("pubmed", "lung and cancer", offset=1, limit=10)
print(result["data"]["count"], result["data"]["ids"][:3])
```
