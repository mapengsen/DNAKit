# cooler and cooltools: Hi-C contact matrices

This group of functions in `dnakit.epigenome` is a set of adapters over [cooler](https://github.com/open2c/cooler) (BSD-3-Clause), the HDF5 format and command-line tool for Hi-C contact matrices, and [cooltools](https://github.com/open2c/cooltools) (MIT), the analysis toolkit built on it. They build a cooler from read pairs, pixel tables or hiclib files, bin, sort and index contact lists, balance, coarsen, merge, zoom and dump matrices, inspect and reorganise cooler files, draw contact-matrix regions to image files, and list and fetch cooltools example data or draw a saddle plot. They do not call loops, domains or compartments.

Most functions run exactly one `cooler` subcommand with a typed keyword list that maps to fixed flags; no shell is used and no other option can be passed. The `cooler` program is looked up as an explicit `cooler_path`, then on `PATH`, then next to the running Python interpreter. These functions return a dictionary with `status` (`"completed"` or `"failed"`), `backend`, `backend_version`, `subcommand`, `command`, `return_code`, `outputs` (the existing output files) and `log_tail` (the end of cooler's log), plus the extra fields named in each section. Invalid arguments raise `ConfigurationError`; a missing program or package raises `BackendUnavailableError`. A few functions run in process instead (`cooler_csort_flip_pairs`, `cooler_pairs_header`, `cooler_show_view`, `cooler_resolve_bins` and the cooltools functions). Output files must be new paths unless a section says otherwise.

The functions are imported with `from dnakit.epigenome import ...`. The dependencies are optional and are imported only when a function is called. If a package is missing, the function raises an error that includes the install command. The cooler functions need the `dnakit[bio-cooler]` extra, the cooltools functions the `dnakit[bio-cooltools]` extra, and the functions that draw images also `dnakit[bio-plots]`.

## Example data

All examples use a two-chromosome genome of 160 bp that is binned at 20 bp (8 bins). Each example builds its own cooler in a temporary folder from the first two files, so nothing is written next to the example files.

**Example file:** `chrom.sizes`

```text
chr1	100
chr2	60
```

**Example file:** `pairs.txt`

```text
## pairs format v1.0
#columns: readID chrom1 pos1 chrom2 pos2
r1	chr1	5	chr1	45
r2	chr1	10	chr1	30
r3	chr1	25	chr1	85
r4	chr1	50	chr2	10
r5	chr2	5	chr2	50
r6	chr1	95	chr1	3
r7	chr1	12	chr2	45
r8	chr1	5	chr3	7
```

The last read pair is on `chr3`, which is not in `chrom.sizes`, so binning drops it. The sixth pair has its mates in the wrong order and is flipped.

**Example file:** `pixels.coo`

```text
0	1	3
0	2	1
3	3	4
```

**Example file:** `genome.fa`

```text
>chr1
AAAAAAAAAAAAAAAAAAAAAAGCTTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCAAGCTTGGGGGGGGGGGGGGGGGGGGGGGGGGGG
>chr2
TTTTTTTTTTTTTTTTTTTTTTTTTTTTAAGCTTCCCCCCCCCCCCCCCCCCCCCCCCCC
```

**Example file:** `datasets.tsv`

```text
# key	comment	link	filename	checksum
tiny	A tiny text file	https://example.org/tiny.txt	tiny.txt	b1946ac92492d2347c6235b4d2611184
small	Another file	https://example.org/small.txt	small.txt	00000000000000000000000000000000
```

## 1) `dnakit.epigenome.cooler_attrs()`

- **Purpose:** Return the HDF5 attributes of a cooler file or group as YAML text and as a dictionary.
- **Method:** Run `cooler attrs` with an optional `--level`. The text cooler prints (`cooler.fileops.pprint_attr_tree`) is returned in `attributes_yaml` and, when PyYAML is installed, parsed with `yaml.safe_load` into `attributes`.
- **API:** `cooler_attrs(uri: 'str', *, level: 'int | None' = None, cooler_path: 'str | None' = None, timeout_seconds: 'float' = 600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package cooler (extra `dnakit[bio-cooler]`); PyYAML is a cooler dependency.
- **Input:** `uri` is an existing cooler file or `file::/group`. `level` is the maximum depth (an integer of at least 0).
- **Output:** The common fields plus `uri`, `attributes_yaml` and `attributes` (`None` if the run failed). The attributes of the root group are under the key `"@attrs"`.
- **Example code:**

```python
import tempfile

from dnakit.epigenome import cooler_attrs, cooler_cload_pairs

with tempfile.TemporaryDirectory() as folder:
    cool = f"{folder}/a.cool"
    cooler_cload_pairs("chrom.sizes:20", "pairs.txt", cool, chrom1=2, pos1=3, chrom2=4, pos2=5, verbose=False)
    result = cooler_attrs(cool, verbose=False)
    root = result["attributes"]["@attrs"]
    print(result["status"], root["bin-size"], root["nbins"], root["sum"])
    print(sorted(key for key in result["attributes"] if key != "@attrs"))
```

- **Example output:**

```text
completed 20 8 7
['bins', 'chroms', 'indexes', 'pixels']
```

## 2) `dnakit.epigenome.cooler_balance()`

- **Purpose:** Balance a cooler by iterative correction and store the weights in its bin table.
- **Method:** Run `cooler balance` on the file, which is modified in place (a `weight` column, or the column called `name`, is added to `bins`). When the run completes, DNAKit calls `cooler_balance_check()` and reports the answer as `weights_stored`. Options left as `None` keep cooler's defaults (ignore_diags 2, mad_max 5, min_nnz 10, min_count 0, tol 1e-5, max_iters 200, convergence_policy store_final, 8 processes). `cooler balance --stdout` is not exposed: read the weights with `cooler_dump(table="bins", columns=["weight"])`.
- **API:** `cooler_balance(cool_uri: 'str', *, cis_only: 'bool' = False, trans_only: 'bool' = False, ignore_diags: 'int | None' = None, ignore_dist: 'int | None' = None, mad_max: 'int | None' = None, min_nnz: 'int | None' = None, min_count: 'int | None' = None, blacklist: 'str | None' = None, nproc: 'int | None' = None, chunksize: 'int | None' = None, tol: 'float | None' = None, max_iters: 'int | None' = None, name: 'str | None' = None, force: 'bool' = False, convergence_policy: "Literal['store_final', 'store_nan', 'discard', 'error'] | None" = None, cooler_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package cooler (extra `dnakit[bio-cooler]`).
- **Input:** `cool_uri` is an existing cooler. `cis_only` and `trans_only` cannot both be true. The integer options must be at least 0 (`nproc`, `chunksize` and `max_iters` at least 1). `blacklist` is an existing BED file with at least two lines (cooler 0.10.4 reads a one-line file as a header and fails). `force` recomputes existing weights.
- **Output:** The common fields plus `cool_uri`, `weights_column` and `weights_stored`. A non-converging run is still stored under the default `store_final` policy, and the reason is in `log_tail`.
- **Example code:**

```python
import tempfile

from dnakit.epigenome import cooler_balance, cooler_cload_pairs

with tempfile.TemporaryDirectory() as folder:
    cool = f"{folder}/a.cool"
    cooler_cload_pairs("chrom.sizes:20", "pairs.txt", cool, chrom1=2, pos1=3, chrom2=4, pos2=5, verbose=False)
    result = cooler_balance(cool, ignore_diags=0, min_nnz=0, mad_max=0, cis_only=True, nproc=1, verbose=False)
    print(result["status"], result["weights_column"], result["weights_stored"])
```

- **Example output:**

```text
completed weight True
```

## 3) `dnakit.epigenome.cooler_balance_check()`

- **Purpose:** Report whether a cooler already has a balancing weight column.
- **Method:** Run `cooler balance --check`. Cooler exits with code 1 both for a cooler without weights and for an unreadable file, so DNAKit reads its message: "is balanced" gives `True`, "column found" gives `False`, anything else is reported as a failed run with `balanced` set to `None`.
- **API:** `cooler_balance_check(cool_uri: 'str', *, name: 'str | None' = None, cooler_path: 'str | None' = None, timeout_seconds: 'float' = 600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package cooler (extra `dnakit[bio-cooler]`).
- **Input:** `cool_uri` is an existing cooler; `name` is the weight column to look for (default `weight`).
- **Output:** The common fields plus `cool_uri`, `weights_column` and `balanced`.
- **Example code:**

```python
import tempfile

from dnakit.epigenome import cooler_balance, cooler_balance_check, cooler_cload_pairs

with tempfile.TemporaryDirectory() as folder:
    cool = f"{folder}/a.cool"
    cooler_cload_pairs("chrom.sizes:20", "pairs.txt", cool, chrom1=2, pos1=3, chrom2=4, pos2=5, verbose=False)
    print(cooler_balance_check(cool, verbose=False)["balanced"])
    cooler_balance(cool, verbose=False)
    print(cooler_balance_check(cool, verbose=False)["balanced"])
```

- **Example output:**

```text
False
True
```

## 4) `dnakit.epigenome.cooler_cload_hiclib()`

- **Purpose:** Bin a hiclib HDF5 contact list into a cooler.
- **Method:** Run `cooler cload hiclib`. The HDF5 file holds the datasets `chrms1`, `cuts1`, `chrms2` and `cuts2` (zero-based chromosome numbers in the order of the bins, zero-based positions). The pairs must be sorted on the first chromosome and upper triangular for the chromosome order of `bins`.
- **API:** `cooler_cload_hiclib(bins: 'str', pairs_path: 'str', output_path: 'str', *, chunksize: 'int | None' = None, metadata_path: 'str | None' = None, assembly: 'str | None' = None, cooler_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package cooler (extra `dnakit[bio-cooler]`); h5py is a cooler dependency.
- **Input:** `bins` is `<chromsizes file>:<bin size in bp>` or the path of a BED file with the bins. `pairs_path` is an existing HDF5 file; `output_path` must not exist. `chunksize` (default 100,000,000 in cooler) is a positive integer. `metadata_path` is an existing JSON file and `assembly` a genome name such as `hg38`.
- **Output:** The common fields (`subcommand` is `"cload hiclib"`) plus `output_uri`.
- **Example code:**

```python
import tempfile

import h5py
import numpy as np

from dnakit.epigenome import cooler_cload_hiclib, cooler_info

with tempfile.TemporaryDirectory() as folder:
    with h5py.File(f"{folder}/contacts.h5", "w") as handle:
        handle["chrms1"] = np.array([0, 0, 0, 1], dtype="int32")
        handle["cuts1"] = np.array([5, 10, 50, 5], dtype="int64")
        handle["chrms2"] = np.array([0, 0, 1, 1], dtype="int32")
        handle["cuts2"] = np.array([45, 30, 10, 50], dtype="int64")
    result = cooler_cload_hiclib("chrom.sizes:20", f"{folder}/contacts.h5", f"{folder}/h.cool", verbose=False)
    info = cooler_info(f"{folder}/h.cool", verbose=False)["info"]
    print(result["status"], result["subcommand"], info["nbins"], info["sum"])
```

- **Example output:**

```text
completed cload hiclib 8 4
```

## 5) `dnakit.epigenome.cooler_cload_pairix()`

- **Purpose:** Bin a bgzip-compressed, pairix-indexed pairs file into a cooler.
- **Method:** Check that the `.px2` index exists next to the file, then run `cooler cload pairix`. Cooler needs the Python package `pypairix` for this mode; without it the run is reported as failed with the import error in `log_tail`.
- **API:** `cooler_cload_pairix(bins: 'str', pairs_path: 'str', output_path: 'str', *, zero_based: 'bool' = False, nproc: 'int | None' = None, max_split: 'int | None' = None, block_char: 'str | None' = None, metadata_path: 'str | None' = None, assembly: 'str | None' = None, cooler_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package cooler (extra `dnakit[bio-cooler]`) and the Python package `pypairix` (`python -m pip install pypairix`); the pairix program is needed to build the index (`conda install -c bioconda pairix`). The example uses placeholder files: it passes the index check and then shows what cooler reports when `pypairix` is not installed; with the package installed, use the call shown below the output.
- **Input:** `bins` as in `cooler_cload_hiclib()`. `pairs_path` is an existing bgzip file with a `.px2` index next to it (otherwise `ConfigurationError`). `nproc` (cooler default 8) and `max_split` (default 2) are positive integers; `block_char` (default `|`) is a short string without quotes or backslashes.
- **Output:** The common fields (`subcommand` is `"cload pairix"`) plus `output_uri`.
- **Example code:**

```python
import tempfile

from dnakit.epigenome import cooler_cload_pairix

with tempfile.TemporaryDirectory() as folder:
    pairs = f"{folder}/pairs.gz"
    for path in (pairs, pairs + ".px2"):
        with open(path, "wb") as handle:
            handle.write(b"x")
    result = cooler_cload_pairix("chrom.sizes:20", pairs, f"{folder}/p.cool", verbose=False)
    print(result["status"], result["subcommand"], "pypairix is required" in result["log_tail"])
```

- **Example output:**

```text
failed cload pairix True
```

Usage with the dependency installed:

```python
from dnakit.epigenome import cooler_csort, cooler_cload_pairix

cooler_csort("pairs.txt", "chrom.sizes", "pairs.sorted.gz", chrom1=2, pos1=3, chrom2=4, pos2=5, index="pairix")
cooler_cload_pairix("chrom.sizes:20", "pairs.sorted.gz", "pairix.cool")
```

## 6) `dnakit.epigenome.cooler_cload_pairs()`

- **Purpose:** Bin an unsorted text file of read pairs into a cooler.
- **Method:** Run `cooler cload pairs` with the four field numbers. Cooler sorts and merges the pairs itself, orients every contact upper triangular, and drops pairs on chromosomes that are not in the bins.
- **API:** `cooler_cload_pairs(bins: 'str', pairs_path: 'str', output_path: 'str', *, chrom1: 'int', pos1: 'int', chrom2: 'int', pos2: 'int', zero_based: 'bool' = False, comment_char: 'str | None' = None, no_symmetric_upper: 'bool' = False, input_copy_status: "Literal['unique', 'duplex'] | None" = None, fields: 'Sequence[str] | None' = None, chunksize: 'int | None' = None, mergebuf: 'int | None' = None, max_merge: 'int | None' = None, temp_dir: 'str | None' = None, no_delete_temp: 'bool' = False, storage_options: 'str | None' = None, metadata_path: 'str | None' = None, assembly: 'str | None' = None, append: 'bool' = False, cooler_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package cooler (extra `dnakit[bio-cooler]`).
- **Input:** `bins` is `<chromsizes file>:<bin size in bp>` or a BED file. `pairs_path` is an existing text file, optionally gzipped; standard input is not supported. `chrom1`, `pos1`, `chrom2` and `pos2` are four different one-based column numbers. Positions are one-based unless `zero_based` is true. Comment lines start with `comment_char` (default `#`). `fields` entries look like `score=7:dtype=float,agg=mean`. `storage_options` looks like `compression=gzip,compression_opts=6`. With `append=True` the cooler is added to an existing file and `output_path` may be `file.mcool::resolutions/10000`; otherwise `output_path` must not exist.
- **Output:** The common fields (`subcommand` is `"cload pairs"`) plus `output_uri`.
- **Example code:**

```python
import tempfile

from dnakit.epigenome import cooler_cload_pairs, cooler_info

with tempfile.TemporaryDirectory() as folder:
    cool = f"{folder}/a.cool"
    result = cooler_cload_pairs("chrom.sizes:20", "pairs.txt", cool, chrom1=2, pos1=3, chrom2=4, pos2=5, verbose=False)
    info = cooler_info(cool, verbose=False)["info"]
    print(result["status"], result["subcommand"], result["return_code"])
    print(info["bin-size"], info["nbins"], info["nnz"], info["sum"])
```

- **Example output:**

```text
completed cload pairs 0
20 8 7 7
```

## 7) `dnakit.epigenome.cooler_cload_tabix()`

- **Purpose:** Bin a bgzip-compressed, tabix-indexed pairs file into a cooler.
- **Method:** Check that a `.tbi` or `.csi` index exists next to the file, then run `cooler cload tabix`. The index must be built on chrom1 and pos1, and the chromosome order of `bins` must agree with the order in the file.
- **API:** `cooler_cload_tabix(bins: 'str', pairs_path: 'str', output_path: 'str', *, chrom2: 'int | None' = None, pos2: 'int | None' = None, zero_based: 'bool' = False, nproc: 'int | None' = None, max_split: 'int | None' = None, metadata_path: 'str | None' = None, assembly: 'str | None' = None, cooler_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package cooler (extra `dnakit[bio-cooler]`) and pysam, which cooler uses to read the file. The example builds its index with pysam; the `tabix` and `bgzip` programs (htslib) are only needed for `cooler_csort()`.
- **Input:** `bins` as above. `pairs_path` is an existing bgzip file with an index next to it (otherwise `ConfigurationError`). `chrom2` (cooler default 4) and `pos2` (default 5) are one-based column numbers; `nproc` (default 8) and `max_split` (default 2) are positive integers.
- **Output:** The common fields (`subcommand` is `"cload tabix"`) plus `output_uri`.
- **Example code:**

```python
import tempfile

import pysam

from dnakit.epigenome import cooler_cload_tabix, cooler_csort_flip_pairs, cooler_csort_sort_pairs, cooler_info

with tempfile.TemporaryDirectory() as folder:
    fields = dict(chrom1=2, pos1=3, chrom2=4, pos2=5)
    cooler_csort_flip_pairs("pairs.txt", "chrom.sizes", f"{folder}/flipped.tsv", verbose=False, **fields)
    cooler_csort_sort_pairs(f"{folder}/flipped.tsv", f"{folder}/sorted.tsv", index="tabix", verbose=False, **fields)
    pysam.tabix_compress(f"{folder}/sorted.tsv", f"{folder}/sorted.tsv.gz")
    pysam.tabix_index(f"{folder}/sorted.tsv.gz", seq_col=1, start_col=2, end_col=2, zerobased=False, force=True)
    result = cooler_cload_tabix("chrom.sizes:20", f"{folder}/sorted.tsv.gz", f"{folder}/t.cool", chrom2=4, pos2=5, verbose=False)
    print(result["status"], result["subcommand"], cooler_info(f"{folder}/t.cool", verbose=False)["info"]["sum"])
```

- **Example output:**

```text
completed cload tabix 7
```

## 8) `dnakit.epigenome.cooler_coarsen()`

- **Purpose:** Coarsen a cooler to a lower resolution by pooling k-by-k pixel tiles.
- **Method:** Run `cooler coarsen --out`. The value columns named in `fields` are aggregated (default: the `count` column, summed).
- **API:** `cooler_coarsen(cool_uri: 'str', output_path: 'str', *, factor: 'int | None' = None, nproc: 'int | None' = None, chunksize: 'int | None' = None, fields: 'Sequence[str] | None' = None, append: 'bool' = False, cooler_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package cooler (extra `dnakit[bio-cooler]`).
- **Input:** `cool_uri` is an existing cooler. `factor` (cooler default 2) is an integer of at least 2; `nproc` (default 1) and `chunksize` (default 10,000,000) are positive integers. `fields` entries look like `score:dtype=float,agg=mean`. With `append=True` the result is added to an existing file and `output_path` may carry `::group`; otherwise it must not exist.
- **Output:** The common fields plus `input_uri` and `output_uri`.
- **Example code:**

```python
import tempfile

from dnakit.epigenome import cooler_cload_pairs, cooler_coarsen, cooler_info

with tempfile.TemporaryDirectory() as folder:
    cool = f"{folder}/a.cool"
    cooler_cload_pairs("chrom.sizes:20", "pairs.txt", cool, chrom1=2, pos1=3, chrom2=4, pos2=5, verbose=False)
    result = cooler_coarsen(cool, f"{folder}/c.cool", factor=2, verbose=False)
    info = cooler_info(f"{folder}/c.cool", verbose=False)["info"]
    print(result["status"], info["bin-size"], info["nbins"], info["nnz"], info["sum"])
```

- **Example output:**

```text
completed 40 5 6 7
```

## 9) `dnakit.epigenome.cooler_cp()`

- **Purpose:** Copy a cooler to another file or to another group of a file.
- **Method:** Run `cooler cp`. With `overwrite=True` DNAKit passes `-w`, which truncates and replaces the whole destination file if it exists; when source and destination are in the same file, cooler 0.10.4 cannot honour it and fails with "unable to truncate a file which is already open".
- **API:** `cooler_cp(source_uri: 'str', destination_uri: 'str', *, overwrite: 'bool' = False, cooler_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package cooler (extra `dnakit[bio-cooler]`).
- **Input:** `source_uri` is an existing cooler (`file` or `file::/group`); `destination_uri` is `file` or `file::/group` in an existing folder.
- **Output:** The common fields plus `source_uri`, `destination_uri` and `source_file`; `outputs` lists the destination file when it exists.
- **Example code:**

```python
import tempfile

from dnakit.epigenome import cooler_cload_pairs, cooler_cp, cooler_info

with tempfile.TemporaryDirectory() as folder:
    cool = f"{folder}/a.cool"
    cooler_cload_pairs("chrom.sizes:20", "pairs.txt", cool, chrom1=2, pos1=3, chrom2=4, pos2=5, verbose=False)
    result = cooler_cp(cool, f"{folder}/b.cool", verbose=False)
    print(result["status"], result["return_code"], cooler_info(f"{folder}/b.cool", field="sum", verbose=False)["value"])
```

- **Example output:**

```text
completed 0 7
```

## 10) `dnakit.epigenome.cooler_csort()`

- **Purpose:** Orient, sort, bgzip and index a contact list so that it can be loaded with `cooler_cload_tabix()` or `cooler_cload_pairix()`.
- **Method:** Run `cooler csort`. Mates are ordered so that every contact is upper triangular for the chromosome order of `chromosomes_path` (pairs on other chromosomes are dropped), rows are sorted (positional for tabix, block for pairix), compressed with bgzip and indexed. Before running, DNAKit checks that `sort`, `bgzip` and `tabix` or `pairix` are on `PATH` (or next to cooler) and raises `BackendUnavailableError` otherwise. The result lists the sorted file and its index (`.tbi` or `.px2`).
- **API:** `cooler_csort(pairs_path: 'str', chromosomes_path: 'str', output_path: 'str', *, chrom1: 'int', pos1: 'int', chrom2: 'int', pos2: 'int', index: "Literal['tabix', 'pairix']" = 'pairix', nproc: 'int | None' = None, zero_based: 'bool' = False, sep: 'str | None' = None, comment_char: 'str | None' = None, cooler_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package cooler (extra `dnakit[bio-cooler]`) and the external programs `sort`, `bgzip` and `tabix` (htslib, `conda install -c bioconda htslib`) or `pairix` on `PATH`. Reading standard input and the free `--sort-options` text are not supported; see `cooler_csort_sort_pairs()` for sort tuning. The programs `bgzip` and `tabix` are not installed in the environment that produced this page, so the example shows the error that is raised in that case.
- **Input:** `pairs_path` is an existing text file; `chromosomes_path` lists the chromosomes in the wanted order (a chromsizes file works); `output_path` must not exist and differ from `pairs_path`. `chrom1`, `pos1`, `chrom2` and `pos2` are four different one-based column numbers. `index` is `"tabix"` or `"pairix"` (default). `sep` defaults to a tab and `comment_char` to `#`.
- **Output:** The common fields (`outputs` holds the sorted file and its index) plus `pairs_path`, `index` and `index_file`.
- **Example code:**

```python
import tempfile

from dnakit.epigenome import cooler_csort
from dnakit.exceptions import BackendUnavailableError

with tempfile.TemporaryDirectory() as folder:
    try:
        cooler_csort("pairs.txt", "chrom.sizes", f"{folder}/pairs.sorted.gz", chrom1=2, pos1=3, chrom2=4, pos2=5, index="tabix", verbose=False)
    except BackendUnavailableError as error:
        print(type(error).__name__, "cooler csort needs bgzip" in str(error))
```

- **Example output:**

```text
BackendUnavailableError True
```

Usage with the dependency installed:

```python
from dnakit.epigenome import cooler_csort

result = cooler_csort("pairs.txt", "chrom.sizes", "pairs.sorted.gz", chrom1=2, pos1=3, chrom2=4, pos2=5, index="tabix")
print(result["status"], result["index_file"])
```

## 11) `dnakit.epigenome.cooler_csort_flip_pairs()`

- **Purpose:** Orient each contact upper triangular for a chromosome order (the flip step of csort).
- **Method:** Run in process with the algorithm that `cooler csort` runs as a `python -c` pipeline stage. Contacts on chromosomes that are not listed are dropped, lines starting with `comment_char` are skipped, and the mates swap when chrom1 ranks after chrom2, or on the same chromosome when pos1 is larger than pos2. A `.gz` input is read through gzip. Output columns are joined with tabs. No cooler or other program is needed.
- **API:** `cooler_csort_flip_pairs(pairs_path: 'str', chromosomes_path: 'str', output_path: 'str', *, chrom1: 'int', pos1: 'int', chrom2: 'int', pos2: 'int', sep: 'str' = '\t', comment_char: 'str' = '#', verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Works after plain pip install.
- **Input:** `pairs_path` and `chromosomes_path` are existing files (the first tab-separated column of `chromosomes_path` gives the chromosome order). `output_path` must not exist. The four field numbers are different one-based integers; a line with too few columns or a non-integer position raises `ConfigurationError` and removes the partial output.
- **Output:** `status`, `backend` (`"python"`), `outputs`, `chromosomes` and the counts `lines_read`, `comment_lines`, `dropped`, `flipped` and `written`.
- **Example code:**

```python
import tempfile

from dnakit.epigenome import cooler_csort_flip_pairs

with tempfile.TemporaryDirectory() as folder:
    out = f"{folder}/flipped.tsv"
    result = cooler_csort_flip_pairs("pairs.txt", "chrom.sizes", out, chrom1=2, pos1=3, chrom2=4, pos2=5, verbose=False)
    print({key: value for key, value in result.items() if key not in ("outputs", "backend")})
    lines = open(out).read().splitlines()
    print(" ".join(lines[5].split()))
```

- **Example output:**

```text
{'status': 'completed', 'chromosomes': 2, 'lines_read': 10, 'comment_lines': 2, 'dropped': 1, 'flipped': 1, 'written': 7}
r6 chr1 3 chr1 95
```

## 12) `dnakit.epigenome.cooler_csort_index_pairs()`

- **Purpose:** Index a sorted, bgzip-compressed contact list with tabix or pairix (the index step of csort).
- **Method:** Check that the file starts with the gzip magic bytes, find the `tabix` or `pairix` program, and run it with the same arguments as `cooler csort`. Tabix indexes chrom1 and pos1 and writes `<pairs_path>.tbi`; pairix indexes chrom1, chrom2, pos1 and pos2 and writes `<pairs_path>.px2`.
- **API:** `cooler_csort_index_pairs(pairs_path: 'str', *, chrom1: 'int', pos1: 'int', chrom2: 'int', pos2: 'int', index: "Literal['tabix', 'pairix']" = 'pairix', zero_based: 'bool' = False, index_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `tabix` (htslib, `conda install -c bioconda htslib`) or `pairix` (`conda install -c bioconda pairix`) on `PATH`, next to the Python interpreter, or given as `index_path`. The example gives an `index_path` that does not exist, so that its result does not depend on the machine.
- **Input:** `pairs_path` is an existing bgzip file sorted in the order that `cooler_csort_sort_pairs()` writes for the same `index` (otherwise `ConfigurationError`). The four field numbers are different one-based integers. `index_path` is the executable to use.
- **Output:** `status`, `backend`, `command`, `return_code`, `outputs`, `index_file` and `log_tail`.
- **Example code:**

```python
import gzip
import tempfile

from dnakit.epigenome import cooler_csort_index_pairs
from dnakit.exceptions import BackendUnavailableError

with tempfile.TemporaryDirectory() as folder:
    pairs = f"{folder}/pairs.sorted.gz"
    with gzip.open(pairs, "wt") as handle:
        handle.write("r1\tchr1\t5\tchr1\t45\n")
    try:
        cooler_csort_index_pairs(pairs, chrom1=2, pos1=3, chrom2=4, pos2=5, index="tabix", index_path="/nonexistent/tabix", verbose=False)
    except BackendUnavailableError as error:
        print(type(error).__name__, "tabix_path does not exist" in str(error))
```

- **Example output:**

```text
BackendUnavailableError True
```

Usage with the dependency installed:

```python
from dnakit.epigenome import cooler_csort_index_pairs

result = cooler_csort_index_pairs("pairs.sorted.gz", chrom1=2, pos1=3, chrom2=4, pos2=5, index="tabix")
print(result["status"], result["index_file"])
```

## 13) `dnakit.epigenome.cooler_csort_read_pairs()`

- **Purpose:** Write a contact list as plain text, decompressing `.gz` input (the read step of csort).
- **Method:** As in `cooler csort`, a file ending in `.gz` is read with `pigz -dc` when pigz is available and with `gzip -dc` otherwise; the program's output goes straight to `output_path`. Any other file is copied unchanged (cooler pipes it through `cat`) and `command` is `None`.
- **API:** `cooler_csort_read_pairs(pairs_path: 'str', output_path: 'str', *, pigz_path: 'str | None' = None, gzip_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Plain files need nothing. A `.gz` file needs the external program `pigz` (`conda install -c conda-forge pigz`) or `gzip` on `PATH`; with neither the function raises `BackendUnavailableError`.
- **Input:** `pairs_path` is an existing file; `output_path` must not exist. `pigz_path` and `gzip_path` are optional explicit program paths and must exist when given.
- **Output:** `status`, `backend` (`"copy"`, `"pigz"` or `"gzip"`), `command`, `return_code`, `outputs`, `bytes_written` and `log_tail`.
- **Example code:**

```python
import gzip
import os
import shutil
import tempfile

from dnakit.epigenome import cooler_csort_read_pairs

with tempfile.TemporaryDirectory() as folder:
    with open("pairs.txt", "rb") as source, gzip.open(f"{folder}/pairs.txt.gz", "wb") as target:
        shutil.copyfileobj(source, target)
    result = cooler_csort_read_pairs(f"{folder}/pairs.txt.gz", f"{folder}/plain.txt", verbose=False)
    print(result["status"], result["return_code"], result["bytes_written"] == os.path.getsize("pairs.txt"))
    copied = cooler_csort_read_pairs("pairs.txt", f"{folder}/copy.txt", verbose=False)
    print(copied["backend"], copied["command"], copied["bytes_written"])
```

- **Example output:**

```text
completed 0 True
copy None 209
```

## 14) `dnakit.epigenome.cooler_csort_sort_pairs()`

- **Purpose:** Sort a contact list with cooler csort's keys (the sort step of csort).
- **Method:** Run the system `sort` with `LC_ALL=C` (set for the child process only) and write the result with `sort -o`. `index="tabix"` sorts by chrom1, pos1, chrom2, pos2 (positional); `index="pairix"` by chrom1, chrom2, pos1, pos2 (block). Positions sort numerically, names byte-wise. `nproc` adds `--parallel`, `buffer_size` adds `--buffer-size` and `temp_dir` adds `-T`.
- **API:** `cooler_csort_sort_pairs(pairs_path: 'str', output_path: 'str', *, chrom1: 'int', pos1: 'int', chrom2: 'int', pos2: 'int', index: "Literal['tabix', 'pairix']" = 'pairix', nproc: 'int | None' = None, buffer_size: 'str | None' = None, temp_dir: 'str | None' = None, sort_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program `sort` (GNU coreutils; present on Linux and macOS).
- **Input:** `pairs_path` is an existing tab-separated file that should already be flipped (see `cooler_csort_flip_pairs()`); `output_path` must not exist. The four field numbers are different one-based integers. `buffer_size` looks like `50%` or `2G`; `temp_dir` is an existing folder; `sort_path` is an explicit program path.
- **Output:** `status`, `backend` (`"sort"`), `command`, `return_code`, `outputs` and `log_tail`.
- **Example code:**

```python
import tempfile

from dnakit.epigenome import cooler_csort_flip_pairs, cooler_csort_sort_pairs

with tempfile.TemporaryDirectory() as folder:
    fields = dict(chrom1=2, pos1=3, chrom2=4, pos2=5)
    cooler_csort_flip_pairs("pairs.txt", "chrom.sizes", f"{folder}/flipped.tsv", verbose=False, **fields)
    result = cooler_csort_sort_pairs(f"{folder}/flipped.tsv", f"{folder}/sorted.tsv", index="tabix", verbose=False, **fields)
    print(result["status"], result["return_code"])
    for line in open(f"{folder}/sorted.tsv").read().splitlines()[:3]:
        print(" ".join(line.split()))
```

- **Example output:**

```text
completed 0
r6 chr1 3 chr1 95
r1 chr1 5 chr1 45
r2 chr1 10 chr1 30
```

## 15) `dnakit.epigenome.cooler_digest()`

- **Purpose:** Write the restriction-fragment bins of a genome as a BED file.
- **Method:** Run `cooler digest --out`. Cooler finds the cut sites of a Biopython restriction enzyme in each chromosome and writes one bin per fragment, in the chromosome order of the chromsizes file. The result adds the line count and the first lines of the BED file.
- **API:** `cooler_digest(chromsizes_path: 'str', fasta_path: 'str', enzyme: 'str', output_path: 'str', *, header: 'bool' = False, rel_ids: 'Literal[0, 1] | None' = None, cooler_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package cooler (extra `dnakit[bio-cooler]`) and, for this subcommand, Biopython and pyfaidx; pyfaidx may create a `.fai` index next to the FASTA file. The example prints the fragment boundaries that cooler reports for the two `AAGCTT` sites of `chr1`.
- **Input:** `chromsizes_path` is an existing chromsizes file. `fasta_path` is an uncompressed genome FASTA file or a folder of `*.fa` / `*.fasta` files. `enzyme` is a Biopython enzyme name such as `HindIII` or `DpnII`. `header` adds a header line; `rel_ids` (0 or 1) adds a column of per-chromosome fragment numbers. `output_path` must not exist.
- **Output:** The common fields plus `line_count` and `preview` (up to five lines, as a list of strings).
- **Example code:**

```python
import tempfile

from dnakit.epigenome import cooler_digest

with tempfile.TemporaryDirectory() as folder:
    result = cooler_digest("chrom.sizes", "genome.fa", "HindIII", f"{folder}/fragments.bed", verbose=False)
    print(result["status"], result["line_count"])
    for line in result["preview"]:
        print(" ".join(line.split()))
```

- **Example output:**

```text
completed 5
chr1 0 23
chr1 23 69
chr1 69 100
chr2 0 31
chr2 31 60
```

## 16) `dnakit.epigenome.cooler_dump()`

- **Purpose:** Dump a cooler table (chroms, bins or pixels) to tab-separated text.
- **Method:** Run `cooler dump --out`. For pixels, `region` and `region2` select a UCSC-style query box, `balanced` adds a `balanced` column (weight1 * weight2 * count), `join` replaces bin ids with chromosome coordinates, `annotate` joins extra bin-table columns (suffixed 1 and 2) and `fill_lower` generates the implicit lower triangle of a symmetric-upper cooler. `one_based_ids` and `one_based_starts` are applied by cooler only together with `join`, `balanced` or `annotate`. The result adds the line count and the first lines of the file.
- **API:** `cooler_dump(cool_uri: 'str', output_path: 'str', *, table: "Literal['chroms', 'bins', 'pixels']" = 'pixels', columns: 'Sequence[str] | None' = None, header: 'bool' = False, na_rep: 'str | None' = None, float_format: 'str | None' = None, region: 'str | None' = None, region2: 'str | None' = None, fill_lower: 'bool' = False, balanced: 'bool' = False, join: 'bool' = False, annotate: 'Sequence[str] | None' = None, one_based_ids: 'bool' = False, one_based_starts: 'bool' = False, chunksize: 'int | None' = None, cooler_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package cooler (extra `dnakit[bio-cooler]`).
- **Input:** `cool_uri` is an existing cooler; `output_path` must not exist. `table` selects what to dump. `columns` (plain column names) applies to chroms and bins. `region2` needs `region`. `float_format` looks like `g`, `.12g` or `03.2f` (default `g`). `na_rep` must be non-empty when given (cooler's default is empty). `chunksize` is a positive integer.
- **Output:** The common fields plus `input_uri`, `table`, `line_count` and `preview`.
- **Example code:**

```python
import tempfile

from dnakit.epigenome import cooler_cload_pairs, cooler_dump

with tempfile.TemporaryDirectory() as folder:
    cool = f"{folder}/a.cool"
    cooler_cload_pairs("chrom.sizes:20", "pairs.txt", cool, chrom1=2, pos1=3, chrom2=4, pos2=5, verbose=False)
    result = cooler_dump(cool, f"{folder}/pixels.tsv", join=True, header=True, verbose=False)
    print(result["status"], result["table"], result["line_count"])
    for line in result["preview"][:3]:
        print(" ".join(line.split()))
    bins = cooler_dump(cool, f"{folder}/bins.tsv", table="bins", columns=["chrom", "end"], verbose=False)
    print(bins["line_count"], bins["preview"])
```

- **Example output:**

```text
completed pixels 8
chrom1 start1 end1 chrom2 start2 end2 count
chr1 0 20 chr1 20 40 1
chr1 0 20 chr1 40 60 1
8 ['chr1\t20', 'chr1\t40', 'chr1\t60', 'chr1\t80', 'chr1\t100']
```

## 17) `dnakit.epigenome.cooler_info()`

- **Purpose:** Return the info block of a cooler: bin size, bin and chromosome counts, pixel sum.
- **Method:** Run `cooler info --out` into a temporary file and parse the JSON. `field` returns one value (parsed as JSON when possible, otherwise as text); `metadata=True` returns the user metadata instead.
- **API:** `cooler_info(cool_uri: 'str', *, field: 'str | None' = None, metadata: 'bool' = False, cooler_path: 'str | None' = None, timeout_seconds: 'float' = 600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package cooler (extra `dnakit[bio-cooler]`).
- **Input:** `cool_uri` is an existing cooler or `file::/group`. `field` is an info key such as `nbins`; it cannot be combined with `metadata`.
- **Output:** The common fields plus `cool_uri`, `info` (all fields), `value` (the one requested field) and `metadata`; the ones that were not requested are `None`.
- **Example code:**

```python
import tempfile

from dnakit.epigenome import cooler_cload_pairs, cooler_info

with tempfile.TemporaryDirectory() as folder:
    cool = f"{folder}/a.cool"
    cooler_cload_pairs("chrom.sizes:20", "pairs.txt", cool, chrom1=2, pos1=3, chrom2=4, pos2=5, verbose=False)
    info = cooler_info(cool, verbose=False)["info"]
    print({key: info[key] for key in ("bin-size", "bin-type", "nbins", "nchroms", "nnz", "storage-mode", "sum")})
    print(cooler_info(cool, field="nbins", verbose=False)["value"], cooler_info(cool, metadata=True, verbose=False)["metadata"])
```

- **Example output:**

```text
{'bin-size': 20, 'bin-type': 'fixed', 'nbins': 8, 'nchroms': 2, 'nnz': 7, 'storage-mode': 'symmetric-upper', 'sum': 7}
8 {}
```

## 18) `dnakit.epigenome.cooler_ln()`

- **Purpose:** Link a cooler instead of copying it.
- **Method:** Run `cooler ln`: a hard link within one file by default; with `soft=True` a soft link in the same file, or an external link when the destination is another file. `overwrite=True` passes `-w`, which truncates the destination file (and fails inside one file, see `cooler_cp()`). Do not link the root group `file::/` into one of its own sub-groups: that makes a cycle, and cooler then fails to list the file with a `RecursionError`.
- **API:** `cooler_ln(source_uri: 'str', destination_uri: 'str', *, overwrite: 'bool' = False, soft: 'bool' = False, cooler_path: 'str | None' = None, timeout_seconds: 'float' = 600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package cooler (extra `dnakit[bio-cooler]`).
- **Input:** `source_uri` is an existing cooler (`file` or `file::/group`); `destination_uri` is `file` or `file::/group` in an existing folder.
- **Output:** The common fields plus `source_uri`, `destination_uri` and `source_file`.
- **Example code:**

```python
import tempfile

from dnakit.epigenome import cooler_cload_pairs, cooler_cp, cooler_ln, cooler_ls

with tempfile.TemporaryDirectory() as folder:
    cool = f"{folder}/a.cool"
    multi = f"{folder}/m.mcool"
    cooler_cload_pairs("chrom.sizes:20", "pairs.txt", cool, chrom1=2, pos1=3, chrom2=4, pos2=5, verbose=False)
    cooler_cp(cool, f"{multi}::resolutions/20", verbose=False)
    result = cooler_ln(f"{multi}::resolutions/20", f"{multi}::alias", verbose=False)
    print(result["status"], [uri.split("::")[1] for uri in cooler_ls(multi, verbose=False)["coolers"]])
```

- **Example output:**

```text
completed ['/alias', '/resolutions/20']
```

## 19) `dnakit.epigenome.cooler_load()`

- **Purpose:** Create a cooler from a pre-binned matrix in COO or BG2 text format.
- **Method:** Run `cooler load --format`. `coo` rows are `bin1_id`, `bin2_id`, `count`; `bg2` rows are `chrom1`, `start1`, `end1`, `chrom2`, `start2`, `end2`, `count`. Once any `fields` entry is given, the count column is loaded only if it is listed too (`count=7` or `count:dtype=float`).
- **API:** `cooler_load(bins: 'str', pixels_path: 'str', output_path: 'str', *, input_format: "Literal['coo', 'bg2']", metadata_path: 'str | None' = None, assembly: 'str | None' = None, fields: 'Sequence[str] | None' = None, count_as_float: 'bool' = False, one_based: 'bool' = False, comment_char: 'str | None' = None, no_symmetric_upper: 'bool' = False, input_copy_status: "Literal['unique', 'duplex'] | None" = None, chunksize: 'int | None' = None, mergebuf: 'int | None' = None, temp_dir: 'str | None' = None, max_merge: 'int | None' = None, no_delete_temp: 'bool' = False, storage_options: 'str | None' = None, append: 'bool' = False, cooler_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package cooler (extra `dnakit[bio-cooler]`).
- **Input:** `bins` is `<chromsizes file>:<bin size in bp>` or a BED file. `pixels_path` is an existing text file, optionally gzipped (standard input is not supported). `input_format` is `"coo"` or `"bg2"`. `fields` entries look like `score=4:dtype=float`. `count_as_float` stores counts as floats; `one_based` says that the bin ids of a COO file are one-based. `output_path` must not exist unless `append=True`.
- **Output:** The common fields plus `output_uri`.
- **Example code:**

```python
import tempfile

from dnakit.epigenome import cooler_info, cooler_load

with tempfile.TemporaryDirectory() as folder:
    result = cooler_load("chrom.sizes:20", "pixels.coo", f"{folder}/l.cool", input_format="coo", verbose=False)
    info = cooler_info(f"{folder}/l.cool", verbose=False)["info"]
    print(result["status"], info["nbins"], info["nnz"], info["sum"])
```

- **Example output:**

```text
completed 8 3 8
```

## 20) `dnakit.epigenome.cooler_ls()`

- **Purpose:** List all coolers inside an HDF5 file, for example the resolutions of a `.mcool` file.
- **Method:** Run `cooler ls` (with `--long` when asked) and parse each line. With `long=True` the bin size of every cooler is reported, `None` for variable-size bins.
- **API:** `cooler_ls(cool_path: 'str', *, long: 'bool' = False, cooler_path: 'str | None' = None, timeout_seconds: 'float' = 600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package cooler (extra `dnakit[bio-cooler]`).
- **Input:** `cool_path` is an existing cooler or multi-resolution file (a plain path, not a `::group` URI).
- **Output:** The common fields plus `cool_path`, `coolers` (a list of `file::/group` URIs) and, with `long=True`, `bin_sizes` (URI to bin size).
- **Example code:**

```python
import tempfile

from dnakit.epigenome import cooler_cload_pairs, cooler_coarsen, cooler_cp, cooler_ls

with tempfile.TemporaryDirectory() as folder:
    cool = f"{folder}/a.cool"
    multi = f"{folder}/m.mcool"
    cooler_cload_pairs("chrom.sizes:20", "pairs.txt", cool, chrom1=2, pos1=3, chrom2=4, pos2=5, verbose=False)
    cooler_cp(cool, f"{multi}::resolutions/20", verbose=False)
    cooler_coarsen(cool, f"{multi}::resolutions/40", append=True, verbose=False)
    result = cooler_ls(multi, long=True, verbose=False)
    print({uri.split("::")[1]: size for uri, size in result["bin_sizes"].items()})
```

- **Example output:**

```text
{'/resolutions/20': 20, '/resolutions/40': 40}
```

## 21) `dnakit.epigenome.cooler_makebins()`

- **Purpose:** Write fixed-width genomic bins as a BED file from a chromsizes file.
- **Method:** Run `cooler makebins --out`. `rel_ids` adds a column of per-chromosome bin numbers, zero-based (0) or one-based (1). The result adds the line count and the first lines of the file.
- **API:** `cooler_makebins(chromsizes_path: 'str', bin_size: 'int', output_path: 'str', *, header: 'bool' = False, rel_ids: 'Literal[0, 1] | None' = None, cooler_path: 'str | None' = None, timeout_seconds: 'float' = 600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package cooler (extra `dnakit[bio-cooler]`).
- **Input:** `chromsizes_path` is an existing chromsizes file; `bin_size` is a positive integer in bp; `output_path` must not exist; `header` adds a header line.
- **Output:** The common fields plus `line_count` and `preview` (up to five lines).
- **Example code:**

```python
import tempfile

from dnakit.epigenome import cooler_makebins

with tempfile.TemporaryDirectory() as folder:
    result = cooler_makebins("chrom.sizes", 40, f"{folder}/bins.bed", rel_ids=1, verbose=False)
    print(result["status"], result["line_count"])
    for line in result["preview"]:
        print(" ".join(line.split()))
```

- **Example output:**

```text
completed 5
chr1 0 40 1
chr1 40 80 2
chr1 80 100 3
chr2 0 40 1
chr2 40 60 2
```

## 22) `dnakit.epigenome.cooler_merge()`

- **Purpose:** Merge coolers that share identical bins into one cooler by summing their pixels.
- **Method:** Run `cooler merge` with the output first and the inputs after it. `fields` are as in `cooler_coarsen()`.
- **API:** `cooler_merge(output_path: 'str', input_paths: 'Sequence[str]', *, chunksize: 'int | None' = None, fields: 'Sequence[str] | None' = None, append: 'bool' = False, cooler_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package cooler (extra `dnakit[bio-cooler]`).
- **Input:** `input_paths` is a list of at least two existing coolers (files or URIs) with identical chromosomes and bins. `chunksize` (merge buffer in pixel rows, cooler default 20,000,000) is a positive integer. `output_path` must not exist unless `append=True`.
- **Output:** The common fields plus `output_uri` and `input_uris`.
- **Example code:**

```python
import tempfile

from dnakit.epigenome import cooler_cload_pairs, cooler_cp, cooler_info, cooler_merge

with tempfile.TemporaryDirectory() as folder:
    cool = f"{folder}/a.cool"
    cooler_cload_pairs("chrom.sizes:20", "pairs.txt", cool, chrom1=2, pos1=3, chrom2=4, pos2=5, verbose=False)
    cooler_cp(cool, f"{folder}/b.cool", verbose=False)
    result = cooler_merge(f"{folder}/merged.cool", [cool, f"{folder}/b.cool"], verbose=False)
    info = cooler_info(f"{folder}/merged.cool", verbose=False)["info"]
    print(result["status"], len(result["input_uris"]), info["nnz"], info["sum"])
```

- **Example output:**

```text
completed 2 7 14
```

## 23) `dnakit.epigenome.cooler_mv()`

- **Purpose:** Rename a cooler group inside one HDF5 file.
- **Method:** Run `cooler mv`. `overwrite=True` passes `-w`, which cooler 0.10.4 cannot honour inside one file (see `cooler_cp()`); it is accepted for completeness.
- **API:** `cooler_mv(source_uri: 'str', destination_uri: 'str', *, overwrite: 'bool' = False, cooler_path: 'str | None' = None, timeout_seconds: 'float' = 600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package cooler (extra `dnakit[bio-cooler]`).
- **Input:** `source_uri` is an existing `file::/group`; `destination_uri` is a `file::/group` in the same file.
- **Output:** The common fields plus `source_uri`, `destination_uri` and `source_file`.
- **Example code:**

```python
import tempfile

from dnakit.epigenome import cooler_cload_pairs, cooler_cp, cooler_ls, cooler_mv

with tempfile.TemporaryDirectory() as folder:
    cool = f"{folder}/a.cool"
    multi = f"{folder}/m.mcool"
    cooler_cload_pairs("chrom.sizes:20", "pairs.txt", cool, chrom1=2, pos1=3, chrom2=4, pos2=5, verbose=False)
    cooler_cp(cool, f"{multi}::resolutions/20", verbose=False)
    result = cooler_mv(f"{multi}::resolutions/20", f"{multi}::resolutions/renamed", verbose=False)
    print(result["status"], [uri.split("::")[1] for uri in cooler_ls(multi, verbose=False)["coolers"]])
```

- **Example output:**

```text
completed ['/resolutions/renamed']
```

## 24) `dnakit.epigenome.cooler_pairs_header()`

- **Purpose:** Return the leading comment lines of a contact list, plain or gzipped.
- **Method:** Run in process, as `cooler cload` does when it reads the header. Lines are read until the first one that does not start with `comment_char`; they are stripped of trailing spaces and newlines and keep their comment character. A `.gz` file is opened with gzip. No cooler or other program is needed.
- **API:** `cooler_pairs_header(pairs_path: 'str', *, comment_char: 'str' = '#', verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Works after plain pip install.
- **Input:** `pairs_path` is an existing text or `.gz` file; `comment_char` is a short string without quotes or backslashes. An unreadable or non-UTF-8 file raises `ConfigurationError`.
- **Output:** `status`, `pairs_path`, `header` (list of lines), `header_line_count` and `has_data_lines`.
- **Example code:**

```python
from dnakit.epigenome import cooler_pairs_header

result = cooler_pairs_header("pairs.txt", verbose=False)
print(result["header_line_count"], result["has_data_lines"])
for line in result["header"]:
    print(line)
```

- **Example output:**

```text
2 True
## pairs format v1.0
#columns: readID chrom1 pos1 chrom2 pos2
```

## 25) `dnakit.epigenome.cooler_show()`

- **Purpose:** Draw a contact-matrix region of a cooler to an image file without a display.
- **Method:** Run `cooler show --out` with `MPLBACKEND=Agg`. When the cooler package is importable, DNAKit then reads the same region with `Cooler.matrix(...).fetch()` and returns the values behind the image: a summary and, for at most 40,000 cells, the matrix itself (non-finite values become `None`). Above 4,000,000 cells the values are omitted. The colour scale is applied to the image only; the returned matrix is not transformed.
- **API:** `cooler_show(cool_uri: 'str', region: 'str', output_path: 'str', *, region2: 'str | None' = None, balanced: 'bool' = False, dpi: 'int | None' = None, scale: "Literal['linear', 'log2', 'log10'] | None" = None, force: 'bool' = False, zmin: 'float | None' = None, zmax: 'float | None' = None, cmap: 'str | None' = None, field: 'str | None' = None, cooler_path: 'str | None' = None, timeout_seconds: 'float' = 1800.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package cooler (extra `dnakit[bio-cooler]`) and matplotlib (extra `dnakit[bio-plots]`).
- **Input:** `region` is a UCSC-style region without whitespace (`chr1:10,000,000-11,000,000`); `region2` gives the column range for trans or asymmetric views. `output_path` must not exist and must end in `.png`, `.jpg`, `.jpeg`, `.svg`, `.pdf`, `.ps` or `.eps`. `scale` is `linear`, `log2` or `log10` (cooler default); `dpi` is at least 10; `cmap` is a matplotlib colormap name (cooler default `YlOrRd`); `field` is the pixel column (default `count`); `zmin` and `zmax` are finite numbers.
- **Output:** The common fields plus `image_path`, `image_bytes`, `cool_uri`, `region`, `region2`, `balanced`, `scale`, `field` and `data` (`binsize`, `summary` with `shape`, `finite_cells`, `nan_cells`, `min`, `max`, `sum`, `mean`, and `matrix`).
- **Example code:**

```python
import tempfile

from dnakit.epigenome import cooler_cload_pairs, cooler_show

with tempfile.TemporaryDirectory() as folder:
    cool = f"{folder}/a.cool"
    cooler_cload_pairs("chrom.sizes:20", "pairs.txt", cool, chrom1=2, pos1=3, chrom2=4, pos2=5, verbose=False)
    result = cooler_show(cool, "chr1", f"{folder}/chr1.png", verbose=False)
    print(result["status"], result["image_bytes"] > 0, result["scale"], result["data"]["summary"]["shape"])
    print(result["data"]["matrix"][0])
    trans = cooler_show(cool, "chr1", f"{folder}/trans.png", region2="chr2", scale="linear", verbose=False)
    print(trans["data"]["summary"]["shape"], trans["data"]["summary"]["sum"])
```

- **Example output:**

```text
completed True log10 [5, 5]
[0.0, 1.0, 1.0, 0.0, 1.0]
[5, 3] 2.0
```

## 26) `dnakit.epigenome.cooler_show_view()`

- **Purpose:** Draw the heat map that the interactive `cooler show` window displays after a zoom or pan, without a display.
- **Method:** Run in process with the cooler and matplotlib packages (no `cooler` program is started). Like the interactive window, it rounds the window outwards to whole bins, trims it to the chromosome (the last, partial bin is excluded), and draws a placeholder text instead of the matrix when the window holds `max_matrix_size` cells or more. Without `view` the whole of `region` is drawn.
- **API:** `cooler_show_view(cool_uri: 'str', region: 'str', output_path: 'str', *, region2: 'str | None' = None, view: 'Sequence[float] | None' = None, balanced: 'bool' = False, field: 'str' = 'count', scale: "Literal['linear', 'log2', 'log10']" = 'log10', cmap: 'str' = 'YlOrRd', zmin: 'float | None' = None, zmax: 'float | None' = None, dpi: 'int | None' = None, max_matrix_size: 'int' = 10000000, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages cooler (extra `dnakit[bio-cooler]`), numpy (`dnakit[bio-analysis]`) and matplotlib (`dnakit[bio-plots]`). The cooler must have fixed-size bins.
- **Input:** `region` and `region2` as in `cooler_show()`. `view` is `[column start, column end, row start, row end]` in bp on the chromosomes of `region2` (default `region`) and `region`, with start below end. `balanced=True` needs a `weight` column (run `cooler_balance()` first). `output_path` must not exist and needs an image suffix. `max_matrix_size` is a positive integer.
- **Output:** `status`, `backend`, `image_path`, `image_bytes`, `cool_uri`, `row_region` and `column_region` (`[chrom, start, end]` after rounding), `view_applied`, `placeholder`, `binsize`, `balanced`, `scale`, `field` and `data` (`summary`, `matrix`, `matrix_omitted`).
- **Example code:**

```python
import tempfile

from dnakit.epigenome import cooler_cload_pairs, cooler_show_view

with tempfile.TemporaryDirectory() as folder:
    cool = f"{folder}/a.cool"
    cooler_cload_pairs("chrom.sizes:20", "pairs.txt", cool, chrom1=2, pos1=3, chrom2=4, pos2=5, verbose=False)
    result = cooler_show_view(cool, "chr1", f"{folder}/view.png", view=[10, 70, 10, 70], verbose=False)
    print(result["row_region"], result["column_region"], result["view_applied"], result["placeholder"])
    print(result["data"]["matrix"])
```

- **Example output:**

```text
['chr1', 0, 80] ['chr1', 0, 80] True False
[[0.0, 1.0, 1.0, 0.0], [1.0, 0.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0]]
```

## 27) `dnakit.epigenome.cooler_tree()`

- **Purpose:** Return the HDF5 data hierarchy of a cooler file or group as text.
- **Method:** Run `cooler tree` with an optional `--level`. The text cooler prints (`cooler.fileops.pprint_data_tree`) is returned in `tree`.
- **API:** `cooler_tree(uri: 'str', *, level: 'int | None' = None, cooler_path: 'str | None' = None, timeout_seconds: 'float' = 600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package cooler (extra `dnakit[bio-cooler]`).
- **Input:** `uri` is an existing cooler file or `file::/group`. `level` is the maximum depth (an integer of at least 0).
- **Output:** The common fields plus `uri` and `tree` (`None` if the run failed).
- **Example code:**

```python
import tempfile

from dnakit.epigenome import cooler_cload_pairs, cooler_tree

with tempfile.TemporaryDirectory() as folder:
    cool = f"{folder}/a.cool"
    cooler_cload_pairs("chrom.sizes:20", "pairs.txt", cool, chrom1=2, pos1=3, chrom2=4, pos2=5, verbose=False)
    print(cooler_tree(cool, level=1, verbose=False)["tree"])
```

- **Example output:**

```text
/
 ├── bins
 ├── chroms
 ├── indexes
 └── pixels
```

## 28) `dnakit.epigenome.cooler_zoomify()`

- **Purpose:** Build a multi-resolution cooler (`.mcool`) by recursive coarsening.
- **Method:** Run `cooler zoomify --out`. `resolutions` is a comma-separated list such as `1000,5000,25000` or a progression: `B` (binary, steps of 2), `N` (nice, 1-2-5 steps), `<start>N` or `4DN`; the default is `B`. The progression `<start>B` (for example `1000B`) is documented by cooler but fails in 0.10.4, so use `B` or an explicit list. `balance=True` balances every zoom level and `balance_options` passes cooler balance settings.
- **API:** `cooler_zoomify(cool_uri: 'str', output_path: 'str', *, resolutions: 'str | None' = None, nproc: 'int | None' = None, chunksize: 'int | None' = None, balance: 'bool' = False, balance_options: 'Mapping[str, Any] | None' = None, base_uris: 'Sequence[str] | None' = None, fields: 'Sequence[str] | None' = None, legacy: 'bool' = False, cooler_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package cooler (extra `dnakit[bio-cooler]`).
- **Input:** `cool_uri` is an existing cooler; `output_path` must not exist. `balance_options` needs `balance=True` and may hold `cis_only`, `trans_only`, `ignore_diags`, `ignore_dist`, `mad_max`, `min_nnz`, `min_count`, `nproc`, `max_iters` and `tol`. `base_uris` are extra base coolers to aggregate from. `legacy=True` writes the integer-labelled zoom levels of the legacy layout. `nproc` and `chunksize` are positive integers.
- **Output:** The common fields plus `input_uri`, `output_path`, `legacy` and `balanced`.
- **Example code:**

```python
import tempfile

from dnakit.epigenome import cooler_cload_pairs, cooler_ls, cooler_zoomify

with tempfile.TemporaryDirectory() as folder:
    cool = f"{folder}/a.cool"
    cooler_cload_pairs("chrom.sizes:20", "pairs.txt", cool, chrom1=2, pos1=3, chrom2=4, pos2=5, verbose=False)
    result = cooler_zoomify(cool, f"{folder}/z.mcool", resolutions="20,40,80", verbose=False)
    listing = cooler_ls(f"{folder}/z.mcool", long=True, verbose=False)
    print(result["status"], result["legacy"], result["balanced"])
    print({uri.split("::")[1]: size for uri, size in listing["bin_sizes"].items()})
```

- **Example output:**

```text
completed False False
{'/resolutions/20': 20, '/resolutions/40': 40, '/resolutions/80': 80}
```

## 29) `dnakit.epigenome.cooler_resolve_bins()`

- **Purpose:** Resolve a cooler `BINS` argument and write the bin table to a new BED file.
- **Method:** Follows `cooler.cli._util.parse_bins` (cooler 0.10.4). A string with a colon is `<chromsizes file>:<binsize in bp>` and gives fixed-size bins from `cooler.util.binnify`; otherwise it is the path of a three-column BED file (chrom, start, end) that defines the bins. A BED file that fails to parse raises `ConfigurationError` (upstream raises `AttributeError` there with pandas 2). The table is written without a header.
- **API:** `cooler_resolve_bins(bins: 'str', output_bed: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package cooler (extra `dnakit[bio-cooler]`); pandas is a cooler dependency.
- **Input:** `bins` is a non-empty string as above, with an existing chromsizes or BED file and a positive integer bin size. `output_bed` must not exist.
- **Output:** `operation`, `status`, `backend`, `package`, `source` (`"chromsizes"` or `"bed"`), `binsize` (`None` for a BED source), `bin_count`, `chromosomes` (in order of first appearance), `output_bed` and `coordinates` (`"0-based half-open BED"`).
- **Example code:**

```python
import tempfile

from dnakit.epigenome import cooler_resolve_bins

with tempfile.TemporaryDirectory() as folder:
    result = cooler_resolve_bins("chrom.sizes:40", f"{folder}/bins.bed", verbose=False)
    print(result["source"], result["binsize"], result["bin_count"], result["chromosomes"])
    print(open(f"{folder}/bins.bed").read().split("\n")[0].split("\t"))
    again = cooler_resolve_bins(f"{folder}/bins.bed", f"{folder}/bins2.bed", verbose=False)
    print(again["source"], again["binsize"], again["bin_count"])
```

- **Example output:**

```text
chromsizes 40 5 ['chr1', 'chr2']
['chr1', '0', '40']
bed None 5
```

## 30) `dnakit.epigenome.download_cooltools_example_data()`

- **Purpose:** Download one cooltools example dataset, or all of them, into a named folder.
- **Method:** Call `cooltools.lib.io.download_data()`. With `cache=True` a file that already exists in `data_dir` is reused without any download, after its md5 checksum is compared with the one in the dataset table (unless `ignore_checksum` is true). Otherwise the file is fetched over HTTP(S).
- **API:** `download_cooltools_example_data(name: 'str' = 'all', *, data_dir: 'str', cache: 'bool' = True, ignore_checksum: 'bool' = False, info_url: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package cooltools (extra `dnakit[bio-cooltools]`) and internet access for files that are not cached; the example uses a local table and a cached file, so it needs no network.
- **Input:** `name` is a dataset key or `"all"`; `data_dir` is a folder (created if missing). `cache` and `ignore_checksum` are booleans. `info_url` is an http(s) URL or the path of a local dataset table (by default cooltools' published table is used). A table has a header line starting with `# ` and tab-separated rows of key, comment, link, filename and md5 checksum.
- **Output:** `status` (`"complete"`), `operation`, `service`, `backend`, `version`, `request`, `retrieved_utc` and `data` with `files` (paths), `file_count` and `data_dir`. An unknown dataset or a checksum mismatch raises `BackendExecutionError`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import download_cooltools_example_data
from dnakit.exceptions import BackendExecutionError

with tempfile.TemporaryDirectory() as folder:
    Path(folder, "tiny.txt").write_text("hello\n")
    result = download_cooltools_example_data("tiny", data_dir=folder, info_url="datasets.tsv", verbose=False)
    print(result["status"], result["data"]["file_count"], [Path(path).name for path in result["data"]["files"]])
    try:
        download_cooltools_example_data("missing", data_dir=folder, info_url="datasets.tsv", verbose=False)
    except BackendExecutionError as error:
        print(type(error).__name__, "is not available" in str(error))
```

- **Example output:**

```text
complete 1 ['tiny.txt']
BackendExecutionError True
```

Usage with internet access:

```python
from dnakit.epigenome import download_cooltools_example_data

result = download_cooltools_example_data("HFF_MicroC", data_dir="cooltools_data")
print(result["data"]["files"])
```

## 31) `dnakit.epigenome.download_cooltools_file()`

- **Purpose:** Download one file from an http(s) URL to a new local path.
- **Method:** Call `cooltools.lib.io.download_file()` (an HTTP GET streamed to the file), then report the size and md5 checksum of the file. A run that writes no data raises `BackendExecutionError`.
- **API:** `download_cooltools_file(url: 'str', output_path: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package cooltools (extra `dnakit[bio-cooltools]`) and internet access. The example needs no network: a URL with another scheme is rejected before anything is downloaded.
- **Input:** `url` is an http or https URL (at most 2,048 characters); `output_path` must not exist (its folder is created if missing).
- **Output:** The same envelope as above (`status`, `operation`, `service`, `backend`, `version`, `request`, `retrieved_utc`) and `data` with `output_path`, `bytes` and `md5`.
- **Example code:**

```python
from dnakit.epigenome import download_cooltools_file
from dnakit.exceptions import ConfigurationError

try:
    download_cooltools_file("ftp://example.org/data.cool", "data.cool", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, "url must be an http or https URL" in str(error))
```

- **Example output:**

```text
ConfigurationError True
```

Usage with internet access:

```python
from dnakit.epigenome import download_cooltools_file

result = download_cooltools_file("https://example.org/data.txt", "data.txt")
print(result["data"]["bytes"], result["data"]["md5"])
```

## 32) `dnakit.epigenome.list_cooltools_example_datasets()`

- **Purpose:** Return the cooltools example-dataset table that `print_available_datasets` prints.
- **Method:** Read the table with cooltools' own parser (`cooltools.lib.io._get_datasets_info`) from `url`, or from cooltools' published table when `url` is `None`. Each row holds the dataset key, comment, source link, local file name and md5 checksum, plus the row `index`.
- **API:** `list_cooltools_example_datasets(*, url: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package cooltools (extra `dnakit[bio-cooltools]`). Without `url` it needs internet access; with a local table it works offline.
- **Input:** `url` is an http(s) URL or the path of an existing local table with a `# ` header line and tab-separated rows.
- **Output:** `status` (`"complete"`), `operation`, `service`, `backend`, `version`, `request`, `retrieved_utc` and `data` with `datasets` (list of row dictionaries) and `count`.
- **Example code:**

```python
from dnakit.epigenome import list_cooltools_example_datasets

result = list_cooltools_example_datasets(url="datasets.tsv", verbose=False)
print(result["status"], result["operation"], result["data"]["count"])
for row in result["data"]["datasets"]:
    print(row["index"], row["key"], row["filename"], row["checksum"][:8])
```

- **Example output:**

```text
complete list_cooltools_example_datasets 2
1 tiny tiny.txt b1946ac9
2 small small.txt 00000000
```

## 33) `dnakit.epigenome.plot_cooltools_saddle()`

- **Purpose:** Draw a cooltools saddle plot of a track against a saddle matrix into a new PNG file.
- **Method:** Digitise the track with `cooltools.api.saddle.digitize`, call `cooltools.api.saddle.saddleplot` (which draws on pyplot's current figure), save that figure at 150 dpi and close every figure the call opened. The returned counts come from the digitised track; the heatmap values are the saddle matrix without its two flank rows and columns.
- **API:** `plot_cooltools_saddle(track_rows: 'Sequence[Sequence[object]]', saddle_matrix: 'Sequence[Sequence[object]]', n_bins: 'int', output_png: 'str', *, vrange: 'Sequence[float] | None' = None, qrange: 'Sequence[float] | None' = (0.0, 1.0), scale: 'str' = 'log', vmin: 'float' = 0.5, vmax: 'float' = 2.0, cmap: 'str' = 'coolwarm', title: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages cooltools (extra `dnakit[bio-cooltools]`) and matplotlib (extra `dnakit[bio-plots]`).
- **Input:** `track_rows` is a non-empty list of `[chrom, start, end, value]` rows (non-negative integer positions with end above start; `None` for a missing value). `saddle_matrix` is a square list of `n_bins + 2` rows, including the two outlier flanks that `cooltools.api.saddle` returns. Exactly one of `vrange` (value limits) and `qrange` (quantile limits, default `(0.0, 1.0)`) must be given, so pass `qrange=None` when using `vrange`. `scale` is `"log"` or `"linear"`; `vmin` must be below `vmax` (and above 0 for `log`). `output_png` must not exist and end in `.png`.
- **Output:** `status`, `operation`, `service`, `backend`, `version`, `image_path`, `data` and `upstream_warnings`. `data` holds `n_bins`, `vrange`, `qrange`, `scale`, `vmin`, `vmax`, `cmap`, `track_rows`, `bin_edges`, `bin_counts`, `left_outlier_count`, `right_outlier_count`, `missing_count`, `heatmap_shape` and `heatmap_values`.
- **Example code:**

```python
import tempfile

from dnakit.epigenome import plot_cooltools_saddle

track = [["chr1", 0, 20, 0.1], ["chr1", 20, 40, 0.5], ["chr1", 40, 60, -0.3],
         ["chr1", 60, 80, 0.9], ["chr1", 80, 100, -0.8], ["chr2", 0, 20, 0.2]]
saddle = [[0] * 5, [0, 1.2, 0.8, 0.5, 0], [0, 0.8, 1.0, 0.9, 0], [0, 0.5, 0.9, 1.5, 0], [0] * 5]
with tempfile.TemporaryDirectory() as folder:
    result = plot_cooltools_saddle(track, saddle, 3, f"{folder}/saddle.png", vrange=[-1, 1], qrange=None, verbose=False)
    data = result["data"]
    print(result["status"], data["bin_counts"], data["left_outlier_count"], data["right_outlier_count"])
    print([round(edge, 3) for edge in data["bin_edges"]], data["heatmap_shape"])
    print(data["heatmap_values"][0])
```

- **Example output:**

```text
completed [1, 3, 2] 0 0
[-1.0, -0.333, 0.333, 1.0] [3, 3]
[1.2, 0.8, 0.5]
```
