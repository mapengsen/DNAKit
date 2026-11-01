# SeqKit and screed: FASTA and FASTQ command-line operations

This group of functions in `dnakit.fasta_fastq` wraps two upstream tools for sequence files: the external command-line program [SeqKit](https://github.com/shenwei356/seqkit) (21 functions, in `seqkit_cli_functions`) and the Python package [screed](https://github.com/dib-lab/screed) (3 functions, in `screed_gap_functions`). The SeqKit functions run one whitelisted subcommand (`subseq`, `stats`, `sum`, `watch`, `fa2fq`, `fx2tab`, `fish`, `grep`, `duplicate`, `range`, `head-genome`, `split`, `concat`, `rename`, `replace`, `bam`, `merge-slides`) without a shell, with a timeout and an output size limit, and write results to explicit new paths; 2.14.0 is the tested SeqKit version. The screed functions index a FASTA or FASTQ file into a screed database and dump it back to text. Import them with `from dnakit.fasta_fastq import ...`.

The dependencies are optional and are imported or looked up only when a function is called. SeqKit is an external program: install it (for example `conda install -c bioconda seqkit` or a release binary), then put it on PATH or pass `seqkit_path`. If it cannot be found, `BackendUnavailableError` (code `SEQKIT_UNAVAILABLE`) is raised with the install hint; if SeqKit exits with an error, `BackendExecutionError` (code `SEQKIT_EXECUTION_FAILED`) is raised with the tail of its log. screed needs the extra `dnakit[bio-screed]`; its package metadata lists the licence as BSD-3-Clause. DNAKit does not bundle SeqKit.

Every SeqKit function validates its arguments first, never overwrites an existing output and returns a dictionary. All SeqKit results share the same fields: `status` (`completed`), `backend` (`seqkit`), `operation`, `version`, `command`, `outputs`, `log_tail` and `provenance`; the section for each function lists only the extra fields. SeqKit is not installed in the environment used to produce this page, so the 21 SeqKit examples show the validation or missing-executable behaviour that is real on any machine; each of them is followed by a call to run with SeqKit installed. The three screed examples run for real.

## Example data

The SeqKit examples create their own tiny inputs in a temporary directory. The screed examples read these two files.

**Example file:** `reads.fa`

```text
>r1 first read
ACGTACGT
>r2 second read
GGCCTTAA
```

**Example file:** `reads.fq`

```text
@r1 first read
ACGTACGT
+
IIIIHHHH
@r2 second read
GGCCTTAA
+
########
```

## 1) `dnakit.fasta_fastq.seqkit_bam_alignment_stats()`

- **Purpose:** Count mapped, unmapped, primary, secondary and supplementary alignments of a BAM file.
- **Method:** Run `seqkit bam --stat` (full scan) or `seqkit bam --idx-stat` (when `index_based` is true) with `--quiet-mode`, and parse the tab-separated table that SeqKit prints to stderr; the first column is `PrimAlnPerc` in scan mode and `AlnPerc` in index mode.
- **API:** `seqkit_bam_alignment_stats(bam_path: 'str', *, index_based: 'bool' = False, seqkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program `seqkit` (2.14.0 is the tested version) on PATH, or an explicit `seqkit_path`; it is not installed in the environment that produced this page, and no Python package is required. The example below therefore shows behaviour that does not need SeqKit: a `seqkit_path` that does not exist raises `BackendUnavailableError` (code `SEQKIT_UNAVAILABLE`).
- **Input:** `bam_path` is an existing BAM file. `index_based` must be a bool; True reads the BAM index and needs a `.bai` next to the BAM.
- **Output:** The common result dictionary plus `records` (the parsed table rows, numbers converted, the `File` column kept as text) and `index_based`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.fasta_fastq import seqkit_bam_alignment_stats

with tempfile.TemporaryDirectory() as tmp:
    bam = Path(tmp) / "reads.bam"
    bam.write_bytes(b"not a real BAM")
    try:
        seqkit_bam_alignment_stats(str(bam), index_based=True, seqkit_path="/nonexistent/seqkit", verbose=False)
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError SEQKIT_UNAVAILABLE
```

**Usage with the dependency installed:**

```python
from dnakit.fasta_fastq import seqkit_bam_alignment_stats

result = seqkit_bam_alignment_stats("reads.bam", verbose=False)
print(result["records"])
```

## 2) `dnakit.fasta_fastq.seqkit_bam_bundles()`

- **Purpose:** Partition a BAM file into loci or into bundles of at least a given number of reads.
- **Method:** Create `output_dir`, run `seqkit bam --bundle=<min_reads> --quiet-mode` with that directory as the working directory, parse the `Bundle` table from SeqKit's log and read the "Written N BAM records to L loci and B bundles" summary line. SeqKit writes the bundles to `<bam name>_bundles/` inside `output_dir`.
- **API:** `seqkit_bam_bundles(bam_path: 'str', output_dir: 'str', *, min_reads: 'int' = -1, seqkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program `seqkit` (2.14.0 is the tested version) on PATH, or an explicit `seqkit_path`; it is not installed in the environment that produced this page, and no Python package is required. The example below therefore shows behaviour that does not need SeqKit: a `seqkit_path` that does not exist raises `BackendUnavailableError` (code `SEQKIT_UNAVAILABLE`), and the output directory has already been created by then.
- **Input:** `bam_path` is an existing BAM file. `output_dir` must not exist yet and its parent must exist. `min_reads` is an integer of at least -1; -1 (default) means one bundle per locus.
- **Output:** The common result dictionary plus `records` (the bundle table rows, `Chrom` kept as text), `min_reads`, and, when the summary line is found, `record_count`, `locus_count` and `bundle_count`. `outputs` holds `directory` and `bundles`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.fasta_fastq import seqkit_bam_bundles

with tempfile.TemporaryDirectory() as tmp:
    bam = Path(tmp) / "reads.bam"
    bam.write_bytes(b"not a real BAM")
    out = Path(tmp) / "bundles"
    try:
        seqkit_bam_bundles(str(bam), str(out), min_reads=10, seqkit_path="/nonexistent/seqkit", verbose=False)
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
    print("output directory created before the failure:", out.is_dir())
    try:
        seqkit_bam_bundles(str(bam), str(out), seqkit_path="/nonexistent/seqkit", verbose=False)
    except Exception as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError SEQKIT_UNAVAILABLE
output directory created before the failure: True
ConfigurationError CONFIGURATION_ERROR
```

**Usage with the dependency installed:**

```python
from dnakit.fasta_fastq import seqkit_bam_bundles

result = seqkit_bam_bundles("reads.bam", "bundles_out", min_reads=100, verbose=False)
print(result["record_count"], result["locus_count"], result["bundle_count"])
```

## 3) `dnakit.fasta_fastq.seqkit_bam_feature_histogram()`

- **Purpose:** Build a histogram of one per-record BAM feature such as `ReadLen`, `MapQual` or `MeanQual`.
- **Method:** Run `seqkit bam --field=<field> --bins=<bins> --dump --quiet-mode` (plus `--log` and `--img=<image_path>` when requested) and parse the `Bin` table that SeqKit prints. The field name is only checked for shape; SeqKit decides whether it is known (`seqkit bam -H` lists them). `Acc` is not passed because SeqKit 2.14.0 panics on a BAM without an NM tag, per the source comment.
- **API:** `seqkit_bam_feature_histogram(bam_path: 'str', field: 'str', *, bins: 'int' = -1, log_transform: 'bool' = False, image_path: 'str | None' = None, seqkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program `seqkit` (2.14.0 is the tested version) on PATH, or an explicit `seqkit_path`; it is not installed in the environment that produced this page, and no Python package is required. The example below therefore shows behaviour that does not need SeqKit: an invalid argument raises `ConfigurationError` before SeqKit is looked up.
- **Input:** `bam_path` is an existing BAM file. `field` is a name made of letters and digits that starts with a letter. `bins` is an integer of at least -1 (-1 lets SeqKit choose). `log_transform` is a bool. `image_path` is optional and must be a new file in an existing folder.
- **Output:** The common result dictionary plus `records` (the histogram bins) and `field`. `outputs` holds `image` when an image was requested.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.fasta_fastq import seqkit_bam_feature_histogram

with tempfile.TemporaryDirectory() as tmp:
    bam = Path(tmp) / "reads.bam"
    bam.write_bytes(b"not a real BAM")
    try:
        seqkit_bam_feature_histogram(str(bam), "Read Len", verbose=False)
    except ConfigurationError as error:
        print(type(error).__name__, error)
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] field must be a SeqKit BAM feature name such as ReadLen.
```

**Usage with the dependency installed:**

```python
from dnakit.fasta_fastq import seqkit_bam_feature_histogram

result = seqkit_bam_feature_histogram("reads.bam", "ReadLen", bins=10, verbose=False)
print(result["field"], len(result["records"]))
```

## 4) `dnakit.fasta_fastq.seqkit_bam_top_records()`

- **Purpose:** Save the BAM records with the largest value of a feature to a new BAM file.
- **Method:** Run `seqkit bam --field=<field> --top-bam=<output_bam> --top-size=<top_size> --quiet-mode`. The facade checks that the new BAM exists afterwards.
- **API:** `seqkit_bam_top_records(bam_path: 'str', output_bam: 'str', *, field: 'str', top_size: 'int' = 100, seqkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program `seqkit` (2.14.0 is the tested version) on PATH, or an explicit `seqkit_path`; it is not installed in the environment that produced this page, and no Python package is required. The example below therefore shows behaviour that does not need SeqKit: an invalid argument raises `ConfigurationError` before SeqKit is looked up.
- **Input:** `bam_path` is an existing BAM file. `output_bam` must not exist and its folder must exist. `field` is a feature name (letters and digits, starting with a letter) and is required. `top_size` is an integer of at least 1 (default 100).
- **Output:** The common result dictionary plus `field` and `top_size`. `outputs` holds `bam`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.fasta_fastq import seqkit_bam_top_records

with tempfile.TemporaryDirectory() as tmp:
    bam = Path(tmp) / "reads.bam"
    bam.write_bytes(b"not a real BAM")
    try:
        seqkit_bam_top_records(str(bam), str(Path(tmp) / "top.bam"), field="ReadLen", top_size=0, verbose=False)
    except ConfigurationError as error:
        print(type(error).__name__, error)
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] top_size must be at least 1.
```

**Usage with the dependency installed:**

```python
from dnakit.fasta_fastq import seqkit_bam_top_records

result = seqkit_bam_top_records("reads.bam", "longest.bam", field="ReadLen", top_size=10, verbose=False)
print(result["outputs"]["bam"].endswith("longest.bam"))
```

## 5) `dnakit.fasta_fastq.seqkit_concat_by_id()`

- **Purpose:** Concatenate the sequences that share an ID across several files, in file order.
- **Method:** Run `seqkit concat` with `--separator=<separator>` and, when `full_join` is true, `--full` (and `--fill=<fill>`), writing to `output_path` with `-o`. The record count of the output is counted afterwards.
- **API:** `seqkit_concat_by_id(input_paths: 'Sequence[str]', output_path: 'str', *, full_join: 'bool' = False, fill: 'str | None' = None, separator: 'str' = '|', seqkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program `seqkit` (2.14.0 is the tested version) on PATH, or an explicit `seqkit_path`; it is not installed in the environment that produced this page, and no Python package is required. The example below therefore shows behaviour that does not need SeqKit: an invalid argument raises `ConfigurationError` before SeqKit is looked up.
- **Input:** `input_paths` is a non-empty list of existing files. `output_path` must not exist. `full_join` is a bool. `fill` is optional, allowed only with `full_join=True`, and may contain only letters or `-`. `separator` is a string used to join the descriptions (default `|`).
- **Output:** The common result dictionary plus `record_count`. `outputs` holds `output`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.fasta_fastq import seqkit_concat_by_id

with tempfile.TemporaryDirectory() as tmp:
    a = Path(tmp) / "a.fa"
    b = Path(tmp) / "b.fa"
    a.write_text(">s1\nAAAA\n")
    b.write_text(">s1\nCCCC\n")
    try:
        seqkit_concat_by_id([str(a), str(b)], str(Path(tmp) / "joined.fa"), fill="N", verbose=False)
    except ConfigurationError as error:
        print(type(error).__name__, error)
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] fill only applies with full_join=True.
```

**Usage with the dependency installed:**

```python
from dnakit.fasta_fastq import seqkit_concat_by_id

result = seqkit_concat_by_id(["part1.fa", "part2.fa"], "joined.fa", full_join=True, fill="N", verbose=False)
print(result["record_count"])
```

## 6) `dnakit.fasta_fastq.seqkit_count_matching_records()`

- **Purpose:** Count the records whose ID, full name or sequence matches any of the given patterns.
- **Method:** Run `seqkit grep --count` with one `--pattern=` per pattern, plus `--by-name`, `--by-seq`, `--use-regexp`, `--invert-match` and `--ignore-case` as requested, and read the last all-digit line of the output as the count. Nothing is written to disk.
- **API:** `seqkit_count_matching_records(input_path: 'str', *, patterns: 'Sequence[str]', by_name: 'bool' = False, by_sequence: 'bool' = False, use_regexp: 'bool' = False, invert_match: 'bool' = False, ignore_case: 'bool' = False, seqkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program `seqkit` (2.14.0 is the tested version) on PATH, or an explicit `seqkit_path`; it is not installed in the environment that produced this page, and no Python package is required. The example below therefore shows behaviour that does not need SeqKit: an invalid argument raises `ConfigurationError` before SeqKit is looked up.
- **Input:** `input_path` is an existing FASTA/FASTQ file. `patterns` is a non-empty list of non-empty strings without commas (SeqKit splits pattern lists on commas). The other options are bools.
- **Output:** The common result dictionary plus `count` (an int). `outputs` is empty. If SeqKit prints no count, `BackendExecutionError` with code `SEQKIT_OUTPUT_UNREADABLE` is raised.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.fasta_fastq import seqkit_count_matching_records

with tempfile.TemporaryDirectory() as tmp:
    fasta = Path(tmp) / "reads.fa"
    fasta.write_text(">s1\nAAAA\n>s2\nCCCC\n")
    try:
        seqkit_count_matching_records(str(fasta), patterns=["s1,s2"], verbose=False)
    except ConfigurationError as error:
        print(type(error).__name__, error)
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] each pattern must be a non-empty string without commas (SeqKit splits pattern lists on commas).
```

**Usage with the dependency installed:**

```python
from dnakit.fasta_fastq import seqkit_count_matching_records

result = seqkit_count_matching_records("reads.fa", patterns=["r1", "r2"], verbose=False)
print(result["count"])
```

## 7) `dnakit.fasta_fastq.seqkit_duplicate_records()`

- **Purpose:** Write every record several times, next to each other, with unchanged IDs.
- **Method:** Run `seqkit duplicate --times=<times>` and write to `output_path` with `-o`; the record count of the output is counted afterwards.
- **API:** `seqkit_duplicate_records(input_path: 'str', output_path: 'str', *, times: 'int', seqkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program `seqkit` (2.14.0 is the tested version) on PATH, or an explicit `seqkit_path`; it is not installed in the environment that produced this page, and no Python package is required. The example below therefore shows behaviour that does not need SeqKit: an invalid argument raises `ConfigurationError` before SeqKit is looked up.
- **Input:** `input_path` is an existing FASTA/FASTQ file. `output_path` must not exist. `times` is a required integer of at least 1.
- **Output:** The common result dictionary plus `record_count`. `outputs` holds `output`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.fasta_fastq import seqkit_duplicate_records

with tempfile.TemporaryDirectory() as tmp:
    fasta = Path(tmp) / "reads.fa"
    fasta.write_text(">s1\nAAAA\n")
    try:
        seqkit_duplicate_records(str(fasta), str(Path(tmp) / "dup.fa"), times=0, verbose=False)
    except ConfigurationError as error:
        print(type(error).__name__, error)
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] times must be at least 1.
```

**Usage with the dependency installed:**

```python
from dnakit.fasta_fastq import seqkit_duplicate_records

result = seqkit_duplicate_records("reads.fa", "tripled.fa", times=3, verbose=False)
print(result["record_count"])
```

## 8) `dnakit.fasta_fastq.seqkit_extended_stats()`

- **Purpose:** Compute extended statistics for each file: length quartiles, N50 and Nx, Q20/Q30, average quality and GC.
- **Method:** Run `seqkit stats --all --tabular --fq-encoding=<fastq_encoding>` (plus one `--N=<value>` per entry of `n_values`) and write the table to `output_path`; then re-read that file, take the first line as the header and convert numeric cells.
- **API:** `seqkit_extended_stats(input_paths: 'Sequence[str]', output_path: 'str', *, fastq_encoding: 'str' = 'sanger', n_values: 'Sequence[int]' = (), seqkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program `seqkit` (2.14.0 is the tested version) on PATH, or an explicit `seqkit_path`; it is not installed in the environment that produced this page, and no Python package is required. The example below therefore shows behaviour that does not need SeqKit: an invalid argument raises `ConfigurationError` before SeqKit is looked up.
- **Input:** `input_paths` is a non-empty list of existing files. `output_path` must not exist. `fastq_encoding` is one of `sanger`, `solexa`, `illumina-1.3+`, `illumina-1.5+`, `illumina-1.8+`. `n_values` is a sequence of integers from 0 to 100.
- **Output:** The common result dictionary plus `records`, one dictionary per input file with the table columns (`file`, `format` and `type` are kept as text). `outputs` holds `output`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.fasta_fastq import seqkit_extended_stats

with tempfile.TemporaryDirectory() as tmp:
    fastq = Path(tmp) / "reads.fq"
    fastq.write_text("@r1\nACGT\n+\nIIII\n")
    try:
        seqkit_extended_stats([str(fastq)], str(Path(tmp) / "stats.tsv"), fastq_encoding="phred64", verbose=False)
    except ConfigurationError as error:
        print(type(error).__name__, error)
    try:
        seqkit_extended_stats([str(fastq)], str(Path(tmp) / "stats.tsv"), n_values=[101], verbose=False)
    except ConfigurationError as error:
        print(type(error).__name__, error)
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] fastq_encoding must be one of ['sanger', 'solexa', 'illumina-1.3+', 'illumina-1.5+', 'illumina-1.8+']; got 'phred64'.
ConfigurationError [CONFIGURATION_ERROR] n_values must be at most 100.
```

**Usage with the dependency installed:**

```python
from dnakit.fasta_fastq import seqkit_extended_stats

result = seqkit_extended_stats(["reads.fq"], "stats.tsv", n_values=[90], verbose=False)
print(result["records"][0]["num_seqs"])
```

## 9) `dnakit.fasta_fastq.seqkit_fish_search()`

- **Purpose:** Locate short query sequences inside longer target sequences by local alignment.
- **Method:** Run `seqkit fish --query-fastx=<query_path> --min-qual=<min_mapping_quality>` (plus `--aln-params`, `--all`, `--stranded` as requested) on the target. Because SeqKit 2.14.0 ignores `-o` for `fish`, the facade takes the table that starts with the `Ref` header from SeqKit's output and writes it to `output_path` itself, then parses the rows (`Ref`, `Query` and `Strand` stay text).
- **API:** `seqkit_fish_search(target_path: 'str', query_path: 'str', output_path: 'str', *, min_mapping_quality: 'float' = 5.0, search_all: 'bool' = False, stranded: 'bool' = False, alignment_params: 'str | None' = None, seqkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program `seqkit` (2.14.0 is the tested version) on PATH, or an explicit `seqkit_path`; it is not installed in the environment that produced this page, and no Python package is required. The example below therefore shows behaviour that does not need SeqKit: an invalid argument raises `ConfigurationError` before SeqKit is looked up.
- **Input:** `target_path` and `query_path` are existing FASTA/FASTQ files. `output_path` must not exist. `min_mapping_quality` is a number (default 5.0). `alignment_params` is None or `match,mismatch,gap_open,gap_extend` as four integers (SeqKit default 4,-4,-2,-1). `search_all` and `stranded` are bools. Per the source docstring, targets containing N make SeqKit 2.14.0 exit with code 2, which is reported as `BackendExecutionError`.
- **Output:** The common result dictionary plus `records` (one row per hit). `outputs` holds `output`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.fasta_fastq import seqkit_fish_search

with tempfile.TemporaryDirectory() as tmp:
    target = Path(tmp) / "target.fa"
    query = Path(tmp) / "query.fa"
    target.write_text(">t1\nACGTACGTACGTAAGGCCTT\n")
    query.write_text(">q1\nAAGGCC\n")
    try:
        seqkit_fish_search(str(target), str(query), str(Path(tmp) / "hits.tsv"), alignment_params="4,-4,-2", verbose=False)
    except ConfigurationError as error:
        print(type(error).__name__, error)
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] alignment_params must look like '4,-4,-2,-1'.
```

**Usage with the dependency installed:**

```python
from dnakit.fasta_fastq import seqkit_fish_search

result = seqkit_fish_search("targets.fa", "queries.fa", "hits.tsv", search_all=True, verbose=False)
print(len(result["records"]))
```

## 10) `dnakit.fasta_fastq.seqkit_fastq_records_for_fasta()`

- **Purpose:** Return the FASTQ records that have a FASTA record with the same ID and the same sequence.
- **Method:** Run `seqkit fa2fq --fasta-file=<fasta_path>` on the FASTQ file (plus `--only-positive-strand` when requested) and write to `output_path` with `-o`. Per the source docstring, SeqKit 2.14.0 compares the full sequence on both strands unless `only_positive_strand` is True; a FASTQ record that is absent from the FASTA or has a different sequence is not returned.
- **API:** `seqkit_fastq_records_for_fasta(fastq_path: 'str', fasta_path: 'str', output_path: 'str', *, only_positive_strand: 'bool' = False, seqkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program `seqkit` (2.14.0 is the tested version) on PATH, or an explicit `seqkit_path`; it is not installed in the environment that produced this page, and no Python package is required. The example below therefore shows behaviour that does not need SeqKit: a `seqkit_path` that does not exist raises `BackendUnavailableError` (code `SEQKIT_UNAVAILABLE`), and no output file is written.
- **Input:** `fastq_path` and `fasta_path` are existing files. `output_path` must not exist. `only_positive_strand` is a bool.
- **Output:** The common result dictionary plus `record_count`. `outputs` holds `output`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.fasta_fastq import seqkit_fastq_records_for_fasta

with tempfile.TemporaryDirectory() as tmp:
    fastq = Path(tmp) / "reads.fq"
    fasta = Path(tmp) / "reads.fa"
    fastq.write_text("@r1\nACGT\n+\nIIII\n")
    fasta.write_text(">r1\nACGT\n")
    out = Path(tmp) / "matched.fq"
    try:
        seqkit_fastq_records_for_fasta(str(fastq), str(fasta), str(out), seqkit_path="/nonexistent/seqkit", verbose=False)
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
    print("output file written:", out.exists())
```

- **Example output:**

```text
BackendUnavailableError SEQKIT_UNAVAILABLE
output file written: False
```

**Usage with the dependency installed:**

```python
from dnakit.fasta_fastq import seqkit_fastq_records_for_fasta

result = seqkit_fastq_records_for_fasta("reads.fq", "reads.fa", "matched.fq", verbose=False)
print(result["record_count"])
```

## 11) `dnakit.fasta_fastq.seqkit_flanking_sequences()`

- **Purpose:** Extract BED or GTF features together with their upstream and downstream flanks.
- **Method:** Run `seqkit subseq` with `--bed=<feature_path>` or `--gtf=<feature_path>`, plus `--feature=<name>` for each GTF feature type, `--up-stream`, `--down-stream` and `--only-flank` as requested, and write to `output_path` with `-o`. SeqKit may leave a `<sequence>.seqkit.fai` index next to the FASTA; if it exists its path is returned as `index_file`.
- **API:** `seqkit_flanking_sequences(sequence_path: 'str', feature_path: 'str', output_path: 'str', *, feature_format: 'str' = 'bed', up_stream: 'int' = 0, down_stream: 'int' = 0, only_flank: 'bool' = False, feature_types: 'Sequence[str] | None' = None, seqkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program `seqkit` (2.14.0 is the tested version) on PATH, or an explicit `seqkit_path`; it is not installed in the environment that produced this page, and no Python package is required. The example below therefore shows behaviour that does not need SeqKit: an invalid argument raises `ConfigurationError` before SeqKit is looked up.
- **Input:** `sequence_path` and `feature_path` are existing files. `output_path` must not exist. `feature_format` is `bed` or `gtf`. `up_stream` and `down_stream` are integers of at least 0. With `only_flank=True` exactly one of them must be above 0. `feature_types` (a list of names) is allowed only for `gtf`.
- **Output:** The common result dictionary plus `record_count` and `index_file` (a path string or None). `outputs` holds `output`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.fasta_fastq import seqkit_flanking_sequences

with tempfile.TemporaryDirectory() as tmp:
    genome = Path(tmp) / "genome.fa"
    bed = Path(tmp) / "genes.bed"
    genome.write_text(">chr1\nACGTACGTACGTACGTACGT\n")
    bed.write_text("chr1\t5\t10\tgeneA\n")
    try:
        seqkit_flanking_sequences(str(genome), str(bed), str(Path(tmp) / "out.fa"), up_stream=2, down_stream=2, only_flank=True, verbose=False)
    except ConfigurationError as error:
        print(type(error).__name__, error)
    try:
        seqkit_flanking_sequences(str(genome), str(bed), str(Path(tmp) / "out.fa"), feature_types=["gene"], verbose=False)
    except ConfigurationError as error:
        print(type(error).__name__, error)
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] only_flank needs exactly one of up_stream or down_stream above 0.
ConfigurationError [CONFIGURATION_ERROR] feature_types only applies to feature_format='gtf'.
```

**Usage with the dependency installed:**

```python
from dnakit.fasta_fastq import seqkit_flanking_sequences

result = seqkit_flanking_sequences("genome.fa", "genes.bed", "genes_flanked.fa", up_stream=100, down_stream=100, verbose=False)
print(result["record_count"])
```

## 12) `dnakit.fasta_fastq.seqkit_grep_motif()`

- **Purpose:** Keep the records whose sequence contains any of the given motifs.
- **Method:** Run `seqkit grep --by-seq` with one `--pattern=` per motif, plus `--degenerate`, `--max-mismatch=<n>`, `--only-positive-strand`, `--ignore-case` and `--invert-match` as requested, and write to `output_path` with `-o`. Matching is literal by default and covers both strands unless `only_positive_strand` is True.
- **API:** `seqkit_grep_motif(input_path: 'str', output_path: 'str', *, motifs: 'Sequence[str]', degenerate: 'bool' = False, max_mismatch: 'int' = 0, only_positive_strand: 'bool' = False, ignore_case: 'bool' = False, invert_match: 'bool' = False, seqkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program `seqkit` (2.14.0 is the tested version) on PATH, or an explicit `seqkit_path`; it is not installed in the environment that produced this page, and no Python package is required. The example below therefore shows behaviour that does not need SeqKit: an invalid argument raises `ConfigurationError` before SeqKit is looked up.
- **Input:** `input_path` is an existing FASTA/FASTQ file. `output_path` must not exist. `motifs` is a non-empty list of non-empty strings without commas. `max_mismatch` is an integer of at least 0. The other options are bools.
- **Output:** The common result dictionary plus `record_count`. `outputs` holds `output`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.fasta_fastq import seqkit_grep_motif

with tempfile.TemporaryDirectory() as tmp:
    fasta = Path(tmp) / "reads.fa"
    fasta.write_text(">s1\nAAGGTCTCAA\n>s2\nCCCCCCCC\n")
    try:
        seqkit_grep_motif(str(fasta), str(Path(tmp) / "hits.fa"), motifs=["GGTCTC,GAGACC"], verbose=False)
    except ConfigurationError as error:
        print(type(error).__name__, error)
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] each pattern must be a non-empty string without commas (SeqKit splits pattern lists on commas).
```

**Usage with the dependency installed:**

```python
from dnakit.fasta_fastq import seqkit_grep_motif

result = seqkit_grep_motif("reads.fa", "with_site.fa", motifs=["GGTCTC"], verbose=False)
print(result["record_count"])
```

## 13) `dnakit.fasta_fastq.seqkit_head_genome()`

- **Purpose:** Return the records of the first genome in a file, where a genome is a run of records that share leading name words.
- **Method:** Run `seqkit head-genome --mini-common-words=<minimal_common_words>` and write to `output_path` with `-o`; the record count of the output is counted afterwards.
- **API:** `seqkit_head_genome(input_path: 'str', output_path: 'str', *, minimal_common_words: 'int' = 1, seqkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program `seqkit` (2.14.0 is the tested version) on PATH, or an explicit `seqkit_path`; it is not installed in the environment that produced this page, and no Python package is required. The example below therefore shows behaviour that does not need SeqKit: an invalid argument raises `ConfigurationError` before SeqKit is looked up.
- **Input:** `input_path` is an existing FASTA/FASTQ file. `output_path` must not exist. `minimal_common_words` is an integer of at least 1 (default 1).
- **Output:** The common result dictionary plus `record_count`. `outputs` holds `output`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.fasta_fastq import seqkit_head_genome

with tempfile.TemporaryDirectory() as tmp:
    fasta = Path(tmp) / "genomes.fa"
    fasta.write_text(">Escherichia coli contig1\nACGT\n>Escherichia coli contig2\nGGCC\n")
    try:
        seqkit_head_genome(str(fasta), str(Path(tmp) / "first.fa"), minimal_common_words=0, verbose=False)
    except ConfigurationError as error:
        print(type(error).__name__, error)
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] minimal_common_words must be at least 1.
```

**Usage with the dependency installed:**

```python
from dnakit.fasta_fastq import seqkit_head_genome

result = seqkit_head_genome("genomes.fa", "first_genome.fa", minimal_common_words=2, verbose=False)
print(result["record_count"])
```

## 14) `dnakit.fasta_fastq.seqkit_mean_read_quality()`

- **Purpose:** Report the average quality of every FASTQ read.
- **Method:** Run `seqkit fx2tab --name --avg-qual --qual-ascii-base=<quality_ascii_base>`, write the two-column table to `output_path` with `-o`, and parse it into rows. The values are SeqKit's own averages.
- **API:** `seqkit_mean_read_quality(fastq_path: 'str', output_path: 'str', *, quality_ascii_base: 'int' = 33, seqkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program `seqkit` (2.14.0 is the tested version) on PATH, or an explicit `seqkit_path`; it is not installed in the environment that produced this page, and no Python package is required. The example below therefore shows behaviour that does not need SeqKit: an invalid argument raises `ConfigurationError` before SeqKit is looked up.
- **Input:** `fastq_path` is an existing FASTQ file. `output_path` must not exist. `quality_ascii_base` is an integer from 0 to 127 (33 for Phred+33, the default).
- **Output:** The common result dictionary plus `records`, one `{"id", "avg_qual"}` row per read. `outputs` holds `output`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.fasta_fastq import seqkit_mean_read_quality

with tempfile.TemporaryDirectory() as tmp:
    fastq = Path(tmp) / "reads.fq"
    fastq.write_text("@r1\nACGT\n+\nIIII\n")
    try:
        seqkit_mean_read_quality(str(fastq), str(Path(tmp) / "quality.tsv"), quality_ascii_base=200, verbose=False)
    except ConfigurationError as error:
        print(type(error).__name__, error)
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] quality_ascii_base must be at most 127.
```

**Usage with the dependency installed:**

```python
from dnakit.fasta_fastq import seqkit_mean_read_quality

result = seqkit_mean_read_quality("reads.fq", "quality.tsv", verbose=False)
print(result["records"])
```

## 15) `dnakit.fasta_fastq.seqkit_merge_slides_to_bed()`

- **Purpose:** Merge overlapping sliding-window IDs into BED intervals.
- **Method:** Run `seqkit merge-slides --min-overlap=<min_overlap> --max-gap=<max_gap>` on a text file of window IDs (for example `chr1_sliding:1-1000`, one per line) and write the BED output to `output_path` with `-o`; then parse it into rows.
- **API:** `seqkit_merge_slides_to_bed(window_names_path: 'str', output_path: 'str', *, min_overlap: 'int' = 1, max_gap: 'int' = 0, seqkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program `seqkit` (2.14.0 is the tested version) on PATH, or an explicit `seqkit_path`; it is not installed in the environment that produced this page, and no Python package is required. The example below therefore shows behaviour that does not need SeqKit: an invalid argument raises `ConfigurationError` before SeqKit is looked up.
- **Input:** `window_names_path` is an existing text file. `output_path` must not exist. `min_overlap` is an integer of at least 1 (default 1). `max_gap` is an integer of at least 0 (default 0, meaning no limit on the distance between window starts).
- **Output:** The common result dictionary plus `records`, one `{"chrom", "start", "end"}` row per merged interval (0-based start). `outputs` holds `output`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.fasta_fastq import seqkit_merge_slides_to_bed

with tempfile.TemporaryDirectory() as tmp:
    names = Path(tmp) / "windows.txt"
    names.write_text("chr1_sliding:1-10\nchr1_sliding:6-15\n")
    try:
        seqkit_merge_slides_to_bed(str(names), str(Path(tmp) / "merged.bed"), min_overlap=0, verbose=False)
    except ConfigurationError as error:
        print(type(error).__name__, error)
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] min_overlap must be at least 1.
```

**Usage with the dependency installed:**

```python
from dnakit.fasta_fastq import seqkit_merge_slides_to_bed

result = seqkit_merge_slides_to_bed("windows.txt", "merged.bed", verbose=False)
print(result["records"])
```

## 16) `dnakit.fasta_fastq.seqkit_range_records()`

- **Purpose:** Select the records from one position to another, counting from the start or the end of the file.
- **Method:** Run `seqkit range --range=<start>:<end>` and write to `output_path` with `-o`. Positions are 1-based and inclusive; negative numbers count from the end, so (1, 12) is the first 12 records and (-3, -1) the last three.
- **API:** `seqkit_range_records(input_path: 'str', output_path: 'str', *, start: 'int', end: 'int', seqkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program `seqkit` (2.14.0 is the tested version) on PATH, or an explicit `seqkit_path`; it is not installed in the environment that produced this page, and no Python package is required. The example below therefore shows behaviour that does not need SeqKit: an invalid argument raises `ConfigurationError` before SeqKit is looked up.
- **Input:** `input_path` is an existing FASTA/FASTQ file. `output_path` must not exist. `start` and `end` are required integers (bools are rejected).
- **Output:** The common result dictionary plus `record_count`. `outputs` holds `output`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.fasta_fastq import seqkit_range_records

with tempfile.TemporaryDirectory() as tmp:
    fasta = Path(tmp) / "reads.fa"
    fasta.write_text(">s1\nAAAA\n>s2\nCCCC\n>s3\nGGGG\n")
    try:
        seqkit_range_records(str(fasta), str(Path(tmp) / "range.fa"), start="1", end=2, verbose=False)
    except ConfigurationError as error:
        print(type(error).__name__, error)
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] start must be an integer.
```

**Usage with the dependency installed:**

```python
from dnakit.fasta_fastq import seqkit_range_records

result = seqkit_range_records("reads.fa", "last_two.fa", start=-2, end=-1, verbose=False)
print(result["record_count"])
```

## 17) `dnakit.fasta_fastq.seqkit_read_feature_histogram()`

- **Purpose:** Build a histogram of one read feature (length, mean quality, GC or GC skew) for a FASTA/FASTQ file.
- **Method:** Run `seqkit watch --fields=<field> --bins=<bins> --qual-ascii-base=<quality_ascii_base> --dump --quiet-mode` (plus `--log`, `--validate-seq` and `--img=<image_path>` as requested) and parse the `Bin` table from SeqKit's output. The report is produced once at the end of the file; online monitoring options are not offered.
- **API:** `seqkit_read_feature_histogram(input_path: 'str', field: 'str', *, bins: 'int' = -1, log_transform: 'bool' = False, quality_ascii_base: 'int' = 33, validate_seq: 'bool' = False, image_path: 'str | None' = None, seqkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program `seqkit` (2.14.0 is the tested version) on PATH, or an explicit `seqkit_path`; it is not installed in the environment that produced this page, and no Python package is required. The example below therefore shows behaviour that does not need SeqKit: an invalid argument raises `ConfigurationError` before SeqKit is looked up.
- **Input:** `input_path` is an existing FASTA/FASTQ file. `field` is one of `ReadLen`, `MeanQual`, `GC`, `GCSkew`. `bins` is an integer of at least -1. `quality_ascii_base` is an integer from 0 to 127. `log_transform` and `validate_seq` are bools. `image_path` is optional and must be a new file in an existing folder.
- **Output:** The common result dictionary plus `records` (the histogram bins) and `field`. `outputs` holds `image` when an image was requested.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.fasta_fastq import seqkit_read_feature_histogram

with tempfile.TemporaryDirectory() as tmp:
    fastq = Path(tmp) / "reads.fq"
    fastq.write_text("@r1\nACGT\n+\nIIII\n")
    try:
        seqkit_read_feature_histogram(str(fastq), "Length", verbose=False)
    except ConfigurationError as error:
        print(type(error).__name__, error)
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] field must be one of ['ReadLen', 'MeanQual', 'GC', 'GCSkew']; got 'Length'.
```

**Usage with the dependency installed:**

```python
from dnakit.fasta_fastq import seqkit_read_feature_histogram

result = seqkit_read_feature_histogram("reads.fq", "MeanQual", bins=5, verbose=False)
print(result["field"], len(result["records"]))
```

## 18) `dnakit.fasta_fastq.seqkit_rename_duplicate_ids()`

- **Purpose:** Make repeated IDs unique by adding a counter suffix.
- **Method:** Run `seqkit rename --start-num=<start_number> --separator=<separator>` (plus `--by-name` and `--rename-1st-rec` as requested) and write to `output_path` with `-o`. The first occurrence keeps its ID unless `rename_first` is true, so `id, id, id` becomes `id, id_2, id_3` with the defaults.
- **API:** `seqkit_rename_duplicate_ids(input_path: 'str', output_path: 'str', *, by_name: 'bool' = False, start_number: 'int' = 2, separator: 'str' = '_', rename_first: 'bool' = False, seqkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program `seqkit` (2.14.0 is the tested version) on PATH, or an explicit `seqkit_path`; it is not installed in the environment that produced this page, and no Python package is required. The example below therefore shows behaviour that does not need SeqKit: an invalid argument raises `ConfigurationError` before SeqKit is looked up.
- **Input:** `input_path` is an existing FASTA/FASTQ file. `output_path` must not exist. `by_name` and `rename_first` are bools. `start_number` is an integer of at least 1 (default 2). `separator` is a non-empty string (default `_`).
- **Output:** The common result dictionary plus `record_count`. `outputs` holds `output`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.fasta_fastq import seqkit_rename_duplicate_ids

with tempfile.TemporaryDirectory() as tmp:
    fasta = Path(tmp) / "dups.fa"
    fasta.write_text(">s1\nAAAA\n>s1\nCCCC\n")
    try:
        seqkit_rename_duplicate_ids(str(fasta), str(Path(tmp) / "unique.fa"), start_number=0, verbose=False)
    except ConfigurationError as error:
        print(type(error).__name__, error)
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] start_number must be at least 1.
```

**Usage with the dependency installed:**

```python
from dnakit.fasta_fastq import seqkit_rename_duplicate_ids

result = seqkit_rename_duplicate_ids("dups.fa", "unique.fa", verbose=False)
print(result["record_count"])
```

## 19) `dnakit.fasta_fastq.seqkit_replace_sequence_regex()`

- **Purpose:** Apply a regular-expression substitution to the sequences of a file.
- **Method:** Run `seqkit replace --by-seq --pattern=<pattern> --replacement=<replacement>` (plus `--ignore-case` when requested) and write to `output_path` with `-o`. The replacement may use capture groups such as `$1`. For FASTQ the sequence length must not change, because the quality string keeps its length.
- **API:** `seqkit_replace_sequence_regex(input_path: 'str', output_path: 'str', *, pattern: 'str', replacement: 'str', ignore_case: 'bool' = False, seqkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program `seqkit` (2.14.0 is the tested version) on PATH, or an explicit `seqkit_path`; it is not installed in the environment that produced this page, and no Python package is required. The example below therefore shows behaviour that does not need SeqKit: an invalid argument raises `ConfigurationError` before SeqKit is looked up.
- **Input:** `input_path` is an existing FASTA/FASTQ file. `output_path` must not exist. `pattern` is a non-empty string. `replacement` is a string (it may be empty). `ignore_case` is a bool.
- **Output:** The common result dictionary plus `record_count`. `outputs` holds `output`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.fasta_fastq import seqkit_replace_sequence_regex

with tempfile.TemporaryDirectory() as tmp:
    fasta = Path(tmp) / "reads.fa"
    fasta.write_text(">s1\nAANNAA\n")
    try:
        seqkit_replace_sequence_regex(str(fasta), str(Path(tmp) / "fixed.fa"), pattern="N", replacement=None, verbose=False)
    except ConfigurationError as error:
        print(type(error).__name__, error)
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] replacement must be a string.
```

**Usage with the dependency installed:**

```python
from dnakit.fasta_fastq import seqkit_replace_sequence_regex

result = seqkit_replace_sequence_regex("reads.fa", "no_n.fa", pattern="N", replacement="A", verbose=False)
print(result["record_count"])
```

## 20) `dnakit.fasta_fastq.seqkit_sequence_digest()`

- **Purpose:** Compute an order-aware message digest of the sequences of each file.
- **Method:** Run `seqkit sum` (plus `--include-id`, `--remove-gaps` and `--rna2dna` as requested) and write the output to `output_path` with `-o`; each line is parsed into a digest and a file name. The digest has the form `seqkit.v0.1_...` and is not a SHA-256 of the file.
- **API:** `seqkit_sequence_digest(input_paths: 'Sequence[str]', output_path: 'str', *, include_ids: 'bool' = False, remove_gaps: 'bool' = False, rna_to_dna: 'bool' = False, seqkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program `seqkit` (2.14.0 is the tested version) on PATH, or an explicit `seqkit_path`; it is not installed in the environment that produced this page, and no Python package is required. The example below therefore shows behaviour that does not need SeqKit: an existing output file is refused and left untouched.
- **Input:** `input_paths` is a non-empty list of existing files. `output_path` must not exist (an existing file is never overwritten). `include_ids`, `remove_gaps` and `rna_to_dna` are bools.
- **Output:** The common result dictionary plus `records`, one `{"digest", "file"}` row per input file. `outputs` holds `output`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.fasta_fastq import seqkit_sequence_digest

with tempfile.TemporaryDirectory() as tmp:
    fasta = Path(tmp) / "reads.fa"
    fasta.write_text(">s1\nAAAA\n")
    existing = Path(tmp) / "digest.tsv"
    existing.write_text("old\n")
    try:
        seqkit_sequence_digest([str(fasta)], str(existing), verbose=False)
    except ConfigurationError as error:
        print(type(error).__name__, error.code)
    print("existing file kept:", existing.read_text().strip())
```

- **Example output:**

```text
ConfigurationError CONFIGURATION_ERROR
existing file kept: old
```

**Usage with the dependency installed:**

```python
from dnakit.fasta_fastq import seqkit_sequence_digest

result = seqkit_sequence_digest(["reads.fa"], "digest.tsv", verbose=False)
print(result["records"][0]["digest"].startswith("seqkit.v0.1_"))
```

## 21) `dnakit.fasta_fastq.seqkit_split_by_region()`

- **Purpose:** Split the records into files according to the subsequence found in a region, for example a barcode.
- **Method:** Run `seqkit split --by-region=<region> --out-dir=<output_dir>`; SeqKit creates the directory and one file per distinct subsequence. The facade then lists the files in the directory and counts the records of each.
- **API:** `seqkit_split_by_region(input_path: 'str', output_dir: 'str', *, region: 'str', seqkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program `seqkit` (2.14.0 is the tested version) on PATH, or an explicit `seqkit_path`; it is not installed in the environment that produced this page, and no Python package is required. The example below therefore shows behaviour that does not need SeqKit: an invalid argument raises `ConfigurationError` before SeqKit is looked up.
- **Input:** `input_path` is an existing FASTA/FASTQ file. `output_dir` must not exist and its parent must exist. `region` is a string of the form `start:end` with optional minus signs, such as `1:5` or `-12:-1`.
- **Output:** The common result dictionary plus `records`, one `{"file", "record_count"}` row per file written. `outputs` holds `directory`. If SeqKit does not create the directory, `BackendExecutionError` with code `SEQKIT_OUTPUT_MISSING` is raised.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.fasta_fastq import seqkit_split_by_region

with tempfile.TemporaryDirectory() as tmp:
    fastq = Path(tmp) / "reads.fq"
    fastq.write_text("@r1\nACGTAC\n+\nIIIIII\n")
    try:
        seqkit_split_by_region(str(fastq), str(Path(tmp) / "split"), region="first5", verbose=False)
    except ConfigurationError as error:
        print(type(error).__name__, error)
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] region must look like '1:5' or '-12:-1'.
```

**Usage with the dependency installed:**

```python
from dnakit.fasta_fastq import seqkit_split_by_region

result = seqkit_split_by_region("reads.fq", "by_barcode", region="1:4", verbose=False)
print(result["records"])
```

## 22) `dnakit.fasta_fastq.create_screed_database()`

- **Purpose:** Index a FASTA or FASTQ file into a screed database at a path of your choice.
- **Method:** Open the file with `screed.openscreed.Open(..., parse_description=True)` to detect the format from the first record, choose the FASTA or FASTQ field types, and build the SQLite database with `screed.createscreed.create_db`. If the output name does not end with screed's `_screed` suffix, the suffix is added. The record count is then read back from the database with SQL.
- **API:** `create_screed_database(sequence_file: 'str', output_database: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package `screed` (extra `dnakit[bio-screed]`); it is installed in the environment that produced this page and the example runs for real.
- **Input:** `sequence_file` is an existing FASTA or FASTQ file. `output_database` must not exist (the name with the `_screed` suffix must not exist either) and its folder must exist.
- **Output:** A dictionary with `status`, `backend` (`screed`), `input_format` (`FASTA` or `FASTQ`), `database` (the absolute path, with the suffix) and `record_count`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.fasta_fastq import create_screed_database

with tempfile.TemporaryDirectory() as tmp:
    fasta_result = create_screed_database("reads.fa", str(Path(tmp) / "fa_db"), verbose=False)
    print(fasta_result["input_format"], fasta_result["record_count"], Path(fasta_result["database"]).name)
    fastq_result = create_screed_database("reads.fq", str(Path(tmp) / "fq_db_screed"), verbose=False)
    print(fastq_result["input_format"], fastq_result["record_count"], Path(fastq_result["database"]).name)
```

- **Example output:**

```text
FASTA 2 fa_db_screed
FASTQ 2 fq_db_screed
```

## 23) `dnakit.fasta_fastq.dump_screed_fasta()`

- **Purpose:** Write every record of a screed database as FASTA text.
- **Method:** Check that the database has records and a `name` and `sequence` field, then call `screed.conversion.ToFasta(database, output)`. The header line is the record name followed by its description; for a database built from FASTQ the quality values are dropped.
- **API:** `dump_screed_fasta(database: 'str', output_fasta: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package `screed` (extra `dnakit[bio-screed]`); it is installed in the environment that produced this page and the example runs for real.
- **Input:** `database` is an existing screed database file. `output_fasta` must not exist and its folder must exist. An empty database raises `ConfigurationError`.
- **Output:** A dictionary with `status`, `backend`, `database`, `format` (`FASTA`), `output` (the absolute path) and `record_count`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.fasta_fastq import create_screed_database, dump_screed_fasta

with tempfile.TemporaryDirectory() as tmp:
    database = create_screed_database("reads.fq", str(Path(tmp) / "reads_db"), verbose=False)["database"]
    out = Path(tmp) / "reads_dump.fa"
    result = dump_screed_fasta(database, str(out), verbose=False)
    print(result["format"], result["record_count"])
    print(out.read_text(), end="")
```

- **Example output:**

```text
FASTA 2
>r1 first read
ACGTACGT
>r2 second read
GGCCTTAA
```

## 24) `dnakit.fasta_fastq.dump_screed_fastq()`

- **Purpose:** Write every record of a screed database as FASTQ text.
- **Method:** Check that the database has records and a `name`, `sequence` and `quality` field, then call `screed.conversion.ToFastq(database, output)`. A database built from FASTA has no quality field and is rejected with `ConfigurationError`.
- **API:** `dump_screed_fastq(database: 'str', output_fastq: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the optional package `screed` (extra `dnakit[bio-screed]`); it is installed in the environment that produced this page and the example runs for real.
- **Input:** `database` is an existing screed database built from FASTQ. `output_fastq` must not exist and its folder must exist.
- **Output:** A dictionary with `status`, `backend`, `database`, `format` (`FASTQ`), `output` (the absolute path) and `record_count`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.fasta_fastq import create_screed_database, dump_screed_fastq

with tempfile.TemporaryDirectory() as tmp:
    fq_db = create_screed_database("reads.fq", str(Path(tmp) / "fq_db"), verbose=False)["database"]
    out = Path(tmp) / "reads_dump.fq"
    result = dump_screed_fastq(fq_db, str(out), verbose=False)
    print(result["format"], result["record_count"])
    print(out.read_text(), end="")
    fa_db = create_screed_database("reads.fa", str(Path(tmp) / "fa_db"), verbose=False)["database"]
    try:
        dump_screed_fastq(fa_db, str(Path(tmp) / "none.fq"), verbose=False)
    except ConfigurationError as error:
        print(type(error).__name__, "no quality field" in str(error))
```

- **Example output:**

```text
FASTQ 2
@r1 first read
ACGTACGT
+
IIIIHHHH
@r2 second read
GGCCTTAA
+
########
ConfigurationError True
```
