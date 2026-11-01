# seqtk and NanoPlot: sequence-file and long-read command-line operations

This group of functions wraps the external command-line programs [seqtk](https://github.com/lh3/seqtk) (16 functions, module `seqtk_p12_cli_functions` in `dnakit.fasta_fastq`), [NanoPlot](https://github.com/wdecoster/NanoPlot) (1 functions, module `nanoplot_p12_cli_functions` in `dnakit.sequencing_qc`). Every function runs one fixed program or sub-command with an explicit, validated argument list (no shell, timeout and output-size limits), writes results to explicit new paths and never overwrites an existing output. Tested versions: seqtk 1.5-r133, NanoPlot 1.48.0. Import them with `from <area> import ...` once the area package exports them, or directly from the module.

The programs are not bundled with DNAKit: install them (seqtk: `conda install -c bioconda seqtk`; NanoPlot: `conda install -c bioconda nanoplot` or `pip install NanoPlot`), then put them on PATH or pass the path option of the function. If a program cannot be found, `BackendUnavailableError` is raised with an install hint; if it exits with an error, `BackendExecutionError` is raised with the tail of its log. Results are JSON-serialisable dictionaries with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and `operation`; each section lists only the extra keys. The programs are not installed in the environment used to produce this page, so every checked example shows behaviour that is real on any machine (a missing program); each is followed by the call to run with the program installed.

## 1) `dnakit.fasta_fastq.seqtk_cutn()`

- **Purpose:** Cut sequences at long N tracts with ``seqtk cutN``, or list the tracts only.
- **Method:** ``min_n_tract`` is ``-n`` (default 1000) and ``non_n_penalty`` the penalty for a non-N base inside a tract (``-p``, default 10). With ``gaps_only=True`` (``-g``) no sequence is written; ``output_path`` then holds the gaps as ``seqid``, ``start``, ``end`` rows, which are also returned in ``records``. Otherwise the pieces are written as FASTA records named ``seqid:start-end`` and ``record_count`` is their number.
- **API:** `seqtk_cutn(input_path: 'str', output_path: 'str', *, min_n_tract: 'int' = 1000, non_n_penalty: 'int' = 10, gaps_only: 'bool' = False, max_records: 'int' = 10000, seqtk_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program(s) of seqtk (1.5-r133 is the tested version) on PATH or through `seqtk_path`; it is not installed in the environment that produced this page. The example below shows behaviour that does not need it: the error raised for a missing program.
- **Input:** required: `input_path`, `output_path`; optional keyword-only: `min_n_tract`, `non_n_penalty`, `gaps_only`, `max_records`, `seqtk_path` (meanings and defaults are in the Method text and the signature). Input files must exist; output paths must not exist yet.
- **Output:** The common result dictionary. No extra keys beyond `outputs` and `log_tail`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_cutn

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'in.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        seqtk_cutn(
            input_path=str(root / 'in.fa'),
            output_path=str(root / 'result.txt'),
            seqtk_path='/nonexistent/seqtk',
            verbose=False,
        )
    except DNAKitError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError BACKEND_UNAVAILABLE
```

**Usage with the program installed:**

```python
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_cutn

result = seqtk_cutn(
    input_path='in.fa',
    output_path='result.txt',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 2) `dnakit.fasta_fastq.seqtk_dropse()`

- **Purpose:** Drop the unpaired reads from an interleaved paired-end FASTA or FASTQ file with ``seqtk dropse``.
- **Method:** Two consecutive records are kept only when their names are the same except for a trailing ``/1`` ``/2``. ``output_path`` receives the paired reads.
- **API:** `seqtk_dropse(input_path: 'str', output_path: 'str', *, seqtk_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program(s) of seqtk (1.5-r133 is the tested version) on PATH or through `seqtk_path`; it is not installed in the environment that produced this page. The example below shows behaviour that does not need it: the error raised for a missing program.
- **Input:** required: `input_path`, `output_path`; optional keyword-only: `seqtk_path` (meanings and defaults are in the Method text and the signature). Input files must exist; output paths must not exist yet.
- **Output:** The common result dictionary. Result keys: ``output_file`` and ``record_count``.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_dropse

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'in.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        seqtk_dropse(
            input_path=str(root / 'in.fa'),
            output_path=str(root / 'result.txt'),
            seqtk_path='/nonexistent/seqtk',
            verbose=False,
        )
    except DNAKitError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError BACKEND_UNAVAILABLE
```

**Usage with the program installed:**

```python
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_dropse

result = seqtk_dropse(
    input_path='in.fa',
    output_path='result.txt',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 3) `dnakit.fasta_fastq.seqtk_famask()`

- **Purpose:** Apply an X-coded mask FASTA to a source FASTA with ``seqtk famask``.
- **Method:** Positions that are lower case or ``X`` in ``mask_path`` become lower case in the copy of ``source_path`` written to ``output_path``.
- **API:** `seqtk_famask(source_path: 'str', mask_path: 'str', output_path: 'str', *, seqtk_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program(s) of seqtk (1.5-r133 is the tested version) on PATH or through `seqtk_path`; it is not installed in the environment that produced this page. The example below shows behaviour that does not need it: the error raised for a missing program.
- **Input:** required: `source_path`, `mask_path`, `output_path`; optional keyword-only: `seqtk_path` (meanings and defaults are in the Method text and the signature). Input files must exist; output paths must not exist yet.
- **Output:** The common result dictionary. Result keys: ``output_file`` and ``record_count``.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_famask

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'mask.fa').write_text('x\n')
    (root / 'src.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        seqtk_famask(
            source_path=str(root / 'src.fa'),
            mask_path=str(root / 'mask.fa'),
            output_path=str(root / 'result.txt'),
            seqtk_path='/nonexistent/seqtk',
            verbose=False,
        )
    except DNAKitError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError BACKEND_UNAVAILABLE
```

**Usage with the program installed:**

```python
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_famask

result = seqtk_famask(
    source_path='src.fa',
    mask_path='mask.fa',
    output_path='result.txt',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 4) `dnakit.fasta_fastq.seqtk_fqchk()`

- **Purpose:** Summarise base composition and base quality per read position with ``seqtk fqchk``.
- **Method:** ``quality_threshold`` is ``-q`` (default 20): the table then ends with the percentage of bases below and at or above that quality (``%low`` and ``%high``); 0 gives the percentage of every observed quality value instead. The table is saved in ``output_path``.
- **API:** `seqtk_fqchk(input_path: 'str', output_path: 'str', *, quality_threshold: 'int' = 20, max_records: 'int' = 10000, seqtk_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program(s) of seqtk (1.5-r133 is the tested version) on PATH or through `seqtk_path`; it is not installed in the environment that produced this page. The example below shows behaviour that does not need it: the error raised for a missing program.
- **Input:** required: `input_path`, `output_path`; optional keyword-only: `quality_threshold`, `max_records`, `seqtk_path` (meanings and defaults are in the Method text and the signature). Input files must exist; output paths must not exist yet.
- **Output:** The common result dictionary. Result keys: ``summary`` (``min_len``, ``max_len``, ``avg_len``, ``distinct_quality_values``), ``columns``, ``rows`` (the row ``ALL`` first, then one row per read position, at most ``max_records``) and ``row_count``.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_fqchk

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'in.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        seqtk_fqchk(
            input_path=str(root / 'in.fa'),
            output_path=str(root / 'result.txt'),
            seqtk_path='/nonexistent/seqtk',
            verbose=False,
        )
    except DNAKitError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError BACKEND_UNAVAILABLE
```

**Usage with the program installed:**

```python
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_fqchk

result = seqtk_fqchk(
    input_path='in.fa',
    output_path='result.txt',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 5) `dnakit.fasta_fastq.seqtk_gap()`

- **Purpose:** List the N-gap locations (runs of non-ACGT bases) of at least ``min_length`` with ``seqtk gap``.
- **Method:** ``min_length`` is ``-l`` (default 50).
- **API:** `seqtk_gap(input_path: 'str', output_path: 'str', *, min_length: 'int' = 50, max_records: 'int' = 10000, seqtk_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program(s) of seqtk (1.5-r133 is the tested version) on PATH or through `seqtk_path`; it is not installed in the environment that produced this page. The example below shows behaviour that does not need it: the error raised for a missing program.
- **Input:** required: `input_path`, `output_path`; optional keyword-only: `min_length`, `max_records`, `seqtk_path` (meanings and defaults are in the Method text and the signature). Input files must exist; output paths must not exist yet.
- **Output:** The common result dictionary. Result keys: ``records`` (``seqid``, 0-based ``start``, ``end``; at most ``max_records``) and ``record_count``; the same table is saved in ``output_path``.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_gap

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'in.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        seqtk_gap(
            input_path=str(root / 'in.fa'),
            output_path=str(root / 'result.txt'),
            seqtk_path='/nonexistent/seqtk',
            verbose=False,
        )
    except DNAKitError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError BACKEND_UNAVAILABLE
```

**Usage with the program installed:**

```python
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_gap

result = seqtk_gap(
    input_path='in.fa',
    output_path='result.txt',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 6) `dnakit.fasta_fastq.seqtk_gc()`

- **Purpose:** Identify GC-rich (default) or AT-rich regions of a FASTA with the X-dropoff algorithm of ``seqtk gc``.
- **Method:** ``at_rich=True`` looks for high-AT regions (``-w``) and then ``min_fraction`` is an AT fraction. ``min_fraction`` is ``-f`` (default 0.60), ``min_length`` the shortest region reported (``-l``, default 20) and ``x_dropoff`` the X-dropoff (``-x``, default 10.0).
- **API:** `seqtk_gc(input_path: 'str', output_path: 'str', *, at_rich: 'bool' = False, min_fraction: 'float' = 0.6, min_length: 'int' = 20, x_dropoff: 'float' = 10.0, max_records: 'int' = 10000, seqtk_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program(s) of seqtk (1.5-r133 is the tested version) on PATH or through `seqtk_path`; it is not installed in the environment that produced this page. The example below shows behaviour that does not need it: the error raised for a missing program.
- **Input:** required: `input_path`, `output_path`; optional keyword-only: `at_rich`, `min_fraction`, `min_length`, `x_dropoff`, `max_records`, `seqtk_path` (meanings and defaults are in the Method text and the signature). Input files must exist; output paths must not exist yet.
- **Output:** The common result dictionary. Result keys: ``records`` (``seqid``, ``start`` 0-based, ``end``, ``hit_bases``; at most ``max_records``) and ``record_count``.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_gc

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'in.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        seqtk_gc(
            input_path=str(root / 'in.fa'),
            output_path=str(root / 'result.txt'),
            seqtk_path='/nonexistent/seqtk',
            verbose=False,
        )
    except DNAKitError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError BACKEND_UNAVAILABLE
```

**Usage with the program installed:**

```python
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_gc

result = seqtk_gc(
    input_path='in.fa',
    output_path='result.txt',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 7) `dnakit.fasta_fastq.seqtk_hety()`

- **Purpose:** Measure the regional heterozygosity of a FASTA with IUPAC ambiguity codes with ``seqtk hety``.
- **Method:** ``window_size`` is ``-w`` (default 50000), ``start_positions`` the number of window start positions per window (``-t``, default 5; the window slides by ``window_size / start_positions``) and ``mask_lowercase`` treats lower-case bases as masked (``-m``).
- **API:** `seqtk_hety(input_path: 'str', output_path: 'str', *, window_size: 'int' = 50000, start_positions: 'int' = 5, mask_lowercase: 'bool' = False, max_records: 'int' = 10000, seqtk_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program(s) of seqtk (1.5-r133 is the tested version) on PATH or through `seqtk_path`; it is not installed in the environment that produced this page. The example below shows behaviour that does not need it: the error raised for a missing program.
- **Input:** required: `input_path`, `output_path`; optional keyword-only: `window_size`, `start_positions`, `mask_lowercase`, `max_records`, `seqtk_path` (meanings and defaults are in the Method text and the signature). Input files must exist; output paths must not exist yet.
- **Output:** The common result dictionary. Result keys: ``records`` (``seqid``, ``start``, ``end``, ``hets_per_window``, ``called_sites``, ``het_sites``; at most ``max_records``) and ``record_count``.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_hety

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'in.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        seqtk_hety(
            input_path=str(root / 'in.fa'),
            output_path=str(root / 'result.txt'),
            seqtk_path='/nonexistent/seqtk',
            verbose=False,
        )
    except DNAKitError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError BACKEND_UNAVAILABLE
```

**Usage with the program installed:**

```python
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_hety

result = seqtk_hety(
    input_path='in.fa',
    output_path='result.txt',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 8) `dnakit.fasta_fastq.seqtk_hpc()`

- **Purpose:** Homopolymer-compress every sequence of a FASTA or FASTQ file with ``seqtk hpc``.
- **Method:** Each run of identical bases is replaced by a single base.
- **API:** `seqtk_hpc(input_path: 'str', output_path: 'str', *, seqtk_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program(s) of seqtk (1.5-r133 is the tested version) on PATH or through `seqtk_path`; it is not installed in the environment that produced this page. The example below shows behaviour that does not need it: the error raised for a missing program.
- **Input:** required: `input_path`, `output_path`; optional keyword-only: `seqtk_path` (meanings and defaults are in the Method text and the signature). Input files must exist; output paths must not exist yet.
- **Output:** The common result dictionary. Result keys: ``output_file`` and ``record_count``.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_hpc

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'in.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        seqtk_hpc(
            input_path=str(root / 'in.fa'),
            output_path=str(root / 'result.txt'),
            seqtk_path='/nonexistent/seqtk',
            verbose=False,
        )
    except DNAKitError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError BACKEND_UNAVAILABLE
```

**Usage with the program installed:**

```python
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_hpc

result = seqtk_hpc(
    input_path='in.fa',
    output_path='result.txt',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 9) `dnakit.fasta_fastq.seqtk_kfreq()`

- **Purpose:** Count the occurrences of a k-mer and of its single-mismatch neighbours in every sequence with ``seqtk kfreq``.
- **Method:** ``kmer`` is a DNA word of 1 to 12 bases (A, C, G, T). Both strands are scanned; the reported strand is the one with more neighbour hits.
- **API:** `seqtk_kfreq(kmer: 'str', input_path: 'str', output_path: 'str', *, max_records: 'int' = 10000, seqtk_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program(s) of seqtk (1.5-r133 is the tested version) on PATH or through `seqtk_path`; it is not installed in the environment that produced this page. The example below shows behaviour that does not need it: the error raised for a missing program.
- **Input:** required: `kmer`, `input_path`, `output_path`; optional keyword-only: `max_records`, `seqtk_path` (meanings and defaults are in the Method text and the signature). Input files must exist; output paths must not exist yet.
- **Output:** The common result dictionary. Result keys: ``records`` (``seqid``, ``length``, ``strand``, ``neighbour_count`` and ``exact_count``; at most ``max_records``) and ``record_count``; the same table is saved in ``output_path``.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_kfreq

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'in.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        seqtk_kfreq(
            kmer='ACGT',
            input_path=str(root / 'in.fa'),
            output_path=str(root / 'result.txt'),
            seqtk_path='/nonexistent/seqtk',
            verbose=False,
        )
    except DNAKitError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError BACKEND_UNAVAILABLE
```

**Usage with the program installed:**

```python
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_kfreq

result = seqtk_kfreq(
    kmer='ACGT',
    input_path='in.fa',
    output_path='result.txt',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 10) `dnakit.fasta_fastq.seqtk_listhet()`

- **Purpose:** List the position of every heterozygous (two-allele IUPAC) base of a FASTA with ``seqtk listhet``.
- **Method:** List the position of every heterozygous (two-allele IUPAC) base of a FASTA with ``seqtk listhet``.
- **API:** `seqtk_listhet(input_path: 'str', output_path: 'str', *, max_records: 'int' = 10000, seqtk_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program(s) of seqtk (1.5-r133 is the tested version) on PATH or through `seqtk_path`; it is not installed in the environment that produced this page. The example below shows behaviour that does not need it: the error raised for a missing program.
- **Input:** required: `input_path`, `output_path`; optional keyword-only: `max_records`, `seqtk_path` (meanings and defaults are in the Method text and the signature). Input files must exist; output paths must not exist yet.
- **Output:** The common result dictionary. Result keys: ``records`` (``seqid``, 1-based ``position``, ``base``; at most ``max_records``) and ``record_count``; the same table is saved in ``output_path``.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_listhet

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'in.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        seqtk_listhet(
            input_path=str(root / 'in.fa'),
            output_path=str(root / 'result.txt'),
            seqtk_path='/nonexistent/seqtk',
            verbose=False,
        )
    except DNAKitError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError BACKEND_UNAVAILABLE
```

**Usage with the program installed:**

```python
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_listhet

result = seqtk_listhet(
    input_path='in.fa',
    output_path='result.txt',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 11) `dnakit.fasta_fastq.seqtk_mergefa()`

- **Purpose:** Merge two FASTA/FASTQ allele sets (for example the two haplotypes of one individual) with ``seqtk mergefa``.
- **Method:** Bases that differ become IUPAC heterozygous codes. ``quality_threshold`` is ``-q`` (default 0), ``intersection`` takes the intersection of the inputs (``-i``), ``lowercase_on_n`` converts to lower case where one input base is N (``-m``), ``random_allele`` picks a random allele from heterozygous input bases (``-r``) and ``suppress_hets`` removes heterozygous bases from the input (``-h``).
- **API:** `seqtk_mergefa(first_path: 'str', second_path: 'str', output_path: 'str', *, quality_threshold: 'int' = 0, intersection: 'bool' = False, lowercase_on_n: 'bool' = False, random_allele: 'bool' = False, suppress_hets: 'bool' = False, seqtk_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program(s) of seqtk (1.5-r133 is the tested version) on PATH or through `seqtk_path`; it is not installed in the environment that produced this page. The example below shows behaviour that does not need it: the error raised for a missing program.
- **Input:** required: `first_path`, `second_path`, `output_path`; optional keyword-only: `quality_threshold`, `intersection`, `lowercase_on_n`, `random_allele`, `suppress_hets`, `seqtk_path` (meanings and defaults are in the Method text and the signature). Input files must exist; output paths must not exist yet.
- **Output:** The common result dictionary. Result keys: ``output_file`` and ``record_count``.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_mergefa

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'a.fa').write_text('x\n')
    (root / 'b.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        seqtk_mergefa(
            first_path=str(root / 'a.fa'),
            second_path=str(root / 'b.fa'),
            output_path=str(root / 'result.txt'),
            seqtk_path='/nonexistent/seqtk',
            verbose=False,
        )
    except DNAKitError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError BACKEND_UNAVAILABLE
```

**Usage with the program installed:**

```python
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_mergefa

result = seqtk_mergefa(
    first_path='a.fa',
    second_path='b.fa',
    output_path='result.txt',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 12) `dnakit.fasta_fastq.seqtk_mutfa()`

- **Purpose:** Introduce point mutations into a FASTA at listed positions with ``seqtk mutfa``.
- **Method:** ``snp_path`` is a text file whose lines have at least four columns: sequence name, 1-based position, any value (ignored) and the base to change to. The mutated FASTA is written to ``output_path``.
- **API:** `seqtk_mutfa(input_path: 'str', snp_path: 'str', output_path: 'str', *, seqtk_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program(s) of seqtk (1.5-r133 is the tested version) on PATH or through `seqtk_path`; it is not installed in the environment that produced this page. The example below shows behaviour that does not need it: the error raised for a missing program.
- **Input:** required: `input_path`, `snp_path`, `output_path`; optional keyword-only: `seqtk_path` (meanings and defaults are in the Method text and the signature). Input files must exist; output paths must not exist yet.
- **Output:** The common result dictionary. Result keys: ``output_file`` and ``record_count``.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_mutfa

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'in.fa').write_text('x\n')
    (root / 'snps.txt').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        seqtk_mutfa(
            input_path=str(root / 'in.fa'),
            snp_path=str(root / 'snps.txt'),
            output_path=str(root / 'result.txt'),
            seqtk_path='/nonexistent/seqtk',
            verbose=False,
        )
    except DNAKitError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError BACKEND_UNAVAILABLE
```

**Usage with the program installed:**

```python
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_mutfa

result = seqtk_mutfa(
    input_path='in.fa',
    snp_path='snps.txt',
    output_path='result.txt',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 13) `dnakit.fasta_fastq.seqtk_randbase()`

- **Purpose:** Resolve every heterozygous IUPAC base of a FASTA to one of its two alleles at random with ``seqtk randbase``.
- **Method:** seqtk draws with the C ``drand48`` generator without seeding, so a given file always gives the same output. Sequences are written 60 bases per line.
- **API:** `seqtk_randbase(input_path: 'str', output_path: 'str', *, seqtk_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program(s) of seqtk (1.5-r133 is the tested version) on PATH or through `seqtk_path`; it is not installed in the environment that produced this page. The example below shows behaviour that does not need it: the error raised for a missing program.
- **Input:** required: `input_path`, `output_path`; optional keyword-only: `seqtk_path` (meanings and defaults are in the Method text and the signature). Input files must exist; output paths must not exist yet.
- **Output:** The common result dictionary. Result keys: ``output_file`` and ``record_count``.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_randbase

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'in.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        seqtk_randbase(
            input_path=str(root / 'in.fa'),
            output_path=str(root / 'result.txt'),
            seqtk_path='/nonexistent/seqtk',
            verbose=False,
        )
    except DNAKitError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError BACKEND_UNAVAILABLE
```

**Usage with the program installed:**

```python
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_randbase

result = seqtk_randbase(
    input_path='in.fa',
    output_path='result.txt',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 14) `dnakit.fasta_fastq.seqtk_seq()`

- **Purpose:** Transform FASTA or FASTQ records with ``seqtk seq`` (comment removal, quality masking, upper-casing, mate selection).
- **Method:** This single function covers four operations of the work list through its parameters: ``drop_comments=True`` removes the comment after the first white space of every header (``seqtk seq -C``); ``mask_quality_below=Q`` masks bases whose Phred quality is lower than Q (``-q Q``, the bases become lower case unless ``mask_character`` gives the replacement character, ``-n``); ``uppercase=True`` converts all bases to upper case (``-U``); ``mate=1`` keeps the odd records (reads 1, ``-1``) and ``mate=2`` the even records (reads 2, ``-2``) of an interleaved paired-end file. ``quality_shift`` is the ASCII offset of the quality characters (``-Q``, default 33). ``line_length`` re-wraps sequences (``-l``; None keeps every sequence on one line). With no option the records are copied unchanged. ``input_path`` may be gzip-compressed.
- **API:** `seqtk_seq(input_path: 'str', output_path: 'str', *, drop_comments: 'bool' = False, mask_quality_below: 'int | None' = None, mask_character: 'str | None' = None, quality_shift: 'int' = 33, uppercase: 'bool' = False, mate: 'int | None' = None, line_length: 'int | None' = None, seqtk_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program(s) of seqtk (1.5-r133 is the tested version) on PATH or through `seqtk_path`; it is not installed in the environment that produced this page. The example below shows behaviour that does not need it: the error raised for a missing program.
- **Input:** required: `input_path`, `output_path`; optional keyword-only: `drop_comments`, `mask_quality_below`, `mask_character`, `quality_shift`, `uppercase`, `mate`, `line_length`, `seqtk_path` (meanings and defaults are in the Method text and the signature). Input files must exist; output paths must not exist yet.
- **Output:** The common result dictionary. Result keys: ``output_file`` and ``record_count``.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_seq

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'in.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        seqtk_seq(
            input_path=str(root / 'in.fa'),
            output_path=str(root / 'result.txt'),
            seqtk_path='/nonexistent/seqtk',
            verbose=False,
        )
    except DNAKitError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError BACKEND_UNAVAILABLE
```

**Usage with the program installed:**

```python
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_seq

result = seqtk_seq(
    input_path='in.fa',
    output_path='result.txt',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 15) `dnakit.fasta_fastq.seqtk_size()`

- **Purpose:** Count the sequences and the bases of a FASTA or FASTQ file with ``seqtk size``.
- **Method:** ``input_path`` may be gzip-compressed. The tab-separated line written by seqtk is stored in ``output_path`` and parsed into ``sequence_count`` and ``base_count``.
- **API:** `seqtk_size(input_path: 'str', output_path: 'str', *, seqtk_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program(s) of seqtk (1.5-r133 is the tested version) on PATH or through `seqtk_path`; it is not installed in the environment that produced this page. The example below shows behaviour that does not need it: the error raised for a missing program.
- **Input:** required: `input_path`, `output_path`; optional keyword-only: `seqtk_path` (meanings and defaults are in the Method text and the signature). Input files must exist; output paths must not exist yet.
- **Output:** The common result dictionary. No extra keys beyond `outputs` and `log_tail`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_size

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'in.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        seqtk_size(
            input_path=str(root / 'in.fa'),
            output_path=str(root / 'result.txt'),
            seqtk_path='/nonexistent/seqtk',
            verbose=False,
        )
    except DNAKitError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError BACKEND_UNAVAILABLE
```

**Usage with the program installed:**

```python
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_size

result = seqtk_size(
    input_path='in.fa',
    output_path='result.txt',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 16) `dnakit.fasta_fastq.seqtk_telo()`

- **Purpose:** Identify telomere repeats at the ends of assembly sequences or long reads with ``seqtk telo``.
- **Method:** ``motif`` is the repeat unit as written on the 5' end (``-m``, default ``CCCTAA``; the reverse complement is searched on the 3' end), ``penalty`` the score subtracted for a non-motif base (``-p``, default 1), ``max_drop`` the maximum score drop before extension stops (``-d``, default 2000) and ``min_score`` the minimum score to report (``-s``, default 300).
- **API:** `seqtk_telo(input_path: 'str', output_path: 'str', *, motif: 'str' = 'CCCTAA', penalty: 'int' = 1, max_drop: 'int' = 2000, min_score: 'int' = 300, max_records: 'int' = 10000, seqtk_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program(s) of seqtk (1.5-r133 is the tested version) on PATH or through `seqtk_path`; it is not installed in the environment that produced this page. The example below shows behaviour that does not need it: the error raised for a missing program.
- **Input:** required: `input_path`, `output_path`; optional keyword-only: `motif`, `penalty`, `max_drop`, `min_score`, `max_records`, `seqtk_path` (meanings and defaults are in the Method text and the signature). Input files must exist; output paths must not exist yet.
- **Output:** The common result dictionary. Result keys: ``records`` (``seqid``, ``start``, ``end``, ``length``; at most ``max_records``), ``record_count``, ``telomere_bases`` and ``total_bases`` (the totals seqtk prints to standard error).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_telo

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'in.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        seqtk_telo(
            input_path=str(root / 'in.fa'),
            output_path=str(root / 'result.txt'),
            seqtk_path='/nonexistent/seqtk',
            verbose=False,
        )
    except DNAKitError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError BACKEND_UNAVAILABLE
```

**Usage with the program installed:**

```python
from dnakit.fasta_fastq.seqtk_p12_cli_functions import seqtk_telo

result = seqtk_telo(
    input_path='in.fa',
    output_path='result.txt',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 17) `dnakit.sequencing_qc.nanoplot_plots()`

- **Purpose:** Draw the NanoPlot long-read plots and statistics of one dataset and return the plot files by name (``NanoPlot``).
- **Method:** ``input_paths`` are one or more files of ``input_type``: ``fastq``, ``fasta``, ``fastq_rich`` and ``fastq_minimal`` (reads), ``summary`` (a Guppy/MinKNOW sequencing summary, which also gives the time and channel plots), ``bam``, ``ubam``, ``cram`` (alignments, which also give the alignment plots) or ``feather``. ``output_dir`` must not exist. NanoPlot writes every plot the input allows, as interactive HTML and, when ``static_images`` is true (the NanoPlot default), as static images in ``image_formats`` (default PNG); with ``static_images=False`` only HTML is written, which is much faster. The result groups the files under these plot names, which cover the plots of the work list: ``read_length_histogram`` (weighted and non-weighted, plus log-transformed with ``log_length``), ``yield_by_minimal_length``, ``read_length_vs_quality_scatter``, ``aligned_vs_sequenced_length_scatter``, ``mapq_vs_base_quality_scatter``, ``mapq_vs_read_length_scatter``, ``percent_identity_vs_base_quality_scatter``, ``percent_identity_vs_aligned_length_scatter``, ``percent_identity_histogram``, ``reads_per_channel_activity_map``, ``reads_over_time``, ``active_pores_over_time``, ``cumulative_yield_gigabases``, ``cumulative_yield_reads``, ``read_length_over_time_violin``, ``quality_over_time_violin`` and ``sequencing_speed_over_time_violin``. ``plot_names`` limits what is reported (NanoPlot still writes all of them). ``bivariate_plots`` selects the bivariate plot kinds among ``kde``, ``hex`` and ``dot`` (``--plots``). Filters and transforms: ``min_length``/``max_length``, ``min_quality``, ``downsample``, ``drop_outliers``, ``no_supplementary``, ``log_length`` (``--loglength``), ``percent_quality`` (``--percentqual``), ``aligned_length`` (``--alength``) and ``show_n50`` (True ``--N50``, False ``--no-N50``). Plot style: ``color``, ``colormap``, ``title`` and ``dpi``. ``raw_table=True`` also writes the per-read table (``--raw``: read id, quality, length, and for alignments aligned length, mapping quality) as ``NanoPlot-data.tsv.gz``.
- **API:** `nanoplot_plots(input_paths: 'Sequence[str]', output_dir: 'str', *, input_type: 'str' = 'fastq', plot_names: 'Sequence[str] | None' = None, bivariate_plots: 'Sequence[str] | None' = None, static_images: 'bool' = True, image_formats: 'Sequence[str] | None' = None, log_length: 'bool' = False, percent_quality: 'bool' = False, aligned_length: 'bool' = False, min_length: 'int | None' = None, max_length: 'int | None' = None, min_quality: 'float | None' = None, downsample: 'int | None' = None, drop_outliers: 'bool' = False, no_supplementary: 'bool' = False, show_n50: 'bool | None' = None, color: 'str | None' = None, colormap: 'str | None' = None, title: 'str | None' = None, dpi: 'int | None' = None, raw_table: 'bool' = False, threads: 'int' = 1, timeout_seconds: 'float' = 86400.0, nanoplot_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program(s) of NanoPlot (1.48.0 is the tested version) on PATH or through `nanoplot_path`; it is not installed in the environment that produced this page. The example below shows behaviour that does not need it: the error raised for a missing program.
- **Input:** required: `input_paths`, `output_dir`; optional keyword-only: `input_type`, `plot_names`, `bivariate_plots`, `static_images`, `image_formats`, `log_length`, `percent_quality`, `aligned_length`, `min_length`, `max_length`, `min_quality`, `downsample`, `drop_outliers`, `no_supplementary`, `show_n50`, `color`, `colormap`, `title`, `dpi`, `raw_table`, `threads`, `timeout_seconds`, `nanoplot_path` (meanings and defaults are in the Method text and the signature). Input files must exist; output paths must not exist yet.
- **Output:** The common result dictionary. Result keys: ``output_dir``, ``plots`` (plot name to the list of files written), ``missing_plots`` (requested names for which NanoPlot wrote nothing, for example the time plots of a FASTQ input or the mapping-quality plots when the mapping qualities do not vary), ``statistics`` (the ``NanoStats`` table as a dictionary: ``number_of_reads``, ``number_of_bases``, ``median_read_length``, ``mean_read_length``, ``n50``, ``mean_qual``, ``median_qual``, ...), ``stats_file``, ``report_file`` and, with ``raw_table``, ``raw_table_file``, ``raw_table_columns`` and ``raw_table_rows``.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.sequencing_qc.nanoplot_p12_cli_functions import nanoplot_plots

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'in.fastq').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        nanoplot_plots(
            input_paths=[str(root / 'in.fastq')],
            output_dir=str(root / 'out'),
            nanoplot_path='/nonexistent/NanoPlot',
            verbose=False,
        )
    except DNAKitError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError BACKEND_UNAVAILABLE
```

**Usage with the program installed:**

```python
from dnakit.sequencing_qc.nanoplot_p12_cli_functions import nanoplot_plots

result = nanoplot_plots(
    input_paths=['in.fastq'],
    output_dir='out',
    verbose=False,
)
print(result['status'], result['outputs'])
```
