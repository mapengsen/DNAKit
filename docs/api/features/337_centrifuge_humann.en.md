# Centrifuge and HUMAnN: metagenomic reference and functional-profile command-line operations

This group of functions wraps the external command-line programs [Centrifuge](https://github.com/DaehwanKimLab/centrifuge) (9 functions, module `centrifuge_p12_cli_functions` in `dnakit.microbiome`), [HUMAnN](https://github.com/biobakery/humann) (13 functions, module `humann_p12_cli_functions` in `dnakit.microbiome`). Every function runs one fixed program or sub-command with an explicit, validated argument list (no shell, timeout and output-size limits), writes results to explicit new paths and never overwrites an existing output. Tested versions: Centrifuge 1.0.4, HUMAnN 3.9. Import them with `from <area> import ...` once the area package exports them, or directly from the module.

The programs are not bundled with DNAKit: install them (Centrifuge: `conda install -c bioconda centrifuge`; HUMAnN: `conda install -c bioconda humann` or `pip install humann`), then put them on PATH or pass the path option of the function. If a program cannot be found, `BackendUnavailableError` is raised with an install hint; if it exits with an error, `BackendExecutionError` is raised with the tail of its log. Results are JSON-serialisable dictionaries with `status`, `backend`, `version`, `command`, `return_code`, `outputs`, `log_tail`, `provenance` and `operation`; each section lists only the extra keys. The programs are not installed in the environment used to produce this page, so every checked example shows behaviour that is real on any machine (a missing program); each is followed by the call to run with the program installed.

## 1) `dnakit.microbiome.centrifuge_build()`

- **Purpose:** Build a Centrifuge classification index from reference FASTA files and a taxonomy (``centrifuge-build``).
- **Method:** ``reference_files`` are DNA FASTA files; ``conversion_table`` maps every sequence id to a taxonomy id (two tab separated columns), ``taxonomy_tree`` is the NCBI ``nodes.dmp`` file and ``name_table`` (optional) the NCBI ``names.dmp`` file; ``size_table`` is an optional table of genome sizes. ``index_prefix`` is the path prefix of the four index files ``<prefix>.1.cf`` to ``<prefix>.4.cf`` (none may exist yet). ``threads`` is ``-p``, ``offrate`` ``-o`` (the suffix array is sampled every 2^offrate characters), ``ftab_chars`` ``-t``, ``bmax``, ``bmax_divisor`` (``--bmaxdivn``), ``diff_cover_period`` (``--dcv``) and ``seed`` are passed on unchanged. ``kmer_count=K`` also counts the distinct k-mers of the reference for k from 1 to K (``--kmer-count K``, K at most 64) and returns them in ``distinct_kmers`` (a mapping from k, as text, to the count).
- **API:** `centrifuge_build(reference_files: 'Sequence[str]', index_prefix: 'str', conversion_table: 'str', taxonomy_tree: 'str', *, name_table: 'str | None' = None, size_table: 'str | None' = None, threads: 'int' = 1, offrate: 'int | None' = None, ftab_chars: 'int | None' = None, bmax: 'int | None' = None, bmax_divisor: 'int | None' = None, diff_cover_period: 'int | None' = None, seed: 'int | None' = None, kmer_count: 'int | None' = None, centrifuge_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program(s) of Centrifuge (1.0.4 is the tested version) on PATH or through `centrifuge_dir`; it is not installed in the environment that produced this page. The example below shows behaviour that does not need it: the error raised for a missing program.
- **Input:** required: `reference_files`, `index_prefix`, `conversion_table`, `taxonomy_tree`; optional keyword-only: `name_table`, `size_table`, `threads`, `offrate`, `ftab_chars`, `bmax`, `bmax_divisor`, `diff_cover_period`, `seed`, `kmer_count`, `centrifuge_dir` (meanings and defaults are in the Method text and the signature). Input files must exist; output paths must not exist yet.
- **Output:** The common result dictionary. Result keys: ``index_prefix`` and ``index_files``.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.microbiome.centrifuge_p12_cli_functions import centrifuge_build

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'conv.tsv').write_text('x\n')
    (root / 'nodes.dmp').write_text('x\n')
    (root / 'ref.fa').write_text('x\n')
    for suffix in ('.1.cf', '.2.cf', '.3.cf', '.4.cf'):
        (root / f'idx{suffix}').write_text('x')
    (root / 'bin').mkdir()
    try:
        centrifuge_build(
            reference_files=[str(root / 'ref.fa')],
            index_prefix=str(root / 'idx'),
            conversion_table=str(root / 'conv.tsv'),
            taxonomy_tree=str(root / 'nodes.dmp'),
            centrifuge_dir=str(root / 'bin'),
            verbose=False,
        )
    except DNAKitError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
ConfigurationError CONFIGURATION_ERROR
```

**Usage with the program installed:**

```python
from dnakit.microbiome.centrifuge_p12_cli_functions import centrifuge_build

result = centrifuge_build(
    reference_files=['ref.fa'],
    index_prefix='idx',
    conversion_table='conv.tsv',
    taxonomy_tree='nodes.dmp',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 2) `dnakit.microbiome.centrifuge_compress_genomes()`

- **Purpose:** Compress a reference genome library by merging near-identical genomes of a taxonomic level (``centrifuge-compress.pl``).
- **Method:** ``download_dir`` holds the genome FASTA files (``.fna``, ``.fa``, ``.ffn``, ``.fasta``), ``taxonomy_dir`` the NCBI ``nodes.dmp`` and ``names.dmp`` (and ``gi_taxid_nucl.dmp`` unless ``mapping_file``, a two-column file mapping sequence ids to taxonomy ids, is given). ``level`` is the taxonomic level at which genomes are merged (``-l``, default species), ``max_genome_size`` the largest genome that is compressed (``--maxG``), ``no_compress`` only collects and cleans the sequences (``-noCompress``) and ``no_dustmasker`` skips low-complexity masking (``-noDustmasker``). The Perl program needs Jellyfish 2 and MUMmer (``nucmer``) on PATH (and ``dustmasker`` for masking). ``output_dir`` must not exist; the program runs inside it and writes ``<output_name>.fa``, ``<output_name>.map`` and ``<output_name>.size`` there.
- **API:** `centrifuge_compress_genomes(download_dir: 'str', taxonomy_dir: 'str', output_dir: 'str', *, output_name: 'str' = 'compressed', level: 'str' = 'species', mapping_file: 'str | None' = None, threads: 'int' = 1, max_genome_size: 'int' = 50000000, no_compress: 'bool' = False, no_dustmasker: 'bool' = False, centrifuge_dir: 'str | None' = None, timeout_seconds: 'float' = 86400.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program(s) of Centrifuge (1.0.4 is the tested version) on PATH or through `centrifuge_dir`; it is not installed in the environment that produced this page. The example below shows behaviour that does not need it: the error raised for a missing program.
- **Input:** required: `download_dir`, `taxonomy_dir`, `output_dir`; optional keyword-only: `output_name`, `level`, `mapping_file`, `threads`, `max_genome_size`, `no_compress`, `no_dustmasker`, `centrifuge_dir`, `timeout_seconds` (meanings and defaults are in the Method text and the signature). Input files must exist; output paths must not exist yet.
- **Output:** The common result dictionary. Result keys: ``output_dir``, ``fasta``, ``taxid_map``, ``size_table``, ``sequence_count`` and ``base_count``.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.microbiome.centrifuge_p12_cli_functions import centrifuge_compress_genomes

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'download').mkdir()
    (root / 'taxonomy').mkdir()
    (root / 'bin').mkdir()
    try:
        centrifuge_compress_genomes(
            download_dir=str(root / 'download'),
            taxonomy_dir=str(root / 'taxonomy'),
            output_dir=str(root / 'out'),
            centrifuge_dir=str(root / 'bin'),
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
from dnakit.microbiome.centrifuge_p12_cli_functions import centrifuge_compress_genomes

result = centrifuge_compress_genomes(
    download_dir='download',
    taxonomy_dir='taxonomy',
    output_dir='out',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 3) `dnakit.microbiome.centrifuge_download()`

- **Purpose:** Download reference sequences or the NCBI taxonomy for a Centrifuge index (``centrifuge-download``; needs the network).
- **Method:** ``database`` selects the source: ``refseq`` or ``genbank`` download genome FASTA files of the listed ``domains`` (one or more of bacteria, viral, archaea, fungi, protozoa, invertebrate, plant, vertebrate_mammalian and vertebrate_other); ``contaminants`` downloads the UniVec and EmVec vector sequences, and ``taxonomy`` the NCBI ``nodes.dmp`` and ``names.dmp``. The options of the genome databases are ``assembly_level`` (``-a``; one of "Complete Genome", "Chromosome", "Scaffold", "Contig" or "Any"), ``refseq_category`` (``-c``), ``taxids`` (``-t``), ``include_rna`` to also download RNA sequences (``-r``), ``filter_unplaced`` (``-u``), ``mask_low_complexity`` to mask low-complexity regions with ``dustmasker`` (``-m``, needs BLAST+), ``taxid_in_headers`` to put the taxonomy id into every header (``-l``), ``threads`` (``-P``) and ``download_program`` (``-g``: rsync, curl or wget). ``output_dir`` must not exist; the program writes into it. Centrifuge prints the sequence-id to taxonomy-id table to standard output; it is saved as ``seqid_to_taxid.map`` in ``output_dir`` (``mapping_file`` may name another new path).
- **API:** `centrifuge_download(output_dir: 'str', database: 'str', *, domains: 'Sequence[str] | None' = None, assembly_level: 'str' = 'Complete Genome', refseq_category: 'str | None' = None, taxids: 'Sequence[int] | None' = None, include_rna: 'bool' = False, filter_unplaced: 'bool' = False, mask_low_complexity: 'bool' = False, taxid_in_headers: 'bool' = False, threads: 'int' = 1, download_program: 'str | None' = None, mapping_file: 'str | None' = None, timeout_seconds: 'float' = 86400.0, centrifuge_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program(s) of Centrifuge (1.0.4 is the tested version) on PATH or through `centrifuge_dir`; it is not installed in the environment that produced this page. The example below shows behaviour that does not need it: the error raised for a missing program.
- **Input:** required: `output_dir`, `database`; optional keyword-only: `domains`, `assembly_level`, `refseq_category`, `taxids`, `include_rna`, `filter_unplaced`, `mask_low_complexity`, `taxid_in_headers`, `threads`, `download_program`, `mapping_file`, `timeout_seconds`, `centrifuge_dir` (meanings and defaults are in the Method text and the signature). Input files must exist; output paths must not exist yet.
- **Output:** The common result dictionary. Result keys: ``output_dir``, ``mapping_file``, ``mapping_rows`` (number of rows) and ``files`` (the downloaded files, at most 1000, as paths relative to ``output_dir``).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.microbiome.centrifuge_p12_cli_functions import centrifuge_download

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'bin').mkdir()
    try:
        centrifuge_download(
            output_dir=str(root / 'out'),
            database='taxonomy',
            centrifuge_dir=str(root / 'bin'),
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
from dnakit.microbiome.centrifuge_p12_cli_functions import centrifuge_download

result = centrifuge_download(
    output_dir='out',
    database='taxonomy',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 4) `dnakit.microbiome.centrifuge_inspect()`

- **Purpose:** Extract information from a Centrifuge index (``centrifuge-inspect``).
- **Method:** ``output`` selects what is written to ``output_path`` and parsed: ``fasta`` the indexed reference sequences as FASTA with ``line_width`` bases per line (``-a``, default 60; ``record_count`` and ``base_count`` are returned), ``names`` the reference sequence names (``-n``), ``summary`` the index parameters and the sequences with their lengths (``-s``; returned as ``parameters`` and ``records``), ``conversion_table`` the sequence-to-taxon table (``--conversion-table``), ``taxonomy_tree`` (``--taxonomy-tree``), ``name_table`` (``--name-table``), ``size_table`` the summed reference length per taxon (``--size-table``) and ``estimate_n_kmers`` the approximate number of distinct k-mers of the reference (``--estimate-n-kmers``; returned as ``estimated_kmers``). Tables are returned as ``records`` (at most ``max_records``) with ``record_count``.
- **API:** `centrifuge_inspect(index_prefix: 'str', output_path: 'str', *, output: 'str' = 'fasta', line_width: 'int' = 60, max_records: 'int' = 10000, centrifuge_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program(s) of Centrifuge (1.0.4 is the tested version) on PATH or through `centrifuge_dir`; it is not installed in the environment that produced this page. The example below shows behaviour that does not need it: the error raised for a missing program.
- **Input:** required: `index_prefix`, `output_path`; optional keyword-only: `output`, `line_width`, `max_records`, `centrifuge_dir` (meanings and defaults are in the Method text and the signature). Input files must exist; output paths must not exist yet.
- **Output:** The common result dictionary. No extra keys beyond `outputs` and `log_tail`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.microbiome.centrifuge_p12_cli_functions import centrifuge_inspect

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    for suffix in ('.1.cf', '.2.cf', '.3.cf', '.4.cf'):
        (root / f'idx{suffix}').write_text('x')
    (root / 'bin').mkdir()
    try:
        centrifuge_inspect(
            index_prefix=str(root / 'idx'),
            output_path=str(root / 'result.txt'),
            centrifuge_dir=str(root / 'bin'),
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
from dnakit.microbiome.centrifuge_p12_cli_functions import centrifuge_inspect

result = centrifuge_inspect(
    index_prefix='idx',
    output_path='result.txt',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 5) `dnakit.microbiome.centrifuge_kreport()`

- **Purpose:** Create a Kraken-style taxonomic report from Centrifuge classification output (``centrifuge-kreport``).
- **Method:** ``index_prefix`` is the index used for the classification (its taxonomy is read), ``classification_files`` the Centrifuge ``-S`` output files. ``min_score`` and ``min_length`` keep only hits with at least that score or hit length; ``no_lca`` reports count fractions at the taxa instead of the lowest common ancestor of multiple assignments; ``show_zeros`` also lists clades without reads; ``is_count_table`` reads files with the two columns ``taxid`` and ``count`` instead. The report is saved in ``output_path``.
- **API:** `centrifuge_kreport(index_prefix: 'str', classification_files: 'Sequence[str]', output_path: 'str', *, min_score: 'int | None' = None, min_length: 'int | None' = None, no_lca: 'bool' = False, show_zeros: 'bool' = False, is_count_table: 'bool' = False, max_records: 'int' = 10000, centrifuge_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program(s) of Centrifuge (1.0.4 is the tested version) on PATH or through `centrifuge_dir`; it is not installed in the environment that produced this page. The example below shows behaviour that does not need it: the error raised for a missing program.
- **Input:** required: `index_prefix`, `classification_files`, `output_path`; optional keyword-only: `min_score`, `min_length`, `no_lca`, `show_zeros`, `is_count_table`, `max_records`, `centrifuge_dir` (meanings and defaults are in the Method text and the signature). Input files must exist; output paths must not exist yet.
- **Output:** The common result dictionary. Result keys: ``records`` (``percent``, ``clade_reads``, ``taxon_reads``, ``rank`` code, ``taxid``, ``name`` and ``depth``; at most ``max_records``) and ``record_count``.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.microbiome.centrifuge_p12_cli_functions import centrifuge_kreport

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'cls.tsv').write_text('x\n')
    for suffix in ('.1.cf', '.2.cf', '.3.cf', '.4.cf'):
        (root / f'idx{suffix}').write_text('x')
    (root / 'bin').mkdir()
    try:
        centrifuge_kreport(
            index_prefix=str(root / 'idx'),
            classification_files=[str(root / 'cls.tsv')],
            output_path=str(root / 'result.txt'),
            centrifuge_dir=str(root / 'bin'),
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
from dnakit.microbiome.centrifuge_p12_cli_functions import centrifuge_kreport

result = centrifuge_kreport(
    index_prefix='idx',
    classification_files=['cls.tsv'],
    output_path='result.txt',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 6) `dnakit.microbiome.centrifuge_remove_empty_sequences()`

- **Purpose:** Remove the FASTA records that have an empty sequence, as left behind by masking (``centrifuge-RemoveEmptySequence.pl``).
- **Method:** Remove the FASTA records that have an empty sequence, as left behind by masking (``centrifuge-RemoveEmptySequence.pl``).
- **API:** `centrifuge_remove_empty_sequences(input_path: 'str', output_path: 'str', *, centrifuge_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program(s) of Centrifuge (1.0.4 is the tested version) on PATH or through `centrifuge_dir`; it is not installed in the environment that produced this page. The example below shows behaviour that does not need it: the error raised for a missing program.
- **Input:** required: `input_path`, `output_path`; optional keyword-only: `centrifuge_dir` (meanings and defaults are in the Method text and the signature). Input files must exist; output paths must not exist yet.
- **Output:** The common result dictionary. Result keys: ``output_file``, ``record_count`` and ``base_count`` of the cleaned file.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.microbiome.centrifuge_p12_cli_functions import centrifuge_remove_empty_sequences

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'in.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        centrifuge_remove_empty_sequences(
            input_path=str(root / 'in.fa'),
            output_path=str(root / 'result.txt'),
            centrifuge_dir=str(root / 'bin'),
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
from dnakit.microbiome.centrifuge_p12_cli_functions import centrifuge_remove_empty_sequences

result = centrifuge_remove_empty_sequences(
    input_path='in.fa',
    output_path='result.txt',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 7) `dnakit.microbiome.centrifuge_remove_n()`

- **Purpose:** Delete the N bases of every FASTA record and re-wrap the sequence to 80 columns (``centrifuge-RemoveN.pl``).
- **Method:** Both ``N`` and ``n`` are removed; other bases keep their case. A record that consists only of Ns keeps its header with no sequence line.
- **API:** `centrifuge_remove_n(input_path: 'str', output_path: 'str', *, centrifuge_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program(s) of Centrifuge (1.0.4 is the tested version) on PATH or through `centrifuge_dir`; it is not installed in the environment that produced this page. The example below shows behaviour that does not need it: the error raised for a missing program.
- **Input:** required: `input_path`, `output_path`; optional keyword-only: `centrifuge_dir` (meanings and defaults are in the Method text and the signature). Input files must exist; output paths must not exist yet.
- **Output:** The common result dictionary. Result keys: ``output_file``, ``record_count`` and ``base_count``.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.microbiome.centrifuge_p12_cli_functions import centrifuge_remove_n

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'in.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        centrifuge_remove_n(
            input_path=str(root / 'in.fa'),
            output_path=str(root / 'result.txt'),
            centrifuge_dir=str(root / 'bin'),
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
from dnakit.microbiome.centrifuge_p12_cli_functions import centrifuge_remove_n

result = centrifuge_remove_n(
    input_path='in.fa',
    output_path='result.txt',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 8) `dnakit.microbiome.centrifuge_shared_sequence()`

- **Purpose:** Detect sequence segments shared by several genomes and merge them (``centrifuge-BuildSharedSequence.pl``).
- **Method:** ``fasta_files`` are the genomes of one taxon (at least two). The program counts k-mers with Jellyfish, picks the pair of genomes that share most k-mers and aligns them with MUMmer ``nucmer``; both must be on PATH. ``kmer_size`` is ``-kmerSize`` (default 53), ``kmer_portion`` the fraction of k-mers sampled (``-kmerPortion``, default 0.01), ``nucmer_identity`` the minimum percent identity of shared segments (``-nucmerIdy``, default 99), ``overlap`` the bases kept at the ends of a shared segment (``-overlap``, default 250) and ``fragment`` keeps the shared segments as separate fragments instead of one chunk (``-fragment``). ``output_dir`` must not exist; the program runs inside it.
- **API:** `centrifuge_shared_sequence(fasta_files: 'Sequence[str]', output_dir: 'str', *, prefix: 'str' = 'shared', kmer_size: 'int' = 53, kmer_portion: 'float' = 0.01, nucmer_identity: 'int' = 99, overlap: 'int' = 250, fragment: 'bool' = False, centrifuge_dir: 'str | None' = None, timeout_seconds: 'float' = 86400.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program(s) of Centrifuge (1.0.4 is the tested version) on PATH or through `centrifuge_dir`; it is not installed in the environment that produced this page. The example below shows behaviour that does not need it: the error raised for a missing program.
- **Input:** required: `fasta_files`, `output_dir`; optional keyword-only: `prefix`, `kmer_size`, `kmer_portion`, `nucmer_identity`, `overlap`, `fragment`, `centrifuge_dir`, `timeout_seconds` (meanings and defaults are in the Method text and the signature). Input files must exist; output paths must not exist yet.
- **Output:** The common result dictionary. Result keys: ``output_dir``, ``files`` (the ``<prefix>_<n>.fa`` FASTA files that were written, with their ``records`` and ``bases``) and ``program_output`` (the progress lines).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.microbiome.centrifuge_p12_cli_functions import centrifuge_shared_sequence

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'g1.fa').write_text('x\n')
    (root / 'g2.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        centrifuge_shared_sequence(
            fasta_files=[str(root / 'g1.fa'), str(root / 'g2.fa')],
            output_dir=str(root / 'out'),
            centrifuge_dir=str(root / 'bin'),
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
from dnakit.microbiome.centrifuge_p12_cli_functions import centrifuge_shared_sequence

result = centrifuge_shared_sequence(
    fasta_files=['g1.fa', 'g2.fa'],
    output_dir='out',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 9) `dnakit.microbiome.centrifuge_sort_nt()`

- **Purpose:** Sort the records of a nucleotide FASTA by taxonomy id (``centrifuge-sort-nt.pl``).
- **Method:** ``sequence_file`` is a nucleotide FASTA whose headers start with an accession with version (for example the NCBI nt database), ``mapping_files`` are one or more NCBI ``accession2taxid`` files (plain or ``.gz``; four columns: accession, accession.version, taxid, gi). The sorted FASTA goes to ``output_path``. ``taxid_map_output`` (``-m``) receives the accession-to-taxid table of the sequences that were found and ``unmapped_output`` (``-a``) the file the upstream script writes when some accession has no mapping (Centrifuge 1.0.4 writes the accessions that do have a mapping into it, a defect of the upstream script that this facade does not alter; the file is not required to exist afterwards).
- **API:** `centrifuge_sort_nt(sequence_file: 'str', mapping_files: 'Sequence[str]', output_path: 'str', *, taxid_map_output: 'str | None' = None, unmapped_output: 'str | None' = None, centrifuge_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program(s) of Centrifuge (1.0.4 is the tested version) on PATH or through `centrifuge_dir`; it is not installed in the environment that produced this page. The example below shows behaviour that does not need it: the error raised for a missing program.
- **Input:** required: `sequence_file`, `mapping_files`, `output_path`; optional keyword-only: `taxid_map_output`, `unmapped_output`, `centrifuge_dir` (meanings and defaults are in the Method text and the signature). Input files must exist; output paths must not exist yet.
- **Output:** The common result dictionary. Result keys: ``output_file``, ``record_count`` and ``base_count``.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.microbiome.centrifuge_p12_cli_functions import centrifuge_sort_nt

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'a2t.txt').write_text('x\n')
    (root / 'nt.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        centrifuge_sort_nt(
            sequence_file=str(root / 'nt.fa'),
            mapping_files=[str(root / 'a2t.txt')],
            output_path=str(root / 'result.txt'),
            centrifuge_dir=str(root / 'bin'),
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
from dnakit.microbiome.centrifuge_p12_cli_functions import centrifuge_sort_nt

result = centrifuge_sort_nt(
    sequence_file='nt.fa',
    mapping_files=['a2t.txt'],
    output_path='result.txt',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 10) `dnakit.microbiome.humann_barplot()`

- **Purpose:** Plot the taxon-stratified contributions of one gene family or pathway across samples (``humann_barplot``).
- **Method:** ``input_path`` is a stratified HUMAnN table (TSV or BIOM) and ``focal_feature`` the feature id to plot (the id, not the full name). The figure is written to ``output_image``; its extension selects the format (png, pdf, svg). ``top_taxa`` is the largest number of taxa highlighted (default 18), ``as_genera`` collapses species to genera, ``exclude_unclassified`` drops the unclassified stratum, ``remove_zeros`` skips samples with a zero total for the feature. ``sort`` lists sample sorting methods in order (none, sum, dominant, braycurtis, braycurtis_w, metadata); ``focal_metadata`` and ``last_metadata`` use metadata rows of the table; ``scaling`` is ``original``, ``logstack`` or ``totalsum``; ``units`` names the y axis and ``dimensions`` is the width and height in inches. ``sample_order_output`` saves the sample order.
- **API:** `humann_barplot(input_path: 'str', focal_feature: 'str', output_image: 'str', *, top_taxa: 'int' = 18, as_genera: 'bool' = False, exclude_unclassified: 'bool' = False, remove_zeros: 'bool' = False, sort: 'Sequence[str] | None' = None, focal_metadata: 'str | None' = None, last_metadata: 'str | None' = None, scaling: 'str' = 'original', units: 'str' = 'generic', dimensions: 'Sequence[float]' = (11.0, 6.0), sample_order_output: 'str | None' = None, humann_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program(s) of HUMAnN (3.9 is the tested version) on PATH or through `humann_dir`; it is not installed in the environment that produced this page. The example below shows behaviour that does not need it: the error raised for a missing program.
- **Input:** required: `input_path`, `focal_feature`, `output_image`; optional keyword-only: `top_taxa`, `as_genera`, `exclude_unclassified`, `remove_zeros`, `sort`, `focal_metadata`, `last_metadata`, `scaling`, `units`, `dimensions`, `sample_order_output`, `humann_dir` (meanings and defaults are in the Method text and the signature). Input files must exist; output paths must not exist yet.
- **Output:** The common result dictionary. Result keys: ``image`` (the figure path, with its size in ``image_bytes``) and ``sample_order_file``. The plotted values are in the input table.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.microbiome.humann_p12_cli_functions import humann_barplot

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'in.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        humann_barplot(
            input_path=str(root / 'in.fa'),
            focal_feature='UniRef90_X',
            output_image=str(root / 'plot.png'),
            humann_dir=str(root / 'bin'),
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
from dnakit.microbiome.humann_p12_cli_functions import humann_barplot

result = humann_barplot(
    input_path='in.fa',
    focal_feature='UniRef90_X',
    output_image='plot.png',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 11) `dnakit.microbiome.humann_build_custom_database()`

- **Purpose:** Build a custom nucleotide or protein search database from gene sequences, optionally limited to genera (``humann_build_custom_database``).
- **Method:** ``fasta_path`` is a FASTA file of gene sequences. ``id_mapping`` maps FASTA ids to genus names and ``taxonomic_profile`` is a MetaPhlAn profile; with both, only the genera whose abundance reaches ``genus_abundance_threshold`` (``--genus-abundance-threshold``) are kept. ``database_format`` is ``fasta`` (a filtered FASTA file, the default) or ``diamond`` (also builds a DIAMOND index, which needs ``diamond`` on PATH). ``output_dir`` must not exist.
- **API:** `humann_build_custom_database(fasta_path: 'str', output_dir: 'str', *, id_mapping: 'str | None' = None, taxonomic_profile: 'str | None' = None, database_format: 'str' = 'fasta', genus_abundance_threshold: 'float | None' = None, humann_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program(s) of HUMAnN (3.9 is the tested version) on PATH or through `humann_dir`; it is not installed in the environment that produced this page. The example below shows behaviour that does not need it: the error raised for a missing program.
- **Input:** required: `fasta_path`, `output_dir`; optional keyword-only: `id_mapping`, `taxonomic_profile`, `database_format`, `genus_abundance_threshold`, `humann_dir` (meanings and defaults are in the Method text and the signature). Input files must exist; output paths must not exist yet.
- **Output:** The common result dictionary. Result keys: ``output_dir`` and ``database_files`` (the files written, sorted).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.microbiome.humann_p12_cli_functions import humann_build_custom_database

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'genes.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        humann_build_custom_database(
            fasta_path=str(root / 'genes.fa'),
            output_dir=str(root / 'out'),
            humann_dir=str(root / 'bin'),
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
from dnakit.microbiome.humann_p12_cli_functions import humann_build_custom_database

result = humann_build_custom_database(
    fasta_path='genes.fa',
    output_dir='out',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 12) `dnakit.microbiome.humann_expand_cluster()`

- **Purpose:** Expand one UniRef90 gene family to all UniRef90 families in its UniRef50 cluster (``humann_expand_cluster``).
- **Method:** ``input_path`` is a UniRef90 gene family table (TSV) and ``uniref90`` the family of interest (with or without the ``UniRef90_`` prefix). The output table keeps the input rows of every UniRef90 that clusters with it. The program reads the UniRef90-to-UniRef50 mapping of the HUMAnN utility mapping data.
- **API:** `humann_expand_cluster(input_path: 'str', uniref90: 'str', output_path: 'str', *, humann_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program(s) of HUMAnN (3.9 is the tested version) on PATH or through `humann_dir`; it is not installed in the environment that produced this page. The example below shows behaviour that does not need it: the error raised for a missing program.
- **Input:** required: `input_path`, `uniref90`, `output_path`; optional keyword-only: `humann_dir` (meanings and defaults are in the Method text and the signature). Input files must exist; output paths must not exist yet.
- **Output:** The common result dictionary. Result keys: ``table`` (``file``, ``columns``, ``row_count`` and the first ``rows``).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.microbiome.humann_p12_cli_functions import humann_expand_cluster

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'in.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        humann_expand_cluster(
            input_path=str(root / 'in.fa'),
            uniref90='UniRef90_X',
            output_path=str(root / 'result.txt'),
            humann_dir=str(root / 'bin'),
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
from dnakit.microbiome.humann_p12_cli_functions import humann_expand_cluster

result = humann_expand_cluster(
    input_path='in.fa',
    uniref90='UniRef90_X',
    output_path='result.txt',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 13) `dnakit.microbiome.humann_infer_taxonomy()`

- **Purpose:** Infer the taxonomy of unclassified gene families from the lowest common ancestor of their UniRef cluster (``humann_infer_taxonomy``).
- **Method:** ``input_path`` is a HUMAnN gene family table. The LCA data come from ``database`` (the name of a tree-of-life/LCA file of the HUMAnN utility mapping data, for example ``uniref90`` or ``uniref50``) or from ``database_file`` (a path, ``--devdb``); give exactly one. ``level`` is the target taxonomic level (Kingdom, Phylum, Class, Order, Family or Genus), ``mode`` selects the rows used (``totals``, ``unclassified`` or ``stratified``), ``lca_choice`` the per-gene annotation (``source_tax``, ``uniref_lca`` or ``humann_lca``) and ``threshold`` the minimum frequency of a new taxon.
- **API:** `humann_infer_taxonomy(input_path: 'str', output_path: 'str', *, level: 'str' = 'Family', database: 'str | None' = None, database_file: 'str | None' = None, mode: 'str' = 'totals', lca_choice: 'str' = 'humann_lca', threshold: 'float' = 0.001, humann_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program(s) of HUMAnN (3.9 is the tested version) on PATH or through `humann_dir`; it is not installed in the environment that produced this page. The example below shows behaviour that does not need it: the error raised for a missing program.
- **Input:** required: `input_path`, `output_path`; optional keyword-only: `level`, `database`, `database_file`, `mode`, `lca_choice`, `threshold`, `humann_dir` (meanings and defaults are in the Method text and the signature). Input files must exist; output paths must not exist yet.
- **Output:** The common result dictionary. Result keys: ``table`` (``file``, ``columns``, ``row_count`` and the first ``rows``) and ``summary`` (the stratification summary line of the program, when present).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.microbiome.humann_p12_cli_functions import humann_infer_taxonomy

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'in.fa').write_text('x\n')
    (root / 'tol.tsv').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        humann_infer_taxonomy(
            input_path=str(root / 'in.fa'),
            output_path=str(root / 'result.txt'),
            database_file=str(root / 'tol.tsv'),
            humann_dir=str(root / 'bin'),
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
from dnakit.microbiome.humann_p12_cli_functions import humann_infer_taxonomy

result = humann_infer_taxonomy(
    input_path='in.fa',
    output_path='result.txt',
    database_file='tol.tsv',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 14) `dnakit.microbiome.humann_join_tables()`

- **Purpose:** Join the per-sample gene, pathway or taxonomy tables of a folder into one multi-sample table (``humann_join_tables``).
- **Method:** Every table in ``input_dir`` (or only those whose file name contains ``file_name``, ``--file_name``) becomes columns of ``output_path``; ``search_subdirectories`` also reads sub-folders (``-s``).
- **API:** `humann_join_tables(input_dir: 'str', output_path: 'str', *, file_name: 'str | None' = None, search_subdirectories: 'bool' = False, humann_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program(s) of HUMAnN (3.9 is the tested version) on PATH or through `humann_dir`; it is not installed in the environment that produced this page. The example below shows behaviour that does not need it: the error raised for a missing program.
- **Input:** required: `input_dir`, `output_path`; optional keyword-only: `file_name`, `search_subdirectories`, `humann_dir` (meanings and defaults are in the Method text and the signature). Input files must exist; output paths must not exist yet.
- **Output:** The common result dictionary. Result keys: ``table`` (``file``, ``columns``, ``row_count`` and the first ``rows``).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.microbiome.humann_p12_cli_functions import humann_join_tables

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'input').mkdir()
    (root / 'bin').mkdir()
    try:
        humann_join_tables(
            input_dir=str(root / 'input'),
            output_path=str(root / 'result.txt'),
            humann_dir=str(root / 'bin'),
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
from dnakit.microbiome.humann_p12_cli_functions import humann_join_tables

result = humann_join_tables(
    input_dir='input',
    output_path='result.txt',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 15) `dnakit.microbiome.humann_quantify_alignments()`

- **Purpose:** Quantify gene families and pathways from alignments or from a gene family table with ``humann``.
- **Method:** This function covers three operations of the work list through ``input_format``: ``sam`` or ``bam`` (read alignments to a UniRef/ChocoPhlAn reference, ``--input-format sam|bam``), ``blastm8`` (a BLAST-like tabular alignment file such as DIAMOND output, ``--input-format blastm8``) and ``genetable`` or ``biom`` (an existing gene family table, from which only pathway abundance and coverage are computed). The XIPE option (``--xipe on``) is not offered: in HUMAnN 3.9 it fails on Python 3 (``xipe.py`` uses the Python 2 name ``file``). The search steps of HUMAnN are not run. ``minpath`` and ``gap_fill`` switch the MinPath and gap-fill steps (both on by default as in HUMAnN; MinPath must be installed for ``minpath=True``); ``pathways`` is ``metacyc`` or ``unipathway``; ``pathways_database`` is one or two custom pathway mapping files (``--pathways-database``); ``id_mapping`` maps alignment reference ids to gene families (``--id-mapping``) and ``annotation_gene_index`` is the column of the gene id in the reference annotation (``--annotation-gene-index``). ``output_basename`` names the outputs (default: the input file name before the first dot); ``output_format`` is ``tsv`` or ``biom``; ``remove_stratified_output`` and ``remove_column_description_output`` simplify the tables. ``output_dir`` must not exist yet.
- **API:** `humann_quantify_alignments(input_path: 'str', output_dir: 'str', *, input_format: 'str' = 'blastm8', minpath: 'bool' = True, gap_fill: 'bool' = True, pathways: 'str' = 'metacyc', pathways_database: 'Sequence[str] | None' = None, id_mapping: 'str | None' = None, annotation_gene_index: 'int | None' = None, output_basename: 'str | None' = None, output_format: 'str' = 'tsv', output_max_decimals: 'int | None' = None, remove_stratified_output: 'bool' = False, remove_column_description_output: 'bool' = False, remove_temp_output: 'bool' = False, nucleotide_database: 'str | None' = None, protein_database: 'str | None' = None, threads: 'int' = 1, timeout_seconds: 'float' = 86400.0, humann_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program(s) of HUMAnN (3.9 is the tested version) on PATH or through `humann_dir`; it is not installed in the environment that produced this page. The example below shows behaviour that does not need it: the error raised for a missing program.
- **Input:** required: `input_path`, `output_dir`; optional keyword-only: `input_format`, `minpath`, `gap_fill`, `pathways`, `pathways_database`, `id_mapping`, `annotation_gene_index`, `output_basename`, `output_format`, `output_max_decimals`, `remove_stratified_output`, `remove_column_description_output`, `remove_temp_output`, `nucleotide_database`, `protein_database`, `threads`, `timeout_seconds`, `humann_dir` (meanings and defaults are in the Method text and the signature). Input files must exist; output paths must not exist yet.
- **Output:** The common result dictionary. Result keys: ``output_dir``, ``gene_families`` (absent for ``genetable``/``biom`` input), ``pathway_abundance`` and ``pathway_coverage`` (each with ``file``, ``columns``, ``row_count`` and the first ``rows``; for ``biom`` output only ``file``).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.microbiome.humann_p12_cli_functions import humann_quantify_alignments

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'in.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        humann_quantify_alignments(
            input_path=str(root / 'in.fa'),
            output_dir=str(root / 'out'),
            humann_dir=str(root / 'bin'),
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
from dnakit.microbiome.humann_p12_cli_functions import humann_quantify_alignments

result = humann_quantify_alignments(
    input_path='in.fa',
    output_dir='out',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 16) `dnakit.microbiome.humann_reduce_table()`

- **Purpose:** Reduce a gene or pathway table to one value per feature with a summary function (``humann_reduce_table``).
- **Method:** ``function`` is ``max`` (the HUMAnN default), ``sum``, ``min`` or ``mean`` and is applied across the sample columns of every feature; ``sort_by`` (``level``, ``name`` or ``value``) sorts the output.
- **API:** `humann_reduce_table(input_path: 'str', output_path: 'str', *, function: 'str' = 'max', sort_by: 'str | None' = None, humann_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program(s) of HUMAnN (3.9 is the tested version) on PATH or through `humann_dir`; it is not installed in the environment that produced this page. The example below shows behaviour that does not need it: the error raised for a missing program.
- **Input:** required: `input_path`, `output_path`; optional keyword-only: `function`, `sort_by`, `humann_dir` (meanings and defaults are in the Method text and the signature). Input files must exist; output paths must not exist yet.
- **Output:** The common result dictionary. Result keys: ``table`` (``file``, ``columns``, ``row_count`` and the first ``rows``).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.microbiome.humann_p12_cli_functions import humann_reduce_table

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'in.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        humann_reduce_table(
            input_path=str(root / 'in.fa'),
            output_path=str(root / 'result.txt'),
            humann_dir=str(root / 'bin'),
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
from dnakit.microbiome.humann_p12_cli_functions import humann_reduce_table

result = humann_reduce_table(
    input_path='in.fa',
    output_path='result.txt',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 17) `dnakit.microbiome.humann_regroup_table()`

- **Purpose:** Regroup gene families into functional groups such as GO, KO, EC, Pfam or MetaCyc reactions (``humann_regroup_table``).
- **Method:** Give exactly one of ``groups`` (a built-in mapping, ``uniref90_rxn`` or ``uniref50_rxn``, which needs the HUMAnN utility mapping data; other mappings such as GO, KO, EC and Pfam are supplied with ``custom_groups``) and ``custom_groups`` (a TSV or gzipped TSV with one group per line followed by its member features; with ``reversed_mapping=True`` one feature per line followed by its groups). ``function`` combines grouped features (``sum`` or ``mean``), ``precision`` rounds to that many decimals, ``include_ungrouped`` adds an UNGROUPED group (``--ungrouped Y``) and ``keep_protected`` carries UNMAPPED through (``--protected Y``).
- **API:** `humann_regroup_table(input_path: 'str', output_path: 'str', *, groups: 'str | None' = None, custom_groups: 'str | None' = None, reversed_mapping: 'bool' = False, function: 'str' = 'sum', precision: 'int | None' = None, include_ungrouped: 'bool' = True, keep_protected: 'bool' = True, humann_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program(s) of HUMAnN (3.9 is the tested version) on PATH or through `humann_dir`; it is not installed in the environment that produced this page. The example below shows behaviour that does not need it: the error raised for a missing program.
- **Input:** required: `input_path`, `output_path`; optional keyword-only: `groups`, `custom_groups`, `reversed_mapping`, `function`, `precision`, `include_ungrouped`, `keep_protected`, `humann_dir` (meanings and defaults are in the Method text and the signature). Input files must exist; output paths must not exist yet.
- **Output:** The common result dictionary. Result keys: ``table`` (``file``, ``columns``, ``row_count`` and the first ``rows``).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.microbiome.humann_p12_cli_functions import humann_regroup_table

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'groups.tsv').write_text('x\n')
    (root / 'in.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        humann_regroup_table(
            input_path=str(root / 'in.fa'),
            output_path=str(root / 'result.txt'),
            custom_groups=str(root / 'groups.tsv'),
            humann_dir=str(root / 'bin'),
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
from dnakit.microbiome.humann_p12_cli_functions import humann_regroup_table

result = humann_regroup_table(
    input_path='in.fa',
    output_path='result.txt',
    custom_groups='groups.tsv',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 18) `dnakit.microbiome.humann_rename_table()`

- **Purpose:** Attach feature names (glosses) to the gene or pathway ids of a table (``humann_rename_table``).
- **Method:** Give exactly one of ``names`` (a built-in mapping: kegg-orthology, kegg-pathway, kegg-module, ec, metacyc-rxn, metacyc-pwy, pfam, eggnog, go, infogo1000 or uniref50; these need the HUMAnN utility mapping data) and ``custom_mapping`` (a two-column TSV, optionally gzipped, from feature id to full name). ``simplify`` removes non-alphanumeric characters from the names.
- **API:** `humann_rename_table(input_path: 'str', output_path: 'str', *, names: 'str | None' = None, custom_mapping: 'str | None' = None, simplify: 'bool' = False, humann_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program(s) of HUMAnN (3.9 is the tested version) on PATH or through `humann_dir`; it is not installed in the environment that produced this page. The example below shows behaviour that does not need it: the error raised for a missing program.
- **Input:** required: `input_path`, `output_path`; optional keyword-only: `names`, `custom_mapping`, `simplify`, `humann_dir` (meanings and defaults are in the Method text and the signature). Input files must exist; output paths must not exist yet.
- **Output:** The common result dictionary. Result keys: ``table`` (``file``, ``columns``, ``row_count`` and the first ``rows``).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.microbiome.humann_p12_cli_functions import humann_rename_table

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'in.fa').write_text('x\n')
    (root / 'names.tsv').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        humann_rename_table(
            input_path=str(root / 'in.fa'),
            output_path=str(root / 'result.txt'),
            custom_mapping=str(root / 'names.tsv'),
            humann_dir=str(root / 'bin'),
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
from dnakit.microbiome.humann_p12_cli_functions import humann_rename_table

result = humann_rename_table(
    input_path='in.fa',
    output_path='result.txt',
    custom_mapping='names.tsv',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 19) `dnakit.microbiome.humann_renorm_table()`

- **Purpose:** Renormalise gene or pathway abundances to copies per million or relative abundance (``humann_renorm_table``).
- **Method:** ``units`` is ``cpm`` (copies per million, the HUMAnN default) or ``relab`` (relative abundance); ``mode`` normalises all levels of a stratified table by the ``community`` total (default) or by their own totals (``levelwise``); ``include_special`` keeps the special features UNMAPPED, UNINTEGRATED and UNGROUPED (``--special y``, the default); ``update_sample_names`` replaces ``-RPK`` in the sample names by the new suffix (``--update-snames``).
- **API:** `humann_renorm_table(input_path: 'str', output_path: 'str', *, units: 'str' = 'cpm', mode: 'str' = 'community', include_special: 'bool' = True, update_sample_names: 'bool' = False, humann_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program(s) of HUMAnN (3.9 is the tested version) on PATH or through `humann_dir`; it is not installed in the environment that produced this page. The example below shows behaviour that does not need it: the error raised for a missing program.
- **Input:** required: `input_path`, `output_path`; optional keyword-only: `units`, `mode`, `include_special`, `update_sample_names`, `humann_dir` (meanings and defaults are in the Method text and the signature). Input files must exist; output paths must not exist yet.
- **Output:** The common result dictionary. Result keys: ``table`` (``file``, ``columns``, ``row_count`` and the first ``rows``).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.microbiome.humann_p12_cli_functions import humann_renorm_table

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'in.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        humann_renorm_table(
            input_path=str(root / 'in.fa'),
            output_path=str(root / 'result.txt'),
            humann_dir=str(root / 'bin'),
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
from dnakit.microbiome.humann_p12_cli_functions import humann_renorm_table

result = humann_renorm_table(
    input_path='in.fa',
    output_path='result.txt',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 20) `dnakit.microbiome.humann_split_stratified_table()`

- **Purpose:** Split a stratified table into a stratified table and an unstratified (community totals) table (``humann_split_stratified_table``).
- **Method:** ``input_path`` is a TSV (plain, gzip or bzip2) or BIOM HUMAnN table; ``output_dir`` must not exist and receives ``<name>_stratified.tsv`` and ``<name>_unstratified.tsv``.
- **API:** `humann_split_stratified_table(input_path: 'str', output_dir: 'str', *, humann_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program(s) of HUMAnN (3.9 is the tested version) on PATH or through `humann_dir`; it is not installed in the environment that produced this page. The example below shows behaviour that does not need it: the error raised for a missing program.
- **Input:** required: `input_path`, `output_dir`; optional keyword-only: `humann_dir` (meanings and defaults are in the Method text and the signature). Input files must exist; output paths must not exist yet.
- **Output:** The common result dictionary. Result keys: ``stratified`` and ``unstratified`` (each with ``file``, ``columns``, ``row_count`` and the first ``rows``).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.microbiome.humann_p12_cli_functions import humann_split_stratified_table

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'in.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        humann_split_stratified_table(
            input_path=str(root / 'in.fa'),
            output_dir=str(root / 'out'),
            humann_dir=str(root / 'bin'),
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
from dnakit.microbiome.humann_p12_cli_functions import humann_split_stratified_table

result = humann_split_stratified_table(
    input_path='in.fa',
    output_dir='out',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 21) `dnakit.microbiome.humann_split_table()`

- **Purpose:** Split a merged gene or pathway table into one table per sample (``humann_split_table``).
- **Method:** ``output_dir`` must not exist; one table per sample column is written there. ``taxonomy_index`` is the index of the gene in the taxonomy data and ``taxonomy_level`` (Kingdom, Phylum, Class, Order, Family, Genus or Species) the level of taxonomy for the output when the input comes from PICRUSt ``metagenome_contributions.py``.
- **API:** `humann_split_table(input_path: 'str', output_dir: 'str', *, taxonomy_index: 'int | None' = None, taxonomy_level: 'str | None' = None, humann_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program(s) of HUMAnN (3.9 is the tested version) on PATH or through `humann_dir`; it is not installed in the environment that produced this page. The example below shows behaviour that does not need it: the error raised for a missing program.
- **Input:** required: `input_path`, `output_dir`; optional keyword-only: `taxonomy_index`, `taxonomy_level`, `humann_dir` (meanings and defaults are in the Method text and the signature). Input files must exist; output paths must not exist yet.
- **Output:** The common result dictionary. Result keys: ``tables`` (the files written, sorted) and ``table_count``.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.microbiome.humann_p12_cli_functions import humann_split_table

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'in.fa').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        humann_split_table(
            input_path=str(root / 'in.fa'),
            output_dir=str(root / 'out'),
            humann_dir=str(root / 'bin'),
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
from dnakit.microbiome.humann_p12_cli_functions import humann_split_table

result = humann_split_table(
    input_path='in.fa',
    output_dir='out',
    verbose=False,
)
print(result['status'], result['outputs'])
```

## 22) `dnakit.microbiome.humann_unpack_pathways()`

- **Purpose:** Unpack pathway abundances into the abundances of the genes that make up each pathway (``humann_unpack_pathways``).
- **Method:** ``gene_table`` is the gene family (or EC) abundance table and ``pathway_table`` the pathway abundance table of the same samples. ``gene_mapping`` (gene family to reaction) and ``pathway_mapping`` (reaction to pathway) default to the MetaCyc files shipped with HUMAnN. ``remove_taxonomy`` drops the taxonomy strata (``--remove-taxonomy``).
- **API:** `humann_unpack_pathways(gene_table: 'str', pathway_table: 'str', output_path: 'str', *, gene_mapping: 'str | None' = None, pathway_mapping: 'str | None' = None, remove_taxonomy: 'bool' = False, humann_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external program(s) of HUMAnN (3.9 is the tested version) on PATH or through `humann_dir`; it is not installed in the environment that produced this page. The example below shows behaviour that does not need it: the error raised for a missing program.
- **Input:** required: `gene_table`, `pathway_table`, `output_path`; optional keyword-only: `gene_mapping`, `pathway_mapping`, `remove_taxonomy`, `humann_dir` (meanings and defaults are in the Method text and the signature). Input files must exist; output paths must not exist yet.
- **Output:** The common result dictionary. Result keys: ``table`` (``file``, ``columns``, ``row_count`` and the first ``rows``).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import DNAKitError
from dnakit.microbiome.humann_p12_cli_functions import humann_unpack_pathways

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'genes.tsv').write_text('x\n')
    (root / 'paths.tsv').write_text('x\n')
    (root / 'bin').mkdir()
    try:
        humann_unpack_pathways(
            gene_table=str(root / 'genes.tsv'),
            pathway_table=str(root / 'paths.tsv'),
            output_path=str(root / 'result.txt'),
            humann_dir=str(root / 'bin'),
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
from dnakit.microbiome.humann_p12_cli_functions import humann_unpack_pathways

result = humann_unpack_pathways(
    gene_table='genes.tsv',
    pathway_table='paths.tsv',
    output_path='result.txt',
    verbose=False,
)
print(result['status'], result['outputs'])
```
