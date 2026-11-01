# GATK reference, FASTA and interval tools (GATK, pack P01)

This group of 15 functions in `dnakit.sequence.gatk_p01_reference_cli_functions` runs tools of the Broad Institute's [GATK](https://github.com/broadinstitute/gatk) 4.7.0.0 command line (Picard tools bundled in GATK included). Each function runs one fixed tool, or a family of tools chosen by a mode argument, without a shell and with a timeout and an output size limit. Only the long options that `gatk <Tool> --help` lists are accepted through `options`; options that a parameter already sets cannot be repeated, and GATK's script-executing `JAVASCRIPT_FILE` is not reachable. Outputs are written to explicit new paths and never overwritten.

GATK is an external program: install it (for example `conda install -c bioconda gatk4` or the release zip) with Java 17, then put the `gatk` launcher on PATH, set `DNAKIT_GATK` or pass `gatk_path`. If it cannot be found, `BackendUnavailableError` (code `GATK_UNAVAILABLE`) is raised with the install hint; a failing tool raises `BackendExecutionError` (code `GATK_EXECUTION_FAILED`) with the tail of its log; a tool that exits 0 without writing a requested output raises `GATK_OUTPUT_MISSING`. DNAKit does not bundle GATK. Spark tools run in local mode (`spark-master` can be set through `options`). Reference FASTA files need `.fai` and `.dict` files next to them, BAM/VCF inputs their indexes, as for GATK itself. Results are dictionaries; the section for each function lists the extra fields.

## 1) `dnakit.sequence.gatk_p01_reference_cli_functions.gatk_index_feature_file()`

- **Purpose:** Create a tribble index (.idx, or .tbi for bgzipped VCF) for a feature file.
- **Method:** Run `gatk IndexFeatureFile`; without `index_path` the index file next to the input is detected after the run and reported as `outputs.index`.
- **API:** `gatk_index_feature_file(feature_path: 'str', *, index_path: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): IndexFeatureFile.
- **Input:** `feature_path` (required): Feature file (VCF, BED, GFF ...) in a tribble-supported format.; `index_path` (optional): Explicit index file; by default the index is written next to the feature file.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs.index` is the index. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.sequence.gatk_p01_reference_cli_functions import gatk_index_feature_file

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    feature_path = folder / "feature_path.dat"; feature_path.write_bytes(b"x")
    try:
        gatk_index_feature_file(
            feature_path=str(feature_path),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 2) `dnakit.sequence.gatk_p01_reference_cli_functions.gatk_sort_gff()`

- **Purpose:** Sort a GFF3 file and add the flush (`###`) directives.
- **Method:** Run Picard `gatk SortGff`; contigs are ordered by `sequence_dictionary` when given.
- **API:** `gatk_sort_gff(input_path: 'str', output_path: 'str', *, sequence_dictionary: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): SortGff.
- **Input:** `input_path` (required): GFF3 file.; `output_path` (required): New sorted GFF3 file (with flush directives).; `sequence_dictionary` (optional): Sequence dictionary that defines the contig order.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs.output_path` is the sorted GFF3. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.sequence.gatk_p01_reference_cli_functions import gatk_sort_gff

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    try:
        gatk_sort_gff(
            input_path=str(input_path),
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 3) `dnakit.sequence.gatk_p01_reference_cli_functions.gatk_gtf_to_bed()`

- **Purpose:** Convert a Gencode GTF file to a BED file.
- **Method:** Run `gatk GtfToBed`.
- **API:** `gatk_gtf_to_bed(gtf_path: 'str', output_path: 'str', sequence_dictionary: 'str', *, sort_by_transcript: 'bool | None' = None, use_basic_transcript: 'bool | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): GtfToBed.
- **Input:** `gtf_path` (required): Gencode GTF file.; `output_path` (required): New BED file.; `sequence_dictionary` (required): Sequence dictionary (.dict) of the reference, used to order contigs.; `sort_by_transcript` (optional): Sort the output by transcript.; `use_basic_transcript` (optional): Keep only transcripts tagged `basic`.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs.output_path` is the BED file. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.sequence.gatk_p01_reference_cli_functions import gatk_gtf_to_bed

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    gtf_path = folder / "gtf_path.dat"; gtf_path.write_bytes(b"x")
    sequence_dictionary = folder / "sequence_dictionary.dat"; sequence_dictionary.write_bytes(b"x")
    try:
        gatk_gtf_to_bed(
            gtf_path=str(gtf_path),
            output_path=str(folder / "output_path.out"),
            sequence_dictionary=str(sequence_dictionary),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 4) `dnakit.sequence.gatk_p01_reference_cli_functions.gatk_bait_designer()`

- **Purpose:** Design oligonucleotide baits for hybrid selection (exome capture) from target intervals.
- **Method:** Run Picard `gatk BaitDesigner`; the files written to `output_dir` are listed in `files`.
- **API:** `gatk_bait_designer(targets: 'str', reference_path: 'str', design_name: 'str', output_dir: 'str', *, bait_size: 'int | None' = None, bait_offset: 'int | None' = None, design_strategy: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): BaitDesigner.
- **Input:** `targets` (required): Interval list of the target regions.; `reference_path` (required): Reference FASTA (with .fai and .dict).; `design_name` (required): Name of the bait design.; `output_dir` (required): New directory for the design files.; `bait_size` (optional): Length of each bait.; `bait_offset` (optional): Offset between consecutive baits.; `design_strategy` (optional): How baits are placed over a target.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields and `files` (file names inside `output_dir`). Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.sequence.gatk_p01_reference_cli_functions import gatk_bait_designer

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    targets = folder / "targets.dat"; targets.write_bytes(b"x")
    reference_path = folder / "reference_path.dat"; reference_path.write_bytes(b"x")
    try:
        gatk_bait_designer(
            targets=str(targets),
            reference_path=str(reference_path),
            design_name="design1",
            output_dir=str(folder / "output_dir"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 5) `dnakit.sequence.gatk_p01_reference_cli_functions.gatk_bwa_mem_index_image()`

- **Purpose:** Create a BWA-MEM index image (`.img`) of a reference for the GATK BWA tools.
- **Method:** Run `gatk BwaMemIndexImageCreator`; the image `<fasta>.img` is written next to the FASTA and must not exist yet; it is reported as `outputs.image`.
- **API:** `gatk_bwa_mem_index_image(reference_path: 'str', *, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): BwaMemIndexImageCreator.
- **Input:** `reference_path` (required): Reference FASTA.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs.image` is the index image. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.sequence.gatk_p01_reference_cli_functions import gatk_bwa_mem_index_image

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    reference_path = folder / "reference_path.dat"; reference_path.write_bytes(b"x")
    try:
        gatk_bwa_mem_index_image(
            reference_path=str(reference_path),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 6) `dnakit.sequence.gatk_p01_reference_cli_functions.gatk_check_reference_compatibility()`

- **Purpose:** Check a BAM/CRAM or VCF for compatibility with one or more reference FASTA files.
- **Method:** Run `gatk CheckReferenceCompatibility` (experimental); at least one of `alignment_path` and `variant_path` is needed.
- **API:** `gatk_check_reference_compatibility(references_to_compare: 'Sequence[str]', *, alignment_path: 'str | None' = None, variant_path: 'str | None' = None, output_path: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): CheckReferenceCompatibility.
- **Input:** `references_to_compare` (required): Reference FASTA files to test against.; `alignment_path` (optional): BAM/CRAM whose sequence dictionary is checked.; `variant_path` (optional): VCF whose contig header is checked.; `output_path` (optional): New report file (default: standard output).; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields (the table is in `log_tail` unless `output_path` is given). Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.sequence.gatk_p01_reference_cli_functions import gatk_check_reference_compatibility

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    references_to_compare = folder / "references_to_compare.dat"; references_to_compare.write_bytes(b"x")
    alignment_path = folder / "alignment_path.dat"; alignment_path.write_bytes(b"x")
    try:
        gatk_check_reference_compatibility(
            references_to_compare=[str(references_to_compare)],
            alignment_path=str(alignment_path),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 7) `dnakit.sequence.gatk_p01_reference_cli_functions.gatk_compare_references()`

- **Purpose:** Compare reference FASTA files contig by contig (MD5, length, base differences) and summarise the differences.
- **Method:** Run `gatk CompareReferences` (experimental) and parse the table it prints.
- **API:** `gatk_compare_references(reference_path: 'str', references_to_compare: 'Sequence[str]', *, display_only_differing_sequences: 'bool | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): CompareReferences.
- **Input:** `reference_path` (required): Main reference FASTA.; `references_to_compare` (required): Reference FASTA files to compare with.; `display_only_differing_sequences` (optional): List only sequences that differ.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields plus `comparison_table` (rows of the printed table). Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Not run to completion here: the arguments and option whitelist are checked against the 4.7.0.0 help, but the tool needs inputs (GATK-SV evidence, Funcotator data sources, trios or similar) that were not available as toy data. (arguments checked; upstream raised an error when a reference is compared with itself, and the two-reference run was not completed)
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.sequence.gatk_p01_reference_cli_functions import gatk_compare_references

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    reference_path = folder / "reference_path.dat"; reference_path.write_bytes(b"x")
    references_to_compare = folder / "references_to_compare.dat"; references_to_compare.write_bytes(b"x")
    try:
        gatk_compare_references(
            reference_path=str(reference_path),
            references_to_compare=[str(references_to_compare)],
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 8) `dnakit.sequence.gatk_p01_reference_cli_functions.gatk_compose_str_table()`

- **Purpose:** Compose a genome-wide short-tandem-repeat (STR) location table for DRAGstr model calibration.
- **Method:** Run `gatk ComposeSTRTableFile`.
- **API:** `gatk_compose_str_table(reference_path: 'str', output_path: 'str', *, max_period: 'int | None' = None, max_repeats: 'int | None' = None, decimation: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): ComposeSTRTableFile.
- **Input:** `reference_path` (required): Reference FASTA (with .fai and .dict).; `output_path` (required): New STR table (binary).; `max_period` (optional): Maximum repeat unit length.; `max_repeats` (optional): Maximum repeat count.; `decimation` (optional): Decimation table (`period:factor,...`).; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs.output_path` is the STR table. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.sequence.gatk_p01_reference_cli_functions import gatk_compose_str_table

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    reference_path = folder / "reference_path.dat"; reference_path.write_bytes(b"x")
    try:
        gatk_compose_str_table(
            reference_path=str(reference_path),
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 9) `dnakit.sequence.gatk_p01_reference_cli_functions.gatk_count_reference_bases()`

- **Purpose:** Count the bases of a FASTA file: the number of A, C, G, T and N, or the number of non-N bases.
- **Method:** Run `gatk CountBasesInReference` (`mode=count_bases`, parsed into `base_counts`) or Picard `NonNFastaSize` (`mode=non_n_size`, result in the output file, optionally limited to an interval list). `mode`: `count_bases` = CountBasesInReference, `non_n_size` = NonNFastaSize.
- **API:** `gatk_count_reference_bases(fasta_path: 'str', *, mode: 'str' = 'count_bases', output_path: 'str | None' = None, intervals: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): CountBasesInReference, NonNFastaSize.
- **Input:** `fasta_path` (required): FASTA file.; `output_path` (required in some modes): New file for the non-N size.; `intervals` (optional): Interval list limiting the count.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields and `mode`; `base_counts` and `total_bases` for `count_bases`, `outputs.output_path` for `non_n_size`. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.sequence.gatk_p01_reference_cli_functions import gatk_count_reference_bases

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    fasta_path = folder / "fasta_path.dat"; fasta_path.write_bytes(b"x")
    try:
        gatk_count_reference_bases(
            fasta_path=str(fasta_path),
            mode='count_bases',
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 10) `dnakit.sequence.gatk_p01_reference_cli_functions.gatk_extract_reference_sequences()`

- **Purpose:** Extract the sequence of intervals from a reference into a new FASTA file.
- **Method:** Run Picard `ExtractSequences` (`source=interval_list`, intervals from an interval-list file) or `gatk FastaReferenceMaker` (`source=intervals`, intervals written as `chr:start-end`). `source`: `interval_list` = ExtractSequences, `intervals` = FastaReferenceMaker.
- **API:** `gatk_extract_reference_sequences(reference_path: 'str', output_path: 'str', *, source: 'str' = 'interval_list', interval_list: 'str | None' = None, intervals: 'Sequence[str] | None' = None, line_length: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): ExtractSequences, FastaReferenceMaker.
- **Input:** `reference_path` (required): Reference FASTA (with .fai and .dict).; `output_path` (required): New FASTA file.; `interval_list` (required in some modes): Interval list to extract (Picard format).; `intervals` (required in some modes): Intervals to extract, for example `chr1:100-200`.; `line_length` (optional): Line width of the output FASTA.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields and `mode`; `outputs.output_path` is the FASTA. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.sequence.gatk_p01_reference_cli_functions import gatk_extract_reference_sequences

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    reference_path = folder / "reference_path.dat"; reference_path.write_bytes(b"x")
    interval_list = folder / "interval_list.dat"; interval_list.write_bytes(b"x")
    try:
        gatk_extract_reference_sequences(
            reference_path=str(reference_path),
            output_path=str(folder / "output_path.out"),
            interval_list=str(interval_list),
            source='interval_list',
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 11) `dnakit.sequence.gatk_p01_reference_cli_functions.gatk_fasta_alternate_reference_maker()`

- **Purpose:** Create an alternative reference by applying the variants of a VCF to a reference FASTA.
- **Method:** Run `gatk FastaAlternateReferenceMaker`; restrict the regions with the `intervals` option.
- **API:** `gatk_fasta_alternate_reference_maker(reference_path: 'str', variant_path: 'str', output_path: 'str', *, snp_mask: 'str | None' = None, use_iupac_sample: 'str | None' = None, line_width: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): FastaAlternateReferenceMaker.
- **Input:** `reference_path` (required): Reference FASTA (with .fai and .dict).; `variant_path` (required): Indexed VCF whose alleles are applied.; `output_path` (required): New alternate-reference FASTA.; `snp_mask` (optional): VCF of positions to mask with N.; `use_iupac_sample` (optional): Sample whose heterozygous genotypes are written as IUPAC codes.; `line_width` (optional): Line width of the output FASTA.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs.output_path` is the FASTA. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.sequence.gatk_p01_reference_cli_functions import gatk_fasta_alternate_reference_maker

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    reference_path = folder / "reference_path.dat"; reference_path.write_bytes(b"x")
    variant_path = folder / "variant_path.dat"; variant_path.write_bytes(b"x")
    try:
        gatk_fasta_alternate_reference_maker(
            reference_path=str(reference_path),
            variant_path=str(variant_path),
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 12) `dnakit.sequence.gatk_p01_reference_cli_functions.gatk_find_bad_genomic_kmers_spark()`

- **Purpose:** Find the k-mers that occur at high frequency in a reference (for example to mask them in SV discovery).
- **Method:** Run `gatk FindBadGenomicKmersSpark` (beta, Spark local mode).
- **API:** `gatk_find_bad_genomic_kmers_spark(reference_path: 'str', output_path: 'str', *, k_size: 'int | None' = None, kmer_max_dust_score: 'int | None' = None, high_copy_fasta: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): FindBadGenomicKmersSpark.
- **Input:** `reference_path` (required): Reference FASTA (a .2bit file also works).; `output_path` (required): New text file of over-represented k-mers.; `k_size` (optional): Length of the k-mers.; `kmer_max_dust_score` (optional): Maximum DUST score of a reported k-mer.; `high_copy_fasta` (optional): FASTA of known high-copy sequences to include.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs.output_path` is the k-mer list. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Not run to completion here: the arguments and option whitelist are checked against the 4.7.0.0 help, but the tool needs inputs (GATK-SV evidence, Funcotator data sources, trios or similar) that were not available as toy data. (the Spark job did not finish on the toy genome and was stopped)
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.sequence.gatk_p01_reference_cli_functions import gatk_find_bad_genomic_kmers_spark

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    reference_path = folder / "reference_path.dat"; reference_path.write_bytes(b"x")
    try:
        gatk_find_bad_genomic_kmers_spark(
            reference_path=str(reference_path),
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 13) `dnakit.sequence.gatk_p01_reference_cli_functions.gatk_normalize_fasta()`

- **Purpose:** Rewrite a FASTA file so that all sequence lines have the same length.
- **Method:** Run Picard `gatk NormalizeFasta`.
- **API:** `gatk_normalize_fasta(input_path: 'str', output_path: 'str', *, line_length: 'int | None' = None, truncate_sequence_names_at_whitespace: 'bool | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): NormalizeFasta.
- **Input:** `input_path` (required): FASTA file.; `output_path` (required): New normalised FASTA.; `line_length` (optional): Sequence line length (default 100).; `truncate_sequence_names_at_whitespace` (optional): Cut sequence names at the first whitespace.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs.output_path` is the new FASTA. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.sequence.gatk_p01_reference_cli_functions import gatk_normalize_fasta

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    try:
        gatk_normalize_fasta(
            input_path=str(input_path),
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 14) `dnakit.sequence.gatk_p01_reference_cli_functions.gatk_scatter_intervals_by_ns()`

- **Purpose:** Write an interval list made by splitting a reference at its runs of N.
- **Method:** Run Picard `gatk ScatterIntervalsByNs`.
- **API:** `gatk_scatter_intervals_by_ns(reference_path: 'str', output_path: 'str', *, output_type: 'str | None' = None, max_to_merge: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): ScatterIntervalsByNs.
- **Input:** `reference_path` (required): Reference FASTA (with .fai and .dict).; `output_path` (required): New interval list.; `output_type` (optional): Which intervals to write.; `max_to_merge` (optional): Merge N-runs of at most this length into the neighbouring ACGT interval.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs.output_path` is the interval list. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.sequence.gatk_p01_reference_cli_functions import gatk_scatter_intervals_by_ns

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    reference_path = folder / "reference_path.dat"; reference_path.write_bytes(b"x")
    try:
        gatk_scatter_intervals_by_ns(
            reference_path=str(reference_path),
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 15) `dnakit.sequence.gatk_p01_reference_cli_functions.gatk_shift_fasta()`

- **Purpose:** Create a shifted FASTA (circular contigs such as chrM rotated by an offset) and the file that shifts coordinates back.
- **Method:** Run `gatk ShiftFasta` (beta).
- **API:** `gatk_shift_fasta(reference_path: 'str', output_path: 'str', shift_back_output: 'str', *, shift_offset_list: 'Sequence[str] | None' = None, interval_file_name: 'str | None' = None, line_width: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): ShiftFasta.
- **Input:** `reference_path` (required): Reference FASTA (with .fai and .dict).; `output_path` (required): New shifted FASTA.; `shift_back_output` (required): New chain file that maps shifted coordinates back.; `shift_offset_list` (optional): Offset per contig, in the order of the dictionary.; `interval_file_name` (optional): New interval file for the shifted contigs.; `line_width` (optional): Line width of the output FASTA.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs` lists the FASTA and the shift-back file. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.sequence.gatk_p01_reference_cli_functions import gatk_shift_fasta

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    reference_path = folder / "reference_path.dat"; reference_path.write_bytes(b"x")
    try:
        gatk_shift_fasta(
            reference_path=str(reference_path),
            output_path=str(folder / "output_path.out"),
            shift_back_output=str(folder / "shift_back_output.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **Example output:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```
