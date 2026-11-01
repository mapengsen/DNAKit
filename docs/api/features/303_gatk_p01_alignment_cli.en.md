# GATK alignment and read-file tools (GATK, pack P01)

This group of 35 functions in `dnakit.alignment.gatk_p01_alignment_cli_functions` runs tools of the Broad Institute's [GATK](https://github.com/broadinstitute/gatk) 4.7.0.0 command line (Picard tools bundled in GATK included). Each function runs one fixed tool, or a family of tools chosen by a mode argument, without a shell and with a timeout and an output size limit. Only the long options that `gatk <Tool> --help` lists are accepted through `options`; options that a parameter already sets cannot be repeated, and GATK's script-executing `JAVASCRIPT_FILE` is not reachable. Outputs are written to explicit new paths and never overwritten.

GATK is an external program: install it (for example `conda install -c bioconda gatk4` or the release zip) with Java 17, then put the `gatk` launcher on PATH, set `DNAKIT_GATK` or pass `gatk_path`. If it cannot be found, `BackendUnavailableError` (code `GATK_UNAVAILABLE`) is raised with the install hint; a failing tool raises `BackendExecutionError` (code `GATK_EXECUTION_FAILED`) with the tail of its log; a tool that exits 0 without writing a requested output raises `GATK_OUTPUT_MISSING`. DNAKit does not bundle GATK. Spark tools run in local mode (`spark-master` can be set through `options`). Reference FASTA files need `.fai` and `.dict` files next to them, BAM/VCF inputs their indexes, as for GATK itself. Results are dictionaries; the section for each function lists the extra fields.

## 1) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_pathseq()`

- **Purpose:** Run the PathSeq pathogen-discovery pipeline: host filtering, microbe alignment, abundance scoring, or all steps at once.
- **Method:** Run `gatk PathSeqFilterSpark`, `PathSeqBwaSpark`, `PathSeqScoreSpark` or `PathSeqPipelineSpark` (Spark local mode) depending on `step`. Host and microbe resources (k-mer file, host BWA image, hash files) are passed as GATK options such as `filter-bwa-image` or `kmer-file`. `step`: `filter` = host-read filtering (PathSeqFilterSpark), `align` = microbe alignment (PathSeqBwaSpark), `score` = abundance scoring (PathSeqScoreSpark), `pipeline` = all steps (PathSeqPipelineSpark).
- **API:** `gatk_pathseq(*, step: 'str' = 'filter', input_path: 'str | None' = None, paired_input: 'str | None' = None, unpaired_input: 'str | None' = None, microbe_bwa_image: 'str | None' = None, microbe_dict: 'str | None' = None, taxonomy_file: 'str | None' = None, paired_output: 'str | None' = None, unpaired_output: 'str | None' = None, scores_output: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): PathSeqBwaSpark, PathSeqFilterSpark, PathSeqPipelineSpark, PathSeqScoreSpark.
- **Input:** `input_path` (required in some modes): Query BAM/CRAM/SAM of reads (host-filtering input).; `paired_input` (optional): Paired reads BAM from the filter step.; `unpaired_input` (optional): Unpaired reads BAM from the filter step.; `microbe_bwa_image` (required in some modes): BWA-MEM index image of the microbe reference (see `gatk_bwa_mem_index_image`).; `microbe_dict` (required in some modes): Sequence dictionary of the microbe reference.; `taxonomy_file` (required in some modes): Taxonomy database built by PathSeqBuildReferenceTaxonomy.; `paired_output` (optional): Output BAM for paired reads.; `unpaired_output` (optional): Output BAM for unpaired reads.; `scores_output` (required in some modes): Pathogen abundance score table (tab-separated).; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields, `mode` and the output paths that were requested. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_pathseq

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    try:
        gatk_pathseq(
            input_path=str(input_path),
            paired_output=str(folder / "paired_output.out"),
            step='filter',
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

## 2) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_detect_cram_issue_8768()`

- **Purpose:** Check a CRAM file for base corruption caused by GATK issue 8768.
- **Method:** Run `gatk CRAMIssue8768Detector`; the mismatch-rate threshold and TSV output are GATK options (`mismatch-rate-threshold`, `output-tsv`).
- **API:** `gatk_detect_cram_issue_8768(cram_path: 'str', reference_path: 'str', output_path: 'str', *, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): CRAMIssue8768Detector.
- **Input:** `cram_path` (required): CRAM file to check.; `reference_path` (required): Reference FASTA used for the CRAM (with .fai and .dict).; `output_path` (required): New report file.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs.output_path` is the report. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_detect_cram_issue_8768

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    cram_path = folder / "cram_path.dat"; cram_path.write_bytes(b"x")
    reference_path = folder / "reference_path.dat"; reference_path.write_bytes(b"x")
    try:
        gatk_detect_cram_issue_8768(
            cram_path=str(cram_path),
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

## 3) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_create_bam_splitting_index()`

- **Purpose:** Create a Hadoop-BAM splitting index (.sbi) for a BAM file.
- **Method:** Run `gatk CreateHadoopBamSplittingIndex`; the index `<bam>.sbi` is written next to the BAM and must not exist yet.
- **API:** `gatk_create_bam_splitting_index(bam_path: 'str', *, create_bai: 'bool | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): CreateHadoopBamSplittingIndex.
- **Input:** `bam_path` (required): Coordinate-sorted BAM file.; `create_bai` (optional): Also write a .bai index.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs.sbi` is the splitting index path. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_create_bam_splitting_index

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    bam_path = folder / "bam_path.dat"; bam_path.write_bytes(b"x")
    try:
        gatk_create_bam_splitting_index(
            bam_path=str(bam_path),
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

## 4) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_print_bgzf_block_information()`

- **Purpose:** Report the compressed and uncompressed size and offset of every BGZF block of a file.
- **Method:** Run `gatk PrintBGZFBlockInformation` and parse its `Block #n at file offset ...` listing (the first 10,000 blocks are returned; `block_count` is the total).
- **API:** `gatk_print_bgzf_block_information(bgzf_path: 'str', *, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): PrintBGZFBlockInformation.
- **Input:** `bgzf_path` (required): BGZF-compressed file (BAM, bgzipped VCF ...).; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields plus `block_count` and `blocks` (block, offset, compressed_size, uncompressed_size). Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_print_bgzf_block_information

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    bgzf_path = folder / "bgzf_path.dat"; bgzf_path.write_bytes(b"x")
    try:
        gatk_print_bgzf_block_information(
            bgzf_path=str(bgzf_path),
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

## 5) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_print_file_diagnostics()`

- **Purpose:** Write diagnostic information (record counts, header and format details) about a genomics file.
- **Method:** Run `gatk PrintFileDiagnostics`; `count-limit` limits the number of records examined.
- **API:** `gatk_print_file_diagnostics(input_path: 'str', output_path: 'str', *, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): PrintFileDiagnostics.
- **Input:** `input_path` (required): Genomics file in a format the tool analyses (CRAM; GATK 4.7.0.0 rejects BAM).; `output_path` (required): New local diagnostics file.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs.output_path` is the diagnostics file. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_print_file_diagnostics

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    try:
        gatk_print_file_diagnostics(
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

## 6) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_read_anonymizer()`

- **Purpose:** Replace the bases of reads with the reference bases they align to, so that no personal sequence remains.
- **Method:** Run `gatk ReadAnonymizer`.
- **API:** `gatk_read_anonymizer(input_path: 'str', reference_path: 'str', output_path: 'str', *, ref_base_quality: 'int | None' = None, use_simple_cigar: 'bool | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): ReadAnonymizer.
- **Input:** `input_path` (required): Input BAM/CRAM/SAM.; `reference_path` (required): Reference FASTA (with .fai and .dict).; `output_path` (required): New output alignment file.; `ref_base_quality` (optional): Base quality assigned to the replacing reference bases.; `use_simple_cigar` (optional): Rewrite the CIGAR as a simple full-length match.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs.output_path` is the anonymised file. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_read_anonymizer

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    reference_path = folder / "reference_path.dat"; reference_path.write_bytes(b"x")
    try:
        gatk_read_anonymizer(
            input_path=str(input_path),
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

## 7) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_add_comments_to_bam()`

- **Purpose:** Add comment (@CO) lines to the header of a SAM/BAM/CRAM file.
- **Method:** Run Picard-based `gatk AddCommentsToBam`; one `--COMMENT` per entry of `comments`.
- **API:** `gatk_add_comments_to_bam(input_path: 'str', output_path: 'str', comments: 'Sequence[str]', *, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): AddCommentsToBam.
- **Input:** `input_path` (required): Input SAM/BAM/CRAM.; `output_path` (required): New output file.; `comments` (required): Comment texts; each becomes a @CO header line.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs.output_path` is the new file. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_add_comments_to_bam

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    try:
        gatk_add_comments_to_bam(
            input_path=str(input_path),
            output_path=str(folder / "output_path.out"),
            comments=["processed by dnakit"],
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

## 8) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_add_original_alignment_tags()`

- **Purpose:** Record the current alignment of every read in the OA tag (and the original mate contig).
- **Method:** Run Picard `AddOATag` (`mode=oa_tag`) or GATK `AddOriginalAlignmentTags` (`mode=original_alignment`, which also adds the original mate contig tag). `mode`: `oa_tag` = AddOATag, `original_alignment` = AddOriginalAlignmentTags.
- **API:** `gatk_add_original_alignment_tags(input_path: 'str', output_path: 'str', *, mode: 'str' = 'oa_tag', interval_list: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): AddOATag, AddOriginalAlignmentTags.
- **Input:** `input_path` (required): Input SAM/BAM/CRAM.; `output_path` (required): New output file.; `interval_list` (optional): Only reads overlapping these intervals receive the OA tag.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields and `mode`; `outputs.output_path` is the new file. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_add_original_alignment_tags

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    try:
        gatk_add_original_alignment_tags(
            input_path=str(input_path),
            output_path=str(folder / "output_path.out"),
            mode='oa_tag',
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

## 9) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_base_recalibrator()`

- **Purpose:** Build the base quality score recalibration (BQSR) table.
- **Method:** Run `gatk BaseRecalibrator`; covariates and quality settings are GATK options such as `low-quality-tail` or `maximum-cycle-value`.
- **API:** `gatk_base_recalibrator(input_path: 'str', reference_path: 'str', known_sites: 'Sequence[str]', output_path: 'str', *, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): BaseRecalibrator.
- **Input:** `input_path` (required): Input BAM/CRAM.; `reference_path` (required): Reference FASTA (with .fai and .dict).; `known_sites` (required): Known-variant VCF files (indexed) to mask.; `output_path` (required): New recalibration table.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs.output_path` is the recalibration table. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_base_recalibrator

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    reference_path = folder / "reference_path.dat"; reference_path.write_bytes(b"x")
    known_sites = folder / "known_sites.dat"; known_sites.write_bytes(b"x")
    try:
        gatk_base_recalibrator(
            input_path=str(input_path),
            reference_path=str(reference_path),
            known_sites=[str(known_sites)],
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

## 10) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_apply_bqsr()`

- **Purpose:** Apply a base quality score recalibration table to the reads.
- **Method:** Run `gatk ApplyBQSR`.
- **API:** `gatk_apply_bqsr(input_path: 'str', recalibration_table: 'str', output_path: 'str', *, reference_path: 'str | None' = None, static_quantized_quals: 'Sequence[str] | None' = None, use_original_qualities: 'bool | None' = None, emit_original_quals: 'bool | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): ApplyBQSR.
- **Input:** `input_path` (required): Input BAM/CRAM.; `recalibration_table` (required): Table from `gatk_base_recalibrator`.; `output_path` (required): New recalibrated BAM/CRAM.; `reference_path` (optional): Reference FASTA (needed for CRAM).; `static_quantized_quals` (optional): Static quality levels to quantize to.; `use_original_qualities` (optional): Use the OQ tag as the input qualities.; `emit_original_quals` (optional): Keep the original qualities in the OQ tag.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs.output_path` is the recalibrated file. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_apply_bqsr

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    recalibration_table = folder / "recalibration_table.dat"; recalibration_table.write_bytes(b"x")
    try:
        gatk_apply_bqsr(
            input_path=str(input_path),
            recalibration_table=str(recalibration_table),
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

## 11) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_bqsr_pipeline_spark()`

- **Purpose:** Run both BQSR steps (BaseRecalibrator and ApplyBQSR) in one Spark-parallel pass.
- **Method:** Run `gatk BQSRPipelineSpark` (beta); pass `spark-master` (default `local[*]`) through `options`.
- **API:** `gatk_bqsr_pipeline_spark(input_path: 'str', reference_path: 'str', known_sites: 'Sequence[str]', output_path: 'str', *, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): BQSRPipelineSpark.
- **Input:** `input_path` (required): Input BAM/CRAM.; `reference_path` (required): Reference FASTA (a .2bit file also works for Spark).; `known_sites` (required): Known-variant VCF files (indexed).; `output_path` (required): New recalibrated BAM.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs.output_path` is the recalibrated BAM. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_bqsr_pipeline_spark

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    reference_path = folder / "reference_path.dat"; reference_path.write_bytes(b"x")
    known_sites = folder / "known_sites.dat"; known_sites.write_bytes(b"x")
    try:
        gatk_bqsr_pipeline_spark(
            input_path=str(input_path),
            reference_path=str(reference_path),
            known_sites=[str(known_sites)],
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

## 12) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_bam_to_bfq()`

- **Purpose:** Convert a BAM file to BFQ (binary FASTQ) files for the Maq aligner.
- **Method:** Run Picard `gatk BamToBfq`; the files written to `analysis_dir` are listed in `files`.
- **API:** `gatk_bam_to_bfq(input_path: 'str', analysis_dir: 'str', paired_run: 'bool', *, output_file_prefix: 'str | None' = None, flowcell_barcode: 'str | None' = None, lane: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): BamToBfq.
- **Input:** `input_path` (required): Input BAM/SAM.; `analysis_dir` (required): New directory for the .bfq files.; `output_file_prefix` (optional): Prefix of the .bfq file names (use this or `flowcell_barcode`).; `flowcell_barcode` (optional): Flowcell barcode from which the file names are built (use this or `output_file_prefix`).; `paired_run` (required): Whether the run is paired-end.; `lane` (optional): Lane number to write.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields and `files` (file names inside `analysis_dir`). Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_bam_to_bfq

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    try:
        gatk_bam_to_bfq(
            input_path=str(input_path),
            analysis_dir=str(folder / "analysis_dir"),
            output_file_prefix="reads",
            paired_run=False,
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

## 13) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_bwa_spark()`

- **Purpose:** Align reads to a reference with BWA-MEM on Spark, optionally followed by duplicate marking.
- **Method:** Run `gatk BwaSpark` (`pipeline=align`) or `BwaAndMarkDuplicatesPipelineSpark` (`pipeline=align_mark_duplicates`, name-sorted input); both are beta tools. `pipeline`: `align` = BwaSpark, `align_mark_duplicates` = BwaAndMarkDuplicatesPipelineSpark.
- **API:** `gatk_bwa_spark(input_path: 'str', reference_path: 'str', output_path: 'str', *, pipeline: 'str' = 'align', bwa_mem_index_image: 'str | None' = None, single_end_alignment: 'bool | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): BwaAndMarkDuplicatesPipelineSpark, BwaSpark.
- **Input:** `input_path` (required): Unaligned (or name-sorted) BAM/CRAM.; `reference_path` (required): Reference FASTA (with .fai and .dict).; `output_path` (required): New aligned BAM.; `bwa_mem_index_image` (optional): Prebuilt BWA-MEM index image (see `gatk_bwa_mem_index_image`).; `single_end_alignment` (optional): Align reads as single-end.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields and `mode`; `outputs.output_path` is the aligned BAM. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_bwa_spark

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    reference_path = folder / "reference_path.dat"; reference_path.write_bytes(b"x")
    try:
        gatk_bwa_spark(
            input_path=str(input_path),
            reference_path=str(reference_path),
            output_path=str(folder / "output_path.out"),
            pipeline='align',
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

## 14) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_clean_sam()`

- **Purpose:** Clean a SAM/BAM/CRAM file: soft-clip alignments that hang off the end of the reference and set MAPQ to 0 for unmapped reads.
- **Method:** Run Picard `gatk CleanSam`.
- **API:** `gatk_clean_sam(input_path: 'str', output_path: 'str', *, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): CleanSam.
- **Input:** `input_path` (required): Input SAM/BAM/CRAM.; `output_path` (required): New cleaned file.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs.output_path` is the cleaned file. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_clean_sam

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    try:
        gatk_clean_sam(
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

## 15) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_clip_reads()`

- **Purpose:** Clip reads by sequence, cycle or quality and report the clipping statistics.
- **Method:** Run `gatk ClipReads`; `clip_representation` takes GATK's `WRITE_NS` (default), `WRITE_Q0S`, `WRITE_NS_Q0S`, `SOFTCLIP_BASES`, `HARDCLIP_BASES` or `REVERT_SOFTCLIPPED_BASES`.
- **API:** `gatk_clip_reads(input_path: 'str', output_path: 'str', *, clip_sequences: 'Sequence[str] | None' = None, cycles_to_trim: 'str | None' = None, q_trimming_threshold: 'int | None' = None, clip_representation: 'str | None' = None, statistics_path: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): ClipReads.
- **Input:** `input_path` (required): Input BAM/CRAM/SAM.; `output_path` (required): New clipped file.; `clip_sequences` (optional): Sequences to clip wherever they occur in a read.; `cycles_to_trim` (optional): Cycle ranges to trim, for example `1-5,10-12`.; `q_trimming_threshold` (optional): Trim ends below this base quality.; `clip_representation` (optional): How clipped bases are written.; `statistics_path` (optional): New file for the clipping statistics.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs.output_path` is the clipped file (and `statistics_path` when requested). Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_clip_reads

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    try:
        gatk_clip_reads(
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

## 16) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_convert_headerless_bam_shard()`

- **Purpose:** Turn a headerless BAM shard (Hadoop-BAM output) into a readable BAM by adding a header.
- **Method:** Run `gatk ConvertHeaderlessHadoopBamShardToBam` (beta).
- **API:** `gatk_convert_headerless_bam_shard(shard_path: 'str', header_bam_path: 'str', output_path: 'str', *, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): ConvertHeaderlessHadoopBamShardToBam.
- **Input:** `shard_path` (required): Headerless BAM shard written by a Spark tool.; `header_bam_path` (required): BAM file whose header is used.; `output_path` (required): New readable BAM.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs.output_path` is the BAM. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_convert_headerless_bam_shard

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    shard_path = folder / "shard_path.dat"; shard_path.write_bytes(b"x")
    header_bam_path = folder / "header_bam_path.dat"; header_bam_path.write_bytes(b"x")
    try:
        gatk_convert_headerless_bam_shard(
            shard_path=str(shard_path),
            header_bam_path=str(header_bam_path),
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

## 17) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_downsample_by_duplicate_set()`

- **Purpose:** Discard a set fraction of the duplicate sets of a UMI-grouped BAM.
- **Method:** Run `gatk DownsampleByDuplicateSet` (beta).
- **API:** `gatk_downsample_by_duplicate_set(input_path: 'str', output_path: 'str', *, fraction_to_keep: 'float | None' = None, min_reads: 'int | None' = None, min_per_strand_reads: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): DownsampleByDuplicateSet.
- **Input:** `input_path` (required): UMI-grouped BAM carrying MI tags (`<id>/A`, `<id>/B`) and sorted by MI.; `output_path` (required): New downsampled BAM.; `fraction_to_keep` (optional): Fraction of duplicate sets to keep.; `min_reads` (optional): Minimum reads in a set to be kept.; `min_per_strand_reads` (optional): Minimum reads per strand in a set.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs.output_path` is the downsampled BAM. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_downsample_by_duplicate_set

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    try:
        gatk_downsample_by_duplicate_set(
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

## 18) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_downsample_sam()`

- **Purpose:** Randomly keep a fraction of the reads of a SAM/BAM/CRAM file.
- **Method:** Run Picard `gatk DownsampleSam`; strategy defaults to GATK's own default (`ConstantMemory`).
- **API:** `gatk_downsample_sam(input_path: 'str', output_path: 'str', *, probability: 'float | None' = None, strategy: 'str | None' = None, accuracy: 'float | None' = None, random_seed: 'int | None' = None, metrics_path: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): DownsampleSam.
- **Input:** `input_path` (required): Input SAM/BAM/CRAM.; `output_path` (required): New downsampled file.; `probability` (optional): Probability of keeping each read pair (0-1).; `strategy` (optional): Downsampling strategy.; `accuracy` (optional): Accuracy of the kept fraction (not for ConstantMemory).; `random_seed` (optional): Random seed (default 1).; `metrics_path` (optional): New metrics file.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs.output_path` is the downsampled file. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_downsample_sam

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    try:
        gatk_downsample_sam(
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

## 19) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_position_based_downsample_sam()`

- **Purpose:** Downsample reads by their physical position on the flowcell tile (keeps whole neighbourhoods of clusters, for duplicate studies).
- **Method:** Run Picard `gatk PositionBasedDownsampleSam`.
- **API:** `gatk_position_based_downsample_sam(input_path: 'str', output_path: 'str', fraction: 'float', *, stop_after: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): PositionBasedDownsampleSam.
- **Input:** `input_path` (required): Coordinate-sorted input with Illumina-style read names.; `output_path` (required): New downsampled file.; `fraction` (required): Fraction of the read positions (tile areas) to keep (0-1).; `stop_after` (optional): Stop after this many reads (0 = all).; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs.output_path` is the downsampled file. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_position_based_downsample_sam

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    try:
        gatk_position_based_downsample_sam(
            input_path=str(input_path),
            output_path=str(folder / "output_path.out"),
            fraction=0.5,
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

## 20) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_extract_original_alignment_records_by_name()`

- **Purpose:** Subset a BAM to the reads whose names are listed in a file (with their original alignments).
- **Method:** Run `gatk ExtractOriginalAlignmentRecordsByNameSpark` (beta).
- **API:** `gatk_extract_original_alignment_records_by_name(input_path: 'str', read_name_file: 'str', output_path: 'str', *, invert_match: 'bool | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): ExtractOriginalAlignmentRecordsByNameSpark.
- **Input:** `input_path` (required): Input BAM.; `read_name_file` (required): Text file with one read name per line.; `output_path` (required): New BAM with the selected reads.; `invert_match` (optional): Keep the reads that are NOT listed.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs.output_path` is the subset. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_extract_original_alignment_records_by_name

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    read_name_file = folder / "read_name_file.dat"; read_name_file.write_bytes(b"x")
    try:
        gatk_extract_original_alignment_records_by_name(
            input_path=str(input_path),
            read_name_file=str(read_name_file),
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

## 21) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_filter_sam_reads()`

- **Purpose:** Subset the reads of a SAM/BAM/CRAM file by aligned status, read name, interval or tag value.
- **Method:** Run Picard `gatk FilterSamReads`. The `includeJavascript` filter is not offered because it would execute a user script.
- **API:** `gatk_filter_sam_reads(input_path: 'str', output_path: 'str', filter: 'str', *, read_list_file: 'str | None' = None, interval_list: 'str | None' = None, tag: 'str | None' = None, tag_values: 'Sequence[str] | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): FilterSamReads.
- **Input:** `input_path` (required): Input SAM/BAM/CRAM.; `output_path` (required): New filtered file.; `filter` (required): Which filter to apply.; `read_list_file` (optional): Read names for the read-list filters.; `interval_list` (optional): Intervals for `includePairedIntervals`.; `tag` (optional): Tag name for the tag-value filters.; `tag_values` (optional): Tag values for the tag-value filters.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs.output_path` is the filtered file. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_filter_sam_reads

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    try:
        gatk_filter_sam_reads(
            input_path=str(input_path),
            output_path=str(folder / "output_path.out"),
            filter="includeAligned",
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

## 22) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_fix_misencoded_base_quality_reads()`

- **Purpose:** Fix Illumina base quality scores that were stored with the old Phred+64 encoding.
- **Method:** Run `gatk FixMisencodedBaseQualityReads`.
- **API:** `gatk_fix_misencoded_base_quality_reads(input_path: 'str', output_path: 'str', *, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): FixMisencodedBaseQualityReads.
- **Input:** `input_path` (required): Input BAM/CRAM/SAM with Illumina-encoded (Phred+64) qualities.; `output_path` (required): New file with Phred+33 qualities.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs.output_path` is the fixed file. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_fix_misencoded_base_quality_reads

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    try:
        gatk_fix_misencoded_base_quality_reads(
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

## 23) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_left_align_indels()`

- **Purpose:** Left-align the indels in the CIGARs of the reads of a SAM/BAM/CRAM file.
- **Method:** Run `gatk LeftAlignIndels`.
- **API:** `gatk_left_align_indels(input_path: 'str', reference_path: 'str', output_path: 'str', *, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): LeftAlignIndels.
- **Input:** `input_path` (required): Input BAM/CRAM/SAM.; `reference_path` (required): Reference FASTA (with .fai and .dict).; `output_path` (required): New file with left-aligned indels.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs.output_path` is the new file. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_left_align_indels

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    reference_path = folder / "reference_path.dat"; reference_path.write_bytes(b"x")
    try:
        gatk_left_align_indels(
            input_path=str(input_path),
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

## 24) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_mark_duplicates_with_mate_cigar()`

- **Purpose:** Mark duplicate reads using the mate CIGAR (standard, simplified or UMI-aware).
- **Method:** Run Picard `MarkDuplicatesWithMateCigar`, `SimpleMarkDuplicatesWithMateCigar` or `UmiAwareMarkDuplicatesWithMateCigar` depending on `variant`; the metrics file is parsed into `tables.metrics_path`. `variant`: `standard` = MarkDuplicatesWithMateCigar, `simple` = SimpleMarkDuplicatesWithMateCigar (experimental), `umi_aware` = UmiAwareMarkDuplicatesWithMateCigar (experimental).
- **API:** `gatk_mark_duplicates_with_mate_cigar(input_path: 'str', output_path: 'str', metrics_path: 'str', *, variant: 'str' = 'standard', umi_metrics_path: 'str | None' = None, remove_duplicates: 'bool | None' = None, umi_tag_name: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): MarkDuplicatesWithMateCigar, SimpleMarkDuplicatesWithMateCigar, UmiAwareMarkDuplicatesWithMateCigar.
- **Input:** `input_path` (required): Coordinate-sorted BAM with MC (mate CIGAR) tags.; `output_path` (required): New BAM with duplicates marked.; `metrics_path` (required): New duplication metrics file.; `umi_metrics_path` (required in some modes): New UMI metrics file.; `remove_duplicates` (optional): Drop duplicates instead of flagging them.; `umi_tag_name` (optional): Read tag that holds the UMI (default RX).; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields, `mode` and `tables.metrics_path` (the Picard metrics rows). Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_mark_duplicates_with_mate_cigar

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    try:
        gatk_mark_duplicates_with_mate_cigar(
            input_path=str(input_path),
            output_path=str(folder / "output_path.out"),
            metrics_path=str(folder / "metrics_path.out"),
            variant='standard',
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

## 25) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_merge_bam_alignment()`

- **Purpose:** Merge the alignment data of an aligner BAM with the original data of the unmapped BAM.
- **Method:** Run Picard `gatk MergeBamAlignment`.
- **API:** `gatk_merge_bam_alignment(unmapped_bam: 'str', reference_path: 'str', output_path: 'str', *, aligned_bams: 'Sequence[str] | None' = None, paired_run: 'bool | None' = None, sort_order: 'str | None' = None, add_mate_cigar: 'bool | None' = None, clip_adapters: 'bool | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): MergeBamAlignment.
- **Input:** `unmapped_bam` (required): Original unmapped BAM (query-name sorted).; `reference_path` (required): Reference FASTA (with .fai and .dict).; `output_path` (required): New merged BAM.; `aligned_bams` (optional): Aligner output files (query-name sorted); omit for unmapped-only output.; `paired_run` (optional): Whether the run is paired-end.; `sort_order` (optional): Sort order of the output.; `add_mate_cigar` (optional): Add the MC tag.; `clip_adapters` (optional): Clip adapters marked by MarkIlluminaAdapters.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs.output_path` is the merged BAM. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_merge_bam_alignment

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    unmapped_bam = folder / "unmapped_bam.dat"; unmapped_bam.write_bytes(b"x")
    reference_path = folder / "reference_path.dat"; reference_path.write_bytes(b"x")
    try:
        gatk_merge_bam_alignment(
            unmapped_bam=str(unmapped_bam),
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

## 26) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_print_distant_mates()`

- **Purpose:** Unmap reads whose mates are mapped far away and write them out.
- **Method:** Run `gatk PrintDistantMates`.
- **API:** `gatk_print_distant_mates(input_path: 'str', output_path: 'str', *, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): PrintDistantMates.
- **Input:** `input_path` (required): Input BAM/CRAM/SAM.; `output_path` (required): New file.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs.output_path` is the new file. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_print_distant_mates

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    try:
        gatk_print_distant_mates(
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

## 27) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_reorder_sam()`

- **Purpose:** Reorder the reads of a SAM/BAM/CRAM file to match the contig order of another reference.
- **Method:** Run Picard `gatk ReorderSam`.
- **API:** `gatk_reorder_sam(input_path: 'str', output_path: 'str', sequence_dictionary: 'str', *, allow_incomplete_dict_concordance: 'bool | None' = None, allow_contig_length_discordance: 'bool | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): ReorderSam.
- **Input:** `input_path` (required): Input SAM/BAM/CRAM.; `output_path` (required): New reordered file.; `sequence_dictionary` (required): SAM/BAM/CRAM/dict/fasta whose contig order is the target.; `allow_incomplete_dict_concordance` (optional): Allow contigs missing from the target.; `allow_contig_length_discordance` (optional): Allow contigs of different length.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs.output_path` is the reordered file. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_reorder_sam

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    sequence_dictionary = folder / "sequence_dictionary.dat"; sequence_dictionary.write_bytes(b"x")
    try:
        gatk_reorder_sam(
            input_path=str(input_path),
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

## 28) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_revert_base_quality_scores()`

- **Purpose:** Revert base qualities to the original values stored in the OQ tag (optionally adding the mate CIGAR tag).
- **Method:** Run `gatk RevertBaseQualityScores` (`mode=original_qualities`) or Picard `RevertOriginalBaseQualitiesAndAddMateCigar` (`mode=with_mate_cigar`). `mode`: `original_qualities` = RevertBaseQualityScores, `with_mate_cigar` = RevertOriginalBaseQualitiesAndAddMateCigar.
- **API:** `gatk_revert_base_quality_scores(input_path: 'str', output_path: 'str', *, mode: 'str' = 'original_qualities', restore_original_qualities: 'bool | None' = None, sort_order: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): RevertBaseQualityScores, RevertOriginalBaseQualitiesAndAddMateCigar.
- **Input:** `input_path` (required): Input BAM/CRAM/SAM that carries OQ tags.; `output_path` (required): New file.; `restore_original_qualities` (optional): Restore the OQ values.; `sort_order` (optional): Sort order of the output.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields and `mode`; `outputs.output_path` is the new file. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_revert_base_quality_scores

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    try:
        gatk_revert_base_quality_scores(
            input_path=str(input_path),
            output_path=str(folder / "output_path.out"),
            mode='original_qualities',
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

## 29) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_revert_sam()`

- **Purpose:** Revert a SAM/BAM/CRAM file to an earlier state (unmapped, original qualities, no duplicate flags) for re-alignment.
- **Method:** Run Picard `gatk RevertSam` with a single output file; `OUTPUT_BY_READGROUP` is not offered by this function.
- **API:** `gatk_revert_sam(input_path: 'str', output_path: 'str', *, restore_original_qualities: 'bool | None' = None, remove_duplicate_information: 'bool | None' = None, remove_alignment_information: 'bool | None' = None, sanitize: 'bool | None' = None, sample_alias: 'str | None' = None, library_name: 'str | None' = None, sort_order: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): RevertSam.
- **Input:** `input_path` (required): Input SAM/BAM/CRAM.; `output_path` (required): New reverted file.; `restore_original_qualities` (optional): Restore OQ qualities (Picard default true).; `remove_duplicate_information` (optional): Clear the duplicate flag.; `remove_alignment_information` (optional): Remove alignment information.; `sanitize` (optional): Remove reads that cannot be reverted consistently.; `sample_alias` (optional): Sample name for the read groups.; `library_name` (optional): Library name for the read groups.; `sort_order` (optional): Sort order of the output.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs.output_path` is the reverted file. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_revert_sam

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    try:
        gatk_revert_sam(
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

## 30) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_convert_sam_format()`

- **Purpose:** Convert a SAM file to BAM or a BAM file to SAM (or CRAM).
- **Method:** Run Picard `gatk SamFormatConverter`.
- **API:** `gatk_convert_sam_format(input_path: 'str', output_path: 'str', *, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): SamFormatConverter.
- **Input:** `input_path` (required): Input SAM/BAM/CRAM.; `output_path` (required): New file; the format follows its extension (.sam, .bam, .cram).; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs.output_path` is the converted file. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_convert_sam_format

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    try:
        gatk_convert_sam_format(
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

## 31) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_sam_to_fastq_with_tags()`

- **Purpose:** Convert a SAM/BAM/CRAM file to FASTQ and also write FASTQ files built from read tags (for example UMIs).
- **Method:** Run Picard `gatk SamToFastqWithTags`.
- **API:** `gatk_sam_to_fastq_with_tags(input_path: 'str', fastq_path: 'str', sequence_tag_groups: 'Sequence[str]', *, quality_tag_groups: 'Sequence[str] | None' = None, second_end_fastq: 'str | None' = None, unpaired_fastq: 'str | None' = None, interleave: 'bool | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): SamToFastqWithTags.
- **Input:** `input_path` (required): Input SAM/BAM/CRAM.; `fastq_path` (required): New FASTQ for the first read (or all reads when interleaved).; `sequence_tag_groups` (required): Tags whose values are written as extra FASTQ files; join several tags with commas.; `quality_tag_groups` (optional): Quality tags that go with the sequence tag groups.; `second_end_fastq` (optional): New FASTQ for the second read.; `unpaired_fastq` (optional): New FASTQ for unpaired reads.; `interleave` (optional): Write both reads of a pair to one FASTQ.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs` lists the FASTQ files that were requested. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_sam_to_fastq_with_tags

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    try:
        gatk_sam_to_fastq_with_tags(
            input_path=str(input_path),
            fastq_path=str(folder / "fastq_path.out"),
            sequence_tag_groups=["RX"],
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

## 32) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_split_n_cigar_reads()`

- **Purpose:** Split RNA-seq reads at the N operators of their CIGARs into exon segments and clip overhangs.
- **Method:** Run `gatk SplitNCigarReads`.
- **API:** `gatk_split_n_cigar_reads(input_path: 'str', reference_path: 'str', output_path: 'str', *, max_bases_in_overhang: 'int | None' = None, max_mismatches_in_overhang: 'int | None' = None, skip_mapping_quality_transform: 'bool | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): SplitNCigarReads.
- **Input:** `input_path` (required): RNA-seq BAM/CRAM/SAM with N operators in CIGARs.; `reference_path` (required): Reference FASTA (with .fai and .dict).; `output_path` (required): New file.; `max_bases_in_overhang` (optional): Maximum overhang bases before a split piece is hard-clipped.; `max_mismatches_in_overhang` (optional): Maximum mismatches allowed in an overhang.; `skip_mapping_quality_transform` (optional): Do not change MAPQ 255 to 60.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs.output_path` is the split file. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_split_n_cigar_reads

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    reference_path = folder / "reference_path.dat"; reference_path.write_bytes(b"x")
    try:
        gatk_split_n_cigar_reads(
            input_path=str(input_path),
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

## 33) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_split_sam()`

- **Purpose:** Split a SAM/BAM/CRAM file into several files, by library name or by number of reads.
- **Method:** Run Picard `SplitSamByLibrary` (`by=library`, reads without a library go to `unknown`) or `SplitSamByNumberOfReads` (`by=number_of_reads`); the files written are listed in `files`. `by`: `library` = SplitSamByLibrary, `number_of_reads` = SplitSamByNumberOfReads.
- **API:** `gatk_split_sam(input_path: 'str', output_dir: 'str', *, by: 'str' = 'library', split_to_n_files: 'int | None' = None, split_to_n_reads: 'int | None' = None, total_reads_in_input: 'int | None' = None, out_prefix: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): SplitSamByLibrary, SplitSamByNumberOfReads.
- **Input:** `input_path` (required): Input SAM/BAM/CRAM (query-grouped for `number_of_reads`).; `output_dir` (required): New directory for the split files.; `split_to_n_files` (optional): Number of output files.; `split_to_n_reads` (optional): Reads per output file.; `total_reads_in_input` (optional): Total reads, which saves a counting pass.; `out_prefix` (optional): File name prefix (default `shard`).; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields, `mode` and `files` (file names inside `output_dir`). Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_split_sam

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    try:
        gatk_split_sam(
            input_path=str(input_path),
            output_dir=str(folder / "output_dir"),
            by='library',
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

## 34) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_transfer_read_tags()`

- **Purpose:** Copy read tags (for example UMIs) from an unmapped SAM/BAM onto the matching reads of an aligned file.
- **Method:** Run `gatk TransferReadTags` (experimental).
- **API:** `gatk_transfer_read_tags(input_path: 'str', unmapped_sam: 'str', output_path: 'str', read_tags: 'Sequence[str]', *, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): TransferReadTags.
- **Input:** `input_path` (required): Aligned BAM/SAM (query-name sorted) that receives the tags.; `unmapped_sam` (required): Query-name sorted unmapped SAM/BAM that holds the tags.; `output_path` (required): New output file.; `read_tags` (required): Tag names to transfer.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs.output_path` is the new file. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_transfer_read_tags

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    unmapped_sam = folder / "unmapped_sam.dat"; unmapped_sam.write_bytes(b"x")
    try:
        gatk_transfer_read_tags(
            input_path=str(input_path),
            unmapped_sam=str(unmapped_sam),
            output_path=str(folder / "output_path.out"),
            read_tags=["RX"],
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

## 35) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_unmark_duplicates()`

- **Purpose:** Clear the duplicate flag (0x400) of all reads.
- **Method:** Run `gatk UnmarkDuplicates`.
- **API:** `gatk_unmark_duplicates(input_path: 'str', output_path: 'str', *, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the external launcher `gatk` (GATK 4.7.0.0 tested, Java 17) on PATH, in `DNAKIT_GATK` or passed as `gatk_path`; GATK is not installed in the environment that produced this page, so the runnable example shows the validated missing-launcher behaviour (`BackendUnavailableError`, code `GATK_UNAVAILABLE`). GATK tool(s): UnmarkDuplicates.
- **Input:** `input_path` (required): Input BAM/CRAM/SAM.; `output_path` (required): New file.; `options`: GATK long options of this tool as `{name: value}` (only options listed by `gatk <Tool> --help`; `True`/`False` become `true`/`false`, a list repeats the option); `gatk_path`, `java_max_memory` (for example `'4g'`), `verbose`..
- **Output:** The common result fields; `outputs.output_path` is the new file. Common fields: `status`, `backend`, `operation`, `tool`, `version`, `command`, `outputs`, `log_tail`, `provenance`.
- **Verification:** Verified: run to completion with GATK 4.7.0.0 on toy data (see the verification folder of the pack).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_unmark_duplicates

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    try:
        gatk_unmark_duplicates(
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
