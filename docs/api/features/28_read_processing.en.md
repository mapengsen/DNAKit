# Sequencing read processing

## 1) Pair reads by ID {#bio-10}

- **Purpose:** Match R1/R2 records by name and identify orphans.
- **Method:** Match after removing trailing `/1` or `/2`; a missing mate is `None`.
- **API:** `dnakit.bioinformatics.reads.pair_reads(first, second, max_pending=100000, include_orphans=True)`
- **Input:** Two iterators of `Read` objects with unique normalized IDs in each mate stream.
- **Output:** Returns an iterator of `ReadPair` entries with a pair ID and both reads.
- **Example code:**

```python
from dnakit.bioinformatics.reads import Read, pair_reads

pairs = pair_reads([Read("r1/1", "ACGT")], [Read("r1/2", "TGCA")])
print([(pair.id, pair.first.sequence, pair.second.sequence) for pair in pairs])
```

- **Example output:**

```text
[('r1', 'ACGT', 'TGCA')]
```

[Dependencies and implementation notes](../../faq/bioinformatics.md#reads).

## 6) Long-read filtering and summaries {#bio-long-read-qc}

DNAKit provides one bounded external boundary for NanoFilt, Filtlong and
NanoPlot. The caller installs each executable separately; DNAKit does not copy
the GPL/MIT backends or download databases.

- **NanoFilt API:** `dnakit.sequencing_qc.filter_nanopore_reads(source, output_directory, min_quality=None, min_length=None, headcrop=None, tailcrop=None, ...)`.
- **Filtlong API:** `dnakit.sequencing_qc.filter_long_reads(source, output_directory, min_length=None, keep_percent=None, target_bases=None, ...)`.
- **NanoPlot API:** `dnakit.sequencing_qc.summarize_long_reads(source, output_directory, threads=1, min_length=None, min_quality=None, ...)`.
- **Porechop API:** `dnakit.sequencing_qc.trim_porechop(source, output_directory, threads=1, no_split=False, ...)`.
- **Trimmomatic API:** `dnakit.sequencing_qc.trim_trimmomatic(read1, output_directory, read2=None, adapter_fasta=None, leading=None, trailing=None, min_length=None, ...)`.
- **AfterQC API:** `dnakit.sequencing_qc.run_afterqc(read1, output_directory, read2=None, qc_only=False, ...)`.
- **Input and output:** Both filters read a local FASTQ and capture `filtered.fastq` plus stderr. NanoPlot runs in FASTQ mode with static plots disabled and checks `NanoStats.txt` and `NanoPlot-report.html`. The output directory must be new; paths, numeric options, timeouts and artifact sizes are validated.
- **Agent/MCP:** All three functions are catalogued as `external_process`; callers must explicitly allow the external side effect. The boundary covers fixed quality, length, cropping, retention, target-base and thread options.
- **Citation and license:** NanoFilt: De Coster et al. (2018), [10.1093/bioinformatics/bty149](https://doi.org/10.1093/bioinformatics/bty149), GPL-3.0; Filtlong: no DOI supplied by the official project, GPL-3.0; NanoPlot: De Coster & Rademakers (2023), [10.1093/bioinformatics/btad311](https://doi.org/10.1093/bioinformatics/btad311), MIT. The NanoFilt repository states that the project is no longer updated.
- **Porechop citation and license:** Wick et al. (2017), [10.1099/mgen.0.000132](https://doi.org/10.1099/mgen.0.000132), GPL-3.0; the official repository states that the project is unsupported.
- **Trimmomatic citation and license:** Bolger, Lohse & Usadel (2014), [10.1093/bioinformatics/btu170](https://doi.org/10.1093/bioinformatics/btu170), GPL-3.0; Illumina supplies the adapter sequences under separate upstream terms.
- **AfterQC citation and license:** Chen et al. (2017), [10.1186/s12859-017-1469-3](https://doi.org/10.1186/s12859-017-1469-3), MIT. AfterQC's official README recommends fastp for a maintained C++ implementation; this facade remains available for reproducibility of published AfterQC workflows.

These facades do not claim the complete upstream command surfaces. Use the
upstream CLI directly for Filtlong reference inputs or unit suffixes and for
NanoPlot input formats beyond FASTQ, while retaining the upstream citation.

## 2) Recover valid records from damaged FASTQ {#bio-11}

- **Purpose:** Recover recognizable four-line FASTQ records from a damaged file.
- **Method:** Scan for headers, sequences, separators and matching quality lengths, counting discarded lines.
- **API:** `dnakit.bioinformatics.reads.repair_fastq(source, destination, phred_offset=33)`
- **Input:** A damaged FASTQ and a new output path; use strict `read_fastx` for multiline FASTQ.
- **Output:** Returns `RepairResult` with the repaired FASTQ path, recovered-read count, and discarded-line count.
- **Example code:**

```python
from tempfile import TemporaryDirectory
from pathlib import Path
from dnakit.bioinformatics.reads import repair_fastq

with TemporaryDirectory() as folder:
    source = Path(folder) / "broken.fq"
    source.write_text("invalid\n@r1\nACGT\n+\nIIII\n")
    result = repair_fastq(source, Path(folder) / "repaired.fq")
    print(result.recovered_reads, result.discarded_lines)
```

- **Example output:**

```text
1 1
```

[Dependencies and implementation notes](../../faq/bioinformatics.md#reads).

## 3) Sort, intersect, rename and fingerprint reads {#bio-12}

- **Purpose:** Organize records, extract shared sequences and compute a sequence digest.
- **Method:** Sort or shuffle with a seed, intersect by ID/sequence, rename with regex, and hash length-delimited sequences.
- **API:** `dnakit.bioinformatics.reads.arrange_reads(reads, operation='sort')`; `dnakit.bioinformatics.reads.common_reads(collections)`; `dnakit.bioinformatics.reads.rename_reads(reads, pattern, replacement)`; `dnakit.bioinformatics.reads.sequence_digest(reads)`
- **Input:** Collections of `Read` objects; the digest depends on record order, not IDs.
- **Output:** Returns sorted, intersected, or renamed `Read` records; the check function returns a comparison result.
- **Example code:**

```python
from dnakit.bioinformatics.reads import Read, arrange_reads, common_reads, rename_reads, sequence_digest

reads = [Read("b", "TGCA"), Read("a", "ACGT")]
ordered = arrange_reads(reads)
shared = common_reads([ordered, [Read("a", "ACGT")]])
renamed = list(rename_reads(shared, "^a$", "sample_a"))
print([read.id for read in ordered])
print([read.id for read in renamed])
print(sequence_digest(shared) == sequence_digest(renamed))
```

- **Example output:**

```text
['a', 'b']
['sample_a']
True
```

[Dependencies and implementation notes](../../faq/bioinformatics.md#reads).

## 4) Prepare, mark duplicates and index BAM {#bio-16}

- **Purpose:** Prepare SAM/BAM for variant calling.
- **Method:** Run samtools name-sort, fixmate, coordinate sort, markdup, index and stats.
- **API:** `dnakit.bioinformatics.ngs.prepare_bam(source, output_directory, mark_duplicates=True, config=None)`
- **Input:** Existing `aligned.sam` or BAM, samtools on PATH, and a new output directory.
- **Output:** Returns `ToolResult` with `backend`, artifact paths in `artifacts`, and execution `log`.
- **Agent/MCP API:** `dnakit.ngs.prepare_bam(source, output_directory, mark_duplicates=True, ...)`; returns JSON-serializable artifact paths, backend metadata, and diagnostics.
- **Citation and license:** cite SAMtools as Danecek et al. (2021), DOI `10.1093/gigascience/giab008` (MIT). The caller installs the backend.
- **Example code:**

```python
from dnakit.bioinformatics.ngs import prepare_bam

result = prepare_bam("aligned.sam", "new_bam")
print(result.artifacts["aligned.bam"])
```

- **Example output:**

Returns `ToolResult` containing `aligned.bam`, its index and statistics files. This external-program example was not executed.

[Dependencies and implementation notes](../../faq/bioinformatics.md#reads).

## 5) External sequencing QC tools {#bio-external-qc}

DNAKit provides one Python/Agent/MCP boundary for Cutadapt, fastp, FastQC and MultiQC. These are real bounded executable adapters: the caller installs each backend, while DNAKit does not bundle the executable, Java, plugins or databases.

- **Shared method:** DNAKit constructs an allowlisted argv vector and uses a no-shell bounded executor. It checks the executable, input paths, timeout, output size, exit status and declared artifacts, and refuses to overwrite existing outputs. Arbitrary command-line flags are rejected.
- **Shared output:** Success returns a JSON-serializable envelope with `status`, `backend`, exact argv tokens, `outputs`, bounded backend output and `provenance`. Missing backends, non-zero exits, missing/empty artifacts, invalid JSON reports and resource-limit violations become structured errors.
- **Agent/MCP:** All four functions are catalogued with an `external_process` effect; MCP callers still need to explicitly allow the external side effect. The boundary covers a fixed DNAKit option subset and does not claim to re-export each upstream CLI in full.

### Cutadapt

- **API:** `dnakit.sequencing_qc.run_cutadapt(input, *, options=None, verbose=True)`.
- **Input:** `input` is a FASTQ path or `{"read1": path, "read2": path}`. `options` requires `executable_path` and `output`/`output1`; paired input also requires `output2`, plus at least one adapter such as `adapters=["AGATCGGAAG"]`. Fixed options include `cores`, `minimum_overlap`, `minimum_length`, `error_rate`, `quality_cutoff`, `discard_untrimmed` and `json_report`.
- **Output:** One or two trimmed FASTQ files and an optional JSON report. All outputs must be in an existing directory and must not overwrite existing files.
- **Citation and license:** Martin (2011), [10.14806/ej.17.1.200](https://doi.org/10.14806/ej.17.1.200); upstream Cutadapt license: MIT.

### fastp

- **API:** `dnakit.sequencing_qc.run_fastp(input, *, options=None, verbose=True)`.
- **Input:** Single-end or paired `{"read1": path, "read2": path}` FASTQ files. `options` requires `executable_path` and `output`/`output1`, with `output2` for paired input; fixed quality, length, thread, adapter and filtering options are supported.
- **Output:** One or two filtered FASTQ files plus JSON and HTML reports. When report paths are omitted, DNAKit chooses non-conflicting defaults beside the output FASTQ.
- **Citation and license:** Chen (2025), [10.1002/imt2.70078](https://doi.org/10.1002/imt2.70078); upstream fastp license: MIT.

### FastQC

- **API:** `dnakit.sequencing_qc.run_fastqc(input, *, options=None, verbose=True)`.
- **Input:** One to eight local FASTQ/BAM/SAM-style paths. `options` requires `executable_path` and an existing empty `output_directory`. The adapter follows the FastQC 0.12.x CLI for `threads`, `noextract`, `nogroup`, `quiet`, `format`, `contaminant`, `adapter` and a 2–10 `kmer` value.
- **Output:** `<sample>_fastqc.html` and `<sample>_fastqc.zip` for every input, with checks that both files were created. One call monitors at most eight inputs (sixteen reports) to preserve artifact limits.
- **Citation and license:** [FastQC project page](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/); no DOI was supplied by the survey; upstream license: GPL-3.0. The option contract targets FastQC 0.12.x and should be checked separately for other versions.

### MultiQC

- **API:** `dnakit.sequencing_qc.run_multiqc(input, *, options=None, verbose=True)`.
- **Input:** One or more local files/directories. `options` requires `executable_path` and an existing empty `output_directory`, with `filename` selecting a `.html` report. Fixed options include `title`, `interactive`, `zip_data`, `flat`, `dirs`, `dirs_depth` and a local `config`; `zip_data` maps to MultiQC `--zip-data-dir`.
- **Output:** The declared HTML summary report. An empty output directory lets DNAKit detect a zero-exit backend that produced no report instead of treating it as success.
- **Citation and license:** Ewels et al. (2016), [10.1093/bioinformatics/btw354](https://doi.org/10.1093/bioinformatics/btw354); upstream MultiQC license: GPL-3.0-or-later.

**Example code (Cutadapt):**

```python
from pathlib import Path
from dnakit.sequencing_qc import run_cutadapt

result = run_cutadapt(
    {"read1": Path("reads_R1.fastq"), "read2": Path("reads_R2.fastq")},
    options={
        "executable_path": "/usr/local/bin/cutadapt",
        "output1": Path("results/trimmed_R1.fastq"),
        "output2": Path("results/trimmed_R2.fastq"),
        "adapters": ["AGATCGGAAG"],
        "json_report": Path("results/cutadapt.json"),
    },
    verbose=False,
)
print(result["outputs"])
```

[Dependencies and implementation notes](../../faq/bioinformatics.md#reads).
