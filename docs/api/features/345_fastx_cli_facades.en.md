# FASTA/FASTQ command-line facades: seqmagick, SeqKit, AfterQC and Trimmomatic

This group of functions runs four caller-installed programs through fixed argument whitelists, without a shell: [seqmagick](https://github.com/fhcrc/seqmagick) (residue slicing, alignment clean-up, record selection, quality filter, file summary), [SeqKit](https://github.com/shenwei356/seqkit) (FASTA index, sanitising, streaming, common records, header replacement, BAM toolbox), [AfterQC](https://github.com/OpenGene/AfterQC) (overlap analysis, debubble, barcodes) and [Trimmomatic](https://github.com/usadellab/Trimmomatic) (trimming and filtering steps, mate pairing). The programs are not bundled: a missing program raises an error that carries the install command. Output files must not exist yet; nothing is overwritten. Every result is a dictionary with `status`, `backend`, `version`, `command`, `outputs`, parsed tables and `log_tail`.

Tested versions: seqmagick 0.8.6, SeqKit 2.14.0, Trimmomatic 0.39. AfterQC 0.9.7 is a Python 2 program that was not available in the test environment, so its three functions are tested for argument handling only. The examples run seqmagick, SeqKit and Trimmomatic.

## Example data

All inputs are tiny text files; the examples write their results to temporary folders.

**Example file:** `aln.fa`

```text
>s1
AC--GTACGT
>s2
AC--GTTCGT
>s3
ACGTGTACGA
>s4
----------
```

**Example file:** `x.fa`

```text
>a first
ACGTACGTAC
>b second
GGGGCCCCAA
>c third
TTTTTTTTTT
```

**Example file:** `other.fa`

```text
>a other
ACGTACGTAC
>x only
CCCCCCCCCC
>c third
AAAAAAAAAA
```

**Example file:** `reads.fq`

```text
@r1
CAGATGACGTACGTACGT
+
IIIIIIIIIIIIIIIIII
@r2
TTTTACGTACGTACGTAC
+
IIIIIIIIIIIIIIIIII
@r3
CAGATGACNTACGTACGT
+
IIIIIIIIIIIIIIIIII
```

**Example file:** `broken.fq`

```text
@r1
ACGT
+
IIII
@bad
AC
+
IIII
@r3
ACGTN
+
IIIII
nonsense
```

**Example file:** `r1.fq`

```text
@p1/1
ACGTACGTACGTACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@p2/1
TTGCATTGCATTGCATTGCATTGCATTGCA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII
```

**Example file:** `r2.fq`

```text
@p1/2
GGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@p2/2
ACACAACACAACACAACACAACACAACACA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII
```

**Example file:** `m1.fq`

```text
@m1/1
ACGTAC
+
IIIIII
@m2/1
ACGTAC
+
IIIIII
@m3/1
ACGTAC
+
IIIIII
```

**Example file:** `m2.fq`

```text
@m1/2
GGGGGG
+
IIIIII
@m3/2
TTTTTT
+
IIIIII
```

## 1) `dnakit.fasta_fastq.seqmagick_p15_cli_functions.seqmagick_slice_residues()`

- **Purpose:** Cut, delete or mask residues of every sequence by position range.
- **Method:** Runs `seqmagick convert --cut/--drop/--mask`. Ranges are 1-based, inclusive `start:end`, several separated by commas; `relative_to` takes numbering from a gapped sequence.
- **API:** `dnakit.fasta_fastq.seqmagick_p15_cli_functions.seqmagick_slice_residues(source, destination, *, cut=None, drop=None, mask=None, relative_to=None, input_format=None, output_format=None, alphabet=None, line_wrap=None, seqmagick_path=None)`.
- **Input:** A sequence file, a new destination path and at least one of `cut`, `drop`, `mask`.
- **Output:** `output`, `output_record_count`, command and version.
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.fasta_fastq.seqmagick_p15_cli_functions import seqmagick_slice_residues

result = seqmagick_slice_residues("aln.fa", str(out / "cut.fa"), cut="3:6")
print(result["output_record_count"], result["version"])
print((out / "cut.fa").read_text().splitlines()[1])
```

- **Example output:**

```text
4 0.8.6
--GT
```

## 2) `dnakit.fasta_fastq.seqmagick_p15_cli_functions.seqmagick_clean_alignment()`

- **Purpose:** Normalise gap characters, remove gap-dominated columns and drop all-gap sequences.
- **Method:** Runs `seqmagick convert --dash-gap`, `--squeeze` / `--squeeze-threshold` and `--prune-empty`.
- **API:** `dnakit.fasta_fastq.seqmagick_p15_cli_functions.seqmagick_clean_alignment(source, destination, *, dash_gap=False, squeeze=False, squeeze_threshold=None, prune_empty=False, input_format=None, output_format=None, alphabet=None, line_wrap=None, seqmagick_path=None)`.
- **Input:** An alignment file, a new destination and at least one option.
- **Output:** `output`, `output_record_count`.
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.fasta_fastq.seqmagick_p15_cli_functions import seqmagick_clean_alignment

result = seqmagick_clean_alignment("aln.fa", str(out / "clean.fa"), squeeze_threshold=0.75, prune_empty=True)
print(result["output_record_count"])
print((out / "clean.fa").read_text().splitlines()[1])
```

- **Example output:**

```text
3
ACGTACGT
```

## 3) `dnakit.fasta_fastq.seqmagick_p15_cli_functions.seqmagick_select_sequences()`

- **Purpose:** Keep or drop sequences by an ID list or by a regular expression on the sequence.
- **Method:** Runs `seqmagick convert` with `--include-from-file`, `--exclude-from-file`, `--seq-pattern-include` or `--seq-pattern-exclude`.
- **API:** `dnakit.fasta_fastq.seqmagick_p15_cli_functions.seqmagick_select_sequences(source, destination, *, include_ids_file=None, exclude_ids_file=None, seq_pattern_include=None, seq_pattern_exclude=None, input_format=None, output_format=None, alphabet=None, line_wrap=None, seqmagick_path=None)`.
- **Input:** A sequence file, a new destination and at least one selector.
- **Output:** `output`, `output_record_count`.
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.fasta_fastq.seqmagick_p15_cli_functions import seqmagick_select_sequences

ids = out / "ids.txt"
ids.write_text("s1\ns3\n")
result = seqmagick_select_sequences("aln.fa", str(out / "kept.fa"), include_ids_file=str(ids))
print(result["output_record_count"])
result = seqmagick_select_sequences("aln.fa", str(out / "pattern.fa"), seq_pattern_include="GTAC")
print(result["output_record_count"])
```

- **Example output:**

```text
2
2
```

## 4) `dnakit.fasta_fastq.seqmagick_p15_cli_functions.seqmagick_quality_filter()`

- **Purpose:** Filter reads by quality, length, ambiguous bases and a leading IUPAC primer.
- **Method:** Runs `seqmagick quality-filter` with its defaults (mean quality 25, length 200 to 1000). `ambiguous_action` truncates or drops reads at `N`. A primer is checked through seqmagick's barcode machinery, so a one-row barcode table with an empty barcode is supplied; passing reads lose the primer.
- **API:** `dnakit.fasta_fastq.seqmagick_p15_cli_functions.seqmagick_quality_filter(sequence_file, output_file, *, input_qual=None, min_mean_quality=25.0, min_length=200, max_length=1000, ambiguous_action=None, max_ambiguous=None, pct_ambiguous=None, primer=None, barcode_file=None, barcode_header=False, sample_map_out=None, quality_window=0, quality_window_mean_qual=None, quality_window_prop=1.0, details_out=None, seqmagick_path=None)`.
- **Input:** A FASTQ (or FASTA plus `input_qual`) file, a new output file, thresholds and optional `primer`, `barcode_file`.
- **Output:** `output`, `output_record_count`, `report` (one row per filter with passed/failed counts).
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.fasta_fastq.seqmagick_p15_cli_functions import seqmagick_quality_filter

result = seqmagick_quality_filter("reads.fq", str(out / "q.fq"), min_length=5, max_length=100, ambiguous_action="truncate")
print(result["output_record_count"])
for row in result["report"]:
    print(row["name"], row["failed"])
primed = seqmagick_quality_filter("reads.fq", str(out / "p.fq"), min_length=5, max_length=100, primer="CAGATG")
print(primed["output_record_count"], (out / "p.fq").read_text().splitlines()[1])
```

- **Example output:**

```text
3
Quality Score [min_mean: 25.0] 0
Maximum Length [100] 0
Minimum Length [5] 0
Ambiguous Base [truncate] 0
2 ACGTACGTACGT
```

## 5) `dnakit.fasta_fastq.seqmagick_p15_cli_functions.seqmagick_info()`

- **Purpose:** Summarise sequence files: record count, length range, mean length and alignment detection.
- **Method:** Runs `seqmagick info`; a file is an alignment when all sequences have the same length.
- **API:** `dnakit.fasta_fastq.seqmagick_p15_cli_functions.seqmagick_info(sequence_files, output_file, *, input_format=None, table_format='tab', threads=1, seqmagick_path=None)`.
- **Input:** A list of files, a new output file and `table_format` (`tab`, `csv`, `align`).
- **Output:** `files`: one row per file with `alignment`, `min_len`, `max_len`, `avg_len`, `num_seqs`.
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.fasta_fastq.seqmagick_p15_cli_functions import seqmagick_info

result = seqmagick_info(["aln.fa", "reads.fq"], str(out / "info.tsv"))
for row in result["files"]:
    print(Path(row["name"]).name, row["alignment"], row["min_len"], row["max_len"], row["num_seqs"])
```

- **Example output:**

```text
aln.fa TRUE 10 10 4
reads.fq TRUE 18 18 3
```

## 6) `dnakit.fasta_fastq.seqkit_p15_cli_functions.seqkit_faidx()`

- **Purpose:** Build the FASTA index (`.fai`) and optionally extract regions.
- **Method:** Runs `seqkit faidx`. Regions are 1-based and may use negative positions; `full_head` writes `.seqkit.fai` and prints full headers.
- **API:** `dnakit.fasta_fastq.seqkit_p15_cli_functions.seqkit_faidx(fasta_path, *, regions=None, output_path=None, full_head=False, use_regexp=False, ignore_case=False, update_index=True, seqkit_path=None)`.
- **Input:** A FASTA file; optional `regions` with a new `output_path`.
- **Output:** `index_file`, `sequence_count`, `index` rows (name, length, offset, line_bases, line_width), `extracted_records`.
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
import shutil
from dnakit.fasta_fastq.seqkit_p15_cli_functions import seqkit_faidx

fasta = shutil.copy("x.fa", out / "x.fa")
result = seqkit_faidx(str(fasta))
print(result["sequence_count"], result["index"][1])
region = seqkit_faidx(str(fasta), regions=["a:2-4", "b"], output_path=str(out / "regions.fa"))
print(region["extracted_records"], (out / "regions.fa").read_text().splitlines()[1])
```

- **Example output:**

```text
3 {'name': 'b', 'length': 10, 'offset': 30, 'line_bases': 10, 'line_width': 11}
2 CGT
```

## 7) `dnakit.fasta_fastq.seqkit_p15_cli_functions.seqkit_sanitize_fastx()`

- **Purpose:** Repair a broken single-line FASTQ or FASTA file by skipping malformed records.
- **Method:** Runs `seqkit sana`, which keeps parsing after an error and drops records whose sequence and quality lengths differ.
- **API:** `dnakit.fasta_fastq.seqkit_p15_cli_functions.seqkit_sanitize_fastx(input_path, output_path, *, format='fastq', input_format=None, output_format=None, allow_gaps=False, quality_ascii_base=33, seqkit_path=None)`.
- **Input:** A damaged file and a new output path; `format`, `input_format`, `output_format`.
- **Output:** `passed_records`, `discarded_lines`, `output_record_count`.
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.fasta_fastq.seqkit_p15_cli_functions import seqkit_sanitize_fastx

result = seqkit_sanitize_fastx("broken.fq", str(out / "clean.fq"))
print(result["passed_records"], result["discarded_lines"], result["output_record_count"])
```

- **Example output:**

```text
2 5 2
```

## 8) `dnakit.fasta_fastq.seqkit_p15_cli_functions.seqkit_stream_directory()`

- **Purpose:** Concatenate and stream the FASTX files of a directory as they appear.
- **Method:** Runs `seqkit scat`. `find_only` joins the files that exist and stops; otherwise `time_limit` is required so the call ends.
- **API:** `dnakit.fasta_fastq.seqkit_p15_cli_functions.seqkit_stream_directory(directory, output_path, *, format='fastq', find_only=False, time_limit=None, gz_only=False, file_pattern=None, minimum_size_increase_kb=5, quality_ascii_base=33, seqkit_path=None)`.
- **Input:** A directory, a new output path, `format`, `find_only` or `time_limit`.
- **Output:** `passed_records`, `output_record_count`.
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
import shutil
from dnakit.fasta_fastq.seqkit_p15_cli_functions import seqkit_stream_directory

run = out / "run"
run.mkdir()
shutil.copy("r1.fq", run / "one.fq")
shutil.copy("r2.fq", run / "two.fastq")
result = seqkit_stream_directory(str(run), str(out / "all.fq"), find_only=True)
print(result["passed_records"], result["output_record_count"])
```

- **Example output:**

```text
4 4
```

## 9) `dnakit.fasta_fastq.seqkit_p15_cli_functions.seqkit_common_records()`

- **Purpose:** Records shared by two or more files, matched by ID, full name or sequence.
- **Method:** Runs `seqkit common`. Sequence matching compares both strands unless `only_positive_strand`.
- **API:** `dnakit.fasta_fastq.seqkit_p15_cli_functions.seqkit_common_records(input_paths, output_path, *, by='id', ignore_case=False, only_positive_strand=False, check_embedded_seqs=False, threads=1, seqkit_path=None)`.
- **Input:** At least two files and a new output path; `by` is `id`, `name` or `sequence`.
- **Output:** `shared_record_count`, `key`.
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.fasta_fastq.seqkit_p15_cli_functions import seqkit_common_records

for by in ("id", "name", "sequence"):
    result = seqkit_common_records(["x.fa", "other.fa"], str(out / f"{by}.fa"), by=by)
    print(by, result["shared_record_count"])
result = seqkit_common_records(["x.fa", "other.fa"], str(out / "pos.fa"), by="sequence", only_positive_strand=True)
print("positive strand", result["shared_record_count"])
```

- **Example output:**

```text
id 2
name 1
sequence 2
positive strand 1
```

## 10) `dnakit.fasta_fastq.seqkit_p15_cli_functions.seqkit_replace_header()`

- **Purpose:** Replace sequence names with a regular expression.
- **Method:** Runs `seqkit replace` on names. The replacement may use `$1`, `{nr}` (record number, `number_width` zero-pads), `{fn}` and `{kv}` (value from a key-value file).
- **API:** `dnakit.fasta_fastq.seqkit_p15_cli_functions.seqkit_replace_header(input_path, output_path, pattern, replacement, *, ignore_case=False, key_value_file=None, key_capture_index=1, keep_key=False, keep_untouched=False, missing_key_replacement=None, number_width=1, only_records_matching=None, only_records_matching_regexp=False, threads=1, seqkit_path=None)`.
- **Input:** A sequence file, a new output path, `pattern` and `replacement`.
- **Output:** `output`, `output_record_count`.
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.fasta_fastq.seqkit_p15_cli_functions import seqkit_replace_header

seqkit_replace_header("x.fa", str(out / "renamed.fa"), "^(\\w)", "seq{nr}_$1", number_width=2)
print([l for l in (out / "renamed.fa").read_text().splitlines() if l.startswith(">")])
```

- **Example output:**

```text
['>seq01_a first', '>seq02_b second', '>seq03_c third']
```

## 11) `dnakit.fasta_fastq.seqkit_p15_cli_functions.seqkit_bam_dump_records()`

- **Purpose:** Dump per-record alignment features of a BAM file as a table.
- **Method:** Runs `seqkit bam -T Dump`. Features that need the `NM` tag (`Acc`, `Match`, `Mismatch`) require it in the BAM.
- **API:** `dnakit.fasta_fastq.seqkit_p15_cli_functions.seqkit_bam_dump_records(bam_path, tsv_path, *, fields=None, seqkit_path=None)`.
- **Input:** A BAM file, a new TSV path and an optional list of `fields`.
- **Output:** `rows`, `row_count`, `fields`.
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
import pysam
from dnakit.fasta_fastq.seqkit_p15_cli_functions import seqkit_bam_dump_records

header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 100}]}
with pysam.AlignmentFile(str(out / "x.bam"), "wb", header=header) as handle:
    for index, (start, nm) in enumerate(((5, 0), (10, 1))):
        read = pysam.AlignedSegment()
        read.query_name, read.query_sequence, read.flag = f"read{index}", "ACGT" * 5, 0
        read.reference_id, read.reference_start, read.mapping_quality, read.cigarstring = 0, start, 60, "20M"
        read.query_qualities = pysam.qualitystring_to_array("I" * 20)
        read.set_tag("NM", nm)
        handle.write(read)
result = seqkit_bam_dump_records(str(out / "x.bam"), str(out / "dump.tsv"), fields=["Read", "Pos", "Acc", "MapQual"])
for row in result["rows"]:
    print(row)
```

- **Example output:**

```text
{'Read': 'read0', 'Pos': 5, 'Acc': 100.0, 'MapQual': 60}
{'Read': 'read1', 'Pos': 10, 'Acc': 95.0, 'MapQual': 60}
```

## 12) `dnakit.fasta_fastq.seqkit_p15_cli_functions.seqkit_bam_filter_by_alignment_context()`

- **Purpose:** Keep BAM records by the reference sequence at their start and end.
- **Method:** Runs `seqkit bam -T AlnContext`: the reference windows `[pos+left_shift, pos+right_shift)` around both ends are tested with Go regular expressions; `stranded` reverse-complements reverse-strand reads and `invert` flips the test.
- **API:** `dnakit.fasta_fastq.seqkit_p15_cli_functions.seqkit_bam_filter_by_alignment_context(bam_path, reference_fasta, output_bam, *, regex_start=None, regex_end=None, left_shift=-10, right_shift=10, stranded=False, invert=False, tsv_path=None, seqkit_path=None)`.
- **Input:** A BAM, an indexed reference FASTA, a new output BAM and `regex_start` and/or `regex_end`.
- **Output:** `output_bam`, `context_rows` (examined windows and match flags).
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
import pysam
from dnakit.fasta_fastq.seqkit_p15_cli_functions import seqkit_bam_filter_by_alignment_context

sequence = "ACGTTGCAAGGCTTAACCGGATATCGCGATTTAAACCCGGGTTTACGTACGT"
reference = out / "ref.fa"
reference.write_text(f">chr1\n{sequence}\n")
pysam.faidx(str(reference))
header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": len(sequence)}]}
with pysam.AlignmentFile(str(out / "x.bam"), "wb", header=header) as handle:
    for index, start in enumerate((5, 10)):
        read = pysam.AlignedSegment()
        read.query_name, read.query_sequence, read.flag = f"read{index}", sequence[start:start + 20], 0
        read.reference_id, read.reference_start, read.mapping_quality, read.cigarstring = 0, start, 60, "20M"
        read.query_qualities = pysam.qualitystring_to_array("I" * 20)
        handle.write(read)
result = seqkit_bam_filter_by_alignment_context(str(out / "x.bam"), str(reference), str(out / "kept.bam"),
                                                regex_start="^GC", left_shift=0, right_shift=3, tsv_path=str(out / "ctx.tsv"))
for row in result["context_rows"]:
    print(row["Read"], row["StartSeq"], row["StartMatch"])
print(pysam.AlignmentFile(str(out / "kept.bam")).count(until_eof=True))
```

- **Example output:**

```text
read0 GCA -1
read1 GCT -1
2
```

## 13) `dnakit.fasta_fastq.seqkit_p15_cli_functions.seqkit_bam_accuracy_stats()`

- **Purpose:** Mean alignment accuracy of a BAM file, plain and weighted by alignment length.
- **Method:** Runs `seqkit bam -T AccStats`; accuracy of a record is `(1 - NM / alignment length) * 100`.
- **API:** `dnakit.fasta_fastq.seqkit_p15_cli_functions.seqkit_bam_accuracy_stats(bam_path, *, tsv_path, seqkit_path=None)`.
- **Input:** A BAM carrying the `NM` tag and a new TSV path.
- **Output:** `accuracy_mean`, `accuracy_weighted_mean`.
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
import pysam
from dnakit.fasta_fastq.seqkit_p15_cli_functions import seqkit_bam_accuracy_stats

header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 100}]}
with pysam.AlignmentFile(str(out / "x.bam"), "wb", header=header) as handle:
    for index, (start, nm) in enumerate(((5, 0), (10, 1))):
        read = pysam.AlignedSegment()
        read.query_name, read.query_sequence, read.flag = f"read{index}", "ACGT" * 5, 0
        read.reference_id, read.reference_start, read.mapping_quality, read.cigarstring = 0, start, 60, "20M"
        read.query_qualities = pysam.qualitystring_to_array("I" * 20)
        read.set_tag("NM", nm)
        handle.write(read)
result = seqkit_bam_accuracy_stats(str(out / "x.bam"), tsv_path=str(out / "acc.tsv"))
print(result["accuracy_mean"], result["accuracy_weighted_mean"])
```

- **Example output:**

```text
97.5 97.5
```

## 14) `dnakit.sequencing_qc.afterqc_p15_cli_functions.afterqc_overlap_analysis()`

- **Purpose:** Trim adapters, correct mismatched bases and store the overlap of read pairs by overlap analysis.
- **Method:** Runs `after.py` with its overlap analysis on; `analyse_overlap=False` is `--no_overlap`, `correct_mismatches=False` is `--no_correction`, `mask_mismatch` is `--mask_mismatch`, `store_overlap` is `--store_overlap on`.
- **API:** `dnakit.sequencing_qc.afterqc_p15_cli_functions.afterqc_overlap_analysis(read1, read2, output_directory, *, correct_mismatches=True, mask_mismatch=False, analyse_overlap=True, store_overlap=False, trim_front=-1, trim_tail=-1, qualified_quality_phred=15, unqualified_base_limit=60, poly_size_limit=35, n_base_limit=5, seq_len_req=35, gzip_output=False, afterqc_path=None)`.
- **Input:** Two FASTQ files of one library and a new output folder.
- **Output:** `good_files`, `bad_files`, `qc_files`, `overlap_files`.
- **Example call (needs the external program; not executed in the documentation build):**

```python
from dnakit.sequencing_qc.afterqc_p15_cli_functions import afterqc_overlap_analysis

result = afterqc_overlap_analysis("sample_R1.fq", "sample_R2.fq", "afterqc_out", store_overlap=True)
print(result["good_files"])
```

## 15) `dnakit.sequencing_qc.afterqc_p15_cli_functions.afterqc_debubble()`

- **Purpose:** Remove reads from sequencer bubbles before filtering.
- **Method:** Runs `after.py --debubble` on a directory of Illumina FASTQ files; reads inside circular regions with abnormal bad-read density are moved to `bad`.
- **API:** `dnakit.sequencing_qc.afterqc_p15_cli_functions.afterqc_debubble(input_directory, output_directory, *, draw=True, trim_front=-1, trim_tail=-1, qualified_quality_phred=15, unqualified_base_limit=60, poly_size_limit=35, n_base_limit=5, seq_len_req=35, gzip_output=False, afterqc_path=None)`.
- **Input:** An input directory and a new output folder.
- **Output:** `good_files`, `bad_files`, `debubble_files`, `circles`.
- **Example call (needs the external program; not executed in the documentation build):**

```python
from dnakit.sequencing_qc.afterqc_p15_cli_functions import afterqc_debubble

result = afterqc_debubble("lane_fastq_dir", "afterqc_debubble_out", draw=False)
print(result["circles"])
```

## 16) `dnakit.sequencing_qc.afterqc_p15_cli_functions.afterqc_move_barcodes_to_read_names()`

- **Purpose:** Move single-molecule barcodes from read starts into read names.
- **Method:** Runs `after.py --barcode on`: the first `barcode_length` bases are the barcode, the next bases must equal `barcode_verify`. AfterQC treats a file as barcoded when its name contains `barcode_flag`.
- **API:** `dnakit.sequencing_qc.afterqc_p15_cli_functions.afterqc_move_barcodes_to_read_names(read1, output_directory, *, barcode_length=12, barcode_verify='CAGTA', barcode_flag='barcode', read2=None, qualified_quality_phred=15, unqualified_base_limit=60, poly_size_limit=35, n_base_limit=5, seq_len_req=35, gzip_output=False, afterqc_path=None)`.
- **Input:** A FASTQ file whose name contains `barcode` and a new output folder.
- **Output:** `good_files`, `bad_files`, `qc_files`.
- **Example call (needs the external program; not executed in the documentation build):**

```python
from dnakit.sequencing_qc.afterqc_p15_cli_functions import afterqc_move_barcodes_to_read_names

result = afterqc_move_barcodes_to_read_names("run_barcode_R1.fq", "afterqc_barcode_out", barcode_length=12, barcode_verify="CAGTA")
print(result["good_files"])
```

## 17) `dnakit.sequencing_qc.trimmomatic_p15_cli_functions.trimmomatic_process_reads()`

- **Purpose:** Trim and filter single-end or paired-end reads with Trimmomatic steps (sliding window, MAXINFO, AVGQUAL, BASECOUNT, MINLEN, MAXLEN, TOPHRED33/64, CROP, HEADCROP, TAILCROP, LEADING, TRAILING).
- **Method:** Builds the step list in the order Trimmomatic's manual recommends and runs `SE` or `PE` with `-summary`. Trimmomatic is a Java program: pass `trimmomatic_jar` (the `.jar` file), or put the `trimmomatic` launcher on PATH.
- **API:** `dnakit.sequencing_qc.trimmomatic_p15_cli_functions.trimmomatic_process_reads(read1, output_directory, *, read2=None, input_phred=None, convert_quality=None, crop=None, head_crop=None, tail_crop=None, leading=None, trailing=None, sliding_window=None, max_info=None, average_quality=None, base_count=None, min_length=None, max_length=None, threads=1, trimmomatic_jar=None)`.
- **Input:** `read1` (and `read2`), a new output folder and at least one step parameter; `input_phred` sets the input encoding.
- **Output:** `steps`, `survival` (input, surviving and dropped reads), `read_outputs`.
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
import os
from dnakit.sequencing_qc.trimmomatic_p15_cli_functions import trimmomatic_process_reads

jar = os.environ.get("TRIMMOMATIC_JAR")
single = trimmomatic_process_reads("r1.fq", str(out / "se"), sliding_window=[4, 20], min_length=25, max_length=100,
                                   average_quality=20, base_count=["N", 0, 2], input_phred=33, trimmomatic_jar=jar)
print(single["steps"])
print(single["survival"]["input_reads"], single["survival"]["surviving_reads"])
paired = trimmomatic_process_reads("r1.fq", str(out / "pe"), read2="r2.fq", convert_quality=64, input_phred=33, trimmomatic_jar=jar)
print(paired["survival"]["both_surviving_reads"])
```

- **Example output:**

```text
['SLIDINGWINDOW:4:20', 'AVGQUAL:20', 'BASECOUNT:N:0:2', 'MINLEN:25', 'MAXLEN:100']
2 2
2
```

## 18) `dnakit.sequencing_qc.trimmomatic_p15_cli_functions.trimmomatic_pair_reads()`

- **Purpose:** Pair FASTQ mates by read name and split off reads without a mate.
- **Method:** Runs Trimmomatic's `Pairomatic` class from the jar. Names are compared whole unless `delimiter` is given (`/` for `name/1` and `name/2`); reads present in both files must be in the same order.
- **API:** `dnakit.sequencing_qc.trimmomatic_p15_cli_functions.trimmomatic_pair_reads(read1, read2, output_directory, *, delimiter=None, trimmomatic_jar=None)`.
- **Input:** Two FASTQ files, a new output folder and `trimmomatic_jar` (required).
- **Output:** `read_counts` per output file and `read_outputs`.
- **Example code:**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
import os
from dnakit.sequencing_qc.trimmomatic_p15_cli_functions import trimmomatic_pair_reads

result = trimmomatic_pair_reads("m1.fq", "m2.fq", str(out / "paired"), delimiter="/", trimmomatic_jar=os.environ["TRIMMOMATIC_JAR"])
print(result["read_counts"])
```

- **Example output:**

```text
{'R1.paired.fastq': 2, 'R1.unpaired.fastq': 1, 'R2.paired.fastq': 2, 'R2.unpaired.fastq': 0}
```
