# BEDOPS Starch archives and convert2bed

These functions in `dnakit.annotation` (module `bedops_p05_cli_functions`) wrap the Starch archive programs of [BEDOPS](https://github.com/bedops/bedops) 2.4.42 (`starch`, `unstarch`, `starchcat`, `starchstrip`) and `convert2bed`, which turns BAM, SAM, GFF3, GTF, GVF, PSL, RepeatMasker and VCF/WIG files into sorted BED or Starch. Starch is BEDOPS' compressed BED format with per-chromosome streams and metadata. `starch-diff` is not wrapped: in the 2.4.42 release it is a Python 2 script (`unicode`, `iteritems`) that fails on Python 3. Programs write to standard output; the functions stream that to an explicit new file and never overwrite.

BEDOPS is an external program: install it (for example `conda install -c bioconda bedops`, or the release archive from https://github.com/bedops/bedops/releases), then put its `bin` folder on PATH or pass `bedops_bin_dir`. If a program cannot be found, `BackendUnavailableError` is raised with the install hint; if it exits with an error, `BackendExecutionError` (code `BEDOPS_EXECUTION_FAILED`) is raised with the tail of its log and no partial output is left behind. BEDOPS is released under the GPL-2.0 licence and is not bundled with DNAKit.

All results share `status` (`completed`), `backend` (`bedops`), `operation`, `program`, `version` (read from the program), `command`, `outputs` (paths of the files written), `log_tail` and `provenance` (tool, version, tested version 2.4.42, executable, citation); the sections below list only the extra fields.

## Example data

The examples read these files from the working directory and write their outputs to a temporary directory.

**Example file:** `ref.bed`

```text
chr1	100	200	r1
chr1	300	400	r2
chr1	1000	1100	r3
chr2	10	60	r4
```

**Example file:** `map.bed`

```text
chr1	150	250	m1	2	+
chr1	180	350	m2	4	-
chr1	320	330	m3	6	+
chr2	0	20	m4	10	+
```

**Example file:** `variants.vcf`

```text
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	101	rs1	A	G	50	PASS	.
chr1	51	rs2	C	T	60	PASS	.
```

**Example file:** `genes.gtf`

```text
chr1	demo	exon	101	200	.	+	.	gene_id "g1"; transcript_id "t1";
chr1	demo	exon	301	400	.	-	.	gene_id "g2"; transcript_id "t2";
```

**Example file:** `repeats.out`

```text
   SW  perc perc perc  query      position in query           matching       repeat              position in  repeat
score  div. del. ins.  sequence    begin     end    (left)    repeat         class/family         begin  end (left)   ID

  463  10.2  0.0  1.2  chr1        101     200  (900)  +  L1PA2         LINE/L1                 1   100  (10)    1
```

## 1) `dnakit.annotation.starch_compress()`

- **Purpose:** Compress a sorted BED file into a BEDOPS Starch archive.
- **Method:** Run `starch [--bzip2|--gzip] [--note=...] [--omit-signature] [--header] [<unique-tag>] <bed>` and stream the archive to `output_path`; the element count is then read with `unstarch --elements`.
- **API:** `starch_compress(input_path: 'str', output_path: 'str', *, compression: 'str' = 'bzip2', note: 'str | None' = None, omit_signature: 'bool' = False, header: 'bool' = False, unique_tag: 'str | None' = None, bedops_bin_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the BEDOPS program `starch` (2.4.42 is the tested version) on PATH, or `bedops_bin_dir` pointing to the BEDOPS `bin` folder; no Python package is required. The example output below was produced with BEDOPS 2.4.42 installed. A missing program raises `BackendUnavailableError` with an install hint.
- **Input:** `input_path` is an existing BED file sorted as `sort-bed` sorts it. `output_path` must not exist yet. `compression` is `bzip2` (default, smaller) or `gzip` (extracts faster). `note` is a single-line string stored in the metadata. `omit_signature` skips the per-chromosome signature. `header` accepts UCSC track, SAM or VCF headers. `unique_tag` is an optional identifier.
- **Output:** The common result dictionary plus `element_count` and `size_bytes`. `outputs` holds `archive`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.annotation.bedops_p05_cli_functions import starch_compress, unstarch_to_bed

with tempfile.TemporaryDirectory() as tmp:
    archive = str(Path(tmp) / "ref.starch")
    result = starch_compress("ref.bed", archive, compression="gzip", note="documentation example", verbose=False)
    print(result["status"], result["element_count"], list(result["outputs"]))
    unstarch_to_bed(archive, str(Path(tmp) / "back.bed"), verbose=False)
    print("round trip equals input:", (Path(tmp) / "back.bed").read_text() == Path("ref.bed").read_text())
```

- **Example output:**

```text
completed 4 ['archive']
round trip equals input: True
```

## 2) `dnakit.annotation.unstarch_to_bed()`

- **Purpose:** Decompress a Starch archive to BED, optionally one chromosome only.
- **Method:** Run `unstarch [<chromosome>] <archive>` and stream the BED rows to `output_path`.
- **API:** `unstarch_to_bed(input_path: 'str', output_path: 'str', *, chromosome: 'str | None' = None, bedops_bin_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the BEDOPS program `unstarch` (2.4.42 is the tested version) on PATH, or `bedops_bin_dir` pointing to the BEDOPS `bin` folder; no Python package is required. The example output below was produced with BEDOPS 2.4.42 installed. A missing program raises `BackendUnavailableError` with an install hint.
- **Input:** `input_path` is an existing Starch archive. `output_path` must not exist yet. `chromosome` is optional.
- **Output:** The common result dictionary plus `line_count` and `preview`. `outputs` holds `output`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.annotation.bedops_p05_cli_functions import starch_compress, unstarch_to_bed

with tempfile.TemporaryDirectory() as tmp:
    archive = str(Path(tmp) / "ref.starch")
    starch_compress("ref.bed", archive, verbose=False)
    whole = unstarch_to_bed(archive, str(Path(tmp) / "all.bed"), verbose=False)
    one = unstarch_to_bed(archive, str(Path(tmp) / "chr1.bed"), chromosome="chr1", verbose=False)
    print(whole["line_count"], one["line_count"], one["preview"][0].replace("\t", " "))
```

- **Example output:**

```text
4 3 chr1 100 200 r1
```

## 3) `dnakit.annotation.unstarch_summary()`

- **Purpose:** Read one metadata value of a Starch archive without decompressing it.
- **Method:** Run `unstarch [<chromosome>] --<statistic> <archive>` and parse the output. `statistic` is `elements` (element count), `bases` (total bases), `bases-uniq` (distinct bases), `has-duplicate` and `has-nested` (returned as bool), `list-chromosomes` (returned as a list), `elements-max-string-length` or `list-json` (the archive metadata as a dictionary). The function covers the work-list rows element count, total bases, unique bases, detect duplicate elements, detect nested elements and list chromosomes of a Starch archive.
- **API:** `unstarch_summary(input_path: 'str', *, statistic: 'str' = 'elements', chromosome: 'str | None' = None, bedops_bin_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the BEDOPS program `unstarch` (2.4.42 is the tested version) on PATH, or `bedops_bin_dir` pointing to the BEDOPS `bin` folder; no Python package is required. The example output below was produced with BEDOPS 2.4.42 installed. A missing program raises `BackendUnavailableError` with an install hint.
- **Input:** `input_path` is an existing Starch archive. `statistic` is one of the values above. `chromosome` (optional) restricts the value to that chromosome.
- **Output:** The common result dictionary plus `statistic`, `chromosome` and `value`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.annotation.bedops_p05_cli_functions import starch_compress, unstarch_summary

with tempfile.TemporaryDirectory() as tmp:
    archive = str(Path(tmp) / "ref.starch")
    starch_compress("ref.bed", archive, verbose=False)
    for statistic in ("elements", "bases", "bases-uniq", "has-duplicate", "has-nested", "list-chromosomes"):
        print(statistic, unstarch_summary(archive, statistic=statistic, verbose=False)["value"])
    print("chr1 elements", unstarch_summary(archive, statistic="elements", chromosome="chr1", verbose=False)["value"])
```

- **Example output:**

```text
elements 4
bases 350
bases-uniq 350
has-duplicate False
has-nested False
list-chromosomes ['chr1', 'chr2']
chr1 elements 3
```

## 4) `dnakit.annotation.starchcat_archives()`

- **Purpose:** Merge Starch archives into one multiset-union archive, or re-compress one archive.
- **Method:** Run `starchcat [--bzip2|--gzip] [--note=...] [--omit-signature] <archives>` and stream the new archive to `output_path`; the element count is read with `unstarch --elements`. With one input it only re-compresses the archive or updates its note.
- **API:** `starchcat_archives(input_paths: 'Sequence[str]', output_path: 'str', *, compression: 'str' = 'bzip2', note: 'str | None' = None, omit_signature: 'bool' = False, bedops_bin_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the BEDOPS program `starchcat` (2.4.42 is the tested version) on PATH, or `bedops_bin_dir` pointing to the BEDOPS `bin` folder; no Python package is required. The example output below was produced with BEDOPS 2.4.42 installed. A missing program raises `BackendUnavailableError` with an install hint.
- **Input:** `input_paths` is a non-empty list of existing sorted, header-less Starch archives. `output_path` must not exist yet. `compression`, `note` and `omit_signature` are as in `starch_compress`.
- **Output:** The common result dictionary plus `element_count` and `size_bytes`. `outputs` holds `archive`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.annotation.bedops_p05_cli_functions import starch_compress, starchcat_archives, unstarch_summary

with tempfile.TemporaryDirectory() as tmp:
    first, second, merged = (str(Path(tmp) / name) for name in ("ref.starch", "map.starch", "merged.starch"))
    starch_compress("ref.bed", first, verbose=False)
    starch_compress("map.bed", second, verbose=False)
    result = starchcat_archives([first, second], merged, verbose=False)
    print(result["element_count"], unstarch_summary(merged, statistic="list-chromosomes", verbose=False)["value"])
```

- **Example output:**

```text
8 ['chr1', 'chr2']
```

## 5) `dnakit.annotation.starchstrip_chromosomes()`

- **Purpose:** Keep or drop chromosomes of a Starch archive without recompressing it.
- **Method:** Run `starchstrip --include <list> <archive>` and stream the result to `output_path`; chromosomes missing from the archive are ignored. BEDOPS 2.4.42 `starchstrip --exclude` writes an archive that `unstarch` cannot read when the first chromosome is removed, so for `mode="exclude"` the function lists the archive's chromosomes and passes the complementary `--include` list, which gives the intended archive. Excluding every chromosome is refused.
- **API:** `starchstrip_chromosomes(input_path: 'str', output_path: 'str', chromosomes: 'Sequence[str]', *, mode: 'str' = 'include', bedops_bin_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the BEDOPS program `starchstrip` (2.4.42 is the tested version) on PATH, or `bedops_bin_dir` pointing to the BEDOPS `bin` folder; no Python package is required. The example output below was produced with BEDOPS 2.4.42 installed. A missing program raises `BackendUnavailableError` with an install hint.
- **Input:** `input_path` is an existing Starch archive. `output_path` must not exist yet. `chromosomes` is a non-empty list of names without commas. `mode` is `include` or `exclude`.
- **Output:** The common result dictionary plus `mode`, `strategy`, `chromosomes_kept` and `element_count`. `outputs` holds `archive`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.annotation.bedops_p05_cli_functions import starch_compress, starchstrip_chromosomes, unstarch_to_bed

with tempfile.TemporaryDirectory() as tmp:
    archive = str(Path(tmp) / "ref.starch")
    starch_compress("ref.bed", archive, verbose=False)
    result = starchstrip_chromosomes(archive, str(Path(tmp) / "no_chr1.starch"), ["chr1"], mode="exclude", verbose=False)
    print(result["chromosomes_kept"], result["element_count"])
    rows = unstarch_to_bed(result["outputs"]["archive"], str(Path(tmp) / "no_chr1.bed"), verbose=False)
    print(rows["preview"][0].replace("\t", " "))
```

- **Example output:**

```text
['chr2'] 1
chr2 10 60 r4
```

## 6) `dnakit.annotation.convert_to_bed()`

- **Purpose:** Convert BAM, SAM, GFF3, GTF, GVF, PSL, RepeatMasker or VCF/WIG files to sorted BED or to a Starch archive.
- **Method:** Run `convert2bed --input=<format> --output=<bed|starch> [format options] < input` with the input file attached to standard input, and stream the result to `output_path`. `input_format` is `bam`, `sam`, `gff` (GFF3), `gtf`, `gvf`, `psl`, `rmsk` (RepeatMasker .out), `vcf` or `wig`; `output_format` is `bed` or `starch`. One function covers `bam2bed`, `sam2bed`, `gff2bed`, `gtf2bed`, `gvf2bed`, `psl2bed`, `rmsk2bed`, `vcf2bed`, `wig2bed` and the matching `*2starch` programs, including the cluster wrappers `bam2starch_sge`, `bam2starch_slurm` and `bam2starch_gnuParallel`, which only distribute the same conversion. `format_options` are boolean switches valid for the input format (`keep-header`; `split`, `all-reads`, `reduced`, `split-with-deletions` for sam/bam; `split` for psl; `snvs`, `insertions`, `deletions`, `do-not-split` for vcf; `zero-indexed` for wig); `attribute_key` selects the GENCODE attribute copied to the BED name column for `gff` and `gtf`. `sam` and `bam` input needs `samtools` on PATH, otherwise `BackendUnavailableError` is raised. The function covers the work-list rows converting BAM, SAM, GFF3, GTF, GVF, PSL, RepeatMasker OUT, VCF and WIG to Starch and SAM, GVF, PSL, RepeatMasker OUT, VCF and WIG to BED.
- **API:** `convert_to_bed(input_path: 'str', output_path: 'str', *, input_format: 'str', output_format: 'str' = 'bed', format_options: 'Sequence[str] | None' = None, attribute_key: 'str | None' = None, sort_output: 'bool' = True, max_memory: 'str | None' = None, temp_dir: 'str | None' = None, starch_compression: 'str' = 'bzip2', starch_note: 'str | None' = None, bedops_bin_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** needs the BEDOPS program `convert2bed` (2.4.42 is the tested version) on PATH, or `bedops_bin_dir` pointing to the BEDOPS `bin` folder; no Python package is required. The example output below was produced with BEDOPS 2.4.42 installed. A missing program raises `BackendUnavailableError` with an install hint.
- **Input:** `input_path` is an existing file in `input_format`. `output_path` must not exist yet. `format_options` is a list of option names allowed for that format; `attribute_key` is one of the GENCODE keys listed by `convert2bed --help-gff` / `--help-gtf` and only applies to gff and gtf. `sort_output=False` adds `--do-not-sort` (BED only). `max_memory` (`8G`) and `temp_dir` set `--max-mem` and `--sort-tmpdir`. `starch_compression` (`bzip2` or `gzip`) and `starch_note` apply to `output_format="starch"` only.
- **Output:** The common result dictionary plus `input_format`, `output_format`; for BED output `line_count` and `preview` (`outputs` holds `output`), for Starch output `element_count` and `size_bytes` (`outputs` holds `archive`).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.annotation.bedops_p05_cli_functions import convert_to_bed

with tempfile.TemporaryDirectory() as tmp:
    for name, fmt in (("variants.vcf", "vcf"), ("genes.gtf", "gtf"), ("repeats.out", "rmsk")):
        result = convert_to_bed(name, str(Path(tmp) / f"{fmt}.bed"), input_format=fmt, verbose=False)
        print(fmt, "->", result["preview"][0].replace("\t", " ")[:60])
    snvs = convert_to_bed("variants.vcf", str(Path(tmp) / "v.starch"), input_format="vcf", output_format="starch", verbose=False)
    print("vcf -> starch elements:", snvs["element_count"])
    named = convert_to_bed(
        "genes.gtf", str(Path(tmp) / "named.bed"), input_format="gtf", attribute_key="transcript_id", sort_output=False, verbose=False
    )
    print([line.split("\t")[3] for line in named["preview"]])
```

- **Example output:**

```text
vcf -> chr1 50 51 rs2 60 C T PASS .
gtf -> chr1 100 200 g1 . + demo exon . gene_id "g1"; transcript_id 
rmsk -> chr1 100 200 L1PA2 463 + 10.2 0.0 1.2 (900) LINE/L1 1 100 (1
vcf -> starch elements: 2
['t1', 't2']
```
