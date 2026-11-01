# GimmeMotifs: motif scanning, clustering, background sets and de novo discovery

This group of 40 functions in `dnakit.sequence` is a set of adapters over [GimmeMotifs](https://github.com/vanheeringen-lab/gimmemotifs), a suite for transcription factor motif analysis. They scan sequences with position frequency matrices (PFM), match, cluster and draw motifs, build background sequence sets (random, genomic and GC-matched), combine peak summits, count reads in peaks, compare motif frequencies between sequence sets, write ROC and enrichment reports, run the de novo motif tools, and build motif-to-factor tables from ortholog groups. Some functions also use [pybedtools](https://github.com/daler/pybedtools) with the [bedtools](https://github.com/arq5x/bedtools2) program, [ete3](https://github.com/etetoolkit/ete), [gffread](https://github.com/gpertea/gffread) and [OrthoFinder](https://github.com/davidemms/OrthoFinder). All functions work on local files. None of them downloads a genome or a database, and a genome is always given as a FASTA path, never as a genomepy name. The only function that talks to the network is `gimmemotifs_jaspar_motif_info()`.

GimmeMotifs is MIT-licensed (the licence file in the installed package is the MIT licence). The functions are imported from `dnakit.sequence`, for example `from dnakit.sequence import gimmemotifs_pfm_scan`. The dependencies are optional and are imported only when a function is called: install them with `pip install "dnakit[bio-gimmemotifs]"` (GimmeMotifs 0.18 and ete3). External programs (`gimme`, `bedtools`, `gffread`, `orthofinder` and the de novo tools) are looked up when the function runs, and a missing one raises `BackendUnavailableError`. Every function except `gimmemotifs_jaspar_motif_info()` returns a dictionary with `status`, `operation`, `backend`, `version` and `outputs` (output names mapped to path strings), plus the fields named in each section; `gimmemotifs_jaspar_motif_info()` returns only `motif_id` and `info`. Output files and directories must be new: a path or directory that already exists is refused with `ConfigurationError`. The one exception is `gimmemotifs_maelstrom_html_report()`, which writes its report into the existing run folder.

## Example data

Most examples use three small motifs in `motifs.pfm`, six short sequences in `seqs.fa` (each with a `TGACAC` or `GCATAC` site in some of them) and one database motif in `db.pfm`. Examples that need more data (genomes, peaks, read files, background sets) generate it with a seeded random generator and write it to a temporary directory. Nothing is written outside that directory.

**Example file:** `motifs.pfm`

```text
>MA0001.1 AR1
0.2 0.2 0.2 0.4
0.1 0.1 0.7 0.1
0.8 0.1 0.05 0.05
0.05 0.8 0.1 0.05
0.6 0.2 0.1 0.1
0.1 0.7 0.1 0.1
>MA0002.1 GATA
0.1 0.1 0.7 0.1
0.1 0.7 0.1 0.1
0.7 0.1 0.1 0.1
0.1 0.1 0.1 0.7
0.4 0.2 0.2 0.2
0.2 0.4 0.2 0.2
>MA0001.2 AR1b
0.2 0.2 0.2 0.4
0.1 0.1 0.7 0.1
0.7 0.1 0.1 0.1
0.05 0.8 0.1 0.05
0.6 0.2 0.1 0.1
0.1 0.6 0.2 0.1
```

**Example file:** `seqs.fa`

```text
>s1
ACGTACGTAATTCCGGTGACACAGTCCATGACTTAGC
>s2
GGATCCATGACCTGACACTTAGGCATACCGATCAGTA
>s3
TTAGCGCATACGGATCCTGACACGTACTGAGCATCGA
>s4
CCGATTAGCGTACGTTAGCCATGCATGGATCAGTTCA
>s5
AGCATACTTGACACGGTACCGATAGCTTAGGCATCCA
>s6
TACCGTAGCTAGGATCCGATAGCTTACGGATCCGTAA
```

**Example file:** `db.pfm`

```text
>DB1 long
0.25 0.25 0.25 0.25
0.1 0.1 0.7 0.1
0.8 0.1 0.05 0.05
0.05 0.8 0.1 0.05
0.6 0.2 0.1 0.1
0.1 0.7 0.1 0.1
0.25 0.25 0.25 0.25
0.25 0.25 0.25 0.25
```

## 1) `dnakit.sequence.gimmemotifs_cluster_motif_report()`

- **Purpose:** Cluster motifs by similarity and write an HTML report with logos.
- **Method:** Create the output directory (with an `images` subdirectory), then call GimmeMotifs `cluster_motifs_with_report()` on the PFM file with the similarity threshold. The clustered PFM is written to `clustered_motifs.pfm` in that directory.
- **API:** `gimmemotifs_cluster_motif_report(pfm_file: 'str', *, output_dir: 'str', threshold: 'float', title: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package `gimmemotifs` (extra `dnakit[bio-gimmemotifs]`).
- **Input:** `pfm_file` is an existing PFM file. `output_dir` must not exist. `threshold` is a number in (0, 1]. `title` is an optional report title.
- **Output:** `clusters` (a list of `{"representative", "members"}` with motif ids), `cluster_count`, and `outputs` with `clustered_pfm` and, when GimmeMotifs wrote it, `report`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import gimmemotifs_cluster_motif_report

out = Path(tempfile.mkdtemp()) / "report"
result = gimmemotifs_cluster_motif_report(
    "motifs.pfm", output_dir=out, threshold=0.9, verbose=False
)
print(result["cluster_count"], sorted(result["outputs"]))
for cluster in sorted(result["clusters"], key=lambda c: c["representative"]):
    print(cluster["representative"], sorted(cluster["members"]))
```

- **Example output:**

```text
2 ['clustered_pfm', 'report']
Average_1 ['MA0001.1 AR1', 'MA0001.2 AR1b']
MA0002.1 GATA ['MA0002.1 GATA']
```

## 2) `dnakit.sequence.gimmemotifs_cluster_pfm()`

- **Purpose:** Cluster PFM motifs by similarity and write the clustered PFM and an HTML report.
- **Method:** Run the `gimme cluster` command with a fixed argument list and no shell. `gimme` creates the new output directory, writes `clustered_motifs.pfm`, `gimme.clustered.html` and one logo per cluster, and the function checks that the first two exist.
- **API:** `gimmemotifs_cluster_pfm(input_pfm: 'str', *, output_dir: 'str', threshold: 'float | None' = None, compare_reverse_complement: 'bool' = True, threads: 'int' = 1, gimme_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package `gimmemotifs` (extra `dnakit[bio-gimmemotifs]`) and its `gimme` program. The program is taken from `gimme_path`, then from PATH, then from the directory of the running Python interpreter.
- **Input:** `input_pfm` is an existing PFM file. `output_dir` must not exist. `threshold` is None (the GimmeMotifs default) or a number in (0, 1]. `compare_reverse_complement=False` passes `-s` (no reverse complement comparison). `threads` is a positive integer.
- **Output:** `outputs` with `clustered_pfm` and `report`, and `log_tail`, the last 2000 characters of the program output.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import gimmemotifs_cluster_pfm

out = Path(tempfile.mkdtemp()) / "clusters"
result = gimmemotifs_cluster_pfm("motifs.pfm", output_dir=out, threshold=0.9, verbose=False)
print(sorted(result["outputs"]))
print((out / "clustered_motifs.pfm").read_text().count(">"))
```

- **Example output:**

```text
['clustered_pfm', 'report']
2
```

## 3) `dnakit.sequence.gimmemotifs_combine_peak_summits()`

- **Purpose:** Combine MACS2 summit or narrowPeak files into one set of non-overlapping summits.
- **Method:** Call GimmeMotifs `combine_peaks()` on the files with the window width. Within each window only the summit with the highest score is kept. GimmeMotifs indexes the genome FASTA, and the index files (`.fai`, `.sizes` and `.gaps.bed`) are written next to it. The rows are returned, and optionally written as a headerless tab-separated file.
- **API:** `gimmemotifs_combine_peak_summits(peak_files: 'Sequence[str]', *, genome_fasta: 'str', window: 'int' = 200, scale_value: 'bool' = False, output_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package `gimmemotifs` (extra `dnakit[bio-gimmemotifs]`).
- **Input:** `peak_files` is a non-empty list of existing summit or narrowPeak files. `genome_fasta` is an existing FASTA file. `window` is a positive integer. `scale_value` rescales the scores per file. `output_path`, when given, must not exist.
- **Output:** `columns`, `rows` (one dictionary per kept summit), `row_count`, and `outputs` with `table` when `output_path` is given. Values in `rows` are strings for the text columns read from the files.
- **Example code:**

```python
import random
import tempfile
from pathlib import Path

from dnakit.sequence import gimmemotifs_combine_peak_summits

t = Path(tempfile.mkdtemp())
rng = random.Random(1)
genome = "".join(
    f">{c}\n{''.join(rng.choice('ACGT') for _ in range(3000))}\n" for c in ("chr1", "chr2")
)
(t / "genome.fa").write_text(genome)
(t / "a.bed").write_text("chr1\t500\t501\ta1\t12.5\nchr1\t1500\t1501\ta2\t8.0\nchr2\t700\t701\ta3\t5.0\n")
(t / "b.bed").write_text("chr1\t520\t521\tb1\t20.0\nchr2\t2000\t2001\tb2\t9.0\n")
result = gimmemotifs_combine_peak_summits(
    [str(t / "a.bed"), str(t / "b.bed")], genome_fasta=str(t / "genome.fa"), verbose=False
)
print(result["columns"], result["row_count"])
for row in result["rows"]:
    print(row["chrom"], row["start"], row["end"])
```

- **Example output:**

```text
['chrom', 'start', 'end'] 4
chr1 520 521
chr1 1500 1501
chr2 700 701
chr2 2000 2001
```

## 4) `dnakit.sequence.gimmemotifs_coverage_table()`

- **Purpose:** Count the reads of BAM or BED files in each peak region.
- **Method:** Copy the peaks to a temporary BED file with a `+` strand column, because upstream reverses the window of any other strand value and then fails; only the peak centre is used, so the result does not change. Then call GimmeMotifs `coverage_table()` with the window, log transform, normalization, top-region selection, duplicate and repeat filters. BAM files must already have a `.bai` index, because the upstream automatic indexing refers to an undefined name.
- **API:** `gimmemotifs_coverage_table(peak_file: 'str', data_files: 'Sequence[str]', *, output_path: 'str | None' = None, window: 'int' = 200, log_transform: 'bool' = True, normalization: 'str | None' = None, top: 'int' = -1, top_method: 'str' = 'var', remove_duplicates: 'bool' = True, remove_repeats: 'bool' = True, threads: 'int' = 1, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package `gimmemotifs` (extra `dnakit[bio-gimmemotifs]`), which depends on pysam and pybedtools.
- **Input:** `peak_file` is a BED file with at least chrom, start and end. `data_files` is a non-empty list of BAM (with index) or BED files. `normalization` is None, `"quantile"` or `"scale"`. `top_method` is `"var"`, `"std"`, `"mean"` or `"random"`; `top=-1` keeps all regions. `output_path`, when given, must not exist.
- **Output:** `columns` (one per data file), `rows` (one dictionary per region with an `index` such as `chr1:400-600` and one value per file), `row_count`, and `outputs` with `table` when `output_path` is given.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import gimmemotifs_coverage_table

t = Path(tempfile.mkdtemp())
(t / "peaks.bed").write_text("chr1\t400\t600\tp1\t1\nchr1\t1400\t1600\tp2\t1\nchr2\t600\t800\tp3\t1\n")
reads1 = "".join(f"chr1\t{450 + i * 10}\t{480 + i * 10}\tr{i}\t0\t+\n" for i in range(8))
reads2 = "".join(f"chr1\t{1450 + i * 10}\t{1480 + i * 10}\tr{i}\t0\t+\n" for i in range(5))
(t / "reads1.bed").write_text(reads1 + "chr2\t650\t680\tq1\t0\t+\n")
(t / "reads2.bed").write_text(reads2 + "chr2\t650\t680\tq1\t0\t+\nchr2\t660\t690\tq2\t0\t+\n")
result = gimmemotifs_coverage_table(
    str(t / "peaks.bed"),
    [str(t / "reads1.bed"), str(t / "reads2.bed")],
    log_transform=False,
    remove_duplicates=False,
    remove_repeats=False,
    verbose=False,
)
print(result["columns"], result["row_count"])
for row in sorted(result["rows"], key=lambda r: r["index"]):
    print(row["index"], row["reads1"], row["reads2"])
```

- **Example output:**

```text
['reads1', 'reads2'] 3
chr1:1400-1600 0.0 5.0
chr1:400-600 8.0 0.0
chr2:600-800 1.0 2.0
```

## 5) `dnakit.sequence.gimmemotifs_create_background()`

- **Purpose:** Create a set of background sequences (or regions) for motif analysis.
- **Method:** Check the arguments first, because upstream exits the process on invalid ones, then call GimmeMotifs `create_background_file()`. `random` builds a first-order Markov model of `input_file`; `genomic` draws random genome regions of `size` (or the median length of `input_file`); `gc` draws genome regions that match the GC content of `input_file` and needs the `bedtools` program for GC windows. The `promoter` type needs a genomepy annotation and is rejected.
- **API:** `gimmemotifs_create_background(output_path: 'str', *, background_type: 'str', fmt: 'str' = 'fasta', size: 'int | None' = None, genome: 'str | None' = None, input_file: 'str | None' = None, number: 'int' = 10000, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package `gimmemotifs` (extra `dnakit[bio-gimmemotifs]`). The `gc` type also needs `bedtools` on PATH.
- **Input:** `output_path` must not exist. `background_type` is `random`, `genomic` or `gc`. `fmt` is `fasta` or `bed` (a random background can only be FASTA). `random` needs `input_file`; `genomic` and `gc` need `genome`, a FASTA path; `genomic` also needs `size` or `input_file`; `gc` also needs `input_file`. `number` is a positive integer.
- **Output:** `outputs` with `background` and `sequence_count`, which is the requested `number` (the function does not recount the file).
- **Example code:**

```python
import random
import tempfile
from pathlib import Path

from dnakit.sequence import gimmemotifs_create_background

t = Path(tempfile.mkdtemp())
rng = random.Random(1)
genome = "".join(
    f">{c}\n{''.join(rng.choice('ACGT') for _ in range(3000))}\n" for c in ("chr1", "chr2")
)
(t / "genome.fa").write_text(genome)
result = gimmemotifs_create_background(
    t / "background.fa",
    background_type="genomic",
    genome=str(t / "genome.fa"),
    size=100,
    number=20,
    verbose=False,
)
records = (t / "background.fa").read_text().split(">")[1:]
print(result["sequence_count"], len(records), {len(r.split("\n", 1)[1].strip()) for r in records})
```

- **Example output:**

```text
20 20 {100}
```

## 6) `dnakit.sequence.gimmemotifs_denovo_motifs()`

- **Purpose:** Predict de novo motifs from a FASTA, BED or narrowPeak file with an ensemble of external tools, then filter, cluster and report them.
- **Method:** Check that the input exists and that at least one de novo program known to GimmeMotifs is configured and executable, otherwise raise `BackendUnavailableError` before GimmeMotifs is called. Then call GimmeMotifs `gimme_motifs()` with the output directory, the parameter mapping and the three switches. Every file written below the output directory is listed as `file_0`, `file_1`, and so on.
- **API:** `gimmemotifs_denovo_motifs(input_file: 'str', *, output_dir: 'str', params: 'Mapping[str, Any] | None' = None, filter_significant: 'bool' = True, cluster: 'bool' = True, create_report: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package `gimmemotifs` (extra `dnakit[bio-gimmemotifs]`) and at least one configured de novo program (for example MEME or HOMER on PATH, registered in `~/.config/gimmemotifs/gimmemotifs.cfg`). The default parameters ask for a `gc` background, which needs a genome (`params["genome"]`); without one GimmeMotifs logs "Need a genome to create background" and exits the whole process. A full de novo run needs tools and a genome that cannot be assumed here, so the example below shows the input check: a missing input file is refused with `ConfigurationError`.
- **Input:** `input_file` is an existing file. `output_dir` must not exist. `params` follows the GimmeMotifs parameter mapping (`tools`, `analysis`, `background`, `genome`, `fraction`, and others).
- **Output:** `outputs` with the written files, in sorted order of their relative paths.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.sequence import gimmemotifs_denovo_motifs

t = Path(tempfile.mkdtemp())
try:
    gimmemotifs_denovo_motifs(str(t / "missing.fa"), output_dir=t / "denovo", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(": ")[0])
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] input_file does not exist
```

Usage with the dependency installed:

```python
from dnakit.sequence import gimmemotifs_denovo_motifs

result = gimmemotifs_denovo_motifs(
    "peaks.fa",
    output_dir="denovo_out",
    params={"tools": "MEME,Homer", "background": "random", "analysis": "small"},
)
print(result["outputs"].keys())
```

## 7) `dnakit.sequence.gimmemotifs_denovo_predict()`

- **Purpose:** Run the configured de novo tools on a FASTA file against a background FASTA file and write the predicted motifs.
- **Method:** Check the two input files and the new output path, check that at least one de novo program is configured and executable, then call GimmeMotifs `predict_motifs()`, which runs the tools named in `params["tools"]`. The function raises `BackendExecutionError` if no motif file is written.
- **API:** `gimmemotifs_denovo_predict(input_fasta: 'str', background_fasta: 'str', output_path: 'str', *, params: 'Mapping[str, Any] | None' = None, stats_fg: 'Mapping[str, Any] | None' = None, stats_bg: 'Mapping[str, Any] | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package `gimmemotifs` (extra `dnakit[bio-gimmemotifs]`) and at least one configured de novo program. The tools cannot be assumed here, so the example below shows the input check: a missing background file is refused with `ConfigurationError`.
- **Input:** `input_fasta` and `background_fasta` are existing FASTA files. `output_path` must not exist. `params` is the GimmeMotifs parameter mapping; `stats_fg` and `stats_bg` are optional statistics inputs that are passed through.
- **Output:** `outputs` with `motifs`, the path of the PFM file written by GimmeMotifs.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.sequence import gimmemotifs_denovo_predict

t = Path(tempfile.mkdtemp())
try:
    gimmemotifs_denovo_predict("seqs.fa", str(t / "missing.fa"), t / "motifs.pfm", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(": ")[0])
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] background_fasta does not exist
```

Usage with the dependency installed:

```python
from dnakit.sequence import gimmemotifs_denovo_predict

result = gimmemotifs_denovo_predict(
    "peaks.fa",
    "background.fa",
    "predicted.pfm",
    params={"tools": "MEME,Homer", "analysis": "small", "genome": "hg38"},
)
print(result["outputs"]["motifs"])
```

## 8) `dnakit.sequence.gimmemotifs_gc_bin_bed()`

- **Purpose:** Sample BED regions spread evenly over GC-content bins of a genome.
- **Method:** Call GimmeMotifs `gc_bin_bedfile()`. It computes the GC content of genome windows with pybedtools and caches that index in `<cache dir>/gimmemotifs/<genome file name>.gcfreq.<min_bin_size>.feather`, keyed by the file name and not by the path, so a different genome with the same file name reuses the cached index. The cache directory follows `XDG_CACHE_HOME`; the example sets it to a temporary directory.
- **API:** `gimmemotifs_gc_bin_bed(output_path: 'str', *, genome: 'str', number: 'int', length: 'int' = 200, bins: 'Sequence[float] | None' = None, random_state: 'int | None' = None, min_bin_size: 'int' = 100, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package `gimmemotifs` (extra `dnakit[bio-gimmemotifs]`) and the `bedtools` program on PATH (used through pybedtools).
- **Input:** `output_path` must not exist. `genome` is an existing FASTA file. `number` is a positive integer and must be at least the number of bins; `number // bins` regions are sampled per bin, with replacement. `length` is the region length. `bins` is a list of `(low, high)` GC ranges; the default is 0-0.2, 0.2-0.8 in steps of 0.05, and 0.8-1. `random_state` makes the sample repeatable.
- **Output:** `outputs` with `bed`, and `region_count`, which is the requested `number`. The BED file has chrom, start, end and the bin label.
- **Example code:**

```python
import os
import random
import tempfile
from pathlib import Path

t = Path(tempfile.mkdtemp())
os.environ["XDG_CACHE_HOME"] = str(t / "cache")

from dnakit.sequence import gimmemotifs_gc_bin_bed

rng = random.Random(1)
with (t / "genome.fa").open("w") as handle:
    for chrom in ("chr1", "chr2"):
        seq = ""
        for block in range(24):
            gc = (0.3, 0.4, 0.5, 0.6)[block % 4]
            seq += "".join(
                rng.choice("GC") if rng.random() < gc else rng.choice("AT") for _ in range(500)
            )
        handle.write(f">{chrom}\n{seq}\n")
result = gimmemotifs_gc_bin_bed(
    t / "bins.bed",
    genome=str(t / "genome.fa"),
    number=12,
    length=100,
    bins=[(0.2, 0.4), (0.4, 0.6)],
    random_state=1,
    verbose=False,
)
rows = [line.split("\t") for line in (t / "bins.bed").read_text().splitlines()]
print(result["region_count"], len(rows), sorted({row[3] for row in rows}))
```

- **Example output:**

```text
12 12 ['0.20-0.40', '0.40-0.60']
```

## 9) `dnakit.sequence.gimmemotifs_gc_bin_index()`

- **Purpose:** Write an index of the GC content of fixed-size genome windows, for GC-matched sampling.
- **Method:** Call GimmeMotifs `create_gc_bin_index()`. It cuts the genome into windows of `min_bin_size` with pybedtools, computes their GC content and N count, adds averages over 2 and 5 windows, and writes a Feather file. The function checks that the file was created.
- **API:** `gimmemotifs_gc_bin_index(genome: 'str', output_path: 'str', *, min_bin_size: 'int' = 100, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package `gimmemotifs` (extra `dnakit[bio-gimmemotifs]`) and the `bedtools` program on PATH (used through pybedtools).
- **Input:** `genome` is an existing FASTA file. `output_path` must not exist. `min_bin_size` is a positive integer; a small value makes a very large index.
- **Output:** `outputs` with `index`, the Feather file.
- **Example code:**

```python
import random
import tempfile
from pathlib import Path

import pandas

from dnakit.sequence import gimmemotifs_gc_bin_index

t = Path(tempfile.mkdtemp())
rng = random.Random(1)
(t / "genome.fa").write_text(f">chr1\n{''.join(rng.choice('ACGT') for _ in range(2000))}\n")
result = gimmemotifs_gc_bin_index(str(t / "genome.fa"), t / "gc.feather", verbose=False)
table = pandas.read_feather(t / "gc.feather")
print(sorted(result["outputs"]), len(table), list(table.columns))
```

- **Example output:**

```text
['index'] 20 ['chrom', 'start', 'end', 'w100', 'n100', 'w200', 'n200', 'w500', 'n500']
```

## 10) `dnakit.sequence.gimmemotifs_histogram_plot()`

- **Purpose:** Draw a histogram of numbers as an SVG file.
- **Method:** Convert the values to floats, append `.svg` to the output path when it is missing, and call GimmeMotifs `plot_histogram()` with the bin count, range, title and axis label.
- **API:** `gimmemotifs_histogram_plot(values: 'Sequence[float]', output_path: 'str', *, x_range: 'tuple[float, float] | None' = None, bins: 'int' = 10, title: 'str | None' = None, x_label: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package `gimmemotifs` (extra `dnakit[bio-gimmemotifs]`), which uses matplotlib.
- **Input:** `values` is a non-empty list of numbers. `output_path` must not exist after `.svg` is appended. `x_range` is None or `(low, high)` with `low < high`. `bins` is a positive integer.
- **Output:** `outputs` with `image` (the SVG path that was written), `count`, `minimum` and `maximum` of the values.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import gimmemotifs_histogram_plot

t = Path(tempfile.mkdtemp())
result = gimmemotifs_histogram_plot(
    [0.1, 0.2, 0.2, 0.3, 0.5, 0.9], t / "hist", bins=5, title="scores", verbose=False
)
image = Path(result["outputs"]["image"])
print(result["count"], result["minimum"], result["maximum"], image.name, image.read_text()[:4])
```

- **Example output:**

```text
6 0.1 0.9 hist.svg <?xm
```

## 11) `dnakit.sequence.gimmemotifs_jaspar_motif_info()`

- **Purpose:** Fetch the metadata of a JASPAR matrix from the JASPAR web service.
- **Method:** Check that the id is a non-empty string, then call GimmeMotifs `get_jaspar_motif_info()`, which requests `http://jaspar.genereg.net/api/v1/matrix/<id>?format=json` and returns the JSON object. This is the only function of the group that uses the network.
- **API:** `gimmemotifs_jaspar_motif_info(motif_id: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package `gimmemotifs` (extra `dnakit[bio-gimmemotifs]`) and internet access to `jaspar.genereg.net`; no key is needed. The example below does not touch the network, so it shows the input check: an empty id is refused with `ConfigurationError`.
- **Input:** `motif_id` is a JASPAR matrix id such as `MA0106.3`, a non-empty string.
- **Output:** `motif_id` and `info`, the JSON object that JASPAR sends.
- **Example code:**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.sequence import gimmemotifs_jaspar_motif_info

try:
    gimmemotifs_jaspar_motif_info("  ", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error))
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] motif_id must be a non-empty string.
```

Usage with the dependency installed:

```python
from dnakit.sequence import gimmemotifs_jaspar_motif_info

result = gimmemotifs_jaspar_motif_info("MA0106.3")
print(result["info"]["name"])
```

## 12) `dnakit.sequence.gimmemotifs_logo_images()`

- **Purpose:** Draw sequence logos for the motifs of a PFM file as PNG files.
- **Method:** Run the `gimme logo` command in the new output directory, with `-i` for the selected ids, `-k` for the logo kind and `--notitle` when titles are off. Every PNG found below the directory is returned; if there is none the function raises `BackendExecutionError`.
- **API:** `gimmemotifs_logo_images(pfm_file: 'str', *, output_dir: 'str', motif_ids: 'Sequence[str] | None' = None, kind: 'str | None' = None, show_title: 'bool' = True, gimme_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package `gimmemotifs` (extra `dnakit[bio-gimmemotifs]`) and its `gimme` program (from `gimme_path`, PATH, or the directory of the running interpreter).
- **Input:** `pfm_file` is an existing PFM file. `output_dir` must not exist. `motif_ids` selects motifs (all by default). `kind` is None, `information`, `frequency`, `energy` or `ensembl`.
- **Output:** `outputs` with `image_0`, `image_1`, and so on (sorted PNG paths), and `log_tail`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import gimmemotifs_logo_images

out = Path(tempfile.mkdtemp()) / "logos"
result = gimmemotifs_logo_images(
    "motifs.pfm", output_dir=out, motif_ids=["MA0001.1 AR1", "MA0002.1 GATA"], verbose=False
)
print(sorted(Path(path).name for path in result["outputs"].values()))
```

- **Example output:**

```text
['MA0001.1 AR1.png', 'MA0002.1 GATA.png']
```

## 13) `dnakit.sequence.gimmemotifs_match_motifs()`

- **Purpose:** Find the best matching database motif for each query motif.
- **Method:** Run `gimme match -d <database> -n <nreport> <query>` and stream its table to the new output file, then parse the table into rows. The graphical report of `gimme match -o` is not offered, because it fails in this GimmeMotifs release with a numpy broadcasting error.
- **API:** `gimmemotifs_match_motifs(pfm_file: 'str', *, database_pfm: 'str', output_path: 'str', nreport: 'int' = 1, gimme_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package `gimmemotifs` (extra `dnakit[bio-gimmemotifs]`) and its `gimme` program (from `gimme_path`, PATH, or the directory of the running interpreter).
- **Input:** `pfm_file` is the query PFM file and `database_pfm` the database PFM file, both existing. `output_path` must not exist. `nreport` is the number of matches per query (a positive integer).
- **Output:** `rows` (one dictionary per match with the table columns `Motif`, `Match`, `Score` and `P-value`, all strings), `row_count`, `outputs` with `table`, and `log_tail`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import gimmemotifs_match_motifs

out = Path(tempfile.mkdtemp()) / "match.tsv"
result = gimmemotifs_match_motifs("motifs.pfm", database_pfm="db.pfm", output_path=out, verbose=False)
print(result["row_count"])
for row in result["rows"]:
    print(row["Motif"], "->", row["Match"], row["Score"])
```

- **Example output:**

```text
3
MA0001.1 AR1 -> DB1 long 0.99
MA0002.1 GATA -> DB1 long 0.23
MA0001.2 AR1b -> DB1 long 0.96
```

## 14) `dnakit.sequence.gimmemotifs_match_plot()`

- **Purpose:** Plot each query motif next to its best database match with the p-value.
- **Method:** Read every query and database motif by id from its PFM file, then call GimmeMotifs `match_plot()` with the list of `(query, match, p-value)` triples and write the image to the new output path.
- **API:** `gimmemotifs_match_plot(matches: 'Sequence[Mapping[str, Any]]', output_path: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package `gimmemotifs` (extra `dnakit[bio-gimmemotifs]`), which uses matplotlib.
- **Input:** `matches` is a non-empty list of mappings with `query_pfm`, `query_id`, `database_pfm`, `database_id` (existing files and motif ids in them) and a numeric `p_value`. `output_path` must not exist.
- **Output:** `outputs` with `image`, and `pairs`, the number of plotted pairs.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import gimmemotifs_match_plot

out = Path(tempfile.mkdtemp()) / "match.png"
result = gimmemotifs_match_plot(
    [
        {
            "query_pfm": "motifs.pfm",
            "query_id": "MA0001.1 AR1",
            "database_pfm": "db.pfm",
            "database_id": "DB1 long",
            "p_value": 0.01,
        }
    ],
    out,
    verbose=False,
)
print(result["pairs"], out.read_bytes()[:4])
```

- **Example output:**

```text
1 b'\x89PNG'
```

## 15) `dnakit.sequence.gimmemotifs_matched_gc_bed()`

- **Purpose:** Sample BED regions whose GC content follows that of the regions of another BED file.
- **Method:** Call GimmeMotifs `matched_gc_bedfile()`. It reads the GC content of the regions in `match_bed` (a BED or FASTA file), draws genome windows from a GC index of the genome and writes the sampled regions. The index is cached by genome file name in the cache directory, as described in section 8.
- **API:** `gimmemotifs_matched_gc_bed(output_path: 'str', *, match_bed: 'str', genome: 'str', number: 'int', size: 'int | None' = None, min_bin_size: 'int' = 100, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package `gimmemotifs` (extra `dnakit[bio-gimmemotifs]`) and the `bedtools` program on PATH (used through pybedtools).
- **Input:** `output_path` must not exist. `match_bed` and `genome` are existing files; the genome is a FASTA file. `number` is a positive integer. `size` is the region length (the length of the input regions by default). `min_bin_size` is a positive integer.
- **Output:** `outputs` with `bed`. The function does not count the regions.
- **Example code:**

```python
import os
import random
import tempfile
from pathlib import Path

t = Path(tempfile.mkdtemp())
os.environ["XDG_CACHE_HOME"] = str(t / "cache")

from dnakit.sequence import gimmemotifs_matched_gc_bed

rng = random.Random(1)
with (t / "genome.fa").open("w") as handle:
    for chrom in ("chr1", "chr2"):
        seq = ""
        for block in range(24):
            gc = (0.3, 0.4, 0.5, 0.6)[block % 4]
            seq += "".join(
                rng.choice("GC") if rng.random() < gc else rng.choice("AT") for _ in range(500)
            )
        handle.write(f">{chrom}\n{seq}\n")
(t / "match.bed").write_text("chr1\t0\t100\nchr1\t500\t600\nchr2\t1000\t1100\n")
result = gimmemotifs_matched_gc_bed(
    t / "matched.bed",
    match_bed=str(t / "match.bed"),
    genome=str(t / "genome.fa"),
    number=10,
    verbose=False,
)
rows = [line.split("\t") for line in (t / "matched.bed").read_text().splitlines()]
print(sorted(result["outputs"]), len(rows), {int(r[2]) - int(r[1]) for r in rows})
```

- **Example output:**

```text
['bed'] 10 {100}
```

## 16) `dnakit.sequence.gimmemotifs_motif_image_series()`

- **Purpose:** Write one logo PNG for each listed motif id.
- **Method:** Read the motifs of the PFM file, check that every id is present, then call GimmeMotifs `motif_to_img_series()`, which draws a PNG per motif under `<output_dir>/<subdir>`. Characters other than letters, digits and `-` in the id are replaced by `_` in the file name. The function checks that every image exists.
- **API:** `gimmemotifs_motif_image_series(motif_ids: 'Sequence[str]', *, pfm_file: 'str', output_dir: 'str', subdir: 'str' = 'logos', verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package `gimmemotifs` (extra `dnakit[bio-gimmemotifs]`) and pandas.
- **Input:** `motif_ids` is a non-empty list of ids that are in `pfm_file`. `output_dir` must not exist. `subdir` is the subdirectory for the images.
- **Output:** `images` (motif id mapped to the PNG path) and `image_count`. The `outputs` mapping is empty.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import gimmemotifs_motif_image_series

out = Path(tempfile.mkdtemp()) / "images"
result = gimmemotifs_motif_image_series(
    ["MA0001.1 AR1", "MA0002.1 GATA"], pfm_file="motifs.pfm", output_dir=out, verbose=False
)
print(result["image_count"])
for motif_id, path in result["images"].items():
    print(motif_id, "->", Path(path).relative_to(out))
```

- **Example output:**

```text
2
MA0001.1 AR1 -> logos/MA0001_1_AR1.png
MA0002.1 GATA -> logos/MA0002_1_GATA.png
```

## 17) `dnakit.sequence.gimmemotifs_motif_localization()`

- **Purpose:** Plot where a motif occurs along equal-length sequences and test whether the positions are uniform.
- **Method:** Read the motif by id and call GimmeMotifs `motif_localization()`, which scans the FASTA file, draws a histogram of the match positions as SVG and returns a Kolmogorov-Smirnov p-value against a uniform distribution. Upstream writes the SVG only when the motif has matches; otherwise `outputs` has no `image`.
- **API:** `gimmemotifs_motif_localization(fasta_file: 'str', *, pfm_file: 'str', motif_id: 'str', output_path: 'str', size: 'int', cutoff: 'float' = 0.9, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package `gimmemotifs` (extra `dnakit[bio-gimmemotifs]`), which uses matplotlib.
- **Input:** `fasta_file` and `pfm_file` are existing files; `motif_id` is an id in the PFM file. `output_path` must not exist. `size` is the sequence length used for the histogram and must be larger than the motif length. `cutoff` is the score cutoff, a number in (0, 1].
- **Output:** `motif_id`, `p_value`, and `outputs` with `image` when the SVG was written.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import gimmemotifs_motif_localization

out = Path(tempfile.mkdtemp()) / "location.svg"
result = gimmemotifs_motif_localization(
    "seqs.fa",
    pfm_file="motifs.pfm",
    motif_id="MA0001.1 AR1",
    output_path=out,
    size=37,
    cutoff=0.9,
    verbose=False,
)
print(result["motif_id"], round(result["p_value"], 3), sorted(result["outputs"]), out.exists())
```

- **Example output:**

```text
MA0001.1 AR1 0.244 ['image'] True
```

## 18) `dnakit.sequence.gimmemotifs_motif_location_plots()`

- **Purpose:** Plot where each motif occurs relative to the centre of the sequences.
- **Method:** Run `gimme location` in the new output directory with the score cutoff, optional size and optional motif ids. Every file found below the directory is returned.
- **API:** `gimmemotifs_motif_location_plots(fasta_file: 'str', *, pfm_file: 'str', output_dir: 'str', size: 'int | None' = None, motif_ids: 'Sequence[str] | None' = None, cutoff: 'float' = 0.95, gimme_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package `gimmemotifs` (extra `dnakit[bio-gimmemotifs]`) and its `gimme` program (from `gimme_path`, PATH, or the directory of the running interpreter).
- **Input:** `fasta_file` and `pfm_file` are existing files. `output_dir` must not exist. `size` is a positive integer (the sequence length used for the plot). `motif_ids` selects motifs. `cutoff` is a number in (0, 1].
- **Output:** `outputs` with `file_0`, `file_1`, and so on (sorted paths; one histogram SVG per motif), and `log_tail`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import gimmemotifs_motif_location_plots

out = Path(tempfile.mkdtemp()) / "location"
result = gimmemotifs_motif_location_plots(
    "seqs.fa", pfm_file="motifs.pfm", output_dir=out, size=36, verbose=False
)
print(sorted(Path(path).name for path in result["outputs"].values()))
```

- **Example output:**

```text
['MA0001.1 AR1_histogram.svg', 'MA0001.2 AR1b_histogram.svg', 'MA0002.1 GATA_histogram.svg']
```

## 19) `dnakit.sequence.gimmemotifs_parallel_denovo_predict()`

- **Purpose:** Predict motifs in parallel with several external de novo tools.
- **Method:** Check the FASTA file and the new output path, check that at least one de novo program is configured and executable, then call GimmeMotifs `pp_predict_motifs()` with the analysis size, organism, tool list, number of threads and time limit. The function raises `BackendExecutionError` if no motif file is written.
- **API:** `gimmemotifs_parallel_denovo_predict(fasta_file: 'str', output_path: 'str', *, analysis: 'str' = 'small', organism: 'str' = 'hg19', single: 'bool' = False, background: 'str' = '', tools: 'str | None' = None, threads: 'int' = 8, max_time: 'int' = -1, stats_fg: 'Mapping[str, Any] | None' = None, stats_bg: 'Mapping[str, Any] | None' = None, gc: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package `gimmemotifs` (extra `dnakit[bio-gimmemotifs]`) and at least one configured de novo program. The example below does not run a tool: a missing input file is refused with `ConfigurationError` before any program starts.
- **Input:** `fasta_file` is an existing FASTA file. `output_path` must not exist. `analysis` is the GimmeMotifs analysis size (`xs`, `small`, `medium`, `large`, `xl`). `tools` is a comma-separated list of tool names (the GimmeMotifs default set when None). `threads` is a positive integer. `max_time` is a time limit in seconds (-1 for none).
- **Output:** `outputs` with `motifs`, the PFM file written by GimmeMotifs.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.sequence import gimmemotifs_parallel_denovo_predict

t = Path(tempfile.mkdtemp())
try:
    gimmemotifs_parallel_denovo_predict(str(t / "missing.fa"), t / "motifs.pfm", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(": ")[0])
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] fasta_file does not exist
```

Usage with the dependency installed:

```python
from dnakit.sequence import gimmemotifs_parallel_denovo_predict

result = gimmemotifs_parallel_denovo_predict(
    "peaks.fa", "predicted.pfm", analysis="small", tools="MEME,Homer", threads=4
)
print(result["outputs"]["motifs"])
```

## 20) `dnakit.sequence.gimmemotifs_pfm_scan()`

- **Purpose:** Scan FASTA or BED sequences with PFM motifs using a fixed score cutoff.
- **Method:** Run `gimme scan -p <pfm> -c <cutoff> -n <nreport> -N <threads>` with `-t` (table), `-b` (BED) or `-T` (score table), plus `-r` when reverse complements are off, and stream the result to `scan.tsv` or `scan.bed` in the new output directory. The cutoff is required because the FPR mode of `gimme scan` needs a background set or a genome, which this wrapper does not take. Tabular output is parsed into rows.
- **API:** `gimmemotifs_pfm_scan(input_file: 'str', *, pfm_file: 'str', output_dir: 'str', cutoff: 'float', nreport: 'int' = 1, reverse_complement: 'bool' = True, output_format: 'str' = 'table', threads: 'int' = 1, gimme_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package `gimmemotifs` (extra `dnakit[bio-gimmemotifs]`) and its `gimme` program (from `gimme_path`, PATH, or the directory of the running interpreter).
- **Input:** `input_file` is a FASTA or BED file and `pfm_file` a PFM file, both existing. `output_dir` must not exist. `cutoff` is a number in (0, 1] (the motif score cutoff). `output_format` is `table`, `bed` or `score_table`. `nreport` and `threads` are positive integers.
- **Output:** `outputs` with `scan`, `rows` (one dictionary per sequence; for the table format the count of matches per motif as a string), `row_count` (0 for BED output) and `log_tail`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import gimmemotifs_pfm_scan

out = Path(tempfile.mkdtemp()) / "scan"
result = gimmemotifs_pfm_scan(
    "seqs.fa", pfm_file="motifs.pfm", output_dir=out, cutoff=0.9, verbose=False
)
print(result["row_count"], Path(result["outputs"]["scan"]).name)
for row in result["rows"]:
    print(row["sequence"], row["MA0001.1 AR1"], row["MA0002.1 GATA"], row["MA0001.2 AR1b"])
```

- **Example output:**

```text
6 scan.tsv
s1 1 0 1
s2 1 1 1
s3 1 1 1
s4 0 0 0
s5 1 1 1
s6 0 0 0
```

## 21) `dnakit.sequence.gimmemotifs_roc_curve_plot()`

- **Purpose:** Plot one or several ROC curves with 1 - specificity on the x axis.
- **Method:** Convert the inputs to float arrays, add `.png` when the path has no image extension, and call GimmeMotifs `roc_plot()` with the curve ids as legend. The accepted extensions are those of GimmeMotifs (`.png`, `.pdf`, `.svg`, `.ps`).
- **API:** `gimmemotifs_roc_curve_plot(output_path: 'str', *, plot_x: 'Sequence[Any]', plot_y: 'Sequence[Any]', curve_ids: 'Sequence[str] | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package `gimmemotifs` (extra `dnakit[bio-gimmemotifs]`), which uses matplotlib.
- **Input:** `plot_x` and `plot_y` are non-empty and of equal length: either one curve (a list of numbers) or a list of curves (a list of lists). `curve_ids` names the curves. `output_path` must not exist after the extension is added.
- **Output:** `outputs` with `image` (the path that was written) and `curves`, the number of curves.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import gimmemotifs_roc_curve_plot

out = Path(tempfile.mkdtemp()) / "roc"
result = gimmemotifs_roc_curve_plot(
    out,
    plot_x=[[0, 0.2, 1], [0, 0.4, 1]],
    plot_y=[[0, 0.8, 1], [0, 0.6, 1]],
    curve_ids=["motif A", "motif B"],
    verbose=False,
)
image = Path(result["outputs"]["image"])
print(result["curves"], image.name, image.read_bytes()[1:4])
```

- **Example output:**

```text
2 roc.png b'PNG'
```

## 22) `dnakit.sequence.gimmemotifs_gimme_motifs()`

- **Purpose:** Run `gimme motifs` on a sample FASTA file: scan known motifs (and optionally de novo motifs), compute enrichment statistics against a background and write the reports.
- **Method:** Check the arguments and the de novo tools (each must be configured and installed, otherwise `BackendUnavailableError` is raised before any work starts), then call the GimmeMotifs `motifs` command function with an explicit argument set. Known motifs of `pfm_file` are scanned on the sample and on the background; ROC, precision-recall and enrichment statistics go to `gimme.roc.report.txt`; a non-redundant subset is chosen and its scan results are written as BED files in `motif_scan_results`; the HTML reports `gimme.motifs.html` and `gimme.motifs.redundant.html` are written when `report` is true. The genome is not an input, so z-score and GC normalisation are off. When `background_fasta` is None, GimmeMotifs generates a random background from the sample. In the verified environment (xxhash 4.0.1) a run in which at least one non-redundant motif is selected fails with `TypeError: Strings must be encoded before hashing`, because the later BED scan hashes the motif with a `str`; the example uses data in which no motif is selected.
- **API:** `gimmemotifs_gimme_motifs(sample_fasta: 'str', pfm_file: 'str', output_dir: 'str', *, background_fasta: 'str | None' = None, known: 'bool' = True, denovo_tools: 'Sequence[str] | None' = None, analysis: 'str' = 'xl', fraction: 'float' = 0.2, single: 'bool' = False, report: 'bool' = True, ncpus: 'int' = 1, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package `gimmemotifs` (extra `dnakit[bio-gimmemotifs]`). De novo tools named in `denovo_tools` must be configured and installed; the default (known motifs only) needs none.
- **Input:** `sample_fasta`, `pfm_file` and the optional `background_fasta` are existing files. `output_dir` must not exist. `analysis` is `small`, `medium`, `large` or `xl`. `fraction` is a number in (0, 1]. `known=False` needs `denovo_tools`. `denovo_tools` is a list of tool names such as `["MEME", "Homer", "BioProspector"]`. `ncpus` is at least 1.
- **Output:** `known`, `denovo_tools`, `motifs_scored` (the number of motifs in the ROC report), and `outputs` with `roc_report`, `motif_scan_results` and, when present, `gimme.motifs.html`, `gimme.motifs.redundant.html` and `combined.pfm`.
- **Example code:**

```python
import os
import random
import tempfile
from pathlib import Path

t = Path(tempfile.mkdtemp())
os.environ["XDG_CACHE_HOME"] = str(t / "cache")

from dnakit.sequence import gimmemotifs_gimme_motifs

rng = random.Random(1)


def seq(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


sample = [seq(200) for _ in range(40)]
sample = [s[:50] + "TGACGTCA" + s[58:] if i % 2 == 0 else s for i, s in enumerate(sample)]
(t / "sample.fa").write_text("".join(f">s{i}\n{s}\n" for i, s in enumerate(sample)))
(t / "background.fa").write_text("".join(f">b{i}\n{seq(200)}\n" for i in range(80)))
result = gimmemotifs_gimme_motifs(
    str(t / "sample.fa"),
    "motifs.pfm",
    t / "run",
    background_fasta=str(t / "background.fa"),
    verbose=False,
)
print(result["known"], result["denovo_tools"], result["motifs_scored"])
print(sorted(result["outputs"]))
print((t / "run" / "gimme.roc.report.txt").read_text().splitlines()[0].split("\t")[:5])
```

- **Example output:**

```text
True [] 3
['gimme.motifs.html', 'gimme.motifs.redundant.html', 'motif_scan_results', 'roc_report']
['Motif', '# matches', '% matches input', '# matches background', '% matches background']
```

## 23) `dnakit.sequence.gimmemotifs_maelstrom_html_report()`

- **Purpose:** Write the HTML report of a finished maelstrom run folder.
- **Method:** Read `final.out.txt` (the motif-by-method table) from the run folder, count the motifs that have an absolute value at or above the threshold in at least one method column (columns whose names contain `corr` or `% with motif` are not method columns), then call GimmeMotifs `maelstrom_html_report()`. The report `gimme.maelstrom.report.html` is written in the run folder and an existing report there is replaced, as upstream does.
- **API:** `gimmemotifs_maelstrom_html_report(run_folder: 'str', pfm_file: 'str', *, threshold: 'float' = 3, plot_all_motifs: 'bool' = False, plot_no_motifs: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package `gimmemotifs` (extra `dnakit[bio-gimmemotifs]`) and pandas.
- **Input:** `run_folder` is an existing folder with `final.out.txt`; `pfm_file` is the motif database of the run. `threshold` is a number of at least 0. `plot_all_motifs` draws a logo for every motif; `plot_no_motifs` leaves the logo column out.
- **Output:** `threshold`, `motifs_in_report`, `methods` (the number of method columns), and `outputs` with `report` and, when the folder exists, `logos`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import gimmemotifs_maelstrom_html_report

run = Path(tempfile.mkdtemp())
(run / "final.out.txt").write_text(
    "\tHypergeom\tMean\tHypergeom % with motif\tcorr Hypergeom\n"
    "MA0001.1 AR1\t5.2\t3.1\t80.0\t0.9\n"
    "MA0002.1 GATA\t-1.2\t0.2\t20.0\t0.1\n"
)
result = gimmemotifs_maelstrom_html_report(str(run), "motifs.pfm", threshold=3, verbose=False)
print(result["threshold"], result["motifs_in_report"], result["methods"], sorted(result["outputs"]))
html = (run / "gimme.maelstrom.report.html").read_text()
print("MA0001.1 AR1" in html, "MA0002.1 GATA" in html)
```

- **Example output:**

```text
3.0 1 2 ['logos', 'report']
True False
```

## 24) `dnakit.sequence.gimmemotifs_predict_motifs()`

- **Purpose:** Run one de novo motif tool on a FASTA file and write the predicted motifs as a PFM file.
- **Method:** Look up the tool in GimmeMotifs (it must be configured in `gimmemotifs.cfg` and executable), run its `run()` method with the FASTA file and the parameter mapping that `gimme prediction` reads from its parameter file, and write every predicted motif with `to_pfm()`, one after another, to the new output file.
- **API:** `gimmemotifs_predict_motifs(tool: 'str', fasta: 'str', output_pfm: 'str', *, params: 'Mapping[str, Any] | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package `gimmemotifs` (extra `dnakit[bio-gimmemotifs]`) and the named tool configured and installed (for example MEME or Homer; GimmeMotifs also bundles some tools). A tool that is known but not configured or not executable raises `BackendUnavailableError`. The example below uses an unknown tool name, which is refused before any program starts.
- **Input:** `tool` is a GimmeMotifs tool name (case-insensitive, one of its list of tools). `fasta` is an existing FASTA file. `output_pfm` must not exist. `params` is a mapping of tool parameters (for example `background`, `width`, `number`, `analysis`, `use_strand`) or None.
- **Output:** `tool`, `motifs` (the number of predicted motifs), `stdout_tail` and `stderr_tail` (bounded tails of the tool output), and `outputs` with `pfm`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.sequence import gimmemotifs_predict_motifs

out = Path(tempfile.mkdtemp()) / "predicted.pfm"
try:
    gimmemotifs_predict_motifs("nosuchtool", "seqs.fa", out, verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(";")[0])
print(out.exists())
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] Unknown gimmemotifs tool 'nosuchtool'
False
```

Usage with the dependency installed:

```python
from dnakit.sequence import gimmemotifs_predict_motifs

result = gimmemotifs_predict_motifs(
    "BioProspector",
    "peaks.fa",
    "predicted.pfm",
    params={"background": "background.fa", "width": 8, "number": 3, "analysis": "small", "use_strand": False},
)
print(result["motifs"])
```

## 25) `dnakit.sequence.gimmemotifs_annot2primpep()`

- **Purpose:** Write the primary peptide FASTA of an assembly: the longest protein of each gene.
- **Method:** Check that `gffread` is on PATH, copy the genome FASTA and the GTF annotation into `<output_dir>/<name>/`, and call GimmeMotifs `annot2primpep()`. It runs gffread to translate the annotation and keeps the longest protein per gene. Records are named `gene_name|gene_id` (`.` when the name is missing). The result is `<output_dir>/prim_genes/<name>.pep.fa`.
- **API:** `gimmemotifs_annot2primpep(genome_fasta: 'str', annotation_gtf: 'str', output_dir: 'str', *, name: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package `gimmemotifs` (extra `dnakit[bio-gimmemotifs]`) and the external program `gffread` on PATH (for example `conda install -c bioconda gffread`). The example below runs where `gffread` is missing and shows the error; the program is not installed in the environment used to verify this page.
- **Input:** `genome_fasta` and `annotation_gtf` are existing files. `output_dir` must not exist. `name` is a plain identifier without `/` or `'` (the genome file stem by default).
- **Output:** `assembly` (the name), `protein_count` (the number of records in the peptide file) and `outputs` with `peptides`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.sequence import gimmemotifs_annot2primpep

out = Path(tempfile.mkdtemp()) / "assembly"
try:
    gimmemotifs_annot2primpep("seqs.fa", "motifs.pfm", out, verbose=False)
except BackendUnavailableError as error:
    print(type(error).__name__, str(error).split(" Hint")[0])
print(out.exists())
```

- **Example output:**

```text
BackendUnavailableError [BACKEND_UNAVAILABLE] The gffread program was not found on PATH.
False
```

Usage with the dependency installed:

```python
from dnakit.sequence import gimmemotifs_annot2primpep

result = gimmemotifs_annot2primpep("genome.fa", "genes.gtf", "primary_peptides", name="mygenome")
print(result["assembly"], result["protein_count"], result["outputs"]["peptides"])
```

## 26) `dnakit.sequence.gimmemotifs_best_motif_in_cluster()`

- **Purpose:** Return the best motif of each de novo motif cluster.
- **Method:** Read the single and cluster motifs, resolve the ids of each cluster, compute statistics for every involved motif against each background FASTA (`calc_stats()` with the z-score and GC options off), then call GimmeMotifs `best_motif_in_cluster()` with those statistics, which ranks the motifs of each cluster by the given metrics. The ranking metrics are rank statistics, so they do not change under a per-motif z-score.
- **API:** `gimmemotifs_best_motif_in_cluster(single_motif_file: 'str', cluster_motif_file: 'str', clusters: 'Sequence[Mapping[str, Any]]', foreground_fasta: 'str', background_fastas: 'Mapping[str, str]', genome_fasta: 'str', *, metrics: 'Sequence[str]' = ('roc_auc', 'recall_at_fdr'), verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package `gimmemotifs` (extra `dnakit[bio-gimmemotifs]`).
- **Input:** `single_motif_file` and `cluster_motif_file` are PFM files. `clusters` is a non-empty list of mappings with `representative` (an id in the cluster file) and `members` (ids in the single motif file). `foreground_fasta`, `genome_fasta` and every path in `background_fastas` (a mapping from a background name to a FASTA path) exist. `metrics` are GimmeMotifs statistic names.
- **Output:** `best_motifs`, the ids of the best motif of each cluster, from best to worst.
- **Example code:**

```python
import os
import random
import tempfile
from pathlib import Path

t = Path(tempfile.mkdtemp())
os.environ["XDG_CACHE_HOME"] = str(t / "cache")

from dnakit.sequence import gimmemotifs_best_motif_in_cluster

rng = random.Random(1)


def seq(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


foreground = [seq(200) for _ in range(40)]
foreground = [s[:50] + "TGACAC" + s[56:] if i % 2 == 0 else s for i, s in enumerate(foreground)]
(t / "fg.fa").write_text("".join(f">s{i}\n{s}\n" for i, s in enumerate(foreground)))
(t / "bg.fa").write_text("".join(f">b{i}\n{seq(200)}\n" for i in range(80)))
(t / "genome.fa").write_text("".join(f">{c}\n{seq(5000)}\n" for c in ("chr1", "chr2")))
result = gimmemotifs_best_motif_in_cluster(
    "motifs.pfm",
    "motifs.pfm",
    [{"representative": "MA0001.1 AR1", "members": ["MA0001.1 AR1", "MA0001.2 AR1b"]}],
    str(t / "fg.fa"),
    {"random": str(t / "bg.fa")},
    str(t / "genome.fa"),
    verbose=False,
)
print(result["best_motifs"])
```

- **Example output:**

```text
['MA0001.1 AR1']
```

## 27) `dnakit.sequence.gimmemotifs_denovo_create_background()`

- **Purpose:** Create one background FASTA file for de novo motif analysis.
- **Method:** For `random`, `genomic`, `gc` and `custom` call GimmeMotifs `create_background()` of the de novo module. For `promoter`, build random promoters with `PromoterFasta` from a BED file of genes and write them. `nr_times` background sequences are made per foreground sequence. The genome is always a FASTA path.
- **API:** `gimmemotifs_denovo_create_background(background_type: 'str', fasta_file: 'str', output_path: 'str', *, genome_fasta: 'str | None' = None, annotation_bed: 'str | None' = None, size: 'int' = 200, nr_times: 'int' = 10, custom_background: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package `gimmemotifs` (extra `dnakit[bio-gimmemotifs]`). The `gc` type also needs `bedtools` on PATH.
- **Input:** `background_type` is `random` (first-order Markov model of `fasta_file`), `genomic` (random regions of `genome_fasta`), `gc` (regions of `genome_fasta` matched to the GC content of `fasta_file`), `promoter` (needs `genome_fasta` and `annotation_bed`) or `custom` (`custom_background` is copied as it is). `output_path` must not exist. `size` and `nr_times` are positive integers.
- **Output:** `background_type`, `sequence_count` and `outputs` with `background`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import gimmemotifs_denovo_create_background

out = Path(tempfile.mkdtemp()) / "background.fa"
result = gimmemotifs_denovo_create_background(
    "random", "seqs.fa", out, nr_times=3, verbose=False
)
records = out.read_text().split(">")[1:]
print(result["background_type"], result["sequence_count"], len(records))
```

- **Example output:**

```text
random 18 18
```

## 28) `dnakit.sequence.gimmemotifs_denovo_create_backgrounds()`

- **Purpose:** Create the prediction background and the validation backgrounds of a de novo run.
- **Method:** Copy the prediction and validation FASTA files into the new output directory as `prediction.fa` and `validation.fa`, then call GimmeMotifs `create_backgrounds()`. The prediction background uses the first type (or `gc` when it is listed); one validation background file per type is written.
- **API:** `gimmemotifs_denovo_create_backgrounds(output_dir: 'str', *, prediction_fasta: 'str', validation_fasta: 'str', background_types: 'Sequence[str]' = ('random',), genome_fasta: 'str | None' = None, size: 'int' = 200, custom_background: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package `gimmemotifs` (extra `dnakit[bio-gimmemotifs]`). The `gc` type also needs `bedtools` on PATH.
- **Input:** `output_dir` must not exist. `prediction_fasta` and `validation_fasta` are existing files. `background_types` is a list from `random`, `genomic`, `gc` and `custom` (`promoter` is available through `gimmemotifs_denovo_create_background`). `genomic` and `gc` need `genome_fasta`; `custom` needs `custom_background`.
- **Output:** `background_types`, `background_files` (type mapped to the validation file path) and `outputs` with `prediction_background` and one `validation_<type>` entry per type.
- **Example code:**

```python
import random
import tempfile
from pathlib import Path

from dnakit.sequence import gimmemotifs_denovo_create_backgrounds

t = Path(tempfile.mkdtemp())
rng = random.Random(1)


def seq(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


(t / "prediction.fa").write_text("".join(f">p{i}\n{seq(100)}\n" for i in range(10)))
(t / "validation.fa").write_text("".join(f">v{i}\n{seq(100)}\n" for i in range(10)))
(t / "genome.fa").write_text("".join(f">{c}\n{seq(3000)}\n" for c in ("chr1", "chr2")))
result = gimmemotifs_denovo_create_backgrounds(
    t / "run",
    prediction_fasta=str(t / "prediction.fa"),
    validation_fasta=str(t / "validation.fa"),
    background_types=["random", "genomic"],
    genome_fasta=str(t / "genome.fa"),
    size=100,
    verbose=False,
)
print(sorted(result["outputs"]))
print(sorted(path.name for path in (t / "run").iterdir()))
```

- **Example output:**

```text
['prediction_background', 'validation_genomic', 'validation_random']
['bg.genomic.fa', 'bg.random.fa', 'prediction.bg.fa', 'prediction.fa', 'validation.fa']
```

## 29) `dnakit.sequence.gimmemotifs_denovo_motif_report()`

- **Purpose:** Write the text and HTML report of a set of de novo motifs.
- **Method:** Follow the steps of the upstream report: draw a ROC plot per motif and background, find the closest database motif with `MotifComparer`, compute statistics (`calc_stats()` with z-score and GC off, then `add_star()`), draw a localization histogram per motif from `location_fasta`, then write the text report and the HTML report. Every file written is returned.
- **API:** `gimmemotifs_denovo_motif_report(input_file: 'str', pfm_file: 'str', foreground_fasta: 'str', background_fastas: 'Mapping[str, str]', location_fasta: 'str', output_dir: 'str', *, genome_fasta: 'str', cutoff_fpr: 'float' = 0.9, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package `gimmemotifs` (extra `dnakit[bio-gimmemotifs]`), which uses matplotlib.
- **Input:** `input_file` is an existing file whose name is shown in the report. `pfm_file` holds the motifs. `foreground_fasta`, `location_fasta` (sequences centred on the peaks) and `genome_fasta` are existing FASTA files. `background_fastas` maps background names to FASTA paths. `output_dir` must not exist. `cutoff_fpr` is the score cutoff of the histograms (default 0.9).
- **Output:** `report_folder` and `outputs`, a mapping from relative file path to path: `gimme.denovo.html`, `stats.<background>.txt`, and for each motif a logo, `<id>_roc.<background>.png` and `<id>_histogram.svg`. The closest-match logos of database motifs are also listed.
- **Example code:**

```python
import os
import random
import tempfile
from pathlib import Path

t = Path(tempfile.mkdtemp())
os.environ["XDG_CACHE_HOME"] = str(t / "cache")

from dnakit.sequence import gimmemotifs_denovo_motif_report

rng = random.Random(1)


def seq(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


foreground = [seq(200) for _ in range(40)]
foreground = [s[:50] + "TGACAC" + s[56:] if i % 2 == 0 else s for i, s in enumerate(foreground)]
(t / "fg.fa").write_text("".join(f">s{i}\n{s}\n" for i, s in enumerate(foreground)))
(t / "bg.fa").write_text("".join(f">b{i}\n{seq(200)}\n" for i in range(80)))
(t / "genome.fa").write_text("".join(f">{c}\n{seq(5000)}\n" for c in ("chr1", "chr2")))
result = gimmemotifs_denovo_motif_report(
    str(t / "fg.fa"),
    "motifs.pfm",
    str(t / "fg.fa"),
    {"random": str(t / "bg.fa")},
    str(t / "fg.fa"),
    t / "report",
    genome_fasta=str(t / "genome.fa"),
    verbose=False,
)
names = sorted(Path(path).name for path in result["outputs"].values())
print("gimme.denovo.html" in names, "stats.random.txt" in names)
print([n for n in names if n.endswith("_roc.random.png")])
print([n for n in names if n.endswith("_histogram.svg")])
```

- **Example output:**

```text
True True
['MA0001.1 AR1_roc.random.png', 'MA0001.2 AR1b_roc.random.png', 'MA0002.1 GATA_roc.random.png']
['MA0001.1 AR1_histogram.svg', 'MA0001.2 AR1b_histogram.svg', 'MA0002.1 GATA_histogram.svg']
```

## 30) `dnakit.sequence.gimmemotifs_diff_motifs()`

- **Purpose:** Compare motif frequencies of several FASTA sets with a background FASTA and plot the differences.
- **Method:** Follow `gimme diff`: scan every motif with the cutoff and count the sequences with a hit in the background and in each foreground set; frequencies are `(count + 0.01) / n` for a set of `n` sequences. Set names are the FASTA file stems. Then call `gimmemotifs_diff_plot()` with the tables. Upstream compares `log2(frequency / background frequency)` with `min_enrichment`, so the default 3 means an eight-fold enrichment, and it prints the kept motifs to standard output.
- **API:** `gimmemotifs_diff_motifs(foreground_fastas: 'Sequence[str]', background_fasta: 'str', pfm_file: 'str', output_path: 'str', *, cutoff: 'float' = 0.9, min_enrichment: 'float' = 3, min_frequency: 'float' = 0.01, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package `gimmemotifs` (extra `dnakit[bio-gimmemotifs]`), which uses matplotlib and scipy.
- **Input:** `foreground_fastas` is a non-empty list of FASTA paths with distinct file stems. `background_fasta` and `pfm_file` are existing files. `output_path` must not exist. `cutoff` is a score cutoff (at least 0). `min_enrichment` and `min_frequency` are the plot thresholds.
- **Output:** The fields of `gimmemotifs_diff_plot()` (`plot_written`, `motif_ids`, `sample_names`, `frequencies`, `counts`, `background_frequencies`, `background_counts`) plus `cutoff`, `background_size` and `foreground_sizes`.
- **Example code:**

```python
import contextlib
import io
import os
import random
import tempfile
from pathlib import Path

t = Path(tempfile.mkdtemp())
os.environ["XDG_CACHE_HOME"] = str(t / "cache")

from dnakit.sequence import gimmemotifs_diff_motifs

rng = random.Random(1)


def seq(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def plant(site, n):
    out = []
    for _ in range(n):
        s = seq(100)
        out.append(s[:40] + site + s[40 + len(site):])
    return out


sets = {"setA": plant("TGACAC", 30), "setB": plant("GCATAC", 30), "bg": [seq(100) for _ in range(60)]}
for name, seqs in sets.items():
    (t / f"{name}.fa").write_text("".join(f">{name}{i}\n{s}\n" for i, s in enumerate(seqs)))
with contextlib.redirect_stdout(io.StringIO()):
    result = gimmemotifs_diff_motifs(
        [str(t / "setA.fa"), str(t / "setB.fa")],
        str(t / "bg.fa"),
        "motifs.pfm",
        t / "diff.png",
        min_enrichment=1,
        verbose=False,
    )
print(result["plot_written"], result["background_size"], result["foreground_sizes"])
for motif, row in zip(result["motif_ids"], result["counts"]):
    print(motif, row)
```

- **Example output:**

```text
True 60 {'setA': 30, 'setB': 30}
MA0001.1 AR1 [30.0, 4.0]
MA0002.1 GATA [6.0, 30.0]
MA0001.2 AR1b [30.0, 3.0]
```

## 31) `dnakit.sequence.gimmemotifs_diff_plot()`

- **Purpose:** Plot motif frequencies of several sequence sets against a background.
- **Method:** Check the table shapes and call GimmeMotifs `diff_plot()`. Upstream keeps a motif only if `log2(frequency / background frequency)` exceeds `min_enrichment` in some set, a frequency exceeds `min_frequency`, the spread of the log enrichment over the sets exceeds `min_difference`, and some count is above 2. It draws nothing when no motif passes (`plot_written` is False), and with one set the spread is 0, so no plot is made at the default `min_difference`. Upstream also prints the kept motifs to standard output.
- **API:** `gimmemotifs_diff_plot(motif_ids: 'Sequence[str]', pfm_file: 'str', sample_names: 'Sequence[str]', frequencies: 'Sequence[Sequence[float]]', counts: 'Sequence[Sequence[float]]', background_frequencies: 'Sequence[float]', background_counts: 'Sequence[float]', output_path: 'str', *, min_difference: 'float' = 0, min_enrichment: 'float' = 3, min_frequency: 'float' = 0.01, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package `gimmemotifs` (extra `dnakit[bio-gimmemotifs]`), which uses matplotlib.
- **Input:** `motif_ids` are ids in `pfm_file`. `frequencies` and `counts` have one row per motif and one column per entry of `sample_names`; `background_frequencies` and `background_counts` have one value per motif. `output_path` must not exist. The thresholds are numbers (`min_frequency` at least 0).
- **Output:** `plot_written`, the input tables as given (`motif_ids`, `sample_names`, `frequencies`, `counts`, `background_frequencies`, `background_counts`) and `outputs` with `plot` when it was written.
- **Example code:**

```python
import contextlib
import io
import tempfile
from pathlib import Path

from dnakit.sequence import gimmemotifs_diff_plot

t = Path(tempfile.mkdtemp())
with contextlib.redirect_stdout(io.StringIO()):
    result = gimmemotifs_diff_plot(
        ["MA0001.1 AR1", "MA0002.1 GATA"],
        "motifs.pfm",
        ["setA", "setB"],
        [[0.8, 0.1], [0.1, 0.7]],
        [[24, 3], [3, 21]],
        [0.1, 0.1],
        [6, 6],
        t / "diff.png",
        min_enrichment=1,
        verbose=False,
    )
    single = gimmemotifs_diff_plot(
        ["MA0001.1 AR1"], "motifs.pfm", ["setA"], [[0.8]], [[24]], [0.1], [6], t / "one.png",
        min_enrichment=1, verbose=False,
    )
print(result["plot_written"], sorted(result["outputs"]), single["plot_written"], single["outputs"])
```

- **Example output:**

```text
True ['plot'] False {}
```

## 32) `dnakit.sequence.gimmemotifs_factor2orthogroups()`

- **Purpose:** Return the orthogroups of a transcription factor in an ortholog database.
- **Method:** Call GimmeMotifs `factor2orthogroups()` with the `strict` strategy, which matches the factor name or id (case-insensitively) against the genes of the reference assemblies in the SQLite database. The medium and lenient strategies look names up on mygene.info, so they are refused. Upstream caches results per factor, references and database path for the life of the process.
- **API:** `gimmemotifs_factor2orthogroups(factor: 'str', references: 'Sequence[str]', database_path: 'str', *, strategy: 'str' = 'strict', verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package `gimmemotifs` (extra `dnakit[bio-gimmemotifs]`) and a database made by `gimmemotifs_load_orthogroups_in_db()`.
- **Input:** `factor` is a non-empty name without a quote. `references` is a non-empty list of assembly names in the database. `database_path` is an existing SQLite file. `strategy` must be `"strict"`.
- **Output:** `factor` and `orthogroups`, a sorted list of orthogroup ids (empty when the factor is not found or is on the upstream blacklist).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import gimmemotifs_factor2orthogroups, gimmemotifs_load_orthogroups_in_db

t = Path(tempfile.mkdtemp())
(t / "results" / "Orthogroups").mkdir(parents=True)
(t / "results" / "Orthogroups" / "Orthogroups.tsv").write_text(
    "Orthogroup\thuman.pep\tfish.pep\n"
    "OG0000000\tGATA1|ENSG1\tgata1a|ENSDARG1, gata1b|ENSDARG2\n"
    "OG0000001\tTP53|ENSG2\ttp53|ENSDARG3\n"
)
(t / "results" / "Orthogroups" / "Orthogroups_UnassignedGenes.tsv").write_text(
    "Orthogroup\thuman.pep\tfish.pep\nOG0000002\t.|ENSG9\t\n"
)
database = t / "orthologs.sqlite"
gimmemotifs_load_orthogroups_in_db(database, ["human", "fish"], str(t / "results"), verbose=False)
for factor in ("gata1", "ENSG2", "unknown"):
    result = gimmemotifs_factor2orthogroups(factor, ["human"], str(database), verbose=False)
    print(result["factor"], result["orthogroups"])
```

- **Example output:**

```text
gata1 ['OG0000000']
ENSG2 ['OG0000001']
unknown []
```

## 33) `dnakit.sequence.gimmemotifs_filter_significant_motifs()`

- **Purpose:** Keep the de novo motifs whose statistics pass minimum values.
- **Method:** Read the motifs, key the statistics of each motif as `<id>_<consensus>` (the key GimmeMotifs uses), and call GimmeMotifs `filter_significant_motifs()` for the chosen background type. The significant motifs are written as PFM text. A motif is kept only if every statistic in `metrics` reaches its minimum.
- **API:** `gimmemotifs_filter_significant_motifs(motif_file: 'str', output_path: 'str', *, background_type: 'str', stats: 'Mapping[str, Mapping[str, Mapping[str, float]]]', metrics: 'Mapping[str, float] | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package `gimmemotifs` (extra `dnakit[bio-gimmemotifs]`).
- **Input:** `motif_file` is an existing PFM file. `output_path` must not exist. `background_type` is a non-empty string such as `random`. `stats` maps a motif id to a background type to statistics (`max_enrichment`, `roc_auc`, `enr_at_fpr`, and others); motifs without an entry are dropped. `metrics` maps statistic names to minimum values; None uses the upstream defaults (`max_enrichment` 3, `roc_auc` 0.55, `enr_at_fpr` 0.55).
- **Output:** `significant_motifs` (the kept ids), `background_type` and `outputs` with `motifs`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import gimmemotifs_filter_significant_motifs

out = Path(tempfile.mkdtemp()) / "significant.pfm"
stats = {
    "MA0001.1 AR1": {"random": {"max_enrichment": 5.0, "roc_auc": 0.8, "enr_at_fpr": 0.9}},
    "MA0002.1 GATA": {"random": {"max_enrichment": 1.0, "roc_auc": 0.5, "enr_at_fpr": 0.2}},
}
result = gimmemotifs_filter_significant_motifs(
    "motifs.pfm", out, background_type="random", stats=stats, verbose=False
)
print(result["significant_motifs"], out.read_text().count(">"))
```

- **Example output:**

```text
['MA0001.1 AR1'] 1
```

## 34) `dnakit.sequence.gimmemotifs_load_orthogroups_in_db()`

- **Purpose:** Store OrthoFinder orthogroups in a new SQLite database.
- **Method:** Check that the two OrthoFinder tables exist and have a `<assembly>.pep` column for every assembly, then call GimmeMotifs `load_orthogroups_in_db()`. It creates the tables `orthogroups`, `assemblies` and `genes` and indexes on gene names and ids. Genes are written as `gene_name|gene_id` (`.` for a missing name). The function counts the rows afterwards.
- **API:** `gimmemotifs_load_orthogroups_in_db(database_path: 'str', assemblies: 'Sequence[str]', orthofinder_results_dir: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package `gimmemotifs` (extra `dnakit[bio-gimmemotifs]`) and pandas. OrthoFinder itself is not run.
- **Input:** `database_path` must not exist (upstream would replace it). `assemblies` is a non-empty list of names. `orthofinder_results_dir` contains `Orthogroups/Orthogroups.tsv` and `Orthogroups/Orthogroups_UnassignedGenes.tsv`.
- **Output:** `orthogroup_count` (including the `UNASSIGNED` group), `gene_counts` per assembly and `outputs` with `database`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import gimmemotifs_load_orthogroups_in_db

t = Path(tempfile.mkdtemp())
(t / "results" / "Orthogroups").mkdir(parents=True)
(t / "results" / "Orthogroups" / "Orthogroups.tsv").write_text(
    "Orthogroup\thuman.pep\tfish.pep\n"
    "OG0000000\tGATA1|ENSG1\tgata1a|ENSDARG1, gata1b|ENSDARG2\n"
    "OG0000001\tTP53|ENSG2\ttp53|ENSDARG3\n"
)
(t / "results" / "Orthogroups" / "Orthogroups_UnassignedGenes.tsv").write_text(
    "Orthogroup\thuman.pep\tfish.pep\nOG0000002\t.|ENSG9\t\n"
)
result = gimmemotifs_load_orthogroups_in_db(
    t / "orthologs.sqlite", ["human", "fish"], str(t / "results"), verbose=False
)
print(result["orthogroup_count"], result["gene_counts"], sorted(result["outputs"]))
```

- **Example output:**

```text
3 {'human': 3, 'fish': 3} ['database']
```

## 35) `dnakit.sequence.gimmemotifs_make_motif2factors()`

- **Purpose:** Write a motif-to-factor table for a new assembly from an ortholog database.
- **Method:** For each motif, find the orthogroups of its factors in the `database_references` assemblies, and write every gene of `new_reference` in those orthogroups as a factor of the motif (GimmeMotifs `make_motif2factors()`). A motif without any ortholog gets the row `NO ORTHOLOGS FOUND`. The PFM file of the motifs is written next to the table. Only the `strict` strategy is allowed because the others use the network.
- **API:** `gimmemotifs_make_motif2factors(output_prefix: 'str', *, motif_file: 'str', motif_factors: 'Mapping[str, Sequence[str]]', database_path: 'str', new_reference: 'str', database_references: 'Sequence[str]', strategy: 'str' = 'strict', verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package `gimmemotifs` (extra `dnakit[bio-gimmemotifs]`) and an ortholog database from `gimmemotifs_load_orthogroups_in_db()`.
- **Input:** `output_prefix` is a path prefix; `<prefix>.motif2factors.txt` and `<prefix>.pfm` must not exist. `motif_file` is a PFM file. `motif_factors` maps motif ids to lists of factor names. `database_path` is an existing database; `new_reference` and `database_references` are assembly names in it. `strategy` must be `"strict"`.
- **Output:** `motif_count` and `outputs` with `motif2factors` and `pfm`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import gimmemotifs_load_orthogroups_in_db, gimmemotifs_make_motif2factors

t = Path(tempfile.mkdtemp())
(t / "results" / "Orthogroups").mkdir(parents=True)
(t / "results" / "Orthogroups" / "Orthogroups.tsv").write_text(
    "Orthogroup\thuman.pep\tfish.pep\n"
    "OG0000000\tGATA1|ENSG1\tgata1a|ENSDARG1, gata1b|ENSDARG2\n"
    "OG0000001\tTP53|ENSG2\ttp53|ENSDARG3\n"
)
(t / "results" / "Orthogroups" / "Orthogroups_UnassignedGenes.tsv").write_text(
    "Orthogroup\thuman.pep\tfish.pep\nOG0000002\t.|ENSG9\t\n"
)
database = t / "orthologs.sqlite"
gimmemotifs_load_orthogroups_in_db(database, ["human", "fish"], str(t / "results"), verbose=False)
result = gimmemotifs_make_motif2factors(
    str(t / "fish"),
    motif_file="motifs.pfm",
    motif_factors={"MA0001.1 AR1": ["GATA1"], "MA0002.1 GATA": ["TP53"], "MA0001.2 AR1b": ["XYZ"]},
    database_path=str(database),
    new_reference="fish",
    database_references=["human"],
    verbose=False,
)
print(result["motif_count"], sorted(result["outputs"]))
for line in sorted((t / "fish.motif2factors.txt").read_text().splitlines()[1:]):
    print(line.replace("\t", " | "))
```

- **Example output:**

```text
3 ['motif2factors', 'pfm']
MA0001.1 AR1 | gata1a | Orthologs | N
MA0001.1 AR1 | gata1b | Orthologs | N
MA0001.2 AR1b | NO ORTHOLOGS FOUND | Orthologs | N
MA0002.1 GATA | tp53 | Orthologs | N
```

## 36) `dnakit.sequence.gimmemotifs_motif2factor_from_orthologs()`

- **Purpose:** Build motif-to-factor tables for new assemblies from local genomes and annotations by ortholog inference.
- **Method:** Check that `gffread` and `orthofinder` are on PATH, then, for every assembly, copy the genome and GTF into `intermediate/annotated/<name>/` and run GimmeMotifs `annot2primpep()`; run OrthoFinder on the peptide folder; load the orthogroups into an SQLite database; and call `make_motif2factors()` for each name in `new_references` with the `strict` strategy. Nothing is downloaded.
- **API:** `gimmemotifs_motif2factor_from_orthologs(genomes: 'Mapping[str, Mapping[str, str]]', database_references: 'Sequence[str]', new_references: 'Sequence[str]', motif_file: 'str', motif_factors: 'Mapping[str, Sequence[str]]', output_dir: 'str', *, threads: 'int' = 1, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package `gimmemotifs` (extra `dnakit[bio-gimmemotifs]`) and the external programs `gffread` and `orthofinder` on PATH (for example `conda install -c bioconda gffread orthofinder`). The example below runs where they are missing and shows the error; neither program is installed in the environment used to verify this page.
- **Input:** `genomes` maps each assembly name to `{"fasta": path, "annotation": path}`. `database_references` and `new_references` are lists of names that are keys of `genomes`. `motif_file` is a PFM file and `motif_factors` maps motif ids to factor names. `output_dir` must not exist. `threads` is a positive integer.
- **Output:** `orthofinder_results` (the OrthoFinder results folder) and `outputs` with `<reference>_motif2factors` for each new reference, in `output_dir/<reference>.<motif file name>.motif2factors.txt`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.sequence import gimmemotifs_motif2factor_from_orthologs

out = Path(tempfile.mkdtemp()) / "motif2factor"
genomes = {name: {"fasta": "seqs.fa", "annotation": "motifs.pfm"} for name in ("human", "fish")}
try:
    gimmemotifs_motif2factor_from_orthologs(
        genomes, ["human"], ["fish"], "motifs.pfm", {"MA0001.1 AR1": ["GATA1"]}, out, verbose=False
    )
except BackendUnavailableError as error:
    print(type(error).__name__, str(error).split(" Hint")[0])
print(out.exists())
```

- **Example output:**

```text
BackendUnavailableError [BACKEND_UNAVAILABLE] The gffread program was not found on PATH.
False
```

Usage with the dependency installed:

```python
from dnakit.sequence import gimmemotifs_motif2factor_from_orthologs

result = gimmemotifs_motif2factor_from_orthologs(
    {
        "human": {"fasta": "human.fa", "annotation": "human.gtf"},
        "mouse": {"fasta": "mouse.fa", "annotation": "mouse.gtf"},
    },
    database_references=["human"],
    new_references=["mouse"],
    motif_file="motifs.pfm",
    motif_factors={"MA0001.1 AR1": ["GATA1"]},
    output_dir="motif2factor_out",
    threads=4,
)
print(result["outputs"])
```

## 37) `dnakit.sequence.gimmemotifs_motif_tree_plot()`

- **Purpose:** Draw a motif clustering tree whose branches are coloured by values per motif.
- **Method:** Check that ete3 can be imported (before anything is written), validate the Newick text and the per-motif values, build a pandas table with one column per motif, and call GimmeMotifs `motif_tree_plot()`. The colour of a branch is the mean of the values of the leaves below it.
- **API:** `gimmemotifs_motif_tree_plot(output_path: 'str', newick_tree: 'str', leaf_values: 'Mapping[str, Sequence[float]]', *, circle: 'bool' = True, vmin: 'float | None' = None, vmax: 'float | None' = None, dpi: 'int' = 300, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages `gimmemotifs` and `ete3` (both in extra `dnakit[bio-gimmemotifs]`). ete3 is not installed in the environment used to verify this page, so the example shows the error that is raised before any file is written.
- **Input:** `output_path` must not exist. `newick_tree` is Newick text ending with `;` whose leaves are motif ids. `leaf_values` maps each motif id to a list of numbers, all of the same length. `circle` selects a circular layout. `vmin` and `vmax` set the colour range. `dpi` is a positive integer.
- **Output:** `leaf_count` and `outputs` with `plot`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.sequence import gimmemotifs_motif_tree_plot

out = Path(tempfile.mkdtemp()) / "tree.png"
try:
    gimmemotifs_motif_tree_plot(
        out, "((m1,m2),m3);", {"m1": [1.0, 2.0], "m2": [1.5, 2.5], "m3": [-1.0, 0.0]}, verbose=False
    )
except BackendUnavailableError as error:
    print(type(error).__name__, str(error).split(" Hint")[0])
print(out.exists())
```

- **Example output:**

```text
BackendUnavailableError [BACKEND_UNAVAILABLE] ete3 could not be imported: No module named 'ete3'
False
```

Usage with the dependency installed:

```python
from dnakit.sequence import gimmemotifs_motif_tree_plot

result = gimmemotifs_motif_tree_plot(
    "tree.png",
    "((m1,m2),m3);",
    {"m1": [1.0, 2.0], "m2": [1.5, 2.5], "m3": [-1.0, 0.0]},
    circle=True,
    dpi=150,
)
print(result["leaf_count"], result["outputs"]["plot"])
```

## 38) `dnakit.sequence.gimmemotifs_roc_html_report()`

- **Purpose:** Write the HTML report of a motif enrichment table.
- **Method:** Read the table, add the Benjamini-Hochberg corrected P-value and its negative log10, keep the motifs with a corrected P-value at or below `threshold`, draw their logos into `<output_dir>/logos`, and write the HTML file (GimmeMotifs `roc_html_report()`).
- **API:** `gimmemotifs_roc_html_report(output_dir: 'str', input_table: 'str', pfm_file: 'str', output_name: 'str', *, threshold: 'float' = 0.01, use_motifs: 'Sequence[str] | None' = None, link_matches: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package `gimmemotifs` (extra `dnakit[bio-gimmemotifs]`) and pandas.
- **Input:** `output_dir` must not exist. `input_table` is the tab-separated enrichment table of a finished run (motif ids in the first column and the columns `# matches`, `% matches input`, `# matches background`, `% matches background`, `P-value`, `log10 P-value`, `ROC AUC`, `PR AUC`, `Enr. at 1% FPR`, `Recall at 10% FDR`). `pfm_file` holds the motifs. `output_name` is a plain file name. `threshold` is a number in [0, 1]. `use_motifs` restricts the ids. `link_matches` links the match counts to BED files of a scan.
- **Output:** `threshold` and `outputs`: `report`, and every file below the directory (the HTML file and the logos) by relative path.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import gimmemotifs_roc_html_report

t = Path(tempfile.mkdtemp())
header = (
    "Motif\t# matches\t% matches input\t# matches background\t% matches background\t"
    "P-value\tlog10 P-value\tROC AUC\tPR AUC\tEnr. at 1% FPR\tRecall at 10% FDR\n"
)
(t / "roc.txt").write_text(
    header
    + "MA0001.1 AR1\t30\t75.0\t8\t10.0\t1e-12\t12.0\t0.90\t0.85\t6.0\t0.50\n"
    + "MA0002.1 GATA\t20\t50.0\t10\t12.5\t1e-04\t4.0\t0.70\t0.65\t3.0\t0.20\n"
    + "MA0001.2 AR1b\t28\t70.0\t8\t10.0\t1e-10\t10.0\t0.88\t0.80\t5.0\t0.45\n"
)
result = gimmemotifs_roc_html_report(
    t / "report", str(t / "roc.txt"), "motifs.pfm", "report.html", threshold=1e-6, verbose=False
)
print(sorted(result["outputs"]))
html = (t / "report" / "report.html").read_text()
print("MA0001.1 AR1" in html, "MA0002.1 GATA" in html)
```

- **Example output:**

```text
['logos/MA0001_1_AR1.png', 'logos/MA0001_2_AR1b.png', 'report', 'report.html']
True False
```

## 39) `dnakit.sequence.gimmemotifs_roc_plots()`

- **Purpose:** Draw one ROC curve per motif and background set.
- **Method:** For every motif and background FASTA, compute the ROC values of the motif on the foreground and the background (`calc_stats()` with the z-score and GC options off, because the default options hash the motif with a `str`, which xxhash 3 and later reject) and draw the curve with the GimmeMotifs report function. The function checks that one file per pair was written.
- **API:** `gimmemotifs_roc_plots(pfm_file: 'str', foreground_fasta: 'str', background_fastas: 'Mapping[str, str]', output_dir: 'str', *, genome_fasta: 'str', verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package `gimmemotifs` (extra `dnakit[bio-gimmemotifs]`), which uses matplotlib.
- **Input:** `pfm_file`, `foreground_fasta` and `genome_fasta` are existing files. `background_fastas` maps background names to FASTA paths. `output_dir` must not exist.
- **Output:** `motif_count`, `background_count`, and `outputs`, a mapping from `<motif id>_roc.<background>.png` to its path (files are in `output_dir/images`).
- **Example code:**

```python
import os
import random
import tempfile
from pathlib import Path

t = Path(tempfile.mkdtemp())
os.environ["XDG_CACHE_HOME"] = str(t / "cache")

from dnakit.sequence import gimmemotifs_roc_plots

rng = random.Random(1)


def seq(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


foreground = [seq(200) for _ in range(40)]
foreground = [s[:50] + "TGACAC" + s[56:] if i % 2 == 0 else s for i, s in enumerate(foreground)]
(t / "fg.fa").write_text("".join(f">s{i}\n{s}\n" for i, s in enumerate(foreground)))
(t / "bg.fa").write_text("".join(f">b{i}\n{seq(200)}\n" for i in range(80)))
(t / "genome.fa").write_text("".join(f">{c}\n{seq(5000)}\n" for c in ("chr1", "chr2")))
result = gimmemotifs_roc_plots(
    "motifs.pfm",
    str(t / "fg.fa"),
    {"random": str(t / "bg.fa")},
    t / "roc",
    genome_fasta=str(t / "genome.fa"),
    verbose=False,
)
print(result["motif_count"], result["background_count"])
print(sorted(result["outputs"]))
```

- **Example output:**

```text
3 1
['MA0001.1 AR1_roc.random.png', 'MA0001.2 AR1b_roc.random.png', 'MA0002.1 GATA_roc.random.png']
```

## 40) `dnakit.sequence.gimmemotifs_select_nonredundant_motifs()`

- **Purpose:** Select a non-redundant set of motifs from a ROC report and two motif score tables.
- **Method:** Call GimmeMotifs `select_nonredundant_motifs()`. Motifs with `Enr. at 1% FPR` below 2 are dropped, the rest are ranked by the mean rank of four report columns, motifs similar to a better one (seqcor of at least 0.7) are removed, and, if at least two motifs remain, recursive feature elimination with an L1 logistic regression keeps the motifs that improve the average precision.
- **API:** `gimmemotifs_select_nonredundant_motifs(roc_report: 'str', pfm_file: 'str', foreground_table: 'str', background_table: 'str', *, tolerance: 'float' = 0.001, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package `gimmemotifs` (extra `dnakit[bio-gimmemotifs]`), which uses scikit-learn.
- **Input:** `roc_report` is a tab-separated table with motif ids in the first column and the columns `ROC AUC`, `PR AUC`, `Enr. at 1% FPR` and `Recall at 10% FDR`. `pfm_file` holds the motifs. `foreground_table` and `background_table` are tab-separated score tables with one row per sequence and one column per motif id. `tolerance` is a number of at least 0.
- **Output:** `selected_motifs`, the ids of the selected motifs.
- **Example code:**

```python
import random
import tempfile
from pathlib import Path

from dnakit.sequence import gimmemotifs_select_nonredundant_motifs

t = Path(tempfile.mkdtemp())
header = (
    "Motif\t# matches\t% matches input\t# matches background\t% matches background\t"
    "P-value\tlog10 P-value\tROC AUC\tPR AUC\tEnr. at 1% FPR\tRecall at 10% FDR\n"
)
(t / "roc.txt").write_text(
    header
    + "MA0001.1 AR1\t30\t75.0\t8\t10.0\t1e-12\t12.0\t0.90\t0.85\t6.0\t0.50\n"
    + "MA0002.1 GATA\t20\t50.0\t10\t12.5\t1e-04\t4.0\t0.70\t0.65\t3.0\t0.20\n"
    + "MA0001.2 AR1b\t28\t70.0\t8\t10.0\t1e-10\t10.0\t0.88\t0.80\t5.0\t0.45\n"
)
rng = random.Random(3)


def table(n, shift):
    rows = ["\tMA0001.1 AR1\tMA0002.1 GATA\tMA0001.2 AR1b"]
    for i in range(n):
        a = rng.gauss(shift, 1)
        rows.append(f"s{i}\t{a:.3f}\t{rng.gauss(shift / 2, 1):.3f}\t{a + rng.gauss(0, 0.1):.3f}")
    return "\n".join(rows) + "\n"


(t / "fg.tsv").write_text(table(60, 2))
(t / "bg.tsv").write_text(table(60, 0))
result = gimmemotifs_select_nonredundant_motifs(
    str(t / "roc.txt"), "motifs.pfm", str(t / "fg.tsv"), str(t / "bg.tsv"), verbose=False
)
print(result["selected_motifs"])
```

- **Example output:**

```text
['MA0001.1 AR1']
```
