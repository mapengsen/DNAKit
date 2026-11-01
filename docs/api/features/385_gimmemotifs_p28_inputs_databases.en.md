# GimmeMotifs FASTA tools, statistics, thresholds, match counts and motif databases

These functions in `dnakit.regulation` wrap the FASTA, peak-input, statistics, threshold, match-counting and motif-database routines of [GimmeMotifs](https://github.com/vanheeringen-lab/gimmemotifs) 0.18.4. Each one calls the GimmeMotifs routine on local files or explicit lists and returns a JSON-serialisable dictionary with `status`, the result and a `provenance` entry. Genome arguments are FASTA paths (no genomepy download); output files must not exist. The two motif-database functions use the network.

GimmeMotifs is an optional dependency (`pip install "dnakit[bio-gimmemotifs]"`; the BED and narrowPeak input preparation also needs `genomepy`). The examples import from `dnakit.regulation.gimmemotifs_p28_data_functions`; after the area package is wired the same names are available from `dnakit.regulation`.

They cover `Fasta.hardmask`, `Fasta.get_random`, `utils.as_seqdict`, `stats.rank_motifs`, `stats.write_stats`, `Scanner.set_threshold` (and the `gimme threshold` command), `Scanner.count` and `Scanner.total_count`, `utils.narrowpeak_to_bed`, `motif.denovo.prepare_denovo_input_fa`, `_bed` and `_narrowpeak`, and the `download` and `annotate_factors` methods of the JASPAR, HOMER, HOCOMOCO, ENCODE, SwissRegulon and IMAGE motif databases. Rows already served by older functions are `rocmetrics` (`motif_score_metric`), `preprocessing.peakfile2bedfile` (`normalize_peak_widths`), and `utils.get_jaspar_motif_info` (`gimmemotifs_jaspar_motif_info`).

## Example data

The first example reads the two-record soft-masked FASTA file `soft.fa`. The other examples write their inputs into a temporary directory with a small deterministic generator (`make_fasta`) or pass them directly.

**Example file:** `soft.fa`

```text
>a
ACGTacgtnNAC
>b
ggggCCCC
```

## 1) `dnakit.regulation.gimmemotifs_hardmask_fasta()`

- **Purpose:** Replace the lower-case (soft-masked) bases of a FASTA file by N.
- **Method:** Read the file with `gimmemotifs.fasta.Fasta`, call `Fasta.hardmask` (lower-case a, c, g, t and n become `N`, upper-case letters are kept) and write the records with `Fasta.writefasta`, without line wrapping. The output path must be new.
- **API:** `dnakit.regulation.gimmemotifs_hardmask_fasta(input_fasta, output_fasta)`.
- **Input:** A FASTA file that starts with `>` and a new output path.
- **Output:** `masked_bases`, `n_sequences` and the written file in `outputs`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.regulation.gimmemotifs_p28_data_functions import gimmemotifs_hardmask_fasta

out = Path(tempfile.mkdtemp()) / "hard.fa"
result = gimmemotifs_hardmask_fasta("soft.fa", out)
print(result["masked_bases"], result["n_sequences"])
print(out.read_text().strip())
```

- **Example output:**

```text
9 2
>a
ACGTNNNNNNAC
>b
NNNNCCCC
```

## 2) `dnakit.regulation.gimmemotifs_random_fasta_subset()`

- **Purpose:** Draw random whole records or random fixed-length windows from a FASTA file.
- **Method:** Call `Fasta.get_random`. Without `length` it samples `n_sequences` records without replacement and keeps their names. With `length` it visits the records in random order and takes one random window of `length` bases from each record that is long enough, naming them `random1`, `random2`, ...; it fails when fewer than `n_sequences` records are long enough. `seed` seeds Python's random generator for the call (the previous state is restored). The subset can be written to a new FASTA file.
- **API:** `dnakit.regulation.gimmemotifs_random_fasta_subset(input_fasta, n_sequences, *, length=None, seed=None, output_fasta=None)`.
- **Input:** A FASTA file, the number of sequences, an optional window `length`, `seed` and `output_fasta`.
- **Output:** `sequences` (name to sequence), `n_sequences` and the optional file in `outputs`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.regulation.gimmemotifs_p28_data_functions import gimmemotifs_random_fasta_subset

def make_fasta(path, count, length, seed):
    state, rows = seed, []
    for index in range(count):
        bases = []
        for _ in range(length):
            state = (state * 1103515245 + 12345) % 2**31
            bases.append("ACGT"[(state >> 16) % 4])
        rows.append(f">r{index}\n{''.join(bases)}\n")
    path.write_text("".join(rows))
    return path

fasta = make_fasta(Path(tempfile.mkdtemp()) / "r.fa", 20, 60, 7)
records = gimmemotifs_random_fasta_subset(fasta, 4, seed=3)
print(len(records["sequences"]), all(len(s) == 60 for s in records["sequences"].values()))
again = gimmemotifs_random_fasta_subset(fasta, 4, seed=3)
print(records["sequences"] == again["sequences"])
windows = gimmemotifs_random_fasta_subset(fasta, 3, length=25, seed=3)
print(list(windows["sequences"]), [len(s) for s in windows["sequences"].values()])
```

- **Example output:**

```text
4 True
True
['random1', 'random2', 'random3'] [25, 25, 25]
```

## 3) `dnakit.regulation.gimmemotifs_sequence_dict()`

- **Purpose:** Read sequences into a dictionary keyed by sequence name.
- **Method:** Call `gimmemotifs.utils.as_seqdict`. A FASTA file gives the complete header lines as names; a list of DNA strings gives the names `1`, `2`, ...; a BED file or a list of `chrom:start-end` regions is cut out of `genome_fasta` (read as a genomepy genome from its FASTA path). With `minsize` every sequence must be at least that long, otherwise the call fails.
- **API:** `dnakit.regulation.gimmemotifs_sequence_dict(source, *, genome_fasta=None, minsize=None)`.
- **Input:** A FASTA, BED or region file path, or a list of DNA strings or regions; optional `genome_fasta` and `minsize`.
- **Output:** `sequences` (name to sequence) and `n_sequences`.
- **Example code:**

```python
from dnakit.regulation.gimmemotifs_p28_data_functions import gimmemotifs_sequence_dict

print(gimmemotifs_sequence_dict("soft.fa")["sequences"])
print(gimmemotifs_sequence_dict(["ACGT", "GGAA"])["sequences"])
try:
    gimmemotifs_sequence_dict("soft.fa", minsize=10)
except Exception as error:
    print(type(error).__name__, "shorter than 10" in str(error))
```

- **Example output:**

```text
{'a': 'ACGTacgtnNAC', 'b': 'ggggCCCC'}
{'1': 'ACGT', '2': 'GGAA'}
BackendExecutionError True
```

## 4) `dnakit.regulation.gimmemotifs_rank_motif_stats()`

- **Purpose:** Rank motifs by their mean rank over several enrichment metrics.
- **Method:** Call `gimmemotifs.stats.rank_motifs` on a table `{motif: {background: {metric: value}}}`, the shape that the motif-enrichment statistics produce. For each metric the motif values are averaged over the backgrounds and ranked with `scipy.stats.rankdata` (rank 1 is the lowest value, ties get the average rank); the result is the mean of these ranks, so a higher number is a better motif.
- **API:** `dnakit.regulation.gimmemotifs_rank_motif_stats(stats, *, metrics=("roc_auc", "recall_at_fdr"))`.
- **Input:** The nested `stats` mapping with the same backgrounds and metrics for every motif, and the list of `metrics` to use.
- **Output:** `mean_rank` per motif, `order` (best first), `metrics` and `backgrounds`.
- **Example code:**

```python
from dnakit.regulation.gimmemotifs_p28_data_functions import gimmemotifs_rank_motif_stats

stats = {
    "m_1": {"bg": {"roc_auc": 0.9, "recall_at_fdr": 0.5}},
    "m_2": {"bg": {"roc_auc": 0.6, "recall_at_fdr": 0.7}},
    "m_3": {"bg": {"roc_auc": 0.7, "recall_at_fdr": 0.1}},
}
result = gimmemotifs_rank_motif_stats(stats)
print(result["mean_rank"], result["order"])
print(gimmemotifs_rank_motif_stats(stats, metrics=["recall_at_fdr"])["order"])
```

- **Example output:**

```text
{'m_1': 2.5, 'm_2': 2.0, 'm_3': 1.5} ['m_1', 'm_2', 'm_3']
['m_2', 'm_1', 'm_3']
```

## 5) `dnakit.regulation.gimmemotifs_write_motif_stats()`

- **Purpose:** Write the motif statistics table, one tab-separated file per background.
- **Method:** Call `gimmemotifs.stats.write_stats`. `output_template` is a path with one `{}` that is replaced by each background name. Each file has the optional `header` text, a `Motif` column and one column per metric in alphabetical order. As upstream does, the motif name is written without its last `_` part (`m_1` becomes `m`). Existing files are never overwritten.
- **API:** `dnakit.regulation.gimmemotifs_write_motif_stats(stats, output_template, *, header=None)`.
- **Input:** The `stats` mapping as for `gimmemotifs_rank_motif_stats`, a template with one `{}` and an optional header.
- **Output:** `outputs` (background to file path), `backgrounds` and `n_motifs`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.regulation.gimmemotifs_p28_data_functions import gimmemotifs_write_motif_stats

stats = {
    "m_1": {"bg": {"roc_auc": 0.9, "recall_at_fdr": 0.5}},
    "m_2": {"bg": {"roc_auc": 0.6, "recall_at_fdr": 0.7}},
}
folder = Path(tempfile.mkdtemp())
result = gimmemotifs_write_motif_stats(stats, str(folder / "stats.{}.txt"), header="# demo\n")
print(list(result["outputs"]))
print((folder / "stats.bg.txt").read_text().strip())
```

- **Example output:**

```text
['bg']
# demo
Motif	recall_at_fdr	roc_auc
m	0.5	0.9
m	0.7	0.6
```

## 6) `dnakit.regulation.gimmemotifs_motif_thresholds()`

- **Purpose:** Motif score thresholds at a target false-positive rate on background sequences.
- **Method:** Scan every motif of the PFM file over the background FASTA (best hit per sequence, both strands) with the GimmeMotifs scanner and take, per motif, the `1 - fpr` quantile of those best scores as the threshold (`Scanner.get_gc_thresholds`, the step behind `Scanner.set_threshold` and `gimme threshold`). GimmeMotifs resamples the background scores, so pass `seed` for reproducible values. No GC correction is applied. A threshold of `None` means the quantile reaches the maximum score of the motif, so no sequence can pass it. The thresholds are also given as the fraction of the distance between the motif's minimum and maximum score.
- **API:** `dnakit.regulation.gimmemotifs_motif_thresholds(motifs_pfm, background_fasta, *, fpr=0.01, seed=None)`.
- **Input:** A PFM file, a background FASTA file, `fpr` strictly between 0 and 1 and an optional integer `seed`.
- **Output:** `thresholds` per motif (`threshold`, `fraction`, `max_score`) and `n_background_sequences`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.regulation.gimmemotifs_p28_data_functions import gimmemotifs_motif_thresholds

def make_fasta(path, count, length, seed):
    state, rows = seed, []
    for index in range(count):
        bases = []
        for _ in range(length):
            state = (state * 1103515245 + 12345) % 2**31
            bases.append("ACGT"[(state >> 16) % 4])
        rows.append(f">r{index}\n{''.join(bases)}\n")
    path.write_text("".join(rows))
    return path

folder = Path(tempfile.mkdtemp())
background = make_fasta(folder / "bg.fa", 60, 100, 11)
(folder / "m.pfm").write_text((
    ">weak\n6 2 1 1\n1 6 2 1\n1 1 6 2\n2 1 1 6\n3 3 2 2\n"
    ">acgt\n10 0 0 0\n0 10 0 0\n0 0 10 0\n0 0 0 10\n"
))
result = gimmemotifs_motif_thresholds(folder / "m.pfm", background, fpr=0.5, seed=1)
print(result["n_background_sequences"], list(result["thresholds"]))
weak = result["thresholds"]["weak"]
print(round(weak["threshold"], 4), round(weak["fraction"], 4), round(weak["max_score"], 4))
print(result["thresholds"]["acgt"]["threshold"])
```

- **Example output:**

```text
60 ['weak', 'acgt']
2.6183 0.8544 3.7085
-0.40047363884467657
```

## 7) `dnakit.regulation.gimmemotifs_count_motif_matches()`

- **Purpose:** Count motif matches per sequence and in total.
- **Method:** Use the GimmeMotifs `Scanner`: `Scanner.count` gives the number of matches of every motif in every sequence and `Scanner.total_count` the sum over all sequences. A match is a hit above the motif's cutoff, on both strands when `scan_rc` is true, and at most `nreport` hits are counted per motif and sequence. By default the cutoff is `threshold`, the fraction of the distance between the motif's minimum and maximum score (0.95 is the GimmeMotifs default without a genome). With `fpr` and `background_fasta` the cutoff is instead the score that a fraction `fpr` of the background sequences reaches (see `gimmemotifs_motif_thresholds`).
- **API:** `dnakit.regulation.gimmemotifs_count_motif_matches(sequences_fasta, motifs_pfm, *, threshold=0.95, background_fasta=None, fpr=None, nreport=100, scan_rc=True, seed=None)`.
- **Input:** A FASTA file, a PFM file and either `threshold` or the pair `fpr` plus `background_fasta`.
- **Output:** `counts` (sequence to motif to count), `total` (motif to count), `motifs` and the cutoff mode.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.regulation.gimmemotifs_p28_data_functions import gimmemotifs_count_motif_matches

def make_fasta(path, count, length, seed):
    state, rows = seed, []
    for index in range(count):
        bases = []
        for _ in range(length):
            state = (state * 1103515245 + 12345) % 2**31
            bases.append("ACGT"[(state >> 16) % 4])
        rows.append(f">r{index}\n{''.join(bases)}\n")
    path.write_text("".join(rows))
    return path

folder = Path(tempfile.mkdtemp())
sequences = make_fasta(folder / "seq.fa", 10, 100, 21)
(folder / "m.pfm").write_text((
    ">weak\n6 2 1 1\n1 6 2 1\n1 1 6 2\n2 1 1 6\n3 3 2 2\n"
    ">acgt\n10 0 0 0\n0 10 0 0\n0 0 10 0\n0 0 0 10\n"
))
result = gimmemotifs_count_motif_matches(sequences, folder / "m.pfm", threshold=0.8)
print(result["motifs"], result["total"])
print(result["counts"]["r0"])
print(sum(sum(row.values()) for row in result["counts"].values()) == sum(result["total"].values()))
```

- **Example output:**

```text
['weak', 'acgt'] {'weak': 41, 'acgt': 14}
{'weak': 4, 'acgt': 0}
True
```

## 8) `dnakit.regulation.gimmemotifs_narrowpeak_to_bed()`

- **Purpose:** Convert a narrowPeak file to a four-column BED file, optionally as fixed-size windows around the peak summits.
- **Method:** Call `gimmemotifs.utils.narrowpeak_to_bed`. With `size=0` every peak keeps its coordinates. With a positive `size` the interval is `size` bases long and starts at `peak start + summit - size // 2`; the summit is narrowPeak column 10, and -1 (no summit) falls back to `(end - start) // 2`. The BED columns are chromosome, start, end and the signal value (narrowPeak column 7). Comment, `track` and `browser` lines are skipped. The output file must not exist.
- **API:** `dnakit.regulation.gimmemotifs_narrowpeak_to_bed(input_file, output_bed, *, size=0)`.
- **Input:** A narrowPeak file with 10 columns, a new BED path and `size`.
- **Output:** `intervals` as `[chromosome, start, end, signal]`, `n_intervals` and the file in `outputs`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.regulation.gimmemotifs_p28_data_functions import gimmemotifs_narrowpeak_to_bed

folder = Path(tempfile.mkdtemp())
peaks = folder / "p.narrowPeak"
peaks.write_text(
    "chr1\t100\t300\tp0\t0\t.\t12.5\t-1\t-1\t80\n"
    "chr2\t500\t900\tp1\t0\t.\t7.0\t-1\t-1\t-1\n"
)
print(gimmemotifs_narrowpeak_to_bed(peaks, folder / "plain.bed")["intervals"])
print(gimmemotifs_narrowpeak_to_bed(peaks, folder / "summit.bed", size=60)["intervals"])
```

- **Example output:**

```text
[['chr1', 100, 300, '12.5'], ['chr2', 500, 900, '7.0']]
[['chr1', 150, 210, '12.5'], ['chr2', 670, 730, '7.0']]
```

## 9) `dnakit.regulation.gimmemotifs_prepare_denovo_input()`

- **Purpose:** Prepare the input files of a de novo motif run: prediction, validation and localization sets.
- **Method:** `input_type="fasta"` calls `prepare_denovo_input_fa`: the records are split into `prediction.fa` (a fraction `fraction` of them, at most `abs_max`) and `validation.fa`, and `localization.fa` is a copy of the validation set. `"bed"` calls `prepare_denovo_input_bed`: the regions are resized to `size` bases around their midpoint (`input.bed`), split into `prediction.bed` and `validation.bed`, and cut out of `genome_fasta` into the matching FASTA files; the localization FASTA is extended to `lsize` bases. `"narrowpeak"` calls `prepare_denovo_input_narrowpeak`, which first converts the peaks to a BED of `size` bases around the summit (`<name>_to.bed`) and then proceeds as for BED. Unset `size`, `fraction`, `abs_max` and `lsize` take the GimmeMotifs defaults (200, 0.2, 1000 and 500). `seed` seeds the random split.
- **API:** `dnakit.regulation.gimmemotifs_prepare_denovo_input(input_file, output_dir, *, input_type="fasta", genome_fasta=None, size=None, fraction=None, abs_max=None, lsize=None, use_strand=False, seed=None)`.
- **Input:** An input file, a new or empty `output_dir`, the `input_type` and, for BED and narrowPeak, `genome_fasta`.
- **Output:** `outputs` (file name to path), `records_per_file` and the `parameters` used.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.regulation.gimmemotifs_p28_data_functions import gimmemotifs_prepare_denovo_input

def make_fasta(path, count, length, seed):
    state, rows = seed, []
    for index in range(count):
        bases = []
        for _ in range(length):
            state = (state * 1103515245 + 12345) % 2**31
            bases.append("ACGT"[(state >> 16) % 4])
        rows.append(f">r{index}\n{''.join(bases)}\n")
    path.write_text("".join(rows))
    return path

folder = Path(tempfile.mkdtemp())
genome = make_fasta(folder / "genome.fa", 1, 5000, 3)
genome.write_text(genome.read_text().replace(">r0", ">chr1"))
peaks = folder / "peaks.bed"
peaks.write_text("".join(f"chr1\t{100 + 150 * i}\t{300 + 150 * i}\tp{i}\t{100 - i}\n" for i in range(10)))
fa = make_fasta(folder / "in.fa", 20, 60, 5)

result = gimmemotifs_prepare_denovo_input(fa, folder / "from_fasta", fraction=0.5, seed=1)
print(result["records_per_file"])
result = gimmemotifs_prepare_denovo_input(
    peaks, folder / "from_bed", input_type="bed", genome_fasta=genome, size=100, fraction=0.5, lsize=140, seed=1
)
print(result["records_per_file"])
print(len((folder / "from_bed" / "prediction.fa").read_text().splitlines()[1]),
      len((folder / "from_bed" / "localization.fa").read_text().splitlines()[1]))
```

- **Example output:**

```text
{'localization.fa': 10, 'prediction.fa': 10, 'validation.fa': 10}
{'input.bed': 10, 'localization.fa': 5, 'prediction.bed': 5, 'prediction.fa': 5, 'validation.bed': 5, 'validation.fa': 5}
100 140
```

## 10) `dnakit.regulation.gimmemotifs_download_motif_database()`

- **Purpose:** Download a public motif database and write it in GimmeMotifs PFM format with its factor annotation (network).
- **Method:** Create the GimmeMotifs `MotifDb` of the chosen database and call its `download` method. `database` is `jaspar` (all CORE groups of `jaspar_version`), `homer`, `hocomoco` (v11 human and mouse), `encode`, `swissregulon` or `image`. Each PFM file (with `# Retrieved from` and `# Date` header lines) is written together with a `.motif2factors.txt` file made by `annotate_factors`. The directory is created if needed and files are overwritten as upstream does. The upstream downloads of Factorbook, CIS-BP and RSAT are not offered because their source addresses answer 404 (checked 2026-10-20). The example replaces the web request by a local text so that it runs offline.
- **API:** `dnakit.regulation.gimmemotifs_download_motif_database(database, output_dir, *, jaspar_version="2022")`.
- **Input:** A database name, an output directory and, for JASPAR, a four-digit release year.
- **Output:** `outputs` (file name to path of the files written), `date` and `database`.
- **Example code:**

```python
import tempfile
from pathlib import Path

import gimmemotifs.db as db

from dnakit.regulation.gimmemotifs_p28_data_functions import gimmemotifs_download_motif_database

HOMER = (
    ">ATGACTCATC_PU.1(ETS)/ThioMac-PU.1-ChIP-Seq(GSE21512)/Homer\t6.04\t-1782\t0\t9805\n"
    "0.419\t0.275\t0.277\t0.028\n0.001\t0.001\t0.001\t0.997\n0.010\t0.002\t0.965\t0.023\n"
    ">GGAAGTGGGG_Gfi1b(Zf)/Bone-Gfi1b-ChIP-Seq(GSE23725)/Homer\t6.1\t-12\t0\t9\n"
    "0.1\t0.1\t0.7\t0.1\n0.1\t0.7\t0.1\t0.1\n0.7\t0.1\t0.1\t0.1\n"
)


class Reply(list):
    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


# a local text stands in for the HOMER web server so that the example runs offline;
# remove the next line to download the real database
db.urlopen = lambda url: Reply(HOMER.encode().splitlines(keepends=True))

folder = Path(tempfile.mkdtemp())
result = gimmemotifs_download_motif_database("homer", folder)
print(sorted(result["outputs"]))
print((folder / "HOMER.pfm").read_text().count(">"))
print((folder / "HOMER.motif2factors.txt").read_text().strip())
```

- **Example output:**

```text
['HOMER.motif2factors.txt', 'HOMER.pfm']
2
Motif	Factor	Evidence	Curated
ATGACTCATC_PU.1(ETS)/ThioMac-PU.1-ChIP-Seq(GSE21512)/Homer_6.04	PU.1	ChIP-seq	Y
GGAAGTGGGG_Gfi1b(Zf)/Bone-Gfi1b-ChIP-Seq(GSE23725)/Homer_6.1	Gfi1b	ChIP-seq	Y
```

## 11) `dnakit.regulation.gimmemotifs_annotate_motif_factors()`

- **Purpose:** Annotate the motifs of a PFM file with transcription factors from a motif database (network).
- **Method:** Read the motifs and call `annotate_factors` of the chosen GimmeMotifs `MotifDb`. `database` is `jaspar` (looks every `MA` matrix up at JASPAR), `homer`, `hocomoco` (`hocomoco_group` `HUMAN` or `MOUSE`, which selects the annotation table), `encode`, `swissregulon` or `image`. The motif IDs must follow that database's own naming, as in the PFM file written by `gimmemotifs_download_motif_database`; motifs without a match are left out. With `output_dir` the table is also written as `<pfm name>.motif2factors.txt` (columns Motif, Factor, Evidence, Curated). The RSAT, CIS-BP and Factorbook annotations are not offered: the first is an empty stub upstream and the other two read addresses that answer 404 (checked 2026-10-20). The example replaces the web request by a local text.
- **API:** `dnakit.regulation.gimmemotifs_annotate_motif_factors(database, motifs_pfm, *, hocomoco_group="HUMAN", output_dir=None)`.
- **Input:** A database name, a PFM file, `hocomoco_group` and an optional output directory.
- **Output:** `annotation` (motif to a list of `[factor, evidence, curated]`), `n_motifs`, `n_annotated` and the optional file in `outputs`.
- **Example code:**

```python
import tempfile
from pathlib import Path

import gimmemotifs.db as db

from dnakit.regulation.gimmemotifs_p28_data_functions import gimmemotifs_annotate_motif_factors

HOMER = (
    ">ATGACTCATC_PU.1(ETS)/ThioMac-PU.1-ChIP-Seq(GSE21512)/Homer\t6.04\t-1782\t0\t9805\n"
    "0.419\t0.275\t0.277\t0.028\n0.001\t0.001\t0.001\t0.997\n0.010\t0.002\t0.965\t0.023\n"
    ">GGAAGTGGGG_Gfi1b(Zf)/Bone-Gfi1b-ChIP-Seq(GSE23725)/Homer\t6.1\t-12\t0\t9\n"
    "0.1\t0.1\t0.7\t0.1\n0.1\t0.7\t0.1\t0.1\n0.7\t0.1\t0.1\t0.1\n"
)


class Reply(list):
    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


# a local text stands in for the HOMER web server so that the example runs offline;
# remove the next line to download the real database
db.urlopen = lambda url: Reply(HOMER.encode().splitlines(keepends=True))

folder = Path(tempfile.mkdtemp())
(folder / "homer.pfm").write_text(
    ">ATGACTCATC_PU.1(ETS)/ThioMac-PU.1-ChIP-Seq(GSE21512)/Homer\n0.4 0.3 0.2 0.1\n0.1 0.1 0.1 0.7\n"
    ">NNNNNNNN_AP-1(bZIP)/ThioMac-PU.1-ChIP-Seq(GSE21512)/Homer\n0.4 0.3 0.2 0.1\n0.1 0.1 0.1 0.7\n"
)
result = gimmemotifs_annotate_motif_factors("homer", folder / "homer.pfm", output_dir=folder / "out")
print(result["n_motifs"], result["n_annotated"])
print(result["annotation"])
print(sorted(p.name for p in (folder / "out").iterdir()))
```

- **Example output:**

```text
2 1
{'ATGACTCATC_PU.1(ETS)/ThioMac-PU.1-ChIP-Seq(GSE21512)/Homer': [['PU.1', 'ChIP-seq', 'Y']]}
['homer.motif2factors.txt']
```
