# Merqury: k-mer based assembly evaluation, trio hap-mers and k-mer database tools

This group of functions in `dnakit.assembly` runs the scripts of [Merqury](https://github.com/marbl/merqury) (k-mer based evaluation of genome assemblies, built on Meryl) as typed facades: per-sequence QV, k-mer multiplicity tracks, false-duplication counts and tracks, Mash QV, k-mer size choice, 10x read splitting, solid-k-mer filtering, k-mer database set operations, histogram export, trio hap-mer derivation, hap-mer spectra and blob counts, Hamming error, read assignment and exclusion, phase blocks with switch errors, N*/NG* statistics, bed-to-bigWig conversion and the three R plots. Whole-assembly QV and completeness with copy-number spectra (`merqury.sh`) is already provided by `dnakit.assembly.evaluate_merqury`. Each function runs one script as `bash <MERQURY>/<script> <arguments>` in a new output directory, after linking the inputs into it, so the input folders are never written to and an existing output directory is refused. Results are dictionaries with `status`, `backend`, `version` (Meryl), `command`, `return_code`, `outputs`, `log_tail`, `operation`, `merqury_home` and `provenance` plus parsed values. `status` is `completed` only when the exit code is 0 and every expected file exists.

Merqury 1.4.1 (Meryl 1.4.2) is the tested version; both are public-domain software of the NIH/NHGRI (Merqury is distributed under the terms in its repository) and DNAKit starts them as separate programs without copying code. Install them with `conda install -c bioconda merqury`, which also brings Meryl, R with ggplot2, Java, bedtools and samtools, and let DNAKit find the scripts through the `MERQURY` environment variable, `merqury.sh` on PATH, or the `merqury_home` argument. A missing script or program raises `BackendUnavailableError` with an install hint before any output is created. Known limits of the upstream scripts: `count_10x.sh` pipes into `meryl ... -`, which Meryl 1.4.2 rejects, and `split_ont_fa.sh` hard-codes a cluster scratch folder, so neither is wrapped; some trailing R plots fail on very small inputs.

## Example data

The examples create their own small data in the preparation step: a random 4,000 bp genome split into two contigs, a mutated second haplotype (3 percent differences), read sets of 100 bp reads, two parental read sets and a child read set, all counted into 15-mer Meryl databases. Every other output goes to a fresh temporary folder.

**Example code:**

```python
import random
import shutil
import subprocess
from pathlib import Path

rng = random.Random(5)
genome = "".join(rng.choice("ACGT") for _ in range(4000))
other = "".join(base if rng.random() > 0.03 else rng.choice("ACGT") for base in genome)


def write_reads(path, sequence, count, error=0.0):
    records = []
    for index in range(count):
        start = rng.randint(0, len(sequence) - 100)
        read = list(sequence[start:start + 100])
        if rng.random() < error:
            read[rng.randint(0, 99)] = "A"
        records.append(f">r{index}\n{''.join(read)}\n")
    Path(path).write_text("".join(records))


Path("asm.fa").write_text(f">c1\n{genome[:2000]}\n>c2\n{genome[2000:]}\n")
Path("hap2.fa").write_text(f">d1\n{other[:2000]}\n>d2\n{other[2000:]}\n")
write_reads("reads.fa", genome, 600, error=0.1)
write_reads("mat.fa", genome, 1500)
write_reads("pat.fa", other, 1500)
write_reads("child.fa", genome, 700)
with open("child.fa", "a") as handle:
    handle.write(Path("child.fa").read_text().replace(">r", ">p")[:0])
write_reads("child2.fa", other, 700)
Path("child.fa").write_text(Path("child.fa").read_text() + Path("child2.fa").read_text().replace(">r", ">q"))
for name in ("asm", "reads", "mat", "pat", "child"):
    shutil.rmtree(f"{name}.meryl", ignore_errors=True)
    subprocess.run(["meryl", "k=15", "count", "output", f"{name}.meryl", f"{name}.fa"], check=True, capture_output=True)
    print(name, "counted")
```

**Example output:**

```text
asm counted
reads counted
mat counted
pat counted
child counted
```

Preparation: FASTA files and 15-mer Meryl databases (`asm.meryl`, `reads.meryl`, `mat.meryl`, `pat.meryl`, `child.meryl`).

## 1) `dnakit.assembly.merqury_assign_reads()`

- **Purpose:** Assign reads to haplotype 1, haplotype 2 or unknown by their hap-mer content (``trio/assign_reads.sh``).
- **Method:** One read file (FASTQ/FASTA, gzipped or not). For each read the output ``<read file>.assigned`` has the read id, k-mer count, hap-1 k-mers and how many were found, the same for hap 2, and the assignment (``hap1_name``, ``hap2_name`` or ``unk``). Result keys: ``assignments`` (file), ``counts`` (reads per assignment).
- **API:** `merqury_assign_reads(read_path: 'str', hap1_name: 'str', hap2_name: 'str', hap1_hapmer_meryl: 'str', hap2_hapmer_meryl: 'str', output_dir: 'str', *, merqury_home: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs Merqury and the programs of its scripts (`conda install -c bioconda merqury`; Meryl, and depending on the script Java, R, samtools, bedtools, Mash, pigz). DNAKit finds the scripts through `merqury_home`, `MERQURY` or `merqury.sh` on PATH and raises `BackendUnavailableError` before creating anything if something is missing.
- **Input:** See the parameter descriptions above and the API signature.
- **Output:** The common result dictionary plus the keys named above.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.assembly import merqury_assign_reads

work = Path(tempfile.mkdtemp())
r = merqury_assign_reads("child.fa", "mat", "pat", "mat.meryl", "pat.meryl", str(work / "assign"), verbose=False)
print(r["status"], sum(r["counts"].values()))
```

- **Example output:**

```text
completed 1400
```

## 2) `dnakit.assembly.merqury_bed_to_bigwig()`

- **Purpose:** Turn a BED file of intervals into a bigWig coverage track (``util/bed_to_bigwig.sh``).
- **Method:** ``fasta_index`` is the ``.fai`` of the assembly the BED refers to. Intervals are converted to a bedGraph with ``bedtools genomecov``, sorted per sequence and converted with ``bedGraphToBigWig``. ``convert_pipes=True`` (Merqury's default) replaces ``|`` by ``_`` in sequence names. The result is ``<bed>.bigwig`` in the new ``output_dir``. Needs bedtools and bedGraphToBigWig.
- **API:** `merqury_bed_to_bigwig(bed_path: 'str', fasta_index: 'str', output_dir: 'str', *, convert_pipes: 'bool' = True, merqury_home: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs Merqury and the programs of its scripts (`conda install -c bioconda merqury`; Meryl, and depending on the script Java, R, samtools, bedtools, Mash, pigz). DNAKit finds the scripts through `merqury_home`, `MERQURY` or `merqury.sh` on PATH and raises `BackendUnavailableError` before creating anything if something is missing.
- **Input:** See the parameter descriptions above and the API signature.
- **Output:** The common result dictionary plus the keys named above.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.assembly import merqury_bed_to_bigwig

work = Path(tempfile.mkdtemp())
__import__('subprocess').run(["samtools", "faidx", "asm.fa"], check=True)
Path("regions.bed").write_text("c1\t10\t200\nc2\t5\t90\n")
r = merqury_bed_to_bigwig("regions.bed", "asm.fa.fai", str(work / "bw"), verbose=False)
print(r["status"], Path(r["bigwig_file"]).name)
```

- **Example output:**

```text
completed regions.bigwig
```

## 3) `dnakit.assembly.merqury_best_k()`

- **Purpose:** Choose the k-mer size for a genome size and a tolerable k-mer collision rate (``best_k.sh``; Fofanov et al. 2004).
- **Method:** ``genome_size`` is the haploid (or diploid, depending on what is evaluated) genome size in bp. ``tolerable_collision_rate`` is the error rate of the read set (Merqury's default 0.001 suits Illumina WGS). k = log4(genome_size * (1 - rate) / rate). Result keys: ``k`` (the real value printed by Merqury) and ``k_rounded_up``.
- **API:** `merqury_best_k(genome_size: 'int', *, tolerable_collision_rate: 'float' = 0.001, merqury_home: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs Merqury and the programs of its scripts (`conda install -c bioconda merqury`; Meryl, and depending on the script Java, R, samtools, bedtools, Mash, pigz). DNAKit finds the scripts through `merqury_home`, `MERQURY` or `merqury.sh` on PATH and raises `BackendUnavailableError` before creating anything if something is missing.
- **Input:** See the parameter descriptions above and the API signature.
- **Output:** The common result dictionary plus the keys named above.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.assembly import merqury_best_k

work = Path(tempfile.mkdtemp())
r = merqury_best_k(3_000_000_000, verbose=False)
print(r["status"], r["k"], r["k_rounded_up"])
```

- **Example output:**

```text
completed 20.7233 21
```

## 4) `dnakit.assembly.merqury_block_n_stats()`

- **Purpose:** Contig and phase-block N* / NG* statistics and the N* curves (``trio/block_n_stats.sh``).
- **Method:** ``assemblies`` (one or two FASTA files) and ``block_beds`` (the matching ``*.phased_block.bed`` files from ``merqury_phase_blocks``, in the same order). ``genome_size`` (bp) adds NG* values. Writes the size tables and N* plots (R) with the prefix. Result key: ``files``.
- **API:** `merqury_block_n_stats(assemblies: 'Sequence[str]', block_beds: 'Sequence[str]', output_dir: 'str', output_prefix: 'str', *, genome_size: 'int | None' = None, merqury_home: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs Merqury and the programs of its scripts (`conda install -c bioconda merqury`; Meryl, and depending on the script Java, R, samtools, bedtools, Mash, pigz). DNAKit finds the scripts through `merqury_home`, `MERQURY` or `merqury.sh` on PATH and raises `BackendUnavailableError` before creating anything if something is missing.
- **Input:** See the parameter descriptions above and the API signature.
- **Output:** The common result dictionary plus the keys named above.
This function has no runnable example here because it was not exercised against the real scripts on the example data.

## 5) `dnakit.assembly.merqury_combine_kmer_dbs()`

- **Purpose:** Combine Meryl k-mer databases (``build/diff.sh``, ``build/intersect.sh`` and ``build/union_sum.sh``).
- **Method:** ``operation`` selects the script. ``"difference"`` (``diff.sh``) takes exactly two databases and keeps the k-mers of the first that are not in the second (``<db1>_not_<db2>.meryl``); with ``filter_solid=True`` (Merqury's default) the result is then filtered with ``merqury_filter_solid_kmers``' rule and linked as ``<output_name>.meryl``, otherwise (``-no-filt``) only linked. ``"intersect"`` (``intersect.sh``) takes two databases and keeps the k-mers of both with the counts of the first (``<db1>_and_<db2>.meryl``, linked as ``<output_name>.meryl``). ``"union_sum"`` (``union_sum.sh``) takes two or more databases and the k-mer size ``k`` and writes the union with summed counts as ``<output_name>.k<k>.meryl``. Result key: ``database``.
- **API:** `merqury_combine_kmer_dbs(operation: 'str', meryl_dbs: 'Sequence[str]', output_dir: 'str', output_name: 'str', *, k: 'int | None' = None, filter_solid: 'bool' = True, merqury_home: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs Merqury and the programs of its scripts (`conda install -c bioconda merqury`; Meryl, and depending on the script Java, R, samtools, bedtools, Mash, pigz). DNAKit finds the scripts through `merqury_home`, `MERQURY` or `merqury.sh` on PATH and raises `BackendUnavailableError` before creating anything if something is missing.
- **Input:** See the parameter descriptions above and the API signature.
- **Output:** The common result dictionary plus the keys named above.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.assembly import merqury_combine_kmer_dbs

work = Path(tempfile.mkdtemp())
r = merqury_combine_kmer_dbs("intersect", ["asm.meryl", "reads.meryl"], str(work / "both"), "shared", verbose=False)
print(r["status"], Path(r["database"]).name)
u = merqury_combine_kmer_dbs("union_sum", ["asm.meryl", "reads.meryl"], str(work / "union"), "all", k=15, verbose=False)
print(u["status"], Path(u["database"]).name)
```

- **Example output:**

```text
completed asm_and_reads.meryl
completed all.k15.meryl
```

## 6) `dnakit.assembly.merqury_exclude_reads()`

- **Purpose:** Remove every read (and its mate) that contains a k-mer of a Meryl database (``trio/exclude_reads.sh``).
- **Method:** ``read_pairs`` is a list of ``[R1]`` or ``[R1, R2]`` FASTQ paths; the filtered reads are written under the new ``output_dir`` with the input file names. ``tenx=True`` ignores the first 23 bases of each R1 (10x barcodes). Result key: ``filtered_reads``.
- **API:** `merqury_exclude_reads(meryl_db: 'str', read_pairs: 'Sequence[Sequence[str]]', output_dir: 'str', *, tenx: 'bool' = False, merqury_home: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs Merqury and the programs of its scripts (`conda install -c bioconda merqury`; Meryl, and depending on the script Java, R, samtools, bedtools, Mash, pigz). DNAKit finds the scripts through `merqury_home`, `MERQURY` or `merqury.sh` on PATH and raises `BackendUnavailableError` before creating anything if something is missing.
- **Input:** See the parameter descriptions above and the API signature.
- **Output:** The common result dictionary plus the keys named above.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.assembly import merqury_exclude_reads

work = Path(tempfile.mkdtemp())
r = merqury_exclude_reads("mat.meryl", [["child.fa"]], str(work / "excl"), verbose=False)
print(r["status"], len(r["filtered_reads"]))
```

- **Example output:**

```text
completed 1
```

## 7) `dnakit.assembly.merqury_false_duplications()`

- **Purpose:** Count assembly k-mers that are extra copies in the single- and two-copy read peaks (``eval/false_duplications.sh``).
- **Method:** ``spectra_cn_hist`` is the ``*.spectra-cn.hist`` table of a pseudo-haplotype or mixed-haplotype assembly (columns Copies, kmer_multiplicity, Count). The cutoff is 1.5 times the read multiplicity at the 1-copy peak. Nothing is written. Result keys: ``cutoff``, ``copies`` (k-mer counts for 1, 2, 3, 4 and >4 copies), ``duplicated_kmers`` (copies above 1), ``all_kmers`` and ``duplicated_percent``.
- **API:** `merqury_false_duplications(spectra_cn_hist: 'str', *, merqury_home: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs Merqury and the programs of its scripts (`conda install -c bioconda merqury`; Meryl, and depending on the script Java, R, samtools, bedtools, Mash, pigz). DNAKit finds the scripts through `merqury_home`, `MERQURY` or `merqury.sh` on PATH and raises `BackendUnavailableError` before creating anything if something is missing.
- **Input:** See the parameter descriptions above and the API signature.
- **Output:** The common result dictionary plus the keys named above.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.assembly import merqury_false_duplications

work = Path(tempfile.mkdtemp())
Path("spectra.hist").write_text("Copies\tkmer_multiplicity\tCount\n1\t1\t50\n1\t2\t10\n1\t5\t400\n2\t5\t30\n")
r = merqury_false_duplications("spectra.hist", verbose=False)
print(r["status"], r["cutoff"], r["copies"]["1"], r["duplicated_percent"])
```

- **Example output:**

```text
completed 8 460 6.12245
```

## 8) `dnakit.assembly.merqury_false_duplications_track()`

- **Purpose:** Locate the extra-copy k-mers of ``merqury_false_duplications`` along the assembly as BED and WIG tracks (``eval/false_duplications_track.sh``).
- **Method:** ``assembly_meryl`` is the assembly's Meryl database (the ``<assembly>.meryl`` that Merqury builds) and ``reads_meryl`` the read database. Outputs in the new ``output_dir``: ``<assembly>.fd.bed`` (merged intervals) and ``<assembly>.fd.wig``. Result keys: ``bed_file``, ``wig_file``, ``interval_count``.
- **API:** `merqury_false_duplications_track(spectra_cn_hist: 'str', assembly_fasta: 'str', assembly_meryl: 'str', reads_meryl: 'str', output_dir: 'str', *, merqury_home: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs Merqury and the programs of its scripts (`conda install -c bioconda merqury`; Meryl, and depending on the script Java, R, samtools, bedtools, Mash, pigz). DNAKit finds the scripts through `merqury_home`, `MERQURY` or `merqury.sh` on PATH and raises `BackendUnavailableError` before creating anything if something is missing.
- **Input:** See the parameter descriptions above and the API signature.
- **Output:** The common result dictionary plus the keys named above.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.assembly import merqury_false_duplications_track

work = Path(tempfile.mkdtemp())
r = merqury_false_duplications_track("spectra.hist", "asm.fa", "asm.meryl", "reads.meryl", str(work / "fd"), verbose=False)
print(r["status"], r["interval_count"])
```

- **Example output:**

```text
completed 0
```

## 9) `dnakit.assembly.merqury_filter_solid_kmers()`

- **Purpose:** Drop the erroneous low-count k-mers of a read Meryl database using the ploidy-derived threshold (``build/filt.sh``).
- **Method:** The k-mer histogram is written (``<db>.hist``), ``kmerHistToPloidyDepth.jar`` derives the ploidy and the count threshold (``<db>.hist.ploidy``) and ``meryl greater-than`` keeps the k-mers above it as ``<db>.gt<threshold>.meryl``. Result keys: ``threshold``, ``solid_database``, ``histogram_file``, ``ploidy_file``.
- **API:** `merqury_filter_solid_kmers(meryl_db: 'str', output_dir: 'str', *, merqury_home: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs Merqury and the programs of its scripts (`conda install -c bioconda merqury`; Meryl, and depending on the script Java, R, samtools, bedtools, Mash, pigz). DNAKit finds the scripts through `merqury_home`, `MERQURY` or `merqury.sh` on PATH and raises `BackendUnavailableError` before creating anything if something is missing.
- **Input:** See the parameter descriptions above and the API signature.
- **Output:** The common result dictionary plus the keys named above.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.assembly import merqury_filter_solid_kmers

work = Path(tempfile.mkdtemp())
r = merqury_filter_solid_kmers("reads.meryl", str(work / "filt"), verbose=False)
print(r["status"], r["threshold"], Path(r["solid_database"]).name)
```

- **Example output:**

```text
completed 2 reads.gt2.meryl
```

## 10) `dnakit.assembly.merqury_hamming_error()`

- **Purpose:** Hamming (phasing) error of assemblies from the hap-mer counts per contig (``trio/hamming_error.sh``).
- **Method:** ``hapmer_count_table`` is the ``*.hapmers.count`` file of ``merqury_hapmer_counts`` and ``assemblies`` the one or two FASTA files (their names, without suffix, match the first column of the table). Writes ``<table>.ham`` (per contig) and ``phase_hammingerror.stats`` in the new ``output_dir``. Result keys: ``rows`` (Assembly, Hap1, Hap2, Hap1Error(%), Hap2Error(%), MajorHapError(%)).
- **API:** `merqury_hamming_error(hapmer_count_table: 'str', assemblies: 'Sequence[str]', output_dir: 'str', *, merqury_home: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs Merqury and the programs of its scripts (`conda install -c bioconda merqury`; Meryl, and depending on the script Java, R, samtools, bedtools, Mash, pigz). DNAKit finds the scripts through `merqury_home`, `MERQURY` or `merqury.sh` on PATH and raises `BackendUnavailableError` before creating anything if something is missing.
- **Input:** See the parameter descriptions above and the API signature.
- **Output:** The common result dictionary plus the keys named above.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.assembly import merqury_hamming_error

work = Path(tempfile.mkdtemp())
from dnakit.assembly import merqury_hapmers

derived = merqury_hapmers("mat.meryl", "pat.meryl", str(work / "hapmers"), child_meryl="child.meryl", verbose=False)
hapmers = {Path(p).name: p for p in derived["hapmer_databases"]}
from dnakit.assembly import merqury_hapmer_counts

hc = merqury_hapmer_counts(
    hapmers["mat.hapmer.meryl"], hapmers["pat.hapmer.meryl"], ["asm.fa", "hap2.fa"], str(work / "counts"), "demo", verbose=False
)
r = merqury_hamming_error(hc["count_file"], ["asm.fa", "hap2.fa"], str(work / "ham"), verbose=False)
print(r["status"], r["columns"])
for row in r["rows"]:
    print(row)
```

- **Example output:**

```text
completed ['Assembly', 'Hap1', 'Hap2', 'Hap1Error(%)', 'Hap2Error(%)', 'MajorHapError(%)']
['asm', '1025', '0', '0', '100', '0']
['hap2', '0', '1009', '100', '0', '0']
```

## 11) `dnakit.assembly.merqury_hapmer_counts()`

- **Purpose:** Count the hap-mers of each haplotype in every assembly contig and draw the blob plot (``trio/hap_blob.sh``).
- **Method:** ``hap1_meryl`` and ``hap2_meryl`` are the hap-mer databases; ``assemblies`` one or two FASTA files. Writes the table ``<prefix>.hapmers.count`` (Assembly, Contig, hap1 count, hap2 count, Size) and the plot ``<prefix>.hapmers.blob.png`` (R). Result keys: ``count_file``, ``columns``, ``row_count``, ``rows``.
- **API:** `merqury_hapmer_counts(hap1_meryl: 'str', hap2_meryl: 'str', assemblies: 'Sequence[str]', output_dir: 'str', output_prefix: 'str', *, merqury_home: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs Merqury and the programs of its scripts (`conda install -c bioconda merqury`; Meryl, and depending on the script Java, R, samtools, bedtools, Mash, pigz). DNAKit finds the scripts through `merqury_home`, `MERQURY` or `merqury.sh` on PATH and raises `BackendUnavailableError` before creating anything if something is missing.
- **Input:** See the parameter descriptions above and the API signature.
- **Output:** The common result dictionary plus the keys named above.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.assembly import merqury_hapmer_counts

work = Path(tempfile.mkdtemp())
from dnakit.assembly import merqury_hapmers

derived = merqury_hapmers("mat.meryl", "pat.meryl", str(work / "hapmers"), child_meryl="child.meryl", verbose=False)
hapmers = {Path(p).name: p for p in derived["hapmer_databases"]}
hc = merqury_hapmer_counts(
    hapmers["mat.hapmer.meryl"], hapmers["pat.hapmer.meryl"], ["asm.fa", "hap2.fa"], str(work / "counts"), "demo", verbose=False
)
print(hc["status"], hc["columns"])
for row in hc["rows"]:
    print(row)
```

- **Example output:**

```text
completed ['Assembly', 'Contig', 'mat.inherited.gt1', 'pat.inherited.gt6', 'Size']
['asm', 'c1', '482', '0', '1986']
['asm', 'c2', '543', '0', '1986']
['hap2', 'd1', '0', '472', '1986']
['hap2', 'd2', '0', '537', '1986']
```

## 12) `dnakit.assembly.merqury_hapmer_spectra()`

- **Purpose:** Hap-mer copy-number spectra and hap-mer completeness of one or two assemblies (``trio/spectra-hap.sh``).
- **Method:** ``reads_meryl`` holds all read k-mers; ``hap1_meryl`` and ``hap2_meryl`` the haplotype-specific k-mers with their read counts (for example ``mat.hapmer.meryl`` and ``pat.hapmer.meryl``); ``assemblies`` is one or two FASTA files. Per assembly and haplotype it writes ``<prefix>.<asm>.<hap>.spectra-cn.hist`` and plots (R), and ``<prefix>.completeness.stats`` (assembly, haplotype, hap-mers found, hap-mers total, percent). ``homopolymer_compressed`` is the script's ``-c`` flag. Result keys: ``completeness`` (rows) and ``histograms``.
- **API:** `merqury_hapmer_spectra(reads_meryl: 'str', hap1_meryl: 'str', hap2_meryl: 'str', assemblies: 'Sequence[str]', output_dir: 'str', output_prefix: 'str', *, homopolymer_compressed: 'bool' = False, merqury_home: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs Merqury and the programs of its scripts (`conda install -c bioconda merqury`; Meryl, and depending on the script Java, R, samtools, bedtools, Mash, pigz). DNAKit finds the scripts through `merqury_home`, `MERQURY` or `merqury.sh` on PATH and raises `BackendUnavailableError` before creating anything if something is missing.
- **Input:** See the parameter descriptions above and the API signature.
- **Output:** The common result dictionary plus the keys named above.
This function has no runnable example here because it was not exercised against the real scripts on the example data.

## 13) `dnakit.assembly.merqury_hapmers()`

- **Purpose:** Derive haplotype-specific k-mers (hap-mers) from parental k-mer databases, and the inherited ones from a child (``trio/hapmers.sh``).
- **Method:** ``hap1_meryl`` and ``hap2_meryl`` hold the k-mers of both parents (all built with the same k). The script writes ``<hap1>.only.meryl``, ``<hap2>.only.meryl`` and ``shrd.meryl``; with ``child_meryl`` it adds the inherited hap-mers ``<hap>.hapmer.meryl``, ``shrd.inherited.meryl``, ``cutoffs.txt`` and the ``inherited_hapmers`` spectra plot. ``no_filt=True`` skips the solid k-mer filtering of the parental databases (use only when they do not come from regular sequencing). Needs Java and, with a child, R. Result keys: ``hapmer_databases`` (the .meryl directories that exist), ``cutoffs`` (parsed ``cutoffs.txt``).
- **API:** `merqury_hapmers(hap1_meryl: 'str', hap2_meryl: 'str', output_dir: 'str', *, child_meryl: 'str | None' = None, no_filt: 'bool' = False, merqury_home: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs Merqury and the programs of its scripts (`conda install -c bioconda merqury`; Meryl, and depending on the script Java, R, samtools, bedtools, Mash, pigz). DNAKit finds the scripts through `merqury_home`, `MERQURY` or `merqury.sh` on PATH and raises `BackendUnavailableError` before creating anything if something is missing.
- **Input:** See the parameter descriptions above and the API signature.
- **Output:** The common result dictionary plus the keys named above.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.assembly import merqury_hapmers

work = Path(tempfile.mkdtemp())
hapmers_result = merqury_hapmers("mat.meryl", "pat.meryl", str(work / "hapmers"), child_meryl="child.meryl", verbose=False)
print(hapmers_result["status"], hapmers_result["cutoffs"])
hapmers = {Path(p).name: p for p in hapmers_result["hapmer_databases"]}
print(sorted(name for name in hapmers if name.endswith(".hapmer.meryl")))
```

- **Example output:**

```text
completed {'mat': 1, 'pat': 6, 'shared': 2}
['mat.hapmer.meryl', 'pat.hapmer.meryl']
```

## 14) `dnakit.assembly.merqury_kmer_histogram_table()`

- **Purpose:** Export the k-mer multiplicity histograms of Meryl databases as the table that Merqury's plot script reads (``plot/to_hist_for_plotting.sh``).
- **Method:** ``meryl_dbs`` maps a label to a Meryl database. The TSV ``output_path`` has the columns ``kmer`` (the label), ``kmer_multiplicity`` and ``Count``. Result keys: ``columns``, ``row_count``, ``rows``.
- **API:** `merqury_kmer_histogram_table(meryl_dbs: 'Mapping[str, str]', output_path: 'str', *, merqury_home: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs Merqury and the programs of its scripts (`conda install -c bioconda merqury`; Meryl, and depending on the script Java, R, samtools, bedtools, Mash, pigz). DNAKit finds the scripts through `merqury_home`, `MERQURY` or `merqury.sh` on PATH and raises `BackendUnavailableError` before creating anything if something is missing.
- **Input:** See the parameter descriptions above and the API signature.
- **Output:** The common result dictionary plus the keys named above.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.assembly import merqury_kmer_histogram_table

work = Path(tempfile.mkdtemp())
r = merqury_kmer_histogram_table({"reads": "reads.meryl"}, str(work / "hist.tsv"), verbose=False)
print(r["status"], r["columns"], r["rows"][:2])
```

- **Example output:**

```text
completed ['kmer', 'kmer_multiplicity', 'Count'] [['reads', '1', '487'], ['reads', '2', '4']]
```

## 15) `dnakit.assembly.merqury_kmer_multiplicity_track()`

- **Purpose:** Write a k-mer count track over the assembly sequences as WIG and bigWig (``eval/asm_multiplicity.sh`` and ``eval/read_multiplicity.sh``).
- **Method:** ``source="assembly"`` (``asm_multiplicity.sh``) counts how often each assembly k-mer occurs in the assembly itself (copy number; ``meryl_db`` is the Meryl database of the assembly); ``source="reads"`` (``read_multiplicity.sh``) reports the read multiplicity of each assembly k-mer (``meryl_db`` is the read database). The files are ``<prefix>.copies.wig`` / ``<prefix>.asm_multiplicity.bigWig`` or ``<prefix>.read.wig`` / ``<prefix>.read_multiplicity.bigwig``. Needs samtools and wigToBigWig.
- **API:** `merqury_kmer_multiplicity_track(assembly_fasta: 'str', meryl_db: 'str', output_dir: 'str', *, source: 'str' = 'assembly', prefix: 'str' = 'track', merqury_home: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs Merqury and the programs of its scripts (`conda install -c bioconda merqury`; Meryl, and depending on the script Java, R, samtools, bedtools, Mash, pigz). DNAKit finds the scripts through `merqury_home`, `MERQURY` or `merqury.sh` on PATH and raises `BackendUnavailableError` before creating anything if something is missing.
- **Input:** See the parameter descriptions above and the API signature.
- **Output:** The common result dictionary plus the keys named above.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.assembly import merqury_kmer_multiplicity_track

work = Path(tempfile.mkdtemp())
for source, db in (("assembly", "asm.meryl"), ("reads", "reads.meryl")):
    r = merqury_kmer_multiplicity_track("asm.fa", db, str(work / source), source=source, verbose=False)
    print(source, r["status"], Path(r["wig_file"]).name, Path(r["bigwig_file"]).name)
```

- **Example output:**

```text
assembly completed track.copies.wig track.asm_multiplicity.bigWig
reads completed track.read.wig track.read_multiplicity.bigwig
```

## 16) `dnakit.assembly.merqury_mash_qv()`

- **Purpose:** Quick assembly QV estimate from a Mash sketch of the assembly screened against the reads (``merqury-mash.sh``).
- **Method:** ``k`` is the k-mer size and ``read_paths`` the read files (FASTA/FASTQ, possibly gzipped); they are listed in a file of file names as the script expects. Mash screen identities are converted to QV; the values do not match the full Merqury run exactly (the script says so). Output files ``<assembly>.msh.idy`` and ``<assembly>.msh.qv`` are written next to the linked assembly. Result keys: ``qv_file``, ``columns``, ``rows``.
- **API:** `merqury_mash_qv(assembly_fasta: 'str', read_paths: 'Sequence[str]', output_dir: 'str', *, k: 'int' = 21, cpus: 'int' = 1, merqury_home: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs Merqury and the programs of its scripts (`conda install -c bioconda merqury`; Meryl, and depending on the script Java, R, samtools, bedtools, Mash, pigz). DNAKit finds the scripts through `merqury_home`, `MERQURY` or `merqury.sh` on PATH and raises `BackendUnavailableError` before creating anything if something is missing.
- **Input:** See the parameter descriptions above and the API signature.
- **Output:** The common result dictionary plus the keys named above.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.assembly import merqury_mash_qv

work = Path(tempfile.mkdtemp())
r = merqury_mash_qv("asm.fa", ["reads.fa"], str(work / "mash"), k=15, verbose=False)
print(r["status"], r["columns"], r["rows"][0][:2])
```

- **Example output:**

```text
completed ['assembly', 'reads', 'qv', 'error_rate'] ['asm', '3965']
```

## 17) `dnakit.assembly.merqury_per_sequence_qv()`

- **Purpose:** Estimate the consensus quality value (QV) of every sequence of an assembly from a read k-mer database (``eval/per_seq_qv.sh``).
- **Method:** ``read_meryl`` is the Meryl database of the reads. The assembly is counted into k-mers of the same size, the k-mers found only in the assembly are looked up per sequence, and QV = -10 log10(error rate) with the error rate 1 - (1 - e/n)^(1/k). Output ``<prefix>.qv`` in the new ``output_dir`` has the columns sequence, assembly-only k-mers, k-mers, QV and error rate. Result keys: ``qv_file``, ``columns``, ``row_count``, ``rows``.
- **API:** `merqury_per_sequence_qv(assembly_fasta: 'str', read_meryl: 'str', output_dir: 'str', *, prefix: 'str' = 'per_seq', merqury_home: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs Merqury and the programs of its scripts (`conda install -c bioconda merqury`; Meryl, and depending on the script Java, R, samtools, bedtools, Mash, pigz). DNAKit finds the scripts through `merqury_home`, `MERQURY` or `merqury.sh` on PATH and raises `BackendUnavailableError` before creating anything if something is missing.
- **Input:** See the parameter descriptions above and the API signature.
- **Output:** The common result dictionary plus the keys named above.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.assembly import merqury_per_sequence_qv

work = Path(tempfile.mkdtemp())
r = merqury_per_sequence_qv("asm.fa", "reads.meryl", str(work / "qv"), verbose=False)
print(r["status"], r["columns"])
for row in r["rows"]:
    print(row)
```

- **Example output:**

```text
completed ['sequence', 'error_kmers', 'total_kmers', 'qv', 'error_rate']
['c1', '0', '1986', '+inf', '0']
['c2', '7', '1986', '36.2826', '0.000235366']
```

## 18) `dnakit.assembly.merqury_phase_blocks()`

- **Purpose:** Build phase blocks from hap-mer markers along an assembly and count switch errors (``trio/phase_block.sh`` and ``trio/switch_error.sh``).
- **Method:** ``phase_block.sh`` finds the gaps, looks up the hap-mers of both haplotypes along the assembly (BED and WIG tracks) and builds blocks with its fixed default (100 switches within 20,000 bp). When ``num_switch`` and ``short_range`` differ from that, ``switch_error.sh`` is run again on ``<prefix>.sort.bed`` with your values (``include_gaps=False`` restricts blocks to within contigs). Needs Java. Result keys: ``block_files`` (``*.phased_block.bed`` per setting), ``switch_error`` (rows of ``*.switches.txt``), ``stats_files``.
- **API:** `merqury_phase_blocks(assembly_fasta: 'str', hap1_meryl: 'str', hap2_meryl: 'str', output_dir: 'str', output_prefix: 'str', *, num_switch: 'int' = 100, short_range: 'int' = 20000, include_gaps: 'bool' = True, merqury_home: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs Merqury and the programs of its scripts (`conda install -c bioconda merqury`; Meryl, and depending on the script Java, R, samtools, bedtools, Mash, pigz). DNAKit finds the scripts through `merqury_home`, `MERQURY` or `merqury.sh` on PATH and raises `BackendUnavailableError` before creating anything if something is missing.
- **Input:** See the parameter descriptions above and the API signature.
- **Output:** The common result dictionary plus the keys named above.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.assembly import merqury_phase_blocks

work = Path(tempfile.mkdtemp())
from dnakit.assembly import merqury_hapmers

derived = merqury_hapmers("mat.meryl", "pat.meryl", str(work / "hapmers"), child_meryl="child.meryl", verbose=False)
hapmers = {Path(p).name: p for p in derived["hapmer_databases"]}
r = merqury_phase_blocks(
    "asm.fa", hapmers["mat.hapmer.meryl"], hapmers["pat.hapmer.meryl"], str(work / "phase"), "demo",
    num_switch=50, short_range=10000, verbose=False,
)
print(r["status"], [Path(p).name for p in r["block_files"]])
```

- **Example output:**

```text
completed ['demo.100_20000.phased_block.bed', 'demo.50_10000.phased_block.bed']
```

## 19) `dnakit.assembly.merqury_plot_blob()`

- **Purpose:** Plot hap-mer (or phased-block) counts per contig as a blob plot (``plot/plot_blob.R``).
- **Method:** ``count_table`` is a tab-separated table with a header and the columns category, sequence id, hap1 count, hap2 count and sequence size (the ``*.hapmers.count`` of ``merqury_hapmer_counts``). Writes ``<prefix>.png`` (or ``.pdf``). Needs Rscript with ggplot2, argparse and scales. Result key: ``plot``.
- **API:** `merqury_plot_blob(count_table: 'str', output_dir: 'str', output_prefix: 'str', *, width: 'float' = 6.5, height: 'float' = 6.0, pdf: 'bool' = False, merqury_home: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs Merqury and the programs of its scripts (`conda install -c bioconda merqury`; Meryl, and depending on the script Java, R, samtools, bedtools, Mash, pigz). DNAKit finds the scripts through `merqury_home`, `MERQURY` or `merqury.sh` on PATH and raises `BackendUnavailableError` before creating anything if something is missing.
- **Input:** See the parameter descriptions above and the API signature.
- **Output:** The common result dictionary plus the keys named above.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.assembly import merqury_plot_blob

work = Path(tempfile.mkdtemp())
Path("c.count").write_text("category\tseq\thap1\thap2\tsize\nhap1\tc1\t10\t3\t2000\nhap1\tc2\t4\t5\t1000\n")
r = merqury_plot_blob("c.count", str(work / "blob"), "blob", verbose=False)
print(r["status"], Path(r["plot"]).name)
```

- **Example output:**

```text
completed blob.png
```

## 20) `dnakit.assembly.merqury_plot_block_n()`

- **Purpose:** Plot N* (or NG* with ``genome_size``) continuity curves of phase blocks, scaffolds and contigs (``plot/plot_block_N.R``).
- **Method:** Give at least one of ``block_sizes``, ``scaffold_sizes`` and ``contig_sizes``: ``.sizes`` tables with three tab-separated columns (type, group, size; longest first). Writes ``<prefix>.block.N.png`` and ``<prefix>.continuity.N.png`` (``NG`` instead of ``N`` with ``genome_size``; ``.pdf`` with ``pdf=True``) in the new ``output_dir``. Needs Rscript with ggplot2, argparse and scales. Result key: ``plots``.
- **API:** `merqury_plot_block_n(output_dir: 'str', output_prefix: 'str', *, block_sizes: 'str | None' = None, scaffold_sizes: 'str | None' = None, contig_sizes: 'str | None' = None, genome_size: 'int' = 0, width: 'float' = 6.0, height: 'float' = 5.0, pdf: 'bool' = False, merqury_home: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs Merqury and the programs of its scripts (`conda install -c bioconda merqury`; Meryl, and depending on the script Java, R, samtools, bedtools, Mash, pigz). DNAKit finds the scripts through `merqury_home`, `MERQURY` or `merqury.sh` on PATH and raises `BackendUnavailableError` before creating anything if something is missing.
- **Input:** See the parameter descriptions above and the API signature.
- **Output:** The common result dictionary plus the keys named above.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.assembly import merqury_plot_block_n

work = Path(tempfile.mkdtemp())
Path("c.sizes").write_text("contig\tasm\t2000\ncontig\tasm\t1000\n")
r = merqury_plot_block_n(str(work / "bn"), "bn", contig_sizes="c.sizes", verbose=False)
print(r["status"], [Path(p).name for p in r["plots"]])
```

- **Example output:**

```text
completed ['bn.block.N.png', 'bn.continuity.N.png']
```

## 21) `dnakit.assembly.merqury_plot_spectra_cn()`

- **Purpose:** Plot k-mer copy-number spectra from a ``*.spectra-cn.hist`` table (``plot/plot_spectra_cn.R``).
- **Method:** ``histogram`` has the columns ``Assembly``/``Copies``, ``kmer_multiplicity`` and ``Count`` (see ``merqury_kmer_histogram_table``). ``only_histogram`` is the optional ``*.only.hist`` (assembly-only k-mers) and ``cutoff_file`` the optional solid k-mer cutoffs. ``plot_type`` is ``line``, ``fill``, ``stack`` or ``all`` (one file each: ``<prefix>.ln|fl|st.png``); ``xmax``/``ymax`` of 0 let the script choose. ``pdf=True`` writes PDF. Needs Rscript with ggplot2, argparse and scales. Result key: ``plots``.
- **API:** `merqury_plot_spectra_cn(histogram: 'str', output_dir: 'str', output_prefix: 'str', *, only_histogram: 'str | None' = None, cutoff_file: 'str | None' = None, plot_type: 'str' = 'all', width: 'float' = 6.0, height: 'float' = 5.0, xmax: 'int' = 0, ymax: 'int' = 0, pdf: 'bool' = False, merqury_home: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs Merqury and the programs of its scripts (`conda install -c bioconda merqury`; Meryl, and depending on the script Java, R, samtools, bedtools, Mash, pigz). DNAKit finds the scripts through `merqury_home`, `MERQURY` or `merqury.sh` on PATH and raises `BackendUnavailableError` before creating anything if something is missing.
- **Input:** See the parameter descriptions above and the API signature.
- **Output:** The common result dictionary plus the keys named above.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.assembly import merqury_plot_spectra_cn

work = Path(tempfile.mkdtemp())
from dnakit.assembly import merqury_kmer_histogram_table

merqury_kmer_histogram_table({"reads": "reads.meryl"}, str(work / "h.tsv"), verbose=False)
r = merqury_plot_spectra_cn(str(work / "h.tsv"), str(work / "plots"), "spectra", plot_type="line", verbose=False)
print(r["status"], [Path(p).name for p in r["plots"]])
```

- **Example output:**

```text
completed ['spectra.ln.png']
```

## 22) `dnakit.assembly.merqury_split_10x_fastq()`

- **Purpose:** Trim the first 23 bases of 10x R1 reads and split each file into pieces of 300 million lines, gzip-compressed (``build/split_10x.sh``).
- **Method:** ``read_paths`` are FASTQ R1 files (gzipped or not). The pieces go to ``<output_dir>/split/<name>.NNNN.fq.gz``; ``<fofn>.<n>`` lists the pieces of the n-th input. Needs pigz. Result keys: ``pieces`` (per input file).
- **API:** `merqury_split_10x_fastq(read_paths: 'Sequence[str]', output_dir: 'str', *, threads: 'int' = 8, merqury_home: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs Merqury and the programs of its scripts (`conda install -c bioconda merqury`; Meryl, and depending on the script Java, R, samtools, bedtools, Mash, pigz). DNAKit finds the scripts through `merqury_home`, `MERQURY` or `merqury.sh` on PATH and raises `BackendUnavailableError` before creating anything if something is missing.
- **Input:** See the parameter descriptions above and the API signature.
- **Output:** The common result dictionary plus the keys named above.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.assembly import merqury_split_10x_fastq

work = Path(tempfile.mkdtemp())
r = merqury_split_10x_fastq(["reads.fa"], str(work / "split"), verbose=False)
print(r["status"], [v for v in r["pieces"].values()])
```

- **Example output:**

```text
completed [['split/reads.fa.0000.fq.gz']]
```
