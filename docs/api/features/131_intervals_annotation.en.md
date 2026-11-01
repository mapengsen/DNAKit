# Genomic intervals and annotation: pybedtools, bioframe, pyranges, genomicranges, gtars and pyGenomeViz

This group of functions is a set of adapters over optional packages and external programs: [pybedtools](https://github.com/daler/pybedtools) (with bedtools, samtools and the UCSC bigBed/bigWig tools), [bioframe](https://github.com/open2c/bioframe), [pyranges1](https://github.com/pyranges/pyranges1) (with [pysam](https://github.com/pysam-developers/pysam)), [genomicranges](https://github.com/BiocPy/genomicranges), [gtars](https://github.com/databio/gtars) and [pyGenomeViz](https://github.com/moshi4/pyGenomeViz) (with MMseqs2, MUMmer and progressiveMauve). They get chromosome sizes and centromeres, convert between BED, bedGraph, wiggle, bigBed and bigWig, count mapped reads, build and draw interval sets (Venn sets, binary heat maps, interval tracks, itemRgb colours), read BAM alignments as intervals, read UCSC gene-model GTF files, map an HGVS variant to a GA4GH VRS identifier, and draw genome feature maps and genome comparisons. The 26 interval and annotation functions are imported with `from dnakit.annotation import ...`; the six pyGenomeViz functions are imported with `from dnakit.visualization import ...`.

The installed metadata of pybedtools, bioframe, pyranges1, genomicranges, pyGenomeViz and pysam states the MIT licence; the installed gtars package states no licence in its metadata. DNAKit calls these packages and copies none of their code. The UCSC tools, bedtools, samtools, MMseqs2, MUMmer and progressiveMauve are separate programs with their own licences that you install yourself; DNAKit only starts them with fixed argument lists and never through a shell.

The dependencies are optional and are imported only when the matching function is called. If a package is missing, the function raises `BackendUnavailableError` with the install command; if a program is missing it raises `BackendUnavailableError` with an install hint. Output files must be new: the functions refuse to overwrite an existing file. Several examples below show that missing-program error, because the UCSC tools, samtools, pypairix, MMseqs2, MUMmer and progressiveMauve are not installed in the environment where the outputs were produced; for these, a second block labelled "Usage with the dependency installed" shows a realistic call that was not executed.

## Example data

The examples use these tiny files. Each example writes its outputs to a temporary directory.

**Example file:** `peaks_a.bed`

```text
chr1	100	200	a1	10	+
chr1	300	400	a2	50	-
chr2	50	150	a3	90	+
```

**Example file:** `peaks_b.bed`

```text
chr1	150	250	b1	20	+
chr1	700	800	b2	60	+
chr2	100	200	b3	70	-
```

**Example file:** `chrom.sizes`

```text
chr1	1000
chr2	800
```

**Example file:** `signal.bedgraph`

```text
chr1	0	100	1.5
chr1	100	200	3
```

**Example file:** `signal.wig`

```text
variableStep chrom=chr1 span=10
101	2.0
111	3.0
```

**Example file:** `reads.sam`

```text
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:1000
r1	0	chr1	101	60	10M	*	0	0	ACGTACGTAC	IIIIIIIIII
r2	16	chr1	201	30	8M	*	0	0	ACGTACGT	IIIIIIII
r3	1024	chr1	301	60	8M	*	0	0	ACGTACGT	IIIIIIII
r4	4	*	0	0	*	*	0	0	ACGTACGT	IIIIIIII
```

**Example file:** `pairs.txt`

```text
r1	chr1	100	chr1	500	+	-
r2	chr1	120	chr2	300	+	+
```

**Example file:** `genes.gtf`

```text
chr1	refGene	exon	101	200	.	+	.	gene_id "GENE1"; transcript_id "TX1";
chr1	refGene	exon	301	400	.	+	.	gene_id "GENE1"; transcript_id "TX1";
chr2	refGene	exon	51	150	.	-	.	gene_id "GENE2"; transcript_id "TX2";
```

**Example file:** `tiny.fa`

```text
>chr1
ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
```

**Example file:** `genomeA.gbk`

```text
LOCUS       genomeA                   84 bp    DNA              UNK 01-JAN-1980
DEFINITION  demo genome.
ACCESSION   genomeA
VERSION     genomeA
KEYWORDS    .
SOURCE      .
  ORGANISM  .
            .
FEATURES             Location/Qualifiers
     CDS             11..37
                     /product="protein 10"
                     /protein_id="genomeA_10"
                     /translation="MAKGLTEV"
     CDS             complement(46..72)
                     /product="protein 45"
                     /protein_id="genomeA_45"
                     /translation="MKRLGTLD"
ORIGIN
        1 actgactgac atggctaaag gtctgaccga agtttaaact gactgttaat ccagggtacc
       61 cagacgtttc atactgactg actg
//
```

**Example file:** `genomeB.gbk`

```text
LOCUS       genomeB                   86 bp    DNA              UNK 01-JAN-1980
DEFINITION  demo genome.
ACCESSION   genomeB
VERSION     genomeB
KEYWORDS    .
SOURCE      .
  ORGANISM  .
            .
FEATURES             Location/Qualifiers
     CDS             15..41
                     /product="protein 14"
                     /protein_id="genomeB_14"
                     /translation="MAKGLTEV"
     CDS             complement(48..74)
                     /product="protein 47"
                     /protein_id="genomeB_47"
                     /translation="MKRLGTLD"
ORIGIN
        1 tgcatgcatg catgatggct aaaggtctga ccgaagttta atgcatgtta atccagggta
       61 cccagacgtt tcattgcatg catgca
//
```

## 1) `dnakit.annotation.add_color()`

- **Purpose:** Write a BED file whose itemRgb colour follows each feature's score through a matplotlib colormap.
- **Method:** Read the scores from the fifth column, scale them linearly between the minimum and maximum (or `score_min` and `score_max`) with `matplotlib.colors.Normalize`, and apply `pybedtools.featurefuncs.add_color` to every feature with the chosen colormap. The result is saved with `saveas()`. The returned colours are read back from the ninth column of the new file.
- **API:** `add_color(bed: 'str', output_bed: 'str', *, colormap: 'str' = 'viridis', score_min: 'float | None' = None, score_max: 'float | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages pybedtools (extra `dnakit[bio-pybedtools]`) and matplotlib (extra `dnakit[bio-plots]`).
- **Input:** `bed`, an existing BED file in which every feature has a numeric score in column 5; `output_bed`, a new file in an existing folder; `colormap`, the name of a matplotlib colormap (an unknown name raises `ConfigurationError`); `score_min` and `score_max`, optional finite bounds with `score_min <= score_max` (the file's minimum and maximum by default). At most 100,000 features.
- **Output:** A dictionary with `status`, `backend`, `output_path` and `data`, where `data` has `feature_count`, `colormap`, `score_min`, `score_max` and `rgb`, the list of itemRgb strings in file order.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.annotation import add_color

out = Path(tempfile.mkdtemp()) / "colored.bed"
result = add_color("peaks_a.bed", str(out), colormap="viridis", verbose=False)
data = result["data"]
print(data["feature_count"], data["score_min"], data["score_max"])
print(data["rgb"])
```

- **Example output:**

```text
3 10.0 90.0
['68,1,84', '32,144,140', '253,231,36']
```

## 2) `dnakit.annotation.bam_to_bigwig()`

- **Purpose:** Write a bigWig of read coverage from a BAM file, optionally scaled to reads per million mapped reads.
- **Method:** Run `BedTool(bam).genome_coverage(bg=True, split=True)` (bedtools genomecov, as pybedtools does) into a temporary bedGraph, then convert it with UCSC `bedGraphToBigWig`. With `scale=True` the factor `1e6 / mapped reads` is passed to genomecov, where the read count comes from `mapped_read_count()` (samtools). If bedGraphToBigWig complains about order, the bedGraph is sorted with bedtools and the conversion is retried. The temporary bedGraph files are removed afterwards. The `chromsizes` file goes to bedGraphToBigWig only; bedtools ignores `-g` for BAM input.
- **API:** `bam_to_bigwig(bam: 'str', chromsizes: 'str', output_bigwig: 'str', *, scale: 'bool' = False, samtools_path: 'str | None' = None, bedgraphtobigwig_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package pybedtools (extra `dnakit[bio-pybedtools]`), the external programs bedtools and UCSC bedGraphToBigWig on PATH (or `bedgraphtobigwig_path`), and samtools on PATH when `scale=True`.
- **Input:** `bam`, an existing coordinate-sorted BAM; `chromsizes`, a two-column chromosome sizes file; `output_bigwig`, a new path; `scale`, a boolean. A BAM with no mapped reads and `scale=True` raises `ConfigurationError`.
- **Output:** A dictionary with `status` (`completed` or `failed`), `backend`, `executable`, `command`, `return_code`, `outputs`, `log_tail`, `mapped_reads` (None without scaling), `scale_factor` and, if the retry was used, `sorted_before_conversion`.
- **Example code:**

```python
import tempfile
from pathlib import Path

import pysam

from dnakit.annotation import bam_to_bigwig
from dnakit.exceptions import BackendUnavailableError

work = Path(tempfile.mkdtemp())
pysam.view("-b", "-o", str(work / "reads.bam"), "reads.sam", catch_stdout=False)
try:
    bam_to_bigwig(
        str(work / "reads.bam"),
        "chrom.sizes",
        str(work / "coverage.bw"),
        bedgraphtobigwig_path="/nonexistent/bedGraphToBigWig",
        verbose=False,
    )
except BackendUnavailableError as error:
    print(type(error).__name__, "bedGraphToBigWig" in str(error))
```

- **Example output:**

```text
BackendUnavailableError True
```

The BAM is made from `reads.sam` with pysam. In this environment bedGraphToBigWig is not installed, so the example points `bedgraphtobigwig_path` at a program that does not exist and shows the error. (The check for bedtools comes first, so on a machine without bedtools the error names bedtools instead.)

Usage with the dependency installed:

```python
from dnakit.annotation import bam_to_bigwig

result = bam_to_bigwig("reads.bam", "chrom.sizes", "coverage.bw", scale=True)
print(result["status"], result["mapped_reads"], result["outputs"])
```

## 3) `dnakit.annotation.bedgraph_to_bigwig()`

- **Purpose:** Convert a sorted bedGraph file to an indexed bigWig file.
- **Method:** Run `bedGraphToBigWig <bedgraph> <chromsizes> <output>` through DNAKit's bounded command runner, without a shell. Success means a zero exit code and an output file.
- **API:** `bedgraph_to_bigwig(bedgraph: 'str', chromsizes: 'str', output_bigwig: 'str', *, bedgraphtobigwig_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program UCSC bedGraphToBigWig on PATH (or `bedgraphtobigwig_path`). It needs no Python package.
- **Input:** `bedgraph`, an existing bedGraph sorted by chromosome and start; `chromsizes`, a chromosome sizes file; `output_bigwig`, a new path in an existing folder; optional `bedgraphtobigwig_path`.
- **Output:** A dictionary with `status` (`completed` or `failed`), `backend`, `executable`, `command`, `return_code`, `outputs` (the bigWig path when it was produced) and `log_tail` (the last 4000 characters of the program output).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.annotation import bedgraph_to_bigwig
from dnakit.exceptions import BackendUnavailableError

out = Path(tempfile.mkdtemp()) / "signal.bw"
try:
    bedgraph_to_bigwig(
        "signal.bedgraph",
        "chrom.sizes",
        str(out),
        bedgraphtobigwig_path="/nonexistent/bedGraphToBigWig",
        verbose=False,
    )
except BackendUnavailableError as error:
    print(type(error).__name__, "bedGraphToBigWig path does not exist" in str(error))
```

- **Example output:**

```text
BackendUnavailableError True
```

The program is not installed here, so the example shows the error for a program path that does not exist.

Usage with the dependency installed:

```python
from dnakit.annotation import bedgraph_to_bigwig

result = bedgraph_to_bigwig("signal.bedgraph", "chrom.sizes", "signal.bw")
print(result["status"], result["outputs"])
```

## 4) `dnakit.annotation.bigbed()`

- **Purpose:** Convert a BED file to an indexed bigBed file.
- **Method:** Run UCSC `bedToBigBed` with `-blockSize`, `-itemsPerSlot` and `-type` (and `-unc`, `-tab` and `-as=` when requested). When `bed_type` is not given, the type is `bed<N>` where N is the number of tab-separated columns of the first record.
- **API:** `bigbed(bed: 'str', chromsizes: 'str', output_bigbed: 'str', *, bed_type: 'str | None' = None, block_size: 'int' = 256, items_per_slot: 'int' = 512, autosql: 'str | None' = None, unc: 'bool' = False, tab: 'bool' = False, bedtobigbed_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program UCSC bedToBigBed on PATH (or `bedtobigbed_path`). It needs no Python package.
- **Input:** `bed`, an existing BED file sorted by chromosome and start; `chromsizes`, a chromosome sizes file; `output_bigbed`, a new path; `bed_type` such as `bed3`, `bed6` or `bed6+4`; `block_size` and `items_per_slot`, positive integers (at most 1,000,000); `autosql`, an optional existing `.as` file; `unc` and `tab`, booleans.
- **Output:** The same dictionary as `bedgraph_to_bigwig()`: `status`, `backend`, `executable`, `command`, `return_code`, `outputs` and `log_tail`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.annotation import bigbed
from dnakit.exceptions import BackendUnavailableError

out = Path(tempfile.mkdtemp()) / "peaks.bb"
try:
    bigbed(
        "peaks_a.bed",
        "chrom.sizes",
        str(out),
        bed_type="bed6",
        bedtobigbed_path="/nonexistent/bedToBigBed",
        verbose=False,
    )
except BackendUnavailableError as error:
    print(type(error).__name__, "bedToBigBed path does not exist" in str(error))
```

- **Example output:**

```text
BackendUnavailableError True
```

The program is not installed here, so the example shows the error for a program path that does not exist.

Usage with the dependency installed:

```python
from dnakit.annotation import bigbed

result = bigbed("peaks_a.bed", "chrom.sizes", "peaks.bb", bed_type="bed6")
print(result["status"], result["command"])
```

## 5) `dnakit.annotation.bigbed_to_bed()`

- **Purpose:** Extract the intervals of an indexed bigBed file to a BED file.
- **Method:** Run UCSC `bigBedToBed` with optional `-chrom`, `-start`, `-end` and `-maxItems`, followed by the input and output paths.
- **API:** `bigbed_to_bed(bigbed: 'str', output_bed: 'str', *, chrom: 'str | None' = None, start: 'int | None' = None, end: 'int | None' = None, max_items: 'int | None' = None, bigbedtobed_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program UCSC bigBedToBed on PATH (or `bigbedtobed_path`). It needs no Python package.
- **Input:** `bigbed`, an existing bigBed file; `output_bed`, a new path; `chrom`, `start` and `end` to restrict the region (`start` or `end` without `chrom` raises `ConfigurationError`); `max_items`, a positive integer.
- **Output:** A dictionary with `status`, `backend`, `executable`, `command`, `return_code`, `outputs` and `log_tail`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.annotation import bigbed_to_bed
from dnakit.exceptions import BackendUnavailableError

work = Path(tempfile.mkdtemp())
(work / "peaks.bb").write_bytes(b"")
try:
    bigbed_to_bed(
        str(work / "peaks.bb"),
        str(work / "peaks.bed"),
        chrom="chr1",
        bigbedtobed_path="/nonexistent/bigBedToBed",
        verbose=False,
    )
except BackendUnavailableError as error:
    print(type(error).__name__, "bigBedToBed path does not exist" in str(error))
```

- **Example output:**

```text
BackendUnavailableError True
```

The program is not installed here, so the example shows the error for a program path that does not exist. An empty placeholder file stands in for the bigBed, because the program is looked up before the file is opened.

Usage with the dependency installed:

```python
from dnakit.annotation import bigbed_to_bed

result = bigbed_to_bed("peaks.bb", "peaks_chr1.bed", chrom="chr1", start=0, end=500)
print(result["status"], result["outputs"])
```

## 6) `dnakit.annotation.bigwig_to_bedgraph()`

- **Purpose:** Extract the values of an indexed bigWig file to a bedGraph file.
- **Method:** Run UCSC `bigWigToBedGraph` with optional `-chrom`, `-start`, `-end` and `-udcDir`, followed by the input and output paths.
- **API:** `bigwig_to_bedgraph(bigwig: 'str', output_bedgraph: 'str', *, chrom: 'str | None' = None, start: 'int | None' = None, end: 'int | None' = None, udc_dir: 'str | None' = None, bigwigtobedgraph_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program UCSC bigWigToBedGraph on PATH (or `bigwigtobedgraph_path`). It needs no Python package.
- **Input:** `bigwig`, an existing bigWig file; `output_bedgraph`, a new path; `chrom`, `start` and `end` to restrict the region (`start` or `end` without `chrom` raises `ConfigurationError`); `udc_dir`, a cache folder name passed to `-udcDir`.
- **Output:** A dictionary with `status`, `backend`, `executable`, `command`, `return_code`, `outputs` and `log_tail`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.annotation import bigwig_to_bedgraph
from dnakit.exceptions import BackendUnavailableError

work = Path(tempfile.mkdtemp())
(work / "signal.bw").write_bytes(b"")
try:
    bigwig_to_bedgraph(
        str(work / "signal.bw"),
        str(work / "signal.bedgraph"),
        bigwigtobedgraph_path="/nonexistent/bigWigToBedGraph",
        verbose=False,
    )
except BackendUnavailableError as error:
    print(type(error).__name__, "bigWigToBedGraph path does not exist" in str(error))
```

- **Example output:**

```text
BackendUnavailableError True
```

The program is not installed here, so the example shows the error for a program path that does not exist. An empty placeholder file stands in for the bigWig, because the program is looked up before the file is opened.

Usage with the dependency installed:

```python
from dnakit.annotation import bigwig_to_bedgraph

result = bigwig_to_bedgraph("signal.bw", "signal_chr1.bedgraph", chrom="chr1", start=0, end=300)
print(result["status"], result["outputs"])
```

## 7) `dnakit.annotation.binary_heatmap()`

- **Purpose:** Draw a binary heat map of which BED files contain each interval and return the class counts.
- **Method:** Run `BedTool().multi_intersect(i=beds, names=names, cluster=True)` (bedtools multiinter), count the class label of each row (for example `a` or `a,b`), build the 0/1 presence matrix from the per-file columns, order the rows with the same column-weighted score as `sort_binary_matrix()`, and draw the matrix with matplotlib `imshow` (colormap `binary`) into a PNG at 120 dpi.
- **API:** `binary_heatmap(beds: 'Sequence[str]', names: 'Sequence[str]', output_png: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages pybedtools (extra `dnakit[bio-pybedtools]`) and matplotlib (extra `dnakit[bio-plots]`), and the external program bedtools on PATH.
- **Input:** `beds`, 2 to 30 existing BED files; `names`, one unique name per file (letters, digits, `_`, `.`, `-` and spaces, at most 40 characters); `output_png`, a new path. If the files share no intervals a `ConfigurationError` is raised.
- **Output:** A dictionary with `status`, `backend`, `image_path` and `data`, where `data` has `names`, `interval_count`, `class_counts` (sorted by count, then name), `matrix` and `row_order`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.annotation import binary_heatmap

out = Path(tempfile.mkdtemp()) / "heatmap.png"
result = binary_heatmap(["peaks_a.bed", "peaks_b.bed"], ["a", "b"], str(out), verbose=False)
data = result["data"]
print(data["interval_count"], data["class_counts"])
print(data["matrix"], data["row_order"])
print(out.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n")
```

- **Example output:**

```text
4 {'a,b': 2, 'a': 1, 'b': 1}
[[1, 1], [1, 0], [0, 1], [1, 1]] [3, 0, 1, 2]
True
```

## 8) `dnakit.annotation.chromsizes()`

- **Purpose:** Return the chromosome sizes of a genome assembly.
- **Method:** Call `pybedtools.chromsizes(genome)`, which uses its built-in registry of assemblies first and otherwise genomepy and then the UCSC MySQL server. The `service` field says which path applies: the registry is offline, any other name needs genomepy or the network. Sizes are returned as a name-to-length map; with `output_path` they are also written as a chromosome sizes file.
- **API:** `chromsizes(genome: 'str', *, output_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package pybedtools (extra `dnakit[bio-pybedtools]`). Assemblies in the built-in registry (for example `hg19`) work offline; other assemblies need genomepy or internet access to UCSC.
- **Input:** `genome`, an assembly name made of letters, digits, `_`, `.` and `-` (at most 64 characters); `output_path`, an optional new file path. A failed lookup raises `BackendExecutionError`.
- **Output:** A dictionary with `status`, `service`, `request`, `retrieved_at` (UTC time), `pybedtools_version`, `row_count`, `data` (chromosome name to length) and `output_path`.
- **Example code:**

```python
from dnakit.annotation import chromsizes

result = chromsizes("hg19", verbose=False)
print(result["status"], result["service"])
print(result["row_count"], result["data"]["chr1"], result["data"]["chrM"])
```

- **Example output:**

```text
completed pybedtools built-in genome registry (offline)
93 249250621 16571
```

## 9) `dnakit.annotation.cleaned_intersect()`

- **Purpose:** Build the interval-string sets of 2 to 4 BED files in which overlapping intervals are made identical.
- **Method:** Call `pybedtools.contrib.venn_maker.cleaned_intersect()` on the files, which intersects every overlapping pair so that overlapping intervals of different files become the same interval. Each interval is then turned into a string with tabs replaced by `|`, and each set is sorted and de-duplicated.
- **API:** `cleaned_intersect(beds: 'Sequence[str]', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package pybedtools (extra `dnakit[bio-pybedtools]`) and the external program bedtools on PATH.
- **Input:** `beds`, a list of 2 to 4 existing BED files.
- **Output:** A dictionary with `status`, `backend`, `pybedtools_version`, `files`, `sets` (one sorted list of strings per file) and `set_sizes`.
- **Example code:**

```python
from dnakit.annotation import cleaned_intersect

result = cleaned_intersect(["peaks_a.bed", "peaks_b.bed"], verbose=False)
print(result["set_sizes"])
for interval_set in result["sets"]:
    print(interval_set)
```

- **Example output:**

```text
[3, 3]
['chr1|100|200', 'chr1|300|400', 'chr2|50|150']
['chr1|100|200', 'chr1|700|800', 'chr2|50|150']
```

## 10) `dnakit.annotation.get_chromsizes_from_genomepy()`

- **Purpose:** Read the chromosome sizes of a genome that is installed with genomepy.
- **Method:** Call `pybedtools.helpers.get_chromsizes_from_genomepy(genome)`. When genomepy has no sizes file for the genome the function does not raise an error: it returns `status` `not_found` with a message and `data` set to None.
- **API:** `get_chromsizes_from_genomepy(genome: 'str', *, output_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages genomepy (extra `dnakit[bio-genomepy]`) and pybedtools (extra `dnakit[bio-pybedtools]`). The genome must already be installed with genomepy; nothing is downloaded by this function.
- **Input:** `genome`, an assembly name; `output_path`, an optional new file path (written only when sizes are found). Any other genomepy failure raises `BackendExecutionError`.
- **Output:** A dictionary with `status` (`completed` or `not_found`), `service`, `request`, `retrieved_at`, `pybedtools_version`, `row_count`, `data` (chromosome name to length, or None), `output_path` and, for `not_found`, a `message`.
- **Example code:**

```python
from dnakit.annotation import get_chromsizes_from_genomepy

result = get_chromsizes_from_genomepy("no_such_genome_demo", verbose=False)
print(result["status"], result["row_count"], result["data"])
print(result["message"])
```

- **Example output:**

```text
not_found 0 None
genomepy has no local sizes file for this genome; install it with genomepy install.
```

No genome of this name is installed, so the result is `not_found`.

Usage with the dependency installed:

```python
from dnakit.annotation import get_chromsizes_from_genomepy

# after: genomepy install hg38
result = get_chromsizes_from_genomepy("hg38", output_path="hg38.chrom.sizes")
print(result["status"], result["row_count"])
```

## 11) `dnakit.annotation.get_chromsizes_from_ucsc()`

- **Purpose:** Download the chromosome sizes of a UCSC assembly from the UCSC MySQL server.
- **Method:** Call `pybedtools.get_chromsizes_from_ucsc(genome, saveas=..., timeout=..., host_url=...)`, which runs the `mysql` and `fetchChromSizes` programs against the UCSC server. Network and program failures are re-raised as `BackendExecutionError`.
- **API:** `get_chromsizes_from_ucsc(genome: 'str', *, host_url: 'str' = 'genome-mysql.cse.ucsc.edu', timeout: 'float | None' = None, output_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package pybedtools (extra `dnakit[bio-pybedtools]`), internet access to the UCSC MySQL server (no key), and the programs `mysql` and `fetchChromSizes` on PATH, which pybedtools starts for the lookup.
- **Input:** `genome`, an assembly name such as `hg38` (letters, digits, `_`, `.` and `-`); `host_url`, a host name; `timeout`, seconds in (0, 3600] or None; `output_path`, an optional new file path.
- **Output:** A dictionary with `status`, `service`, `request` (with `genome` and `host_url`), `retrieved_at`, `pybedtools_version`, `row_count`, `data` and `output_path`.
- **Example code:**

```python
from dnakit.annotation import get_chromsizes_from_ucsc
from dnakit.exceptions import ConfigurationError

try:
    get_chromsizes_from_ucsc("hg38; drop", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, "genome must be an assembly name" in str(error))
```

- **Example output:**

```text
ConfigurationError True
```

The example does not touch the network. It shows how an invalid assembly name is rejected before any request is made.

Usage with the dependency installed:

```python
from dnakit.annotation import get_chromsizes_from_ucsc

result = get_chromsizes_from_ucsc("hg38", timeout=60, output_path="hg38.chrom.sizes")
print(result["row_count"], result["data"]["chr1"])
```

## 12) `dnakit.annotation.mapped_read_count()`

- **Purpose:** Count the mapped reads of a BAM file.
- **Method:** Run `samtools view -c -F 0x4 <bam>` and read the last output line as an integer. Unlike pybedtools, which caches the count in a `<bam>.scale` file next to the BAM, this function does not read or write that cache, so the input folder stays unchanged.
- **API:** `mapped_read_count(bam: 'str', *, samtools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program samtools on PATH (or `samtools_path`). It needs no Python package.
- **Input:** `bam`, an existing BAM file; optional `samtools_path`.
- **Output:** A dictionary with `status` (`completed` or `failed`), `backend`, `executable`, `command`, `return_code`, `mapped_reads` (an integer, or None when the count could not be read), `outputs` (empty) and `log_tail`.
- **Example code:**

```python
import tempfile
from pathlib import Path

import pysam

from dnakit.annotation import mapped_read_count
from dnakit.exceptions import BackendUnavailableError

work = Path(tempfile.mkdtemp())
pysam.view("-b", "-o", str(work / "reads.bam"), "reads.sam", catch_stdout=False)
try:
    mapped_read_count(str(work / "reads.bam"), samtools_path="/nonexistent/samtools", verbose=False)
except BackendUnavailableError as error:
    print(type(error).__name__, "samtools path does not exist" in str(error))
```

- **Example output:**

```text
BackendUnavailableError True
```

samtools is not installed here, so the example shows the error for a program path that does not exist. The BAM is made from `reads.sam` with pysam.

Usage with the dependency installed:

```python
from dnakit.annotation import mapped_read_count

result = mapped_read_count("reads.bam")
print(result["status"], result["mapped_reads"])
```

## 13) `dnakit.annotation.sort_binary_matrix()`

- **Purpose:** Return the row order of a 0/1 presence matrix, sorted by a column-weighted score with the highest first.
- **Method:** Give column j (counting from 1 at the left, with C columns) the weight `2 ** (C - j + 1)`, so the first column weighs most. Sum the weighted values of every row and return the row indices with `numpy.argsort(...)[::-1]`, as pybedtools `sort_binary_matrix` does. An empty matrix returns an empty list.
- **API:** `sort_binary_matrix(matrix: 'Sequence[Sequence[int]]', *, verbose: bool = True) -> 'list[int]'`.
- **Requirements:** Needs numpy (installed with the `dnakit[bio-pybedtools]` environment; the function reports a missing numpy with that extra).
- **Input:** `matrix`, a list of equal-length, non-empty rows that contain only the integers 0 and 1 (booleans are rejected); at most 100,000 rows. Anything else raises `ConfigurationError`.
- **Output:** A list of row indices (Python integers).
- **Example code:**

```python
from dnakit.annotation import sort_binary_matrix

matrix = [
    [0, 1, 1],
    [1, 0, 0],
    [1, 1, 0],
    [0, 0, 1],
]
print(sort_binary_matrix(matrix, verbose=False))
```

- **Example output:**

```text
[2, 1, 0, 3]
```

## 14) `dnakit.annotation.venn_maker()`

- **Purpose:** Write an R script that draws a Venn diagram of 2 to 4 BED files with the VennDiagram package, and render it when Rscript is available.
- **Method:** Build the interval sets with `cleaned_intersect()`, call `pybedtools.contrib.venn_maker.venn_maker(..., run=False)` to write the R script, and count the items in each exact combination of sets (for example `a` or `a&b`). If `render` is true, a `figure_path` is given and `Rscript` is on PATH, run the script; otherwise only the script is written and `render_note` says why.
- **API:** `venn_maker(beds: 'Sequence[str]', output_script: 'str', *, names: 'Sequence[str] | None' = None, figure_path: 'str | None' = None, additional_args: 'Sequence[str] | None' = None, render: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package pybedtools (extra `dnakit[bio-pybedtools]`) and the external program bedtools on PATH. Rendering the picture additionally needs R with the VennDiagram package and `Rscript` on PATH; writing the script and the counts does not.
- **Input:** `beds`, 2 to 4 existing BED files; `output_script`, a new path for the R script; `names`, unique labels (default `a`, `b`, `c`, `d`); `figure_path`, an optional new path without quotes or backslashes; `additional_args`, R named arguments such as `euler.d=TRUE` (no semicolons or line breaks); `render`, a boolean.
- **Output:** A dictionary with `status`, `backend`, `script_path`, `figure_path`, `rendered`, `render_note`, `render` (the Rscript return code and log tail, or None) and `data` with `names`, `set_sizes` and `regions`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.annotation import venn_maker

script = Path(tempfile.mkdtemp()) / "venn.R"
result = venn_maker(
    ["peaks_a.bed", "peaks_b.bed"],
    str(script),
    names=["A", "B"],
    render=False,
    verbose=False,
)
print(result["data"])
print(result["rendered"], result["render_note"])
print(script.is_file(), "venn.diagram" in script.read_text())
```

- **Example output:**

```text
{'names': ['A', 'B'], 'set_sizes': {'A': 3, 'B': 3}, 'regions': {'A': 1, 'A&B': 2, 'B': 1}}
False Rendering was not requested; the R script was written only.
True True
```

The example uses `render=False`, so the output is the same with or without R.

Usage with the dependency installed:

```python
from dnakit.annotation import venn_maker

result = venn_maker(["peaks_a.bed", "peaks_b.bed"], "venn.R", names=["A", "B"], figure_path="venn.png")
print(result["rendered"], result["render_note"])
```

## 15) `dnakit.annotation.wig_to_bigwig()`

- **Purpose:** Convert a wiggle file to an indexed bigWig file.
- **Method:** Run UCSC `wigToBigWig <wig> <chromsizes> <output>` through DNAKit's bounded command runner, without a shell.
- **API:** `wig_to_bigwig(wig: 'str', chromsizes: 'str', output_bigwig: 'str', *, wigtobigwig_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program UCSC wigToBigWig on PATH (or `wigtobigwig_path`). It needs no Python package.
- **Input:** `wig`, an existing wiggle file; `chromsizes`, a chromosome sizes file; `output_bigwig`, a new path; optional `wigtobigwig_path`.
- **Output:** A dictionary with `status`, `backend`, `executable`, `command`, `return_code`, `outputs` and `log_tail`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.annotation import wig_to_bigwig
from dnakit.exceptions import BackendUnavailableError

out = Path(tempfile.mkdtemp()) / "signal.bw"
try:
    wig_to_bigwig(
        "signal.wig",
        "chrom.sizes",
        str(out),
        wigtobigwig_path="/nonexistent/wigToBigWig",
        verbose=False,
    )
except BackendUnavailableError as error:
    print(type(error).__name__, "wigToBigWig path does not exist" in str(error))
```

- **Example output:**

```text
BackendUnavailableError True
```

The program is not installed here, so the example shows the error for a program path that does not exist.

Usage with the dependency installed:

```python
from dnakit.annotation import wig_to_bigwig

result = wig_to_bigwig("signal.wig", "chrom.sizes", "signal.bw")
print(result["status"], result["outputs"])
```

## 16) `dnakit.annotation.fetch_centromeres()`

- **Purpose:** Fetch the centromere position of every chromosome of an assembly.
- **Method:** Call `bioframe.fetch_centromeres(db, provider=provider)`. With `provider="ucsc"` the tables come from the UCSC Genome Browser; with `provider="local"` bioframe's bundled tables are used and no network is needed. The rows are returned as records, and with `output_path` the table is also written as a tab-separated file with a header.
- **API:** `fetch_centromeres(db: 'str', *, provider: 'str' = 'ucsc', output_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package bioframe (extra `dnakit[bio-bioframe]`). `provider="local"` works offline for the assemblies that bioframe bundles; `provider="ucsc"` (the default) needs internet access.
- **Input:** `db`, an assembly name such as `hg38`; `provider`, `ucsc` or `local`; `output_path`, an optional new file path. A lookup failure (`OSError`, `ValueError` or `KeyError` from bioframe) raises `BackendExecutionError`.
- **Output:** A dictionary with `status`, `service`, `request`, `retrieved_at`, `bioframe_version`, `row_count`, `data` (a list of records with the columns bioframe returns, for example `chrom`, `start`, `end`, `mid`) and `output_path`.
- **Example code:**

```python
from dnakit.annotation import fetch_centromeres

result = fetch_centromeres("hg38", provider="local", verbose=False)
print(result["service"])
print(result["row_count"], result["data"][0])
```

- **Example output:**

```text
bioframe bundled assembly tables (no network access)
24 {'chrom': 'chr1', 'start': 121700000, 'end': 125100000, 'mid': 123400000}
```

## 17) `dnakit.annotation.fetch_chromsizes()`

- **Purpose:** Fetch the chromosome lengths of an assembly as a name-to-length map or as a BED3 table.
- **Method:** Call `bioframe.fetch_chromsizes(db, provider=..., as_bed=..., filter_chroms=...)`. `filter_chroms=True` keeps the primary chromosomes only. With `output_path` the table is written without a header (a two-column sizes file, or BED3 when `as_bed` is true).
- **API:** `fetch_chromsizes(db: 'str', *, provider: 'str' = 'ucsc', as_bed: 'bool' = False, filter_chroms: 'bool' = True, output_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package bioframe (extra `dnakit[bio-bioframe]`). `provider="local"` works offline for the assemblies that bioframe bundles; `provider="ucsc"` (the default) needs internet access.
- **Input:** `db`, an assembly name; `provider`, `ucsc` or `local`; `as_bed` and `filter_chroms`, booleans; `output_path`, an optional new file path. A lookup failure raises `BackendExecutionError`.
- **Output:** A dictionary with `status`, `service`, `request`, `retrieved_at`, `bioframe_version`, `row_count`, `data` (a chromosome-to-length map, or a list of `chrom`/`start`/`end` records when `as_bed` is true) and `output_path`.
- **Example code:**

```python
from dnakit.annotation import fetch_chromsizes

sizes = fetch_chromsizes("hg38", provider="local", verbose=False)
print(sizes["row_count"], sizes["data"]["chr1"], sizes["data"]["chrX"])
bed = fetch_chromsizes("hg38", provider="local", as_bed=True, verbose=False)
print(bed["data"][0])
```

- **Example output:**

```text
25 248956422 156040895
{'chrom': 'chr1', 'start': 0, 'end': 248956422}
```

## 18) `dnakit.annotation.plot_intervals()`

- **Purpose:** Plot a BED-like interval file as one stacked track per chromosome in a PNG.
- **Method:** Read the tab-separated file without a header (columns 0 to 2 are chromosome, start and end), group the rows by chromosome, assign each interval a row with `bioframe.core.arrops.stack_intervals` so that overlapping intervals do not collide, and draw each interval as a rectangle on a matplotlib Figure (not pyplot), with an optional text label from a chosen column.
- **API:** `plot_intervals(intervals: 'str', output_png: 'str', *, label_column: 'int | None' = None, color: 'str' = 'skyblue', show_coords: 'bool' = False, figsize: 'Sequence[float]' = (10.0, 2.0), verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages bioframe (extra `dnakit[bio-bioframe]`) and matplotlib (extra `dnakit[bio-plots]`), plus pandas and numpy, which bioframe installs.
- **Input:** `intervals`, an existing tab-separated file with at least three columns, integer start and end, `0 <= start <= end`, and 1 to 100,000 rows; `output_png`, a new path; `label_column`, a column index or None; `color`, a matplotlib colour name; `show_coords`, a boolean; `figsize`, two positive numbers (width, height per track).
- **Output:** A dictionary with `status`, `backend`, `image_path` and `data`, where `data` has `source`, `interval_count`, `chromosomes`, `tracks` (chromosome, interval count and highest stack level) and `intervals` (chromosome, start, end, level and label for each interval).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.annotation import plot_intervals

out = Path(tempfile.mkdtemp()) / "intervals.png"
result = plot_intervals("peaks_a.bed", str(out), label_column=3, verbose=False)
data = result["data"]
print(data["interval_count"], data["chromosomes"], data["tracks"])
print(data["intervals"][0])
print(out.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n")
```

- **Example output:**

```text
3 ['chr1', 'chr2'] [{'chrom': 'chr1', 'interval_count': 2, 'max_level': 0}, {'chrom': 'chr2', 'interval_count': 1, 'max_level': 0}]
{'chrom': 'chr1', 'start': 100, 'end': 200, 'level': 0, 'label': 'a1'}
True
```

## 19) `dnakit.annotation.plot_intervals_arr()`

- **Purpose:** Plot one set of intervals given as start and end arrays in a PNG.
- **Method:** Use the stack levels you give, or compute them with `bioframe.core.arrops.stack_intervals`, and draw one rectangle per interval with an optional centred label on a matplotlib Figure. The x range is `-0.5` to `max(ends) + 0.5` unless `xlim` is given.
- **API:** `plot_intervals_arr(starts: 'Sequence[float]', ends: 'Sequence[float]', output_png: 'str', *, levels: 'Sequence[float] | None' = None, labels: 'Sequence[str | None] | None' = None, colors: 'str | Sequence[str] | None' = None, xlim: 'Sequence[float] | None' = None, show_coords: 'bool' = False, figsize: 'Sequence[float]' = (10.0, 2.0), verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages bioframe (extra `dnakit[bio-bioframe]`) and matplotlib (extra `dnakit[bio-plots]`), plus numpy, which bioframe installs.
- **Input:** `starts` and `ends`, equal-length non-empty lists of finite numbers with `start <= end`; `output_png`, a new path; optional `levels` (one per interval), `labels` (one string or None per interval), `colors` (one colour name or one per interval), `xlim` (two numbers, lower first), `show_coords` and `figsize`.
- **Output:** A dictionary with `status`, `backend`, `image_path` and `data`, where `data` has `interval_count`, `starts`, `ends`, `levels`, `labels`, `colors`, `xlim` and `ylim`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.annotation import plot_intervals_arr

out = Path(tempfile.mkdtemp()) / "arr.png"
result = plot_intervals_arr(
    [0, 10, 25],
    [20, 30, 40],
    str(out),
    labels=["x", "y", None],
    colors="salmon",
    verbose=False,
)
data = result["data"]
print(data["levels"], data["labels"], data["xlim"], data["ylim"])
print(out.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n")
```

- **Example output:**

```text
[0, 1, 0] ['x', 'y', None] [-0.5, 40.5] [-0.5, 1.5]
True
```

## 20) `dnakit.annotation.read_pairix()`

- **Purpose:** Read the pairs of a pairix-indexed file that fall in one region, or between two regions.
- **Method:** Call `bioframe.read_pairix(path, region1, region2=..., chromsizes=..., usecols=...)`, which queries the index with pypairix, and return the table as records.
- **API:** `read_pairix(path: 'str', region1: 'str', *, region2: 'str | None' = None, chromsizes: 'Mapping[str, int] | None' = None, usecols: 'Sequence[str] | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages pypairix and cytoolz (extra `dnakit[bio-pairix]`) and bioframe (extra `dnakit[bio-bioframe]`). The file must be bgzip-compressed and pairix-indexed.
- **Input:** `path`, an existing indexed pairs file; `region1`, a region string such as `chr1:1-1000`; `region2`, an optional second region; `chromsizes`, an optional map of chromosome name to positive integer length; `usecols`, an optional non-empty list of column names. pypairix and read failures raise `BackendExecutionError`.
- **Output:** A dictionary with `status`, `backend`, `bioframe_version`, `request`, `row_count`, `columns` and `data` (a list of records).
- **Example code:**

```python
from dnakit.annotation import read_pairix
from dnakit.exceptions import BackendUnavailableError

try:
    read_pairix("pairs.txt", "chr1:1-1000", verbose=False)
except BackendUnavailableError as error:
    print(type(error).__name__, "pypairix" in str(error))
```

- **Example output:**

```text
BackendUnavailableError True
```

pypairix is not installed here, so the example shows `BackendUnavailableError`. The plain text file stands in for an indexed file, because the check for pypairix happens before the file is read.

Usage with the dependency installed:

```python
from dnakit.annotation import read_pairix

# pairs.sorted.gz: bgzip-compressed and indexed with `pairix -p pairs`
result = read_pairix("pairs.sorted.gz", "chr1:1-1000", region2="chr2:1-1000")
print(result["row_count"], result["columns"])
```

## 21) `dnakit.annotation.to_bigbed()`

- **Purpose:** Write a sorted BED file as an indexed bigBed file with an explicit chromosome sizes file.
- **Method:** Run UCSC `bedToBigBed [-type=<bed_format>] <bed> <chromsizes> <output>` through DNAKit's bounded command runner, without a shell. Unlike `bigbed()`, it passes no block-size or autoSql options and does not guess the type.
- **API:** `to_bigbed(bed: 'str', chromsizes: 'str', output_bigbed: 'str', *, bed_format: 'str | None' = None, bedtobigbed_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program UCSC bedToBigBed on PATH (or `bedtobigbed_path`). It needs no Python package.
- **Input:** `bed`, an existing BED file sorted by chromosome and start; `chromsizes`, a chromosome sizes file; `output_bigbed`, a new path; `bed_format`, an optional type such as `bed3`, `bed6` or `bed6+4`.
- **Output:** A dictionary with `status`, `backend`, `executable`, `command`, `return_code`, `outputs` and `log_tail`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.annotation import to_bigbed
from dnakit.exceptions import BackendUnavailableError

out = Path(tempfile.mkdtemp()) / "peaks.bb"
try:
    to_bigbed(
        "peaks_a.bed",
        "chrom.sizes",
        str(out),
        bed_format="bed6",
        bedtobigbed_path="/nonexistent/bedToBigBed",
        verbose=False,
    )
except BackendUnavailableError as error:
    print(type(error).__name__, "bedToBigBed path does not exist" in str(error))
```

- **Example output:**

```text
BackendUnavailableError True
```

The program is not installed here, so the example shows the error for a program path that does not exist.

Usage with the dependency installed:

```python
from dnakit.annotation import to_bigbed

result = to_bigbed("peaks_a.bed", "chrom.sizes", "peaks.bb", bed_format="bed6")
print(result["status"], result["command"])
```

## 22) `dnakit.annotation.to_bigwig()`

- **Purpose:** Write a sorted bedGraph file as an indexed bigWig file with an explicit chromosome sizes file.
- **Method:** Run UCSC `bedGraphToBigWig <bedgraph> <chromsizes> <output>` through DNAKit's bounded command runner, without a shell. It does the same as `bedgraph_to_bigwig()`; it is the bioframe-named counterpart of that function.
- **API:** `to_bigwig(bedgraph: 'str', chromsizes: 'str', output_bigwig: 'str', *, bedgraphtobigwig_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the external program UCSC bedGraphToBigWig on PATH (or `bedgraphtobigwig_path`). It needs no Python package.
- **Input:** `bedgraph`, an existing bedGraph sorted by chromosome and start; `chromsizes`, a chromosome sizes file; `output_bigwig`, a new path; optional `bedgraphtobigwig_path`.
- **Output:** A dictionary with `status`, `backend`, `executable`, `command`, `return_code`, `outputs` and `log_tail`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.annotation import to_bigwig
from dnakit.exceptions import BackendUnavailableError

out = Path(tempfile.mkdtemp()) / "signal.bw"
try:
    to_bigwig(
        "signal.bedgraph",
        "chrom.sizes",
        str(out),
        bedgraphtobigwig_path="/nonexistent/bedGraphToBigWig",
        verbose=False,
    )
except BackendUnavailableError as error:
    print(type(error).__name__, "bedGraphToBigWig path does not exist" in str(error))
```

- **Example output:**

```text
BackendUnavailableError True
```

The program is not installed here, so the example shows the error for a program path that does not exist.

Usage with the dependency installed:

```python
from dnakit.annotation import to_bigwig

result = to_bigwig("signal.bedgraph", "chrom.sizes", "signal.bw")
print(result["status"], result["outputs"])
```

## 23) `dnakit.annotation.read_pyranges1_bam()`

- **Purpose:** Read BAM alignments as pyranges1 intervals with the columns `Chromosome`, `Start`, `End`, `Strand` and `Flag`.
- **Method:** Read the BAM with pysam (`pyranges1.read_bam` needs a separate bamread package, so it is not called), apply the filters, take `Start` from the 0-based reference start and `End` from the CIGAR reference span (end exclusive), set `Strand` to `-` for flag 16, and build a `pyranges1.PyRanges` from the table. The returned records are the first `max_rows` intervals.
- **API:** `read_pyranges1_bam(path: 'str | Path', *, mapq: 'int' = 0, required_flag: 'int' = 0, filter_flag: 'int' = 1540, sparse: 'bool' = True, max_rows: 'int' = 1000, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages pysam (extra `dnakit[bio-reads]`), pyranges1 and pandas (extra `dnakit[bio-pyranges1]`).
- **Input:** `path`, an existing BAM file; `mapq` (0 to 255), keep reads with at least this mapping quality; `required_flag`, keep reads that have all of these SAM flag bits; `filter_flag`, drop reads that have any of these bits (default 1540: unmapped, QC fail, duplicate); `sparse`, a boolean (`False` adds `QueryStart`, `QueryEnd`, `QuerySequence`, `Name`, `Cigar` and `Quality`); `max_rows`, 1 to 1,000,000. The whole file is read; reads without a reference span are skipped.
- **Output:** A dictionary with `status`, `operation`, `provenance` (package, backend and the coordinate convention), `columns`, `record_count` (all kept intervals), `truncated` (true when more than `max_rows` were kept) and `records`.
- **Example code:**

```python
import tempfile
from pathlib import Path

import pysam

from dnakit.annotation import read_pyranges1_bam

work = Path(tempfile.mkdtemp())
pysam.view("-b", "-o", str(work / "reads.bam"), "reads.sam", catch_stdout=False)
result = read_pyranges1_bam(str(work / "reads.bam"), verbose=False)
print(result["columns"], result["record_count"], result["truncated"])
for record in result["records"]:
    print(record)
strict = read_pyranges1_bam(str(work / "reads.bam"), mapq=50, verbose=False)
print(strict["record_count"])
```

- **Example output:**

```text
['Chromosome', 'Start', 'End', 'Strand', 'Flag'] 2 False
{'Chromosome': 'chr1', 'Start': 100, 'End': 110, 'Strand': '+', 'Flag': 0}
{'Chromosome': 'chr1', 'Start': 200, 'End': 208, 'Strand': '-', 'Flag': 16}
1
```

The BAM is made from `reads.sam` with pysam. Read `r3` is a duplicate and `r4` is unmapped, so the default filter drops both.

## 24) `dnakit.annotation.read_ucsc_gtf()`

- **Purpose:** Read a UCSC gene-model GTF into a table and write it as a TSV file.
- **Method:** Parse the GTF with `genomicranges.io.gtf.parse_gtf`, write the table with `to_csv(sep="\t", index=False)` and summarise it. Without `local_gtf` the file is downloaded from the address given by `ucsc_gtf_url()`; with `local_gtf` (plain or `.gz`) the local file is read, nothing is downloaded, and `genome` and `annotation` only label the result. `source_date` is the file's modification date for a local file and today's date for a download.
- **API:** `read_ucsc_gtf(genome: 'str', output_table: 'str', *, annotation: 'str' = 'refGene', local_gtf: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional packages genomicranges and pandas (extra `dnakit[bio-genomicranges]`). Without `local_gtf` it also needs internet access to UCSC.
- **Input:** `genome`, a UCSC genome code such as `hg38` (no `/`, not starting with `-`); `output_table`, a new file in an existing folder; `annotation`, one of `refGene`, `ensGene`, `knownGene` and `ncbiRefSeq`; `local_gtf`, an optional existing file.
- **Output:** A dictionary with `status`, `backend`, `genome`, `annotation`, `source`, `source_date`, `row_count`, `columns`, `rows_per_chromosome` and `output_table`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.annotation import read_ucsc_gtf

out = Path(tempfile.mkdtemp()) / "genes.tsv"
result = read_ucsc_gtf("hg38", str(out), local_gtf="genes.gtf", verbose=False)
print(result["row_count"], result["rows_per_chromosome"])
print(result["columns"])
print(len(out.read_text().splitlines()))
```

- **Example output:**

```text
3 {'chr1': 2, 'chr2': 1}
['seqnames', 'source', 'feature', 'starts', 'ends', 'score', 'strand', 'frame', 'gene_id', 'transcript_id']
4
```

The example reads the local `genes.gtf`, so nothing is downloaded.

## 25) `dnakit.annotation.ucsc_gtf_url()`

- **Purpose:** Return the UCSC download address of a genome's gene-model GTF.
- **Method:** Build the string `http://hgdownload.cse.ucsc.edu/goldenPath/<genome>/bigZips/genes/<genome>.<annotation>.gtf.gz`. The address that genomicranges builds has a doubled slash before the file name; this function writes a single slash. No request is made.
- **API:** `ucsc_gtf_url(genome: 'str', annotation: 'str' = 'refGene', *, verbose: bool = True) -> 'str'`.
- **Requirements:** Works after plain pip install (it builds a string and imports no optional package).
- **Input:** `genome`, a UCSC genome code (not empty, no `/`, not starting with `-`); `annotation`, one of `refGene`, `ensGene`, `knownGene` and `ncbiRefSeq`. Other values raise `ConfigurationError`.
- **Output:** The address as a string.
- **Example code:**

```python
from dnakit.annotation import ucsc_gtf_url

print(ucsc_gtf_url("hg38", verbose=False))
print(ucsc_gtf_url("mm10", "ncbiRefSeq", verbose=False))
```

- **Example output:**

```text
http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz
http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/genes/mm10.ncbiRefSeq.gtf.gz
```

## 26) `dnakit.annotation.hgvs_to_vrs_allele_id()`

- **Purpose:** Map an HGVS variant to its GA4GH VRS allele identifier.
- **Method:** Open the transcript store (`gtars.reftx.ReadonlyTxStore.open`) and wrap it in a `ReftxProvider`, open the refget store directory with `gtars.refget.RefgetStore.open_local`, and call `gtars.vrs.hgvs.hgvs_to_vrs_id`. A reference mismatch (the HGVS reference base differs from the sequence) or any other mapping failure is raised as `ConfigurationError`; nothing is downloaded.
- **API:** `hgvs_to_vrs_allele_id(hgvs: 'str', transcript_store: 'str', refget_store: 'str', collection_digest: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package gtars (extra `dnakit[bio-gtars]`) and two local inputs that you build once: a transcript store file and a refget store directory.
- **Input:** `hgvs`, an HGVS string such as `NM_TEST.1:c.4G>A`; `transcript_store`, an existing store file; `refget_store`, an existing store directory; `collection_digest`, the digest of the reference collection in that store.
- **Output:** A dictionary with `status`, `backend`, `version`, `hgvs`, `vrs_allele_id` (`ga4gh:VA.<digest>`), `transcript_store`, `refget_store` and `collection_digest`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from gtars import refget, reftx

from dnakit.annotation import hgvs_to_vrs_allele_id
from dnakit.exceptions import ConfigurationError

work = Path(tempfile.mkdtemp())
store = refget.RefgetStore.on_disk(str(work / "refget"))
store.set_quiet(True)
meta, _ = store.add_sequence_collection_from_fasta("tiny.fa")
sequence = store.list_sequences()[0]
accession = "SQ." + sequence.sha512t24u

builder = reftx.TxStoreBuilder()
builder.add_chrom_mapping("chr1", accession)
builder.add_transcript(
    {
        "accession": "NM_TEST.1",
        "gene": "TEST",
        "chrom": accession,
        "strand": "+",
        "cds_start": 0,
        "cds_end": 39,
        "exons": [{"start": 0, "end": 39}],
        "mane": None,
    }
)
builder.build(str(work / "tx.reftx"))

args = (str(work / "tx.reftx"), str(work / "refget"), meta.digest)
result = hgvs_to_vrs_allele_id("NM_TEST.1:c.4G>A", *args, verbose=False)
print(result["status"], result["vrs_allele_id"])
try:
    hgvs_to_vrs_allele_id("NM_TEST.1:c.4C>A", *args, verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, "REF mismatch" in str(error))
```

- **Example output:**

```text
completed ga4gh:VA.931BbxClFtmksLkac20ozAYn01o1uv_f
ConfigurationError True
```

The example builds a one-sequence refget store from `tiny.fa` and a one-transcript store with gtars, then maps `c.4G>A` (the fourth base of the coding sequence, which is G) and a deliberately wrong reference base.

## 27) `dnakit.visualization.compare_genomes_mmseqs()`

- **Purpose:** Link neighbouring genomes by MMseqs2 reciprocal best protein hits and draw the links in `result.png`.
- **Method:** Read the GenBank files with pyGenomeViz, write the translated CDS proteins of each genome to `<name>.faa`, and for every pair of consecutive genomes run `mmseqs easy-rbh <a.faa> <b.faa> <result.tsv> <tmp> --threads N -e EVALUE -v 0`. The hits become alignment coordinates with their identity in percent; all of them are written to `align_coords.tsv`, and those with identity of at least `min_identity` are drawn as links between CDS tracks with an identity colour bar. The temporary MMseqs2 folders are removed.
- **API:** `compare_genomes_mmseqs(gbk_resources: 'Sequence[str]', output_dir: 'str', *, evalue: 'float' = 0.001, min_identity: 'float' = 0.0, threads: 'int' = 1, fig_width: 'float' = 15.0, dpi: 'int' = 300, timeout_seconds: 'float' = 3600.0, mmseqs_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package pygenomeviz (extra `dnakit[bio-pygenomeviz]`) and the external program MMseqs2 (`mmseqs`) on PATH or `mmseqs_path`.
- **Input:** `gbk_resources`, at least two GenBank files (the syntax `file:start-end` selects a range; genome names, taken from the records, may contain only letters, digits, `.`, `_` and `-` and must be unique); `output_dir`, a new or empty folder; `evalue` and `fig_width`, positive numbers; `min_identity`, a percentage from 0 to 100; `threads` and `dpi`, positive integers; `timeout_seconds`, a limit per program call. A genome without any CDS that has a translation raises `ConfigurationError`.
- **Output:** A dictionary with `status` (`completed`, or `failed` when a program call fails), `backend`, `commands`, `return_codes`, `outputs` (the hit tables, `align_coords.tsv` and `result.png`), `image_path`, `log_tail` and `data` (`genomes`, `alignment_count` and `alignments`).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.visualization import compare_genomes_mmseqs

out = Path(tempfile.mkdtemp()) / "mmseqs_out"
try:
    compare_genomes_mmseqs(
        ["genomeA.gbk", "genomeB.gbk"],
        str(out),
        mmseqs_path="/nonexistent/mmseqs",
        verbose=False,
    )
except BackendUnavailableError as error:
    print(type(error).__name__, "mmseqs was not found" in str(error), out.exists())
```

- **Example output:**

```text
BackendUnavailableError True False
```

MMseqs2 is not installed here, so the example shows the error for a program path that does not exist. The GenBank files are read first, so they must be valid.

Usage with the dependency installed:

```python
from dnakit.visualization import compare_genomes_mmseqs

result = compare_genomes_mmseqs(["genomeA.gbk", "genomeB.gbk"], "mmseqs_out", min_identity=30, dpi=150)
print(result["status"], result["data"]["alignment_count"], result["image_path"])
```

## 28) `dnakit.visualization.compare_genomes_mummer()`

- **Purpose:** Align neighbouring genomes with MUMmer and draw the alignment links in `result.png`.
- **Method:** Write each genome as a FASTA file and, for every pair of consecutive genomes, run `promer` (for `seqtype="protein"`) or `nucmer` (for `"nucleotide"`) with `--mum`, then `delta-filter -m` (`-1` with `one_to_one`), then `show-coords -H -T -r -k`. The coordinates are parsed by pyGenomeViz, shifted to genome coordinates, written to `align_coords.tsv`, and the alignments with at least `min_length` and `min_identity` percent are drawn as links. The temporary delta and coordinate files are removed.
- **API:** `compare_genomes_mummer(gbk_resources: 'Sequence[str]', output_dir: 'str', *, seqtype: 'str' = 'protein', min_length: 'int' = 0, min_identity: 'float' = 0.0, one_to_one: 'bool' = False, fig_width: 'float' = 15.0, dpi: 'int' = 300, timeout_seconds: 'float' = 3600.0, mummer_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package pygenomeviz (extra `dnakit[bio-pygenomeviz]`) and the MUMmer 4 programs `promer` or `nucmer`, `delta-filter` and `show-coords` on PATH, or in the folder `mummer_dir`.
- **Input:** `gbk_resources`, at least two GenBank files; `output_dir`, a new or empty folder; `seqtype`, `protein` or `nucleotide`; `min_length`, a non-negative integer; `min_identity`, a percentage from 0 to 100; `one_to_one`, a boolean; `fig_width`, `dpi`, `timeout_seconds`; `mummer_dir`, an optional folder with the programs.
- **Output:** A dictionary with `status` (`completed` or `failed`), `backend`, `commands`, `return_codes`, `outputs` (`align_coords.tsv` and `result.png`), `image_path`, `log_tail` and `data` (`genomes`, `alignment_count` and `alignments`).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.visualization import compare_genomes_mummer

out = Path(tempfile.mkdtemp()) / "mummer_out"
try:
    compare_genomes_mummer(
        ["genomeA.gbk", "genomeB.gbk"],
        str(out),
        seqtype="nucleotide",
        mummer_dir="/nonexistent/mummer",
        verbose=False,
    )
except BackendUnavailableError as error:
    print(type(error).__name__, "nucmer was not found" in str(error), out.exists())
```

- **Example output:**

```text
BackendUnavailableError True False
```

MUMmer is not installed here, so the example points `mummer_dir` at a folder that does not exist and shows the error. The output folder is not created when a program is missing.

Usage with the dependency installed:

```python
from dnakit.visualization import compare_genomes_mummer

result = compare_genomes_mummer(["genomeA.gbk", "genomeB.gbk"], "mummer_out", seqtype="nucleotide", min_identity=70)
print(result["status"], result["data"]["alignment_count"], result["image_path"])
```

## 29) `dnakit.visualization.compare_genomes_progressive_mauve()`

- **Purpose:** Align sequence files with progressiveMauve and draw the conserved blocks and their links in `result.png`.
- **Method:** Copy the input files into a `seqfiles` folder, run `progressiveMauve --output=mauve.xmfa --backbone-output=mauve_bbone.tsv <files>`, parse the backbone file into block coordinates (blocks shorter than 20 bp and blocks missing from at least one genome are dropped, and the blocks are ordered by the reference genome `refid`), write them to `align_coords.tsv`, and draw one track per genome with a coloured box for each block and links between neighbouring genomes.
- **API:** `compare_genomes_progressive_mauve(seq_files: 'Sequence[str]', output_dir: 'str', *, refid: 'int' = 0, cmap: 'str' = 'hsv', fig_width: 'float' = 15.0, dpi: 'int' = 300, timeout_seconds: 'float' = 3600.0, progressive_mauve_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package pygenomeviz (extra `dnakit[bio-pygenomeviz]`) and the external program progressiveMauve on PATH or `progressive_mauve_path`.
- **Input:** `seq_files`, at least two FASTA or GenBank files (genomes are named after the file name without its suffix; the names must be unique and use only letters, digits, `.`, `_` and `-`); `output_dir`, a new or empty folder; `refid`, the index of the reference file; `cmap`, a matplotlib colormap name; `fig_width`, `dpi`, `timeout_seconds`.
- **Output:** A dictionary with `status` (`completed` or `failed`), `backend`, `commands`, `return_codes`, `outputs` (`mauve.xmfa`, `mauve_bbone.tsv`, `align_coords.tsv` and `result.png`), `image_path`, `log_tail` and `data` (`genomes` with the sequence file of each, `alignment_count` and `alignments`).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.visualization import compare_genomes_progressive_mauve

out = Path(tempfile.mkdtemp()) / "mauve_out"
try:
    compare_genomes_progressive_mauve(
        ["genomeA.gbk", "genomeB.gbk"],
        str(out),
        progressive_mauve_path="/nonexistent/progressiveMauve",
        verbose=False,
    )
except BackendUnavailableError as error:
    print(type(error).__name__, "progressiveMauve was not found" in str(error), out.exists())
```

- **Example output:**

```text
BackendUnavailableError True False
```

progressiveMauve is not installed here, so the example shows the error for a program path that does not exist.

Usage with the dependency installed:

```python
from dnakit.visualization import compare_genomes_progressive_mauve

result = compare_genomes_progressive_mauve(["genomeA.gbk", "genomeB.gbk"], "mauve_out", refid=0)
print(result["status"], result["data"]["alignment_count"], result["image_path"])
```

## 30) `dnakit.visualization.load_example_dataset()`

- **Purpose:** Download a pyGenomeViz example GenBank dataset from GitHub and list its files and links.
- **Method:** Check the name against `pygenomeviz.config.DATASETS`, then call `pygenomeviz.utils.load_example_dataset(name, cache_dir=..., overwrite_cache=...)`, which downloads the GenBank files from the pygenomeviz-data repository into a cache folder (or reuses the cache) and returns their paths and the link definitions.
- **API:** `load_example_dataset(name: 'str', *, cache_dir: 'str | None' = None, overwrite_cache: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package pygenomeviz (extra `dnakit[bio-pygenomeviz]`) and internet access to GitHub (no key); a cached dataset can be reused offline.
- **Input:** `name`, one of the dataset names in pyGenomeViz (an unknown name raises `ConfigurationError` before any download); `cache_dir`, an optional folder; `overwrite_cache`, a boolean.
- **Output:** A dictionary with `status`, `service`, `source_url`, `dataset`, `request`, `retrieved_at_utc`, `gbk_files` (paths), `link_count` and `links`.
- **Example code:**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.visualization import load_example_dataset

try:
    load_example_dataset("no_such_dataset", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, "Unknown pyGenomeViz example dataset" in str(error))
```

- **Example output:**

```text
ConfigurationError True
```

The example does not touch the network. It shows that an unknown dataset name is rejected first.

Usage with the dependency installed:

```python
from dnakit.visualization import load_example_dataset

result = load_example_dataset("escherichia_coli", cache_dir="pgv_cache")
print(len(result["gbk_files"]), result["link_count"])
```

## 31) `dnakit.visualization.load_example_gff()`

- **Purpose:** Download a pyGenomeViz example GFF file from GitHub into a cache folder and return its path.
- **Method:** Check the file name against `pygenomeviz.config.EXAMPLE_GFF_FILES`, then call `pygenomeviz.utils.load_example_gff(filename, cache_dir=..., overwrite_cache=...)`.
- **API:** `load_example_gff(filename: 'str', *, cache_dir: 'str | None' = None, overwrite_cache: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package pygenomeviz (extra `dnakit[bio-pygenomeviz]`) and internet access to GitHub (no key); a cached file can be reused offline.
- **Input:** `filename`, one of the example GFF names in pyGenomeViz (an unknown name raises `ConfigurationError` before any download); `cache_dir`, an optional folder; `overwrite_cache`, a boolean.
- **Output:** A dictionary with `status`, `service`, `source_url`, `filename`, `request`, `retrieved_at_utc` and `gff_path` (the absolute path of the cached file).
- **Example code:**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.visualization import load_example_gff

try:
    load_example_gff("no_such_file.gff", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, "Unknown pyGenomeViz example GFF" in str(error))
```

- **Example output:**

```text
ConfigurationError True
```

The example does not touch the network. It shows that an unknown file name is rejected first.

Usage with the dependency installed:

```python
from dnakit.visualization import load_example_gff

result = load_example_gff("enterobacteria_phage.gff", cache_dir="pgv_cache")
print(result["gff_path"])
```

## 32) `dnakit.visualization.plot_genome_features()`

- **Purpose:** Draw the CDS features of GenBank genomes as stacked pyGenomeViz tracks in one PNG.
- **Method:** Read the GenBank files with pyGenomeViz, add one feature track per genome, and draw every CDS that has a `translation` qualifier and lies inside the genome range as an arrow (`bigarrow` or `arrow`) pointing along its strand. The strand is read from `feature.location.strand`, because upstream code reads `feature.strand`, which Biopython 1.88 no longer defines. The figure is saved with `bbox_inches="tight"`.
- **API:** `plot_genome_features(gbk_resources: 'Sequence[str]', image_path: 'str', *, fig_width: 'float' = 15.0, fig_track_height: 'float' = 1.0, feature_color: 'str' = 'orange', feature_plotstyle: 'str' = 'bigarrow', align_type: 'str' = 'left', tick_style: 'str | None' = None, dpi: 'int' = 300, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package pygenomeviz (extra `dnakit[bio-pygenomeviz]`), which brings matplotlib and Biopython. No external program is needed.
- **Input:** `gbk_resources`, at least one GenBank file (`file:start-end` selects a range); `image_path`, a new path that ends in `.png` in an existing folder; `fig_width` and `fig_track_height`, positive numbers; `feature_color`, a matplotlib colour; `feature_plotstyle`, `bigarrow` or `arrow`; `align_type`, `left`, `center` or `right`; `tick_style`, None, `bar` or `axis`; `dpi`, a positive integer.
- **Output:** A dictionary with `status`, `image_path` and `data`, where `data["genomes"]` has, for each genome, `name`, `min_range`, `max_range`, `range_size`, `reverse`, `full_genome_length` and `cds_count`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.visualization import plot_genome_features

image = Path(tempfile.mkdtemp()) / "features.png"
result = plot_genome_features(
    ["genomeA.gbk", "genomeB.gbk"], str(image), dpi=80, verbose=False
)
for genome in result["data"]["genomes"]:
    print(genome)
print(image.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n")
```

- **Example output:**

```text
{'name': 'genomeA', 'min_range': 0, 'max_range': 84, 'range_size': 84, 'reverse': False, 'full_genome_length': 84, 'cds_count': 2}
{'name': 'genomeB', 'min_range': 0, 'max_range': 86, 'range_size': 86, 'reverse': False, 'full_genome_length': 86, 'cds_count': 2}
True
```
