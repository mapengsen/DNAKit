# Long-read QC, assembly, structural variants, depth, methylation and FastTree facades

These functions run [FreeBayes](https://github.com/freebayes/freebayes) helper tools, [Filtlong](https://github.com/rrwick/Filtlong), [NanoFilt](https://github.com/wdecoster/nanofilt), [Porechop](https://github.com/rrwick/Porechop), [mosdepth](https://github.com/brentp/mosdepth), [Sniffles2](https://github.com/fritzsedlazeck/Sniffles), [cuteSV](https://github.com/tjiangHIT/cuteSV), [vcfanno](https://github.com/brentp/vcfanno), [Flye](https://github.com/mikolmogorov/Flye), [Liftoff](https://github.com/agshumate/Liftoff), [miniprot](https://github.com/lh3/miniprot), [MethylDackel](https://github.com/dpryan79/MethylDackel) and [FastTree](http://www.microbesonline.org/fasttree/) with fixed argument whitelists and no shell. The FreeBayes region/VCF scripts, the NanoFilt time filter are re-implemented in Python with the exact semantics of the original scripts (checked against them) and need no program.

The programs are not installed by DNAKit; put them on `PATH` or pass the `*_path` argument. Outputs go to new paths. Tested versions: mosdepth 0.3.14, FastTree 2.2.0, miniprot 0.18, Filtlong 0.3.1, NanoFilt 2.8.0, Sniffles2 2.6.2, cuteSV 2.1.4, vcfanno 0.3.9, Flye 2.9.6, Liftoff 1.6.3, MethylDackel 0.6.1, bamleftalign (FreeBayes 1.3.10). Porechop is archived and its conda build could not be run here; its command assembly is covered by a unit test with a recording stand-in.

## Example data

The declared input files are written before the examples run; the other examples create their data in code.

**Example file:** `ref.fa.fai`

```text
chr1	1000	6	60	61
chr2	450	1100	60	61
tiny	7	2000	60	61
```

**Example file:** `a.vcf`

```text
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	10	.	A	T	.	.	.
```

**Example file:** `b.vcf`

```text
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
2	20	.	G	C	.	.	.
```

**Example file:** `aln.fa`

```text
>t1
ACGTACGTAGCTAGCTAGGATCGATCGTAGCTAGCATCGATCGAT
>t2
ACGTACGTAGCTAGCTAGGATCGATCGTAGCTAGCATCGATCGAT
>t3
ACGTACGTAGCAAGCTAGGATCGATCGTAGCTAGCATCGTTCGAT
>t4
ACGTTCGTAGCAAGCTAGGATCGATCGAAGCTAGCATCGTTCGAT
>t5
TCGTTCGTAGCAAGCTAGGATCGTTCGAAGCTAGCATCGTTCGAC
```

**Example file:** `reads.fastq`

```text
@low runid=x start_time=2024-01-01T10:00:00Z
ATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATAT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@mid runid=x start_time=2024-01-02T10:00:00Z
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@high runid=x start_time=2024-01-03T10:00:00Z
GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
```

## 1) `dnakit.variants.freebayes_p16_cli_functions.keep_first_vcf_header()`

- **Purpose:** Concatenate VCF files keeping only the first header block (FreeBayes `vcffirstheader`).
- **Method:** The files are read in order as one stream (as `cat a.vcf b.vcf | vcffirstheader` would). Meta lines (`##`) and the `#CHROM` line are written until the first `#CHROM` line has been seen; every later header line is dropped. Data lines are copied with trailing white space removed. `.gz` inputs are read transparently.
- **API:** `dnakit.variants.freebayes_p16_cli_functions.keep_first_vcf_header(inputs, output_path, verbose=True)`.
- **Input:** Existing input files, a new output path (existing outputs are never overwritten) and the keyword options shown in the API line.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (files written), `log_tail` and `provenance`; also `header_lines_written`, `header_lines_dropped`, `records_written`.
- **See the example under:** `partition_reference_into_regions()`

## 2) `dnakit.variants.freebayes_p16_cli_functions.left_align_bam_indels()`

- **Purpose:** Left-align and merge indels in BAM alignments (FreeBayes `bamleftalign`).
- **Method:** `reference` is the FASTA the reads were mapped to (required by the program, `--fasta-reference`). The BAM is streamed to the program's standard input and the realigned BAM is written to `output_path` (a new file). `compressed=True` writes a compressed BAM (`--compressed`; the program's own default is uncompressed). `max_iterations` limits the left-realignment passes (`--max-iterations`).
- **API:** `dnakit.variants.freebayes_p16_cli_functions.left_align_bam_indels(alignment, output_path, *, reference, compressed=True, max_iterations=None, bamleftalign_path=None, verbose=True)`.
- **Input:** Existing input files, a new output path (existing outputs are never overwritten) and the keyword options shown in the API line.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (files written), `log_tail` and `provenance`; also `alignment`, `output`.
- **See the example under:** `compute_mosdepth_depth()`

## 3) `dnakit.variants.freebayes_p16_cli_functions.partition_reference_by_bam_size()`

- **Purpose:** Split a reference into regions holding equal amounts of BAM data (FreeBayes `split_ref_by_bai_datasize.py`).
- **Method:** The data amount per 16 kb interval is read from the linear index of every `<bam>.bai` (the difference of the BGZF virtual file offsets, used as a size exactly as the original script does) and summed over the BAM files; each sequence is cut where the cumulated amount reaches multiples of `target_data_size` (the script's `-s`, default `100e6`, multiplied by 1,000,000 as in the script; values below 1 are refused). Sequences with fewer than two index intervals are returned whole. Needs numpy and scipy (`dnakit[bio-analysis]`). Rows are `(chrom, start, end)`; `output_path` receives them tab-separated as the original script prints them, and `regions` in the result is cut at 100,000 rows.
- **API:** `dnakit.variants.freebayes_p16_cli_functions.partition_reference_by_bam_size(reference, bam_files, *, output_path=None, target_data_size=100000000.0, bai_interval_size=16384, verbose=True)`.
- **Input:** Existing input files, a new output path (existing outputs are never overwritten) and the keyword options shown in the API line.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (files written), `log_tail` and `provenance`; also `region_count`, `regions`.

## 4) `dnakit.variants.freebayes_p16_cli_functions.partition_reference_by_coverage()`

- **Purpose:** Make regions of about equal total coverage from per-position depth (FreeBayes `coverage_to_regions.py`).
- **Method:** `coverage` is the three-column `chrom<TAB>pos<TAB>depth` output of `bamtools coverage` (one line per covered position, sorted); `reference` is the FASTA or its `.fai`. A new region is closed when the running depth sum exceeds `total_depth / num_regions`; the last region of each sequence ends at its length. Returns `chrom:start-end` regions as the original script prints them.
- **API:** `dnakit.variants.freebayes_p16_cli_functions.partition_reference_by_coverage(reference, coverage, num_regions, *, output_path=None, verbose=True)`.
- **Input:** Existing input files, a new output path (existing outputs are never overwritten) and the keyword options shown in the API line.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (files written), `log_tail` and `provenance`; also `region_count`, `regions`, `total_depth`.

## 5) `dnakit.variants.freebayes_p16_cli_functions.partition_reference_into_regions()`

- **Purpose:** Split a reference into fixed-size `chrom:start-end` regions (FreeBayes `fasta_generate_regions.py`).
- **Method:** `reference` is a FASTA path or its `.fai` index; the index must exist. With `chunks=False` each region is `region_size` bases long; with `chunks=True` every sequence is divided into `region_size` equal chunks. `chromosomes` limits the sequences. Regions are 0-based half-open, exactly as the original script prints them. `output_path` (a new file) receives one region per line; `bed_prefix` instead writes one BED file per region named `<bed_prefix>.<chrom>.region.<n>.bed` (the Snakemake mode of the script). Returns the (bounded) region list.
- **API:** `dnakit.variants.freebayes_p16_cli_functions.partition_reference_into_regions(reference, region_size, *, output_path=None, chunks=False, chromosomes=None, bed_prefix=None, verbose=True)`.
- **Input:** Existing input files, a new output path (existing outputs are never overwritten) and the keyword options shown in the API line.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (files written), `log_tail` and `provenance`; also `region_count`, `regions`.
- **Example code:**

```python
from dnakit.variants.freebayes_p16_cli_functions import keep_first_vcf_header, partition_reference_into_regions

print(partition_reference_into_regions("ref.fa", 300, verbose=False)["regions"])
print(partition_reference_into_regions("ref.fa", 3, chunks=True, chromosomes=["chr1"], verbose=False)["regions"])
result = keep_first_vcf_header(["a.vcf", "b.vcf"], "merged.vcf", verbose=False)
print(result["header_lines_dropped"], result["records_written"])
```

- **Example output:**

```text
['chr1:0-300', 'chr1:300-600', 'chr1:600-900', 'chr1:900-1000', 'chr2:0-300', 'chr2:300-450', 'tiny:0-7']
['chr1:0-334', 'chr1:334-668', 'chr1:668-1000']
2 2
```

## 6) `dnakit.sequencing_qc.p16_readqc_cli_functions.demultiplex_nanopore_barcodes()`

- **Purpose:** Bin reads by native/rapid barcode while trimming adapters (`porechop --barcode_dir`).
- **Method:** One FASTQ per barcode (`BC01.fastq` ...) and `none.fastq` for unassigned reads is written to the new `output_directory`. `require_two_barcodes` (`--require_two_barcodes`) demands a barcode at both ends; `barcode_threshold` is `--barcode_threshold` (percent identity, default 75); `discard_unassigned` is `--discard_unassigned`.
- **API:** `dnakit.sequencing_qc.p16_readqc_cli_functions.demultiplex_nanopore_barcodes(reads, output_directory, *, require_two_barcodes=False, barcode_threshold=None, discard_unassigned=False, threads=1, porechop_path=None, verbose=True)`.
- **Input:** Existing input files, a new output path (existing outputs are never overwritten) and the keyword options shown in the API line.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (files written), `log_tail` and `provenance`; also `reads_per_bin`.

## 7) `dnakit.sequencing_qc.p16_readqc_cli_functions.filter_long_reads_with_reference()`

- **Purpose:** Score long reads by k-mer matches to a reference, then filter, trim or split them (`filtlong`).
- **Method:** Exactly one kind of reference may be given: Illumina reads (`illumina_reads_1` and optionally `illumina_reads_2`, `--short_1`/`--short_2`) or an `assembly` (`--assembly`); read quality is then judged by k-mer identity to it instead of the Phred scores. `trim` removes non-matching bases from both read ends (`--trim`); `split` cuts reads at runs of at least that many non-matching bases (`--split`). Both need a reference. Thresholds are the Filtlong ones: `min_length`, `max_length`, `min_mean_quality`, `target_bases`, `keep_percent`. The kept reads are written to `output_path` (FASTQ).
- **API:** `dnakit.sequencing_qc.p16_readqc_cli_functions.filter_long_reads_with_reference(reads, output_path, *, illumina_reads_1=None, illumina_reads_2=None, assembly=None, trim=False, split=None, min_length=None, max_length=None, min_mean_quality=None, target_bases=None, keep_percent=None, filtlong_path=None, verbose=True)`.
- **Input:** Existing input files, a new output path (existing outputs are never overwritten) and the keyword options shown in the API line.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (files written), `log_tail` and `provenance`; also `reads_kept`.

## 8) `dnakit.sequencing_qc.p16_readqc_cli_functions.filter_nanopore_reads_by_gc_summary()`

- **Purpose:** Filter long reads by GC content or by the quality in a sequencing-summary file (`NanoFilt`).
- **Method:** `min_gc`/`max_gc` are fractions in 0-1 (`--minGC`/`--maxGC`; NanoFilt ignores them when a summary is used). `sequencing_summary` (`--summary`) takes read quality from the ONT `sequencing_summary.txt` instead of the FASTQ; `read_type` is 1D, 2D or 1D2. `min_quality`, `min_length` and `max_length` are the NanoFilt thresholds. The uncompressed FASTQ is read from `reads` and the kept reads are written to `output_path`.
- **API:** `dnakit.sequencing_qc.p16_readqc_cli_functions.filter_nanopore_reads_by_gc_summary(reads, output_path, *, min_gc=None, max_gc=None, sequencing_summary=None, min_quality=None, min_length=None, max_length=None, read_type=None, nanofilt_path=None, verbose=True)`.
- **Input:** Existing input files, a new output path (existing outputs are never overwritten) and the keyword options shown in the API line.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (files written), `log_tail` and `provenance`; also `reads_kept`.
- **Example code:**

```python
import gzip
from dnakit.sequencing_qc.p16_readqc_cli_functions import filter_nanopore_reads_by_gc_summary, filter_nanopore_reads_by_start_time

result = filter_nanopore_reads_by_gc_summary("reads.fastq", "gc.fastq", min_gc=0.3, max_gc=0.7, verbose=False)
print(result["reads_kept"], open("gc.fastq").readline().split()[0])
with gzip.open("reads.fastq.gz", "wt") as handle:
    handle.write(open("reads.fastq").read())
timed = filter_nanopore_reads_by_start_time("reads.fastq.gz", "timed.fastq", time_from="2024-01-02T00:00:00Z", time_to="2024-01-03T00:00:00Z", verbose=False)
print(timed["reads_total"], timed["reads_kept"])
```

- **Example output:**

```text
1 @mid
3 1
```

## 9) `dnakit.sequencing_qc.p16_readqc_cli_functions.filter_nanopore_reads_by_start_time()`

- **Purpose:** Keep reads whose sequencing start time lies in a window (NanoFilt `scripts/timefilt.py`).
- **Method:** `reads_gz` is a gzip-compressed FASTQ whose headers carry `start_time=<ISO time>` (as written by MinKNOW/Guppy); `time_from`/`time_to` are inclusive bounds in any format `dateutil` parses; give a time zone (for example `2017-10-13T00:00:00Z`) because read times carry one and naive and aware times cannot be compared. Reads are written unchanged (the original rewrites them through Biopython, which gives the same four lines) to `output_path`; reads without `start_time` raise an error, as in the script. Needs python-dateutil (extra `nanopore-time`).
- **API:** `dnakit.sequencing_qc.p16_readqc_cli_functions.filter_nanopore_reads_by_start_time(reads_gz, output_path, *, time_from, time_to, verbose=True)`.
- **Input:** Existing input files, a new output path (existing outputs are never overwritten) and the keyword options shown in the API line.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (files written), `log_tail` and `provenance`; also `reads_total`, `reads_kept`.
- **See the example under:** `filter_nanopore_reads_by_gc_summary()`

## 10) `dnakit.sequencing_qc.p16_readqc_cli_functions.trim_nanopore_adapters_discard_middle()`

- **Purpose:** Trim adapters and discard reads with an adapter in the middle (`porechop --discard_middle`).
- **Method:** Reads with an internal adapter are chimeras; Porechop's default splits them, `discard_middle=True` (this function's default) throws them away instead. `False` reproduces Porechop's default behaviour (split).
- **API:** `dnakit.sequencing_qc.p16_readqc_cli_functions.trim_nanopore_adapters_discard_middle(reads, output_path, *, discard_middle=True, threads=1, porechop_path=None, verbose=True)`.
- **Input:** Existing input files, a new output path (existing outputs are never overwritten) and the keyword options shown in the API line.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (files written), `log_tail` and `provenance`; also `discard_middle`.

## 11) `dnakit.variants.p16_variant_cli_functions.annotate_vcf_from_bam_vcfanno()`

- **Purpose:** Annotate VCF variants with read depth, mapping quality and read bases from a BAM (`vcfanno`).
- **Method:** A vcfanno configuration is written next to the output (`<output_path>.toml`) with one `[[annotation]]` block for the BAM and the chosen `operations` (vcfanno BAM fields): `depth` (read count, INFO `<bam>_depth`), `mapq` (mean mapping quality, `_mapq`) and `seq` (the read bases covering each variant, concatenated, `_seq`). The BAM must be indexed. Existing annotate_vcf_vcfanno takes a hand-written configuration; this function builds it for the BAM case.
- **API:** `dnakit.variants.p16_variant_cli_functions.annotate_vcf_from_bam_vcfanno(vcf, alignment, output_path, *, operations=('depth', 'mapq', 'seq'), threads=1, vcfanno_path=None, verbose=True)`.
- **Input:** Existing input files, a new output path (existing outputs are never overwritten) and the keyword options shown in the API line.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (files written), `log_tail` and `provenance`; also `configuration`, `records`.

## 12) `dnakit.variants.p16_variant_cli_functions.call_structural_variants_sniffles2()`

- **Purpose:** Call SVs from an aligned long-read BAM/CRAM with Sniffles2, optionally mosaic, phased, and saving an SNF.
- **Method:** `mosaic=True` (`--mosaic`) detects low-allele-fraction/somatic SVs; `phase=True` (`--phase`) phases calls with the reads' haplotype tags (HP/PS; needs a phased BAM); `snf_path` (`--snf`) saves the per-sample candidates in a new SNF file for later multi-sample calling (:func:`merge_structural_variants_sniffles2_snf`). `min_support` is `--minsupport` (integer or `auto`). The VCF is written to `vcf_path`.
- **API:** `dnakit.variants.p16_variant_cli_functions.call_structural_variants_sniffles2(alignment, vcf_path, *, reference=None, snf_path=None, mosaic=False, phase=False, min_sv_length=None, min_support=None, sample_id=None, threads=1, sniffles_path=None, verbose=True)`.
- **Input:** Existing input files, a new output path (existing outputs are never overwritten) and the keyword options shown in the API line.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (files written), `log_tail` and `provenance`; also `mosaic`, `phase`, `records`.

## 13) `dnakit.variants.p16_variant_cli_functions.compute_mosdepth_depth()`

- **Purpose:** Per-base depth, genome-wide depth distribution and per-chromosome summary of a BAM/CRAM (`mosdepth`).
- **Method:** Writes `<prefix>.per-base.bed.gz` (unless `per_base=False`, `--no-per-base`), `<prefix>.mosdepth.global.dist.txt` (cumulative fraction of bases at each depth) and `<prefix>.mosdepth.summary.txt` (length, bases, mean, min, max per chromosome and total, parsed into `summary`). The remaining parameters are mosdepth's: `--chrom`, `--mapq`, `--min-frag-len`, `--max-frag-len`, `--fast-mode`, `--fragment-mode`, `--flag` (default excludes 1796), `--include-flag`, `--fasta` (CRAM) and `--read-groups`.
- **API:** `dnakit.variants.p16_variant_cli_functions.compute_mosdepth_depth(alignment, prefix, *, per_base=True, chromosome=None, threads=0, mapq=0, min_fragment=None, max_fragment=None, fast_mode=False, fragment_mode=False, exclude_flag=None, include_flag=None, reference=None, read_groups=None, mosdepth_path=None, verbose=True)`.
- **Input:** Existing input files, a new output path (existing outputs are never overwritten) and the keyword options shown in the API line.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (files written), `log_tail` and `provenance`; also `summary`, `distribution_head`.
- **Example code:**

```python
import random
import pysam
from dnakit.variants.p16_variant_cli_functions import compute_mosdepth_depth, compute_mosdepth_regions

rng = random.Random(11)
genome = "".join(rng.choice("ACGT") for _ in range(6000))
open("genome.fa", "w").write(">chrA\n" + genome + "\n")
pysam.faidx("genome.fa")
header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chrA", "LN": 6000}]}
with pysam.AlignmentFile("reads.bam", "wb", header=header) as out:
    for n, start in enumerate(range(0, 5850, 25)):
        read = pysam.AlignedSegment()
        read.query_name, read.query_sequence, read.flag = f"r{n}", genome[start:start + 150], 0
        read.reference_id, read.reference_start, read.mapping_quality, read.cigar = 0, start, 60, [(0, 150)]
        read.query_qualities = pysam.qualitystring_to_array("I" * 150)
        out.write(read)
pysam.index("reads.bam")
depth = compute_mosdepth_depth("reads.bam", "depth", verbose=False)
print(depth["summary"][-1]["chrom"], depth["summary"][-1]["length"], depth["summary"][-1]["mean"])
regions = compute_mosdepth_regions("reads.bam", "win", by="2000", thresholds=[1, 5], verbose=False)
print([(r["start"], r["end"], r["depth"]) for r in regions["regions"]])
```

- **Example output:**

```text
total 6000 5.85
[(0, 2000, 5.81), (2000, 4000, 6.0), (4000, 6000, 5.74)]
```

## 14) `dnakit.variants.p16_variant_cli_functions.compute_mosdepth_quantized_coverage()`

- **Purpose:** Collapse per-base depth into labelled coverage segments (`mosdepth --quantize`).
- **Method:** `quantize` is mosdepth's bin specification, for example `"0:1:4:"` (bins 0, 1-3, 4+) or `"0:1:5:150"`, and writes `<prefix>.quantized.bed.gz` with `chrom, start, end, label` rows (parsed into `segments`, first 1,000).
- **API:** `dnakit.variants.p16_variant_cli_functions.compute_mosdepth_quantized_coverage(alignment, prefix, *, quantize, chromosome=None, threads=0, mapq=0, fast_mode=False, reference=None, mosdepth_path=None, verbose=True)`.
- **Input:** Existing input files, a new output path (existing outputs are never overwritten) and the keyword options shown in the API line.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (files written), `log_tail` and `provenance`; also `segments`.
- **See the example under:** `compute_mosdepth_depth()`

## 15) `dnakit.variants.p16_variant_cli_functions.compute_mosdepth_regions()`

- **Purpose:** Mean (or median) depth, depth histogram and threshold counts per window or BED region (`mosdepth --by`).
- **Method:** `by` is an integer window size (for example `"1000"`) or the path of a BED file. Outputs are `<prefix>.regions.bed.gz` (mean depth per region; median with `use_median`), `<prefix>.mosdepth.region.dist.txt` (cumulative depth distribution per chromosome over the regions) and, with `thresholds` (`--thresholds t1,t2`), `<prefix>.thresholds.bed.gz` with the number of bases covered at least that deeply. Per-base output is off by default (`--no-per-base`); `per_base=True` also writes it. Other parameters as in :func:`compute_mosdepth_depth`.
- **API:** `dnakit.variants.p16_variant_cli_functions.compute_mosdepth_regions(alignment, prefix, *, by, thresholds=None, use_median=False, per_base=False, chromosome=None, threads=0, mapq=0, fast_mode=False, fragment_mode=False, exclude_flag=None, include_flag=None, reference=None, read_groups=None, mosdepth_path=None, verbose=True)`.
- **Input:** Existing input files, a new output path (existing outputs are never overwritten) and the keyword options shown in the API line.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (files written), `log_tail` and `provenance`; also `summary`, `regions`.
- **See the example under:** `compute_mosdepth_depth()`

## 16) `dnakit.variants.p16_variant_cli_functions.genotype_structural_variants_cutesv()`

- **Purpose:** Call SVs from long reads and genotype them (`cuteSV --genotype`).
- **Method:** `genotype=True` adds genotypes and allele depths (`--genotype`). `force_calling_vcf` (`-Ivcf`) genotypes the SVs of that VCF instead of calling de novo. `platform` sets the preset: `ont` (`--max_cluster_bias_INS 100 ...`), `hifi` or `clr` (Pacific Biosciences) as the program's documented parameter sets. `work_directory` must be new.
- **API:** `dnakit.variants.p16_variant_cli_functions.genotype_structural_variants_cutesv(alignment, reference, vcf_path, work_directory, *, genotype=True, force_calling_vcf=None, platform=None, min_support=None, min_size=None, sample=None, threads=1, cutesv_path=None, verbose=True)`.
- **Input:** Existing input files, a new output path (existing outputs are never overwritten) and the keyword options shown in the API line.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (files written), `log_tail` and `provenance`; also `genotype`, `records`.

## 17) `dnakit.variants.p16_variant_cli_functions.merge_structural_variants_sniffles2_snf()`

- **Purpose:** Multi-sample SV calling from per-sample SNF candidate files (`sniffles --input a.snf b.snf ...`).
- **Method:** 
- **API:** `dnakit.variants.p16_variant_cli_functions.merge_structural_variants_sniffles2_snf(snf_files, vcf_path, *, reference=None, min_support=None, threads=1, sniffles_path=None, verbose=True)`.
- **Input:** Existing input files, a new output path (existing outputs are never overwritten) and the keyword options shown in the API line.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (files written), `log_tail` and `provenance`; also `samples`, `records`.

## 18) `dnakit.annotation.p16_assembly_annotation_cli_functions.align_proteins_to_genome_miniprot()`

- **Purpose:** Align protein sequences to a genome, spliced and frameshift-aware (`miniprot`).
- **Method:** `genome_or_index` is a genome FASTA or an index from :func:`build_miniprot_index`. `output_format` is `paf` (default), `gff` (GFF3, `--gff`) or `gtf` (`--gtf`). `translate` adds the translated protein of each alignment (`--trans`, skipping frameshifts) and `include_alignment` the residue alignment (`--aln`); both apply to the PAF/GFF3 output. `no_splicing` is `-S`; `max_intron_size` is `-G`; `splice_model` is `-j` (2 vertebrate/insect, 1 general, 0 none). The result file is written to `output_path`.
- **API:** `dnakit.annotation.p16_assembly_annotation_cli_functions.align_proteins_to_genome_miniprot(genome_or_index, proteins, output_path, *, output_format='paf', translate=False, include_alignment=False, no_splicing=False, max_intron_size=None, splice_model=None, translation_table=1, threads=1, miniprot_path=None, verbose=True)`.
- **Input:** Existing input files, a new output path (existing outputs are never overwritten) and the keyword options shown in the API line.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (files written), `log_tail` and `provenance`; also `output_format`, `lines`, `head`.
- **See the example under:** `build_miniprot_index()`

## 19) `dnakit.annotation.p16_assembly_annotation_cli_functions.assemble_long_reads_flye()`

- **Purpose:** De novo assembly of long reads with Flye, optionally in metagenome mode and with scaffolding.
- **Method:** `read_type` is one of pacbio-raw, pacbio-corr, pacbio-hifi, nano-raw, nano-corr, nano-hq. `metagenome=True` (`--meta`) is metaFlye for uneven coverage; `scaffold=True` (`--scaffold`) joins contigs into gapped scaffolds along the repeat graph. `genome_size` (`--genome-size`, for example `"5m"`) is optional in recent Flye. `iterations` is the number of polishing rounds (`--iterations`). `output_directory` must be new; the assembly is `assembly.fasta` with `assembly_info.txt` (parsed into `contigs`).
- **API:** `dnakit.annotation.p16_assembly_annotation_cli_functions.assemble_long_reads_flye(reads, output_directory, *, read_type='nano-raw', genome_size=None, metagenome=False, scaffold=False, iterations=1, min_overlap=None, threads=1, flye_path=None, verbose=True)`.
- **Input:** Existing input files, a new output path (existing outputs are never overwritten) and the keyword options shown in the API line.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (files written), `log_tail` and `provenance`; also `contigs` (from assembly_info.txt).

## 20) `dnakit.annotation.p16_assembly_annotation_cli_functions.build_miniprot_index()`

- **Purpose:** Build a reusable miniprot genome index (`miniprot -d FILE ref.fna`).
- **Method:** 
- **API:** `dnakit.annotation.p16_assembly_annotation_cli_functions.build_miniprot_index(genome, index_path, *, kmer=6, translation_table=1, threads=1, miniprot_path=None, verbose=True)`.
- **Input:** Existing input files, a new output path (existing outputs are never overwritten) and the keyword options shown in the API line.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (files written), `log_tail` and `provenance`; also `index`.
- **Example code:**

```python
from Bio.Seq import Seq
from dnakit.annotation.p16_assembly_annotation_cli_functions import align_proteins_to_genome_miniprot, build_miniprot_index

genome = open("genome.fa").read().split("\n")[1]
open("protein.faa", "w").write(">prot1\n" + str(Seq(genome[3000:3600]).translate()).replace("*", "A") + "\n")
index = build_miniprot_index("genome.fa", "genome.mpi", verbose=False)
paf = align_proteins_to_genome_miniprot(index["index"], "protein.faa", "aln.paf", verbose=False)
fields = paf["head"][0].split("\t")
print(fields[0], fields[5], fields[7], fields[8])
```

- **Example output:**

```text
prot1 chrA 3000 3600
```

## 21) `dnakit.annotation.p16_assembly_annotation_cli_functions.lift_over_annotation_liftoff()`

- **Purpose:** Lift a GFF/GTF gene annotation from a reference to a target assembly (`liftoff`).
- **Method:** Genes are aligned with minimap2 and mapped exon by exon. `extra_copies=True` (`-copies`) searches for additional gene copies with sequence identity at least `copy_identity` (`-sc`); `polish=True` (`-polish`) re-aligns exons to restore CDS start and stop codons. `min_coverage` (`-a`) and `min_identity` (`-s`) are Liftoff's thresholds; `feature_types` is a feature-type list file (`-f`). Unmapped gene ids go to `unmapped_path` (default `<output_gff>.unmapped.txt`).
- **API:** `dnakit.annotation.p16_assembly_annotation_cli_functions.lift_over_annotation_liftoff(target_fasta, reference_fasta, reference_annotation, output_gff, *, unmapped_path=None, extra_copies=False, polish=False, min_coverage=0.5, min_identity=0.5, copy_identity=0.95, feature_types=None, threads=1, liftoff_path=None, verbose=True)`.
- **Input:** Existing input files, a new output path (existing outputs are never overwritten) and the keyword options shown in the API line.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (files written), `log_tail` and `provenance`; also `lifted_genes`, `unmapped_features`.

## 22) `dnakit.annotation.p16_assembly_annotation_cli_functions.polish_assembly_flye()`

- **Purpose:** Polish a given draft assembly with long reads (standalone Flye polisher, `flye --polish-target`).
- **Method:** `reads` may be FASTA/FASTQ or an aligned BAM; the polished contigs are `polished_<iterations>.fasta`.
- **API:** `dnakit.annotation.p16_assembly_annotation_cli_functions.polish_assembly_flye(draft, reads, output_directory, *, read_type='nano-raw', iterations=1, threads=1, flye_path=None, verbose=True)`.
- **Input:** Existing input files, a new output path (existing outputs are never overwritten) and the keyword options shown in the API line.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (files written), `log_tail` and `provenance`; also `polished`.

## 23) `dnakit.epigenome.p16_methyldackel_cli_functions.convert_mappability_bigwig_to_bbm_methyldackel()`

- **Purpose:** Convert a mappability bigWig to a Binary BisMap (`.bbm`) file (`MethylDackel extract -M ... -N`).
- **Method:** MethylDackel writes the BBM as a by-product of an `extract` run that filters reads by the bigWig; this function runs that command with `-O` (the CpG bedGraph goes to `output_prefix`) and returns the BBM path. `mappability_threshold` is `-t` and `min_mappable_bases` is `-b`. MethylDackel 0.6.1 does not accept `-N <file>`, so the file is produced next to a symlink of the bigWig in a temporary folder and moved to `bbm_path`.
- **API:** `dnakit.epigenome.p16_methyldackel_cli_functions.convert_mappability_bigwig_to_bbm_methyldackel(reference, alignment, bigwig, bbm_path, output_prefix, *, mappability_threshold=0.01, min_mappable_bases=15, region=None, methyldackel_path=None, verbose=True)`.
- **Input:** Existing input files, a new output path (existing outputs are never overwritten) and the keyword options shown in the API line.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (files written), `log_tail` and `provenance`; also `bbm`, `bbm_bytes`.

## 24) `dnakit.epigenome.p16_methyldackel_cli_functions.extract_methylation_report_methyldackel()`

- **Purpose:** Per-cytosine methylation report, or per-CpG/CHG merged metrics, from a bisulfite BAM (`MethylDackel extract`).
- **Method:** `report` selects the output: `cytosine_report` (`--cytosine_report`) writes the exhaustive per-cytosine table with chromosome, 1-based position, strand, methylated count, unmethylated count, CG/CHG/CHH and trinucleotide context for every cytosine including uncovered ones, to `<prefix>.cytosine_report.txt`; `merge_context` (`--mergeContext`) merges the two strands of each CpG/CHG into `<prefix>_CpG.bedGraph` (and `_CHG`/`_CHH`); `fraction` (`--fraction`) writes `<prefix>_CpG.meth.bedGraph`; `per_cytosine` is the plain extract (`<prefix>_CpG.bedGraph` with one line per cytosine, the input of :func:`merge_cytosine_context_methyldackel`). `include_chg`/`include_chh` add contexts and `include_cpg=False` is `--noCpG`. `inclusion_bounds_ot`/`_ob` are `--OT`/`--OB` bounds from :func:`methylation_bias_methyldackel`. The first 20 lines of the main result are returned in `head`.
- **API:** `dnakit.epigenome.p16_methyldackel_cli_functions.extract_methylation_report_methyldackel(reference, alignment, output_prefix, *, report='cytosine_report', include_chg=False, include_chh=False, include_cpg=True, min_depth=1, min_mapq=10, min_base_quality=5, region=None, regions_bed=None, inclusion_bounds_ot=None, inclusion_bounds_ob=None, threads=1, methyldackel_path=None, verbose=True)`.
- **Input:** Existing input files, a new output path (existing outputs are never overwritten) and the keyword options shown in the API line.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (files written), `log_tail` and `provenance`; also `report`, `files`, `head`.

## 25) `dnakit.epigenome.p16_methyldackel_cli_functions.merge_cytosine_context_methyldackel()`

- **Purpose:** Merge per-cytosine metrics from `MethylDackel extract` into per-CpG/CHG metrics (`MethylDackel mergeContext`).
- **Method:** The input must be coordinate sorted (for example `<prefix>_CpG.bedGraph` of a plain extract run).
- **API:** `dnakit.epigenome.p16_methyldackel_cli_functions.merge_cytosine_context_methyldackel(reference, per_cytosine_metrics, output_path, *, methyldackel_path=None, verbose=True)`.
- **Input:** Existing input files, a new output path (existing outputs are never overwritten) and the keyword options shown in the API line.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (files written), `log_tail` and `provenance`; also `sites`.

## 26) `dnakit.epigenome.p16_methyldackel_cli_functions.methylation_bias_methyldackel()`

- **Purpose:** Methylation bias by read position: SVG plots and suggested trim bounds (`MethylDackel mbias`).
- **Method:** SVG files `<prefix>_<strand>_<context>.svg` are written; the suggested inclusion bounds the program prints (for use as `--OT`, `--OB`, `--CTOT`, `--CTOB` in :func:`extract_methylation_report_methyldackel`) are parsed into `suggested_bounds`.
- **API:** `dnakit.epigenome.p16_methyldackel_cli_functions.methylation_bias_methyldackel(reference, alignment, output_prefix, *, include_chg=False, include_chh=False, min_mapq=10, min_base_quality=5, region=None, regions_bed=None, max_depth=2000, threads=1, methyldackel_path=None, verbose=True)`.
- **Input:** Existing input files, a new output path (existing outputs are never overwritten) and the keyword options shown in the API line.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (files written), `log_tail` and `provenance`; also `suggested_bounds`, SVG files.

## 27) `dnakit.epigenome.p16_methyldackel_cli_functions.summarize_per_read_methylation_methyldackel()`

- **Purpose:** Average CpG methylation per read (`MethylDackel perRead`): read, chromosome, position, methylation %, informative bases.
- **Method:** 
- **API:** `dnakit.epigenome.p16_methyldackel_cli_functions.summarize_per_read_methylation_methyldackel(reference, alignment, output_path, *, min_mapq=10, min_base_quality=5, region=None, regions_bed=None, threads=1, methyldackel_path=None, verbose=True)`.
- **Input:** Existing input files, a new output path (existing outputs are never overwritten) and the keyword options shown in the API line.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (files written), `log_tail` and `provenance`; also `reads`, `rows`.

## 28) `dnakit.phylogeny.p16_fasttree_cli_functions.compute_fasttree_distance_matrix()`

- **Purpose:** Pairwise distance matrix of an alignment in PHYLIP format (`FastTree -makematrix`).
- **Method:** `nucleotide=True` uses `-nt` (Jukes-Cantor distances; protein alignments use log-corrected BLOSUM45 distances). `raw_distance` is `-rawdist` (fraction different) and `pseudocounts` is `-pseudo <weight>` (for sequences with little overlap). The matrix is written to `output_path` and returned as `labels` and `matrix`.
- **API:** `dnakit.phylogeny.p16_fasttree_cli_functions.compute_fasttree_distance_matrix(alignment, output_path, *, nucleotide=True, raw_distance=False, pseudocounts=None, fasttree_path=None, verbose=True)`.
- **Input:** Existing input files, a new output path (existing outputs are never overwritten) and the keyword options shown in the API line.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (files written), `log_tail` and `provenance`; also `labels`, `matrix`.
- **Example code:**

```python
from dnakit.phylogeny.p16_fasttree_cli_functions import compute_fasttree_distance_matrix, infer_fasttree_tree

matrix = compute_fasttree_distance_matrix("aln.fa", "dist.phy", verbose=False)
print(matrix["labels"], [round(x, 4) for x in matrix["matrix"][0]])
tree = infer_fasttree_tree("aln.fa", "tree.nwk", refine="noml", bootstrap_replicates=None, verbose=False)
print(tree["tips"], tree["newick"].count("("))
```

- **Example output:**

```text
['t1', 't2', 't3', 't4', 't5'] [0.0, 0.0, 0.0458, 0.0946, 0.1743]
5 3
```

## 29) `dnakit.phylogeny.p16_fasttree_cli_functions.infer_fasttree_tree()`

- **Purpose:** Infer or refine a tree with FastTree modes beyond a plain run.
- **Method:** `starting_tree` is a Newick file used as the starting topology (`-intree`; branch lengths are ignored). `refine` selects what happens next: `None` (default FastTree search), `"mllen"` (`-mllen`: optimise maximum-likelihood branch lengths without topology moves; with a starting tree this is the ML refinement of that tree), `"mllen_fixed_topology"` (`-mllen -nome`: branch lengths on a fixed topology) or `"noml"` (`-noml`: a minimum-evolution tree without maximum-likelihood steps). `constraints` is an alignment of 0/1 split constraints (`-constraints`) with `constraint_weight` (`-constraintWeight`, default 100). `bootstrap_replicates` (`-boot N`) adds resampled support values; `no_support` (`-nosupport`) omits support values. The Newick tree is written to `output_path` (`-out`).
- **API:** `dnakit.phylogeny.p16_fasttree_cli_functions.infer_fasttree_tree(alignment, output_path, *, nucleotide=True, gtr=False, gamma=False, starting_tree=None, refine=None, constraints=None, constraint_weight=None, bootstrap_replicates=None, no_support=False, seed=None, fasttree_path=None, verbose=True)`.
- **Input:** Existing input files, a new output path (existing outputs are never overwritten) and the keyword options shown in the API line.
- **Output:** A dictionary with `status`, `backend`, `version`, `command`, `return_code`, `outputs` (files written), `log_tail` and `provenance`; also `refine`, `newick`, `tips`.
- **See the example under:** `compute_fasttree_distance_matrix()`
