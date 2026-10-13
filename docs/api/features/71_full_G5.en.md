# Variant and interval helpers from pysam, cyvcf2 and vcfpy

These functions in `dnakit.variants` are adapters over optional packages that run offline on VCF files, genomic intervals, reference sequences and genotype arrays. They cover [pysam](https://github.com/pysam-developers/pysam), [cyvcf2](https://github.com/brentp/cyvcf2) and [vcfpy](https://github.com/dahlo/vcfpy). pysam's bcftools functions use the htslib and bcftools build that pysam ships, so no separate bcftools or tabix install is needed. Every function works on explicit files or sequences. None of them uses the network, calls a variant caller or queries a database. Dependencies are optional and are imported only when the matching function is called; a missing package raises an error with the install command.

## Example data

Run the examples in order: later examples use files that earlier ones create. The examples import from `dnakit.variants`. The `sample.vcf` below has 4 sites and 3 samples, with AD and DP fields; `ref.fa` is a 50 bp reference; `annotation.tsv` is a two-row annotation table. The files only demonstrate the interface, and their statistics do not describe a real population.

**Example file:** `sample.vcf`

```text
##fileformat=VCFv4.2
##contig=<ID=chr1,length=50>
##FILTER=<ID=q10,Description="Quality below 10">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="AF">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="DP">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3
chr1	3	rs1	G	T	50	PASS	DP=20;AF=0.5	GT:AD:DP	0/1:10,10:20	1/1:0,12:12	0/0:15,0:15
chr1	10	.	C	CA	40	PASS	DP=10	GT:AD:DP	0/1:5,5:10	0/0:10,0:10	0/1:4,6:10
chr1	22	rs3	C	G,A	60	q10	DP=30;AF=0.3,0.1	GT:AD:DP	1|2:0,10,5:15	0/0:20,0,0:20	0/2:10,0,10:20
chr1	35	rs4	G	C	70	PASS	DP=12;AF=0.4	GT:AD:DP	0/1:6,6:12	0/1:5,5:10	./.:0,0:0
```

**Example file:** `ref.fa`

```text
>chr1
ACGTACGTAC
GTACGTTTGA
CCAAGGTTAC
GATCGATCGA
ACGTACGTAC
```

**Example file:** `annotation.tsv`

```text
chr1	3	G	T	common
chr1	22	C	G	rare
```

## 1) `dnakit.variants.index_tabix_file()`

- **Purpose:** Build a tabix index for a VCF, BED or GFF file so that regions can be read directly.
- **Method:** Call `pysam.tabix_index(..., keep_original=True)`. A plain input is first compressed to `<name>.gz` beside it, and the original is kept. `preset` accepts only `vcf`, `bed` or `gff`. An existing index is not replaced unless `force=True`.
- **API:** `dnakit.variants.index_tabix_file(source, *, preset="vcf", force=False)`.
- **Input:** A plain-text or bgzip-compressed VCF, BED or GFF path.
- **Output:** A dictionary with `compressed` (the `.gz` path), `index` (the `.tbi` path), `preset` and `provenance`.
- **Example code:**

```python
from pathlib import Path

from dnakit.variants import index_tabix_file

result = index_tabix_file("sample.vcf", force=True, verbose=False)
print(Path(result["compressed"]).name, Path(result["index"]).name)
print(result["preset"])
```

- **Example output:**

```text
sample.vcf.gz sample.vcf.gz.tbi
vcf
```

## 2) `dnakit.variants.compress_bgzip_file()`

- **Purpose:** Write a text file as BGZF, the block-gzip format that tabix and htslib read.
- **Method:** Call `pysam.tabix_compress`. The output name must end in `.gz`. An existing file is not replaced unless `force=True`.
- **API:** `dnakit.variants.compress_bgzip_file(source, output_path, *, force=False)`.
- **Input:** A source text file and a new `.gz` output path.
- **Output:** A dictionary with `output`, the two file sizes and `provenance`.
- **Example code:**

```python
from dnakit.variants import compress_bgzip_file

result = compress_bgzip_file("annotation.tsv", "annotation.tsv.gz", force=True, verbose=False)
print(result["operation"], result["uncompressed_bytes"])
```

- **Example output:**

```text
compress_bgzip_file 35
```

## 3) `dnakit.variants.parse_genomic_region()`

- **Purpose:** Convert `chr1:start-end` (1-based, inclusive) to pysam's 0-based half-open interval.
- **Method:** Check the format strictly, then call `pysam.libcutils.parse_region`. A bare contig name gives the whole contig, with `end` set to `None`. Open-ended text such as `chr1:1000` is rejected.
- **API:** `dnakit.variants.parse_genomic_region(region)`.
- **Input:** A region string. Commas in numbers are allowed, for example `chr1:1,001-1,020`.
- **Output:** A dictionary with `reference`, `start` (0-based), `end` (exclusive) and `length`.
- **Example code:**

```python
from dnakit.variants import parse_genomic_region

result = parse_genomic_region("chr1:11-20", verbose=False)
print(result["reference"], result["start"], result["end"], result["length"])
```

- **Example output:**

```text
chr1 10 20 10
```

## 4) `dnakit.variants.parse_genomic_regions()`

- **Purpose:** Parse a comma-separated region list such as `chr1:10-20,chr2`.
- **Method:** Each piece is checked against the strict format, then passed to `pysam.libcvcf.parse_regions`. Upstream reads `chr1:10-20-30` as `chr1:10-20` without an error, so this function rejects such input first.
- **API:** `dnakit.variants.parse_genomic_regions(regions)`.
- **Input:** A comma-separated region string with at most 10000 pieces.
- **Output:** A dictionary whose `regions` items hold `reference`, `start` (0-based), `end` (exclusive; `None` for a whole contig) and `count`.
- **Example code:**

```python
from dnakit.variants import parse_genomic_regions

result = parse_genomic_regions("chr1:10-20,chr2", verbose=False)
for item in result["regions"]:
    print(item["reference"], item["start"], item["end"])
```

- **Example output:**

```text
chr1 9 20
chr2 0 None
```

## 5) `dnakit.variants.fetch_padded_sequence()`

- **Purpose:** Read one window of a reference sequence; positions outside the contig are filled with N.
- **Method:** Open the FASTA with `pysam.FastaFile`, then call `pysam.libcvcf.get_sequence`. The window is 0-based and half-open, and `start` may be negative. The FASTA needs a `.fai` index.
- **API:** `dnakit.variants.fetch_padded_sequence(fasta, reference, start, end)`.
- **Input:** A FASTA path with a `.fai` index, a contig name, `start` and `end`. The window is at most 10000000 bases.
- **Output:** A dictionary with `sequence` (upper case), `reference_length`, `left_padding` and `right_padding`.
- **Example code:**

```python
import pysam

from dnakit.variants import fetch_padded_sequence

pysam.faidx("ref.fa")
right = fetch_padded_sequence("ref.fa", "chr1", 45, 60, verbose=False)
left = fetch_padded_sequence("ref.fa", "chr1", -3, 4, verbose=False)
print(right["sequence"], right["right_padding"])
print(left["sequence"], left["left_padding"])
```

- **Example output:**

```text
CGTACNNNNNNNNNN 10
NNNACGT 3
```

## 6) `dnakit.variants.query_vcf_fields()`

- **Purpose:** Write selected VCF columns as text rows with `bcftools query`.
- **Method:** `fields` can be `CHROM`, `POS`, `ID`, `REF`, `ALT`, `QUAL`, `FILTER` or `INFO/<tag>`. With `genotypes=True`, each sample adds a `GT:<sample>` column. Values stay as text, and missing values remain `.`. A `region` needs an index.
- **API:** `dnakit.variants.query_vcf_fields(source, fields, *, region=None, genotypes=False, max_records=100000)`.
- **Input:** A VCF path, a list of fields, an optional region and `max_records` (1 to 1000000).
- **Output:** A dictionary with `columns`, `rows` (lists of strings), `record_count` and `provenance`.
- **Example code:**

```python
from dnakit.variants import query_vcf_fields

result = query_vcf_fields("sample.vcf", ["POS", "REF", "ALT", "FILTER", "INFO/DP"], verbose=False)
print(result["columns"])
for row in result["rows"]:
    print(row)
```

- **Example output:**

```text
['POS', 'REF', 'ALT', 'FILTER', 'INFO/DP']
['3', 'G', 'T', 'PASS', '20']
['10', 'C', 'CA', 'PASS', '10']
['22', 'C', 'G,A', 'q10', '30']
['35', 'G', 'C', 'PASS', '12']
```

## 7) `dnakit.variants.select_vcf_records()`

- **Purpose:** Write a new uncompressed VCF filtered by samples, region and PASS status.
- **Method:** Call `bcftools view -s -r -f PASS`. Samples keep the order given. A region needs an index. The output is not overwritten by default.
- **API:** `dnakit.variants.select_vcf_records(source, output_path, *, samples=None, region=None, pass_only=False, force=False)`.
- **Input:** A VCF path, a new output path, and optional sample names, a region and `pass_only`.
- **Output:** A dictionary with `output`, the number of kept records `records`, `samples` and `provenance`.
- **Example code:**

```python
from dnakit.variants import select_vcf_records

result = select_vcf_records(
    "sample.vcf.gz",
    "subset.vcf",
    samples=["S2", "S1"],
    region="chr1:9-22",
    pass_only=True,
    force=True,
    verbose=False,
)
print(result["records"], result["samples"])
```

- **Example output:**

```text
1 ['S2', 'S1']
```

## 8) `dnakit.variants.concatenate_vcf_files()`

- **Purpose:** Join several VCFs in the order given; typical use is combining regions of one sample set.
- **Method:** Call `bcftools concat`. Each input must be a bgzip file with an index. The output is an uncompressed VCF.
- **API:** `dnakit.variants.concatenate_vcf_files(sources, output_path, *, force=False)`.
- **Input:** At least two indexed VCF paths and a new output path.
- **Output:** A dictionary with `output`, the number of inputs `inputs`, the record count `records` and `provenance`.
- **Example code:**

```python
from pathlib import Path

import pysam

from dnakit.variants import concatenate_vcf_files

lines = Path("sample.vcf").read_text().splitlines(keepends=True)
header = [line for line in lines if line.startswith("#")]
body = [line for line in lines if not line.startswith("#")]
Path("part1.vcf").write_text("".join(header + body[:2]))
Path("part2.vcf").write_text("".join(header + body[2:]))
for name in ("part1.vcf", "part2.vcf"):
    pysam.tabix_index(name, preset="vcf", force=True, keep_original=True)
result = concatenate_vcf_files(
    ["part1.vcf.gz", "part2.vcf.gz"], "joined.vcf", force=True, verbose=False
)
print(result["inputs"], result["records"])
```

- **Example output:**

```text
2 4
```

## 9) `dnakit.variants.intersect_vcf_sites()`

- **Purpose:** Count sites private to each of two VCFs and sites shared by both.
- **Method:** Call `bcftools isec -p`. Sites match on chromosome, position and alleles, and each site is counted once. `output_dir` must not exist; bcftools writes VCF files and `sites.txt` into it.
- **API:** `dnakit.variants.intersect_vcf_sites(first, second, output_dir)`.
- **Input:** Two indexed VCF paths and a new directory path.
- **Output:** A dictionary with `private_first`, `private_second`, `shared`, `directory`, `files` and `provenance`.
- **Example code:**

```python
import tempfile
from pathlib import Path

import pysam

from dnakit.variants import intersect_vcf_sites

extra = Path("sample.vcf").read_text() + "chr1\t40\t.\tA\tT\t30\tPASS\t.\tGT\t0/1\t0/1\t0/1\n"
Path("more.vcf").write_text(extra)
pysam.tabix_index("more.vcf", preset="vcf", force=True, keep_original=True)
target = Path(tempfile.mkdtemp(), "isec")
result = intersect_vcf_sites("sample.vcf.gz", "more.vcf.gz", target, verbose=False)
print(result["private_first"], result["private_second"], result["shared"])
```

- **Example output:**

```text
0 1 4
```

## 10) `dnakit.variants.sort_vcf_records()`

- **Purpose:** Reorder VCF records by contig and position.
- **Method:** Call `bcftools sort` and write an uncompressed VCF. Contigs follow the order in the input header.
- **API:** `dnakit.variants.sort_vcf_records(source, output_path, *, force=False)`.
- **Input:** A VCF path (plain text is fine) and a new output path.
- **Output:** A dictionary with `output`, the record count `records` and `provenance`.
- **Example code:**

```python
from pathlib import Path

from dnakit.variants import sort_vcf_records

lines = Path("sample.vcf").read_text().splitlines(keepends=True)
header = [line for line in lines if line.startswith("#")]
body = [line for line in lines if not line.startswith("#")]
Path("shuffled.vcf").write_text("".join(header + body[::-1]))
result = sort_vcf_records("shuffled.vcf", "sorted.vcf", force=True, verbose=False)
positions = [line.split("\t")[1] for line in Path("sorted.vcf").read_text().splitlines()
             if not line.startswith("#")]
print(result["records"], positions)
```

- **Example output:**

```text
4 ['3', '10', '22', '35']
```

## 11) `dnakit.variants.annotate_vcf_from_tabix()`

- **Purpose:** Copy INFO tags from an annotation table onto VCF records at the same site.
- **Method:** Call `bcftools annotate -a -c -h`. The annotation table must be bgzip-compressed with a tabix index. `columns` follow the table's column order and must include `CHROM` and `POS` and at least one `INFO/<tag>`. Each new tag needs a matching `##INFO` line in `header_lines`.
- **API:** `dnakit.variants.annotate_vcf_from_tabix(source, annotations, columns, output_path, *, header_lines=(), force=False)`.
- **Input:** A VCF path, an annotation table path, a list of column names, a new output path and optional header lines.
- **Output:** A dictionary with `output`, the record count `records`, the number of annotated records per tag `annotated` and `provenance`.
- **Example code:**

```python
import pysam

from dnakit.variants import annotate_vcf_from_tabix

pysam.tabix_compress("annotation.tsv", "annotation.tsv.gz", force=True)
pysam.tabix_index("annotation.tsv.gz", seq_col=0, start_col=1, end_col=1, force=True)
result = annotate_vcf_from_tabix(
    "sample.vcf",
    "annotation.tsv.gz",
    ["CHROM", "POS", "REF", "ALT", "INFO/NOTE"],
    "annotated.vcf",
    header_lines=['##INFO=<ID=NOTE,Number=1,Type=String,Description="note">'],
    force=True,
    verbose=False,
)
print(result["annotated"], result["records"])
```

- **Example output:**

```text
{'NOTE': 2} 4
```

## 12) `dnakit.variants.build_vcf_consensus()`

- **Purpose:** Apply one sample's variants on one haplotype to a reference, producing a consensus sequence.
- **Method:** Call `bcftools consensus -f -s -H`. The VCF must be bgzip-compressed and indexed, and the FASTA needs a `.fai` index. At a heterozygous site the allele for that haplotype comes from the GT field.
- **API:** `dnakit.variants.build_vcf_consensus(source, fasta, sample, output_path, *, haplotype=1, force=False)`.
- **Input:** An indexed bgzip VCF, a FASTA, a sample name, the haplotype (1 or 2) and a new output path.
- **Output:** A dictionary with `output`, `sample`, `haplotype`, the length of each sequence `sequences` and `provenance`.
- **Example code:**

```python
from pathlib import Path

from dnakit.variants import build_vcf_consensus

result = build_vcf_consensus(
    "sample.vcf.gz", "ref.fa", "S2", "s2_hap1.fa", haplotype=1, force=True, verbose=False
)
print(result["sample"], result["haplotype"], result["sequences"])
print(Path("s2_hap1.fa").read_text().splitlines()[1][:20])
```

- **Example output:**

```text
S2 1 {'chr1': 50}
ACTTACGTACGTACGTTTGA
```

## 13) `dnakit.variants.summarize_sample_heterozygosity()`

- **Purpose:** Compute each sample's heterozygous fraction and allele-balance percentiles at given sites.
- **Method:** Call `cyvcf2.par_het`. The caller supplies the sites, so upstream's thinning rule is not used. Depths come from the FORMAT AD field, and only heterozygous calls with depth above `min_depth` enter the percentiles. Upstream skips sites with allele frequency below 0.01 or call rate below 0.8; `sites_used` lists the sites that were kept.
- **API:** `dnakit.variants.summarize_sample_heterozygosity(source, samples, sites, *, ncpus=1, min_depth=8, percentiles=(10, 90))`.
- **Input:** An indexed bgzip VCF, sample names and a list of sites, each `(chrom, pos)` or `(chrom, pos, ref, alt)`.
- **Output:** A dictionary whose `samples` entries hold `het_count`, `het_ratio`, `p10`, `p90`, `mean_depth` and related values, plus `sites_used` and `warnings`.
- **Example code:**

```python
from dnakit.variants import summarize_sample_heterozygosity

result = summarize_sample_heterozygosity(
    "sample.vcf.gz", ["S1", "S2", "S3"], [("chr1", 3), ("chr1", 10)], verbose=False
)
for name, values in result["samples"].items():
    print(name, values["het_count"], values["het_ratio"], values["p10"])
```

- **Example output:**

```text
S1 2 1.0 0.5
S2 0 0.0 0.0
S3 1 0.5 0.6
```

## 14) `dnakit.variants.summarize_pairwise_relatedness()`

- **Purpose:** Compute pairwise IBS counts and a relatedness coefficient between samples at given sites.
- **Method:** Call `cyvcf2.par_relatedness`, which matches `VCF.site_relatedness`. Each sample pair appears once, sorted by name. A sample with no heterozygous sites raises a warning, and its coefficient is then unreliable.
- **API:** `dnakit.variants.summarize_pairwise_relatedness(source, samples, sites, *, ncpus=1, min_depth=5)`.
- **Input:** An indexed bgzip VCF, sample names and a list of sites.
- **Output:** A dictionary whose `pairs` items hold `sample_a`, `sample_b`, `relatedness`, `ibs0`, `ibs2` and `informative_sites`, plus `warnings`.
- **Example code:**

```python
from dnakit.variants import summarize_pairwise_relatedness

result = summarize_pairwise_relatedness(
    "sample.vcf.gz",
    ["S1", "S2", "S3"],
    [("chr1", 3, "G", "T"), ("chr1", 10, "C", "CA")],
    verbose=False,
)
for pair in result["pairs"]:
    print(pair["sample_a"], pair["sample_b"], pair["relatedness"], pair["ibs2"])
```

- **Example output:**

```text
S1 S2 0.0 0
S1 S3 1.0 1
S2 S3 -2.0 0
```

## 15) `dnakit.variants.pearson_genotype_correlation()`

- **Purpose:** Compute the Pearson correlation of two genotype dosage vectors (0, 1, 2).
- **Method:** Call `cyvcf2.r_unphased`. A position is skipped when either call is missing. The result is `None` when fewer than two usable positions remain or when one side is constant over them. The upstream frequency argument is unused, so 0.0 is passed.
- **API:** `dnakit.variants.pearson_genotype_correlation(first, second)`.
- **Input:** Two equal-length lists whose items are 0, 1, 2 or `None` (missing).
- **Output:** A dictionary with `correlation`, `usable_pairs`, `missing_pairs` and `provenance`.
- **Example code:**

```python
from dnakit.variants import pearson_genotype_correlation

print(pearson_genotype_correlation([0, 1, 2, 1], [0, 2, 1, 1], verbose=False)["correlation"])
print(pearson_genotype_correlation([0, 1, None, 1], [0, 2, 1, 1], verbose=False)["correlation"])
```

- **Example output:**

```text
0.5
0.8660253882408142
```

## 16) `dnakit.variants.classify_vcf_alt()`

- **Purpose:** Classify one ALT allele, such as a substitution, a symbolic allele or a breakend.
- **Method:** Call `vcfpy.parser.process_alt`. Sequence alleles are then classified with `process_sub` as SNV, MNV, INS, DEL or INDEL. Breakend details are returned in `breakend`.
- **API:** `dnakit.variants.classify_vcf_alt(ref, alt)`.
- **Input:** REF bases (A, C, G, T, N) and an ALT string. ALT may not be `.`.
- **Output:** A dictionary with `kind`, `subtype`, `value`, `breakend` and `provenance`.
- **Example code:**

```python
from dnakit.variants import classify_vcf_alt

pairs = [("G", "T"), ("CAT", "C"), ("C", "CA"), ("C", "<DEL>"), ("C", "G]17:198982]")]
for ref, alt in pairs:
    result = classify_vcf_alt(ref, alt, verbose=False)
    print(ref, alt, result["kind"], result["subtype"], result["value"])
```

- **Example output:**

```text
G T substitution SNV T
CAT C substitution DEL C
C CA substitution INS CA
C <DEL> symbolic None DEL
C G]17:198982] breakend None None
```

## 17) `dnakit.variants.classify_vcf_substitution()`

- **Purpose:** Classify a REF/ALT pair as SNV, MNV, INS, DEL or INDEL.
- **Method:** Call `vcfpy.parser.process_sub`. Equal lengths give SNV or MNV; a longer ALT gives INS or INDEL; a shorter ALT gives DEL or INDEL. The helpers `process_sub_grow` and `process_sub_shrink` are called through it.
- **API:** `dnakit.variants.classify_vcf_substitution(ref, alt)`.
- **Input:** Two non-empty base strings, REF and ALT.
- **Output:** A dictionary with `subtype`, `value` and `provenance`.
- **Example code:**

```python
from dnakit.variants import classify_vcf_substitution

for ref, alt in [("AT", "A"), ("AC", "G"), ("ACG", "TTT")]:
    print(ref, alt, classify_vcf_substitution(ref, alt, verbose=False)["subtype"])
```

- **Example output:**

```text
AT A DEL
AC G INDEL
ACG TTT MNV
```

## 18) `dnakit.variants.parse_vcf_breakend()`

- **Purpose:** Parse a bracket-notation breakend ALT such as `G]17:198982]` into its mate position and orientations.
- **Method:** Call `vcfpy.parser.parse_breakend`. Non-bracket text is rejected. `orientation` and `mate_orientation` are `+` or `-`. A contig in angle brackets sets `within_main_assembly` to false.
- **API:** `dnakit.variants.parse_vcf_breakend(alt)`.
- **Input:** A breakend ALT string.
- **Output:** A dictionary with `mate_chrom`, `mate_pos`, `orientation`, `mate_orientation`, `sequence`, `within_main_assembly` and `provenance`.
- **Example code:**

```python
from dnakit.variants import parse_vcf_breakend

result = parse_vcf_breakend("]13:123456]T", verbose=False)
print(result["mate_chrom"], result["mate_pos"], result["orientation"],
      result["mate_orientation"], result["sequence"])
```

- **Example output:**

```text
13 123456 + - T
```

## 19) `dnakit.variants.parse_vcf_field_value()`

- **Purpose:** Convert one INFO or FORMAT value to Python objects using its declared Number and Type.
- **Method:** Call `vcfpy.parser.parse_field_value`. `number` is an integer or one of `A`, `R`, `G`, `.`. A Number of 1 returns one value, and other Numbers return a list; `.` means missing. A value that cannot be converted raises an error instead of staying a string.
- **API:** `dnakit.variants.parse_vcf_field_value(text, *, value_type, number, field_id='')`.
- **Input:** The field text, its type (`Integer`, `Float`, `Flag`, `Character` or `String`) and its Number.
- **Output:** A dictionary with `value`, `value_type`, `number` and `provenance`.
- **Example code:**

```python
from dnakit.variants import parse_vcf_field_value

print(parse_vcf_field_value("0.5,.", value_type="Float", number="A", verbose=False)["value"])
print(parse_vcf_field_value("12", value_type="Integer", number=1, verbose=False)["value"])
```

- **Example output:**

```text
[0.5, None]
12
```

## 20) `dnakit.variants.format_vcf_field_value()`

- **Purpose:** Format one INFO or FORMAT value as VCF text; this is the inverse of the parser.
- **Method:** Call `vcfpy.writer.format_value`. A value that does not match the type is rejected. `Flag` values are not formatted here; write the bare key instead. `None` becomes `.`.
- **API:** `dnakit.variants.format_vcf_field_value(value, *, value_type, number, section='INFO', field_id='')`.
- **Input:** The value, its type, its Number, `section` (`INFO` or `FORMAT`) and an optional field ID.
- **Output:** A dictionary with `text`, `value_type`, `number`, `section` and `provenance`.
- **Example code:**

```python
from dnakit.variants import format_vcf_field_value

print(format_vcf_field_value([0.5, 0.25], value_type="Float", number="A", verbose=False)["text"])
print(format_vcf_field_value(None, value_type="Integer", number=1, verbose=False)["text"])
```

- **Example output:**

```text
0.5,0.25
.
```

## 21) `dnakit.variants.remove_vcf_header_lines()`

- **Purpose:** Read a VCF header and remove the lines named by `(key, ID)` pairs.
- **Method:** Read the header with `vcfpy.Reader`, then call `vcfpy.header_without_lines`. For lines without an ID, give the whole value in place of the ID. The `#CHROM` line is not part of the result.
- **API:** `dnakit.variants.remove_vcf_header_lines(source, remove)`.
- **Input:** A VCF path (a compressed file is fine) and a list of `(key, ID)` pairs.
- **Output:** A dictionary with `header_lines` (lines starting with `##`), `samples`, `removed_lines` and `provenance`.
- **Example code:**

```python
from dnakit.variants import remove_vcf_header_lines

result = remove_vcf_header_lines("sample.vcf", [("FILTER", "q10")], verbose=False)
print(result["removed_lines"], len(result["header_lines"]), result["samples"])
```

- **Example output:**

```text
1 7 ['S1', 'S2', 'S3']
```

## Installing the optional dependencies

```bash
python -m pip install "dnakit[bio-reads]"        # pysam (bcftools and htslib bundled)
python -m pip install "dnakit[bio-cyvcf2]"       # cyvcf2
python -m pip install "dnakit[bio-vcfpy]"        # vcfpy
python -m pip install "dnakit[bio-analysis]"     # numpy, used by pearson_genotype_correlation
```

Project pages: [pysam](https://github.com/pysam-developers/pysam), [cyvcf2](https://github.com/brentp/cyvcf2), [vcfpy](https://github.com/dahlo/vcfpy).

## How these functions are verified

- The test file is `others/tests/unit/test_full_G5.py`, with 24 tests that run in about one second.
- pysam tabix indexes, bgzip output, region parsing and FASTA windows are checked against hand-calculated values. bcftools query, view, concat, isec, sort and annotate results are checked the same way. The consensus test also compares against a direct `pysam.bcftools.consensus` call.
- cyvcf2 heterozygosity is compared with `VCF.het_check`, and relatedness with `VCF.site_relatedness`. `pearson_genotype_correlation` is checked against hand-calculated values.
- vcfpy ALT and substitution classes are compared with `vcfpy.parser.process_alt` and `process_sub`. Breakend parsing and field-value parsing and formatting are checked by hand. Header removal is compared with `vcfpy.header_without_lines`.
- Each module checks that a missing package raises `BackendUnavailableError` with the extra's install hint. Invalid input is rejected, for example an unindexed or uncompressed VCF, an unknown sample, an oversized window or an unsupported field.

## Current scope

- Upstream operations wrapped here: pysam's tabix and bgzip, region parsing and FASTA window reads, and the bcftools query, view, concat, isec, sort, annotate and consensus commands; cyvcf2's `par_het`, `par_relatedness` and `r_unphased`; and vcfpy's ALT and substitution classes, breakend parsing, field-value parsing and formatting, and header-line removal.
- Not wrapped from pysam: `bcftools csq` needs a GFF3 gene model and a phasing policy for heterozygous calls; `roh`, `gtcheck`, `merge`, `cnv`, `convert` and `filter` are outside this batch; `pysam.Pileup` is the legacy `samtools pileup -c` interface, and `tabix_iterator` needs pysam parser objects. BAM and SAM commands from samtools are not variant operations and are not wrapped here.
- Already covered elsewhere: `bcftools norm` and `bcftools stats` duplicate `dnakit.databases.vcf_normalize` and `vcf_summary_stats`; `bcftools call` duplicates `dnakit.ngs.call_variants`; `samtools faidx` duplicates `dnakit.bioinformatics.reads.index_reference`.
- The command-line printers in cyvcf2 (`cyvcf2.cli`) and the bcftools plugin runner are not wrapped.
