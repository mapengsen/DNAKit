# Review G18k: varcode and scikit-allel variant utilities

This page records how the functions with review status `other` or `cap` in varcode and scikit-allel were handled. It adds 10 varcode wrappers and 31 scikit-allel wrappers. All of them run offline and accept only local inputs. The modules are `dnakit.variants` and `dnakit.variants`; they are not yet listed in `dnakit.variants`, which the main line wires later. Dependencies: `bio-varcode` (varcode) and `bio-scikit-allel` (scikit-allel; HDF5 archives also need h5py).

## 1) Test whether one nucleotide is a purine

- **Purpose:** Report whether one nucleotide code is a purine (A or G).
- **Method:** Calls varcode.nucleotides.is_purine; non-standard codes raise unless allow_extended is True.
- **API:** `is_purine_base(nucleotide, *, allow_extended=False, verbose=True)`
- **Input:** One letter `nucleotide`; optional `allow_extended`.
- **Output:** Dict with uppercase `nucleotide`, `is_purine` and `provenance`.
- **Example code:**

```python
from dnakit.variants import is_purine_base

result = is_purine_base("A", verbose=False)
print(result["nucleotide"], result["is_purine"])
```

- **Example output:**

```text
A True
```

## 2) Check that a sequence holds only A, C, G and T

- **Purpose:** Check whether a base string holds only standard bases and list the 1-based positions of the rest.
- **Method:** Calls varcode.nucleotides.all_standard_nucleotides; lists each position outside A, C, G and T.
- **API:** `has_only_standard_nucleotides(nucleotides, *, verbose=True)`
- **Input:** Non-empty letter string `nucleotides` without whitespace, at most 1,000,000 characters.
- **Output:** Returns `length`, `standard` (boolean) and `non_standard_positions` (1-based list).
- **Example code:**

```python
from dnakit.variants import has_only_standard_nucleotides

result = has_only_standard_nucleotides("ACGTNACG", verbose=False)
print(result["standard"], result["non_standard_positions"])
```

- **Example output:**

```text
False [5]
```

## 3) Test whether a contig name is mitochondrial

- **Purpose:** Report whether names such as chrM, MT or mtDNA refer to the mitochondrial genome, which selects a codon table.
- **Method:** Calls varcode.effects.codon_tables.is_mitochondrial_contig; case-insensitive, None returns False.
- **API:** `is_mitochondrial_contig(contig, *, verbose=True)`
- **Input:** `contig` string or None.
- **Output:** Returns `contig` and the boolean `mitochondrial`.
- **Example code:**

```python
from dnakit.variants import is_mitochondrial_contig

print(is_mitochondrial_contig("chrM", verbose=False)["mitochondrial"], is_mitochondrial_contig("chr1", verbose=False)["mitochondrial"])
```

- **Example output:**

```text
True False
```

## 4) Test whether a genome name is a UCSC assembly name

- **Purpose:** Report whether a genome name is one of varcode's UCSC assembly names, such as hg19 or mm10.
- **Method:** Calls varcode.reference.is_ucsc_reference_name; the name is normalized for case and separators first.
- **API:** `is_ucsc_reference_name(name, *, verbose=True)`
- **Input:** Non-empty string `name`.
- **Output:** Returns `name` and the boolean `ucsc`.
- **Example code:**

```python
from dnakit.variants import is_ucsc_reference_name

print(is_ucsc_reference_name("hg38", verbose=False)["ucsc"], is_ucsc_reference_name("GRCh38", verbose=False)["ucsc"])
```

- **Example output:**

```text
True False
```

## 5) Infer the assembly name from a reference name or path

- **Purpose:** Infer the canonical assembly name from a reference name or FASTA path without reading the file.
- **Method:** Calls varcode.reference.infer_reference_name, matching the file name against the assembly table; no match raises ConfigurationError.
- **API:** `infer_reference_assembly(reference, *, verbose=True)`
- **Input:** Reference name or path string `reference`.
- **Output:** Returns `reference` and `assembly`.
- **Example code:**

```python
from dnakit.variants import infer_reference_assembly

print(infer_reference_assembly("/data/refs/GRCh37.fa", verbose=False)["assembly"])
```

- **Example output:**

```text
GRCh37
```

## 6) Span of the DEL, DUP or INV that a breakend junction describes

- **Purpose:** Give the DEL, DUP or INV span described by one same-contig novel junction.
- **Method:** Calls varcode.structural_variant.junction_span_details; both ends must share a contig and state their kept side, otherwise `span` is None.
- **API:** `breakend_junction_span(near, far, *, verbose=True)`
- **Input:** Mappings `near` and `far` with `contig`, 1-based `position` and `keeps` ("left", "right" or None).
- **Output:** `span` is None or a dict with `sv_type`, `start`, `end`, `affected_start` and `affected_end`.
- **Example code:**

```python
from dnakit.variants import breakend_junction_span

span = breakend_junction_span(
    {"contig": "chr1", "position": 100, "keeps": "left"},
    {"contig": "chr1", "position": 200, "keeps": "right"},
    verbose=False,
)["span"]
print(span["sv_type"], span["start"], span["end"], span["affected_start"], span["affected_end"])
```

- **Example output:**

```text
DEL 100 199 101 199
```

## 7) Join two local layouts at breakend sides

- **Purpose:** Keep one breakend side of each reference interval and join them into an oriented layout with the novel junction.
- **Method:** Calls varcode.genomic_layout.join_breakends; segments are reversed when needed; no sequence is read.
- **API:** `join_breakend_layouts(first_segments, first_breakend, second_segments, second_breakend, *, verbose=True)`
- **Input:** Two segment lists (`contig`, `start`, `end`, optional `strand`) and two breakends (`contig`, `position`, `keeps` as "left" or "right").
- **Output:** Returns `segments` (contig, start, end, strand, length) and the total `length`.
- **Example code:**

```python
from dnakit.variants import join_breakend_layouts

result = join_breakend_layouts(
    [{"contig": "chr1", "start": 100, "end": 110}],
    {"contig": "chr1", "position": 105, "keeps": "left"},
    [{"contig": "chr1", "start": 300, "end": 310}],
    {"contig": "chr1", "position": 302, "keeps": "right"},
    verbose=False,
)
print([(item["start"], item["end"]) for item in result["segments"]], result["length"])
```

- **Example output:**

```text
[(100, 105), (302, 310)] 15
```

## 8) Project reference intervals onto a local layout

- **Purpose:** Project reference intervals onto a local layout, one piece per continuous match.
- **Method:** Calls varcode.genomic_layout.project_intervals; intersects coordinates without expanding bases; inverted segments do not match.
- **API:** `project_intervals_onto_layout(segments, intervals, *, verbose=True)`
- **Input:** Layout `segments` and reference `intervals` (`contig`, `start`, `end`; 1-based inclusive).
- **Output:** Returns `count` and `pieces`; each piece has `segments` and `length`.
- **Example code:**

```python
from dnakit.variants import project_intervals_onto_layout

result = project_intervals_onto_layout(
    [{"contig": "chr1", "start": 100, "end": 110}],
    [{"contig": "chr1", "start": 103, "end": 120}],
    verbose=False,
)
piece = result["pieces"][0]["segments"][0]
print(result["count"], piece["start"], piece["end"], result["pieces"][0]["length"])
```

- **Example output:**

```text
1 103 110 8
```

## 9) Merge MATEID-paired BND records of a local VCF

- **Purpose:** Read a local VCF and merge BND records that name each other as MATEID into one combined record per event.
- **Method:** Reads the file with varcode.load_vcf (structural variants parsed), then calls varcode.transforms.pair_breakends; `genome` only names contigs.
- **API:** `pair_breakend_mates(path, *, genome="GRCh38", verbose=True)`
- **Input:** Local VCF `path`; optional `genome` (default GRCh38).
- **Output:** Returns `input_count`, `output_count`, `paired_count` and position-sorted `records` (contig, start, ref, alt, source_starts).
- **Example code:**

```python
from pathlib import Path
from tempfile import TemporaryDirectory

from dnakit.variants import pair_breakend_mates

ROWS = [
    "##fileformat=VCFv4.2",
    "##contig=<ID=chr1,length=1000>",
    '##ALT=<ID=DEL,Description="Deletion">',
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">',
    '##INFO=<ID=MATEID,Number=1,Type=String,Description="Mate">',
    '##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Length">',
    '##INFO=<ID=END,Number=1,Type=Integer,Description="End">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
    "chr1\t100\tbnd_a\tA\tA]chr1:300]\t50\tPASS\tSVTYPE=BND;MATEID=bnd_b\tGT\t0/1",
    "chr1\t300\tbnd_b\tT\t[chr1:100[T\t50\tPASS\tSVTYPE=BND;MATEID=bnd_a\tGT\t0/1",
    "chr1\t500\tdel1\tN\t<DEL>\t50\tPASS\tSVTYPE=DEL;END=600;SVLEN=-100\tGT\t0/1",
]
with TemporaryDirectory() as folder:
    path = Path(folder) / "bnd.vcf"
    path.write_text("\n".join(ROWS) + "\n", encoding="utf-8")
    result = pair_breakend_mates(str(path), verbose=False)
print(result["input_count"], result["output_count"], result["paired_count"])
for row in result["records"]:
    print(row["contig"], row["start"], row["alt"], row["source_starts"])
```

- **Example output:**

```text
3 2 1
chr1 100 A]chr1:300] [100, 300]
chr1 500 <DEL> []
```

## 10) Rewrite a local VCF with varcode rules

- **Purpose:** Read a local VCF and write a new VCF that follows varcode's INFO and FORMAT rules.
- **Method:** Calls varcode.vcf_output.variants_to_vcf; the output must be a new file, and an unfinished file is removed on failure.
- **API:** `rewrite_vcf_records(path, output, *, genome="GRCh38", verbose=True)`
- **Input:** Source VCF `path`, new output `output` (parent directory must exist), optional `genome`.
- **Output:** Returns the `output` path and the `variants` count.
- **Example code:**

```python
from pathlib import Path
from tempfile import TemporaryDirectory

from dnakit.variants import rewrite_vcf_records

ROWS = [
    "##fileformat=VCFv4.2",
    "##contig=<ID=chr1,length=1000>",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3",
    "chr1\t10\trs1\tA\tG\t50\tPASS\tDP=10\tGT\t0|1\t0/0\t1|1",
    "chr1\t20\trs2\tC\tT\t50\tPASS\tDP=12\tGT\t0/1\t1/1\t0/0",
    "chr1\t35\trs3\tG\tA\t50\tPASS\tDP=8\tGT\t1/1\t0/1\t0|0",
]
with TemporaryDirectory() as folder:
    source = Path(folder) / "small.vcf"
    source.write_text("\n".join(ROWS) + "\n", encoding="utf-8")
    output = Path(folder) / "rewritten.vcf"
    result = rewrite_vcf_records(str(source), str(output), verbose=False)
    last = output.read_text(encoding="utf-8").splitlines()[-1]
print(result["variants"], " ".join(last.split("\t")))
```

- **Example output:**

```text
3 chr1 35 rs3 G A 50.0 PASS DP=8 GT 1/1 0/1 0|0
```

## 11) Find alleles present in exactly one population

- **Purpose:** From allele counts of two to fifty populations, find alleles that occur in exactly one population.
- **Method:** Calls allel.locate_private_alleles; an allele is private when exactly one population has a count above 0.
- **API:** `locate_population_private_alleles(allele_counts, *, verbose=True)`
- **Input:** Per-population count matrices (variants, alleles) with the same number of variants.
- **Output:** Returns `private_alleles` (variant_index, allele_index) and `private_variant_indices`.
- **Example code:**

```python
from dnakit.variants import locate_population_private_alleles

counts = [
    [[4, 0], [2, 2], [0, 4]],
    [[0, 4], [2, 2], [4, 0]],
    [[4, 0], [2, 2], [1, 3]],
]
result = locate_population_private_alleles(counts, verbose=False)
print(result["private_variant_indices"], result["private_alleles"])
```

- **Example output:**

```text
[0] [{'variant_index': 0, 'allele_index': 1}]
```

## 12) Find fixed differences between two populations

- **Purpose:** Find variants with no allele shared between two populations (fixed differences).
- **Method:** Calls allel.locate_fixed_differences.
- **API:** `locate_population_fixed_differences(first_counts, second_counts, *, verbose=True)`
- **Input:** Two allele-count matrices (variants, alleles) with the same number of rows.
- **Output:** Returns `fixed_indices` and `fixed_count`.
- **Example code:**

```python
from dnakit.variants import locate_population_fixed_differences

result = locate_population_fixed_differences([[4, 0], [2, 2], [0, 4]], [[0, 4], [2, 2], [3, 1]], verbose=False)
print(result["fixed_indices"], result["fixed_count"])
```

- **Example output:**

```text
[0] 1
```

## 13) Count variants in windows

- **Purpose:** Count the variants in each fixed-length window along one contig.
- **Method:** Calls allel.windowed_count; windows are 1-based inclusive and the default step equals the window size.
- **API:** `count_variants_in_windows(positions, *, window_size, start=None, stop=None, step=None, verbose=True)`
- **Input:** Sorted 1-based `positions`, `window_size`, and optional `start`, `stop`, `step`.
- **Output:** Returns `windows`, each with `start`, `stop` and `variants`.
- **Example code:**

```python
from dnakit.variants import count_variants_in_windows

result = count_variants_in_windows([10, 20, 30, 40, 50, 60], window_size=20, verbose=False)
print([(w["start"], w["stop"], w["variants"]) for w in result["windows"]])
```

- **Example output:**

```text
[(10, 29, 2), (30, 49, 2), (50, 60, 2)]
```

## 14) Divide a windowed statistic by callable bases

- **Purpose:** Divide each window's statistic by its number of bases, or by its accessible bases when a mask is given.
- **Method:** Calls allel.per_base; with `is_accessible`, only bases marked True are counted.
- **API:** `windowed_values_per_base(values, windows, *, is_accessible=None, verbose=True)`
- **Input:** One value per window in `values`; `windows` as 1-based inclusive (start, stop); optional boolean mask from base 1.
- **Output:** Returns `windows`, each with `start`, `stop`, `callable_bases` and `value_per_base`.
- **Example code:**

```python
from dnakit.variants import windowed_values_per_base

result = windowed_values_per_base([1.0, 2.0], [[1, 3], [4, 6]], verbose=False)
print([round(w["value_per_base"], 4) for w in result["windows"]])
```

- **Example output:**

```text
[0.3333, 0.6667]
```

## 15) Build windows with equal accessible bases

- **Purpose:** From an accessibility mask, build windows that each hold the same number of accessible bases.
- **Method:** Calls allel.equally_accessible_windows; `window_size` counts accessible bases.
- **API:** `build_equally_accessible_windows(is_accessible, *, window_size, step=None, verbose=True)`
- **Input:** Boolean mask `is_accessible` from base 1, `window_size` and optional `step`.
- **Output:** Returns `windows`, each with `start`, `stop` and `accessible_bases`.
- **Example code:**

```python
from dnakit.variants import build_equally_accessible_windows

mask = [True, True, False, True, True, True]
result = build_equally_accessible_windows(mask, window_size=2, verbose=False)
print([(w["start"], w["stop"]) for w in result["windows"]])
```

- **Example output:**

```text
[(1, 2), (4, 5)]
```

## 16) Between-population mean pairwise differences

- **Purpose:** Compute the mean number of pairwise differences between two populations for each variant.
- **Method:** Calls allel.mean_pairwise_difference_between from allele counts.
- **API:** `calculate_between_population_pi(first_counts, second_counts, *, verbose=True)`
- **Input:** Two allele-count matrices (variants, alleles).
- **Output:** Returns `values`, one per variant; undefined values are None.
- **Example code:**

```python
from dnakit.variants import calculate_between_population_pi

first = [[4, 0], [2, 2], [3, 1]]
second = [[0, 4], [2, 2], [1, 3]]
result = calculate_between_population_pi(first, second, verbose=False)
print(result["values"])
```

- **Example output:**

```text
[1.0, 0.5, 0.625]
```

## 17) Fixed-difference density in windows

- **Purpose:** Compute the density of fixed differences between two populations in each window.
- **Method:** Calls allel.windowed_df; density is the fixed-difference count divided by callable bases in the window.
- **API:** `calculate_windowed_fixed_difference_density(positions, first_counts, second_counts, *, window_size, start=None, stop=None, step=None, is_accessible=None, verbose=True)`
- **Input:** Positions, two count matrices, window parameters and an optional accessibility mask.
- **Output:** Returns `windows`, each with `start`, `stop`, `sites`, `callable_bases` and `density`.
- **Example code:**

```python
from dnakit.variants import calculate_windowed_fixed_difference_density

positions = [10, 20, 30, 40]
first = [[4, 0], [2, 2], [0, 4], [3, 1]]
second = [[0, 4], [2, 2], [4, 0], [1, 3]]
result = calculate_windowed_fixed_difference_density(positions, first, second, window_size=20, verbose=False)
print([(w["start"], w["stop"], w["sites"], w["callable_bases"], w["density"]) for w in result["windows"]])
```

- **Example output:**

```text
[(10, 29, 2, 20, 0.05), (30, 40, 2, 11, 0.0909090909)]
```

## 18) Patterson Fst between two populations

- **Purpose:** Compute Patterson's Fst between two populations, with per-variant numerator and denominator and their overall ratio.
- **Method:** Calls allel.patterson_fst; the overall Fst is the sum of numerators divided by the sum of denominators.
- **API:** `calculate_patterson_fst(first_counts, second_counts, *, verbose=True)`
- **Input:** Two allele-count matrices (variants, alleles).
- **Output:** Returns `numerator`, `denominator` and `fst`.
- **Example code:**

```python
from dnakit.variants import calculate_patterson_fst

first = [[4, 0], [2, 2], [0, 4], [3, 1]]
second = [[0, 4], [2, 2], [4, 0], [1, 3]]
result = calculate_patterson_fst(first, second, verbose=False)
print(result["numerator"], result["denominator"], result["fst"])
```

- **Example output:**

```text
[1.0, -0.1666666667, 1.0, 0.125] [1.0, 0.5, 1.0, 0.625] 0.6266666667
```

## 19) Patterson Fst in windows

- **Purpose:** Compute Patterson's Fst in each base-pair window.
- **Method:** Calls allel.windowed_patterson_fst.
- **API:** `calculate_windowed_patterson_fst(positions, first_counts, second_counts, *, window_size, start=None, stop=None, step=None, verbose=True)`
- **Input:** Positions, two count matrices and window parameters.
- **Output:** Returns `windows`, each with `start`, `stop`, `sites` and `fst`.
- **Example code:**

```python
from dnakit.variants import calculate_windowed_patterson_fst

positions = [10, 20, 30, 40]
first = [[4, 0], [2, 2], [0, 4], [3, 1]]
second = [[0, 4], [2, 2], [4, 0], [1, 3]]
result = calculate_windowed_patterson_fst(positions, first, second, window_size=20, verbose=False)
print([(w["start"], w["stop"], w["sites"], w["fst"]) for w in result["windows"]])
```

- **Example output:**

```text
[(10, 29, 2, 0.5555555556), (30, 40, 2, 0.6923076923)]
```

## 20) Patterson Fst with block-jackknife error

- **Purpose:** Compute Patterson's Fst and its block-jackknife standard error.
- **Method:** Calls allel.blockwise_patterson_fst; `block_length` counts variants.
- **API:** `calculate_blockwise_patterson_fst(first_counts, second_counts, *, block_length, verbose=True)`
- **Input:** Two count matrices and `block_length`.
- **Output:** Returns `fst`, `standard_error`, `block_values` and `jackknife_values`.
- **Example code:**

```python
from dnakit.variants import calculate_blockwise_patterson_fst

first = [[4, 0], [2, 2], [0, 4], [3, 1]]
second = [[0, 4], [2, 2], [4, 0], [1, 3]]
result = calculate_blockwise_patterson_fst(first, second, block_length=2, verbose=False)
print(result["fst"], result["standard_error"], result["block_values"])
```

- **Example output:**

```text
0.6266666667 0.0683760684 [0.5555555556, 0.6923076923]
```

## 21) Patterson Fst in moving variant windows

- **Purpose:** Compute Patterson's Fst in moving windows of consecutive variants.
- **Method:** Calls allel.moving_patterson_fst; windows count variants, not base positions.
- **API:** `calculate_moving_patterson_fst(first_counts, second_counts, *, window_size, start=0, stop=None, step=None, verbose=True)`
- **Input:** Two count matrices and window parameters.
- **Output:** Returns `window_size`, `step` and `values`.
- **Example code:**

```python
from dnakit.variants import calculate_moving_patterson_fst

first = [[4, 0], [2, 2], [0, 4], [3, 1]]
second = [[0, 4], [2, 2], [4, 0], [1, 3]]
result = calculate_moving_patterson_fst(first, second, window_size=2, verbose=False)
print(result["values"])
```

- **Example output:**

```text
[0.5555555556, 0.6923076923]
```

## 22) Pairwise dxy matrix between samples

- **Purpose:** Compute pairwise nucleotide divergence (dxy) between samples as a symmetric square matrix.
- **Method:** Calls allel.pairwise_dxy for the condensed matrix and expands it to a square matrix.
- **API:** `calculate_pairwise_dxy(positions, genotype_allele_counts, *, start=None, stop=None, is_accessible=None, verbose=True)`
- **Input:** Positions, counts shaped (variants, samples, alleles), optional region and accessibility mask.
- **Output:** Returns `samples` and `matrix`.
- **Example code:**

```python
from dnakit.variants import calculate_pairwise_dxy

positions = [10, 20, 30]
genotypes = [
    [[4, 0], [2, 2], [0, 4]],
    [[3, 1], [1, 3], [2, 2]],
    [[4, 0], [4, 0], [0, 4]],
]
result = calculate_pairwise_dxy(positions, genotypes, verbose=False)
print(result["samples"], result["matrix"][0])
```

- **Example output:**

```text
3 [0.0, 0.0535714286, 0.119047619]
```

## 23) Patterson f2 between two populations

- **Purpose:** Compute Patterson's f2 branch-length estimate per variant, with the mean over defined variants.
- **Method:** Calls allel.patterson_f2; the mean uses only defined variants.
- **API:** `calculate_patterson_f2(first_counts, second_counts, *, verbose=True)`
- **Input:** Two count matrices.
- **Output:** Returns `values` and `mean`.
- **Example code:**

```python
from dnakit.variants import calculate_patterson_f2

first = [[4, 0], [2, 2], [0, 4], [3, 1]]
second = [[0, 4], [2, 2], [4, 0], [1, 3]]
result = calculate_patterson_f2(first, second, verbose=False)
print(result["values"], result["mean"])
```

- **Example output:**

```text
[1.0, -0.1666666667, 1.0, 0.125] 0.4895833333
```

## 24) Patterson f3 and f3* admixture statistic

- **Purpose:** Compute the three-population admixture statistic f3(C; A, B): numerator and heterozygosity sums, and f3*.
- **Method:** Calls allel.patterson_f3; f3* is the numerator sum divided by the heterozygosity sum.
- **API:** `calculate_patterson_f3(target_counts, first_counts, second_counts, *, verbose=True)`
- **Input:** Count matrices for target C and populations A and B.
- **Output:** Returns `numerator`, `heterozygosity`, `numerator_sum`, `heterozygosity_sum` and `f3_star`.
- **Example code:**

```python
from dnakit.variants import calculate_patterson_f3

target = [[4, 0], [2, 2], [1, 3], [3, 1]]
first = [[4, 0], [2, 2], [0, 4], [3, 1]]
second = [[0, 4], [2, 2], [4, 0], [1, 3]]
result = calculate_patterson_f3(target, first, second, verbose=False)
print(result["numerator_sum"], result["heterozygosity_sum"], result["f3_star"])
```

- **Example output:**

```text
-0.3958333333 1.6666666667 -0.2375
```

## 25) Patterson f3 with block-jackknife error

- **Purpose:** Compute f3 with its block-jackknife standard error and z-score.
- **Method:** Calls allel.blockwise_patterson_f3; undefined values are None.
- **API:** `calculate_blockwise_patterson_f3(target_counts, first_counts, second_counts, *, block_length, normed=True, verbose=True)`
- **Input:** Three count matrices, `block_length`, optional `normed`.
- **Output:** Returns `f3`, `standard_error`, `z_score`, `block_values` and `jackknife_values`.
- **Example code:**

```python
from dnakit.variants import calculate_blockwise_patterson_f3

target = [[4, 0], [2, 2], [1, 3], [3, 1]]
first = [[4, 0], [2, 2], [0, 4], [3, 1]]
second = [[0, 4], [2, 2], [4, 0], [1, 3]]
result = calculate_blockwise_patterson_f3(target, first, second, block_length=2, verbose=False)
print(result["f3"], result["standard_error"], result["z_score"])
```

- **Example output:**

```text
-0.2375 0.09375 -2.5333333333
```

## 26) Patterson f3 in moving variant windows

- **Purpose:** Compute f3 in moving windows of consecutive variants.
- **Method:** Calls allel.moving_patterson_f3.
- **API:** `calculate_moving_patterson_f3(target_counts, first_counts, second_counts, *, window_size, start=0, stop=None, step=None, normed=True, verbose=True)`
- **Input:** Three count matrices and window parameters.
- **Output:** Returns `window_size`, `step`, `normed` and `values`.
- **Example code:**

```python
from dnakit.variants import calculate_moving_patterson_f3

target = [[4, 0], [2, 2], [1, 3], [3, 1]]
first = [[4, 0], [2, 2], [0, 4], [3, 1]]
second = [[0, 4], [2, 2], [4, 0], [1, 3]]
result = calculate_moving_patterson_f3(target, first, second, window_size=2, verbose=False)
print(result["values"])
```

- **Example output:**

```text
[-0.125, -0.3125]
```

## 27) Paint haplotypes by shared prefix

- **Purpose:** Sort haplotypes by prefix and give each shared prefix one integer.
- **Method:** Calls allel.voight_painting; accepts only biallelic 0/1 matrices without missing calls.
- **API:** `paint_shared_haplotype_prefixes(haplotypes, *, verbose=True)`
- **Input:** 0/1 matrix shaped (variants, haplotypes).
- **Output:** Returns the `painting` matrix and `haplotype_order`.
- **Example code:**

```python
from dnakit.variants import paint_shared_haplotype_prefixes

haplotypes = [[0, 0, 1, 1], [0, 0, 1, 0], [1, 1, 0, 0], [0, 1, 0, 1]]
result = paint_shared_haplotype_prefixes(haplotypes, verbose=False)
print(result["haplotype_order"], result["painting"][0])
```

- **Example output:**

```text
[0, 1, 3, 2] [1, 1, 2, 2]
```

## 28) Haplotype diversity

- **Purpose:** Compute haplotype diversity: one minus the sum of squared haplotype frequencies, with an n/(n-1) correction.
- **Method:** Calls allel.haplotype_diversity.
- **API:** `calculate_haplotype_diversity(haplotypes, *, verbose=True)`
- **Input:** 0/1 matrix shaped (variants, haplotypes), at least two haplotypes.
- **Output:** Returns `diversity`, `variants` and `haplotypes`.
- **Example code:**

```python
from dnakit.variants import calculate_haplotype_diversity

haplotypes = [[0, 0, 1, 1], [0, 0, 1, 0], [1, 1, 0, 0], [0, 1, 0, 1]]
result = calculate_haplotype_diversity(haplotypes, verbose=False)
print(result["haplotypes"], result["diversity"])
```

- **Example output:**

```text
4 1.0
```

## 29) Haplotype diversity in moving windows

- **Purpose:** Compute haplotype diversity in moving windows of consecutive variants.
- **Method:** Calls allel.moving_haplotype_diversity.
- **API:** `calculate_moving_haplotype_diversity(haplotypes, *, window_size, start=0, stop=None, step=None, verbose=True)`
- **Input:** 0/1 haplotype matrix and window parameters.
- **Output:** Returns `window_size`, `step` and `values`.
- **Example code:**

```python
from dnakit.variants import calculate_moving_haplotype_diversity

haplotypes = [[0, 0, 1, 1], [0, 0, 1, 0], [1, 1, 0, 0], [0, 1, 0, 1], [1, 0, 1, 0], [1, 1, 1, 0]]
result = calculate_moving_haplotype_diversity(haplotypes, window_size=3, verbose=False)
print(result["values"])
```

- **Example output:**

```text
[0.8333333333, 0.8333333333]
```

## 30) Cross-population nSL (XP-nSL)

- **Purpose:** Compute the unstandardized cross-population nSL score between two haplotype sets.
- **Method:** Calls allel.xpnsl with threading disabled so results are deterministic.
- **API:** `calculate_xpnsl(haplotypes_a, haplotypes_b, *, verbose=True)`
- **Input:** Two 0/1 matrices shaped (variants, haplotypes) with the same number of variants.
- **Output:** Returns `variants` and the `xpnsl` list.
- **Example code:**

```python
from dnakit.variants import calculate_xpnsl

haplotypes_a = [[0, 0, 1, 1], [0, 0, 1, 0], [1, 1, 0, 0], [0, 1, 0, 1], [1, 0, 1, 0], [1, 1, 1, 0]]
haplotypes_b = [[0, 1, 1, 0], [0, 0, 1, 1], [1, 1, 0, 0], [1, 0, 0, 1], [0, 1, 0, 1], [1, 1, 1, 0]]
result = calculate_xpnsl(haplotypes_a, haplotypes_b, verbose=False)
print(result["xpnsl"])
```

- **Example output:**

```text
[0.4054651081, 0.2876820725, 0.0, 0.5596157879, 0.3364722366, 0.1335313926]
```

## 31) Standardize a selection score within allele-count bins

- **Purpose:** Standardize selection scores such as iHS or nSL within allele-count bins.
- **Method:** Calls allel.standardize_by_allele_count with plot diagnostics off; `bins` must be strictly increasing.
- **API:** `standardize_score_by_allele_count(score, allele_counts, *, n_bins=None, bins=None, verbose=True)`
- **Input:** Scores `score`, per-site allele counts `allele_counts`, optional `n_bins` or `bins`.
- **Output:** Returns `values` and the `bins` used.
- **Example code:**

```python
from dnakit.variants import standardize_score_by_allele_count

score = [0.1, -0.2, 0.3, 0.5, -0.4, 0.2, 0.0, 1.0]
counts = [1, 2, 3, 4, 5, 6, 7, 8]
result = standardize_score_by_allele_count(score, counts, n_bins=2, verbose=False)
print(result["bins"], [round(v, 4) for v in result["values"][:3]])
```

- **Example output:**

```text
[1, 8] [-0.2163, -0.9581, 0.2781]
```

## 32) Tajima D difference between populations in moving windows

- **Purpose:** Compute the difference in Tajima's D between two populations in moving variant windows.
- **Method:** Calls allel.moving_delta_tajima_d; windows with too few sites return None.
- **API:** `calculate_moving_delta_tajimas_d(first_counts, second_counts, *, window_size, start=0, stop=None, step=None, verbose=True)`
- **Input:** Two count matrices and window parameters.
- **Output:** Returns `window_size`, `step` and `values`.
- **Example code:**

```python
from dnakit.variants import calculate_moving_delta_tajimas_d

first = [[4, 0], [3, 1], [2, 2], [1, 3], [4, 0], [0, 4], [3, 1], [2, 2], [1, 3], [4, 0], [2, 2], [3, 1]]
second = [[0, 4], [2, 2], [1, 3], [3, 1], [2, 2], [4, 0], [1, 3], [3, 1], [2, 2], [0, 4], [4, 0], [2, 2]]
result = calculate_moving_delta_tajimas_d(first, second, window_size=6, verbose=False)
print(result["values"])
```

- **Example output:**

```text
[-1.0, 1.0]
```

## 33) Fold a joint site frequency spectrum

- **Purpose:** Fold a two-population joint site frequency spectrum onto the minor-allele corner.
- **Method:** Calls allel.fold_joint_sfs; each dimension may hold at most sample size plus one classes.
- **API:** `fold_joint_site_frequency_spectrum(spectrum, *, first_sample_size, second_sample_size, verbose=True)`
- **Input:** Two-dimensional integer spectrum and two sample sizes.
- **Output:** Returns the `folded` matrix and both sample sizes.
- **Example code:**

```python
from dnakit.variants import fold_joint_site_frequency_spectrum

spectrum = [[0, 1, 2, 3], [4, 5, 6, 7], [8, 9, 10, 11]]
result = fold_joint_site_frequency_spectrum(spectrum, first_sample_size=2, second_sample_size=3, verbose=False)
print(result["folded"])
```

- **Example output:**

```text
[[14, 14]]
```

## 34) Scale a site frequency spectrum

- **Purpose:** Multiply class i of a spectrum by i, as used in i times xi_i statistics.
- **Method:** Calls allel.scale_sfs.
- **API:** `scale_site_frequency_spectrum(spectrum, *, verbose=True)`
- **Input:** One-dimensional integer spectrum.
- **Output:** Returns `classes` and `scaled`.
- **Example code:**

```python
from dnakit.variants import scale_site_frequency_spectrum

result = scale_site_frequency_spectrum([5, 3, 2, 1], verbose=False)
print(result["scaled"])
```

- **Example output:**

```text
[0, 3, 4, 3]
```

## 35) Scale a folded site frequency spectrum

- **Purpose:** Scale class k of a folded spectrum by k(n−k)/n.
- **Method:** Calls allel.scale_sfs_folded; the class count may not exceed the sample size plus one.
- **API:** `scale_folded_site_frequency_spectrum(spectrum, *, sample_size, verbose=True)`
- **Input:** One-dimensional folded spectrum and `sample_size` (number of chromosomes).
- **Output:** Returns `sample_size` and `scaled`.
- **Example code:**

```python
from dnakit.variants import scale_folded_site_frequency_spectrum

result = scale_folded_site_frequency_spectrum([5, 3, 2], sample_size=6, verbose=False)
print(result["scaled"])
```

- **Example output:**

```text
[0.0, 2.5, 2.6666666667]
```

## 36) Scale a joint site frequency spectrum

- **Purpose:** Multiply each cell (i, j) of a joint spectrum by i times j.
- **Method:** Calls allel.scale_joint_sfs.
- **API:** `scale_joint_site_frequency_spectrum(spectrum, *, verbose=True)`
- **Input:** Two-dimensional integer spectrum.
- **Output:** Returns `shape` and `scaled`.
- **Example code:**

```python
from dnakit.variants import scale_joint_site_frequency_spectrum

result = scale_joint_site_frequency_spectrum([[0, 1, 2], [3, 4, 5]], verbose=False)
print(result["scaled"])
```

- **Example output:**

```text
[[0, 0, 0], [0, 4, 10]]
```

## 37) Scale a folded joint site frequency spectrum

- **Purpose:** Scale each cell of a folded joint spectrum by i·j·(n1−i)·(n2−j).
- **Method:** Calls allel.scale_joint_sfs_folded.
- **API:** `scale_folded_joint_site_frequency_spectrum(spectrum, *, first_sample_size, second_sample_size, verbose=True)`
- **Input:** Two-dimensional spectrum and two sample sizes.
- **Output:** Returns both sample sizes and `scaled`.
- **Example code:**

```python
from dnakit.variants import scale_folded_joint_site_frequency_spectrum

result = scale_folded_joint_site_frequency_spectrum(
    [[0, 1, 2], [3, 4, 5]], first_sample_size=2, second_sample_size=3, verbose=False
)
print(result["scaled"])
```

- **Example output:**

```text
[[0, 0, 0], [0, 8, 10]]
```

## 38) Paint progeny haplotypes by parental transmission

- **Purpose:** Label each progeny allele by the parental haplotype it came from, or as a non-parental allele.
- **Method:** Calls allel.paint_transmission; the `codes` field gives the meaning of each code.
- **API:** `paint_parental_transmission(parent_haplotypes, progeny_haplotypes, *, verbose=True)`
- **Input:** Parental haplotypes (variants, 2) and progeny haplotypes (variants, progeny); -1 marks a missing call.
- **Output:** Returns the `painting` matrix and the `codes` table.
- **Example code:**

```python
from dnakit.variants import paint_parental_transmission

parent = [[0, 1], [1, 1]]
progeny = [[0, 1, 0], [1, 0, 1]]
result = paint_parental_transmission(parent, progeny, verbose=False)
print(result["painting"])
```

- **Example output:**

```text
[[1, 2, 1], [4, 5, 4]]
```

## 39) Summarize the header of a local VCF

- **Purpose:** Read sample names, INFO and FORMAT fields, filters and contig lengths from a local VCF header.
- **Method:** Calls allel.read_vcf_headers and parses the ##contig lines.
- **API:** `read_vcf_header_summary(path, *, verbose=True)`
- **Input:** Local VCF `path`.
- **Output:** Returns `samples`, `info_fields`, `format_fields`, `filters` and `contigs`.
- **Example code:**

```python
from pathlib import Path
from tempfile import TemporaryDirectory

from dnakit.variants import read_vcf_header_summary

ROWS = [
    "##fileformat=VCFv4.2",
    "##contig=<ID=chr1,length=1000>",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3",
    "chr1\t10\trs1\tA\tG\t50\tPASS\tDP=10\tGT\t0|1\t0/0\t1|1",
    "chr1\t20\trs2\tC\tT\t50\tPASS\tDP=12\tGT\t0/1\t1/1\t0/0",
    "chr1\t35\trs3\tG\tA\t50\tPASS\tDP=8\tGT\t1/1\t0/1\t0|0",
]
with TemporaryDirectory() as folder:
    source = Path(folder) / "small.vcf"
    source.write_text("\n".join(ROWS) + "\n", encoding="utf-8")
    result = read_vcf_header_summary(str(source), verbose=False)
print(result["samples"], result["info_fields"], result["contigs"])
```

- **Example output:**

```text
['S1', 'S2', 'S3'] ['DP'] [{'id': 'chr1', 'length': 1000}]
```

## 40) Convert a local VCF to an npz array archive

- **Purpose:** Read a local VCF into NumPy arrays and save them as a new compressed .npz file.
- **Method:** Calls allel.vcf_to_npz with overwrite=False; the output must be a new .npz file.
- **API:** `vcf_to_npz_archive(source, output, *, fields=None, verbose=True)`
- **Input:** Source VCF, new `.npz` output path, optional `fields` list.
- **Output:** Returns the `output` path and each array's shape in `arrays`.
- **Example code:**

```python
from pathlib import Path
from tempfile import TemporaryDirectory

from dnakit.variants import vcf_to_npz_archive

ROWS = [
    "##fileformat=VCFv4.2",
    "##contig=<ID=chr1,length=1000>",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3",
    "chr1\t10\trs1\tA\tG\t50\tPASS\tDP=10\tGT\t0|1\t0/0\t1|1",
    "chr1\t20\trs2\tC\tT\t50\tPASS\tDP=12\tGT\t0/1\t1/1\t0/0",
    "chr1\t35\trs3\tG\tA\t50\tPASS\tDP=8\tGT\t1/1\t0/1\t0|0",
]
with TemporaryDirectory() as folder:
    source = Path(folder) / "small.vcf"
    source.write_text("\n".join(ROWS) + "\n", encoding="utf-8")
    result = vcf_to_npz_archive(
        str(source), str(Path(folder) / "small.npz"),
        fields=["variants/POS", "calldata/GT"], verbose=False,
    )
print(result["arrays"])
```

- **Example output:**

```text
{'calldata/GT': [3, 3, 2], 'variants/POS': [3]}
```

## 41) Convert a local VCF to an HDF5 archive

- **Purpose:** Read a local VCF into arrays and write them to a new HDF5 file.
- **Method:** Calls allel.vcf_to_hdf5 (needs h5py); dataset shapes are read back after writing.
- **API:** `vcf_to_hdf5_archive(source, output, *, fields=None, verbose=True)`
- **Input:** Source VCF, new output path and optional `fields` list.
- **Output:** Returns the `output` path and `datasets` (dataset path to shape).
- **Example code:**

```python
from pathlib import Path
from tempfile import TemporaryDirectory

from dnakit.variants import vcf_to_hdf5_archive

ROWS = [
    "##fileformat=VCFv4.2",
    "##contig=<ID=chr1,length=1000>",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3",
    "chr1\t10\trs1\tA\tG\t50\tPASS\tDP=10\tGT\t0|1\t0/0\t1|1",
    "chr1\t20\trs2\tC\tT\t50\tPASS\tDP=12\tGT\t0/1\t1/1\t0/0",
    "chr1\t35\trs3\tG\tA\t50\tPASS\tDP=8\tGT\t1/1\t0/1\t0|0",
]
with TemporaryDirectory() as folder:
    source = Path(folder) / "small.vcf"
    source.write_text("\n".join(ROWS) + "\n", encoding="utf-8")
    result = vcf_to_hdf5_archive(
        str(source), str(Path(folder) / "small.h5"),
        fields=["variants/POS", "calldata/GT"], verbose=False,
    )
print(result["datasets"])
```

- **Example output:**

```text
{'calldata/GT': [3, 3, 2], 'variants/POS': [3]}
```

## Functions not wrapped

The functions below are not wrapped here. The reason for each is recorded in the `reason` field of the ledger.

- **Duplicates or already covered:** 18 in total. varcode `compare_samples` and `query_in_cis`; scikit-allel `average_patterson_fst`, `average_patterson_f3`, `fold_sfs` (covered by the dadi-based population wrapper) and `write_fasta` (covered by `dnakit.io.write`).
- **Family variants:** 32 in total. DataFrame, CSV and recarray readers, and `allel.write_vcf`; their record-level forms are covered by `dnakit.variants.read_vcf` or `rewrite_vcf_records`.
- **Needs pyensembl transcripts, effects or ExACTO files:** 7 varcode functions; `allel.vcf_to_zarr` needs zarr 2 (this environment has zarr 3.4.1), 1 in total.
- **Internal helpers:** 712 in total, covering registries and global configuration, grouping and sort keys, CSV metadata, window and index constructors, generic moving statistics, HDF5 storage helpers, and functions that only produce intermediate objects.
