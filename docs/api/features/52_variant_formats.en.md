# Variant nomenclature and genotype files: varcode, hgvs, PLINK, and BGEN

This group of functions in `dnakit.variants` is a set of adapters over optional packages for variant nomenclature and genotype files. The nomenclature part uses [varcode](https://github.com/openvax/varcode) to give the minimal representation and kind of a variant, and [hgvs](https://github.com/biocommons/hgvs) to parse, format, and syntax-check HGVS strings. The genotype-file part uses [bed-reader](https://github.com/fastlmm/bed-reader) and [pandas-plink](https://github.com/limix/pandas-plink) to read and write PLINK 1 `.bed`, `.bim`, and `.fam` files, and [bgen-reader](https://github.com/limix/bgen-reader) to read BGEN 1.2 metadata, probabilities, and allele frequencies. These functions do not compute gene-level effects, and they do not run HGVS operations that need reference sequences (sequence fetching, UTA projection, validation, and normalization). The dependencies are optional and are imported only when the matching function is called.

## Example data

The VCF below has 5 records and only demonstrates VCF reading. PLINK and BGEN are binary formats, so the example code writes its own files: bgen-reader has no writer, and its sample files must be downloaded from the network. The BGEN example file is therefore written by code in the BGEN 1.2 layout and read back by bgen-reader to check it.

**Example file:** `records.vcf`
```text
##fileformat=VCFv4.2
##contig=<ID=1,length=248956422>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	101	rs1	A	G	50	PASS	DP=5
1	200	.	ACG	A	40	PASS	DP=6
1	300	.	A	AGT	.	q10	DP=3
1	400	.	AC	TG	.	.	DP=3
1	500	.	A	G	.	PASS	DP=3
```

## 1) `dnakit.variants.describe_variant()`

- **Purpose:** Give the minimal representation and kind of a variant from its contig, 1-based start, REF, and ALT.
- **Method:** Call `varcode.Variant`. It removes the common prefix and common suffix of REF and ALT to give the `trimmed_*` coordinates. The kind is `snv`, `insertion`, `deletion`, `substitution` (equal-length change), or `complex`. Only SNVs get the transition flag. No reference sequence is read.
- **API:** `dnakit.variants.describe_variant(contig, start, ref, alt)`.
- **Input:** Contig name (no whitespace), 1-based start, REF, and ALT (only A, C, G, T, N, case-insensitive). REF and ALT cannot both be empty, and they cannot be equal.
- **Output:** A dictionary with `contig`, `start`, `ref`, `alt` (as given), `trimmed_start`, `trimmed_end`, `trimmed_ref`, `trimmed_alt`, `kind`, `is_transition`, and `provenance`.
- **Example code:**

```python
from dnakit.variants import describe_variant

result = describe_variant("1", 200, "ACG", "A")
print(result["kind"], result["trimmed_start"], result["trimmed_end"], repr(result["trimmed_ref"]), repr(result["trimmed_alt"]))
```

- **Example output:**

```text
deletion 201 202 'CG' ''
```

## 2) `dnakit.variants.load_vcf_variants()`

- **Purpose:** Read VCF records through varcode and give the minimal representation and kind of each record.
- **Method:** Call `varcode.load_vcf`. `distinct=False` keeps every record and does not merge identical sites. With `only_passing=True`, only records whose FILTER is PASS or missing are kept. `genome` is the reference name that varcode requires; this function does not read reference sequences. varcode truncates silently at `max_variants`, so this function requests one extra record and raises an error when the limit is exceeded.
- **API:** `dnakit.variants.load_vcf_variants(path, *, genome="GRCh38", only_passing=True, max_variants=100000)`.
- **Input:** VCF file path. `max_variants` is an integer from 1 to 1000000.
- **Output:** `count` and a `variants` list. Each item has the same fields as `describe_variant()`; `start`, `ref`, and `alt` keep the VCF input form.
- **Example code:**

```python
from dnakit.variants import load_vcf_variants

result = load_vcf_variants("records.vcf")
print(result["count"])
print([(item["start"], item["kind"]) for item in result["variants"]])
```

- **Example output:**

```text
4
[(101, 'snv'), (200, 'deletion'), (400, 'substitution'), (500, 'snv')]
```

## 3) `dnakit.variants.parse_hgvs_variant()`

- **Purpose:** Parse one HGVS string (c., g., n., m., r., or p.) and return its coordinates, edit type, and bases.
- **Method:** Call `hgvs.parser.Parser().parse_hgvs_variant`. Only syntax is parsed; no reference data is accessed. Syntax errors raise an error. An interval whose start is after its end (compared by base and offset) also raises an error, because hgvs does not check this.
- **API:** `dnakit.variants.parse_hgvs_variant(text)`.
- **Input:** One HGVS string without whitespace, at most 1000 characters.
- **Output:** `reference` (reference sequence ID), `gene`, `coordinate_type`, `start` and `end` (with `base` and `offset`; `offset` is None for p. variants), `edit_type`, `reference_bases`, `alternate_bases`, `formatted`, and `provenance`.
- **Example code:**

```python
from dnakit.variants import parse_hgvs_variant

result = parse_hgvs_variant("NM_004006.2:c.4375C>T")
print(result["coordinate_type"], result["start"]["base"], result["edit_type"], result["reference_bases"], result["alternate_bases"])
print(parse_hgvs_variant("NP_003997.1:p.Val600Glu")["formatted"])
```

- **Example output:**

```text
c 4375 sub C T
NP_003997.1:p.Val600Glu
```

## 4) `dnakit.variants.format_hgvs_variant()`

- **Purpose:** Format an HGVS string into hgvs's canonical text and report whether the text changed.
- **Method:** Parse the string and call hgvs `format()` with the default formatting settings. Formatting needs no reference data. For example, `c.123delA` is written as `c.123del`, because hgvs's default output does not write deleted bases.
- **API:** `dnakit.variants.format_hgvs_variant(text)`.
- **Input:** Same as `parse_hgvs_variant()`.
- **Output:** `input`, `formatted`, `changed`, and `provenance`.
- **Example code:**

```python
from dnakit.variants import format_hgvs_variant

result = format_hgvs_variant("NM_004006.2:c.123delA")
print(result["formatted"], result["changed"])
print(format_hgvs_variant("NM_004006.2:c.4375C>T")["changed"])
```

- **Example output:**

```text
NM_004006.2:c.123del True
False
```

## 5) `dnakit.variants.check_hgvs_syntax()`

- **Purpose:** Check the syntax of a batch of HGVS strings and report each entry as valid or invalid, with its canonical text or error.
- **Method:** Each string goes through the same parsing steps as `parse_hgvs_variant()`. Syntax errors and reversed intervals are recorded as invalid entries and do not stop the batch. Non-string entries raise an error.
- **API:** `dnakit.variants.check_hgvs_syntax(texts, *, max_variants=10000)`.
- **Input:** A list of strings. Its length cannot exceed `max_variants` (1 to 10000).
- **Output:** `count`, `valid`, `invalid`, and `items`. Each item has `input`, `valid`, `formatted`, and `error`.
- **Example code:**

```python
from dnakit.variants import check_hgvs_syntax

result = check_hgvs_syntax(["NM_004006.2:c.4375C>T", "NM_004006.2:c.123_122del", "not an hgvs string"])
print(result["count"], result["valid"], result["invalid"])
print([(item["input"], item["valid"]) for item in result["items"]])
```

- **Example output:**

```text
3 1 2
[('NM_004006.2:c.4375C>T', True), ('NM_004006.2:c.123_122del', False), ('not an hgvs string', False)]
```

## 6) `dnakit.variants.read_plink_bed()`

- **Purpose:** Read a PLINK 1 `.bed` genotype matrix, optionally selecting samples and variants by ID.
- **Method:** `bed_reader.open_bed` reads the same-named `.fam` and `.bim`, then genotypes are read by index. `count_allele="allele_1"` counts the first allele (A1), which matches bed-reader's `count_A1=True`. Missing genotypes are returned as `None`. Before reading, the record counts of `.fam` and `.bim` and the matrix cell count are checked against their limits.
- **API:** `dnakit.variants.read_plink_bed(path, *, samples=None, variants=None, count_allele="allele_1", max_records=100000, max_cells=5000000)`.
- **Input:** `.bed` file path; the same folder must hold same-named `.bim` and `.fam` files. `samples` and `variants` are lists of IDs, or None.
- **Output:** `samples`, `variants` (each with `id`, `chromosome`, `position`, `allele_1`, `allele_2`), `genotypes` (rows are samples, columns are variants, values 0, 1, 2, or None), `count_allele`, and `provenance`.
- **Example code:**

```python
from dnakit.variants import read_plink_bed, write_plink_bed

write_plink_bed(
    "panel.bed",
    [[0, 1, 2], [2, None, 1], [1, 1, 0]],
    sample_ids=["S1", "S2", "S3"],
    variant_ids=["rs1", "rs2", "rs3"],
    chromosomes=["1", "1", "2"],
    positions=[100, 200, 300],
    allele_1=["A", "C", "G"],
    allele_2=["G", "T", "A"],
    overwrite=True,
)
result = read_plink_bed("panel.bed")
print(result["samples"], [variant["id"] for variant in result["variants"]])
print(result["genotypes"])
```

- **Example output:**

```text
['S1', 'S2', 'S3'] ['rs1', 'rs2', 'rs3']
[[0, 1, 2], [2, None, 1], [1, 1, 0]]
```

## 7) `dnakit.variants.write_plink_bed()`

- **Purpose:** Write a genotype matrix as PLINK 1 `.bed`, `.bim`, and `.fam` files.
- **Method:** Call `bed_reader.to_bed`. `count_allele` says which allele the values count. Existing files are not overwritten unless `overwrite=True`.
- **API:** `dnakit.variants.write_plink_bed(path, genotypes, *, sample_ids, variant_ids, chromosomes, positions, allele_1, allele_2, count_allele="allele_1", overwrite=False, max_cells=5000000)`.
- **Input:** `genotypes` is a matrix with sample rows and variant columns, values 0, 1, 2, or None. The other lists must match the matrix dimensions; sample IDs and variant IDs must be unique.
- **Output:** `path`, `files`, `samples`, `variants`, `count_allele`, and `provenance`.
- **Example code:**

```python
from pathlib import Path

from dnakit.variants import write_plink_bed

result = write_plink_bed(
    "written.bed",
    [[0, 1], [2, None]],
    sample_ids=["S1", "S2"],
    variant_ids=["rs1", "rs2"],
    chromosomes=["1", "1"],
    positions=[10, 20],
    allele_1=["A", "C"],
    allele_2=["G", "T"],
    overwrite=True,
)
print(result["operation"], result["samples"], result["variants"], [Path(name).name for name in result["files"]])
```

- **Example output:**

```text
write_plink_bed 2 2 ['written.bed', 'written.bim', 'written.fam']
```

## 8) `dnakit.variants.summarize_plink_bed()`

- **Purpose:** For each PLINK variant, count called and missing genotypes and the frequency of the counted allele.
- **Method:** The full matrix is read (limited by `max_cells`). The frequency is `counted allele copies / (2 × called genotypes)`, because PLINK genotypes are diploid.
- **API:** `dnakit.variants.summarize_plink_bed(path, *, count_allele="allele_1", max_records=100000, max_cells=5000000)`.
- **Input:** `.bed` file path; same as `read_plink_bed()`.
- **Output:** `samples`, `variants` (each with `id`, `called`, `missing`, `allele_frequency`), `missing_genotypes`, `count_allele`, and `provenance`.
- **Example code:**

```python
from dnakit.variants import summarize_plink_bed, write_plink_bed

write_plink_bed(
    "summary.bed",
    [[0, 2, 1], [0, 2, None], [0, None, 2]],
    sample_ids=["S1", "S2", "S3"],
    variant_ids=["rs1", "rs2", "rs3"],
    chromosomes=["1", "1", "2"],
    positions=[100, 200, 300],
    allele_1=["A", "C", "G"],
    allele_2=["G", "T", "A"],
    overwrite=True,
)
result = summarize_plink_bed("summary.bed")
print([(item["id"], item["called"], item["missing"], item["allele_frequency"]) for item in result["variants"]])
print(result["missing_genotypes"])
```

- **Example output:**

```text
[('rs1', 3, 0, 0.0), ('rs2', 2, 1, 1.0), ('rs3', 2, 1, 0.75)]
2
```

## 9) `dnakit.variants.read_plink_metadata()`

- **Purpose:** Read the `.fam` sample table and the `.bim` variant table without reading genotypes.
- **Method:** `pandas_plink.read_plink` reads both tables; the genotype part is never accessed. The record limit is checked before reading.
- **API:** `dnakit.variants.read_plink_metadata(path, *, max_records=100000)`.
- **Input:** `.bed` file path (same-named `.bim` and `.fam` must exist in the same folder).
- **Output:** `samples` (with `family_id`, `sample_id`, `father_id`, `mother_id`, `sex`, `phenotype`), `variants` (with `chromosome`, `id`, `genetic_position`, `position`, `allele_1`, `allele_2`), and `provenance`.
- **Example code:**

```python
from dnakit.variants import read_plink_metadata

result = read_plink_metadata("panel.bed")
print(result["samples"][0]["sample_id"], len(result["samples"]))
print(result["variants"][2])
```

- **Example output:**

```text
S1 3
{'chromosome': '2', 'id': 'rs3', 'genetic_position': 0.0, 'position': 300, 'allele_1': 'G', 'allele_2': 'A'}
```

## 10) `dnakit.variants.read_plink_genotypes()`

- **Purpose:** Read a PLINK 1 genotype matrix through pandas-plink and return the same structure as `read_plink_bed()`.
- **Method:** `pandas_plink.read_plink1_bin` returns the counted allele for `ref`: `allele_1` maps to `ref="a0"` and `allele_2` maps to `ref="a1"`. The cell limit is checked before reading. Selecting a subset is not supported by this function.
- **API:** `dnakit.variants.read_plink_genotypes(path, *, count_allele="allele_1", max_records=100000, max_cells=5000000)`.
- **Input:** `.bed` file path; same-named `.bim` and `.fam` files must be in the same folder.
- **Output:** `samples`, `variants`, `genotypes`, `count_allele`, and `provenance`.
- **Example code:**

```python
from dnakit.variants import read_plink_genotypes

print(read_plink_genotypes("panel.bed")["genotypes"])
print(read_plink_genotypes("panel.bed", count_allele="allele_2")["genotypes"])
```

- **Example output:**

```text
[[0, 1, 2], [2, None, 1], [1, 1, 0]]
[[2, 1, 0], [0, None, 1], [1, 1, 2]]
```

## 11) `dnakit.variants.write_plink_genotypes()`

- **Purpose:** Write PLINK 1 files through the pandas-plink writer.
- **Method:** The matrix is wrapped in an xarray `DataArray` and passed to `pandas_plink.write_plink1_bin`. That writer stores counts of the second allele (a1), so with `count_allele="allele_1"` each value is first converted to `2 - value`.
- **API:** `dnakit.variants.write_plink_genotypes(path, genotypes, *, sample_ids, variant_ids, chromosomes, positions, allele_1, allele_2, count_allele="allele_1", overwrite=False, max_cells=5000000)`.
- **Input:** Same as `write_plink_bed()`.
- **Output:** Same as `write_plink_bed()`: `path`, `files`, `samples`, `variants`, `count_allele`, and `provenance`.
- **Example code:**

```python
from dnakit.variants import read_plink_bed, write_plink_genotypes

write_plink_genotypes(
    "pandas.bed",
    [[0, 1, 2], [2, None, 1], [1, 1, 0]],
    sample_ids=["S1", "S2", "S3"],
    variant_ids=["rs1", "rs2", "rs3"],
    chromosomes=["1", "1", "2"],
    positions=[100, 200, 300],
    allele_1=["A", "C", "G"],
    allele_2=["G", "T", "A"],
    overwrite=True,
)
print(read_plink_bed("pandas.bed")["genotypes"])
```

- **Example output:**

```text
[[0, 1, 2], [2, None, 1], [1, 1, 0]]
```

## 12) `dnakit.variants.read_bgen_variants()`

- **Purpose:** Read the sample count and variant metadata of a BGEN file, without reading probabilities.
- **Method:** `bgen_reader.open_bgen` opens the file, and variant IDs, rsIDs, chromosomes, positions, alleles, phasing, and genotype combination counts are read. The first time a file is opened, bgen-reader writes a `.metadata2.mmm` cache file next to it. `max_variants` limits the number of variants. This function was tested only on zlib-compressed BGEN 1.2 files with layout 2; uncompressed files were not tested.
- **API:** `dnakit.variants.read_bgen_variants(path, *, max_variants=100000)`.
- **Input:** BGEN 1.2 file path.
- **Output:** `sample_count`, `variants` (each with `id`, `rsid`, `chromosome`, `position`, `alleles`, `phased`, `genotype_combinations`), and `provenance`.
- **Example code:**

```python
import struct
import zlib

from dnakit.variants import read_bgen_variants


def write_bgen(path, samples, variants):
    n = len(samples)
    flags = 0x80000000 | (2 << 2) | 1
    header = struct.pack("<III", 20, len(variants), n) + b"bgen" + struct.pack("<I", flags)
    names = b"".join(struct.pack("<H", len(s)) + s.encode() for s in samples)
    body = b""
    for ident, rsid, chrom, pos, alleles, probs in variants:
        for text in (ident, rsid, chrom):
            body += struct.pack("<H", len(text)) + text.encode()
        body += struct.pack("<IH", pos, len(alleles))
        body += b"".join(struct.pack("<I", len(a)) + a.encode() for a in alleles)
        raw = struct.pack("<IHBB", n, 2, 2, 2) + bytes([2] * n) + bytes([0, 8])
        raw += bytes(value for row in probs for value in row)
        packed = zlib.compress(raw)
        body += struct.pack("<II", 4 + len(packed), len(raw)) + packed
    sample_block = struct.pack("<II", 8 + len(names), n) + names
    with open(path, "wb") as handle:
        handle.write(struct.pack("<I", len(header) + len(sample_block)) + header + sample_block + body)


write_bgen(
    "two.bgen",
    ["S1", "S2", "S3"],
    [
        ("v1", "rs1", "1", 100, ["A", "G"], [[255, 0], [0, 255], [128, 127]]),
        ("v2", "rs2", "1", 200, ["C", "T"], [[0, 0], [255, 0], [0, 0]]),
    ],
)
result = read_bgen_variants("two.bgen")
print(result["sample_count"], [(v["rsid"], v["position"], v["alleles"]) for v in result["variants"]])
```

- **Example output:**

```text
3 [('rs1', 100, ['A', 'G']), ('rs2', 200, ['C', 'T'])]
```

## 13) `dnakit.variants.read_bgen_probabilities()`

- **Purpose:** Read genotype probabilities for selected variants.
- **Method:** `open_bgen(...).read(..., return_missings=True)` returns probabilities and missing flags; missing samples are returned as `None`. The output is ordered by variant, then sample, then genotype combination. For diploid biallelic sites the combination order is AA, AB, BB.
- **API:** `dnakit.variants.read_bgen_probabilities(path, *, variant_indices=None, max_cells=5000000, max_variants=100000)`.
- **Input:** `variant_indices` is a list of 0-based variant indices; None reads all variants.
- **Output:** `samples`, `variants`, `probabilities` (variant → sample → list of probabilities, or None when missing), and `provenance`.
- **Example code:**

```python
from dnakit.variants import read_bgen_probabilities

result = read_bgen_probabilities("two.bgen", variant_indices=[1])
print(result["variants"][0]["rsid"], result["samples"])
print(result["probabilities"][0])
```

- **Example output:**

```text
rs2 ['S1', 'S2', 'S3']
[[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]
```

## 14) `dnakit.variants.calculate_bgen_allele_frequency()`

- **Purpose:** Compute allele frequencies for unphased BGEN variants.
- **Method:** `allele_expectation` gives each sample's expected allele copies, derived from the probabilities. After missing samples are removed, each allele's frequency is its expected copies divided by the total expected copies. bgen-reader's own `allele_frequency` returns copy sums divided by the allele count (the values sum to the sample count), so this function does not use it. Phased variants raise an error.
- **API:** `dnakit.variants.calculate_bgen_allele_frequency(path, *, variant_indices=None, max_cells=5000000, max_variants=100000)`.
- **Input:** Same as `read_bgen_probabilities()`.
- **Output:** `variants` (each with `id`, `rsid`, `chromosome`, `position`, `alleles`, `called_samples`, `allele_copies`, `allele_frequencies`) and `provenance`.
- **Example code:**

```python
from dnakit.variants import calculate_bgen_allele_frequency

result = calculate_bgen_allele_frequency("two.bgen")
print([(v["rsid"], [round(f, 4) for f in v["allele_frequencies"]]) for v in result["variants"]])
```

- **Example output:**

```text
[('rs1', [0.7503, 0.2497]), ('rs2', [0.3333, 0.6667])]
```

## How these functions are verified

- The test file is `others/tests/unit/test_variant_formats_backends.py`, with 32 tests. Most run against the real installed packages; the missing-package tests simulate an absent package with `sys.modules`.
- PLINK: bed-reader and pandas-plink results are cross-checked against each other for both counting conventions. They are also compared with direct calls to `bed_reader.open_bed` and `pandas_plink.read_plink`. Summary frequencies are compared with hand calculations.
- BGEN: the test writes zlib-compressed BGEN 1.2 layout 2 files following the byte structure of bgen-reader's sample files, and bgen-reader reads them back. Probabilities, missing samples, and allele frequencies are compared with hand calculations; frequencies are also compared with direct `allele_expectation` calls.
- varcode: `describe_variant()` trimming is compared with `varcode.Variant`; `load_vcf_variants()` is compared with `varcode.load_vcf`.
- hgvs: parse and format results are compared with direct calls to `hgvs.parser.Parser()`. Syntax errors, reversed intervals, non-string entries, and batches above the limit are rejected.
- When an optional package is missing, each module raises an install hint (`bio-varcode`, `bio-hgvs`, `bio-bed-reader`, `bio-pandas-plink`, `bio-bgen-reader`).

## Current scope

- Wrapped: varcode minimal representation, variant kind, and VCF reading; hgvs syntax parsing, formatting, and batch syntax checks; PLINK 1 reading, writing, and summaries with bed-reader and pandas-plink; BGEN metadata, probabilities, and allele frequencies with bgen-reader.
- Not wrapped: varcode gene-level effects (`Variant.effects()` and related functions). They need pyensembl to download Ensembl reference data, which was not run here.
- Not wrapped: hgvs reference-based functions (sequence fetching, UTA projection, validator, normalizer). They need seqrepo and a psycopg2 database connection; psycopg2 and seqrepo are missing in this environment, so these functions were not run and are not tested.
- Not wrapped: pandas-plink `read_grm` and `read_rel`. They read relatedness matrices, which are not DNA-level or gene-level operations.
- Not tested: uncompressed BGEN files and BGEN files with `allow_complex=True`.
- Not wrapped: sgkit. It is not installed in this environment, so no sgkit functions were added.
