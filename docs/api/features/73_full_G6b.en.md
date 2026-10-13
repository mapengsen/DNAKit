# G6b: PLINK relationship matrices, bed-reader genotype matrices and BGEN expectations and dosages

These functions in `dnakit.variants` are adapters over [pandas-plink](https://github.com/limix/pandas-plink), [bed-reader](https://github.com/fastlmm/bed-reader) and [bgen-reader](https://github.com/limix/bgen-reader). The relationship-matrix readers load GCTA GRM and PLINK 2 REL files; the bed-reader functions standardize and subset genotype matrices and compute a standardized relationship matrix from a `.bed` file; the bgen-reader functions compute allele expectations and dosages for unphased BGEN variants and create metafiles. Every function works only on explicit local files or lists. None of them uses the network, runs association tests or applies quality filters.

## Example data

The relationship-matrix examples use the same small matrix for three samples. It only demonstrates the interface and does not describe a real population. The BGEN examples write a two-variant, three-sample file at run time with `bgen_fixture.py`.

**Example file:** `study.grm`
```text
1	1	10	1.0
2	1	10	0.5
2	2	10	2.0
3	1	10	0.25
3	2	10	-0.5
3	3	10	1.5
```

**Example file:** `study.grm.id`
```text
fam1	ind1
fam1	ind2
fam2	ind3
```

**Example file:** `study.rel`
```text
1.0
0.5	2.0
0.25	-0.5	1.5
```

**Example file:** `study.rel.id`
```text
fam1	ind1
fam1	ind2
fam2	ind3
```

**Example file:** `bgen_fixture.py`
```python
import struct
import zlib


def write_bgen(path, samples, variants):
    """Write a diallelic, unphased BGEN 1.2 file with zlib-compressed 8-bit probabilities."""
    block = struct.pack("<I", len(samples)) + b"".join(
        struct.pack("<H", len(s)) + s.encode() for s in samples)
    block = struct.pack("<I", len(block) + 4) + block
    head = struct.pack("<III", 20, len(variants), len(samples)) + b"bgen"
    head += struct.pack("<I", 0x80000009)
    out = bytearray(struct.pack("<I", len(head) + len(block)) + head + block)
    for vid, rsid, chrom, pos, (a, b), probs in variants:
        item = b"".join(struct.pack("<H", len(t)) + t.encode() for t in (vid, rsid, chrom))
        item += struct.pack("<IH", pos, 2) + b"".join(struct.pack("<I", len(x)) + x.encode() for x in (a, b))
        data = struct.pack("<IHBB", len(samples), 2, 2, 2) + bytes(0x82 if p is None else 2 for p in probs)
        data += bytes([0, 8]) + bytes(v for p in probs for v in ([0, 0] if p is None else [round(x * 255) for x in p[:2]]))
        packed = zlib.compress(data)
        out += item + struct.pack("<II", len(packed) + 4, len(data)) + packed
    open(path, "wb").write(bytes(out))
```

## 1) `dnakit.variants.read_gcta_relationship_matrix()`

- **Purpose:** Reads a GCTA genetic relationship matrix (GRM) in text or binary form.
- **Method:** Calls pandas-plink's `read_grm`. Samples come from the `.id` file (tab-separated family and individual IDs). If the file lists more samples than `max_cells` allows, it fails before the matrix is read. When `id_path` is omitted, the name follows pandas-plink: a binary file `x.bin` maps to `x.id`; any other file loses a trailing `.gz` (if present) and gains `.id`, so `study.grm` maps to `study.grm.id`.
- **API:** `dnakit.variants.read_gcta_relationship_matrix(path, *, id_path=None, n_snps_path=None, max_cells=5000000, verbose=True)`.
- **Input:** GRM file path; `id_path` and `n_snps_path` are optional; `max_cells` limits the n × n matrix size.
- **Output:** A dictionary with `samples` (`family_id`, `sample_id`), `matrix` (an n × n list), `n_snps` (one count per stored pair, in file order) and `provenance`.
- **Example code:**

```python
from dnakit.variants import read_gcta_relationship_matrix

result = read_gcta_relationship_matrix("study.grm")
print(result["samples"][0])
print(result["matrix"][1])
print(result["n_snps"][:3])
```

- **Example output:**

```text
{'family_id': 'fam1', 'sample_id': 'ind1'}
[0.5, 2.0, -0.5]
[10, 10, 10]
```

## 2) `dnakit.variants.read_plink_relationship_matrix()`

- **Purpose:** Reads a PLINK 2 realized relationship matrix (`.rel` or `.rel.bin`).
- **Method:** Calls pandas-plink's `read_rel`. Text files are read row by row as the lower triangle and then symmetrized; binary files are read as an n × n float64 array. The `.id` name is the file name without `.bin`, `.zst` or `.gz`, plus `.id`. `.zst` files need the `zstandard` package.
- **API:** `dnakit.variants.read_plink_relationship_matrix(path, *, id_path=None, max_cells=5000000, verbose=True)`.
- **Input:** REL file path; `id_path` is optional; `max_cells` limits the n × n matrix size.
- **Output:** A dictionary with `samples`, `matrix` and `provenance`.
- **Example code:**

```python
from dnakit.variants import read_plink_relationship_matrix

result = read_plink_relationship_matrix("study.rel")
print(result["matrix"][2])
print(len(result["samples"]))
```

- **Example output:**

```text
[0.25, -0.5, 1.5]
3
```

## 3) `dnakit.variants.standardize_genotypes()`

- **Purpose:** Standardizes each variant column to mean 0 and population standard deviation 1.
- **Method:** Calls bed-reader's `standardize_f64`; with `dtype="float32"` it calls `standardize_f32`. The standard deviation is the population value (ddof 0). Missing values and zero-variance variants become 0, and a zero-variance variant reports its standard deviation as `None`. A variant with no observed genotype is an error, because bed-reader cannot compute its mean.
- **API:** `dnakit.variants.standardize_genotypes(genotypes, *, dtype="float64", max_cells=5000000, verbose=True)`.
- **Input:** A list of sample rows whose elements are 0, 1, 2 or `None` (missing); `dtype` is `"float64"` or `"float32"`.
- **Output:** A dictionary with `standardized` (the standardized matrix), `variant_statistics` (`mean` and `std` for each column) and `provenance`.
- **Example code:**

```python
from dnakit.variants import standardize_genotypes

result = standardize_genotypes([[0, 1], [1, 1], [2, 1], [None, 2]])
print([round(value, 4) for value in result["standardized"][0]])
print(result["variant_statistics"][0])
```

- **Example output:**

```text
[-1.2247, -0.5774]
{'index': 0, 'mean': 1.0, 'std': 0.816496580927726}
```

## 4) `dnakit.variants.calculate_genotype_relationship_matrix()`

- **Purpose:** Computes the standardized genotype relationship matrix `Z Z' / m` from a PLINK 1 `.bed` file.
- **Method:** Reads all genotypes with bed-reader's `open_bed`, skips variants with no observed genotype, standardizes with `standardize_f64`, and forms `Z Z'` with NumPy. `m` is the number of variants with a positive standard deviation; missing genotypes contribute 0. The matrix product is computed by NumPy, not by bed-reader.
- **API:** `dnakit.variants.calculate_genotype_relationship_matrix(path, *, max_cells=5000000, verbose=True)`.
- **Input:** Path to a `.bed` file; the `.bim` and `.fam` files with the same stem must exist; `max_cells` limits samples times variants.
- **Output:** A dictionary with `samples`, `matrix` (an n × n list), `variant_count`, `observed_variants`, `polymorphic_variants` and `provenance`.
- **Example code:**

```python
from dnakit.variants import calculate_genotype_relationship_matrix
from dnakit.variants.genotype_formats import write_plink_bed

write_plink_bed(
    "panel.bed",
    [[0, 1, 2], [1, 1, 0], [2, None, 1], [0, 0, 2]],
    sample_ids=["a", "b", "c", "d"],
    variant_ids=["s1", "s2", "s3"],
    chromosomes=["1", "1", "1"],
    positions=[1, 2, 3],
    allele_1=["A", "A", "A"],
    allele_2=["G", "G", "G"],
    overwrite=True,
)
result = calculate_genotype_relationship_matrix("panel.bed")
print(result["observed_variants"], result["polymorphic_variants"])
print([round(value, 4) for value in result["matrix"][0]])
```

- **Example output:**

```text
3 3
[0.7121, -0.3788, -0.5455, 0.2121]
```

## 5) `dnakit.variants.subset_genotype_matrix()`

- **Purpose:** Selects a sub-matrix of an in-memory genotype matrix by sample and variant index.
- **Method:** Calls bed-reader's `subset_f64_f64`; with `dtype="float32"` it calls `subset_f32_f32`. Rows and columns follow the order of the indices, indices may repeat, and missing values stay `None`.
- **API:** `dnakit.variants.subset_genotype_matrix(genotypes, *, sample_indices, variant_indices, dtype="float64", max_cells=5000000, verbose=True)`.
- **Input:** The genotype matrix; `sample_indices` and `variant_indices` are non-empty integer lists, and out-of-range values are an error.
- **Output:** A dictionary with `sample_indices`, `variant_indices`, `genotypes` and `provenance`.
- **Example code:**

```python
from dnakit.variants import subset_genotype_matrix

result = subset_genotype_matrix(
    [[0, 1, 2], [None, 1, 0], [2, 2, 2]],
    sample_indices=[0, 2],
    variant_indices=[1, 2],
)
print(result["genotypes"])
```

- **Example output:**

```text
[[1.0, 2.0], [2.0, 2.0]]
```

## 6) `dnakit.variants.read_bgen_allele_expectation()`

- **Purpose:** Computes, for each sample, the expected number of copies of each allele at unphased BGEN variants.
- **Method:** Opens the file with bgen-reader and calls `allele_expectation` for the selected variants. A sample's expectation is the sum of each genotype probability times its allele copy count, so AG contributes 1 to both A and G. Phased variants are rejected, and missing samples are `None`.
- **API:** `dnakit.variants.read_bgen_allele_expectation(path, *, variant_indices=None, max_cells=5000000, max_variants=100000, verbose=True)`.
- **Input:** BGEN file path; `variant_indices` defaults to all variants; `max_variants` and `max_cells` are explicit limits.
- **Output:** A dictionary with `samples`, `variants`, `expectations` (ordered by variant, sample and allele) and `provenance`.
- **Example code:**

```python
from bgen_fixture import write_bgen
from dnakit.variants import read_bgen_allele_expectation

write_bgen(
    "small.bgen",
    ["s1", "s2", "s3"],
    [
        ("v1", "rs1", "1", 100, ("A", "G"), [[51 / 255, 76 / 255, 128 / 255], [1, 0, 0], None]),
        ("v2", "rs2", "1", 200, ("C", "T"), [[0, 0, 1], [64 / 255, 128 / 255, 63 / 255], None]),
    ],
)
result = read_bgen_allele_expectation("small.bgen")
print(result["variants"][0]["rsid"], result["samples"])
print([round(value, 4) for value in result["expectations"][0][0]])
print(result["expectations"][0][2])
```

- **Example output:**

```text
rs1 ['s1', 's2', 's3']
[0.698, 1.302]
None
```

## 7) `dnakit.variants.read_bgen_dosage()`

- **Purpose:** Computes, for each sample, the dosage of one alternative allele: its expected number of copies.
- **Method:** Calls `compute_dosage` on the expectations from `read_bgen_allele_expectation`. When `alt_allele_index` is `None`, bgen-reader's default is used, which is the last allele column. The bgen-reader `compute_dosage` docstring says minor allele (MAF), but the code takes the last column.
- **API:** `dnakit.variants.read_bgen_dosage(path, *, variant_indices=None, alt_allele_index=None, max_cells=5000000, max_variants=100000, verbose=True)`.
- **Input:** The same as `read_bgen_allele_expectation`, plus `alt_allele_index` (a non-negative integer; an index beyond the allele count is an error).
- **Output:** A dictionary with `samples`, `variants` (each with `alt_allele_index` and `alt_allele`), `dosages` and `provenance`.
- **Example code:**

```python
from bgen_fixture import write_bgen
from dnakit.variants import read_bgen_dosage

write_bgen(
    "small.bgen",
    ["s1", "s2", "s3"],
    [
        ("v1", "rs1", "1", 100, ("A", "G"), [[51 / 255, 76 / 255, 128 / 255], [1, 0, 0], None]),
        ("v2", "rs2", "1", 200, ("C", "T"), [[0, 0, 1], [64 / 255, 128 / 255, 63 / 255], None]),
    ],
)
result = read_bgen_dosage("small.bgen", alt_allele_index=0)
print([variant["alt_allele"] for variant in result["variants"]])
print([round(value, 4) for value in result["dosages"][0][:2]])
```

- **Example output:**

```text
['A', 'C']
[0.698, 2.0]
```

## 8) `dnakit.variants.create_bgen_metafile()`

- **Purpose:** Creates a variant metadata file that makes later `read_bgen` calls faster.
- **Method:** Calls bgen-reader's `create_metafile`, which does not overwrite an existing file, so `overwrite=True` removes the old file first. It then reads the variant table back with `read_bgen(metafile_filepath=...)` to check the file and count the variants.
- **API:** `dnakit.variants.create_bgen_metafile(path, metafile_path, *, overwrite=False, verbose=True)`.
- **Input:** BGEN file path and metafile path; the metafile must differ from the BGEN file, and its directory must exist.
- **Output:** A dictionary with `path`, `metafile`, `variants`, `bytes` and `provenance`.
- **Example code:**

```python
from bgen_fixture import write_bgen
from dnakit.variants import create_bgen_metafile

write_bgen(
    "small.bgen",
    ["s1", "s2", "s3"],
    [
        ("v1", "rs1", "1", 100, ("A", "G"), [[51 / 255, 76 / 255, 128 / 255], [1, 0, 0], None]),
        ("v2", "rs2", "1", 200, ("C", "T"), [[0, 0, 1], [64 / 255, 128 / 255, 63 / 255], None]),
    ],
)
result = create_bgen_metafile("small.bgen", "small.metafile", overwrite=True)
print(result["variants"], result["metafile"].endswith("small.metafile"))
```

- **Example output:**

```text
2 True
```

## Installing the optional dependencies

```bash
python -m pip install "dnakit[bio-analysis]"      # numpy, needed by every function here
python -m pip install "dnakit[bio-pandas-plink]"  # pandas-plink and xarray, relationship matrices
python -m pip install "dnakit[bio-bed-reader]"    # bed-reader, standardization, subsets and GRM
python -m pip install "dnakit[bio-bgen-reader]"   # bgen-reader, BGEN expectations, dosages and metafiles
```

Project pages: [pandas-plink](https://github.com/limix/pandas-plink), [bed-reader](https://github.com/fastlmm/bed-reader), [bgen-reader](https://github.com/limix/bgen-reader).

## How these functions are verified

The test file is `others/tests/unit/test_full_G6b.py`, with 20 tests. They check:

- GRM and REL files, in binary and text form, read back to the hand-written values and match direct calls to pandas-plink's `read_grm` and `read_rel`.
- `standardize_genotypes` matches hand-computed means and population standard deviations; the float32 result agrees with float64 within 1e-6; `subset_genotype_matrix` matches a direct call to bed-reader's subset kernel.
- `calculate_genotype_relationship_matrix` matches an entry-by-entry manual formula, including missing and monomorphic variants.
- BGEN expectations and dosages match hand-computed k/255 probabilities and direct calls to `bgen_reader.read_bgen` with `allele_expectation` and `compute_dosage`; phased variants are rejected; metafiles can be read back with `read_bgen`.
- A missing optional package gives an install hint; values above `max_cells` or `max_variants`, invalid indices and out-of-range allele indices are rejected.
- All BGEN and PLINK fixtures are written in `tmp_path`, and no test uses the network.

## Current scope

This batch wraps the GRM and REL readers of pandas-plink; the standardization, subset and standardized-relationship-matrix functions of bed-reader; and the expectation, dosage and metafile functions of bgen-reader. Upstream functions that are not wrapped:

- Low-level Rust kernels (`read_f32`, `read_f64`, `read_i8`, `write_*`, `encode1_*`, `file_*_piece` and so on): reading and writing are already covered by `read_plink_bed` and `write_plink_bed`; the rest are internal linear-algebra or encoding steps. `subset_f32_f64` is a mixed-precision variant that is not exposed.
- Cloud access (`read_cloud_*`, `url_to_bytes`, `check_file_cloud`) and sample-data downloads (`sample_file`, `sample_url`, `example_filepath`): these need the network.
- Upstream functions that `genotype_formats.py` already calls (already wrapped): `read_plink`, `read_plink1_bin`, `write_plink1_bin`, `to_bed`, and the classes `open_bed` and `open_bgen`. Functions that duplicate an existing wrapper: `read_bgen`, `allele_frequency`, and the class `create_bed`.
- Test-suite, conftest, `tmp_path`, `get_num_threads` and the `Chunk` parameter class.

The examples import from `dnakit.variants`.
