# Dependency-free VCF record operations

This page documents eight deterministic record-level operations. They accept explicit VCF-like text or records and do not call a variant caller, reference genome, annotation database, or genotype imputation backend. Every result contains `status`, `operation`, `result`, and `provenance`.

## 1) Parse VCF records

- **Purpose:** Parse VCF text into serializable headers, samples, and records.
- **Method:** Validate the mandatory eight columns, 1-based `POS`, numeric `QUAL`, and `ALT`; parse INFO, FORMAT, samples, and GT without reference validation.
- **API:** `dnakit.variants.parse_vcf_records(input[required], options[optional])`.
- **Input:** VCF text, UTF-8 bytes, VCF lines, or a mapping containing `text`, `vcf`, or `content`.
- **Output:** `result` contains `headers`, `samples`, `records`, `record_count`, and `warnings`.
- **Example code:**

```python
from dnakit.variants import parse_vcf_records

text = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t10\trs1\tA\tG\t50\tPASS\tDP=12"
result = parse_vcf_records(text, verbose=False)
print(result["result"]["records"][0]["alt"])
```

- **Example output:** `['G']`

## 2) Split multiallelic records

- **Purpose:** Convert each record with multiple ALT alleles into biallelic records.
- **Method:** Emit one record per ALT and add `original_alt_index`; remap GT to the selected ALT by default.
- **API:** `dnakit.variants.split_multiallelic(input[required], options[optional])`.
- **Input:** A record sequence or mapping containing `variants`/`records`, with `chrom`, positive `pos`, `ref`, and `alt`.
- **Output:** Split records plus input/output counts.
- **Example code:**

```python
from dnakit.variants import split_multiallelic

result = split_multiallelic({"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": ["G", "T"], "genotypes": ["0/1", "1/2"]}]}, verbose=False)
print([row["alt"] for row in result["result"]["records"]])
```

- **Example output:** `['G', 'T']`

## 3) Count records and classes

- **Purpose:** Count records, ALT alleles, chromosome distribution, and SNP/MNP/indel classes.
- **Method:** Classify explicit `REF`/`ALT` lengths; multi-ALT records contribute one count per ALT.
- **API:** `dnakit.variants.count(input[required], options[optional])`。
- **Input:** A record sequence or mapping containing `variants`/`records`.
- **Output:** `record_count`, `alternate_allele_count`, `by_chromosome`, and `allele_type_counts`.
- **Example code:**

```python
from dnakit.variants import count

result = count({"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G"}]}, verbose=False)
print(result["result"]["record_count"])
```

- **Example output:** `1`

## 4) Compare exact variant identities

- **Purpose:** Compare two explicit variant collections by exact `CHROM/POS/REF/ALT` identity.
- **Method:** Expand ALT alleles and compute intersection, differences, and Jaccard index; no left-alignment or normalization.
- **API:** `dnakit.variants.compare(input[required], options[optional])`.
- **Input:** A mapping with `left` and `right` collections.
- **Output:** Counts, `only_left`, `only_right`, `shared`, and `jaccard`.
- **Example code:**

```python
from dnakit.variants import compare

result = compare({"left": {"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G"}]}, "right": {"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G"}, {"chrom": "chr1", "pos": 20, "ref": "C", "alt": "T"}]}}, verbose=False)
print(result["result"]["jaccard"])
```

- **Example output:** `0.5`

## 5) Calculate minor allele frequency

- **Purpose:** Calculate MAF from explicit genotype calls.
- **Method:** Count called reference and every declared ALT allele, including zero-count ALTs; missing calls are excluded and never imputed. If frequencies tie, the smallest allele index is reported for deterministic metadata while the MAF value is unchanged.
- **API:** `dnakit.variants.calculate_minor_allele_frequency(input[required], options[optional])`.
- **Input:** VCF-like records with explicit `genotypes`.
- **Output:** Per-site allele frequencies, minor allele, MAF, called allele count, and `mean_maf`.
- **Example code:**

```python
from dnakit.variants import calculate_minor_allele_frequency

result = calculate_minor_allele_frequency({"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}, verbose=False)
print(result["result"]["variants"][0]["maf"])
```

- **Example output:** `0.5`

## 6) Filter by QUAL

- **Purpose:** Keep records whose explicit QUAL meets a threshold.
- **Method:** Use `qual`, `QUAL`, or `quality`; missing QUAL is excluded unless `keep_missing=True`.
- **API:** `dnakit.variants.filter_variants_by_quality(input[required], options[required])`.
- **Input:** Records and `options.min_quality`.
- **Output:** Kept records, per-record metrics, predicate, and counts.
- **Example code:**

```python
from dnakit.variants import filter_variants_by_quality

result = filter_variants_by_quality({"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "qual": 50}, {"chrom": "chr1", "pos": 20, "ref": "C", "alt": "T", "qual": 10}]}, options={"min_quality": 20}, verbose=False)
print(result["result"]["output_record_count"])
```

- **Example output:** `1`

## 7) Filter by genotype missingness

- **Purpose:** Remove records whose missing-genotype fraction exceeds a threshold.
- **Method:** Compute missingness from `.` GT alleles or explicit `missing_count/sample_count`; no imputation.
- **API:** `dnakit.variants.filter_variants_by_missingness(input[required], options[required])`.
- **Input:** Records with GT or missing/sample counts and `options.max_missingness` in [0, 1].
- **Output:** Kept records, missingness metrics, predicate, and counts.
- **Example code:**

```python
from dnakit.variants import filter_variants_by_missingness

result = filter_variants_by_missingness({"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "./.", "0/0"]}]}, options={"max_missingness": 0.5}, verbose=False)
print(result["result"]["metrics"][0]["missingness"])
```

- **Example output:** `0.3333333333333333`

## 8) Filter by explicit allele frequency

- **Purpose:** Filter records by a selected ALT allele's explicit frequency range.
- **Method:** Calculate frequency from GT or `alt_count/allele_count`; select an ALT with `allele_index`.
- **API:** `dnakit.variants.filter_variants_by_allele_frequency(input[required], options[optional])`.
- **Input:** Variant records and optional `min_frequency`, `max_frequency`, and `allele_index`.
- **Output:** Kept records, per-site frequencies, predicate, and counts.
- **Example code:**

```python
from dnakit.variants import filter_variants_by_allele_frequency

result = filter_variants_by_allele_frequency({"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1"]}, {"chrom": "chr1", "pos": 20, "ref": "C", "alt": "T", "genotypes": ["0/0", "0/0"]}]}, options={"min_frequency": 0.25}, verbose=False)
print(result["result"]["output_record_count"])
```

- **Example output:** `1`

For publication, record the DNAKit version, input VCF version, coordinate convention, thresholds, and versions of any external caller or annotation software used upstream.
