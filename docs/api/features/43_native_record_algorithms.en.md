# Native FASTA and variant-record algorithms

These functions operate only on explicit FASTA/FASTQ or VCF-like records. They do not align reads, call variants, normalize against a reference genome, or query annotation databases. Coordinate and genotype rules are stated in the method and provenance so the steps can be reproduced.

Method citations can use the [Biopython sequence-processing paper](https://doi.org/10.1093/bioinformatics/btp163) and the [SAMtools/BCFtools/HTSlib paper](https://doi.org/10.1093/gigascience/giab008). These DOIs establish journal sources for the related software methods; they do not replace independent validation of versions, parameters, and data.

## 1) Calculate a FASTA GC profile

- **Purpose:** Calculate length, GC fraction, and GC percentage for each FASTA record.
- **Method:** Parse text or materialized records, normalize DNA case, and count each record independently.
- **API:** `dnakit.fasta_fastq.calculate_gc_profile(source[required])`.
- **Input:** FASTA text, a path, or `(id, sequence)` records.
- **Output:** A list of rows containing `id`, `length`, `gc_fraction`, and `gc_percent`.
- **Example code:**

```python
from dnakit.fasta_fastq import calculate_gc_profile

print(calculate_gc_profile(">a\nACGT\n>b\nGGAA\n", verbose=False))
```

- **Example output:**

```json
[{"id": "a", "length": 4, "gc_fraction": 0.5, "gc_percent": 50.0}, {"id": "b", "length": 4, "gc_fraction": 0.5, "gc_percent": 50.0}]
```

## 2) Summarize FASTA/FASTQ lengths

- **Purpose:** Summarize record count, minimum, maximum, mean, and length counts.
- **Method:** Use only sequence lengths; infer neither quality nor sequencing platform.
- **API:** `dnakit.fasta_fastq.calculate_length_distribution(source[required])`.
- **Input:** FASTA/FASTQ text, a path, or materialized records.
- **Output:** A mapping with `count`, `min`, `max`, `mean`, and `counts`.
- **Example code:**

```python
from dnakit.fasta_fastq import calculate_length_distribution

print(calculate_length_distribution(">a\nACGT\n>b\nGGAA\n", verbose=False))
```

- **Example output:**

```json
{"count": 2, "min": 4, "max": 4, "mean": 4.0, "counts": {"4": 2}}
```

## 3) Remove duplicate FASTA sequences

- **Purpose:** Remove duplicate normalized sequences while keeping the first record ID.
- **Method:** Parse and normalize each sequence, then use the sequence text as the key; reverse complements are not treated as duplicates.
- **API:** `dnakit.fasta_fastq.remove_duplicate_sequences(source[required])`.
- **Input:** FASTA text, a path, or `(id, sequence)` records.
- **Output:** The first `(id, sequence)` record for each unique sequence.
- **Example code:**

```python
from dnakit.fasta_fastq import remove_duplicate_sequences

print(remove_duplicate_sequences(">a\nACGT\n>dup\nACGT\n>b\nGGAA\n", verbose=False))
```

- **Example output:**

```json
[["a", "ACGT"], ["b", "GGAA"]]
```

## 4) Reverse-complement records

- **Purpose:** Generate a reverse-complement sequence for every FASTA record while preserving IDs and order.
- **Method:** Apply the DNA/IUPAC complement table after reversing each sequence.
- **API:** `dnakit.fasta_fastq.reverse_complement_records(source[required])`.
- **Input:** FASTA text, a path, or `(id, sequence)` records.
- **Output:** Reverse-complemented `(id, sequence)` records.
- **Example code:**

```python
from dnakit.fasta_fastq import reverse_complement_records

print(reverse_complement_records(">a\nACGT\n>b\nGGAA\n", verbose=False))
```

- **Example output:**

```json
[["a", "ACGT"], ["b", "TTCC"]]
```

## 5) Validate FASTA sequence IDs

- **Purpose:** Check whether FASTA IDs are non-empty and unique.
- **Method:** Collect header IDs in input order and report duplicates explicitly.
- **API:** `dnakit.fasta_fastq.validate_sequence_ids(source[required])`.
- **Input:** FASTA text, a path, or `(id, sequence)` records.
- **Output:** A mapping with `valid`, `ids`, and `duplicates`.
- **Example code:**

```python
from dnakit.fasta_fastq import validate_sequence_ids

print(validate_sequence_ids(">a\nACGT\n>b\nGGAA\n", verbose=False))
```

- **Example output:**

```json
{"valid": true, "ids": ["a", "b"], "duplicates": []}
```

## 6) Calculate VCF-like allele frequencies

- **Purpose:** Calculate ALT counts, called-allele counts, and frequencies from explicit genotypes.
- **Method:** Parse `0/1` and `1|1` style genotypes, ignore missing alleles, and perform no calling, imputation, or reference normalization.
- **API:** `dnakit.variants.calculate_allele_frequency(input, *, options=None, verbose=True)`.
- **Input:** `{"variants": [{"chrom", "pos", "ref", "alt", "genotypes"}]}` with VCF-like 1-based `pos`.
- **Output:** A completion envelope whose result contains per-variant and aggregate counts.
- **Example code:**

```python
from dnakit.variants import calculate_allele_frequency

result = calculate_allele_frequency(
    {"variants": [{"chrom": "chr1", "pos": 1, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]},
    verbose=False,
)
print(result["result"])
```

- **Example output:**

```json
{"variants": [{"chrom": "chr1", "pos": 1, "ref": "A", "alt": ["G"], "alt_count": 3, "called_alleles": 6, "frequency": 0.5}], "total_alt_alleles": 3, "total_called_alleles": 6, "aggregate_frequency": 0.5}
```

## 7) Calculate transition/transversion statistics

- **Purpose:** Classify explicit SNPs as transitions or transversions and calculate the Ti/Tv ratio.
- **Method:** Accept only single-base REF/ALT pairs; count other records in `excluded` without left-alignment or normalization.
- **API:** `dnakit.variants.calculate_transition_transversion(input, *, options=None, verbose=True)`.
- **Input:** `input.variants` records containing `chrom`, `pos`, `ref`, and `alt`.
- **Output:** A completion envelope containing classified records, transition/transversion counts, exclusions, and `ti_tv_ratio`.
- **Example code:**

```python
from dnakit.variants import calculate_transition_transversion

result = calculate_transition_transversion(
    {"variants": [{"chrom": "chr1", "pos": 1, "ref": "A", "alt": "G"}, {"chrom": "chr1", "pos": 2, "ref": "A", "alt": "C"}]},
    verbose=False,
)
print(result["result"])
```

- **Example output:**

```json
{"variants": [{"chrom": "chr1", "pos": 1, "ref": "A", "alt": "G", "classification": "transition"}, {"chrom": "chr1", "pos": 2, "ref": "A", "alt": "C", "classification": "transversion"}], "transitions": 1, "transversions": 1, "excluded": 0, "ti_tv_ratio": 1.0}
```

## 8) Filter variants by read depth

- **Purpose:** Keep VCF-like records using explicit `DP`, `depth`, or `read_depth` values.
- **Method:** Apply `min_depth`, `max_depth`, and `keep_missing` to each record's depth.
- **API:** `dnakit.variants.filter_variants_by_depth(input, *, options=None)`.
- **Input:** Records in `input.variants` with `chrom`, `pos`, `ref`, `alt`, and an explicit depth field.
- **Output:** A completion envelope with retained records, per-record metrics, counts, and predicate.
- **Example code:** `dnakit.variants.filter_variants_by_depth({"variants": [{"chrom": "chr1", "pos": 1, "ref": "A", "alt": "G", "DP": 30}]}, options={"min_depth": 20})`.
- **Example output:** `{"output_record_count": 1, "removed_record_count": 0}`.

## 9) Filter variants by genotype quality

- **Purpose:** Filter variant records using explicit GQ values.
- **Method:** Read GQ from the record or sample mappings and compare the minimum observed GQ with the threshold.
- **API:** `dnakit.variants.filter_variants_by_genotype_quality(input, *, options=None)`.
- **Input:** Records with `GQ`, `genotype_quality`, `genotype_qualities`, or sample GQ values.
- **Output:** A completion envelope with retained records, GQ metrics, counts, and predicate.
- **Example code:** `dnakit.variants.filter_variants_by_genotype_quality({"variants": [{"chrom": "chr1", "pos": 1, "ref": "A", "alt": "G", "GQ": 45}]}, options={"min_genotype_quality": 20})`.
- **Example output:** `{"output_record_count": 1, "removed_record_count": 0}`.
