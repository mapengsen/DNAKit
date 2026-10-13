# Publication and variant-statistics evidence

This page records how DNAKit cites third-party software and documents the first dependency-free variant statistics. Being citable in a paper does not grant permission to redistribute source code, model weights, databases or hosted services.

## 1) Query a software citation record

- **Purpose:** Query the paper, license and integration boundary for tools such as Biopython, SeqKit, scikit-bio, pydna, DNA Chisel, seqme, Edlib, Cutadapt and MultiQC.
- **Method:** Read the frozen offline evidence registry bundled with DNAKit.
- **API:** `dnakit.integrations.describe_publication_record(name[required])`.
- **Input:** A registered software name or alias.
- **Output:** Citation, license, source URLs, citation status and `reuse_route`; this is not a claim that DNAKit re-exports the full upstream API.
- **Example code:**

```python
from dnakit.integrations import describe_publication_record

record = describe_publication_record("Biopython", verbose=False)
print(record["doi"])
print(record["reuse_route"])
```

- **Example output:**

```text
10.1093/bioinformatics/btp163
optional_pip_adapter
```

## 2) Check manuscript evidence

- **Purpose:** Check whether software, backend and database versions, access date, citation and license evidence are ready for a Methods section.
- **Method:** Check citation evidence, license evidence and runtime-version evidence separately.
- **API:** `dnakit.integrations.publication_readiness(name[required], version[optional], backend_version[optional], database_version[optional], accessed[optional])`.
- **Input:** The software name and the exact version used.
- **Output:** A status, individual checks, citation, license, sources and reuse boundary.
- **Example code:**

```python
from dnakit.integrations import publication_readiness

result = publication_readiness("SeqKit", version="2.10.0", accessed="2026-10-12", verbose=False)
print(result["status"])
```

- **Example output:**

```text
ready_for_manuscript_record
```

## 3) Build a software citation

- **Purpose:** Build a structured citation record for a Methods section or supplement.
- **Method:** Combine software name, version, DOI or project URL, access date, license and backend versions.
- **API:** `dnakit.integrations.build_software_citation(name[required], version[required], accessed[optional], backend_version[optional], database_version[optional])`.
- **Input:** The exact package or executable version used in the experiment.
- **Output:** Citation text, DOI, license, source URLs, reuse route and evidence checks.
- **Example code:**

```python
from dnakit.integrations import build_software_citation

record = build_software_citation("scikit-bio", "0.7.1", accessed="2026-10-12", verbose=False)
print(record["citation"])
```

- **Example output:**

```text
scikit-bio version 0.7.1; Aton et al. (2026), Scikit-bio: a fundamental Python library for biological omic data analysis. DOI: 10.1038/s41592-025-02981-z.
```

## 4) Calculate allele frequency

- **Purpose:** Calculate alternate-allele frequency from explicit VCF-like genotypes.
- **Method:** Parse called alleles from `0/0`, `0/1`, `1/1`, `./.` or equivalent integer arrays; no genotype calling or imputation is performed.
- **API:** `dnakit.variants.calculate_allele_frequency(input[required], options[optional])`.
- **Input:** A mapping containing `variants`, or a sequence of mappings with `chrom`, `pos`, `ref`, `alt` and `genotypes`.
- **Output:** Per-site alternate allele counts, called allele counts and frequencies, plus an aggregate frequency.
- **Example code:**

```python
from dnakit.variants import calculate_allele_frequency

result = calculate_allele_frequency({"variants": [{
    "chrom": "chr1", "pos": 10, "ref": "A", "alt": "G",
    "genotypes": ["0/1", "1/1", "0/0"],
}]})
print(result["result"]["aggregate_frequency"])
```

- **Example output:**

```text
0.5
```

## 5) Calculate observed heterozygosity

- **Purpose:** Calculate observed heterozygosity for explicit diploid genotypes.
- **Method:** Count non-missing samples with two different alleles.
- **API:** `dnakit.variants.calculate_heterozygosity(input[required], options[optional])`.
- **Input:** VCF-like records with diploid `genotypes`.
- **Output:** Per-site and aggregate observed heterozygosity.
- **Example code:**

```python
from dnakit.variants import calculate_heterozygosity

result = calculate_heterozygosity({"variants": [{
    "chrom": "chr1", "pos": 10, "ref": "A", "alt": "G",
    "genotypes": ["0/1", "1/1", "0/0"],
}]})
print(result["result"]["observed_heterozygosity"])
```

- **Example output:**

```text
0.3333333333333333
```

## 6) Calculate missingness

- **Purpose:** Calculate the fraction of samples with missing genotypes.
- **Method:** Genotypes containing `.` are missing; no imputation is performed.
- **API:** `dnakit.variants.calculate_missingness(input[required], options[optional])`.
- **Input:** VCF-like variant records and sample genotypes.
- **Output:** Missing samples, total samples and missingness.
- **Example code:**

```python
from dnakit.variants import calculate_missingness

result = calculate_missingness({"variants": [{
    "chrom": "chr1", "pos": 10, "ref": "A", "alt": "G",
    "genotypes": ["0/1", "./.", "0/0"],
}]})
print(result["result"]["missingness"])
```

- **Example output:**

```text
0.3333333333333333
```

## 7) Calculate transition/transversion ratio

- **Purpose:** Count SNP transitions, transversions and their ratio.
- **Method:** A↔G and C↔T are transitions; other single-base substitutions are transversions; indels are excluded.
- **API:** `dnakit.variants.calculate_transition_transversion(input[required], options[optional])`.
- **Input:** Variant records with `ref` and `alt`.
- **Output:** Counts, excluded records and `ti_tv_ratio`.
- **Example code:**

```python
from dnakit.variants import calculate_transition_transversion

result = calculate_transition_transversion({"variants": [
    {"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G"},
    {"chrom": "chr1", "pos": 20, "ref": "A", "alt": "C"},
]})
print(result["result"]["ti_tv_ratio"])
```

- **Example output:**

```text
1.0
```

## 8) Calculate variant density

- **Purpose:** Calculate variant count and density for an explicitly sized region.
- **Method:** Divide the number of input records by an explicit region length in bases.
- **API:** `dnakit.variants.calculate_variant_density(input[required], options[required])`.
- **Input:** Variant records and `options["region_length"]`.
- **Output:** Region length, variant count, density per base, density per megabase and chromosome counts.
- **Example code:**

```python
from dnakit.variants import calculate_variant_density

result = calculate_variant_density({"variants": [
    {"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G"},
    {"chrom": "chr1", "pos": 20, "ref": "C", "alt": "T"},
]}, options={"region_length": 1_000_000})
print(result["result"]["density_per_megabase"])
```

- **Example output:**

```text
2.0
```

## 9) Manuscript citation boundary

- **Purpose:** Keep paper citation eligibility separate from permission to redistribute code, models, databases or hosted services.
- **Method:** Prefer the registry DOI and record the exact software, Python/CLI backend, model or database version and access date.
- **Input:** The software and resources actually run in the Methods section.
- **Output:** `verified_doi` means DOI evidence is recorded; `preprint_or_project` requires explicit preprint or project wording; license status is checked separately for source, model, database and service reuse.
- **Example:** Prefer SeqKit2 DOI `10.1002/imt2.191` for SeqKit 2.x; cite pyfaidx with `10.7287/peerj.preprints.970v1` as a preprint; cite dnaio with software archive DOI `10.5281/zenodo.10548864`, not as a peer-reviewed article.

See [43) Publication audit for 320 DNA and genomics tools](43_publication_audit.en.md) for the itemized journal-eligibility audit.
