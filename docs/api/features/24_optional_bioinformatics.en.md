# Optional bioinformatics functions

This page documents capabilities that were not previously native to DNAKit and are now exposed
under capability-oriented namespaces: `annotation`, `comparative`, and `molbio`. Users do not
call an `integrations.tooluniverse` namespace.

Install the optional backend with:

```bash
python -m pip install "dnakit[external-tools]"
```

No function trains or fine-tunes a model. Ensembl, ClinVar, dbSNP, and gnomAD calls require network
access; dN/dS and Golden Gate run locally. One backend tool
is loaded per call. Every function returns `ProviderResult`; use `result.data` for the provider
payload and retain `result.parameters`, `result.metadata`, and `result.provenance` for auditing.

## 1) VEP annotation from HGVS

- **Purpose:** Obtain transcript consequences for one HGVS expression from Ensembl VEP.
- **API:** `dnakit.annotation.annotate_variant_vep(hgvs_notation[required], species[optional])`
- **Input:** One HGVS expression; species defaults to `human`.
- **Example:**

```python
from dnakit.annotation import annotate_variant_vep

result = annotate_variant_vep("NM_000546.6:c.215C>G")
print(result.data)
```

- **Result:** Structured Ensembl VEP consequence records.

## 2) VEP annotation from rsID

- **Purpose:** Annotate one dbSNP rsID with Ensembl VEP.
- **API:** `dnakit.annotation.annotate_rsid_vep(rsid[required], species[optional])`
- **Input:** An rsID such as `rs28934578`.
- **Example:**

```python
from dnakit.annotation import annotate_rsid_vep

result = annotate_rsid_vep("rs28934578")
print(result.data)
```

- **Result:** VEP consequences associated with the rsID.

## 3) Variant identifier recoding

- **Purpose:** Resolve an rsID or HGVS-like identifier to equivalent Ensembl representations.
- **API:** `dnakit.annotation.recode_variant(variant_id[required], species[optional])`
- **Input:** One variant identifier and an optional species.
- **Example:**

```python
from dnakit.annotation import recode_variant

result = recode_variant("rs28934578")
print(result.data)
```

- **Result:** Equivalent HGVS, genomic-location, and database identifiers.

## 4) ClinVar variant search

- **Purpose:** Search ClinVar by gene, condition, ClinVar ID, HGVS, or protein change.
- **API:** `dnakit.annotation.search_clinvar_variants(gene[optional], condition[optional], variant_id[optional], variant_name[optional], clinical_significance[optional], limit[optional])`
- **Input:** At least one search field; `limit` must be 1–100.
- **Example:**

```python
from dnakit.annotation import search_clinvar_variants

result = search_clinvar_variants(gene="TP53", variant_name="c.215C>G", limit=10)
print(result.data)
```

- **Result:** ClinVar identifiers and summary records.

## 5) ClinVar variant details

- **Purpose:** Retrieve accession, gene, coordinates, and review status for one ClinVar variation.
- **API:** `dnakit.annotation.get_clinvar_variant(variant_id[required])`
- **Input:** A ClinVar variation ID such as `12345`.
- **Example:**

```python
from dnakit.annotation import get_clinvar_variant

result = get_clinvar_variant("12345")
print(result.data)
```

- **Result:** ClinVar variant details.

## 6) ClinVar clinical significance

- **Purpose:** Retrieve clinical-significance assertions for one ClinVar variation.
- **API:** `dnakit.annotation.get_clinvar_significance(variant_id[required])`
- **Input:** A ClinVar variation ID.
- **Example:**

```python
from dnakit.annotation import get_clinvar_significance

result = get_clinvar_significance("12345")
print(result.data)
```

- **Result:** Pathogenicity classifications and clinical interpretations.

## 7) dbSNP variant information

- **Purpose:** Retrieve alleles and GRCh37/GRCh38 coordinates for one rsID.
- **API:** `dnakit.annotation.get_dbsnp_variant(rsid[required])`
- **Input:** An rsID with or without the `rs` prefix.
- **Example:**

```python
from dnakit.annotation import get_dbsnp_variant

result = get_dbsnp_variant("rs28934578")
print(result.data)
```

- **Result:** dbSNP coordinates, assembly, and allele records.

## 8) dbSNP allele frequencies

- **Purpose:** Retrieve population allele-frequency records for one rsID.
- **API:** `dnakit.annotation.get_dbsnp_frequencies(rsid[required])`
- **Input:** One rsID.
- **Example:**

```python
from dnakit.annotation import get_dbsnp_frequencies

result = get_dbsnp_frequencies("rs28934578")
print(result.data)
```

- **Result:** Study/population allele-count and frequency records.

## 9) gnomAD variant search

- **Purpose:** Resolve an rsID or variant query to canonical `chrom-pos-ref-alt` identifiers.
- **API:** `dnakit.annotation.search_gnomad_variants(query[required], dataset[optional])`
- **Input:** Query text; DNAKit consistently defaults to `gnomad_r4`.
- **Example:**

```python
from dnakit.annotation import search_gnomad_variants

result = search_gnomad_variants("rs28934578")
print(result.data)
```

- **Result:** Matching gnomAD variant IDs and basic metadata.

## 10) gnomAD aggregate variant data

- **Purpose:** Retrieve aggregate frequency data and metadata for one gnomAD variant.
- **API:** `dnakit.annotation.get_gnomad_variant(variant_id[required], dataset[optional])`
- **Input:** A `chrom-pos-ref-alt` identifier; dataset defaults to `gnomad_r4`.
- **Example:**

```python
from dnakit.annotation import get_gnomad_variant

result = get_gnomad_variant("17-7675088-C-G")
print(result.data)
```

- **Result:** Aggregate genome/exome allele counts, allele numbers, and frequencies.

## 11) gnomAD population frequencies

- **Purpose:** Retrieve ancestry- and sex-stratified allele frequencies.
- **API:** `dnakit.annotation.get_gnomad_population_frequencies(variant_id[required], dataset[optional])`
- **Input:** A gnomAD variant ID; dataset defaults to `gnomad_r4`.
- **Example:**

```python
from dnakit.annotation import get_gnomad_population_frequencies

result = get_gnomad_population_frequencies("17-7675088-C-G")
print(result.data)
```

- **Result:** AC, AN, and AF values stratified by ancestry and sex.

## 12) dN/dS selection analysis

- **Purpose:** Estimate Nei–Gojobori dN, dS, and dN/dS with Jukes–Cantor correction.
- **API:** `dnakit.comparative.calculate_dn_ds(sequence_a[required], sequence_b[required])`
- **Input:** Equal-length, codon-aligned, linear canonical DNA values.
- **Example:**

```python
from dnakit import DNASequence
from dnakit.comparative import calculate_dn_ds

result = calculate_dn_ds(DNASequence("ATGGCTGAA"), DNASequence("ATGGCCGAG"))
print(result.data["dN"], result.data["dS"], result.data["dN_dS"])
```

- **Result:** dN, dS, ratio, site/difference counts, and interpretation; the ratio is `None` when dS is zero.

## 13) Golden Gate part design

- **Purpose:** Assign unique non-palindromic junction overhangs and add BsaI/BbsI sites.
- **API:** `dnakit.molbio.design_golden_gate(parts[required], enzyme[optional])`
- **Input:** At least two canonical DNA values; design supports `BsaI` and `BbsI`.
- **Example:**

```python
from dnakit import DNASequence
from dnakit.molbio import design_golden_gate

result = design_golden_gate(
    (DNASequence("ATGGCTGAA"), DNASequence("GCCAAATAA")),
    enzyme="BsaI",
)
print(result.data["parts_with_overhangs"])
```

- **Result:** Left/right overhangs, full synthetic part sequences, and a protocol note.

## 14) Golden Gate reaction assembly

- **Purpose:** Simulate Type IIS digestion and ligation by matching overhangs.
- **API:** `dnakit.molbio.assemble_golden_gate(fragments[required], enzyme[optional], circular[optional], labels[optional])`
- **Input:** At least two canonical DNA values carrying the selected Type IIS sites.
- **Example:**

```python
from dnakit import DNASequence
from dnakit.molbio import assemble_golden_gate

result = assemble_golden_gate(
    (
        DNASequence("GGTCTCAAAACATGGAGGAGCCGCAGTCAGATCCTAGCGTTGAATTCGGATCCCTTTTGAGACC"),
        DNASequence("GGTCTCAAAAGGGATCCAAGCTTACGTACGTATGGAGGAGCCGCAGTCAGATATTTTGAGACC"),
        DNASequence("GGTCTCAAAATCCTAGCGTTGAATTCGGATCCAAGCTTACGTACGTATGGAGCGTTTGAGACC"),
        DNASequence("GGTCTCAAACGATCGATCGATCGATCGATCGGTTTTGAGACC"),
    ),
    labels=("vector", "insert-1", "insert-2", "insert-3"),
)
print(result.data["product_sequence"])
```

- **Result:** Product sequence, length, assembly order, junction overhangs, and unused fragments.

## Boundaries

- Missing optional dependencies raise `BackendUnavailableError` with an installation hint.
- Remote providers may rate-limit requests, enter maintenance, or update database versions; retain provenance.
- ClinVar, VEP, dbSNP, and gnomAD results are research annotations, not clinical diagnoses.
- DNAKit exposes only the allowlisted functions above, not arbitrary ToolUniverse execution.
