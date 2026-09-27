# Variant Calling and Statistics

## 1) Call variants from BAM {#bio-17}

- **Purpose:** Produce normalized, indexed VCF variants.
- **Method:** Run bcftools mpileup, call and norm.
- **API:** `dnakit.bioinformatics.ngs.call_variants(bam, reference, output_directory, ploidy=2, config=None)`
- **Input:** Sorted, indexed BAM, reference FASTA with `.fai`, explicit ploidy, bcftools on PATH, and a new output directory.
- **Output:** Returns `ToolResult` with `backend`, artifact paths in `artifacts`, and execution `log`.
- **Example code:**

```python
from dnakit.bioinformatics.ngs import call_variants

result = call_variants("aligned.bam", "reference.fa", "new_variants", ploidy=2)
print(result.artifacts["variants.vcf.gz"])
```

- **Example output:**

Returns `ToolResult`; `artifacts['variants.vcf.gz']` points to the compressed VCF, with an index also produced. This external-program example was not executed.

Dependencies and details: [FAQ](../../faq/bioinformatics.md#reads).

<!-- BEGIN GENERATED API: vcf -->

<span id="vcf"></span>

## 2) Count filtered variants in a VCF or BCF file {#dnakit.databases.vcf_count_variants}

- **Purpose:** Count filtered variants in a VCF or BCF file.
- **Method:** Filter with bcftools view, then count variant records with bcftools stats.
- **API：** `dnakit.databases.vcf_count_variants(*, operation=None, vcf_path, pass_only=None, regions=None, min_qual=None, include_expr=None, verbose=True)`
- **Input:** Provide `vcf_path`. Requires a local VCF/BCF file and `bcftools` on PATH. Prepare the local files and external programs used by the example, and replace paths with your input or output locations.
- **Output:** Returns a mapping with `status`, `data`, `metadata` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.vcf_count_variants(
    operation='count_variants',
    vcf_path='/tmp/sample.vcf.gz',
    min_qual=30,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "status": "success",
  "data": "<variant counts and filters_applied>",
  "metadata": "<dict>"
}
```

[Full parameters and usage notes](../../faq/api/33_variant_calling.md#dnakit.databases.vcf_count_variants-parameters)

## 3) Normalize VCF variants and multiallelic sites {#dnakit.databases.vcf_normalize}

- **Purpose:** Normalize VCF variants and multiallelic sites.
- **Method:** Use bcftools norm to split or join multiallelic sites, left-align when a reference is supplied, and compare counts.
- **API：** `dnakit.databases.vcf_normalize(*, operation=None, vcf_path, multiallelics=None, reference_fasta=None, verbose=True)`
- **Input:** Provide `vcf_path`. Requires a local VCF/BCF file and `bcftools` on PATH. Prepare the local files and external programs used by the example, and replace paths with your input or output locations.
- **Output:** Returns a mapping with `status`, `data`, `metadata` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.vcf_normalize(
    operation='normalize',
    vcf_path='/tmp/sample.vcf.gz',
    multiallelics='split',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "status": "success",
  "data": {
    "before": "<variant counts>",
    "after": "<variant counts>",
    "left_aligned": "<bool>",
    "multiallelics": "<str>"
  },
  "metadata": "<dict>"
}
```

[Full parameters and usage notes](../../faq/api/33_variant_calling.md#dnakit.databases.vcf_normalize-parameters)

## 4) Summarize VCF variant types and sample genotype statistics {#dnakit.databases.vcf_summary_stats}

- **Purpose:** Summarize VCF variant types and sample genotype statistics.
- **Method:** Run bcftools stats on the filtered VCF and parse variant types and sample statistics.
- **API：** `dnakit.databases.vcf_summary_stats(*, operation=None, vcf_path, pass_only=None, regions=None, min_qual=None, include_expr=None, verbose=True)`
- **Input:** Provide `vcf_path`. Requires a local VCF/BCF file and `bcftools` on PATH. Prepare the local files and external programs used by the example, and replace paths with your input or output locations.
- **Output:** Returns a mapping with `status`, `data`, `metadata` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.vcf_summary_stats(
    operation='summary_stats',
    vcf_path='/tmp/sample.vcf.gz',
    pass_only=True,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "status": "success",
  "data": "<variant and sample statistics; filters_applied>",
  "metadata": "<dict>"
}
```

[Full parameters and usage notes](../../faq/api/33_variant_calling.md#dnakit.databases.vcf_summary_stats-parameters)

<!-- END GENERATED API: vcf -->

<!-- BEGIN GENERATED API: variant-statistics -->

<span id="variant-statistics"></span>

## 5) Calculate the fraction of coding variants matching filters {#dnakit.databases.coding_variant_fraction}

- **Purpose:** Calculate the fraction of coding variants matching filters.
- **Method:** Use coding variants as the denominator and count those meeting the annotation and VAF threshold.
- **API：** `dnakit.databases.coding_variant_fraction(*, file, vaf_threshold=None, annotation=None, header_rows=None, vaf_col=None, so_col=None, verbose=True)`
- **Input:** Provide `file`.
- **Output:** Returns a mapping with `total_variants`, `vaf_filtered`, `coding_variants_below_vaf`, `matching_annotation` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.coding_variant_fraction(
    file='/tmp/tu_variants_demo.csv',
    vaf_threshold=0.3,
    annotation='synonymous_variant',
    header_rows=1,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "total_variants": "<integer>",
  "vaf_filtered": "<integer>",
  "coding_variants_below_vaf": "<integer>",
  "matching_annotation": "<integer>",
  "fraction": "<number>",
  "naive_fraction": "<number>",
  "vaf_column": "<string>",
  "so_column": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/33_variant_calling.md#dnakit.databases.coding_variant_fraction-parameters)

<!-- END GENERATED API: variant-statistics -->
