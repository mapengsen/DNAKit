# Variant Calling and Statistics

## 1) Call small variants with FreeBayes

- **Purpose:** Pass caller-prepared BAM/CRAM files and a reference FASTA to FreeBayes and write a VCF.
- **Method:** Run the fixed `freebayes -f reference -p ploidy ...` command with optional minimum alternate-fraction and count filters; FreeBayes source is not copied.
- **API:** `dnakit.variants.call_freebayes(reference, alignments, output_directory, ploidy=2, ...)`
- **Input/output:** Prepared alignments, a reference FASTA and a new output directory; returns a JSON-serializable stage envelope with `variants.vcf`.
- **Citation and license:** Garrison & Marth (2012), [arXiv:1207.3907](https://arxiv.org/abs/1207.3907); the official repository records MIT. The caller installs FreeBayes separately.

## 2) Filter and recode a VCF with VCFtools

- **Purpose:** Apply bounded indel/SNP, missingness and minimum-allele-frequency filters to a local VCF.
- **API:** `dnakit.variants.filter_vcf_vcftools(source, output_directory, remove_indels=False, remove_snps=False, max_missing=None, min_allele_frequency=None, ...)`
- **Output:** VCFtools `filtered.recode.vcf` and diagnostic artifacts; the input is never overwritten.
- **Citation and license:** Danecek et al. (2011), DOI `10.1093/bioinformatics/btr330`; the official repository records LGPL-3.0. The full VCFtools option surface is not re-exported.

## 3) Annotate a VCF with vcfanno

- **Purpose:** Add INFO annotations from caller-provided VCF, BED or indexed annotation files using a vcfanno TOML configuration.
- **API:** `dnakit.variants.annotate_vcf_vcfanno(source, configuration, output_directory, threads=1, ...)`.
- **Output:** `annotated.vcf`, stderr and a structured execution envelope. The caller supplies the configuration and annotation data; DNAKit does not download them.
- **Citation and license:** Pedersen et al. (2016), [10.1186/s13059-016-0973-5](https://doi.org/10.1186/s13059-016-0973-5); the official repository is MIT. The complete TOML/Lua option surface is not re-exported.

## 4) Merge structural-variant VCFs with SURVIVOR

- **Purpose:** Run SURVIVOR's `merge` mode on at least two SV VCFs with bounded breakpoint distance, support, type, strand and minimum-SV-size settings.
- **API:** `dnakit.variants.merge_structural_variants_survivor(sources, output_directory, max_distance=100, min_support=1, min_sv_size=30, ...)`.
- **Output:** `merged.vcf`, an input manifest and diagnostic log in a fresh output directory.
- **Citation and license:** Jeffares et al. (2017), [10.1038/NCOMMS14061](https://doi.org/10.1038/NCOMMS14061); the official repository is MIT. SURVIVOR simulation, conversion and other filtering commands remain outside this facade.

## 5) Build and query a GEMINI database

- **Purpose:** Load a prepared VCF into a local GEMINI variant database and run a bounded SQL query against it.
- **Method:** Call GEMINI's documented `load` and `query` subcommands; the caller prepares annotation databases, annotation resources and downloads, and DNAKit does not download or fabricate annotation results implicitly.
- **API:** `dnakit.variants.load_gemini_database(source_vcf, output_directory, annotation_type=None, cores=1, ...)`; `dnakit.variants.query_gemini_database(database, query, output_directory, header=True, result_format=None, ...)`.
- **Input/output:** `load` requires a local VCF and a new output directory, returning a JSON-serializable envelope containing `variants.db`; `query` requires a local GEMINI database and SQL text, returning `query.tsv`, query parameters and diagnostics.
- **Example code:**

```python
from dnakit.variants import load_gemini_database, query_gemini_database

loaded = load_gemini_database("prepared.vcf", "gemini_database")
queried = query_gemini_database(
    loaded["result"]["database"],
    "SELECT chrom, start, ref, alt FROM variants LIMIT 10",
    "gemini_query",
)
print(queried["result"]["output"])
```

- **Citation and license:** The official GEMINI repository records MIT; cite Paila et al. (2013), [10.1371/journal.pcbi.1003153](https://doi.org/10.1371/journal.pcbi.1003153). GEMINI, its annotation resources and database-building environment remain caller-installed and caller-prepared.

## 6) Annotate a VCF with Ensembl VEP

- **Purpose:** Run Ensembl Variant Effect Predictor on a caller-prepared VCF and produce functional consequence annotations as VCF, tabular or JSON output.
- **Method:** Use VEP's `--input_file`/`--output_file` boundary with bounded species, assembly, offline cache, reference FASTA, output-format and fork options; the caller supplies caches, transcript resources, plugins and database access.
- **API:** `dnakit.variants.annotate_vcf_ensembl_vep(source_vcf, output_directory, species="homo_sapiens", assembly=None, cache_directory=None, fasta=None, offline=True, output_format="vcf", fork=1, ...)`.
- **Input/output:** Requires a local VCF and a new output directory; by default returns a JSON-serializable envelope containing `annotated.vcf`, with `annotated.tsv` or `annotated.json` available by selecting the output format.
- **Example code:**

```python
from dnakit.variants import annotate_vcf_ensembl_vep

result = annotate_vcf_ensembl_vep(
    "variants.vcf",
    "vep_output",
    assembly="GRCh38",
    cache_directory="/data/vep_cache",
    offline=True,
)
print(result["result"]["output"])
```

- **Citation and license:** The official Ensembl VEP repository records Apache-2.0; cite McLaren et al. (2016), [10.1186/s13059-016-0974-4](https://doi.org/10.1186/s13059-016-0974-4). DNAKit does not download caches, reference sequences or plugins.

## 1) Call variants from BAM {#bio-17}

- **Purpose:** Produce normalized, indexed VCF variants.
- **Method:** Run bcftools mpileup, call and norm.
- **API:** `dnakit.bioinformatics.ngs.call_variants(bam, reference, output_directory, ploidy=2, config=None)`
- **Input:** Sorted, indexed BAM, reference FASTA with `.fai`, explicit ploidy, bcftools on PATH, and a new output directory.
- **Output:** Returns `ToolResult` with `backend`, artifact paths in `artifacts`, and execution `log`.
- **Agent/MCP API:** `dnakit.ngs.call_variants(bam, reference, output_directory, ploidy=2, ...)`; returns JSON-serializable VCF, index, and diagnostic metadata.
- **Citation and license:** cite BCFtools as Danecek et al. (2021), DOI `10.1093/gigascience/giab008` (MIT). The caller installs the backend.
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
