# Variant Calling and Statistics: parameters and notes

See the [tutorial](../../api/features/33_variant_calling.md) for usage and examples.

Each function includes its purpose, API, inputs, outputs, example code and result description. Prepare any required credentials, local files and optional dependencies; remote results vary with database releases and query time.

## Calculate the fraction of coding variants matching filters {#dnakit.databases.coding_variant_fraction-parameters}

`dnakit.databases.coding_variant_fraction(*, file, vaf_threshold=None, annotation=None, header_rows=None, vaf_col=None, so_col=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `file` | string | Required | Path to variant Excel (.xlsx) or CSV file. |
| `vaf_threshold` | number | None | Keep variants with VAF strictly less than this value. Default 0.3. |
| `annotation` | string | None | Sequence Ontology term to count in numerator. Default &#x27;synonymous_variant&#x27;. |
| `header_rows` | integer | None | 1 for flat headers, 2 for MultiIndex (common in Excel variant reports). Default 1. {&quot;enum&quot;: &#91;1, 2&#93;} |
| `vaf_col` | string &#124; null | None | Explicit VAF column name. If empty, auto-detect from &#x27;variant allele freq&#x27; or &#x27;vaf&#x27;. |
| `so_col` | string &#124; null | None | Explicit Sequence Ontology column name. If empty, auto-detect from &#x27;sequence ontology&#x27;. |
| `verbose` | bool | True | Progress flag. |


## Count filtered variants in a VCF or BCF file {#dnakit.databases.vcf_count_variants-parameters}

`dnakit.databases.vcf_count_variants(*, operation=None, vcf_path, pass_only=None, regions=None, min_qual=None, include_expr=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | object &#124; None | None | Operation (fixed) {&quot;const&quot;: &quot;count_variants&quot;} |
| `vcf_path` | string | Required | Path to a local .vcf, .vcf.gz, or .bcf file |
| `pass_only` | boolean &#124; null | None | If true, count only FILTER=PASS or &#x27;.&#x27; records |
| `regions` | string &#124; null | None | Restrict to a region, e.g. &#x27;chr1&#x27; or &#x27;chr1:1-1000000&#x27; |
| `min_qual` | number &#124; null | None | Keep only records with QUAL &gt;= this value |
| `include_expr` | string &#124; null | None | Extra bcftools -i include expression, e.g. &#x27;INFO/AF&gt;0.01&#x27; |
| `verbose` | bool | True | Progress flag. |


## Normalize VCF variants and multiallelic sites {#dnakit.databases.vcf_normalize-parameters}

`dnakit.databases.vcf_normalize(*, operation=None, vcf_path, multiallelics=None, reference_fasta=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | object &#124; None | None | Operation (fixed) {&quot;const&quot;: &quot;normalize&quot;} |
| `vcf_path` | string | Required | Path to a local .vcf, .vcf.gz, or .bcf file |
| `multiallelics` | string &#124; null | None | &#x27;split&#x27; (default) breaks multiallelic records into biallelic; &#x27;join&#x27; merges; &#x27;none&#x27; leaves as-is {&quot;enum&quot;: &#91;&quot;split&quot;, &quot;join&quot;, &quot;none&quot;, null&#93;} |
| `reference_fasta` | string &#124; null | None | Optional reference FASTA path; enables indel left-alignment (-f) |
| `verbose` | bool | True | Progress flag. |


## Summarize VCF variant types and sample genotype statistics {#dnakit.databases.vcf_summary_stats-parameters}

`dnakit.databases.vcf_summary_stats(*, operation=None, vcf_path, pass_only=None, regions=None, min_qual=None, include_expr=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | object &#124; None | None | Operation (fixed) {&quot;const&quot;: &quot;summary_stats&quot;} |
| `vcf_path` | string | Required | Path to a local .vcf, .vcf.gz, or .bcf file |
| `pass_only` | boolean &#124; null | None | If true, count only FILTER=PASS or &#x27;.&#x27; records |
| `regions` | string &#124; null | None | Restrict to a region, e.g. &#x27;chr1&#x27; or &#x27;chr1:1-1000000&#x27; (streamed, no index needed) |
| `min_qual` | number &#124; null | None | Keep only records with QUAL &gt;= this value |
| `include_expr` | string &#124; null | None | Extra bcftools -i include expression, e.g. &#x27;INFO/DP&gt;10&#x27; |
| `verbose` | bool | True | Progress flag. |
