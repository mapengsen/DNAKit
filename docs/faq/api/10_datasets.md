# 数据划分：参数与补充说明

用法和示例见[功能教程](../../api/features/10_datasets.md)。

各功能按作用、API、输入、输出、示例代码和结果说明介绍。运行前按说明准备凭据、本地文件与可选依赖；远程结果随数据库版本与查询时间变化。

ENCODE 接口返回文本：多数是 JSON，`get_citations` 可返回 BibTeX/RIS，`export_data` 可返回 CSV/TSV。仅 JSON 文本使用 `json.loads()`。追踪与凭据保存在 `~/.dnakit-encode/`；持久化凭据功能需要 `dnakit[encode-credentials]`。

## 比较 ENCODE 实验的联合分析兼容性 {#dnakit.encode.compare_experiments-parameters}

`dnakit.encode.compare_experiments(accession1, accession2, *, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `accession1` | str | 必填 | First experiment accession (e.g., &quot;ENCSR133RZO&quot;) |
| `accession2` | str | 必填 | Second experiment accession (e.g., &quot;ENCSR000AKS&quot;) |
| `verbose` | bool | True | 是否显示 DNAKit 运行提示；不会改变返回值。 |


## 将追踪的 ENCODE 实验导出为表格 {#dnakit.encode.export_data-parameters}

`dnakit.encode.export_data(format='csv', assay_title=None, organism=None, organ=None, *, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `format` | Literal&#91;&#x27;csv&#x27;, &#x27;tsv&#x27;, &#x27;json&#x27;&#93; | &#x27;csv&#x27; | Output format: - &quot;csv&quot;: Comma-separated values (default, for Excel/spreadsheets) - &quot;tsv&quot;: Tab-separated values (for R, pandas) - &quot;json&quot;: JSON array (for programmatic use) |
| `assay_title` | str &#124; None | None | Filter by assay type (partial match) |
| `organism` | str &#124; None | None | Filter by organism (partial match) |
| `organ` | str &#124; None | None | Filter by organ (partial match) |
| `verbose` | bool | True | 是否显示 DNAKit 运行提示；不会改变返回值。 |


## 获取追踪实验的论文与参考文献 {#dnakit.encode.get_citations-parameters}

`dnakit.encode.get_citations(accession=None, export_format='json', *, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `accession` | str &#124; None | None | Specific experiment accession. If None, returns all publications. |
| `export_format` | Literal&#91;&#x27;json&#x27;, &#x27;bibtex&#x27;, &#x27;ris&#x27;&#93; | &#x27;json&#x27; | Output format: - &quot;json&quot;: Structured data (default) - &quot;bibtex&quot;: BibTeX format for LaTeX - &quot;ris&quot;: RIS format (Endnote, Zotero, Mendeley) |
| `verbose` | bool | True | 是否显示 DNAKit 运行提示；不会改变返回值。 |


## 查询衍生文件的来源记录 {#dnakit.encode.get_provenance-parameters}

`dnakit.encode.get_provenance(file_path=None, source_accession=None, *, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `file_path` | str &#124; None | None | Get provenance for a specific derived file |
| `source_accession` | str &#124; None | None | List all files derived from a specific ENCODE accession |
| `verbose` | bool | True | 是否显示 DNAKit 运行提示；不会改变返回值。 |


## 查询实验关联的外部参考标识 {#dnakit.encode.get_references-parameters}

`dnakit.encode.get_references(experiment_accession=None, reference_type=None, *, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `experiment_accession` | str &#124; None | None | Filter by specific experiment (optional) |
| `reference_type` | Literal&#91;&#x27;pmid&#x27;, &#x27;doi&#x27;, &#x27;nct_id&#x27;, &#x27;preprint_doi&#x27;, &#x27;geo_accession&#x27;, &#x27;other&#x27;&#93; &#124; None | None | Filter by reference type (optional): &quot;pmid&quot;, &quot;doi&quot;, &quot;nct_id&quot;, &quot;preprint_doi&quot;, &quot;geo_accession&quot;, &quot;other&quot; |
| `verbose` | bool | True | 是否显示 DNAKit 运行提示；不会改变返回值。 |


## 为追踪实验关联外部参考标识 {#dnakit.encode.link_reference-parameters}

`dnakit.encode.link_reference(experiment_accession, reference_type, reference_id, description='', *, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `experiment_accession` | str | 必填 | ENCODE experiment accession (e.g., &quot;ENCSR133RZO&quot;) |
| `reference_type` | Literal&#91;&#x27;pmid&#x27;, &#x27;doi&#x27;, &#x27;nct_id&#x27;, &#x27;preprint_doi&#x27;, &#x27;geo_accession&#x27;, &#x27;other&#x27;&#93; | 必填 | Type of external reference: - &quot;pmid&quot;: PubMed ID (e.g., &quot;32728249&quot;) - &quot;doi&quot;: DOI (e.g., &quot;10.1038/s41586-020-2493-4&quot;) - &quot;nct_id&quot;: ClinicalTrials.gov ID (e.g., &quot;NCT04567890&quot;) - &quot;preprint_doi&quot;: bioRxiv/medRxiv DOI - &quot;geo_accession&quot;: GEO accession (e.g., &quot;GSE123456&quot;) - &quot;other&quot;: Any other identifier |
| `reference_id` | str | 必填 | The actual identifier value |
| `description` | str | &#x27;&#x27; | Optional description of why this reference is linked |
| `verbose` | bool | True | 是否显示 DNAKit 运行提示；不会改变返回值。 |


## 列出本地追踪的 ENCODE 实验 {#dnakit.encode.list_tracked-parameters}

`dnakit.encode.list_tracked(assay_title=None, organism=None, organ=None, *, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `assay_title` | str &#124; None | None | Filter by assay type (partial match) |
| `organism` | str &#124; None | None | Filter by organism (partial match) |
| `organ` | str &#124; None | None | Filter by organ (partial match) |
| `verbose` | bool | True | 是否显示 DNAKit 运行提示；不会改变返回值。 |


## 登记衍生文件与数据来源 {#dnakit.encode.log_derived_file-parameters}

`dnakit.encode.log_derived_file(file_path, source_accessions, description='', file_type='', tool_used='', parameters='', *, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `file_path` | str | 必填 | Path to the derived file you created |
| `source_accessions` | list&#91;str&#93; | 必填 | List of ENCODE accessions this file was derived from (experiment or file accessions, e.g., &#91;&quot;ENCSR133RZO&quot;, &quot;ENCFF635JIA&quot;&#93;) |
| `description` | str | &#x27;&#x27; | What this derived file contains |
| `file_type` | str | &#x27;&#x27; | Type of file (e.g., &quot;filtered_peaks&quot;, &quot;merged_signal&quot;, &quot;differential&quot;) |
| `tool_used` | str | &#x27;&#x27; | Tool/software used to create it (e.g., &quot;bedtools intersect&quot;, &quot;DESeq2&quot;) |
| `parameters` | str | &#x27;&#x27; | Parameters or command used |
| `verbose` | bool | True | 是否显示 DNAKit 运行提示；不会改变返回值。 |


## 按字段汇总追踪的实验集合 {#dnakit.encode.summarize_collection-parameters}

`dnakit.encode.summarize_collection(assay_title=None, organism=None, organ=None, *, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `assay_title` | str &#124; None | None | Filter by assay type (partial match) |
| `organism` | str &#124; None | None | Filter by organism (partial match) |
| `organ` | str &#124; None | None | Filter by organ (partial match) |
| `verbose` | bool | True | 是否显示 DNAKit 运行提示；不会改变返回值。 |


## 在本地追踪 ENCODE 实验及其来源信息 {#dnakit.encode.track_experiment-parameters}

`dnakit.encode.track_experiment(accession, fetch_publications=True, fetch_pipelines=True, notes='', *, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `accession` | str | 必填 | ENCODE experiment accession (e.g., &quot;ENCSR133RZO&quot;) |
| `fetch_publications` | bool | True | Also fetch and store publications/citations (default True) |
| `fetch_pipelines` | bool | True | Also fetch and store pipeline/analysis info (default True) |
| `notes` | str | &#x27;&#x27; | Optional notes to attach to this experiment |
| `verbose` | bool | True | 是否显示 DNAKit 运行提示；不会改变返回值。 |
