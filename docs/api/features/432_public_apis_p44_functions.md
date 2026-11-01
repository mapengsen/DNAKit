# P44 公共数据库操作

通过固定的官方端点执行真实检索和下载。普通测试与下方回放示例不联网；未声称已经逐条完成线上数据验证。

下载继承字节限制、原子写入、校验和与进度条。返回 JSON 可序列化结果，保留服务响应、查询日期和分页信息；分页需按上游令牌显式继续。旧数据集文件名与版本需先通过目录列表确认。

## `datasets_gene_id_download_summary()`

- **作用:** Gene data package download summary (preview) by GeneID. NCBI Datasets v2; one page,
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `datasets_gene_id_download_summary(gene_ids: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `datasets_gene_id_expression_report()`

- **作用:** Gene expression report by GeneID. NCBI Datasets v2; one page, preserving next_page_token.
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `datasets_gene_id_expression_report(gene_ids: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `datasets_gene_id_links()`

- **作用:** Gene external resource links by GeneID. NCBI Datasets v2; one page, preserving
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `datasets_gene_id_links(gene_ids: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `datasets_gene_product_report()`

- **作用:** Gene product report by GeneID, RefSeq accession, symbol+taxon or locus tag. NCBI
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `datasets_gene_product_report(identifier: 'str', *, by: 'str' = 'id', taxon: 'str | None' = None, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `datasets_gene_taxon_annotation_chromosome_summary()`

- **作用:** Gene counts per chromosome by taxon and annotation name. NCBI Datasets v2; one page,
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `datasets_gene_taxon_annotation_chromosome_summary(taxon: 'str', annotation_name: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `datasets_gene_taxon_counts()`

- **作用:** Gene counts by taxon. NCBI Datasets v2; one page, preserving next_page_token.
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `datasets_gene_taxon_counts(taxon: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `datasets_gene_taxon_dataset_report()`

- **作用:** Gene data report by taxon (search). NCBI Datasets v2; one page, preserving
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `datasets_gene_taxon_dataset_report(taxon: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `datasets_gene_taxon_product_report()`

- **作用:** Gene product report by taxon (search). NCBI Datasets v2; one page, preserving
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `datasets_gene_taxon_product_report(taxon: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `datasets_genome_accession_annotation_report()`

- **作用:** Genome annotation report (gene records) by accession. NCBI Datasets v2; one page,
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `datasets_genome_accession_annotation_report(accession: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `datasets_genome_accession_annotation_report_download_summary()`

- **作用:** Genome annotation data package download summary (preview) by accession. NCBI Datasets
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `datasets_genome_accession_annotation_report_download_summary(accession: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `datasets_genome_accession_annotation_summary()`

- **作用:** Genome annotation report summary (facets) by accession. NCBI Datasets v2; one page,
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `datasets_genome_accession_annotation_summary(accession: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `datasets_genome_accession_check()`

- **作用:** Genome assembly accession validity check. NCBI Datasets v2; one page, preserving
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `datasets_genome_accession_check(accessions: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `datasets_genome_accession_download_summary()`

- **作用:** Genome data package download summary (preview) by accession. NCBI Datasets v2; one page,
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `datasets_genome_accession_download_summary(accessions: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `datasets_genome_accession_links()`

- **作用:** Genome assembly external resource links by accession. NCBI Datasets v2; one page,
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `datasets_genome_accession_links(accessions: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `datasets_genome_accession_revision_history()`

- **作用:** Genome assembly revision history by accession. NCBI Datasets v2; one page, preserving
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `datasets_genome_accession_revision_history(accession: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `datasets_genome_report()`

- **作用:** Genome assembly dataset report by identifier (accession, WGS, assembly name). NCBI
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `datasets_genome_report(identifier: 'str', *, by: 'str' = 'accession', parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `datasets_genome_sequence_accession_sequence_assemblies()`

- **作用:** Genome assembly accession for a nucleotide sequence accession. NCBI Datasets v2; one
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `datasets_genome_sequence_accession_sequence_assemblies(accession: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `datasets_genome_taxon_checkm_histogram()`

- **作用:** CheckM quality histogram of genome assemblies by species taxon. NCBI Datasets v2; one
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `datasets_genome_taxon_checkm_histogram(species_taxon: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `datasets_organelle_accession_download()`

- **作用:** Organelle genome data package download by nucleotide accession. NCBI Datasets v2; one
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `datasets_organelle_accession_download(accessions: 'str', *, output_path: 'str', overwrite: 'bool' = False, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `datasets_organelle_accessions_dataset_report()`

- **作用:** Organelle genome data report by nucleotide accession. NCBI Datasets v2; one page,
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `datasets_organelle_accessions_dataset_report(accessions: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `datasets_organelle_taxon_dataset_report()`

- **作用:** Organelle genome data report by taxon (search). NCBI Datasets v2; one page, preserving
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `datasets_organelle_taxon_dataset_report(taxons: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `datasets_protein_accession_download()`

- **作用:** Prokaryote gene data package download by RefSeq protein accession. NCBI Datasets v2; one
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `datasets_protein_accession_download(accessions: 'str', *, output_path: 'str', overwrite: 'bool' = False, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `datasets_virus_accession_annotation_report()`

- **作用:** Virus annotation report by nucleotide accession. NCBI Datasets v2; one page, preserving
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `datasets_virus_accession_annotation_report(accessions: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `datasets_virus_accession_check()`

- **作用:** Virus genome nucleotide accession validity check. NCBI Datasets v2; one page, preserving
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `datasets_virus_accession_check(accessions: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `datasets_virus_download_summary()`

- **作用:** Virus genome data package download summary (preview) by taxon or SARS-CoV-2 protein.
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `datasets_virus_download_summary(identifier: 'str', *, by: 'str' = 'taxon', parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `datasets_virus_taxon_annotation_report()`

- **作用:** Virus annotation report by taxon (search). NCBI Datasets v2; one page, preserving
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `datasets_virus_taxon_annotation_report(taxon: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `ensembl_ga4gh_beacon_query()`

- **作用:** GA4GH Beacon allele presence query. Official GET /ga4gh/beacon/query; upstream defaults
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `ensembl_ga4gh_beacon_query(*, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `ensembl_ga4gh_callsets_record()`

- **作用:** GA4GH callset by ID. Official GET /ga4gh/callsets/:id; upstream defaults retained.
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `ensembl_ga4gh_callsets_record(id: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `ensembl_ga4gh_callsets_search()`

- **作用:** GA4GH genotype callsets search. Official POST /ga4gh/callsets/search; returns the
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `ensembl_ga4gh_callsets_search(query: 'Mapping[str, Any]', *, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `ensembl_ga4gh_datasets_record()`

- **作用:** GA4GH dataset by ID. Official GET /ga4gh/datasets/:id; upstream defaults retained.
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `ensembl_ga4gh_datasets_record(id: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `ensembl_ga4gh_datasets_search()`

- **作用:** GA4GH datasets search. Official POST /ga4gh/datasets/search; returns the service page
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `ensembl_ga4gh_datasets_search(query: 'Mapping[str, Any]', *, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `ensembl_ga4gh_features_record()`

- **作用:** GA4GH feature record by ID. Official GET /ga4gh/features/:id; upstream defaults retained.
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `ensembl_ga4gh_features_record(id: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `ensembl_ga4gh_features_search()`

- **作用:** GA4GH sequence features in a region. Official POST /ga4gh/features/search; returns the
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `ensembl_ga4gh_features_search(query: 'Mapping[str, Any]', *, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `ensembl_ga4gh_featuresets_record()`

- **作用:** GA4GH feature set by ID. Official GET /ga4gh/featuresets/:id; upstream defaults retained.
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `ensembl_ga4gh_featuresets_record(id: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `ensembl_ga4gh_featuresets_search()`

- **作用:** GA4GH feature sets search. Official POST /ga4gh/featuresets/search; returns the service
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `ensembl_ga4gh_featuresets_search(query: 'Mapping[str, Any]', *, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `ensembl_ga4gh_references_record()`

- **作用:** GA4GH reference sequence by ID. Official GET /ga4gh/references/:id; upstream defaults
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `ensembl_ga4gh_references_record(id: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `ensembl_ga4gh_references_search()`

- **作用:** GA4GH reference sequences search. Official POST /ga4gh/references/search; returns the
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `ensembl_ga4gh_references_search(query: 'Mapping[str, Any]', *, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `ensembl_ga4gh_referencesets_record()`

- **作用:** GA4GH reference set by ID. Official GET /ga4gh/referencesets/:id; upstream defaults
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `ensembl_ga4gh_referencesets_record(id: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `ensembl_ga4gh_referencesets_search()`

- **作用:** GA4GH reference sets search. Official POST /ga4gh/referencesets/search; returns the
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `ensembl_ga4gh_referencesets_search(query: 'Mapping[str, Any]', *, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `ensembl_ga4gh_variantannotations_search()`

- **作用:** GA4GH variant annotations in a region. Official POST /ga4gh/variantannotations/search;
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `ensembl_ga4gh_variantannotations_search(query: 'Mapping[str, Any]', *, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `ensembl_ga4gh_variantannotationsets_record()`

- **作用:** GA4GH variant annotation set by ID. Official GET /ga4gh/variantannotationsets/:id;
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `ensembl_ga4gh_variantannotationsets_record(id: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `ensembl_ga4gh_variantannotationsets_search()`

- **作用:** GA4GH variant annotation sets search. Official POST /ga4gh/variantannotationsets/search;
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `ensembl_ga4gh_variantannotationsets_search(query: 'Mapping[str, Any]', *, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `ensembl_ga4gh_variants_record()`

- **作用:** GA4GH variant record by ID. Official GET /ga4gh/variants/:id; upstream defaults retained.
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `ensembl_ga4gh_variants_record(id: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `ensembl_ga4gh_variants_search()`

- **作用:** GA4GH variant calls in a region. Official POST /ga4gh/variants/search; returns the
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `ensembl_ga4gh_variants_search(query: 'Mapping[str, Any]', *, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `ensembl_ga4gh_variantsets_record()`

- **作用:** GA4GH variant set by ID. Official GET /ga4gh/variantsets/:id; upstream defaults retained.
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `ensembl_ga4gh_variantsets_record(id: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `ensembl_ga4gh_variantsets_search()`

- **作用:** GA4GH variant sets search. Official POST /ga4gh/variantsets/search; returns the service
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `ensembl_ga4gh_variantsets_search(query: 'Mapping[str, Any]', *, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `ensembl_info_analysis()`

- **作用:** List annotation analyses of a species. Official GET /info/analysis/:species; upstream
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `ensembl_info_analysis(species: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `ensembl_info_assembly()`

- **作用:** Return top-level assembly sequence metadata from GET info/assembly/:species.
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `ensembl_info_assembly(species: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `ensembl_info_biotypes()`

- **作用:** List gene and transcript biotypes of a species. Official GET /info/biotypes/:species;
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `ensembl_info_biotypes(species: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `ensembl_info_biotypes_groups()`

- **作用:** List biotype groups. Official GET /info/biotypes/groups/:group/:object_type; upstream
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `ensembl_info_biotypes_groups(group: 'str', object_type: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `ensembl_info_biotypes_name()`

- **作用:** Biotype properties by name. Official GET /info/biotypes/name/:name/:object_type;
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `ensembl_info_biotypes_name(name: 'str', object_type: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `ensembl_info_compara_methods()`

- **作用:** List comparative genomics methods. Official GET /info/compara/methods; upstream defaults
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `ensembl_info_compara_methods(*, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `ensembl_info_compara_species_sets()`

- **作用:** List species sets analysed by a comparative method. Official GET
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `ensembl_info_compara_species_sets(method: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `ensembl_info_external_dbs()`

- **作用:** List external databases linked to a species. Official GET /info/external_dbs/:species;
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `ensembl_info_external_dbs(species: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `ensembl_info_genomes_accession()`

- **作用:** Genome record by INSDC accession. Official GET /info/genomes/accession/:accession;
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `ensembl_info_genomes_accession(accession: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `ensembl_info_genomes_assembly()`

- **作用:** Genome record by assembly ID. Official GET /info/genomes/assembly/:assembly_id; upstream
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `ensembl_info_genomes_assembly(assembly_id: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `ensembl_info_genomes_division()`

- **作用:** List genomes in a division. Official GET /info/genomes/division/:division_name; upstream
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `ensembl_info_genomes_division(division_name: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `ensembl_info_genomes_taxonomy()`

- **作用:** List genomes under a taxonomy node. Official GET /info/genomes/taxonomy/:taxon_name;
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `ensembl_info_genomes_taxonomy(taxon_name: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `ensembl_info_variation()`

- **作用:** List variation sources of a species. Official GET /info/variation/:species; upstream
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `ensembl_info_variation(species: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `ensembl_info_variation_consequence_types()`

- **作用:** List variant consequence types. Official GET /info/variation/consequence_types; upstream
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `ensembl_info_variation_consequence_types(*, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `ensembl_info_variation_populations()`

- **作用:** List individuals of a population. Official GET
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `ensembl_info_variation_populations(species: 'str', population_name: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `ensembl_info_variation_populations_list()`

- **作用:** List populations of a species. Official GET /info/variation/populations/:species;
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `ensembl_info_variation_populations_list(species: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `ensembl_map_cds()`

- **作用:** CDS to genomic coordinate conversion. Official GET /map/cds/:id/:region; upstream
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `ensembl_map_cds(id: 'str', region: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `ensembl_map_translation()`

- **作用:** Protein to genomic coordinate conversion. Official GET /map/translation/:id/:region;
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `ensembl_map_translation(id: 'str', region: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `ensembl_overlap_translation()`

- **作用:** Features of a translation (domains, variants). Official GET /overlap/translation/:id;
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `ensembl_overlap_translation(id: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `ensembl_phenotype_accession()`

- **作用:** Phenotype annotations by ontology accession. Official GET
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `ensembl_phenotype_accession(species: 'str', accession: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `ensembl_phenotype_gene()`

- **作用:** Phenotype annotations of a gene. Official GET /phenotype/gene/:species/:gene; upstream
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `ensembl_phenotype_gene(species: 'str', gene: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `ensembl_phenotype_term()`

- **作用:** Phenotype annotations by ontology term. Official GET /phenotype/term/:species/:term;
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `ensembl_phenotype_term(species: 'str', term: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `ensembl_variant_recoder()`

- **作用:** Variant ID, HGVS and SPDI recoding. Official GET /variant_recoder/:species/:id; upstream
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `ensembl_variant_recoder(species: 'str', id: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `ensembl_variation_publication()`

- **作用:** Variants in a publication by PMID or PMCID. Official GET
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `ensembl_variation_publication(species: 'str', publication: 'str', *, publication_type: 'str' = 'pmcid', parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `ensembl_vep_region()`

- **作用:** Variant consequences by genomic region and allele. Official GET
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `ensembl_vep_region(species: 'str', region: 'str', allele: 'str', *, parameters: 'Mapping[str, Any] | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## 离线输入检查示例

**示例代码：**

```python
from dnakit.databases.public_apis_p44_functions import ensembl_info_analysis
from dnakit.exceptions import ConfigurationError
try:
    ensembl_info_analysis("", verbose=False)
except ConfigurationError:
    print("invalid input rejected before network access")
```

**示例输出：**

```text
invalid input rejected before network access
```

## 上游文档

- [官方文档 1](https://rest.ensembl.org/)
- [官方文档 2](https://github.com/ncbi/datasets/blob/master/datasets.openapi.yaml)
- [官方文档 3](https://www.ebi.ac.uk/ena/browser/api/swagger-ui/index.html)
