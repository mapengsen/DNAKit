# 序列、组装与测序数据查询：参数与补充说明

用法和示例见[功能教程](../../api/features/22_database_query_sequences.md)。

各功能按作用、API、输入、输出、示例代码和结果说明介绍。运行前按说明准备凭据、本地文件与可选依赖；远程结果随数据库版本与查询时间变化。

## 按分类与关键词检索 ENA 生物样本 {#dnakit.databases.enaportal_search_samples-parameters}

`dnakit.databases.enaportal_search_samples(*, query, limit=None, fields=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `query` | string | 必填 | ENA search query. Examples: &#x27;tax_tree(562)&#x27; (E. coli), &#x27;description=&quot;liver&quot;&#x27;, &#x27;tax_tree(9606) AND description=&quot;tumor&quot;&#x27;. |
| `limit` | integer &#124; null | None | Maximum results to return (1-100, default 10). |
| `fields` | string &#124; null | None | Comma-separated fields. Default: &#x27;sample_accession,sample_alias,description,tax_id,scientific_name,first_public&#x27;. |
| `verbose` | bool | True | 运行提示开关。 |


## 检索 ENA 测序研究项目 {#dnakit.databases.enaportal_search_studies-parameters}

`dnakit.databases.enaportal_search_studies(*, query, limit=None, fields=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `query` | string | 必填 | ENA search query. IMPORTANT: Plain text keywords (e.g., &#x27;cancer&#x27;) will cause HTTP 400. Always wrap text keywords with description=&quot;...&quot;. Examples: &#x27;description=&quot;cancer&quot;&#x27;, &#x27;tax_tree(9606)&#x27; (human), &#x27;description=&quot;transcriptomics&quot; AND tax_tree(10090)&#x27; (mouse transcriptomics). |
| `limit` | integer &#124; null | None | Maximum results to return (1-100, default 10). |
| `fields` | string &#124; null | None | Comma-separated fields to return. Default: &#x27;study_accession,study_title,center_name,first_public,description&#x27;. Other fields: &#x27;last_updated,tax_id,scientific_name&#x27;. |
| `verbose` | bool | True | 运行提示开关。 |


## 获取基因、转录本或区间序列 {#dnakit.databases.ensembl_get_sequence-parameters}

`dnakit.databases.ensembl_get_sequence(*, id, type=None, species=None, multiple_sequences=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `id` | string | 必填 | Ensembl gene/transcript ID (e.g., &#x27;ENSG00000139618&#x27; or &#x27;ENST00000380152&#x27;) |
| `type` | string | None | Sequence type: &#x27;genomic&#x27; (genomic DNA), &#x27;cds&#x27; (coding sequence), &#x27;cdna&#x27; (cDNA), &#x27;protein&#x27; (protein/peptide sequence -- &#x27;peptide&#x27; is also accepted and mapped to &#x27;protein&#x27;, the name Ensembl&#x27;s own REST API actually expects) {&quot;enum&quot;: &#91;&quot;genomic&quot;, &quot;cds&quot;, &quot;cdna&quot;, &quot;protein&quot;, &quot;peptide&quot;&#93;, &quot;default&quot;: &quot;genomic&quot;} |
| `species` | string | None | Species name (optional, defaults to human) |
| `multiple_sequences` | boolean | None | Return multiple sequences if available (e.g., all transcripts) {&quot;default&quot;: true} |
| `verbose` | bool | True | 运行提示开关。 |


## 按 Ensembl 编号获取核酸或蛋白质序列 {#dnakit.databases.ensembl_seq_get_id_sequence-parameters}

`dnakit.databases.ensembl_seq_get_id_sequence(*, ensembl_id, type=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `ensembl_id` | string | 必填 | Ensembl stable ID. Examples: &#x27;ENSP00000269305&#x27; (TP53 protein, 393 aa), &#x27;ENST00000269305&#x27; (TP53-201 transcript), &#x27;ENSG00000141510&#x27; (TP53 gene). Protein IDs start with ENSP, transcript with ENST, gene with ENSG. |
| `type` | string | None | Sequence type to return. Options: &#x27;protein&#x27; (amino acid), &#x27;cdna&#x27; (spliced cDNA), &#x27;cds&#x27; (coding sequence), &#x27;genomic&#x27; (unspliced DNA). Default: &#x27;protein&#x27; for ENSP IDs, auto-detected for others. {&quot;default&quot;: &quot;protein&quot;} |
| `verbose` | bool | True | 运行提示开关。 |


## 按基因组坐标获取 Ensembl 参考序列 {#dnakit.databases.ensembl_seq_get_region_sequence-parameters}

`dnakit.databases.ensembl_seq_get_region_sequence(*, region, species=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `region` | string | 必填 | Genomic region in format &#x27;chr:start-end&#x27; or &#x27;chr:start..end:strand&#x27;. Examples: &#x27;17:7668421-7668520&#x27; (TP53 locus, 100bp), &#x27;7:140424943-140425043&#x27; (BRAF locus, 100bp), &#x27;X:15560138-15602945&#x27; (MECP2 gene). Strand: 1 for forward, -1 for reverse. |
| `species` | string | None | Ensembl species name. Default: &#x27;homo_sapiens&#x27;. Examples: &#x27;mus_musculus&#x27;, &#x27;danio_rerio&#x27;, &#x27;drosophila_melanogaster&#x27;. {&quot;default&quot;: &quot;homo_sapiens&quot;} |
| `verbose` | bool | True | 运行提示开关。 |


## 将 GenBank UID 转换为序列编号 {#dnakit.databases.ncbi_fetch_accessions-parameters}

`dnakit.databases.ncbi_fetch_accessions(*, operation=None, uids, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | object &#124; None | None | Operation type (fixed: fetch_accession) {&quot;const&quot;: &quot;fetch_accession&quot;, &quot;default&quot;: &quot;fetch_accession&quot;} |
| `uids` | object | 必填 | GenBank UID(s) from NCBI_search_nucleotide. Can be single UID or array of UIDs. {&quot;oneOf&quot;: &#91;{&quot;type&quot;: &quot;array&quot;, &quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}}, {&quot;type&quot;: &quot;string&quot;}&#93;} |
| `verbose` | bool | True | 运行提示开关。 |


## 按 NCBI 序列编号获取核酸记录 {#dnakit.databases.ncbi_get_sequence-parameters}

`dnakit.databases.ncbi_get_sequence(*, operation=None, accession, format=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | object &#124; None | None | Operation type (fixed: fetch_sequence) {&quot;const&quot;: &quot;fetch_sequence&quot;, &quot;default&quot;: &quot;fetch_sequence&quot;} |
| `accession` | string | 必填 | NCBI accession number (e.g., &#x27;U00096&#x27;, &#x27;NC_045512&#x27;, &#x27;NM_000546&#x27;). Works with GenBank, EMBL, and RefSeq accessions. |
| `format` | string | None | Sequence format: &#x27;fasta&#x27; for FASTA, &#x27;gb&#x27; for GenBank, &#x27;embl&#x27; for EMBL, &#x27;gp&#x27; for GenPept protein. Default: fasta {&quot;enum&quot;: &#91;&quot;fasta&quot;, &quot;gb&quot;, &quot;embl&quot;, &quot;gp&quot;&#93;, &quot;default&quot;: &quot;fasta&quot;} |
| `verbose` | bool | True | 运行提示开关。 |


## 按物种、基因或菌株检索核酸记录 {#dnakit.databases.ncbi_search_nucleotide-parameters}

`dnakit.databases.ncbi_search_nucleotide(*, operation=None, organism=None, gene=None, strain=None, keywords=None, seq_type=None, query=None, limit=None, sort=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | object &#124; None | None | Operation type (fixed: search) {&quot;const&quot;: &quot;search&quot;, &quot;default&quot;: &quot;search&quot;} |
| `organism` | string | None | Organism name to search (e.g., &#x27;Escherichia coli&#x27;, &#x27;Homo sapiens&#x27;). Searches in &#91;Organism&#93; field. |
| `gene` | string | None | Gene name to search (e.g., &#x27;BRCA1&#x27;, &#x27;TP53&#x27;). Searches in &#91;Gene&#93; field. |
| `strain` | string | None | Strain name to filter (e.g., &#x27;K-12&#x27;, &#x27;MG1655&#x27;). Searches in &#91;Strain&#93; field. |
| `keywords` | string | None | Keywords to search in titles (e.g., &#x27;complete genome&#x27;, &#x27;mitochondrial&#x27;). Searches in &#91;Title&#93; field. |
| `seq_type` | string | None | Sequence type filter: &#x27;complete_genome&#x27; for complete genomes, &#x27;mrna&#x27; for mRNA sequences, &#x27;refseq&#x27; for RefSeq curated sequences {&quot;enum&quot;: &#91;&quot;complete_genome&quot;, &quot;mrna&quot;, &quot;refseq&quot;&#93;} |
| `query` | string | None | Free-form search query if not using structured filters. Use NCBI query syntax (e.g., &#x27;(Escherichia coli&#91;Organism&#93;) AND (K-12&#91;Strain&#93;)&#x27;) |
| `limit` | integer | None | Maximum number of results to return (default: 20, max: 100) {&quot;default&quot;: 20, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `sort` | string | None | Sort order for results (default: relevance) {&quot;enum&quot;: &#91;&quot;relevance&quot;, &quot;pub_date&quot;, &quot;title&quot;&#93;, &quot;default&quot;: &quot;relevance&quot;} |
| `verbose` | bool | True | 运行提示开关。 |


## 按研究与实验条件检索 SRA 测序运行 {#dnakit.databases.ncbi_sra_search_runs-parameters}

`dnakit.databases.ncbi_sra_search_runs(*, operation=None, study=None, organism=None, strategy=None, platform=None, source=None, query=None, limit=None, sort=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | object &#124; None | None | Operation type (fixed: search) {&quot;const&quot;: &quot;search&quot;} |
| `study` | string | None | Study accession to search (e.g., &#x27;SRP000001&#x27;, &#x27;ERP000001&#x27;, &#x27;DRP000001&#x27;). Searches in &#91;Study&#93; field. |
| `organism` | string | None | Organism name to filter (e.g., &#x27;Homo sapiens&#x27;, &#x27;Mus musculus&#x27;, &#x27;SARS-CoV-2&#x27;). Searches in &#91;Organism&#93; field. |
| `strategy` | string | None | Sequencing strategy to filter. Common values: &#x27;RNA-Seq&#x27; (transcriptomics), &#x27;WGS&#x27; (whole genome), &#x27;WXS&#x27; (exome), &#x27;ChIP-Seq&#x27; (chromatin), &#x27;ATAC-Seq&#x27; (accessibility), &#x27;Bisulfite-Seq&#x27; (methylation), &#x27;AMPLICON&#x27; (targeted), &#x27;Hi-C&#x27; (3D genome) |
| `platform` | string | None | Sequencing platform to filter. Common values: &#x27;ILLUMINA&#x27;, &#x27;OXFORD_NANOPORE&#x27;, &#x27;PACBIO_SMRT&#x27;, &#x27;ION_TORRENT&#x27;, &#x27;BGISEQ&#x27; |
| `source` | string | None | Library source to filter. Common values: &#x27;GENOMIC&#x27; (DNA), &#x27;TRANSCRIPTOMIC&#x27; (RNA), &#x27;METAGENOMIC&#x27;, &#x27;METATRANSCRIPTOMIC&#x27;, &#x27;SYNTHETIC&#x27;, &#x27;VIRAL_RNA&#x27; |
| `query` | string | None | Free-form search query if not using structured filters. Use NCBI query syntax (e.g., &#x27;(Homo sapiens&#91;Organism&#93;) AND (RNA-Seq&#91;Strategy&#93;)&#x27;) |
| `limit` | integer | None | Maximum number of results to return (default: 20, max: 100) {&quot;default&quot;: 20, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `sort` | string | None | Sort order for results (default: relevance) {&quot;enum&quot;: &#91;&quot;relevance&quot;, &quot;pub_date&quot;&#93;, &quot;default&quot;: &quot;relevance&quot;} |
| `verbose` | bool | True | 运行提示开关。 |


## 查询基因组组装的版本与元数据 {#dnakit.databases.ncbidatasets_get_genome_assembly-parameters}

`dnakit.databases.ncbidatasets_get_genome_assembly(*, accession, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `accession` | string | 必填 | Assembly accession, e.g. &#x27;GCF_000005845.2&#x27; (E. coli K-12) or &#x27;GCF_000001405.40&#x27; (human GRCh38.p14). {&quot;pattern&quot;: &quot;^(GCF&#124;GCA)_&#91;0-9&#93;+(?:\\.&#91;0-9&#93;+)?$&quot;} |
| `verbose` | bool | True | 运行提示开关。 |


## 获取组装中染色体与支架的序列报告 {#dnakit.databases.ncbidatasets_get_sequence_reports-parameters}

`dnakit.databases.ncbidatasets_get_sequence_reports(*, accession, page_size=None, page_token=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `accession` | string | 必填 | Assembly accession, e.g. &#x27;GCF_000005845.2&#x27;. {&quot;pattern&quot;: &quot;^(GCF&#124;GCA)_&#91;0-9&#93;+(?:\\.&#91;0-9&#93;+)?$&quot;} |
| `page_size` | integer &#124; null | None | Maximum sequence records in this page (default 100, maximum 1000). {&quot;minimum&quot;: 1, &quot;maximum&quot;: 1000} |
| `page_token` | string &#124; null | None | Opaque next_page_token returned by a previous call. Omit for the first page. {&quot;minLength&quot;: 1, &quot;pattern&quot;: &quot;.*\\S.*&quot;} |
| `verbose` | bool | True | 运行提示开关。 |


## 按分类编号查询物种谱系与统计 {#dnakit.databases.ncbidatasets_get_taxonomy-parameters}

`dnakit.databases.ncbidatasets_get_taxonomy(*, tax_id, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `tax_id` | string | 必填 | NCBI Taxonomy ID (numeric). Examples: &#x27;9606&#x27; (Homo sapiens), &#x27;10090&#x27; (Mus musculus), &#x27;7227&#x27; (Drosophila melanogaster), &#x27;562&#x27; (E. coli). {&quot;pattern&quot;: &quot;^&#91;1-9&#93;&#91;0-9&#93;*$&quot;} |
| `verbose` | bool | True | 运行提示开关。 |


## 列出指定分类单元的基因组组装 {#dnakit.databases.ncbidatasets_list_genomes_by_taxon-parameters}

`dnakit.databases.ncbidatasets_list_genomes_by_taxon(*, taxon, limit=None, reference_only=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `taxon` | string | 必填 | NCBI tax id or scientific name, e.g. &#x27;562&#x27;, &#x27;Escherichia coli&#x27;, &#x27;Saccharomyces cerevisiae&#x27;. {&quot;minLength&quot;: 1, &quot;pattern&quot;: &quot;.*\\S.*&quot;} |
| `limit` | integer &#124; null | None | Max assemblies to return (default 20, max 100). {&quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `reference_only` | boolean &#124; null | None | If true, return only reference/representative genomes. |
| `verbose` | bool | True | 运行提示开关。 |


## 按部分名称查找物种分类候选 {#dnakit.databases.ncbidatasets_suggest_taxonomy-parameters}

`dnakit.databases.ncbidatasets_suggest_taxonomy(*, query, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `query` | string | 必填 | Partial organism name to search. Examples: &#x27;drosophila&#x27;, &#x27;escherichia&#x27;, &#x27;arabidopsis&#x27;, &#x27;saccharomyces&#x27;. {&quot;minLength&quot;: 1, &quot;pattern&quot;: &quot;.*\\S.*&quot;} |
| `verbose` | bool | True | 运行提示开关。 |


## 按实验编号查询 SRA 测序元数据 {#dnakit.databases.sra_get_experiment-parameters}

`dnakit.databases.sra_get_experiment(*, accession, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `accession` | string | 必填 | SRA accession number (e.g., &#x27;SRX26927637&#x27;, &#x27;SRP000001&#x27;, &#x27;SRS000001&#x27;) |
| `verbose` | bool | True | 运行提示开关。 |


## 检索 SRA 测序实验 {#dnakit.databases.sra_search_experiments-parameters}

`dnakit.databases.sra_search_experiments(*, query=None, organism=None, library_strategy=None, platform=None, limit=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `query` | string | None | Search keywords (e.g., &#x27;human gut metagenome&#x27;, &#x27;BRCA1 RNA-Seq&#x27;, &#x27;single cell lung cancer&#x27;) |
| `organism` | string | None | Organism name to filter by (e.g., &#x27;Homo sapiens&#x27;, &#x27;Mus musculus&#x27;, &#x27;human gut metagenome&#x27;) |
| `library_strategy` | string | None | Sequencing strategy filter (e.g., &#x27;RNA-Seq&#x27;, &#x27;WGS&#x27;, &#x27;WXS&#x27;, &#x27;ChIP-Seq&#x27;, &#x27;ATAC-seq&#x27;, &#x27;Bisulfite-Seq&#x27;, &#x27;AMPLICON&#x27;) |
| `platform` | string | None | Sequencing platform filter (e.g., &#x27;ILLUMINA&#x27;, &#x27;PACBIO_SMRT&#x27;, &#x27;OXFORD_NANOPORE&#x27;, &#x27;ION_TORRENT&#x27;) |
| `limit` | integer | None | Maximum results to return (1-50, default 10) {&quot;default&quot;: 10} |
| `verbose` | bool | True | 运行提示开关。 |


## 批量将序列 UID 解析为带版本的编号 {#dnakit.search.ncbi.accessions-parameters}

`dnakit.search.ncbi.accessions(uids, *, database='nuccore', config=None, api_base_url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils', verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `uids` | Sequence&#91;str &#124; int&#93; | 必填 | NCBI 序列 UID 列表，按每批 200 个解析。 |
| `database` | Literal&#91;&#x27;nuccore&#x27;, &#x27;protein&#x27;&#93; | &#x27;nuccore&#x27; | NCBI 序列数据库 nuccore 或 protein。 |
| `config` | SearchConfig &#124; None | None | 查询使用 SearchConfig；下载使用 DownloadConfig。 |
| `api_base_url` | str | &#x27;https://eutils.ncbi.nlm.nih.gov/entrez/eutils&#x27; | 服务 API 根地址；默认值见签名。 |
| `verbose` | bool | True | 是否显示 DNAKit 运行提示；不会改变返回值。 |
