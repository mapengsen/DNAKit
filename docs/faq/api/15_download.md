# 下载：参数与补充说明

用法和示例见[功能教程](../../api/features/15_download.md)。

各功能按作用、API、输入、输出、示例代码和结果说明介绍。运行前按说明准备凭据、本地文件与可选依赖；远程结果随数据库版本与查询时间变化。

ENCODE 接口返回文本：多数是 JSON，`get_citations` 可返回 BibTeX/RIS，`export_data` 可返回 CSV/TSV。仅 JSON 文本使用 `json.loads()`。追踪与凭据保存在 `~/.dnakit-encode/`；持久化凭据功能需要 `dnakit[encode-credentials]`。

各功能按作用、API、输入、输出、示例代码和结果说明介绍。运行前按说明准备凭据、本地文件与可选依赖；远程结果随数据库版本与查询时间变化。

## 通过 HTTP 或 HTTPS 下载文件 {#dnakit.databases.download_file-parameters}

`dnakit.databases.download_file(*, url, output_path=None, timeout=None, return_content=None, chunk_size=None, follow_redirects=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `url` | string | 必填 | HTTP or HTTPS URL to download from (e.g., https://example.com/file.txt) {&quot;format&quot;: &quot;uri&quot;} |
| `output_path` | string | None | Optional path to save the file. If not specified, file will be saved to system temp directory. {&quot;example&quot;: &quot;/tmp/downloaded_file.txt or C:\\Users\\Downloads\\file.txt&quot;} |
| `timeout` | integer | None | Request timeout in seconds {&quot;minimum&quot;: 1, &quot;maximum&quot;: 300, &quot;default&quot;: 30} |
| `return_content` | boolean | None | If true, return file content as text instead of saving to disk (default: false) {&quot;default&quot;: false} |
| `chunk_size` | integer | None | Download chunk size in bytes (default: 8192) {&quot;minimum&quot;: 1024, &quot;maximum&quot;: 10485760, &quot;default&quot;: 8192} |
| `follow_redirects` | boolean | None | Follow HTTP redirects (default: true) {&quot;default&quot;: true} |
| `verbose` | bool | True | 运行提示开关。 |


## 下载 ENA 核酸序列并保存校验清单 {#dnakit.download.ena.sequence-parameters}

`dnakit.download.ena.sequence(accession, output_path, *, format='fasta', config=None, search_config=None, api_base_url='https://www.ebi.ac.uk/ena/browser/api', verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `accession` | str | 必填 | 版本化序列 accession；ENA 使用 INSDC 编号，NCBI 也支持 RefSeq。 |
| `output_path` | str &#124; Path | 必填 | 目标序列文件路径。 |
| `format` | Literal&#91;&#x27;fasta&#x27;, &#x27;embl&#x27;, &#x27;genbank&#x27;&#93; | &#x27;fasta&#x27; | fasta、genbank 或 embl；提供方及转换依赖见序列下载教程。 |
| `config` | DownloadConfig &#124; None | None | 查询使用 SearchConfig；下载使用 DownloadConfig。 |
| `search_config` | SearchConfig &#124; None | None | 控制查询超时、返回大小和 NCBI 凭据的 SearchConfig。 |
| `api_base_url` | str | &#x27;https://www.ebi.ac.uk/ena/browser/api&#x27; | 服务 API 根地址；默认值见签名。 |
| `verbose` | bool | True | 是否显示 DNAKit 运行提示；不会改变返回值。 |


## 下载 NCBI 核酸或蛋白质序列 {#dnakit.download.ncbi.sequence-parameters}

`dnakit.download.ncbi.sequence(accession, output_path, *, database='nuccore', format='fasta', seq_start=None, seq_stop=None, strand=1, config=None, search_config=None, api_base_url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils', verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `accession` | str | 必填 | 版本化序列 accession；ENA 使用 INSDC 编号，NCBI 也支持 RefSeq。 |
| `output_path` | str &#124; Path | 必填 | 目标序列文件路径。 |
| `database` | Literal&#91;&#x27;nuccore&#x27;, &#x27;protein&#x27;&#93; | &#x27;nuccore&#x27; | NCBI 序列数据库 nuccore 或 protein。 |
| `format` | Literal&#91;&#x27;fasta&#x27;, &#x27;genbank&#x27;, &#x27;embl&#x27;&#93; | &#x27;fasta&#x27; | fasta、genbank 或 embl；提供方及转换依赖见序列下载教程。 |
| `seq_start` | int &#124; None | None | NCBI 核酸子区间起点，1-based 闭区间；须同时指定 seq_stop。 |
| `seq_stop` | int &#124; None | None | NCBI 核酸子区间终点，1-based 闭区间；须同时指定 seq_start。 |
| `strand` | Literal&#91;1, 2&#93; | 1 | NCBI 核酸序列链方向，1 正向、2 反向互补。 |
| `config` | DownloadConfig &#124; None | None | 查询使用 SearchConfig；下载使用 DownloadConfig。 |
| `search_config` | SearchConfig &#124; None | None | 控制查询超时、返回大小和 NCBI 凭据的 SearchConfig。 |
| `api_base_url` | str | &#x27;https://eutils.ncbi.nlm.nih.gov/entrez/eutils&#x27; | 服务 API 根地址；默认值见签名。 |
| `verbose` | bool | True | 是否显示 DNAKit 运行提示；不会改变返回值。 |


## 批量检索并下载 ENCODE 文件 {#dnakit.encode.batch_download-parameters}

`dnakit.encode.batch_download(download_dir, file_format=None, output_type=None, output_category=None, assembly=None, assay_title=None, organism='Homo sapiens', organ=None, biosample_type=None, target=None, preferred_default=None, organize_by='experiment', verify_md5=True, limit=100, dry_run=True, *, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `download_dir` | str | 必填 | Local directory to save files |
| `file_format` | str &#124; None | None | File format filter (&quot;fastq&quot;, &quot;bam&quot;, &quot;bed&quot;, &quot;bigWig&quot;, etc.) |
| `output_type` | str &#124; None | None | Output type filter (&quot;reads&quot;, &quot;peaks&quot;, &quot;signal&quot;, etc.) |
| `output_category` | str &#124; None | None | Output category (&quot;raw data&quot;, &quot;alignment&quot;, &quot;annotation&quot;, etc.) |
| `assembly` | str &#124; None | None | Genome assembly (&quot;GRCh38&quot;, &quot;mm10&quot;, etc.) |
| `assay_title` | str &#124; None | None | Assay type (&quot;Histone ChIP-seq&quot;, &quot;ATAC-seq&quot;, &quot;total RNA-seq&quot;, etc.) |
| `organism` | str | &#x27;Homo sapiens&#x27; | Organism (default: &quot;Homo sapiens&quot;) |
| `organ` | str &#124; None | None | Organ/tissue (&quot;pancreas&quot;, &quot;brain&quot;, &quot;liver&quot;, etc.) |
| `biosample_type` | str &#124; None | None | Biosample type (&quot;tissue&quot;, &quot;cell line&quot;, &quot;primary cell&quot;, etc.) |
| `target` | str &#124; None | None | ChIP/CUT&amp;RUN target (&quot;H3K27me3&quot;, &quot;CTCF&quot;, etc.) |
| `preferred_default` | bool &#124; None | None | If True, only download default/recommended files |
| `organize_by` | Literal&#91;&#x27;flat&#x27;, &#x27;experiment&#x27;, &#x27;format&#x27;, &#x27;experiment_format&#x27;&#93; | &#x27;experiment&#x27; | File organization (&quot;flat&quot;, &quot;experiment&quot;, &quot;format&quot;, &quot;experiment_format&quot;) |
| `verify_md5` | bool | True | Verify downloads with MD5 checksums (default True) |
| `limit` | int | 100 | Max files to download (default 100, safety limit) |
| `dry_run` | bool | True | If True (default), only preview what would be downloaded. Set False to download. |
| `verbose` | bool | True | 是否显示 DNAKit 运行提示；不会改变返回值。 |


## 按文件编号下载 ENCODE 数据 {#dnakit.encode.download_files-parameters}

`dnakit.encode.download_files(file_accessions, download_dir, organize_by='flat', verify_md5=True, *, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `file_accessions` | list&#91;str&#93; | 必填 | List of file accessions to download (e.g., &#91;&quot;ENCFF635JIA&quot;, &quot;ENCFF388RZD&quot;&#93;) |
| `download_dir` | str | 必填 | Local directory path to save files (e.g., &quot;./data/encode&quot;) |
| `organize_by` | Literal&#91;&#x27;flat&#x27;, &#x27;experiment&#x27;, &#x27;format&#x27;, &#x27;experiment_format&#x27;&#93; | &#x27;flat&#x27; | How to organize downloaded files: - &quot;flat&quot;: All files in download_dir (default) - &quot;experiment&quot;: download_dir/ENCSR.../filename - &quot;format&quot;: download_dir/bed/filename - &quot;experiment_format&quot;: download_dir/ENCSR.../bed/filename |
| `verify_md5` | bool | True | Verify file integrity with MD5 checksum (default True) |
| `verbose` | bool | True | 是否显示 DNAKit 运行提示；不会改变返回值。 |
