# P43 公共数据库操作

通过固定的官方端点执行真实检索和下载。普通测试与下方回放示例不联网；未声称已经逐条完成线上数据验证。

下载继承字节限制、原子写入、校验和与进度条。返回 JSON 可序列化结果，保留服务响应、查询日期和分页信息；分页需按上游令牌显式继续。旧数据集文件名与版本需先通过目录列表确认。

## `ncbi_fetch_dna_records()`

- **作用:** Fetch nucleotide FASTA/gbc/ft/native/CDS, Gene XML/table/ASN.1, SNP, ClinVar or SRA records.
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `ncbi_fetch_dna_records(database: 'str', ids: 'Sequence[str]', *, rettype: 'str' = 'fasta', retmode: 'str' = 'text', email: 'str | None' = None, api_key: 'str | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `ncbi_link_dna_records()`

- **作用:** ELink one of the 28 curated nuccore/assembly/gene/variant/SRA/GEO/OMIM relationships.
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `ncbi_link_dna_records(linkname: 'str', ids: 'Sequence[str]', *, email: 'str | None' = None, api_key: 'str | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `ncbi_search_dna_database()`

- **作用:** ESearch assembly, annotinfo, snp, dbvar, grasp, geoprofiles, gtr or seqannot; return IDs
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `ncbi_search_dna_database(database: 'str', term: 'str', *, retmax: 'int' = 20, retstart: 'int' = 0, email: 'str | None' = None, api_key: 'str | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `ncbi_summarize_dna_records()`

- **作用:** Return unmodified ESummary JSON records from snp, sra or omim.
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `ncbi_summarize_dna_records(database: 'str', ids: 'Sequence[str]', *, email: 'str | None' = None, api_key: 'str | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## 离线输入检查示例

**示例代码：**

```python
from dnakit.databases.ncbi_eutils_p43_functions import ncbi_fetch_dna_records
from dnakit.exceptions import ConfigurationError
try:
    ncbi_fetch_dna_records("protein", ["1"], verbose=False)
except ConfigurationError:
    print("invalid input rejected before network access")
```

**示例输出：**

```text
invalid input rejected before network access
```

## 上游文档

- [官方文档 1](https://www.ncbi.nlm.nih.gov/books/NBK25499/)
