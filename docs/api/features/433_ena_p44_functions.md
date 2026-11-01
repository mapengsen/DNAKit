# P44 公共数据库操作

通过固定的官方端点执行真实检索和下载。普通测试与下方回放示例不联网；未声称已经逐条完成线上数据验证。

下载继承字节限制、原子写入、校验和与进度条。返回 JSON 可序列化结果，保留服务响应、查询日期和分页信息；分页需按上游令牌显式继续。旧数据集文件名与版本需先通过目录列表确认。

## `ena_count_keyword_results()`

- **作用:** Count ENA keyword-search hits for a record format and domain.
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `ena_count_keyword_results(query: 'str', domain: 'str', result: 'str', *, format: 'str' = 'embl', verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `ena_download_cram_reference()`

- **作用:** Download exact reference sequence bytes by CRAM registry MD5 or SHA1 checksum.
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `ena_download_cram_reference(checksum: 'str', output_path: 'str', *, algorithm: 'str' = 'md5', overwrite: 'bool' = False, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `ena_download_keyword_results()`

- **作用:** Download keyword-search EMBL/FASTA/XML records, including WGS/TSA/TLS master sets.
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `ena_download_keyword_results(query: 'str', domain: 'str', result: 'str', output_path: 'str', *, format: 'str' = 'embl', limit: 'int' = 100000, offset: 'int' = 0, overwrite: 'bool' = False, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `ena_download_linked_records()`

- **作用:** Download EMBL/FASTA/XML records linked to a study, sample or taxon accession.
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `ena_download_linked_records(link_type: 'str', accession: 'str', result: 'str', output_path: 'str', *, format: 'str' = 'embl', subtree: 'bool' = False, overwrite: 'bool' = False, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `ena_download_record()`

- **作用:** Download EMBL contig-set records or XML assembly records by accession.
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `ena_download_record(accession: 'str', output_path: 'str', *, format: 'str' = 'embl', full_contig_set: 'bool' = False, overwrite: 'bool' = False, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `ena_download_search()`

- **作用:** Download advanced-search EMBL/FASTA sequences or XML run/experiment/analysis/assembly
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `ena_download_search(query: 'str', result: 'str', output_path: 'str', *, format: 'str' = 'embl', limit: 'int' = 100000, overwrite: 'bool' = False, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `ena_linked_accessions()`

- **作用:** List portal DNA record accessions linked to a study, sample or taxon.
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `ena_linked_accessions(link_type: 'str', accession: 'str', result: 'str', *, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `ena_search_cross_references()`

- **作用:** Search ENA cross-references to an external database and return original TSV.
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `ena_search_cross_references(source: 'str', *, source_accession: 'str | None' = None, target: 'str | None' = None, target_accession: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## 离线输入检查示例

**示例代码：**

```python
from dnakit.download.ena_p44_functions import ena_download_record
from dnakit.exceptions import ConfigurationError
try:
    ena_download_record("X", "unused", format="bad", verbose=False)
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
