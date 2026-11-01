# P45 公共数据库操作

通过固定的官方端点执行真实检索和下载。普通测试与下方回放示例不联网；未声称已经逐条完成线上数据验证。

下载继承字节限制、原子写入、校验和与进度条。返回 JSON 可序列化结果，保留服务响应、查询日期和分页信息；分页需按上游令牌显式继续。旧数据集文件名与版本需先通过目录列表确认。

## `encode_body_system_summary()`

- **作用:** Aggregate ENCODE Experiment counts by body system and organ.
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `encode_body_system_summary(*, filters: 'Mapping[str, Any] | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `encode_download_urls()`

- **作用:** Retrieve ENCODE batch-download files.txt URL list for a search or cart.
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `encode_download_urls(*, filters: 'Mapping[str, Any] | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `encode_experiment_matrix()`

- **作用:** Aggregate ENCODE Experiment counts in the assay x biosample matrix.
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `encode_experiment_matrix(*, filters: 'Mapping[str, Any] | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `encode_object_metadata()`

- **作用:** Retrieve ENCODE object metadata by accession or @id, using frame=object.
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `encode_object_metadata(accession: 'str', *, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `encode_tabular_report()`

- **作用:** Return ENCODE report.tsv for selected fields and a bounded or explicit all-results search.
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `encode_tabular_report(fields: 'Sequence[str]', *, filters: 'Mapping[str, Any] | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `rcsb_nucleotide_motif_search()`

- **作用:** Search PDB DNA/RNA chain motifs in simple, prosite or regex notation.
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `rcsb_nucleotide_motif_search(pattern: 'str', *, sequence_type: 'str' = 'dna', pattern_type: 'str' = 'simple', start: 'int' = 0, rows: 'int' = 10, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `rcsb_nucleotide_sequence_search()`

- **作用:** Search PDB DNA or RNA chains by sequence similarity, returning polymer-entity hits.
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `rcsb_nucleotide_sequence_search(sequence: 'str', *, sequence_type: 'str' = 'dna', identity_cutoff: 'float' = 0.0, evalue_cutoff: 'float' = 0.1, start: 'int' = 0, rows: 'int' = 10, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `ucsc_search_assembly()`

- **作用:** Search UCSC track descriptions, names and genomic feature data in one assembly.
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `ucsc_search_assembly(genome: 'str', term: 'str', *, categories: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## `ucsc_track_schema()`

- **作用:** Return original UCSC annotation-track column schema.
- **处理方法:** 调用对应上游 API；下载函数写入实际文件，错误直接抛出。
- **API:** `ucsc_track_schema(genome: 'str', track: 'str', *, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **输入:** 必填标识符、数据集选择与可选参数见签名。参数名与含义遵循上游文档。
- **输出:** `status`、`data` 或下载文件清单，以及来源信息。

## 离线输入检查示例

**示例代码：**

```python
from dnakit.databases.public_apis_p45_functions import ucsc_track_schema
from dnakit.exceptions import ConfigurationError
try:
    ucsc_track_schema("", "knownGene", verbose=False)
except ConfigurationError:
    print("invalid input rejected before network access")
```

**示例输出：**

```text
invalid input rejected before network access
```

## 上游文档

- [官方文档 1](https://genome.ucsc.edu/goldenPath/help/api.html)
- [官方文档 2](https://www.encodeproject.org/help/rest-api/)
- [官方文档 3](https://search.rcsb.org/)
