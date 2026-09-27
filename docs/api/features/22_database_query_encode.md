# ENCODE 实验与文件查询

查询 ENCODE 实验、样本和可下载文件。

## 1) 筛选 ENCODE 对象与文件 {#search-encode-search}

- **作用：** 按对象类型、关键词和字段条件查询 ENCODE 记录，例如已发布的 bigWig 文件。
- **计算方法：** 提交 ENCODE 搜索条件，将匹配对象及文件信息整理为 QueryResult。
- **API：** `dnakit.search.encode_search(object_type[可选], search_term[可选], filters[可选], limit[可选], config[可选], verbose[可选])`
- **输入：** 可选对象类型、搜索词、过滤字典、记录上限和 SearchConfig。需要联网。
- **输出：** 返回 `QueryResult`，包含查询记录 `records`、数据源 `provider` 和分页信息。
- **示例代码：**

```python
from dnakit.search import encode_search

result = encode_search(
    object_type="File",
    filters={"file_format": "bigWig", "status": "released"},
    limit=20,
    verbose=False,
)
print(result.query_type, result.provider, len(result.records))
```

- **示例结果：** 以下为返回格式，未执行联网查询；实际记录数由在线数据库决定。

```text
regulation ENCODE <返回记录数>
```

<!-- BEGIN GENERATED API: encode -->

<span id="encode"></span>

## 2) 按编号查询 ENCODE 实验详情 {#dnakit.encode.get_experiment}

- **作用：** 获取指定实验的完整信息。
- **处理方法：** 按指定条件访问 ENCODE，解析实验、文件或下载结果。
- **API：** `dnakit.encode.get_experiment(accession, *, verbose=True)`
- **输入：** 必须提供 `accession`。
- **输出：** 返回 JSON 文本，包含实验详情与文件列表。
- **示例代码：**

```python
import dnakit

result = dnakit.encode.get_experiment(
    accession='ENCSR133RZO',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回 JSON 文本，包含实验详情与文件列表。

[完整参数与使用说明](../../faq/api/22_database_query_encode.md#dnakit.encode.get_experiment-parameters)

## 3) 查询 ENCODE 筛选条件的实时记录数 {#dnakit.encode.get_facets}

- **作用：** 查询可用于筛选的字段取值及数量。
- **处理方法：** 按指定条件访问 ENCODE，解析实验、文件或下载结果。
- **API：** `dnakit.encode.get_facets(search_type='Experiment', assay_title=None, organism=None, organ=None, biosample_type=None, *, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回 JSON 文本，包含筛选字段及各取值的记录数。
- **示例代码：**

```python
import dnakit

result = dnakit.encode.get_facets(
    assay_title='ATAC-seq',
    organism='Homo sapiens',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回 JSON 文本，包含筛选字段及各取值的记录数。

[完整参数与使用说明](../../faq/api/22_database_query_encode.md#dnakit.encode.get_facets-parameters)

## 4) 查询 ENCODE 文件详情 {#dnakit.encode.get_file_info}

- **作用：** 获取指定文件的详细元数据。
- **处理方法：** 按指定条件访问 ENCODE，解析实验、文件或下载结果。
- **API：** `dnakit.encode.get_file_info(accession, *, verbose=True)`
- **输入：** 必须提供 `accession`。 先用文件查询功能取得实际编号，替换 `ENCFF_REPLACE_ME`。
- **输出：** 返回 JSON 文本，包含文件格式、大小、下载地址、MD5 和组装版本。
- **示例代码：**

```python
import dnakit

result = dnakit.encode.get_file_info(
    accession='ENCFF_REPLACE_ME',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回 JSON 文本，包含文件格式、大小、下载地址、MD5 和组装版本。

[完整参数与使用说明](../../faq/api/22_database_query_encode.md#dnakit.encode.get_file_info-parameters)

## 5) 获取 ENCODE 检索可用的筛选值 {#dnakit.encode.get_metadata}

- **作用：** 获取实验的结构化元数据。
- **处理方法：** 按指定条件访问 ENCODE，解析实验、文件或下载结果。
- **API：** `dnakit.encode.get_metadata(metadata_type, *, verbose=True)`
- **输入：** 必须提供 `metadata_type`。
- **输出：** 返回 JSON 文本，列出指定元数据类别的可用取值。
- **示例代码：**

```python
import dnakit

result = dnakit.encode.get_metadata(
    metadata_type='assays',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回 JSON 文本，列出指定元数据类别的可用取值。

[完整参数与使用说明](../../faq/api/22_database_query_encode.md#dnakit.encode.get_metadata-parameters)

## 6) 列出指定 ENCODE 实验的文件 {#dnakit.encode.list_files}

- **作用：** 列出实验关联的文件。
- **处理方法：** 按指定条件访问 ENCODE，解析实验、文件或下载结果。
- **API：** `dnakit.encode.list_files(experiment_accession, file_format=None, file_type=None, output_type=None, output_category=None, assembly=None, status=None, preferred_default=None, limit=200, *, verbose=True)`
- **输入：** 必须提供 `experiment_accession`。
- **输出：** 返回 JSON 文本，列出文件 accession、格式、大小、下载地址及元数据。
- **示例代码：**

```python
import dnakit

result = dnakit.encode.list_files(
    experiment_accession='ENCSR133RZO',
    file_format='bed',
    limit=3,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回 JSON 文本，列出文件 accession、格式、大小、下载地址及元数据。

[完整参数与使用说明](../../faq/api/22_database_query_encode.md#dnakit.encode.list_files-parameters)

## 7) 按实验条件检索 ENCODE 数据 {#dnakit.encode.search_experiments}

- **作用：** 按实验类型、物种、组织等条件检索实验。
- **处理方法：** 按指定条件访问 ENCODE，解析实验、文件或下载结果。
- **API：** `dnakit.encode.search_experiments(assay_title=None, organism='Homo sapiens', organ=None, biosample_type=None, biosample_term_name=None, target=None, status='released', lab=None, award=None, assembly=None, replication_type=None, life_stage=None, sex=None, treatment=None, genetic_modification=None, perturbed=None, search_term=None, date_released_from=None, date_released_to=None, limit=25, offset=0, *, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回 JSON 文本，包含实验结果、总数和分页信息。
- **示例代码：**

```python
import dnakit

result = dnakit.encode.search_experiments(
    assay_title='ATAC-seq',
    organism='Homo sapiens',
    limit=3,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回 JSON 文本，包含实验结果、总数和分页信息。

[完整参数与使用说明](../../faq/api/22_database_query_encode.md#dnakit.encode.search_experiments-parameters)

## 8) 联合实验与文件条件检索 ENCODE 文件 {#dnakit.encode.search_files}

- **作用：** 按文件类型、组装和实验条件检索文件。
- **处理方法：** 按指定条件访问 ENCODE，解析实验、文件或下载结果。
- **API：** `dnakit.encode.search_files(file_format=None, file_type=None, output_type=None, output_category=None, assembly=None, assay_title=None, organism=None, organ=None, biosample_type=None, target=None, status='released', preferred_default=None, search_term=None, limit=25, offset=0, *, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回 JSON 文本，包含文件结果、总数和分页信息。
- **示例代码：**

```python
import dnakit

result = dnakit.encode.search_files(
    file_format='bed',
    assembly='GRCh38',
    assay_title='ATAC-seq',
    limit=3,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回 JSON 文本，包含文件结果、总数和分页信息。

[完整参数与使用说明](../../faq/api/22_database_query_encode.md#dnakit.encode.search_files-parameters)

<!-- END GENERATED API: encode -->
