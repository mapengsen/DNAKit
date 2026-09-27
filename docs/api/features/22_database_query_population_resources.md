# 人群、样本与连锁查询

依赖与返回值见[FAQ](../../faq/tutorial_integrations.md#optional-bioinformatics)。

<!-- BEGIN GENERATED API: population-resources -->

<span id="population-resources"></span>

## 1) 列出 IGSR 与千人基因组数据集合 {#dnakit.databases.igsr_list_data_collections}

- **作用：** 列出 IGSR 与千人基因组数据集合。
- **处理方法：** 通过 `IGSR` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.igsr_list_data_collections(*, operation=None, limit=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `total`, `returned`, `collections` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.igsr_list_data_collections(
    operation='list_data_collections',
    limit=20,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "total": "<integer>",
  "returned": "<integer>",
  "collections": [
    {
      "code": "<string>",
      "title": "<string>",
      "short_title": "<string>",
      "sample_count": "<integer>",
      "population_count": "<integer>",
      "data_types": "<array>",
      "website": "<string>"
    }
  ]
}
```

[完整参数与使用说明](../../faq/api/22_database_query_population_resources.md#dnakit.databases.igsr_list_data_collections-parameters)

## 2) 检索千人基因组人群信息 {#dnakit.databases.igsr_search_populations}

- **作用：** 检索千人基因组人群信息。
- **处理方法：** 通过 `IGSR` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.igsr_search_populations(*, operation=None, superpopulation=None, query=None, limit=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `total`, `returned`, `populations` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.igsr_search_populations(
    operation='search_populations',
    superpopulation='EAS',
    limit=10,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "total": "<integer>",
  "returned": "<integer>",
  "populations": [
    {
      "code": "<string>",
      "name": "<string>",
      "description": "<string>",
      "sample_count": "<integer>",
      "superpopulation_code": "<string>",
      "superpopulation_name": "<string>",
      "latitude": "<number>",
      "longitude": "<number>"
    }
  ]
}
```

[完整参数与使用说明](../../faq/api/22_database_query_population_resources.md#dnakit.databases.igsr_search_populations-parameters)

## 3) 按人群或集合检索千人基因组样本 {#dnakit.databases.igsr_search_samples}

- **作用：** 按人群或集合检索千人基因组样本。
- **处理方法：** 通过 `IGSR` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.igsr_search_samples(*, operation=None, population=None, data_collection=None, sample_name=None, limit=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `total`, `returned`, `samples` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.igsr_search_samples(
    operation='search_samples',
    population='GBR',
    limit=5,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "total": "<integer>",
  "returned": "<integer>",
  "samples": [
    {
      "name": "<string>",
      "sex": "<string>",
      "biosample_id": "<string>",
      "populations": "<array>",
      "data_collections": "<array>"
    }
  ]
}
```

[完整参数与使用说明](../../faq/api/22_database_query_population_resources.md#dnakit.databases.igsr_search_samples-parameters)

## 4) 查询 SNP 的群体特异性连锁不平衡代理变异 {#dnakit.databases.ldlink_get_proxies}

- **作用：** 查询 SNP 的群体特异性连锁不平衡代理变异。
- **处理方法：** 通过 `LDlink` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.ldlink_get_proxies(*, variant, population=None, r2_threshold=None, genome_build=None, limit=None, verbose=True)`
- **输入：** 必须提供 `variant`。
- **输出：** 返回记录列表；每条记录包含 `rsid`, `r2`, `coord`, `alleles` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.ldlink_get_proxies(
    variant='rs7903146',
    population='CEU',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "rsid": "<string>",
    "r2": "<number>",
    "coord": "<string>",
    "alleles": "<string>",
    "maf": "<string>",
    "distance": "<string>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_population_resources.md#dnakit.databases.ldlink_get_proxies-parameters)

<!-- END GENERATED API: population-resources -->
