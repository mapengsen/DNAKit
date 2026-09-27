# 查找与调用数据库接口

按名称或关键词查找并调用数据库接口。

<!-- BEGIN GENERATED API: discovery -->

<span id="discovery"></span>

## 1) 按操作名称与参数调用数据库功能 {#dnakit.databases.call}

- **作用：** 根据操作名和参数字典调用数据库接口。
- **处理方法：** 校验操作名和参数后调用对应接口，保留提供方返回的字段。
- **API：** `dnakit.databases.call(operation, arguments=None, *, max_result_bytes=20000000, verbose=True)`
- **输入：** 必须提供 `operation`。
- **输出：** 返回所选操作的结果；本例为 BRCA1 的 Ensembl 查询记录。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.call(
    operation='ensembl_lookup_gene',
    arguments={'id': 'BRCA1', 'species': 'homo_sapiens'},
    verbose=False,
)
print(result)
```

- **示例结果：** 返回所选操作的结果；本例为 BRCA1 的 Ensembl 查询记录。

[完整参数与使用说明](../../faq/api/22_database_query_discovery.md#dnakit.databases.call-parameters)

## 2) 查看数据库功能的输入约束与来源 {#dnakit.databases.describe}

- **作用：** 查看单个接口的实际输入 Schema、上游别名和实现信息。
- **处理方法：** 从本地接口目录读取参数约束与来源信息，不联网。
- **API：** `dnakit.databases.describe(operation, *, verbose=True)`
- **输入：** 必须提供 `operation`。
- **输出：** 字典；input_schema 为实际输入约束，config 保存提供方配置，upstream_name 为原接口名。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.describe(
    operation='gnomad_get_sv_by_region',
    verbose=False,
)
print(result)
```

- **示例结果：** 字典；input_schema 为实际输入约束，config 保存提供方配置，upstream_name 为原接口名。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "name": "<str>",
  "upstream_name": "<str>",
  "input_schema": "<dict>",
  "config": "<dict>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_discovery.md#dnakit.databases.describe-parameters)

## 3) 使用结构化请求执行数据库操作 {#dnakit.databases.execute}

- **作用：** 使用操作名或结构化请求分发一次调用。
- **处理方法：** 解析操作名或 name/arguments 请求，再调用对应数据库接口。
- **API：** `dnakit.databases.execute(request, arguments=None, *, verbose=True)`
- **输入：** 必须提供 `request`。
- **输出：** 返回请求中所选操作的结果；本例为 BRCA1 的 Ensembl 查询记录。
- **示例代码：**

```python
import dnakit

example_request = {'name': 'ensembl_lookup_gene',
 'arguments': {'id': 'BRCA1', 'species': 'homo_sapiens'}}

result = dnakit.databases.execute(
    request=example_request,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回请求中所选操作的结果；本例为 BRCA1 的 Ensembl 查询记录。

[完整参数与使用说明](../../faq/api/22_database_query_discovery.md#dnakit.databases.execute-parameters)

## 4) 按兼容工具名称调用数据库功能 {#dnakit.databases.execute_tool}

- **作用：** 根据工具名与参数调用本地 DNAKit 数据库操作；禁止递归调用自身。
- **处理方法：** 通过 `Execute` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.execute_tool(*, tool_name, arguments=None, verbose=True)`
- **输入：** 必须提供 `tool_name`。
- **输出：** 返回 tool_name 所选接口的结果，字段由该接口定义。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.execute_tool(
    tool_name='list_tools',
    arguments={'mode': 'names', 'limit': '3'},
    verbose=False,
)
print(result)
```

- **示例结果：** 返回 tool_name 所选接口的结果，字段由该接口定义。

[完整参数与使用说明](../../faq/api/22_database_query_discovery.md#dnakit.databases.execute_tool-parameters)

## 5) 离线列出数据库操作与兼容名称 {#dnakit.databases.list_operations}

- **作用：** 离线列出全部数据库操作名与上游别名。
- **处理方法：** 读取本地登记目录，列出函数名及上游别名。
- **API：** `dnakit.databases.list_operations(*, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 列表，每项包含 name、upstream_name。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.list_operations(
    verbose=False,
)
print(result)
```

- **示例结果：** 列表，每项包含 name、upstream_name。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "name": "<str>",
    "upstream_name": "<str>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_discovery.md#dnakit.databases.list_operations-parameters)

## 6) 按任务关键词离线查找数据库功能 {#dnakit.databases.search}

- **作用：** 按名称、任务关键词或描述离线搜索数据库接口。
- **处理方法：** 在本地接口名称和描述中匹配关键词，按匹配数排序。
- **API：** `dnakit.databases.search(query, *, limit=10, verbose=True)`
- **输入：** 必须提供 `query`。
- **输出：** 列表，每项包含 name、description、upstream_name。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.search(
    query='gnomad',
    limit=5,
    verbose=False,
)
print(result)
```

- **示例结果：** 列表，每项包含 name、description、upstream_name。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "name": "<str>",
    "description": "<str>",
    "upstream_name": "<str>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_discovery.md#dnakit.databases.search-parameters)

<!-- END GENERATED API: discovery -->
