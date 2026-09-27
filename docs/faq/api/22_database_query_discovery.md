# 查找与调用数据库接口：参数与补充说明

用法和示例见[功能教程](../../api/features/22_database_query_discovery.md)。

各功能按作用、API、输入、输出、示例代码和结果说明介绍。运行前按说明准备凭据、本地文件与可选依赖；远程结果随数据库版本与查询时间变化。

scripts、resource_path、describe、search、list_operations 用于本地发现或资源读取。run、run_support 会执行脚本，可能联网或写文件；解释器、科学软件和参考数据需要按具体任务准备。原有严格 YAML 工作流 run_workflow 与这里的脚本执行入口是不同函数。

## 按操作名称与参数调用数据库功能 {#dnakit.databases.call-parameters}

`dnakit.databases.call(operation, arguments=None, *, max_result_bytes=20000000, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | str | 必填 | 本地操作名或登记的上游别名。 |
| `arguments` | Mapping&#91;str, object&#93; &#124; None | None | 数据库调用传参数字典；工作流传 CLI 字符串列表。 |
| `max_result_bytes` | int | 20000000 | 序列化结果的字节上限，超限抛出异常。 |
| `verbose` | bool | True | 是否显示 DNAKit 运行提示；不会改变返回值。 |


## 查看数据库功能的输入约束与来源 {#dnakit.databases.describe-parameters}

`dnakit.databases.describe(operation, *, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | str | 必填 | 本地操作名或登记的上游别名。 |
| `verbose` | bool | True | 是否显示 DNAKit 运行提示；不会改变返回值。 |


## 使用结构化请求执行数据库操作 {#dnakit.databases.execute-parameters}

`dnakit.databases.execute(request, arguments=None, *, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `request` | str &#124; Mapping&#91;str, object&#93; | 必填 | 操作名，或包含 name、arguments 的字典。 |
| `arguments` | Mapping&#91;str, object&#93; &#124; None | None | 数据库调用传参数字典；工作流传 CLI 字符串列表。 |
| `verbose` | bool | True | 是否显示 DNAKit 运行提示；不会改变返回值。 |


## 按兼容工具名称调用数据库功能 {#dnakit.databases.execute_tool-parameters}

`dnakit.databases.execute_tool(*, tool_name, arguments=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `tool_name` | string | 必填 | Name of the tool to execute (e.g., &#x27;example_tool_name&#x27;) |
| `arguments` | object &#124; None | None | Tool arguments as either (1) a JSON object dictionary or (2) a JSON string that decodes to an object. Object example: {&quot;param1&quot;: &quot;value1&quot;, &quot;param2&quot;: 5}. JSON string example: &quot;{\&quot;param1\&quot;: \&quot;value1\&quot;, \&quot;param2\&quot;: 5}&quot;. {&quot;oneOf&quot;: &#91;{&quot;type&quot;: &quot;object&quot;, &quot;additionalProperties&quot;: true}, {&quot;type&quot;: &quot;string&quot;}&#93;} |
| `verbose` | bool | True | 运行提示开关。 |


## 离线列出数据库操作与兼容名称 {#dnakit.databases.list_operations-parameters}

`dnakit.databases.list_operations(*, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `verbose` | bool | True | 是否显示 DNAKit 运行提示；不会改变返回值。 |


## 按任务关键词离线查找数据库功能 {#dnakit.databases.search-parameters}

`dnakit.databases.search(query, *, limit=10, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `query` | str | 必填 | 本地目录搜索词，匹配名称与描述；使用英文功能词或函数名。 |
| `limit` | int | 10 | 最多返回的匹配条数；数据库目录搜索允许 1–100。 |
| `verbose` | bool | True | 是否显示 DNAKit 运行提示；不会改变返回值。 |
