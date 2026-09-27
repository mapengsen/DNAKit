# 工程化与扩展：参数与补充说明

用法和示例见[功能教程](../../api/features/14_engineering.md)。

各功能按作用、API、输入、输出、示例代码和结果说明介绍。运行前按说明准备凭据、本地文件与可选依赖；远程结果随数据库版本与查询时间变化。

ENCODE 接口返回文本：多数是 JSON，`get_citations` 可返回 BibTeX/RIS，`export_data` 可返回 CSV/TSV。仅 JSON 文本使用 `json.loads()`。追踪与凭据保存在 `~/.dnakit-encode/`；持久化凭据功能需要 `dnakit[encode-credentials]`。

各功能按作用、API、输入、输出、示例代码和结果说明介绍。运行前按说明准备凭据、本地文件与可选依赖；远程结果随数据库版本与查询时间变化。

scripts、resource_path、describe、search、list_operations 用于本地发现或资源读取。run、run_support 会执行脚本，可能联网或写文件；解释器、科学软件和参考数据需要按具体任务准备。原有严格 YAML 工作流 run_workflow 与这里的脚本执行入口是不同函数。

## 读取已执行 Notebook 的单元格与输出 {#dnakit.databases.read_executed_notebook-parameters}

`dnakit.databases.read_executed_notebook(*, data_folder=None, notebook_path=None, search=None, regex=None, max_output_chars=None, include_source=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `data_folder` | string &#124; null | None | Directory to search (recursively) for a .ipynb file. Prefers *_executed.ipynb over regular .ipynb. Provide either this OR notebook_path. |
| `notebook_path` | string &#124; null | None | Direct path to a specific .ipynb file. Provide either this OR data_folder. |
| `search` | string | None | Optional: comma-separated terms (or regex if regex=true). Cells whose source or output contains any term are returned in &#x27;matching_cells&#x27;. Case-insensitive. Example: &#x27;upregulated,log2FoldChange&#x27; or &#x27;padj.*0\.05&#x27; |
| `regex` | boolean | None | Treat &#x27;search&#x27; as a Python regex instead of a comma-separated term list (default: false) |
| `max_output_chars` | integer | None | Truncate each cell source/output to this many characters (default 400) |
| `include_source` | boolean | None | Include cell source code in the response (default true). Set false to halve response size when you only need outputs. |
| `verbose` | bool | True | 运行提示开关。 |


## 管理 ENCODE API 凭据 {#dnakit.encode.manage_credentials-parameters}

`dnakit.encode.manage_credentials(action, access_key=None, secret_key=None, *, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `action` | Literal&#91;&#x27;store&#x27;, &#x27;check&#x27;, &#x27;clear&#x27;&#93; | 必填 | What to do: - &quot;store&quot;: Save new credentials (requires access_key and secret_key) - &quot;check&quot;: Check if credentials are configured - &quot;clear&quot;: Remove stored credentials |
| `access_key` | str &#124; None | None | Your ENCODE access key (only for action=&quot;store&quot;) |
| `secret_key` | str &#124; None | None | Your ENCODE secret key (only for action=&quot;store&quot;) |
| `verbose` | bool | True | 是否显示 DNAKit 运行提示；不会改变返回值。 |


## 定位并校验随包资源文件 {#dnakit.workflows.resource_path-parameters}

`dnakit.workflows.resource_path(skill_id, resource, *, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `skill_id` | str | 必填 | 已收录技能名称；兼容登记过的旧名称。 |
| `resource` | str | 必填 | 技能包内相对路径，读取前校验哈希。 |
| `verbose` | bool | True | 是否显示 DNAKit 运行提示；不会改变返回值。 |


## 执行随包工作流并获取退出码与输出 {#dnakit.workflows.run-parameters}

`dnakit.workflows.run(skill_id, script, *, arguments=None, workdir=None, timeout_seconds=300.0, max_output_bytes=1000000, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `skill_id` | str | 必填 | 已收录技能名称；兼容登记过的旧名称。 |
| `script` | str | 必填 | 技能包内 .py、.R 或 .sh 脚本相对路径。 |
| `arguments` | list&#91;str&#93; &#124; None | None | 数据库调用传参数字典；工作流传 CLI 字符串列表。 |
| `workdir` | str &#124; Path &#124; None | None | 子进程工作目录；None 使用当前目录。 |
| `timeout_seconds` | float | 300.0 | 子进程执行超时秒数。 |
| `max_output_bytes` | int | 1000000 | 最多捕获的子进程输出字节数。 |
| `verbose` | bool | True | 是否显示 DNAKit 运行提示；不会改变返回值。 |


## 执行组合流程的支持脚本 {#dnakit.workflows.run_support-parameters}

`dnakit.workflows.run_support(name, *, arguments=None, workdir=None, timeout_seconds=300, max_output_bytes=1000000, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `name` | str | 必填 | 支持流程名称，允许值见本条说明。 |
| `arguments` | list&#91;str&#93; &#124; None | None | 数据库调用传参数字典；工作流传 CLI 字符串列表。 |
| `workdir` | str &#124; Path &#124; None | None | 子进程工作目录；None 使用当前目录。 |
| `timeout_seconds` | float | 300 | 子进程执行超时秒数。 |
| `max_output_bytes` | int | 1000000 | 最多捕获的子进程输出字节数。 |
| `verbose` | bool | True | 是否显示 DNAKit 运行提示；不会改变返回值。 |


## 列出随包 Python、R 与 Shell 脚本 {#dnakit.workflows.scripts-parameters}

`dnakit.workflows.scripts(skill_id, *, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `skill_id` | str | 必填 | 已收录技能名称；兼容登记过的旧名称。 |
| `verbose` | bool | True | 是否显示 DNAKit 运行提示；不会改变返回值。 |


## 将已有 Markdown 报告渲染为 PDF {#dnakit.workflows.write_wes_report-parameters}

`dnakit.workflows.write_wes_report(input_report, output_dir, *, language='en', title=None, logo_left=None, logo_right=None, font_path=None, overwrite=False, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `input_report` | str &#124; Path | 必填 | 已有 UTF-8 Markdown 报告路径，最多 5 MB。 |
| `output_dir` | str &#124; Path | 必填 | 报告输出目录，生成 PDF 和来源清单。 |
| `language` | str | &#x27;en&#x27; | 报告标签语言 en 或 es；不翻译输入正文。 |
| `title` | str &#124; None | None | 自定义报告标题；None 使用语言对应默认标题。 |
| `logo_left` | str &#124; Path &#124; None | None | 左侧 Logo 图片路径。 |
| `logo_right` | str &#124; Path &#124; None | None | 右侧 Logo 图片路径。 |
| `font_path` | str &#124; Path &#124; None | None | 可选 TTF 字体路径。 |
| `overwrite` | bool | False | 是否覆盖已存在的报告文件。 |
| `verbose` | bool | True | 是否显示 DNAKit 运行提示；不会改变返回值。 |
