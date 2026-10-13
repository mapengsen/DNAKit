# 流水线定义工具：WDL、CWL、Snakemake、nf-core 与 Galaxy

这一组函数是 `dnakit.workflows` 中基于可选依赖的适配器，分别使用 [miniwdl](https://github.com/chanzuckerberg/miniwdl)、[cwltool](https://github.com/common-workflow-language/cwltool)、[Snakemake](https://github.com/snakemake/snakemake)、[nf-core tools](https://github.com/nf-core/tools) 和 [bioblend](https://github.com/galaxyproject/bioblend)。它们是流水线定义工具，不是 DNA 算法：只解析和校验 WDL、CWL、Snakefile 和 nf-core 检出目录，列出任务、步骤、规则或进程，做 Snakemake 的 dry-run 计划，并通过 Galaxy API 读取工具、历史和工作流列表。它们不读取或计算任何序列结果，也不执行流水线。

依赖都是可选的，只有调用对应函数时才会导入。缺少包时会报错并给出安装命令，对应的 extra 为 `bio-wdl`、`bio-cwl`、`bio-snakemake`、`bio-nf-core` 和 `bio-bioblend`。

**示例文件：** `pipe.wdl`
```text
version 1.0
task Count {
  input { File bam }
  command <<< wc -l ~{bam} >>>
  output { Int lines = 1 }
}
workflow Pipeline {
  input {
    File reads
    String sample
  }
  call Count { input: bam = reads }
  output { Int lines = Count.lines }
}
```

## 1) `dnakit.workflows.parse_wdl_document()`

- **作用：** 解析 WDL 文档，列出任务、输入、输出和 workflow 中的 call。
- **处理方法：** 用 miniwdl 的 `WDL.load()` 读取文件并做类型检查，不运行任何任务。没有默认值的输入记为必填。任务和 call 按名称排序。
- **API：** `dnakit.workflows.parse_wdl_document(source, *, max_bytes=1000000)`。
- **输入：** WDL 文件路径；`max_bytes` 是文件大小上限，超过时报错。导入的文件由 miniwdl 读取。
- **输出：** 含 `wdl_version`、`imports`、`tasks`、`workflow` 和 `provenance` 的字典。
- **示例代码：**

```python
from dnakit.workflows import parse_wdl_document

result = parse_wdl_document("pipe.wdl")
print([task["name"] for task in result["tasks"]])
print([(item["name"], item["required"]) for item in result["workflow"]["inputs"]])
print(result["workflow"]["calls"])
```

- **示例输出：**

```text
['Count']
[('reads', True), ('sample', True)]
[{'name': 'Count', 'task': 'Count'}]
```

## 2) `dnakit.workflows.check_wdl_inputs()`

- **作用：** 检查一组 JSON 输入是否符合 WDL workflow 的输入声明。
- **处理方法：** 比较输入名称与必填项，再用 miniwdl 的 `Value.from_json()` 检查每个值能否构造成声明的类型。`File` 类型只检查字符串形式，不检查文件是否存在。
- **API：** `dnakit.workflows.check_wdl_inputs(source, inputs, *, max_bytes=1000000)`。
- **输入：** WDL 文件路径；`inputs` 是 `{输入名: 值}` 的映射，名称可以是 `调用名.输入名`。
- **输出：** `valid`、`missing`（缺少的必填输入）、`unknown`（未声明的输入）、`invalid`（类型不符的输入）和 `provenance`。
- **示例代码：**

```python
from dnakit.workflows import check_wdl_inputs

result = check_wdl_inputs("pipe.wdl", {"reads": 3, "sample": "S1", "extra": 1})
print(result["valid"], result["missing"], result["unknown"])
print([item["name"] for item in result["invalid"]])
```

- **示例输出：**

```text
False [] ['extra']
['reads']
```

## 3) `dnakit.workflows.read_cwl_workflow()`

- **作用：** 用 cwltool 加载 CWL Workflow，列出输入、步骤和输出。
- **处理方法：** 调用 cwltool 的 `load_tool()`，由 schema-salad 做校验；不执行任何步骤。步骤按 id 排序，因为 cwltool 的内部顺序依赖 `PYTHONHASHSEED`。每个步骤列出输入的数据来源和输出名。
- **API：** `dnakit.workflows.read_cwl_workflow(source, *, max_bytes=1000000)`。
- **输入：** CWL 工作流文件路径，`class` 必须是 `Workflow`；`run:` 引用的工具文件按工作流文件所在目录解析。
- **输出：** 含 `inputs`、`steps`、`outputs` 和 `provenance` 的字典。`steps` 的每一项包含 `name`、`run`、`run_class`、`inputs` 和 `outputs`。
- **示例文件：** `wf.cwl`
```text
cwlVersion: v1.2
class: Workflow
inputs:
  bam: File
outputs:
  lines:
    type: File
    outputSource: count/lines
steps:
  count:
    run: count.cwl
    in:
      bam: bam
    out: [lines]
```

- **示例文件：** `count.cwl`
```text
cwlVersion: v1.2
class: CommandLineTool
baseCommand: [wc, -l]
inputs:
  bam:
    type: File
    inputBinding: {position: 1}
outputs:
  lines:
    type: stdout
stdout: count.txt
```

- **示例代码：**

```python
from dnakit.workflows import read_cwl_workflow

result = read_cwl_workflow("wf.cwl")
step = result["steps"][0]
print(step["name"], step["run"], step["run_class"])
print(step["inputs"], result["outputs"][0]["sources"])
```

- **示例输出：**

```text
count count.cwl CommandLineTool
[{'name': 'bam', 'sources': ['bam']}] ['count/lines']
```

## 4) `dnakit.workflows.read_snakefile()`

- **作用：** 解析 Snakefile，列出规则及其输入和输出。
- **处理方法：** 通过 Snakemake API 解析文件和 `include`，不运行任何规则。解析在临时目录中进行，因此不会在源目录中写入 `.snakemake/`。规则按名称排序。
- **API：** `dnakit.workflows.read_snakefile(source, *, max_bytes=1000000)`。
- **输入：** Snakefile 路径；`max_bytes` 是文件大小上限。
- **输出：** 含 `default_target`、`rules`（每项为 `name`、`inputs`、`outputs`）和 `provenance` 的字典。
- **示例文件：** `Snakefile`
```text
rule all:
    input:
        expand("counts/{sample}.txt", sample=["A", "B"])

rule align:
    input:
        "reads/{sample}.fastq"
    output:
        "bam/{sample}.bam"
    shell:
        "cp {input} {output}"

rule count:
    input:
        "bam/{sample}.bam"
    output:
        "counts/{sample}.txt"
    shell:
        "wc -l {input} > {output}"
```

- **示例代码：**

```python
from dnakit.workflows import read_snakefile

result = read_snakefile("Snakefile")
print(result["default_target"])
for rule in result["rules"]:
    print(rule["name"], rule["inputs"], rule["outputs"])
```

- **示例输出：**

```text
all
align ['reads/{sample}.fastq'] ['bam/{sample}.bam']
all ['counts/A.txt', 'counts/B.txt'] []
count ['bam/{sample}.bam'] ['counts/{sample}.txt']
```

## 5) `dnakit.workflows.plan_snakemake_jobs()`

- **作用：** 对 Snakefile 做 dry-run，列出为达到目标需要运行的作业。
- **处理方法：** 用 Snakemake API 构建 DAG，只保留需要运行的作业（已存在的输出不再列出）。不运行任何规则，也不获取锁。`workdir` 默认是 Snakefile 所在目录，输入文件相对于它查找；Snakemake 会在该目录中写入 `.snakemake/` 文件夹和日志。作业按规则名、输出和输入排序。
- **API：** `dnakit.workflows.plan_snakemake_jobs(source, *, targets=(), workdir=None, max_bytes=1000000, max_jobs=10000)`。
- **输入：** Snakefile 路径；`targets` 是目标文件或规则名列表，留空时使用默认目标；`max_jobs` 限制计划中的作业数。
- **输出：** 含 `job_count`、`jobs_by_rule`、`jobs`（每项为 `rule`、`inputs`、`outputs`、`wildcards`）和 `provenance` 的字典。
- **示例代码：**

```python
from pathlib import Path
from dnakit.workflows import plan_snakemake_jobs

Path("reads").mkdir(exist_ok=True)
Path("bam").mkdir(exist_ok=True)
for sample in ("A", "B"):
    Path(f"reads/{sample}.fastq").write_text("ACGT\n")
Path("bam/A.bam").write_text("bam\n")

result = plan_snakemake_jobs("Snakefile")
print(result["job_count"], result["jobs_by_rule"])
for job in result["jobs"]:
    print(job["rule"], job["outputs"])
```

- **示例输出：**

```text
4 {'align': 1, 'all': 1, 'count': 2}
align ['bam/B.bam']
all []
count ['counts/A.txt']
count ['counts/B.txt']
```

## 6) `dnakit.workflows.list_nfcore_processes()`

- **作用：** 列出 nf-core 流水线检出目录中的 Nextflow 进程。
- **处理方法：** 用 nf_core 的 `NFCoreComponent` 解析 `modules/nf-core` 和 `modules/local` 中的 `main.nf`，以及这两个目录中直接存放的 `.nf` 文件。全过程离线，不访问 GitHub，也不克隆仓库。只返回进程名和路径；容器和输入输出需要 `nextflow` 命令，因此不读取。
- **API：** `dnakit.workflows.list_nfcore_processes(source, *, max_components=500, max_bytes=1000000)`。
- **输入：** 流水线目录，必须包含 `modules/` 子目录；`max_components` 限制组件数量，`max_bytes` 限制每个 `main.nf` 的大小。
- **输出：** 含 `component_count`、`process_count`、`components`（每项为 `name`、`kind`、`process`、`path`）和 `provenance` 的字典。`kind` 为 `nf-core` 或 `local`。
- **示例代码：**

```python
from pathlib import Path
from dnakit.workflows import list_nfcore_processes

Path("pipe/modules/nf-core/fastqc").mkdir(parents=True, exist_ok=True)
Path("pipe/modules/nf-core/fastqc/main.nf").write_text(
    'process FASTQC {\n    script:\n    """\n    fastqc\n    """\n}\n'
)
Path("pipe/modules/local").mkdir(parents=True, exist_ok=True)
Path("pipe/modules/local/count_lines.nf").write_text(
    'process COUNT_LINES {\n    script:\n    """\n    wc -l\n    """\n}\n'
)

result = list_nfcore_processes("pipe")
print(result["process_count"])
for row in result["components"]:
    print(row["kind"], row["name"], row["process"], row["path"])
```

- **示例输出：**

```text
2
local count_lines COUNT_LINES modules/local/count_lines.nf
nf-core fastqc FASTQC modules/nf-core/fastqc/main.nf
```

## 7) `dnakit.workflows.list_galaxy_tools()`

- **作用：** 列出 Galaxy 服务器上安装的工具。
- **处理方法：** 通过 bioblend 的 `GalaxyInstance.tools.get_tools()` 请求 `/api/tools`，并按 `id` 排序。API key 放在 `x-api-key` 请求头中发送，不会出现在输出里。`url` 不能包含账号信息、查询串或片段。
- **API：** `dnakit.workflows.list_galaxy_tools(url, *, api_key=None, name=None, max_items=1000)`。
- **输入：** Galaxy 服务器地址（`http://` 或 `https://`）；`api_key` 可选；`name` 按工具名过滤；`max_items` 是返回条数上限。
- **输出：** 含 `server`、`count`、`items`（每项为 `id`、`name`、`version`、`section`）和 `provenance` 的字典。
- **示例代码（未本地运行）：**

```python
import os
from dnakit.workflows import list_galaxy_tools

result = list_galaxy_tools(
    "https://usegalaxy.org", api_key=os.environ["GALAXY_API_KEY"], name="Concatenate"
)
print(result["count"], [item["id"] for item in result["items"]])
```

- **示例输出（未本地运行）：**

```text
（未本地运行：需要可访问的 Galaxy 服务器和 API key。输出为匹配的工具数量和 ID 列表。）
```

## 8) `dnakit.workflows.list_galaxy_histories()`

- **作用：** 列出当前用户的未删除历史。
- **处理方法：** 通过 bioblend 的 `GalaxyInstance.histories.get_histories()` 请求 `/api/histories`，并按 `id` 排序。认证方式与 `list_galaxy_tools()` 相同。
- **API：** `dnakit.workflows.list_galaxy_histories(url, *, api_key=None, name=None, max_items=1000)`。
- **输入：** 与 `list_galaxy_tools()` 相同；`name` 按历史名过滤。
- **输出：** 含 `server`、`count`、`items`（每项为 `id`、`name`、`count`、`deleted`、`published`）和 `provenance` 的字典。
- **示例代码（未本地运行）：**

```python
import os
from dnakit.workflows import list_galaxy_histories

result = list_galaxy_histories("https://usegalaxy.org", api_key=os.environ["GALAXY_API_KEY"])
print(result["count"], [item["name"] for item in result["items"]])
```

- **示例输出（未本地运行）：**

```text
（未本地运行：需要可访问的 Galaxy 服务器和 API key。输出为历史数量和名称列表。）
```

## 9) `dnakit.workflows.list_galaxy_workflows()`

- **作用：** 列出 Galaxy 中保存的工作流，即 Galaxy 的流水线定义。
- **处理方法：** 通过 bioblend 的 `GalaxyInstance.workflows.get_workflows()` 请求 `/api/workflows`，并按 `id` 排序。本函数只读取列表，不调用工作流。
- **API：** `dnakit.workflows.list_galaxy_workflows(url, *, api_key=None, name=None, max_items=1000)`。
- **输入：** 与 `list_galaxy_tools()` 相同；`name` 按工作流名过滤。
- **输出：** 含 `server`、`count`、`items`（每项为 `id`、`name`、`published`）和 `provenance` 的字典。
- **示例代码（未本地运行）：**

```python
import os
from dnakit.workflows import list_galaxy_workflows

result = list_galaxy_workflows("https://usegalaxy.org", api_key=os.environ["GALAXY_API_KEY"])
print(result["count"], [item["name"] for item in result["items"]])
```

- **示例输出（未本地运行）：**

```text
（未本地运行：需要可访问的 Galaxy 服务器和 API key。输出为工作流数量和名称列表。）
```

## 验证方式

- 测试文件为 `others/tests/unit/test_workflow_backends.py`，共 25 个用例。在安装了 miniwdl、cwltool、snakemake、nf-core 和 bioblend 的环境中，24 个通过；实时 Galaxy 用例只在设置 `DNAKIT_LIVE_TESTS=1` 时运行，默认跳过。在没有这些包的环境中，依赖上游包的用例会跳过，缺包报错的用例仍然执行。
- WDL 和 CWL 的测试在 `tmp_path` 中写入手写的定义文件，并与手写的任务、输入、步骤和来源对照。CWL 工作流列表与 Snakemake 计划在不同 `PYTHONHASHSEED` 下完全一致。
- Snakemake 的 dry-run 测试用手工计算的作业数（4 个）核对计划，并检查 dry-run 没有生成输出文件；输入缺失和不存在的目标会报错。
- nf-core 测试在临时目录中构造一个流水线检出目录，离线运行。
- Galaxy 测试用假的 `requests.get` 替代服务器，bioblend 的客户端代码照常运行。测试检查 API key 只出现在请求头中，超出 `max_items` 和 HTTP 错误会报错。测试不访问网络。
- 以上结果由手写定义直接核对，没有与上游命令行的输出逐行比较。

## 当前范围

- 已包装：WDL 解析与输入检查（miniwdl）；CWL Workflow 加载与步骤列表（cwltool）；Snakefile 规则列表与 dry-run 作业计划（Snakemake API）；nf-core 检出目录的进程列表（nf_core，离线）；Galaxy 工具、历史和工作流列表（bioblend）。
- 未包装：流水线执行（`miniwdl run`、`cwltool` 执行、Snakemake 执行、Nextflow 运行）。Galaxy 的工作流调用、数据上传和工具运行也未包装，因为它们会改变服务器状态。nf-core 的 `pipelines create`、`lint` 和 `sync` 等命令同样未包装。
- 限制：只支持本地文件，不支持远程 WDL 导入或 URL 形式的 CWL 输入。nf-core 只输出进程名和路径。Snakemake 的 dry-run 需要上游的私有属性 `WorkflowApi._workflow`，因此 `bio-snakemake` 限定了版本范围；dry-run 会在工作目录中写入 `.snakemake/`，但不执行任何规则。
