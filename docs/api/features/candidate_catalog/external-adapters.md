# 外部工具适配

本类收录 `external_adapters` 模块的 **9** 个函数，全部会计算出真实结果。`run` 以参数列表启动已安装的程序（绝不经过 shell，带超时），仅当所请求的程序未安装时才返回 `status: "not_implemented"`。`check` 在 `PATH` 中查找程序；其余七个（`build_command`、`convert_input`、`convert_output`、`make_report`、`parse_log`、`parse_output`、`validate_output`）是纯文本转换。它们本身都不封装任何科学分析工具。示例输出均为真实运行结果；`input_summary` 简写为 `{...}`。

## 1) build_command

- **作用：** 由程序、选项和参数构建外部程序的 argv 列表（绝不生成 shell 字符串）。不执行任何程序。

- **处理方法：** 先是程序，然后是 `options.subcommand`（给出时），接着是选项，最后是参数。选项按映射顺序输出：`true` 输出裸选项，`false`/`null` 省略，列表则每项重复一次该选项，其他值输出 `选项 值`；选项名必须以 `-` 开头，任何词都不能含 NUL 字节或为布尔值。`options.flags` 和 `options.args` 追加在输入自带的之后。由于选项排在参数之前，`{"program": "samtools", "args": ["view", "x.bam"], "flags": {"-b": true}}` 会得到 `samtools -b view x.bam`，这对需要先给子命令的程序是错的。请改用 `options.subcommand` 给出子命令：它紧跟在程序之后、选项之前输出（`samtools view -b x.bam`）。对 argv 列表和命令行字符串形式，子命令插入在第一个词之后。不给 `subcommand` 时顺序不变。

- **API：** `dnakit.external_adapters.build_command(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：程序名、argv 列表、命令行字符串（用 `shlex` 切分，不经过 shell）、path-like，或映射 `{"program": str, "args": [...], "flags": {...}}`。`options`：`subcommand`（字符串如 `"view"`，或列表如 `["index", "build"]`；每项须为非空字符串，否则抛 `ConfigurationError`，`INVALID_OPTION`）、`flags` 和 `args`（形式同上，追加在后）。非法命令会抛 `InputFormatError`（`INVALID_COMMAND`）。

- **输出：** 封装包含 `status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`；`result` 含 `argv`、`program` 和 `command_line`（`shlex.join`，仅供显示）。

- **示例代码：**

```python
import json
from dnakit.external_adapters import build_command

result = build_command(
    {"program": "fastqc", "args": ["reads.fq"], "flags": {"--quiet": True, "--threads": 4, "--nogroup": None, "--adapters": ["a.txt", "b.txt"]}}
)
print(json.dumps(result, indent=2))

# a program with a sub-command: options.subcommand goes right after the program, before the flags
view = build_command(
    {"program": "samtools", "args": ["x.bam"], "flags": {"-b": True}}, options={"subcommand": "view"}
)
print(json.dumps(view["result"], indent=2))
```

- **示例输出：**

```json
{
  "module": "external_adapters",
  "operation": "build_command",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "argv": [
      "fastqc",
      "--quiet",
      "--threads",
      "4",
      "--adapters",
      "a.txt",
      "--adapters",
      "b.txt",
      "reads.fq"
    ],
    "program": "fastqc",
    "command_line": "fastqc --quiet --threads 4 --adapters a.txt --adapters b.txt reads.fq"
  }
}
{
  "argv": [
    "samtools",
    "view",
    "-b",
    "x.bam"
  ],
  "program": "samtools",
  "command_line": "samtools view -b x.bam"
}
```

## 2) check

- **作用：** 检查外部程序是否已安装，可选读取其版本。

- **处理方法：** 用 `shutil.which` 解析每个程序。仅当给出 `version_args` 时才运行 `[program, *version_args]`（不经 shell、无标准输入、带超时），并把其 stdout 与 stderr 合并后的第一行作为 `version`。程序缺失时报告 `available: false`，调用本身仍成功。

- **API：** `dnakit.external_adapters.check(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：程序名/路径，或其列表。`options`：`version_args`（如 `["--version"]` 的列表；默认无，即不执行任何程序）、`timeout`（秒，默认 10；超时时 `version` 为 `null` 并设置 `version_error`）。

- **输出：** 封装包含 `status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`；`result` 含 `all_available` 和 `programs`（每个程序一行：`program`、`available`、`path`、`version`）。示例中的 `path` 和 `version` 取决于运行机器。

- **示例代码：**

```python
import json
from dnakit.external_adapters import check

result = check(["python3", "no_such_tool_xyz"], options={"version_args": ["--version"]})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "external_adapters",
  "operation": "check",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "all_available": false,
    "programs": [
      {
        "program": "python3",
        "available": true,
        "path": "/opt/conda/bin/python3",
        "version": "Python 3.8.12"
      },
      {
        "program": "no_such_tool_xyz",
        "available": false,
        "path": null,
        "version": null
      }
    ]
  }
}
```

## 3) convert_input

- **作用：** 把 DNA 序列记录转换为外部工具所需的文本格式（FASTA、FASTQ 或 TSV）。

- **处理方法：** 序列由 `dnakit.sequence.normalize_dna_sequence` 规范化（转大写、`U` 转 `T`，只允许 IUPAC 字母；其他字符会抛 `InputFormatError`，`INVALID_SEQUENCE`）。没有 id 的记录依次命名为 `seq1`、`seq2`……。FASTQ 对每个碱基使用同一个质量字符（并不存在真实质量值），因此这只是格式转换，不产生质量数据。FASTA 头只保留 id（丢弃描述）。

- **API：** `dnakit.external_adapters.convert_input(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：FASTA 文本（以 `>` 开头）、FASTA 文件的 path-like、单条序列字符串，或序列列表、`(id, sequence)` 对、带 `id`/`sequence` 的映射。`options`：`target_format`（默认 `"fasta"`，或 `"fastq"`、`"tsv"`）、`line_width`（FASTA 换行宽度，默认 80，0 表示不换行）、`quality`（FASTQ 中每个碱基使用的单个字符，默认 `"I"`）。

- **输出：** 封装包含 `status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`；`result` 含 `format`、`text` 和 `record_count`。

- **示例代码：**

```python
import json
from dnakit.external_adapters import convert_input

result = convert_input(">a desc\nacgu\n>b\nGGAA\n", options={"target_format": "fastq"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "external_adapters",
  "operation": "convert_input",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "fastq",
    "text": "@a\nACGT\n+\nIIII\n@b\nGGAA\n+\nIIII\n",
    "record_count": 2
  }
}
```

## 4) convert_output

- **作用：** 把外部工具的输出（TSV、CSV、JSON、JSON Lines、FASTA、纯文本行）转换为 JSON、JSON Lines、TSV 或 CSV 文本。

- **处理方法：** 先按 `parse_output` 的方式解析源文本，再以 `target_format` 写出。表格目标需要字典记录；列为所有键的并集（按首次出现顺序），`null` 写为空单元格。JSON 使用 2 空格缩进。

- **API：** `dnakit.external_adapters.convert_output(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：文本字符串、已存在文件的 `pathlib.Path`，或含 `stdout`（或 `text`/`output`/`path`）的运行结果映射。`options`：`target_format`（默认 `"json"`，或 `"jsonl"`、`"tsv"`、`"csv"`）以及 `parse_output` 的解析选项：`format`（默认 `"auto"`：以 `>` 开头为 FASTA，以 `{` 或 `[` 开头为 JSON，含制表符为 TSV，否则按行；也可指定 `tsv`、`csv`、`json`、`jsonl`、`fasta`、`lines`）、`delimiter`、`header`（默认 true：首行为列名）、`columns`（显式列名，代替表头）、`comment_prefix`（默认 `"#"`；TSV/CSV 中此类行被跳过）、`convert_types`（默认 true：整数和浮点单元格转为数字）。

- **输出：** 封装包含 `status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`；`result` 含 `format`（目标格式）、`text` 和 `record_count`。

- **示例代码：**

```python
import json
from dnakit.external_adapters import convert_output

result = convert_output("#comment\nid\tcount\tname\nr1\t3\talpha\nr2\t4.5\tbeta\n", options={"target_format": "csv"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "external_adapters",
  "operation": "convert_output",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "csv",
    "text": "id,count,name\nr1,3,alpha\nr2,4.5,beta\n",
    "record_count": 2
  }
}
```

## 5) make_report

- **作用：** 把一次外部程序运行（`run` 的结果）汇总为 Markdown 或 JSON 报告。

- **处理方法：** 报告列出命令（`argv` 的 `shlex.join`）、返回码、是否成功（返回码为 0）、耗时、对 stderr 再接 stdout 应用 `parse_log` 关键词规则得到的各级别计数，以及 stdout 和 stderr 的前 `max_lines` 行。不执行任何程序。

- **API：** `dnakit.external_adapters.make_report(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：`run` 的结果映射（需要 `argv` 和 `returncode`；可选 `stdout`、`stderr`、`duration_s`）或完整的 `run` 封装。`options`：`format`（默认 `"markdown"`，或 `"json"`）、`title`（默认 `"External program run"`）、`max_lines`（默认 20）。

- **输出：** 封装包含 `status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`；`result` 含 `format` 和 `report`（Markdown 文本；`format: "json"` 时为含 `command`、`returncode`、`success`、`duration_s`、`log`、`stdout_head`、`stderr_head` 的映射）。

- **示例代码：**

```python
import json
from dnakit.external_adapters import make_report

run_result = {
    "argv": ["/usr/bin/samtools", "flagstat", "x.bam"],
    "returncode": 1,
    "stdout": "10 + 0 in total\n",
    "stderr": "[W::hts_open] warning: low quality\n[E::main] error: truncated file\n",
    "duration_s": 0.5,
}
result = make_report(run_result, options={"title": "samtools flagstat"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "external_adapters",
  "operation": "make_report",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "markdown",
    "report": "# samtools flagstat\n\n- command: `/usr/bin/samtools flagstat x.bam`\n- return code: 1\n- success: no\n- duration: 0.5 s\n- log: 1 error, 1 warning line(s)\n\n## stdout (head)\n\n    10 + 0 in total\n\n## stderr (head)\n\n    [W::hts_open] warning: low quality\n    [E::main] error: truncated file\n"
  }
}
```

## 6) parse_log

- **作用：** 用启发式关键词规则按严重程度对外部程序日志的各行分类。

- **处理方法：** 每行按第一个匹配的整词分类，匹配顺序固定且不区分大小写：error（`fatal`、`critical`、`error`、`exception`、`traceback`、`segmentation fault`）、warning（`warn`、`warning`）、`info`、`debug`；不含关键词的行不分类（只计入 `line_count`）。该规则是启发式的，并非针对某个具体工具日志格式的解析器。超过 200,000 行的日志会抛 `InputFormatError`。

- **API：** `dnakit.external_adapters.parse_log(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：日志文本、文件的 `pathlib.Path`，或运行结果映射（含 `stderr` 时先扫描 stderr 再扫描 stdout，否则扫描 `stdout`）。`options` 不起作用。

- **输出：** 封装包含 `status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`；`result` 含 `line_count`、`levels`（error、warning、info、debug 的计数）、`errors` 和 `warnings`（各最多 100 条样例 `{line, text}`，文本截断为 500 字符）以及 `has_errors`。

- **示例代码：**

```python
import json
from dnakit.external_adapters import parse_log

result = parse_log("INFO start\nWARNING low coverage\nERROR failed to open file\nplain line\n")
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "external_adapters",
  "operation": "parse_log",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "line_count": 4,
    "levels": {
      "error": 1,
      "warning": 1,
      "info": 1,
      "debug": 0
    },
    "errors": [
      {
        "line": 3,
        "text": "ERROR failed to open file"
      }
    ],
    "warnings": [
      {
        "line": 2,
        "text": "WARNING low coverage"
      }
    ],
    "has_errors": true
  }
}
```

## 7) parse_output

- **作用：** 把外部工具的输出解析为结构化记录。

- **处理方法：** TSV/CSV 的每行变成以表头（或 `columns`，无表头时为 `column_1`、`column_2`……）为键的字典；字段数不一致的行会抛 `InputFormatError`（`INVALID_OUTPUT`）。启用 `convert_types` 时单元格先尝试转为整数，再尝试转为浮点数，否则保留文本。FASTA 记录为 `{id, description, sequence}`（序列行已拼接）。JSON 返回列表（单个值包装成列表），JSON Lines 每行一条记录，`lines` 为非空行。超过 50,000,000 字节的输入会被拒绝。

- **API：** `dnakit.external_adapters.parse_output(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：文本字符串、已存在文件的 `pathlib.Path`，或含 `stdout`（或 `text`/`output`/`path`）的运行结果映射。`options`：`format`（默认 `"auto"`：以 `>` 开头为 FASTA，以 `{` 或 `[` 开头为 JSON，含制表符为 TSV，否则按行；也可指定 `tsv`、`csv`、`json`、`jsonl`、`fasta`、`lines`）、`delimiter`、`header`（默认 true：首行为列名）、`columns`（显式列名，代替表头）、`comment_prefix`（默认 `"#"`；TSV/CSV 中此类行被跳过）、`convert_types`（默认 true：整数和浮点单元格转为数字）。

- **输出：** 封装包含 `status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`；`result` 含 `format`（实际使用的格式）、`records` 和 `record_count`。

- **示例代码：**

```python
import json
from dnakit.external_adapters import parse_output

result = parse_output("id\tcount\tname\nr1\t3\talpha\nr2\t4.5\tbeta\n")
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "external_adapters",
  "operation": "parse_output",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "tsv",
    "records": [
      {
        "id": "r1",
        "count": 3,
        "name": "alpha"
      },
      {
        "id": "r2",
        "count": 4.5,
        "name": "beta"
      }
    ],
    "record_count": 2
  }
}
```

## 8) run

- **作用：** 运行外部程序（argv 列表、不经 shell、带超时）并捕获其输出。

- **处理方法：** 命令按 `build_command` 的方式转为 argv 列表，程序用 `shutil.which` 解析，并通过不经 shell 的 `subprocess.run` 启动。程序未安装时返回 `status: "not_implemented"`，`required_dependency` 为 `"the external program '<name>' on PATH"`，`message` 说明在 PATH 中找不到该程序。非零退出码记录在结果中而不抛异常。超时则终止进程并抛 `BackendTimeoutError`（`EXTERNAL_PROGRAM_TIMEOUT`）；无法启动的程序抛 `ConfigurationError`。stdout 和 stderr 各按 `max_output_bytes` 截断，并以 `truncated` 标记且给出警告。

- **API：** `dnakit.external_adapters.run(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：程序名、argv 列表、命令行字符串、path-like 或映射（见 `build_command`）。`options`：`timeout`（秒，默认 60，须在 (0, 3600] 内）、`cwd`（已存在的目录）、`env`（额外环境变量映射）、`stdin`（发送给程序的文本）、`max_output_bytes`（默认 1,000,000）、`flags` 和 `args`（追加）以及 `subcommand`（紧跟在程序之后，与 `build_command` 完全一致，因此 `run` 执行的就是 `build_command` 给出的 argv）。

- **输出：** 封装包含 `status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`（`provenance.implementation` 为 `subprocess`）；`result` 含 `argv`（带解析后的程序路径）、`returncode`、`success`、`stdout`、`stderr`、`duration_s` 和 `truncated`。程序缺失时 `status` 为 `not_implemented`，没有 `result`。示例中 `echo` 的解析路径取决于运行机器，实测的 `duration_s` 被替换为占位符 `"<seconds>"`。

- **示例代码：**

```python
import json
from dnakit.external_adapters import run

result = run(["echo", "hello"], options={"timeout": 10})
result["result"]["duration_s"] = "<seconds>"
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "external_adapters",
  "operation": "run",
  "provenance": {
    "implementation": "subprocess",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "argv": [
      "/usr/bin/echo",
      "hello"
    ],
    "returncode": 0,
    "success": true,
    "stdout": "hello\n",
    "stderr": "",
    "duration_s": "<seconds>",
    "truncated": false
  }
}
```

## 9) validate_output

- **作用：** 检查外部程序是否产生了可用的输出。

- **处理方法：** 每项检查以 `{name, passed, detail}` 记入 `checks`：`not_empty`（总是执行）、`returncode` 是否等于 `expected_returncode`（仅对运行结果）、`parses`（文本能否按给定 `format` 解析）、`min_records`、`required_columns`（针对表格或 JSON 字典记录），以及每个子串一条 `must_contain`。输出未被解析（未给 `format`）时 `min_records` 和 `required_columns` 判为失败。不执行任何程序。

- **API：** `dnakit.external_adapters.validate_output(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：文本字符串、已存在文件的 `pathlib.Path`，或含 `stdout`（或 `text`/`output`/`path`）的运行结果映射。`options`：`format`（启用 `parses` 检查）、`min_records`、`required_columns`、`must_contain`（子串列表）、`expected_returncode`（默认 0）以及 `parse_output` 的解析选项。

- **输出：** 封装包含 `status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`；`result` 含 `valid`（所有检查均通过）和 `checks`。

- **示例代码：**

```python
import json
from dnakit.external_adapters import validate_output

result = validate_output(
    "id\tcount\nr1\t3\nr2\t4\n",
    options={"format": "tsv", "min_records": 3, "required_columns": ["id", "count", "name"], "must_contain": ["r1"]},
)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "external_adapters",
  "operation": "validate_output",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "valid": false,
    "checks": [
      {
        "name": "not_empty",
        "passed": true,
        "detail": "19 characters"
      },
      {
        "name": "parses",
        "passed": true,
        "detail": "2 records"
      },
      {
        "name": "min_records",
        "passed": false,
        "detail": "2 records, need 3"
      },
      {
        "name": "required_columns",
        "passed": false,
        "detail": "missing: ['name']"
      },
      {
        "name": "must_contain",
        "passed": true,
        "detail": "'r1'"
      }
    ]
  }
}
```
