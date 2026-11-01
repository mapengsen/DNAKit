# DendroPy 与 khmer：树格式、PAUP* 与 k-mer 表

这一组函数是两个可选包的适配器。[DendroPy](https://github.com/jeetsukumaran/DendroPy) 是用于系统发育计算的 Python 库；`dnakit.phylogeny` 用它在 FASTA、Newick、NEXUS、NeXML 和 PHYLIP 之间转换序列与树文件，运行它的 `dendropy-format`、`sumtrees`、`sumlabels` 命令行应用，从 NCBI Entrez 获取记录，并驱动外部程序 MUSCLE（比对）和 PAUP*（建树与模型估计、剪枝、简约性得分）。[khmer](https://github.com/dib-lab/khmer) 是近似 k-mer 计数库；`dnakit.microbiome` 用它估算 khmer 表的大小，并从读段构建并保存计数表（countgraph）和存在性表（nodegraph）。函数通过 `from dnakit.phylogeny import ...` 和 `from dnakit.microbiome import ...` 导入。

按包元数据，DendroPy 和 khmer 使用 BSD 许可。MUSCLE 和 PAUP* 是独立的程序，DNAKit 不分发它们，这里也没有核实它们的许可。

依赖都是可选的，只有调用对应函数时才会导入。缺少包时会报错并给出安装命令（`dnakit[bio-dendropy]` 或 `dnakit[bio-khmer]`）。输出文件总是新文件：输出路径已存在时会被拒绝，不会被覆盖。

生成本页的环境有 DendroPy 5.1.0（含 `dendropy-format`、`sumtrees`、`sumlabels` 脚本）和 khmer 2.1.1，但没有 MUSCLE 和 PAUP*，也没有网络。因此有十二个示例运行真实计算（第 2、4 至 8、14 至 19 节）。第 1、3、9 至 13 节展示程序或输入缺失时函数真实抛出的错误，并且每一节后面都有一段未执行的“安装依赖后的用法”代码块。

## 示例数据

示例使用五个很小的已声明文件：一个四条序列的 DNA 比对、一个含模糊编码的两条序列文件、一个含三棵树的 Newick 文件、一个 Newick 目标树，以及一个含五条短读段的 FASTA 文件（两条相同、一条过短、一条含 `N`）。所有输出都写入临时目录。

**示例文件：** `aln.fa`

```text
>t1
ACGTACGTAC
>t2
ACGTACGAAC
>t3
ACCTACGTTC
>t4
TCCTACGTTC
```

**示例文件：** `iupac.fa`

```text
>a
ACGTRYAC
>b
ACGTACNC
```

**示例文件：** `trees.nwk`

```text
((t1:0.1,t2:0.1)0.9:0.2,(t3:0.1,t4:0.1)0.8:0.2);
((t1:0.1,t2:0.1)0.7:0.2,(t3:0.1,t4:0.1)0.6:0.2);
((t1:0.1,t3:0.1)0.5:0.2,(t2:0.1,t4:0.1)0.4:0.2);
```

**示例文件：** `target.nwk`

```text
((t1,t2),(t3,t4));
```

**示例文件：** `reads.fa`

```text
>r1
ACGTACGTACGTAGCT
>r2
ACGTACGTACGTAGCT
>r3
TTGCAGGTCCAT
>r4_short
ACG
>r5_n
ACGTNCGTACGT
```

## 1) `dnakit.phylogeny.align_with_muscle()`

- **作用：** 用外部程序 MUSCLE 对 DNA 序列做比对，并把比对结果保存为 FASTA。
- **处理方法：** 把 FASTA 文件读成 dendropy 的 `DnaCharacterMatrix`，通过 dendropy 的 `muscle_align` 启动 MUSCLE（序列从标准输入送入，比对结果从 MUSCLE 的标准输出读回），再把比对矩阵以 FASTA 格式写入新文件。可执行文件依次从 `muscle_path`、PATH、当前解释器所在目录查找。退出码非零时抛出带有 MUSCLE 消息的 `BackendExecutionError`。
- **API：** `align_with_muscle(fasta: 'str', output: 'str', *, muscle_args: 'Sequence[str] | None' = None, muscle_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 dendropy 5.1（额外依赖 `dnakit[bio-dendropy]`）以及 PATH 上的外部程序 `muscle`（或通过 `muscle_path` 指定）。生成本页的环境没有安装 MUSCLE，因此下面的示例展示找不到程序时真实抛出的错误；其后给出安装程序后的调用（未执行）。
- **输入：** `fasta` 是已存在的 DNA FASTA 文件；`output` 是已存在文件夹中尚不存在的路径；`muscle_args` 是可选的非空字符串列表，会追加到 MUSCLE 命令行（默认为空，即 MUSCLE 3.8 风格）；`muscle_path` 是可选的可执行文件路径。
- **输出：** 字典，包含 `status`、`backend`、`version`、`command`、`sequence_count`、`alignment_length` 和 `outputs`（新建的 FASTA 文件）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.phylogeny import align_with_muscle

with tempfile.TemporaryDirectory() as folder:
    try:
        align_with_muscle("aln.fa", str(Path(folder) / "aligned.fa"), muscle_path="no_such_dir/muscle", verbose=False)
    except BackendUnavailableError as error:
        print(type(error).__name__, str(error).split(":")[0])
```

- **示例输出：**

```text
BackendUnavailableError [BACKEND_UNAVAILABLE] muscle executable does not exist
```

安装依赖后的用法：

```python
from dnakit.phylogeny import align_with_muscle

result = align_with_muscle("aln.fa", "aligned.fa", verbose=False)
print(result["sequence_count"], result["alignment_length"])
```

## 2) `dnakit.phylogeny.convert_format()`

- **作用：** 把 DNA 序列文件或树文件转换成另一种格式。
- **处理方法：** 沿用 `dendropy-format` 应用的步骤，但在当前进程内完成：用 `dendropy.DataSet.get` 读取文件，按需重编码模糊与多态 DNA 编码，再用 `DataSet.write` 写出。它有意不调用该应用的 `convert()` 函数，因为后者会关闭 `sys.stdout`。`recode_uncertain="missing"` 时这些编码变成 `?`，`"gap"` 时变成 `-`。
- **API：** `convert_format(source: 'str', output: 'str', *, input_format: 'str', output_format: 'str | None' = None, data_type: 'str | None' = None, unquoted_underscores: 'bool' = False, recode_uncertain: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 dendropy 5.1（额外依赖 `dnakit[bio-dendropy]`）。安装该额外依赖后即可使用，不需要外部程序。
- **输入：** `source` 是已存在的文件；`output` 是新路径；`input_format` 取 `fasta`、`newick`、`nexml`、`nexus` 或 `phylip-*` 名称（`phylip-strict`、`phylip-relaxed-singlespace`、`phylip-relaxed-multispace` 以及三个 `-interleaved` 变体）；`output_format` 取 `fasta`、`newick`、`nexml`、`nexus`、`phylip`、`phylip-strict`，默认与读取格式相同；读取 FASTA 和 PHYLIP 时必须给 `data_type`（`dna`、`rna` 或 `standard`）；`unquoted_underscores` 作用于 NEXUS 和 Newick 输出；`recode_uncertain` 取 `missing` 或 `gap`。
- **输出：** 字典，包含 `status`、`backend`、`version`、`operation`、`input_format`、`output_format`、`data_type`、`tree_lists` 与 `character_matrices`（读到的数量）以及 `outputs`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.phylogeny import convert_format

with tempfile.TemporaryDirectory() as folder:
    output = Path(folder) / "iupac.nex"
    result = convert_format(
        "iupac.fa", str(output), input_format="fasta", output_format="nexus", data_type="dna",
        recode_uncertain="missing", verbose=False,
    )
    print(result["status"], result["input_format"], result["output_format"], result["character_matrices"])
    for line in output.read_text().splitlines():
        if len(line.split()) == 2 and line.split()[0] in ("a", "b"):
            print(" ".join(line.split()))
```

- **示例输出：**

```text
completed fasta nexus 1
a ACGT??AC
b ACGTAC?C
```

## 3) `dnakit.phylogeny.fetch_entrez_records()`

- **作用：** 按数字形式的 Entrez UID 从 NCBI 数据库获取序列记录。
- **处理方法：** 先检查参数，再调用 `dendropy.interop.entrez.efetch`，它通过 HTTP 向 NCBI 发送 E-utilities efetch 请求，读取完整响应后关闭连接。多个标识符用逗号连接。不做缓存，也不写磁盘。`record_count` 统计以 `>` 开头的标题行，仅在 `rettype` 为 `fasta` 且 `retmode` 为 `text` 时给出。
- **API：** `fetch_entrez_records(database: 'str', ids: 'Sequence[str | int]', *, rettype: 'str' = 'fasta', retmode: 'str' = 'text', email: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要能访问 NCBI 的网络（不传 key，`email` 可选）以及可选包 dendropy 5.1（额外依赖 `dnakit[bio-dendropy]`）。生成本页的环境没有网络，因此示例展示发出请求之前就会执行的真实输入检查；其后给出联网时的调用（未执行）。
- **输入：** `database` 是不以 `-` 开头的非空名称（例如 `nuccore`）；`ids` 是 1 到 200 个数字 UID 的列表（整数或数字字符串；`NM_000546` 这样的登录号会被拒绝）；`rettype` 是非空名称（默认 `fasta`）；`retmode` 取 `text` 或 `xml`；`email` 可选，必须含有 `@`。
- **输出：** 字典，包含 `status`、`backend`、`service`、`database`、`rettype`、`retmode`、`identifier_count`、`record_count`、`access_date`（UTC 日期）、`url`（最终请求地址）和 `records`（响应文本）。
- **示例代码：**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.phylogeny import fetch_entrez_records

try:
    fetch_entrez_records("nuccore", ["NM_000546"], verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(". ")[0])
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] Entrez identifiers are numeric UIDs; got 'NM_000546'
```

安装依赖后的用法：

```python
from dnakit.phylogeny import fetch_entrez_records

# 1234567 stands for any numeric nuccore UID, for example one found with an ESearch query.
result = fetch_entrez_records("nuccore", [1234567], email="you@example.org", verbose=False)
print(result["record_count"], result["access_date"])
print(result["records"][:80])
```

## 4) `dnakit.phylogeny.format_to_fasta()`

- **作用：** 把序列文件（或 NEXUS 文件中的序列）转换成 FASTA。
- **处理方法：** 以 `output_format="fasta"` 调用 `convert_format`，相当于 `dendropy-format --to fasta`。
- **API：** `format_to_fasta(source: 'str', output: 'str', *, input_format: 'str', data_type: 'str | None' = None, recode_uncertain: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 dendropy 5.1（额外依赖 `dnakit[bio-dendropy]`）。安装该额外依赖后即可使用，不需要外部程序。
- **输入：** `source`、`output`、`input_format`、`data_type`、`recode_uncertain` 与 `convert_format` 相同（FASTA 和 PHYLIP 输入需要 `data_type`）。
- **输出：** 与 `convert_format` 相同的字典，其中 `output_format` 为 `fasta`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.phylogeny import format_to_fasta

with tempfile.TemporaryDirectory() as folder:
    output = Path(folder) / "gapped.fa"
    result = format_to_fasta(
        "iupac.fa", str(output), input_format="fasta", data_type="dna", recode_uncertain="gap", verbose=False
    )
    print(result["output_format"])
    print(output.read_text().strip())
```

- **示例输出：**

```text
fasta
>a
ACGT--AC
>b
ACGTAC-C
```

## 5) `dnakit.phylogeny.format_to_newick()`

- **作用：** 把树文件（或 NEXUS 文件中的树）转换成 Newick。
- **处理方法：** 以 `output_format="newick"` 调用 `convert_format`，相当于 `dendropy-format --to newick`。每棵树单独占一行。
- **API：** `format_to_newick(source: 'str', output: 'str', *, input_format: 'str', data_type: 'str | None' = None, unquoted_underscores: 'bool' = False, recode_uncertain: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 dendropy 5.1（额外依赖 `dnakit[bio-dendropy]`）。安装该额外依赖后即可使用，不需要外部程序。
- **输入：** `source`、`output`、`input_format`、`data_type`、`unquoted_underscores`、`recode_uncertain` 与 `convert_format` 相同。
- **输出：** 与 `convert_format` 相同的字典，其中 `output_format` 为 `newick`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.phylogeny import format_to_newick

with tempfile.TemporaryDirectory() as folder:
    output = Path(folder) / "trees_out.nwk"
    result = format_to_newick("trees.nwk", str(output), input_format="newick", verbose=False)
    print(result["output_format"], result["tree_lists"], result["character_matrices"])
    print(output.read_text().splitlines()[0])
```

- **示例输出：**

```text
newick 1 0
((t1:0.1,t2:0.1)0.9:0.2,(t3:0.1,t4:0.1)0.8:0.2);
```

## 6) `dnakit.phylogeny.format_to_nexml()`

- **作用：** 把序列文件或树文件转换成 NeXML。
- **处理方法：** 以 `output_format="nexml"` 调用 `convert_format`，相当于 `dendropy-format --to nexml`。没有 `unquoted_underscores` 选项，因为它只适用于 NEXUS 和 Newick。
- **API：** `format_to_nexml(source: 'str', output: 'str', *, input_format: 'str', data_type: 'str | None' = None, recode_uncertain: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 dendropy 5.1（额外依赖 `dnakit[bio-dendropy]`）。安装该额外依赖后即可使用，不需要外部程序。
- **输入：** `source`、`output`、`input_format`、`data_type`、`recode_uncertain` 与 `convert_format` 相同。
- **输出：** 与 `convert_format` 相同的字典，其中 `output_format` 为 `nexml`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.phylogeny import format_to_nexml

with tempfile.TemporaryDirectory() as folder:
    output = Path(folder) / "aln.xml"
    result = format_to_nexml("aln.fa", str(output), input_format="fasta", data_type="dna", verbose=False)
    text = output.read_text()
    print(result["output_format"], text.startswith("<?xml"), text.count("<otu "), text.count("<row"))
```

- **示例输出：**

```text
nexml True 4 4
```

## 7) `dnakit.phylogeny.format_to_nexus()`

- **作用：** 把序列文件或树文件转换成 NEXUS。
- **处理方法：** 以 `output_format="nexus"` 调用 `convert_format`，相当于 `dendropy-format --to nexus`。输出包含 TAXA 块以及 TREES 或 CHARACTERS 块。
- **API：** `format_to_nexus(source: 'str', output: 'str', *, input_format: 'str', data_type: 'str | None' = None, unquoted_underscores: 'bool' = False, recode_uncertain: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 dendropy 5.1（额外依赖 `dnakit[bio-dendropy]`）。安装该额外依赖后即可使用，不需要外部程序。
- **输入：** `source`、`output`、`input_format`、`data_type`、`unquoted_underscores`、`recode_uncertain` 与 `convert_format` 相同。
- **输出：** 与 `convert_format` 相同的字典，其中 `output_format` 为 `nexus`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.phylogeny import format_to_nexus

with tempfile.TemporaryDirectory() as folder:
    output = Path(folder) / "trees.nex"
    result = format_to_nexus("trees.nwk", str(output), input_format="newick", verbose=False)
    print(result["output_format"], result["tree_lists"])
    for line in output.read_text().splitlines():
        if line.strip().startswith("TREE"):
            print(line.strip())
```

- **示例输出：**

```text
nexus 1
TREE 1 = ((t1:0.1,t2:0.1)0.9:0.2,(t3:0.1,t4:0.1)0.8:0.2);
TREE 2 = ((t1:0.1,t2:0.1)0.7:0.2,(t3:0.1,t4:0.1)0.6:0.2);
TREE 3 = ((t1:0.1,t3:0.1)0.5:0.2,(t2:0.1,t4:0.1)0.4:0.2);
```

## 8) `dnakit.phylogeny.format_to_phylip()`

- **作用：** 把序列文件转换成宽松（relaxed）PHYLIP。
- **处理方法：** 以 `output_format="phylip"`（dendropy 的 relaxed PHYLIP 写出器）调用 `convert_format`，相当于 `dendropy-format --to phylip`。严格格式可在 `convert_format` 中用 `output_format="phylip-strict"` 得到。
- **API：** `format_to_phylip(source: 'str', output: 'str', *, input_format: 'str', data_type: 'str | None' = None, recode_uncertain: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 dendropy 5.1（额外依赖 `dnakit[bio-dendropy]`）。安装该额外依赖后即可使用，不需要外部程序。
- **输入：** `source`、`output`、`input_format`、`data_type`、`recode_uncertain` 与 `convert_format` 相同。
- **输出：** 与 `convert_format` 相同的字典，其中 `output_format` 为 `phylip`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.phylogeny import format_to_phylip

with tempfile.TemporaryDirectory() as folder:
    output = Path(folder) / "aln.phy"
    result = format_to_phylip("aln.fa", str(output), input_format="fasta", data_type="dna", verbose=False)
    print(result["output_format"])
    print(output.read_text().strip())
```

- **示例输出：**

```text
phylip
4 10
t1  ACGTACGTAC
t2  ACGTACGAAC
t3  ACCTACGTTC
t4  TCCTACGTTC
```

## 9) `dnakit.phylogeny.paup_estimate_model()`

- **作用：** 用 PAUP* 估计 DNA 比对的替换模型参数。
- **处理方法：** 读取比对（FASTA），若给出 `tree_file` 则读取其中第一棵树（首个非空行以 `#NEXUS` 开头时按 NEXUS 解析，否则按 Newick 解析）。调用 dendropy 的 `paup.estimate_model`，在该树上拟合模型；未给树时先用 `criterion` 估计一棵树，然后解析 PAUP* 的报告。PAUP* 所用的树以 Newick 写出。
- **API：** `paup_estimate_model(alignment: 'str', output_tree: 'str', *, tree_file: 'str | None' = None, criterion: 'str' = 'likelihood', substitution_types: 'int' = 6, unequal_base_freqs: 'bool' = True, gamma_rates: 'bool' = True, prop_invar: 'bool' = True, user_branch_lengths: 'bool' = True, paup_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 dendropy 5.1（额外依赖 `dnakit[bio-dendropy]`）以及 PATH 上的外部程序 `paup`（PAUP* 4.0），或通过 `paup_path` 指定。生成本页的环境没有安装 PAUP*，因此示例展示找不到程序时真实抛出的错误；其后给出安装 PAUP* 后的调用（未执行）。可执行文件路径只能包含字母、数字和 `._/+:@%-`，因为 dendropy 的部分函数通过 shell 启动 PAUP*。
- **输入：** `alignment` 是已存在的 DNA FASTA 文件；`output_tree` 是新路径；`tree_file` 是可选的 Newick 或 NEXUS 文件；`criterion` 取 `likelihood`、`parsimony`、`nj` 或 `upgma`；`substitution_types` 取 1、2 或 6；`unequal_base_freqs`、`gamma_rates`、`prop_invar`、`user_branch_lengths` 是布尔值；`paup_path` 可选。
- **输出：** 字典，包含 `status`、`backend`、`operation`、`version`、`outputs`（Newick 文件）、`substitution_types` 和 `estimates`（`log_likelihood`、`kappa`、`rAC` 至 `rGT`、`pA` 至 `pT`、`prop_invar`、`alpha`；PAUP* 未报告的为 `None`，其余为浮点数）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.phylogeny import paup_estimate_model

with tempfile.TemporaryDirectory() as folder:
    try:
        paup_estimate_model("aln.fa", str(Path(folder) / "model_tree.nwk"), paup_path="no_such_dir/paup", verbose=False)
    except BackendUnavailableError as error:
        print(type(error).__name__, str(error).split(":")[0])
```

- **示例输出：**

```text
BackendUnavailableError [BACKEND_UNAVAILABLE] paup executable does not exist
```

安装依赖后的用法：

```python
from dnakit.phylogeny import paup_estimate_model

result = paup_estimate_model("aln.fa", "model_tree.nwk", tree_file="start_tree.nwk", substitution_types=2, verbose=False)
print(result["estimates"]["log_likelihood"], result["estimates"]["kappa"])
```

## 10) `dnakit.phylogeny.paup_estimate_tree()`

- **作用：** 用 PAUP* 从 DNA 比对估计一棵树，并保存为 Newick。
- **处理方法：** 读取比对（FASTA），用所给准则和模型选项调用 dendropy 的 `paup.estimate_tree`，再把树以 Newick 写出。模型选项（`substitution_types`、`unequal_base_freqs`、`gamma_rates`、`prop_invar`）只作用于似然估计。程序在读取比对之前就会被查找，缺少 PAUP* 时抛出 `BackendUnavailableError`。
- **API：** `paup_estimate_tree(alignment: 'str', output_tree: 'str', *, criterion: 'str' = 'likelihood', substitution_types: 'int' = 6, unequal_base_freqs: 'bool' = True, gamma_rates: 'bool' = True, prop_invar: 'bool' = True, paup_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 dendropy 5.1（额外依赖 `dnakit[bio-dendropy]`）以及 PATH 上的外部程序 `paup`（PAUP* 4.0），或通过 `paup_path` 指定。生成本页的环境没有安装 PAUP*，因此示例展示找不到程序时真实抛出的错误；其后给出安装 PAUP* 后的调用（未执行）。可执行文件路径只能包含字母、数字和 `._/+:@%-`，因为 dendropy 的部分函数通过 shell 启动 PAUP*。
- **输入：** `alignment` 是已存在的 DNA FASTA 文件；`output_tree` 是新路径；`criterion` 取 `likelihood`（默认）、`parsimony`、`nj` 或 `upgma`；`substitution_types` 取 1、2 或 6（PAUP* 的 `lset nst`）；三个模型开关是布尔值；`paup_path` 可选。
- **输出：** 字典，包含 `status`、`backend`、`operation`、`version`、`outputs`（Newick 文件）、`criterion`、`substitution_types` 和 `taxon_count`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.phylogeny import paup_estimate_tree

with tempfile.TemporaryDirectory() as folder:
    try:
        paup_estimate_tree("aln.fa", str(Path(folder) / "ml_tree.nwk"), paup_path="no_such_dir/paup", verbose=False)
    except BackendUnavailableError as error:
        print(type(error).__name__, str(error).split(":")[0])
```

- **示例输出：**

```text
BackendUnavailableError [BACKEND_UNAVAILABLE] paup executable does not exist
```

安装依赖后的用法：

```python
from dnakit.phylogeny import paup_estimate_tree

result = paup_estimate_tree("aln.fa", "nj_tree.nwk", criterion="nj", verbose=False)
print(result["criterion"], result["taxon_count"], result["outputs"])
```

## 11) `dnakit.phylogeny.paup_estimate_ultrametric_tree()`

- **作用：** 用 PAUP* 估计超度量（受分子钟约束）树。
- **处理方法：** 读取比对（FASTA），若给出 `topology_tree` 则读取其第一棵树。调用 dendropy 的 `paup.estimate_ultrametric_tree`：没有拓扑时先建邻接法树，再估计受分子钟约束的枝长；给出拓扑时保留该拓扑，只估计枝长。树以带枝长的 Newick 写出。
- **API：** `paup_estimate_ultrametric_tree(alignment: 'str', output_tree: 'str', *, topology_tree: 'str | None' = None, paup_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 dendropy 5.1（额外依赖 `dnakit[bio-dendropy]`）以及 PATH 上的外部程序 `paup`（PAUP* 4.0），或通过 `paup_path` 指定。生成本页的环境没有安装 PAUP*，因此示例展示找不到程序时真实抛出的错误；其后给出安装 PAUP* 后的调用（未执行）。可执行文件路径只能包含字母、数字和 `._/+:@%-`，因为 dendropy 的部分函数通过 shell 启动 PAUP*。
- **输入：** `alignment` 是已存在的 DNA FASTA 文件；`output_tree` 是新路径；`topology_tree` 是可选的 Newick 或 NEXUS 文件，取其第一棵树；`paup_path` 可选。
- **输出：** 字典，包含 `status`、`backend`、`operation`、`version`、`outputs`（Newick 文件）、`taxon_count` 和 `topology_given`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.phylogeny import paup_estimate_ultrametric_tree

with tempfile.TemporaryDirectory() as folder:
    try:
        paup_estimate_ultrametric_tree(
            "aln.fa", str(Path(folder) / "clock_tree.nwk"), paup_path="no_such_dir/paup", verbose=False
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, str(error).split(":")[0])
```

- **示例输出：**

```text
BackendUnavailableError [BACKEND_UNAVAILABLE] paup executable does not exist
```

安装依赖后的用法：

```python
from dnakit.phylogeny import paup_estimate_ultrametric_tree

result = paup_estimate_ultrametric_tree("aln.fa", "clock_tree.nwk", topology_tree="target.nwk", verbose=False)
print(result["topology_given"], result["taxon_count"])
```

## 12) `dnakit.phylogeny.paup_prune_taxa()`

- **作用：** 用 PAUP* 从文件的每棵树中删除指定类群，并把剪枝后的树保存为 NEXUS。
- **处理方法：** 把文件中的所有树读入同一个类群命名空间，并逐个查找 `taxa` 中的名称（名称不存在时在查找 PAUP* 之前就抛出 `ConfigurationError`）。然后调用 dendropy 的 `paup.prune_taxa_from_trees`，把剪枝后的树连同类群块以 NEXUS 写出。
- **API：** `paup_prune_taxa(trees_file: 'str', taxa: 'Sequence[str]', output_trees: 'str', *, paup_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 dendropy 5.1（额外依赖 `dnakit[bio-dendropy]`）以及 PATH 上的外部程序 `paup`（PAUP* 4.0），或通过 `paup_path` 指定。生成本页的环境没有安装 PAUP*，因此示例展示找不到程序时真实抛出的错误；其后给出安装 PAUP* 后的调用（未执行）。可执行文件路径只能包含字母、数字和 `._/+:@%-`，因为 dendropy 的部分函数通过 shell 启动 PAUP*。未知类群名则无需 PAUP* 就会被拒绝，示例中也展示了这一点。
- **输入：** `trees_file` 是已存在的 Newick 或 NEXUS 文件；`taxa` 是非空的单行标签列表，且都必须出现在树中；`output_trees` 是新路径；`paup_path` 可选。
- **输出：** 字典，包含 `status`、`backend`、`operation`、`version`、`outputs`（NEXUS 文件）、`removed`（标签）和 `tree_count`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.phylogeny import paup_prune_taxa

with tempfile.TemporaryDirectory() as folder:
    output = str(Path(folder) / "pruned.nex")
    try:
        paup_prune_taxa("trees.nwk", ["t4"], output, paup_path="no_such_dir/paup", verbose=False)
    except BackendUnavailableError as error:
        print(type(error).__name__, str(error).split(":")[0])
    try:
        paup_prune_taxa("trees.nwk", ["no_such_taxon"], output, paup_path="no_such_dir/paup", verbose=False)
    except ConfigurationError as error:
        print(type(error).__name__, str(error).split(":")[0])
```

- **示例输出：**

```text
BackendUnavailableError [BACKEND_UNAVAILABLE] paup executable does not exist
ConfigurationError [CONFIGURATION_ERROR] Taxon not found in the trees
```

安装依赖后的用法：

```python
from dnakit.phylogeny import paup_prune_taxa

result = paup_prune_taxa("trees.nwk", ["t4"], "pruned.nex", verbose=False)
print(result["removed"], result["tree_count"])
```

## 13) `dnakit.phylogeny.paup_pscore_trees()`

- **作用：** 计算文件中每棵树相对于 DNA 比对的简约性得分（PAUP* 的 `pscore`）。
- **处理方法：** 把树和比对读入同一个类群命名空间，调用 dendropy 的 `paup.pscore_trees`，得分转为整数。`pset_options` 的每一项会被写进 `pset` 命令，因此必须匹配简单词元模式 `[A-Za-z0-9_.=:+-]+`（例如 `gapmode=missing`），其余内容在查找 PAUP* 之前就被拒绝。dendropy 自带的 `pscore_option_list` 没有开放，因为它把选项与得分文件路径直接拼接、中间没有分隔符。
- **API：** `paup_pscore_trees(trees_file: 'str', alignment: 'str', *, pset_options: 'Sequence[str] | None' = None, paup_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 dendropy 5.1（额外依赖 `dnakit[bio-dendropy]`）以及 PATH 上的外部程序 `paup`（PAUP* 4.0），或通过 `paup_path` 指定。生成本页的环境没有安装 PAUP*，因此示例展示找不到程序时真实抛出的错误；其后给出安装 PAUP* 后的调用（未执行）。可执行文件路径只能包含字母、数字和 `._/+:@%-`，因为 dendropy 的部分函数通过 shell 启动 PAUP*。非法的 `pset_options` 无需 PAUP* 就会被拒绝，示例中也展示了这一点。
- **输入：** `trees_file` 是已存在的 Newick 或 NEXUS 文件；`alignment` 是已存在的 DNA FASTA 文件；`pset_options` 是可选的简单词元列表；`paup_path` 可选。
- **输出：** 字典，包含 `status`、`backend`、`operation`、`version`、`outputs`（始终为空）、`scores`（每棵树一个整数）和 `tree_count`。
- **示例代码：**

```python
from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.phylogeny import paup_pscore_trees

try:
    paup_pscore_trees("trees.nwk", "aln.fa", pset_options=["gapmode=missing"], paup_path="no_such_dir/paup", verbose=False)
except BackendUnavailableError as error:
    print(type(error).__name__, str(error).split(":")[0])
try:
    paup_pscore_trees("trees.nwk", "aln.fa", pset_options=["gapmode=missing; quit"], paup_path="no_such_dir/paup", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(";")[0])
```

- **示例输出：**

```text
BackendUnavailableError [BACKEND_UNAVAILABLE] paup executable does not exist
ConfigurationError [CONFIGURATION_ERROR] pset_options must be plain PAUP* tokens such as gapmode=missing
```

安装依赖后的用法：

```python
from dnakit.phylogeny import paup_pscore_trees

result = paup_pscore_trees("trees.nwk", "aln.fa", pset_options=["gapmode=missing"], verbose=False)
print(result["tree_count"], result["scores"])
```

## 14) `dnakit.phylogeny.run_dendropy_format()`

- **作用：** 对文件运行 `dendropy-format` 命令行应用，并把它的标准输出保存为新文件。
- **处理方法：** 构造固定的参数列表（`--from`、`--to`、`--data-type`、`--unquoted-underscores`、`--recode-uncertain`，最后是源文件路径），不经过 shell，通过 DNAKit 的受限命令运行器运行 `dendropy-format`。退出码为 0 时把捕获的输出写入新文件；否则不写文件，`status` 为 `failed`。`recode_uncertain="missing"` 会被 `ConfigurationError` 拒绝，因为 `dendropy-format` 5.1 使用它会崩溃；请改用 `"gap"`，或使用会写出 `?` 的 `convert_format`。不启动子进程的同类转换是 `convert_format`。
- **API：** `run_dendropy_format(source: 'str', output: 'str', *, input_format: 'str', output_format: 'str | None' = None, data_type: 'str | None' = None, unquoted_underscores: 'bool' = False, recode_uncertain: 'str | None' = None, dendropy_format_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 dendropy 5.1（额外依赖 `dnakit[bio-dendropy]`），它会安装 `dendropy-format` 命令行脚本；脚本从 PATH、当前解释器所在目录或 `dendropy_format_path` 查找。
- **输入：** `source` 是已存在的文件；`output` 是新路径；`input_format`、`output_format`、`data_type`、`unquoted_underscores` 与 `convert_format` 相同（读取 FASTA 和 PHYLIP 时必须给 `data_type`）；`recode_uncertain` 取 `gap`；`dendropy_format_path` 可选。
- **输出：** 字典，包含 `status`（`completed` 或 `failed`）、`backend`、`version`、`command`、`return_code`、`outputs`（新文件，失败时为空列表）和 `log_tail`（成功时为空，失败时为程序输出的末尾）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.phylogeny import run_dendropy_format

with tempfile.TemporaryDirectory() as folder:
    output = Path(folder) / "trees.nex"
    result = run_dendropy_format("trees.nwk", str(output), input_format="newick", output_format="nexus", verbose=False)
    print(result["status"], result["backend"], result["return_code"], len(result["outputs"]))
    for line in output.read_text().splitlines():
        if line.strip().startswith("TREE"):
            print(line.strip())
```

- **示例输出：**

```text
completed dendropy-format 0 1
TREE 1 = ((t1:0.1,t2:0.1)0.9:0.2,(t3:0.1,t4:0.1)0.8:0.2);
TREE 2 = ((t1:0.1,t2:0.1)0.7:0.2,(t3:0.1,t4:0.1)0.6:0.2);
TREE 3 = ((t1:0.1,t3:0.1)0.5:0.2,(t2:0.1,t4:0.1)0.4:0.2);
```

## 15) `dnakit.phylogeny.summarize_tree_samples()`

- **作用：** 把一组树样本（例如 MCMC 输出）汇总成一棵带支持率的总结树。
- **处理方法：** 为 `sumtrees` 应用构造固定的参数列表（`--quiet`、`--burnin`、`--summary-target` 或 `--target-tree-filepath`、`--labels`、`--decimals`、`--output`，最后是输入文件），不经过 shell 运行。`summary_target` 取 `consensus`（默认）、`mcct` 或 `msct`；`target_tree` 则对给定的一棵树的各个进化枝打分，两者同时给出会被 `ConfigurationError` 拒绝，因为 sumtrees 在这种组合下会崩溃。总结结果由 sumtrees 写入新的输出文件。
- **API：** `summarize_tree_samples(tree_files: 'Sequence[str]', output: 'str', *, input_format: 'str | None' = None, burnin: 'int' = 0, rooted: 'bool | None' = None, weighted_trees: 'bool' = False, summary_target: 'str | None' = None, target_tree: 'str | None' = None, min_clade_freq: 'float | None' = None, edges: 'str | None' = None, labels: 'str' = 'support', percentages: 'bool' = False, decimals: 'int' = 8, output_format: 'str | None' = None, sumtrees_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 dendropy 5.1（额外依赖 `dnakit[bio-dendropy]`），它会安装 `sumtrees` 命令行脚本；脚本从 PATH、当前解释器所在目录或 `sumtrees_path` 查找。
- **输入：** `tree_files` 是已存在树文件的非空列表；`output` 是新路径；`input_format` 取 `nexus/newick`、`nexus`、`newick`、`phylip` 或 `nexml`；`burnin` 是非负整数（每个文件开头丢弃的树数）；`rooted` 取 `True`（强制有根）、`False`（强制无根）或 `None`；`weighted_trees`、`percentages` 是布尔值；`summary_target` 或 `target_tree`；`min_clade_freq` 在 0 到 1 之间；`edges` 取 `mean-length`、`median-length`、`mean-age`、`median-age`、`support`、`keep` 或 `clear`；`labels` 取 `support`、`keep` 或 `clear`；`decimals` 是非负整数；`output_format` 取 `nexus`、`newick`、`phylip` 或 `nexml`；`sumtrees_path` 可选。
- **输出：** 字典，包含 `status`（退出码为 0 且输出文件存在时为 `completed`，否则为 `failed`）、`backend`、`version`、`command`、`return_code`、`outputs` 和 `log_tail`。
- **示例代码：**

```python
import re
import tempfile
from pathlib import Path

from dnakit.phylogeny import summarize_tree_samples

with tempfile.TemporaryDirectory() as folder:
    output = Path(folder) / "summary.nwk"
    result = summarize_tree_samples(
        ["trees.nwk"], str(output), input_format="newick", rooted=True, output_format="newick", verbose=False
    )
    print(result["status"], result["backend"], result["return_code"])
    print(re.findall(r"\)([0-9.]+):", output.read_text()))
```

- **示例输出：**

```text
completed sumtrees 0
['0.66666667', '0.66666667', '1.00000000']
```

## 16) `dnakit.phylogeny.transfer_tree_labels()`

- **作用：** 把源树在匹配分支（split）上携带的标签复制到一棵目标树上。
- **处理方法：** 为 `sumlabels` 应用构造固定的参数列表（`--quiet`、`--target`、`--separator=...`、可选标志、`--output`，最后是源文件），不经过 shell 运行。目标树的每个 split 会得到源树中同一 split 上的标签，标签之间用 `separator` 连接。输出默认是带类群块的 NEXUS，`to_newick=True` 时为 Newick。
- **API：** `transfer_tree_labels(tree_files: 'Sequence[str]', target_tree: 'str', output: 'str', *, separator: 'str' = '/', preserve_target_labels: 'bool' = False, rooted: 'bool | None' = None, include_taxa_block: 'bool' = True, additional_comments: 'str | None' = None, to_newick: 'bool' = False, sumlabels_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 dendropy 5.1（额外依赖 `dnakit[bio-dendropy]`），它会安装 `sumlabels` 命令行脚本；脚本从 PATH、当前解释器所在目录或 `sumlabels_path` 查找。
- **输入：** `tree_files` 是已存在、内部节点带标签的树文件的非空列表；`target_tree` 是已存在的树文件；`output` 是新路径；`separator` 是非空单行字符串（默认 `/`）；`preserve_target_labels`、`include_taxa_block`、`to_newick` 是布尔值；`rooted` 取 `True`、`False` 或 `None`；`additional_comments` 是可选的单行字符串；`sumlabels_path` 可选。
- **输出：** 字典，包含 `status`（退出码为 0 且输出文件存在时为 `completed`，否则为 `failed`）、`backend`、`version`、`command`、`return_code`、`outputs` 和 `log_tail`。`version` 可能是 `None`，因为 `sumlabels --version` 不输出普通的版本号。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.phylogeny import transfer_tree_labels

with tempfile.TemporaryDirectory() as folder:
    output = Path(folder) / "labelled.nwk"
    result = transfer_tree_labels(["trees.nwk"], "target.nwk", str(output), rooted=True, to_newick=True, verbose=False)
    print(result["status"], result["backend"], result["return_code"])
    print(output.read_text().strip())
```

- **示例输出：**

```text
completed sumlabels 0
[&R] ((t1,t2)0.9/0.7,(t3,t4)0.8/0.6);
```

## 17) `dnakit.microbiome.khmer_build_countgraph()`

- **作用：** 把 FASTA 或 FASTQ 读段的规范（canonical）k-mer 计入 khmer countgraph 并保存。
- **处理方法：** 用 `khmer.ReadParser` 读取文件，把每条读段转为大写，将长度不小于 `ksize` 且只含 A、C、G、T 的读段送入 `Countgraph.consume`（`small_count=True` 时用 `SmallCountgraph`）。更短的读段和含其他符号的读段会被跳过，并计入结果。表用 `save()` 保存到新路径。计数是近似的：哈希冲突只会抬高计数，k-mer 与它的反向互补序列共用一个计数。
- **API：** `khmer_build_countgraph(source: 'str', output_path: 'str', *, ksize: 'int' = 32, table_size: 'int' = 1000000, n_tables: 'int' = 4, small_count: 'bool' = False, max_reads: 'int' = 1000000, max_bases: 'int' = 100000000, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 khmer 2.1（额外依赖 `dnakit[bio-khmer]`），不需要外部程序。
- **输入：** `source` 是已存在的 FASTA 或 FASTQ 文件；`output_path` 是已存在文件夹中的新路径；`ksize` 为 1 到 32（默认 32）；`table_size` 是每张表的桶数（默认 1,000,000）；`n_tables` 为 1 到 20（默认 4）；`small_count` 是布尔值；`max_reads`（默认 1,000,000）和 `max_bases`（默认 100,000,000）是显式上限，输入超过时抛出 `ConfigurationError`。
- **输出：** 字典，包含 `status`、`source`、`output`、`kind`（`countgraph` 或 `smallcountgraph`）、`ksize`、`table_size`、`n_tables`、`hashsizes`（每张表的素数桶数）、`occupied_bins`、`kmers_consumed`、`reads`、`used_reads`、`skipped_non_acgt_reads` 和 `skipped_short_reads`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

import khmer

from dnakit.microbiome import khmer_build_countgraph

with tempfile.TemporaryDirectory() as folder:
    output = Path(folder) / "reads.ct"
    result = khmer_build_countgraph("reads.fa", str(output), ksize=8, table_size=1000, n_tables=2, verbose=False)
    print(result["kind"], result["hashsizes"], result["occupied_bins"], result["kmers_consumed"])
    print(result["reads"], result["used_reads"], result["skipped_short_reads"], result["skipped_non_acgt_reads"])
    table = khmer.load_countgraph(str(output))
    print(table.get("ACGTACGT"), table.get("ACGTACGA"))
```

- **示例输出：**

```text
countgraph [997, 991] 11 23
5 3 1 1
4 0
```

## 18) `dnakit.microbiome.khmer_build_nodegraph()`

- **作用：** 把 FASTA 或 FASTQ 读段中出现过的规范 k-mer 记录到 khmer nodegraph（存在性表）并保存。
- **处理方法：** 做法与 `khmer_build_countgraph` 相同，但构建的是每个桶只占一位的 `khmer.Nodegraph`。它回答“是否出现过”的问题，误报率较小，取决于表大小与 k-mer 数量。表用 `save()` 保存到新路径，可用 `khmer.load_nodegraph` 读回。
- **API：** `khmer_build_nodegraph(source: 'str', output_path: 'str', *, ksize: 'int' = 32, table_size: 'int' = 1000000, n_tables: 'int' = 4, max_reads: 'int' = 1000000, max_bases: 'int' = 100000000, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 khmer 2.1（额外依赖 `dnakit[bio-khmer]`），不需要外部程序。
- **输入：** `source`、`output_path`、`ksize`、`table_size`、`n_tables`、`max_reads`、`max_bases` 与 `khmer_build_countgraph` 相同（没有 `small_count`）。
- **输出：** 与 `khmer_build_countgraph` 相同的字典，其中 `kind` 为 `nodegraph`。`kmers_consumed` 统计送入的 k-mer 数，而不是不同 k-mer 的数量。
- **示例代码：**

```python
import tempfile
from pathlib import Path

import khmer

from dnakit.microbiome import khmer_build_nodegraph

with tempfile.TemporaryDirectory() as folder:
    output = Path(folder) / "reads.nodegraph"
    result = khmer_build_nodegraph("reads.fa", str(output), ksize=8, table_size=1000, n_tables=2, verbose=False)
    print(result["kind"], result["hashsizes"], result["occupied_bins"], result["kmers_consumed"])
    table = khmer.load_nodegraph(str(output))
    print(table.get("ACGTACGT"), table.get("ACGTACGA"))
```

- **示例输出：**

```text
nodegraph [997, 991] 11 23
1 0
```

## 19) `dnakit.microbiome.khmer_table_size_from_memory()`

- **作用：** 根据内存预算或表大小，计算 khmer 表每张表的桶数和估计内存。
- **处理方法：** 读取 khmer 的 `_buckets_per_byte`（countgraph 为 1，smallcountgraph 为 2，nodegraph 为 8）。给出 `max_memory_bytes` 时，表大小为 `int(multiplier * buckets_per_byte * max_memory_bytes / n_tables)`，与 khmer 的 `calculate_graphsize` 一致；给出 `max_tablesize` 时直接使用。内存估计为 `n_tables * table_size / buckets_per_byte`。不会构建任何表，也不写文件。
- **API：** `khmer_table_size_from_memory(graph_type: "Literal['countgraph', 'smallcountgraph', 'nodegraph']", *, max_memory_bytes: 'int | None' = None, max_tablesize: 'int | None' = None, n_tables: 'int' = 4, multiplier: 'float' = 1.0, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 khmer 2.1（额外依赖 `dnakit[bio-khmer]`），不需要外部程序。
- **输入：** `graph_type` 取 `countgraph`、`smallcountgraph` 或 `nodegraph`；`max_memory_bytes`（正数）与 `max_tablesize`（1 到 10**12 的整数）必须恰好给出一个；`n_tables` 为 1 到 20（默认 4）；`multiplier` 是缩放由内存推出的表大小的正数。
- **输出：** 字典，包含 `status`、`operation`、`graph_type`、`buckets_per_byte`、`n_tables`、`table_size` 和 `estimated_memory_bytes`（浮点数）。
- **示例代码：**

```python
from dnakit.microbiome import khmer_table_size_from_memory

for graph_type in ("countgraph", "smallcountgraph", "nodegraph"):
    result = khmer_table_size_from_memory(graph_type, max_memory_bytes=1_000_000, verbose=False)
    print(result["graph_type"], result["buckets_per_byte"], result["table_size"], result["estimated_memory_bytes"])
```

- **示例输出：**

```text
countgraph 1 250000 1000000.0
smallcountgraph 2 500000 1000000.0
nodegraph 8 2000000 1000000.0
```
