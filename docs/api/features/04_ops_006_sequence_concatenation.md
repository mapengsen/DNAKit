# 序列拼接与组装

将多个 DNA 序列进行拼接，或根据精确重叠区域进行去重叠拼接。

## 1) 普通序列拼接

- **作用：** 按输入顺序连接多条 DNA 片段，可在连接处插入 linker 或结构化 Gap，并返回各片段在新序列中的坐标，适合构建组合序列。
- **处理方法：** 按输入顺序连接片段，并在相邻片段之间插入指定 linker 或 Gap。
- **API：** `dnakit.ops.concat(sequences[必须], linker[可选], gap[可选])`。
- **输入：** 至少两条 `DNASequence` 或规范化字符串；`linker` 和 `gap` 不能同时使用。
- **输出：** 返回拼接后的 `DNASequence`。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.ops import concat

joined = concat(
    [DNASequence("AA"), DNASequence("CC")],
    linker="T",
)
print(joined.symbols)
```

- **示例结果：**

```text
AATCC
```

## 2) 去重叠拼接

- **作用：** 根据两段 DNA 的末端重叠关系完成连接，使已匹配的重叠区域只保留一次，并报告接头位置和最终序列。
- **处理方法：** 选择最长的“左片段后缀 = 右片段前缀”精确重叠，只保留一次；IUPAC 字符按字面比较。
- **API：** `dnakit.ops.concat_overlap(sequences[必须], min_overlap[可选], max_overlap[可选])`。
- **输入：** 必须是两条线性、无显式 Gap 的 `DNASequence` 或规范化字符串。 `min_overlap` 默认 1，`max_overlap` 可限制重叠搜索长度。
- **输出：** 返回处理后的 `DNASequence`，序列内容可从 `symbols` 读取。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.ops import concat_overlap

joined = concat_overlap(
    [DNASequence("AAACCC"), DNASequence("CCCGG")],
)
print(joined.symbols)
```

- **示例结果：**

```text
AAACCCGG
```

<!-- BEGIN GENERATED API: cloning -->

<span id="cloning"></span>

## 3) 设计 Gibson 无缝组装重叠序列 {#dnakit.databases.dna_gibson_design}

- **作用：** 设计 Gibson 无缝组装重叠序列。
- **处理方法：** 从相邻片段末端提取指定长度的同源区，添加到待组装片段。
- **API：** `dnakit.databases.dna_gibson_design(*, operation=None, fragments, overlap_length=None, verbose=True)`
- **输入：** 必须提供 `fragments`。
- **输出：** 返回字典，包含 `assembly_fragments`, `n_fragments`, `overlap_length` 等字段。
- **示例代码：**

```python
import dnakit

example_fragments = ['ATGGAGGAGCCGCAGTCAGATCCTAGCGTTGAATTCGGATCCAAGCTTACGTACGT',
 'CGTACGTAAGCTTGGATCCGAATTCAAGCGATAGGATCTGACTGCGGCTCCTCCAT',
 'TCCATGGAGGAGCCGCAGTCAGATCCTAGCGTTGAATTCGGATCCAAGCTTACGTACGT']

result = dnakit.databases.dna_gibson_design(
    operation='gibson_design',
    fragments=example_fragments,
    overlap_length=20,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "assembly_fragments": [
    {
      "name": "<string>",
      "original_sequence": "<string>",
      "with_overlaps": "<string>",
      "left_overlap": "<string>",
      "right_overlap": "<string>"
    }
  ],
  "n_fragments": "<integer>",
  "overlap_length": "<integer>"
}
```

[完整参数与使用说明](../../faq/api/04_ops_006_sequence_concatenation.md#dnakit.databases.dna_gibson_design-parameters)

## 4) 推演 Golden Gate 反应的组装产物 {#dnakit.databases.dna_golden_gate_assemble}

- **作用：** 推演 Golden Gate 反应的组装产物。
- **处理方法：** 模拟 Type IIS 酶切，移除仍含识别位点的片段，再按突出端连接。
- **API：** `dnakit.databases.dna_golden_gate_assemble(*, operation=None, fragments, enzyme=None, circular=None, labels=None, verbose=True)`
- **输入：** 必须提供 `fragments`。
- **输出：** 返回字典，包含 `enzyme`, `product_length`, `num_fragments_assembled`, `assembly_order` 等字段。
- **示例代码：**

```python
import dnakit

example_fragments = ['GGTCTCAAAACATGGAGGAGCCGCAGTCAGATCCTAGCGTTGAATTCGGATCCCTTTTGAGACC',
 'GGTCTCAAAAGGGATCCAAGCTTACGTACGTATGGAGGAGCCGCAGTCAGATATTTTGAGACC',
 'GGTCTCAAAATCCTAGCGTTGAATTCGGATCCAAGCTTACGTACGTATGGAGCGTTTGAGACC',
 'GGTCTCAAACGATCGATCGATCGATCGATCGGTTTTGAGACC']

result = dnakit.databases.dna_golden_gate_assemble(
    fragments=example_fragments,
    enzyme='BsaI',
    circular=True,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "enzyme": "<string>",
  "product_length": "<integer>",
  "num_fragments_assembled": "<integer>",
  "assembly_order": [
    "<value>"
  ],
  "unused_fragments": [
    "<value>"
  ],
  "product_sequence": "<string>"
}
```

[完整参数与使用说明](../../faq/api/04_ops_006_sequence_concatenation.md#dnakit.databases.dna_golden_gate_assemble-parameters)

## 5) 设计 Golden Gate 酶切位点与连接末端 {#dnakit.databases.dna_golden_gate_design}

- **作用：** 设计 Golden Gate 酶切位点与连接末端。
- **处理方法：** 分配不回文的唯一四碱基突出端，并添加 Type IIS 识别位点。
- **API：** `dnakit.databases.dna_golden_gate_design(*, operation=None, parts, enzyme=None, verbose=True)`
- **输入：** 必须提供 `parts`。
- **输出：** 返回字典，包含 `parts_with_overhangs`, `overhangs`, `enzyme`, `protocol_note` 等字段。
- **示例代码：**

```python
import dnakit

example_parts = ['ATGGAGGAGCCGCAGTCAGATCCTAGCGTTGAATTCGGATCC',
 'GGATCCAAGCTTACGTACGTATGGAGGAGCCGCAGTCAGAT',
 'CCTAGCGTTGAATTCGGATCCAAGCTTACGTACGTATGGAG']

result = dnakit.databases.dna_golden_gate_design(
    operation='golden_gate_design',
    parts=example_parts,
    enzyme='BsaI',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "parts_with_overhangs": [
    {
      "name": "<string>",
      "sequence": "<string>",
      "left_overhang": "<string>",
      "right_overhang": "<string>",
      "full_sequence": "<string>"
    }
  ],
  "overhangs": [
    "<string>"
  ],
  "enzyme": "<string>",
  "protocol_note": "<string>"
}
```

[完整参数与使用说明](../../faq/api/04_ops_006_sequence_concatenation.md#dnakit.databases.dna_golden_gate_design-parameters)

<!-- END GENERATED API: cloning -->

## 6) 使用 SPAdes 组装 reads {#bio-18}

- **作用：** 将短 reads 或短长读混合数据组装为 contigs。
- **计算方法：** 调用 SPAdes；`second` 提供双端，`nanopore` 提供混合组装长读。
- **API：** `dnakit.bioinformatics.ngs.assemble_reads(first, output_directory, second=None, nanopore=None, memory_gb=8)`
- **输入：** FASTQ 输入；SPAdes 在 PATH 中；新输出目录；可用内存满足 `memory_gb`。
- **输出：** 返回 `ToolResult`，包含后端 `backend`、产物路径 `artifacts` 和运行日志 `log`。
- **Agent/MCP API：** `dnakit.ngs.assemble_reads(read1, output_directory, read2=None, nanopore=None, memory_gb=8, ...)`；返回 JSON 可序列化的组装产物路径与日志。
- **引用与许可：** SPAdes 使用 Bankevich et al. (2012), DOI `10.1089/cmb.2012.0021`（GPL-2.0）；后端由调用者安装。
- **示例代码：**

```python
from dnakit.bioinformatics.ngs import assemble_reads

result = assemble_reads("R1.fq.gz", "new_assembly", second="R2.fq.gz", memory_gb=8)
print(sorted(result.artifacts))
```

- **示例结果：**

返回 `ToolResult`；`artifacts` 列出 SPAdes 生成的 contigs、scaffolds 等文件，名称以实际产物为准。外部程序示例未执行。

[依赖与实现说明](../../faq/bioinformatics.md#reads)。

## 7) 连接双链 DNA 片段 {#bio-52-ligate}

- **作用：** 连接具有兼容末端的双链 DNA 片段。
- **计算方法：** 检查两片段末端的极性和互补性，再合并两条链。
- **API：** `dnakit.bioinformatics.cloning.Duplex.from_sequence(sequence, circular=False)`; `Duplex.ligate(other)`
- **输入：** 两个具有兼容末端的线性 `Duplex`；示例先切割同一分子再重新连接。
- **输出：** 返回连接后的 `Duplex`，包含 Watson 链序列及末端信息。
- **示例代码：**

```python
from dnakit.bioinformatics.cloning import Duplex

molecule = Duplex.from_sequence("AACCGGTT")
left, right = molecule.cut(3, 5)
print(left.ligate(right).watson)
```

- **示例结果：**

```text
AACCGGTT
```

[依赖与实现说明](../../faq/bioinformatics.md#optimization)。

## 8) 保存克隆步骤与中间产物 {#bio-53}

- **作用：** 将克隆过程组织为可追踪的命名步骤。
- **计算方法：** `add` 注册产物，`apply` 调用函数并保存父步骤和结果。
- **API：** `dnakit.bioinformatics.cloning.CloningPlan()`; `CloningPlan.add(name, product)`; `CloningPlan.apply(name, function, parents=(), **parameters)`
- **输入：** 实际 Python 函数、父步骤名称和操作参数。
- **输出：** 返回 `CloningPlan`，示例展示 `add`、`apply`、`steps`。
- **示例代码：**

```python
from dnakit.bioinformatics.cloning import CloningPlan, Duplex

plan = CloningPlan()
plan.add("template", Duplex.from_sequence("AACCGGTT"))
cut = plan.apply("cut", lambda molecule: molecule.cut(3, 5), parents=["template"])
print([step.name for step in plan.steps])
print(len(cut.product))
```

- **示例结果：**

```text
['template', 'cut']
2
```

[依赖与实现说明](../../faq/bioinformatics.md#optimization)。

## 9) 组装同源重叠双链片段 {#bio-54}

- **作用：** 寻找输入片段可形成的线性或环形组装产物。
- **计算方法：** 使用 pydna 搜索满足最小重叠长度的组装路径。
- **API：** `dnakit.bioinformatics.cloning.assemble_duplexes(fragments, circular=False, minimum_overlap=25, max_assemblies=50)`
- **输入：** `Duplex` 片段列表；需要 `bio-assembly`；重叠长度根据实验设计设置。
- **输出：** 返回可行组装产物的 `Duplex` 元组。
- **示例代码：**

```python
from dnakit.bioinformatics.cloning import Duplex, assemble_duplexes

fragments = [Duplex.from_sequence("TTATACGATGCTAGCTTACGGTCA"),
             Duplex.from_sequence("AGCTTACGGTCACCTATGAGTCCA")]
products = assemble_duplexes(fragments, minimum_overlap=12)
longest = max(products, key=lambda product: len(product.watson))
print(longest.watson)
```

- **示例结果：**

```text
TTATACGATGCTAGCTTACGGTCACCTATGAGTCCA
```

[依赖与实现说明](../../faq/bioinformatics.md#optimization)。

## 10) Golden Gate 部件设计 {#13-golden-gate}

- **作用：** 为多个 DNA 部件分配唯一非回文连接 overhang，并添加 BsaI/BbsI 位点。
- **处理方法：** 分配唯一非回文 overhang，并拼接 BsaI/BbsI 识别位点。
- **API：** `dnakit.molbio.design_golden_gate(parts[必须], enzyme[可选])`
- **输入：** 至少两条 canonical DNA；设计阶段支持 `BsaI` 或 `BbsI`。
- **输出：** 返回 `ProviderResult`，查询内容在 `data`，并附有 `provider` 和来源信息。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.molbio import design_golden_gate

result = design_golden_gate(
    (DNASequence("ATGGCTGAA"), DNASequence("GCCAAATAA")),
    enzyme="BsaI",
)
print(result.data["parts_with_overhangs"])
```

- **示例结果：** 每个部件的左右 overhang、完整合成序列和实验说明。

[可选后端依赖与返回值说明](../../faq/tutorial_integrations.md#optional-bioinformatics)

## 11) Golden Gate 反应组装 {#14-golden-gate}

- **作用：** 模拟 Type IIS 酶切、按 overhang 匹配连接，并返回组装产物和片段顺序。
- **处理方法：** 模拟 Type IIS 酶切，再按 overhang 兼容性连接片段。
- **API：** `dnakit.molbio.assemble_golden_gate(fragments[必须], enzyme[可选], circular[可选], labels[可选])`
- **输入：** 至少两条带相应 Type IIS 位点的 canonical DNA；支持 BsaI、BbsI、Esp3I/BsmBI 和 SapI。
- **输出：** 返回 `ProviderResult`，查询内容在 `data`，并附有 `provider` 和来源信息。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.molbio import assemble_golden_gate

result = assemble_golden_gate(
    (
        DNASequence("GGTCTCAAAACATGGAGGAGCCGCAGTCAGATCCTAGCGTTGAATTCGGATCCCTTTTGAGACC"),
        DNASequence("GGTCTCAAAAGGGATCCAAGCTTACGTACGTATGGAGGAGCCGCAGTCAGATATTTTGAGACC"),
        DNASequence("GGTCTCAAAATCCTAGCGTTGAATTCGGATCCAAGCTTACGTACGTATGGAGCGTTTGAGACC"),
        DNASequence("GGTCTCAAACGATCGATCGATCGATCGATCGGTTTTGAGACC"),
    ),
    labels=("vector", "insert-1", "insert-2", "insert-3"),
)
print(result.data["product_sequence"])
```

- **示例结果：** 产物序列、长度、assembly order、junction overhang 和未使用片段。

[可选后端依赖与返回值说明](../../faq/tutorial_integrations.md#optional-bioinformatics)
