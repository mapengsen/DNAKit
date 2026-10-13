# 可选分子实验室后端：DNA 设计、记录与绘图

这一组函数是 `dnakit.molecular_lab` 中基于可选依赖的适配器，覆盖 [DnaChisel](https://github.com/Edinburgh-Genome-Foundry/DnaChisel)、[DnaCauldron](https://github.com/Edinburgh-Genome-Foundry/DnaCauldron)、[pydna](https://github.com/BjornFJohansson/pydna)、[pySBOL3](https://github.com/SynBioDex/pySBOL3)、[snapgene_reader](https://github.com/Edinburgh-Genome-Foundry/SnapGeneReader)、[DNA Features Viewer](https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer)、[pyGenomeViz](https://github.com/moshi4/pyGenomeViz) 和 [primer3-py](https://github.com/libnano/primer3-py) 的 DNA 层面功能。它们计算序列约束评分、优化记录、设计引物与重叠尾巴、查找重组位点、构建和读取 SBOL3 文档、读取 SnapGene 文件、绘制特征图和基因组轨道，并计算引物的 Tm、发夹与二聚体热力学。它们不预测湿实验结果，不访问远程数据库，也不替代原生的 Golden Gate 和 Gibson 模拟器。

primer3-py 使用 GPL-2.0 许可。DNAKit 只在调用 `calculate_primer3_*()` 时导入它，不复制其任何代码。

依赖都是可选的，只有调用对应函数时才会导入。缺少包时会报错并给出安装命令。示例从 `dnakit.molecular_lab` 导入。

## 示例数据

唯一需要声明的输入是一个单条记录的 FASTA 文件，其中 `GGTCTC` 位于第 3 至 9 位（0 起始）。其余示例直接传入输入，或把 GenBank 与 SBOL3 文件写入临时目录。

**示例文件：** `part.fa`

```text
>part1
ATGGGTCTCATGCCCAAAGCTTTAA
```

## 1) `dnakit.molecular_lab.evaluate_dnachisel_specs()`

- **作用：** 不修改序列，只对 DnaChisel 的约束与目标打分。
- **处理方法：** 按请求构建 DnaChisel 规格，用 `constraints_evaluations()` 和 `objectives_evaluations()` 评估。每一行报告该规格是否通过、分数、0 起始的位置和 DnaChisel 的消息。
- **API：** `dnakit.molecular_lab.evaluate_dnachisel_specs(sequence, *, constraints=(), objectives=(), circular=False)`。
- **输入：** DNA 字符串。每个约束或目标是一个映射，其中 `name` 必须属于允许列表（`AvoidPattern`、`AvoidChanges`、`AvoidHairpins`、`AvoidStopCodons`、`CodonOptimize`、`EnforceChanges`、`EnforceGCContent`、`EnforceMeltingTemperature`、`EnforceTranslation`、`MaximizeCAI`），其余键是关键字参数。
- **输出：** `all_constraints_pass`、`failed_constraints`，以及 `constraints` 和 `objectives` 两组行。
- **示例代码：**

```python
from dnakit.molecular_lab import evaluate_dnachisel_specs

result = evaluate_dnachisel_specs(
    "ATGGGTCTCATGCCCAAA",
    constraints=[
        {"name": "AvoidPattern", "pattern": "GGTCTC"},
        {"name": "EnforceGCContent", "mini": 0.3, "maxi": 0.7, "window": 18},
    ],
)
print(result["all_constraints_pass"], result["failed_constraints"])
for row in result["constraints"]:
    print(row["name"], row["passes"], row["locations"])
```

- **示例输出：**

```text
False ['AvoidPattern[0-18](pattern:GGTCTC)']
AvoidPattern False [[3, 9]]
EnforceGCContent True []
```

## 2) `dnakit.molecular_lab.optimize_dnachisel_record()`

- **作用：** 优化本地 FASTA、GenBank 或 SnapGene 记录，并把结果写成带注释的 GenBank 文件。
- **处理方法：** 用 DnaChisel 读取记录，解析约束；若给出目标则优化目标；再用 `to_record()` 写出结果。输出路径必须是新文件。
- **API：** `dnakit.molecular_lab.optimize_dnachisel_record(source, *, output_path, constraints=(), objectives=(), circular=False)`。
- **输入：** 本地记录文件、新的 `.gb`、`.gbk` 或 `.genbank` 输出路径，以及与上相同的规格映射。
- **输出：** 序列长度、编辑次数、是否所有约束都通过，以及输出路径。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.molecular_lab import optimize_dnachisel_record

output = Path(tempfile.mkdtemp()) / "part_optimized.gb"
result = optimize_dnachisel_record(
    "part.fa",
    output_path=output,
    constraints=[{"name": "AvoidPattern", "pattern": "GGTCTC"}],
)
print(result["sequence_length"], result["edits"], result["all_constraints_pass"])
```

- **示例输出：**

```text
25 1 True
```

## 3) `dnakit.molecular_lab.reverse_translate_protein()`

- **作用：** 把蛋白质反向翻译为 DNA。
- **处理方法：** 调用 DnaChisel 的 `reverse_translate()`，每个氨基酸取密码子表中的第一个密码子。结果是确定的。末尾的 `*` 会变成终止密码子。不使用密码子使用偏好。
- **API：** `dnakit.molecular_lab.reverse_translate_protein(protein, *, table="Standard", max_residues=10000)`。
- **输入：** 由 20 种标准氨基酸组成的蛋白质字符串，末尾可以有 `*`；`table` 是 DnaChisel 列表中的密码子表名称。
- **输出：** `sequence`、`sequence_length`、`protein_length` 和 `table`。
- **示例代码：**

```python
from dnakit.molecular_lab import reverse_translate_protein

result = reverse_translate_protein("MAK*")
print(result["sequence"], result["sequence_length"])
```

- **示例输出：**

```text
ATGGCTAAATAA 12
```

## 4) `dnakit.molecular_lab.read_dnacauldron_parts()`

- **作用：** 通过 DnaCauldron 读取本地零件文件，并列出注释中找到的重叠末端。
- **处理方法：** 对每个文件调用 DnaCauldron 的 `load_record()`，再调用 `list_overhangs_from_record_annotations()`。零件 ID 必须唯一，总序列长度有上限。
- **API：** `dnakit.molecular_lab.read_dnacauldron_parts(sources, *, topology="default_to_linear", max_parts=50)`。
- **输入：** GenBank、FASTA 或 SnapGene 文件路径列表。
- **输出：** 每个零件的 `id`、`file`、`length`、`topology`、`overhangs` 和 `sequence`。
- **示例代码：**

```python
from dnakit.molecular_lab import read_dnacauldron_parts

result = read_dnacauldron_parts(["part.fa"])
part = result["parts"][0]
print(part["id"], part["length"], part["topology"])
```

- **示例输出：**

```text
part1 25 linear
```

## 5) `dnakit.molecular_lab.write_dnacauldron_constructs()`

- **作用：** 模拟 DnaCauldron 的 Gibson 或 Type IIS 组装，并把每个构建体写成带注释的 GenBank 文件。
- **处理方法：** 用指定的零件构建 DnaCauldron 仓库，运行 `GibsonAssembly` 或 `Type2sRestrictionAssembly`，再用 `write_record()` 写出每个构建体。输出目录由本函数创建，调用前不能存在。`dnakit.integrations` 中的 `simulate_dnacauldron_assembly()` 只返回序列；本函数会把特征注释保留在文件中。
- **API：** `dnakit.molecular_lab.write_dnacauldron_constructs(parts, *, output_dir, method="gibson", enzyme="auto", max_constructs=20)`。
- **输入：** 从零件名到 DNA 的映射、新的输出目录、`method`（`"gibson"` 或 `"type2s"`）；当 `method="type2s"` 时还需要 Type IIS 酶名。
- **输出：** `constructs` 列表，包含每个文件名、长度和拓扑结构，以及 DnaCauldron 给出的警告和错误。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.molecular_lab import write_dnacauldron_constructs

a = "ATGGCTAAAGGCGAATTGTTTACCGGCGTGGTCCCGATTC"
s1 = "TGGAACTGGATGGCGATGTG"
b = "AACGGCCATAAATTTAGCGTGAGCGGCGAAGGCGAAGGTG"
s2 = "CCGGGTTCAAGTCTCCAACA"
c = "TTAGGCCAACGTTCGATTAGCATCGATCGA"
parts = {"A": a + s1, "B": s1 + b + s2, "C": s2 + c + a[:20]}
result = write_dnacauldron_constructs(
    parts, output_dir=Path(tempfile.mkdtemp()) / "constructs"
)
construct = result["constructs"][0]
print(len(result["constructs"]), construct["length"], construct["topology"])
```

- **示例输出：**

```text
1 150 circular
```

## 6) `dnakit.molecular_lab.design_pcr_primers()`

- **作用：** 为线性 DNA 模板设计一对跨越全长的正向和反向引物。
- **处理方法：** 对模板调用 pydna 的 `primer_design()`。引物位于模板两端，因此产物就是整个模板，`product_length` 等于模板长度。
- **API：** `dnakit.molecular_lab.design_pcr_primers(template, *, target_tm=55.0, min_primer_length=13)`。
- **输入：** 至少含 `2 * min_primer_length` 个碱基的 DNA 字符串。`target_tm` 是 40 至 75 °C 之间的温度。
- **输出：** `forward_primer`、`reverse_primer`、两者的长度以及 `product_length`。
- **示例代码：**

```python
from dnakit.molecular_lab import design_pcr_primers

template = "ATGGCTAAAGGCGAATTGTTTACCGGCGTGGTCCCGATTCTGGTGGAACTGGATGGCGATGTGAACGG"
result = design_pcr_primers(template)
print(result["forward_primer"], result["reverse_primer"], result["product_length"])
```

- **示例输出：**

```text
ATGGCTAAAGGCGAATTG CCGTTCACATCGCC 68
```

## 7) `dnakit.molecular_lab.design_assembly_primers()`

- **作用：** 为多个 PCR 片段的引物加上重叠尾巴，使相邻片段可以通过 Gibson 组装或融合 PCR 连接。
- **处理方法：** 先用 pydna 的 `primer_design()` 为每个模板设计一对引物，再调用 pydna 的 `assembly_fragments()`。第 i 个片段的反向引物的反向互补序列，末尾就是第 i+1 个片段正向引物的前 `overlap` 个碱基。
- **API：** `dnakit.molecular_lab.design_assembly_primers(templates, *, overlap=35, max_link=40, circular=False)`。
- **输入：** 2 至 20 个 DNA 模板，重叠长度为 10 至 100 个碱基，`max_link` 不小于 `overlap`。
- **输出：** 每个片段的 `forward_primer`、`reverse_primer` 和 `product_length`。
- **示例代码：**

```python
from dnakit.molecular_lab import design_assembly_primers

fragment_1 = (
    "ATGGCTAAAGGCGAATTGTTTACCGGCGTGGTCCCGATTCTGGTGGAACTGGATGGCGATGTGAACGGCCATAAATTTAGC"
)
fragment_2 = (
    "GTGAGCGGCGAAGGCGAAGGTGATGCACCGTACCCGGAACTGACCCTGAAATTCATTTGCACCACCGGCAAACTGCCGGTG"
)
result = design_assembly_primers([fragment_1, fragment_2], overlap=20)
for item in result["fragments"]:
    print(item["index"], item["forward_primer"], item["reverse_primer"], item["product_length"])
```

- **示例输出：**

```text
0 ATGGCTAAAGGCGAATTG CGCCGCTCACGCTAAATTTATGGCCGTTCA 91
1 TAAATTTAGCGTGAGCGGCGAAG CACCGGCAGTTTGC 91
```

## 8) `dnakit.molecular_lab.find_recombination_sites()`

- **作用：** 在 DNA 序列中查找 loxP 位点或 Gateway att 位点。
- **处理方法：** `system="loxp"` 时调用 pydna 的 `find_loxP_sites()`，`system="gateway"` 时调用 `find_gateway_sites()`。结果按起始位置排序。
- **API：** `dnakit.molecular_lab.find_recombination_sites(sequence, *, system="loxp", greedy=False)`。
- **输入：** DNA 字符串。`system` 为 `"loxp"` 或 `"gateway"`；`greedy` 传给 Gateway 查找器。
- **输出：** `site_count` 和 `sites`。每个位点包含 `name`、`start`、`end`（0 起始、末端不含）和 `strand`。
- **示例代码：**

```python
from dnakit.molecular_lab import find_recombination_sites

loxp = find_recombination_sites("GGG" + "ATAACTTCGTATAGCATACATTATACGAAGTTAT" + "CCC")
print(loxp["site_count"], loxp["sites"])
attb = find_recombination_sites("AAACAAGTTTGTACAAAAAAGCAGGCTCCCC", system="gateway")
print(attb["site_count"], attb["sites"])
```

- **示例输出：**

```text
1 [{'name': 'loxP', 'start': 3, 'end': 37, 'strand': 1}]
1 [{'name': 'attB1', 'start': 3, 'end': 24, 'strand': 1}]
```

## 9) `dnakit.molecular_lab.build_sbol3_document()`

- **作用：** 构建 SBOL3 文档，每个命名的 DNA 序列对应一个 Component 和一个 Sequence。
- **处理方法：** 为每个显示 ID 创建带 IUPAC DNA 编码的 pySBOL3 `Sequence` 和 DNA 类型的 `Component`。标识符为 `namespace + display_id` 与 `namespace + display_id + "_sequence"`。序列化前先执行文档校验。
- **API：** `dnakit.molecular_lab.build_sbol3_document(components, *, namespace, file_format="turtle", output_path=None)`。
- **输入：** 从合法 SBOL 显示 ID 到 DNA 字符串的映射；以 `/` 或 `#` 结尾的 `http(s)` 命名空间；`file_format` 可取 `turtle`、`ntriples`、`rdfxml` 或 `jsonld`。
- **输出：** `component_count`、`sequence_bases`、`validation`（错误和警告的数量及详情），以及未指定 `output_path` 时的 `document` 文本，或指定时的 `output_path`。
- **示例代码：**

```python
from dnakit.molecular_lab import build_sbol3_document

result = build_sbol3_document(
    {"promoter_1": "TTGACAGCTAGCTCAGTCCTAGGTATAATGCTAGC", "cds_1": "ATGGCTAAAGGCGAATAA"},
    namespace="https://dnakit.example.org/parts/",
)
print(result["component_count"], result["sequence_bases"], result["validation"]["error_count"])
```

- **示例输出：**

```text
2 53 0
```

## 10) `dnakit.molecular_lab.read_sbol3_sequences()`

- **作用：** 读取本地 SBOL3 文档中保存的 DNA 序列。
- **处理方法：** 用 pySBOL3 的 `Document.read()` 读取文件，返回每个顶层 `Sequence` 对象的标识、显示 ID、编码和序列内容。
- **API：** `dnakit.molecular_lab.read_sbol3_sequences(source, *, max_objects=50000, max_sequences=1000)`。
- **输入：** 本地 SBOL3 文件，格式可由 pySBOL3 根据扩展名推断。
- **输出：** `object_count`、`sequence_count` 和 `sequences`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.molecular_lab import build_sbol3_document, read_sbol3_sequences

path = Path(tempfile.mkdtemp()) / "parts.ttl"
build_sbol3_document(
    {"promoter_1": "TTGACAGCTAGC", "cds_1": "ATGGCTAAATAA"},
    namespace="https://dnakit.example.org/parts/",
    output_path=path,
)
for row in read_sbol3_sequences(path)["sequences"]:
    print(row["display_id"], row["length"], row["elements"])
```

- **示例输出：**

```text
cds_1_sequence 12 ATGGCTAAATAA
promoter_1_sequence 12 TTGACAGCTAGC
```

## 11) `dnakit.molecular_lab.validate_sbol3_document()`

- **作用：** 校验本地 SBOL3 文档，并列出其中的错误和警告。
- **处理方法：** 用 pySBOL3 读取文件并调用 `Document.validate()`。每个问题都报告对象、规则和消息。
- **API：** `dnakit.molecular_lab.validate_sbol3_document(source, *, max_issues=200)`。
- **输入：** 本地 SBOL3 文件。
- **输出：** `valid`、`error_count`、`warning_count`，以及最多 `max_issues` 条 `errors` 和 `warnings`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.molecular_lab import build_sbol3_document, validate_sbol3_document

path = Path(tempfile.mkdtemp()) / "parts.ttl"
build_sbol3_document(
    {"promoter_1": "TTGACAGCTAGC"},
    namespace="https://dnakit.example.org/parts/",
    output_path=path,
)
checked = validate_sbol3_document(path)
print(checked["valid"], checked["error_count"], checked["warning_count"])
```

- **示例输出：**

```text
True 0 0
```

## 12) `dnakit.molecular_lab.read_snapgene_sequence()`

- **作用：** 从 SnapGene `.dna` 文件读取序列、拓扑结构和甲基化标记。
- **处理方法：** 调用 snapgene_reader 的 `snapgene_file_to_dict()`，并把字段映射为 DNAKit 的名称。snapgene_reader 只读取文件，DNAKit 不写入 `.dna` 文件。
- **API：** `dnakit.molecular_lab.read_snapgene_sequence(source, *, max_bases=2000000)`。
- **输入：** 由 SnapGene 写出的 `.dna` 文件。
- **输出：** `sequence`、`length`、`topology`、`strandedness`、三个甲基化标记，以及导出和导入版本号。
- **示例代码（未本地运行）：**

```python
from dnakit.molecular_lab import read_snapgene_sequence

result = read_snapgene_sequence("pUC19.dna")
print(result["length"], result["topology"], result["sequence"][:20])
```

- **示例输出（未本地运行）：**

```text
<length> <topology> <first 20 bases>
```

本例未在本地运行，因为仓库中没有由 SnapGene 写出的 `.dna` 文件，检查器无法读取。输出块只是模板，不是实际记录的结果。

## 13) `dnakit.molecular_lab.read_snapgene_features()`

- **作用：** 列出 SnapGene `.dna` 文件中的特征及其坐标。
- **处理方法：** 读取 snapgene_reader 的 `features` 块，把每个起点转换为 0 起始、末端不含的区间。标签来自 `label` 限定符；缺少该限定符时使用特征名。
- **API：** `dnakit.molecular_lab.read_snapgene_features(source, *, max_features=10000)`。
- **输入：** SnapGene `.dna` 文件。
- **输出：** `feature_count` 和 `features`，每个特征包含 `name`、`type`、`start`、`end`、`strand`、`color` 和 `label`。
- **示例代码（未本地运行）：**

```python
from dnakit.molecular_lab import read_snapgene_features

result = read_snapgene_features("pUC19.dna")
for row in result["features"][:3]:
    print(row["label"], row["start"], row["end"], row["strand"])
```

- **示例输出（未本地运行）：**

```text
<label> <start> <end> <strand>
```

本例未在本地运行，原因与上一个函数相同。

## 14) `dnakit.molecular_lab.convert_snapgene_to_genbank()`

- **作用：** 把 SnapGene `.dna` 文件转换为新的 GenBank 文件。
- **处理方法：** 调用 snapgene_reader 的 `snapgene_file_to_gbk()`。GenBank 文本由该转换器生成，其 LOCUS 日期也由它固定。失败时会删除未完成的输出文件。
- **API：** `dnakit.molecular_lab.convert_snapgene_to_genbank(source, output_path)`。
- **输入：** SnapGene `.dna` 文件，以及新的 `.gb`、`.gbk` 或 `.genbank` 输出路径。
- **输出：** `output_path` 和写入的字节数。
- **示例代码（未本地运行）：**

```python
from dnakit.molecular_lab import convert_snapgene_to_genbank

result = convert_snapgene_to_genbank("pUC19.dna", "pUC19.gb")
print(result["output_path"], result["bytes"])
```

- **示例输出（未本地运行）：**

```text
<output path> <bytes>
```

本例未在本地运行，原因与上一个函数相同。

## 15) `dnakit.molecular_lab.build_feature_map()`

- **作用：** 返回线性序列的特征图数据：坐标、颜色，以及每个特征所在的泳道。
- **处理方法：** 创建 DNA Features Viewer 的 `GraphicFeature` 对象，并调用 `compute_features_levels()`。重叠的特征会被移到更高的泳道。泳道从 0.0 开始。
- **API：** `dnakit.molecular_lab.build_feature_map(sequence_length, features, *, max_features=10000)`。
- **输入：** 序列长度，以及特征映射列表。每个映射包含 `start` 和 `end`（0 起始、末端不含），可选 `strand`（1、-1 或 0）、`label` 和 `color`（`#RRGGBB`）。
- **输出：** `lane_count` 和 `features`，每个特征包含坐标、链方向、颜色、标签和 `lane`。
- **示例代码：**

```python
from dnakit.molecular_lab import build_feature_map

result = build_feature_map(
    100,
    [
        {"start": 0, "end": 40, "label": "promoter"},
        {"start": 30, "end": 60, "strand": -1, "label": "tag"},
        {"start": 70, "end": 90, "label": "terminator"},
    ],
)
for row in result["features"]:
    print(row["label"], row["start"], row["end"], row["lane"])
```

- **示例输出：**

```text
promoter 0 40 0.0
tag 30 60 1.0
terminator 70 90 0.0
```

## 16) `dnakit.molecular_lab.render_feature_map_svg()`

- **作用：** 绘制线性特征图，并写入新的 SVG 文件。
- **处理方法：** 在使用 Agg 画布的 Matplotlib 图上绘制 DNA Features Viewer 的 `GraphicRecord`，因此不需要显示设备。SVG 不含日期元数据，元素 ID 使用固定的哈希盐，所以同一输入得到相同字节。
- **API：** `dnakit.molecular_lab.render_feature_map_svg(sequence_length, features, *, output_path, figure_width=8.0, max_features=10000)`。
- **输入：** 与 `build_feature_map()` 相同的特征，以及新的 `.svg` 输出路径。`figure_width` 为 2 至 20 英寸。
- **输出：** `output_path`、`format`、`bytes`、`feature_count` 和 `sequence_length`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.molecular_lab import render_feature_map_svg

path = Path(tempfile.mkdtemp()) / "map.svg"
result = render_feature_map_svg(
    100,
    [
        {"start": 0, "end": 40, "label": "promoter"},
        {"start": 30, "end": 60, "strand": -1, "label": "tag"},
    ],
    output_path=path,
)
print(result["format"], result["feature_count"], path.read_bytes()[:5])
```

- **示例输出：**

```text
svg 2 b'<?xml'
```

## 17) `dnakit.molecular_lab.build_genbank_feature_map()`

- **作用：** 返回本地 GenBank 记录的完整特征图数据，并给出每个特征的泳道。
- **处理方法：** 用 DNA Features Viewer 的 `BiopythonTranslator` 读取记录，再用 `compute_features_levels()` 计算泳道。`dnakit.integrations` 中的 `summarize_dna_features_viewer()` 只返回前几个特征的预览。
- **API：** `dnakit.molecular_lab.build_genbank_feature_map(source, *, max_features=10000)`。
- **输入：** 本地 GenBank 文件。
- **输出：** `sequence_length`、`feature_count`、`lane_count` 和 `features`，每个特征包含 `label`、坐标、`strand`、`color` 和 `lane`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from dnakit.molecular_lab import build_genbank_feature_map

record = SeqRecord(Seq("A" * 100), id="map1", name="map1", description="example")
record.annotations["molecule_type"] = "DNA"
record.features = [
    SeqFeature(FeatureLocation(10, 40, strand=1), type="CDS", qualifiers={"gene": ["geneX"]}),
    SeqFeature(FeatureLocation(30, 60, strand=-1), type="gene", qualifiers={"gene": ["geneY"]}),
]
path = Path(tempfile.mkdtemp()) / "map.gb"
SeqIO.write(record, str(path), "genbank")
result = build_genbank_feature_map(path)
print(result["sequence_length"], result["lane_count"])
print([(row["label"], row["start"], row["end"]) for row in result["features"]])
```

- **示例输出：**

```text
100 2
[('geneX', 10, 40), ('geneY', 30, 60)]
```

## 18) `dnakit.molecular_lab.read_genbank_genome_tracks()`

- **作用：** 把本地 GenBank 的每条记录读取为一个基因组轨道，包含序列长度和所选特征。
- **处理方法：** 用 pyGenomeViz 的 `Genbank` 解析器获取长度，再从解析后的记录中筛选一种特征。筛选由 DNAKit 自己的循环完成，因为 pyGenomeViz 0.4.4 的 `get_seqid2features()` 读取 `SeqFeature.strand`，而 Biopython 1.88 已不再提供该属性。多段位置取整体跨度。
- **API：** `dnakit.molecular_lab.read_genbank_genome_tracks(source, *, feature_type="CDS", max_records=100, max_features=50000)`。
- **输入：** 本地 GenBank 文件，以及一种特征键，例如 `"CDS"` 或 `"gene"`。
- **输出：** `track_count`、`feature_count` 和 `tracks`。每个轨道包含 `name`、`length` 和 `features`，每个特征包含 `type`、`start`、`end`、`strand` 和 `label`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from dnakit.molecular_lab import read_genbank_genome_tracks

record = SeqRecord(Seq("A" * 100), id="map1", name="map1", description="example")
record.annotations["molecule_type"] = "DNA"
record.features = [
    SeqFeature(FeatureLocation(10, 40, strand=1), type="CDS", qualifiers={"gene": ["geneX"]}),
    SeqFeature(FeatureLocation(30, 60, strand=-1), type="gene", qualifiers={"gene": ["geneY"]}),
]
path = Path(tempfile.mkdtemp()) / "map.gb"
SeqIO.write(record, str(path), "genbank")
result = read_genbank_genome_tracks(path, feature_type="CDS")
track = result["tracks"][0]
print(result["track_count"], track["length"], track["features"])
```

- **示例输出：**

```text
1 100 [{'type': 'CDS', 'start': 10, 'end': 40, 'strand': 1, 'label': 'geneX'}]
```

## 19) `dnakit.molecular_lab.render_genome_tracks_svg()`

- **作用：** 绘制特征轨道，可选地在相邻轨道之间绘制连接，并写入新的 SVG 文件。
- **处理方法：** 为每个条目添加一条 pyGenomeViz 特征轨道和其中的特征，再用 `add_link()` 添加连接。连接只允许出现在按给定顺序相邻的两条轨道之间，因为 pyGenomeViz 只在那里绘制连接。pyGenomeViz 会在 SVG 中写入创建时间戳和随机的分组 ID，因此每次运行的字节不同。
- **API：** `dnakit.molecular_lab.render_genome_tracks_svg(tracks, *, output_path, links=(), fig_width=15.0, max_tracks=20)`。
- **输入：** 轨道列表，每条包含 `name`、`length` 和 `features`。可选的 `links` 是包含 `track_a`、`start_a`、`end_a`、`track_b`、`start_b` 和 `end_b` 的映射。
- **输出：** `output_path`、`format`、`bytes`、`track_count`、`feature_count` 和 `link_count`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from dnakit.molecular_lab import read_genbank_genome_tracks, render_genome_tracks_svg

record = SeqRecord(Seq("A" * 100), id="map1", name="map1", description="example")
record.annotations["molecule_type"] = "DNA"
record.features = [
    SeqFeature(FeatureLocation(10, 40, strand=1), type="CDS", qualifiers={"gene": ["geneX"]}),
]
path = Path(tempfile.mkdtemp()) / "map.gb"
SeqIO.write(record, str(path), "genbank")
track = read_genbank_genome_tracks(path)["tracks"][0]
other = {"name": "other", "length": 120, "features": [{"start": 5, "end": 60, "strand": -1}]}
result = render_genome_tracks_svg(
    [track, other],
    output_path=Path(tempfile.mkdtemp()) / "tracks.svg",
    links=[
        {
            "track_a": track["name"],
            "start_a": 10,
            "end_a": 40,
            "track_b": "other",
            "start_b": 5,
            "end_b": 60,
        }
    ],
)
print(result["format"], result["track_count"], result["feature_count"], result["link_count"])
```

- **示例输出：**

```text
svg 2 2 1
```

## 20) `dnakit.molecular_lab.calculate_primer3_tm()`

- **作用：** 用 primer3-py 计算一条引物的熔解温度。
- **处理方法：** 调用 primer3-py 的 `calc_tm()`，采用 SantaLucia 最近邻模型，并传入一价盐、二价盐、dNTP 和引物浓度。primer3-py 为 GPL-2.0；DNAKit 只在运行本函数时导入它。
- **API：** `dnakit.molecular_lab.calculate_primer3_tm(sequence, *, mv_conc=50.0, dv_conc=1.5, dntp_conc=0.6, dna_conc=50.0)`。
- **输入：** 只含 A、C、G、T 的 1 至 200 个碱基引物。浓度中，一价、二价盐和 dNTP 以 mM 为单位，引物以 nM 为单位。
- **输出：** 以摄氏度表示的 `tm_c`、`sequence_length`，以及所用的 `conditions`。
- **示例代码：**

```python
from dnakit.molecular_lab import calculate_primer3_tm

result = calculate_primer3_tm("ACGTTGCAAGCTTGCA")
print(result["sequence_length"], round(result["tm_c"], 2))
```

- **示例输出：**

```text
16 54.92
```

## 21) `dnakit.molecular_lab.calculate_primer3_hairpin()`

- **作用：** 用 primer3-py 报告一条引物的发夹稳定性。
- **处理方法：** 调用 primer3-py 的 `calc_hairpin()`。若 primer3 没有找到发夹，热力学字段为 None。能量以 cal/mol 报告。
- **API：** `dnakit.molecular_lab.calculate_primer3_hairpin(sequence, *, mv_conc=50.0, dv_conc=1.5, dntp_conc=0.6, dna_conc=50.0, temp_c=37.0, max_loop=30)`。
- **输入：** 与 `calculate_primer3_tm()` 相同的引物。`temp_c` 为 0 至 100 °C，`max_loop` 为 3 至 100。
- **输出：** `structure_found`、`tm_c`、`dg_cal_per_mol`、`dh_cal_per_mol` 和 `ds_cal_per_mol_k`。
- **示例代码：**

```python
from dnakit.molecular_lab import calculate_primer3_hairpin

result = calculate_primer3_hairpin("GGGGAAAACCCC")
print(result["structure_found"], round(result["tm_c"], 2), round(result["dg_cal_per_mol"], 1))
```

- **示例输出：**

```text
True 55.6 -1312.3
```

## 22) `dnakit.molecular_lab.calculate_primer3_dimer()`

- **作用：** 用 primer3-py 报告一条引物的自身二聚体；若给出第二条引物，则报告异源二聚体。
- **处理方法：** 未给出 `partner` 时调用 primer3-py 的 `calc_homodimer()`，否则调用 `calc_heterodimer()`。若 primer3 没有找到二聚体，热力学字段为 None。
- **API：** `dnakit.molecular_lab.calculate_primer3_dimer(sequence, *, partner=None, mv_conc=50.0, dv_conc=1.5, dntp_conc=0.6, dna_conc=50.0, temp_c=37.0, max_loop=30)`。
- **输入：** 一条引物，可选的 `partner` 引物，以及与 `calculate_primer3_hairpin()` 相同的条件。
- **输出：** `dimer_type`（`"homodimer"` 或 `"heterodimer"`）、`structure_found`、`tm_c`、`dg_cal_per_mol`、`dh_cal_per_mol`、`ds_cal_per_mol_k` 和 `partner_length`。
- **示例代码：**

```python
from dnakit.molecular_lab import calculate_primer3_dimer

result = calculate_primer3_dimer("ATGCGTACGT", partner="ACGTACGCAT")
print(result["dimer_type"], result["structure_found"], round(result["tm_c"], 2))
```

- **示例输出：**

```text
heterodimer True 31.61
```

## 安装可选依赖

```bash
python -m pip install "dnakit[bio-dnachisel]"            # dnachisel：函数 1 至 3
python -m pip install "dnakit[bio-dnacauldron]"          # dnacauldron：函数 4 和 5
python -m pip install "dnakit[bio-assembly]"             # pydna：函数 6 至 8
python -m pip install "dnakit[bio-pysbol3]"              # sbol3：函数 9 至 11
python -m pip install "dnakit[bio-snapgene]"             # snapgene_reader：函数 12 至 14
python -m pip install "dnakit[bio-dna-features-viewer]" "dnakit[bio-plots]"  # 函数 15 至 17
python -m pip install "dnakit[bio-pygenomeviz]"          # pyGenomeViz：函数 18 和 19
python -m pip install "dnakit[bio-primer3]"              # primer3-py：函数 20 至 22
```

各包的项目地址：[DnaChisel](https://github.com/Edinburgh-Genome-Foundry/DnaChisel)、[DnaCauldron](https://github.com/Edinburgh-Genome-Foundry/DnaCauldron)、[pydna](https://github.com/BjornFJohansson/pydna)、[pySBOL3](https://github.com/SynBioDex/pySBOL3)、[snapgene_reader](https://github.com/Edinburgh-Genome-Foundry/SnapGeneReader)、[DNA Features Viewer](https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer)、[pyGenomeViz](https://github.com/moshi4/pyGenomeViz)、[primer3-py](https://github.com/libnano/primer3-py)。

## 验证方式

测试文件是 `others/tests/unit/test_molecular_lab_backends.py`，共 41 项：40 项运行，1 项在未设置 `DNAKIT_SNAPGENE_TEST_FILE`（指向由 SnapGene 写出的 `.dna` 文件）时跳过。检查内容包括：

- 手工计算：DnaChisel 模式位置 `[3, 9)`、第一密码子反向翻译 `ATGGCTAAATAA`、loxP 位点 `[3, 37)`、三个特征的泳道 `[0, 1, 0]`、370 bp 的环状 Gibson 产物、303 bp 的 GenBank 轨道，以及 `design_assembly_primers()` 相邻引物之间 20 个碱基的重叠。
- 与上游直接调用的结果比较：DnaChisel 约束评估、pydna 引物设计和 Gateway 位点查找、DNA Features Viewer 的泳道、pyGenomeViz 的 GenBank 长度解析，以及 primer3-py 的 Tm、发夹和二聚体函数。
- 通过真实上游完成解析和校验：SBOL3 的读写往返，以及含未标注类型 Component 的无效文档；snapgene_reader 拒绝非 SnapGene 字节。错误的序列、错误的命名空间、pyGenomeViz 中不相邻的连接和已存在的输出路径都会被拒绝。
- 确定性：两次渲染 SVG 特征图得到的字节完全相同。
- 缺少任一可选包时，八个包都会抛出带额外依赖名的 `BackendUnavailableError`。

## 当前范围

这一批包装的是现有适配器尚未覆盖的 DNA 层面操作：DnaChisel 的约束评分和记录输出；DnaCauldron 的文件输入和带注释的 GenBank 输出；pydna 的引物设计、重叠尾巴和位点查找；SBOL3 的构建、序列读取和校验；SnapGene 的序列与特征读取及格式转换；DNA Features Viewer 的特征图数据和 SVG 输出；pyGenomeViz 的 GenBank 轨道和 SVG 输出；以及 primer3 的 Tm、发夹和二聚体数值。

未包装的内容及原因：

- 已在 `dnakit.integrations` 中存在的适配器：基于字符串的 DnaChisel 优化（`optimize_with_dnachisel()`）、基于字符串的 DnaCauldron 模拟（`simulate_dnacauldron_assembly()`）、SBOL3 清单（`summarize_sbol3_document()`）和特征图预览（`summarize_dna_features_viewer()`）。
- pydna 的 Gibson 产物已存在于 `dnakit.bioinformatics.cloning` 的 `assemble_duplexes()`。
- primer3-py 的 `designPrimers()`：引物设计已由 `dnakit.molbio.Primer3CLIDesignAdapter` 提供，它调用 `primer3_core` 程序。本批只包装缺失的热力学函数。
- SnapGene 写入：snapgene_reader 只读，因此不写 `.dna` 文件。
- DnaChisel 的 BLAST 辅助函数（需要网络访问）、pydna 的融合 PCR 辅助函数和凝胶模拟（`simulate_gel()` 已存在于 `dnakit.bioinformatics.cloning`）、DnaCauldron 的 PDF 报告和混合方案规划，以及 pyGenomeViz 的 HTML 导出和 GC 或偏斜曲线。最后一组需要交互式或额外的输出，留待以后的批次处理。
- 基于 SnapGene 的读取示例不在本地运行，因为仓库中没有由 SnapGene 写出的 `.dna` 文件。读取器的错误路径用非 SnapGene 文件测试。
