# 84. 第三方包接口 G16b：SBOL3、SnapGene、DNA Features Viewer、pyGenomeViz、primer3-py

本页包装五个上游包中与 DNA 或基因相关、可离线运行的公开函数，统一放在 `dnakit.molecular_lab` 下，共 12 个函数。每个函数返回 JSON 可序列化的字典，包含 `status`、`operation`、结果字段和 `provenance`。

- 可选依赖：`bio-pysbol3`、`bio-snapgene`、`bio-dna-features-viewer`、`bio-pygenomeviz`、`bio-primer3`。上游包只在函数内部导入，缺失时报 `BackendUnavailableError`。
- 许可证：primer3-py 为 GPL-2.0，仅作为可选导入，dnakit 不复制其代码。
- 网络：`fetch_genbank_by_accession()` 访问 NCBI；测试中替换为模拟客户端，默认测试不发起网络请求。
- SnapGene：本地没有 `.dna` 样例文件，正向读取测试标记为跳过，并在本页说明原因。
- 未包装的函数：绘图与布局函数、命令行入口、示例数据下载、需要外部程序（MMseqs2、MUMmer、progressiveMauve）的工作流、已弃用的别名，以及已由 `backends.py` 或 `dnakit.ops` 等现有接口覆盖的函数。

## copy_sbol3_design() {#84-sbol3-001}

- **作用：** 把本地 SBOL3 文件中的顶层对象复制到新的命名空间，并可把结果写成新文件。
- **计算方法：** 调用 `sbol3.copy` 克隆顶层对象，保留显示 ID 和路径。`sbol3.copy(into_namespace=...)` 在引用对象先于被引用对象处理时，会保留旧命名空间中的引用；本函数在整组对象上再做一次引用更新，结果与输入顺序无关。
- **API：** `dnakit.molecular_lab.copy_sbol3_design(source, *, into_namespace, output=None, verbose=True)`
- **输入：** SBOL 文件路径（`.ttl`、`.nt`、`.rdf`、`.xml` 或 `.jsonld`）。`into_namespace` 必须是绝对 URI；`output` 是新文件路径，后缀决定序列化格式。
- **输出：** `copied_count`，以及按标识排序的 `copied` 列表（类型、原标识、新标识、显示 ID）；还有 `output`、`output_format` 和 `provenance`。
- **示例文件：** `design.ttl`

```turtle
@prefix sbol: <http://sbols.org/v3#> .

<https://example.org/a/part1> a sbol:Component ;
    sbol:displayId "part1" ;
    sbol:hasNamespace <https://example.org/a> ;
    sbol:hasSequence <https://example.org/a/seq1> ;
    sbol:type <https://identifiers.org/SBO:0000251> .

<https://example.org/a/seq1> a sbol:Sequence ;
    sbol:displayId "seq1" ;
    sbol:elements "ATGCGTACGATCG" ;
    sbol:encoding <https://identifiers.org/edam:format_1207> ;
    sbol:hasNamespace <https://example.org/a> .
```

- **示例代码：**

```python
import pathlib

from dnakit.molecular_lab import copy_sbol3_design

pathlib.Path("copied.nt").unlink(missing_ok=True)
result = copy_sbol3_design(
    "design.ttl", into_namespace="https://example.org/b", output="copied.nt"
)
print(result["copied_count"], result["output_format"])
for row in result["copied"]:
    print(row["type"], row["identity"])
```

- **示例输出：**

```text
2 NTRIPLES
Component https://example.org/b/part1
Sequence https://example.org/b/seq1
```

## to_sbol3_display_id() {#84-sbol3-002}

- **作用：** 把任意文本转换为有效的 SBOL3 displayId，便于为 DNA 部件命名。
- **计算方法：** 调用 `sbol3.string_to_display_id` 替换不允许的字符，数字开头时加下划线；结果再用 `sbol3.is_valid_display_id` 检查。空白文本报错。
- **API：** `dnakit.molecular_lab.to_sbol3_display_id(name, *, verbose=True)`
- **输入：** 非空字符串，最长 256 个字符。
- **输出：** `display_id`、`valid`、`rule` 和 `provenance`。

- **示例代码：**

```python
from dnakit.molecular_lab import to_sbol3_display_id

print(to_sbol3_display_id("pTac promoter 1")["display_id"])
result = to_sbol3_display_id("1abc")
print(result["display_id"], result["valid"])
```

- **示例输出：**

```text
pTac_promoter_1
_1abc True
```

## check_sbol3_display_id() {#84-sbol3-003}

- **作用：** 检查字符串是否为有效的 SBOL3 displayId。
- **计算方法：** 调用 `sbol3.is_valid_display_id`。规则：只能包含字母、数字和下划线，且不能以数字开头。
- **API：** `dnakit.molecular_lab.check_sbol3_display_id(display_id, *, verbose=True)`
- **输入：** 字符串（空字符串返回无效）。
- **输出：** `display_id`、`valid`、`rule` 和 `provenance`。

- **示例代码：**

```python
from dnakit.molecular_lab import check_sbol3_display_id

for text in ["pTac_1", "1bad", "a-b"]:
    print(text, check_sbol3_display_id(text)["valid"])
```

- **示例输出：**

```text
pTac_1 True
1bad False
a-b False
```

## extract_sbol3_display_id() {#84-sbol3-004}

- **作用：** 从 SBOL3 URL 标识中取出 displayId，即 URL 路径的最后一段。
- **计算方法：** 调用 `sbol3.extract_display_id`。非 URL 标识返回 `None`；路径末段不是有效 displayId 时抛出 `ConfigurationError`。
- **API：** `dnakit.molecular_lab.extract_sbol3_display_id(identity, *, verbose=True)`
- **输入：** SBOL 对象的标识字符串。
- **输出：** `identity`、`display_id`、`found` 和 `provenance`。

- **示例代码：**

```python
from dnakit.molecular_lab import extract_sbol3_display_id

print(extract_sbol3_display_id("https://example.org/lab/pTac")["display_id"])
print(extract_sbol3_display_id("plain text")["found"])
```

- **示例输出：**

```text
pTac
False
```

## snapgene_seqrecord_summary() {#84-snapgene-001}

- **作用：** 用 snapgene_reader 读取 SnapGene `.dna` 文件，返回 Biopython SeqRecord 视图的摘要（序列和特征）。
- **计算方法：** 调用 `snapgene_file_to_seqrecord`。序列长度和特征数量受 `max_bases`、`max_features` 限制；特征按位置排序。只读，不写文件。缺失文件或无法解析的文件抛出 `ConfigurationError`。
- **API：** `dnakit.molecular_lab.snapgene_seqrecord_summary(source, *, max_bases=2000000, max_features=10000, verbose=True)`
- **输入：** 本地 `.dna` 文件路径。
- **输出：** `record_id`、`name`、`length`、`sequence`、`feature_count` 和按位置排序的 `features`（类型、起止、链、标签）。
- **说明：** 本地环境没有 SnapGene `.dna` 样例文件，因此正向读取测试标记为跳过；缺失文件和损坏文件的测试保留。下面的示例只展示错误路径。

- **示例代码：**

```python
import pathlib

from dnakit.exceptions import ConfigurationError
from dnakit.molecular_lab import snapgene_seqrecord_summary

pathlib.Path("bad.dna").write_bytes(b"not a SnapGene file")
for name in ["missing.dna", "bad.dna"]:
    try:
        snapgene_seqrecord_summary(name)
    except ConfigurationError as exc:
        print(name, exc.code, exc.message.split(":")[0])
```

- **示例输出：**

```text
missing.dna CONFIGURATION_ERROR Input file does not exist
bad.dna CONFIGURATION_ERROR snapgene_reader could not parse bad.dna
```

## annotate_sequence_feature() {#84-dfv-001}

- **作用：** 为一条 DNA 序列添加一个注释特征，并返回这个特征。
- **计算方法：** 构造 Biopython SeqRecord，调用 `dna_features_viewer.annotate_biopython_record`。位置为 0-based、起点包含、终点不包含的区间。
- **API：** `dnakit.molecular_lab.annotate_sequence_feature(sequence, *, start=0, end=None, strand=1, feature_type="misc_feature", qualifiers=None, verbose=True)`
- **输入：** DNA 字符串（支持 IUPAC 简并码）。`end` 缺省为序列末尾；`strand` 取 1、-1 或 0；`qualifiers` 是字符串到字符串的映射。
- **输出：** `sequence_length` 和 `feature`（类型、起止、链、qualifiers），以及 `provenance`。

- **示例代码：**

```python
from dnakit.molecular_lab import annotate_sequence_feature

result = annotate_sequence_feature(
    "ATGCCCGGGTAA", start=3, end=9, feature_type="CDS", qualifiers={"label": "core"}
)
feature = result["feature"]
print(result["sequence_length"], feature["type"], feature["start"], feature["end"])
print(feature["strand"], feature["qualifiers"]["label"])
```

- **示例输出：**

```text
12 CDS 3 9
1 core
```

## graphical_translation() {#84-dfv-002}

- **作用：** 翻译 DNA 子序列，结果按图中从左到右的顺序排列。
- **计算方法：** 调用 `dna_features_viewer.biotools.extract_graphical_translation`。负链先取反向互补再翻译，结果反转，使氨基酸顺序与正向序列一致。
- **API：** `dnakit.molecular_lab.graphical_translation(sequence, *, start=0, end=None, strand=1, long_form=False, verbose=True)`
- **输入：** DNA 字符串；`strand` 取 1 或 -1；`long_form=True` 时返回三字母氨基酸名列表。
- **输出：** `translation`、`residue_count`、`start`、`end`、`strand`、`long_form` 和 `provenance`。

- **示例代码：**

```python
from dnakit.molecular_lab import graphical_translation

print(graphical_translation("ATGGCCTAA")["translation"])
print(graphical_translation("TTAGCCAT", end=8, strand=-1)["translation"])
print(graphical_translation("TTAGCCATG", strand=-1, long_form=True)["translation"])
```

- **示例输出：**

```text
MA*
AM
['*', 'Gly', 'His']
```

## summarize_genbank_records() {#84-pgv-001}

- **作用：** 用 pyGenomeViz 的 GenBank 解析器读取一个或多个 GenBank 文件，汇总基因组长度、范围、GC 百分比和记录。
- **计算方法：** 调用 `pygenomeviz.scripts.gbk_files2gbk_objects`。每个输入可写成 `路径:起点-终点`（取区间）或 `路径::-1`（取反向互补）。记录按序列号排序。
- **API：** `dnakit.molecular_lab.summarize_genbank_records(sources, *, verbose=True)`
- **输入：** GenBank 文件路径列表，1 到 100 个。
- **输出：** `genome_count` 和 `genomes`（名称、全长、范围、区间长度、反向标记、GC 百分比、记录），以及 `provenance`。
- **示例文件：** `genome.gb`

```text
LOCUS       TEST1                     49 bp    DNA              UNK 01-JAN-1980
DEFINITION  Example genome fixture.
ACCESSION   TEST1
VERSION     TEST1.1
KEYWORDS    .
SOURCE      .
  ORGANISM  .
            .
FEATURES             Location/Qualifiers
     CDS             1..12
ORIGIN
        1 atgcgtacga tcgtttaagg ccatgcaagt cgatcgatcg ggatcctag
//
```

- **示例代码：**

```python
from dnakit.molecular_lab import summarize_genbank_records

result = summarize_genbank_records(["genome.gb", "genome.gb:10-30"])
for genome in result["genomes"]:
    print(genome["name"], genome["full_length"], genome["range"], genome["gc_percent"])
```

- **示例输出：**

```text
genome 49 [0, 49] 51.020408
genome 49 [10, 30] 45.0
```

## fetch_genbank_by_accession() {#84-pgv-002}

- **作用：** 按登录号从 NCBI 获取一条核酸 GenBank 记录。这是网络函数。
- **计算方法：** 调用 `pygenomeviz.utils.fetch_genbank_by_accid`（NCBI Entrez efetch），再用 Biopython 把返回文本解析为单条记录。登录号和邮箱在联网前校验；请求失败抛出 `BackendExecutionError`。测试中用模拟的 Entrez 客户端，不发起真实请求。
- **API：** `dnakit.molecular_lab.fetch_genbank_by_accession(accession, *, output=None, email=None, verbose=True)`
- **输入：** 登录号，例如 `NC_000913.3`、`JX128258.1`；`output` 为新文件路径；`email` 为可选的联系邮箱。
- **输出：** `accession`、`record_id`、`name`、`length`、`output`、`bytes` 和 `provenance`。

- **示例代码：**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.molecular_lab import fetch_genbank_by_accession

try:
    fetch_genbank_by_accession("bad accession")
except ConfigurationError as exc:
    print(exc.code, exc.message)
```

- **示例输出：**

```text
CONFIGURATION_ERROR accession must look like NC_000913.3, JX128258.1 or CP000001.
```

## calculate_primer3_end_stability() {#84-p3-001}

- **作用：** 计算引物 3' 端与模板的末端稳定性。
- **计算方法：** 调用 `primer3.calc_end_stability`（SantaLucia 最近邻模型）。返回 Tm、ΔG、ΔH、ΔS；未找到结构时数值为 `None`。primer3-py 为 GPL-2.0，仅作为可选导入。
- **API：** `dnakit.molecular_lab.calculate_primer3_end_stability(primer, template, *, mv_conc=50.0, dv_conc=1.5, dntp_conc=0.6, dna_conc=50.0, temp_c=37.0, max_loop=30, verbose=True)`
- **输入：** 引物和模板各 1 到 200 个碱基（A、C、G、T、N）。`mv_conc`、`dv_conc`、`dntp_conc` 单位为 mM，`dna_conc` 单位为 nM。
- **输出：** `structure_found`、`tm_c`、`dg_cal_per_mol`、`dh_cal_per_mol`、`ds_cal_per_mol_k`、`conditions` 和 `provenance`。

- **示例代码：**

```python
from dnakit.molecular_lab import calculate_primer3_end_stability

result = calculate_primer3_end_stability("GCTAGCTCAGTCCTAGG", "CCTAGGACTGAGCTAGC")
print(result["structure_found"], round(result["tm_c"], 2), round(result["dg_cal_per_mol"], 1))
```

- **示例输出：**

```text
True 52.24 -17299.9
```

## design_primer3_pairs() {#84-p3-002}

- **作用：** 用 primer3 为 DNA 模板设计 PCR 引物对。
- **计算方法：** 调用 `primer3.design_primers`。可用 `target=(起点, 长度)` 指定扩增子必须覆盖的区域（0-based）。返回的引物对按 primer3 罚分从低到高排序。反向引物的 5' 起点由其 3' 端位置和长度推出。
- **API：** `dnakit.molecular_lab.design_primer3_pairs(template, *, target=None, min_product=100, max_product=300, num_return=5, min_size=18, opt_size=20, max_size=25, min_tm=57.0, opt_tm=60.0, max_tm=63.0, min_gc=40.0, max_gc=60.0, verbose=True)`
- **输入：** 50 到 10000 bp 的模板（A、C、G、T、N）。
- **输出：** `pair_count` 和 `pairs`。每对包含左右引物的序列、起点、长度、Tm、GC 百分比，以及产物长度和罚分。
- **说明：** 示例模板含有自互补片段，因此第 2 对的左右引物序列相同。primer3 仍按罚分返回该对；实际使用时应检查引物二聚体和自身互补性。

- **示例代码：**

```python
from dnakit.molecular_lab import design_primer3_pairs

template = "ATGCGTACGATCGATCGTAGCTAGCTAGCATCGATCGTAGCTACGATCGTAGCTAGCTAGCTAGCTGATCGATCGATCGTAGCTAGCTAGCATGCATGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCATCGATGCATGCATGCTAGCTAGCTAGCATCGATCGATGCTAGCTAGCTAGCTAGCATGCTAGCTAGCTAGCTACGATCGATCGATCGTAGCTACGATCGA"
result = design_primer3_pairs(template, min_product=80, max_product=200, num_return=2)
print(result["pair_count"])
for pair in result["pairs"]:
    print(pair["rank"], pair["product_size"], pair["left"]["sequence"], pair["right"]["sequence"])
```

- **示例输出：**

```text
2
1 189 AGCTAGCATCGATCGTAGCT CGATCGATCGTAGCTAGCTAGC
2 142 CGATCGATCGTAGCTAGCTAGC CGATCGATCGTAGCTAGCTAGC
```

## sanitize_primer3_sequence() {#84-p3-003}

- **作用：** 把 IUPAC 简并码替换为 N，使序列符合 primer3 的输入要求。
- **计算方法：** 调用 `primer3.p3helpers.sanitize_sequence`。A、C、G、T、N 保持原样（包括大小写）；X、Z、空位等非 IUPAC 字符报错。
- **API：** `dnakit.molecular_lab.sanitize_primer3_sequence(sequence, *, verbose=True)`
- **输入：** 非空字符串。
- **输出：** `input_length`、`sanitized`、`changed_count` 和 `provenance`。

- **示例代码：**

```python
from dnakit.molecular_lab import sanitize_primer3_sequence

result = sanitize_primer3_sequence("ACGTRYN")
print(result["sanitized"], result["changed_count"])
```

- **示例输出：**

```text
ACGTNNN 2
```
