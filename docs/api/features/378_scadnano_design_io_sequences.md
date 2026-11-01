# scadnano 设计：文件读写、导出、序列与修饰

读写 scadnano 设计（`.sc`），转换为 cadnano v2、oxView 与 IDT 订购格式，生成板图，并处理链的 DNA 序列与修饰的函数。

这些函数是对 [scadnano](https://github.com/UC-Davis-molecular-computing/scadnano-python-package)（MIT 许可，测试版本 0.21.1）的封装，scadnano 是设计 DNA 纳米结构的 Python 包。它是可选依赖，只在函数运行时导入；用 `pip install "dnakit[dna-nanotechnology]"` 安装（Excel 订单还需要 openpyxl，板图表格需要 tabulate，二者由 scadnano 一并安装）。设计是 scadnano 的 JSON 映射（即 `.sc` 文件内容）。编辑类函数不会修改传入的参数，而是在 `design` 键下返回新设计。链用其在 `design["strands"]` 中的 0 起始位置定位，结构域用其在链内的位置定位；螺旋索引与偏移沿用 scadnano 的定义。这些函数不预测折叠、产率或稳定性。若序列是在插入、缺失或结构域改动之前指定的，新设计仍保留旧序列，因此应最后再指定序列。

## 1) `dnakit.genomics.load_scadnano_design()`

- **作用：** 从 `.sc` 文件、JSON 字符串或 JSON 映射读取并校验 scadnano 设计。
- **处理方法：** 按所给参数调用 `Design.from_scadnano_file`、`from_scadnano_json_str` 或 `from_scadnano_json_map`，再把校验过的设计序列化。
- **API：** `dnakit.genomics.load_scadnano_design(design: 'Mapping | None' = None, *, path: 'str | None' = None, json_text: 'str | None' = None) -> 'dict[str, Any]'`。
- **输入：** `design`（映射）、`path`（文件）、`json_text` 三者恰好给出一个。
- **输出：** `design`（规范化 JSON 映射）、网格、螺旋数、链数、支架链数、链名与长度、`source_kind`。
- **示例代码：**

```python
from dnakit.genomics.scadnano_p26_functions import (
    assign_dna_sequence,
    build_strand,
    load_scadnano_design,
)

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]
design = assign_dna_sequence(design, "ACGT" * 12, strand_index=0)["design"]

r = load_scadnano_design(path="design.sc")
print(r["source_kind"], r["helices"], r["strands"], r["strand_names"])
r = load_scadnano_design(json_text=open("design.sc").read())
print(r["source_kind"], r["strand_lengths"])
r = load_scadnano_design(design)
print(r["source_kind"], r["grid"])
```

- **示例输出：**

```text
file 2 2 [None, 'stapA']
json_string [48, 24]
json_map square
```

## 2) `dnakit.genomics.export_scadnano_design()`

- **作用：** 把设计序列化为 scadnano `.sc` JSON 文本，并可写入新文件。
- **处理方法：** 调用 `Design.to_json`（`suppress_indent` 为 scadnano 的紧凑排版），并像 `write_scadnano_file` 一样写出文本；已存在的输出文件会被拒绝。
- **API：** `dnakit.genomics.export_scadnano_design(design: 'Mapping', *, output_path: 'str | None' = None, suppress_indent: 'bool' = True) -> 'dict[str, Any]'`。
- **输入：** 设计映射，可选的新 `output_path`，`suppress_indent`。
- **输出：** `json_text`、`output_path` 与设计计数。
- **示例代码：**

```python
from dnakit.genomics.scadnano_p26_functions import (
    assign_dna_sequence,
    build_strand,
    export_scadnano_design,
)

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]
design = assign_dna_sequence(design, "ACGT" * 12, strand_index=0)["design"]

import os, tempfile
out = os.path.join(tempfile.mkdtemp(), "out.sc")
r = export_scadnano_design(design, output_path=out)
print(r["json_text"].splitlines()[0:3])
print(os.path.getsize(out) == len(r["json_text"].encode()))
```

- **示例输出：**

```text
['{', '  "version": "0.21.1",', '  "grid": "square",']
True
```

## 3) `dnakit.genomics.import_cadnano_v2_design()`

- **作用：** 把 cadnano v2 设计转换为 scadnano 设计。
- **处理方法：** 调用 `Design.from_cadnano_v2`，传入解析后的映射（`json_dict`）或 `.json` 文件。
- **API：** `dnakit.genomics.import_cadnano_v2_design(source: 'str | Mapping', *) -> 'dict[str, Any]'`。
- **输入：** cadnano v2 映射，或 cadnano v2 JSON 文件路径。
- **输出：** 转换后的 `design` 及其计数。
- **示例代码：**

```python
from dnakit.genomics.scadnano_p26_functions import (
    create_origami_rectangle,
    export_cadnano_v2_design,
    import_cadnano_v2_design,
)

design = create_origami_rectangle(num_helices=4, num_cols=10)["design"]

cad = export_cadnano_v2_design(design)["cadnano_v2"]
r = import_cadnano_v2_design(cad)
print(r["grid"], r["helices"], r["strands"])
```

- **示例输出：**

```text
square 4 21
```

## 4) `dnakit.genomics.export_cadnano_v2_design()`

- **作用：** 把设计导出为 cadnano v2 JSON。
- **处理方法：** 调用 `Design.to_cadnano_v2_json`（即 `write_cadnano_v2_file` 写出的内容），返回文本与解析结果，可写入新文件。
- **API：** `dnakit.genomics.export_cadnano_v2_design(design: 'Mapping', *, name: 'str' = '', whitespace: 'bool' = True, output_path: 'str | None' = None) -> 'dict[str, Any]'`。
- **输入：** 设计映射、`name`、`whitespace`、可选的新 `output_path`。
- **输出：** `json_text`、`cadnano_v2`、`output_path`。
- **示例代码：**

```python
from dnakit.genomics.scadnano_p26_functions import (
    create_origami_rectangle,
    export_cadnano_v2_design,
)

design = create_origami_rectangle(num_helices=4, num_cols=10)["design"]

r = export_cadnano_v2_design(design, name="rect")
print(sorted(r["cadnano_v2"])[:4], len(r["cadnano_v2"]["vstrands"]))
```

- **示例输出：**

```text
['name', 'vstrands'] 4
```

## 5) `dnakit.genomics.export_oxview_design()`

- **作用：** 把设计导出为 oxView 文件。
- **处理方法：** 调用 `Design.to_oxview_format`（与 `to_oxview_json`、`write_oxview_file` 内容相同），并统计结果中的链数与核苷酸数。
- **API：** `dnakit.genomics.export_oxview_design(design: 'Mapping', *, use_strand_colors: 'bool' = True, output_path: 'str | None' = None) -> 'dict[str, Any]'`。
- **输入：** 设计映射、`use_strand_colors`、可选的新 `output_path`。
- **输出：** `json_text`、`strands`、`nucleotides`、`output_path`。
- **示例代码：**

```python
from dnakit.genomics.scadnano_p26_functions import (
    assign_dna_sequence,
    build_strand,
    export_oxview_design,
)

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]
design = assign_dna_sequence(design, "ACGT" * 12, strand_index=0)["design"]

r = export_oxview_design(design)
print(r["strands"], r["nucleotides"])
```

- **示例输出：**

```text
2 72
```

## 6) `dnakit.genomics.export_idt_order()`

- **作用：** 把链导出为 IDT 批量订购文本或 96/384 孔板 Excel 工作簿。
- **处理方法：** `order_format="bulk_text"` 调用 `Design.to_idt_bulk_input_format`；`"plate_excel"` 调用 `write_idt_plate_excel_file`，需要 openpyxl 和新的 `.xlsx` 路径。修饰代码写入序列。
- **API：** `dnakit.genomics.export_idt_order(design: 'Mapping', *, order_format: 'str' = 'bulk_text', output_path: 'str | None' = None, delimiter: 'str' = ',', domain_delimiter: 'str' = '', only_strands_with_vendor_fields: 'bool' = False, export_scaffold: 'bool' = False, export_non_modified_strand_version: 'bool' = False, plate_type: 'int' = 96) -> 'dict[str, Any]'`。
- **输入：** 设计映射、`order_format`、分隔符、支架与 vendor 字段开关、`plate_type`（96 或 384）。
- **输出：** 文本模式：`text`、`lines`；孔板模式：`output_path`、`bytes`。
- **示例代码：**

```python
from dnakit.genomics.scadnano_p26_functions import (
    assign_dna_sequence,
    build_strand,
    export_idt_order,
)

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]
design = assign_dna_sequence(design, "ACGT" * 12, strand_index=0)["design"]

r = export_idt_order(design)
print(r["text"])
print(r["lines"])
```

- **示例输出：**

```text
stapA,ACGTACGTACGTACGTACGTACGT,25nm,STD
1
```

## 7) `dnakit.genomics.build_plate_maps()`

- **作用：** 生成 96 或 384 孔合成板图。
- **处理方法：** 调用 `Design.plate_maps` 与 `PlateMap.to_table`。只有 vendor 字段给出板名与孔位的链才会被放入。
- **API：** `dnakit.genomics.build_plate_maps(design: 'Mapping', *, plate_type: 'int' = 96, strand_indices: 'Sequence[int] | None' = None, table_format: 'str' = 'pipe') -> 'dict[str, Any]'`。
- **输入：** 设计映射、`plate_type`、可选 `strand_indices`、`table_format`。
- **输出：** `plates`：每板含 `plate_name`、`wells`（孔位到链名）与可打印 `table`。
- **示例代码：**

```python
from dnakit.genomics.scadnano_p26_functions import (
    assign_dna_sequence,
    build_plate_maps,
    build_strand,
)

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]
design = assign_dna_sequence(design, "ACGT" * 12, strand_index=0)["design"]

r = build_plate_maps(design)
print(r["plates"][0]["plate_name"], r["plates"][0]["wells"])
```

- **示例输出：**

```text
P1 {'A1': 'stapA'}
```

## 8) `dnakit.genomics.get_strand_sequence()`

- **作用：** 返回链或其某个结构域的 DNA 序列，可带结构域分隔符或使用供应商格式。
- **处理方法：** 按参数选择 `Strand.dna_sequence_delimited`、`Strand.vendor_dna_sequence`，以及结构域、环或延伸的 `dna_sequence` / `vendor_dna_sequence`。未指定的碱基为 `?`。
- **API：** `dnakit.genomics.get_strand_sequence(design: 'Mapping', strand_index: 'int', *, domain_index: 'int | None' = None, delimiter: 'str' = '', vendor_format: 'bool' = False) -> 'dict[str, Any]'`。
- **输入：** 设计映射、`strand_index`、可选 `domain_index`、`delimiter`、`vendor_format`。
- **输出：** `sequence`（未指定 DNA 时为 `None`）。
- **示例代码：**

```python
from dnakit.genomics.scadnano_p26_functions import (
    assign_dna_sequence,
    build_strand,
    get_strand_sequence,
)

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]
design = assign_dna_sequence(design, "ACGT" * 12, strand_index=0)["design"]

print(get_strand_sequence(design, 0, delimiter="|")["sequence"])
print(get_strand_sequence(design, 0, domain_index=1)["sequence"])
print(get_strand_sequence(design, 1, vendor_format=True)["sequence"])
```

- **示例输出：**

```text
ACGTACGTACGTACGTACGTACGT|ACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGT
```

## 9) `dnakit.genomics.get_domain_subsequence()`

- **作用：** 返回结构域在两个螺旋偏移之间的 DNA 子序列。
- **处理方法：** 调用 `Domain.dna_sequence_in`；两个偏移都包含端点，并考虑缺失与插入。
- **API：** `dnakit.genomics.get_domain_subsequence(design: 'Mapping', strand_index: 'int', domain_index: 'int', offset_left: 'int', offset_right: 'int', *) -> 'dict[str, Any]'`。
- **输入：** 设计映射、`strand_index`、`domain_index`、`offset_left`、`offset_right`。
- **输出：** `sequence` 与两个偏移。
- **示例代码：**

```python
from dnakit.genomics.scadnano_p26_functions import (
    assign_dna_sequence,
    build_strand,
    get_domain_subsequence,
)

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]
design = assign_dna_sequence(design, "ACGT" * 12, strand_index=0)["design"]

print(get_domain_subsequence(design, 0, 0, 2, 5)["sequence"])
```

- **示例输出：**

```text
GTAC
```

## 10) `dnakit.genomics.assign_dna_sequence()`

- **作用：** 给链或某个结构域指定 DNA，可同时把互补序列写到结合的链上。
- **处理方法：** 调用 `Design.assign_dna`（涵盖 `Strand.set_dna_sequence` 与构建器的 `with_sequence`）。`?` 为未知碱基通配符。
- **API：** `dnakit.genomics.assign_dna_sequence(design: 'Mapping', sequence: 'str', *, strand_index: 'int', domain_index: 'int | None' = None, assign_complement: 'bool' = True, check_length: 'bool' = False) -> 'dict[str, Any]'`。
- **输入：** 设计映射、`sequence`、`strand_index`、可选 `domain_index`、`assign_complement`、`check_length`。
- **输出：** 新的 `design`、计数与 `assigned` 链摘要。
- **示例代码：**

```python
from dnakit.genomics.scadnano_p26_functions import assign_dna_sequence, build_strand

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]

r = assign_dna_sequence(design, "ACGT" * 12, strand_index=0)
print(r["assigned"]["dna_length"], r["strand_lengths"])
print(r["design"]["strands"][1]["sequence"])
r = assign_dna_sequence(design, "AAAA", strand_index=0, domain_index=0, assign_complement=False)
print(r["design"]["strands"][0]["sequence"])
```

- **示例输出：**

```text
48 [48, 24]
ACGTACGTACGTACGTACGTACGT
AAAA????????????????????????????????????????????
```

## 11) `dnakit.genomics.assign_complement_sequence()`

- **作用：** 让一条链获得与其结合的链的 Watson-Crick 互补序列。
- **处理方法：** 调用 `Strand.assign_dna_complement_from`；只有结合区域获得碱基。
- **API：** `dnakit.genomics.assign_complement_sequence(design: 'Mapping', strand_index: 'int', other_strand_index: 'int', *) -> 'dict[str, Any]'`。
- **输入：** 设计映射、`strand_index`（接收链）、`other_strand_index`（来源链）。
- **输出：** 新的 `design` 与 `assigned` 摘要。
- **示例代码：**

```python
from dnakit.genomics.scadnano_p26_functions import (
    assign_complement_sequence,
    assign_dna_sequence,
    build_strand,
)

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]

design = assign_dna_sequence(design, "ACGT" * 12, strand_index=0, assign_complement=False)["design"]
r = assign_complement_sequence(design, 1, 0)
print(r["assigned"]["dna_sequence"])
```

- **示例输出：**

```text
ACGTACGTACGTACGTACGTACGT
```

## 12) `dnakit.genomics.assign_helix_domain_sequences()`

- **作用：** 给一条链指定某个螺旋上结构域 DNA 的拼接序列。
- **处理方法：** 调用 `Design.assign_dna_sequence`：把该螺旋上方向相反的结构域按螺旋顺序拼接并原样写入。
- **API：** `dnakit.genomics.assign_helix_domain_sequences(design: 'Mapping', strand_index: 'int', *, helix: 'int', forward: 'bool') -> 'dict[str, Any]'`。
- **输入：** 设计映射、`strand_index`、`helix`、`forward`（接收链方向）。
- **输出：** 新的 `design` 与 `assigned` 摘要。
- **示例代码：**

```python
from dnakit.genomics.scadnano_p26_functions import (
    assign_dna_sequence,
    assign_helix_domain_sequences,
    build_strand,
)

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]

design = assign_dna_sequence(design, "ACGT" * 12, strand_index=0, assign_complement=False)["design"]
r = assign_helix_domain_sequences(design, 1, helix=0, forward=False)
print(r["assigned"]["dna_sequence"])
```

- **示例输出：**

```text
ACGTACGTACGTACGTACGTACGT
```

## 13) `dnakit.genomics.check_complementarity()`

- **作用：** 判断两个碱基或两条序列是否 Watson-Crick 互补。
- **处理方法：** `mode="bases"` 调用 `bases_complementary`；`mode="sequences"` 调用 `reverse_complementary`，逐位比较两条序列（scadnano 的定义）。
- **API：** `dnakit.genomics.check_complementarity(first: 'str', second: 'str', *, mode: 'str' = 'bases', allow_wildcard: 'bool' = False) -> 'dict[str, Any]'`。
- **输入：** `first`、`second`、`mode`、`allow_wildcard`。
- **输出：** `complementary`（布尔）与 `mode`。
- **示例代码：**

```python
from dnakit.genomics.scadnano_p26_functions import check_complementarity

print(check_complementarity("A", "T")["complementary"])
print(check_complementarity("A", "G")["complementary"])
print(check_complementarity("ACG", "TGC", mode="sequences")["complementary"])
print(check_complementarity("A?", "TG", mode="sequences", allow_wildcard=True)["complementary"])
```

- **示例输出：**

```text
True
False
True
True
```

## 14) `dnakit.genomics.get_m13_sequence()`

- **作用：** 返回 M13mp18 支架序列。
- **处理方法：** 调用 `scadnano.m13`，参数为旋转量（默认 5587）与变体（`p7249`、`p7560`、`p8064`、`p8634`）。
- **API：** `dnakit.genomics.get_m13_sequence(*, rotation: 'int' = 5587, variant: 'str' = 'p7249') -> 'dict[str, Any]'`。
- **输入：** `rotation`、`variant`。
- **输出：** `sequence`、`length`、`variant`、`rotation`。
- **示例代码：**

```python
from dnakit.genomics.scadnano_p26_functions import get_m13_sequence

r = get_m13_sequence()
print(r["length"], r["sequence"][:20])
print(get_m13_sequence(rotation=0, variant="p8064")["length"])
```

- **示例输出：**

```text
7249 TTCCCTTCCTTTCTCGCCAC
8064
```

## 15) `dnakit.genomics.assign_m13_scaffold()`

- **作用：** 把 M13mp18 指定给支架链，并把互补序列指定给订书钉链。
- **处理方法：** 调用 `Design.assign_m13_to_scaffold`。
- **API：** `dnakit.genomics.assign_m13_scaffold(design: 'Mapping', *, rotation: 'int' = 5587, variant: 'str' = 'p7249') -> 'dict[str, Any]'`。
- **输入：** 含一条支架链的设计映射、`rotation`、`variant`。
- **输出：** 新的 `design`。
- **示例代码：**

```python
from dnakit.genomics.scadnano_p26_functions import assign_m13_scaffold, build_strand

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]

r = assign_m13_scaffold(design)
print(r["design"]["strands"][0]["sequence"][:20], r["design"]["strands"][1]["sequence"][:20])
```

- **示例输出：**

```text
TTCCCTTCCTTTCTCGCCAC TAGAGCTTGACGGGGAAAGC
```

## 16) `dnakit.genomics.get_strand_without_modifications()`

- **作用：** 返回去除修饰的链副本，原设计不变。
- **处理方法：** 调用 `Strand.no_modifications_version`。
- **API：** `dnakit.genomics.get_strand_without_modifications(design: 'Mapping', strand_index: 'int', *) -> 'dict[str, Any]'`。
- **输入：** 设计映射与 `strand_index`。
- **输出：** `strand`（JSON 映射）与 `dna_sequence`。
- **示例代码：**

```python
from dnakit.genomics.scadnano_p26_functions import (
    assign_dna_sequence,
    build_strand,
    get_strand_without_modifications,
    set_strand_modification,
)

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]
design = assign_dna_sequence(design, "ACGT" * 12, strand_index=0)["design"]

design = set_strand_modification(design, 1, location="5p", display_text="B", vendor_code="/5Biosg/")["design"]
print(design["strands"][1].get("5prime_modification"))
r = get_strand_without_modifications(design, 1)
print(r["strand"].get("5prime_modification"), r["dna_sequence"])
```

- **示例输出：**

```text
/5Biosg/
None ACGTACGTACGTACGTACGTACGT
```

## 17) `dnakit.genomics.inline_deletions_insertions()`

- **作用：** 把全部缺失与插入折算进结构域长度。
- **处理方法：** 调用 `Design.inline_deletions_insertions`；碱基总数不变。
- **API：** `dnakit.genomics.inline_deletions_insertions(design: 'Mapping', *) -> 'dict[str, Any]'`。
- **输入：** 设计映射。
- **输出：** 新的 `design` 与 `strand_lengths_before`。
- **示例代码：**

```python
from dnakit.genomics.scadnano_p26_functions import (
    add_design_insertion,
    build_strand,
    inline_deletions_insertions,
)

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]

design = add_design_insertion(design, 0, 10, 2)["design"]
print(design["strands"][0]["domains"][0])
r = inline_deletions_insertions(design)
print(r["design"]["strands"][0]["domains"][0], r["strand_lengths"] == r["strand_lengths_before"])
```

- **示例输出：**

```text
{'helix': 0, 'forward': True, 'start': 0, 'end': 24, 'insertions': [[10, 2]]}
{'helix': 0, 'forward': True, 'start': 0, 'end': 26} True
```

## 18) `dnakit.genomics.set_strand_modification()`

- **作用：** 给链加上 5'、3' 或内部 DNA 修饰。
- **处理方法：** 按 `location` 调用 `Strand.set_modification_5p`、`set_modification_3p` 或 `set_modification_internal`。
- **API：** `dnakit.genomics.set_strand_modification(design: 'Mapping', strand_index: 'int', *, location: 'str', display_text: 'str', vendor_code: 'str', connector_length: 'int' = 4, index: 'int | None' = None, allowed_bases: 'Sequence[str] | None' = None) -> 'dict[str, Any]'`。
- **输入：** 设计映射、`strand_index`、`location`（`5p`、`3p`、`internal`）、`display_text`、`vendor_code`、`connector_length`，内部修饰另需 `index` 与 `allowed_bases`。
- **输出：** 新的 `design` 与 `assigned` 摘要。
- **示例代码：**

```python
from dnakit.genomics.scadnano_p26_functions import (
    assign_dna_sequence,
    build_strand,
    set_strand_modification,
)

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]
design = assign_dna_sequence(design, "ACGT" * 12, strand_index=0)["design"]

r = set_strand_modification(design, 1, location="5p", display_text="B", vendor_code="/5Biosg/")
print(r["design"]["strands"][1]["5prime_modification"])
r = set_strand_modification(r["design"], 1, location="internal", display_text="T*", vendor_code="/iBiodT/", index=3, allowed_bases=["T"])
print(r["design"]["strands"][1]["internal_modifications"])
```

- **示例输出：**

```text
/5Biosg/
{'3': '/iBiodT/'}
```

## 19) `dnakit.genomics.remove_strand_modification()`

- **作用：** 移除链的 5'、3' 或内部修饰。
- **处理方法：** 调用 `Strand.remove_modification_5p`、`remove_modification_3p` 或 `remove_modification_internal`。
- **API：** `dnakit.genomics.remove_strand_modification(design: 'Mapping', strand_index: 'int', *, location: 'str', index: 'int | None' = None) -> 'dict[str, Any]'`。
- **输入：** 设计映射、`strand_index`、`location`，内部修饰需 `index`。
- **输出：** 新的 `design`。
- **示例代码：**

```python
from dnakit.genomics.scadnano_p26_functions import (
    assign_dna_sequence,
    build_strand,
    remove_strand_modification,
    set_strand_modification,
)

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]
design = assign_dna_sequence(design, "ACGT" * 12, strand_index=0)["design"]

design = set_strand_modification(design, 1, location="3p", display_text="Q", vendor_code="/3IABkFQ/")["design"]
print("3prime_modification" in design["strands"][1])
r = remove_strand_modification(design, 1, location="3p")
print("3prime_modification" in r["design"]["strands"][1])
```

- **示例输出：**

```text
True
False
```

## 20) `dnakit.genomics.list_design_modifications()`

- **作用：** 列出设计中用到的不同 DNA 修饰。
- **处理方法：** 调用 `Design.modifications`，可限定某一类型。
- **API：** `dnakit.genomics.list_design_modifications(design: 'Mapping', *, location: 'str | None' = None) -> 'dict[str, Any]'`。
- **输入：** 设计映射与可选 `location`。
- **输出：** `modifications`（映射列表）与 `count`。
- **示例代码：**

```python
from dnakit.genomics.scadnano_p26_functions import (
    assign_dna_sequence,
    build_strand,
    list_design_modifications,
    set_strand_modification,
)

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]
design = assign_dna_sequence(design, "ACGT" * 12, strand_index=0)["design"]

design = set_strand_modification(design, 1, location="5p", display_text="B", vendor_code="/5Biosg/")["design"]
r = list_design_modifications(design)
print(r["count"], r["modifications"])
print(list_design_modifications(design, location="3p")["count"])
```

- **示例输出：**

```text
1 [{'display_text': 'B', 'vendor_code': '/5Biosg/', 'location': "5'"}]
0
```

## 21) `dnakit.genomics.compute_base_pairs()`

- **作用：** 计算设计中的碱基对。
- **处理方法：** 调用 `Design.base_pairs`（Watson-Crick 配对；`allow_mismatches=True` 时含错配）。
- **API：** `dnakit.genomics.compute_base_pairs(design: 'Mapping', *, allow_mismatches: 'bool' = False) -> 'dict[str, Any]'`。
- **输入：** 已指定序列的设计映射、`allow_mismatches`。
- **输出：** `base_pairs`（螺旋索引到偏移列表）与 `total`。
- **示例代码：**

```python
from dnakit.genomics.scadnano_p26_functions import (
    assign_dna_sequence,
    build_strand,
    compute_base_pairs,
)

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]
design = assign_dna_sequence(design, "ACGT" * 12, strand_index=0)["design"]

r = compute_base_pairs(design)
print(r["total"], r["base_pairs"]["1"][:5])
```

- **示例输出：**

```text
24 [0, 1, 2, 3, 4]
```
