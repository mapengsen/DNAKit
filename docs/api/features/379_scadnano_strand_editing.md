# scadnano 设计：链与结构域编辑

添加或删除螺旋和链、逐个结构域构建链、编辑结构域列表、添加交叉、缺失、插入、切口与连接，以及整体平移设计的函数。

这些函数是对 [scadnano](https://github.com/UC-Davis-molecular-computing/scadnano-python-package)（MIT 许可，测试版本 0.21.1）的封装，scadnano 是设计 DNA 纳米结构的 Python 包。它是可选依赖，只在函数运行时导入；用 `pip install "dnakit[dna-nanotechnology]"` 安装（Excel 订单还需要 openpyxl，板图表格需要 tabulate，二者由 scadnano 一并安装）。设计是 scadnano 的 JSON 映射（即 `.sc` 文件内容）。编辑类函数不会修改传入的参数，而是在 `design` 键下返回新设计。链用其在 `design["strands"]` 中的 0 起始位置定位，结构域用其在链内的位置定位；螺旋索引与偏移沿用 scadnano 的定义。这些函数不预测折叠、产率或稳定性。若序列是在插入、缺失或结构域改动之前指定的，新设计仍保留旧序列，因此应最后再指定序列。

## 1) `dnakit.genomics.add_design_helix()`

- **作用：** 给设计添加一个螺旋。
- **处理方法：** 调用 `Design.add_helix`，并设置螺旋索引使保存的 JSON 可重新读取（scadnano 会把它留作 -1）。
- **API：** `dnakit.genomics.add_design_helix(design: 'Mapping', *, idx: 'int', max_offset: 'int', min_offset: 'int' = 0, grid_position: 'Sequence[int] | None' = None, roll: 'float' = 0.0, group: 'str' = 'default_group') -> 'dict[str, Any]'`。
- **输入：** 设计映射、`idx`、`max_offset`、`min_offset`、`grid_position`、`roll`、`group`。
- **输出：** 新的 `design` 与计数。
- **示例代码：**

```python
from dnakit.genomics.scadnano_p26_functions import (
    add_design_helix,
    assign_dna_sequence,
    build_strand,
)

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]
design = assign_dna_sequence(design, "ACGT" * 12, strand_index=0)["design"]

r = add_design_helix(design, idx=2, max_offset=64, grid_position=[0, 2])
print(r["helices"], r["design"]["helices"][2])
```

- **示例输出：**

```text
3 {'max_offset': 64, 'grid_position': [0, 2]}
```

## 2) `dnakit.genomics.add_design_strand()`

- **作用：** 把以 scadnano JSON 给出的链添加到设计。
- **处理方法：** 把映射追加到链列表并重新载入设计，由此校验（`Design.add_strand`）。
- **API：** `dnakit.genomics.add_design_strand(design: 'Mapping', strand: 'Mapping', *) -> 'dict[str, Any]'`。
- **输入：** 设计映射与含 `domains` 列表的链映射。
- **输出：** 新的 `design` 与 `added_strand_index`。
- **示例代码：**

```python
from dnakit.genomics.scadnano_p26_functions import add_design_strand, build_strand

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]

r = add_design_strand(design, {"domains": [{"helix": 0, "forward": True, "start": 30, "end": 40}], "name": "extra"})
print(r["added_strand_index"], r["strand_names"])
```

- **示例输出：**

```text
2 [None, 'stapA', 'extra']
```

## 3) `dnakit.genomics.build_strand()`

- **作用：** 逐个结构域构建链（交叉、环、延伸、缺失、插入、序列、支架、环状、vendor 字段）。
- **处理方法：** 在 `Design.draw_strand` 得到的 scadnano `StrandBuilder` 上依次执行各步骤；只允许列出的操作，未知键会被拒绝。
- **API：** `dnakit.genomics.build_strand(design: 'Mapping', *, helix: 'int', offset: 'int', steps: 'Sequence[Mapping]') -> 'dict[str, Any]'`。
- **输入：** 设计映射、起点 `helix` 与 `offset`，以及 `steps`：含 `op` 的映射列表（`to`、`move`、`cross`、`loopout`、`extension_3p`、`extension_5p`、`with_deletions`、`with_insertions`、`with_sequence`、`with_domain_sequence`、`with_name`、`with_label`、`with_domain_name`、`as_scaffold`、`as_circular`、`with_vendor_fields`）。
- **输出：** 新的 `design`、`added_strand_index` 与 `strand` 摘要。
- **示例代码：**

```python
from dnakit.genomics.scadnano_p26_functions import build_strand

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]

steps = [{"op": "to", "offset": 40}, {"op": "loopout", "helix": 1, "length": 3}, {"op": "move", "delta": -16}, {"op": "with_deletions", "deletions": [30]}, {"op": "with_name", "name": "hairpin"}]
r = build_strand(design, helix=0, offset=24, steps=steps)
print(r["added_strand_index"], r["strand"]["name"], r["strand"]["domains"], r["strand"]["dna_length"])
```

- **示例输出：**

```text
2 hairpin 3 34
```

## 4) `dnakit.genomics.remove_design_strand()`

- **作用：** 从设计中删除一条链。
- **处理方法：** 调用 `Design.remove_strand`。
- **API：** `dnakit.genomics.remove_design_strand(design: 'Mapping', strand_index: 'int', *) -> 'dict[str, Any]'`。
- **输入：** 设计映射与 `strand_index`。
- **输出：** 新的 `design`。
- **示例代码：**

```python
from dnakit.genomics.scadnano_p26_functions import (
    assign_dna_sequence,
    build_strand,
    remove_design_strand,
)

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]
design = assign_dna_sequence(design, "ACGT" * 12, strand_index=0)["design"]

r = remove_design_strand(design, 1)
print(r["strands"], r["removed_strand_index"])
```

- **示例输出：**

```text
1 1
```

## 5) `dnakit.genomics.find_strand_by_name()`

- **作用：** 按名称查找链。
- **处理方法：** 像 `Design.strand_with_name` 那样查找；名称不存在时返回 `found=False` 而不是报错。
- **API：** `dnakit.genomics.find_strand_by_name(design: 'Mapping', name: 'str', *) -> 'dict[str, Any]'`。
- **输入：** 设计映射与 `name`。
- **输出：** `found` 与链摘要。
- **示例代码：**

```python
from dnakit.genomics.scadnano_p26_functions import (
    assign_dna_sequence,
    build_strand,
    find_strand_by_name,
)

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]
design = assign_dna_sequence(design, "ACGT" * 12, strand_index=0)["design"]

print(find_strand_by_name(design, "stapA")["strand_index"])
print(find_strand_by_name(design, "nope")["found"])
```

- **示例输出：**

```text
1
False
```

## 6) `dnakit.genomics.list_strands_at_helix_end()`

- **作用：** 列出 5' 端或 3' 端位于某螺旋上的链。
- **处理方法：** 调用 `Design.strands_starting_on_helix`（`end="5p"`）或 `strands_ending_on_helix`（`end="3p"`）。
- **API：** `dnakit.genomics.list_strands_at_helix_end(design: 'Mapping', helix: 'int', *, end: 'str' = '5p') -> 'dict[str, Any]'`。
- **输入：** 设计映射、`helix`、`end`。
- **输出：** `strands`：含 `strand_index` 的链摘要。
- **示例代码：**

```python
from dnakit.genomics.scadnano_p26_functions import (
    assign_dna_sequence,
    build_strand,
    list_strands_at_helix_end,
)

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]
design = assign_dna_sequence(design, "ACGT" * 12, strand_index=0)["design"]

print([s["strand_index"] for s in list_strands_at_helix_end(design, 1, end="5p")["strands"]])
print([s["strand_index"] for s in list_strands_at_helix_end(design, 1, end="3p")["strands"]])
```

- **示例输出：**

```text
[1]
[0, 1]
```

## 7) `dnakit.genomics.find_domains_at_offset()`

- **作用：** 查找覆盖某个螺旋偏移的结构域。
- **处理方法：** 调用 `Design.domains_at`；给出 `forward` 时调用 `Design.domain_at`。
- **API：** `dnakit.genomics.find_domains_at_offset(design: 'Mapping', helix: 'int', offset: 'int', *, forward: 'bool | None' = None) -> 'dict[str, Any]'`。
- **输入：** 设计映射、`helix`、`offset`、可选 `forward`。
- **输出：** `domains`：含链与结构域索引、螺旋、方向、起止的行。
- **示例代码：**

```python
from dnakit.genomics.scadnano_p26_functions import (
    assign_dna_sequence,
    build_strand,
    find_domains_at_offset,
)

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]
design = assign_dna_sequence(design, "ACGT" * 12, strand_index=0)["design"]

rows = find_domains_at_offset(design, 1, 10)["domains"]
print([(d["strand_index"], d["domain_index"], d["forward"]) for d in rows])
print(find_domains_at_offset(design, 1, 10, forward=False)["domains"][0]["strand_index"])
```

- **示例输出：**

```text
[(0, 1, False), (1, 0, True)]
0
```

## 8) `dnakit.genomics.insert_strand_domain()`

- **作用：** 向链中插入（或追加）结构域、环或延伸。
- **处理方法：** 调用 `Design.insert_domain`（省略 `order` 时调用 `Design.append_domain`）。
- **API：** `dnakit.genomics.insert_strand_domain(design: 'Mapping', strand_index: 'int', domain: 'Mapping', *, order: 'int | None' = None) -> 'dict[str, Any]'`。
- **输入：** 设计映射、`strand_index`、结构域映射、可选 `order`。
- **输出：** 新的 `design` 与 `assigned` 摘要。
- **示例代码：**

```python
from dnakit.genomics.scadnano_p26_functions import (
    assign_dna_sequence,
    build_strand,
    insert_strand_domain,
)

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]
design = assign_dna_sequence(design, "ACGT" * 12, strand_index=0)["design"]

r = insert_strand_domain(design, 1, {"helix": 1, "forward": True, "start": 30, "end": 38})
print(r["assigned"]["domains"], r["design"]["strands"][1]["domains"][-1])
r = insert_strand_domain(design, 1, {"loopout": 4}, order=0)
print(r["design"]["strands"][1]["domains"][0])
```

- **示例输出：**

```text
2 {'helix': 1, 'forward': True, 'start': 30, 'end': 38}
{'loopout': 4}
```

## 9) `dnakit.genomics.remove_strand_domain()`

- **作用：** 从链中删除一个结构域。
- **处理方法：** 调用 `Design.remove_domain`。
- **API：** `dnakit.genomics.remove_strand_domain(design: 'Mapping', strand_index: 'int', domain_index: 'int', *) -> 'dict[str, Any]'`。
- **输入：** 设计映射、`strand_index`、`domain_index`。
- **输出：** 新的 `design` 与 `assigned` 摘要。
- **示例代码：**

```python
from dnakit.genomics.scadnano_p26_functions import (
    create_origami_rectangle,
    remove_strand_domain,
)

design = create_origami_rectangle(num_helices=4, num_cols=10)["design"]

r = remove_strand_domain(design, 3, 0)
print(design["strands"][3]["domains"])
print(r["design"]["strands"][3]["domains"])
```

- **示例输出：**

```text
[{'helix': 0, 'forward': False, 'start': 160, 'end': 176}, {'helix': 1, 'forward': True, 'start': 160, 'end': 176}]
[{'helix': 1, 'forward': True, 'start': 160, 'end': 176}]
```

## 10) `dnakit.genomics.set_strand_domains()`

- **作用：** 替换一条链的完整结构域列表。
- **处理方法：** 用映射构造结构域后调用 `Strand.set_domains`。
- **API：** `dnakit.genomics.set_strand_domains(design: 'Mapping', strand_index: 'int', domains: 'Sequence[Mapping]', *) -> 'dict[str, Any]'`。
- **输入：** 设计映射、`strand_index`、结构域映射列表。
- **输出：** 新的 `design` 与 `assigned` 摘要。
- **示例代码：**

```python
from dnakit.genomics.scadnano_p26_functions import (
    create_origami_rectangle,
    set_strand_domains,
)

design = create_origami_rectangle(num_helices=4, num_cols=10)["design"]

r = set_strand_domains(design, 3, [{"helix": 3, "forward": True, "start": 16, "end": 24}])
print(r["assigned"]["domains"], r["assigned"]["dna_length"])
```

- **示例输出：**

```text
1 8
```

## 11) `dnakit.genomics.rotate_strand_domains()`

- **作用：** 旋转链中结构域的顺序。
- **处理方法：** 调用 `Strand.rotate_domains`。
- **API：** `dnakit.genomics.rotate_strand_domains(design: 'Mapping', strand_index: 'int', rotation: 'int', *, forward: 'bool' = True) -> 'dict[str, Any]'`。
- **输入：** 设计映射、`strand_index`、`rotation`、`forward`。
- **输出：** 新的 `design` 与 `assigned` 摘要。
- **示例代码：**

```python
from dnakit.genomics.scadnano_p26_functions import (
    create_origami_rectangle,
    rotate_strand_domains,
)

design = create_origami_rectangle(num_helices=4, num_cols=10)["design"]

r = rotate_strand_domains(design, 3, 1)
print([d["helix"] for d in design["strands"][3]["domains"]], [d["helix"] for d in r["design"]["strands"][3]["domains"]])
```

- **示例输出：**

```text
[0, 1] [1, 0]
```

## 12) `dnakit.genomics.reverse_strands()`

- **作用：** 反转一条链或全部链的极性。
- **处理方法：** 对单链调用 `Strand.reverse`；`strand_index` 为 `None` 时调用 `Design.reverse_all`。
- **API：** `dnakit.genomics.reverse_strands(design: 'Mapping', *, strand_index: 'int | None' = None) -> 'dict[str, Any]'`。
- **输入：** 设计映射与可选 `strand_index`。
- **输出：** 新的 `design`。
- **示例代码：**

```python
from dnakit.genomics.scadnano_p26_functions import (
    create_origami_rectangle,
    reverse_strands,
)

design = create_origami_rectangle(num_helices=4, num_cols=10)["design"]

r = reverse_strands(design, strand_index=3)
print(design["strands"][3]["domains"][0]["helix"], r["design"]["strands"][3]["domains"][0]["helix"])
print(reverse_strands(design)["strands"])
```

- **示例输出：**

```text
0 1
21
```

## 13) `dnakit.genomics.set_strand_properties()`

- **作用：** 把链标为支架、设为环状或线状，或设置 vendor 字段。
- **处理方法：** 调用 `Strand.set_scaffold`、`set_circular` / `set_linear` 并设置 `VendorFields`；为 `None` 的参数保持不变。
- **API：** `dnakit.genomics.set_strand_properties(design: 'Mapping', strand_index: 'int', *, scaffold: 'bool | None' = None, circular: 'bool | None' = None, vendor_fields: 'Mapping[str, str | None] | None' = None) -> 'dict[str, Any]'`。
- **输入：** 设计映射、`strand_index`、`scaffold`、`circular`、`vendor_fields`。
- **输出：** 新的 `design` 与 `assigned` 摘要。
- **示例代码：**

```python
from dnakit.genomics.scadnano_p26_functions import (
    assign_dna_sequence,
    build_strand,
    set_strand_properties,
)

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]
design = assign_dna_sequence(design, "ACGT" * 12, strand_index=0)["design"]

r = set_strand_properties(design, 1, scaffold=True, vendor_fields={"scale": "100nm", "purification": "HPLC"})
print(r["assigned"]["scaffold"], r["design"]["strands"][1]["vendor"])
r = set_strand_properties(design, 0, circular=False)
print(r["assigned"]["circular"])
```

- **示例输出：**

```text
True {'scale': '100nm', 'purification': 'HPLC'}
False
```

## 14) `dnakit.genomics.add_crossover()`

- **作用：** 在两个螺旋间添加半交叉或全交叉。
- **处理方法：** 调用 `Design.add_half_crossover`；`full=True` 时调用 `Design.add_full_crossover`（两个半交叉）。
- **API：** `dnakit.genomics.add_crossover(design: 'Mapping', *, helix: 'int', helix2: 'int', offset: 'int', forward: 'bool', offset2: 'int | None' = None, forward2: 'bool | None' = None, full: 'bool' = False) -> 'dict[str, Any]'`。
- **输入：** 设计映射、`helix`、`helix2`、`offset`、`forward`，可选 `offset2`、`forward2`、`full`。
- **输出：** 新的 `design`。
- **示例代码：**

```python
from dnakit.genomics.scadnano_p26_functions import add_crossover, build_strand

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
design = build_strand(empty, helix=0, offset=0, steps=[{"op": "to", "offset": 24}])["design"]
design = build_strand(design, helix=1, offset=24, steps=[{"op": "move", "delta": -24}])["design"]

r = add_crossover(design, helix=0, helix2=1, offset=23, forward=True)
print(len(design["strands"]), "->", r["strands"], r["design"]["strands"][0]["domains"])
r = add_crossover(design, helix=0, helix2=1, offset=10, forward=True, full=True)
print(r["strands"])
```

- **示例输出：**

```text
2 -> 1 [{'helix': 0, 'forward': True, 'start': 0, 'end': 24}, {'helix': 1, 'forward': False, 'start': 0, 'end': 24}]
2
```

## 15) `dnakit.genomics.add_design_deletion()`

- **作用：** 在螺旋偏移处添加单碱基缺失。
- **处理方法：** 调用 `Design.add_deletion`（构建器的 `with_deletions` 也是同一操作）。
- **API：** `dnakit.genomics.add_design_deletion(design: 'Mapping', helix: 'int', offset: 'int', *) -> 'dict[str, Any]'`。
- **输入：** 设计映射、`helix`、`offset`。
- **输出：** 新的 `design`。
- **示例代码：**

```python
from dnakit.genomics.scadnano_p26_functions import add_design_deletion, build_strand

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]

r = add_design_deletion(design, 0, 10)
print(r["strand_lengths"], r["design"]["strands"][0]["domains"][0])
```

- **示例输出：**

```text
[47, 24] {'helix': 0, 'forward': True, 'start': 0, 'end': 24, 'deletions': [10]}
```

## 16) `dnakit.genomics.add_design_insertion()`

- **作用：** 在螺旋偏移处添加额外碱基插入。
- **处理方法：** 调用 `Design.add_insertion`（构建器的 `with_insertions` 也是同一操作）。
- **API：** `dnakit.genomics.add_design_insertion(design: 'Mapping', helix: 'int', offset: 'int', length: 'int', *) -> 'dict[str, Any]'`。
- **输入：** 设计映射、`helix`、`offset`、`length`。
- **输出：** 新的 `design`。
- **示例代码：**

```python
from dnakit.genomics.scadnano_p26_functions import add_design_insertion, build_strand

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]

r = add_design_insertion(design, 0, 10, 2)
print(r["strand_lengths"], r["design"]["strands"][0]["domains"][0])
```

- **示例输出：**

```text
[50, 24] {'helix': 0, 'forward': True, 'start': 0, 'end': 24, 'insertions': [[10, 2]]}
```

## 17) `dnakit.genomics.add_design_nick()`

- **作用：** 在结构域上切口，把链一分为二。
- **处理方法：** 调用 `Design.add_nick`；切口位于 `offset - 1` 与 `offset` 之间。
- **API：** `dnakit.genomics.add_design_nick(design: 'Mapping', helix: 'int', offset: 'int', forward: 'bool', *, new_color: 'bool' = True) -> 'dict[str, Any]'`。
- **输入：** 设计映射、`helix`、`offset`、`forward`、`new_color`。
- **输出：** 新的 `design`。
- **示例代码：**

```python
from dnakit.genomics.scadnano_p26_functions import add_design_nick, build_strand

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]

r = add_design_nick(design, 0, 12, True)
print(r["strands"], r["strand_lengths"])
```

- **示例输出：**

```text
3 [24, 12, 36]
```

## 18) `dnakit.genomics.ligate_design_nick()`

- **作用：** 连接在切口处相遇的两条链。
- **处理方法：** 调用 `Design.ligate`。
- **API：** `dnakit.genomics.ligate_design_nick(design: 'Mapping', helix: 'int', offset: 'int', forward: 'bool', *) -> 'dict[str, Any]'`。
- **输入：** 设计映射、`helix`、`offset`、`forward`。
- **输出：** 新的 `design`。
- **示例代码：**

```python
from dnakit.genomics.scadnano_p26_functions import (
    add_design_nick,
    build_strand,
    ligate_design_nick,
)

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]

nicked = add_design_nick(design, 0, 12, True)["design"]
r = ligate_design_nick(nicked, 0, 12, True)
print(r["strands"], r["strand_lengths"])
```

- **示例输出：**

```text
2 [24, 48]
```

## 19) `dnakit.genomics.move_strands()`

- **作用：** 沿螺旋或跨螺旋平移所有链。
- **处理方法：** `axis="offset"` 调用 `Design.move_strand_offsets`；`axis="helix"` 调用 `Design.move_strands_on_helices`。
- **API：** `dnakit.genomics.move_strands(design: 'Mapping', delta: 'int', *, axis: 'str' = 'offset') -> 'dict[str, Any]'`。
- **输入：** 设计映射、`delta`、`axis`。
- **输出：** 新的 `design`。
- **示例代码：**

```python
from dnakit.genomics.scadnano_p26_functions import (
    assign_dna_sequence,
    build_strand,
    move_strands,
)

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]
design = assign_dna_sequence(design, "ACGT" * 12, strand_index=0)["design"]

r = move_strands(design, 5)
print(design["strands"][0]["domains"][0], r["design"]["strands"][0]["domains"][0])
r = move_strands(design, 0, axis="helix")
print(r["strands"])
```

- **示例输出：**

```text
{'helix': 0, 'forward': True, 'start': 0, 'end': 24} {'helix': 0, 'forward': True, 'start': 5, 'end': 29}
2
```
