# scadnano 设计：重叠、交叉、螺旋几何与 DNA 折纸

检测结构域重叠、列出交叉与潜在交叉位点、运行 autostaple 与滚转松弛、报告骨架角度，以及生成带扭转修正缺失的 DNA 折纸矩形的函数。

这些函数是对 [scadnano](https://github.com/UC-Davis-molecular-computing/scadnano-python-package)（MIT 许可，测试版本 0.21.1）的封装，scadnano 是设计 DNA 纳米结构的 Python 包。它是可选依赖，只在函数运行时导入；用 `pip install "dnakit[dna-nanotechnology]"` 安装（Excel 订单还需要 openpyxl，板图表格需要 tabulate，二者由 scadnano 一并安装）。设计是 scadnano 的 JSON 映射（即 `.sc` 文件内容）。编辑类函数不会修改传入的参数，而是在 `design` 键下返回新设计。链用其在 `design["strands"]` 中的 0 起始位置定位，结构域用其在链内的位置定位；螺旋索引与偏移沿用 scadnano 的定义。这些函数不预测折叠、产率或稳定性。若序列是在插入、缺失或结构域改动之前指定的，新设计仍保留旧序列，因此应最后再指定序列。

## 1) `dnakit.genomics.check_domain_overlap()`

- **作用：** 判断两个结构域或两条链在螺旋上是否重叠，并给出重叠区间。
- **处理方法：** 调用 `Domain.overlaps`（方向相反）、`Domain.compute_overlap` 与 `Domain.overlaps_illegally`（方向相同），整链则调用 `Strand.overlaps`。
- **API：** `dnakit.genomics.check_domain_overlap(design: 'Mapping', strand_index: 'int', other_strand_index: 'int', *, domain_index: 'int | None' = None, other_domain_index: 'int | None' = None) -> 'dict[str, Any]'`。
- **输入：** 设计映射、两个链索引，可选两个结构域索引。
- **输出：** `overlaps`、`overlap_interval`（含起点、不含终点）、`overlaps_illegally`。
- **示例代码：**

```python
from dnakit.genomics.scadnano_p26_functions import (
    assign_dna_sequence,
    build_strand,
    check_domain_overlap,
)

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]
design = assign_dna_sequence(design, "ACGT" * 12, strand_index=0)["design"]

r = check_domain_overlap(design, 0, 0, domain_index=0, other_domain_index=0)
print(r["overlap_interval"], r["overlaps"], r["overlaps_illegally"])
print(check_domain_overlap(design, 0, 1)["overlaps"])
```

- **示例输出：**

```text
[0, 24] False True
True
```

## 2) `dnakit.genomics.find_overlapping_domains()`

- **作用：** 查找偏移范围与给定结构域重叠的结构域。
- **处理方法：** 对设计中其余结构域调用 `Domain.find_overlapping_ranges`。它只比较范围，不看螺旋与方向。
- **API：** `dnakit.genomics.find_overlapping_domains(design: 'Mapping', strand_index: 'int', domain_index: 'int', *) -> 'dict[str, Any]'`。
- **输入：** 设计映射、`strand_index`、`domain_index`。
- **输出：** `domains`：带位置的行。
- **示例代码：**

```python
from dnakit.genomics.scadnano_p26_functions import (
    assign_dna_sequence,
    build_strand,
    find_overlapping_domains,
)

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]
design = assign_dna_sequence(design, "ACGT" * 12, strand_index=0)["design"]

r = find_overlapping_domains(design, 0, 0)
print([(d["strand_index"], d["domain_index"]) for d in r["domains"]])
```

- **示例输出：**

```text
[(0, 1), (1, 0)]
```

## 3) `dnakit.genomics.find_helix_domain_overlap_pairs()`

- **作用：** 列出在螺旋上重叠的正向/反向结构域对。
- **处理方法：** 对一个或全部螺旋调用 `scadnano.find_overlapping_domains_on_helix`。
- **API：** `dnakit.genomics.find_helix_domain_overlap_pairs(design: 'Mapping', *, helix: 'int | None' = None) -> 'dict[str, Any]'`。
- **输入：** 设计映射与可选 `helix`。
- **输出：** `pairs` 与 `count`。
- **示例代码：**

```python
from dnakit.genomics.scadnano_p26_functions import (
    assign_dna_sequence,
    build_strand,
    find_helix_domain_overlap_pairs,
)

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]
design = assign_dna_sequence(design, "ACGT" * 12, strand_index=0)["design"]

r = find_helix_domain_overlap_pairs(design)
print(r["count"], r["pairs"])
```

- **示例输出：**

```text
1 [{'helix': 1, 'first': {'strand_index': 1, 'domain_index': 0}, 'second': {'strand_index': 0, 'domain_index': 1}}]
```

## 4) `dnakit.genomics.domain_has_crossover()`

- **作用：** 判断结构域在某偏移处是否有交叉。
- **处理方法：** 调用 `Domain.has_crossover_at`。
- **API：** `dnakit.genomics.domain_has_crossover(design: 'Mapping', strand_index: 'int', domain_index: 'int', offset: 'int', *) -> 'dict[str, Any]'`。
- **输入：** 设计映射、`strand_index`、`domain_index`、`offset`。
- **输出：** `has_crossover`。
- **示例代码：**

```python
from dnakit.genomics.scadnano_p26_functions import (
    assign_dna_sequence,
    build_strand,
    domain_has_crossover,
)

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]
design = assign_dna_sequence(design, "ACGT" * 12, strand_index=0)["design"]

print(domain_has_crossover(design, 0, 0, 23)["has_crossover"])
print(domain_has_crossover(design, 0, 0, 10)["has_crossover"])
```

- **示例输出：**

```text
True
False
```

## 5) `dnakit.genomics.map_domain_offset_to_dna_index()`

- **作用：** 把结构域位置映射为链 DNA 序列中的下标。
- **处理方法：** 调用 `Strand.dna_index_start_domain`（等同于 `Domain.get_seq_start_idx` 与 `Loopout.get_seq_start_idx`），给出 `offset` 时再调用 `Domain.domain_offset_to_strand_dna_idx`。
- **API：** `dnakit.genomics.map_domain_offset_to_dna_index(design: 'Mapping', strand_index: 'int', domain_index: 'int', *, offset: 'int | None' = None, offset_closer_to_5p: 'bool' = True) -> 'dict[str, Any]'`。
- **输入：** 设计映射、`strand_index`、`domain_index`、可选 `offset`、`offset_closer_to_5p`。
- **输出：** `start_index`，给出 `offset` 时还有 `dna_index`。
- **示例代码：**

```python
from dnakit.genomics.scadnano_p26_functions import (
    assign_dna_sequence,
    build_strand,
    map_domain_offset_to_dna_index,
)

empty = {"grid": "square", "helices": [{"max_offset": 48}, {"max_offset": 48}], "strands": []}
scaffold = [{"op": "to", "offset": 24}, {"op": "cross", "helix": 1}, {"op": "move", "delta": -24}, {"op": "as_scaffold"}]
staple = [{"op": "to", "offset": 24}, {"op": "with_name", "name": "stapA"}, {"op": "with_vendor_fields", "plate": "P1", "well": "A1"}]
design = build_strand(empty, helix=0, offset=0, steps=scaffold)["design"]
design = build_strand(design, helix=1, offset=0, steps=staple)["design"]
design = assign_dna_sequence(design, "ACGT" * 12, strand_index=0)["design"]

r = map_domain_offset_to_dna_index(design, 0, 1, offset=10)
print(r["start_index"], r["dna_index"])
```

- **示例输出：**

```text
24 37
```

## 6) `dnakit.genomics.list_helix_crossovers()`

- **作用：** 列出与某螺旋相连的交叉。
- **处理方法：** 调用 `Helix.crossover_addresses`。
- **API：** `dnakit.genomics.list_helix_crossovers(design: 'Mapping', helix: 'int', *, allow_intrahelix: 'bool' = True, allow_intergroup: 'bool' = True) -> 'dict[str, Any]'`。
- **输入：** 设计映射、`helix`、`allow_intrahelix`、`allow_intergroup`。
- **输出：** `crossovers`：对方螺旋、偏移、方向；`count`。
- **示例代码：**

```python
from dnakit.genomics.scadnano_p26_functions import (
    create_origami_rectangle,
    list_helix_crossovers,
)

design = create_origami_rectangle(num_helices=4, num_cols=10)["design"]

r = list_helix_crossovers(design, 0)
print(r["count"], r["crossovers"][:3])
```

- **示例输出：**

```text
10 [{'helix': 1, 'offset': 16, 'forward': True}, {'helix': 1, 'offset': 175, 'forward': True}, {'helix': 1, 'offset': 31, 'forward': False}]
```

## 7) `dnakit.genomics.split_helix_strands()`

- **作用：** 把螺旋上的结构域分为订书钉与支架结构域。
- **处理方法：** 调用 `Helix.get_strand_sets`，即 cadnano 的奇偶规则。
- **API：** `dnakit.genomics.split_helix_strands(design: 'Mapping', helix: 'int', *) -> 'dict[str, Any]'`。
- **输入：** 设计映射与 `helix`。
- **输出：** `staple_domains` 与 `scaffold_domains`。
- **示例代码：**

```python
from dnakit.genomics.scadnano_p26_functions import (
    create_origami_rectangle,
    split_helix_strands,
)

design = create_origami_rectangle(num_helices=4, num_cols=10)["design"]

r = split_helix_strands(design, 0)
print(len(r["staple_domains"]), len(r["scaffold_domains"]))
```

- **示例输出：**

```text
9 1
```

## 8) `dnakit.genomics.list_crossover_points()`

- **作用：** 返回设计网格对应的 cadnano 交叉点表。
- **处理方法：** 调用 `Design.get_crossover_points`（方形与蜂窝网格）。
- **API：** `dnakit.genomics.list_crossover_points(design: 'Mapping', *) -> 'dict[str, Any]'`。
- **输入：** 设计映射。
- **输出：** `scaffold_low`、`scaffold_high`、`staple_low`、`staple_high`。
- **示例代码：**

```python
from dnakit.genomics.scadnano_p26_functions import (
    create_origami_rectangle,
    list_crossover_points,
)

design = create_origami_rectangle(num_helices=4, num_cols=10)["design"]

r = list_crossover_points(design)
print(r["grid"], r["scaffold_low"][0], r["staple_low"])
```

- **示例输出：**

```text
square [4, 26, 15] [[31], [23], [15], [7]]
```

## 9) `dnakit.genomics.list_potential_crossovers()`

- **作用：** 列出可与相邻螺旋形成交叉的偏移。
- **处理方法：** 调用 `Design.potential_crossover_list`。
- **API：** `dnakit.genomics.list_potential_crossovers(design: 'Mapping', helix: 'int', *, index: 'int | None' = None) -> 'dict[str, Any]'`。
- **输入：** 设计映射、`helix`、可选 `index`。
- **输出：** `sites`：相邻螺旋、偏移、`scaffold` 或 `staple`、`is_low`。
- **示例代码：**

```python
from dnakit.genomics.scadnano_p26_functions import (
    create_origami_rectangle,
    list_potential_crossovers,
)

design = create_origami_rectangle(num_helices=4, num_cols=10)["design"]

r = list_potential_crossovers(design, 0)
print(r["count"], r["sites"][:2])
```

- **示例输出：**

```text
46 [{'neighbor_helix': 1, 'offset': 18, 'strand_type': 'scaffold', 'is_low': True}, {'neighbor_helix': 1, 'offset': 28, 'strand_type': 'scaffold', 'is_low': True}]
```

## 10) `dnakit.genomics.get_helix_neighbors()`

- **作用：** 返回某螺旋在晶格上的相邻螺旋。
- **处理方法：** 调用 `Design.get_helix_neighbors`；`None` 表示晶格上无螺旋的位置。
- **API：** `dnakit.genomics.get_helix_neighbors(design: 'Mapping', helix: 'int', *) -> 'dict[str, Any]'`。
- **输入：** 设计映射与 `helix`。
- **输出：** `neighbors`。
- **示例代码：**

```python
from dnakit.genomics.scadnano_p26_functions import (
    create_origami_rectangle,
    get_helix_neighbors,
)

design = create_origami_rectangle(num_helices=4, num_cols=10)["design"]

print(get_helix_neighbors(design, 0)["neighbors"])
print(get_helix_neighbors(design, 1)["neighbors"])
```

- **示例输出：**

```text
[None, 1, None, None]
[None, 0, None, 2]
```

## 11) `dnakit.genomics.autostaple_design()`

- **作用：** 用 cadnano 的 autostaple 算法在支架上生成订书钉。
- **处理方法：** 调用 `Design.autostaple`，它会先移除旧订书钉。
- **API：** `dnakit.genomics.autostaple_design(design: 'Mapping', *) -> 'dict[str, Any]'`。
- **输入：** 不含环或延伸的设计映射。
- **输出：** 新的 `design`。
- **示例代码：**

```python
from dnakit.genomics.scadnano_p26_functions import (
    autostaple_design,
    create_origami_rectangle,
)

design = create_origami_rectangle(num_helices=4, num_cols=10)["design"]

r = autostaple_design(design)
print(len(design["strands"]), "->", r["strands"], r["scaffolds"])
```

- **示例输出：**

```text
21 -> 9 1
```

## 12) `dnakit.genomics.relax_helix_rolls()`

- **作用：** 松弛螺旋的滚转角，或仅计算滚转修正量。
- **处理方法：** `compute_only=False` 调用 `Design.relax_helix_rolls`；`True` 只报告每个螺旋的 `Helix.compute_relaxed_roll_delta`，设计保持不变。
- **API：** `dnakit.genomics.relax_helix_rolls(design: 'Mapping', *, compute_only: 'bool' = False) -> 'dict[str, Any]'`。
- **输入：** 设计映射、`compute_only`。
- **输出：** `applied`，以及 `roll_deltas`，或新的 `design` 与 `rolls`。
- **示例代码：**

```python
from dnakit.genomics.scadnano_p26_functions import (
    create_origami_rectangle,
    relax_helix_rolls,
)

design = create_origami_rectangle(num_helices=4, num_cols=10)["design"]

r = relax_helix_rolls(design, compute_only=True)
print(r["applied"], {k: round(v, 2) for k, v in r["roll_deltas"].items()})
r = relax_helix_rolls(design)
print(r["applied"], {k: round(v, 2) for k, v in r["rolls"].items()})
```

- **示例输出：**

```text
False {'0': 349.71, '1': 331.71, '2': 349.71, '3': 331.71}
True {'0': 349.71, '1': 331.71, '2': 349.71, '3': 331.71}
```

## 13) `dnakit.genomics.get_backbone_angle()`

- **作用：** 返回螺旋偏移处的骨架角度。
- **处理方法：** 按设计几何参数调用 `Helix.backbone_angle_at_offset`。
- **API：** `dnakit.genomics.get_backbone_angle(design: 'Mapping', helix: 'int', offset: 'int', forward: 'bool', *) -> 'dict[str, Any]'`。
- **输入：** 设计映射、`helix`、`offset`、`forward`。
- **输出：** `angle_degrees`。
- **示例代码：**

```python
from dnakit.genomics.scadnano_p26_functions import (
    create_origami_rectangle,
    get_backbone_angle,
)

design = create_origami_rectangle(num_helices=4, num_cols=10)["design"]

print(round(get_backbone_angle(design, 0, 10, True)["angle_degrees"], 2))
print(round(get_backbone_angle(design, 0, 10, False)["angle_degrees"], 2))
```

- **示例输出：**

```text
342.86
132.86
```

## 14) `dnakit.genomics.create_origami_rectangle()`

- **作用：** 生成 DNA 折纸矩形设计。
- **处理方法：** 调用 `origami_rectangle.create`，默认指定 M13mp18。
- **API：** `dnakit.genomics.create_origami_rectangle(*, num_helices: 'int', num_cols: 'int', assign_seq: 'bool' = True, seam_left_column: 'int' = -1, nick_pattern: 'str' = 'staggered', twist_correction_deletion_spacing: 'int' = 0, twist_correction_start_col: 'int' = 1, twist_correction_deletion_offset: 'int' = -1, num_flanking_columns: 'int' = 1, num_flanking_helices: 'int' = 0, custom_scaffold: 'str | None' = None, edge_staples: 'bool' = True, scaffold_nick_offset: 'int' = -1) -> 'dict[str, Any]'`。
- **输入：** `num_helices`、`num_cols` 及上游选项（`nick_pattern`、扭转修正间隔、侧翼列数、`custom_scaffold` 等）。
- **输出：** `design` 与其计数。
- **示例代码：**

```python
from dnakit.genomics.scadnano_p26_functions import create_origami_rectangle

r = create_origami_rectangle(num_helices=4, num_cols=10)
print(r["helices"], r["strands"], r["scaffolds"], r["strand_lengths"][:3])
r = create_origami_rectangle(num_helices=4, num_cols=10, twist_correction_deletion_spacing=3)
print(r["strand_lengths"][0])
```

- **示例输出：**

```text
4 21 1 [640, 32, 32]
628
```

## 15) `dnakit.genomics.add_twist_correction_deletions()`

- **作用：** 给折纸矩形添加扭转修正缺失。
- **处理方法：** 调用 `origami_rectangle.add_twist_correction_deletions`：每隔 `deletion_spacing` 列在每个螺旋上加一个缺失。
- **API：** `dnakit.genomics.add_twist_correction_deletions(design: 'Mapping', *, offset_start: 'int', deletion_spacing: 'int', deletion_start_col: 'int', deletion_offset: 'int', num_helices: 'int', num_cols: 'int', num_flanking_helices: 'int' = 0) -> 'dict[str, Any]'`。
- **输入：** 设计映射与矩形布局参数。
- **输出：** 新的 `design`。
- **示例代码：**

```python
from dnakit.genomics.scadnano_p26_functions import (
    add_twist_correction_deletions,
    create_origami_rectangle,
)

design = create_origami_rectangle(num_helices=4, num_cols=10, assign_seq=False)["design"]

r = add_twist_correction_deletions(design, offset_start=0, deletion_spacing=2, deletion_start_col=1, deletion_offset=-1, num_helices=4, num_cols=10)
print(sum(r["strand_lengths"]) - sum(create_origami_rectangle(num_helices=4, num_cols=10, assign_seq=False)["strand_lengths"]))
```

- **示例输出：**

```text
-40
```

## 16) `dnakit.genomics.add_deletion_in_range()`

- **作用：** 在螺旋某范围内的空闲偏移处添加一个缺失。
- **处理方法：** 调用 `origami_rectangle.add_deletion_in_range`；空闲偏移须被两个结构域覆盖且不是结构域端点。
- **API：** `dnakit.genomics.add_deletion_in_range(design: 'Mapping', *, helix: 'int', start: 'int', end: 'int', deletion_offset: 'int' = -1) -> 'dict[str, Any]'`。
- **输入：** 设计映射、`helix`、`start`、`end`、`deletion_offset`。
- **输出：** 新的 `design`。
- **示例代码：**

```python
from dnakit.genomics.scadnano_p26_functions import (
    add_deletion_in_range,
    create_origami_rectangle,
)

design = create_origami_rectangle(num_helices=4, num_cols=10, assign_seq=False)["design"]

r = add_deletion_in_range(design, helix=1, start=17, end=30)
print(sum(r["strand_lengths"]) - sum(create_origami_rectangle(num_helices=4, num_cols=10, assign_seq=False)["strand_lengths"]))
```

- **示例输出：**

```text
-2
```
