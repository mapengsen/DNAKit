# tskit 与 toytree：树序列与树的绘制

这一组函数是基于两个可选依赖的适配器。[tskit](https://github.com/tskit-dev/tskit) 是树序列工具包：`dnakit.population` 用它绘制树序列（SVG、HTML 和文本）、生成随机核苷酸串，并对 `.trees` 文件复现 `python -m tskit` 的命令（汇总信息、列出各树、按表导出文本、导出 FASTA 和 VCF）。[toytree](https://github.com/eaton-lab/toytree)（渲染由 [toyplot](https://github.com/sandialabs/toyplot) 完成）是系统发育树的 Python 库：`dnakit.phylogeny` 用它绘制单棵树、树网格和树云，按数值特征给节点着色或设定大小，导出节点数据表，以及修剪和定根。tskit 的函数从 `dnakit.population` 导入，toytree 的函数从 `dnakit.phylogeny` 导入，例如 `from dnakit.population import tree_file_info`。它们都不预测或模拟数据，只读取、转换并绘制你给定的内容。

tskit 使用 MIT 许可。toytree 使用 GPL-3.0 许可（其包元数据标明 “GNU General Public License v3”）；DNAKit 只在这些函数内部导入它，不复制其任何代码。toyplot 使用 BSD 许可。

依赖都是可选的，只有调用对应函数时才会导入。缺少包时函数会抛出 `BackendUnavailableError` 并给出安装命令（`dnakit[bio-tskit]` 或 `dnakit[bio-toytree]`）。绘图函数从不覆盖文件：输出路径必须是新文件，所在目录必须存在。

## 示例数据

读取 `.trees` 文件的函数（第 6 到 17 节）需要二进制的树序列文件，无法用文本展示。因此它们的示例先用 `tskit.load_text()` 读取下面七张很小的文本表，加入一条溯源记录，再把结果写成临时的 `example.trees`。该树序列有 4 个样本节点（两个二倍体个体）、长度 100、两棵树（`[0, 50)` 和 `[50, 100)`）、两个位点和两个突变。tskit 的绘图示例（第 1 到 5 节和第 18 节）把同一个树序列以表格行的形式直接传入，toytree 的示例直接传入 Newick 文本，所以都不需要文件。

**示例文件：** `nodes.txt`

```text
is_sample time population individual
1 0 0 0
1 0 0 0
1 0 0 1
1 0 0 1
0 1 0 -1
0 1 0 -1
0 2 0 -1
0 1 0 -1
0 1 0 -1
0 3 0 -1
```

**示例文件：** `edges.txt`

```text
left right parent child
0 50 4 0
0 50 4 1
0 50 5 2
0 50 5 3
0 50 6 4
0 50 6 5
50 100 7 0
50 100 7 2
50 100 8 1
50 100 8 3
50 100 9 7
50 100 9 8
```

**示例文件：** `sites.txt`

```text
position ancestral_state
10 A
60 C
```

**示例文件：** `mutations.txt`

```text
site node derived_state
0 4 G
1 7 T
```

**示例文件：** `individuals.txt`

```text
flags location parents
0 0.5 -1,-1
0 1.5 -1,-1
```

**示例文件：** `populations.txt`

```text
metadata
POP1
```

**示例文件：** `migrations.txt`

```text
source dest node time left right
0 0 4 0.5 0 50
```

## 1) `dnakit.population.draw_tree_html()`

- **作用：** 把树序列中的一棵树写成 tskit 的 HTML 视图，保存为新文件。
- **处理方法：** 用 tskit 由显式表格行构建树序列，用 `ts.at_index(tree_index)` 取出这棵树，并把 `tskit.tree_html(tree)` 的文本以 UTF-8 写入文件。目标文件必须不存在，所在目录必须存在。
- **API：** `draw_tree_html(tables: 'Mapping[str, Any]', tree_index: 'int', destination: 'str | Path', *, max_rows: 'int' = 100000, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 tskit（extra `dnakit[bio-tskit]`），不需要外部程序、网络或模型。
- **输入：** `tables` 是映射，含 `sequence_length`、`nodes`（`time`，可选 `sample`）、`edges`（`left`、`right`、`parent`、`child`）、`sites`（`position`、`ancestral_state`）和 `mutations`（`site`、`node`、`derived_state`）；`max_rows` 限制每张表的行数（上限 1,000,000）。`tree_index` 为 0 到 `num_trees - 1` 的整数。`destination` 为新文件路径。
- **输出：** 字典，含 `status`、`operation`、`image_path`（HTML 文件）、`format`（`html`）、`bytes`、`tree_index`、`interval`（`[left, right)`）、`num_roots`、`num_sites` 和 `provenance`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.population import draw_tree_html

tables = {
    "sequence_length": 100,
    "nodes": [{"time": 0, "sample": True}] * 4 + [{"time": t} for t in (1, 1, 2, 1, 1, 3)],
    "edges": [
        {"left": left, "right": right, "parent": parent, "child": child}
        for left, right, parent, child in [
            (0, 50, 4, 0), (0, 50, 4, 1), (0, 50, 5, 2), (0, 50, 5, 3),
            (0, 50, 6, 4), (0, 50, 6, 5), (50, 100, 7, 0), (50, 100, 7, 2),
            (50, 100, 8, 1), (50, 100, 8, 3), (50, 100, 9, 7), (50, 100, 9, 8),
        ]
    ],
    "sites": [
        {"position": 10, "ancestral_state": "A"},
        {"position": 60, "ancestral_state": "C"},
    ],
    "mutations": [
        {"site": 0, "node": 4, "derived_state": "G"},
        {"site": 1, "node": 7, "derived_state": "T"},
    ],
}

with tempfile.TemporaryDirectory() as tmp:
    out = Path(tmp) / "tree1.html"
    result = draw_tree_html(tables, 1, out, verbose=False)
    html = out.read_text(encoding="utf-8")
    print(result["format"], result["tree_index"], result["interval"])
    print(result["num_roots"], result["num_sites"], "<div>" in html)
```

- **示例输出：**

```text
html 1 [50.0, 100.0]
1 1 True
```

## 2) `dnakit.population.draw_tree_sequence_html()`

- **作用：** 把整个树序列的 tskit HTML 摘要写成新文件。
- **处理方法：** 由显式表格行构建树序列，把 `tskit.tree_sequence_html(ts)` 的文本以 UTF-8 写入文件。摘要列出各表大小、溯源记录和元数据。目标文件必须不存在。
- **API：** `draw_tree_sequence_html(tables: 'Mapping[str, Any]', destination: 'str | Path', *, max_rows: 'int' = 100000, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 tskit（extra `dnakit[bio-tskit]`），不需要外部程序、网络或模型。
- **输入：** `tables` 与 `draw_tree_html()` 相同；`destination` 为新文件路径；`max_rows` 限制每张表的行数。
- **输出：** 字典，含 `status`、`operation`、`image_path`、`format`（`html`）、`bytes`、`sequence_length`、`num_trees`、`num_samples`、`num_sites`、`num_mutations`、`num_provenances` 和 `provenance`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.population import draw_tree_sequence_html

tables = {
    "sequence_length": 100,
    "nodes": [{"time": 0, "sample": True}] * 4 + [{"time": t} for t in (1, 1, 2, 1, 1, 3)],
    "edges": [
        {"left": left, "right": right, "parent": parent, "child": child}
        for left, right, parent, child in [
            (0, 50, 4, 0), (0, 50, 4, 1), (0, 50, 5, 2), (0, 50, 5, 3),
            (0, 50, 6, 4), (0, 50, 6, 5), (50, 100, 7, 0), (50, 100, 7, 2),
            (50, 100, 8, 1), (50, 100, 8, 3), (50, 100, 9, 7), (50, 100, 9, 8),
        ]
    ],
    "sites": [
        {"position": 10, "ancestral_state": "A"},
        {"position": 60, "ancestral_state": "C"},
    ],
    "mutations": [
        {"site": 0, "node": 4, "derived_state": "G"},
        {"site": 1, "node": 7, "derived_state": "T"},
    ],
}

with tempfile.TemporaryDirectory() as tmp:
    out = Path(tmp) / "summary.html"
    result = draw_tree_sequence_html(tables, out, verbose=False)
    html = out.read_text(encoding="utf-8")
    print(result["format"], result["sequence_length"], result["num_trees"])
    print(result["num_samples"], result["num_sites"], result["num_mutations"], result["num_provenances"])
    print("Provenances" in html)
```

- **示例输出：**

```text
html 100.0 2
4 2 2 0
True
```

## 3) `dnakit.population.draw_tree_sequence_svg()`

- **作用：** 把树序列或其中一棵树画成 SVG 文件。
- **处理方法：** 由显式表格行构建树序列。`tree_index=None` 时调用 `TreeSequence.draw_svg(path=...)`（绘制全部树，可用 `x_lim` 限定区间）；给出索引时对 `ts.at_index(index)` 调用 `Tree.draw_svg(path=...)`。目标文件必须不存在。
- **API：** `draw_tree_sequence_svg(tables: 'Mapping[str, Any]', destination: 'str | Path', *, tree_index: 'int | None' = None, x_lim: 'Sequence[float] | None' = None, max_rows: 'int' = 100000, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 tskit（extra `dnakit[bio-tskit]`），不需要外部程序、网络或模型。
- **输入：** `tables` 与 `draw_tree_html()` 相同；`destination` 为新文件路径；`tree_index` 为 null 或 0 到 `num_trees - 1` 的整数；`x_lim` 为 null 或两个数 `[left, right]`，要求 `0 <= left < right <= sequence_length`，只适用于整个序列（与 `tree_index` 同时给出会抛出 `ConfigurationError`）。
- **输出：** 字典，含 `status`、`operation`、`image_path`、`format`（`svg`）、`scope`（`tree_sequence` 或 `tree`）、`tree_index`、`x_lim`、`bytes`、`sequence_length`、`num_trees`、`intervals`（每棵树一个 `[left, right]`）、`num_samples`、`num_sites`、`num_mutations` 和 `provenance`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.population import draw_tree_sequence_svg

tables = {
    "sequence_length": 100,
    "nodes": [{"time": 0, "sample": True}] * 4 + [{"time": t} for t in (1, 1, 2, 1, 1, 3)],
    "edges": [
        {"left": left, "right": right, "parent": parent, "child": child}
        for left, right, parent, child in [
            (0, 50, 4, 0), (0, 50, 4, 1), (0, 50, 5, 2), (0, 50, 5, 3),
            (0, 50, 6, 4), (0, 50, 6, 5), (50, 100, 7, 0), (50, 100, 7, 2),
            (50, 100, 8, 1), (50, 100, 8, 3), (50, 100, 9, 7), (50, 100, 9, 8),
        ]
    ],
    "sites": [
        {"position": 10, "ancestral_state": "A"},
        {"position": 60, "ancestral_state": "C"},
    ],
    "mutations": [
        {"site": 0, "node": 4, "derived_state": "G"},
        {"site": 1, "node": 7, "derived_state": "T"},
    ],
}

with tempfile.TemporaryDirectory() as tmp:
    whole = draw_tree_sequence_svg(tables, Path(tmp) / "all.svg", x_lim=[0, 60], verbose=False)
    single = draw_tree_sequence_svg(tables, Path(tmp) / "tree1.svg", tree_index=1, verbose=False)
    print(whole["scope"], whole["x_lim"], whole["intervals"], whole["num_sites"])
    print(single["scope"], single["tree_index"], single["num_trees"])
    print("<svg" in (Path(tmp) / "tree1.svg").read_text(encoding="utf-8"))
```

- **示例输出：**

```text
tree_sequence [0.0, 60.0] [[0.0, 50.0], [50.0, 100.0]] 2
tree 1 2
True
```

## 4) `dnakit.population.draw_variant_html()`

- **作用：** 把一个变异（一个位点）的 tskit HTML 视图写成新文件。
- **处理方法：** 由显式表格行构建树序列，取 `ts.variants()`（位点按位置排序）的第 `site_index` 项，把 `tskit.variant_html(variant)` 的文本以 UTF-8 写入文件。目标文件必须不存在。
- **API：** `draw_variant_html(tables: 'Mapping[str, Any]', site_index: 'int', destination: 'str | Path', *, max_rows: 'int' = 100000, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 tskit（extra `dnakit[bio-tskit]`），不需要外部程序、网络或模型。
- **输入：** `tables` 与 `draw_tree_html()` 相同；`site_index` 为 0 到 `num_sites - 1` 的整数；`destination` 为新文件路径。
- **输出：** 字典，含 `status`、`operation`、`image_path`、`format`（`html`）、`bytes`、`site_index`、`position`、`alleles`（缺失等位基因为 null）、`num_samples` 和 `provenance`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.population import draw_variant_html

tables = {
    "sequence_length": 100,
    "nodes": [{"time": 0, "sample": True}] * 4 + [{"time": t} for t in (1, 1, 2, 1, 1, 3)],
    "edges": [
        {"left": left, "right": right, "parent": parent, "child": child}
        for left, right, parent, child in [
            (0, 50, 4, 0), (0, 50, 4, 1), (0, 50, 5, 2), (0, 50, 5, 3),
            (0, 50, 6, 4), (0, 50, 6, 5), (50, 100, 7, 0), (50, 100, 7, 2),
            (50, 100, 8, 1), (50, 100, 8, 3), (50, 100, 9, 7), (50, 100, 9, 8),
        ]
    ],
    "sites": [
        {"position": 10, "ancestral_state": "A"},
        {"position": 60, "ancestral_state": "C"},
    ],
    "mutations": [
        {"site": 0, "node": 4, "derived_state": "G"},
        {"site": 1, "node": 7, "derived_state": "T"},
    ],
}

with tempfile.TemporaryDirectory() as tmp:
    out = Path(tmp) / "site1.html"
    result = draw_variant_html(tables, 1, out, verbose=False)
    print(result["site_index"], result["position"], result["alleles"], result["num_samples"])
    print("<div>" in out.read_text(encoding="utf-8"))
```

- **示例输出：**

```text
1 60.0 ['C', 'T'] 4
True
```

## 5) `dnakit.population.tskit_random_nucleotides()`

- **作用：** 按给定长度生成随机核苷酸串，可用种子复现。
- **处理方法：** 调用 `tskit.random_nucleotides(length, seed=seed)`，它用 numpy 的 `RandomState` 从 `ACTG` 中均匀抽取。相同种子得到与 tskit 完全相同的字符串；它与 DNAKit 的 `generate_random_dna()` 不是同一随机流。
- **API：** `tskit_random_nucleotides(*, length: 'int', seed: 'int | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 tskit（extra `dnakit[bio-tskit]`），不需要外部程序、网络或模型。
- **输入：** `length` 为 0 到 10,000,000 的整数（仅关键字参数）；`seed` 为 null（每次不同）或 0 到 2^32 - 1 的整数。
- **输出：** 字典，含 `status`、`operation`、`length`、`seed`、`alphabet`（`ACTG`）、`sequence` 和 `provenance`。
- **示例代码：**

```python
from dnakit.population import tskit_random_nucleotides

first = tskit_random_nucleotides(length=16, seed=7, verbose=False)
second = tskit_random_nucleotides(length=16, seed=7, verbose=False)
print(first["sequence"])
print(first["sequence"] == second["sequence"], first["length"], first["alphabet"])
```

- **示例输出：**

```text
GACTGGGGACTGACTT
True 16 ACTG
```

## 6) `dnakit.population.export_tree_file_edges_table()`

- **作用：** 把 `.trees` 文件的边表写成文本。
- **处理方法：** 先检查文件大小，用 `tskit.load()` 读取，再调用 `TreeSequence.dump_text(edges=handle, precision=precision)` 写入新文件。这与 tskit 的 `edges` 命令输出一致。
- **API：** `export_tree_file_edges_table(tree_sequence_path: 'str | Path', destination: 'str | Path', *, precision: 'int' = 6, max_bytes: 'int' = 100000000, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 tskit（extra `dnakit[bio-tskit]`），不需要外部程序、网络或模型。
- **输入：** `tree_sequence_path` 为已存在的 `.trees` 文件，大小不超过 `max_bytes`（默认 100,000,000，上限 1,000,000,000）；`destination` 为新文件路径，所在目录必须存在；`precision` 为 0 到 16 的整数（默认 6）。
- **输出：** 字典，含 `status`、`operation`、`source`、`destination`、`table`（表名）、`row_count`（不含表头的行数）、`bytes` 和 `provenance`。文本文件含一行表头，之后每行一条记录。
- **示例代码：**

```python
import json
import tempfile
from pathlib import Path

import tskit

from dnakit.population import export_tree_file_edges_table

names = ["nodes", "edges", "sites", "mutations", "individuals", "populations", "migrations"]
tables = tskit.load_text(**{name: open(f"{name}.txt") for name in names}, strict=False).dump_tables()
tables.provenances.add_row(json.dumps({"software": "demo"}), "2024-01-01T00:00:00")

with tempfile.TemporaryDirectory() as tmp:
    tree_file = Path(tmp) / "example.trees"
    tables.tree_sequence().dump(str(tree_file))
    out = Path(tmp) / "edges.txt"
    result = export_tree_file_edges_table(tree_file, out, verbose=False)
    print(result["table"], result["row_count"])
    for line in out.read_text(encoding="utf8").splitlines():
        print(line.replace("\t", " | ").rstrip())
```

- **示例输出：**

```text
edges 12
left | right | parent | child | metadata
0.000000 | 50.000000 | 4 | 0 |
0.000000 | 50.000000 | 4 | 1 |
0.000000 | 50.000000 | 5 | 2 |
0.000000 | 50.000000 | 5 | 3 |
50.000000 | 100.000000 | 7 | 0 |
50.000000 | 100.000000 | 7 | 2 |
50.000000 | 100.000000 | 8 | 1 |
50.000000 | 100.000000 | 8 | 3 |
0.000000 | 50.000000 | 6 | 4 |
0.000000 | 50.000000 | 6 | 5 |
50.000000 | 100.000000 | 9 | 7 |
50.000000 | 100.000000 | 9 | 8 |
```

## 7) `dnakit.population.export_tree_file_fasta()`

- **作用：** 把 `.trees` 文件中样本的单倍型序列写成 FASTA。
- **处理方法：** 用 `tskit.load()` 读取文件，调用 `TreeSequence.write_fasta(path, wrap_width=wrap_width)`，即 tskit 的 `run_fasta()` 处理函数所做的操作。每个样本节点对应一条名为 `n<节点编号>` 的记录。文件没有参考序列，因此没有位点的位置写作 `N`。
- **API：** `export_tree_file_fasta(tree_sequence_path: 'str | Path', destination: 'str | Path', *, wrap_width: 'int' = 60, max_bytes: 'int' = 100000000, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 tskit（extra `dnakit[bio-tskit]`），不需要外部程序、网络或模型。
- **输入：** `tree_sequence_path` 为已存在的 `.trees` 文件，大小不超过 `max_bytes`；`destination` 为新文件路径；`wrap_width` 为不小于 0 的整数（默认 60，0 表示每条序列写在一行）。
- **输出：** 字典，含 `status`、`operation`、`source`、`destination`、`file_format`（`fasta`）、`wrap_width`、`bytes`、`num_samples`、`sequence_length`、`num_sites` 和 `provenance`。
- **示例代码：**

```python
import json
import tempfile
from pathlib import Path

import tskit

from dnakit.population import export_tree_file_fasta

names = ["nodes", "edges", "sites", "mutations", "individuals", "populations", "migrations"]
tables = tskit.load_text(**{name: open(f"{name}.txt") for name in names}, strict=False).dump_tables()
tables.provenances.add_row(json.dumps({"software": "demo"}), "2024-01-01T00:00:00")

with tempfile.TemporaryDirectory() as tmp:
    tree_file = Path(tmp) / "example.trees"
    tables.tree_sequence().dump(str(tree_file))
    out = Path(tmp) / "samples.fa"
    result = export_tree_file_fasta(tree_file, out, wrap_width=0, verbose=False)
    print(result["file_format"], result["num_samples"], result["sequence_length"], result["num_sites"])
    for line in out.read_text(encoding="utf8").splitlines():
        print(line[:12], len(line))
```

- **示例输出：**

```text
fasta 4 100.0 2
>n0 3
NNNNNNNNNNGN 100
>n1 3
NNNNNNNNNNGN 100
>n2 3
NNNNNNNNNNAN 100
>n3 3
NNNNNNNNNNAN 100
```

## 8) `dnakit.population.export_tree_file_individuals_table()`

- **作用：** 把 `.trees` 文件的个体表写成文本。
- **处理方法：** 先检查文件大小，用 `tskit.load()` 读取，再调用 `TreeSequence.dump_text(individuals=handle, precision=precision)` 写入新文件。这与 tskit 的 `individuals` 命令输出一致。
- **API：** `export_tree_file_individuals_table(tree_sequence_path: 'str | Path', destination: 'str | Path', *, precision: 'int' = 6, max_bytes: 'int' = 100000000, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 tskit（extra `dnakit[bio-tskit]`），不需要外部程序、网络或模型。
- **输入：** `tree_sequence_path` 为已存在的 `.trees` 文件，大小不超过 `max_bytes`（默认 100,000,000，上限 1,000,000,000）；`destination` 为新文件路径，所在目录必须存在；`precision` 为 0 到 16 的整数（默认 6）。
- **输出：** 字典，含 `status`、`operation`、`source`、`destination`、`table`（表名）、`row_count`（不含表头的行数）、`bytes` 和 `provenance`。文本文件含一行表头，之后每行一条记录。
- **示例代码：**

```python
import json
import tempfile
from pathlib import Path

import tskit

from dnakit.population import export_tree_file_individuals_table

names = ["nodes", "edges", "sites", "mutations", "individuals", "populations", "migrations"]
tables = tskit.load_text(**{name: open(f"{name}.txt") for name in names}, strict=False).dump_tables()
tables.provenances.add_row(json.dumps({"software": "demo"}), "2024-01-01T00:00:00")

with tempfile.TemporaryDirectory() as tmp:
    tree_file = Path(tmp) / "example.trees"
    tables.tree_sequence().dump(str(tree_file))
    out = Path(tmp) / "individuals.txt"
    result = export_tree_file_individuals_table(tree_file, out, verbose=False)
    print(result["table"], result["row_count"])
    for line in out.read_text(encoding="utf8").splitlines():
        print(line.replace("\t", " | ").rstrip())
```

- **示例输出：**

```text
individuals 2
id | flags | location | parents | metadata
0 | 0 | 0.5 | -1,-1 |
1 | 0 | 1.5 | -1,-1 |
```

## 9) `dnakit.population.export_tree_file_migrations_table()`

- **作用：** 把 `.trees` 文件的迁移表写成文本。
- **处理方法：** 先检查文件大小，用 `tskit.load()` 读取，再调用 `TreeSequence.dump_text(migrations=handle, precision=precision)` 写入新文件。这与 tskit 的 `migrations` 命令输出一致。
- **API：** `export_tree_file_migrations_table(tree_sequence_path: 'str | Path', destination: 'str | Path', *, precision: 'int' = 6, max_bytes: 'int' = 100000000, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 tskit（extra `dnakit[bio-tskit]`），不需要外部程序、网络或模型。
- **输入：** `tree_sequence_path` 为已存在的 `.trees` 文件，大小不超过 `max_bytes`（默认 100,000,000，上限 1,000,000,000）；`destination` 为新文件路径，所在目录必须存在；`precision` 为 0 到 16 的整数（默认 6）。
- **输出：** 字典，含 `status`、`operation`、`source`、`destination`、`table`（表名）、`row_count`（不含表头的行数）、`bytes` 和 `provenance`。文本文件含一行表头，之后每行一条记录。
- **示例代码：**

```python
import json
import tempfile
from pathlib import Path

import tskit

from dnakit.population import export_tree_file_migrations_table

names = ["nodes", "edges", "sites", "mutations", "individuals", "populations", "migrations"]
tables = tskit.load_text(**{name: open(f"{name}.txt") for name in names}, strict=False).dump_tables()
tables.provenances.add_row(json.dumps({"software": "demo"}), "2024-01-01T00:00:00")

with tempfile.TemporaryDirectory() as tmp:
    tree_file = Path(tmp) / "example.trees"
    tables.tree_sequence().dump(str(tree_file))
    out = Path(tmp) / "migrations.txt"
    result = export_tree_file_migrations_table(tree_file, out, verbose=False)
    print(result["table"], result["row_count"])
    for line in out.read_text(encoding="utf8").splitlines():
        print(line.replace("\t", " | ").rstrip())
```

- **示例输出：**

```text
migrations 1
left | right | node | source | dest | time | metadata
0.0 | 50.0 | 4 | 0 | 0 | 0.5 |  |
```

## 10) `dnakit.population.export_tree_file_mutations_table()`

- **作用：** 把 `.trees` 文件的突变表写成文本。
- **处理方法：** 先检查文件大小，用 `tskit.load()` 读取，再调用 `TreeSequence.dump_text(mutations=handle, precision=precision)` 写入新文件。这与 tskit 的 `mutations` 命令输出一致。
- **API：** `export_tree_file_mutations_table(tree_sequence_path: 'str | Path', destination: 'str | Path', *, precision: 'int' = 6, max_bytes: 'int' = 100000000, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 tskit（extra `dnakit[bio-tskit]`），不需要外部程序、网络或模型。
- **输入：** `tree_sequence_path` 为已存在的 `.trees` 文件，大小不超过 `max_bytes`（默认 100,000,000，上限 1,000,000,000）；`destination` 为新文件路径，所在目录必须存在；`precision` 为 0 到 16 的整数（默认 6）。
- **输出：** 字典，含 `status`、`operation`、`source`、`destination`、`table`（表名）、`row_count`（不含表头的行数）、`bytes` 和 `provenance`。文本文件含一行表头，之后每行一条记录。
- **示例代码：**

```python
import json
import tempfile
from pathlib import Path

import tskit

from dnakit.population import export_tree_file_mutations_table

names = ["nodes", "edges", "sites", "mutations", "individuals", "populations", "migrations"]
tables = tskit.load_text(**{name: open(f"{name}.txt") for name in names}, strict=False).dump_tables()
tables.provenances.add_row(json.dumps({"software": "demo"}), "2024-01-01T00:00:00")

with tempfile.TemporaryDirectory() as tmp:
    tree_file = Path(tmp) / "example.trees"
    tables.tree_sequence().dump(str(tree_file))
    out = Path(tmp) / "mutations.txt"
    result = export_tree_file_mutations_table(tree_file, out, verbose=False)
    print(result["table"], result["row_count"])
    for line in out.read_text(encoding="utf8").splitlines():
        print(line.replace("\t", " | ").rstrip())
```

- **示例输出：**

```text
mutations 2
site | node | time | derived_state | parent | metadata
0 | 4 | unknown | G | -1 |
1 | 7 | unknown | T | -1 |
```

## 11) `dnakit.population.export_tree_file_nodes_table()`

- **作用：** 把 `.trees` 文件的节点表写成文本。
- **处理方法：** 先检查文件大小，用 `tskit.load()` 读取，再调用 `TreeSequence.dump_text(nodes=handle, precision=precision)` 写入新文件。这与 tskit 的 `nodes` 命令输出一致。
- **API：** `export_tree_file_nodes_table(tree_sequence_path: 'str | Path', destination: 'str | Path', *, precision: 'int' = 6, max_bytes: 'int' = 100000000, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 tskit（extra `dnakit[bio-tskit]`），不需要外部程序、网络或模型。
- **输入：** `tree_sequence_path` 为已存在的 `.trees` 文件，大小不超过 `max_bytes`（默认 100,000,000，上限 1,000,000,000）；`destination` 为新文件路径，所在目录必须存在；`precision` 为 0 到 16 的整数（默认 6）。
- **输出：** 字典，含 `status`、`operation`、`source`、`destination`、`table`（表名）、`row_count`（不含表头的行数）、`bytes` 和 `provenance`。文本文件含一行表头，之后每行一条记录。
- **示例代码：**

```python
import json
import tempfile
from pathlib import Path

import tskit

from dnakit.population import export_tree_file_nodes_table

names = ["nodes", "edges", "sites", "mutations", "individuals", "populations", "migrations"]
tables = tskit.load_text(**{name: open(f"{name}.txt") for name in names}, strict=False).dump_tables()
tables.provenances.add_row(json.dumps({"software": "demo"}), "2024-01-01T00:00:00")

with tempfile.TemporaryDirectory() as tmp:
    tree_file = Path(tmp) / "example.trees"
    tables.tree_sequence().dump(str(tree_file))
    out = Path(tmp) / "nodes.txt"
    result = export_tree_file_nodes_table(tree_file, out, verbose=False)
    print(result["table"], result["row_count"])
    for line in out.read_text(encoding="utf8").splitlines():
        print(line.replace("\t", " | ").rstrip())
```

- **示例输出：**

```text
nodes 10
id | is_sample | time | population | individual | metadata
0 | 1 | 0.000000 | 0 | 0 |
1 | 1 | 0.000000 | 0 | 0 |
2 | 1 | 0.000000 | 0 | 1 |
3 | 1 | 0.000000 | 0 | 1 |
4 | 0 | 1.000000 | 0 | -1 |
5 | 0 | 1.000000 | 0 | -1 |
6 | 0 | 2.000000 | 0 | -1 |
7 | 0 | 1.000000 | 0 | -1 |
8 | 0 | 1.000000 | 0 | -1 |
9 | 0 | 3.000000 | 0 | -1 |
```

## 12) `dnakit.population.export_tree_file_populations_table()`

- **作用：** 把 `.trees` 文件的群体表写成文本。
- **处理方法：** 先检查文件大小，用 `tskit.load()` 读取，再调用 `TreeSequence.dump_text(populations=handle)` 写入新文件。这与 tskit 的 `populations` 命令输出一致。
- **API：** `export_tree_file_populations_table(tree_sequence_path: 'str | Path', destination: 'str | Path', *, max_bytes: 'int' = 100000000, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 tskit（extra `dnakit[bio-tskit]`），不需要外部程序、网络或模型。
- **输入：** `tree_sequence_path` 为已存在的 `.trees` 文件，大小不超过 `max_bytes`（默认 100,000,000，上限 1,000,000,000）；`destination` 为新文件路径，所在目录必须存在；此函数没有 `precision` 参数。
- **输出：** 字典，含 `status`、`operation`、`source`、`destination`、`table`（表名）、`row_count`（不含表头的行数）、`bytes` 和 `provenance`。文本文件含一行表头，之后每行一条记录。
- **示例代码：**

```python
import json
import tempfile
from pathlib import Path

import tskit

from dnakit.population import export_tree_file_populations_table

names = ["nodes", "edges", "sites", "mutations", "individuals", "populations", "migrations"]
tables = tskit.load_text(**{name: open(f"{name}.txt") for name in names}, strict=False).dump_tables()
tables.provenances.add_row(json.dumps({"software": "demo"}), "2024-01-01T00:00:00")

with tempfile.TemporaryDirectory() as tmp:
    tree_file = Path(tmp) / "example.trees"
    tables.tree_sequence().dump(str(tree_file))
    out = Path(tmp) / "populations.txt"
    result = export_tree_file_populations_table(tree_file, out, verbose=False)
    print(result["table"], result["row_count"])
    for line in out.read_text(encoding="utf8").splitlines():
        print(line.replace("\t", " | ").rstrip())
```

- **示例输出：**

```text
populations 1
id | metadata
0 | POP1
```

## 13) `dnakit.population.export_tree_file_provenances_table()`

- **作用：** 把 `.trees` 文件的溯源记录写成表格文本或可读的分块文本。
- **处理方法：** 用 `tskit.load()` 读取文件。默认调用 `TreeSequence.dump_text(provenances=handle)` 写出表格文本。`human_readable=True` 时每条记录写成一块，格式为 `id=..., timestamp=..., record=<缩进 4 个空格的 JSON>`，与 tskit 的 `-H` 选项一致；记录不是合法 JSON 时抛出 `ConfigurationError`。
- **API：** `export_tree_file_provenances_table(tree_sequence_path: 'str | Path', destination: 'str | Path', *, human_readable: 'bool' = False, max_bytes: 'int' = 100000000, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 tskit（extra `dnakit[bio-tskit]`），不需要外部程序、网络或模型。
- **输入：** `tree_sequence_path` 为已存在的 `.trees` 文件，大小不超过 `max_bytes`；`destination` 为新文件路径；`human_readable` 为布尔值（默认 false）。
- **输出：** 字典，含 `status`、`operation`、`source`、`destination`、`table`（`provenances`）、`row_count`（记录数）、`bytes` 和 `provenance`；可读形式还有值为 true 的 `human_readable`。
- **示例代码：**

```python
import json
import tempfile
from pathlib import Path

import tskit

from dnakit.population import export_tree_file_provenances_table

names = ["nodes", "edges", "sites", "mutations", "individuals", "populations", "migrations"]
tables = tskit.load_text(**{name: open(f"{name}.txt") for name in names}, strict=False).dump_tables()
tables.provenances.add_row(json.dumps({"software": "demo"}), "2024-01-01T00:00:00")

with tempfile.TemporaryDirectory() as tmp:
    tree_file = Path(tmp) / "example.trees"
    tables.tree_sequence().dump(str(tree_file))
    table = export_tree_file_provenances_table(tree_file, Path(tmp) / "table.txt", verbose=False)
    blocks = export_tree_file_provenances_table(
        tree_file, Path(tmp) / "blocks.txt", human_readable=True, verbose=False
    )
    print(table["row_count"], blocks["row_count"], blocks["human_readable"])
    for line in (Path(tmp) / "table.txt").read_text(encoding="utf8").splitlines():
        print(line.replace("\t", " | ").rstrip())
    print((Path(tmp) / "blocks.txt").read_text(encoding="utf8"), end="")
```

- **示例输出：**

```text
1 1 True
id | timestamp | record
0 | 2024-01-01T00:00:00 | {"software": "demo"} |
id=0, timestamp=2024-01-01T00:00:00, record={
    "software": "demo"
}
```

## 14) `dnakit.population.export_tree_file_sites_table()`

- **作用：** 把 `.trees` 文件的位点表写成文本。
- **处理方法：** 先检查文件大小，用 `tskit.load()` 读取，再调用 `TreeSequence.dump_text(sites=handle, precision=precision)` 写入新文件。这与 tskit 的 `sites` 命令输出一致。
- **API：** `export_tree_file_sites_table(tree_sequence_path: 'str | Path', destination: 'str | Path', *, precision: 'int' = 6, max_bytes: 'int' = 100000000, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 tskit（extra `dnakit[bio-tskit]`），不需要外部程序、网络或模型。
- **输入：** `tree_sequence_path` 为已存在的 `.trees` 文件，大小不超过 `max_bytes`（默认 100,000,000，上限 1,000,000,000）；`destination` 为新文件路径，所在目录必须存在；`precision` 为 0 到 16 的整数（默认 6）。
- **输出：** 字典，含 `status`、`operation`、`source`、`destination`、`table`（表名）、`row_count`（不含表头的行数）、`bytes` 和 `provenance`。文本文件含一行表头，之后每行一条记录。
- **示例代码：**

```python
import json
import tempfile
from pathlib import Path

import tskit

from dnakit.population import export_tree_file_sites_table

names = ["nodes", "edges", "sites", "mutations", "individuals", "populations", "migrations"]
tables = tskit.load_text(**{name: open(f"{name}.txt") for name in names}, strict=False).dump_tables()
tables.provenances.add_row(json.dumps({"software": "demo"}), "2024-01-01T00:00:00")

with tempfile.TemporaryDirectory() as tmp:
    tree_file = Path(tmp) / "example.trees"
    tables.tree_sequence().dump(str(tree_file))
    out = Path(tmp) / "sites.txt"
    result = export_tree_file_sites_table(tree_file, out, verbose=False)
    print(result["table"], result["row_count"])
    for line in out.read_text(encoding="utf8").splitlines():
        print(line.replace("\t", " | ").rstrip())
```

- **示例输出：**

```text
sites 2
position | ancestral_state | metadata
10.000000 | A |
60.000000 | C |
```

## 15) `dnakit.population.export_tree_file_vcf()`

- **作用：** 把 `.trees` 文件中的变异写成 VCF 文件。
- **处理方法：** 用 `tskit.load()` 读取文件，调用 `TreeSequence.write_vcf(handle, ploidy=..., contig_id=..., allow_position_zero=...)`，即 tskit 的 `run_vcf()` 处理函数所做的操作。tskit 按个体各写一列（文件没有个体信息时按每 `ploidy` 个样本一列）。tskit 出错时转为 `ConfigurationError` 并删除写了一半的文件。
- **API：** `export_tree_file_vcf(tree_sequence_path: 'str | Path', destination: 'str | Path', *, ploidy: 'int | None' = None, contig_id: 'str' = '1', allow_position_zero: 'bool' = False, max_bytes: 'int' = 100000000, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 tskit（extra `dnakit[bio-tskit]`），不需要外部程序、网络或模型。
- **输入：** `tree_sequence_path` 为已存在的 `.trees` 文件，大小不超过 `max_bytes`；`destination` 为新文件路径；`ploidy` 为 null（用 tskit 默认值）或能整除样本数的不小于 1 的整数（文件已含个体信息时 tskit 会拒绝）；`contig_id` 为 CHROM 列文本（默认 `1`）；`allow_position_zero` 为布尔值。
- **输出：** 字典，含 `status`、`operation`、`source`、`destination`、`file_format`（`vcf`）、`ploidy`、`contig_id`、`allow_position_zero`、`bytes`、`num_samples`、`num_records`（非表头行数）和 `provenance`。
- **示例代码：**

```python
import json
import tempfile
from pathlib import Path

import tskit

from dnakit.population import export_tree_file_vcf

names = ["nodes", "edges", "sites", "mutations", "individuals", "populations", "migrations"]
tables = tskit.load_text(**{name: open(f"{name}.txt") for name in names}, strict=False).dump_tables()
tables.provenances.add_row(json.dumps({"software": "demo"}), "2024-01-01T00:00:00")

with tempfile.TemporaryDirectory() as tmp:
    tree_file = Path(tmp) / "example.trees"
    tables.tree_sequence().dump(str(tree_file))
    out = Path(tmp) / "variants.vcf"
    result = export_tree_file_vcf(tree_file, out, contig_id="chr1", verbose=False)
    print(result["file_format"], result["contig_id"], result["num_samples"], result["num_records"])
    lines = out.read_text(encoding="utf8").splitlines()
    for line in lines:
        if not line.startswith("##"):
            print(line.replace("\t", " | ").rstrip())
```

- **示例输出：**

```text
vcf chr1 4 2
#CHROM | POS | ID | REF | ALT | QUAL | FILTER | INFO | FORMAT | tsk_0 | tsk_1
chr1 | 10 | 0 | A | G | . | PASS | . | GT | 1|1 | 0|0
chr1 | 60 | 1 | C | T | . | PASS | . | GT | 1|0 | 1|0
```

## 16) `dnakit.population.tree_file_info()`

- **作用：** 汇总 `.trees` 文件：tskit 的摘要表和各表计数。
- **处理方法：** 先检查文件大小，用 `tskit.load()` 读取，返回 `str(TreeSequence)` 以及从树序列读出的各项计数，与 tskit 的 `info` 命令一致。
- **API：** `tree_file_info(tree_sequence_path: 'str | Path', *, max_bytes: 'int' = 100000000, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 tskit（extra `dnakit[bio-tskit]`），不需要外部程序、网络或模型。
- **输入：** `tree_sequence_path` 为已存在的 `.trees` 文件，大小不超过 `max_bytes`（默认 100,000,000，上限 1,000,000,000）。
- **输出：** 字典，含 `status`、`operation`、`source`、`summary_text`（多行 tskit 表格，其中包含文件的内存大小）、`sequence_length`、`num_trees`、`num_samples`、`num_nodes`、`num_edges`、`num_sites`、`num_mutations`、`num_individuals`、`num_populations`、`num_migrations`、`num_provenances` 和 `provenance`。
- **示例代码：**

```python
import json
import tempfile
from pathlib import Path

import tskit

from dnakit.population import tree_file_info

names = ["nodes", "edges", "sites", "mutations", "individuals", "populations", "migrations"]
tables = tskit.load_text(**{name: open(f"{name}.txt") for name in names}, strict=False).dump_tables()
tables.provenances.add_row(json.dumps({"software": "demo"}), "2024-01-01T00:00:00")

with tempfile.TemporaryDirectory() as tmp:
    tree_file = Path(tmp) / "example.trees"
    tables.tree_sequence().dump(str(tree_file))
    result = tree_file_info(tree_file, verbose=False)
    print(result["sequence_length"], result["num_trees"], result["num_samples"], result["num_nodes"], result["num_edges"])
    print(result["num_sites"], result["num_mutations"], result["num_individuals"], result["num_populations"])
    print(result["num_migrations"], result["num_provenances"])
    print(result["summary_text"].splitlines()[3])
```

- **示例输出：**

```text
100.0 2 4 10 12
2 2 2 1
1 1
║Trees          │      2║
```

## 17) `dnakit.population.tree_file_trees()`

- **作用：** 列出 `.trees` 文件中每棵树的序号、区间和位点数。
- **处理方法：** 用 `tskit.load()` 读取文件并遍历 `TreeSequence.trees()`。区间按 `precision` 位小数取整，与 tskit 的 `trees` 命令一致。`draw=True` 时每条记录还包含 `Tree.draw(format="unicode")` 的文本。
- **API：** `tree_file_trees(tree_sequence_path: 'str | Path', *, precision: 'int' = 6, draw: 'bool' = False, max_bytes: 'int' = 100000000, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 tskit（extra `dnakit[bio-tskit]`），不需要外部程序、网络或模型。
- **输入：** `tree_sequence_path` 为已存在的 `.trees` 文件，大小不超过 `max_bytes`；`precision` 为 0 到 16 的整数（默认 6）；`draw` 为布尔值（默认 false）。
- **输出：** 字典，含 `status`、`operation`、`source`、`precision`、`num_trees`、`trees`（每项含 `index`、`num_sites`、`interval`，`draw` 时还有 `drawing`）和 `provenance`。
- **示例代码：**

```python
import json
import tempfile
from pathlib import Path

import tskit

from dnakit.population import tree_file_trees

names = ["nodes", "edges", "sites", "mutations", "individuals", "populations", "migrations"]
tables = tskit.load_text(**{name: open(f"{name}.txt") for name in names}, strict=False).dump_tables()
tables.provenances.add_row(json.dumps({"software": "demo"}), "2024-01-01T00:00:00")

with tempfile.TemporaryDirectory() as tmp:
    tree_file = Path(tmp) / "example.trees"
    tables.tree_sequence().dump(str(tree_file))
    result = tree_file_trees(tree_file, precision=2, draw=True, verbose=False)
    print(result["num_trees"], result["precision"])
    for tree in result["trees"]:
        print(tree["index"], tree["interval"], tree["num_sites"])
        print(tree["drawing"])
```

- **示例输出：**

```text
2 2
0 [0.0, 50.0] 1
   6   
 ┏━┻━┓ 
 4   5 
┏┻┓ ┏┻┓
0 1 2 3
1 [50.0, 100.0] 1
   9   
 ┏━┻━┓ 
 7   8 
┏┻┓ ┏┻┓
0 2 1 3
```

## 18) `dnakit.population.draw_tree_text()`

- **作用：** 把树序列中的一棵树写成 tskit 的文本图，保存为新文件。
- **处理方法：** 由显式表格行构建树序列，用 `ts.at_index(tree_index)` 取出这棵树，用 `Tree.draw(format=text_format)` 绘制并以 UTF-8 写入文件。结果同时返回图文本、该树的 Newick 字符串和区间。
- **API：** `draw_tree_text(tables: 'Mapping[str, Any]', tree_index: 'int', destination: 'str | Path', *, text_format: 'str' = 'unicode', max_rows: 'int' = 100000, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 tskit（extra `dnakit[bio-tskit]`），不需要外部程序、网络或模型。
- **输入：** `tables` 与 `draw_tree_html()` 相同；`tree_index` 为 0 到 `num_trees - 1` 的整数；`destination` 为新文件路径；`text_format` 为 `unicode`（默认，框线字符）或 `ascii`，其他值会抛出 `ConfigurationError`。
- **输出：** 字典，含 `status`、`operation`、`image_path`、`format`、`bytes`、`tree_index`、`interval`、`newick`、`num_roots`、`num_sites`、`text`（图文本）和 `provenance`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.population import draw_tree_text

tables = {
    "sequence_length": 100,
    "nodes": [{"time": 0, "sample": True}] * 4 + [{"time": t} for t in (1, 1, 2, 1, 1, 3)],
    "edges": [
        {"left": left, "right": right, "parent": parent, "child": child}
        for left, right, parent, child in [
            (0, 50, 4, 0), (0, 50, 4, 1), (0, 50, 5, 2), (0, 50, 5, 3),
            (0, 50, 6, 4), (0, 50, 6, 5), (50, 100, 7, 0), (50, 100, 7, 2),
            (50, 100, 8, 1), (50, 100, 8, 3), (50, 100, 9, 7), (50, 100, 9, 8),
        ]
    ],
    "sites": [
        {"position": 10, "ancestral_state": "A"},
        {"position": 60, "ancestral_state": "C"},
    ],
    "mutations": [
        {"site": 0, "node": 4, "derived_state": "G"},
        {"site": 1, "node": 7, "derived_state": "T"},
    ],
}

with tempfile.TemporaryDirectory() as tmp:
    unicode_result = draw_tree_text(tables, 0, Path(tmp) / "tree0.txt", verbose=False)
    ascii_result = draw_tree_text(
        tables, 1, Path(tmp) / "tree1.txt", text_format="ascii", verbose=False
    )
    print(unicode_result["interval"], unicode_result["newick"])
    print(unicode_result["text"])
    print(ascii_result["format"], ascii_result["interval"])
    print(ascii_result["text"])
```

- **示例输出：**

```text
[0.0, 50.0] ((1:1.00000000000000,2:1.00000000000000):1.00000000000000,(3:1.00000000000000,4:1.00000000000000):1.00000000000000);
   6   
 ┏━┻━┓ 
 4   5 
┏┻┓ ┏┻┓
0 1 2 3
ascii [50.0, 100.0]
   9   
 +-+-+ 
 7   8 
+++ +++
0 2 1 3
```

## 19) `dnakit.phylogeny.draw_tree_cloud_figure()`

- **作用：** 把多棵具有相同叶标签的树叠加成一朵树云并保存为图像。
- **处理方法：** 用 `toytree.tree()` 解析每棵树，检查它们的叶标签完全相同，构建 `toytree.mtree(trees)`，用 `MultiTree.draw_cloud_tree()` 绘制，再用 `toytree.save()` 保存画布。只有第一棵树显示叶标签，叶的顺序跟随一致树。
- **API：** `draw_tree_cloud_figure(trees: 'Sequence[str]', output_path: 'str | Path', *, image_format: 'str' = 'svg', width: 'int | None' = None, height: 'int | None' = None, max_tips: 'int' = 10000, max_file_bytes: 'int' = 10000000, verbose: bool = True) -> 'dict[str, object]'`。
- **运行要求：** 需要可选包 toytree（extra `dnakit[bio-toytree]`，会带入 toyplot、numpy 和 pandas），不需要外部程序、网络或模型。
- **输入：** `trees` 为 1 到 100 个 Newick 字符串或 Newick 文件路径组成的列表，叶标签必须唯一且彼此相同；`output_path` 必须是新文件，后缀与 `image_format`（`svg`、`png`、`html` 或 `pdf`）一致；`width`、`height` 为 null 或 50 到 10,000 的整数；`max_tips` 与 `max_file_bytes` 限制每棵树的规模。
- **输出：** 字典，含 `status`、`operation`、`result` 和 `provenance`。`result` 包含 `image_path`、`image_format`、`image_bytes`、`tree_count`、`tip_count`、`tip_labels`（已排序）、`newick`（每棵输入树一个字符串）和 `consensus_newick`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.phylogeny import draw_tree_cloud_figure

trees = [
    "((A:1,B:1):1,(C:2,D:2):1);",
    "((A:1,C:1):1,(B:2,D:2):1);",
    "((A:1,B:1):1,(C:1,D:1):1);",
]

with tempfile.TemporaryDirectory() as tmp:
    out = Path(tmp) / "cloud.svg"
    result = draw_tree_cloud_figure(trees, out, verbose=False)["result"]
    print(result["image_format"], result["tree_count"], result["tip_count"], result["tip_labels"])
    print(result["newick"][1])
    print(result["consensus_newick"])
    print(out.read_text(encoding="utf-8").startswith("<svg"))
```

- **示例输出：**

```text
svg 3 4 ['A', 'B', 'C', 'D']
((A:1,C:1):1,(B:2,D:2):1);
((A:2,B:2.66666666667)0.666666666667:1,C:2.66666666667,D:3.33333333333);
True
```

## 20) `dnakit.phylogeny.draw_tree_figure()`

- **作用：** 用 toytree 把一棵树画成 SVG、PNG、HTML 或 PDF 文件，并返回绘图所用的数据。
- **处理方法：** 用 `toytree.tree()` 读取 Newick 文本（或文件），按所选布局、样式和边类型调用 `ToyTree.draw()`，再把请求的 `toytree.annotate` 或 `toytree.utils` 注释加到坐标轴上，最后用 `toytree.save()` 保存画布。注释中含 `tip_labels` 时，叶标签由 `add_tip_labels()` 添加（按 `feature` 或叶名标注），而不是由 `tree.draw` 绘制。节点表给出每个节点的布局坐标 `x` 和 `y`。
- **API：** `draw_tree_figure(tree: 'object', output_path: 'str | Path', *, image_format: 'str' = 'svg', layout: 'str' = 'r', tree_style: 'str' = 'n', edge_type: 'str | None' = None, width: 'int | None' = None, height: 'int | None' = None, tip_labels_align: 'bool' = False, annotations: 'Sequence[str]' = (), feature: 'str | None' = None, pie_proportions: 'Sequence[Sequence[float]] | None' = None, max_tips: 'int' = 10000, max_file_bytes: 'int' = 10000000, verbose: bool = True) -> 'dict[str, object]'`。
- **运行要求：** 需要可选包 toytree（extra `dnakit[bio-toytree]`，会带入 toyplot、numpy 和 pandas），不需要外部程序、网络或模型。
- **输入：** `tree` 为 Newick 文本（必须含 `;`）或单棵树的 Newick 文件路径；`output_path` 为新文件，后缀与 `image_format`（`svg`、`png`、`html` 或 `pdf`）一致；`layout` 取 `r`、`l`、`u`、`d`、`c`、`c0-180`、`un` 之一；`tree_style` 取 `n`、`r`、`c`、`s`、`o`、`b` 之一；`edge_type` 为 null、`p`、`c` 或 `b`；`annotations` 为不重复的名称列表，可选 `tip_labels`、`tip_markers`、`node_markers`、`edge_markers`、`node_labels`、`edge_labels`、`node_bars`、`scale_bar`、`box_outline`、`box_outline_spines`、`external_ticks`、`node_pie_charts`、`edge_pie_charts`（`node_bars` 需要 `feature`；两种饼图需要 `pie_proportions`，每个节点索引一行、每行之和为 1）；`width`、`height` 为 null 或 50 到 10,000 的整数。
- **输出：** 字典，含 `status`、`operation`、`result` 和 `provenance`。`result` 包含 `image_path`、`image_format`、`image_bytes`、`layout`、`tree_style`、`edge_type`、`annotations`、`tip_count`、`tip_labels`（已排序）、`newick`、`nodes`（每项含 `idx`、`name`、`is_tip`、`height`、`dist`、`support`、`x`、`y`）和 `tip_label_angles`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.phylogeny import draw_tree_figure

tree = "((A:1,B:1)n1:1,(C:2,D:2)n2:1)root;"

with tempfile.TemporaryDirectory() as tmp:
    out = Path(tmp) / "tree.svg"
    result = draw_tree_figure(
        tree, out, annotations=["tip_labels", "node_markers", "scale_bar"], verbose=False
    )["result"]
    print(result["image_format"], result["layout"], result["annotations"])
    print(result["tip_count"], result["tip_labels"], result["newick"])
    for node in result["nodes"][:3]:
        print(node["idx"], node["name"], node["is_tip"], node["x"], node["y"])
    print(out.read_text(encoding="utf-8").startswith("<svg"))
```

- **示例输出：**

```text
svg r ['tip_labels', 'node_markers', 'scale_bar']
4 ['A', 'B', 'C', 'D'] ((A:1,B:1):1,(C:2,D:2):1);
0 A True -1.0 0.0
1 B True -1.0 1.0
2 C True 0.0 2.0
True
```

## 21) `dnakit.phylogeny.draw_tree_grid_figure()`

- **作用：** 把多棵树按行列网格并排绘制并保存为图像。
- **处理方法：** 用 `toytree.tree()` 解析每棵树，检查叶标签完全相同，构建 `toytree.mtree(trees)`，用 `MultiTree.draw(shape=(rows, columns))` 绘制（每棵树占一个格子），再用 `toytree.save()` 保存画布。
- **API：** `draw_tree_grid_figure(trees: 'Sequence[str]', output_path: 'str | Path', *, image_format: 'str' = 'svg', rows: 'int' = 1, columns: 'int' = 4, width: 'int | None' = None, height: 'int | None' = None, max_tips: 'int' = 10000, max_file_bytes: 'int' = 10000000, verbose: bool = True) -> 'dict[str, object]'`。
- **运行要求：** 需要可选包 toytree（extra `dnakit[bio-toytree]`，会带入 toyplot、numpy 和 pandas），不需要外部程序、网络或模型。
- **输入：** `trees` 为 1 到 100 个 Newick 字符串或文件路径，叶标签必须唯一且彼此相同，数量不超过 `rows * columns`；`rows`、`columns` 为 1 到 20 的整数（默认 1 和 4）；`output_path`、`image_format`、`width`、`height`、`max_tips` 与 `max_file_bytes` 同 `draw_tree_cloud_figure()`。
- **输出：** 字典，含 `status`、`operation`、`result` 和 `provenance`。`result` 包含 `image_path`、`image_format`、`image_bytes`、`shape`（`[rows, columns]`）、`tree_count` 和 `trees`（每棵树的 Newick 字符串和已排序的叶标签）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.phylogeny import draw_tree_grid_figure

trees = [
    "((A:1,B:1):1,(C:2,D:2):1);",
    "((A:1,C:1):1,(B:2,D:2):1);",
]

with tempfile.TemporaryDirectory() as tmp:
    out = Path(tmp) / "grid.svg"
    result = draw_tree_grid_figure(trees, out, rows=1, columns=2, verbose=False)["result"]
    print(result["shape"], result["tree_count"])
    for tree in result["trees"]:
        print(tree["newick"], tree["tip_labels"])
    print(out.read_text(encoding="utf-8").startswith("<svg"))
```

- **示例输出：**

```text
[1, 2] 2
((A:1,B:1):1,(C:2,D:2):1); ['A', 'B', 'C', 'D']
((A:1,C:1):1,(B:2,D:2):1); ['A', 'B', 'C', 'D']
True
```

## 22) `dnakit.phylogeny.draw_tree_node_feature_figure()`

- **作用：** 按某个数值型节点特征给节点着色或设定大小来绘制树，并返回映射后的取值。
- **处理方法：** 用 `toytree.tree()` 读取树，从 `tree.get_node_data()` 取特征列。`encoding="color"` 时用 `toytree.style.get_color_mapped_feature()` 映射；`encoding="size"` 时用 `get_range_mapped_feature()` 映射到 `min_size` 与 `max_size` 之间（缺失值为 0，除非设置了 `nan_value`）。再用这些节点颜色或大小调用 `ToyTree.draw()` 绘制，并用 `toytree.save()` 保存。
- **API：** `draw_tree_node_feature_figure(tree: 'object', feature: 'str', output_path: 'str | Path', *, encoding: 'str' = 'color', image_format: 'str' = 'svg', layout: 'str' = 'r', tree_style: 'str' = 'n', tips_only: 'bool' = False, cmap: 'str | None' = None, domain_min: 'float | None' = None, domain_max: 'float | None' = None, nan_value: 'float | None' = None, min_size: 'float' = 5.0, max_size: 'float' = 15.0, reverse: 'bool' = False, max_tips: 'int' = 10000, max_file_bytes: 'int' = 10000000, verbose: bool = True) -> 'dict[str, object]'`。
- **运行要求：** 需要可选包 toytree（extra `dnakit[bio-toytree]`，会带入 toyplot、numpy 和 pandas），不需要外部程序、网络或模型。
- **输入：** `tree` 为 Newick 文本或单棵树的文件路径；`feature` 必须是节点数据中的列（如 `dist`、`height` 或 `support`），否则抛出 `ConfigurationError`；`output_path` 为新文件，后缀与 `image_format` 一致；`encoding` 为 `color` 或 `size`；`layout`、`tree_style` 同 `draw_tree_figure()`；着色时可用 `cmap`（toyplot 的色板名，名称未知会抛出 `ConfigurationError`）、`domain_min`、`domain_max`、`nan_value`、`tips_only` 和 `reverse`；按大小时可用 `min_size` 与 `max_size`（0 到 1000，且 `max_size >= min_size`）。
- **输出：** 字典，含 `status`、`operation`、`result` 和 `provenance`。`result` 包含 `image_path`、`image_format`、`image_bytes`、`feature`、`encoding`、`tips_only`、`tip_labels`、`newick`、`nodes`（同 `draw_tree_figure()`）、`feature_values`（每个节点索引一个值）和 `mapped`（十六进制颜色或大小，每个节点索引一个）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.phylogeny import draw_tree_node_feature_figure

tree = "((A:1,B:1)n1:1,(C:2,D:2)n2:1)root;"

with tempfile.TemporaryDirectory() as tmp:
    colored = draw_tree_node_feature_figure(tree, "height", Path(tmp) / "color.svg", verbose=False)["result"]
    sized = draw_tree_node_feature_figure(
        tree, "height", Path(tmp) / "size.svg", encoding="size", verbose=False
    )["result"]
    print(colored["feature_values"])
    print(colored["mapped"])
    print([round(value, 2) for value in sized["mapped"]])
```

- **示例输出：**

```text
[1.0, 1.0, 0.0, 0.0, 2.0, 2.0, 3.0]
['#bfe5a0', '#bfe5a0', '#5e4fa2', '#5e4fa2', '#fdbf6f', '#fdbf6f', '#9e0142']
[8.33, 8.33, 5.0, 5.0, 11.67, 11.67, 15.0]
```

## 23) `dnakit.phylogeny.export_tree_node_data()`

- **作用：** 返回树的节点数据表，并可选地写成分隔符文本文件。
- **处理方法：** 用 `toytree.tree()` 读取树，取 `tree.get_node_data(features)`；`tips_only` 为 true 时取 `tree.get_tip_data(features)`。指定的 `nodes` 按给定顺序返回，不在树中的名称会抛出 `ConfigurationError`。给出 `output_path` 时用 `to_csv(sep=separator, index=False)` 写出，`idx` 列在最前，缺失值为空。这对应 toytree 的 `get-node-data` 命令。这里可以取到内部节点名，尽管 `tree.write()` 不输出它们。
- **API：** `export_tree_node_data(tree: 'object', *, features: 'Sequence[str] | None' = None, nodes: 'Sequence[str] | None' = None, tips_only: 'bool' = False, separator: 'str' = '\t', output_path: 'str | Path | None' = None, max_tips: 'int' = 10000, max_file_bytes: 'int' = 10000000, verbose: bool = True) -> 'dict[str, object]'`。
- **运行要求：** 需要可选包 toytree（extra `dnakit[bio-toytree]`，会带入 toyplot、numpy 和 pandas），不需要外部程序、网络或模型。
- **输入：** `tree` 为 Newick 文本或单棵树的文件路径；`features` 为 null（全部列：`idx`、`name`、`height`、`dist`、`support` 以及其他特征）或列名列表；`nodes` 为 null 或互不相同的节点名列表；`tips_only` 为布尔值；`separator` 为单个字符（默认制表符）；`output_path` 为 null 或新文件。
- **输出：** 字典，含 `status`、`operation`、`result` 和 `provenance`。`result` 包含 `tips_only`、`columns`、`row_count`、`rows`（字典列表，缺失数值为 null）、`output_path`、`output_bytes` 和 `separator`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.phylogeny import export_tree_node_data

tree = "((A:1,B:1)n1:1,(C:2,D:2)n2:1)root;"

with tempfile.TemporaryDirectory() as tmp:
    out = Path(tmp) / "nodes.csv"
    result = export_tree_node_data(
        tree, features=["name", "height"], nodes=["n2", "A"], separator=",", output_path=out, verbose=False
    )["result"]
    print(result["columns"], result["row_count"])
    for row in result["rows"]:
        print(row)
    print(out.read_text(encoding="utf-8"), end="")
    tips = export_tree_node_data(tree, tips_only=True, verbose=False)["result"]
    print(tips["row_count"], [row["name"] for row in tips["rows"]])
```

- **示例输出：**

```text
['idx', 'name', 'height'] 2
{'idx': 5, 'name': 'n2', 'height': 2.0}
{'idx': 0, 'name': 'A', 'height': 1.0}
idx,name,height
5,n2,2.0
0,A,1.0
4 ['A', 'B', 'C', 'D']
```

## 24) `dnakit.phylogeny.map_node_feature_values()`

- **作用：** 用 toytree 的取值映射器把一组数字映射为颜色或大小，不绘图。
- **处理方法：** 先把 `None` 转为 NaN，再调用 `toytree.style.get_color_mapped_values()`（十六进制颜色，`cmap` 为 null 时用 toytree 默认色板）或 `get_range_mapped_values()`（大小在 `min_size` 到 `max_size` 之间）。NaN 用 `nan_value` 填充；按大小映射时默认填 0.0，着色时缺失值若未设置 `nan_value` 则得到 `#000000`。
- **API：** `map_node_feature_values(values: 'Sequence[float | None]', *, encoding: 'str' = 'color', cmap: 'str | None' = None, domain_min: 'float | None' = None, domain_max: 'float | None' = None, nan_value: 'float | None' = None, min_size: 'float' = 5.0, max_size: 'float' = 15.0, reverse: 'bool' = False, verbose: bool = True) -> 'dict[str, object]'`。
- **运行要求：** 需要可选包 toytree（extra `dnakit[bio-toytree]`，会带入 toyplot、numpy 和 pandas），不需要外部程序、网络或模型。
- **输入：** `values` 为非空列表，最多 1,000,000 项，元素为数字或 null；`encoding` 为 `color` 或 `size`；`cmap` 为 null 或 toyplot 色板名（名称未知会抛出 `ConfigurationError`）；`domain_min`、`domain_max` 和 `nan_value` 为 null 或有限数；`min_size`、`max_size` 为 0 到 1000 的数，且 `max_size >= min_size`；`reverse` 为布尔值。
- **输出：** 字典，含 `status`、`operation`、`result` 和 `provenance`。`result` 包含 `encoding`、`count` 和 `mapped`（十六进制颜色字符串或浮点数，顺序与输入一致）。
- **示例代码：**

```python
from dnakit.phylogeny import map_node_feature_values

values = [0, 1, None, 2]
colors = map_node_feature_values(values, verbose=False)["result"]
filled = map_node_feature_values(values, nan_value=1, domain_min=0, domain_max=2, verbose=False)["result"]
sizes = map_node_feature_values(values, encoding="size", verbose=False)["result"]
print(colors["count"], colors["mapped"])
print(filled["mapped"])
print(sizes["encoding"], sizes["mapped"])
```

- **示例输出：**

```text
4 ['#5e4fa2', '#ffffbf', '#000000', '#9e0142']
['#5e4fa2', '#ffffbf', '#ffffbf', '#9e0142']
size [5.0, 10.0, 0.0, 15.0]
```

## 25) `dnakit.phylogeny.prune_tree_tips()`

- **作用：** 保留连接指定节点的最小子树，并以 Newick 返回。
- **处理方法：** 用 `toytree.tree()` 读取树，检查每个名称都是树中的节点名，然后调用 `toytree.mod.prune(tree, *names, preserve_dists=..., require_root=...)`。其余叶被删除，变成单子节点的节点也被去掉。这对应 toytree 的 `prune` 命令。
- **API：** `prune_tree_tips(tree: 'object', names: 'Sequence[str]', *, preserve_dists: 'bool' = True, require_root: 'bool' = False, max_tips: 'int' = 10000, max_file_bytes: 'int' = 10000000, verbose: bool = True) -> 'dict[str, object]'`。
- **运行要求：** 需要可选包 toytree（extra `dnakit[bio-toytree]`，会带入 toyplot、numpy 和 pandas），不需要外部程序、网络或模型。
- **输入：** `tree` 为 Newick 文本或单棵树的文件路径；`names` 为至少两个互不相同的叶名或有名内部节点名；`preserve_dists`（默认 true）保留保留节点之间的路径长度，false 则各保留的边沿用自己的长度；`require_root` 让原树的根作为结果中的一个节点保留下来。
- **输出：** 字典，含 `status`、`operation`、`result` 和 `provenance`。`result` 包含 `newick`、`tip_count`、`tips`、`retained_names`、`preserve_dists` 和 `require_root`。
- **示例代码：**

```python
from dnakit.phylogeny import prune_tree_tips

tree = "((A:1,B:1)n1:1,(C:2,D:2)n2:1)root;"
merged = prune_tree_tips(tree, ["A", "C"], verbose=False)["result"]
own = prune_tree_tips(tree, ["A", "C"], preserve_dists=False, verbose=False)["result"]
print(merged["newick"], merged["tip_count"], merged["retained_names"])
print(own["newick"], own["preserve_dists"])
```

- **示例输出：**

```text
(A:2,C:3); 2 ['A', 'C']
(A:1,C:2); False
```

## 26) `dnakit.phylogeny.root_toytree_on_outgroup()`

- **作用：** 按节点名指定外类群，在其上方的边上给树定根，并以 Newick 返回。
- **处理方法：** 用 `toytree.tree()` 读取树，检查每个外类群名称都是节点名，然后调用 `toytree.mod.root(tree, *outgroup, edge_features=...)`。给一个名称时根放在该节点的边上；给多个名称时根放在它们最近共同祖先上方的边上。源码模块中这个函数叫 `root_tree_on_outgroup`，但导出名是 `root_toytree_on_outgroup`，因为 `dnakit.phylogeny` 已有另一个会写文件的同名 `root_tree_on_outgroup`。
- **API：** `root_toytree_on_outgroup(tree: 'object', outgroup: 'Sequence[str]', *, edge_features: 'Sequence[str] | None' = None, max_tips: 'int' = 10000, max_file_bytes: 'int' = 10000000, verbose: bool = True) -> 'dict[str, object]'`。
- **运行要求：** 需要可选包 toytree（extra `dnakit[bio-toytree]`，会带入 toyplot、numpy 和 pandas），不需要外部程序、网络或模型。
- **输入：** `tree` 为 Newick 文本或单棵树的文件路径；`outgroup` 为一个或多个互不相同的叶名或有名内部节点名组成的列表；`edge_features` 为 null 或节点特征名列表，toytree 会把这些特征在新根的两条边之间分配。
- **输出：** 字典，含 `status`、`operation`、`result` 和 `provenance`。`result` 包含 `newick`（定根后的树）、`tip_count`、`tips`、`outgroup` 和 `edge_features`。
- **示例代码：**

```python
from dnakit.phylogeny import root_toytree_on_outgroup

tree = "((A:1,B:1)n1:1,(C:2,D:2)n2:1)root;"
single = root_toytree_on_outgroup(tree, ["D"], verbose=False)["result"]
pair = root_toytree_on_outgroup(tree, ["A", "B"], verbose=False)["result"]
print(single["newick"], single["outgroup"])
print(pair["newick"], pair["tips"])
```

- **示例输出：**

```text
(D:1,(C:2,(A:1,B:1):2):1); ['D']
((A:1,B:1):1,(C:2,D:2):1); ['A', 'B', 'C', 'D']
```
