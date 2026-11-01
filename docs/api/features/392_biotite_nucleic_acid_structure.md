# Biotite 序列与核酸结构操作

Biotite >=2,<3. 示例使用源码仓库的本地测试文件，不联网。

## base_pairs_to_dot_bracket

**作用：** Write base pairs as dot-bracket notation (``biotite.structure.dot_bracket``).

**处理方法：** Write base pairs as dot-bracket notation (``biotite.structure.dot_bracket``).

``base_pairs`` are 0-based position pairs ``[i, j]`` of a strand of ``length`` bases.
Pseudoknots use bracket types ``()``, ``[]``, ``{}``, ``<>`` and letters. ``scores`` rank
pairs when the pseudoknot assignment is ambiguous; ``max_pseudoknot_order`` caps the
pseudoknot order. More than one equivalent notation can be returned.

**API：** `dnakit.secondary_structure.biotite_p30_structure_functions.base_pairs_to_dot_bracket(base_pairs: 'Sequence[Sequence[int]]', length: 'int', *, scores: 'Sequence[int] | None' = None, max_pseudoknot_order: 'int | None' = None, verbose: bool = True) -> 'dict[str, Any]'`

**输入：** 参数类型与默认值见 API；序列字符串、对齐行和结构文件均在调用时验证。

**输出：** 可序列化字典，包含状态、实际计算结果或图片路径及来源信息。

**示例代码：**

```python
from dnakit.secondary_structure.biotite_p30_structure_functions import base_pairs_to_dot_bracket
result = base_pairs_to_dot_bracket(base_pairs=[[0, 10], [1, 9], [2, 8], [4, 14], [5, 13]], length=15, verbose=False)
print(result["status"])
```

**示例输出：**

```text
completed
```

## nucleic_acid_base_pairs

**作用：** Find the base pairs of a nucleic-acid 3D structure (``biotite.structure.base_pairs``).

**处理方法：** Find the base pairs of a nucleic-acid 3D structure (``biotite.structure.base_pairs``).

Pairs follow the DSSR criteria. ``annotate_edges=True`` adds the Leontis-Westhof edge of
each base (``base_pairs_edge``: watson_crick, hoogsteen, sugar); ``annotate_glycosidic_bond``
adds the cis/trans orientation (``base_pairs_glycosidic_bond``). The edge annotation needs
hydrogen atoms in the structure (for example NMR models); without them Biotite raises an
error that is reported as a configuration error.

**API：** `dnakit.secondary_structure.biotite_p30_structure_functions.nucleic_acid_base_pairs(structure_file: 'str', *, model: 'int' = 1, chain_id: 'str | None' = None, min_atoms_per_base: 'int' = 3, unique: 'bool' = True, annotate_edges: 'bool' = False, annotate_glycosidic_bond: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`

**输入：** 参数类型与默认值见 API；序列字符串、对齐行和结构文件均在调用时验证。

**输出：** 可序列化字典，包含状态、实际计算结果或图片路径及来源信息。

**示例代码：**

```python
from pathlib import Path
import dnakit
from dnakit.secondary_structure.biotite_p30_structure_functions import nucleic_acid_base_pairs
result = nucleic_acid_base_pairs(structure_file=str(Path(dnakit.__file__).resolve().parents[2] / 'others/tests/data/biotite_p30/rna_hairpin.pdb'), verbose=False)
print(result["status"])
```

**示例输出：**

```text
completed
```

## nucleic_acid_base_stacking

**作用：** Find stacked neighbouring bases of a nucleic-acid 3D structure (``base_stacking``).

**处理方法：** Find stacked neighbouring bases of a nucleic-acid 3D structure (``base_stacking``).

**API：** `dnakit.secondary_structure.biotite_p30_structure_functions.nucleic_acid_base_stacking(structure_file: 'str', *, model: 'int' = 1, chain_id: 'str | None' = None, min_atoms_per_base: 'int' = 3, verbose: bool = True) -> 'dict[str, Any]'`

**输入：** 参数类型与默认值见 API；序列字符串、对齐行和结构文件均在调用时验证。

**输出：** 可序列化字典，包含状态、实际计算结果或图片路径及来源信息。

**示例代码：**

```python
from pathlib import Path
import dnakit
from dnakit.secondary_structure.biotite_p30_structure_functions import nucleic_acid_base_stacking
result = nucleic_acid_base_stacking(structure_file=str(Path(dnakit.__file__).resolve().parents[2] / 'others/tests/data/biotite_p30/rna_hairpin.pdb'), verbose=False)
print(result["status"])
```

**示例输出：**

```text
completed
```

## nucleic_acid_dot_bracket

**作用：** Dot-bracket secondary structure of a nucleic-acid 3D structure.

**处理方法：** Dot-bracket secondary structure of a nucleic-acid 3D structure.

Finds the base pairs with ``base_pairs`` and writes them with
``dot_bracket_from_structure``. Chains of a duplex are concatenated in file order;
``chain_id`` restricts the strand. Returns the residue name of every nucleotide used.

**API：** `dnakit.secondary_structure.biotite_p30_structure_functions.nucleic_acid_dot_bracket(structure_file: 'str', *, model: 'int' = 1, chain_id: 'str | None' = None, scores: 'Sequence[int] | None' = None, max_pseudoknot_order: 'int | None' = None, verbose: bool = True) -> 'dict[str, Any]'`

**输入：** 参数类型与默认值见 API；序列字符串、对齐行和结构文件均在调用时验证。

**输出：** 可序列化字典，包含状态、实际计算结果或图片路径及来源信息。

**示例代码：**

```python
from pathlib import Path
import dnakit
from dnakit.secondary_structure.biotite_p30_structure_functions import nucleic_acid_dot_bracket
result = nucleic_acid_dot_bracket(structure_file=str(Path(dnakit.__file__).resolve().parents[2] / 'others/tests/data/biotite_p30/rna_hairpin.pdb'), verbose=False)
print(result["status"])
```

**示例输出：**

```text
completed
```

## plot_nucleotide_secondary_structure

**作用：** Draw a nucleic-acid secondary structure to an image file.

**处理方法：** Draw a nucleic-acid secondary structure to an image file.

Uses ``plot_nucleotide_secondary_structure``: the 2D layout comes from ViennaRNA's
``RNAplot`` (``layout`` is ``radial``, ``naview``, ``circular``, ``rnaturtle`` or
``rnapuzzler``). ``bases`` is the sequence and ``dot_bracket`` the structure of the same
length. Returns the image path, the base pairs and the layout.

**API：** `dnakit.secondary_structure.biotite_p30_structure_functions.plot_nucleotide_secondary_structure(bases: 'str', dot_bracket: 'str', output_path: 'str', *, layout: 'str' = 'naview', draw_pseudoknots: 'bool' = True, angle: 'float' = 0.0, rnaplot_path: 'str' = 'RNAplot', dpi: 'int' = 150, verbose: bool = True) -> 'dict[str, Any]'`

**输入：** 参数类型与默认值见 API；序列字符串、对齐行和结构文件均在调用时验证。

**输出：** 可序列化字典，包含状态、实际计算结果或图片路径及来源信息。

此绘图接口需要外部 RNAplot；当前环境未提供该程序，因此未编造成功输出。

## pseudoknot_order

**作用：** Assign each base pair its pseudoknot order (``biotite.structure.pseudoknots``).

**处理方法：** Assign each base pair its pseudoknot order (``biotite.structure.pseudoknots``).

Order 0 pairs are nested; order 1 pairs cross them, and so on. Returns one assignment per
optimal solution (``orders[solution][pair]``); -1 marks pairs above ``max_pseudoknot_order``.
``scores`` give an integer weight per pair for choosing among equal solutions.

**API：** `dnakit.secondary_structure.biotite_p30_structure_functions.pseudoknot_order(base_pairs: 'Sequence[Sequence[int]]', *, scores: 'Sequence[int] | None' = None, max_pseudoknot_order: 'int | None' = None, verbose: bool = True) -> 'dict[str, Any]'`

**输入：** 参数类型与默认值见 API；序列字符串、对齐行和结构文件均在调用时验证。

**输出：** 可序列化字典，包含状态、实际计算结果或图片路径及来源信息。

**示例代码：**

```python
from dnakit.secondary_structure.biotite_p30_structure_functions import pseudoknot_order
result = pseudoknot_order(base_pairs=[[0, 10], [1, 9], [2, 8], [4, 14], [5, 13]], verbose=False)
print(result["status"])
```

**示例输出：**

```text
completed
```

