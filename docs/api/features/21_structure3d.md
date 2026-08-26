# 三维结构与力学性质

读取 DNA 三维结构文件并分析坐标几何、NMR 多模型柔性以及 3DNA/DSSR 结构参数。

三维性质必须来自显式坐标或轨迹，不会从普通序列伪造唯一三维结构。下列示例从项目根目录运行，结果均由当前本地源码和已下载样例实际计算得到。

## 1) 读取单个 PDB 模型

- **作用：** 读取 PDB 中指定的一个 DNA 模型，提取原子坐标、残基、链和推导序列，生成后续三维几何计算使用的结构对象。
- **API：** `dnakit.structure3d.load_pdb(path[必须], model_index[可选])`
- **输入：** 必填可读的 legacy PDB 路径；可选正整数模型编号。
- **示例代码：**

```python
from pathlib import Path

from dnakit.structure3d import load_pdb

structure = load_pdb(Path("temp/dna_structures/1BNA.pdb"))
print(structure.pdb_id, structure.model_index)
print(len(structure.atoms), len(structure.residues), structure.chain_ids)
print(structure.sequence_by_chain)
```

- **示例结果：**

```text
1BNA 1
486 24 ('A', 'B')
('CGCGAATTCGCG', 'CGCGAATTCGCG')
```

- **限制：** 只解析 ASCII legacy PDB 和 canonical DNA `ATOM` 记录，不读取 mmCIF、蛋白质或非标准残基。当前上限为 100 MB、1000000 个 DNA 原子和 10000 个模型。

## 2) 读取 PDB 多模型集合

- **作用：** 从同一 PDB 或模型集合读取多个 DNA 构象，统一原子对应关系，作为 RMSF 和构象变化比较的输入。
- **API：** `dnakit.structure3d.load_pdb_ensemble(path[必须])`
- **输入：** 必填包含一个或多个 DNA 模型的 legacy PDB 路径。
- **示例代码：**

```python
from pathlib import Path

from dnakit.structure3d import load_pdb_ensemble

models = load_pdb_ensemble(Path("temp/dna_structures/1AC7.pdb"))
print(len(models), models[0].pdb_id)
print(len(models[0].residues), models[0].sequence_by_chain)
```

- **示例结果：**

```text
10 1AC7
16 ('ATCCTAGTTATAGGAT',)
```

- **限制：** 该 API 只读取模型，不做旋转/平移对齐、轨迹采样评估或力场校验。

## 3) 显式坐标几何分析

- **作用：** 从显式 DNA 原子坐标计算回转半径、SASA、体积、形状、骨架二面角和几何氢键等指标，用于定量描述一个三维构象。
- **API：** `dnakit.structure3d.analyze_structure(structure[必须], sasa_probe_radius_angstrom[可选], sasa_points_per_atom[可选], volume_grid_spacing_angstrom[可选], progress[可选])`
- **输入：** 必填 `DNA3DStructure`；可选 0–5 Å SASA probe、24–4096 个每原子采样点、0.25–5 Å 体素间距和进度回调。
- **示例代码：**

```python
from pathlib import Path

from dnakit.structure3d import analyze_structure, load_pdb

structure = load_pdb(Path("temp/dna_structures/1BNA.pdb"))
result = analyze_structure(
    structure,
    sasa_points_per_atom=24,
    volume_grid_spacing_angstrom=1.5,
)
assert result.helix is not None
print(result.atom_count, result.residue_count)
print(
    round(result.shape.radius_of_gyration_angstrom, 3),
    round(result.solvent_accessible_surface_area_angstrom2, 3),
    round(result.molecular_volume_angstrom3, 3),
)
print(
    round(result.helix.mean_rise_angstrom, 3),
    round(result.helix.mean_twist_degree, 3),
)
```

- **示例结果：**

```text
486 24
13.228 4650.051 4890.375
3.351 35.797
```

- **限制：** SASA 使用 Fibonacci sphere 近似，体积使用体素近似；只有显式氢原子存在时才报告真正 D–H···A 几何计数。原生螺旋参数是两条等长反向互补链的全局轴近似，不是 3DNA/DSSR 局部参考系参数。

## 4) NMR/多模型 RMSF 柔性

- **作用：** 在原子对应一致的多个三维模型间计算逐原子和逐残基 RMSF，量化不同位置在构象集合中的波动程度。
- **API：** `dnakit.structure3d.analyze_ensemble_flexibility(structures[必须])`
- **输入：** 必填 2–10000 个 `DNA3DStructure` 模型，且至少有 3 个公共原子。
- **示例代码：**

```python
from pathlib import Path

from dnakit.structure3d import analyze_ensemble_flexibility, load_pdb_ensemble

models = load_pdb_ensemble(Path("temp/dna_structures/1AC7.pdb"))
result = analyze_ensemble_flexibility(models)
print(result.model_count, result.common_atom_count)
print(
    round(result.mean_atomic_rmsf_angstrom, 3),
    round(result.max_atomic_rmsf_angstrom, 3),
)
```

- **示例结果：**

```text
10 510
0.986 2.281
```

- **限制：** 实现只做平移中心化，不做 Kabsch 旋转拟合；输入模型必须已在共同旋转参考系中。RMSF 不等同于持久长度或弯曲/扭转/伸展刚度。

## 5) 3DNA `bp_step.par` 标准参数解析

- **作用：** 解析外部 3DNA 参数文件，汇总 shift、slide、rise、twist 等碱基对步及螺旋参数，便于结构统计和模型间比较。
- **API：** `dnakit.structure3d.read_3dna_bp_step(path[必须])`
- **输入：** 必填已存在、UTF-8 可解码且符合 3DNA `bp_step.par` 列结构的文件。
- **示例代码：**

```python
from pathlib import Path
from tempfile import TemporaryDirectory

from dnakit.structure3d import read_3dna_bp_step

with TemporaryDirectory() as directory:
    path = Path(directory) / "bp_step.par"
    path.write_text(
        "2 # base-pairs\n"
        "# Shear Stretch Stagger Buckle Prop-T Opening Shift Slide Rise Tilt Roll Twist\n"
        "A-T 0.1 0.2 0.3 1 2 3 0.4 0.5 3.4 4 5 36\n"
        "T-A 0.2 0.3 0.4 2 3 4 0.5 0.6 3.5 5 6 35\n",
        encoding="utf-8",
    )
    result = read_3dna_bp_step(path)
    print(result.base_pair_count)
    print(result.mean_rise_angstrom, result.mean_twist_degree)
    print(round(result.base_pairs_per_turn or 0.0, 4))
```

- **示例结果：**

```text
2
3.5 35.0
10.2857
```

- **限制：** 示例自建的文件只用于演示解析器格式，不是 3DNA 对某个结构的计算结果。DNAKit 只解析已有标准输出，不重新实现 3DNA 局部参考框架。

## 6) DSSR JSON 摘要解析

- **作用：** 解析外部 DSSR JSON，提取碱基、配对、螺旋、茎环和氢键注释，转换为可查询、汇总和报告的结构化结果。
- **API：** `dnakit.structure3d.read_dssr_json(path[必须])`
- **输入：** 必填已存在的 DSSR JSON 输出文件。
- **示例代码：**

```python
from pathlib import Path

from dnakit.structure3d import read_dssr_json

result = read_dssr_json(Path("temp/dna_structures/1ehz-dssr-example.json"))
print(
    result.nucleotide_count,
    result.pair_count,
    result.helix_count,
    result.stem_count,
    result.hairpin_count,
)
print(result.hydrogen_bond_count, result.backbone_torsion_record_count)
```

- **示例结果：**

```text
76 34 2 4 3
116 76
```

- **限制：** `1ehz-dssr-example.json` 是 DSSR 官方 RNA JSON，只用于验证解析器字段架构；它不是 DNA 结果，也不是由本项目的 DNA PDB 样例生成的输出。该 API 只读摘要计数，不执行 DSSR。

<span id="7"></span>**当前条件功能**

以下能力尚无独立的 native 公开计算 API，因此不伪造示例数值：

| 功能 | 所需附加输入/后端 |
| --- | --- |
| 大沟/小沟宽度与深度 | 3DNA/DSSR 或 Curves+ 的局部碱基对框架 |
| 碱基堆积面积 | 标准碱基平面和原子分类 |
| 静电势与电荷分布 | PQR/力场电荷及 APBS 等求解器 |
| 持久长度、弯曲/扭转/伸展刚度 | 已对齐的轨迹/构象集合及明确力学模型 |
| 局部柔性和变形能力 | 足够采样的实验/MD 构象集合 |

<span id="8"></span>**已下载并实测的结构**

| 文件 | 类型 | 模型数 | 当前解析结果摘要 |
| --- | --- | ---: | --- |
| `1BNA.pdb` | B-DNA 晶体 12-mer 双链 | 1 | 24 residues；Rg 13.228 Å；近似 rise 3.351 Å、twist 35.797° |
| `1AC7.pdb` | 16-mer DNA hairpin NMR | 10 | 16 residues；平移中心化 mean atomic RMSF 0.986 Å |
| `139D.pdb` | 平行 DNA G-quadruplex NMR | 4 | 28 residues；4 chains；mean atomic RMSF 1.138 Å |

来源记录分别为 [RCSB 1BNA](https://www.rcsb.org/structure/1BNA)、[RCSB 1AC7](https://www.rcsb.org/structure/1AC7) 和 [RCSB 139D](https://www.rcsb.org/structure/139D)。下载 URL、大小和 SHA-256 在 `temp/dna_structures/manifest.json`；完整本地运行结果在 `temp/dna_structures/analysis_results.json`。

一键重算会显示每个结构的进度条：

```bash
PYTHONPATH=src python examples/analyze_dna_structures.py \
  --input-dir temp/dna_structures \
  --output temp/dna_structures/analysis_results.json \
  --sasa-points 96 \
  --volume-grid 0.75
```

3DNA 参数定义可查阅 [3DNA 官方示意](https://x3dna.org/highlights/schematic-diagrams-of-base-pair-parameters)，JSON 结构见 [DSSR JSON 文档](https://x3dna.org/highlights/dssr-output-in-json-format)。
