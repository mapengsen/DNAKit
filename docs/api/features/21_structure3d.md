# 三维结构与力学性质

读取 DNA 三维结构文件并分析坐标几何、NMR 多模型柔性以及 3DNA/DSSR 结构参数。

三维性质必须来自显式坐标或轨迹，不会从普通序列伪造唯一三维结构。下列示例从项目根目录运行，结果均由当前本地源码和已下载样例实际计算得到。

各项计算的论文来源和内部公式见 [FAQ：三维结构与力学性质的计算依据和参考文献](../../faq.md#structure3d-references)。

## 1) 读取单个 PDB 模型

- **作用：** 读取 PDB 中指定的一个 DNA 模型，提取原子坐标、残基、链和推导序列，生成后续三维几何计算使用的结构对象。
- **计算方法：** 按 legacy PDB 固定列规格解析 `HEADER`、`TITLE`、`MODEL` 和 `ATOM`，只保留 canonical DNA 残基；备选构象只接受空白或 `A`，重复原子按主构象优先、再按 occupancy 选择，最后按残基顺序映射 A/C/G/T 序列。这是文件解析，不是结构预测。
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

## 2) 读取 PDB 多模型集合

- **作用：** 从同一 PDB 或模型集合读取多个 DNA 构象，统一原子对应关系，作为 RMSF 和构象变化比较的输入。
- **计算方法：** 使用与单模型相同的固定列解析规则，按 `MODEL` 编号分组并返回全部模型。此步只建立链/残基/原子标识，不做旋转、平移对齐或 RMSF 计算。
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

## 3) 显式坐标几何分析

- **作用：** 从显式 DNA 原子坐标计算回转半径、SASA、体积、形状、骨架二面角和几何氢键等指标，用于定量描述一个三维构象。
- **计算方法：** 由原子质量加权坐标计算质心、回转半径和张量特征值；SASA 用 `vdW 半径 + probe` 球面的确定性采样点暴露比例求和，体积用 vdW 球占据的体素并集估算；骨架角由四原子二面角计算，氢键用明确的距离/角度阈值筛选。双链螺旋指标是配对 `C1′` 中心和全局轴的 DNAKit 近似值。
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

## 4) NMR/多模型 RMSF 柔性

- **作用：** 在原子对应一致的多个三维模型间计算逐原子和逐残基 RMSF，量化不同位置在构象集合中的波动程度。
- **计算方法：** 先取所有模型共有的原子键，对每个模型减去共有原子的几何中心以去除平移；对原子 `i` 计算 `RMSF_i = sqrt[Σm |r_im − <r_i>|²/M]`，残基 RMSF 再对其原子 RMSF 做均方根聚合。实现不做 Kabsch 旋转拟合。
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

## 5) 3DNA `bp_step.par` 标准参数解析

- **作用：** 解析外部 3DNA 参数文件，汇总 shift、slide、rise、twist 等碱基对步及螺旋参数，便于结构统计和模型间比较。
- **计算方法：** 按 3DNA `bp_step.par` 的 12 个刚体参数列解析每行，校验行数与声明碱基对数；对第二行起的 step 求平均 rise/twist，再计算 `bp/turn = 360/平均twist` 和 `pitch = 平均rise × bp/turn`。DNAKit 不重新计算 3DNA 局部参考框架。
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

## 6) DSSR JSON 摘要解析

- **作用：** 解析外部 DSSR JSON，提取碱基、配对、螺旋、茎环和氢键注释，转换为可查询、汇总和报告的结构化结果。
- **计算方法：** 读取 DSSR JSON 中的程序版本以及 `nts`、`pairs`、`helices`、`stems`、`hairpins`、`hbonds` 和 `ntParams/ntPars` 字段，优先校验声明计数与数组长度后返回摘要。结构识别由生成该 JSON 的外部 DSSR 完成，DNAKit 不执行 DSSR 算法。
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

<span id="7"></span>
