# FAQ

## 序列功能搜索采用什么匹配策略？ {#pattern-matching-strategy}

[序列功能搜索](api/features/06_patterns.md)中的 10 类功能都是**确定性的规则算法**：在算法版本、输入序列和参数相同的情况下，结果一定相同。它们不使用机器学习模型、随机抽样或实验数据推断。PWM 分数也是按固定公式计算的匹配分数，不是模型预测概率。

这些功能返回的是符合规则的**候选序列**。命中启动子 motif、TF motif、PAM 或倒置重复，不代表该位置一定具有真实的启动子活性、转录因子结合、CRISPR 编辑效率或发卡结构。

### 通用匹配约定

- **链方向：** 正链直接扫描；反链先生成反向互补序列再扫描；`strand="both"` 扫描两条链。
- **IUPAC 匹配：** 每个符号表示一组可能的碱基，例如 `N=ACGT`、`R=AG`、`Y=CT`、`W=AT`。目标符号和规则符号表示的碱基集合有交集时，该位置兼容；所有位置都兼容时，整个窗口才命中。
- **Gap：** 显式 Gap 是扫描边界，候选不会跨越 Gap。
- **环状序列：** exact/IUPAC motif、PWM、启动子 motif、限制酶位点和 PAM 可以在支持的条件下跨越环状原点；正则 motif、密码子和重复结构扫描不能跨越环状原点。
- **资源上限：** `max_matches` 限制返回数量，达到上限时结果会标记为截断；比较次数或扫描单元超过对应上限时会报错，不会静默省略计算。

### 1. 通用 motif 与 PWM

- **exact：** 用与 motif 等长的窗口逐位扫描，窗口字符串与 motif 完全相同才命中。
- **IUPAC：** 同样逐窗扫描，但每个位置按上面的 IUPAC 集合兼容规则判断。
- **regex：** 使用受限的 DNA 正则表达式查找；是否允许重叠命中由 `overlapping` 控制。
- **PWM：** 对每个与矩阵等长的窗口计算固定的 log2-odds 分数：`各位置 log2(该碱基的 PWM 概率 / 该碱基的背景概率)之和`。分数大于或等于 `threshold` 才命中；背景默认 A、C、G、T 各为 0.25，包含模糊碱基的目标窗口会跳过。

motif、正则、PWM 和阈值由调用方提供，不是内置的生物功能预测模型。

### 2. 起止密码子

程序分别从每条所选链的第 0、1、2 个碱基开始，每次前进 3 bp，形成三个阅读框；扫描两条链时共得到六个阅读框。每个三联体只要属于起始或终止密码子集合就命中。

- NCBI 遗传密码表 1：起始密码子为 `ATG`，终止密码子为 `TAA`、`TAG`、`TGA`。
- NCBI 遗传密码表 11：起始密码子为 `ATG`、`GTG`、`TTG`，终止密码子为 `TAA`、`TAG`、`TGA`。
- 传入 `start_codons` 或 `stop_codons` 时，使用调用方提供的集合代替对应默认集合。

该功能只报告单独的密码子位置，不要求起始密码子后面存在同阅读框的终止密码子，也不等同于完整 ORF 检测。

### 3. 启动子 motif

程序将每条启动子共识序列作为 IUPAC motif，在所选链上逐窗匹配。内置规则只有：

- 真核 TATA box：`TATAWAWR`；
- 细菌 -10 区：`TATAAT`；
- 细菌 -35 区：`TTGACA`。

传入 `motifs` 时使用调用方提供的 motif 映射代替内置映射。程序不检查 -10/-35 区之间的距离、转录起始位点或其他调控元件，因此只返回启动子 motif 候选，不预测启动子活性。

### 4. TF motif

调用方必须提供转录因子的 PWM 和阈值。程序使用与通用 PWM 相同的 log2-odds 公式逐窗打分，只返回分数大于或等于阈值的窗口。DNAKit 不内置 JASPAR 等 TF motif 数据库，也不把分数解释为真实结合强度或结合概率。

### 5. 限制酶位点

程序在正反两条链上按 IUPAC 规则匹配限制酶识别序列，命中后根据该酶定义中的上下链切割偏移量计算切点。回文识别位点的正反链重复结果会合并。

内置小型目录包含 `BamHI`、`EcoRI`、`HaeIII`、`HindIII`、`NotI` 和 `SmaI`；也可以传入自定义 `RestrictionEnzyme`。该规则不包含完整 REBASE，也不考虑甲基化对酶切的影响。

### 6. CRISPR PAM

程序先在所选链上按 IUPAC 规则逐位寻找 PAM，再根据 PAM 位于 3' 侧还是 5' 侧，从相邻位置提取固定长度的 guide：

- `SpCas9`：`NGG`，PAM 位于 guide 的 3' 侧，默认 guide 长度 20 bp；
- `SaCas9`：`NNGRRT`，PAM 位于 guide 的 3' 侧，默认 guide 长度 21 bp；
- `AsCas12a`：`TTTV`，PAM 位于 guide 的 5' 侧，默认 guide 长度 20 bp。

guide 超出序列边界时不命中；默认排除含模糊碱基的 guide。随后按 `min_gc`、`max_gc` 和 `exclude_motifs` 做确定性过滤。也可以提供自定义 `PAMRule`。该功能不预测编辑效率和脱靶风险。

### 7. 回文序列

程序枚举每个起点以及 `min_length` 到 `max_length` 范围内的每个长度，将候选片段与它自己的反向互补序列逐位比较；全部位置满足 IUPAC 兼容规则才命中。`maximal_per_start=True` 时，每个起点只保留最长命中。

### 8. 倒置重复

程序枚举左臂起点、臂长和 loop 长度，再计算左臂的反向互补序列；右臂与该序列逐位满足 IUPAC 兼容规则时命中。当前不允许 mismatch 或 indel，因此它只是潜在发卡结构的序列候选，不进行折叠或自由能预测。

### 9. 串联重复

程序从每个起点按单元长度从小到大枚举重复单元，检查后续相邻单元是否与第一个单元字符串完全相同，并一直延伸到首次不同的位置。达到对应的最少重复次数后，报告该起点最小的合格单元及其最长连续重复区间。

这里使用字符完全相等规则，IUPAC 模糊符号也必须是同一个字符；不允许 mismatch、indel 或中断。`overlapping=False` 时，命中后从该重复区间末尾继续扫描；设为 `True` 时从下一个碱基继续扫描。

### 10. 微卫星

微卫星调用同一套串联重复算法，但把重复单元长度固定为 1～6 bp。默认最少重复次数为：

- 1 bp 单元至少重复 6 次；
- 2～6 bp 单元至少重复 3 次。

匹配必须连续且完全相同，不允许中断、mismatch 或 indel。`min_repeats_by_unit` 可以修改阈值，但必须为 1～6 bp 的每种单元长度都提供阈值。

## Diversity 和 Novelty 有哪些计算方法与参考文献？ {#diversity-novelty-references}

DNAKit 保留原有归一化相似度方法，并新增论文中的原始 Levenshtein 距离方法：

| 指标 | 默认方法 | 第二种方法 |
| --- | --- | --- |
| Diversity | `calculation="similarity"`：距离定义为 `1 - pair_similarity`，`score` 为平均最近邻距离，并同时返回平均两两距离和阈值 cluster 摘要。 | `calculation="levenshtein"`：`Σ(i≠j) Levenshtein(xᵢ,xⱼ) / [n(n−1)]`，等价于所有无序序列对原始编辑距离的平均值。 |
| Novelty | `novelty_calculation="similarity"`：每条查询为 `1 - nearest_reference_similarity`。 | `novelty_calculation="levenshtein"`：`meanᵢ minₛ Levenshtein(queryᵢ, referenceₛ)`。 |

第二种方法依据 Cherednichenko & Poptsova, *Data augmentation with generative models improves detection of Non-B DNA structures*, **Computers in Biology and Medicine** 184 (2025) 109440，[DOI 10.1016/j.compbiomed.2024.109440](https://doi.org/10.1016/j.compbiomed.2024.109440)，其中第 2.8 节公式 (20) 和 (21) 将距离说明为 Levenshtein 距离。该文沿用了 Jain et al., *Biological Sequence Design with GFlowNets*, ICML 2022，[PMLR 论文页](https://proceedings.mlr.press/v162/jain22a.html)中的术语。

文章有[官方 GitHub 仓库](https://github.com/powidla/nonB-DNA-structures-generation)，相关代码位于 [`seq_analysis.ipynb`](https://github.com/powidla/nonB-DNA-structures-generation/blob/ea61a37f95c5a1effe64324af366c781755fe4c8/notebooks/seq_analysis.ipynb)。截至 2026-09-02，该仓库没有声明开源许可证，并且 notebook 使用固定 100 bp、展平 one-hot/KDTree 及分块计算，与正文公式并不完全一致。因此 DNAKit 没有复制该代码，而是用自身有界 Levenshtein 实现正文公式；不会隐式填充或截断序列，结果也不保证复现论文 Table 2 的 notebook 数值。

Levenshtein 结果单位为“编辑操作数”，没有归一化；数值越大表示越多样或越远离参考库。比较不同数据集时应尽量保证序列长度分布一致。

## 理化性质的计算依据和参考文献是什么？ {#physicochemical-references}

[理化性质](api/features/07_physicochemical.md)中的功能不使用机器学习模型。它们分为理论公式、公开经验参数模型、DNAKit 内部透明规则和外部 Primer3 热力学结构预测。下表逐项列出实际依据；没有论文来源的内部规则会明确标注，不用不存在的引用补充包装。

| 功能 | 计算依据 | 参考文献或来源 |
| --- | --- | --- |
| THERMO-001 分子量 | 无水 DNA 残基质量求和、末端及 5′ 磷酸修正 | 当前实现记录的是标准寡核苷酸质量公式，尚未绑定原始论文；本地验证包含 Biopython 数值对照，但该对照不能替代科学来源引用。 |
| THERMO-014 260 nm 消光系数 | 相邻二核苷酸 hypochromicity 参数之和减内部单碱基参数 | Warshaw & Tinoco, 1966，[DOI 10.1016/0022-2836(66)90115-X](https://doi.org/10.1016/0022-2836(66)90115-X)；Cantor, Warshaw & Shapiro, 1970，[DOI 10.1002/bip.1970.360090909](https://doi.org/10.1002/bip.1970.360090909)。 |
| THERMO-002 Tm（Wallace） | `2 × (A+T) + 4 × (G+C)` 短寡核苷酸经验规则 | Wallace et al., 1979，[DOI 10.1093/nar/6.11.3543](https://doi.org/10.1093/nar/6.11.3543)。当前实现尚未把该 DOI 写入结果 provenance，FAQ 在此补充文献关系。 |
| THERMO-002 Tm（nearest-neighbor） | 相邻堆积、末端、对称性、盐浓度和链浓度模型 | SantaLucia, 1998，[DOI 10.1073/pnas.95.4.1460](https://doi.org/10.1073/pnas.95.4.1460)。 |
| THERMO-003 盐浓度修正 | SantaLucia 单价盐熵修正式 | SantaLucia, 1998，[DOI 10.1073/pnas.95.4.1460](https://doi.org/10.1073/pnas.95.4.1460)。 |
| THERMO-012 局部 Tm | 滑动窗口重复调用 Wallace 或 nearest-neighbor Tm | 不引入新的科学模型；引用继承所选择的 Wallace 1979 或 SantaLucia 1998 方法。 |
| EVAL-013 合成风险 | GC、同碱基连续、串联重复和倒置重复五项分量等权平均 | DNAKit 内部透明启发式规则；没有外部论文或供应商规则集，分数不是实验合成成功概率。 |
| THERMO-004 热力学参数 | 完全互补双链的 ΔH、ΔS、ΔG 和 Tm | SantaLucia, 1998，[DOI 10.1073/pnas.95.4.1460](https://doi.org/10.1073/pnas.95.4.1460)。 |
| THERMO-005 Nearest-neighbor | 与 THERMO-004 相同，同时返回逐堆积步骤明细 | SantaLucia, 1998，[DOI 10.1073/pnas.95.4.1460](https://doi.org/10.1073/pnas.95.4.1460)。 |
| THERMO-006 Duplex stability | 原生完整互补路径使用 SantaLucia；可选 mismatch/dangling 路径调用 Primer3 | SantaLucia, 1998，[DOI 10.1073/pnas.95.4.1460](https://doi.org/10.1073/pnas.95.4.1460)；Untergasser et al., 2012，[DOI 10.1093/nar/gks596](https://doi.org/10.1093/nar/gks596)。其中 `Tm > 设置温度` 的稳定布尔判据是 DNAKit 的结果解释规则。 |
| THERMO-007 碱基堆积 | 逐相邻二核苷酸查询 ΔH/ΔS 参数并计算 ΔG | SantaLucia, 1998，[DOI 10.1073/pnas.95.4.1460](https://doi.org/10.1073/pnas.95.4.1460)。 |
| THERMO-008 Hairpin | 用户安装的 Primer3 `ntthal` 热力学发卡结构预测 | Untergasser et al., 2012，[DOI 10.1093/nar/gks596](https://doi.org/10.1093/nar/gks596)及 [Primer3 官方手册](https://primer3.org/manual.html)；具体数值还取决于实际 Primer3 版本和参数目录。 |
| THERMO-009 Self-dimer | 用户安装的 Primer3 `ntthal` 热力学自二聚体预测 | Untergasser et al., 2012，[DOI 10.1093/nar/gks596](https://doi.org/10.1093/nar/gks596)及 [Primer3 官方手册](https://primer3.org/manual.html)；具体数值还取决于实际 Primer3 版本和参数目录。 |
| THERMO-010 Heterodimer | 用户安装的 Primer3 `ntthal` 热力学异二聚体预测 | Untergasser et al., 2012，[DOI 10.1093/nar/gks596](https://doi.org/10.1093/nar/gks596)及 [Primer3 官方手册](https://primer3.org/manual.html)；具体数值还取决于实际 Primer3 版本和参数目录。 |

这些引用说明算法和参数的来源，不表示任何计算结果等同于实验测量。Tm、ΔG 和结构结果仍受盐浓度、链浓度、温度、序列长度、化学修饰以及实际后端版本影响。

## 双链热力学扩展的计算依据和参考文献是什么？ {#duplex-thermodynamics-references}

[双链热力学扩展](api/features/19_duplex_thermodynamics.md)不使用机器学习。原生双链能量来自 SantaLucia 公开参数，其他结果是标准热力学关系、DNAKit 透明组合规则或外部 Primer3 计算。

| 功能 | 实际计算依据 | 参考文献或来源 |
| --- | --- | --- |
| 1. 完全互补双链 `ΔH/ΔS/ΔG/Tm` | 相邻堆积、末端、对称、单价盐和链浓度公式 | SantaLucia, 1998，[DOI 10.1073/pnas.95.4.1460](https://doi.org/10.1073/pnas.95.4.1460)。该 DOI 已写入原生结果 provenance。 |
| 2. 统一双链稳定性 | `native` 沿用第 1 项；`primer3-cli` 调用 `ntthal` 的 mismatch/dangling-end 结构模型 | 原生路径引用 SantaLucia 1998；Primer3 路径引用 Untergasser et al., 2012，[DOI 10.1093/nar/gks596](https://doi.org/10.1093/nar/gks596) 及 [Primer3 官方手册](https://primer3.org/manual.html)。`Tm > 设置温度` 是 DNAKit 的结果解释规则，没有单独论文。 |
| 3. 相邻碱基对步骤贡献 | 查相邻堆积 `ΔH/ΔS` 表并计算 `ΔG=ΔH−TΔS/1000` | SantaLucia, 1998，[DOI 10.1073/pnas.95.4.1460](https://doi.org/10.1073/pnas.95.4.1460)。不另引用新模型。 |
| 4. 条件与 Na⁺/K⁺ | 条件对象作数值校验和记录；原生模型用 `0.368(N−1)ln([Na⁺]+[K⁺])` 修正 `ΔS` | 盐修正引用 SantaLucia 1998。对象校验和 `Na⁺+K⁺` 字段合并是 DNAKit 实现规则，没有单独论文。 |
| 5. `Ka`、`Kd` 和双链比例 | `Ka=exp(−ΔG/RT)`、`Kd=1/Ka`，再解理想两态质量守恒方程 | `ΔG` 引用 SantaLucia 1998；平衡常数关系和二次方程是标准热力学与 DNAKit 透明代数实现，未绑定单独论文或拟合数据。 |
| 6. 理论熔解曲线 | 在每个温度重复第 1、5 项，对双链比例 0.5 交点做线性插值 | 引用继承 SantaLucia 1998；温度扫描和线性插值是 DNAKit 内部确定性规则，没有新论文。 |
| 7. 5′/3′ 末端稳定性 | 对两个等长末端窗口分别计算最近邻 `ΔG` 后比较 | 能量引用 SantaLucia 1998；窗口截取和“较高 `ΔG` 为较不稳定端”是 DNAKit 内部规则。 |
| 8. DMSO/甲酰胺修正 | Primer3 手册的线性经验加和式 | [Primer3 官方手册](https://primer3.org/manual.html)。默认 DMSO 因子 0.6 由手册引自 Musielski et al., 1981；甲酰胺公式引自 Blake & Delcourt, 1996，[DOI 10.1093/nar/24.11.2095](https://doi.org/10.1093/nar/24.11.2095)。该功能不是机理性自由能模型。 |
| 9. Primer3 CLI 扩展 | `oligotm` 计算 Tm，`ntthal` 计算发卡和二聚体；DNAKit 只做 adapter | Untergasser et al., 2012，[DOI 10.1093/nar/gks596](https://doi.org/10.1093/nar/gks596) 及 [Primer3 官方手册](https://primer3.org/manual.html)。还应记录实际 Primer3 版本和热力学参数目录。 |

## 二级结构性质的计算依据和参考文献是什么？ {#secondary-structure-references}

[二级结构性质](api/features/20_secondary_structure.md)包含 DNAKit 内部解析/公式和条件性 NUPACK adapter。只有显式运行第 6、7 项时，结构能量和平衡结果才是 NUPACK 计算值。

| 功能 | 实际计算依据 | 参考文献或来源 |
| --- | --- | --- |
| 1. Dot-bracket 解析 | 栈配对、连续嵌套配对合并和固定结构分类 | DNAKit 内部确定性解析器，没有用论文参数或训练模型。Dot-parens-plus 符号约定可参考 Zadeh et al., 2011，[DOI 10.1002/jcc.21596](https://doi.org/10.1002/jcc.21596)；`() [] {} <>` 的扩展括号解析是 DNAKit 实现。 |
| 2. 配对概率与窗口可接近性 | NUPACK 风格稠密概率矩阵的边缘概率派生量 | 矩阵语义可参考 [NUPACK 分析文档](https://docs.nupack.org/analysis/) 和 Zadeh et al. 2011。对称性/行和校验及“窗口边缘概率算术平均”是 DNAKit 内部规则，不是 NUPACK 联合开放概率。 |
| 3. 归一化 ensemble defect | 每个碱基未采用目标配对状态的期望比例 | Zadeh, Wolfe & Pierce, 2011，[DOI 10.1002/jcc.21633](https://doi.org/10.1002/jcc.21633)。DNAKit 只对已提供的概率矩阵计算该定义。 |
| 4. 目标结构概率 | `exp[−(Gtarget−Gensemble)/(RT)]` Boltzmann 关系 | 标准统计热力学公式，DNAKit 未绑定单独论文或拟合参数。自由能若来自 NUPACK，应再引用实际 NUPACK 模型与版本。 |
| 5. NUPACK 被动探测 | Python 包位置和元数据检查 | 工程环境检查，不执行 NUPACK，因此没有科学论文引用；该结果不得写成 NUPACK 预测。 |
| 6. NUPACK 单复合物 | 外部 NUPACK 4 的配分函数、MFE、配对概率、次优结构和 Boltzmann 抽样 | DNAKit provenance 当前记录 Zadeh et al. 2011，[DOI 10.1002/jcc.21596](https://doi.org/10.1002/jcc.21596)。NUPACK 4 官方引用页还指定 Fornace, Porubsky & Pierce, 2020，[DOI 10.1021/acssynbio.9b00523](https://doi.org/10.1021/acssynbio.9b00523)；应按实际调用功能与版本选择引用，见 [NUPACK 官方引用页](https://docs.nupack.org/)。 |
| 7. NUPACK tube 平衡 | 外部 NUPACK 4 的多链复合物与试管平衡分析 | Dirks et al., 2007，[DOI 10.1137/060651100](https://doi.org/10.1137/060651100)；Fornace, Porubsky & Pierce, 2020，[DOI 10.1021/acssynbio.9b00523](https://doi.org/10.1021/acssynbio.9b00523)；以及 [NUPACK 官方引用页](https://docs.nupack.org/)。DNAKit 的目标/非目标比例是对 adapter 返回浓度的内部汇总。 |

## 三维结构与力学性质的计算依据和参考文献是什么？ {#structure3d-references}

[三维结构与力学性质](api/features/21_structure3d.md)不从 DNA 序列预测三维结构。第 1–4 项是显式坐标的确定性解析/几何计算，第 5、6 项是对外部 3DNA/DSSR 结果的 adapter。

| 功能 | 实际计算依据 | 参考文献或来源 |
| --- | --- | --- |
| 1. 读取单个 PDB 模型 | legacy PDB `MODEL/ATOM` 固定列、DNA 残基表和备选构象选择 | [wwPDB PDB v3.3 坐标记录规格](https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html)。这是文件格式来源，不是科学预测论文。 |
| 2. 读取 PDB 多模型 | 与第 1 项相同，再按 `MODEL` 序号分组 | [wwPDB PDB v3.3 坐标记录规格](https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html)。多模型组织是格式规则，没有单独预测模型。 |
| 3. 显式坐标几何 | 质心、回转半径、张量特征值、球面点采样 SASA、体素体积、二面角、距离/角度氢键筛选和全局螺旋轴近似 | 当前结果 provenance 将整个原生几何实现标为 DNAKit 内部方法，未绑定统一论文。其中 SASA 是 Shrake–Rupley 风格的球面点暴露法，可参考 Shrake & Rupley, 1973，[DOI 10.1016/0022-2836(73)90011-9](https://doi.org/10.1016/0022-2836(73)90011-9)；体素体积、氢键阈值及 `C1′` 全局轴是 DNAKit 具体规则，不应写成 3DNA/DSSR 标准参数。 |
| 4. NMR/多模型 RMSF | 公共原子平移中心化后的均方位移平方根 | 标准 RMSF 数学定义的 DNAKit 内部实现，没有绑定单独论文。实现不做 Kabsch 旋转拟合，也不推断持久长度或力学模量。 |
| 5. 3DNA `bp_step.par` 解析 | 读取 3DNA 已计算的 12 个刚体参数，DNAKit 只求平均值、`360/twist` 和 pitch | Lu & Olson, 2003，[DOI 10.1093/nar/gkg680](https://doi.org/10.1093/nar/gkg680) 及 [3DNA 参数官方说明](https://x3dna.org/highlights/schematic-diagrams-of-base-pair-parameters)。 |
| 6. DSSR JSON 摘要 | 读取 DSSR 已识别的核苷酸、配对、螺旋、stem、hairpin 和氢键计数 | Lu, Bussemaker & Olson, 2015，[DOI 10.1093/nar/gkv716](https://doi.org/10.1093/nar/gkv716) 及 [DSSR JSON 官方说明](https://x3dna.org/highlights/dssr-output-in-json-format)。DSSR 完成结构识别，DNAKit 只解析已有 JSON。 |

## 光学与浓度换算的计算依据和参考文献是什么？ {#optics-concentration-references}

[光学与浓度换算](api/features/18_optics_concentration.md)不使用机器学习。它使用公开的单链 `ε260` 经验参数、Beer–Lambert 定律、分子量与单位换算，以及调用方显式提供的双链/修饰参数。

| 功能 | 实际计算依据 | 参考文献或来源 |
| --- | --- | --- |
| 1. 260 nm 单链摩尔消光系数 | `Σ相邻二核苷酸系数 − Σ内部单碱基系数` | Warshaw & Tinoco, 1966，[DOI 10.1016/0022-2836(66)90115-X](https://doi.org/10.1016/0022-2836(66)90115-X)；Cantor, Warshaw & Shapiro, 1970，[DOI 10.1002/bip.1970.360090909](https://doi.org/10.1002/bip.1970.360090909)。两篇已写入该结果 provenance。 |
| 2. 单链/双链理论光学性质 | 单链沿用第 1 项；双链使用 `13200 × bp` 平均式或 `(ε1+ε2)(1−h)`；分子量用无水残基求和 | 单链 `ε` 引用上述两篇论文。`13200 M⁻¹·cm⁻¹/bp` 是实现采用的传统 `1 OD260 ≈ 50 µg/mL` 双链平均换算，当前源码未绑定原始论文；`h`、`Δε` 和 `ΔMW` 由调用方提供，应引用其实际实验、文献或厂商数据表。无水残基分子量公式尚未绑定原始论文。 |
| 3. 1 OD260 对应 nmol/质量 | `10⁶/ε` 和 `1000×MW/ε` | 由 OD260 定义和单位换算直接派生，没有新论文；引用应继承所选 `ε260` 和 `MW` 来源。 |
| 4. A260 到浓度和总量 | `c=A/(εl)`，再加稀释倍数、`m=nMW`、`n=cV` 及显式染料减除 | Beer–Lambert 标准物理定律与 DNAKit 透明单位换算，当前实现未绑定特定现代论文。`ε`、`MW` 和染料因子应分别引用它们的实际来源。 |
| 5. 浓度/物质的量/质量互换 | `m=nMW`、`n=m/MW`、`c=n/V`、`n=cV` | DNAKit 内部量纲换算，不是经验模型，没有单独论文。 |
| 6. 染料与修饰校正 | `Σcount×Δε`、`Σcount×ΔMW` 和 `ΣA标记峰值×校正因子` | DNAKit 不内置染料参数表，只做确定性加减法，因此没有统一论文。调用方必须为每个参数记录真实实验、文献或厂商来源。 |

上述引用只说明算法、参数或外部工具的来源。对于调用方提供的自由能、hypochromicity、染料修正或 3DNA/DSSR/NUPACK 输出，还必须记录实际参数表、软件版本、输入条件和数据来源。

## 深度学习性质预测有哪些参考文献？ {#deep-learning-property-prediction-references}

下表覆盖[深度学习性质预测](api/features/23_deep_learning_property_prediction.md)实际集成的 54 个功能所使用的主要模型论文、任务数据集和训练协议。多个输出头共用同一篇论文，因此不按功能重复列出。

| 对应功能 | 参考文献 |
| --- | --- |
| RNA-seq、CAGE、PRO-cap、ATAC-seq、DNase-seq、ChIP-seq、剪接和接触图 | Avsec et al., *Advancing regulatory variant effect prediction with AlphaGenome*, **Nature** 649, 1206–1218 (2026), [DOI 10.1038/s41586-025-10014-0](https://doi.org/10.1038/s41586-025-10014-0)。 |
| 人类和小鼠调控轨道 | Avsec et al., *Effective gene expression prediction from sequence by integrating long-range interactions*, **Nature Methods** 18, 1196–1203 (2021), [DOI 10.1038/s41592-021-01252-x](https://doi.org/10.1038/s41592-021-01252-x)。 |
| NT Revised 18 个分类任务 | 主干模型：Avsec et al., *Effective gene expression prediction from sequence by integrating long-range interactions*, **Nature Methods** 18, 1196–1203 (2021), [DOI 10.1038/s41592-021-01252-x](https://doi.org/10.1038/s41592-021-01252-x)；任务定义：Dalla-Torre et al., *Nucleotide Transformer: building and evaluating robust foundation models for human genomics*, **Nature Methods** 22, 287–297 (2025), [DOI 10.1038/s41592-024-02523-z](https://doi.org/10.1038/s41592-024-02523-z)，以及[修订数据集说明](https://huggingface.co/datasets/InstaDeepAI/nucleotide_transformer_downstream_tasks_revised)；完整微调协议：Wu et al., *GENERator: A Long-Context Generative Genomic Foundation Model*, arXiv (2025), [arXiv:2502.07272](https://arxiv.org/abs/2502.07272)，附录 C.4。 |
| Genomic Benchmarks 9 个分类任务 | 主干模型：Avsec et al., *Effective gene expression prediction from sequence by integrating long-range interactions*, **Nature Methods** 18, 1196–1203 (2021), [DOI 10.1038/s41592-021-01252-x](https://doi.org/10.1038/s41592-021-01252-x)；任务数据集：Grešová et al., *Genomic benchmarks: a collection of datasets for genomic sequence classification*, **BMC Genomic Data** 24, 25 (2023), [DOI 10.1186/s12863-023-01123-8](https://doi.org/10.1186/s12863-023-01123-8)；完整微调协议：Wu et al., *GENERator: A Long-Context Generative Genomic Foundation Model*, arXiv (2025), [arXiv:2502.07272](https://arxiv.org/abs/2502.07272)，附录 C.4。 |
| 14 类单碱基基因组分割的基础编码器 | Dalla-Torre et al., *Nucleotide Transformer: building and evaluating robust foundation models for human genomics*, **Nature Methods** 22, 287–297 (2025), [DOI 10.1038/s41592-024-02523-z](https://doi.org/10.1038/s41592-024-02523-z)。 |
| 14 类单碱基基因组分割头 | de Almeida et al., *Annotating the genome at single-nucleotide resolution with DNA foundation models*, **Nature Methods** (2025), [DOI 10.1038/s41592-025-02881-2](https://doi.org/10.1038/s41592-025-02881-2)。 |
| 长上下文零样本变异效应和外显子概率 | Brixi et al., *Genome modelling and design across all domains of life with Evo 2*, **Nature** (2026), [DOI 10.1038/s41586-026-10176-5](https://doi.org/10.1038/s41586-026-10176-5)。 |
| 等位基因条件概率变异效应 | Wu et al., *GENERator: A Long-Context Generative Genomic Foundation Model*, arXiv (2025), [arXiv:2502.07272](https://arxiv.org/abs/2502.07272)；Li et al., *GENERator-v2: Reconciling Coarse Tokenization with Single-Nucleotide Resolution in Genomic Language Modeling*, **bioRxiv** (2026), [DOI 10.64898/2026.01.27.702015](https://doi.org/10.64898/2026.01.27.702015)。 |
| 中心法则、分类学、物种、蛋白定位/稳定性、ncRNA 家族和序列对互作 | He et al., *Generalized biological foundation model with unified nucleic acid and protein language*, **Nature Machine Intelligence** 7, 942–953 (2025), [DOI 10.1038/s42256-025-01044-4](https://doi.org/10.1038/s42256-025-01044-4)。 |

直接推理还依赖官方代码或已训练权重，其来源如下：

- RNA-seq 等 11 类轨道：[AlphaGenome research](https://github.com/google-deepmind/alphagenome_research) 和 [`google/alphagenome-all-folds`](https://huggingface.co/google/alphagenome-all-folds)。
- 人类/小鼠调控轨道：[Enformer 官方实现](https://github.com/google-deepmind/deepmind-research/tree/master/enformer) 和 [`EleutherAI/enformer-official-rough`](https://huggingface.co/EleutherAI/enformer-official-rough)。
- 27 个任务分类 checkpoint：从[统一的 Google Drive checkpoint 文件夹](https://drive.google.com/drive/folders/1lrZXzkrgAJMqM0wAmnIeZ4DEp0XFNIRI?usp=sharing)下载；任务定义来自 [NT Revised 数据集](https://huggingface.co/datasets/InstaDeepAI/nucleotide_transformer_downstream_tasks_revised)和 [Genomic Benchmarks 仓库](https://github.com/ML-Bioinfo-CEITEC/genomic_benchmarks)，DNAKit 不随包重新分发这些权重。
- 单碱基基因组分割：[Nucleotide Transformer 仓库](https://github.com/instadeepai/nucleotide-transformer) 和 [`InstaDeepAI/segment_nt`](https://huggingface.co/InstaDeepAI/segment_nt)。
- 零样本变异效应/外显子概率：[Evo 2 仓库](https://github.com/ArcInstitute/evo2)、[`arcinstitute/evo2_7b`](https://huggingface.co/arcinstitute/evo2_7b)、[`arcinstitute/evo2_7b_base`](https://huggingface.co/arcinstitute/evo2_7b_base) 和 [`schmojo/evo2-exon-classifier`](https://huggingface.co/schmojo/evo2-exon-classifier)。
- 等位基因条件概率变异效应：[GENERator 仓库](https://github.com/GenerTeam/GENERator) 和 [`GenerTeam/GENERator-v2-eukaryote-1.2b-base`](https://huggingface.co/GenerTeam/GENERator-v2-eukaryote-1.2b-base)。
- 核酸/蛋白下游任务：[LucaOne](https://github.com/LucaOne/LucaOne)、[LucaOneTasks](https://github.com/LucaOne/LucaOneTasks) 和 [Zenodo 10.5281/zenodo.15171943](https://doi.org/10.5281/zenodo.15171943)。

这些引用说明模型、训练任务头和权重的来源，不表示预测结果等同于实验测量或临床结论。
