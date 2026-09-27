# 基因组与 DNA 工程：参数与补充说明

作用、计算方法和示例见[功能教程](../../api/features/38_genomics_extensions.md)。以下接口默认显示进度，可使用 `verbose=False` 关闭。

## 基因组关系矩阵与育种预测 {#breeding}

安装数值计算依赖：

```bash
python -m pip install -e ".[bio-analysis]"
```

关系矩阵使用 VanRaden 方法 1 的加性关系定义；常见亲缘系数矩阵为其二分之一。预测使用固定惩罚系数的 RR-BLUP／岭回归，`ridge_lambda` 必须显式提供；示例值不代表最优值。按育种周期、群体或家系划分验证集，避免亲缘泄漏。

缺失剂量可以使用 `None`，默认以训练样本均值填补；`missing="error"` 可要求完整输入。全缺失和等位基因单态训练标记会排除，返回其索引。预测样本不参与频率估计或参数拟合。当前只接受二倍体；不包含多倍体、多性状、环境协变量、REML 估计或亲本配组优化。

## 原核基因识别 {#genes}

安装真实 Pyrodigal 后端：

```bash
python -m pip install -e ".[genomics-genes]"
```

教程中的人工序列仅演示接口，不保证模型检出基因。真实分析每次传入一条细菌／古菌线性 contig。`meta=True` 使用上游预训练模型；`meta=False` 需显式提供适于训练的 `training_sequence`，短序列训练失败不会悄悄切换模式。`translation_table` 只用于训练模式。

结果包括 CDS、蛋白、模型评分、翻译表与后端版本。坐标为输入序列上的 **0-based、左闭右开**；负链 CDS 已反向互补。模糊 IUPAC 字母转为 N 并保留原位置。该接口识别原核 CDS，不预测真核外显子／内含子或赋予基因功能名称。环状 DNA 需要先明确线性化位置。

## 甲基化计数与汇总 {#methylation}

接受 modkit 的 18 列扩展 bedMethyl，包括 gzip／BGZF。第 11 列使用百分数单位（0–100）。结果按染色体、链、修饰代码和 motif 分组，比例为修饰读数总和／有效覆盖读数总和；零覆盖返回空比例，不混合 5mC 与 5hmC。默认 `max_records=1_000_000`、`max_bytes=268_435_456`（解压后字节数）、`max_line_bytes=65_536`；超限会报错，不返回截断结果。较大数据应分区或显式调整限制。

该统计不执行差异甲基化检验或疾病预测。从带修饰信息的 BAM 生成 bedMethyl，可使用 [外部分析](../../api/features/39_external_genomics.md) 的 `modkit-pileup`。

## DNA 数据存储编码与解码 {#storage}

安装编码器：

```bash
python -m pip install -e ".[dna-storage]"
```

通过 Chamaeleo 的固定 Huffman Goldman 编码器转换字节。必须保存完整 JSON packet：索引、长度和 SHA-256 存在 JSON 中，未编码到 DNA 链内。支持完整片段乱序重排；丢段、损坏或摘要不符会报错。摘要用于检错，不提供纠错或认证。该接口不包含原论文完整的重叠冗余存档协议，也不执行合成或测序。

## DNA 纳米结构与 oxDNA 导出 {#nanotechnology}

安装结构设计依赖：

```bash
python -m pip install -e ".[dna-nanotechnology]"
```

`design` 是可保存为 `.sc` 文件的 scadnano JSON。导出结果包含 `configuration`、`topology`，可供后续 oxDNA 工作流保存使用。也可传入自行设计的 scadnano JSON；所有链必须完整分配 A/C/G/T。

当前仅导出线性链：验证的 scadnano 0.21.1 会丢失环状链闭合键，因此接口拒绝环状设计。无效或非有限几何参数、上游无法处理的偏移会返回明确错误。结构构造和文件导出不等于折叠优化、能量松弛或分子动力学模拟。

## 基因型 PCA {#pca}

`genotype_pca()` 使用 NumPy，安装 `.[genomics]` 或较小的 `.[bio-analysis]` 即可。输入是样本×位点的二倍体双等位剂量矩阵，允许 0–2 的小数剂量，以 `None` 或 NaN 表示缺失。样本名称可通过 `sample_ids` 传入，必须非空、唯一且数量匹配。

先用有观测样本估计各位点均值 μ 和频率 p=μ/2，删除全缺失或观测值不变的位点，再以 μ 填补其余缺失值，按 `(dosage−μ)/sqrt(2p(1−p))` 标准化并计算 SVD。返回样本得分 `U*S`、保留位点载荷 `V` 和协方差特征值 `S²/(n_samples−1)`；方差比例的分母包含全部主成分。`components` 必须在 1 到 `min(n_samples−1, n_markers_used)` 之间。

`marker_indices` 映射载荷行至原始位点；均值和频率保留原始列数，全缺失列为 null。填补计数只统计保留位点。最大绝对载荷取正以固定符号；相等特征值对应的子空间仍可能随数值库旋转。输入不被修改，结果可 JSON 序列化。

这是内存中的稠密矩阵 PCA；不读取 BED/VCF、不剪枝 LD，也不自动排除亲缘样本。它与 PLINK 的关系矩阵特征值采用不同尺度。应先统一位点／等位基因方向并完成 QC；大规模队列继续使用专业后端。

## 基因组 ANI {#ani}

`average_nucleotide_identity(query, reference)` 使用 `.[genomics]` 中的 PyFastANI；无需单独安装 FastANI 可执行程序。它接受 DNA 字符串或本地 FASTA 路径（支持 gzip）。仅由 IUPAC 字母组成的字符串按序列解释；这种名称的文件使用 `./文件名` 或 `Path` 对象。FASTA 内所有记录属于同一个基因组，contig 边界保持不变；不把多条记录当作不同样本。

默认使用 FastANI 的 3,000 bp 查询片段等参数，适合近缘原核基因组；短 contig 与不完整末端片段不能贡献完整查询片段。当前实现将序列读入内存；模糊碱基转换为 N 并记录数量。

只计算 query→reference 方向；需要双向结果时交换输入再调用。`ani_percent` 是同源片段的百分比一致性；`query_aligned_fraction` 为匹配片段数／完整查询片段数，不是全基因组碱基覆盖率。无可报告命中时返回 `status="no_match"`，ANI 和匹配比例为 null，不能解释为 0% 一致性。

结果记录后端版本、参数、输入摘要、片段数和方法引用。默认进度条显示加载、建索引和查询阶段；`verbose=False` 关闭进度，不屏蔽上游对短参考序列等情况的警告。依赖的预编译 wheel 是否可用取决于平台。

## 完整参数

::: dnakit.genomics
    options:
      members:
        - genotype_pca
        - average_nucleotide_identity
        - genomic_relationship_matrix
        - genomic_prediction
        - predict_prokaryotic_genes
        - read_bedmethyl
        - summarize_methylation
        - encode_dna_storage
        - decode_dna_storage
        - design_dna_duplex
        - export_oxdna_design
      show_source: false
