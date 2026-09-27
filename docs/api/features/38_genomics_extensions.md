# 基因组与 DNA 工程

提供基因型 PCA、ANI 相似度、育种预测、原核基因识别、甲基化统计、DNA 存储和纳米结构接口。安装依赖、完整参数与使用限制见 [FAQ：参数与补充说明](../../faq/api/38_genomics_extensions.md)，方法作者、论文与许可证见[来源 FAQ](../../faq/genomics_sources.md)。

<span id="pip-genomics"></span>
当前源码可通过 `python -m pip install -e ".[genomics]"` 一并安装 NumPy、Pyrodigal、pysam 和 PyFastANI。下面的基因型 PCA、ANI、关系矩阵、育种预测和原核 CDS 预测，以及[基于 pysam 的 BAM 统计](32_sequencing_analysis.md#bio-14)，均在 Python 进程内计算，无需另外安装专业命令行程序。

该安装选项不包含 DNA 存储、纳米结构依赖或全部外部分析。43 项外部分析的迁移状态见[FAQ：外部分析如何转为 pip 可用功能](../../faq/external_analysis_migration.md)。新增能力尚未发布到 PyPI。

## 1) 基因组关系矩阵与育种预测

- **作用：** 计算样本间的加性基因组关系，并使用训练样本的基因型和表型预测新样本的表型值。
- **计算方法：** 关系矩阵使用 VanRaden 方法 1；预测使用带截距、固定惩罚系数的 RR-BLUP／岭回归。预测的缺失值填补、标记筛选和中心化仅使用训练样本统计量。
- **API：** `dnakit.genomics.genomic_relationship_matrix(genotypes[必须], missing[可选])`；`dnakit.genomics.genomic_prediction(train_genotypes[必须], phenotypes[必须], predict_genotypes[必须], ridge_lambda[必须], missing[可选])`
- **输入：** 行为样本、列为标记的二倍体剂量矩阵，数值范围为 0–2；训练和预测矩阵的标记顺序及等位基因方向必须一致。每个训练样本需有一个有限表型值，惩罚系数必须为正数。
- **输出：** 返回关系矩阵和基因组预测字典，包含样本数与预测值。
- **示例代码：**

```python
from dnakit.genomics import genomic_relationship_matrix, genomic_prediction

genotypes = [[0, 0, 2], [1, 2, 1], [2, 1, 0], [0, 1, 1]]
relationship = genomic_relationship_matrix(genotypes)
prediction = genomic_prediction(
    genotypes,
    [2.0, 5.0, 4.0, 3.0],
    [[1, 1, 1]],
    ridge_lambda=1.0,
)
print(relationship["n_samples"])
print([round(value, 3) for value in prediction["predictions"]])
```

- **示例结果：**

```text
4
[3.592]
```

[参数、缺失值处理与验证约定](../../faq/api/38_genomics_extensions.md#breeding)

## 2) 原核基因识别

- **作用：** 从细菌或古菌的线性 contig 识别 CDS，返回坐标、编码序列、蛋白和模型评分。
- **计算方法：** 调用本地 Pyrodigal／Prodigal；默认使用上游预训练模型，DNAKit 将输出坐标转换为 0-based、左闭右开区间，并保留后端版本。
- **API：** `dnakit.genomics.predict_prokaryotic_genes(sequence[必须], sequence_id[可选], meta[可选], training_sequence[可选])`
- **输入：** 一条非空、无缺口的 IUPAC DNA 序列；模糊碱基转为 N 并保留位置。示例人工序列仅用于演示接口，不保证检出基因。
- **输出：** 返回 `GenePredictionResult`，示例展示 `sequence_id`、`sequence_length`、`coordinate_system`。
- **示例代码：**

```python
from dnakit.genomics import predict_prokaryotic_genes

sequence = "ATG" + "GCT" * 120 + "TAA"
result = predict_prokaryotic_genes(sequence, sequence_id="example", meta=True)
print(result.sequence_id, result.sequence_length)
print(result.coordinate_system)
```

- **示例结果：** 返回结果的序列标识、长度和坐标约定；基因预测保存在 `result.genes`。

```text
example 366
0-based half-open
```

[参数、训练模式与结果字段](../../faq/api/38_genomics_extensions.md#genes)

## 3) 甲基化计数与汇总

- **作用：** 读取 bedMethyl 位点计数，按染色体、链、修饰代码和 motif 汇总甲基化比例。
- **计算方法：** 解析 modkit 的 18 列扩展 bedMethyl，筛选记录后计算各组的修饰读数总和／有效覆盖读数总和；5mC 与 5hmC 分别统计。
- **API：** `dnakit.genomics.read_bedmethyl(path[必须])`；`dnakit.genomics.summarize_methylation(records[必须], modification_codes[可选], min_coverage[可选])`
- **输入：** bedMethyl 文件路径，支持 gzip／BGZF；第 11 列为 0–100 的百分数。下面在临时目录创建一条覆盖数为 10、修饰读数为 5 的演示记录。
- **输出：** `read_bedmethyl` 返回记录列表；`summarize_methylation` 返回修饰比例及分组统计。
- **示例代码：**

```python
from pathlib import Path
from tempfile import TemporaryDirectory
from dnakit.genomics import read_bedmethyl, summarize_methylation

with TemporaryDirectory() as directory:
    path = Path(directory) / "example.bedmethyl"
    path.write_text(
        "chr1\t0\t1\tm\t10\t+\t0\t1\t255,0,0\t10\t50\t5\t5\t0\t0\t0\t0\t0\n",
        encoding="utf-8",
    )
    records = read_bedmethyl(path)
summary = summarize_methylation(records, modification_codes=["m"], min_coverage=10)
group = summary["groups"][0]
print(summary["records_selected"])
print(group["chromosome"], group["modification_code"], group["fraction_modified"])
```

- **示例结果：**

```text
1
chr1 m 0.5
```

[参数、读取限制与统计边界](../../faq/api/38_genomics_extensions.md#methylation)

## 4) DNA 数据存储编码与解码

- **作用：** 将字节数据编码为 DNA 片段，并从完整数据包恢复原始字节。
- **计算方法：** 使用 Chamaeleo 的固定 Huffman Goldman 编码器；解码时按片段索引重排，并检查长度和 SHA-256 摘要。
- **API：** `dnakit.genomics.encode_dna_storage(data_hex[必须], segment_bytes[可选])`；`dnakit.genomics.decode_dna_storage(packet[必须])`
- **输入：** 编码接收字节数据的十六进制字符串；解码接收完整 JSON packet，包含 DNA 片段、索引、长度和摘要。
- **输出：** 编码返回 DNA 存储包字典；解码返回恢复数据的 `data_hex`。
- **示例代码：**

```python
from dnakit.genomics import encode_dna_storage, decode_dna_storage

packet = encode_dna_storage("DNAKit".encode().hex(), segment_bytes=32)
restored = decode_dna_storage(packet)
print(bytes.fromhex(restored["data_hex"]).decode())
```

- **示例结果：**

```text
DNAKit
```

[参数、数据包保存与检错边界](../../faq/api/38_genomics_extensions.md#storage)

## 5) DNA 纳米结构与 oxDNA 导出

- **作用：** 构造互补 DNA 双链的 scadnano 几何设计，并导出 oxDNA 配置和拓扑文本。
- **计算方法：** 由 scadnano 创建两条反向互补链、分配碱基并转换坐标；导出前校验序列与几何数据，不执行折叠优化或分子动力学模拟。
- **API：** `dnakit.genomics.design_dna_duplex(sequence[必须])`；`dnakit.genomics.export_oxdna_design(design[必须])`
- **输入：** 完整 A/C/G/T 序列，或所有链均已分配完整序列的 scadnano JSON；当前仅支持线性链导出。
- **输出：** 设计返回双链结构字典；oxDNA 导出返回拓扑和坐标文本。
- **示例代码：**

```python
from dnakit.genomics import design_dna_duplex, export_oxdna_design

design = design_dna_duplex("AACCGGTT")
export = export_oxdna_design(design)
print(export["strand_count"], export["nucleotide_count"])
print(export["topology"].splitlines()[0])
```

- **示例结果：** 第一行为链数和核苷酸数；第二行为 oxDNA 拓扑头部的核苷酸数和链数。

```text
2 16
16 2
```

[参数、文件保存与导出限制](../../faq/api/38_genomics_extensions.md#nanotechnology)

## 6) 基因型 PCA {#pca}

- **作用：** 从已整理的二倍体基因型剂量矩阵计算样本主成分和位点载荷。
- **计算方法：** 位点均值填补缺失值，删除全缺失／不变位点，以等位基因频率标准化后计算 SVD。
- **API：** `dnakit.genomics.genotype_pca(genotypes, components=2, sample_ids=None, verbose=True)`。
- **输入：** 样本×位点矩阵，剂量范围 0–2；缺失用 `None`。通过 `.[genomics]` 或 `.[bio-analysis]` 安装 NumPy。
- **输出：** 返回字典，示例展示 `eigenvalues`、`marker_indices`。
- **示例代码：**

```python
from dnakit.genomics import genotype_pca

result = genotype_pca(
    [[0, 0], [1, 2], [2, 1]],
    components=2,
    sample_ids=["sample_A", "sample_B", "sample_C"],
)
print([round(value, 6) for value in result["eigenvalues"]])
print(result["marker_indices"])
```

- **示例结果：**

```text
[3.0, 1.0]
[0, 1]
```

还返回 `scores`、`loadings`、方差比例和预处理记录。此接口为稠密矩阵 PCA，不读取 BED/VCF，也不执行 LD 剪枝或完整 PLINK 流程。[参数与特征值尺度](../../faq/api/38_genomics_extensions.md#pca)。

## 7) 基因组 ANI 相似度 {#ani}

- **作用：** 估计两个近缘原核基因组同源片段的平均核苷酸一致性。
- **计算方法：** 直接调用 PyFastANI 中的 FastANI 算法，保持 contig 边界，计算 query→reference 方向。
- **API：** `dnakit.genomics.average_nucleotide_identity(query, reference, threads=1, verbose=True)`。
- **输入：** 两个组装 FASTA 路径（支持 gzip），或两个 DNA 字符串；安装 `.[genomics]`，无需 FastANI 可执行文件。
- **输出：** 返回字典，示例展示 `status`、`ani_percent`、`matching_fragments`、`query_fragments`。
- **示例代码：** 准备自己的组装文件后运行：

```python
from dnakit.genomics import average_nucleotide_identity

result = average_nucleotide_identity("query.fna", "reference.fna", threads=4)
print(result["status"], result["ani_percent"])
print(result["matching_fragments"], result["query_fragments"])
```

- **示例结果：** 有可报告命中时 `status` 为 `matched`，ANI 单位为百分数；无可报告命中时为 `no_match`，ANI 为 `None`。短于查询片段的输入不能产生有效估计。实际数值取决于输入；人工同序列验证结果接近 100%。[输入约定、方向性与完整参数](../../faq/api/38_genomics_extensions.md#ani)。
