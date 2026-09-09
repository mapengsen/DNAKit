# 原始 DNA 的基因名称与区间注释

`dnakit.annotation.annotate_genes()` 输入一段基因组 DNA 和物种，返回具体基因名称、ID、方向与输入序列中的区间。`assembly` 可选：首次未填写时，查询该物种当前的 NCBI RefSeq 参考组装；成功后把准确的组装和注释文件固定在项目缓存中。

这是基于参考序列比对的已有注释映射。所选组装必须有基因注释；无法获得正式名称时保留基因 ID 或 locus tag，不推测名称。本接口面向基因组 DNA，不用于拼接后的 cDNA、转录本或从头基因预测。

## 安装和命令行

需要安装 DNAKit 和外部程序 minimap2。导入模块不会下载参考文件或执行外部程序。

安装当前工作区的 DNAKit：

```bash
python -m pip install -e .
```

通过 Conda 安装 minimap2：

```bash
conda install -c conda-forge -c bioconda minimap2
```

仅提供物种，自动选择并缓存默认组装；`query.fasta` 应含一条基因组 DNA 序列：

```bash
python -m dnakit.cli annotate-genes \
  --input query.fasta \
  --species "Homo sapiens" \
  --cache-dir .dnakit/gene_annotation \
  --output genes.json
```

明确指定组装版本：

```bash
python -m dnakit.cli annotate-genes \
  --input query.fasta \
  --species "Homo sapiens" \
  --assembly GCF_000001405.40 \
  --cache-dir .dnakit/gene_annotation \
  --output genes.json
```

使用已准备的同版本 FASTA 和 GFF3，完全离线运行：

```bash
python -m dnakit.cli annotate-genes \
  --input query.fasta \
  --species "Homo sapiens" \
  --reference-fasta reference.fna \
  --annotation-gff annotation.gff3 \
  --cache-dir .dnakit/gene_annotation \
  --output genes.json \
  --no-download
```

可通过 `--minimap2 /path/to/minimap2` 指定可执行文件，`--threads 4` 设置线程数。已有输出默认不覆盖，需要覆盖时显式增加 `--overwrite`。输入 FASTA 和本地参考文件支持 gzip；一次命令处理一条序列，多个样本可通过 Python 循环调用并复用同一缓存。

## Python 接口

```python
from dnakit.annotation import GeneAnnotationConfig, annotate_genes
from dnakit.io import read

with read("query.fasta", format="fasta") as source:
    record = next(iter(source))

result = annotate_genes(
    record.sequence.symbols,
    species="Homo sapiens",  # 也接受 NCBI TaxID，例如 9606
    assembly=None,
    sequence_id=record.id,
    config=GeneAnnotationConfig(
        cache_dir=".dnakit/gene_annotation",
        threads=4,
    ),
)
print(result.reference["assembly"])
for gene in result.gene_records:
    print(gene.to_dict())
for segment in result.segments:
    print(segment.to_dict())
```

`assembly` 接受版本化 `GCA_/GCF_` 编号、现有固定别名，以及能在指定物种下唯一解析的组装名称。如果名称同时匹配 RefSeq/GenBank 两套记录，会提示改用明确编号。显式指定旧版时不自动替换为新版。

## 输出和坐标

- `reference`：实际组装编号、组装名、注释元数据、文件 SHA256、来源和选择方式。自动下载同时验证物种与组装的归属；菌株可匹配其所属物种。
- `gene_records`：一个基因实例一条记录，允许重叠。`gene_name` 缺失时为 `null`，`display_name` 回退到 locus tag 或 ID。重复位点的映射身份由 `mapping_id` 保留。
- `segments`：完整覆盖输入序列、彼此不重叠的分段；类型为 `gene`、`intergenic` 或 `unknown`。重叠区间保留多个基因。
- `alignments`：比对位置、方向、MAPQ、相似度和是否达到接受阈值。

Python 对象中的 `Interval` 为 **0-based 半开区间**；`to_dict()` 和 CLI 的 JSON 均为 **1-based 闭区间**。例如输入位置 1–12000，在 Python 内部是 `Interval(0, 12000)`。

对示意基因 A:1–12000、B:18001–45000、C:43000–52000，`gene_records` 保留三条记录；`segments` 为：

```text
1–12000      gene         A
12001–18000  intergenic
18001–42999  gene         B
43000–45000  gene         B + C
45001–52000  gene         C
```

`gene_records.start/end` 是输入中的外包区间，`parts` 才是可靠映射到该基因的具体片段；存在插入或歧义时，中间可能有空洞。正反链和插入/缺失通过 CIGAR 逐块映射，不进行简单的整体坐标平移。GFF3 明确声明环状染色体时，支持跨原点的基因区间。

## 默认组装与可靠性

未提供 `assembly` 时，查询该物种及其菌株下的当前 RefSeq 参考组装。如果找不到参考组装、存在多个候选、物种不匹配或没有 GFF3 基因注释，会返回明确错误；不会随意挑选其他组装。尚未收录的物种或组装可以提供本地文件。

同一缓存中的相同请求会复用第一次成功获取的组装和注释，不再联网解析默认值。参考 FASTA、GFF3 与 minimap2 索引都校验 SHA256；索引按参考内容、preset 和 minimap2 版本隔离。要更新默认组装或注释，请使用新的 `cache_dir`。`--no-download` 可复用已固定的缓存，或者使用一对本地参考文件。

本地文件检查序列 ID、长度、基因边界及可获得的组装声明。离线时物种由调用者声明，结果会记录 `species_validated=false`；FASTA/GFF3 的物种和版本配套关系仍由提供者负责。GFF3 必须包含 gene、pseudogene 或免疫基因片段这类基因级行，不能只含 CDS/转录本。染色体 ID 必须完全匹配，不自动猜测 `chr1` 与 accession 的对应关系。

默认 minimap2 preset 为 `asm5`，适合较长且准确的基因组 DNA；可根据数据修改 `preset`。默认只接受主比对、MAPQ ≥ 20、相似度 ≥ 0.90、实际比对长度 ≥ 100 bp；MAPQ=255 视为未知。阈值是可配置筛选条件，不代表经过校准的基因识别准确率。

**`intergenic` 仅表示可靠映射位置不被本次所选 GFF3 的基因特征覆盖。** 未映射、低置信度、多条可信主比对重叠、输入中的插入，以及参考 contig 没有注释覆盖信息时均输出 `unknown`。缺少注释不等于生物学上不存在基因。

首次下载大型基因组和建立索引可能需要较多磁盘、内存和时间。配置可限制参考长度、下载量、GFF 行数、索引/比对输出大小和运行时间；单批索引避免多批索引导致 MAPQ 不可靠。

## 验证与来源

运行坐标映射、默认组装解析、缓存完整性和 CLI 测试；安装 minimap2 后还会运行真实比对测试：

```bash
PYTHONPATH=src python -m pytest -q others/tests/unit/test_gene_annotation.py
```

参考：[NCBI 参考组装选择规则](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/data-processing/policies-annotation/genome-processing/reference-selection/)、[NCBI 组装数据包](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-packages/genome/)、[NCBI GFF3 说明](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/file-formats/annotation-files/about-ncbi-gff3/)、[minimap2 手册](https://lh3.github.io/minimap2/minimap2.html)。

本版新增：多物种、可选 assembly 的参考基因注释流程、重叠区间输出及可复现缓存。
