# 快速入门

本页只调用当前已经实现并测试的公开 API。内部坐标统一为 0-based 半开区间；标准化修改会保留审计记录，序列计算不会默认跨越 Gap。

## 1. 一个核心对象完成单条和多条输入 { #core-standardization }

```python
import dnakit
from dnakit.ops import delete, reverse_complement

dna = dnakit.DNA(
    " acgu\n",
    id="seq-1",
    metadata={"species": "synthetic"},
    features=[{"type": "motif", "start": 0, "end": 2}],
)
dataset = dnakit.DNA(["ACGT", "TTAA"])

edited = delete(dna, 1, 2)
reversed_dna = reverse_complement(edited)

print(dnakit.__version__)  # 0.1.0.dev0
print(dna.symbols)         # ACG
print(dataset.ids)         # ('sequence_1', 'sequence_2')
print(reversed_dna.symbols)
```

`DNA(...)` 自动标准化原始文本并把审计保存在 `dna.normalization`。单条、多条、ID、拓扑、metadata 和 feature 都使用同一个构造入口；下标和切片仍返回 `DNA`。`DNASequence`、`DNARecord`、`DNASet` 只作为高级模型和旧代码兼容接口。

## 2. 读取、描述符、模式与热力学 { #reading-analysis }

```python
from io import StringIO

import dnakit
from dnakit.descriptors import gc_at_content, linguistic_complexity
from dnakit.patterns import scan_motif, scan_orfs, scan_restriction_sites
from dnakit.thermodynamics import ThermodynamicConditions, melting_temperature

records = dnakit.read(
    StringIO(">a\nATGGAATTCTAA\n>b\nATGGAATTTTAA\n"),
    format="fasta",
    mode="dna",
)
gc = gc_at_content(records[0])
complexity = linguistic_complexity(records[0], max_word_size=3)
motifs = scan_motif(records[0], "GAATTC", mode="exact")
orfs = scan_orfs(records[0], require_complete=True)
sites = scan_restriction_sites(records[0], ("EcoRI",))
tm = melting_temperature(records[0], conditions=ThermodynamicConditions(sodium_molar=0.05))

print(gc.gc_fraction)
print(complexity.score)
print(len(motifs.hits), len(orfs.hits), len(sites.hits))
print(tm.tm_celsius, tm.algorithm_version)
```

内部复现的 nearest-neighbor 热力学限定为线性、无 Gap、标准 A/C/G/T、完整 Watson–Crick 互补、2–60 nt 和 Na⁺+K⁺ 总单价盐。Mg²+、dNTP、mismatch、修饰与 dangling end 不会被隐式近似。

## 3. 指纹、搜索与聚类 { #fingerprint-search-clustering }

```python
from dnakit.datasets import ClusterConfig, cluster_sequences, deduplicate_approximate
from dnakit.fingerprints import hashed_kmer_fingerprint, minhash, panel_fingerprint
from dnakit.similarity import build_sketch_index, nearest_neighbors

sketches = tuple(minhash(record, k=3, num_hashes=32) for record in records)
hashed = hashed_kmer_fingerprint(records[0], k=3)
panel = panel_fingerprint(records[0], {"start": "ATG", "EcoRI": "GAATTC"})
index = build_sketch_index(records, k=3, num_hashes=32)
neighbors = nearest_neighbors(records[0], index, top_k=2)

config = ClusterConfig(method="identity", threshold=0.8, representative_policy="medoid")
clusters = cluster_sequences(records, config=config)
near_duplicates = deduplicate_approximate(records, config=config)

print(len(sketches), hashed.dimension, panel.dense_values())
print(tuple(hit.record_id for hit in neighbors.hits))
print(clusters.labels, near_duplicates.representatives.ids)
```

阈值聚类是有界、穷举 pairwise 的连通分量语义，不是 CD-HIT/MMseqs2 的别名。大型库应使用经验证的外部索引后端；不要把 `max_records` 提高后直接做无界全矩阵。

## 4. 划分、泄漏与参考库评价 { #split-reference-evaluation }

```python
from dnakit.datasets import LeakageConfig, SplitConfig, detect_leakage, split
from dnakit.evaluation import (
    ReferenceSearchConfig,
    create_reference_library,
    evaluate_memorization,
    evaluate_novelty,
    evaluate_synthesis_risk,
)

partitions = split(
    records,
    config=SplitConfig(method="random", ratios={"train": 0.5, "test": 0.5}, seed=7),
)
# 如果输入顺序可能变化，可改用 method="hash"；它按稳定的 record.id 划分。
split_sets = {subset.name: subset.records for subset in partitions.subsets}
leakage = detect_leakage(
    split_sets,
    config=LeakageConfig(method="identity", threshold=0.9),
)

reference = create_reference_library(
    records[0],
    name="fixed-local-reference",
    version="1",
    source="quickstart fixture",
    filters={"scope": "one record"},
)
search_config = ReferenceSearchConfig(method="identity", copy_threshold=0.9)
novelty = evaluate_novelty(records, reference, config=search_config)
memorization = evaluate_memorization(records, reference, config=search_config)
risk = evaluate_synthesis_risk(records)

print(dict(partitions.counts))
print(len(leakage.events))
print(novelty.metrics["novel_fraction"])
print(memorization.metrics["memorization_fraction"])
print(risk.metrics["disclaimer"])
```

novelty 和 memorization 始终相对于调用方提供的、带名称/版本/摘要/过滤条件的本地参考库定义。synthesis-risk 只是透明序列规则，不是实验合成成功率。

## 5. 可视化 { #visualization }

```python
from pathlib import Path
from tempfile import TemporaryDirectory

from dnakit.visualization import (
    ImageExportConfig,
    plot_sequence,
    save_image,
    save_svg,
)

sequence_svg = plot_sequence(records[0])

with TemporaryDirectory() as directory:
    root = Path(directory)
    svg_result = save_svg(sequence_svg, root / "sequence.svg")
    png_result = save_image(
        sequence_svg,
        root / "sequence.png",
        config=ImageExportConfig(dpi=600, width=640),
    )
    print(svg_result.target_artifact.sha256)
    print(png_result.target_artifact.byte_size, png_result.format)
```

PNG/TIFF/PDF 需要 `viz` extra；SVG 不需要额外图形依赖。

## 6. CLI

```bash
dnakit normalize " acgu "
dnakit describe ACGTACGT
dnakit search ACGTACGT ACG --mode exact
dnakit compare ACGT ACGA --method hamming
```

文件型转换、去重、划分、报告和 YAML/JSON pipeline 见[示例](examples/index.md)。workflow 的统一入口是 `dnakit workflow CONFIG`。所有命令都对覆盖、输入规模和格式错误给出非零退出码与结构化提示。

## 7. 大文件与复现规则

- 普通读取使用 `read(path, mode="dna")`，无论一条还是多条都返回 `DNA`；大文件流式处理使用 `read(path, mode="stream")`。旧的 `read_one()`、`read_set()` 只保留兼容。
- `iter_chunks()` 只保留一个 chunk；未压缩普通 FASTA 和严格四行 FASTQ 可建立带源文件摘要的索引并按 ID/坐标读取。
- `run_batch()` 提供稳定输入顺序、显式 seed、线程模式、错误收集和 `resume_completed_ids`。
- `CacheKey` 自动纳入 key schema 和 DNAKit 版本；输入、参数及实际算法/后端版本仍须由调用方放入 components。`JSONCache` 校验 payload 完整性。
- YAML/JSON workflow 只允许 8 个白名单操作，写入配置中声明的专用目录并生成带 SHA-256 的运行清单；未知步骤/字段会被拒绝。
- 所有随机 API 都要求或记录 seed；外部 backend、参考库和参数版本应随结果/manifest 保存。

更完整的端到端调用见仓库文件 `notebooks/01_advanced_workflow.ipynb`。
