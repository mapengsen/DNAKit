# 注释坐标与基因模型

这些函数只处理显式提供的注释区间，不读取参考基因组或注释数据库。每个函数都返回可序列化结果，并在输出中记录坐标约定，适合放入可复现的预处理流程。

区间处理的方法来源可引用 [pybedtools](https://doi.org/10.1093/bioinformatics/btr539) 的期刊论文；本页实现是无依赖的确定性区间算术，仍需在论文中报告坐标约定、版本和输入数据。

## 1) `convert_annotation_coordinates`

- **Purpose:** 在 `1-based-closed` 与 `0-based-half-open` 之间转换注释区间。
- **Method:** 对每条记录执行确定性的起止坐标变换；不会修改其他字段。
- **API:** `dnakit.annotation.convert_annotation_coordinates(input, *, options=None, verbose=True)`
- **Input:** `{"records": [{"chrom": "chr1", "start": 2, "end": 5}]}`；`options.coordinate_system` 指定输入约定，`target_coordinate_system` 指定输出约定。
- **Output:** 包含 `records`、源/目标坐标约定和 `dnakit.native_annotation_intervals` provenance 的完成信封。
- **Example code:**

```python
from dnakit.annotation import convert_annotation_coordinates

result = convert_annotation_coordinates(
    {"records": [{"chrom": "chr1", "start": 2, "end": 5, "id": "g1"}]},
    options={
        "coordinate_system": "1-based-closed",
        "target_coordinate_system": "0-based-half-open",
    },
    verbose=False,
)
print(result["result"]["records"])
```

- **Example output:**

```json
[{"chrom": "chr1", "start": 1, "end": 5, "id": "g1", "coordinate_system": "0-based-half-open"}]
```

## 2) `make_intergenic_regions`

- **Purpose:** 根据基因区间和染色体长度生成染色体边界内的基因间区。
- **Method:** 先按染色体排序并合并重叠或相邻区间，再计算边界与区间之间的空隙。
- **API:** `dnakit.annotation.make_intergenic_regions(input, *, options=None, verbose=True)`
- **Input:** `input.genes` 为带 `chrom/start/end` 的区间，`input.contig_lengths` 为染色体长度；默认使用 `1-based-closed`。
- **Output:** `result.regions` 列表以及基因数、区域数、染色体长度和合并规则。
- **Example code:**

```python
from dnakit.annotation import make_intergenic_regions

result = make_intergenic_regions(
    {
        "genes": [{"chrom": "chr1", "start": 2, "end": 4}],
        "contig_lengths": {"chr1": 10},
    },
    options={"coordinate_system": "1-based-closed"},
    verbose=False,
)
print(result["result"]["regions"])
```

- **Example output:**

```json
[
  {"id": "intergenic_chr1_1", "chrom": "chr1", "start": 1, "end": 1, "feature_type": "intergenic", "coordinate_system": "1-based-closed"},
  {"id": "intergenic_chr1_2", "chrom": "chr1", "start": 5, "end": 10, "feature_type": "intergenic", "coordinate_system": "1-based-closed"}
]
```

## 3) `make_promoter_regions`

- **Purpose:** 按基因链方向生成启动子窗口。
- **Method:** 正链以 `start` 为转录起始点，负链以 `end` 为转录起始点，并按坐标约定裁剪到染色体边界。
- **API:** `dnakit.annotation.make_promoter_regions(input, *, options=None, verbose=True)`
- **Input:** `input.genes` 至少包含 `chrom/start/end/strand`；`options.upstream` 和 `options.downstream` 为非负窗口长度，可选 `contig_lengths` 用于边界裁剪。
- **Output:** `result.promoters` 列表，包含 `gene_id`、链方向、窗口长度和坐标约定。
- **Example code:**

```python
from dnakit.annotation import make_promoter_regions

result = make_promoter_regions(
    {
        "genes": [{"id": "g1", "chrom": "chr1", "start": 4, "end": 7, "strand": "-"}],
        "contig_lengths": {"chr1": 12},
    },
    options={"coordinate_system": "1-based-closed", "upstream": 2, "downstream": 1},
    verbose=False,
)
print(result["result"]["promoters"])
```

- **Example output:**

```json
[{"id": "promoter_g1", "gene_id": "g1", "chrom": "chr1", "start": 7, "end": 9, "strand": "-", "feature_type": "promoter", "coordinate_system": "1-based-closed", "upstream": 2, "downstream": 1}]
```

## 4) `merge_gene_models`

- **Purpose:** 按基因 ID、染色体和链方向合并显式 exon/segment 区间。
- **Method:** 对同一基因的区间排序，合并重叠区间；`merge_adjacent=True` 时同时合并相邻区间。
- **API:** `dnakit.annotation.merge_gene_models(input, *, options=None, verbose=True)`
- **Input:** `input.segments` 中每条记录包含 `gene_id`（或 `gene`/`id`）、`chrom/start/end` 和可选 `strand`。
- **Output:** `result.models` 列表，包含模型边界、合并后的 `segments` 和区间数量。
- **Example code:**

```python
from dnakit.annotation import merge_gene_models

result = merge_gene_models(
    {"segments": [
        {"gene_id": "g1", "chrom": "chr1", "start": 1, "end": 3},
        {"gene_id": "g1", "chrom": "chr1", "start": 4, "end": 5},
    ]},
    options={"coordinate_system": "1-based-closed", "merge_adjacent": True},
    verbose=False,
)
print(result["result"]["models"])
```

- **Example output:**

```json
[{"gene_id": "g1", "gene_name": "g1", "chrom": "chr1", "strand": "+", "start": 1, "end": 5, "segments": [{"start": 1, "end": 5}], "segment_count": 1, "coordinate_system": "1-based-closed"}]
```
