# G18q 复核：selene-sdk、enformer-pytorch 与 alphagenome 变异评分

本页记录 `predictions` 区域第 G18q 组对台账中 `other` 状态的复核结果。包装函数位于 `src/dnakit/predictions/alphagenome_review_functions.py`，导入路径为 `dnakit.predictions`；这些函数尚未加入 `predictions/__init__.py`，示例因此直接从模块导入。

- 包装的 7 个函数都是 alphagenome 的变异评分工具，输入均为本地数据：protobuf 字节、AnnData 与轨道元数据记录。它们不需要网络或模型权重。
- selene-sdk 与 enformer-pytorch 在本组没有可离线运行的 DNA 级函数，对应的 `selene_sdk_review_functions.py` 与 `enformer_pytorch_review_functions.py` 不导出任何函数。
- 测试位于 `others/tests/unit/test_review_G18q.py`。上游调用与手工计算的结果都在测试中核对。
- 示例中的输出由实际运行得到，文档检查器应报告 0 个不一致。

## 复核结果

| 包 | 复核前 other / cap | 包装 | 复核后的台账状态 |
| --- | --- | --- | --- |
| selene-sdk | 27 / 0 | 0 | internal 22；not-dna 2；duplicate-of 2；needs-data 1 |
| enformer-pytorch | 10 / 0 | 0 | internal 6；not-dna 3；needs-data 1 |
| alphagenome | 8 / 0 | 7 | wrapped 7；internal 1 |

复核后三个台账中已没有 `other` 或 `cap` 行。合计 45 行 `other`：7 个被包装，其余 38 个有明确状态与理由。

- selene-sdk 的 `get_reverse_encoding`（两个模块路径）对 one-hot 输入抛出 IndexError，原因是上游用浮点行向量作为索引。它要完成的反向互补与已有的 `reverse_complement_one_hot` 相同，因此标记为 duplicate-of，不另行包装。
- selene-sdk 的 `predict` 与 enformer-pytorch 的 `get_enformer_embeddings` 都需要训练好的模型，标记为 needs-data，不在本组包装。
- enformer-pytorch 的 Poisson 损失与 Pearson 相关只作用于轨道预测值，不是 DNA 操作，标记为 not-dna。
- 参数冻结、日志、配置读取与运行目录等辅助函数标记为 internal。

## 通用约定

- 每个公开函数的 `verbose` 都是仅关键字参数，默认 `True`。
- 返回值是可 JSON 序列化的字典，包含 `status`、`operation`、结果字段与 `provenance`。
- NaN 输出为 `null`；Variant 与 Interval 等基因组对象输出为字符串，例如 `chr22:36201698:A>C`。
- `alphagenome`、`anndata` 与 `pandas` 只在函数运行时导入；缺少依赖时抛出 `BackendUnavailableError`。
- 输入不合法时抛出 `ConfigurationError`。

## 1) `convert_alphagenome_variant_scores`

- **作用：** 把序列化的 `DenseVariantScores` protobuf 消息转换为每个评分器一份 AnnData 摘要。
- **处理方法：** 逐条解析消息，检查 `track_metadata` 覆盖了消息中出现的全部评分器，然后调用上游 `atlas.convert_variant_scores_to_anndata`。解析错误转为 ConfigurationError。
- **API：** `convert_alphagenome_variant_scores(variant_scores[必须], *, track_metadata[必须], verbose[可选])`
- **输入：** `variant_scores` 是 `DenseVariantScores.SerializeToString()` 得到的 bytes 列表，最多 10000 条；`track_metadata` 是评分器名到轨道记录列表的字典，每条记录是含 `name` 与 `strand` 的映射。
- **输出：** 字典，含 `variant_count`、`scorers` 与 `provenance`。`scorers` 按评分器名排序，每项含 `scorer` 与 `anndata`；`anndata` 含 `shape`、`obs_names`、`var_names`、`obs`、`var`、`scores`（即 X）与 `layers`。

**示例代码：**

```python
import numpy as np
from alphagenome.protos import atlas_service_pb2 as atlas
from alphagenome.protos import dna_model_pb2 as model
from dnakit.predictions import convert_alphagenome_variant_scores

message = atlas.DenseVariantScores()
message.variant.chromosome = "chr22"
message.variant.position = 36201698
message.variant.reference_bases = "A"
message.variant.alternate_bases = "C"
message.interval.chromosome = "chr22"
message.interval.start = 36100000
message.interval.end = 36300000
message.interval.strand = model.STRAND_POSITIVE
score = message.scores.add()
score.variant_scorer.name = "ATAC"
score.shape.extend([1, 2])
score.scores = np.array([0.5, -1.0], dtype=np.float32).tobytes()
tracks = {"ATAC": [{"name": "T1", "strand": "."}, {"name": "T2", "strand": "."}]}
result = convert_alphagenome_variant_scores(
    [message.SerializeToString()], track_metadata=tracks, verbose=False
)
item = result["scorers"][0]
print(item["scorer"], item["anndata"]["shape"], item["anndata"]["scores"])
```

**示例输出：**

```text
ATAC [1, 2] [[0.5, -1.0]]
```

## 2) `tidy_alphagenome_anndata`

- **作用：** 把一个变异评分 AnnData 整理为每行一个分数的长表。
- **处理方法：** 调用上游 `variant_scorers.tidy_anndata`。`match_gene_strand=True` 时删除基因链与轨道链不一致的行；`include_extended_metadata=True` 时保留存在的生物样本、检测方法与本体等列。上游要求非基因评分器的 `obs` 为空表（无列），并要求 `uns` 中有 `interval`、`variant` 与 `variant_scorer`。
- **API：** `tidy_alphagenome_anndata(adata[必须], *, match_gene_strand[可选], include_extended_metadata[可选], verbose[可选])`
- **输入：** `anndata.AnnData` 对象。
- **输出：** 字典，含 `row_count`、`columns`、`records`（每行一个字典）与 `provenance`。

**示例代码：**

```python
import anndata
import numpy as np
import pandas as pd
from alphagenome.data import genome
from alphagenome.models.v1 import variant_scorers
from dnakit.predictions import tidy_alphagenome_anndata

var = pd.DataFrame({"name": ["T1", "T2"], "strand": [".", "."]}, index=["0", "1"])
adata = anndata.AnnData(X=np.array([[0.5, -1.0]], dtype=np.float32), var=var)
adata.uns["interval"] = genome.Interval("chr22", 36100000, 36300000, "+")
adata.uns["variant"] = genome.Variant("chr22", 36201698, "A", "C")
adata.uns["variant_scorer"] = variant_scorers.RECOMMENDED_VARIANT_SCORERS["ATAC"]
result = tidy_alphagenome_anndata(adata, verbose=False)
for row in result["records"]:
    print(row["variant_id"], row["track_name"], row["raw_score"])
```

**示例输出：**

```text
chr22:36201698:A>C T1 0.5
chr22:36201698:A>C T2 -1.0
```

## 3) `tidy_alphagenome_scores`

- **作用：** 把多个变异评分 AnnData 合并为一张长表，并保持输入顺序。
- **处理方法：** 输入列表中每个 AnnData 单独作为一组，嵌套列表中每个子列表作为一组，然后调用上游 `variant_scorers.tidy_scores`。上游返回 None 时输出 0 行。
- **API：** `tidy_alphagenome_scores(scores[必须], *, match_gene_strand[可选], include_extended_metadata[可选], verbose[可选])`
- **输入：** AnnData 列表，或 AnnData 列表的列表，最多 10000 组。
- **输出：** 字典，含 `group_count`、`row_count`、`columns`、`records` 与 `provenance`。

**示例代码：**

```python
import anndata
import numpy as np
import pandas as pd
from alphagenome.data import genome
from alphagenome.models.v1 import variant_scorers
from dnakit.predictions import tidy_alphagenome_scores


def score_table(scorer, values):
    var = pd.DataFrame({"name": ["T1", "T2"], "strand": [".", "."]}, index=["0", "1"])
    adata = anndata.AnnData(X=np.array([values], dtype=np.float32), var=var)
    adata.uns["interval"] = genome.Interval("chr22", 36100000, 36300000, "+")
    adata.uns["variant"] = genome.Variant("chr22", 36201698, "A", "C")
    adata.uns["variant_scorer"] = variant_scorers.RECOMMENDED_VARIANT_SCORERS[scorer]
    return adata


result = tidy_alphagenome_scores(
    [[score_table("ATAC", [0.5, -1.0])], [score_table("DNASE", [2.0, 3.0])]],
    verbose=False,
)
print(result["group_count"], result["row_count"], sorted({row["output_type"] for row in result["records"]}))
```

**示例输出：**

```text
2 4 ['ATAC', 'DNASE']
```

## 4) `merge_alphagenome_track_metadata`

- **作用：** 把 `+` 与 `-` 链的轨道记录合并为一个无链（`.`）轨道。
- **处理方法：** 调用上游 `variant_scoring_utils.merge_stranded_track_metadata`。与合并轨道同名的无链轨道会被删除，避免重复。`validate_tracks=True` 时要求正负轨道的名称顺序一致，否则报错。
- **API：** `merge_alphagenome_track_metadata(track_metadata[必须], *, validate_tracks[可选], verbose[可选])`
- **输入：** 记录列表，每条须含 `name` 与 `strand`，`strand` 只能是 `+`、`-` 或 `.`。
- **输出：** 字典，含 `track_count`、`tracks` 与 `provenance`。

**示例代码：**

```python
from dnakit.predictions import merge_alphagenome_track_metadata

rows = [
    {"name": "A", "strand": "+"},
    {"name": "A", "strand": "-"},
    {"name": "B", "strand": "+"},
    {"name": "B", "strand": "-"},
    {"name": "C", "strand": "."},
]
result = merge_alphagenome_track_metadata(rows, verbose=False)
print(result["track_count"], [(track["name"], track["strand"]) for track in result["tracks"]])
```

**示例输出：**

```text
3 [('A', '.'), ('B', '.'), ('C', '.')]
```

## 5) `merge_alphagenome_stranded_gene_tracks`

- **作用：** 把基因分数从 `+` 与 `-` 两条轨道合并为每个名称一条无链轨道。
- **处理方法：** 调用上游 `variant_scoring_utils.merge_stranded_gene_tracks`。`+` 基因取正链轨道的分数，`-` 基因取负链轨道的分数。若 `obs` 中没有 `gene_id` 与 `strand` 列，上游原样返回输入。
- **API：** `merge_alphagenome_stranded_gene_tracks(scores[必须], *, validate_tracks[可选], verbose[可选])`
- **输入：** AnnData。`obs` 含 `gene_id` 与 `strand`，`var` 含 `name` 与 `strand`，X 为浮点数。
- **输出：** 字典，含 `anndata`（结构与第 1 节的 `anndata` 相同）与 `provenance`。

**示例代码：**

```python
import anndata
import numpy as np
import pandas as pd
from dnakit.predictions import merge_alphagenome_stranded_gene_tracks

obs = pd.DataFrame({"gene_id": ["ENSG1", "ENSG2"], "strand": ["+", "-"]}, index=["0", "1"])
var = pd.DataFrame(
    {"name": ["A", "A", "B", "B"], "strand": ["+", "-", "+", "-"]},
    index=["0", "1", "2", "3"],
)
scores = anndata.AnnData(
    X=np.array([[1, 2, 3, 4], [5, 6, 7, 8]], dtype=float), obs=obs, var=var
)
result = merge_alphagenome_stranded_gene_tracks(scores, verbose=False)
print(result["anndata"]["scores"], result["anndata"]["var_names"])
```

**示例输出：**

```text
[[1.0, 3.0], [6.0, 8.0]] ['0', '1']
```

## 6) `unmerge_alphagenome_stranded_gene_tracks`

- **作用：** 用原始轨道元数据把无链基因分数拆回 `+` 与 `-` 轨道。
- **处理方法：** 调用上游 `variant_scoring_utils.unmerge_stranded_gene_tracks`。每个基因保留与其链一致的轨道；链不一致的列置为 NaN，输出为 `null`。
- **API：** `unmerge_alphagenome_stranded_gene_tracks(scores[必须], *, track_metadata[必须], validate_tracks[可选], verbose[可选])`
- **输入：** AnnData，其 `var` 的 `strand` 须全部为 `.`。`track_metadata` 是含 `name` 与 `strand` 的记录列表，名称须全部出现在 `var` 中。
- **输出：** 字典，含 `anndata` 与 `provenance`。

**示例代码：**

```python
import anndata
import numpy as np
import pandas as pd
from dnakit.predictions import unmerge_alphagenome_stranded_gene_tracks

obs = pd.DataFrame({"gene_id": ["ENSG1", "ENSG2"], "strand": ["+", "-"]}, index=["0", "1"])
merged = anndata.AnnData(
    X=np.array([[1.0, 3.0], [6.0, 8.0]]),
    obs=obs,
    var=pd.DataFrame({"name": ["A", "B"], "strand": [".", "."]}, index=["0", "1"]),
)
tracks = [
    {"name": "A", "strand": "+"},
    {"name": "A", "strand": "-"},
    {"name": "B", "strand": "+"},
    {"name": "B", "strand": "-"},
]
result = unmerge_alphagenome_stranded_gene_tracks(merged, track_metadata=tracks, verbose=False)
print(result["anndata"]["scores"])
```

**示例输出：**

```text
[[1.0, None, 3.0, None], [None, 6.0, None, 8.0]]
```

## 7) `list_alphagenome_recommended_scorers`

- **作用：** 列出某个物种的推荐变异评分器。
- **处理方法：** 调用上游 `variant_scorers.get_recommended_scorers`，再按注册表名称排序输出。
- **API：** `list_alphagenome_recommended_scorers(organism[可选], *, verbose[可选])`
- **输入：** `organism` 为 `homo_sapiens`（默认）或 `mus_musculus`。
- **输出：** 字典，含 `organism`、`scorer_count`、`scorers`（每项含 `name` 与 `scorer` 描述）与 `provenance`。人类 19 个，小鼠 18 个。

**示例代码：**

```python
from dnakit.predictions import list_alphagenome_recommended_scorers

for organism in ["homo_sapiens", "mus_musculus"]:
    result = list_alphagenome_recommended_scorers(organism, verbose=False)
    print(organism, result["scorer_count"], [item["name"] for item in result["scorers"][:3]])
```

**示例输出：**

```text
homo_sapiens 19 ['ATAC', 'ATAC_ACTIVE', 'CAGE']
mus_musculus 18 ['ATAC', 'ATAC_ACTIVE', 'CAGE']
```

## 未包装的函数

- **selene-sdk（27 行，全部未包装）**
  - internal 22：日志初始化 5 个；YAML 读取、配置实例化与模型构建 17 个。其中 `set_random_seeds` 会修改全局 torch 与 numpy 随机状态，因此不包装。
  - not-dna 2：`write_to_hdf5_file` 与 `write_to_tsv_file`，是预测结果的输出写入器。
  - duplicate-of 2：`get_reverse_encoding` 的两个模块路径，对应 `reverse_complement_one_hot`。
  - needs-data 1：`predict`，需要训练好的 selene 模型。
- **enformer-pytorch（10 行，全部未包装）**
  - internal 6：六个参数冻结辅助函数，会修改模块的 `requires_grad`。
  - not-dna 3：两个 `poisson_loss` 与 `pearson_corr_coef`，是对轨道预测值的通用计算。
  - needs-data 1：`get_enformer_embeddings`，需要训练好的 Enformer 模型。
- **alphagenome（1 行未包装）**：`get_all_folds` 返回固定的 `fold0` 至 `fold7` 名称，标记为 internal；区间版本是已包装的 `get_alphagenome_fold_intervals`。
