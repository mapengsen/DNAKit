# scikit-bio 样本元数据与区间元数据

这一组函数位于 `dnakit.microbiome.skbio_p21_functions`，用 [scikit-bio](https://scikit.bio/) 合并并查询样本元数据表，以及查询和连接序列上的区间注释（IntervalMetadata）。依赖为可选的 `dnakit[bio-skibio]`。

## 1) `dnakit.microbiome.skbio_p21_functions.merge_sample_metadata()`

- **作用：** 按共有样本 ID 合并多个样本元数据表。
- **处理方法：** 调用 `skbio.SampleMetadata.merge`：只保留所有表共有的样本，各表的列名不能重复。
- **API：** `merge_sample_metadata(tables)`。
- **输入：** 2 到 20 个映射，每个把样本 ID 映射到 {列: 值}。
- **输出：** `ids`、`columns` 和合并后的 `metadata`。
- **示例代码：**

```python
from dnakit.microbiome.skbio_p21_functions import merge_sample_metadata
a = {"s1": {"x": 1, "c": "u"}, "s2": {"x": 5, "c": "v"}, "s3": {"x": 7, "c": "u"}}
b = {"s1": {"y": 2.5}, "s2": {"y": 3.5}, "s4": {"y": 1.0}}
result = merge_sample_metadata([a, b], verbose=False)
print(result["ids"], result["columns"])
print(result["metadata"]["s2"])
```

- **示例输出：**

```text
['s1', 's2'] ['x', 'c', 'y']
{'x': 5.0, 'c': 'v', 'y': 3.5}
```

## 2) `dnakit.microbiome.skbio_p21_functions.query_sample_ids()`

- **作用：** 用 SQL WHERE 子句筛选样本 ID。
- **处理方法：** 调用 `skbio.SampleMetadata.get_ids`：在内存 SQLite 表上执行 `where`。列名需要时用双引号。省略 `where` 时返回全部 ID。
- **API：** `query_sample_ids(metadata, *, where=None)`。
- **输入：** 样本 ID 到 {列: 值} 的映射，以及可选的 WHERE 子句。
- **输出：** 排序后的 `ids`。
- **示例代码：**

```python
from dnakit.microbiome.skbio_p21_functions import query_sample_ids
metadata = {"s1": {"x": 1, "c": "u"}, "s2": {"x": 5, "c": "v"}, "s3": {"x": 7, "c": "u"}}
print(query_sample_ids(metadata, where='"x" > 2', verbose=False)["ids"])
print(query_sample_ids(metadata, where='"c" = \'u\' AND "x" < 5', verbose=False)["ids"])
print(query_sample_ids(metadata, verbose=False)["ids"])
```

- **示例输出：**

```text
['s2', 's3']
['s1']
['s1', 's2', 's3']
```

## 3) `dnakit.microbiome.skbio_p21_functions.query_interval_metadata()`

- **作用：** 按坐标范围和属性查询序列上的区间特征。
- **处理方法：** 调用 `skbio.metadata.IntervalMetadata.query`：特征只需与任一给定范围重叠，并包含全部给定属性。坐标为 0 起始的半开区间。
- **API：** `query_interval_metadata(interval_metadata, *, bounds=None, attributes=None)`。
- **输入：** `{"upper_bound": int, "intervals": [{"bounds": [[起, 止]], "metadata": {...}}]}`，可选的范围列表和属性映射。
- **输出：** `n_matches` 与 `matches`。
- **示例代码：**

```python
from dnakit.microbiome.skbio_p21_functions import query_interval_metadata
features = {"upper_bound": 100, "intervals": [
    {"bounds": [[0, 20]], "metadata": {"gene": "recA"}},
    {"bounds": [[30, 60]], "metadata": {"gene": "sagB"}},
]}
print(query_interval_metadata(features, bounds=[[10, 35]], verbose=False)["n_matches"])
print(query_interval_metadata(features, attributes={"gene": "sagB"}, verbose=False)["matches"])
```

- **示例输出：**

```text
2
[{'bounds': [[30, 60]], 'metadata': {'gene': 'sagB'}}]
```

## 4) `dnakit.microbiome.skbio_p21_functions.concat_interval_metadata()`

- **作用：** 把多条序列的区间注释连接到同一个坐标空间。
- **处理方法：** 调用 `skbio.metadata.IntervalMetadata.concat`：后面每个对象的坐标平移前面各对象上界之和。
- **API：** `concat_interval_metadata(interval_metadata)`。
- **输入：** 1 到 1000 个与 `query_interval_metadata` 格式相同的对象。
- **输出：** 合并后的 `upper_bound` 与 `intervals`。
- **示例代码：**

```python
from dnakit.microbiome.skbio_p21_functions import concat_interval_metadata
first = {"upper_bound": 3, "intervals": [{"bounds": [[0, 2]], "metadata": {"gene": "sagA"}}]}
second = {"upper_bound": 4, "intervals": [{"bounds": [[1, 4]], "metadata": {"gene": "sagB"}}]}
result = concat_interval_metadata([first, second], verbose=False)
print(result["upper_bound"], result["intervals"])
```

- **示例输出：**

```text
7 [{'bounds': [[0, 2]], 'metadata': {'gene': 'sagA'}}, {'bounds': [[4, 7]], 'metadata': {'gene': 'sagB'}}]
```
