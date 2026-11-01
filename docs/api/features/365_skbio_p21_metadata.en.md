# scikit-bio sample metadata and interval metadata

These functions in `dnakit.microbiome.skbio_p21_functions` use [scikit-bio](https://scikit.bio/) to merge and query sample-metadata tables and to query and concatenate interval annotations (IntervalMetadata) on sequences. The dependency is the optional `dnakit[bio-skibio]`.

## 1) `dnakit.microbiome.skbio_p21_functions.merge_sample_metadata()`

- **Purpose:** Merge several sample-metadata tables on their shared sample ids.
- **Method:** Calls `skbio.SampleMetadata.merge`: only samples present in every table are kept and column names must differ between tables.
- **API:** `merge_sample_metadata(tables)`.
- **Input:** 2 to 20 mappings, each mapping sample id to {column: value}.
- **Output:** `ids`, `columns` and the merged `metadata`.
- **Example code:**

```python
from dnakit.microbiome.skbio_p21_functions import merge_sample_metadata
a = {"s1": {"x": 1, "c": "u"}, "s2": {"x": 5, "c": "v"}, "s3": {"x": 7, "c": "u"}}
b = {"s1": {"y": 2.5}, "s2": {"y": 3.5}, "s4": {"y": 1.0}}
result = merge_sample_metadata([a, b], verbose=False)
print(result["ids"], result["columns"])
print(result["metadata"]["s2"])
```

- **Example output:**

```text
['s1', 's2'] ['x', 'c', 'y']
{'x': 5.0, 'c': 'v', 'y': 3.5}
```

## 2) `dnakit.microbiome.skbio_p21_functions.query_sample_ids()`

- **Purpose:** Select sample ids with an SQL WHERE clause.
- **Method:** Calls `skbio.SampleMetadata.get_ids`: `where` runs on an in-memory SQLite table. Quote column names with double quotes when needed. Without `where` all ids are returned.
- **API:** `query_sample_ids(metadata, *, where=None)`.
- **Input:** A mapping of sample id to {column: value} and an optional WHERE clause.
- **Output:** The sorted `ids`.
- **Example code:**

```python
from dnakit.microbiome.skbio_p21_functions import query_sample_ids
metadata = {"s1": {"x": 1, "c": "u"}, "s2": {"x": 5, "c": "v"}, "s3": {"x": 7, "c": "u"}}
print(query_sample_ids(metadata, where='"x" > 2', verbose=False)["ids"])
print(query_sample_ids(metadata, where='"c" = \'u\' AND "x" < 5', verbose=False)["ids"])
print(query_sample_ids(metadata, verbose=False)["ids"])
```

- **Example output:**

```text
['s2', 's3']
['s1']
['s1', 's2', 's3']
```

## 3) `dnakit.microbiome.skbio_p21_functions.query_interval_metadata()`

- **Purpose:** Query interval features of a sequence by bounds and attributes.
- **Method:** Calls `skbio.metadata.IntervalMetadata.query`: a feature must overlap any given span and carry every given attribute. Coordinates are 0-based half-open.
- **API:** `query_interval_metadata(interval_metadata, *, bounds=None, attributes=None)`.
- **Input:** `{"upper_bound": int, "intervals": [{"bounds": [[start, end]], "metadata": {...}}]}`, an optional list of spans and an attribute mapping.
- **Output:** `n_matches` and `matches`.
- **Example code:**

```python
from dnakit.microbiome.skbio_p21_functions import query_interval_metadata
features = {"upper_bound": 100, "intervals": [
    {"bounds": [[0, 20]], "metadata": {"gene": "recA"}},
    {"bounds": [[30, 60]], "metadata": {"gene": "sagB"}},
]}
print(query_interval_metadata(features, bounds=[[10, 35]], verbose=False)["n_matches"])
print(query_interval_metadata(features, attributes={"gene": "sagB"}, verbose=False)["matches"])
```

- **Example output:**

```text
2
[{'bounds': [[30, 60]], 'metadata': {'gene': 'sagB'}}]
```

## 4) `dnakit.microbiome.skbio_p21_functions.concat_interval_metadata()`

- **Purpose:** Concatenate interval annotations of several sequences into one coordinate space.
- **Method:** Calls `skbio.metadata.IntervalMetadata.concat`: the coordinates of each later object are shifted by the summed upper bounds of the earlier ones.
- **API:** `concat_interval_metadata(interval_metadata)`.
- **Input:** 1 to 1000 objects in the format of `query_interval_metadata`.
- **Output:** The combined `upper_bound` and `intervals`.
- **Example code:**

```python
from dnakit.microbiome.skbio_p21_functions import concat_interval_metadata
first = {"upper_bound": 3, "intervals": [{"bounds": [[0, 2]], "metadata": {"gene": "sagA"}}]}
second = {"upper_bound": 4, "intervals": [{"bounds": [[1, 4]], "metadata": {"gene": "sagB"}}]}
result = concat_interval_metadata([first, second], verbose=False)
print(result["upper_bound"], result["intervals"])
```

- **Example output:**

```text
7 [{'bounds': [[0, 2]], 'metadata': {'gene': 'sagA'}}, {'bounds': [[4, 7]], 'metadata': {'gene': 'sagB'}}]
```
