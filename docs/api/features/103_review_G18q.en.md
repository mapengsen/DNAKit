# G18q review: selene-sdk, enformer-pytorch and alphagenome variant scoring

This page records the review of the `other` ledger rows in review group G18q of the `predictions` area. The wrappers live in `src/dnakit/predictions/alphagenome_review_functions.py` and import from `dnakit.predictions`. They are not yet exported from `predictions/__init__.py`, so the examples import them from the module directly.

- The 7 wrapped functions are alphagenome variant-scoring utilities. Their inputs are all local: protobuf bytes, AnnData objects and track-metadata records. None needs the network or model weights.
- selene-sdk and enformer-pytorch have no offline DNA-level function in this group. Their modules, `selene_sdk_review_functions.py` and `enformer_pytorch_review_functions.py`, export nothing.
- The tests are in `others/tests/unit/test_review_G18q.py`. Each result is checked against either the upstream call or a hand calculation.
- The example outputs come from running the code, so the documentation checker should report 0 mismatches.

## Review results

| Package | Other / cap before review | Wrapped | Ledger status after review |
| --- | --- | --- | --- |
| selene-sdk | 27 / 0 | 0 | internal 22; not-dna 2; duplicate-of 2; needs-data 1 |
| enformer-pytorch | 10 / 0 | 0 | internal 6; not-dna 3; needs-data 1 |
| alphagenome | 8 / 0 | 7 | wrapped 7; internal 1 |

After the review, none of the three ledgers has an `other` or `cap` row. Of the 45 `other` rows, 7 are wrapped and the other 38 have an explicit status and reason.

- selene-sdk `get_reverse_encoding` (both module paths) raises IndexError on one-hot input, because upstream uses float row values as indices. The reverse complement it is meant to compute is the same as `reverse_complement_one_hot`, so it is marked duplicate-of and not wrapped again.
- selene-sdk `predict` and enformer-pytorch `get_enformer_embeddings` both need a trained model. They are marked needs-data and are not wrapped in this group.
- enformer-pytorch Poisson loss and Pearson correlation act only on track predictions. They are not DNA operations and are marked not-dna.
- Parameter-freezing, logging, configuration-loading and run-directory helpers are marked internal.

## Shared conventions

- The `verbose` argument of every public function is keyword-only and defaults to `True`.
- Return values are JSON-serializable dicts with `status`, `operation`, the result fields and `provenance`.
- NaN is returned as `null`. Genomic objects such as Variant and Interval are returned as strings, for example `chr22:36201698:A>C`.
- `alphagenome`, `anndata` and `pandas` are imported only when a function runs. A missing dependency raises `BackendUnavailableError`.
- Invalid input raises `ConfigurationError`.

## 1) `convert_alphagenome_variant_scores`

- **Purpose:** Convert serialized `DenseVariantScores` protobuf messages into one AnnData summary per variant scorer.
- **Method:** Each message is parsed, the scorer names used in the messages are checked against `track_metadata`, and the upstream `atlas.convert_variant_scores_to_anndata` is called. Parse errors are raised as ConfigurationError.
- **API:** `convert_alphagenome_variant_scores(variant_scores[required], *, track_metadata[required], verbose[optional])`
- **Input:** `variant_scores` is a list of bytes from `DenseVariantScores.SerializeToString()`, at most 10000 entries. `track_metadata` maps each scorer name to a list of track records, where each record is a mapping with `name` and `strand`.
- **Output:** A dict with `variant_count`, `scorers` and `provenance`. `scorers` is sorted by scorer name, and each item has `scorer` and `anndata`. `anndata` contains `shape`, `obs_names`, `var_names`, `obs`, `var`, `scores` (the X matrix) and `layers`.

**Example code:**

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

**Example output:**

```text
ATAC [1, 2] [[0.5, -1.0]]
```

## 2) `tidy_alphagenome_anndata`

- **Purpose:** Turn one variant-score AnnData into a long table with one score per row.
- **Method:** Calls the upstream `variant_scorers.tidy_anndata`. With `match_gene_strand=True`, rows whose gene strand and track strand disagree are removed. With `include_extended_metadata=True`, biosample, assay and ontology columns are kept when present. The upstream code requires an empty `obs` table (no columns) for non-gene scorers, and requires `interval`, `variant` and `variant_scorer` in `uns`.
- **API:** `tidy_alphagenome_anndata(adata[required], *, match_gene_strand[optional], include_extended_metadata[optional], verbose[optional])`
- **Input:** An `anndata.AnnData` object.
- **Output:** A dict with `row_count`, `columns`, `records` (one dict per row) and `provenance`.

**Example code:**

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

**Example output:**

```text
chr22:36201698:A>C T1 0.5
chr22:36201698:A>C T2 -1.0
```

## 3) `tidy_alphagenome_scores`

- **Purpose:** Combine several variant-score AnnData objects into one long table, keeping the input order.
- **Method:** In a flat input list, each AnnData is one group. In a nested list, each sub-list is one group. The groups are then passed to the upstream `variant_scorers.tidy_scores`. If upstream returns None, the output has 0 rows.
- **API:** `tidy_alphagenome_scores(scores[required], *, match_gene_strand[optional], include_extended_metadata[optional], verbose[optional])`
- **Input:** A list of AnnData objects, or a list of lists of AnnData objects, at most 10000 groups.
- **Output:** A dict with `group_count`, `row_count`, `columns`, `records` and `provenance`.

**Example code:**

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

**Example output:**

```text
2 4 ['ATAC', 'DNASE']
```

## 4) `merge_alphagenome_track_metadata`

- **Purpose:** Merge the `+` and `-` strand track records into one unstranded (`.`) track per name.
- **Method:** Calls the upstream `variant_scoring_utils.merge_stranded_track_metadata`. An unstranded track with the same name as a merged track is dropped to avoid duplication. With `validate_tracks=True`, the positive and negative tracks must have the same ordered names, or an error is raised.
- **API:** `merge_alphagenome_track_metadata(track_metadata[required], *, validate_tracks[optional], verbose[optional])`
- **Input:** A list of records. Each record needs `name` and `strand`, where `strand` is `+`, `-` or `.`.
- **Output:** A dict with `track_count`, `tracks` and `provenance`.

**Example code:**

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

**Example output:**

```text
3 [('A', '.'), ('B', '.'), ('C', '.')]
```

## 5) `merge_alphagenome_stranded_gene_tracks`

- **Purpose:** Merge gene scores from the `+` and `-` tracks into one unstranded track per name.
- **Method:** Calls the upstream `variant_scoring_utils.merge_stranded_gene_tracks`. A `+` gene takes the positive-track scores, and a `-` gene takes the negative-track scores. If `obs` has no `gene_id` and `strand` columns, the upstream returns the input unchanged.
- **API:** `merge_alphagenome_stranded_gene_tracks(scores[required], *, validate_tracks[optional], verbose[optional])`
- **Input:** An AnnData object. `obs` has `gene_id` and `strand`, `var` has `name` and `strand`, and X is floating point.
- **Output:** A dict with `anndata` (same structure as the `anndata` item in section 1) and `provenance`.

**Example code:**

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

**Example output:**

```text
[[1.0, 3.0], [6.0, 8.0]] ['0', '1']
```

## 6) `unmerge_alphagenome_stranded_gene_tracks`

- **Purpose:** Split unstranded gene scores back into `+` and `-` tracks, using the original track metadata.
- **Method:** Calls the upstream `variant_scoring_utils.unmerge_stranded_gene_tracks`. Each gene keeps the track that matches its strand. Columns on the mismatched strand are set to NaN, which is output as `null`.
- **API:** `unmerge_alphagenome_stranded_gene_tracks(scores[required], *, track_metadata[required], validate_tracks[optional], verbose[optional])`
- **Input:** An AnnData object whose `var` `strand` values are all `.`. `track_metadata` is a list of records with `name` and `strand`, and every name must appear in `var`.
- **Output:** A dict with `anndata` and `provenance`.

**Example code:**

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

**Example output:**

```text
[[1.0, None, 3.0, None], [None, 6.0, None, 8.0]]
```

## 7) `list_alphagenome_recommended_scorers`

- **Purpose:** List the recommended variant scorers for one organism.
- **Method:** Calls the upstream `variant_scorers.get_recommended_scorers`, then outputs the results sorted by registry name.
- **API:** `list_alphagenome_recommended_scorers(organism[optional], *, verbose[optional])`
- **Input:** `organism` is `homo_sapiens` (the default) or `mus_musculus`.
- **Output:** A dict with `organism`, `scorer_count`, `scorers` (each item has `name` and a `scorer` description) and `provenance`. There are 19 scorers for human and 18 for mouse.

**Example code:**

```python
from dnakit.predictions import list_alphagenome_recommended_scorers

for organism in ["homo_sapiens", "mus_musculus"]:
    result = list_alphagenome_recommended_scorers(organism, verbose=False)
    print(organism, result["scorer_count"], [item["name"] for item in result["scorers"][:3]])
```

**Example output:**

```text
homo_sapiens 19 ['ATAC', 'ATAC_ACTIVE', 'CAGE']
mus_musculus 18 ['ATAC', 'ATAC_ACTIVE', 'CAGE']
```

## Functions not wrapped

- **selene-sdk (27 rows, none wrapped)**
  - internal 22: 5 logging setup functions; 17 YAML loading, configuration instantiation and model-construction helpers. `set_random_seeds` changes the global torch and numpy random state, so it is not wrapped.
  - not-dna 2: `write_to_hdf5_file` and `write_to_tsv_file`, the output writers for prediction results.
  - duplicate-of 2: the two module paths of `get_reverse_encoding`, which correspond to `reverse_complement_one_hot`.
  - needs-data 1: `predict`, which needs a trained selene model.
- **enformer-pytorch (10 rows, none wrapped)**
  - internal 6: six parameter-freezing helpers that change module `requires_grad`.
  - not-dna 3: the two `poisson_loss` functions and `pearson_corr_coef`, generic computations on track predictions.
  - needs-data 1: `get_enformer_embeddings`, which needs a trained Enformer model.
- **alphagenome (1 row, not wrapped)**: `get_all_folds` returns the fixed names `fold0` to `fold7` and is marked internal. The interval version, `get_alphagenome_fold_intervals`, is already wrapped.
