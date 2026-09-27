# clustering

Sequence groups and representative sequences are obtained by clustering sequences based on sequence identity, k-mer, fingerprint similarity, or DNA base model rep.

## 1) Identity clustering {#1-data-007-identity-clustering}

- **Purpose:** Calculate the pairwise similarity of sequences and group them according to the connected components of the threshold graph, and return cluster labels, members and representative records, which are used to organize similar sequences and observe data structures.
- **Method:** Compute pairwise identity, connect records meeting the threshold and take graph components as clusters.


- **API:** `dnakit.datasets.cluster_sequences(records[required], config[optional])`; `config` uses `dnakit.datasets.ClusterConfig`, this item is set to `method="identity"`.

- **Input:** Required record set; optional identity threshold, representation policy, and comparison resource limit.
- **Output:** Returns `ClusteringResult` with `clusters`, `labels`, and representative records.

- **Example code:**

```python
from dnakit import DNARecord, DNASequence
from dnakit.datasets import ClusterConfig, cluster_sequences

records = [
    DNARecord(DNASequence("AAAA"), "a"),
    DNARecord(DNASequence("AAAT"), "b"),
    DNARecord(DNASequence("CCCC"), "c"),
]
result = cluster_sequences(
    records, config=ClusterConfig(method="identity", threshold=0.7)
)
print(result.labels)
```

- **Example output:**

```text
(0, 0, 1)
```

## 2) k-mer clustering {#2-data-008-k-mer-clustering}

- **Purpose:** Establish groups based on k-mer Jaccard, Containment or Cosine similarity between sequences and return cluster members and labels, suitable for constitutive clustering without base-level alignment.
- **Method:** Compute k-mer similarities, build a threshold graph and group its connected components.


- **API:** `dnakit.datasets.cluster_sequences(records[required], config[optional])`; `config` uses `dnakit.datasets.ClusterConfig`, this item is set to `method="kmer"`.

- **Input:** Required record set; optional k, canonical, threshold, representative strategy and pairwise comparison upper limit.
- **Output:** Returns `ClusteringResult` with `clusters`, `labels`, and representative records.

- **Example code:**

```python
from dnakit import DNARecord, DNASequence
from dnakit.datasets import ClusterConfig, cluster_sequences

records = [
    DNARecord(DNASequence("AAAA"), "a"),
    DNARecord(DNASequence("AAAT"), "b"),
    DNARecord(DNASequence("CCCC"), "c"),
]
result = cluster_sequences(
    records, config=ClusterConfig(method="kmer", threshold=0.5, k=2)
)
print(result.labels)
```

- **Example output:**

```text
(0, 0, 1)
```

## 3) Fingerprint clustering {#3-data-009-fingerprint-clustering}

- **Purpose:** Group sequences using Tanimoto/Jaccard similarity of fixed schema DNA fingerprints, returning cluster labels and members for clustering by local patterns or Panel features.
- **Method:** Compare compatible fingerprints, build a threshold graph and group connected components.


- **API:** `dnakit.datasets.cluster_sequences(records[required], config[optional])`; `config` uses `dnakit.datasets.ClusterConfig`, this item is set to `method="fingerprint"`.

- **Input:** Required set of records; optional thresholds, representative policies, and resource caps.
- **Output:** Returns `ClusteringResult` with `clusters`, `labels`, and representative records.

- **Example code:**

```python
from dnakit import DNARecord, DNASequence
from dnakit.datasets import ClusterConfig, cluster_sequences

records = [
    DNARecord(DNASequence("AAAA"), "a"),
    DNARecord(DNASequence("AAAT"), "b"),
    DNARecord(DNASequence("CCCC"), "c"),
]
result = cluster_sequences(
    records, config=ClusterConfig(method="fingerprint", threshold=0.5)
)
print(result.labels, result.method)
```

- **Example output:**

```text
(0, 0, 1) fingerprint
```

## 4) Hierarchical clustering {#4-data-010-hierarchical-clustering}

- **Purpose:** Gradually merge the nearest clusters based on the pre-calculated sequence distance, and return linkage, hierarchical relationships and segmentation labels, which are used to draw dendrograms and observe groupings at different distance scales.
- **Method:** Start from singleton clusters and iteratively merge the closest pair using single/complete/average linkage.


- **API:** `dnakit.datasets.hierarchical_cluster(records[required], config[optional])`; `config` uses `dnakit.datasets.HierarchicalClusteringConfig`.

- **Input:** Required record collection; optional identity/edit/k-mer/fingerprint methods and single/complete/average linkage.
- **Output:** Returns a `HierarchicalClusteringResult`; the example shows `linkage`.

- **Example code:**

```python
from dnakit import DNARecord, DNASequence
from dnakit.datasets import HierarchicalClusteringConfig, hierarchical_cluster

records = [
    DNARecord(DNASequence("AAAA"), "a"),
    DNARecord(DNASequence("AAAT"), "b"),
    DNARecord(DNASequence("CCCC"), "c"),
]
result = hierarchical_cluster(
    records, config=HierarchicalClusteringConfig(linkage="average")
)
print(len(result.linkage), result.linkage[-1].member_count)
```

- **Example output:**

```text
2 3
```

## 5) Representative sequence selection {#5-data-011-representative-sequence-selection}

- **Purpose:** Select a representative sequence from each cluster according to clear strategies such as first, longest, medoid, etc., and retain member mapping for compressed data sets and subsequent manual inspection.
- **Method:** Group records by labels and select one representative using order, length or medoid distances.


- **API:** `dnakit.datasets.select_representatives(records[required], labels[required], policy[optional], medoid_method[optional], k[optional], canonical[optional], max_records[optional], max_pairwise_comparisons[optional], max_alignment_cells[optional])`.

- **Input:** Required records and equal-length labels; optional selection strategy, medoid method, and resource cap.
- **Output:** Returns a `RepresentativeSelectionResult`; the example shows `representative_ids`.

- **Example code:**

```python
from dnakit import DNARecord, DNASequence
from dnakit.datasets import select_representatives

records = [
    DNARecord(DNASequence("AAAA"), "long"),
    DNARecord(DNASequence("AA"), "short"),
    DNARecord(DNASequence("CCCC"), "other"),
]
result = select_representatives(records, [7, 7, 2], policy="shortest")
print(result.representative_ids)
```

- **Example output:**

```text
('other', 'short')
```

## 6) Neural Network Clustering {#data-027-neural-clustering}

- **Purpose:** Extract a fixed-length rep for each sequence using the selected DNA base model, then perform normalization, optional PCA, and k-means with a fixed seed on the vector, returning the label, center, score, and center nearest representative sequence for grouping by model representation.
- **Method:** Pool model vectors, apply L2 normalization and optional PCA, then run seeded k-means++/Lloyd clustering. Select the record nearest each center as its representative.
- **API:** `dnakit.datasets.neural_cluster_sequences(records, config=NeuralClusteringConfig(...))`.
- **Input:** Nonempty DNA records without explicit gaps, the cluster count, seed and model configuration. Default LucaOne requires its dependencies and `allow_remote_code=True`; IUPAC symbols map to `N` by default.
- **Output:** Returns a `NeuralClusteringResult`; the example shows `labels`, `representatives`, `checkpoint_path`.
- **Example code:**

```python
from dnakit import DNARecord, DNASequence
from dnakit.datasets import NeuralClusteringConfig, neural_cluster_sequences
from dnakit.representations import RepresentationConfig

records = [
    DNARecord(DNASequence("AAAAAA"), "a"),
    DNARecord(DNASequence("AAAAAT"), "b"),
    DNARecord(DNASequence("CCCCCC"), "c"),
    DNARecord(DNASequence("CCCCCG"), "d"),
]
result = neural_cluster_sequences(
    records,
    config=NeuralClusteringConfig(
        representation=RepresentationConfig(
            batch_size=2,
            allow_remote_code=True,
        ),
        n_clusters=2,
        seed=7,
    ),
)
print(result.labels)
print(result.representatives.ids)
print(result.checkpoint_path)
```

- **Example output:** Four cluster labels, two representative IDs and the checkpoint path. The exact labels depend on the model and environment; model inference was not run here.

<span id="checkpoint"></span>[FAQ: details](../../faq/tutorial_basics.md#neural-clustering)

## 7) Cluster samples {#bio-25}

- **Purpose:** Group similar samples by their features.
- **Method:** Apply the selected clustering method; MCL iterates expansion/inflation on a nonnegative affinity matrix.
- **API:** `dnakit.bioinformatics.features.cluster_features(values, method='kmeans', seed=0, parameters=None)`
- **Input:** A numeric feature matrix; scale it beforehand. MCL requires symmetric nonnegative affinity; requires `bio-analysis`; see [available methods](../../faq/bioinformatics.md#encoding).
- **Output:** Returns a clustering model; `labels_` gives the cluster label for each sample.
- **Example code:**

```python
from dnakit.bioinformatics.features import cluster_features

model = cluster_features([[0, 0], [0, 1], [9, 9], [9, 10]],
                         parameters={"n_clusters": 2}, seed=0)
print(model.labels_[0] == model.labels_[1])
print(model.labels_[0] != model.labels_[2])
```

- **Example output:**

```text
True
True
```

Dependencies and details: [FAQ](../../faq/bioinformatics.md#encoding).

## 8) Label-stratified diversity sampling {#diversity-sampling}

`dnakit.datasets.diversity_sample` selects a fixed number of original candidate records.
It defaults to frozen DNABERT2 mean embeddings, L2 normalization, and seeded
k-means++/Lloyd within each class. The model runs inference only; no training is performed.
The `diversity` argument selects high or low diversity:

| Argument | Cluster count | Selection |
| --- | --- | --- |
| `diversity="high"` (default) | `n_samples` in total, allocated proportionally | One nearest actual record per centroid |
| `diversity="low"` | Five per class with a nonzero sample quota by default | The tightest fixed-size neighborhood among the centers |

`n_samples` always sets the **final sample size**; it equals the total cluster count
only in `high` mode. `low_n_clusters` sets the per-class cluster count for `low` mode.
Explicit `n_samples`, `diversity`, and `low_n_clusters` override the corresponding
configuration values without modifying the configuration object. Omitted arguments
use the configuration values, defaulting to 400, `"high"`, and 5 respectively.

Class quotas follow the candidate proportions using the **largest-remainder method**:
floor each proportional quota, then assign remaining slots in descending remainder order.
For 840 candidates with class counts 173 and 667, selecting 400 allocates 82 and 318
records respectively. High mode fits 82 and 318 clusters; low mode fits five per class,
then selects 82 and 318 records. Quota and distance ties follow input order.

Install model dependencies with `python -m pip install "dnakit[neural]"`.
This example expects `sequence` and `label` columns in `candidates.csv`:

```python
from dnakit.datasets import DiversitySamplingConfig, diversity_sample
from dnakit.representations import RepresentationConfig

result = diversity_sample(
    "candidates.csv",
    n_samples=400,
    diversity="high",
    output_path="diversity_400",
    config=DiversitySamplingConfig(
        label_key="label",
        seed=42,
        representation=RepresentationConfig(
            model="dnabert2",
            pooling="mean",
            allow_remote_code=True,  # Load model code from the official checkpoint
            # checkpoint_path="/path/to/dnabert2-117m",  # Reuse local weights
        ),
    ),
)
print(len(result.records))
for allocation in result.allocations:
    print(allocation.label, allocation.input_count, allocation.n_samples)
```

To select a low-diversity subset:

```python
low_result = diversity_sample(
    "candidates.csv",
    n_samples=400,
    diversity="low",
    low_n_clusters=5,
    output_path="low_diversity_400",
    config=DiversitySamplingConfig(
        seed=42,
        representation=RepresentationConfig(model="dnabert2", allow_remote_code=True),
    ),
)
```

Low mode clusters each class separately. Each center proposes the class quota of
nearest neighbors from **all candidates in that class**, yielding up to five proposals.
Score each proposal by its mean Euclidean distance to **its own mean embedding** and
choose the smallest score, breaking score ties by sorted original row indices.
Neighborhoods may cross cluster boundaries, so a small cluster does not cause an
underfilled sample. This is a compactness heuristic, not a global optimum.

File input produces `diversity_400/diversity_sample.csv` and `diversity_400/result.json`.
The CSV retains original columns and selected rows in input order; the JSON records
allocations, source indices, cluster membership, centers, and model provenance.
Alternatively pass `Iterable[DNARecord]` with labels in `record.metadata["label"]`.
`result.records` contains the original records, preserving sequences, labels, annotations,
and metadata.

- `selected_indices` uses zero-based positions in the input; `cluster_labels` maps each
  input record to its cluster. In low mode, use `records` and `selected_indices` for
  the selected subset; cluster representatives alone do not constitute the sample.
  Results record `diversity`, `low_n_clusters`, and `sampling_strategy`; inspect
  `len(result.clusters)` for the actual total cluster count.
- `label_key=None` samples one unlabeled pool. Otherwise every record requires a non-null
  string, numeric, or boolean label. Types are distinct: `1` and `"1"` are separate classes.
- Very small classes can receive zero slots; their `cluster_labels` entries are `None`.
- Set `normalize=False` to cluster raw embeddings using Euclidean distance. No PCA is
  performed. Other DNAKit models can be selected with `RepresentationConfig(model=...)`.
- `n_samples` must be positive and cannot exceed the input count. High mode rejects
  insufficient distinct embeddings. Low mode caps each class's cluster count at its
  number of distinct embeddings and permits distinct rows with identical embeddings.
  Neither mode repeats an input row. Empty fitted clusters raise an error.
  Both sampling and representation `max_records` limits apply, defaulting to 10,000.
- This method targets representative or compact subsets in model feature space. It does
  not guarantee maximum or minimum pairwise sequence distance or deduplicate across labels. Reproducibility
  also depends on input order, model weights, library versions, and the compute environment.

Tests cover both modes with synthetic 840-to-400 embeddings, 82/318 quotas, compact
region selection, duplicate embeddings, and CSV export. They do not establish biological diversity gains
or constitute DNABERT2 inference on a user's candidate pool.
References: [official DNABERT-2 implementation](https://github.com/MAGICS-LAB/DNABERT_2)
and [scikit-learn KMeans](https://scikit-learn.org/stable/modules/generated/sklearn.cluster.KMeans.html).
