# clustering

Sequence groups and representative sequences are obtained by clustering sequences based on sequence identity, k-mer, fingerprint similarity, or DNA base model rep.

## 1) DATA-007 · Identity clustering

**Function:** Calculate the pairwise similarity of sequences and group them according to the connected components of the threshold graph, and return cluster labels, members and representative records, which are used to organize similar sequences and observe data structures.


**API:** `dnakit.datasets.cluster_sequences(records[required], config[optional])`; `config` uses `dnakit.datasets.ClusterConfig`, this item is set to `method="identity"`.

**Input:** Required record set; optional identity threshold, representation policy, and comparison resource limit.

**Sample code:**

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

**Example results:**

```text
(0, 0, 1)
```

## 2) DATA-008 · k-mer clustering

**Function:** Establish groups based on k-mer Jaccard, Containment or Cosine similarity between sequences and return cluster members and labels, suitable for constitutive clustering without base-level alignment.


**API:** `dnakit.datasets.cluster_sequences(records[required], config[optional])`; `config` uses `dnakit.datasets.ClusterConfig`, this item is set to `method="kmer"`.

**Input:** Required record set; optional k, canonical, threshold, representative strategy and pairwise comparison upper limit.

**Sample code:**

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

**Example results:**

```text
(0, 0, 1)
```

## 3) DATA-009 · Fingerprint clustering

**Function:** Group sequences using Tanimoto/Jaccard similarity of fixed schema DNA fingerprints, returning cluster labels and members for clustering by local patterns or Panel features.


**API:** `dnakit.datasets.cluster_sequences(records[required], config[optional])`; `config` uses `dnakit.datasets.ClusterConfig`, this item is set to `method="fingerprint"`.

**Input:** Required set of records; optional thresholds, representative policies, and resource caps.

**Sample code:**

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

**Example results:**

```text
(0, 0, 1) fingerprint
```

## 4) DATA-010 · Hierarchical clustering

**Function:** Gradually merge the nearest clusters based on the pre-calculated sequence distance, and return linkage, hierarchical relationships and segmentation labels, which are used to draw dendrograms and observe groupings at different distance scales.


**API:** `dnakit.datasets.hierarchical_cluster(records[required], config[optional])`; `config` uses `dnakit.datasets.HierarchicalClusteringConfig`.

**Input:** Required record collection; optional identity/edit/k-mer/fingerprint methods and single/complete/average linkage.

**Sample code:**

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

**Example results:**

```text
2 3
```

## 5) DATA-011 · Representative sequence selection

**Function:** Select a representative sequence from each cluster according to clear strategies such as first, longest, medoid, etc., and retain member mapping for compressed data sets and subsequent manual inspection.


**API:** `dnakit.datasets.select_representatives(records[required], labels[required], policy[optional], medoid_method[optional], k[optional], canonical[optional], max_records[optional], max_pairwise_comparisons[optional], max_alignment_cells[optional])`.

**Input:** Required records and equal-length labels; optional selection strategy, medoid method, and resource cap.

**Sample code:**

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

**Example results:**

```text
('other', 'short')
```

## 6) DATA-027 · Neural Network Clustering {#data-027-neural-clustering}

**Function:** Extract a fixed-length rep for each sequence using the selected DNA base model, then perform normalization, optional PCA, and k-means with a fixed seed on the vector, returning the label, center, score, and center nearest representative sequence for grouping by model representation.

What is used here is **k-means**, not k-mer similarity clustering. If you need not load the neural network
k-mer clustering, use `DATA-008` on this page.

**API:**

- `dnakit.representations.extract_representations(records[required], config[optional])`: Extract only rep;
- `dnakit.datasets.neural_cluster_sequences(records[required], config[optional])`: Execute k-means after extracting rep;
- Use `dnakit.representations.RepresentationConfig` for characterization configuration; use `dnakit.representations.RepresentationConfig` for clustering configuration
  `dnakit.datasets.NeuralClusteringConfig`.

Standalone rep API, all 11 models and their checkpoint/dependency boundaries See also
[The neural network characterization column ](08_fingerprints.md#neural-representations) in the sequence characterization page.

**Input:** Required `DNARecord` collection; optional model, checkpoint position, mean/cls/max/last
Pooling, device, dtype, batch size, L2 normalization, PCA dimensionality, number of clusters, and random seeds. explicit
`Gap` is not supported; the default mapping of IUPAC fuzzy bases is `N`, and can also be set to report an error directly.

**Clustering process:** Each sequence is divided into blocks according to the model context window. The pooling vector of each block is first obtained, and then the blocks are
Vector averaging; then L2 normalization is performed by default, PCA is optional, and finally a fixed seed is performed
`k-means++`/Lloyd clustering. Each group of representative sequences takes the input record closest to the cluster center; the center is located
Space after normalization and optional PCA.

**Sample code:**

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

The default model is `lucaone`, corresponding to checkpoint
`LucaGroup/LucaOne-gene-step36.8M`. The first time you call it, it will be downloaded by default to the directory where the command is run.
`ckpt/lucaone-gene-step36-8m/`; If the complete checkpoint already exists, it will be directly reused without repeated downloading.
Each directory will hold a list of `.dnakit-checkpoint.json` sources. Available `checkpoint_dir` to modify the cache
root directory, or use `checkpoint_path` to point to an existing checkpoint directory and skip downloading entirely. Download and
Sequence-by-sequence extraction displays a progress bar by default.

LucaOne checkpoint contains custom Transformers code. DNAKit does not do this implicitly; use
Default model's standard backend must still be set explicitly in `RepresentationConfig`
`allow_remote_code=True`.

Extract only rep:

```python
from dnakit.representations import RepresentationConfig, extract_representations

reps = extract_representations(
    records,
    config=RepresentationConfig(allow_remote_code=True),
)
print(reps.representations.shape)
```

<span id="checkpoint"></span>**Optional models and official checkpoint**

| `model` | Default official checkpoint | Loading conditions |
| --- | --- | --- |
| `grover` | [PoetschLab/GROVER](https://huggingface.co/PoetschLab/GROVER) | `neural` extra; the local independent model environment has completed the rep + k-means smoke of the real checkpoint |
| `dnabert2` | [zhihan1996/DNABERT-2-117M](https://huggingface.co/zhihan1996/DNABERT-2-117M) | `neural` extra; requires explicit `allow_remote_code=True` |
| `ntv2` | [InstaDeepAI/nucleotide-transformer-v2-500m-multi-species](https://huggingface.co/InstaDeepAI/nucleotide-transformer-v2-500m-multi-species) | `neural` extra; requires explicit `allow_remote_code=True` |
| `hyenadna` | [LongSafari/hyenadna-medium-450k-seqlen-hf](https://huggingface.co/LongSafari/hyenadna-medium-450k-seqlen-hf) | `neural` extra; requires explicit `allow_remote_code=True` |
| `caduceus` | [kuleshov-group/caduceus-ph 131k](https://huggingface.co/kuleshov-group/caduceus-ph_seqlen-131k_d_model-256_n_layer-16) | `neural`, `neural-caduceus` extras; need to explicitly `allow_remote_code=True` |
| `lucaone` (default) | [LucaGroup/LucaOne-gene-step36.8M](https://huggingface.co/LucaGroup/LucaOne-gene-step36.8M) | `neural` extra; requires explicit `allow_remote_code=True` |
| `generator` | [GenerTeam/GENERator-v2-eukaryote-1.2b-base](https://huggingface.co/GenerTeam/GENERator-v2-eukaryote-1.2b-base) | `neural` extra; requires explicit `allow_remote_code=True` |
| `enformer` | [EleutherAI/enformer-official-rough](https://huggingface.co/EleutherAI/enformer-official-rough) | `neural`, `neural-enformer` extras |
| `alphagenome` | [google/alphagenome-all-folds](https://huggingface.co/google/alphagenome-all-folds) | Must accept the model terms and install [official research code](https://github.com/google-deepmind/alphagenome_research); high hardware and JAX/Orbax environment requirements |
| `janusdna` | [Harvard Dataverse DOI 10.7910/DVN/HDT0RN](https://doi.org/10.7910/DVN/HDT0RN) | checkpoint automatically verifies the official MD5; you must also obtain the [official source code](https://github.com/Qihao-Duan/JanusDNA) and set up `model_source_path` |
| `evo2` | [arcinstitute/evo2_7b](https://huggingface.co/arcinstitute/evo2_7b) | `neural`, `neural-evo2` extras; follow [ official Evo 2 environment ](https://github.com/ArcInstitute/evo2) and GPU requirements |

Except for GROVER's real smoke, the other models have completed registration, checkpoint analysis and correspondence in this function.
adapter is connected, but full checkpoint value verification has not been completed one by one on the same machine. Therefore
`DATA-027` The status is `conditional`, "Connected" cannot be interpreted as meaning that all models, hardware and versions
Verified.

**Remote code boundaries:** Models noted in the table load Python code from checkpoint. DNAKit
This type of execution is denied by default and can only be set explicitly after the caller reviews the official repository and checkpoint.
`allow_remote_code=True` is loaded. checkpoint is not distributed with wheel/sdist; model license,
Access terms, video memory and storage requirements are separately confirmed by the user.

**Result:** `NeuralClusteringResult` Save labels, cluster members, center nearest representative sequence,
Model name, checkpoint path, original/clustering dimension, inertia, calculated silhouette, PCA
Explained variance, seed, and actual number of iterations.
